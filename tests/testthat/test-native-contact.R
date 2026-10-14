native_contact_pairs <- function() {
  data.frame(
    variant_id = c("observed", "zero", "diagonal"),
    gene_id = c("gene_observed", "gene_zero", "gene_diagonal"),
    chrom = c("1", "1", "1"),
    pos1 = c(2500001, 1, 1),
    tss_chrom = c("1", "2", "1"),
    tss_pos1 = c(5000001, 2500001, 1),
    stringsAsFactors = FALSE
  )
}

test_that("native query planner bounds and deduplicates exact pixel tiles", {
  query <- data.frame(
    chrom1 = c("chr1", "chr1", "chr2", "chr1"),
    start1_0 = c(0, 10, 20, 80),
    chrom2 = c("chr1", "chr1", "chr1", "chr1"),
    start2_0 = c(10, 80, 0, 10),
    stringsAsFactors = FALSE
  )
  query$.pixel_key <- ACV2G:::.contact_pixel_key(
    query$chrom1, query$start1_0, query$chrom2, query$start2_0
  )
  planned <- ACV2G:::.contact_native_query_plan(
    query,
    resolution = 10,
    available_chromosomes = c("chr1", "chr2"),
    chrom_lengths = c(chr1 = 95, chr2 = 50),
    tile_n_bins = 2
  )

  expect_lte(max(planned$plan$tile1_end0 - planned$plan$tile1_start0), 20)
  expect_lte(max(planned$plan$tile2_end0 - planned$plan$tile2_start0), 20)
  expect_true(any(grepl("chr1:1-20\\|chr1:81-95", planned$plan$focus)))
  expect_true(any(grepl("chr1:1-20\\|chr2:21-40", planned$plan$focus)))
  assigned <- unlist(planned$plan$required_keys, use.names = FALSE)
  expect_setequal(assigned, unique(query$.pixel_key))
  expect_false(anyDuplicated(assigned) > 0L)
})

test_that("native pixel import filters extras and rejects conflicting duplicates", {
  required <- ACV2G:::.contact_pixel_key("chr1", 0, "chr1", 10)
  frame <- data.frame(
    seqnames1 = "chr1", start1 = 1, end1 = 10,
    seqnames2 = "chr1", start2 = c(11, 21), end2 = c(20, 30),
    count = c(3, 999), stringsAsFactors = FALSE
  )
  filtered <- ACV2G:::.contact_validate_hce_frame(
    frame, "count", resolution = 10,
    chrom_lengths = c(chr1 = 100),
    required_keys = required, focus = "chr1:1-30"
  )
  expect_equal(nrow(filtered), 1L)
  expect_equal(filtered$count, 3)
  expect_identical(filtered$.pixel_key, required)

  conflicting <- rbind(filtered, transform(
    filtered, count = 4, .query_focus = "chr1:1-20|chr1:11-30"
  ))
  expect_error(
    ACV2G:::.contact_collapse_native_pixels(list(conflicting), "count"),
    "Conflicting native contact values",
    class = "acv2g_validation_error"
  )
})

test_that("declared native bins retain the fixed-width identity at chromosome ends", {
  pairs <- data.frame(
    pair_id = "p", variant_anchor_id = "v", variant_id = "v",
    gene_id_normalized = "g", gene_id = "g",
    chrom = "chr1", pos1 = 95, tss_chrom = "chr1", tss_pos1 = 95,
    tss_id = "t", stringsAsFactors = FALSE
  )
  seed <- ACV2G:::.contact_make_tss_bin_seed(
    pairs, resolution = 10, chrom_lengths = c(chr1 = 95)
  )
  expect_equal(seed$contact_variant_bin_start0, 90)
  expect_equal(seed$contact_variant_bin_end0, 100)
  expect_equal(seed$contact_promoter_bin_end0, 100)
})

test_that("native implicit-zero semantics require named finite chromosome lengths", {
  expect_error(
    ACV2G:::.contact_validate_native_chromosomes("chr1", NULL, "x.cool"),
    "cannot safely distinguish structural zero",
    class = "acv2g_validation_error"
  )
  info <- ACV2G:::.contact_validate_native_chromosomes(
    c("chr1", "chr2"), c(chr1 = 100, chr2 = 80), "x.cool"
  )
  expect_identical(info$chromosomes, c("chr1", "chr2"))
  expect_equal(info$lengths, c(chr1 = 100, chr2 = 80))
})

test_that("Cooler layer roles fall back from absent weights to raw", {
  roles <- ACV2G:::.contact_native_layer_roles(
    c("count", "balanced"), "cool", weights_available = FALSE
  )
  expect_identical(roles$layer, "count")
  expect_warning(
    selected <- ACV2G:::.contact_resolve_layer(
      roles, contact_value = "balanced", contact_layer = NULL
    ),
    class = "acv2g_contact_raw_fallback"
  )
  expect_identical(selected$type, "raw")
})

test_that("real hic extraction preserves raw zero and conservative balanced missingness", {
  skip_if_not_installed("strawr")
  path <- system.file("extdata", "test.hic", package = "strawr")
  skip_if(!nzchar(path) || !file.exists(path), "strawr test.hic is unavailable")
  pairs <- native_contact_pairs()

  raw <- extract_v2g_contacts(
    pairs, path,
    source_type = "hic",
    contact_value = "raw",
    diagonal_policy = "neighbor_max",
    contact_transform = "identity"
  )
  raw_pairs <- raw$pairs[match(pairs$variant_id, raw$pairs$variant_id), ]
  expect_equal(raw_pairs$contact_raw, c(256, 0, 122))
  expect_identical(
    raw_pairs$contact_status,
    c("observed", "observed_zero", "observed")
  )
  expect_identical(raw$qc$contact_layer_resolved, "count")
  expect_identical(raw$qc$contact_normalization, "NONE")
  expect_match(raw$qc$contact_available_layers, "KR")

  expect_error(
    extract_v2g_contacts(
      pairs[1L, ], path,
      source_type = "hic",
      contact_value = "balanced",
      diagonal_policy = "observed",
      contact_transform = "identity"
    ),
    "Multiple balanced contact layers",
    class = "acv2g_validation_error"
  )

  balanced <- extract_v2g_contacts(
    pairs[1:2, ], path,
    source_type = "hic",
    contact_value = "balanced",
    contact_layer = "KR",
    diagonal_policy = "neighbor_max",
    contact_transform = "identity"
  )
  balanced_pairs <- balanced$pairs[
    match(pairs$variant_id[1:2], balanced$pairs$variant_id),
  ]
  expect_equal(balanced_pairs$contact_raw[[1L]], 251.4966, tolerance = 1e-4)
  expect_identical(balanced$qc$contact_normalization, "KR")
  expect_identical(unique(balanced$pairs$contact_normalization), "KR")
  expect_true(is.na(balanced_pairs$contact_raw[[2L]]))
  expect_identical(
    balanced_pairs$contact_status,
    c("observed", "missing_contact_all_tss_bins")
  )
  zero_detail <- balanced$tss_bin_contacts[
    balanced$tss_bin_contacts$variant_id == "zero", , drop = FALSE
  ]
  expect_identical(zero_detail$contact_status, "missing_balance_weight")
})

test_that("failed native tiles remain missing while successful tiles are retained", {
  skip_if_not_installed("strawr")
  path <- system.file("extdata", "test.hic", package = "strawr")
  skip_if(!nzchar(path) || !file.exists(path), "strawr test.hic is unavailable")
  pairs <- native_contact_pairs()[c(2L, 1L), , drop = FALSE]

  testthat::local_mocked_bindings(
    .contact_straw_query = function(
        normalization, path, chr1loc, chr2loc, resolution) {
      if (grepl("^2:", chr2loc)) stop("synthetic tile failure")
      data.frame(x = 2500000, y = 5000000, counts = 7)
    },
    .package = "ACV2G"
  )
  expect_warning(
    result <- extract_v2g_contacts(
      pairs, path,
      source_type = "hic", contact_value = "raw",
      diagonal_policy = "observed", contact_transform = "identity"
    ),
    class = "acv2g_contact_query_warning"
  )
  detail <- result$tss_bin_contacts[
    match(pairs$variant_id, result$tss_bin_contacts$variant_id),
  ]
  expect_true(is.na(detail$contact_raw[[1L]]))
  expect_identical(detail$contact_status, c("missing_contact_query", "observed"))
  expect_equal(detail$contact_raw[[2L]], 7)
  expect_equal(result$qc$contact_query_failure_count, 1L)
  expect_match(result$qc$contact_query_failure_messages, "synthetic tile failure")
})

test_that("direct R hic paths expand a tilde before every backend call", {
  skip_if_not_installed("strawr")
  source <- system.file("extdata", "test.hic", package = "strawr")
  skip_if(!nzchar(source) || !file.exists(source), "strawr test.hic is unavailable")
  fake_home <- tempfile("acv2g-home-")
  dir.create(fake_home)
  on.exit(unlink(fake_home, recursive = TRUE), add = TRUE)
  target <- file.path(fake_home, "test.hic")
  expect_true(file.copy(source, target))
  withr::local_envvar(HOME = fake_home)

  result <- extract_v2g_contacts(
    native_contact_pairs()[1L, ], "~/test.hic",
    source_type = "hic", contact_value = "raw",
    diagonal_policy = "observed", contact_transform = "identity"
  )
  expect_equal(result$pairs$contact_raw, 256)
  expect_identical(result$qc$contact_source_input, "~/test.hic")
  expect_identical(result$qc$contact_source, normalizePath(target))
})

test_that("generated cool smoke resolves raw and source-balanced layers", {
  skip_if_not_installed("strawr")
  skip_if_not_installed("HiCExperiment")
  hic_path <- system.file("extdata", "test.hic", package = "strawr")
  skip_if(!nzchar(hic_path) || !file.exists(hic_path), "test.hic unavailable")
  prefix <- tempfile("acv2g-native-smoke-")
  cool_path <- paste0(prefix, ".cool")
  on.exit(unlink(cool_path), add = TRUE)

  hce <- NULL
  invisible(utils::capture.output(
    hce <- suppressMessages(HiCExperiment::import(
      HiCExperiment::HicFile(hic_path),
      focus = "1:1-10000000", resolution = 2500000
    ))
  ))
  invisible(utils::capture.output(
    suppressMessages(HiCExperiment::export(
      hce, prefix = prefix, format = "cool"
    ))
  ))
  pairs <- native_contact_pairs()[1L, , drop = FALSE]
  raw <- extract_v2g_contacts(
    pairs, cool_path, source_type = "cool", contact_value = "raw",
    diagonal_policy = "observed", contact_transform = "identity"
  )
  balanced <- extract_v2g_contacts(
    pairs, cool_path, source_type = "cool", contact_value = "balanced",
    diagonal_policy = "observed", contact_transform = "identity"
  )
  expect_equal(raw$pairs$contact_raw, 256)
  expect_equal(balanced$pairs$contact_raw, 256)
  expect_true(balanced$qc$contact_weights_source_provided)
  expect_match(balanced$qc$contact_available_layers, "count;balanced")
})
