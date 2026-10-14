.structural_fixture <- function() {
  loci <- read_v2g_loci(data.frame(
    locus_id = c("locus_1", "locus_2"),
    chrom = c("chr1", "chr1"),
    start0 = c(0L, 200L),
    end0 = c(100L, 300L)
  ), assembly = "hg19")
  variants <- read_v2g_variants(data.frame(
    variant_id = c("v1", "v2"),
    chrom = c("chr1", "chr1"),
    pos1 = c(50L, 250L)
  ), assembly = "hg19")
  genes <- read_v2g_genes(data.frame(
    gene_id = c("G1", "G2"),
    chrom = c("chr1", "chr1"),
    tss_pos1 = c(121L, 180L),
    strand = c("+", "+"),
    gene_biotype = c("protein_coding", "protein_coding")
  ), assembly = "hg19")
  list(loci = loci, variants = variants, genes = genes)
}

.structural_example_config_path <- function() {
  source_path <- testthat::test_path(
    "..", "..", "inst", "extdata", "example-config.yaml"
  )
  if (file.exists(source_path)) {
    return(normalizePath(source_path, mustWork = TRUE))
  }
  system.file(
    "extdata", "example-config.yaml",
    package = "ACV2G",
    mustWork = TRUE
  )
}

.tiny_chain <- function(
    source_chrom,
    target_chrom,
    source_start,
    source_end,
    offset,
    reversed = FALSE) {
  block <- methods::new(
    "ChainBlock",
    ranges = IRanges::IRanges(source_start, source_end),
    offset = as.integer(offset),
    score = 1L,
    space = target_chrom,
    reversed = reversed,
    length = as.integer(source_end - source_start + 1L)
  )
  methods::new(
    "Chain",
    listData = stats::setNames(list(block), source_chrom),
    elementType = "ChainBlock",
    elementMetadata = NULL,
    metadata = list()
  )
}

test_that("failed target anchors and gene TSSs fail their whole source loci", {
  fixture <- .structural_fixture()
  failed_variant <- fixture$variants$anchors$variant_anchor_id[[1L]]
  failed_tss <- fixture$genes$tss$tss_id[
    fixture$genes$tss$gene_id_normalized == "G2"
  ]
  attributed <- ACV2G:::.attribute_structural_failures(
    fixture$loci,
    fixture$variants,
    fixture$genes,
    failed_variant_ids = failed_variant,
    failed_tss_ids = failed_tss,
    gene_flank_bp = 25L
  )

  expect_setequal(attributed$failed_loci, c("locus_1", "locus_2"))
  expect_true(any(
    attributed$record_locus_map$input_role == "target_variants" &
      attributed$record_locus_map$locus_id == "locus_1"
  ))
  expect_true(any(
    attributed$record_locus_map$input_role == "genes" &
      attributed$record_locus_map$locus_id == "locus_2"
  ))
})

test_that("unsafe failed-record attribution is a validation error", {
  fixture <- .structural_fixture()
  fixture$loci$locus_assembly <- "hg18"
  failed_variant <- fixture$variants$anchors$variant_anchor_id[[1L]]
  expect_error(
    ACV2G:::.attribute_structural_failures(
      fixture$loci,
      fixture$variants,
      fixture$genes,
      failed_variant_ids = failed_variant,
      gene_flank_bp = 25L
    ),
    "cannot be safely attributed",
    class = "acv2g_validation_error"
  )

  explicit <- data.frame(variant_id = "v1", locus_id = "locus_1")
  attributed <- ACV2G:::.attribute_structural_failures(
    fixture$loci,
    fixture$variants,
    fixture$genes,
    failed_variant_ids = failed_variant,
    variant_locus_map = explicit,
    gene_flank_bp = 25L
  )
  expect_identical(attributed$failed_loci, "locus_1")

  failed_tss <- fixture$genes$tss$tss_id[[1L]]
  fixture$genes$tss$assembly <- "hg17"
  expect_error(
    ACV2G:::.attribute_structural_failures(
      fixture$loci,
      fixture$variants,
      fixture$genes,
      failed_tss_ids = failed_tss,
      gene_flank_bp = 25L
    ),
    "cannot be safely attributed",
    class = "acv2g_validation_error"
  )
})

test_that("fail_locus removes a locus affected by a failed target anchor", {
  config <- validate_v2g_config(.structural_example_config_path())
  config$liftover$structural_failure_policy <- "fail_locus"
  original_liftover <- ACV2G:::liftover_v2g_coordinates
  testthat::local_mocked_bindings(
    liftover_v2g_coordinates = function(x, ..., id_col = NULL) {
      result <- original_liftover(x, ..., id_col = id_col)
      if (identical(id_col, "variant_anchor_id")) {
        result$data$liftover_status[[1L]] <- "unmapped"
        result$data$chrom[[1L]] <- NA_character_
        result$data$pos1[[1L]] <- NA_integer_
        result$qc$status[[1L]] <- "unmapped"
        result$qc$target_chrom[[1L]] <- NA_character_
        result$qc$target_start1[[1L]] <- NA_integer_
        result$qc$target_end1[[1L]] <- NA_integer_
        result$qc$mapping_count[[1L]] <- 0L
      }
      result
    },
    .package = "ACV2G"
  )

  structural <- ACV2G:::.prepare_structural_inputs(config)
  expect_identical(structural$failed_loci, "locus_1")
  expect_equal(nrow(structural$loci), 0L)
  variant_failure <- structural$liftover_qc[
    structural$liftover_qc$input_role == "target_variants" &
      structural$liftover_qc$status == "unmapped", , drop = FALSE
  ]
  expect_equal(nrow(variant_failure), 1L)
  expect_identical(variant_failure$affected_locus_ids, "locus_1")
  expect_true(all(structural$liftover_qc$liftover_role == "structural"))
})

test_that("fail_locus returns an auditable zero-pair partial result when all loci fail", {
  config <- validate_v2g_config(.structural_example_config_path())
  config$liftover$structural_failure_policy <- "fail_locus"
  original_liftover <- ACV2G:::liftover_v2g_coordinates
  testthat::local_mocked_bindings(
    liftover_v2g_coordinates = function(x, ..., id_col = NULL) {
      result <- original_liftover(x, ..., id_col = id_col)
      if (identical(id_col, "variant_anchor_id")) {
        result$data$liftover_status[[1L]] <- "unmapped"
        result$data$chrom[[1L]] <- NA_character_
        result$data$pos1[[1L]] <- NA_integer_
        result$qc$status[[1L]] <- "unmapped"
        result$qc$target_chrom[[1L]] <- NA_character_
        result$qc$target_start1[[1L]] <- NA_integer_
        result$qc$target_end1[[1L]] <- NA_integer_
        result$qc$mapping_count[[1L]] <- 0L
      }
      result
    },
    .package = "ACV2G"
  )

  expect_warning(
    result <- predict_v2g(config),
    class = "acv2g_all_loci_failed_warning"
  )
  expect_identical(result$run_status$status, "partial")
  expect_equal(nrow(v2g_table(result, "pairs")), 0L)
  expect_true(all(
    v2g_table(result, "locus_status")$status ==
      "failed_structural_liftover"
  ))
  expect_true(all(audit_v2g_result(result)$status == "pass"))
})

test_that("cross-build activity requires correctly directed query chains", {
  config <- read_v2g_config(.structural_example_config_path())
  config$profiles$example_profile$activity$assembly <- "hg19"
  expect_error(
    validate_v2g_config(config),
    "hg38 [(]analysis_build[)] -> hg19 [(]activity assembly[)]",
    class = "acv2g_validation_error"
  )
  chain <- .tiny_chain(
    "chr1", "chrActivity", 1L, 200000L, -100L, reversed = FALSE
  )
  config$profiles$example_profile$activity$chain <- chain
  config$profiles$example_profile$activity$chrom_sizes <- data.frame(
    chrom = "chrActivity", length = 1000000L
  )
  resolved <- validate_v2g_config(config)
  expect_s4_class(resolved$profiles$example_profile$activity$chain, "Chain")

  plan <- validate_v2g_assemblies(
    "hg38",
    c(contact = "hg38", activity = "hg19"),
    chains = list(activity = chain)
  )
  activity <- plan[plan$input_role == "activity", , drop = FALSE]
  expect_identical(activity$source_build, "hg38")
  expect_identical(activity$target_build, "hg19")
})

test_that("liftover updates strand across inversions and records reciprocity", {
  forward <- .tiny_chain(
    "chrSource", "chrTarget", 1L, 100L, -100L, reversed = TRUE
  )
  reverse <- .tiny_chain(
    "chrTarget", "chrSource", 101L, 200L, 100L, reversed = TRUE
  )
  lifted <- liftover_v2g_coordinates(
    data.frame(id = "tss1", chrom = "chrSource", pos1 = 11L, strand = "+"),
    chain = forward,
    source_build = "source_build",
    target_build = "target_build",
    id_col = "id",
    strand_col = "strand",
    reciprocal_policy = "require",
    reverse_chain = reverse
  )
  expect_identical(lifted$data$pos1, 190L)
  expect_identical(lifted$data$strand, "-")
  expect_identical(lifted$data$strand_source, "+")
  expect_identical(lifted$qc$reciprocal_status, "passed")
  expect_true(lifted$qc$reciprocal_ok)
  expect_identical(lifted$qc$chain_direction,
                   "source_build->target_build")
  expect_identical(lifted$qc$reverse_chain_direction,
                   "target_build->source_build")
  expect_true(nzchar(lifted$qc$chain_sha256))
  expect_true(nzchar(lifted$qc$reverse_chain_sha256))
})

test_that("activity query liftover preserves universe and makes failures NA", {
  activity_path <- withr::local_tempfile(fileext = ".tsv")
  data.table::fwrite(data.frame(
    chrom = "chrTarget", pos1 = 190L, activity_raw = 4
  ), activity_path, sep = "\t")
  candidates <- list(pairs = data.frame(
    variant_anchor_id = c("anchor_mapped", "anchor_failed"),
    chrom = "chrSource",
    pos1 = c(11L, 150L),
    stringsAsFactors = FALSE
  ))
  profile <- list(
    profile_id = "profile_1",
    activity = list(
      path = activity_path,
      type = "table",
      table_mode = "variant",
      value_col = "activity_raw",
      assembly = "activity_build",
      chain = .tiny_chain(
        "chrSource", "chrTarget", 1L, 100L, -100L, reversed = TRUE
      ),
      chrom_sizes = data.frame(chrom = "chrTarget", length = 1000L)
    )
  )
  config <- list(
    analysis_build = "analysis_build",
    activity = list(
      activity_window_bp = 500L,
      aggregation = "window_mean",
      bam_mapq_min = 30L
    ),
    method = list(proxy_transform = "log1p"),
    liftover = list(
      structural_failure_policy = "error",
      reciprocal_liftover_policy = "none"
    )
  )
  expect_warning(
    result <- ACV2G:::.activity_for_profile(
      candidates, profile, config, chrom_lengths = NULL
    ),
    class = "acv2g_warning"
  )
  activity <- result$value
  expect_identical(activity$variant_anchor_id,
                   c("anchor_mapped", "anchor_failed"))
  expect_equal(activity$activity_raw[[1L]], 4)
  expect_true(is.na(activity$activity_raw[[2L]]))
  expect_identical(
    activity$activity_status,
    c("ok_precomputed", "liftover_unmapped")
  )
  expect_identical(activity$chrom, rep("chrSource", 2L))
  expect_identical(activity$activity_query_chrom,
                   c("chrTarget", NA_character_))
  expect_identical(result$qc$input_role, rep("activity_query", 2L))
  expect_identical(result$qc$source_build, rep("analysis_build", 2L))
  expect_identical(result$qc$target_build, rep("activity_build", 2L))
})
