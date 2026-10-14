test_that("result bundles retain full schema and verify checksums", {
  result <- new_v2g_result(data.frame(
    profile_id = "p1",
    pair_id = "pair1",
    variant_anchor_id = "v1",
    gene_id_normalized = "g1",
    locus_id = "l1",
    chrom = "chr1",
    locus_start0 = 0L,
    locus_pair_z = 0,
    method = "locus_pair_z",
    stringsAsFactors = FALSE
  ))
  output <- tempfile("acv2g-bundle-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  write_v2g_result(result, output)

  expect_true(all(audit_v2g_result(output)$status == "pass"))
  expected <- c(
    "manifest.json", "config.resolved.yaml", "run_status.json",
    "profiles.tsv.gz", "pairs.tsv.gz", "tss_bin_contacts.tsv.gz",
    "variant_record_anchor_map.tsv.gz", "variant_locus_qc.tsv.gz",
    "locus_status.tsv.gz",
    "normalization_qc.tsv.gz", "activity_qc.tsv.gz", "contact_qc.tsv.gz",
    "annotation_qc.tsv.gz", "liftover_qc.tsv.gz", "warnings.tsv.gz",
    "checksums.sha256"
  )
  expect_true(all(file.exists(file.path(output, expected))))

  lazy <- read_v2g_result(output, lazy = TRUE)
  expect_true(lazy$file_backed)
  expect_equal(nrow(v2g_table(lazy, "pairs")), 1L)
  expect_error(write_v2g_result(result, output), class = "acv2g_validation_error")
})

contract_result_tables <- function() {
  list(
    pairs = data.frame(
      profile_id = "p1", pair_id = "pair1", variant_anchor_id = "v1",
      gene_id_normalized = "g1", stringsAsFactors = FALSE
    ),
    profiles = data.frame(
      profile_id = "p1", activity_source = "activity.bw",
      contact_source = "contact.cool", stringsAsFactors = FALSE
    ),
    tss_bin_contacts = data.frame(
      profile_id = "p1", tss_bin_id = "bin1", pair_id = "pair1",
      variant_anchor_id = "v1", gene_id_normalized = "g1",
      stringsAsFactors = FALSE
    ),
    variant_record_anchor_map = data.frame(
      input_role = "target_variants", record_id = "record1",
      variant_anchor_id = "v1",
      stringsAsFactors = FALSE
    ),
    variant_locus_qc = data.frame(
      variant_anchor_id = "v1", variant_key = "hg38:chr1:1",
      variant_locus_status = "assigned_unique_locus",
      locus_membership_count = 1L, locus_assignment_mode = "coordinate",
      stringsAsFactors = FALSE
    ),
    locus_status = data.frame(
      profile_id = "p1", locus_id = "l1", status = "ok",
      stringsAsFactors = FALSE
    ),
    normalization_qc = data.frame(
      profile_id = "p1", method = "locus_pair_z", group_id = "l1",
      status = "ok", stringsAsFactors = FALSE
    )
  )
}

test_that("result constructor enforces required columns and stable keys", {
  tables <- contract_result_tables()
  result <- do.call(new_v2g_result, tables)
  expect_true(all(audit_v2g_result(result)$status == "pass"))

  required_to_remove <- c(
    profiles = "activity_source",
    pairs = "pair_id",
    tss_bin_contacts = "tss_bin_id",
    variant_record_anchor_map = "record_id",
    variant_locus_qc = "variant_key",
    locus_status = "locus_id",
    normalization_qc = "group_id"
  )
  for (name in names(required_to_remove)) {
    broken <- tables
    broken[[name]][[required_to_remove[[name]]]] <- NULL
    expect_error(
      do.call(new_v2g_result, broken),
      regexp = paste0(name, ".*missing required column"),
      class = "acv2g_validation_error"
    )
  }

  for (name in names(required_to_remove)) {
    broken <- tables
    broken[[name]] <- rbind(broken[[name]], broken[[name]])
    expect_error(
      do.call(new_v2g_result, broken),
      regexp = paste0(name, ".*duplicate composite key"),
      class = "acv2g_validation_error"
    )
  }
})

test_that("object and bundle audits expose table schema and key failures", {
  result <- do.call(new_v2g_result, contract_result_tables())

  missing_column <- result
  missing_column$tables$pairs$pair_id <- NULL
  schema_audit <- audit_v2g_result(missing_column)
  expect_equal(
    schema_audit$status[schema_audit$check == "table_schema:pairs"],
    "fail"
  )
  expect_match(
    schema_audit$detail[schema_audit$check == "table_schema:pairs"],
    "pair_id"
  )

  duplicate_key <- result
  duplicate_key$tables$pairs <- rbind(
    duplicate_key$tables$pairs, duplicate_key$tables$pairs
  )
  key_audit <- audit_v2g_result(duplicate_key)
  expect_equal(
    key_audit$status[key_audit$check == "table_key:pairs"],
    "fail"
  )
  expect_error(
    write_v2g_result(duplicate_key, tempfile("invalid-result-")),
    class = "acv2g_validation_error"
  )

  output <- tempfile("acv2g-contract-bundle-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  write_v2g_result(result, output)
  pair_path <- file.path(output, "pairs.tsv.gz")
  pair_table <- data.table::fread(pair_path)
  pair_table <- rbind(pair_table, pair_table)
  data.table::fwrite(
    pair_table, pair_path, sep = "\t", quote = TRUE, na = "NA",
    compress = "gzip"
  )
  bundle_audit <- audit_v2g_result(output, verify_checksums = FALSE)
  expect_equal(
    bundle_audit$status[bundle_audit$check == "table_key:pairs"],
    "fail"
  )

  pair_table <- pair_table[1L, setdiff(names(pair_table), "pair_id"), with = FALSE]
  data.table::fwrite(
    pair_table, pair_path, sep = "\t", quote = TRUE, na = "NA",
    compress = "gzip"
  )
  bundle_schema_audit <- suppressWarnings(
    audit_v2g_result(output, verify_checksums = FALSE)
  )
  expect_equal(
    bundle_schema_audit$status[
      bundle_schema_audit$check == "table_schema:pairs"
    ],
    "fail"
  )
  expect_error(
    suppressWarnings(read_v2g_result(output, verify_checksums = FALSE)),
    class = "acv2g_validation_error"
  )
})

test_that("empty typed result tables retain their types through bundles", {
  result <- new_v2g_result(data.frame(
    profile_id = "p1", pair_id = "pair1", variant_anchor_id = "v1",
    gene_id_normalized = "g1", integer_value = 1L,
    logical_value = TRUE, double_value = 1.5,
    stringsAsFactors = FALSE
  ))
  expect_type(
    v2g_table(result, "variant_locus_qc")$locus_membership_count,
    "integer"
  )

  output <- tempfile("acv2g-typed-bundle-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  write_v2g_result(result, output)
  eager <- read_v2g_result(output)
  lazy <- read_v2g_result(output, lazy = TRUE)
  expect_type(v2g_table(eager, "pairs")$integer_value, "integer")
  expect_type(v2g_table(eager, "pairs")$logical_value, "logical")
  expect_type(v2g_table(eager, "pairs")$double_value, "double")
  expect_type(
    v2g_table(lazy, "variant_locus_qc")$locus_membership_count,
    "integer"
  )
})

test_that("bundle reader preserves long multiline warning messages", {
  multiline_message <- paste(
    c(
      "Five groups exceeded the warning threshold.",
      "Representative groups:",
      paste0("profile=p1|locus=l", seq_len(5L))
    ),
    collapse = "\n"
  )
  warnings <- data.frame(
    schema_version = "1.0.0",
    warning_id = "warning_000001",
    class = "acv2g_warning;warning;condition",
    message = multiline_message,
    profile_id = "p1",
    locus_id = NA_character_,
    condition_data_json = NA_character_,
    timestamp_utc = "2026-08-06 16:00:00 UTC",
    stringsAsFactors = FALSE
  )
  result <- new_v2g_result(
    data.frame(
      profile_id = "p1", pair_id = "pair1", variant_anchor_id = "v1",
      gene_id_normalized = "g1", stringsAsFactors = FALSE
    ),
    warnings = warnings
  )
  output <- tempfile("acv2g-multiline-warning-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)

  write_v2g_result(result, output)
  restored <- read_v2g_result(output)

  expect_true(all(audit_v2g_result(output)$status == "pass"))
  expect_identical(v2g_table(restored, "warnings")$message, multiline_message)
})

test_that("bundle audit requires complete checksum coverage", {
  result <- new_v2g_result(data.frame(
    profile_id = "p1", pair_id = "pair1", variant_anchor_id = "v1",
    gene_id_normalized = "g1", stringsAsFactors = FALSE
  ))
  output <- tempfile("acv2g-checksum-coverage-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  write_v2g_result(result, output)
  checksum_path <- file.path(output, "checksums.sha256")
  lines <- readLines(checksum_path)
  writeLines(lines[-1L], checksum_path)
  audit <- audit_v2g_result(output, verify_checksums = TRUE)
  expect_equal(audit$status[audit$check == "checksums"], "fail")
  expect_match(
    audit$detail[audit$check == "checksums"],
    "cover each payload exactly once"
  )
})

test_that("profile model requires exact explicit pairing", {
  profiles <- validate_v2g_profiles(data.frame(
    profile_id = c("a", "b"),
    activity_source = c("a1.bw", "a2.bw"),
    contact_source = c("c1.cool", "c2.cool")
  ))
  expect_equal(profiles$profile_id, c("a", "b"))
  expect_error(
    validate_v2g_profiles(data.frame(
      profile_id = c("a", "a"), activity_source = "x", contact_source = "y"
    )),
    class = "acv2g_validation_error"
  )
})

test_that("resolved configuration fixes defaults without guessing assemblies", {
  files <- vapply(c("loci.tsv", "variants.tsv", "genes.tsv", "activity.tsv",
                    "contact.tsv", "chrom_sizes.tsv"), function(name) {
    path <- tempfile(name)
    writeLines("placeholder", path)
    path
  }, character(1L))
  on.exit(unlink(files), add = TRUE)
  config <- list(
    analysis_build = "hg38",
    inputs = list(
      loci = list(path = files[[1L]], assembly = "hg38"),
      variants = list(path = files[[2L]], assembly = "hg38"),
      genes = list(path = files[[3L]], assembly = "hg38"),
      chrom_sizes = list(path = files[[6L]], assembly = "hg38")
    ),
    profiles = list(
      neuron = list(
        activity = list(
          path = files[[4L]], type = "bedgraph", assembly = "hg38"
        ),
        contact = list(
          path = files[[5L]], type = "cool", assembly = "hg38"
        )
      )
    )
  )
  resolved <- validate_v2g_config(config)
  expect_equal(resolved$method$name, "locus_pair_z")
  expect_equal(resolved$method$gene_flank_bp, 1000000)
  expect_equal(resolved$activity$activity_window_bp, 500)
  expect_equal(resolved$contact$distance_policy, "none")
  expect_equal(names(resolved$profiles), "neuron")

  config$inputs$variants$assembly <- NULL
  expect_error(validate_v2g_config(config), class = "acv2g_validation_error")
})

test_that("locus plot has no implicit title or threshold", {
  pairs <- data.frame(
    profile_id = "p1",
    locus_id = rep(c("l1", "l2"), each = 3),
    chrom = rep(c("chr1", "chr2"), each = 3),
    locus_start0 = rep(c(10L, 20L), each = 3),
    locus_pair_z = c(-1, 0, 1, 0, 1.5, 2),
    method = "locus_pair_z"
  )
  plot <- plot_locus_z(pairs)
  expect_s3_class(plot, "ggplot")
  expect_null(plot$labels$title)
  expect_null(attr(plot, "acv2g_plot_provenance")$display_threshold)

  highlighted <- plot_locus_z(pairs, display_threshold = 1.5,
                              max_points_per_locus = 2, seed = 7)
  provenance <- attr(highlighted, "acv2g_plot_provenance")
  expect_equal(provenance$display_threshold, 1.5)
  expect_equal(provenance$plotted_points, 4L)
})
