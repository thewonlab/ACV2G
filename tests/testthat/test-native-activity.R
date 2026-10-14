test_that("BigWig activity extraction reads an exact nonzero signal", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomeInfoDb")
  skip_if_not_installed("IRanges")
  skip_if_not_installed("S4Vectors")

  bigwig_path <- tempfile("acv2g-native-activity-", fileext = ".bw")
  on.exit(unlink(bigwig_path), add = TRUE)

  signal <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 101L, end = 110L)
  )
  GenomeInfoDb::seqlengths(signal) <- c(chr1 = 1000L)
  S4Vectors::mcols(signal)$score <- 4
  rtracklayer::export(signal, bigwig_path, format = "BigWig")

  observed <- ACV2G::extract_v2g_activity(
    variants = data.frame(
      variant_id = "bigwig_variant",
      chrom = "chr1",
      pos1 = 106L,
      stringsAsFactors = FALSE
    ),
    activity = bigwig_path,
    source_type = "auto",
    activity_window_bp = 10L,
    activity_transform = "identity"
  )

  expect_identical(observed$activity_source_type, "bigwig")
  expect_identical(observed$activity_status, "ok")
  expect_equal(observed$activity_raw, 4, tolerance = 0)
  expect_equal(observed$activity_covered_bp, 10, tolerance = 0)
  expect_equal(observed$activity_coverage_fraction, 1, tolerance = 0)
})

test_that("indexed BAM activity extraction reads exact fragment coverage CPM", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("GenomicAlignments")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomeInfoDb")
  skip_if_not_installed("IRanges")

  sam_path <- tempfile("acv2g-native-activity-", fileext = ".sam")
  bam_stem <- tempfile("acv2g-native-activity-")
  cleanup_paths <- c(
    sam_path,
    paste0(bam_stem, ".bam"),
    paste0(bam_stem, ".bam.bai"),
    paste0(bam_stem, ".bai")
  )
  on.exit(unlink(unique(cleanup_paths)), add = TRUE)

  writeLines(
    c(
      "@HD\tVN:1.6\tSO:coordinate",
      "@SQ\tSN:chr1\tLN:1000",
      "read1\t0\tchr1\t101\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII"
    ),
    sam_path,
    useBytes = TRUE
  )
  bam_path <- Rsamtools::asBam(
    sam_path,
    destination = bam_stem,
    overwrite = TRUE,
    indexDestination = TRUE
  )
  cleanup_paths <- c(
    cleanup_paths,
    bam_path,
    paste0(bam_path, ".bai"),
    sub("\\.bam$", ".bai", bam_path, ignore.case = TRUE)
  )

  observed <- ACV2G::extract_v2g_activity(
    variants = data.frame(
      variant_id = "bam_variant",
      chrom = "chr1",
      pos1 = 106L,
      stringsAsFactors = FALSE
    ),
    activity = bam_path,
    source_type = "auto",
    activity_window_bp = 10L,
    activity_transform = "identity"
  )

  expect_identical(observed$activity_source_type, "bam")
  expect_identical(observed$activity_status, "ok")
  expect_equal(observed$activity_raw, 1e6, tolerance = 0)
  expect_equal(observed$activity_covered_bp, 10, tolerance = 0)
  expect_equal(observed$activity_coverage_fraction, 1, tolerance = 0)
  expect_identical(observed$activity_units, "fragment_coverage_cpm")
  expect_identical(observed$activity_bam_mode, "fragment_coverage_cpm")
  expect_identical(
    observed$activity_bam_filter_policy,
    "mapped_primary_qcpass_nonduplicate_nonsupplementary_v1"
  )
  expect_equal(observed$activity_bam_library_size, 1, tolerance = 0)
})
