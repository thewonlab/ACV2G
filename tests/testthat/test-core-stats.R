test_that("compute_v2g_proxy gives exact transformed products and preserves missingness", {
  pairs <- data.frame(
    pair_id = c("p0", "p1", "pm", "pz"),
    activity_value = c(0, 3, NA, 0),
    contact_value = c(8, 8, 8, NA)
  )

  out <- compute_v2g_proxy(pairs)

  expect_equal(out$activity_transformed, c(0, log(4), NA, 0))
  expect_equal(out$contact_transformed, c(log(9), log(9), log(9), NA))
  expect_equal(out$proxy_score[1:2], c(0, log(4) * log(9)))
  expect_true(is.na(out$proxy_score[[3L]]))
  expect_true(is.na(out$proxy_score[[4L]]))
  expect_identical(out$pair_id, pairs$pair_id)

  identity <- compute_v2g_proxy(
    pairs,
    activity_transform = "identity",
    contact_transform = "identity"
  )
  expect_equal(identity$proxy_score[1:2], c(0, 24))
})

test_that("compute_v2g_proxy rejects invalid input instead of silently transforming it", {
  negative <- data.frame(
    pair_id = c("ok", "bad"),
    activity_value = c(1, -0.1),
    contact_value = c(1, 1)
  )
  expect_error(
    compute_v2g_proxy(negative),
    "finite negative value",
    class = "acv2g_validation_error"
  )
  expect_error(
    compute_v2g_proxy(
      transform(negative, activity_value = c(1, Inf)),
      activity_transform = "identity"
    ),
    "non-finite, non-missing",
    class = "acv2g_validation_error"
  )
  expect_error(
    compute_v2g_proxy(
      data.frame(activity_value = c(.Machine$double.xmax), contact_value = 2),
      activity_transform = "identity",
      contact_transform = "identity"
    ),
    "multiplication produced",
    class = "acv2g_validation_error"
  )
})

test_that("locus_pair_z uses finite values and sample SD while retaining all rows", {
  pairs <- data.frame(
    row_id = 1:6,
    profile_id = "profile_a",
    locus_id = c(rep("locus_1", 4), rep("locus_2", 2)),
    proxy_score = c(0, 1, 2, NA, 5, 5)
  )

  expect_warning(
    result <- compute_locus_pair_z(pairs),
    class = "acv2g_missing_fraction_warning"
  )

  expect_s3_class(result, "acv2g_normalization_result")
  expect_identical(result$pairs$row_id, pairs$row_id)
  expect_equal(result$pairs$locus_pair_z[1:4], c(-1, 0, 1, NA))
  expect_true(all(is.na(result$pairs$locus_pair_z[5:6])))

  qc1 <- subset(result$normalization_qc, locus_id == "locus_1")
  expect_equal(qc1$normalization_n_total, 4L)
  expect_equal(qc1$normalization_n_valid, 3L)
  expect_equal(qc1$normalization_missing_fraction, 0.25)
  expect_equal(qc1$normalization_mean, 1)
  expect_equal(qc1$normalization_sd, 1)
  expect_identical(qc1$normalization_status, "ok")
  expect_identical(qc1$normalization_warnings, "high_missing_fraction")
  expect_identical(qc1$group_id, qc1$normalization_group)
  expect_identical(qc1$status, qc1$normalization_status)

  qc2 <- subset(result$normalization_qc, locus_id == "locus_2")
  expect_identical(qc2$normalization_status, "constant_proxy")
  expect_true(qc2$normalization_degenerate)
})

test_that("missing-fraction warning comparison is strictly greater than the threshold", {
  pairs <- data.frame(
    profile_id = "p",
    locus_id = "l",
    proxy_score = c(0, 1, 2, 3, NA)
  )
  expect_no_warning(result <- compute_locus_pair_z(pairs))
  expect_equal(result$normalization_qc$normalization_missing_fraction, 0.2)
  expect_identical(result$normalization_qc$normalization_warnings, "")
  expect_true(all(is.finite(result$pairs$locus_pair_z[1:4])))
})

test_that("locus_pair_z reports insufficient and all-missing groups", {
  pairs <- data.frame(
    profile_id = "p",
    locus_id = c("one", "missing", "missing"),
    proxy_score = c(4, NA, NA)
  )
  expect_warning(result <- compute_locus_pair_z(pairs))
  status <- setNames(
    result$normalization_qc$normalization_status,
    result$normalization_qc$locus_id
  )
  expect_identical(status[["one"]], "insufficient_finite_proxy")
  expect_identical(status[["missing"]], "all_proxy_missing")
  expect_true(all(is.na(result$pairs$locus_pair_z)))
})

test_that("gene_wise_z requires an explicit reference mode", {
  pairs <- data.frame(
    profile_id = "p",
    gene_id = "g",
    variant_key = "v",
    chrom = "chr1",
    pos1 = 100L,
    tss_chrom = "chr1",
    tss_pos1 = 100L,
    proxy_score = 1
  )
  expect_error(
    compute_gene_wise_z(pairs),
    "must be explicitly set",
    class = "acv2g_validation_error"
  )
  expect_error(
    compute_gene_wise_z(pairs, reference_mode = "background_variants"),
    "background_pairs.*required",
    class = "acv2g_validation_error"
  )
})

test_that("gene_wise_z distance universe is inclusive and does not drop target rows", {
  pairs <- data.frame(
    row_id = 1:3,
    profile_id = "p",
    gene_id = "g",
    variant_key = c("left", "right", "outside"),
    chrom = "chr1",
    pos1 = c(1L, 2000001L, 2000002L),
    tss_chrom = "chr1",
    tss_pos1 = 1000001L,
    proxy_score = c(0, 2, 100)
  )

  result <- compute_gene_wise_z(pairs, reference_mode = "input_variants")
  expect_identical(result$pairs$row_id, 1:3)
  expect_equal(
    result$pairs$gene_wise_z,
    c(-1 / sqrt(2), 1 / sqrt(2), 99 / sqrt(2))
  )
  expect_equal(result$normalization_qc$normalization_n_total, 2L)
  expect_equal(result$normalization_qc$normalization_mean, 1)
  expect_equal(result$normalization_qc$normalization_sd, sqrt(2))
  expect_identical(
    result$normalization_qc$reference_universe_rule,
    "distance"
  )
  expect_equal(result$normalization_qc$reference_gene_radius_bp, 1000000)
})

test_that("gene_wise_z uses any supplied TSS for distance eligibility", {
  pairs <- data.frame(
    profile_id = "p",
    gene_id = "g",
    variant_key = c("v1", "v2", "v3"),
    chrom = "chr1",
    pos1 = c(100L, 1000L, 5000L),
    proxy_score = c(0, 2, 100)
  )
  tss <- data.frame(
    gene_id = c("g", "g"),
    tss_chrom = c("chr1", "chr1"),
    tss_pos1 = c(100L, 1000L)
  )
  result <- compute_gene_wise_z(
    pairs,
    reference_mode = "input_variants",
    gene_tss = tss,
    reference_gene_radius_bp = 0L
  )
  expect_equal(result$normalization_qc$normalization_n_total, 2L)
  expect_equal(result$normalization_qc$normalization_mean, 1)
})

test_that("explicit mapping completely overrides distance fields and deduplicates anchors", {
  pairs <- data.frame(
    row_id = 1:4,
    profile_id = "p",
    gene_id = "g",
    variant_key = c("v1", "v1", "v2", "v3"),
    proxy_score = c(0, 0, 100, 4)
  )
  mapping <- data.frame(
    gene_id = c("g", "g"),
    variant_key = c("v1", "v3")
  )

  result <- compute_gene_wise_z(
    pairs,
    reference_mode = "input_variants",
    reference_gene_variant_map = mapping,
    reference_gene_variant_map_id = "map_v1"
  )
  expect_equal(result$normalization_qc$normalization_n_total, 2L)
  expect_equal(result$normalization_qc$normalization_mean, 2)
  expect_equal(result$normalization_qc$normalization_sd, sqrt(8))
  expect_identical(
    result$normalization_qc$reference_universe_rule,
    "explicit_mapping"
  )
  expect_identical(
    result$normalization_qc$reference_gene_variant_map_id,
    "map_v1"
  )
  expect_true(is.na(result$normalization_qc$reference_gene_radius_bp))
  expect_identical(result$pairs$row_id, 1:4)
})

test_that("background reference statistics are applied to full target rows", {
  targets <- data.frame(
    profile_id = "p",
    gene_id = "g",
    variant_key = c("target_1", "target_2"),
    proxy_score = c(2, NA)
  )
  background <- data.frame(
    profile_id = "p",
    gene_id = "g",
    variant_key = c("b1", "b2"),
    chrom = "chr1",
    pos1 = c(90L, 110L),
    tss_chrom = "chr1",
    tss_pos1 = 100L,
    proxy_score = c(0, 4)
  )

  result <- compute_gene_wise_z(
    targets,
    reference_mode = "background_variants",
    background_pairs = background,
    reference_set_id = "control_v1"
  )
  expect_equal(result$pairs$gene_wise_z, c(0, NA))
  expect_identical(result$pairs$reference_mode, rep("background_variants", 2))
  expect_identical(result$normalization_qc$reference_set_id, "control_v1")
})

test_that("gene_wise_z retains target groups with no eligible reference", {
  targets <- data.frame(
    profile_id = "p",
    gene_id = "g",
    variant_key = "target",
    proxy_score = 2
  )
  background <- data.frame(
    profile_id = "p",
    gene_id = "g",
    variant_key = "far",
    chrom = "chr1",
    pos1 = 3000000L,
    tss_chrom = "chr1",
    tss_pos1 = 1L,
    proxy_score = 4
  )
  result <- compute_gene_wise_z(
    targets,
    reference_mode = "background_variants",
    background_pairs = background
  )
  expect_equal(nrow(result$pairs), 1L)
  expect_true(is.na(result$pairs$gene_wise_z))
  expect_identical(
    result$normalization_qc$normalization_status,
    "no_reference_variants"
  )
})

test_that("call_v2g_links defaults to no calls", {
  pairs <- data.frame(locus_pair_z = c(-1, 1.5, 2, NA))
  out <- call_v2g_links(pairs)
  expect_true(all(is.na(out$link_call)))
  expect_identical(unique(out$call_rule), "none")
  expect_true(all(is.na(out$call_threshold)))
  expect_identical(unique(out$call_comparator), "none")
  expect_error(
    call_v2g_links(pairs, z_threshold = 1.5),
    "must be `NULL`",
    class = "acv2g_validation_error"
  )
})

test_that("z-threshold calling uses the conditional 1.5 default and inclusive comparator", {
  pairs <- data.frame(gene_wise_z = c(-1, 1.499, 1.5, 2, NA))
  out <- call_v2g_links(
    pairs,
    call_rule = "z_threshold",
    z_col = "gene_wise_z"
  )
  expect_identical(out$link_call, c(FALSE, FALSE, TRUE, TRUE, NA))
  expect_equal(unique(out$call_threshold), 1.5)
  expect_identical(unique(out$call_threshold_source), "conditional_default")
  expect_identical(unique(out$call_comparator), ">=")

  custom <- call_v2g_links(
    pairs,
    call_rule = "z_threshold",
    z_col = "gene_wise_z",
    z_threshold = 2
  )
  expect_identical(custom$link_call, c(FALSE, FALSE, FALSE, TRUE, NA))
  expect_identical(unique(custom$call_threshold_source), "user")
})
