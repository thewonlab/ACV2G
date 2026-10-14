test_that("activity windows use exact half-open geometry and clip without shifting", {
  variants <- data.frame(
    variant_id = c("edge", "middle"),
    chrom = "chr1",
    pos1 = c(1, 11)
  )
  empty_signal <- data.frame(
    chrom = character(), start0 = numeric(), end0 = numeric(), score = numeric()
  )

  even <- extract_v2g_activity(
    variants, empty_signal,
    chrom_lengths = c(chr1 = 100),
    activity_window_bp = 4,
    activity_transform = "identity"
  )
  expect_equal(even$activity_window_requested_start0, c(-2, 8))
  expect_equal(even$activity_window_requested_end0, c(2, 12))
  expect_equal(even$activity_window_effective_start0, c(0, 8))
  expect_equal(even$activity_window_effective_end0, c(2, 12))
  expect_equal(even$activity_effective_window_bp, c(2, 4))
  expect_equal(even$activity_window_clipped, c(TRUE, FALSE))

  odd <- extract_v2g_activity(
    variants[2, , drop = FALSE], empty_signal,
    chrom_lengths = c(chr1 = 100),
    activity_window_bp = 5,
    activity_transform = "identity"
  )
  expect_equal(odd$activity_window_requested_start0, 8)
  expect_equal(odd$activity_window_requested_end0, 13)
})

test_that("activity window_mean and covered_mean have distinct denominators", {
  variant <- data.frame(variant_id = "v", chrom = "chr1", pos1 = 101)
  signal <- data.frame(
    chrom = "chr1", start0 = 95, end0 = 100, score = 4
  )
  window <- extract_v2g_activity(
    variant, signal,
    chrom_lengths = c(chr1 = 1000),
    activity_window_bp = 10,
    activity_aggregate = "window_mean",
    activity_transform = "identity"
  )
  covered <- extract_v2g_activity(
    variant, signal,
    chrom_lengths = c(chr1 = 1000),
    activity_window_bp = 10,
    activity_aggregate = "covered_mean",
    activity_transform = "identity"
  )
  expect_equal(window$activity_raw, 2)
  expect_equal(covered$activity_raw, 4)
  expect_equal(window$activity_covered_bp, 5)
  expect_equal(window$activity_coverage_fraction, 0.5)
  expect_equal(window$activity_uncovered_policy, "zero")
  expect_equal(covered$activity_uncovered_policy, "exclude")
})

test_that("valid uncovered activity is zero only for window_mean", {
  variant <- data.frame(variant_id = "v", chrom = "chr1", pos1 = 101)
  signal <- data.frame(
    chrom = "chr1", start0 = 0, end0 = 10, score = 8
  )
  window <- extract_v2g_activity(
    variant, signal, chrom_lengths = c(chr1 = 1000),
    activity_window_bp = 10, activity_transform = "identity"
  )
  covered <- extract_v2g_activity(
    variant, signal, chrom_lengths = c(chr1 = 1000),
    activity_window_bp = 10, activity_aggregate = "covered_mean",
    activity_transform = "identity"
  )
  expect_identical(window$activity_raw, 0)
  expect_identical(window$activity_status, "ok_uncovered_zero")
  expect_true(is.na(covered$activity_raw))
  expect_identical(covered$activity_status, "no_covered_bases")
})

test_that("absent activity seqlevels remain missing rather than zero", {
  variants <- data.frame(
    variant_id = c("present", "absent"),
    chrom = c("chr1", "chr2"),
    pos1 = c(101, 101)
  )
  signal <- data.frame(
    chrom = "chr1", start0 = 0, end0 = 10, score = 8
  )
  out <- extract_v2g_activity(
    variants, signal,
    chrom_lengths = c(chr1 = 1000, chr2 = 1000),
    activity_window_bp = 10,
    activity_transform = "identity"
  )
  expect_equal(out$activity_raw[1], 0)
  expect_identical(out$activity_status[1], "ok_uncovered_zero")
  expect_true(is.na(out$activity_raw[2]))
  expect_true(is.na(out$activity_covered_bp[2]))
  expect_true(is.na(out$activity_coverage_fraction[2]))
  expect_identical(out$activity_status[2], "missing_activity_seqlevel")
})

test_that("precomputed activity retains missing values and coverage uncertainty", {
  variants <- data.frame(
    variant_id = c("v1", "v2", "v3"), chrom = "chr1", pos1 = c(10, 20, 30)
  )
  values <- data.frame(
    variant_id = c("v1", "v2"), activity_raw = c(0, NA_real_)
  )
  out <- extract_v2g_activity(
    variants, values, table_mode = "variant",
    chrom_lengths = c(chr1 = 100), activity_transform = "identity"
  )
  expect_equal(out$activity_raw[1], 0)
  expect_true(all(is.na(out$activity_raw[2:3])))
  expect_identical(
    out$activity_status,
    c("ok_precomputed", "missing_precomputed_value", "missing_precomputed_variant_id")
  )
  expect_true(all(is.na(out$activity_coverage_fraction)))
})

test_that("precomputed activity prefers explicit coordinates over ambiguous IDs", {
  variants <- data.frame(
    variant_id = c("internal_a", "internal_b"),
    chrom = "chr1",
    pos1 = c(10, 20)
  )
  values <- data.frame(
    variant_id = c("submitted_b", "submitted_a"),
    chrom = "chr1",
    pos1 = c(20, 10),
    activity_raw = c(8, 3)
  )
  out <- extract_v2g_activity(
    variants, values, table_mode = "variant",
    chrom_lengths = c(chr1 = 100), activity_transform = "identity"
  )
  expect_equal(out$activity_raw, c(3, 8))
  expect_identical(out$activity_status, rep("ok_precomputed", 2L))
})

test_that("activity log1p rejects finite negative values", {
  variant <- data.frame(variant_id = "v", chrom = "chr1", pos1 = 10)
  signal <- data.frame(chrom = "chr1", start0 = 0, end0 = 20, score = -0.1)
  expect_error(
    extract_v2g_activity(variant, signal, chrom_lengths = c(chr1 = 100)),
    class = "acv2g_validation_error"
  )
  identity <- extract_v2g_activity(
    variant, signal, chrom_lengths = c(chr1 = 100),
    activity_transform = "identity"
  )
  expect_lt(identity$activity_raw, 0)
})

test_that("activity log1p rejects negative source intervals before aggregation", {
  variant <- data.frame(variant_id = "v", chrom = "chr1", pos1 = 2)
  cancelling <- data.frame(
    chrom = c("chr1", "chr1"),
    start0 = c(0, 1), end0 = c(1, 2),
    score = c(-1, 1)
  )
  expect_error(
    extract_v2g_activity(
      variant, cancelling, table_mode = "interval",
      chrom_lengths = c(chr1 = 100), activity_window_bp = 2,
      activity_transform = "log1p"
    ),
    "before window aggregation",
    class = "acv2g_validation_error"
  )
  identity <- extract_v2g_activity(
    variant, cancelling, table_mode = "interval",
    chrom_lengths = c(chr1 = 100), activity_window_bp = 2,
    activity_transform = "identity"
  )
  expect_equal(identity$activity_raw, 0)
})

test_that("activity query failures warn and remain explicitly missing", {
  measurement <- ACV2G:::.activity_empty_measurement(2L, "window_mean")
  expect_warning(
    failed <- ACV2G:::.activity_mark_query_failure(
      measurement, c(TRUE, FALSE), "bigwig", "signal.bw",
      simpleError("synthetic query failure")
    ),
    class = "acv2g_activity_query_warning"
  )
  expect_true(is.na(failed$value[[1L]]))
  expect_identical(
    failed$status,
    c("missing_activity_query", "ok_uncovered_zero")
  )
  expect_identical(failed$query_failed, c(TRUE, FALSE))
  expect_match(failed$query_failure_message[[1L]], "synthetic query failure")
})

test_that("activity validates bounds, overlaps, and CRAM exclusion", {
  variant <- data.frame(variant_id = "v", chrom = "chr1", pos1 = 101)
  expect_error(
    extract_v2g_activity(
      variant, data.frame(chrom = "chr1", start0 = 0, end0 = 1, score = 1),
      chrom_lengths = c(chr1 = 100)
    ),
    class = "acv2g_validation_error"
  )
  overlapping <- data.frame(
    chrom = "chr1", start0 = c(90, 95), end0 = c(100, 105), score = 1
  )
  expect_error(
    extract_v2g_activity(
      variant, overlapping, chrom_lengths = c(chr1 = 1000),
      activity_window_bp = 20
    ),
    "overlapping records",
    class = "acv2g_validation_error"
  )
  expect_error(
    extract_v2g_activity(
      variant, "not-present.cram", source_type = "bam",
      chrom_lengths = c(chr1 = 1000)
    ),
    "CRAM",
    class = "acv2g_validation_error"
  )
})

make_contact_pairs <- function(
    tss_pos1 = 20001, tss_id = paste0("t", seq_along(tss_pos1)),
    variant_pos1 = 10001, gene_id = "g1") {
  data.frame(
    locus_id = "L1",
    variant_id = "v1",
    gene_id = gene_id,
    chrom = "chr1",
    pos1 = variant_pos1,
    tss_id = tss_id,
    tss_chrom = "chr1",
    tss_pos1 = tss_pos1,
    stringsAsFactors = FALSE
  )
}

make_contact_table <- function(start1, start2, value, column = "raw") {
  out <- data.frame(
    chrom1 = "chr1",
    start1_0 = start1,
    end1_0 = start1 + 10000,
    chrom2 = "chr1",
    start2_0 = start2,
    end2_0 = start2 + 10000,
    stringsAsFactors = FALSE
  )
  out[[column]] <- value
  out
}

test_that("contact bins are exact, symmetric, and preserve observed zero", {
  pairs <- make_contact_pairs(tss_pos1 = 10001, variant_pos1 = 10000)
  contacts <- make_contact_table(10000, 0, 0)
  result <- extract_v2g_contacts(
    pairs, contacts,
    contact_value = "raw",
    contact_score_col = "raw",
    contact_score_type = "raw",
    diagonal_policy = "observed"
  )
  expect_equal(result$tss_bin_contacts$contact_variant_bin_start0, 0)
  expect_equal(result$tss_bin_contacts$contact_promoter_bin_start0, 10000)
  expect_equal(result$pairs$contact_raw, 0)
  expect_identical(result$pairs$contact_status, "observed_zero")
})

test_that("missing table pixel remains missing rather than becoming zero", {
  pairs <- make_contact_pairs(tss_pos1 = c(20001, 30001))
  contacts <- make_contact_table(10000, 20000, 0)
  result <- extract_v2g_contacts(
    pairs, contacts,
    contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw",
    diagonal_policy = "observed"
  )
  expect_equal(result$tss_bin_contacts$contact_raw, c(0, NA_real_))
  expect_identical(
    result$tss_bin_contacts$contact_status,
    c("observed_zero", "missing_contact_pixel")
  )
  expect_equal(result$pairs$contact_raw, 0)
  expect_equal(result$pairs$contact_tss_bin_n_missing, 1)
  expect_identical(
    result$pairs$contact_aggregation_status,
    "finite_max_with_partial_missing"
  )
})

test_that("balanced default falls back only when exactly one raw layer exists", {
  pairs <- make_contact_pairs()
  contacts <- make_contact_table(10000, 20000, 3)
  expect_warning(
    result <- extract_v2g_contacts(
      pairs, contacts,
      contact_score_col = "raw", contact_score_type = "raw",
      diagonal_policy = "observed"
    ),
    class = "acv2g_contact_raw_fallback"
  )
  expect_true(result$qc$contact_fallback_used)
  expect_identical(result$qc$contact_value_requested, "balanced")
  expect_identical(result$qc$contact_value_resolved, "raw")
  expect_identical(result$qc$contact_fallback_reason, "balanced_layer_absent")

  named <- make_contact_table(10000, 20000, 3, column = "custom")
  expect_error(
    extract_v2g_contacts(
      pairs, named,
      contact_score_col = "custom", contact_score_type = "named_score"
    ),
    "raw fallback cannot be resolved",
    class = "acv2g_validation_error"
  )
})

test_that("ambiguous balanced and missing explicit layers fail with availability", {
  pairs <- make_contact_pairs()
  contacts <- make_contact_table(10000, 20000, 3, column = "KR")
  contacts$VC <- 4
  contacts$count <- 5
  expect_error(
    extract_v2g_contacts(
      pairs, contacts,
      contact_layers = list(balanced = c("KR", "VC"), raw = "count")
    ),
    "Multiple balanced",
    class = "acv2g_validation_error"
  )
  expect_error(
    extract_v2g_contacts(
      pairs, contacts,
      contact_layers = list(balanced = c("KR", "VC"), raw = "count"),
      contact_layer = "ICE"
    ),
    "Available layers: KR, VC, count",
    class = "acv2g_validation_error"
  )
  selected <- extract_v2g_contacts(
    pairs, contacts,
    contact_layers = list(balanced = c("KR", "VC"), raw = "count"),
    contact_layer = "VC", diagonal_policy = "observed"
  )
  expect_equal(selected$pairs$contact_raw, 4)
  expect_false(selected$qc$contact_fallback_used)
})

test_that("diagonal neighbor_max uses two, one, zero, or no finite neighbors", {
  pairs <- make_contact_pairs(tss_pos1 = 10001, variant_pos1 = 10001)
  diagonal <- make_contact_table(10000, 10000, 999)
  left <- make_contact_table(10000, 0, 2)
  right <- make_contact_table(10000, 20000, 5)

  two <- extract_v2g_contacts(
    pairs, rbind(diagonal, left, right), contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw"
  )
  expect_equal(two$pairs$contact_raw, 5)
  expect_equal(two$tss_bin_contacts$contact_diagonal_observed, 999)
  expect_equal(two$tss_bin_contacts$contact_diagonal_replacement, 5)
  expect_identical(
    two$tss_bin_contacts$contact_diagonal_status,
    "neighbor_max_two_finite"
  )

  one <- extract_v2g_contacts(
    pairs, rbind(diagonal, left), contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw"
  )
  expect_equal(one$pairs$contact_raw, 2)
  expect_identical(
    one$tss_bin_contacts$contact_diagonal_status,
    "neighbor_max_one_finite"
  )

  zeros <- rbind(
    diagonal,
    make_contact_table(10000, 0, 0),
    make_contact_table(10000, 20000, 0)
  )
  zero <- extract_v2g_contacts(
    pairs, zeros, contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw"
  )
  expect_equal(zero$pairs$contact_raw, 0)

  none <- extract_v2g_contacts(
    pairs, diagonal, contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw"
  )
  expect_true(is.na(none$pairs$contact_raw))
  expect_identical(none$pairs$contact_status, "missing_contact_all_tss_bins")
})

test_that("observed and keep_na diagonal policies remain auditable", {
  pairs <- make_contact_pairs(tss_pos1 = 10001, variant_pos1 = 10001)
  diagonal <- make_contact_table(10000, 10000, 7)
  observed <- extract_v2g_contacts(
    pairs, diagonal, contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw",
    diagonal_policy = "observed"
  )
  kept <- extract_v2g_contacts(
    pairs, diagonal, contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw",
    diagonal_policy = "keep_na"
  )
  expect_equal(observed$pairs$contact_raw, 7)
  expect_true(is.na(kept$pairs$contact_raw))
  expect_equal(kept$tss_bin_contacts$contact_diagonal_observed, 7)
  expect_identical(kept$tss_bin_contacts$contact_diagonal_status, "kept_na")
})

test_that("TSS coordinates and bins deduplicate before finite max aggregation", {
  pairs <- make_contact_pairs(
    tss_pos1 = c(20001, 20001, 30001, 35001),
    tss_id = c("t1", "t1_alias", "t2", "t3")
  )
  contacts <- rbind(
    make_contact_table(10000, 20000, 4),
    make_contact_table(10000, 30000, 7)
  )
  result <- extract_v2g_contacts(
    pairs, contacts, contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw",
    diagonal_policy = "observed"
  )
  expect_equal(nrow(result$tss_bin_contacts), 2)
  expect_equal(result$pairs$contact_tss_n_eligible, 4)
  expect_equal(result$pairs$contact_tss_bin_n_total, 2)
  expect_equal(result$pairs$contact_raw, 7)
  expect_identical(result$pairs$contact_selected_tss_id, "t2;t3")
  expect_identical(
    result$tss_bin_contacts$contact_tss_ids,
    c("t1;t1_alias", "t2;t3")
  )
  expect_identical(result$tss_bin_contacts$contact_is_max, c(FALSE, TRUE))

  contacts$raw <- 7
  tied <- extract_v2g_contacts(
    pairs, contacts, contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw",
    diagonal_policy = "observed"
  )
  expect_equal(tied$pairs$contact_max_tie_n, 2)
  expect_true(all(tied$tss_bin_contacts$contact_is_max))
  expect_identical(tied$pairs$contact_selected_tss_id, "t1;t1_alias")
})

test_that("TSS-bin detail preserves canonical identifiers without internal keys", {
  pairs <- make_contact_pairs(
    tss_pos1 = c(20001, 30001), tss_id = c("tx1", "tx2")
  )
  pairs$pair_id <- "pair_canonical"
  pairs$variant_anchor_id <- "anchor_canonical"
  pairs$gene_id_normalized <- "ENSG_CANONICAL"
  contacts <- rbind(
    make_contact_table(10000, 20000, 4),
    make_contact_table(10000, 30000, 7)
  )
  result <- extract_v2g_contacts(
    pairs, contacts, contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw",
    diagonal_policy = "observed"
  )
  detail <- result$tss_bin_contacts
  expect_true(all(c(
    "pair_id", "variant_anchor_id", "gene_id_normalized", "tss_id",
    "tss_bin_id"
  ) %in% names(detail)))
  expect_false(".pair_key" %in% names(detail))
  expect_identical(detail$pair_id, rep("pair_canonical", 2L))
  expect_identical(detail$variant_anchor_id, rep("anchor_canonical", 2L))
  expect_identical(detail$gene_id_normalized, rep("ENSG_CANONICAL", 2L))
  expect_identical(detail$tss_id, c("tx1", "tx2"))
  expect_false(anyDuplicated(detail$tss_bin_id) > 0L)
  expect_identical(result$pairs$pair_id, "pair_canonical")

  reordered <- extract_v2g_contacts(
    pairs[2:1, , drop = FALSE], contacts[2:1, , drop = FALSE],
    contact_value = "raw", contact_score_col = "raw",
    contact_score_type = "raw", diagonal_policy = "observed"
  )
  expect_identical(
    detail$tss_bin_id,
    reordered$tss_bin_contacts$tss_bin_id
  )
})

test_that("contact resolution, transforms, and supported policy scope validate", {
  pairs <- make_contact_pairs()
  negative <- make_contact_table(10000, 20000, -0.1)
  expect_error(
    extract_v2g_contacts(
      pairs, negative, contact_value = "raw",
      contact_score_col = "raw", contact_score_type = "raw"
    ),
    class = "acv2g_validation_error"
  )
  identity <- extract_v2g_contacts(
    pairs, negative, contact_value = "raw",
    contact_score_col = "raw", contact_score_type = "raw",
    contact_transform = "identity", diagonal_policy = "observed"
  )
  expect_equal(identity$pairs$contact_raw, -0.1)
  expect_error(
    extract_v2g_contacts(
      pairs, make_contact_table(10000, 20000, 1),
      contact_resolution_bp = 5000, contact_value = "raw",
      contact_score_col = "raw", contact_score_type = "raw"
    ),
    "conflicts",
    class = "acv2g_validation_error"
  )
  expect_error(
    extract_v2g_contacts(
      pairs, make_contact_table(10000, 20000, 1),
      contact_value = "raw", contact_score_col = "raw",
      contact_score_type = "raw", distance_policy = "observed_expected"
    ),
    class = "acv2g_validation_error"
  )
  expect_error(
    extract_v2g_contacts(
      pairs, make_contact_table(10000, 20000, 1),
      contact_value = "raw", contact_score_col = "raw",
      contact_score_type = "raw", table_representation = "sparse"
    ),
    "Sparse called-loop",
    class = "acv2g_validation_error"
  )
})

test_that("native contact paths fail through a clear validation boundary", {
  expect_error(
    extract_v2g_contacts(
      make_contact_pairs(), "not-present.cool", source_type = "cool"
    ),
    "does not exist",
    class = "acv2g_validation_error"
  )
})
