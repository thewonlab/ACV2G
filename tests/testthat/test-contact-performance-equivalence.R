reference_contact_make_tss_bin_seed <- function(pairs, resolution) {
  pair_columns <- ACV2G:::.contact_pair_columns(pairs)
  pair_key <- ACV2G:::.contact_pair_key(pairs)
  variant_bin <- floor((pairs$pos1 - 1) / resolution) * resolution
  promoter_bin <- floor((pairs$tss_pos1 - 1) / resolution) * resolution
  x <- pairs
  bin_key <- paste(pair_key, x$tss_chrom, promoter_bin, sep = "\u001f")
  groups <- split(seq_len(nrow(x)), bin_key)
  rows <- lapply(groups, function(index) {
    first <- index[1L]
    out <- x[first, pair_columns, drop = FALSE]
    out$chrom <- x$chrom[first]
    out$tss_chrom <- x$tss_chrom[first]
    tss_ids <- paste(sort(unique(x$tss_id[index])), collapse = ";")
    out$tss_id <- tss_ids
    out$contact_tss_ids <- tss_ids
    out$contact_tss_positions1 <- paste(
      sort(unique(x$tss_pos1[index])), collapse = ";"
    )
    out$contact_tss_n_in_bin <- length(unique(x$tss_id[index]))
    out$contact_variant_bin_start0 <- variant_bin[first]
    out$contact_promoter_bin_start0 <- promoter_bin[first]
    out$contact_variant_bin_end0 <- variant_bin[first] + resolution
    out$contact_promoter_bin_end0 <- promoter_bin[first] + resolution
    out$contact_distance_bp <- min(abs(x$tss_pos1[index] - x$pos1[index]))
    out$tss_bin_id <- ACV2G:::.stable_row_id(
      out$pair_id, out$tss_chrom, out$contact_promoter_bin_start0,
      resolution, prefix = "tss_bin"
    )
    out
  })
  seed <- do.call(rbind, rows)
  seed$.pair_key <- ACV2G:::.contact_pair_key(seed)
  seed <- seed[order(
    seed$.pair_key, seed$tss_chrom, seed$contact_promoter_bin_start0,
    seed$contact_tss_ids
  ), , drop = FALSE]
  rownames(seed) <- NULL
  seed
}

reference_contact_apply_diagonal <- function(
    seed, direct, prepared, diagonal_policy, resolution) {
  detail <- seed
  is_diagonal <- detail$chrom == detail$tss_chrom &
    detail$contact_variant_bin_start0 ==
    detail$contact_promoter_bin_start0
  detail$contact_is_diagonal <- is_diagonal
  detail$contact_diagonal_observed <- ifelse(
    is_diagonal, direct$value, NA_real_
  )
  detail$contact_diagonal_replacement <- NA_real_
  detail$contact_diagonal_status <- ifelse(
    is_diagonal, NA_character_, "not_diagonal"
  )
  detail$contact_raw <- direct$value
  detail$contact_status <- direct$status

  if (diagonal_policy == "observed") {
    detail$contact_diagonal_status[is_diagonal] <- "used_observed_diagonal"
    return(detail)
  }
  if (diagonal_policy == "keep_na") {
    detail$contact_raw[is_diagonal] <- NA_real_
    detail$contact_status[is_diagonal] <- "diagonal_kept_na"
    detail$contact_diagonal_status[is_diagonal] <- "kept_na"
    return(detail)
  }
  if (!any(is_diagonal)) return(detail)

  index <- which(is_diagonal)
  left_start <- detail$contact_promoter_bin_start0[index] - resolution
  left <- data.frame(
    value = rep(NA_real_, length(index)),
    status = rep("outside_chromosome", length(index)),
    stringsAsFactors = FALSE
  )
  valid_left <- left_start >= 0
  if (any(valid_left)) {
    left[valid_left, ] <- ACV2G:::.contact_lookup_pixels(
      detail$tss_chrom[index[valid_left]],
      detail$contact_variant_bin_start0[index[valid_left]],
      detail$tss_chrom[index[valid_left]], left_start[valid_left], prepared
    )
  }
  right <- ACV2G:::.contact_lookup_pixels(
    detail$tss_chrom[index], detail$contact_variant_bin_start0[index],
    detail$tss_chrom[index],
    detail$contact_promoter_bin_start0[index] + resolution, prepared
  )
  for (k in seq_along(index)) {
    values <- c(left$value[k], right$value[k])
    finite <- is.finite(values)
    i <- index[k]
    if (any(finite)) {
      replacement <- max(values[finite])
      detail$contact_raw[i] <- replacement
      detail$contact_diagonal_replacement[i] <- replacement
      detail$contact_status[i] <- if (replacement == 0) {
        "diagonal_neighbor_zero"
      } else {
        "diagonal_neighbor_observed"
      }
      detail$contact_diagonal_status[i] <- if (sum(finite) == 2L) {
        "neighbor_max_two_finite"
      } else {
        "neighbor_max_one_finite"
      }
    } else {
      detail$contact_raw[i] <- NA_real_
      detail$contact_status[i] <- "missing_diagonal_neighbors"
      detail$contact_diagonal_status[i] <- "neighbor_max_no_finite"
    }
  }
  detail
}

reference_contact_aggregate_tss_bins <- function(
    detail, original_pairs, transform) {
  pair_columns <- ACV2G:::.contact_pair_columns(original_pairs)
  groups <- split(seq_len(nrow(detail)), detail$.pair_key)
  original_pair_key <- ACV2G:::.contact_pair_key(original_pairs)
  eligible_tss <- vapply(
    split(seq_len(nrow(original_pairs)), original_pair_key),
    function(i) length(unique(paste(
      original_pairs$tss_id[i], original_pairs$tss_chrom[i],
      original_pairs$tss_pos1[i], sep = "\u001f"
    ))),
    integer(1L)
  )
  rows <- lapply(names(groups), function(key) {
    index <- groups[[key]]
    values <- detail$contact_raw[index]
    finite <- is.finite(values)
    zero <- finite & values == 0
    missing <- !finite
    out <- detail[index[1L], pair_columns, drop = FALSE]
    out$contact_raw <- if (any(finite)) max(values[finite]) else NA_real_
    out$contact_transformed <- if (transform == "log1p") {
      log1p(out$contact_raw)
    } else {
      out$contact_raw
    }
    out$contact_tss_n_eligible <- unname(eligible_tss[key])
    out$contact_tss_bin_n_total <- length(index)
    out$contact_tss_bin_n_finite <- sum(finite)
    out$contact_tss_bin_n_zero <- sum(zero)
    out$contact_tss_bin_n_missing <- sum(missing)
    out$contact_tss_bin_missing_fraction <- mean(missing)
    variant_bins <- unique(detail$contact_variant_bin_start0[index])
    if (length(variant_bins) != 1L) {
      stop("One variant-gene pair resolved to multiple variant contact bins.")
    }
    out$contact_variant_bin_start0 <- variant_bins
    if (any(finite)) {
      ties <- index[finite & values == out$contact_raw]
      ordered_ties <- ties[order(
        detail$tss_chrom[ties], detail$contact_promoter_bin_start0[ties],
        detail$contact_tss_ids[ties]
      )]
      selected <- ordered_ties[1L]
      out$contact_max_tie_n <- length(ties)
      out$contact_selected_tss_id <- detail$contact_tss_ids[selected]
      out$contact_selected_tss_bin_start0 <-
        detail$contact_promoter_bin_start0[selected]
      out$contact_promoter_bin_start0 <-
        detail$contact_promoter_bin_start0[selected]
      out$contact_distance_bp <- detail$contact_distance_bp[selected]
      out$contact_is_diagonal <- detail$contact_is_diagonal[selected]
      out$contact_diagonal_observed <-
        detail$contact_diagonal_observed[selected]
      out$contact_diagonal_replacement <-
        detail$contact_diagonal_replacement[selected]
      out$contact_diagonal_status <- detail$contact_diagonal_status[selected]
      out$contact_aggregation_status <- if (any(missing)) {
        "finite_max_with_partial_missing"
      } else {
        "finite_max_complete"
      }
      detail$contact_is_max[index] <<- FALSE
      detail$contact_is_max[ties] <<- TRUE
    } else {
      out$contact_max_tie_n <- 0L
      out$contact_selected_tss_id <- NA_character_
      out$contact_selected_tss_bin_start0 <- NA_real_
      out$contact_promoter_bin_start0 <- NA_real_
      out$contact_distance_bp <- NA_real_
      out$contact_is_diagonal <- NA
      out$contact_diagonal_observed <- NA_real_
      out$contact_diagonal_replacement <- NA_real_
      out$contact_diagonal_status <- "all_tss_bins_missing"
      out$contact_aggregation_status <- "all_tss_bins_missing"
      detail$contact_is_max[index] <<- FALSE
    }
    out$contact_status <- if (any(finite)) {
      if (out$contact_raw == 0) "observed_zero" else "observed"
    } else {
      "missing_contact_all_tss_bins"
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "detail_contact_is_max") <- detail$contact_is_max
  out
}

make_equivalence_pairs <- function() {
  pairs <- data.frame(
    pair_id = c(
      "pair_d", "pair_b", "pair_a", "pair_c", "pair_a", "pair_b",
      "pair_d", "pair_a", "pair_d", "pair_a"
    ),
    locus_id = c(
      "L4", "L2", "L1", "L3", "L1", "L2", "L4", "L1", "L4", "L1"
    ),
    variant_anchor_id = c(
      "v4", "v2", "v1", "v3", "v1", "v2", "v4", "v1", "v4", "v1"
    ),
    variant_id = c(
      "v4", "v2", "v1", "v3", "v1", "v2", "v4", "v1", "v4", "v1"
    ),
    gene_id_normalized = c(
      "g4", "g2", "g1", "g3", "g1", "g2", "g4", "g1", "g4", "g1"
    ),
    gene_id = c(
      "g4", "g2", "g1", "g3", "g1", "g2", "g4", "g1", "g4", "g1"
    ),
    chrom = "chr1",
    pos1 = c(90001, 30001, 50001, 10001, 50001, 30001,
             90001, 50001, 90001, 50001),
    tss_id = c("d2", "z", "a2", "c1", "a1", "a", "d1", "a3", "d3", "a1_alias"),
    tss_chrom = "chr1",
    tss_pos1 = c(110001, 25000, 71001, 10001, 61001, 21001,
                 100001, 75000, 115000, 61001),
    stringsAsFactors = FALSE
  )
  ACV2G:::.contact_validate_pairs(pairs)
}

test_that("vectorized TSS-bin seeds are identical to the frozen implementation", {
  pairs <- make_equivalence_pairs()
  expected <- reference_contact_make_tss_bin_seed(pairs, 10000)
  observed <- ACV2G:::.contact_make_tss_bin_seed(pairs, 10000)

  expect_identical(observed, expected)
  expect_identical(
    ACV2G:::.contact_vector_stable_row_id(
      observed$pair_id, observed$tss_chrom,
      observed$contact_promoter_bin_start0, 10000,
      prefix = "tss_bin"
    ),
    ACV2G:::.stable_row_id(
      observed$pair_id, observed$tss_chrom,
      observed$contact_promoter_bin_start0, 10000,
      prefix = "tss_bin"
    )
  )
})

make_diagonal_fixture <- function() {
  seed <- data.frame(
    pair_id = c("two", "one_zero", "none", "edge", "off_diagonal"),
    chrom = "chr1",
    tss_chrom = "chr1",
    contact_variant_bin_start0 = c(10000, 40000, 70000, 0, 100000),
    contact_promoter_bin_start0 = c(10000, 40000, 70000, 0, 110000),
    stringsAsFactors = FALSE
  )
  direct <- data.frame(
    value = c(999, 888, 777, 666, 9),
    status = rep("observed", 5L),
    stringsAsFactors = FALSE
  )
  starts1 <- c(10000, 10000, 40000, 40000, 0)
  starts2 <- c(0, 20000, 30000, 50000, 10000)
  values <- c(2, 5, 0, NA_real_, 4)
  pixels <- data.frame(
    .pixel_key = ACV2G:::.contact_pixel_key(
      "chr1", starts1, "chr1", starts2
    ),
    .selected_contact_value = values,
    stringsAsFactors = FALSE
  )
  prepared <- list(
    pixels = pixels,
    available_chromosomes = "chr1",
    chrom_lengths = c(chr1 = 200000),
    implicit_zero = FALSE,
    selected_type = "raw",
    bin_weights = NULL
  )
  list(seed = seed, direct = direct, prepared = prepared)
}

test_that("vectorized diagonal resolution is reference-identical", {
  fixture <- make_diagonal_fixture()
  for (policy in c("neighbor_max", "observed", "keep_na")) {
    expected <- reference_contact_apply_diagonal(
      fixture$seed, fixture$direct, fixture$prepared, policy, 10000
    )
    observed <- ACV2G:::.contact_apply_diagonal(
      fixture$seed, fixture$direct, fixture$prepared, policy, 10000
    )
    expect_identical(observed, expected, info = policy)
  }
})

test_that("vectorized pair aggregation preserves ties, zeros, and missingness", {
  pairs <- make_equivalence_pairs()
  detail <- reference_contact_make_tss_bin_seed(pairs, 10000)
  value_key <- paste(detail$.pair_key, detail$contact_promoter_bin_start0)
  values <- c(
    "pair_a 60000" = 7,
    "pair_a 70000" = 7,
    "pair_b 20000" = 0,
    "pair_c 10000" = NA_real_,
    "pair_d 100000" = 4,
    "pair_d 110000" = NA_real_
  )
  detail$contact_raw <- unname(values[value_key])
  detail$contact_is_diagonal <-
    detail$contact_variant_bin_start0 == detail$contact_promoter_bin_start0
  detail$contact_diagonal_observed <- ifelse(
    detail$contact_is_diagonal, detail$contact_raw, NA_real_
  )
  detail$contact_diagonal_replacement <- NA_real_
  detail$contact_diagonal_status <- ifelse(
    detail$contact_is_diagonal, "used_observed_diagonal", "not_diagonal"
  )
  detail <- detail[c(6, 2, 1, 5, 3, 4), , drop = FALSE]
  rownames(detail) <- NULL

  for (transform in c("identity", "log1p")) {
    expected <- reference_contact_aggregate_tss_bins(
      detail, pairs, transform
    )
    observed <- ACV2G:::.contact_aggregate_tss_bins(
      detail, pairs, transform
    )
    expect_identical(observed, expected, info = transform)
  }
})

test_that("implicit structural zero stays distinct from a missing balance weight", {
  pairs <- ACV2G:::.contact_validate_pairs(data.frame(
    pair_id = c("zero", "missing_weight"),
    locus_id = c("L1", "L2"),
    variant_id = c("v1", "v2"),
    gene_id = c("g1", "g2"),
    chrom = "chr1",
    pos1 = c(10001, 30001),
    tss_id = c("t1", "t2"),
    tss_chrom = "chr1",
    tss_pos1 = c(20001, 40001),
    stringsAsFactors = FALSE
  ))
  seed <- ACV2G:::.contact_make_tss_bin_seed(pairs, 10000)
  prepared <- list(
    pixels = data.frame(
      .pixel_key = character(),
      .selected_contact_value = numeric(),
      stringsAsFactors = FALSE
    ),
    available_chromosomes = "chr1",
    chrom_lengths = c(chr1 = 100000),
    implicit_zero = TRUE,
    selected_type = "balanced",
    bin_weights = data.frame(
      chrom = "chr1",
      start0 = c(10000, 20000, 30000, 40000),
      weight = c(1, 1, NA_real_, 1),
      stringsAsFactors = FALSE
    )
  )
  direct <- ACV2G:::.contact_lookup_pixels(
    seed$chrom, seed$contact_variant_bin_start0,
    seed$tss_chrom, seed$contact_promoter_bin_start0, prepared
  )
  expect_identical(
    direct$status,
    c("missing_balance_weight", "implicit_structural_zero")
  )
  expect_equal(direct$value, c(NA_real_, 0))

  detail <- ACV2G:::.contact_apply_diagonal(
    seed, direct, prepared, "observed", 10000
  )
  expected <- reference_contact_aggregate_tss_bins(
    detail, pairs, "identity"
  )
  observed <- ACV2G:::.contact_aggregate_tss_bins(
    detail, pairs, "identity"
  )
  expect_identical(observed, expected)
})
