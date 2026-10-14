plot_contract_pairs <- function(profile_id = "p1", call_rule = "none") {
  z <- c(-1, 0, 2, NA_real_)
  out <- data.frame(
    profile_id = profile_id,
    pair_id = paste0(profile_id, "_pair", seq_along(z)),
    variant_anchor_id = paste0(profile_id, "_v", seq_along(z)),
    gene_id_normalized = paste0(profile_id, "_g", seq_along(z)),
    locus_id = rep(c("locus_1", "locus_2"), each = 2),
    chrom = rep(c("chr1", "chr2"), each = 2),
    locus_start0 = rep(c(100L, 200L), each = 2),
    locus_pair_z = z,
    method = "locus_pair_z",
    call_rule = call_rule,
    stringsAsFactors = FALSE
  )
  if (identical(call_rule, "z_threshold")) {
    out$call_threshold <- 1.5
    # Deliberately not recomputed from Z: the plot contract uses stored calls.
    out$link_call <- c(TRUE, FALSE, FALSE, NA)
  } else {
    out$call_threshold <- NA_real_
    out$link_call <- rep(NA, nrow(out))
  }
  out
}

plot_contract_profiles <- function(ids) {
  data.frame(
    profile_id = ids,
    activity_source = paste0(ids, ".bw"),
    contact_source = paste0(ids, ".cool"),
    stringsAsFactors = FALSE
  )
}

test_that("canonical nested profiles resolve one activity-contact pair", {
  profiles <- validate_v2g_profiles(list(
    neuron = list(
      activity = list(path = "neuron.bw", type = "bigwig"),
      contact = list(path = "neuron.cool", type = "cool")
    ),
    glia = list(
      activity = list(path = "glia.bedgraph"),
      contact = list(path = "glia.hic")
    )
  ))
  expect_equal(profiles$profile_id, c("neuron", "glia"))
  expect_equal(profiles$activity_source, c("neuron.bw", "glia.bedgraph"))
  expect_equal(profiles$contact_source, c("neuron.cool", "glia.hic"))
  expect_equal(profiles$activity_type, c("bigwig", NA_character_))
  expect_equal(profiles$contact_type, c("cool", NA_character_))

  flat <- validate_v2g_profiles(data.frame(
    profile_id = "p1", activity_source = "a.bw", contact_source = "c.cool"
  ))
  expect_equal(flat$activity_source, "a.bw")

  duplicated <- structure(
    list(
      list(activity = list(path = "a"), contact = list(path = "c")),
      list(activity = list(path = "b"), contact = list(path = "d"))
    ),
    names = c("same", "same")
  )
  expect_error(validate_v2g_profiles(duplicated), class = "acv2g_validation_error")
  expect_error(
    validate_v2g_profiles(list(
      bad = list(activity = list(path = "a"), contact = list(type = "cool"))
    )),
    class = "acv2g_validation_error"
  )
})

test_that("multi-profile results require explicit profile isolation", {
  pairs <- rbind(plot_contract_pairs("p1"), plot_contract_pairs("p2"))
  result <- new_v2g_result(
    pairs,
    profiles = plot_contract_profiles(c("p1", "p2")),
    config = list(method = list(name = "locus_pair_z", call_rule = "none"))
  )
  expect_error(plot_locus_z(result), class = "acv2g_validation_error")

  plot <- plot_locus_z(result, profile_id = "p2")
  provenance <- attr(plot, "acv2g_plot_provenance")
  expect_equal(provenance$profile_id, "p2")
  expect_equal(provenance$input_pair_count, 4L)
  expect_equal(provenance$locus_count, 2L)
  expect_match(plot$labels$caption, "profile: p2", fixed = TRUE)
  expect_error(
    plot_locus_z(result, profile_id = "absent"),
    class = "acv2g_validation_error"
  )
})

test_that("stored z-threshold calls drive highlights and threshold line", {
  pairs <- plot_contract_pairs(call_rule = "z_threshold")
  result <- new_v2g_result(
    pairs,
    profiles = plot_contract_profiles("p1"),
    config = list(
      method = list(
        name = "locus_pair_z", call_rule = "z_threshold", z_threshold = 1.5
      )
    )
  )
  plot <- plot_locus_z(result)
  provenance <- attr(plot, "acv2g_plot_provenance")
  expect_equal(provenance$call_rule, "z_threshold")
  expect_equal(provenance$call_threshold, 1.5)
  expect_equal(provenance$displayed_threshold, 1.5)
  expect_equal(provenance$highlight_source, "link_call")
  expect_identical(plot$layers[[2L]]$data$.acv2g_highlight, c(TRUE, FALSE, FALSE))
  expect_true(any(vapply(
    plot$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1L)
  )))

  override <- plot_locus_z(result, display_threshold = 2.5)
  override_provenance <- attr(override, "acv2g_plot_provenance")
  expect_equal(override_provenance$call_threshold, 1.5)
  expect_equal(override_provenance$displayed_threshold, 2.5)
  expect_equal(override_provenance$threshold_source, "display_threshold_override")
  expect_equal(override_provenance$highlight_source, "link_call")
  expect_identical(
    override$layers[[2L]]$data$.acv2g_highlight,
    plot$layers[[2L]]$data$.acv2g_highlight
  )
})

test_that("call_rule none has no implicit line or calls", {
  pairs <- plot_contract_pairs()
  plot <- plot_locus_z(pairs)
  provenance <- attr(plot, "acv2g_plot_provenance")
  expect_equal(provenance$call_rule, "none")
  expect_null(provenance$displayed_threshold)
  expect_equal(provenance$highlight_source, "none")
  expect_false(any(vapply(
    plot$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1L)
  )))
  expect_true(all(!plot$layers[[2L]]$data$.acv2g_highlight))

  visual <- plot_locus_z(pairs, display_threshold = 1.5)
  visual_provenance <- attr(visual, "acv2g_plot_provenance")
  expect_equal(visual_provenance$call_rule, "none")
  expect_equal(visual_provenance$highlight_source, "none")
  expect_equal(visual_provenance$displayed_threshold, 1.5)
  expect_equal(visual_provenance$threshold_label, "display threshold = 1.5")
  expect_true(all(!visual$layers[[2L]]$data$.acv2g_highlight))
})

test_that("point modes and offsets are deterministic", {
  pairs <- rbind(
    plot_contract_pairs(),
    transform(
      plot_contract_pairs(),
      pair_id = paste0(pair_id, "_extra"),
      variant_anchor_id = paste0(variant_anchor_id, "_extra")
    )
  )
  all_points <- plot_locus_z(pairs, point_mode = "all", seed = 19L)
  expect_equal(attr(all_points, "acv2g_plot_provenance")$plotted_points, 6L)

  no_points <- plot_locus_z(pairs, point_mode = "none")
  expect_equal(attr(no_points, "acv2g_plot_provenance")$plotted_points, 0L)
  expect_equal(nrow(no_points$layers[[2L]]$data), 0L)

  sampled <- plot_locus_z(
    pairs, point_mode = "sample", max_points_per_locus = 2, seed = 19L
  )
  expect_equal(attr(sampled, "acv2g_plot_provenance")$plotted_points, 4L)
  shuffled <- plot_locus_z(
    pairs[rev(seq_len(nrow(pairs))), ],
    point_mode = "sample", max_points_per_locus = 2, seed = 19L
  )
  selected <- sampled$layers[[2L]]$data
  selected <- selected[order(selected$pair_id), c("pair_id", ".acv2g_plot_x")]
  shuffled_selected <- shuffled$layers[[2L]]$data
  shuffled_selected <- shuffled_selected[
    order(shuffled_selected$pair_id), c("pair_id", ".acv2g_plot_x")
  ]
  rownames(selected) <- NULL
  rownames(shuffled_selected) <- NULL
  expect_equal(selected, shuffled_selected)

  compatibility <- plot_locus_z(pairs, max_points_per_locus = 2, seed = 19L)
  expect_equal(attr(compatibility, "acv2g_plot_provenance")$point_mode, "sample")
})

test_that("plot validates locus metadata and records stable result provenance", {
  pairs <- plot_contract_pairs()
  qc <- data.frame(
    profile_id = rep("p1", 2),
    method = rep("locus_pair_z", 2),
    group_id = c("locus_1", "locus_2"),
    locus_id = c("locus_1", "locus_2"),
    normalization_status = c("ok", "constant_proxy"),
    status = c("ok", "constant_proxy"),
    stringsAsFactors = FALSE
  )
  result <- new_v2g_result(
    pairs,
    profiles = plot_contract_profiles("p1"),
    normalization_qc = qc,
    config = list(method = list(name = "locus_pair_z", call_rule = "none"))
  )
  plot <- plot_locus_z(result)
  provenance <- attr(plot, "acv2g_plot_provenance")
  expect_equal(provenance$nonfinite_pair_count, 1L)
  expect_equal(provenance$nonnormalizable_pair_count, 2L)
  expect_equal(provenance$nonnormalizable_locus_count, 1L)
  expect_match(provenance$source_result_hash, "^[0-9a-f]{64}$")
  expect_match(provenance$source_result_identity, "^sha256:[0-9a-f]{64}$")
  expect_true(nzchar(provenance$package_version))
  expect_match(plot$labels$caption, "missing/non-finite pairs: 1/4", fixed = TRUE)
  expect_match(plot$labels$caption, "non-normalizable loci: 1/2", fixed = TRUE)

  reordered <- new_v2g_result(
    pairs[rev(seq_len(nrow(pairs))), ],
    profiles = plot_contract_profiles("p1"),
    normalization_qc = qc[rev(seq_len(nrow(qc))), ],
    config = list(method = list(name = "locus_pair_z", call_rule = "none"))
  )
  expect_equal(
    attr(plot_locus_z(reordered), "acv2g_plot_provenance")$source_result_hash,
    provenance$source_result_hash
  )

  ambiguous <- pairs
  ambiguous$chrom[[2L]] <- "chr9"
  expect_error(plot_locus_z(ambiguous), class = "acv2g_validation_error")
})

test_that("locus plotting rejects other estimands", {
  pairs <- plot_contract_pairs()
  pairs$method <- "gene_wise_z"
  expect_error(plot_locus_z(pairs), class = "acv2g_validation_error")
  expect_error(
    plot_locus_z(transform(pairs, method = "locus_pair_z"), z_col = "other_z"),
    class = "acv2g_validation_error"
  )
})
