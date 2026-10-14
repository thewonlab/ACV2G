.acv2g_plot_scalar_text <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    .acv2g_validation_error("{.arg {arg}} must be one non-missing, non-empty character value.")
  }
  x
}

.acv2g_plot_profile_ids <- function(data, x = NULL) {
  ids <- character()
  if ("profile_id" %in% names(data)) {
    pair_ids <- as.character(data$profile_id)
    if (anyNA(pair_ids) || any(!nzchar(trimws(pair_ids)))) {
      .acv2g_validation_error("Pair data contain missing or empty {.field profile_id} values.")
    }
    ids <- unique(pair_ids)
  }
  if (inherits(x, "v2g_result") && "profiles" %in% names(x$tables)) {
    profiles <- v2g_table(x, "profiles")
    if (nrow(profiles) && "profile_id" %in% names(profiles)) {
      result_ids <- as.character(profiles$profile_id)
      result_ids <- result_ids[!is.na(result_ids) & nzchar(trimws(result_ids))]
      ids <- unique(c(ids, result_ids))
    }
  }
  ids
}

.acv2g_select_plot_profile <- function(data, profile_id = NULL, x = NULL) {
  available <- .acv2g_plot_profile_ids(data, x)
  if (!is.null(profile_id)) {
    profile_id <- .acv2g_plot_scalar_text(profile_id, "profile_id")
    if (!length(available) || !profile_id %in% available) {
      .acv2g_validation_error(
        paste0(
          "Unknown plot profile ", profile_id, ". Available profiles: ",
          if (length(available)) paste(available, collapse = ", ") else "none",
          "."
        )
      )
    }
  } else if (length(available) > 1L) {
    .acv2g_validation_error(
      paste0(
        "The result contains multiple profiles; {.arg profile_id} is required. ",
        "Available profiles: ", paste(available, collapse = ", "), "."
      )
    )
  } else if (length(available) == 1L) {
    profile_id <- available[[1L]]
  }

  if (!is.null(profile_id)) {
    if (!"profile_id" %in% names(data)) {
      .acv2g_validation_error(
        "Pair data must contain {.field profile_id} when {.arg profile_id} is supplied."
      )
    }
    data <- data[as.character(data$profile_id) == profile_id, , drop = FALSE]
  }
  if (!nrow(data)) {
    .acv2g_validation_error("The selected profile contains no variant-gene pairs to plot.")
  }
  list(data = data, profile_id = profile_id, available = available)
}

.acv2g_validate_plot_method <- function(data, x, z_col) {
  if (!identical(z_col, "locus_pair_z")) {
    .acv2g_validation_error(
      "{.fn plot_locus_z} only plots the canonical {.field locus_pair_z} column."
    )
  }
  methods <- character()
  if ("method" %in% names(data)) {
    methods <- unique(stats::na.omit(as.character(data$method)))
    methods <- methods[nzchar(trimws(methods))]
  }
  configured <- NULL
  if (inherits(x, "v2g_result") && is.list(x$config$method)) {
    configured <- x$config$method$name %||% NULL
    if (!is.null(configured)) {
      configured <- .acv2g_plot_scalar_text(configured, "x$config$method$name")
      methods <- unique(c(methods, configured))
    }
  }
  if (length(methods) && any(methods != "locus_pair_z")) {
    .acv2g_validation_error(
      paste0(
        "plot_locus_z() only accepts locus_pair_z results; observed method(s): ",
        paste(methods, collapse = ", "), "."
      )
    )
  }
  invisible("locus_pair_z")
}

.acv2g_validate_locus_metadata <- function(
    data,
    locus_col,
    chromosome_col,
    locus_start_col) {
  locus <- as.character(data[[locus_col]])
  chromosome <- as.character(data[[chromosome_col]])
  if (anyNA(locus) || any(!nzchar(trimws(locus)))) {
    .acv2g_validation_error("Locus identifiers must be non-missing and non-empty.")
  }
  if (anyNA(chromosome) || any(!nzchar(trimws(chromosome)))) {
    .acv2g_validation_error("Chromosome labels must be non-missing and non-empty.")
  }
  data[[locus_col]] <- locus
  data[[chromosome_col]] <- chromosome

  fields <- chromosome_col
  if (locus_start_col %in% names(data)) {
    start <- data[[locus_start_col]]
    if (!is.numeric(start) || anyNA(start) || any(!is.finite(start)) ||
        any(start < 0) || any(start != floor(start))) {
      .acv2g_validation_error(
        "{.field {locus_start_col}} must contain non-negative integer coordinates."
      )
    }
    fields <- c(fields, locus_start_col)
  }

  groups <- split(seq_len(nrow(data)), locus, drop = TRUE)
  ambiguous <- vapply(groups, function(index) {
    any(vapply(fields, function(field) {
      length(unique(data[[field]][index])) != 1L
    }, logical(1L)))
  }, logical(1L))
  if (any(ambiguous)) {
    examples <- paste(utils::head(names(groups)[ambiguous], 5L), collapse = ", ")
    .acv2g_validation_error(
      paste0(
        "Each locus must resolve to exactly one chromosome and locus start. ",
        "Ambiguous locus metadata: ", examples, "."
      )
    )
  }

  first <- vapply(groups, `[[`, integer(1L), 1L)
  columns <- c(locus_col, chromosome_col)
  if (locus_start_col %in% names(data)) columns <- c(columns, locus_start_col)
  info <- data[first, columns, drop = FALSE]
  if (locus_start_col %in% names(info)) {
    info <- info[order(
      .chromosome_order(info[[chromosome_col]]),
      info[[locus_start_col]], info[[locus_col]]
    ), , drop = FALSE]
  } else {
    info <- info[order(
      .chromosome_order(info[[chromosome_col]]), info[[locus_col]]
    ), , drop = FALSE]
  }
  rownames(info) <- NULL
  list(data = data, locus_info = info)
}

.acv2g_unique_plot_value <- function(data, column, type = c("character", "number")) {
  type <- match.arg(type)
  if (!column %in% names(data)) return(NULL)
  value <- data[[column]]
  if (type == "character") {
    value <- as.character(value)
    value <- value[!is.na(value) & nzchar(trimws(value))]
  } else {
    if (!is.numeric(value)) {
      .acv2g_validation_error("{.field {column}} must be numeric.")
    }
    if (any(!is.na(value) & !is.finite(value))) {
      .acv2g_validation_error("{.field {column}} cannot contain infinite values.")
    }
    value <- value[is.finite(value)]
  }
  value <- unique(value)
  if (length(value) > 1L) {
    .acv2g_validation_error(
      "The selected profile must have one consistent {.field {column}} value."
    )
  }
  if (length(value)) value[[1L]] else NULL
}

.acv2g_plot_call_contract <- function(data, x, display_threshold) {
  pair_rule <- .acv2g_unique_plot_value(data, "call_rule", "character")
  config_rule <- NULL
  config_threshold <- NULL
  if (inherits(x, "v2g_result") && is.list(x$config$method)) {
    config_rule <- x$config$method$call_rule %||% NULL
    config_threshold <- x$config$method$z_threshold %||% NULL
    if (!is.null(config_rule)) {
      config_rule <- .acv2g_plot_scalar_text(config_rule, "x$config$method$call_rule")
    }
  }
  if (!is.null(pair_rule) && !is.null(config_rule) && !identical(pair_rule, config_rule)) {
    .acv2g_validation_error(
      "Pair-level {.field call_rule} conflicts with the resolved result configuration."
    )
  }
  rule <- pair_rule %||% config_rule %||% "none"
  if (!rule %in% c("none", "z_threshold")) {
    .acv2g_validation_error(
      "Plotting supports {.field call_rule} values {.val none} and {.val z_threshold} only."
    )
  }

  threshold <- NULL
  highlight <- rep(FALSE, nrow(data))
  highlight_source <- "none"
  if (identical(rule, "z_threshold")) {
    threshold <- .acv2g_unique_plot_value(data, "call_threshold", "number")
    if (is.null(threshold) && !is.null(config_threshold)) {
      threshold <- .assert_scalar_number(config_threshold, "x$config$method$z_threshold")
    }
    if (is.null(threshold)) {
      .acv2g_validation_error(
        "A z_threshold result must retain one finite {.field call_threshold}."
      )
    }
    if (!"link_call" %in% names(data) || !is.logical(data$link_call)) {
      .acv2g_validation_error(
        "A z_threshold result must retain logical {.field link_call} values."
      )
    }
    finite_z <- is.finite(data$locus_pair_z)
    if (any(finite_z & is.na(data$link_call))) {
      .acv2g_validation_error(
        "Finite locus-pair Z-scores cannot have missing {.field link_call} under {.val z_threshold}."
      )
    }
    highlight <- data$link_call %in% TRUE
    highlight_source <- "link_call"
  } else if (!is.null(pair_rule) && "link_call" %in% names(data) &&
             any(!is.na(data$link_call))) {
    .acv2g_validation_error(
      "Rows declaring {.val call_rule = none} cannot contain non-missing {.field link_call} values."
    )
  }

  line_threshold <- threshold
  threshold_source <- if (is.null(threshold)) "none" else "result_call_threshold"
  if (!is.null(display_threshold)) {
    line_threshold <- display_threshold
    threshold_source <- if (is.null(threshold)) {
      "display_threshold"
    } else {
      "display_threshold_override"
    }
  }
  list(
    rule = rule,
    call_threshold = threshold,
    line_threshold = line_threshold,
    threshold_source = threshold_source,
    highlight = highlight,
    highlight_source = highlight_source
  )
}

.acv2g_plot_row_identity <- function(data, locus_col, chromosome_col, z_col) {
  preferred <- intersect(
    c(
      "profile_id", "pair_id", "variant_anchor_id", "variant_id",
      "gene_id_normalized", "gene_id", locus_col, chromosome_col, z_col
    ),
    names(data)
  )
  encoded <- lapply(data[preferred], function(value) {
    value <- if (is.numeric(value)) {
      ifelse(is.na(value), "<NA>", format(value, digits = 17, scientific = TRUE))
    } else {
      value <- as.character(value)
      ifelse(is.na(value), "<NA>", value)
    }
    paste0(nchar(value, type = "bytes"), ":", value)
  })
  identity <- do.call(paste, c(encoded, sep = "\u001f"))
  occurrence <- ave(seq_along(identity), identity, FUN = seq_along)
  paste(identity, occurrence, sep = "\u001e")
}

.acv2g_hash_int <- function(key, seed) {
  vapply(key, function(value) {
    as.double(digest::digest2int(paste(seed, value, sep = "\u001f")))
  }, numeric(1L), USE.NAMES = FALSE)
}

.deterministic_within_group_sample <- function(
    data,
    group_col,
    max_n,
    seed,
    identity_col = ".acv2g_identity") {
  .assert_scalar_number(max_n, "max_points_per_locus", lower = 1, integer = TRUE)
  .assert_scalar_number(seed, "seed", integer = TRUE)
  groups <- split(seq_len(nrow(data)), data[[group_col]], drop = TRUE)
  selected <- unlist(lapply(names(groups), function(group) {
    index <- groups[[group]]
    if (length(index) <= max_n) return(index)
    score <- .acv2g_hash_int(
      paste(group, data[[identity_col]][index], sep = "\u001f"), seed
    )
    index[order(score, data[[identity_col]][index], method = "radix")][seq_len(max_n)]
  }), use.names = FALSE)
  data[sort(selected), , drop = FALSE]
}

.acv2g_canonical_plot_table <- function(data, columns) {
  columns <- sort(intersect(columns, names(data)))
  data <- data[, columns, drop = FALSE]
  for (column in names(data)) {
    if (is.factor(data[[column]])) data[[column]] <- as.character(data[[column]])
  }
  if (nrow(data) && ncol(data)) {
    keys <- lapply(data, function(value) {
      value <- if (is.numeric(value)) {
        ifelse(is.na(value), "<NA>", format(value, digits = 17, scientific = TRUE))
      } else {
        value <- as.character(value)
        ifelse(is.na(value), "<NA>", value)
      }
      paste0(nchar(value, type = "bytes"), ":", value)
    })
    row_key <- do.call(paste, c(keys, sep = "\u001f"))
    data <- data[order(row_key, method = "radix"), , drop = FALSE]
  }
  rownames(data) <- NULL
  data
}

.acv2g_plot_result_hash <- function(x, pairs, normalization_qc) {
  if (!inherits(x, "v2g_result")) return(NULL)
  pair_columns <- c(
    "schema_version", "profile_id", "pair_id", "variant_anchor_id",
    "variant_id", "gene_id_normalized", "gene_id", "locus_id", "chrom",
    "locus_start0", "method", "proxy_score", "locus_pair_z", "call_rule",
    "call_threshold", "link_call", "normalization_status"
  )
  qc_columns <- c(
    "schema_version", "profile_id", "method", "locus_id", "group_id",
    "normalization_status", "status", "normalization_n_total",
    "normalization_n_valid", "normalization_missing_fraction"
  )
  content <- list(
    schema_version = x$schema_version %||% .acv2g_schema_version,
    pairs = .acv2g_canonical_plot_table(pairs, pair_columns),
    normalization_qc = .acv2g_canonical_plot_table(normalization_qc, qc_columns)
  )
  digest::digest(content, algo = "sha256", serialize = TRUE)
}

.acv2g_plot_normalization_counts <- function(
    data,
    normalization_qc,
    locus_col,
    profile_id) {
  loci <- unique(as.character(data[[locus_col]]))
  bad <- character()
  qc <- normalization_qc
  if (nrow(qc)) {
    if (!is.null(profile_id) && "profile_id" %in% names(qc)) {
      qc <- qc[as.character(qc$profile_id) == profile_id, , drop = FALSE]
    }
    if ("method" %in% names(qc)) {
      qc <- qc[is.na(qc$method) | as.character(qc$method) == "locus_pair_z", , drop = FALSE]
    }
    qc_locus_col <- if (locus_col %in% names(qc)) locus_col else
      if ("locus_id" %in% names(qc)) "locus_id" else NULL
    status_col <- if ("normalization_status" %in% names(qc)) {
      "normalization_status"
    } else if ("status" %in% names(qc)) {
      "status"
    } else {
      NULL
    }
    if (!is.null(qc_locus_col) && !is.null(status_col)) {
      status <- as.character(qc[[status_col]])
      bad <- as.character(qc[[qc_locus_col]][is.na(status) | status != "ok"])
    }
  }
  if (!length(bad) && "normalization_status" %in% names(data)) {
    status <- as.character(data$normalization_status)
    bad <- as.character(data[[locus_col]][is.na(status) | status != "ok"])
  }
  all_nonfinite <- vapply(split(data$locus_pair_z, data[[locus_col]]), function(z) {
    !any(is.finite(z))
  }, logical(1L))
  bad <- unique(c(bad, names(all_nonfinite)[all_nonfinite]))
  bad <- intersect(loci, bad)
  list(
    pair_count = sum(as.character(data[[locus_col]]) %in% bad),
    locus_count = length(bad),
    locus_ids = bad
  )
}

.acv2g_chromosome_runs <- function(strip) {
  run <- rle(as.character(strip$.acv2g_chrom))
  end <- cumsum(run$lengths)
  start <- end - run$lengths + 1L
  data.frame(
    .acv2g_x = (start + end) / 2,
    .acv2g_label = run$values,
    stringsAsFactors = FALSE
  )
}

#' Plot locus-level pair Z-scores
#'
#' `plot_locus_z()` only plots the default `locus_pair_z` estimand. A
#' multi-profile result must be reduced explicitly with `profile_id`; profiles
#' are never pooled. Stored `z_threshold` calls are highlighted from
#' `link_call`, while `display_threshold` changes the displayed reference line
#' only and never changes stored calls.
#'
#' @param x A `v2g_result` or pair data frame.
#' @param profile_id Profile to plot. Required for multi-profile results and
#'   inferred only when exactly one profile is available.
#' @param z_col,locus_col,chromosome_col,locus_start_col Column mappings.
#'   `z_col` must remain `"locus_pair_z"` in this release.
#' @param display_threshold Optional plotting-only reference threshold. For a
#'   result with stored threshold calls, this overrides the displayed line but
#'   does not replace `link_call` highlighting.
#' @param title Optional user-supplied title; `NULL` produces no title.
#' @param point_mode One of `"all"`, `"sample"`, or `"none"`. The default is
#'   `"all"`. For compatibility, supplying a finite `max_points_per_locus`
#'   without `point_mode` selects `"sample"`.
#' @param max_points_per_locus Positive finite per-locus limit used only by
#'   sampled points.
#' @param seed Integer used for deterministic sampling and horizontal offsets.
#' @param point_alpha Point alpha.
#'
#' @return A `ggplot` object with `acv2g_plot_provenance` metadata.
#' @export
plot_locus_z <- function(
    x,
    profile_id = NULL,
    z_col = "locus_pair_z",
    locus_col = "locus_id",
    chromosome_col = "chrom",
    locus_start_col = "locus_start0",
    display_threshold = NULL,
    title = NULL,
    point_mode = c("all", "sample", "none"),
    max_points_per_locus = Inf,
    seed = 1L,
    point_alpha = 0.35) {
  is_result <- inherits(x, "v2g_result")
  all_data <- if (is_result) v2g_table(x, "pairs") else
    .as_plain_data_frame(x, "x")
  normalization_qc <- if (is_result && "normalization_qc" %in% names(x$tables)) {
    v2g_table(x, "normalization_qc")
  } else {
    data.frame()
  }
  source_hash <- .acv2g_plot_result_hash(x, all_data, normalization_qc)

  .require_columns(all_data, c(z_col, locus_col, chromosome_col), "pair data")
  selected <- .acv2g_select_plot_profile(all_data, profile_id, x)
  data <- selected$data
  profile_id <- selected$profile_id
  .acv2g_validate_plot_method(data, x, z_col)
  if (!is.numeric(data[[z_col]])) {
    .acv2g_validation_error("Z-score column {.field {z_col}} must be numeric.")
  }
  if (!is.null(display_threshold)) {
    display_threshold <- .assert_scalar_number(
      display_threshold, "display_threshold"
    )
  }
  if (!is.null(title)) title <- .acv2g_plot_scalar_text(title, "title")
  .assert_scalar_number(seed, "seed", integer = TRUE)
  .assert_scalar_number(point_alpha, "point_alpha", lower = 0, upper = 1)

  point_mode_was_missing <- missing(point_mode)
  point_mode <- match.arg(point_mode)
  if (point_mode_was_missing && is.finite(max_points_per_locus)) {
    point_mode <- "sample"
  }
  if (identical(point_mode, "sample")) {
    if (!is.finite(max_points_per_locus)) {
      .acv2g_validation_error(
        "{.arg max_points_per_locus} must be finite when {.arg point_mode = 'sample'}."
      )
    }
    max_points_per_locus <- .assert_scalar_number(
      max_points_per_locus, "max_points_per_locus", lower = 1, integer = TRUE
    )
  } else if (!is.numeric(max_points_per_locus) ||
             length(max_points_per_locus) != 1L ||
             is.na(max_points_per_locus) || max_points_per_locus <= 0 ||
             (is.finite(max_points_per_locus) &&
              max_points_per_locus != floor(max_points_per_locus))) {
    .acv2g_validation_error(
      "{.arg max_points_per_locus} must be a positive integer or positive infinity."
    )
  }

  metadata <- .acv2g_validate_locus_metadata(
    data, locus_col, chromosome_col, locus_start_col
  )
  data <- metadata$data
  locus_info <- metadata$locus_info
  data$locus_pair_z <- data[[z_col]]
  call_contract <- .acv2g_plot_call_contract(data, x, display_threshold)

  locus_levels <- as.character(locus_info[[locus_col]])
  data$.acv2g_locus <- factor(data[[locus_col]], levels = locus_levels)
  data$.acv2g_x <- match(data[[locus_col]], locus_levels)
  data$.acv2g_z <- data[[z_col]]
  data$.acv2g_identity <- .acv2g_plot_row_identity(
    data, locus_col, chromosome_col, z_col
  )
  data$.acv2g_highlight <- call_contract$highlight

  finite <- data[is.finite(data$.acv2g_z), , drop = FALSE]
  if (identical(point_mode, "sample")) {
    plotted <- .deterministic_within_group_sample(
      finite, ".acv2g_locus", max_points_per_locus, seed
    )
  } else if (identical(point_mode, "none")) {
    plotted <- finite[FALSE, , drop = FALSE]
  } else {
    plotted <- finite
  }
  offset_score <- .acv2g_hash_int(plotted$.acv2g_identity, seed)
  offset_unit <- (offset_score + 2147483648) / 4294967296
  plotted$.acv2g_plot_x <- plotted$.acv2g_x + (offset_unit - 0.5) * 0.36

  chromosome_levels <- unique(as.character(locus_info[[chromosome_col]]))
  palette <- grDevices::hcl.colors(max(3L, length(chromosome_levels)), "Dark 3")
  chrom_colors <- setNames(palette[seq_along(chromosome_levels)], chromosome_levels)

  range_values <- c(
    data$.acv2g_z[is.finite(data$.acv2g_z)], call_contract$line_threshold
  )
  range_values <- range_values[is.finite(range_values)]
  z_min <- if (length(range_values)) min(range_values) else -1
  z_max <- if (length(range_values)) max(range_values) else 1
  span <- max(z_max - z_min, 1)
  strip_y <- z_min - 0.11 * span
  strip_height <- 0.075 * span
  strip <- locus_info
  strip$.acv2g_x <- seq_len(nrow(strip))
  strip$.acv2g_chrom <- as.character(strip[[chromosome_col]])
  strip$.acv2g_y <- strip_y
  strip_labels <- .acv2g_chromosome_runs(strip)
  strip_labels$.acv2g_y <- strip_y

  input_pair_count <- nrow(data)
  finite_pair_count <- sum(is.finite(data$.acv2g_z))
  nonfinite_pair_count <- input_pair_count - finite_pair_count
  locus_count <- length(locus_levels)
  normalization_counts <- .acv2g_plot_normalization_counts(
    data, normalization_qc, locus_col, profile_id
  )
  nonnormalizable_pair_count <- normalization_counts$pair_count
  nonnormalizable_locus_count <- normalization_counts$locus_count
  point_text <- switch(
    point_mode,
    all = sprintf("points: all %d finite (offset seed %d)", nrow(plotted), seed),
    sample = sprintf(
      "points: sampled %d/%d finite (max %d/locus; seed %d)",
      nrow(plotted), finite_pair_count, max_points_per_locus, seed
    ),
    none = "points: none"
  )
  caption_parts <- c(
    if (!is.null(profile_id)) paste0("profile: ", profile_id) else NULL,
    sprintf("missing/non-finite pairs: %d/%d", nonfinite_pair_count, input_pair_count),
    sprintf(
      "non-normalizable pairs: %d/%d",
      nonnormalizable_pair_count, input_pair_count
    ),
    sprintf(
      "non-normalizable loci: %d/%d",
      nonnormalizable_locus_count, locus_count
    ),
    point_text
  )
  caption <- paste(caption_parts, collapse = " | ")

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$.acv2g_x, y = .data$.acv2g_z,
      group = .data$.acv2g_locus
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.65, outlier.shape = NA, linewidth = 0.25,
      fill = "grey85", na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = plotted,
      ggplot2::aes(
        x = .data$.acv2g_plot_x, y = .data$.acv2g_z,
        color = .data$.acv2g_highlight
      ),
      inherit.aes = FALSE, size = 0.65, alpha = point_alpha,
      na.rm = TRUE, show.legend = FALSE
    ) +
    ggplot2::geom_tile(
      data = strip,
      ggplot2::aes(
        x = .data$.acv2g_x, y = .data$.acv2g_y,
        fill = .data$.acv2g_chrom
      ),
      inherit.aes = FALSE, height = strip_height, width = 0.98,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = strip_labels,
      ggplot2::aes(
        x = .data$.acv2g_x, y = .data$.acv2g_y,
        label = .data$.acv2g_label
      ),
      inherit.aes = FALSE, color = "white", fontface = "bold", size = 2.5
    ) +
    ggplot2::scale_fill_manual(values = chrom_colors) +
    ggplot2::scale_x_continuous(
      limits = c(0.45, length(locus_levels) + 0.55), expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(strip_y - 0.7 * strip_height, z_max + 0.05 * span),
      clip = "off"
    ) +
    ggplot2::labs(x = NULL, y = "Z-score", title = title, caption = caption) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, color = "grey30")
    )
  if (nrow(plotted)) {
    p <- p + ggplot2::scale_color_manual(
      values = c(`FALSE` = "grey40", `TRUE` = "#D44835"), drop = FALSE
    )
  }
  if (!is.null(call_contract$line_threshold)) {
    threshold_label <- paste0(
      if (grepl("^display_threshold", call_contract$threshold_source)) {
        "display threshold"
      } else {
        "call threshold"
      },
      " = ", format(call_contract$line_threshold, digits = 6, trim = TRUE)
    )
    p <- p + ggplot2::geom_hline(
      yintercept = call_contract$line_threshold,
      linetype = "dashed", color = "#B22222", linewidth = 0.4
    ) + ggplot2::annotate(
      "text", x = Inf, y = call_contract$line_threshold,
      label = threshold_label, hjust = 1.05, vjust = -0.45,
      color = "#B22222", size = 3
    )
  } else {
    threshold_label <- NULL
  }

  package_version <- tryCatch(
    as.character(utils::packageVersion("ACV2G")),
    error = function(e) "0.0.0.9000"
  )
  attr(p, "acv2g_plot_provenance") <- list(
    schema_version = .acv2g_schema_version,
    package_version = package_version,
    source_result_identity = if (is.null(source_hash)) NULL else
      paste0("sha256:", source_hash),
    source_result_id = if (is.null(source_hash)) NULL else
      paste0("sha256:", source_hash),
    source_result_hash = source_hash,
    profile_id = profile_id,
    method = "locus_pair_z",
    z_col = z_col,
    call_rule = call_contract$rule,
    call_threshold = call_contract$call_threshold,
    resolved_threshold = call_contract$line_threshold,
    displayed_threshold = call_contract$line_threshold,
    display_threshold = display_threshold,
    threshold_label = threshold_label,
    threshold_source = call_contract$threshold_source,
    highlight_source = call_contract$highlight_source,
    point_mode = point_mode,
    max_points_per_locus = if (identical(point_mode, "sample"))
      max_points_per_locus else Inf,
    seed = if (identical(point_mode, "none")) NULL else seed,
    deterministic_offset_width = 0.36,
    input_rows = input_pair_count,
    input_pair_count = input_pair_count,
    finite_pair_count = finite_pair_count,
    nonfinite_pair_count = nonfinite_pair_count,
    nonnormalizable_pair_count = nonnormalizable_pair_count,
    locus_count = locus_count,
    nonnormalizable_locus_count = nonnormalizable_locus_count,
    nonnormalizable_locus_ids = normalization_counts$locus_ids,
    ordering = list(
      chromosome = "natural",
      within_chromosome = if (locus_start_col %in% names(data))
        c(locus_start_col, locus_col) else locus_col,
      locus_key = c(if (!is.null(profile_id)) "profile_id", locus_col)
    ),
    chromosome_order = unique(as.character(locus_info[[chromosome_col]])),
    locus_order = locus_levels,
    plotted_points = nrow(plotted),
    caption = caption
  )
  p
}

.chromosome_order <- function(chromosome) {
  core <- sub("^chr", "", as.character(chromosome), ignore.case = TRUE)
  numeric_core <- suppressWarnings(as.integer(core))
  special <- match(toupper(core), c("X", "Y", "M", "MT"))
  ifelse(!is.na(numeric_core), numeric_core,
         ifelse(!is.na(special), 1000L + special,
                2000L + rank(core, ties.method = "min")))
}

#' Save a locus Z-score plot
#'
#' @param plot A plot returned by [plot_locus_z()].
#' @param path Output path ending in `.pdf`, `.png`, or `.svg`.
#' @param width,height Plot dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return `path`, invisibly.
#' @export
save_locus_z_plot <- function(plot, path, width = 14, height = 4.5, dpi = 300) {
  extension <- tolower(tools::file_ext(path))
  if (!extension %in% c("pdf", "png", "svg")) {
    .acv2g_validation_error("Plot output must end in .pdf, .png, or .svg.")
  }
  device <- if (identical(extension, "svg")) svglite::svglite else extension
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height,
                  units = "in", dpi = dpi, device = device)
  invisible(path)
}
