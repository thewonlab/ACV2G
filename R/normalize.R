.acv2g_validate_group_columns <- function(x, group_cols, object_name = "pairs") {
  if (!is.character(group_cols) || !length(group_cols) || anyNA(group_cols) ||
      any(!nzchar(group_cols)) || anyDuplicated(group_cols)) {
    .acv2g_validation_error(
      "`group_cols` must contain one or more distinct, non-empty column names."
    )
  }
  .require_columns(x, group_cols, object_name)

  for (column in group_cols) {
    value <- x[[column]]
    if (!is.atomic(value) || is.list(value)) {
      .acv2g_validation_error(
        sprintf("Grouping column `%s` must be an atomic vector.", column)
      )
    }
    invalid <- is.na(value)
    if (is.character(value) || is.factor(value)) {
      invalid <- invalid | !nzchar(as.character(value))
    }
    if (is.numeric(value)) {
      invalid <- invalid | !is.finite(value)
    }
    if (any(invalid)) {
      .acv2g_validation_error(
        sprintf(
          "Grouping column `%s` contains %d missing or invalid value(s).",
          column,
          sum(invalid)
        )
      )
    }
  }
  invisible(group_cols)
}

.acv2g_group_key <- function(x, group_cols) {
  if (!nrow(x)) return(character())
  encode <- function(value, column) {
    original_class <- paste(class(value), collapse = ",")
    original_type <- typeof(value)
    if (is.factor(value)) value <- as.character(value)
    text <- switch(
      typeof(value),
      double = sprintf("%a", as.numeric(value)),
      integer = as.character(value),
      logical = ifelse(value, "TRUE", "FALSE"),
      complex = paste0(
        sprintf("%a", Re(value)), ",", sprintf("%a", Im(value))
      ),
      character = enc2utf8(value),
      raw = as.character(as.integer(value)),
      as.character(value)
    )
    tag <- paste(column, original_type, original_class, sep = ":")
    paste0(
      nchar(tag, type = "bytes"), ":", tag,
      nchar(text, type = "bytes"), ":", text
    )
  }
  pieces <- Map(
    function(column) encode(x[[column]], column),
    group_cols
  )
  Reduce(function(left, right) paste0(left, right), pieces)
}

.acv2g_group_label <- function(x, group_cols) {
  if (!nrow(x)) return(character())
  pieces <- lapply(group_cols, function(column) {
    paste0(column, "=", as.character(x[[column]]))
  })
  Reduce(function(left, right) paste0(left, "|", right), pieces)
}

.acv2g_match_groups <- function(x, lookup, group_cols) {
  if (!nrow(x)) return(integer())
  if (!nrow(lookup)) return(rep(NA_integer_, nrow(x)))
  lookup_table <- data.table::as.data.table(lookup[group_cols])
  query_table <- data.table::as.data.table(x[group_cols])
  as.integer(lookup_table[
    query_table,
    on = group_cols,
    which = TRUE
  ])
}

.acv2g_validate_probability <- function(x, arg) {
  .assert_scalar_number(x, arg = arg, lower = 0, upper = 1)
}

.acv2g_empty_normalization_qc <- function(data, group_cols) {
  out <- data[FALSE, group_cols, drop = FALSE]
  out$method <- character()
  out$normalization_group <- character()
  out$normalization_universe_id <- character()
  out$normalization_n_total <- integer()
  out$normalization_n_valid <- integer()
  out$normalization_missing_fraction <- numeric()
  out$normalization_warn_missing_fraction <- numeric()
  out$normalization_mean <- numeric()
  out$normalization_sd <- numeric()
  out$normalization_degenerate <- logical()
  out$normalization_status <- character()
  out$normalization_warnings <- character()
  out$group_id <- character()
  out$status <- character()
  out
}

.acv2g_compute_group_statistics <- function(
    data,
    group_cols,
    proxy_col,
    method,
    warn_missing_fraction) {
  .acv2g_validate_group_columns(data, group_cols, "reference pairs")
  .require_columns(data, proxy_col, "reference pairs")
  .acv2g_validate_numeric_evidence(
    data[[proxy_col]], proxy_col, data, transform = NULL
  )

  if (!nrow(data)) {
    return(.acv2g_empty_normalization_qc(data, group_cols))
  }

  grouped <- data.table::as.data.table(data[c(group_cols, proxy_col)])[
    , {
      values <- .SD[[1L]]
      finite <- is.finite(values)
      n_valid <- sum(finite)
      list(
        normalization_n_total = .N,
        normalization_n_valid = n_valid,
        normalization_mean = if (n_valid) mean(values[finite]) else NA_real_,
        normalization_sd = if (n_valid >= 2L) {
          stats::sd(values[finite])
        } else {
          NA_real_
        }
      )
    },
    by = group_cols,
    .SDcols = proxy_col
  ]
  group_table <- as.data.frame(
    grouped[, group_cols, with = FALSE], stringsAsFactors = FALSE
  )
  n_total <- as.integer(grouped$normalization_n_total)
  n_valid <- as.integer(grouped$normalization_n_valid)
  missing_fraction <- (n_total - n_valid) / n_total
  proxy_mean <- as.numeric(grouped$normalization_mean)
  proxy_sd <- as.numeric(grouped$normalization_sd)
  status <- ifelse(
    n_valid == 0L,
    "all_proxy_missing",
    ifelse(
      n_valid < 2L,
      "insufficient_finite_proxy",
      ifelse(
        !is.finite(proxy_sd),
        "invalid_proxy_sd",
        ifelse(proxy_sd == 0, "constant_proxy", "ok")
      )
    )
  )

  high_missing <- missing_fraction > warn_missing_fraction
  out <- group_table
  out$method <- method
  out$normalization_group <- .acv2g_group_label(group_table, group_cols)
  out$normalization_universe_id <- out$normalization_group
  out$normalization_n_total <- n_total
  out$normalization_n_valid <- n_valid
  out$normalization_missing_fraction <- missing_fraction
  out$normalization_warn_missing_fraction <- warn_missing_fraction
  out$normalization_mean <- proxy_mean
  out$normalization_sd <- proxy_sd
  out$normalization_degenerate <- status != "ok"
  out$normalization_status <- status
  out$normalization_warnings <- ifelse(high_missing, "high_missing_fraction", "")
  # Stable bundle aliases; detailed names remain the authoritative API fields.
  out$group_id <- out$normalization_group
  out$status <- out$normalization_status
  out
}

.acv2g_emit_missing_warning <- function(qc) {
  high <- qc$normalization_warnings == "high_missing_fraction"
  if (!any(high)) {
    return(invisible(qc))
  }
  examples <- paste(utils::head(qc$normalization_group[high], 5L), collapse = ", ")
  .acv2g_warn(
    sprintf(
      paste0(
        "%d normalization group(s) exceed `warn_missing_fraction`; ",
        "scores were retained where normalization was otherwise valid. ",
        "Representative group(s): %s."
      ),
      sum(high),
      examples
    ),
    class = "acv2g_missing_fraction_warning"
  )
  invisible(qc)
}

.acv2g_attach_group_statistics <- function(
    pairs,
    qc,
    group_cols,
    proxy_col,
    z_col,
    method) {
  index <- .acv2g_match_groups(pairs, qc, group_cols)

  metadata <- c(
    "normalization_group",
    "normalization_universe_id",
    "normalization_n_total",
    "normalization_n_valid",
    "normalization_missing_fraction",
    "normalization_warn_missing_fraction",
    "normalization_mean",
    "normalization_sd",
    "normalization_degenerate",
    "normalization_status",
    "normalization_warnings"
  )
  for (column in metadata) {
    pairs[[column]] <- qc[[column]][index]
  }
  pairs$method <- rep(method, nrow(pairs))

  z <- rep(NA_real_, nrow(pairs))
  normalizable <- !is.na(index) & qc$normalization_status[index] == "ok"
  finite_proxy <- is.finite(pairs[[proxy_col]])
  calculate <- normalizable & finite_proxy
  z[calculate] <- (
    pairs[[proxy_col]][calculate] - qc$normalization_mean[index[calculate]]
  ) / qc$normalization_sd[index[calculate]]
  pairs[[z_col]] <- z
  pairs
}

.acv2g_normalization_result <- function(pairs, normalization_qc) {
  structure(
    list(pairs = pairs, normalization_qc = normalization_qc),
    class = c("acv2g_normalization_result", "list")
  )
}

#' Standardize pair proxies within fixed loci
#'
#' Standardizes finite proxy values using the sample standard deviation within
#' each fixed grouping combination. By default, the groups are profile by
#' user-supplied locus. Every input row is retained, and missing proxy values
#' receive a missing Z-score. A high missing fraction produces a warning only;
#' it never suppresses an otherwise calculable score.
#'
#' @param pairs A data frame containing the complete variant-gene pair universe.
#' @param proxy_col Numeric proxy column.
#' @param group_cols Columns defining the fixed normalization group.
#' @param warn_missing_fraction Warn when a group's missing fraction is strictly
#'   greater than this value.
#' @param z_col Output column for locus-pair Z-scores.
#'
#' @return An `acv2g_normalization_result` containing `pairs` and one-row-per-
#'   group `normalization_qc` tables.
#' @export
compute_locus_pair_z <- function(
    pairs,
    proxy_col = "proxy_score",
    group_cols = c("profile_id", "locus_id"),
    warn_missing_fraction = 0.20,
    z_col = "locus_pair_z") {
  pairs <- .as_plain_data_frame(pairs, "pairs")
  .acv2g_validate_probability(warn_missing_fraction, "warn_missing_fraction")
  .acv2g_validate_group_columns(pairs, group_cols, "pairs")
  .require_columns(pairs, proxy_col, "pairs")
  .acv2g_validate_numeric_evidence(
    pairs[[proxy_col]], proxy_col, pairs, transform = NULL
  )
  if (!is.character(z_col) || length(z_col) != 1L || is.na(z_col) || !nzchar(z_col)) {
    .acv2g_validation_error("`z_col` must be one non-empty column name.")
  }

  qc <- .acv2g_compute_group_statistics(
    pairs,
    group_cols,
    proxy_col,
    method = "locus_pair_z",
    warn_missing_fraction = warn_missing_fraction
  )
  .acv2g_emit_missing_warning(qc)
  scored <- .acv2g_attach_group_statistics(
    pairs, qc, group_cols, proxy_col, z_col, method = "locus_pair_z"
  )
  .acv2g_normalization_result(scored, qc)
}

.acv2g_validate_id_columns <- function(x, columns, object_name) {
  .require_columns(x, columns, object_name)
  for (column in columns) {
    value <- as.character(x[[column]])
    invalid <- is.na(value) | !nzchar(value)
    if (any(invalid)) {
      .acv2g_validation_error(
        sprintf(
          "Column `%s` in %s contains %d missing or empty identifier(s).",
          column,
          object_name,
          sum(invalid)
        )
      )
    }
  }
  invisible(x)
}

.acv2g_gene_variant_key <- function(x, gene_col, variant_col) {
  gene <- as.character(x[[gene_col]])
  variant <- as.character(x[[variant_col]])
  paste0(nchar(gene), ":", gene, "|", nchar(variant), ":", variant)
}

.acv2g_validate_genomic_coordinates <- function(
    x,
    chrom_col,
    pos_col,
    object_name) {
  .require_columns(x, c(chrom_col, pos_col), object_name)
  chrom <- as.character(x[[chrom_col]])
  pos <- x[[pos_col]]
  if (!is.numeric(pos)) {
    .acv2g_validation_error(
      sprintf("Coordinate column `%s` in %s must be numeric.", pos_col, object_name)
    )
  }
  invalid <- is.na(chrom) | !nzchar(chrom) | is.na(pos) | !is.finite(pos) |
    pos < 1 | pos != floor(pos)
  if (any(invalid)) {
    .acv2g_validation_error(
      sprintf(
        "%s contains %d invalid chromosome/1-based-position row(s).",
        object_name,
        sum(invalid)
      )
    )
  }
  invisible(x)
}

.acv2g_distance_reference_rows <- function(
    source,
    gene_tss,
    radius,
    gene_col,
    variant_chrom_col,
    variant_pos_col,
    gene_chrom_col,
    gene_tss_col) {
  .acv2g_validate_genomic_coordinates(
    source, variant_chrom_col, variant_pos_col, "reference pairs"
  )

  if (is.null(gene_tss)) {
    .require_columns(source, c(gene_col, gene_chrom_col, gene_tss_col), "reference pairs")
    .acv2g_validate_genomic_coordinates(
      source, gene_chrom_col, gene_tss_col, "reference pairs"
    )
    return(
      as.character(source[[variant_chrom_col]]) == as.character(source[[gene_chrom_col]]) &
        abs(source[[variant_pos_col]] - source[[gene_tss_col]]) <= radius
    )
  }

  gene_tss <- .as_plain_data_frame(gene_tss, "gene_tss")
  .acv2g_validate_id_columns(gene_tss, gene_col, "gene_tss")
  .acv2g_validate_genomic_coordinates(
    gene_tss, gene_chrom_col, gene_tss_col, "gene_tss"
  )
  source_gene <- as.character(source[[gene_col]])
  tss_gene <- as.character(gene_tss[[gene_col]])
  missing_gene <- setdiff(unique(source_gene), unique(tss_gene))
  if (length(missing_gene)) {
    .acv2g_validation_error(
      sprintf(
        "No TSS was supplied for %d reference gene(s); representative ID(s): %s.",
        length(missing_gene),
        paste(utils::head(missing_gene, 5L), collapse = ", ")
      )
    )
  }

  source_ranges <- data.table::data.table(
    .gene = source_gene,
    .chrom = as.character(source[[variant_chrom_col]]),
    .start = as.numeric(source[[variant_pos_col]]),
    .end = as.numeric(source[[variant_pos_col]])
  )
  tss_pos <- as.numeric(gene_tss[[gene_tss_col]])
  tss_ranges <- data.table::data.table(
    .gene = tss_gene,
    .chrom = as.character(gene_tss[[gene_chrom_col]]),
    .start = pmax(1, tss_pos - radius),
    .end = tss_pos + radius
  )
  data.table::setkeyv(
    tss_ranges, c(".gene", ".chrom", ".start", ".end")
  )
  hits <- data.table::foverlaps(
    source_ranges,
    tss_ranges,
    by.x = c(".gene", ".chrom", ".start", ".end"),
    by.y = c(".gene", ".chrom", ".start", ".end"),
    type = "within",
    nomatch = 0L,
    which = TRUE
  )
  eligible <- logical(nrow(source))
  if (nrow(hits)) eligible[unique(hits$xid)] <- TRUE
  eligible
}

.acv2g_explicit_reference_rows <- function(
    source,
    mapping,
    gene_col,
    variant_col) {
  mapping <- .as_plain_data_frame(mapping, "reference_gene_variant_map")
  .acv2g_validate_id_columns(
    mapping, c(gene_col, variant_col), "reference_gene_variant_map"
  )
  mapping <- unique(mapping[c(gene_col, variant_col)])
  source_key <- .acv2g_gene_variant_key(source, gene_col, variant_col)
  mapping_key <- .acv2g_gene_variant_key(mapping, gene_col, variant_col)

  absent <- setdiff(mapping_key, unique(source_key))
  if (length(absent)) {
    absent_rows <- match(absent, mapping_key)
    examples <- paste0(
      as.character(mapping[[gene_col]][absent_rows]),
      "/",
      as.character(mapping[[variant_col]][absent_rows])
    )
    .acv2g_validation_error(
      sprintf(
        paste0(
          "The explicit reference mapping contains %d gene-variant pair(s) ",
          "absent from the selected reference pair table; representative ",
          "pair(s): %s."
        ),
        length(absent),
        paste(utils::head(examples, 5L), collapse = ", ")
      )
    )
  }
  list(eligible = source_key %in% mapping_key, mapping = mapping)
}

.acv2g_deduplicate_reference_pairs <- function(
    source,
    group_cols,
    variant_col,
    proxy_col) {
  if (!nrow(source)) {
    return(source)
  }
  key_cols <- unique(c(group_cols, variant_col))
  source_table <- data.table::as.data.table(source)
  conflict <- source_table[
    , {
      values <- get(proxy_col)
      value_labels <- ifelse(
        is.na(values), "<missing>", format(values, digits = 17L)
      )
      list(.acv2g_conflicting_proxy = length(unique(value_labels)) > 1L)
    },
    by = key_cols
  ]
  if (any(conflict$.acv2g_conflicting_proxy)) {
    .acv2g_validation_error(
      paste0(
        "A reference gene-variant anchor has conflicting proxy values; ",
        "reference anchors must be unique or exact duplicates."
      )
    )
  }
  keep <- !duplicated(source_table, by = key_cols)
  source[keep, , drop = FALSE]
}

.acv2g_add_absent_target_groups <- function(
    qc,
    pairs,
    group_cols,
    method,
    warn_missing_fraction) {
  target_groups <- unique(pairs[group_cols])
  absent <- is.na(.acv2g_match_groups(target_groups, qc, group_cols))
  if (!any(absent)) {
    return(qc)
  }

  extra <- target_groups[absent, , drop = FALSE]
  n <- nrow(extra)
  extra$method <- method
  extra$normalization_group <- .acv2g_group_label(extra, group_cols)
  extra$normalization_universe_id <- extra$normalization_group
  extra$normalization_n_total <- rep.int(0L, n)
  extra$normalization_n_valid <- rep.int(0L, n)
  extra$normalization_missing_fraction <- rep(NA_real_, n)
  extra$normalization_warn_missing_fraction <- rep(warn_missing_fraction, n)
  extra$normalization_mean <- rep(NA_real_, n)
  extra$normalization_sd <- rep(NA_real_, n)
  extra$normalization_degenerate <- rep(TRUE, n)
  extra$normalization_status <- rep("no_reference_variants", n)
  extra$normalization_warnings <- rep("", n)
  extra$group_id <- extra$normalization_group
  extra$status <- extra$normalization_status
  for (column in setdiff(names(qc), names(extra))) {
    extra[[column]] <- qc[[column]][rep(NA_integer_, n)]
  }
  for (column in setdiff(names(extra), names(qc))) {
    qc[[column]] <- extra[[column]][rep(NA_integer_, nrow(qc))]
  }
  extra <- extra[names(qc)]
  rbind(qc, extra)
}

#' Standardize pair proxies using an explicit per-gene reference universe
#'
#' `compute_gene_wise_z()` is independent of locus-pair normalization. The
#' caller must explicitly choose whether the reference rows come from the input
#' pairs or a separate background pair table. Within that source, reference
#' variants are selected by an inclusive TSS radius, or by a supplied mapping
#' that completely overrides distance selection. Target rows are never dropped.
#'
#' @param pairs Target variant-gene pairs to score.
#' @param reference_mode Explicitly either `"input_variants"` or
#'   `"background_variants"`; there is no default or fallback.
#' @param background_pairs Pair-level evidence for a separate background. It is
#'   required only for `reference_mode = "background_variants"`.
#' @param gene_tss Optional data frame of one or more TSSs per gene. When
#'   omitted, TSS coordinates are read from each reference-pair row.
#' @param reference_gene_radius_bp Inclusive distance around any declared TSS.
#' @param reference_gene_variant_map Optional gene-variant mapping that fully
#'   replaces distance selection.
#' @param reference_set_id Optional stable identifier for the selected source.
#' @param reference_gene_variant_map_id Optional mapping provenance identifier.
#' @param proxy_col Numeric proxy column present in target and reference pairs.
#' @param group_cols Per-gene normalization columns; this must include
#'   `gene_col` and normally also includes `profile_id`.
#' @param gene_col Gene identifier column.
#' @param variant_col Coordinate-level variant-anchor identifier column.
#' @param variant_chrom_col Reference-variant chromosome column.
#' @param variant_pos_col Reference-variant 1-based position column.
#' @param gene_chrom_col TSS chromosome column.
#' @param gene_tss_col TSS 1-based position column.
#' @param warn_missing_fraction Warn when a reference group's missing fraction
#'   is strictly greater than this value.
#' @param z_col Output column for gene-wise Z-scores.
#'
#' @return An `acv2g_normalization_result` containing full target `pairs` and
#'   detailed per-profile, per-gene `normalization_qc`.
#' @export
compute_gene_wise_z <- function(
    pairs,
    reference_mode,
    background_pairs = NULL,
    gene_tss = NULL,
    reference_gene_radius_bp = 1000000L,
    reference_gene_variant_map = NULL,
    reference_set_id = NULL,
    reference_gene_variant_map_id = NULL,
    proxy_col = "proxy_score",
    group_cols = c("profile_id", "gene_id"),
    gene_col = "gene_id",
    variant_col = "variant_key",
    variant_chrom_col = "chrom",
    variant_pos_col = "pos1",
    gene_chrom_col = "tss_chrom",
    gene_tss_col = "tss_pos1",
    warn_missing_fraction = 0.20,
    z_col = "gene_wise_z") {
  pairs <- .as_plain_data_frame(pairs, "pairs")
  if (missing(reference_mode)) {
    .acv2g_validation_error(
      "`reference_mode` must be explicitly set to `input_variants` or `background_variants`."
    )
  }
  reference_mode <- .assert_scalar_choice(
    reference_mode,
    c("input_variants", "background_variants"),
    "reference_mode"
  )
  .acv2g_validate_probability(warn_missing_fraction, "warn_missing_fraction")
  .assert_scalar_number(
    reference_gene_radius_bp,
    arg = "reference_gene_radius_bp",
    lower = 0,
    integer = TRUE
  )
  .acv2g_validate_group_columns(pairs, group_cols, "pairs")
  if (!gene_col %in% group_cols) {
    .acv2g_validation_error("`group_cols` must include `gene_col` for gene-wise normalization.")
  }
  .acv2g_validate_id_columns(pairs, c(gene_col, variant_col), "pairs")
  .require_columns(pairs, proxy_col, "pairs")
  .acv2g_validate_numeric_evidence(
    pairs[[proxy_col]], proxy_col, pairs, transform = NULL
  )
  if (!is.character(z_col) || length(z_col) != 1L || is.na(z_col) || !nzchar(z_col)) {
    .acv2g_validation_error("`z_col` must be one non-empty column name.")
  }

  if (identical(reference_mode, "input_variants")) {
    if (!is.null(background_pairs)) {
      .acv2g_validation_error(
        "`background_pairs` must be `NULL` when `reference_mode = \"input_variants\"`."
      )
    }
    source <- pairs
  } else {
    if (is.null(background_pairs)) {
      .acv2g_validation_error(
        "`background_pairs` is required when `reference_mode = \"background_variants\"`."
      )
    }
    source <- .as_plain_data_frame(background_pairs, "background_pairs")
  }

  .acv2g_validate_group_columns(source, group_cols, "reference pairs")
  .acv2g_validate_id_columns(source, c(gene_col, variant_col), "reference pairs")
  .require_columns(source, proxy_col, "reference pairs")
  .acv2g_validate_numeric_evidence(
    source[[proxy_col]], proxy_col, source, transform = NULL
  )

  if (is.null(reference_gene_variant_map)) {
    eligible <- .acv2g_distance_reference_rows(
      source,
      gene_tss,
      reference_gene_radius_bp,
      gene_col,
      variant_chrom_col,
      variant_pos_col,
      gene_chrom_col,
      gene_tss_col
    )
    universe_rule <- "distance"
    mapping_id <- NA_character_
    radius_record <- as.numeric(reference_gene_radius_bp)
  } else {
    mapped <- .acv2g_explicit_reference_rows(
      source,
      reference_gene_variant_map,
      gene_col,
      variant_col
    )
    eligible <- mapped$eligible
    universe_rule <- "explicit_mapping"
    mapping_id <- reference_gene_variant_map_id
    if (is.null(mapping_id)) {
      mapping_key <- sort(.acv2g_gene_variant_key(mapped$mapping, gene_col, variant_col))
      mapping_id <- paste0(
        "mapping_",
        digest::digest(mapping_key, algo = "xxhash64", serialize = TRUE)
      )
    }
    if (!is.character(mapping_id) || length(mapping_id) != 1L ||
        is.na(mapping_id) || !nzchar(mapping_id)) {
      .acv2g_validation_error(
        "`reference_gene_variant_map_id` must be one non-empty identifier."
      )
    }
    radius_record <- NA_real_
  }

  reference <- source[eligible, , drop = FALSE]
  reference <- .acv2g_deduplicate_reference_pairs(
    reference, group_cols, variant_col, proxy_col
  )
  qc <- .acv2g_compute_group_statistics(
    reference,
    group_cols,
    proxy_col,
    method = "gene_wise_z",
    warn_missing_fraction = warn_missing_fraction
  )
  qc <- .acv2g_add_absent_target_groups(
    qc,
    pairs,
    group_cols,
    method = "gene_wise_z",
    warn_missing_fraction = warn_missing_fraction
  )

  if (is.null(reference_set_id)) {
    reference_set_id <- reference_mode
  }
  if (!is.character(reference_set_id) || length(reference_set_id) != 1L ||
      is.na(reference_set_id) || !nzchar(reference_set_id)) {
    .acv2g_validation_error("`reference_set_id` must be one non-empty identifier.")
  }
  qc$normalization_universe_id <- paste0(
    reference_set_id, "::", qc$normalization_group
  )
  qc$reference_mode <- rep(reference_mode, nrow(qc))
  qc$reference_universe_rule <- rep(universe_rule, nrow(qc))
  qc$reference_set_id <- rep(reference_set_id, nrow(qc))
  qc$reference_gene_radius_bp <- rep(radius_record, nrow(qc))
  qc$reference_gene_variant_map_id <- rep(mapping_id, nrow(qc))
  qc$reference_variant_count <- qc$normalization_n_total
  qc$reference_pair_count <- qc$normalization_n_total

  .acv2g_emit_missing_warning(qc)
  scored <- .acv2g_attach_group_statistics(
    pairs, qc, group_cols, proxy_col, z_col, method = "gene_wise_z"
  )
  scored$reference_mode <- rep(reference_mode, nrow(scored))
  scored$reference_universe_rule <- rep(universe_rule, nrow(scored))
  scored$reference_set_id <- rep(reference_set_id, nrow(scored))
  scored$reference_gene_radius_bp <- rep(radius_record, nrow(scored))
  scored$reference_gene_variant_map_id <- rep(mapping_id, nrow(scored))
  .acv2g_normalization_result(scored, qc)
}
