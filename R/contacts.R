#' Extract and aggregate variant-to-TSS contact evidence
#'
#' `extract_v2g_contacts()` maps variant and TSS anchors to exact fixed-width,
#' 0-based half-open contact bins, resolves one source-provided contact layer,
#' applies the requested diagonal policy at each unique TSS bin, and then takes
#' the finite maximum across a gene's eligible TSS bins.  Explicit numeric zero
#' is retained.  Unavailable table pixels and unavailable balancing weights are
#' retained as `NA` with machine-readable status values. A failed native tile
#' query warns, leaves only its affected keys as `missing_contact_query`, and
#' does not discard successfully queried tiles.
#'
#' Quantitative tables must be declared complete (not called-loop/sparse-list
#' inputs) and must explicitly map score columns to `raw`, `balanced`, or
#' `named_score` semantics.  Native `.cool`, `.mcool`, and `.hic` files are
#' opened with [HiCExperiment::CoolFile()] or [HiCExperiment::HicFile()].
#' Cooler tiles are imported with [HiCExperiment::import()], while bounded
#' `.hic` tiles are queried with [strawr::straw()] so raw fallback does not
#' depend on a KR layer. ACV2G does not hand-parse these formats.
#'
#' @param pairs A data frame containing `variant_id`, `gene_id`, `chrom`,
#'   1-based `pos1`, `tss_chrom`, and 1-based `tss_pos1`.  Supply one row per
#'   eligible TSS.  `tss_id` and `locus_id` are retained when present.
#' @param contacts A quantitative contact table or one `.cool`, `.mcool`, or
#'   `.hic` path.
#' @param source_type `"auto"`, `"table"`, `"cool"`, `"mcool"`, or
#'   `"hic"`.
#' @param contact_resolution_bp Positive fixed bin width.  A single-resolution
#'   native/table source supplies its intrinsic width.  A multi-resolution
#'   native source uses exactly 10 kb when available unless this is supplied.
#' @param contact_value Requested score semantics: `"balanced"` (default),
#'   `"raw"`, or `"named_score"`.
#' @param contact_layer Optional exact layer name.  A missing explicitly named
#'   layer is always an error and never triggers raw fallback.
#' @param contact_score_col,contact_score_type For a one-score table, the
#'   numeric score column and its explicit semantics (`"balanced"`, `"raw"`,
#'   or `"named_score"`).
#' @param contact_layers Optional named character vector or named list mapping
#'   score semantics to table columns, for example
#'   `list(raw = "count", balanced = "KR")`.  Multiple balanced columns are
#'   allowed only when `contact_layer` explicitly selects one.
#' @param table_representation Must be `"complete"`; sparse called-loop lists
#'   are not a valid sole source in this release.
#' @param diagonal_policy `"neighbor_max"` (default), `"observed"`, or
#'   `"keep_na"`.  The default replaces a diagonal value with the finite
#'   maximum of the two immediate off-diagonal pixels.
#' @param distance_policy Must be `"none"` in this release.
#' @param contact_aggregate Must be `"max"` in this release.
#' @param contact_transform `"log1p"` (default) or `"identity"`.  `log1p`
#'   rejects finite negative resolved contact values.
#' @param contact_assay_id,contact_assembly,contact_units,contact_normalization
#'   Provenance labels.
#' @param contact_weights_source_provided Whether a tabular balanced score used
#'   source-provided weights.  Native balanced layers record `TRUE`.
#'
#' @return A classed list with `pairs` (one row per variant-gene pair),
#'   `tss_bin_contacts` (one row per unique same-gene TSS bin), and `qc`
#'   (resolved source/layer/normalization/resolution and query-failure
#'   provenance).
#' @export
extract_v2g_contacts <- function(
    pairs,
    contacts,
    source_type = c("auto", "table", "cool", "mcool", "hic"),
    contact_resolution_bp = NULL,
    contact_value = c("balanced", "raw", "named_score"),
    contact_layer = NULL,
    contact_score_col = NULL,
    contact_score_type = NULL,
    contact_layers = NULL,
    table_representation = "complete",
    diagonal_policy = c("neighbor_max", "observed", "keep_na"),
    distance_policy = "none",
    contact_aggregate = "max",
    contact_transform = c("log1p", "identity"),
    contact_assay_id = "contact",
    contact_assembly = NA_character_,
    contact_units = NA_character_,
    contact_normalization = NA_character_,
    contact_weights_source_provided = NA) {
  pairs <- .as_plain_data_frame(pairs, "pairs")
  .require_columns(
    pairs,
    c("chrom", "pos1", "tss_chrom", "tss_pos1"),
    "pairs"
  )
  if (!any(c("variant_anchor_id", "variant_id") %in% names(pairs))) {
    .acv2g_validation_error(
      "{.arg pairs} must contain {.field variant_anchor_id} or {.field variant_id}."
    )
  }
  if (!any(c("gene_id_normalized", "gene_id") %in% names(pairs))) {
    .acv2g_validation_error(
      "{.arg pairs} must contain {.field gene_id_normalized} or {.field gene_id}."
    )
  }
  if (!nrow(pairs)) {
    .acv2g_validation_error("{.arg pairs} must contain at least one eligible pair/TSS row.")
  }
  source_type <- match.arg(source_type)
  contact_value <- match.arg(contact_value)
  diagonal_policy <- match.arg(diagonal_policy)
  contact_transform <- match.arg(contact_transform)
  if (!identical(distance_policy, "none")) {
    .acv2g_validation_error(
      "Only {.arg distance_policy = 'none'} is implemented in this release."
    )
  }
  if (!identical(contact_aggregate, "max")) {
    .acv2g_validation_error(
      "Only {.arg contact_aggregate = 'max'} is implemented in this release."
    )
  }
  if (!identical(table_representation, "complete")) {
    .acv2g_validation_error(
      "Sparse called-loop lists cannot be the sole contact source; provide a complete quantitative contact table or native matrix."
    )
  }
  if (!is.null(contact_layer) &&
      (!is.character(contact_layer) || length(contact_layer) != 1L ||
       is.na(contact_layer) || !nzchar(contact_layer))) {
    .acv2g_validation_error("{.arg contact_layer} must be one non-empty layer name.")
  }

  pairs <- .contact_validate_pairs(pairs)
  source_type <- .contact_source_type(contacts, source_type)
  prepared <- if (source_type == "table") {
    .contact_prepare_table(
      contacts = contacts,
      requested_resolution = contact_resolution_bp,
      contact_score_col = contact_score_col,
      contact_score_type = contact_score_type,
      contact_layers = contact_layers
    )
  } else {
    .contact_prepare_native(
      pairs = pairs,
      path = contacts,
      source_type = source_type,
      requested_resolution = contact_resolution_bp,
      diagonal_policy = diagonal_policy,
      contact_value = contact_value,
      contact_layer = contact_layer
    )
  }
  resolved_resolution <- prepared$resolution
  selected <- prepared$selected %||% .contact_resolve_layer(
    layer_roles = prepared$layer_roles,
    contact_value = contact_value,
    contact_layer = contact_layer
  )

  selected_values <- suppressWarnings(as.numeric(
    prepared$pixels[[selected$layer]]
  ))
  if (any(is.infinite(selected_values))) {
    .acv2g_validation_error(
      "Resolved contact layer {.field {selected$layer}} contains infinite values."
    )
  }
  if (contact_transform == "log1p" &&
      any(is.finite(selected_values) & selected_values < 0)) {
    .acv2g_validation_error(
      "{.arg contact_transform = 'log1p'} does not allow negative values in resolved contact layer {.field {selected$layer}}."
    )
  }
  prepared$pixels$.selected_contact_value <- selected_values
  prepared$selected_type <- selected$type

  detail_seed <- .contact_make_tss_bin_seed(
    pairs, resolved_resolution, prepared$chrom_lengths
  )
  direct <- .contact_lookup_pixels(
    detail_seed$chrom,
    detail_seed$contact_variant_bin_start0,
    detail_seed$tss_chrom,
    detail_seed$contact_promoter_bin_start0,
    prepared
  )
  detail <- .contact_apply_diagonal(
    detail_seed, direct, prepared, diagonal_policy, resolved_resolution
  )

  if (contact_transform == "log1p" &&
      any(is.finite(detail$contact_raw) & detail$contact_raw < 0)) {
    .acv2g_validation_error(
      "{.arg contact_transform = 'log1p'} does not allow negative contact evidence after diagonal resolution."
    )
  }
  detail$contact_transform <- contact_transform
  detail$contact_transformed <- if (contact_transform == "log1p") {
    log1p(detail$contact_raw)
  } else {
    detail$contact_raw
  }

  requested_resolution <- if (is.null(contact_resolution_bp)) {
    NA_real_
  } else {
    as.numeric(contact_resolution_bp)
  }
  weights_provided <- if (source_type == "table") {
    if (selected$type == "balanced") {
      declared <- as.logical(contact_weights_source_provided)[1L]
      if (is.na(declared)) TRUE else declared
    } else {
      FALSE
    }
  } else {
    identical(selected$type, "balanced")
  }
  fallback_reason <- if (selected$fallback_used) {
    "balanced_layer_absent"
  } else {
    NA_character_
  }
  resolved_normalization <- if (!is.na(contact_normalization)[1L]) {
    as.character(contact_normalization)
  } else if (identical(source_type, "hic") &&
             identical(selected$layer, "count")) {
    "NONE"
  } else {
    as.character(selected$layer)
  }
  detail$contact_assay_id <- as.character(contact_assay_id)
  detail$contact_assembly <- as.character(contact_assembly)
  detail$contact_value_requested <- contact_value
  detail$contact_value_resolved <- selected$type
  detail$contact_layer_requested <- if (is.null(contact_layer)) NA_character_ else contact_layer
  detail$contact_layer_resolved <- selected$layer
  detail$contact_fallback_used <- selected$fallback_used
  detail$contact_fallback_reason <- fallback_reason
  detail$contact_weights_source_provided <- weights_provided
  detail$contact_units <- as.character(contact_units)
  detail$contact_resolution_requested <- requested_resolution
  detail$contact_resolution_resolved <- resolved_resolution
  detail$contact_score_type <- selected$type
  detail$contact_normalization <- resolved_normalization
  detail$contact_aggregate <- contact_aggregate
  detail$distance_policy <- distance_policy
  detail$diagonal_policy <- diagonal_policy

  aggregated <- .contact_aggregate_tss_bins(
    detail, pairs, contact_transform
  )
  detail$contact_is_max <- attr(aggregated, "detail_contact_is_max")
  attr(aggregated, "detail_contact_is_max") <- NULL
  detail$.pair_key <- NULL
  provenance_columns <- c(
    "contact_assay_id", "contact_assembly", "contact_value_requested",
    "contact_value_resolved", "contact_layer_requested",
    "contact_layer_resolved", "contact_fallback_used",
    "contact_fallback_reason", "contact_weights_source_provided",
    "contact_units", "contact_resolution_requested",
    "contact_resolution_resolved", "contact_score_type",
    "contact_normalization", "contact_aggregate", "distance_policy",
    "diagonal_policy", "contact_transform"
  )
  for (column in provenance_columns) {
    aggregated[[column]] <- detail[[column]][1L]
  }

  qc <- data.frame(
    contact_assay_id = as.character(contact_assay_id),
    contact_source_type = source_type,
    contact_source = if (source_type == "table") {
      "<in-memory-table>"
    } else {
      prepared$source_path
    },
    contact_source_input = if (source_type == "table") {
      "<in-memory-table>"
    } else {
      prepared$source_path_input
    },
    contact_resolution_requested = requested_resolution,
    contact_resolution_resolved = resolved_resolution,
    contact_value_requested = contact_value,
    contact_value_resolved = selected$type,
    contact_layer_requested = if (is.null(contact_layer)) NA_character_ else contact_layer,
    contact_layer_resolved = selected$layer,
    contact_fallback_used = selected$fallback_used,
    contact_fallback_reason = fallback_reason,
    contact_available_layers = paste(prepared$layer_roles$layer, collapse = ";"),
    contact_weights_source_provided = weights_provided,
    contact_score_type = selected$type,
    contact_normalization = resolved_normalization,
    contact_query_failure_count = nrow(
      prepared$query_failures %||% .contact_empty_query_failures()
    ),
    contact_query_failed_foci = paste(
      unique((prepared$query_failures %||%
        .contact_empty_query_failures())$focus),
      collapse = ";"
    ),
    contact_query_failure_messages = paste(
      unique(utils::head(
        (prepared$query_failures %||%
          .contact_empty_query_failures())$message,
        5L
      )),
      collapse = ";"
    ),
    distance_policy = distance_policy,
    diagonal_policy = diagonal_policy,
    contact_aggregate = contact_aggregate,
    stringsAsFactors = FALSE
  )
  result <- list(
    pairs = aggregated,
    tss_bin_contacts = detail,
    qc = qc
  )
  class(result) <- c("acv2g_contact_result", "list")
  result
}

.contact_source_type <- function(contacts, source_type) {
  if (source_type != "auto") {
    if (source_type == "table" && !is.data.frame(contacts)) {
      .acv2g_validation_error(
        "{.arg contacts} must be a data frame when {.arg source_type = 'table'}."
      )
    }
    if (source_type != "table" &&
        (!is.character(contacts) || length(contacts) != 1L || is.na(contacts))) {
      .acv2g_validation_error(
        "{.arg contacts} must be one native contact file path."
      )
    }
    return(source_type)
  }
  if (is.data.frame(contacts)) return("table")
  if (!is.character(contacts) || length(contacts) != 1L || is.na(contacts)) {
    .acv2g_validation_error(
      "{.arg contacts} must be a quantitative table or one native contact path."
    )
  }
  lower <- tolower(contacts)
  if (grepl("\\.mcool$", lower)) return("mcool")
  if (grepl("\\.cool$", lower)) return("cool")
  if (grepl("\\.hic$", lower)) return("hic")
  .acv2g_validation_error(
    "Cannot infer contact format from {.file {contacts}}; set {.arg source_type}."
  )
}

.contact_validate_pairs <- function(pairs) {
  if (!"variant_anchor_id" %in% names(pairs)) {
    pairs$variant_anchor_id <- as.character(pairs$variant_id)
  }
  if (!"variant_id" %in% names(pairs)) {
    pairs$variant_id <- as.character(pairs$variant_anchor_id)
  }
  if (!"gene_id_normalized" %in% names(pairs)) {
    pairs$gene_id_normalized <- as.character(pairs$gene_id)
  }
  if (!"gene_id" %in% names(pairs)) {
    pairs$gene_id <- as.character(pairs$gene_id_normalized)
  }
  for (column in c("variant_id", "gene_id", "chrom", "tss_chrom")) {
    pairs[[column]] <- as.character(pairs[[column]])
    if (anyNA(pairs[[column]]) || any(!nzchar(pairs[[column]]))) {
      .acv2g_validation_error("{.field {column}} cannot be missing or empty.")
    }
  }
  for (column in c("pos1", "tss_pos1")) {
    value <- suppressWarnings(as.numeric(pairs[[column]]))
    if (anyNA(value) || any(!is.finite(value)) || any(value < 1) ||
        any(value != floor(value))) {
      .acv2g_validation_error(
        "{.field {column}} must contain finite, positive, 1-based integer positions."
      )
    }
    pairs[[column]] <- value
  }
  if (!"tss_id" %in% names(pairs)) {
    pairs$tss_id <- paste0(pairs$tss_chrom, ":", pairs$tss_pos1)
  }
  pairs$tss_id <- as.character(pairs$tss_id)
  missing_id <- is.na(pairs$tss_id) | !nzchar(pairs$tss_id)
  pairs$tss_id[missing_id] <- paste0(
    pairs$tss_chrom[missing_id], ":", pairs$tss_pos1[missing_id]
  )
  for (column in c("variant_anchor_id", "gene_id_normalized")) {
    pairs[[column]] <- as.character(pairs[[column]])
    if (anyNA(pairs[[column]]) || any(!nzchar(pairs[[column]]))) {
      .acv2g_validation_error("{.field {column}} cannot be missing or empty.")
    }
  }
  if (!"pair_id" %in% names(pairs)) {
    locus <- if ("locus_id" %in% names(pairs)) pairs$locus_id else ""
    pairs$pair_id <- .stable_row_id(
      locus, pairs$variant_anchor_id, pairs$gene_id_normalized,
      prefix = "pair"
    )
  }
  pairs$pair_id <- as.character(pairs$pair_id)
  if (anyNA(pairs$pair_id) || any(!nzchar(pairs$pair_id))) {
    .acv2g_validation_error("{.field pair_id} cannot be missing or empty.")
  }
  pairs
}

.contact_pair_columns <- function(pairs) {
  candidates <- c(
    "pair_id", "locus_id", "variant_anchor_id", "variant_id",
    "gene_id_normalized", "gene_id"
  )
  intersect(candidates, names(pairs))
}

.contact_pair_key <- function(pairs) {
  if ("pair_id" %in% names(pairs)) return(as.character(pairs$pair_id))
  columns <- .contact_pair_columns(pairs)
  do.call(paste, c(pairs[columns], sep = "\u001f"))
}

.contact_resolution <- function(available, requested, native = TRUE) {
  available <- sort(unique(suppressWarnings(as.numeric(available))))
  available <- available[is.finite(available) & available > 0 &
                           available == floor(available)]
  if (!length(available)) {
    .acv2g_validation_error("Contact source exposes no valid fixed resolution.")
  }
  if (!is.null(requested)) {
    requested <- .assert_scalar_number(
      requested, "contact_resolution_bp", lower = 1, integer = TRUE
    )
  }
  if (length(available) == 1L) {
    if (!is.null(requested) && requested != available) {
      .acv2g_validation_error(
        "Requested contact resolution {requested} conflicts with the source's intrinsic resolution {available}."
      )
    }
    return(available)
  }
  selected <- if (is.null(requested)) 10000 else requested
  if (!selected %in% available) {
    .acv2g_validation_error(
      "Contact resolution {selected} is unavailable. Available resolutions: {paste(available, collapse = ', ')}."
    )
  }
  selected
}

.contact_layer_map <- function(
    available_columns, contact_score_col, contact_score_type, contact_layers) {
  accepted <- c("balanced", "raw", "named_score")
  records <- list()
  if (!is.null(contact_score_col)) {
    if (is.null(contact_score_type) ||
        !is.character(contact_score_type) || length(contact_score_type) != 1L ||
        !contact_score_type %in% accepted) {
      .acv2g_validation_error(
        "{.arg contact_score_type} must explicitly be balanced, raw, or named_score when {.arg contact_score_col} is used."
      )
    }
    records[[length(records) + 1L]] <- data.frame(
      layer = contact_score_col, type = contact_score_type,
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(contact_layers)) {
    if (is.character(contact_layers) && !is.null(names(contact_layers))) {
      contact_layers <- split(unname(contact_layers), names(contact_layers))
    }
    if (!is.list(contact_layers) || is.null(names(contact_layers)) ||
        any(!names(contact_layers) %in% accepted)) {
      .acv2g_validation_error(
        "{.arg contact_layers} must be a named character vector/list keyed by balanced, raw, or named_score."
      )
    }
    for (type in names(contact_layers)) {
      layers <- as.character(contact_layers[[type]])
      if (!length(layers) || anyNA(layers) || any(!nzchar(layers))) {
        .acv2g_validation_error(
          "Every {.arg contact_layers} entry must name at least one non-empty score column."
        )
      }
      records[[length(records) + 1L]] <- data.frame(
        layer = layers, type = rep(type, length(layers)),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(records)) {
    .acv2g_validation_error(
      "A quantitative contact table requires {.arg contact_score_col} plus {.arg contact_score_type}, or an explicit {.arg contact_layers} mapping."
    )
  }
  map <- unique(do.call(rbind, records))
  if (anyDuplicated(map$layer)) {
    duplicated_layers <- unique(map$layer[duplicated(map$layer) |
      duplicated(map$layer, fromLast = TRUE)])
    .acv2g_validation_error(
      "Contact layer{?s} mapped to more than one semantic type: {paste(duplicated_layers, collapse = ', ')}."
    )
  }
  missing_layers <- setdiff(map$layer, available_columns)
  if (length(missing_layers)) {
    .acv2g_validation_error(
      "Mapped contact layer{?s} unavailable: {paste(missing_layers, collapse = ', ')}. Available columns: {paste(available_columns, collapse = ', ')}."
    )
  }
  map
}

.contact_prepare_table <- function(
    contacts, requested_resolution, contact_score_col, contact_score_type,
    contact_layers) {
  contacts <- .as_plain_data_frame(contacts, "contacts")
  .require_columns(
    contacts,
    c("chrom1", "start1_0", "end1_0", "chrom2", "start2_0", "end2_0"),
    "contacts"
  )
  layer_roles <- .contact_layer_map(
    names(contacts), contact_score_col, contact_score_type, contact_layers
  )
  for (column in c("start1_0", "end1_0", "start2_0", "end2_0")) {
    contacts[[column]] <- suppressWarnings(as.numeric(contacts[[column]]))
  }
  bad <- !is.finite(contacts$start1_0) | !is.finite(contacts$end1_0) |
    !is.finite(contacts$start2_0) | !is.finite(contacts$end2_0) |
    contacts$start1_0 < 0 | contacts$start2_0 < 0 |
    contacts$end1_0 <= contacts$start1_0 |
    contacts$end2_0 <= contacts$start2_0 |
    contacts$start1_0 != floor(contacts$start1_0) |
    contacts$end1_0 != floor(contacts$end1_0) |
    contacts$start2_0 != floor(contacts$start2_0) |
    contacts$end2_0 != floor(contacts$end2_0)
  if (any(bad)) {
    .acv2g_validation_error(
      "Contact table coordinates must be valid 0-based half-open intervals."
    )
  }
  widths <- unique(c(
    contacts$end1_0 - contacts$start1_0,
    contacts$end2_0 - contacts$start2_0
  ))
  resolution <- .contact_resolution(widths, requested_resolution, native = FALSE)
  aligned <- contacts$start1_0 %% resolution == 0 &
    contacts$start2_0 %% resolution == 0 &
    contacts$end1_0 - contacts$start1_0 == resolution &
    contacts$end2_0 - contacts$start2_0 == resolution
  if (any(!aligned)) {
    .acv2g_validation_error(
      "Every contact-table interval must align exactly to the resolved fixed-width bins."
    )
  }
  contacts$chrom1 <- as.character(contacts$chrom1)
  contacts$chrom2 <- as.character(contacts$chrom2)
  key <- .contact_pixel_key(
    contacts$chrom1, contacts$start1_0,
    contacts$chrom2, contacts$start2_0
  )
  if (anyDuplicated(key)) {
    .acv2g_validation_error(
      "Contact table contains duplicate pixels after symmetric orientation normalization."
    )
  }
  pixels <- data.frame(
    .pixel_key = key,
    chrom1 = contacts$chrom1,
    start1_0 = contacts$start1_0,
    chrom2 = contacts$chrom2,
    start2_0 = contacts$start2_0,
    stringsAsFactors = FALSE
  )
  for (layer in layer_roles$layer) {
    values <- suppressWarnings(as.numeric(contacts[[layer]]))
    if (any(is.na(values) & !is.na(contacts[[layer]]))) {
      .acv2g_validation_error(
        "Contact layer {.field {layer}} must be numeric."
      )
    }
    pixels[[layer]] <- values
  }
  list(
    pixels = pixels,
    layer_roles = layer_roles,
    resolution = resolution,
    implicit_zero = FALSE,
    available_chromosomes = unique(c(contacts$chrom1, contacts$chrom2)),
    chrom_lengths = NULL,
    bin_weights = NULL,
    native_type = NULL
  )
}

.contact_quiet_native_call <- function(expr) {
  value <- NULL
  invisible(utils::capture.output(
    value <- suppressMessages(expr)
  ))
  value
}

.contact_empty_query_failures <- function() {
  data.frame(
    backend = character(),
    focus = character(),
    layer = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

.contact_query_failure_row <- function(backend, focus, layer, condition) {
  data.frame(
    backend = as.character(backend),
    focus = as.character(focus),
    layer = if (is.null(layer)) NA_character_ else as.character(layer),
    message = conditionMessage(condition),
    stringsAsFactors = FALSE
  )
}

.contact_native_file <- function(path, source_type) {
  if (!file.exists(path)) {
    .acv2g_validation_error("Contact file does not exist: {.file {path}}.")
  }
  tryCatch(
    .contact_quiet_native_call(
      if (source_type %in% c("cool", "mcool")) {
        HiCExperiment::CoolFile(path)
      } else {
        HiCExperiment::HicFile(path)
      }
    ),
    error = function(e) {
      .acv2g_validation_error(
        "Could not open native contact file {.file {path}} with HiCExperiment: {conditionMessage(e)}"
      )
    }
  )
}

.contact_native_query_rows <- function(
    pairs, resolution, diagonal_policy, chrom_lengths = NULL) {
  seed <- .contact_make_tss_bin_seed(pairs, resolution, chrom_lengths)
  query <- data.frame(
    chrom1 = seed$chrom,
    start1_0 = seed$contact_variant_bin_start0,
    chrom2 = seed$tss_chrom,
    start2_0 = seed$contact_promoter_bin_start0,
    stringsAsFactors = FALSE
  )
  diagonal <- query$chrom1 == query$chrom2 &
    query$start1_0 == query$start2_0
  if (identical(diagonal_policy, "neighbor_max") && any(diagonal)) {
    left <- query[diagonal, , drop = FALSE]
    left$start2_0 <- left$start2_0 - resolution
    right <- query[diagonal, , drop = FALSE]
    right$start2_0 <- right$start2_0 + resolution
    query <- rbind(query, left, right)
  }
  query$.pixel_key <- .contact_pixel_key(
    query$chrom1, query$start1_0, query$chrom2, query$start2_0
  )
  query[!duplicated(query$.pixel_key), , drop = FALSE]
}

.contact_native_query_plan <- function(
    query,
    resolution,
    available_chromosomes,
    chrom_lengths = NULL,
    tile_n_bins = 256L) {
  tile_n_bins <- .assert_scalar_number(
    tile_n_bins, "tile_n_bins", lower = 1, integer = TRUE
  )
  tile_span <- as.numeric(resolution) * tile_n_bins
  chromosome_ok <- if (length(available_chromosomes)) {
    query$chrom1 %in% available_chromosomes &
      query$chrom2 %in% available_chromosomes
  } else {
    rep(TRUE, nrow(query))
  }
  coordinate_ok <- query$start1_0 >= 0 & query$start2_0 >= 0 &
    query$start1_0 %% resolution == 0 & query$start2_0 %% resolution == 0
  if (!is.null(chrom_lengths)) {
    length1 <- unname(chrom_lengths[query$chrom1])
    length2 <- unname(chrom_lengths[query$chrom2])
    coordinate_ok <- coordinate_ok & is.finite(length1) & is.finite(length2) &
      query$start1_0 < length1 & query$start2_0 < length2
  }
  query$query_valid <- chromosome_ok & coordinate_ok
  valid <- query[query$query_valid, , drop = FALSE]
  empty_plan <- data.frame(
    chrom1 = character(), tile1_start0 = numeric(), tile1_end0 = numeric(),
    chrom2 = character(), tile2_start0 = numeric(), tile2_end0 = numeric(),
    focus = character(), hic_loc1 = character(), hic_loc2 = character(),
    stringsAsFactors = FALSE
  )
  empty_plan$required_keys <- I(list())
  if (!nrow(valid)) {
    return(list(required = query, plan = empty_plan))
  }

  seqlevel_order <- unique(c(
    as.character(available_chromosomes),
    sort(unique(c(valid$chrom1, valid$chrom2)))
  ))
  rank1 <- match(valid$chrom1, seqlevel_order)
  rank2 <- match(valid$chrom2, seqlevel_order)
  tile1 <- floor(valid$start1_0 / tile_span) * tile_span
  tile2 <- floor(valid$start2_0 / tile_span) * tile_span
  swap <- rank1 > rank2 | (rank1 == rank2 & tile1 > tile2)
  axis1_chrom <- ifelse(swap, valid$chrom2, valid$chrom1)
  axis2_chrom <- ifelse(swap, valid$chrom1, valid$chrom2)
  axis1_start <- ifelse(swap, tile2, tile1)
  axis2_start <- ifelse(swap, tile1, tile2)

  tile_end <- function(chrom, start0) {
    end0 <- start0 + tile_span
    if (!is.null(chrom_lengths)) {
      chromosome_length <- unname(chrom_lengths[chrom])
      finite <- is.finite(chromosome_length)
      end0[finite] <- pmin(end0[finite], chromosome_length[finite])
    }
    end0
  }
  axis1_end <- tile_end(axis1_chrom, axis1_start)
  axis2_end <- tile_end(axis2_chrom, axis2_start)
  group_key <- paste(
    axis1_chrom, axis1_start, axis2_chrom, axis2_start, sep = "\u001f"
  )
  groups <- split(seq_len(nrow(valid)), group_key)
  rows <- lapply(groups, function(index) {
    chrom1 <- axis1_chrom[index[[1L]]]
    chrom2 <- axis2_chrom[index[[1L]]]
    start1 <- axis1_start[index[[1L]]]
    start2 <- axis2_start[index[[1L]]]
    end1 <- axis1_end[index[[1L]]]
    end2 <- axis2_end[index[[1L]]]
    range1 <- paste0(chrom1, ":", sprintf("%.0f", start1 + 1), "-",
                     sprintf("%.0f", end1))
    range2 <- paste0(chrom2, ":", sprintf("%.0f", start2 + 1), "-",
                     sprintf("%.0f", end2))
    focus <- if (identical(chrom1, chrom2) && identical(start1, start2)) {
      range1
    } else {
      paste(range1, range2, sep = "|")
    }
    out <- data.frame(
      chrom1 = chrom1,
      tile1_start0 = start1,
      tile1_end0 = end1,
      chrom2 = chrom2,
      tile2_start0 = start2,
      tile2_end0 = end2,
      focus = focus,
      hic_loc1 = paste0(
        chrom1, ":", sprintf("%.0f", start1), ":",
        sprintf("%.0f", end1 - 1)
      ),
      hic_loc2 = paste0(
        chrom2, ":", sprintf("%.0f", start2), ":",
        sprintf("%.0f", end2 - 1)
      ),
      stringsAsFactors = FALSE
    )
    out$required_keys <- I(list(unique(valid$.pixel_key[index])))
    out
  })
  plan <- do.call(rbind, rows)
  rownames(plan) <- NULL
  list(required = query, plan = plan)
}

.contact_canonicalize_pixels <- function(chrom1, start1_0, chrom2, start2_0) {
  first <- paste(chrom1, sprintf("%.0f", start1_0), sep = ":")
  second <- paste(chrom2, sprintf("%.0f", start2_0), sep = ":")
  swap <- first > second
  data.frame(
    chrom1 = ifelse(swap, chrom2, chrom1),
    start1_0 = ifelse(swap, start2_0, start1_0),
    chrom2 = ifelse(swap, chrom1, chrom2),
    start2_0 = ifelse(swap, start1_0, start2_0),
    stringsAsFactors = FALSE
  )
}

.contact_validate_hce_frame <- function(
    frame, score_names, resolution, chrom_lengths, required_keys, focus) {
  if (!nrow(frame)) {
    out <- data.frame(
      .pixel_key = character(), chrom1 = character(), start1_0 = numeric(),
      chrom2 = character(), start2_0 = numeric(),
      .query_focus = character(), stringsAsFactors = FALSE
    )
    for (layer in score_names) out[[layer]] <- numeric()
    return(out[c(
      ".pixel_key", "chrom1", "start1_0", "chrom2", "start2_0",
      score_names, ".query_focus"
    )])
  }
  .require_columns(
    frame,
    c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2",
      score_names),
    "native contact query"
  )
  numeric_fields <- c("start1", "end1", "start2", "end2")
  for (field in numeric_fields) frame[[field]] <- as.numeric(frame[[field]])
  start1_0 <- frame$start1 - 1
  start2_0 <- frame$start2 - 1
  width1 <- frame$end1 - frame$start1 + 1
  width2 <- frame$end2 - frame$start2 + 1
  bad <- !is.finite(start1_0) | !is.finite(start2_0) |
    !is.finite(width1) | !is.finite(width2) |
    start1_0 < 0 | start2_0 < 0 |
    start1_0 %% resolution != 0 | start2_0 %% resolution != 0 |
    width1 <= 0 | width2 <= 0 | width1 > resolution | width2 > resolution
  if (!is.null(chrom_lengths)) {
    length1 <- unname(chrom_lengths[as.character(frame$seqnames1)])
    length2 <- unname(chrom_lengths[as.character(frame$seqnames2)])
    bad <- bad | !is.finite(length1) | !is.finite(length2) |
      frame$end1 > length1 | frame$end2 > length2 |
      (width1 < resolution & frame$end1 != length1) |
      (width2 < resolution & frame$end2 != length2)
  }
  if (any(bad)) {
    .acv2g_validation_error(paste0(
      "Native contact query returned invalid or non-aligned fixed bins at ",
      focus, "."
    ))
  }
  coords <- .contact_canonicalize_pixels(
    as.character(frame$seqnames1), start1_0,
    as.character(frame$seqnames2), start2_0
  )
  out <- coords
  out$.pixel_key <- .contact_pixel_key(
    out$chrom1, out$start1_0, out$chrom2, out$start2_0
  )
  keep <- out$.pixel_key %in% required_keys
  out <- out[keep, , drop = FALSE]
  for (layer in score_names) {
    raw <- frame[[layer]][keep]
    value <- suppressWarnings(as.numeric(raw))
    if (any(is.na(value) & !is.na(raw)) || any(is.infinite(value))) {
      .acv2g_validation_error(paste0(
        "Native contact layer ", layer,
        " contains a non-numeric or infinite value at ", focus, "."
      ))
    }
    out[[layer]] <- value
  }
  out$.query_focus <- rep(focus, nrow(out))
  out[c(
    ".pixel_key", "chrom1", "start1_0", "chrom2", "start2_0",
    score_names, ".query_focus"
  )]
}

.contact_values_agree <- function(x) {
  if (all(is.na(x))) return(TRUE)
  if (anyNA(x)) return(FALSE)
  length(unique(x)) == 1L
}

.contact_collapse_native_pixels <- function(frames, score_names) {
  empty <- data.frame(
    .pixel_key = character(), chrom1 = character(), start1_0 = numeric(),
    chrom2 = character(), start2_0 = numeric(), stringsAsFactors = FALSE
  )
  for (layer in score_names) empty[[layer]] <- numeric()
  if (!length(frames)) return(empty)
  pixels <- do.call(rbind, frames)
  duplicated_key <- duplicated(pixels$.pixel_key) |
    duplicated(pixels$.pixel_key, fromLast = TRUE)
  if (any(duplicated_key)) {
    groups <- split(
      which(duplicated_key), pixels$.pixel_key[duplicated_key]
    )
    for (key in names(groups)) {
      index <- which(pixels$.pixel_key == key)
      for (layer in score_names) {
        if (!.contact_values_agree(pixels[[layer]][index])) {
          foci <- paste(unique(pixels$.query_focus[index]), collapse = ", ")
          .acv2g_validation_error(paste0(
            "Conflicting native contact values for pixel ", key,
            ", layer ", layer, ", across query focus/foci: ", foci, "."
          ))
        }
      }
    }
  }
  pixels <- pixels[!duplicated(pixels$.pixel_key), , drop = FALSE]
  pixels$.query_focus <- NULL
  rownames(pixels) <- NULL
  pixels[c(
    ".pixel_key", "chrom1", "start1_0", "chrom2", "start2_0",
    score_names
  )]
}

.contact_extract_hce_weights <- function(hce, resolution) {
  bins <- tryCatch(
    HiCExperiment::bins(hce),
    error = function(e) {
      .acv2g_validation_error(paste0(
        "Could not inspect native Cooler bin weights: ", conditionMessage(e)
      ))
    }
  )
  if (!"weight" %in% names(S4Vectors::mcols(bins))) {
    return(list(available = FALSE, weights = NULL))
  }
  out <- data.frame(
    chrom = as.character(GenomeInfoDb::seqnames(bins)),
    start0 = as.numeric(GenomicRanges::start(bins)) - 1,
    weight = as.numeric(S4Vectors::mcols(bins)$weight),
    stringsAsFactors = FALSE
  )
  bad <- !is.finite(out$start0) | out$start0 < 0 |
    out$start0 %% resolution != 0 | is.infinite(out$weight)
  if (any(bad)) {
    .acv2g_validation_error(
      "Native Cooler contains invalid bin coordinates or infinite weights."
    )
  }
  list(available = TRUE, weights = out)
}

.contact_collapse_weights <- function(frames) {
  if (!length(frames)) return(NULL)
  weights <- do.call(rbind, frames)
  key <- paste(weights$chrom, weights$start0, sep = "\u001f")
  duplicate <- duplicated(key) | duplicated(key, fromLast = TRUE)
  if (any(duplicate)) {
    for (value_key in unique(key[duplicate])) {
      if (!.contact_values_agree(weights$weight[key == value_key])) {
        .acv2g_validation_error(paste0(
          "Conflicting native balancing weights for bin ", value_key, "."
        ))
      }
    }
  }
  weights <- weights[!duplicated(key), , drop = FALSE]
  rownames(weights) <- NULL
  weights
}

.contact_import_hce_focus <- function(native, source_type, focus, resolution) {
  imported <- NULL
  invisible(utils::capture.output(
    imported <- suppressMessages(
      if (identical(source_type, "cool")) {
        HiCExperiment::import(native, focus = focus)
      } else {
        HiCExperiment::import(native, focus = focus, resolution = resolution)
      }
    )
  ))
  imported
}

.contact_import_cool_tiles <- function(
    native, path, source_type, plan, resolution,
    available_chromosomes, chrom_lengths) {
  work_plan <- plan$plan
  inspection_only <- !nrow(work_plan)
  if (inspection_only) {
    if (!length(available_chromosomes)) {
      .acv2g_validation_error(
        "No query chromosome is available in the native contact source."
      )
    }
    chrom <- available_chromosomes[[1L]]
    end0 <- resolution
    if (!is.null(chrom_lengths) && is.finite(chrom_lengths[chrom])) {
      end0 <- min(end0, unname(chrom_lengths[chrom]))
    }
    work_plan <- data.frame(
      chrom1 = chrom, tile1_start0 = 0, tile1_end0 = end0,
      chrom2 = chrom, tile2_start0 = 0, tile2_end0 = end0,
      focus = paste0(chrom, ":1-", sprintf("%.0f", end0)),
      hic_loc1 = NA_character_, hic_loc2 = NA_character_,
      stringsAsFactors = FALSE
    )
    work_plan$required_keys <- I(list(character()))
  }

  needed_bins <- unique(c(
    paste(
      plan$required$chrom1[plan$required$query_valid],
      plan$required$start1_0[plan$required$query_valid], sep = "\u001f"
    ),
    paste(
      plan$required$chrom2[plan$required$query_valid],
      plan$required$start2_0[plan$required$query_valid], sep = "\u001f"
    )
  ))
  frames <- list()
  weight_frames <- list()
  score_names <- NULL
  weights_checked <- FALSE
  weights_available <- FALSE
  weight_coverage_complete <- FALSE
  queried_keys <- character()
  query_failures <- list()

  for (i in seq_len(nrow(work_plan))) {
    focus <- work_plan$focus[[i]]
    hce <- tryCatch(
      .contact_import_hce_focus(native, source_type, focus, resolution),
      error = function(e) e
    )
    if (inherits(hce, "error")) {
      query_failures[[length(query_failures) + 1L]] <-
        .contact_query_failure_row(
          "HiCExperiment", focus, NULL, hce
        )
      next
    }
    current_names <- names(HiCExperiment::scores(hce))
    if (is.null(current_names)) current_names <- character()
    if (is.null(score_names)) {
      score_names <- current_names
    } else if (!setequal(score_names, current_names)) {
      .acv2g_validation_error(
        "Native contact layers differ across queried regions."
      )
    }
    required_keys <- work_plan$required_keys[[i]]
    frame <- .contact_validate_hce_frame(
      as.data.frame(HiCExperiment::interactions(hce)),
      score_names, resolution, chrom_lengths,
      required_keys, focus
    )
    if (nrow(frame)) frames[[length(frames) + 1L]] <- frame
    if (!weights_checked || (weights_available && !weight_coverage_complete)) {
      extracted <- .contact_extract_hce_weights(hce, resolution)
      if (!weights_checked) {
        weights_available <- extracted$available
        weights_checked <- TRUE
      } else if (!identical(weights_available, extracted$available)) {
        .acv2g_validation_error(
          "Native Cooler weight availability differs across queried regions."
        )
      }
      if (isTRUE(extracted$available)) {
        weight_frames[[length(weight_frames) + 1L]] <- extracted$weights
        collapsed <- .contact_collapse_weights(weight_frames)
        observed <- paste(collapsed$chrom, collapsed$start0, sep = "\u001f")
        weight_coverage_complete <- all(needed_bins %in% observed)
      }
    }
    if (!inspection_only) {
      queried_keys <- c(queried_keys, required_keys)
    }
  }
  failures <- if (length(query_failures)) {
    do.call(rbind, query_failures)
  } else {
    .contact_empty_query_failures()
  }
  if (is.null(score_names)) {
    examples <- paste(utils::head(failures$focus, 3L), collapse = "; ")
    .acv2g_validation_error(paste0(
      "No native Cooler query succeeded, so contact layers could not be ",
      "validated. Failed focus/foci: ", examples, "."
    ), data = failures)
  }
  list(
    pixels = .contact_collapse_native_pixels(frames, score_names),
    score_names = score_names,
    bin_weights = if (weights_available) {
      .contact_collapse_weights(weight_frames)
    } else {
      NULL
    },
    weights_available = weights_available,
    queried_keys = unique(queried_keys),
    query_failures = failures
  )
}

.contact_straw_query <- function(
    normalization, path, chr1loc, chr2loc, resolution) {
  strawr::straw(
    normalization, path, chr1loc, chr2loc, "BP", resolution,
    matrix = "observed"
  )
}

.contact_validate_straw_layer <- function(
    frame, layer, chrom1, chrom2, resolution,
    chrom_lengths, required_keys, focus) {
  empty <- data.frame(
    .pixel_key = character(), chrom1 = character(), start1_0 = numeric(),
    chrom2 = character(), start2_0 = numeric(), value = numeric(),
    stringsAsFactors = FALSE
  )
  if (!nrow(frame)) return(empty)
  .require_columns(frame, c("x", "y", "counts"), "straw contact query")
  start1 <- suppressWarnings(as.numeric(frame$x))
  start2 <- suppressWarnings(as.numeric(frame$y))
  value <- suppressWarnings(as.numeric(frame$counts))
  bad <- !is.finite(start1) | !is.finite(start2) |
    start1 < 0 | start2 < 0 |
    start1 %% resolution != 0 | start2 %% resolution != 0 |
    is.infinite(value)
  if (!is.null(chrom_lengths)) {
    bad <- bad | start1 >= unname(chrom_lengths[chrom1]) |
      start2 >= unname(chrom_lengths[chrom2])
  }
  if (any(bad)) {
    .acv2g_validation_error(paste0(
      "straw returned invalid bins/values for layer ", layer,
      " at ", focus, "."
    ))
  }
  coords <- .contact_canonicalize_pixels(
    rep(chrom1, length(start1)), start1,
    rep(chrom2, length(start2)), start2
  )
  out <- coords
  out$.pixel_key <- .contact_pixel_key(
    out$chrom1, out$start1_0, out$chrom2, out$start2_0
  )
  out$value <- value
  out <- out[out$.pixel_key %in% required_keys, , drop = FALSE]
  duplicate <- duplicated(out$.pixel_key) |
    duplicated(out$.pixel_key, fromLast = TRUE)
  if (any(duplicate)) {
    for (key in unique(out$.pixel_key[duplicate])) {
      if (!.contact_values_agree(out$value[out$.pixel_key == key])) {
        .acv2g_validation_error(paste0(
          "straw returned conflicting values for pixel ", key,
          ", layer ", layer, "."
        ))
      }
    }
    out <- out[!duplicated(out$.pixel_key), , drop = FALSE]
  }
  rownames(out) <- NULL
  out[c(
    ".pixel_key", "chrom1", "start1_0", "chrom2", "start2_0", "value"
  )]
}

.contact_import_hic_tiles <- function(
    path, plan, resolution, chrom_lengths, contact_value, contact_layer) {
  norm_types <- tryCatch(
    as.character(strawr::readHicNormTypes(path)),
    error = function(e) {
      .acv2g_validation_error(paste0(
        "Could not inspect .hic normalization layers in ", path, ": ",
        conditionMessage(e)
      ))
    }
  )
  norm_types <- unique(toupper(norm_types))
  balanced_norms <- setdiff(norm_types, "NONE")
  available_score_names <- c("count", balanced_norms)
  layer_roles <- .contact_native_layer_roles(
    available_score_names, "hic", weights_available = FALSE
  )
  selected <- .contact_resolve_layer(
    layer_roles, contact_value = contact_value, contact_layer = contact_layer
  )
  score_names <- selected$layer
  frames <- list()
  queried_keys <- character()
  query_failures <- list()
  if (!nrow(plan$plan)) {
    return(list(
      pixels = .contact_collapse_native_pixels(frames, score_names),
      score_names = score_names,
      bin_weights = NULL,
      weights_available = FALSE,
      queried_keys = queried_keys,
      layer_roles = layer_roles,
      selected = selected,
      query_failures = .contact_empty_query_failures()
    ))
  }

  for (i in seq_len(nrow(plan$plan))) {
    row <- plan$plan[i, , drop = FALSE]
    required_keys <- row$required_keys[[1L]]
    layer_frames <- vector("list", length(score_names))
    names(layer_frames) <- score_names
    for (layer in score_names) {
      normalization <- if (identical(layer, "count")) "NONE" else layer
      raw <- tryCatch(
        .contact_straw_query(
          normalization, path, row$hic_loc1[[1L]], row$hic_loc2[[1L]],
          resolution
        ),
        error = function(e) e
      )
      if (inherits(raw, "error")) {
        query_failures[[length(query_failures) + 1L]] <-
          .contact_query_failure_row(
            "strawr", row$focus[[1L]], normalization, raw
          )
        layer_frames <- NULL
        break
      }
      layer_frames[[layer]] <- .contact_validate_straw_layer(
        raw, layer, row$chrom1[[1L]], row$chrom2[[1L]], resolution,
        chrom_lengths, required_keys, row$focus[[1L]]
      )
    }
    if (is.null(layer_frames)) next
    nonempty <- layer_frames[vapply(layer_frames, nrow, integer(1L)) > 0L]
    if (length(nonempty)) {
      coordinates <- do.call(rbind, lapply(nonempty, function(x) {
        x[c(".pixel_key", "chrom1", "start1_0", "chrom2", "start2_0")]
      }))
      coordinates <- coordinates[!duplicated(coordinates$.pixel_key), , drop = FALSE]
      tile <- coordinates
      for (layer in score_names) {
        source <- layer_frames[[layer]]
        tile[[layer]] <- source$value[match(tile$.pixel_key, source$.pixel_key)]
      }
      tile$.query_focus <- rep(row$focus[[1L]], nrow(tile))
      frames[[length(frames) + 1L]] <- tile
    }
    queried_keys <- c(queried_keys, required_keys)
  }
  list(
    pixels = .contact_collapse_native_pixels(frames, score_names),
    score_names = score_names,
    bin_weights = NULL,
    weights_available = FALSE,
    queried_keys = unique(queried_keys),
    layer_roles = layer_roles,
    selected = selected,
    query_failures = if (length(query_failures)) {
      do.call(rbind, query_failures)
    } else {
      .contact_empty_query_failures()
    }
  )
}

.contact_native_layer_roles <- function(
    score_names, source_type, weights_available) {
  lower <- tolower(score_names)
  raw <- lower %in% c("count", "counts", "raw")
  balanced <- if (identical(source_type, "hic")) {
    !raw
  } else {
    lower == "balanced" & isTRUE(weights_available)
  }
  named <- !raw & !balanced
  if (source_type %in% c("cool", "mcool") && !isTRUE(weights_available)) {
    named <- named & lower != "balanced"
  }
  rbind(
    data.frame(
      layer = score_names[raw], type = rep("raw", sum(raw)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      layer = score_names[balanced], type = rep("balanced", sum(balanced)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      layer = score_names[named], type = rep("named_score", sum(named)),
      stringsAsFactors = FALSE
    )
  )
}

.contact_validate_native_chromosomes <- function(
    chromosomes, lengths, path) {
  chromosomes <- as.character(chromosomes)
  bad_chromosome <- is.na(chromosomes) | !nzchar(chromosomes)
  if (!length(chromosomes) || any(bad_chromosome) ||
      anyDuplicated(chromosomes)) {
    .acv2g_validation_error(paste0(
      "Native contact source ", path,
      " does not expose a unique, non-empty chromosome list."
    ))
  }
  if (is.null(lengths) || is.null(names(lengths))) {
    .acv2g_validation_error(paste0(
      "Native contact source ", path,
      " does not expose named chromosome lengths; ACV2G cannot safely ",
      "distinguish structural zero from an out-of-bounds query."
    ))
  }
  lengths <- suppressWarnings(as.numeric(lengths[chromosomes]))
  bad_length <- !is.finite(lengths) | lengths <= 0 |
    lengths != floor(lengths)
  if (any(bad_length)) {
    .acv2g_validation_error(paste0(
      "Native contact source ", path,
      " has missing or invalid chromosome lengths for: ",
      paste(chromosomes[bad_length], collapse = ", "), "."
    ))
  }
  names(lengths) <- chromosomes
  list(chromosomes = chromosomes, lengths = lengths)
}

.contact_prepare_native <- function(
    pairs, path, source_type, requested_resolution, diagonal_policy,
    contact_value, contact_layer) {
  input_path <- path
  expanded_path <- path.expand(path)
  native <- .contact_native_file(expanded_path, source_type)
  path <- normalizePath(expanded_path, mustWork = TRUE)
  if (identical(source_type, "hic")) {
    available_resolutions <- tryCatch(
      strawr::readHicBpResolutions(path),
      error = function(e) {
        .acv2g_validation_error(paste0(
          "Could not inspect native contact resolutions in ", path, ": ",
          conditionMessage(e)
        ))
      }
    )
    chromosome_table <- tryCatch(
      strawr::readHicChroms(path),
      error = function(e) {
        .acv2g_validation_error(paste0(
          "Could not inspect native contact chromosomes in ", path, ": ",
          conditionMessage(e)
        ))
      }
    )
    chromosome_table <- chromosome_table[
      toupper(as.character(chromosome_table$name)) != "ALL", , drop = FALSE
    ]
    available_chromosomes <- as.character(chromosome_table$name)
    chrom_lengths <- as.numeric(chromosome_table$length)
    names(chrom_lengths) <- available_chromosomes
  } else {
    available_resolutions <- tryCatch(
      .contact_quiet_native_call(
        HiCExperiment::availableResolutions(native)
      ),
      error = function(e) {
        .acv2g_validation_error(paste0(
          "Could not inspect native contact resolutions in ", path, ": ",
          conditionMessage(e)
        ))
      }
    )
    available_info <- tryCatch({
      available <- .contact_quiet_native_call(
        HiCExperiment::availableChromosomes(native)
      )
      if (inherits(available, "Seqinfo")) {
        list(
          chromosomes = GenomeInfoDb::seqlevels(available),
          lengths = GenomeInfoDb::seqlengths(available)
        )
      } else {
        list(chromosomes = as.character(available), lengths = NULL)
      }
    }, error = function(e) {
      .acv2g_validation_error(paste0(
        "Could not inspect native contact chromosomes in ", path, ": ",
        conditionMessage(e)
      ))
    })
    available_chromosomes <- available_info$chromosomes
    chrom_lengths <- available_info$lengths
  }
  chromosome_info <- .contact_validate_native_chromosomes(
    available_chromosomes, chrom_lengths, path
  )
  available_chromosomes <- chromosome_info$chromosomes
  chrom_lengths <- chromosome_info$lengths
  resolution <- .contact_resolution(
    available_resolutions, requested_resolution, native = TRUE
  )
  query <- .contact_native_query_rows(
    pairs, resolution, diagonal_policy, chrom_lengths
  )
  plan <- .contact_native_query_plan(
    query, resolution, available_chromosomes, chrom_lengths
  )
  imported <- if (identical(source_type, "hic")) {
    .contact_import_hic_tiles(
      path, plan, resolution, chrom_lengths, contact_value, contact_layer
    )
  } else {
    .contact_import_cool_tiles(
      native, path, source_type, plan, resolution,
      available_chromosomes, chrom_lengths
    )
  }
  layer_roles <- imported$layer_roles %||% .contact_native_layer_roles(
    imported$score_names, source_type, imported$weights_available
  )
  query_failures <- imported$query_failures %||%
    .contact_empty_query_failures()
  if (nrow(query_failures)) {
    .acv2g_warn(
      paste0(
        nrow(query_failures),
        " native contact tile query/queries failed; affected pixels are ",
        "retained as missing_contact_query."
      ),
      data = list(failures = utils::head(query_failures, 5L)),
      class = "acv2g_contact_query_warning"
    )
  }
  list(
    pixels = imported$pixels,
    layer_roles = layer_roles,
    selected = imported$selected %||% NULL,
    resolution = resolution,
    implicit_zero = TRUE,
    query_succeeded_keys = imported$queried_keys,
    available_chromosomes = available_chromosomes,
    chrom_lengths = chrom_lengths,
    bin_weights = imported$bin_weights,
    native_type = source_type,
    source_path_input = input_path,
    source_path = path,
    query_failures = query_failures
  )
}

.contact_resolve_layer <- function(layer_roles, contact_value, contact_layer) {
  available <- as.character(layer_roles$layer)
  available_text <- if (length(available)) paste(available, collapse = ", ") else "<none>"
  if (!is.null(contact_layer)) {
    if (!contact_layer %in% available) {
      .acv2g_validation_error(
        "Explicit contact layer {.field {contact_layer}} is unavailable. Available layers: {available_text}."
      )
    }
    type <- layer_roles$type[match(contact_layer, layer_roles$layer)]
    if (contact_value != "named_score" && type != contact_value) {
      .acv2g_validation_error(
        "Explicit layer {.field {contact_layer}} is declared {.val {type}}, not requested {.val {contact_value}}."
      )
    }
    return(list(layer = contact_layer, type = type, fallback_used = FALSE))
  }
  if (contact_value == "named_score") {
    .acv2g_validation_error(
      "{.arg contact_layer} is required when {.arg contact_value = 'named_score'}. Available layers: {available_text}."
    )
  }
  candidates <- layer_roles$layer[layer_roles$type == contact_value]
  if (length(candidates) > 1L) {
    .acv2g_validation_error(
      "Multiple {contact_value} contact layers are available: {paste(candidates, collapse = ', ')}. Select one with {.arg contact_layer}."
    )
  }
  if (length(candidates) == 1L) {
    return(list(layer = candidates, type = contact_value, fallback_used = FALSE))
  }
  if (contact_value != "balanced") {
    .acv2g_validation_error(
      "No raw contact layer is available. Available layers: {available_text}."
    )
  }
  raw <- layer_roles$layer[layer_roles$type == "raw"]
  if (length(raw) != 1L) {
    reason <- if (!length(raw)) "No raw layer is available" else
      "Multiple raw layers are available"
    .acv2g_validation_error(
      "No balanced contact layer is available and raw fallback cannot be resolved. {reason}. Available layers: {available_text}."
    )
  }
  .acv2g_warn(
    "No balanced contact layer is available; using raw layer {.field {raw}}.",
    class = "acv2g_contact_raw_fallback"
  )
  list(layer = raw, type = "raw", fallback_used = TRUE)
}

.contact_pixel_key <- function(chrom1, start1_0, chrom2, start2_0) {
  first <- paste(chrom1, sprintf("%.0f", start1_0), sep = ":")
  second <- paste(chrom2, sprintf("%.0f", start2_0), sep = ":")
  swap <- first > second
  paste(ifelse(swap, second, first), ifelse(swap, first, second), sep = "\u001e")
}

.contact_vector_stable_row_id <- function(..., prefix) {
  pieces <- list(...)
  n <- max(vapply(pieces, length, integer(1L)), 0L)
  pieces <- lapply(pieces, rep_len, length.out = n)
  key <- do.call(paste, c(pieces, sep = "\u001f"))
  hashes <- if ("getVDigest" %in% getNamespaceExports("digest")) {
    digest::getVDigest("xxhash64")(key, serialize = FALSE)
  } else {
    vapply(
      key, digest::digest, character(1L),
      algo = "xxhash64", serialize = FALSE
    )
  }
  paste0(prefix, "_", unname(hashes))
}

.contact_make_tss_bin_seed <- function(
    pairs, resolution, chrom_lengths = NULL) {
  pair_columns <- .contact_pair_columns(pairs)
  pair_key <- .contact_pair_key(pairs)
  variant_bin <- floor((pairs$pos1 - 1) / resolution) * resolution
  promoter_bin <- floor((pairs$tss_pos1 - 1) / resolution) * resolution
  x <- pairs
  bin_key <- paste(pair_key, x$tss_chrom, promoter_bin, sep = "\u001f")

  # Most pair/TSS-bin groups contain a single row. Resolve those groups with
  # direct vector indexing and use data.table grouping only for bins that
  # actually contain duplicate or alternative TSS records.
  bin_levels <- sort(unique(bin_key))
  group_index <- match(bin_key, bin_levels)
  group_size <- tabulate(group_index, nbins = length(bin_levels))
  first_index <- match(seq_along(bin_levels), group_index)

  contact_tss_ids <- as.character(x$tss_id[first_index])
  contact_tss_positions1 <- vapply(
    first_index,
    function(i) paste(sort(unique(x$tss_pos1[i])), collapse = ";"),
    character(1L)
  )
  contact_tss_n_in_bin <- rep.int(1L, length(bin_levels))
  contact_distance_bp <- abs(
    x$tss_pos1[first_index] - x$pos1[first_index]
  )

  duplicate_groups <- which(group_size > 1L)
  if (length(duplicate_groups)) {
    duplicate_row <- group_index %in% duplicate_groups
    duplicate_data <- data.table::as.data.table(
      x[duplicate_row, c("tss_id", "tss_pos1", "pos1"), drop = FALSE]
    )
    data.table::set(
      duplicate_data,
      j = ".contact_group",
      value = group_index[duplicate_row]
    )
    duplicate_summary <- duplicate_data[, {
      ids <- get("tss_id")
      positions <- get("tss_pos1")
      variant_positions <- get("pos1")
      list(
        contact_tss_ids = paste(sort(unique(ids)), collapse = ";"),
        contact_tss_positions1 = paste(
          sort(unique(positions)), collapse = ";"
        ),
        contact_tss_n_in_bin = length(unique(ids)),
        contact_distance_bp = min(abs(positions - variant_positions))
      )
    }, by = ".contact_group"]
    replace_index <- duplicate_summary[[".contact_group"]]
    contact_tss_ids[replace_index] <-
      duplicate_summary[["contact_tss_ids"]]
    contact_tss_positions1[replace_index] <-
      duplicate_summary[["contact_tss_positions1"]]
    contact_tss_n_in_bin[replace_index] <- as.integer(
      duplicate_summary[["contact_tss_n_in_bin"]]
    )
    contact_distance_bp[replace_index] <-
      duplicate_summary[["contact_distance_bp"]]
  }

  seed <- x[first_index, pair_columns, drop = FALSE]
  seed$chrom <- x$chrom[first_index]
  seed$tss_chrom <- x$tss_chrom[first_index]
  seed$tss_id <- contact_tss_ids
  seed$contact_tss_ids <- contact_tss_ids
  seed$contact_tss_positions1 <- contact_tss_positions1
  seed$contact_tss_n_in_bin <- contact_tss_n_in_bin
  seed$contact_variant_bin_start0 <- variant_bin[first_index]
  seed$contact_promoter_bin_start0 <- promoter_bin[first_index]
  seed$contact_variant_bin_end0 <- variant_bin[first_index] + resolution
  seed$contact_promoter_bin_end0 <- promoter_bin[first_index] + resolution
  seed$contact_distance_bp <- contact_distance_bp
  seed$tss_bin_id <- .contact_vector_stable_row_id(
    seed$pair_id, seed$tss_chrom, seed$contact_promoter_bin_start0,
    resolution, prefix = "tss_bin"
  )
  seed$.pair_key <- .contact_pair_key(seed)
  seed <- seed[order(
    seed$.pair_key, seed$tss_chrom, seed$contact_promoter_bin_start0,
    seed$contact_tss_ids
  ), , drop = FALSE]
  rownames(seed) <- NULL
  seed
}

.contact_lookup_pixels <- function(
    chrom1, start1_0, chrom2, start2_0, prepared) {
  key <- .contact_pixel_key(chrom1, start1_0, chrom2, start2_0)
  match_index <- match(key, prepared$pixels$.pixel_key)
  value <- rep(NA_real_, length(key))
  found <- !is.na(match_index)
  value[found] <- prepared$pixels$.selected_contact_value[match_index[found]]
  chromosome_ok <- chrom1 %in% prepared$available_chromosomes &
    chrom2 %in% prepared$available_chromosomes
  coordinate_ok <- rep(TRUE, length(key))
  if (!is.null(prepared$chrom_lengths)) {
    length1 <- unname(prepared$chrom_lengths[chrom1])
    length2 <- unname(prepared$chrom_lengths[chrom2])
    coordinate_ok <- is.finite(length1) & is.finite(length2) &
      start1_0 >= 0 & start2_0 >= 0 &
      start1_0 < length1 & start2_0 < length2
  }
  status <- rep("missing_contact_pixel", length(key))
  status[!chromosome_ok] <- "missing_contact_seqlevel"
  status[chromosome_ok & !coordinate_ok] <- "outside_contact_bounds"
  status[found & is.finite(value)] <- ifelse(
    value[found & is.finite(value)] == 0, "observed_zero", "observed"
  )
  status[found & is.na(value)] <- "missing_contact_value"

  if (prepared$implicit_zero) {
    absent_valid <- !found & chromosome_ok & coordinate_ok
    query_succeeded <- rep(TRUE, length(key))
    if (!is.null(prepared$query_succeeded_keys)) {
      query_succeeded <- key %in% prepared$query_succeeded_keys
      status[absent_valid & !query_succeeded] <- "missing_contact_query"
      absent_valid <- absent_valid & query_succeeded
    }
    if (identical(prepared$selected_type, "balanced") &&
        !is.null(prepared$bin_weights)) {
      weight_key <- paste(prepared$bin_weights$chrom,
                          prepared$bin_weights$start0, sep = "\u001f")
      i1 <- match(paste(chrom1, start1_0, sep = "\u001f"), weight_key)
      i2 <- match(paste(chrom2, start2_0, sep = "\u001f"), weight_key)
      missing_weight <- is.na(i1) | is.na(i2) |
        is.na(prepared$bin_weights$weight[i1]) |
        is.na(prepared$bin_weights$weight[i2])
      missing_weight <- missing_weight & query_succeeded
      value[missing_weight & chromosome_ok & coordinate_ok] <- NA_real_
      status[missing_weight & chromosome_ok & coordinate_ok] <-
        "missing_balance_weight"
      absent_valid <- absent_valid & !missing_weight
    }
    if (identical(prepared$selected_type, "balanced") &&
        identical(prepared$native_type, "hic")) {
      uncertain <- chromosome_ok & coordinate_ok &
        ((!found & absent_valid) | (found & is.na(value)))
      value[uncertain] <- NA_real_
      status[uncertain] <- "missing_balance_weight"
      absent_valid <- absent_valid & !uncertain
    }
    value[absent_valid] <- 0
    status[absent_valid] <- "implicit_structural_zero"
  }
  data.frame(value = value, status = status, stringsAsFactors = FALSE)
}

.contact_apply_diagonal <- function(
    seed, direct, prepared, diagonal_policy, resolution) {
  detail <- seed
  is_diagonal <- detail$chrom == detail$tss_chrom &
    detail$contact_variant_bin_start0 ==
    detail$contact_promoter_bin_start0
  detail$contact_is_diagonal <- is_diagonal
  detail$contact_diagonal_observed <- ifelse(is_diagonal, direct$value, NA_real_)
  detail$contact_diagonal_replacement <- NA_real_
  detail$contact_diagonal_status <- ifelse(is_diagonal, NA_character_, "not_diagonal")
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
    left[valid_left, ] <- .contact_lookup_pixels(
      detail$tss_chrom[index[valid_left]],
      detail$contact_variant_bin_start0[index[valid_left]],
      detail$tss_chrom[index[valid_left]], left_start[valid_left], prepared
    )
  }
  right <- .contact_lookup_pixels(
    detail$tss_chrom[index], detail$contact_variant_bin_start0[index],
    detail$tss_chrom[index],
    detail$contact_promoter_bin_start0[index] + resolution, prepared
  )
  left_finite <- is.finite(left$value)
  right_finite <- is.finite(right$value)
  finite_n <- as.integer(left_finite) + as.integer(right_finite)
  has_finite <- finite_n > 0L
  replacement <- rep(NA_real_, length(index))
  replacement[left_finite & right_finite] <- pmax(
    left$value[left_finite & right_finite],
    right$value[left_finite & right_finite]
  )
  replacement[left_finite & !right_finite] <-
    left$value[left_finite & !right_finite]
  replacement[!left_finite & right_finite] <-
    right$value[!left_finite & right_finite]

  detail$contact_raw[index] <- replacement
  detail$contact_diagonal_replacement[index[has_finite]] <-
    replacement[has_finite]
  detail$contact_status[index[has_finite]] <- ifelse(
    replacement[has_finite] == 0,
    "diagonal_neighbor_zero",
    "diagonal_neighbor_observed"
  )
  detail$contact_status[index[!has_finite]] <- "missing_diagonal_neighbors"
  detail$contact_diagonal_status[index[finite_n == 2L]] <-
    "neighbor_max_two_finite"
  detail$contact_diagonal_status[index[finite_n == 1L]] <-
    "neighbor_max_one_finite"
  detail$contact_diagonal_status[index[finite_n == 0L]] <-
    "neighbor_max_no_finite"
  detail
}

.contact_aggregate_tss_bins <- function(detail, original_pairs, transform) {
  pair_columns <- .contact_pair_columns(original_pairs)
  pair_levels <- sort(unique(detail$.pair_key))
  group_index <- match(detail$.pair_key, pair_levels)
  group_n <- length(pair_levels)
  first_index <- match(seq_len(group_n), group_index)

  original_pair_key <- .contact_pair_key(original_pairs)
  eligible_identity <- paste(
    original_pairs$tss_id, original_pairs$tss_chrom,
    original_pairs$tss_pos1, sep = "\u001f"
  )
  eligible_unique <- unique(data.table::data.table(
    .pair_key = original_pair_key,
    .tss_identity = eligible_identity
  ))
  eligible_count <- tabulate(
    match(eligible_unique[[".pair_key"]], pair_levels),
    nbins = group_n
  )

  # Retain the original validation that every pair resolves to exactly one
  # variant bin, including the treatment of repeated missing values.
  unique_variant_bins <- unique(data.table::data.table(
    .contact_group = group_index,
    .variant_bin = detail$contact_variant_bin_start0
  ))
  variant_bin_n <- tabulate(
    unique_variant_bins[[".contact_group"]], nbins = group_n
  )
  if (any(variant_bin_n != 1L)) {
    .acv2g_validation_error(
      "One variant-gene pair resolved to multiple variant contact bins."
    )
  }

  values <- detail$contact_raw
  finite <- is.finite(values)
  zero <- finite & values == 0
  n_total <- tabulate(group_index, nbins = group_n)
  n_finite <- tabulate(group_index[finite], nbins = group_n)
  n_zero <- tabulate(group_index[zero], nbins = group_n)
  n_missing <- n_total - n_finite
  has_finite <- n_finite > 0L

  contact_raw <- rep(NA_real_, group_n)
  if (any(finite)) {
    finite_data <- data.table::data.table(
      .contact_group = group_index[finite],
      .contact_value = values[finite]
    )
    maxima <- finite_data[, list(
      contact_raw = max(get(".contact_value"))
    ), by = ".contact_group"]
    contact_raw[maxima[[".contact_group"]]] <- maxima[["contact_raw"]]
  }

  contact_is_max <- finite & values == contact_raw[group_index]
  max_tie_n <- tabulate(group_index[contact_is_max], nbins = group_n)
  selected_index <- rep(NA_integer_, group_n)
  ties <- which(contact_is_max)
  if (length(ties)) {
    tie_order <- order(
      group_index[ties], detail$tss_chrom[ties],
      detail$contact_promoter_bin_start0[ties],
      detail$contact_tss_ids[ties]
    )
    ordered_ties <- ties[tie_order]
    ordered_groups <- group_index[ordered_ties]
    selected_index[ordered_groups[!duplicated(ordered_groups)]] <-
      ordered_ties[!duplicated(ordered_groups)]
  }

  out <- detail[first_index, pair_columns, drop = FALSE]
  out$contact_raw <- contact_raw
  out$contact_transformed <- if (transform == "log1p") {
    log1p(contact_raw)
  } else {
    contact_raw
  }
  out$contact_tss_n_eligible <- eligible_count
  out$contact_tss_bin_n_total <- n_total
  out$contact_tss_bin_n_finite <- n_finite
  out$contact_tss_bin_n_zero <- n_zero
  out$contact_tss_bin_n_missing <- n_missing
  out$contact_tss_bin_missing_fraction <- n_missing / n_total
  out$contact_variant_bin_start0 <-
    detail$contact_variant_bin_start0[first_index]
  out$contact_max_tie_n <- max_tie_n
  out$contact_selected_tss_id <- detail$contact_tss_ids[selected_index]
  out$contact_selected_tss_bin_start0 <-
    detail$contact_promoter_bin_start0[selected_index]
  out$contact_promoter_bin_start0 <-
    detail$contact_promoter_bin_start0[selected_index]
  out$contact_distance_bp <- detail$contact_distance_bp[selected_index]
  out$contact_is_diagonal <- detail$contact_is_diagonal[selected_index]
  out$contact_diagonal_observed <-
    detail$contact_diagonal_observed[selected_index]
  out$contact_diagonal_replacement <-
    detail$contact_diagonal_replacement[selected_index]
  out$contact_diagonal_status <- detail$contact_diagonal_status[selected_index]
  out$contact_diagonal_status[!has_finite] <- "all_tss_bins_missing"
  out$contact_aggregation_status <- ifelse(
    !has_finite,
    "all_tss_bins_missing",
    ifelse(n_missing > 0L, "finite_max_with_partial_missing", "finite_max_complete")
  )
  out$contact_status <- ifelse(
    !has_finite,
    "missing_contact_all_tss_bins",
    ifelse(contact_raw == 0, "observed_zero", "observed")
  )
  rownames(out) <- NULL
  attr(out, "detail_contact_is_max") <- contact_is_max
  out
}

#' @export
print.acv2g_contact_result <- function(x, ...) {
  cat(
    "<acv2g_contact_result>", nrow(x$pairs), "variant-gene pairs;",
    nrow(x$tss_bin_contacts), "unique TSS-bin contacts\n"
  )
  invisible(x)
}
