#' Extract variant-centred regulatory activity
#'
#' `extract_v2g_activity()` measures signal in an exact, 0-based half-open
#' window around each variant anchor.  The default mean uses the complete
#' effective window as its denominator, so bases that are not represented in a
#' valid BigWig or bedGraph track contribute zero.  A missing chromosome,
#' unreadable value, or region-query failure is retained as missing evidence
#' rather than being changed to zero. Query failures emit a warning and are
#' recorded in row-level status/message fields; invalid source metadata or
#' schema still raises a validation error.
#'
#' BigWig and bedGraph files are read with `rtracklayer::import()`.  Indexed,
#' coordinate-sorted BAM files use a single fixed assay: filtered fragment/read
#' coverage in counts per million (CPM).  Properly paired templates are counted
#' once and unpaired alignments are counted as reads.  Secondary,
#' supplementary, QC-failed, duplicate, unmapped, and low-MAPQ records are
#' excluded.  CRAM is deliberately not supported in the first release.
#'
#' @param variants A data frame with `chrom` and 1-based `pos1` columns.
#'   `variant_id` is optional and is generated deterministically when absent.
#' @param activity A BigWig, bedGraph, or BAM path, or a data frame.  A table
#'   can contain genomic intervals (`chrom`, `start0`, `end0`, and `value_col`)
#'   or one precomputed value per variant.
#' @param source_type One of `"auto"`, `"table"`, `"bigwig"`,
#'   `"bedgraph"`, or `"bam"`.
#' @param chrom_lengths Named chromosome lengths, or a data frame with `chrom`
#'   and `length`.  Native BigWig/BAM lengths are used when this is `NULL`.
#' @param activity_window_bp Positive total window width; defaults to 500 bp.
#' @param activity_aggregate Either `"window_mean"` (the default) or
#'   `"covered_mean"`.
#' @param activity_uncovered_policy Must be `"zero"` for `window_mean` and
#'   `"exclude"` for `covered_mean`.  `NULL` selects the corresponding policy.
#' @param activity_transform `"log1p"` (default) or `"identity"`.  `log1p`
#'   rejects every finite negative source value before interval aggregation as
#'   well as every negative precomputed/extracted value.
#' @param activity_assay_id,activity_assembly,activity_units Provenance labels.
#' @param value_col Numeric signal column in a table.  It defaults to `score`
#'   for interval tables and `activity_raw` for precomputed tables when those
#'   columns exist.
#' @param table_mode `"auto"`, `"interval"`, or `"variant"`.
#' @param bam_mapq Minimum BAM mapping quality; defaults to 30.
#' @param bam_library_size Optional positive filtered fragment/read library
#'   size.  When omitted it is counted by streaming through the BAM.
#'
#' @return A data frame containing the input variant columns and activity
#'   value, transformation, exact window geometry, coverage, query-failure,
#'   BAM-policy (when applicable), and status fields.
#' @export
extract_v2g_activity <- function(
    variants,
    activity,
    source_type = c("auto", "table", "bigwig", "bedgraph", "bam"),
    chrom_lengths = NULL,
    activity_window_bp = 500L,
    activity_aggregate = c("window_mean", "covered_mean"),
    activity_uncovered_policy = NULL,
    activity_transform = c("log1p", "identity"),
    activity_assay_id = "activity",
    activity_assembly = NA_character_,
    activity_units = NA_character_,
    value_col = NULL,
    table_mode = c("auto", "interval", "variant"),
    bam_mapq = 30L,
    bam_library_size = NULL) {
  variants <- .as_plain_data_frame(variants, "variants")
  .require_columns(variants, c("chrom", "pos1"), "variants")
  source_type <- match.arg(source_type)
  activity_aggregate <- match.arg(activity_aggregate)
  activity_transform <- match.arg(activity_transform)
  table_mode <- match.arg(table_mode)
  activity_window_bp <- .assert_scalar_number(
    activity_window_bp, "activity_window_bp", lower = 1, integer = TRUE
  )
  bam_mapq <- .assert_scalar_number(
    bam_mapq, "bam_mapq", lower = 0, integer = TRUE
  )

  if (is.null(activity_uncovered_policy)) {
    activity_uncovered_policy <- if (activity_aggregate == "window_mean") {
      "zero"
    } else {
      "exclude"
    }
  }
  expected_policy <- if (activity_aggregate == "window_mean") "zero" else "exclude"
  if (!identical(activity_uncovered_policy, expected_policy)) {
    .acv2g_validation_error(
      "{.arg activity_uncovered_policy} must be {.val {expected_policy}} when {.arg activity_aggregate} is {.val {activity_aggregate}}."
    )
  }

  if (!"variant_id" %in% names(variants)) {
    variants$variant_id <- .stable_row_id(
      variants$chrom, variants$pos1, seq_len(nrow(variants)), prefix = "variant"
    )
  }
  if (anyNA(variants$variant_id) || any(!nzchar(as.character(variants$variant_id)))) {
    .acv2g_validation_error("{.field variant_id} cannot be missing or empty.")
  }
  if (anyDuplicated(variants$variant_id)) {
    .acv2g_validation_error(
      "{.field variant_id} must be unique for activity extraction."
    )
  }
  variants$chrom <- as.character(variants$chrom)
  pos <- suppressWarnings(as.numeric(variants$pos1))
  if (anyNA(pos) || any(!is.finite(pos)) || any(pos != floor(pos)) || any(pos < 1)) {
    .acv2g_validation_error(
      "{.field pos1} must contain finite, positive, 1-based integer positions."
    )
  }
  variants$pos1 <- pos

  source_type <- .activity_source_type(activity, source_type)
  if (source_type == "bam" && grepl("\\.cram$", activity, ignore.case = TRUE)) {
    .acv2g_validation_error(
      "CRAM input is not supported; provide an indexed coordinate-sorted BAM."
    )
  }

  native_lengths <- NULL
  if (source_type == "bigwig") {
    native_lengths <- .activity_bigwig_lengths(activity)
  } else if (source_type == "bam") {
    native_lengths <- .activity_bam_header(activity)$targets
  }
  lengths <- .activity_resolve_lengths(chrom_lengths, native_lengths)
  missing_length <- !variants$chrom %in% names(lengths) |
    is.na(lengths[variants$chrom])
  if (any(missing_length)) {
    .acv2g_validation_error(
      "Chromosome length is unavailable for: {paste(unique(variants$chrom[missing_length]), collapse = ', ')}."
    )
  }
  variant_lengths <- unname(lengths[variants$chrom])
  outside <- (variants$pos1 - 1) >= variant_lengths
  if (any(outside)) {
    bad <- unique(paste0(variants$chrom[outside], ":", variants$pos1[outside]))
    .acv2g_validation_error(
      "Variant anchor{?s} outside chromosome bounds: {paste(bad, collapse = ', ')}."
    )
  }

  anchor0 <- variants$pos1 - 1
  requested_start0 <- anchor0 - floor(activity_window_bp / 2)
  requested_end0 <- anchor0 + ceiling(activity_window_bp / 2)
  effective_start0 <- pmax(0, requested_start0)
  effective_end0 <- pmin(variant_lengths, requested_end0)
  effective_width <- effective_end0 - effective_start0
  windows <- GenomicRanges::GRanges(
    seqnames = variants$chrom,
    ranges = IRanges::IRanges(
      start = effective_start0 + 1,
      end = effective_end0
    )
  )
  names(windows) <- as.character(variants$variant_id)

  measured <- switch(
    source_type,
    table = .activity_from_table(
      variants = variants,
      windows = windows,
      activity = activity,
      value_col = value_col,
      table_mode = table_mode,
      aggregate = activity_aggregate,
      reject_negative = identical(activity_transform, "log1p")
    ),
    bigwig = .activity_from_bigwig(
      windows, activity, activity_aggregate,
      reject_negative = identical(activity_transform, "log1p")
    ),
    bedgraph = .activity_from_bedgraph(
      windows, activity, value_col, activity_aggregate,
      reject_negative = identical(activity_transform, "log1p")
    ),
    bam = .activity_from_bam(
      windows, activity, activity_aggregate, bam_mapq, bam_library_size,
      reject_negative = identical(activity_transform, "log1p")
    )
  )
  if (any(is.infinite(measured$value))) {
    .acv2g_validation_error("Extracted activity contains infinite values.")
  }

  finite_negative <- is.finite(measured$value) & measured$value < 0
  if (activity_transform == "log1p" && any(finite_negative)) {
    ids <- variants$variant_id[finite_negative]
    .acv2g_validation_error(
      "{.arg activity_transform = 'log1p'} does not allow negative activity; offending variant_id{?s}: {paste(head(ids, 20L), collapse = ', ')}."
    )
  }
  transformed <- if (activity_transform == "log1p") {
    log1p(measured$value)
  } else {
    measured$value
  }

  out <- variants
  out$activity_assay_id <- rep_len(as.character(activity_assay_id), nrow(out))
  out$activity_assembly <- rep_len(as.character(activity_assembly), nrow(out))
  out$activity_units <- if (source_type == "bam") {
    "fragment_coverage_cpm"
  } else {
    rep_len(as.character(activity_units), nrow(out))
  }
  out$activity_source_type <- source_type
  out$activity_window_bp <- activity_window_bp
  out$activity_window_requested_start0 <- requested_start0
  out$activity_window_requested_end0 <- requested_end0
  out$activity_window_effective_start0 <- effective_start0
  out$activity_window_effective_end0 <- effective_end0
  out$activity_effective_window_bp <- effective_width
  out$activity_window_clipped <- requested_start0 != effective_start0 |
    requested_end0 != effective_end0
  out$activity_aggregate <- activity_aggregate
  out$activity_uncovered_policy <- activity_uncovered_policy
  out$activity_transform <- activity_transform
  out$activity_raw <- measured$value
  out$activity_transformed <- transformed
  out$activity_covered_bp <- measured$covered_bp
  out$activity_coverage_fraction <- measured$covered_bp / effective_width
  out$activity_status <- measured$status
  out$activity_query_failed <- measured$query_failed
  out$activity_query_failure_message <- measured$query_failure_message
  out$activity_bam_mapq <- if (source_type == "bam") bam_mapq else NA_integer_
  out$activity_bam_mode <- if (source_type == "bam") {
    "fragment_coverage_cpm"
  } else {
    NA_character_
  }
  out$activity_bam_filter_policy <- if (source_type == "bam") {
    "mapped_primary_qcpass_nonduplicate_nonsupplementary_v1"
  } else {
    NA_character_
  }
  out$activity_bam_library_size <- if (source_type == "bam") {
    attr(measured, "library_size")
  } else {
    NA_real_
  }
  rownames(out) <- NULL
  out
}

.activity_source_type <- function(activity, source_type) {
  if (source_type != "auto") {
    if (source_type == "table" && !is.data.frame(activity)) {
      .acv2g_validation_error(
        "{.arg activity} must be a data frame when {.arg source_type = 'table'}."
      )
    }
    if (source_type != "table" &&
        (!is.character(activity) || length(activity) != 1L || is.na(activity))) {
      .acv2g_validation_error(
        "{.arg activity} must be one file path for {.arg source_type = {source_type}}."
      )
    }
    return(source_type)
  }
  if (is.data.frame(activity)) {
    return("table")
  }
  if (!is.character(activity) || length(activity) != 1L || is.na(activity)) {
    .acv2g_validation_error(
      "{.arg activity} must be a data frame or one supported file path."
    )
  }
  lower <- tolower(activity)
  if (grepl("\\.(bw|bigwig)$", lower)) return("bigwig")
  if (grepl("\\.(bedgraph|bdg)(\\.gz)?$", lower)) return("bedgraph")
  if (grepl("\\.bam$", lower)) return("bam")
  if (grepl("\\.cram$", lower)) return("bam")
  .acv2g_validation_error(
    "Cannot infer activity format from {.file {activity}}; set {.arg source_type}."
  )
}

.activity_resolve_lengths <- function(chrom_lengths, native_lengths = NULL) {
  normalise <- function(x, arg) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x)) {
      .require_columns(x, c("chrom", "length"), arg)
      if (anyDuplicated(x$chrom)) {
        .acv2g_validation_error("{.arg {arg}} contains duplicate chromosomes.")
      }
      ans <- suppressWarnings(as.numeric(x$length))
      names(ans) <- as.character(x$chrom)
    } else {
      ans <- suppressWarnings(as.numeric(x))
      names(ans) <- names(x)
    }
    if (is.null(names(ans)) || any(!nzchar(names(ans))) || anyNA(ans) ||
        any(!is.finite(ans)) || any(ans <= 0) || any(ans != floor(ans))) {
      .acv2g_validation_error(
        "{.arg {arg}} must provide named, finite, positive integer lengths."
      )
    }
    ans
  }
  supplied <- normalise(chrom_lengths, "chrom_lengths")
  native <- normalise(native_lengths, "native chromosome lengths")
  if (is.null(supplied) && is.null(native)) {
    .acv2g_validation_error(
      "{.arg chrom_lengths} is required when the activity source does not provide chromosome lengths."
    )
  }
  if (!is.null(supplied) && !is.null(native)) {
    shared <- intersect(names(supplied), names(native))
    conflict <- shared[supplied[shared] != native[shared]]
    if (length(conflict)) {
      .acv2g_validation_error(
        "Supplied and native chromosome lengths disagree for: {paste(conflict, collapse = ', ')}."
      )
    }
  }
  if (!is.null(supplied)) supplied else native
}

.activity_bigwig_lengths <- function(path) {
  if (!file.exists(path)) {
    .acv2g_validation_error("Activity file does not exist: {.file {path}}.")
  }
  info <- tryCatch(
    GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(path)),
    error = function(e) {
      .acv2g_validation_error(
        "Could not read BigWig metadata from {.file {path}}: {conditionMessage(e)}"
      )
    }
  )
  GenomeInfoDb::seqlengths(info)
}

.activity_from_bigwig <- function(
    windows, path, aggregate, reject_negative = FALSE) {
  source_info <- GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(path))
  available <- GenomeInfoDb::seqlevels(source_info)
  compatible <- as.character(GenomeInfoDb::seqnames(windows)) %in% available
  result <- .activity_empty_measurement(length(windows), aggregate)
  if (any(compatible)) {
    signal <- tryCatch(
      rtracklayer::import(
        rtracklayer::BigWigFile(path), which = windows[compatible]
      ),
      error = function(e) e
    )
    if (inherits(signal, "error")) {
      result <- .activity_mark_query_failure(
        result, compatible, "bigwig", path, signal
      )
      result$value[!compatible] <- NA_real_
      result$covered_bp[!compatible] <- NA_real_
      result$status[!compatible] <- "missing_activity_seqlevel"
      return(result)
    }
    scores <- as.numeric(S4Vectors::mcols(signal)$score)
    measured <- .activity_integrate_intervals(
      windows[compatible], signal, scores, aggregate, allow_overlap = FALSE,
      available_seqlevels = available,
      reject_negative = reject_negative
    )
    result[compatible, ] <- measured
  }
  result$value[!compatible] <- NA_real_
  result$covered_bp[!compatible] <- NA_real_
  result$status[!compatible] <- "missing_activity_seqlevel"
  result
}

.activity_from_bedgraph <- function(
    windows, path, value_col, aggregate, reject_negative = FALSE) {
  if (!file.exists(path)) {
    .acv2g_validation_error("Activity file does not exist: {.file {path}}.")
  }
  signal <- tryCatch(
    rtracklayer::import(path, format = "bedGraph"),
    error = function(e) e
  )
  if (inherits(signal, "error")) {
    return(.activity_mark_query_failure(
      .activity_empty_measurement(length(windows), aggregate),
      rep(TRUE, length(windows)), "bedgraph", path, signal
    ))
  }
  m <- as.data.frame(S4Vectors::mcols(signal), stringsAsFactors = FALSE)
  if (is.null(value_col)) {
    value_col <- if ("score" %in% names(m)) "score" else names(m)[1L]
  }
  if (!value_col %in% names(m)) {
    .acv2g_validation_error(
      "bedGraph value column {.field {value_col}} is unavailable."
    )
  }
  .activity_integrate_intervals(
    windows, signal, suppressWarnings(as.numeric(m[[value_col]])), aggregate,
    allow_overlap = FALSE,
    available_seqlevels = GenomeInfoDb::seqlevels(signal),
    reject_negative = reject_negative
  )
}

.activity_from_table <- function(
    variants, windows, activity, value_col, table_mode, aggregate,
    reject_negative = FALSE) {
  activity <- .as_plain_data_frame(activity, "activity")
  if (table_mode == "auto") {
    table_mode <- if (all(c("chrom", "start0", "end0") %in% names(activity))) {
      "interval"
    } else {
      "variant"
    }
  }
  if (table_mode == "interval") {
    .require_columns(activity, c("chrom", "start0", "end0"), "activity")
    if (is.null(value_col)) {
      value_col <- if ("score" %in% names(activity)) "score" else NULL
    }
    if (is.null(value_col) || !value_col %in% names(activity)) {
      .acv2g_validation_error(
        "Interval activity tables require an explicit numeric {.arg value_col}."
      )
    }
    start0 <- suppressWarnings(as.numeric(activity$start0))
    end0 <- suppressWarnings(as.numeric(activity$end0))
    raw_values <- activity[[value_col]]
    values <- suppressWarnings(as.numeric(raw_values))
    if (any(is.na(values) & !is.na(raw_values))) {
      .acv2g_validation_error(
        "Activity value column {.field {value_col}} must be numeric."
      )
    }
    bad <- is.na(start0) | is.na(end0) | start0 != floor(start0) |
      end0 != floor(end0) | start0 < 0 | end0 <= start0
    if (any(bad)) {
      .acv2g_validation_error(
        "Activity intervals must be valid 0-based half-open intervals."
      )
    }
    signal <- GenomicRanges::GRanges(
      as.character(activity$chrom),
      IRanges::IRanges(start0 + 1, end0)
    )
    return(.activity_integrate_intervals(
      windows, signal, values, aggregate, allow_overlap = FALSE,
      available_seqlevels = unique(as.character(activity$chrom)),
      reject_negative = reject_negative
    ))
  }

  if (is.null(value_col)) {
    value_col <- if ("activity_raw" %in% names(activity)) {
      "activity_raw"
    } else if ("score" %in% names(activity)) {
      "score"
    } else {
      NULL
    }
  }
  if (is.null(value_col) || !value_col %in% names(activity)) {
    .acv2g_validation_error(
      "Precomputed activity tables require an explicit numeric {.arg value_col}."
    )
  }
  if (all(c("chrom", "pos1") %in% names(activity))) {
    key <- "chrom + pos1"
    activity_key <- paste(activity$chrom, activity$pos1, sep = "\u001f")
    if (anyDuplicated(activity_key)) {
      .acv2g_validation_error(
        "Precomputed activity has duplicate {.field chrom + pos1} anchors."
      )
    }
    lookup <- match(
      paste(variants$chrom, variants$pos1, sep = "\u001f"), activity_key
    )
  } else if ("variant_id" %in% names(activity)) {
    key <- "variant_id"
    lookup <- match(variants$variant_id, activity$variant_id)
    if (anyDuplicated(activity$variant_id)) {
      .acv2g_validation_error(
        "Precomputed activity has duplicate {.field variant_id} values."
      )
    }
  } else {
    .acv2g_validation_error(
      paste0(
        "Precomputed activity requires chrom + pos1 coordinates or a ",
        "variant_id matching the extraction input."
      )
    )
  }
  table_values <- suppressWarnings(as.numeric(activity[[value_col]]))
  if (any(is.na(table_values) & !is.na(activity[[value_col]]))) {
    .acv2g_validation_error(
      "Precomputed activity column {.field {value_col}} must be numeric."
    )
  }
  values <- rep(NA_real_, nrow(variants))
  values[!is.na(lookup)] <- table_values[lookup[!is.na(lookup)]]
  covered <- rep(NA_real_, nrow(variants))
  if ("activity_covered_bp" %in% names(activity)) {
    covered[!is.na(lookup)] <- suppressWarnings(as.numeric(
      activity$activity_covered_bp[lookup[!is.na(lookup)]]
    ))
    effective_width <- GenomicRanges::width(windows)
    invalid_coverage <- is.finite(covered) &
      (covered < 0 | covered > effective_width | covered != floor(covered))
    if (any(invalid_coverage)) {
      .acv2g_validation_error(
        "Precomputed {.field activity_covered_bp} must be an integer between zero and the effective window width."
      )
    }
  }
  status <- ifelse(
    is.na(lookup), paste0("missing_precomputed_", gsub(" \\+ ", "_", key)),
    ifelse(is.finite(values), "ok_precomputed", "missing_precomputed_value")
  )
  if ("chrom" %in% names(activity)) {
    available_seqlevels <- unique(as.character(activity$chrom))
    absent_seqlevel <- !variants$chrom %in% available_seqlevels
    values[absent_seqlevel] <- NA_real_
    covered[absent_seqlevel] <- NA_real_
    status[absent_seqlevel] <- "missing_activity_seqlevel"
  }
  result <- .activity_empty_measurement(nrow(variants), aggregate)
  result$value <- values
  result$covered_bp <- covered
  result$status <- status
  result
}

.activity_integrate_intervals <- function(
    windows, signal, scores, aggregate, allow_overlap = FALSE,
    available_seqlevels = GenomeInfoDb::seqlevels(signal),
    reject_negative = FALSE) {
  if (length(signal) != length(scores)) {
    .acv2g_validation_error(
      "Activity signal ranges and values have different lengths."
    )
  }
  if (any(is.infinite(scores))) {
    .acv2g_validation_error("Activity signal contains infinite values.")
  }
  if (isTRUE(reject_negative) && any(is.finite(scores) & scores < 0)) {
    .acv2g_validation_error(
      paste0(
        "activity_transform = 'log1p' does not allow negative source ",
        "activity values before window aggregation."
      )
    )
  }
  if (!allow_overlap && length(signal) > 1L) {
    ord <- order(
      as.character(GenomeInfoDb::seqnames(signal)),
      GenomicRanges::start(signal), GenomicRanges::end(signal)
    )
    chr <- as.character(GenomeInfoDb::seqnames(signal))[ord]
    s <- GenomicRanges::start(signal)[ord]
    e <- GenomicRanges::end(signal)[ord]
    overlap <- chr[-1L] == chr[-length(chr)] & s[-1L] <= e[-length(e)]
    if (any(overlap)) {
      .acv2g_validation_error(
        "Activity interval signal contains overlapping records; their combination is undefined."
      )
    }
  }
  result <- .activity_empty_measurement(length(windows), aggregate)
  window_seqlevels <- as.character(GenomeInfoDb::seqnames(windows))
  absent_seqlevel <- !window_seqlevels %in% available_seqlevels
  result$value[absent_seqlevel] <- NA_real_
  result$covered_bp[absent_seqlevel] <- NA_real_
  result$status[absent_seqlevel] <- "missing_activity_seqlevel"
  common_seqlevels <- intersect(
    GenomeInfoDb::seqlevels(windows), available_seqlevels
  )
  if (!length(common_seqlevels)) {
    return(result)
  }
  hits <- GenomicRanges::findOverlaps(windows, signal, ignore.strand = TRUE)
  if (!length(hits)) return(result)
  query_index <- S4Vectors::queryHits(hits)
  subject_index <- S4Vectors::subjectHits(hits)
  overlap_start <- pmax(
    GenomicRanges::start(windows)[query_index],
    GenomicRanges::start(signal)[subject_index]
  )
  overlap_end <- pmin(
    GenomicRanges::end(windows)[query_index],
    GenomicRanges::end(signal)[subject_index]
  )
  overlap_width <- overlap_end - overlap_start + 1

  # All production adapters require non-overlapping signal records. In that
  # case each overlap width contributes to exactly one covered segment, so the
  # integration can be reduced by window in one vectorized pass. Keep the
  # general reduction path below for explicit internal callers that allow
  # overlapping records.
  if (!isTRUE(allow_overlap)) {
    hit_scores <- scores[subject_index]
    missing_score <- is.na(hit_scores)
    weighted <- hit_scores * overlap_width
    weighted[missing_score] <- 0
    aggregate_matrix <- rowsum(
      cbind(
        covered_bp = as.numeric(overlap_width),
        weighted_signal = as.numeric(weighted),
        missing_value_count = as.numeric(missing_score)
      ),
      group = query_index,
      reorder = FALSE
    )
    window_index <- as.integer(rownames(aggregate_matrix))
    missing_window <- aggregate_matrix[, "missing_value_count"] > 0
    if (any(missing_window)) {
      result$status[window_index[missing_window]] <- "missing_activity_value"
    }
    valid <- !missing_window
    if (any(valid)) {
      valid_window <- window_index[valid]
      covered <- aggregate_matrix[valid, "covered_bp"]
      denominator <- if (aggregate == "window_mean") {
        GenomicRanges::width(windows)[valid_window]
      } else {
        covered
      }
      value <- aggregate_matrix[valid, "weighted_signal"] / denominator
      result$value[valid_window] <- value
      result$covered_bp[valid_window] <- covered
      result$status[valid_window] <- ifelse(
        is.finite(value), "ok", "missing_activity_value"
      )
    }
    return(result)
  }

  by_query <- split(S4Vectors::subjectHits(hits), S4Vectors::queryHits(hits))
  for (i_chr in names(by_query)) {
    i <- as.integer(i_chr)
    j <- by_query[[i_chr]]
    overlap_start <- pmax(GenomicRanges::start(windows)[i],
                          GenomicRanges::start(signal)[j])
    overlap_end <- pmin(GenomicRanges::end(windows)[i],
                        GenomicRanges::end(signal)[j])
    overlap_width <- overlap_end - overlap_start + 1
    if (anyNA(scores[j])) {
      result$status[i] <- "missing_activity_value"
      next
    }
    covered_ranges <- IRanges::reduce(IRanges::IRanges(overlap_start, overlap_end))
    covered <- sum(IRanges::width(covered_ranges))
    numerator <- sum(scores[j] * overlap_width)
    denominator <- if (aggregate == "window_mean") {
      GenomicRanges::width(windows)[i]
    } else {
      covered
    }
    result$value[i] <- if (denominator > 0) numerator / denominator else NA_real_
    result$covered_bp[i] <- covered
    result$status[i] <- if (is.finite(result$value[i])) "ok" else "missing_activity_value"
  }
  result
}

.activity_empty_measurement <- function(n, aggregate) {
  data.frame(
    value = if (aggregate == "window_mean") rep(0, n) else rep(NA_real_, n),
    covered_bp = rep(0, n),
    status = if (aggregate == "window_mean") {
      rep("ok_uncovered_zero", n)
    } else {
      rep("no_covered_bases", n)
    },
    query_failed = rep(FALSE, n),
    query_failure_message = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
}

.activity_mark_query_failure <- function(
    result, affected, source_type, path, condition) {
  affected <- as.logical(affected)
  if (length(affected) != nrow(result)) {
    .acv2g_validation_error(
      "Internal activity query-failure mask has the wrong length."
    )
  }
  message <- conditionMessage(condition)
  result$value[affected] <- NA_real_
  result$covered_bp[affected] <- NA_real_
  result$status[affected] <- "missing_activity_query"
  result$query_failed[affected] <- TRUE
  result$query_failure_message[affected] <- message
  .acv2g_warn(
    paste0(
      source_type, " activity query failed for ", sum(affected),
      " variant window(s); affected activity values are retained as missing."
    ),
    data = list(
      source_type = source_type,
      path = path,
      affected_window_count = sum(affected),
      message = message
    ),
    class = "acv2g_activity_query_warning"
  )
  result
}

.activity_bam_index_candidates <- function(path) {
  unique(c(
    paste0(path, ".bai"), sub("\\.bam$", ".bai", path, ignore.case = TRUE),
    paste0(path, ".csi"), sub("\\.bam$", ".csi", path, ignore.case = TRUE)
  ))
}

.activity_bam_index_paths <- function(path) {
  candidates <- .activity_bam_index_candidates(path)
  candidates[file.exists(candidates) & !dir.exists(candidates)]
}

.activity_bam_header <- function(path) {
  if (!file.exists(path)) {
    .acv2g_validation_error("BAM file does not exist: {.file {path}}.")
  }
  index_paths <- .activity_bam_index_paths(path)
  if (!length(index_paths)) {
    .acv2g_validation_error(
      "BAM must have a sibling BAI or CSI index: {.file {path}}."
    )
  }
  header <- tryCatch(
    Rsamtools::scanBamHeader(path)[[1L]],
    error = function(e) {
      .acv2g_validation_error(
        "Could not read BAM header from {.file {path}}: {conditionMessage(e)}"
      )
    }
  )
  text <- header$text
  hd <- text[["@HD"]]
  sort_order <- NA_character_
  if (length(hd)) {
    if (is.list(hd)) {
      first <- hd[[1L]]
      if (is.list(first) && !is.null(first$SO)) sort_order <- first$SO
      if (!is.list(first) && "SO" %in% names(hd)) sort_order <- hd[["SO"]]
    }
    fields <- unlist(hd, recursive = TRUE, use.names = FALSE)
    so_field <- fields[grepl("^SO:", fields)]
    if (length(so_field)) sort_order <- sub("^SO:", "", so_field[1L])
  }
  if (!identical(as.character(sort_order), "coordinate")) {
    .acv2g_validation_error(
      "BAM header must explicitly declare coordinate sorting (SO:coordinate)."
    )
  }
  list(
    targets = header$targets,
    sort_order = sort_order,
    index_paths = normalizePath(index_paths, mustWork = TRUE)
  )
}

.activity_bam_flag <- function() {
  Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isNotPassingQualityControls = FALSE,
    isDuplicate = FALSE,
    isSupplementaryAlignment = FALSE
  )
}

.activity_bam_library_size <- function(path, mapq) {
  bf <- Rsamtools::BamFile(path, yieldSize = 1000000L)
  Rsamtools::open.BamFile(bf)
  on.exit(Rsamtools::close.BamFile(bf), add = TRUE)
  total <- 0
  repeat {
    chunk <- Rsamtools::scanBam(
      bf,
      param = Rsamtools::ScanBamParam(
        flag = .activity_bam_flag(), what = c("flag", "mapq", "isize")
      )
    )[[1L]]
    if (!length(chunk$flag)) break
    keep_mapq <- !is.na(chunk$mapq) & chunk$mapq >= mapq
    paired <- bitwAnd(chunk$flag, 0x1L) != 0L
    proper <- bitwAnd(chunk$flag, 0x2L) != 0L
    mate_mapped <- bitwAnd(chunk$flag, 0x8L) == 0L
    total <- total + sum(keep_mapq & !paired) +
      sum(keep_mapq & paired & proper & mate_mapped & chunk$isize > 0L)
  }
  if (!is.finite(total) || total <= 0) {
    .acv2g_validation_error(
      "No BAM fragments/reads pass the configured filters; CPM is undefined."
    )
  }
  as.numeric(total)
}

.activity_from_bam <- function(
    windows, path, aggregate, mapq, library_size = NULL,
    reject_negative = FALSE) {
  header <- .activity_bam_header(path)
  available_seqlevels <- names(header$targets)
  compatible <- as.character(GenomeInfoDb::seqnames(windows)) %in%
    available_seqlevels
  if (is.null(library_size)) {
    library_size <- tryCatch(
      .activity_bam_library_size(path, mapq),
      error = function(e) e
    )
    if (inherits(library_size, "error")) {
      result <- .activity_empty_measurement(length(windows), aggregate)
      result$value[!compatible] <- NA_real_
      result$covered_bp[!compatible] <- NA_real_
      result$status[!compatible] <- "missing_activity_seqlevel"
      result <- .activity_mark_query_failure(
        result, compatible, "bam_library_scan", path, library_size
      )
      attr(result, "library_size") <- NA_real_
      return(result)
    }
  } else {
    library_size <- .assert_scalar_number(
      library_size, "bam_library_size", lower = 0, inclusive = FALSE
    )
  }
  if (!any(compatible)) {
    result <- .activity_empty_measurement(length(windows), aggregate)
    result$value[] <- NA_real_
    result$covered_bp[] <- NA_real_
    result$status[] <- "missing_activity_seqlevel"
    attr(result, "library_size") <- library_size
    return(result)
  }
  query <- GenomicRanges::reduce(windows[compatible], ignore.strand = TRUE)
  records <- tryCatch(
    Rsamtools::scanBam(
      path,
      param = Rsamtools::ScanBamParam(
        which = query,
        flag = .activity_bam_flag(),
        what = c("qname", "flag", "mapq", "rname", "pos", "cigar", "isize")
      )
    ),
    error = function(e) e
  )
  if (inherits(records, "error")) {
    result <- .activity_empty_measurement(length(windows), aggregate)
    result$value[!compatible] <- NA_real_
    result$covered_bp[!compatible] <- NA_real_
    result$status[!compatible] <- "missing_activity_seqlevel"
    result <- .activity_mark_query_failure(
      result, compatible, "bam", path, records
    )
    attr(result, "library_size") <- library_size
    return(result)
  }
  records <- lapply(records, as.data.frame, stringsAsFactors = FALSE)
  records <- records[vapply(records, nrow, integer(1L)) > 0L]
  if (!length(records)) {
    result <- .activity_empty_measurement(length(windows), aggregate)
    result$value[!compatible] <- NA_real_
    result$covered_bp[!compatible] <- NA_real_
    result$status[!compatible] <- "missing_activity_seqlevel"
    attr(result, "library_size") <- library_size
    return(result)
  }
  reads <- do.call(rbind, records)
  reads <- reads[!is.na(reads$mapq) & reads$mapq >= mapq, , drop = FALSE]
  if (!nrow(reads)) {
    result <- .activity_empty_measurement(length(windows), aggregate)
    result$value[!compatible] <- NA_real_
    result$covered_bp[!compatible] <- NA_real_
    result$status[!compatible] <- "missing_activity_seqlevel"
    attr(result, "library_size") <- library_size
    return(result)
  }
  reads <- unique(reads)
  paired <- bitwAnd(reads$flag, 0x1L) != 0L
  proper <- bitwAnd(reads$flag, 0x2L) != 0L
  mate_mapped <- bitwAnd(reads$flag, 0x8L) == 0L
  reference_width <- GenomicAlignments::cigarWidthAlongReferenceSpace(
    reads$cigar
  )
  read_end <- reads$pos + reference_width - 1L

  single <- reads[!paired, , drop = FALSE]
  single_width <- reference_width[!paired]
  single_ranges <- if (nrow(single)) {
    GenomicRanges::GRanges(
      as.character(single$rname),
      IRanges::IRanges(single$pos, single$pos + single_width - 1L)
    )
  } else {
    GenomicRanges::GRanges()
  }

  pair_keep <- paired & proper & mate_mapped & reads$isize != 0L
  pair_reads <- reads[pair_keep, , drop = FALSE]
  pair_end <- read_end[pair_keep]
  fragment_start <- ifelse(
    pair_reads$isize > 0L,
    pair_reads$pos,
    pair_end - abs(pair_reads$isize) + 1L
  )
  fragment_end <- ifelse(
    pair_reads$isize > 0L,
    pair_reads$pos + pair_reads$isize - 1L,
    pair_end
  )
  fragment_key <- paste(
    pair_reads$qname, pair_reads$rname, fragment_start, fragment_end, sep = "\u001f"
  )
  keep_fragment <- !duplicated(fragment_key)
  pair_ranges <- if (any(keep_fragment)) {
    GenomicRanges::GRanges(
      as.character(pair_reads$rname[keep_fragment]),
      IRanges::IRanges(
        fragment_start[keep_fragment], fragment_end[keep_fragment]
      )
    )
  } else {
    GenomicRanges::GRanges()
  }
  fragments <- c(single_ranges, pair_ranges)
  if (!length(fragments)) {
    result <- .activity_empty_measurement(length(windows), aggregate)
    result$value[!compatible] <- NA_real_
    result$covered_bp[!compatible] <- NA_real_
    result$status[!compatible] <- "missing_activity_seqlevel"
  } else {
    result <- .activity_integrate_intervals(
      windows,
      fragments,
      rep(1e6 / library_size, length(fragments)),
      aggregate,
      allow_overlap = TRUE,
      available_seqlevels = available_seqlevels,
      reject_negative = reject_negative
    )
  }
  attr(result, "library_size") <- library_size
  result
}
