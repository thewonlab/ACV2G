# Coordinate helpers -------------------------------------------------------

.validate_assembly <- function(assembly, arg = "assembly") {
  if (!is.character(assembly) || length(assembly) != 1L || is.na(assembly) ||
      !nzchar(trimws(assembly))) {
    .acv2g_validation_error(
      paste0(arg, " must be one non-empty, explicitly declared assembly name.")
    )
  }
  trimws(assembly)
}

.validate_chrom <- function(chrom, arg = "chrom") {
  chrom <- as.character(chrom)
  bad <- is.na(chrom) | !nzchar(trimws(chrom))
  if (any(bad)) {
    .acv2g_validation_error(
      paste0(arg, " contains ", sum(bad), " missing or empty chromosome name(s).")
    )
  }
  trimws(chrom)
}

.as_integer_coordinate <- function(x, arg, lower = 0) {
  numeric_x <- suppressWarnings(as.numeric(x))
  bad <- is.na(numeric_x) | !is.finite(numeric_x) |
    numeric_x != floor(numeric_x) | numeric_x < lower
  if (any(bad)) {
    .acv2g_validation_error(
      paste0(
        arg, " must contain finite integer coordinates >= ", lower, "; ",
        sum(bad), " value(s) failed validation."
      )
    )
  }
  if (any(numeric_x > .Machine$integer.max)) {
    # R's integer type cannot hold some valid genome coordinates. Retain an
    # exactly represented double rather than overflowing to NA.
    return(numeric_x)
  }
  as.integer(numeric_x)
}

#' Convert 1-based genomic positions to 0-based anchors
#'
#' This is an internal, vectorized coordinate primitive. ACV2G represents a
#' point position as the one-base half-open interval `[anchor0, anchor0 + 1)`.
#'
#' @param pos1 Finite positive integer positions.
#' @return Numeric 0-based anchors.
#' @keywords internal
.pos1_to_anchor0 <- function(pos1) {
  pos1 <- .as_integer_coordinate(pos1, "pos1", lower = 1)
  as.numeric(pos1) - 1
}

#' Convert 0-based anchors to 1-based genomic positions
#'
#' @param anchor0 Finite non-negative integer anchors.
#' @return Numeric 1-based positions.
#' @keywords internal
.anchor0_to_pos1 <- function(anchor0) {
  anchor0 <- .as_integer_coordinate(anchor0, "anchor0", lower = 0)
  as.numeric(anchor0) + 1
}

#' Validate 0-based half-open intervals
#'
#' @param start0,end0 Interval boundaries.
#' @param object_name Label used in validation messages.
#' @return A data frame with validated `start0` and `end0` columns.
#' @keywords internal
.validate_half_open <- function(start0, end0, object_name = "intervals") {
  start0 <- .as_integer_coordinate(start0, paste0(object_name, "$start0"), 0)
  end0 <- .as_integer_coordinate(end0, paste0(object_name, "$end0"), 0)
  if (length(start0) != length(end0)) {
    .acv2g_validation_error(
      paste0(object_name, " start and end vectors must have equal lengths.")
    )
  }
  bad <- end0 <= start0
  if (any(bad)) {
    .acv2g_validation_error(
      paste0(
        object_name, " contains ", sum(bad),
        " empty or reversed 0-based half-open interval(s)."
      )
    )
  }
  data.frame(start0 = start0, end0 = end0, stringsAsFactors = FALSE)
}

#' Map positions to fixed-width 0-based contact bins
#'
#' @param pos1 Finite positive integer genomic positions.
#' @param resolution_bp A finite positive integer resolution.
#' @return A data frame with `anchor0`, `bin_start0`, and `bin_end0`.
#' @keywords internal
.fixed_contact_bins <- function(pos1, resolution_bp) {
  resolution_bp <- .assert_scalar_number(
    resolution_bp,
    arg = "resolution_bp",
    lower = 1,
    integer = TRUE
  )
  anchor0 <- .pos1_to_anchor0(pos1)
  bin_start0 <- floor(anchor0 / resolution_bp) * resolution_bp
  data.frame(
    anchor0 = anchor0,
    bin_start0 = bin_start0,
    bin_end0 = bin_start0 + resolution_bp,
    stringsAsFactors = FALSE
  )
}

#' Construct deterministic variant-centered activity windows
#'
#' The requested 0-based half-open window is
#' `[anchor0 - floor(width / 2), anchor0 + ceiling(width / 2))`. When
#' chromosome lengths are supplied, valid near-edge windows are clipped by
#' intersection and are never shifted.
#'
#' @param chrom Chromosome names.
#' @param pos1 1-based variant positions.
#' @param activity_window_bp Positive integer total window width.
#' @param chromosome_lengths Optional named numeric vector of chromosome
#'   lengths. If supplied, every chromosome must be present and variant
#'   anchors must fall inside their chromosome.
#' @return A data frame of requested and effective half-open coordinates.
#' @keywords internal
.activity_window_geometry <- function(
    chrom,
    pos1,
    activity_window_bp = 500L,
    chromosome_lengths = NULL) {
  activity_window_bp <- .assert_scalar_number(
    activity_window_bp,
    arg = "activity_window_bp",
    lower = 1,
    integer = TRUE
  )
  chrom <- .validate_chrom(chrom)
  if (length(chrom) != length(pos1)) {
    .acv2g_validation_error(
      "{.arg chrom} and {.arg pos1} must have equal lengths."
    )
  }
  anchor0 <- .pos1_to_anchor0(pos1)
  requested_start0 <- anchor0 - floor(activity_window_bp / 2)
  requested_end0 <- anchor0 + ceiling(activity_window_bp / 2)

  if (is.null(chromosome_lengths)) {
    effective_start0 <- requested_start0
    effective_end0 <- requested_end0
    clipped <- rep(FALSE, length(anchor0))
  } else {
    if (is.null(names(chromosome_lengths)) ||
        anyNA(names(chromosome_lengths)) ||
        any(!nzchar(names(chromosome_lengths)))) {
      .acv2g_validation_error(
        "{.arg chromosome_lengths} must be a named vector."
      )
    }
    lengths_numeric <- suppressWarnings(as.numeric(chromosome_lengths))
    bad_lengths <- is.na(lengths_numeric) | !is.finite(lengths_numeric) |
      lengths_numeric != floor(lengths_numeric) | lengths_numeric <= 0
    if (any(bad_lengths)) {
      .acv2g_validation_error(
        "{.arg chromosome_lengths} must contain positive finite integers."
      )
    }
    names(lengths_numeric) <- names(chromosome_lengths)
    unresolved <- setdiff(unique(chrom), names(lengths_numeric))
    if (length(unresolved)) {
      .acv2g_validation_error(
        paste0(
          "Chromosome length(s) unavailable for: ",
          paste(unresolved, collapse = ", "), "."
        )
      )
    }
    chromosome_length <- unname(lengths_numeric[chrom])
    out_of_bounds <- anchor0 >= chromosome_length
    if (any(out_of_bounds)) {
      examples <- unique(paste0(chrom[out_of_bounds], ":", pos1[out_of_bounds]))
      .acv2g_validation_error(
        paste0(
          sum(out_of_bounds),
          " variant anchor(s) fall outside declared chromosome bounds: ",
          paste(head(examples, 5L), collapse = ", "), "."
        )
      )
    }
    effective_start0 <- pmax(0, requested_start0)
    effective_end0 <- pmin(chromosome_length, requested_end0)
    clipped <- effective_start0 != requested_start0 |
      effective_end0 != requested_end0
  }

  data.frame(
    chrom = chrom,
    pos1 = as.numeric(pos1),
    anchor0 = anchor0,
    activity_window_requested_start0 = requested_start0,
    activity_window_requested_end0 = requested_end0,
    activity_window_effective_start0 = effective_start0,
    activity_window_effective_end0 = effective_end0,
    activity_effective_window_bp = effective_end0 - effective_start0,
    activity_window_clipped = clipped,
    stringsAsFactors = FALSE
  )
}
