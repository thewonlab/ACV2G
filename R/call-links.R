#' Apply an optional descriptive Z-score threshold
#'
#' ACV2G Z-scores are dataset-relative prioritization scores, not calibrated
#' null-test statistics. Consequently, no link calls are made by default. When
#' explicitly selected, `z_threshold` calls rows whose finite Z-score is greater
#' than or equal to the threshold; a missing Z-score always gives a missing
#' call rather than `FALSE`.
#'
#' @param pairs A scored pair data frame.
#' @param call_rule Either `"none"` (the overall default) or
#'   `"z_threshold"`.
#' @param z_col Name of the Z-score column to inspect. This can be
#'   `locus_pair_z`, `gene_wise_z`, or a caller-specified score column.
#' @param z_threshold Threshold used only for `call_rule = "z_threshold"`.
#'   When omitted under that rule, the conditional default is 1.5.
#' @param output_col Name of the logical output column.
#' @param call_scope Provenance label describing the scope of each call.
#'
#' @return `pairs`, in its original order, with the logical call and complete
#'   call-rule provenance appended.
#' @export
call_v2g_links <- function(
    pairs,
    call_rule = c("none", "z_threshold"),
    z_col = "locus_pair_z",
    z_threshold = NULL,
    output_col = "link_call",
    call_scope = "pair") {
  pairs <- .as_plain_data_frame(pairs, "pairs")
  call_rule <- match.arg(call_rule)
  if (!is.character(z_col) || length(z_col) != 1L || is.na(z_col) || !nzchar(z_col)) {
    .acv2g_validation_error("`z_col` must be one non-empty column name.")
  }
  if (!is.character(output_col) || length(output_col) != 1L ||
      is.na(output_col) || !nzchar(output_col)) {
    .acv2g_validation_error("`output_col` must be one non-empty column name.")
  }
  if (!is.character(call_scope) || length(call_scope) != 1L ||
      is.na(call_scope) || !nzchar(call_scope)) {
    .acv2g_validation_error("`call_scope` must be one non-empty label.")
  }

  .require_columns(pairs, z_col, "pairs")
  .acv2g_validate_numeric_evidence(
    pairs[[z_col]], z_col, pairs, transform = NULL
  )

  if (identical(call_rule, "none")) {
    if (!is.null(z_threshold)) {
      .acv2g_validation_error(
        "`z_threshold` must be `NULL` when `call_rule = \"none\"`."
      )
    }
    threshold <- NA_real_
    threshold_source <- "not_applicable"
    comparator <- "none"
    link_call <- rep(NA, nrow(pairs))
  } else {
    if (is.null(z_threshold)) {
      threshold <- 1.5
      threshold_source <- "conditional_default"
    } else {
      threshold <- .assert_scalar_number(z_threshold, arg = "z_threshold")
      threshold_source <- "user"
    }
    comparator <- ">="
    z <- pairs[[z_col]]
    link_call <- rep(NA, length(z))
    finite <- is.finite(z)
    link_call[finite] <- z[finite] >= threshold
  }

  pairs[[output_col]] <- as.logical(link_call)
  pairs$call_rule <- rep(call_rule, nrow(pairs))
  pairs$call_threshold <- rep(threshold, nrow(pairs))
  pairs$call_threshold_source <- rep(threshold_source, nrow(pairs))
  pairs$call_comparator <- rep(comparator, nrow(pairs))
  pairs$call_direction <- rep("higher", nrow(pairs))
  pairs$call_scope <- rep(call_scope, nrow(pairs))
  pairs
}
