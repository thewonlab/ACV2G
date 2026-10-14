#' List ACV2G scoring methods
#'
#' @return A data frame describing the two implemented normalization methods.
#' @export
v2g_methods <- function() {
  data.frame(
    method = c("locus_pair_z", "gene_wise_z"),
    default = c(TRUE, FALSE),
    estimand = c(
      "variant-gene pair proxy standardized within a fixed user-supplied locus",
      "variant-gene pair proxy standardized within an explicit per-gene reference universe"
    ),
    reference_requirement = c(
      "fixed locus membership",
      "explicit input_variants or background_variants mode"
    ),
    output_column = c("locus_pair_z", "gene_wise_z"),
    stringsAsFactors = FALSE
  )
}

#' Summarize an ACV2G result
#'
#' @param x A `v2g_result` or bundle path.
#'
#' @return A named list of compact run counts and status.
#' @export
summarize_v2g_result <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    x <- read_v2g_result(x, lazy = TRUE)
  }
  if (!inherits(x, "v2g_result")) {
    .acv2g_validation_error("{.arg x} must be a result object or bundle path.")
  }
  pairs <- v2g_table(x, "pairs")
  warnings <- v2g_table(x, "warnings")
  list(
    schema_version = x$schema_version,
    run_status = x$run_status$status,
    profiles = if ("profile_id" %in% names(pairs))
      length(unique(pairs$profile_id)) else 0L,
    loci = if ("locus_id" %in% names(pairs))
      length(unique(stats::na.omit(pairs$locus_id))) else 0L,
    variants = if ("variant_anchor_id" %in% names(pairs))
      length(unique(pairs$variant_anchor_id)) else 0L,
    genes = if ("gene_id_normalized" %in% names(pairs))
      length(unique(pairs$gene_id_normalized)) else 0L,
    pairs = nrow(pairs),
    finite_scores = {
      score_col <- intersect(c("locus_pair_z", "gene_wise_z"), names(pairs))
      if (length(score_col)) sum(is.finite(pairs[[score_col[[1L]]]])) else 0L
    },
    called_pairs = if ("link_call" %in% names(pairs))
      sum(pairs$link_call %in% TRUE, na.rm = TRUE) else NA_integer_,
    warnings = nrow(warnings)
  )
}
