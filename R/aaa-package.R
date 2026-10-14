#' ACV2G: Activity–Contact Variant-to-Gene Prioritization
#'
#' ACV2G combines regulatory activity and chromatin contact to prioritize
#' variant-gene pairs. It is inspired by the activity-by-contact idea but does
#' not implement the canonical enhancer-centric ABC model.
#'
#' @import data.table
#' @importFrom stats aggregate ave setNames
#' @importFrom utils head
#' @keywords internal
"_PACKAGE"

.acv2g_schema_version <- "1.0.0"

#' Return the ACV2G result schema version
#'
#' @return A length-one character vector.
#' @export
acv2g_schema_version <- function() {
  .acv2g_schema_version
}

utils::globalVariables(".data")
