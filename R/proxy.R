.acv2g_proxy_row_labels <- function(x, rows, row_id_col = NULL, n = 5L) {
  rows <- utils::head(rows, n)
  if (!is.null(row_id_col)) {
    .require_columns(x, row_id_col, "pairs")
    ids <- as.character(x[[row_id_col]][rows])
    ids[is.na(ids) | !nzchar(ids)] <- paste0("row_", rows[is.na(ids) | !nzchar(ids)])
    return(ids)
  }

  candidate <- intersect(
    c("pair_id", "variant_gene_pair_id", "variant_key", "variant_id"),
    names(x)
  )
  if (length(candidate)) {
    ids <- as.character(x[[candidate[[1L]]]][rows])
    ids[is.na(ids) | !nzchar(ids)] <- paste0("row_", rows[is.na(ids) | !nzchar(ids)])
    return(ids)
  }
  paste0("row_", rows)
}

.acv2g_validate_numeric_evidence <- function(
    x,
    column,
    data,
    transform = NULL,
    row_id_col = NULL) {
  if (!is.numeric(x)) {
    .acv2g_validation_error(
      sprintf("Column `%s` must be numeric.", column)
    )
  }

  invalid <- which(!is.na(x) & !is.finite(x))
  if (length(invalid)) {
    labels <- .acv2g_proxy_row_labels(data, invalid, row_id_col)
    .acv2g_validation_error(
      sprintf(
        paste0(
          "Column `%s` contains %d non-finite, non-missing value(s); ",
          "representative row identifier(s): %s."
        ),
        column,
        length(invalid),
        paste(labels, collapse = ", ")
      ),
      data = list(column = column, rows = invalid, row_ids = labels)
    )
  }

  if (identical(transform, "log1p")) {
    negative <- which(is.finite(x) & x < 0)
    if (length(negative)) {
      labels <- .acv2g_proxy_row_labels(data, negative, row_id_col)
      .acv2g_validation_error(
        sprintf(
          paste0(
            "Column `%s` contains %d finite negative value(s), which are ",
            "not allowed with the `log1p` transform; representative row ",
            "identifier(s): %s."
          ),
          column,
          length(negative),
          paste(labels, collapse = ", ")
        ),
        data = list(column = column, rows = negative, row_ids = labels)
      )
    }
  }
  invisible(x)
}

#' Compute the activity-contact proxy
#'
#' Computes a pair-level activity-contact proxy without changing the pair
#' universe. Missing activity or contact remains missing, while observed zero
#' evidence remains numeric zero. The default is
#' `log1p(activity) * log1p(contact)`; `identity` is the only alternative
#' transform in the first ACV2G release.
#'
#' @param pairs A data frame containing one row per variant-gene pair.
#' @param activity_col Name of the numeric activity column.
#' @param contact_col Name of the numeric contact column.
#' @param activity_transform Either `"log1p"` or `"identity"`.
#' @param contact_transform Either `"log1p"` or `"identity"`.
#' @param activity_output_col Name for the transformed activity column. Use
#'   `NULL` to omit this component column.
#' @param contact_output_col Name for the transformed contact column. Use
#'   `NULL` to omit this component column.
#' @param output_col Name for the computed proxy column.
#' @param row_id_col Optional column used to identify invalid rows in errors.
#'
#' @return `pairs`, in its original order, with transformed components (unless
#'   omitted) and the proxy column appended or replaced.
#' @export
compute_v2g_proxy <- function(
    pairs,
    activity_col = "activity_value",
    contact_col = "contact_value",
    activity_transform = c("log1p", "identity"),
    contact_transform = c("log1p", "identity"),
    activity_output_col = "activity_transformed",
    contact_output_col = "contact_transformed",
    output_col = "proxy_score",
    row_id_col = NULL) {
  pairs <- .as_plain_data_frame(pairs, "pairs")
  activity_transform <- match.arg(activity_transform)
  contact_transform <- match.arg(contact_transform)
  .require_columns(pairs, c(activity_col, contact_col), "pairs")

  output_names <- c(activity_output_col, contact_output_col, output_col)
  output_names <- output_names[!vapply(output_names, is.null, logical(1L))]
  if (anyNA(output_names) || any(!nzchar(output_names)) || anyDuplicated(output_names)) {
    .acv2g_validation_error(
      "Requested proxy output column names must be non-empty and distinct."
    )
  }

  activity <- pairs[[activity_col]]
  contact <- pairs[[contact_col]]
  .acv2g_validate_numeric_evidence(
    activity, activity_col, pairs, activity_transform, row_id_col
  )
  .acv2g_validate_numeric_evidence(
    contact, contact_col, pairs, contact_transform, row_id_col
  )

  transformed_activity <- if (identical(activity_transform, "log1p")) {
    log1p(activity)
  } else {
    activity
  }
  transformed_contact <- if (identical(contact_transform, "log1p")) {
    log1p(contact)
  } else {
    contact
  }
  proxy <- transformed_activity * transformed_contact

  invalid_proxy <- which(!is.na(proxy) & !is.finite(proxy))
  if (length(invalid_proxy)) {
    labels <- .acv2g_proxy_row_labels(pairs, invalid_proxy, row_id_col)
    .acv2g_validation_error(
      sprintf(
        paste0(
          "Proxy multiplication produced %d non-finite value(s); ",
          "representative row identifier(s): %s."
        ),
        length(invalid_proxy),
        paste(labels, collapse = ", ")
      ),
      data = list(rows = invalid_proxy, row_ids = labels)
    )
  }

  if (!is.null(activity_output_col)) {
    pairs[[activity_output_col]] <- transformed_activity
  }
  if (!is.null(contact_output_col)) {
    pairs[[contact_output_col]] <- transformed_contact
  }
  pairs[[output_col]] <- proxy
  pairs
}
