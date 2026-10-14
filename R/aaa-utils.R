.acv2g_abort <- function(message, ..., class = "acv2g_error") {
  cli::cli_abort(
    message,
    ...,
    class = c(class, "acv2g_error"),
    .envir = parent.frame()
  )
}

.acv2g_validation_error <- function(message, ..., data = NULL) {
  cli::cli_abort(
    message,
    ...,
    class = c("acv2g_validation_error", "acv2g_error"),
    data = data,
    .envir = parent.frame()
  )
}

.acv2g_warn <- function(message, ..., class = "acv2g_warning") {
  cli::cli_warn(
    message,
    ...,
    class = c(class, "acv2g_warning"),
    .envir = parent.frame()
  )
}

.require_columns <- function(x, columns, object_name = deparse(substitute(x))) {
  missing_columns <- setdiff(columns, names(x))
  if (length(missing_columns)) {
    .acv2g_validation_error(
      "{.arg {object_name}} is missing required column{?s}: {.field {missing_columns}}."
    )
  }
  invisible(x)
}

.as_plain_data_frame <- function(x, object_name = deparse(substitute(x))) {
  if (inherits(x, "DataFrame")) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) {
    .acv2g_validation_error("{.arg {object_name}} must be a data frame.")
  }
  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

.assert_scalar_choice <- function(x, choices, arg = deparse(substitute(x))) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    .acv2g_validation_error(
      "{.arg {arg}} must be exactly one of {.val {choices}}; received {.val {x}}."
    )
  }
  x
}

.assert_scalar_number <- function(
    x,
    arg = deparse(substitute(x)),
    lower = -Inf,
    upper = Inf,
    integer = FALSE,
    inclusive = TRUE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (ok && integer) {
    ok <- x == floor(x)
  }
  if (ok && inclusive) {
    ok <- x >= lower && x <= upper
  } else if (ok) {
    ok <- x > lower && x < upper
  }
  if (!ok) {
    .acv2g_validation_error("{.arg {arg}} is outside its accepted domain.")
  }
  x
}

.empty_table <- function(columns) {
  out <- setNames(vector("list", length(columns)), columns)
  out[] <- lapply(out, function(.) character())
  as.data.frame(out, stringsAsFactors = FALSE)
}

.utc_now <- function() {
  format(Sys.time(), tz = "UTC", usetz = TRUE)
}

.stable_row_id <- function(..., prefix = "row") {
  pieces <- list(...)
  n <- max(vapply(pieces, length, integer(1L)), 0L)
  pieces <- lapply(pieces, rep_len, length.out = n)
  key <- do.call(paste, c(pieces, sep = "\u001f"))
  vector_digest <- tryCatch(
    digest::getVDigest("xxhash64"),
    error = function(e) NULL
  )
  hashes <- if (is.function(vector_digest)) {
    vector_digest(key, serialize = FALSE)
  } else {
    vapply(
      key, digest::digest, character(1L),
      algo = "xxhash64", serialize = FALSE
    )
  }
  paste0(prefix, "_", hashes)
}
