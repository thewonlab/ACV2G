.default_reference_registry <- function() {
  path <- system.file("extdata", "reference-registry.yaml", package = "ACV2G")
  if (!nzchar(path)) {
    .acv2g_validation_error(
      "The bundled reference registry is unavailable; install ACV2G or supply {.arg registry} explicitly."
    )
  }
  path
}

.read_reference_registry <- function(registry = NULL) {
  if (is.null(registry)) registry <- .default_reference_registry()
  if (is.character(registry) && length(registry) == 1L) {
    if (!file.exists(registry)) {
      .acv2g_validation_error("Reference registry does not exist: {.path {registry}}.")
    }
    registry_path <- normalizePath(registry, mustWork = TRUE)
    registry <- yaml::read_yaml(registry_path)
    attr(registry, "registry_path") <- registry_path
  }
  if (!is.list(registry)) {
    .acv2g_validation_error("Reference registry must be a YAML path or list.")
  }
  entries <- if ("references" %in% names(registry)) registry$references else registry
  if (is.null(entries)) entries <- list()
  if (!is.list(entries)) {
    .acv2g_validation_error("Reference registry entries must be a named mapping.")
  }
  if (length(entries) && (is.null(names(entries)) || any(!nzchar(names(entries))) ||
                          anyDuplicated(names(entries)))) {
    .acv2g_validation_error("Reference registry IDs must be unique and non-empty.")
  }
  entries
}

.validate_reference_id <- function(managed_id) {
  if (!is.character(managed_id) || length(managed_id) != 1L ||
      is.na(managed_id) || !nzchar(managed_id) || managed_id == "latest") {
    .acv2g_validation_error(
      "{.arg managed_id} must be one exact pinned identifier and cannot be {.val latest}."
    )
  }
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", managed_id)) {
    .acv2g_validation_error("Managed reference ID contains unsafe characters.")
  }
  managed_id
}

.validate_reference_entry <- function(entry, managed_id) {
  required <- c("url", "sha256", "filename", "assembly", "format")
  missing <- setdiff(required, names(entry))
  if (length(missing)) {
    .acv2g_validation_error(
      "Reference {.val {managed_id}} is missing registry field{?s}: {.field {missing}}."
    )
  }
  if (!grepl("^[0-9a-fA-F]{64}$", entry$sha256)) {
    .acv2g_validation_error(
      "Reference {.val {managed_id}} must declare an exact SHA-256 checksum."
    )
  }
  if (basename(entry$filename) != entry$filename ||
      entry$filename %in% c("", ".", "..")) {
    .acv2g_validation_error("Reference registry filename must be a safe basename.")
  }
  entry
}

#' List pinned managed annotation references
#'
#' @param registry Optional registry YAML path or list. The bundled registry is
#'   used by default.
#'
#' @return A data frame of exact managed IDs and metadata.
#' @export
list_v2g_references <- function(registry = NULL) {
  entries <- .read_reference_registry(registry)
  if (!length(entries)) {
    return(data.frame(
      managed_id = character(), assembly = character(), format = character(),
      filename = character(), sha256 = character(), url = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(names(entries), function(id) {
    entry <- .validate_reference_entry(entries[[id]], id)
    data.frame(
      managed_id = id,
      assembly = as.character(entry$assembly),
      format = as.character(entry$format),
      filename = as.character(entry$filename),
      sha256 = tolower(as.character(entry$sha256)),
      url = as.character(entry$url),
      stringsAsFactors = FALSE
    )
  }))
}

.reference_cache_dir <- function(cache_dir = NULL) {
  if (is.null(cache_dir)) cache_dir <- tools::R_user_dir("ACV2G", "cache")
  if (!is.character(cache_dir) || length(cache_dir) != 1L ||
      is.na(cache_dir) || !nzchar(cache_dir)) {
    .acv2g_validation_error("{.arg cache_dir} must be one non-empty path.")
  }
  path.expand(cache_dir)
}

#' Fetch one exact managed annotation into a portable cache
#'
#' This is the only ACV2G operation that downloads annotation data. Prediction
#' resolves cached resources and never invokes this function implicitly.
#'
#' @param managed_id Exact pinned registry ID.
#' @param registry Optional registry YAML path or list.
#' @param cache_dir Cache root; defaults to the platform-standard R user cache.
#' @param force Replace an existing verified cache entry.
#'
#' @return The cached annotation path, invisibly.
#' @export
fetch_v2g_reference <- function(managed_id, registry = NULL, cache_dir = NULL,
                                force = FALSE) {
  managed_id <- .validate_reference_id(managed_id)
  entries <- .read_reference_registry(registry)
  if (!managed_id %in% names(entries)) {
    .acv2g_validation_error(
      "Unknown managed reference {.val {managed_id}}. Available IDs: {.val {names(entries)}}."
    )
  }
  entry <- .validate_reference_entry(entries[[managed_id]], managed_id)
  cache_dir <- .reference_cache_dir(cache_dir)
  entry_dir <- file.path(cache_dir, managed_id)
  target <- file.path(entry_dir, entry$filename)
  metadata_path <- file.path(entry_dir, "metadata.json")
  if (file.exists(target) && !isTRUE(force)) {
    observed <- digest::digest(file = target, algo = "sha256", serialize = FALSE)
    if (identical(tolower(observed), tolower(entry$sha256))) {
      return(invisible(normalizePath(target, mustWork = TRUE)))
    }
    .acv2g_validation_error(
      "Cached reference checksum mismatch for {.val {managed_id}}; use {.arg force = TRUE} only after reviewing the registry."
    )
  }
  dir.create(entry_dir, recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", entry$filename, "."),
                        tmpdir = entry_dir)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::download.file(entry$url, temporary, mode = "wb", quiet = FALSE)
  observed <- digest::digest(file = temporary, algo = "sha256", serialize = FALSE)
  if (!identical(tolower(observed), tolower(entry$sha256))) {
    .acv2g_validation_error(
      "Downloaded reference checksum mismatch for {.val {managed_id}}."
    )
  }
  if (file.exists(target)) unlink(target)
  if (!file.rename(temporary, target)) {
    .acv2g_abort("Could not publish cached reference {.path {target}}.")
  }
  metadata <- c(
    list(
      managed_id = managed_id,
      fetched_at_utc = .utc_now(),
      schema_version = .acv2g_schema_version
    ),
    entry
  )
  jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE,
                       pretty = TRUE, na = "null")
  invisible(normalizePath(target, mustWork = TRUE))
}

#' Resolve a cached managed reference without network access
#'
#' @param managed_id Exact pinned registry ID.
#' @param registry Optional registry YAML path or list.
#' @param cache_dir Cache root.
#'
#' @return A list containing the verified path and registry metadata.
#' @export
resolve_v2g_reference <- function(managed_id, registry = NULL, cache_dir = NULL) {
  managed_id <- .validate_reference_id(managed_id)
  entries <- .read_reference_registry(registry)
  if (!managed_id %in% names(entries)) {
    .acv2g_validation_error("Unknown managed reference {.val {managed_id}}.")
  }
  entry <- .validate_reference_entry(entries[[managed_id]], managed_id)
  target <- file.path(.reference_cache_dir(cache_dir), managed_id, entry$filename)
  if (!file.exists(target)) {
    .acv2g_validation_error(
      "Managed reference {.val {managed_id}} is not cached. Run the separate reference fetch command first."
    )
  }
  observed <- digest::digest(file = target, algo = "sha256", serialize = FALSE)
  if (!identical(tolower(observed), tolower(entry$sha256))) {
    .acv2g_validation_error("Cached reference checksum mismatch for {.val {managed_id}}.")
  }
  c(list(managed_id = managed_id, path = normalizePath(target, mustWork = TRUE)), entry)
}
