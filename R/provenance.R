#' Hash a local ACV2G input resource
#'
#' @param path Existing local file.
#' @param role Optional input-role label.
#' @param assembly Optional declared assembly.
#'
#' @return A one-row data frame.
#' @export
hash_v2g_resource <- function(path, role = NA_character_, assembly = NA_character_) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path) || dir.exists(path)) {
    .acv2g_validation_error("Resource to hash must be one existing file.")
  }
  normalized <- normalizePath(path, mustWork = TRUE)
  info <- file.info(normalized)
  data.frame(
    role = as.character(role),
    path = normalized,
    assembly = as.character(assembly),
    size_bytes = as.numeric(info$size),
    modified_at_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    sha256 = digest::digest(file = normalized, algo = "sha256", serialize = FALSE),
    stringsAsFactors = FALSE
  )
}

.installed_version_or_na <- function(package) {
  tryCatch(as.character(utils::packageVersion(package)),
           error = function(e) NA_character_)
}

.acv2g_software_versions <- function() {
  packages <- c(
    "ACV2G", "R", "cli", "data.table", "digest", "GenomeInfoDb",
    "GenomicAlignments", "GenomicRanges", "ggplot2", "HiCExperiment",
    "IRanges", "jsonlite", "optparse", "Rsamtools", "rtracklayer",
    "S4Vectors", "strawr", "svglite", "stats", "utils",
    "VariantAnnotation", "yaml"
  )
  versions <- vapply(packages, function(package) {
    if (package == "R") as.character(getRversion()) else
      .installed_version_or_na(package)
  }, character(1L))
  as.list(versions)
}

.collect_config_resources <- function(config) {
  rows <- list()
  add <- function(path, role, assembly) {
    if (is.character(path) && length(path) == 1L && file.exists(path) &&
        !dir.exists(path)) {
      rows[[length(rows) + 1L]] <<- hash_v2g_resource(path, role, assembly)
    }
  }
  for (role in names(config$inputs %||% list())) {
    spec <- config$inputs[[role]]
    if (is.list(spec)) {
      add(spec$path, role, spec$assembly %||% NA_character_)
      add(spec$chain, paste0(role, ":chain"), spec$assembly %||% NA_character_)
      add(spec$reverse_chain, paste0(role, ":reverse_chain"),
          config$analysis_build %||% NA_character_)
      if (identical(role, "genes") && is.null(spec$path) &&
          !is.null(spec$managed_id)) {
        resolved <- resolve_v2g_reference(
          spec$managed_id,
          registry = config$reference_registry %||% NULL,
          cache_dir = config$reference_cache_dir %||% NULL
        )
        add(resolved$path, paste0(role, ":managed:", spec$managed_id),
            resolved$assembly %||% spec$assembly %||% NA_character_)
      }
    }
  }
  for (profile_id in names(config$profiles %||% list())) {
    profile <- config$profiles[[profile_id]]
    add(profile$activity$path, paste0("profile:", profile_id, ":activity"),
        profile$activity$assembly)
    activity_path <- profile$activity$path %||% NULL
    activity_type <- profile$activity$type %||% "auto"
    is_bam <- identical(activity_type, "bam") ||
      (identical(activity_type, "auto") && is.character(activity_path) &&
       length(activity_path) == 1L &&
       grepl("\\.bam$", activity_path, ignore.case = TRUE))
    if (isTRUE(is_bam) && is.character(activity_path) &&
        length(activity_path) == 1L) {
      for (index_path in .activity_bam_index_paths(activity_path)) {
        add(
          index_path,
          paste0(
            "profile:", profile_id, ":activity_index:", basename(index_path)
          ),
          profile$activity$assembly
        )
      }
    }
    add(profile$activity$chain,
        paste0("profile:", profile_id, ":activity:chain"),
        profile$activity$assembly)
    add(profile$activity$reverse_chain,
        paste0("profile:", profile_id, ":activity:reverse_chain"),
        config$analysis_build %||% NA_character_)
    activity_chrom_sizes <- profile$activity$chrom_sizes %||% NULL
    if (is.list(activity_chrom_sizes) &&
        !is.data.frame(activity_chrom_sizes)) {
      add(
        activity_chrom_sizes$path,
        paste0("profile:", profile_id, ":activity:chrom_sizes"),
        activity_chrom_sizes$assembly %||% profile$activity$assembly
      )
    }
    add(profile$contact$path, paste0("profile:", profile_id, ":contact"),
        profile$contact$assembly)
  }
  add(config$reference_registry, "reference_registry", NA_character_)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

#' Build an ACV2G provenance manifest
#'
#' @param config Fully resolved configuration.
#' @param profiles Resolved profile table.
#' @param extra Additional named manifest fields.
#'
#' @return A JSON-serializable list.
#' @export
build_v2g_manifest <- function(config, profiles = NULL, extra = list()) {
  resources <- .collect_config_resources(config)
  profile_ids <- if (is.null(profiles)) names(config$profiles %||% list()) else
    as.character(profiles$profile_id)
  base <- list(
    schema_version = .acv2g_schema_version,
    package = "ACV2G",
    package_version = .installed_version_or_na("ACV2G"),
    created_at_utc = .utc_now(),
    analysis_build = config$analysis_build %||% NA_character_,
    profile_ids = profile_ids,
    resolved_method = config$method %||% list(),
    resolved_activity = config$activity %||% list(),
    resolved_contact = config$contact %||% list(),
    resolved_liftover = config$liftover %||% list(),
    input_resources = if (nrow(resources)) split(resources, seq_len(nrow(resources))) else list(),
    software_versions = .acv2g_software_versions(),
    platform = list(
      os = Sys.info()[["sysname"]],
      release = Sys.info()[["release"]],
      machine = Sys.info()[["machine"]]
    )
  )
  utils::modifyList(base, extra, keep.null = TRUE)
}

.warning_record <- function(condition, index) {
  data <- condition$data
  data_text <- if (is.null(data)) NA_character_ else {
    tryCatch(jsonlite::toJSON(data, auto_unbox = TRUE, na = "null"),
             error = function(e) NA_character_)
  }
  data.frame(
    schema_version = .acv2g_schema_version,
    warning_id = paste0("warning_", sprintf("%06d", index)),
    class = paste(class(condition), collapse = ";"),
    message = conditionMessage(condition),
    profile_id = condition$profile_id %||% NA_character_,
    locus_id = condition$locus_id %||% NA_character_,
    condition_data_json = data_text,
    timestamp_utc = .utc_now(),
    stringsAsFactors = FALSE
  )
}

.capture_warnings <- function(expr) {
  records <- list()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      records[[length(records) + 1L]] <<- .warning_record(w, length(records) + 1L)
    }
  )
  warnings <- if (length(records)) do.call(rbind, records) else
    .default_result_table("warnings")
  list(value = value, warnings = warnings)
}
