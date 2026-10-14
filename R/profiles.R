.acv2g_profile_required_path <- function(profile, role, profile_id) {
  profile_names <- names(profile) %||% character()
  locations <- which(profile_names == role)
  if (length(locations) != 1L) {
    .acv2g_validation_error(
      paste0(
        "Profile ", profile_id, " must contain exactly one nested ", role,
        " specification."
      )
    )
  }
  specification <- profile[[locations]]
  if (!is.list(specification) || is.data.frame(specification)) {
    .acv2g_validation_error(
      paste0("Profile ", profile_id, " ", role, " must be a named list.")
    )
  }
  specification_names <- names(specification) %||% character()
  path_locations <- which(specification_names == "path")
  if (length(path_locations) != 1L) {
    .acv2g_validation_error(
      paste0(
        "Profile ", profile_id, " ", role,
        " must contain exactly one path."
      )
    )
  }
  path <- specification[[path_locations]]
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(trimws(path))) {
    .acv2g_validation_error(
      paste0("Profile ", profile_id, " ", role, " path must be one non-empty string.")
    )
  }
  type_locations <- which(specification_names == "type")
  if (length(type_locations) > 1L) {
    .acv2g_validation_error(
      paste0("Profile ", profile_id, " ", role, " contains duplicate type fields.")
    )
  }
  type <- if (length(type_locations)) specification[[type_locations]] else
    NA_character_
  if (!length(type)) type <- NA_character_
  if (length(type) != 1L || (!is.na(type) &&
      (!is.character(type) || !nzchar(trimws(type))))) {
    .acv2g_validation_error(
      paste0("Profile ", profile_id, " ", role, " type must be one string when supplied.")
    )
  }
  list(path = path, type = as.character(type))
}

.acv2g_profile_flat_scalar <- function(profile, field, profile_id) {
  locations <- which((names(profile) %||% character()) == field)
  if (length(locations) != 1L) {
    .acv2g_validation_error(
      paste0("Profile ", profile_id, " must declare exactly one ", field, ".")
    )
  }
  value <- profile[[locations]]
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    .acv2g_validation_error(
      paste0("Profile ", profile_id, " ", field, " must be one non-empty string.")
    )
  }
  value
}

.acv2g_profile_optional_flat_scalar <- function(profile, field, profile_id) {
  locations <- which((names(profile) %||% character()) == field)
  if (!length(locations)) return(NA_character_)
  if (length(locations) != 1L) {
    .acv2g_validation_error(
      paste0("Profile ", profile_id, " contains duplicate ", field, " fields.")
    )
  }
  value <- profile[[locations]]
  if (!length(value)) return(NA_character_)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    .acv2g_validation_error(
      paste0("Profile ", profile_id, " ", field, " must be one string when supplied.")
    )
  }
  value
}

.acv2g_profile_list_row <- function(profile, profile_id) {
  if (!is.list(profile) || is.data.frame(profile)) {
    .acv2g_validation_error("Each profile specification must be a named list.")
  }
  profile_names <- names(profile) %||% character()
  has_nested <- any(profile_names %in% c("activity", "contact"))
  has_flat <- any(profile_names %in% c("activity_source", "contact_source"))
  if (has_nested && has_flat) {
    .acv2g_validation_error(
      paste0(
        "Profile ", profile_id,
        " mixes nested and flat evidence specifications; use one representation."
      )
    )
  }

  if (has_nested) {
    activity <- .acv2g_profile_required_path(profile, "activity", profile_id)
    contact <- .acv2g_profile_required_path(profile, "contact", profile_id)
  } else {
    activity <- list(
      path = .acv2g_profile_flat_scalar(profile, "activity_source", profile_id),
      type = .acv2g_profile_optional_flat_scalar(
        profile, "activity_type", profile_id
      )
    )
    contact <- list(
      path = .acv2g_profile_flat_scalar(profile, "contact_source", profile_id),
      type = .acv2g_profile_optional_flat_scalar(
        profile, "contact_type", profile_id
      )
    )
  }
  data.frame(
    profile_id = profile_id,
    activity_source = activity$path,
    contact_source = contact$path,
    activity_type = activity$type,
    contact_type = contact$type,
    stringsAsFactors = FALSE
  )
}

#' Validate explicit ACV2G scoring profiles
#'
#' One profile is exactly one activity source plus one contact source. ACV2G
#' never constructs replicate averages or activity-by-contact cross-products.
#' Canonical list input uses named profile entries containing
#' `activity = list(path = ...)` and `contact = list(path = ...)`.
#'
#' @param profiles A data frame with `profile_id`, `activity_source`, and
#'   `contact_source`, or a named list of canonical nested profile
#'   specifications.
#'
#' @return A standardized data frame.
#' @export
validate_v2g_profiles <- function(profiles) {
  if (is.list(profiles) && !is.data.frame(profiles)) {
    profile_ids <- names(profiles)
    if (is.null(profile_ids) || anyNA(profile_ids) ||
        any(!nzchar(trimws(profile_ids))) || anyDuplicated(profile_ids)) {
      .acv2g_validation_error("A profile list must have unique non-empty names.")
    }
    rows <- lapply(seq_along(profiles), function(index) {
      .acv2g_profile_list_row(profiles[[index]], profile_ids[[index]])
    })
    profiles <- if (length(rows)) do.call(rbind, rows) else data.frame()
  }
  profiles <- .as_plain_data_frame(profiles, "profiles")
  .require_columns(profiles, c("profile_id", "activity_source", "contact_source"),
                   "profiles")
  profiles$profile_id <- as.character(profiles$profile_id)
  if (anyNA(profiles$profile_id) || any(!nzchar(trimws(profiles$profile_id))) ||
      anyDuplicated(profiles$profile_id)) {
    .acv2g_validation_error("Profile IDs must be unique, non-missing, and non-empty.")
  }
  for (field in c("activity_source", "contact_source")) {
    value <- profiles[[field]]
    if (is.list(value)) {
      .acv2g_validation_error(
        "Every profile must declare exactly one scalar {.field {field}}."
      )
    }
    value <- as.character(value)
    invalid <- is.na(value) | !nzchar(trimws(value))
    if (any(invalid)) {
      .acv2g_validation_error(
        "Every profile must declare exactly one non-missing {.field {field}}."
      )
    }
    profiles[[field]] <- value
  }
  profiles$schema_version <- .acv2g_schema_version
  profiles[, c("schema_version", setdiff(names(profiles), "schema_version")),
           drop = FALSE]
}
