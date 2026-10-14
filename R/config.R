.config_path_fields <- c(
  "path", "chain", "reverse_chain", "background_variants",
  "reference_gene_variant_map", "variant_locus_map", "chrom_sizes",
  "reference_registry", "reference_cache_dir"
)

.is_absolute_config_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\|//)", path)
}

.resolve_config_paths <- function(x, base_dir, parent_name = NULL) {
  if (!is.list(x) || is.data.frame(x)) return(x)
  for (name in names(x)) {
    value <- x[[name]]
    if (is.list(value) && !is.data.frame(value)) {
      x[[name]] <- .resolve_config_paths(value, base_dir, name)
    } else if ((name %in% .config_path_fields ||
                (name == "protein_coding_gene_ids" &&
                 is.character(value) && length(value) == 1L &&
                 grepl("[/\\\\]|\\.(txt|tsv|csv)(\\.gz)?$", value,
                       ignore.case = TRUE))) && is.character(value) &&
               length(value) == 1L && !is.na(value) && nzchar(value) &&
               !grepl("^[A-Za-z][A-Za-z0-9+.-]*://", value) &&
               !grepl("^<", value)) {
      x[[name]] <- if (.is_absolute_config_path(value)) {
        path.expand(value)
      } else {
        file.path(base_dir, value)
      }
    }
  }
  x
}

.nearest_config_key <- function(key, allowed) {
  if (!length(allowed)) return(NULL)
  distances <- as.numeric(utils::adist(tolower(key), tolower(allowed)))
  best <- which.min(distances)
  threshold <- max(1L, min(3L, floor(nchar(key) / 3L)))
  if (!length(best) || distances[[best]] > threshold) return(NULL)
  allowed[[best]]
}

.assert_known_config_keys <- function(x, allowed, field) {
  if (is.null(x)) return(invisible(x))
  if (!is.list(x) || is.data.frame(x)) {
    .acv2g_validation_error(paste0(field, " must be a mapping/object."))
  }
  keys <- names(x)
  if (length(x) && (is.null(keys) || anyNA(keys) || any(!nzchar(keys)))) {
    .acv2g_validation_error(paste0(field, " must use non-empty named fields."))
  }
  if (anyDuplicated(keys)) {
    duplicates <- unique(keys[duplicated(keys)])
    .acv2g_validation_error(paste0(
      field, " contains duplicate field(s): ",
      paste(duplicates, collapse = ", "), "."
    ))
  }
  unknown <- setdiff(keys, allowed)
  if (length(unknown)) {
    details <- vapply(unknown, function(key) {
      nearest <- .nearest_config_key(key, allowed)
      if (is.null(nearest)) key else paste0(key, " (did you mean ", nearest, "?)")
    }, character(1L))
    .acv2g_validation_error(paste0(
      "Unknown ", field, " field(s): ", paste(details, collapse = ", "), "."
    ))
  }
  invisible(x)
}

.config_input_allowed_keys <- function(name) {
  common <- c("path", "assembly", "chain", "reverse_chain")
  switch(
    name,
    loci = c(common, "format", "column_map"),
    variants = c(common, "format", "column_map", "bed_anchor_column"),
    background_variants = c(
      common, "format", "column_map", "bed_anchor_column"
    ),
    genes = c(
      common, "managed_id", "format", "column_map", "gene_biotype_field",
      "transcript_biotype_field", "protein_coding_values",
      "protein_coding_gene_ids", "id_version_policy"
    ),
    chrom_sizes = c("path", "assembly"),
    variant_locus_map = c("path", "assembly"),
    reference_gene_variant_map = c("path", "assembly"),
    character()
  )
}

.validate_config_key_contract <- function(config) {
  .assert_known_config_keys(
    config,
    c(
      "schema_version", "analysis_build", "inputs", "profiles", "method",
      "activity", "contact", "liftover", "reference_registry",
      "reference_cache_dir"
    ),
    "top-level configuration"
  )
  .assert_known_config_keys(
    config$method,
    c(
      "name", "gene_flank_bp", "reference_gene_radius_bp",
      "proxy_transform", "activity_transform", "contact_transform",
      "max_missing_fraction", "call_rule", "z_threshold", "reference_mode",
      "reference_set_id"
    ),
    "method"
  )
  .assert_known_config_keys(
    config$activity,
    c("activity_window_bp", "aggregation", "bam_mapq_min"),
    "activity"
  )
  .assert_known_config_keys(
    config$contact,
    c(
      "contact_value", "contact_resolution_bp", "diagonal_policy",
      "distance_policy"
    ),
    "contact"
  )
  .assert_known_config_keys(
    config$liftover,
    c("structural_failure_policy", "reciprocal_liftover_policy"),
    "liftover"
  )
  if (!is.null(config$inputs)) {
    .assert_known_config_keys(
      config$inputs,
      c(
        "loci", "variants", "genes", "background_variants", "chrom_sizes",
        "variant_locus_map", "reference_gene_variant_map"
      ),
      "inputs"
    )
    for (name in intersect(names(config$inputs), c(
      "loci", "variants", "genes", "background_variants", "chrom_sizes",
      "variant_locus_map", "reference_gene_variant_map"
    ))) {
      spec <- config$inputs[[name]]
      if (is.character(spec) && length(spec) == 1L) next
      .assert_known_config_keys(
        spec, .config_input_allowed_keys(name), paste0("inputs.", name)
      )
    }
  }
  if (!is.null(config$profiles)) {
    .assert_known_config_keys(config$profiles, names(config$profiles), "profiles")
    for (profile_id in names(config$profiles)) {
      profile <- config$profiles[[profile_id]]
      .assert_known_config_keys(
        profile, c("activity", "contact", "profile_id"),
        paste0("profiles.", profile_id)
      )
      if (!is.null(profile$profile_id) &&
          (!is.character(profile$profile_id) ||
           length(profile$profile_id) != 1L || is.na(profile$profile_id) ||
           !identical(profile$profile_id, profile_id))) {
        .acv2g_validation_error(paste0(
          "profiles.", profile_id,
          ".profile_id must equal its enclosing profile name."
        ))
      }
      if (!is.null(profile$activity) && !is.character(profile$activity)) {
        .assert_known_config_keys(
          profile$activity,
          c(
            "path", "assembly", "type", "chain", "reverse_chain",
            "chrom_sizes", "activity_window_bp", "aggregation",
            "uncovered_policy", "value_col", "table_mode", "bam_mapq",
            "bam_library_size", "assay_id", "units"
          ),
          paste0("profiles.", profile_id, ".activity")
        )
        activity_sizes <- profile$activity$chrom_sizes
        if (!is.null(activity_sizes) && !is.character(activity_sizes) &&
            !is.data.frame(activity_sizes) &&
            !(is.numeric(activity_sizes) && !is.null(names(activity_sizes)))) {
          .assert_known_config_keys(
            activity_sizes, c("path", "assembly"),
            paste0("profiles.", profile_id, ".activity.chrom_sizes")
          )
        }
      }
      if (!is.null(profile$contact) && !is.character(profile$contact)) {
        .assert_known_config_keys(
          profile$contact,
          c(
            "path", "assembly", "type", "resolution_bp", "contact_value",
            "layer", "score_col", "score_type", "layers",
            "table_representation", "diagonal_policy",
            "weights_source_provided", "assay_id", "units", "normalization"
          ),
          paste0("profiles.", profile_id, ".contact")
        )
      }
    }
  }
  invisible(config)
}

.resolve_method_transforms <- function(config) {
  method <- config$method %||% list()
  has_proxy <- "proxy_transform" %in% names(method) &&
    !is.null(method$proxy_transform)
  has_activity <- "activity_transform" %in% names(method) &&
    !is.null(method$activity_transform)
  has_contact <- "contact_transform" %in% names(method) &&
    !is.null(method$contact_transform)
  if (has_proxy && (has_activity || has_contact)) {
    .acv2g_validation_error(paste0(
      "method.proxy_transform is a backward-compatible shorthand only when ",
      "neither method.activity_transform nor method.contact_transform is set."
    ))
  }
  if (has_proxy) {
    transform <- .assert_scalar_choice(
      method$proxy_transform, c("log1p", "identity"),
      "method.proxy_transform"
    )
    method$activity_transform <- transform
    method$contact_transform <- transform
  }
  method$activity_transform <- method$activity_transform %||% "log1p"
  method$contact_transform <- method$contact_transform %||% "log1p"
  method$proxy_transform <- NULL
  config$method <- method
  config
}

#' Read an ACV2G YAML configuration
#'
#' Relative file paths are resolved relative to the configuration file, never
#' to the process working directory.
#'
#' @param path YAML configuration path.
#'
#' @return A configuration list with source-path metadata.
#' @export
read_v2g_config <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  config <- yaml::read_yaml(path)
  if (!is.list(config)) {
    .acv2g_validation_error("The YAML root must be a mapping/object.")
  }
  config <- .resolve_config_paths(config, dirname(path))
  attr(config, "config_path") <- path
  config
}

.input_spec <- function(config, name, required = TRUE) {
  spec <- config$inputs[[name]]
  if (is.null(spec)) {
    if (required) {
      .acv2g_validation_error("Configuration is missing {.field inputs.{name}}.")
    }
    return(NULL)
  }
  if (is.character(spec) && length(spec) == 1L) spec <- list(path = spec)
  if (!is.list(spec)) {
    .acv2g_validation_error("{.field inputs.{name}} must be a mapping.")
  }
  spec
}

.validate_local_resource <- function(
    spec,
    field,
    allow_managed = FALSE,
    check_files = TRUE) {
  if (allow_managed && !is.null(spec$managed_id)) {
    if (!is.character(spec$managed_id) || length(spec$managed_id) != 1L ||
        !nzchar(spec$managed_id) || identical(spec$managed_id, "latest")) {
      .acv2g_validation_error(
        "Managed annotation IDs must be exact pinned identifiers, never {.val latest}."
      )
    }
    if (!is.null(spec$path)) {
      .acv2g_validation_error(
        "Choose exactly one annotation source: {.field path} or {.field managed_id}."
      )
    }
    return(spec)
  }
  if (is.null(spec$path) || !is.character(spec$path) ||
      length(spec$path) != 1L || !nzchar(spec$path)) {
    .acv2g_validation_error("{.field {field}.path} must name one local resource.")
  }
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", spec$path)) {
    .acv2g_validation_error(
      "Prediction never downloads inputs; fetch and pin {.field {field}} separately."
    )
  }
  if (isTRUE(check_files) && !file.exists(spec$path)) {
    .acv2g_validation_error("Input resource does not exist: {.path {spec$path}}.")
  }
  spec
}

.validate_liftover_chain <- function(
    chain,
    field,
    direction,
    check_files = TRUE,
    resource_label = "Liftover chain") {
  if (inherits(chain, "Chain")) return(chain)
  if (!is.character(chain) || length(chain) != 1L || is.na(chain) ||
      !nzchar(chain)) {
    .acv2g_validation_error(
      paste0(field, " is required and must map ", direction, ".")
    )
  }
  if (isTRUE(check_files) && !file.exists(chain)) {
    .acv2g_validation_error(
      paste0(resource_label, " does not exist: ", chain, ".")
    )
  }
  chain
}

.validate_activity_query_chain <- function(
    chain,
    field,
    direction,
    check_files = TRUE) {
  .validate_liftover_chain(
    chain = chain,
    field = field,
    direction = direction,
    check_files = check_files,
    resource_label = "Activity query chain"
  )
}

.validate_structural_input_spec <- function(
    spec,
    name,
    analysis_build,
    reciprocal_policy,
    check_files = TRUE) {
  field <- paste0("inputs.", name)
  spec$assembly <- .assert_build(spec$assembly, paste0(field, ".assembly"))
  if (identical(spec$assembly, analysis_build)) return(spec)

  forward_direction <- paste0(
    spec$assembly, " (", field, " assembly) -> ",
    analysis_build, " (analysis_build)"
  )
  spec$chain <- .validate_liftover_chain(
    chain = spec$chain,
    field = paste0(field, ".chain"),
    direction = forward_direction,
    check_files = check_files,
    resource_label = paste0(field, " liftover chain")
  )
  if (identical(reciprocal_policy, "require")) {
    reverse_direction <- paste0(
      analysis_build, " (analysis_build) -> ",
      spec$assembly, " (", field, " assembly)"
    )
    spec$reverse_chain <- .validate_liftover_chain(
      chain = spec$reverse_chain,
      field = paste0(field, ".reverse_chain"),
      direction = reverse_direction,
      check_files = check_files,
      resource_label = paste0(field, " reverse liftover chain")
    )
  }
  spec
}

.validate_analysis_chrom_sizes <- function(
    spec,
    analysis_build,
    check_files = TRUE) {
  if (!is.list(spec) || is.data.frame(spec)) {
    .acv2g_validation_error(
      paste0(
        "{.field inputs.chrom_sizes} must be a resource mapping with explicit ",
        "{.field path} and {.field assembly}."
      )
    )
  }
  if (is.null(spec$assembly)) {
    .acv2g_validation_error(
      "{.field inputs.chrom_sizes.assembly} must explicitly equal analysis_build."
    )
  }
  spec$assembly <- .assert_build(
    spec$assembly, "inputs.chrom_sizes.assembly"
  )
  if (!identical(spec$assembly, analysis_build)) {
    .acv2g_validation_error(
      paste0(
        "inputs.chrom_sizes assembly must equal analysis_build; observed ",
        spec$assembly, " versus ", analysis_build, "."
      )
    )
  }
  .validate_local_resource(
    spec,
    "inputs.chrom_sizes",
    check_files = check_files
  )
}

.validate_activity_chrom_sizes <- function(spec, profile_id, check_files) {
  sizes <- spec$chrom_sizes
  if (is.null(sizes) || is.data.frame(sizes) ||
      (is.atomic(sizes) && !is.character(sizes))) {
    return(spec)
  }
  if (is.character(sizes) && length(sizes) == 1L) {
    sizes <- list(path = sizes, assembly = spec$assembly)
  }
  if (!is.list(sizes)) {
    .acv2g_validation_error(
      paste0(
        "profiles.", profile_id,
        ".activity.chrom_sizes must be a chromosome-length table or local resource."
      )
    )
  }
  sizes$assembly <- sizes$assembly %||% spec$assembly
  sizes$assembly <- .assert_build(
    sizes$assembly,
    paste0("profiles.", profile_id, ".activity.chrom_sizes.assembly")
  )
  if (!identical(sizes$assembly, spec$assembly)) {
    .acv2g_validation_error(
      paste0(
        "profiles.", profile_id,
        ".activity.chrom_sizes assembly must equal the activity assembly."
      )
    )
  }
  sizes <- .validate_local_resource(
    sizes,
    paste0("profiles.", profile_id, ".activity.chrom_sizes"),
    check_files = check_files
  )
  spec$chrom_sizes <- sizes
  spec
}

.normalize_profile_specs <- function(profiles) {
  if (is.null(profiles) || !is.list(profiles) || !length(profiles)) {
    .acv2g_validation_error("Configuration must define at least one named profile.")
  }
  if (is.null(names(profiles)) || any(!nzchar(names(profiles))) ||
      anyDuplicated(names(profiles))) {
    .acv2g_validation_error("Profiles must have unique non-empty IDs.")
  }
  lapply(names(profiles), function(id) {
    profile <- profiles[[id]]
    if (!is.list(profile) || is.null(profile$activity) || is.null(profile$contact)) {
      .acv2g_validation_error(
        "Profile {.val {id}} must contain exactly one activity and one contact specification."
      )
    }
    if (!is.list(profile$activity)) profile$activity <- list(path = profile$activity)
    if (!is.list(profile$contact)) profile$contact <- list(path = profile$contact)
    profile$profile_id <- id
    profile
  })
}

.resolved_defaults <- function(config) {
  defaults <- list(
    schema_version = .acv2g_schema_version,
    method = list(
      name = "locus_pair_z",
      gene_flank_bp = 1000000L,
      reference_gene_radius_bp = 1000000L,
      activity_transform = "log1p",
      contact_transform = "log1p",
      max_missing_fraction = 0.20,
      call_rule = "none",
      z_threshold = 1.5
    ),
    activity = list(
      activity_window_bp = 500L,
      aggregation = "window_mean",
      bam_mapq_min = 30L
    ),
    contact = list(
      contact_value = "balanced",
      contact_resolution_bp = NULL,
      diagonal_policy = "neighbor_max",
      distance_policy = "none"
    ),
    liftover = list(
      structural_failure_policy = "error",
      reciprocal_liftover_policy = "none"
    )
  )
  utils::modifyList(defaults, config, keep.null = TRUE)
}

#' Validate and resolve an ACV2G configuration
#'
#' @param config Configuration list or YAML path.
#' @param check_files Whether local resources must already exist.
#'
#' @return A fully resolved configuration list.
#' @export
validate_v2g_config <- function(config, check_files = TRUE) {
  if (is.character(config) && length(config) == 1L) config <- read_v2g_config(config)
  if (!is.list(config)) {
    .acv2g_validation_error("{.arg config} must be a list or YAML path.")
  }
  if (!is.logical(check_files) || length(check_files) != 1L ||
      is.na(check_files)) {
    .acv2g_validation_error("{.arg check_files} must be one non-missing logical value.")
  }
  .validate_config_key_contract(config)
  if (!is.null(config$schema_version) &&
      (!is.character(config$schema_version) ||
       length(config$schema_version) != 1L || is.na(config$schema_version) ||
       !identical(config$schema_version, .acv2g_schema_version))) {
    .acv2g_validation_error(paste0(
      "Configuration schema_version must exactly equal ",
      .acv2g_schema_version, "; received ",
      paste(config$schema_version, collapse = ", "), "."
    ))
  }
  config <- .resolve_method_transforms(config)
  config <- .resolved_defaults(config)
  config$analysis_build <- .assert_build(config$analysis_build, "analysis_build")
  config$liftover$structural_failure_policy <- .assert_scalar_choice(
    config$liftover$structural_failure_policy,
    c("error", "fail_locus"),
    "liftover.structural_failure_policy"
  )
  config$liftover$reciprocal_liftover_policy <- .assert_scalar_choice(
    config$liftover$reciprocal_liftover_policy,
    c("none", "require"),
    "liftover.reciprocal_liftover_policy"
  )
  if (is.null(config$inputs) || !is.list(config$inputs)) {
    .acv2g_validation_error("Configuration must contain an {.field inputs} mapping.")
  }
  structural_inputs <- c(
    "loci", "variants", "genes",
    if (!is.null(config$inputs$background_variants)) "background_variants"
  )
  for (name in structural_inputs) {
    spec <- .input_spec(config, name, required = name != "background_variants")
    spec <- .validate_local_resource(
      spec,
      paste0("inputs.", name),
      allow_managed = name == "genes",
      check_files = check_files
    )
    spec <- .validate_structural_input_spec(
      spec = spec,
      name = name,
      analysis_build = config$analysis_build,
      reciprocal_policy = config$liftover$reciprocal_liftover_policy,
      check_files = check_files
    )
    config$inputs[[name]] <- spec
  }
  if (!is.null(config$inputs$chrom_sizes)) {
    config$inputs$chrom_sizes <- .validate_analysis_chrom_sizes(
      config$inputs$chrom_sizes,
      analysis_build = config$analysis_build,
      check_files = check_files
    )
  }

  profiles <- .normalize_profile_specs(config$profiles)
  for (i in seq_along(profiles)) {
    profile <- profiles[[i]]
    profile$activity <- .validate_local_resource(
      profile$activity,
      paste0("profiles.", profile$profile_id, ".activity"),
      check_files = check_files
    )
    profile$contact <- .validate_local_resource(
      profile$contact,
      paste0("profiles.", profile$profile_id, ".contact"),
      check_files = check_files
    )
    profile$activity$assembly <- .assert_build(
      profile$activity$assembly,
      paste0("profiles.", profile$profile_id, ".activity.assembly")
    )
    profile$contact$assembly <- .assert_build(
      profile$contact$assembly,
      paste0("profiles.", profile$profile_id, ".contact.assembly")
    )
    if (!identical(profile$contact$assembly, config$analysis_build)) {
      .acv2g_validation_error(
        "Profile {.val {profile$profile_id}} contact assembly must equal analysis_build; contacts are never lifted."
      )
    }
    if (!identical(profile$activity$assembly, config$analysis_build)) {
      direction <- paste0(
        config$analysis_build, " (analysis_build) -> ",
        profile$activity$assembly, " (activity assembly)"
      )
      profile$activity$chain <- .validate_activity_query_chain(
        profile$activity$chain,
        paste0("profiles.", profile$profile_id, ".activity.chain"),
        direction,
        check_files = check_files
      )
      if (identical(
        config$liftover$reciprocal_liftover_policy %||% "none", "require"
      )) {
        reverse_direction <- paste0(
          profile$activity$assembly, " (activity assembly) -> ",
          config$analysis_build, " (analysis_build)"
        )
        profile$activity$reverse_chain <- .validate_activity_query_chain(
          profile$activity$reverse_chain,
          paste0("profiles.", profile$profile_id,
                 ".activity.reverse_chain"),
          reverse_direction,
          check_files = check_files
        )
      }
    }
    profile$activity <- .validate_activity_chrom_sizes(
      profile$activity, profile$profile_id, check_files
    )
    profiles[[i]] <- profile
  }
  names(profiles) <- vapply(profiles, `[[`, character(1L), "profile_id")
  config$profiles <- profiles

  config$method$name <- .assert_scalar_choice(
    config$method$name, c("locus_pair_z", "gene_wise_z"), "method.name"
  )
  config$method$gene_flank_bp <- .assert_scalar_number(
    config$method$gene_flank_bp, "method.gene_flank_bp", lower = 0,
    integer = TRUE
  )
  config$method$reference_gene_radius_bp <- .assert_scalar_number(
    config$method$reference_gene_radius_bp,
    "method.reference_gene_radius_bp", lower = 0, integer = TRUE
  )
  config$method$activity_transform <- .assert_scalar_choice(
    config$method$activity_transform, c("log1p", "identity"),
    "method.activity_transform"
  )
  config$method$contact_transform <- .assert_scalar_choice(
    config$method$contact_transform, c("log1p", "identity"),
    "method.contact_transform"
  )
  config$method$max_missing_fraction <- .assert_scalar_number(
    config$method$max_missing_fraction, "method.max_missing_fraction",
    lower = 0, upper = 1
  )
  config$method$call_rule <- .assert_scalar_choice(
    config$method$call_rule, c("none", "z_threshold"), "method.call_rule"
  )
  if (config$method$call_rule == "z_threshold") {
    config$method$z_threshold <- .assert_scalar_number(
      config$method$z_threshold, "method.z_threshold"
    )
  }
  if (config$method$name == "gene_wise_z") {
    config$method$reference_mode <- .assert_scalar_choice(
      config$method$reference_mode,
      c("input_variants", "background_variants"),
      "method.reference_mode"
    )
    if (config$method$reference_mode == "background_variants" &&
        is.null(config$inputs$background_variants)) {
      .acv2g_validation_error(
        "{.field inputs.background_variants} is required for background reference mode."
      )
    }
    if (!is.null(config$method$reference_set_id) &&
        (!is.character(config$method$reference_set_id) ||
         length(config$method$reference_set_id) != 1L ||
         is.na(config$method$reference_set_id) ||
         !nzchar(config$method$reference_set_id))) {
      .acv2g_validation_error(
        "{.field method.reference_set_id} must be one non-empty identifier."
      )
    }
  }
  config$activity$activity_window_bp <- .assert_scalar_number(
    config$activity$activity_window_bp, "activity.activity_window_bp",
    lower = 1, integer = TRUE
  )
  config$activity$aggregation <- .assert_scalar_choice(
    config$activity$aggregation, c("window_mean", "covered_mean"),
    "activity.aggregation"
  )
  config$activity$bam_mapq_min <- .assert_scalar_number(
    config$activity$bam_mapq_min,
    "activity.bam_mapq_min",
    lower = 0,
    integer = TRUE
  )
  config$contact$contact_value <- .assert_scalar_choice(
    config$contact$contact_value,
    c("balanced", "raw", "named_score"),
    "contact.contact_value"
  )
  if (!is.null(config$contact$contact_resolution_bp)) {
    config$contact$contact_resolution_bp <- .assert_scalar_number(
      config$contact$contact_resolution_bp,
      "contact.contact_resolution_bp", lower = 1, integer = TRUE
    )
  }
  config$contact$diagonal_policy <- .assert_scalar_choice(
    config$contact$diagonal_policy, c("neighbor_max", "observed", "keep_na"),
    "contact.diagonal_policy"
  )
  config$contact$distance_policy <- .assert_scalar_choice(
    config$contact$distance_policy, "none", "contact.distance_policy"
  )
  config$schema_version <- .acv2g_schema_version
  .validate_evidence_source_specs(config, check_files = check_files)
}

#' Write a resolved ACV2G configuration
#'
#' @param config A configuration list or YAML path.
#' @param path Destination YAML path.
#' @param check_files Validate local inputs before writing.
#'
#' @return `path`, invisibly.
#' @export
write_resolved_v2g_config <- function(config, path, check_files = TRUE) {
  resolved <- validate_v2g_config(config, check_files = check_files)
  yaml::write_yaml(resolved, path)
  invisible(path)
}
