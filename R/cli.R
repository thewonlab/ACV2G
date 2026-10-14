.acv2g_cli_top_help <- function() {
  paste(
    "ACV2G: Activity\u2013Contact Variant-to-Gene Prioritization",
    "",
    "Usage:",
    "  acv2g <command> [options]",
    "",
    "Commands:",
    "  methods          List the implemented scoring methods",
    "  validate         Validate and resolve a YAML configuration",
    "  predict          Run the method declared by a YAML configuration",
    "  gene-wise-z      Run a configuration declaring method.name: gene_wise_z",
    "  summarize        Summarize a result bundle",
    "  plot-locus-z     Save a locus_pair_z result plot",
    "  reference        List, fetch, or inspect pinned references",
    "",
    "Global options:",
    "  -h, --help       Show this help",
    "  --version        Print package and schema versions",
    "",
    "Run 'acv2g <command> --help' for command-specific options.",
    sep = "\n"
  )
}

.acv2g_cli_reference_help <- function() {
  paste(
    "Usage:",
    "  acv2g reference list [options]",
    "  acv2g reference fetch --id ID [options]",
    "  acv2g reference status [--id ID] [options]",
    "",
    "Reference fetching is explicit. Prediction never downloads references.",
    sep = "\n"
  )
}

.acv2g_cli_usage_error <- function(message) {
  condition <- simpleError(message, call = NULL)
  class(condition) <- c("acv2g_cli_usage_error", class(condition))
  stop(condition)
}

.acv2g_cli_json <- function(x, stream = stdout()) {
  value <- jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null",
    na = "null",
    digits = NA
  )
  cat(value, "\n", sep = "", file = stream)
  invisible(x)
}

.acv2g_cli_error <- function(condition, exit_code, error_type) {
  .acv2g_cli_json(
    list(
      status = "error",
      exit_code = as.integer(exit_code),
      error_type = error_type,
      message = conditionMessage(condition)
    ),
    stream = stderr()
  )
  as.integer(exit_code)
}

.acv2g_cli_parser <- function(usage, description = NULL, options = list()) {
  if (is.null(description)) description <- ""
  optparse::OptionParser(
    usage = usage,
    description = description,
    option_list = options,
    add_help_option = FALSE
  )
}

.acv2g_cli_parse <- function(parser, args) {
  if (any(args %in% c("-h", "--help"))) {
    optparse::print_help(parser)
    return(structure(list(), class = "acv2g_cli_help"))
  }
  tryCatch(
    optparse::parse_args(
      parser,
      args = args,
      positional_arguments = FALSE,
      print_help_and_exit = FALSE,
      convert_hyphens_to_underscores = TRUE
    ),
    error = function(e) .acv2g_cli_usage_error(conditionMessage(e))
  )
}

.acv2g_cli_option <- function(..., dest, type = "character", default = NULL,
                               action = NULL, help = NULL, metavar = NULL) {
  values <- list(
    ..., dest = dest, type = type, default = default, action = action,
    help = help, metavar = metavar
  )
  values <- values[!vapply(values, is.null, logical(1L))]
  do.call(optparse::make_option, values)
}

.acv2g_cli_required <- function(options, name, display = paste0("--", name)) {
  value <- options[[name]]
  if (is.null(value) || length(value) != 1L || is.na(value) ||
      (is.character(value) && !nzchar(value))) {
    .acv2g_cli_usage_error(paste0("Required option ", display, " was not provided."))
  }
  value
}

.acv2g_cli_existing_path <- function(path, display) {
  expanded <- path.expand(path)
  if (!file.exists(expanded)) {
    .acv2g_cli_usage_error(paste0(display, " does not exist: ", expanded, "."))
  }
  normalizePath(expanded, mustWork = TRUE)
}

.acv2g_cli_file_option <- function(flags, dest, help) {
  .acv2g_cli_option(
    flags, dest = dest, metavar = "FILE", help = help
  )
}

.acv2g_cli_config_parser <- function(command, include_output = FALSE) {
  options <- list(
    .acv2g_cli_file_option(
      c("-c", "--config"), "config", "YAML configuration file"
    )
  )
  if (isTRUE(include_output)) {
    options <- c(
      options,
      list(
        .acv2g_cli_file_option(
          c("-o", "--output"), "output", "Result bundle directory"
        ),
        .acv2g_cli_option(
          "--overwrite", dest = "overwrite", action = "store_true",
          default = FALSE, type = "logical",
          help = "Explicitly replace an existing output bundle"
        )
      )
    )
  }
  .acv2g_cli_parser(
    paste0("acv2g ", command, " [options]"),
    options = options
  )
}

.acv2g_cli_read_config <- function(path) {
  path <- .acv2g_cli_existing_path(path, "Configuration file")
  config <- tryCatch(
    read_v2g_config(path),
    acv2g_validation_error = function(e) stop(e),
    error = function(e) .acv2g_cli_usage_error(paste0(
      "Could not read configuration ", path, ": ", conditionMessage(e)
    ))
  )
  config <- validate_v2g_config(config, check_files = TRUE)
  list(path = path, config = config)
}

.acv2g_cli_optional_input <- function(config, name, check_files = TRUE) {
  spec <- config$inputs[[name]]
  if (is.null(spec)) return(config)
  if (is.character(spec) && length(spec) == 1L) spec <- list(path = spec)
  spec <- .validate_local_resource(
    spec,
    paste0("inputs.", name),
    check_files = check_files
  )
  if (!is.null(spec$assembly)) {
    spec$assembly <- .assert_build(
      spec$assembly, paste0("inputs.", name, ".assembly")
    )
  }
  config$inputs[[name]] <- spec
  config
}

.acv2g_cli_profile_source_type <- function(spec, profile_id, role) {
  field <- paste0("profiles.", profile_id, ".", role, ".type")
  allowed <- if (role == "activity") {
    c("auto", "table", "bigwig", "bedgraph", "bam")
  } else {
    c("auto", "table", "cool", "mcool", "hic")
  }
  declared <- .assert_scalar_choice(spec$type %||% "auto", allowed, field)
  if (declared != "auto") return(declared)
  tryCatch(
    if (role == "activity") {
      .activity_source_type(spec$path, "auto")
    } else {
      .contact_source_type(spec$path, "auto")
    },
    acv2g_validation_error = function(e) {
      .acv2g_validation_error(paste0(
        field, " is required because the source type cannot be inferred from ",
        spec$path, ". Set one of: ",
        paste(setdiff(allowed, "auto"), collapse = ", "), "."
      ))
    }
  )
}

.acv2g_cli_read_profile_table <- function(spec, field) {
  extension <- tolower(tools::file_ext(
    sub("\\.gz$", "", spec$path, ignore.case = TRUE)
  ))
  separator <- if (extension == "csv") "," else "\t"
  tryCatch(
    as.data.frame(
      data.table::fread(
        spec$path,
        sep = separator,
        nrows = 1000L,
        data.table = FALSE,
        check.names = FALSE,
        na.strings = c("NA", "")
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    acv2g_validation_error = function(e) stop(e),
    error = function(e) .acv2g_validation_error(paste0(
      "Could not read ", field, " table ", spec$path, ": ",
      conditionMessage(e)
    ))
  )
}

.acv2g_cli_validate_numeric_column <- function(table, column, field) {
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column)) {
    .acv2g_validation_error(paste0(field, " must name one non-empty column."))
  }
  if (!column %in% names(table)) {
    .acv2g_validation_error(paste0(
      field, " names unavailable column ", column, ". Available columns: ",
      paste(names(table), collapse = ", "), "."
    ))
  }
  raw <- table[[column]]
  numeric <- suppressWarnings(as.numeric(raw))
  if (any(is.na(numeric) & !is.na(raw)) || any(is.infinite(numeric))) {
    .acv2g_validation_error(paste0(
      field, " must map to a numeric column without infinite values."
    ))
  }
  invisible(column)
}

.acv2g_validate_column_name <- function(column, field) {
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column)) {
    .acv2g_validation_error(paste0(field, " must name one non-empty column."))
  }
  invisible(column)
}

.acv2g_cli_validate_activity_table <- function(spec, profile_id) {
  field <- paste0("profiles.", profile_id, ".activity")
  table <- .acv2g_cli_read_profile_table(spec, field)
  table_mode <- .assert_scalar_choice(
    spec$table_mode %||% "auto", c("auto", "interval", "variant"),
    paste0(field, ".table_mode")
  )
  if (table_mode == "auto") {
    table_mode <- if (all(c("chrom", "start0", "end0") %in% names(table))) {
      "interval"
    } else {
      "variant"
    }
  }
  if (table_mode == "interval") {
    .require_columns(table, c("chrom", "start0", "end0"), field)
    value_col <- spec$value_col %||%
      if ("score" %in% names(table)) "score" else NULL
  } else {
    .require_columns(table, c("chrom", "pos1"), field)
    value_col <- spec$value_col %||% if ("activity_raw" %in% names(table)) {
      "activity_raw"
    } else if ("score" %in% names(table)) {
      "score"
    } else {
      NULL
    }
  }
  if (is.null(value_col)) {
    .acv2g_validation_error(paste0(
      field, " requires value_col because no default score/activity_raw column is available."
    ))
  }
  .acv2g_cli_validate_numeric_column(
    table, value_col, paste0(field, ".value_col")
  )
  invisible(table_mode)
}

.acv2g_cli_validate_activity_spec <- function(
    spec,
    profile_id,
    config,
    check_files = TRUE) {
  field <- paste0("profiles.", profile_id, ".activity")
  source_type <- .acv2g_cli_profile_source_type(spec, profile_id, "activity")
  if (source_type == "bam" && grepl("\\.cram$", spec$path, ignore.case = TRUE)) {
    .acv2g_validation_error(paste0(
      field, " points to CRAM input, which is not supported in this release."
    ))
  }
  if (!is.null(spec$activity_window_bp)) {
    .assert_scalar_number(
      spec$activity_window_bp, paste0(field, ".activity_window_bp"),
      lower = 1, integer = TRUE
    )
  }
  aggregation <- .assert_scalar_choice(
    spec$aggregation %||% config$activity$aggregation,
    c("window_mean", "covered_mean"), paste0(field, ".aggregation")
  )
  if (!is.null(spec$uncovered_policy)) {
    expected <- if (aggregation == "window_mean") "zero" else "exclude"
    if (!identical(spec$uncovered_policy, expected)) {
      .acv2g_validation_error(paste0(
        field, ".uncovered_policy must be ", expected,
        " when aggregation is ", aggregation, "."
      ))
    }
  }
  if (!is.null(spec$table_mode) && source_type != "table") {
    .acv2g_validation_error(paste0(
      field, ".table_mode is valid only when activity type is table."
    ))
  }
  if (!is.null(spec$value_col) &&
      !source_type %in% c("table", "bedgraph")) {
    .acv2g_validation_error(paste0(
      field, ".value_col is valid only for table or bedgraph activity."
    ))
  }
  if (!is.null(spec$value_col)) {
    .acv2g_validate_column_name(spec$value_col, paste0(field, ".value_col"))
  }
  bam_options <- c("bam_mapq", "bam_library_size")
  supplied_bam_options <- bam_options[vapply(
    bam_options, function(name) !is.null(spec[[name]]), logical(1L)
  )]
  if (length(supplied_bam_options) && source_type != "bam") {
    .acv2g_validation_error(paste0(
      field, " option(s) ", paste(supplied_bam_options, collapse = ", "),
      " are valid only for BAM activity."
    ))
  }
  if (!is.null(spec$bam_mapq)) {
    .assert_scalar_number(
      spec$bam_mapq, paste0(field, ".bam_mapq"), lower = 0, integer = TRUE
    )
  }
  if (!is.null(spec$bam_library_size)) {
    .assert_scalar_number(
      spec$bam_library_size, paste0(field, ".bam_library_size"),
      lower = 0, inclusive = FALSE
    )
  }
  if (source_type == "table") {
    .assert_scalar_choice(
      spec$table_mode %||% "auto", c("auto", "interval", "variant"),
      paste0(field, ".table_mode")
    )
  }
  if (source_type == "table" && isTRUE(check_files)) {
    .acv2g_cli_validate_activity_table(spec, profile_id)
  }
  spec$type <- source_type
  spec
}

.acv2g_validate_contact_table_declaration <- function(spec, profile_id) {
  field <- paste0("profiles.", profile_id, ".contact")
  if (!is.null(spec$score_type) && is.null(spec$score_col)) {
    .acv2g_validation_error(paste0(
      field, ".score_type requires a corresponding score_col."
    ))
  }
  if (!is.null(spec$score_col)) {
    .acv2g_validate_column_name(spec$score_col, paste0(field, ".score_col"))
  }
  if (!is.null(spec$layers)) {
    valid_layers <- (is.character(spec$layers) || is.list(spec$layers)) &&
      !is.null(names(spec$layers)) && length(spec$layers) > 0L &&
      all(nzchar(names(spec$layers))) &&
      all(names(spec$layers) %in% c("balanced", "raw", "named_score"))
    if (!valid_layers) {
      .acv2g_validation_error(paste0(
        field, ".layers must be a named character vector/list keyed by ",
        "balanced, raw, or named_score."
      ))
    }
    layer_values <- unlist(lapply(spec$layers, as.character), use.names = FALSE)
    if (!length(layer_values) || anyNA(layer_values) ||
        any(!nzchar(layer_values))) {
      .acv2g_validation_error(paste0(
        field, ".layers entries must name non-empty score columns."
      ))
    }
  }
  declared_layers <- unique(c(
    spec$score_col %||% character(),
    if (is.list(spec$layers) || is.character(spec$layers)) {
      unlist(spec$layers, use.names = FALSE)
    } else {
      character()
    }
  ))
  .contact_layer_map(
    declared_layers,
    contact_score_col = spec$score_col %||% NULL,
    contact_score_type = spec$score_type %||% NULL,
    contact_layers = spec$layers %||% NULL
  )
}

.acv2g_cli_validate_contact_table <- function(spec, profile_id, config) {
  field <- paste0("profiles.", profile_id, ".contact")
  table <- .acv2g_cli_read_profile_table(spec, field)
  .require_columns(
    table,
    c("chrom1", "start1_0", "end1_0", "chrom2", "start2_0", "end2_0"),
    field
  )
  layer_roles <- .acv2g_validate_contact_table_declaration(spec, profile_id)
  declared_layers <- layer_roles$layer
  if (length(declared_layers) && is.character(declared_layers)) {
    missing_layers <- setdiff(declared_layers, names(table))
    if (length(missing_layers)) {
      .acv2g_validation_error(paste0(
        field, " maps unavailable contact score column(s): ",
        paste(missing_layers, collapse = ", "), ". Available columns: ",
        paste(names(table), collapse = ", "), "."
      ))
    }
  }
  for (layer in layer_roles$layer) {
    .acv2g_cli_validate_numeric_column(
      table, layer, paste0(field, ".score mapping")
    )
  }
  contact_value <- .assert_scalar_choice(
    spec$contact_value %||% config$contact$contact_value,
    c("balanced", "raw", "named_score"), paste0(field, ".contact_value")
  )
  withCallingHandlers(
    .contact_resolve_layer(
      layer_roles,
      contact_value = contact_value,
      contact_layer = spec[["layer"]] %||% NULL
    ),
    warning = function(w) invokeRestart("muffleWarning")
  )
  invisible(layer_roles)
}

.acv2g_cli_validate_contact_spec <- function(
    spec,
    profile_id,
    config,
    check_files = TRUE) {
  field <- paste0("profiles.", profile_id, ".contact")
  source_type <- .acv2g_cli_profile_source_type(spec, profile_id, "contact")
  .assert_scalar_choice(
    spec$contact_value %||% config$contact$contact_value,
    c("balanced", "raw", "named_score"), paste0(field, ".contact_value")
  )
  if (!is.null(spec$resolution_bp)) {
    .assert_scalar_number(
      spec$resolution_bp, paste0(field, ".resolution_bp"),
      lower = 1, integer = TRUE
    )
  }
  .assert_scalar_choice(
    spec$diagonal_policy %||% config$contact$diagonal_policy,
    c("neighbor_max", "observed", "keep_na"),
    paste0(field, ".diagonal_policy")
  )
  if (source_type == "table" && is.null(spec$table_representation)) {
    .acv2g_validation_error(paste0(
      field,
      ".table_representation must be explicitly set to complete for tabular contact input."
    ))
  }
  if (!is.null(spec$table_representation) &&
      !identical(spec$table_representation, "complete")) {
    .acv2g_validation_error(paste0(
      field, ".table_representation must be complete; sparse loop lists are unsupported."
    ))
  }
  if (!is.null(spec[["layer"]]) &&
      (!is.character(spec[["layer"]]) || length(spec[["layer"]]) != 1L ||
       is.na(spec[["layer"]]) || !nzchar(spec[["layer"]]))) {
    .acv2g_validation_error(paste0(field, ".layer must be one non-empty name."))
  }
  table_only <- c(
    "score_col", "score_type", "layers", "weights_source_provided",
    "table_representation"
  )
  supplied_table_only <- table_only[vapply(
    table_only, function(name) !is.null(spec[[name]]), logical(1L)
  )]
  if (length(supplied_table_only) && source_type != "table") {
    .acv2g_validation_error(paste0(
      field, " option(s) ", paste(supplied_table_only, collapse = ", "),
      " are valid only for quantitative contact tables."
    ))
  }
  if (!is.null(spec$weights_source_provided) &&
      (!is.logical(spec$weights_source_provided) ||
       length(spec$weights_source_provided) != 1L ||
       is.na(spec$weights_source_provided))) {
    .acv2g_validation_error(paste0(
      field, ".weights_source_provided must be one non-missing logical value."
    ))
  }
  if (source_type == "table") {
    layer_roles <- .acv2g_validate_contact_table_declaration(spec, profile_id)
    withCallingHandlers(
      .contact_resolve_layer(
        layer_roles,
        contact_value = spec$contact_value %||% config$contact$contact_value,
        contact_layer = spec[["layer"]] %||% NULL
      ),
      warning = function(w) invokeRestart("muffleWarning")
    )
  }
  if (source_type == "table" && isTRUE(check_files)) {
    .acv2g_cli_validate_contact_table(spec, profile_id, config)
  }
  spec$type <- source_type
  spec
}

.validate_evidence_source_specs <- function(config, check_files = TRUE) {
  for (name in c("variant_locus_map", "reference_gene_variant_map")) {
    config <- .acv2g_cli_optional_input(
      config, name, check_files = check_files
    )
  }
  mapping_assemblies <- c(
    variant_locus_map = config$inputs$variants$assembly,
    reference_gene_variant_map = config$analysis_build
  )
  for (name in names(mapping_assemblies)) {
    spec <- config$inputs[[name]]
    if (is.null(spec)) next
    expected_assembly <- unname(mapping_assemblies[[name]])
    if (!is.null(spec$assembly) &&
        !identical(spec$assembly, expected_assembly)) {
      .acv2g_validation_error(paste0(
        "inputs.", name, " assembly must equal ", expected_assembly,
        "; observed ", spec$assembly, "."
      ))
    }
    if (isTRUE(check_files)) {
      .read_coordinate_mapping(
        spec,
        expected_assembly = expected_assembly,
        name = paste0("inputs.", name)
      )
    }
  }
  for (profile_id in names(config$profiles)) {
    profile <- config$profiles[[profile_id]]
    profile$activity <- .acv2g_cli_validate_activity_spec(
      profile$activity, profile_id, config, check_files = check_files
    )
    activity_needs_lengths <- profile$activity$type %in% c("table", "bedgraph")
    has_profile_lengths <- !is.null(profile$activity$chrom_sizes)
    has_analysis_lengths <- identical(
      profile$activity$assembly, config$analysis_build
    ) && !is.null(config$inputs$chrom_sizes)
    if (activity_needs_lengths &&
        !has_profile_lengths && !has_analysis_lengths) {
      .acv2g_validation_error(paste0(
        "profiles.", profile_id, ".activity requires assembly-matched ",
        "chromosome lengths. Supply profiles.", profile_id,
        ".activity.chrom_sizes, or inputs.chrom_sizes when activity assembly ",
        "equals analysis_build."
      ))
    }
    profile$contact <- .acv2g_cli_validate_contact_spec(
      profile$contact, profile_id, config, check_files = check_files
    )
    config$profiles[[profile_id]] <- profile
  }
  config
}

.acv2g_cli_predict_api <- function(config, output, overwrite) {
  predict_v2g_from_config(
    config = config,
    output = output,
    overwrite = overwrite
  )
}

.acv2g_cli_methods <- function(args) {
  parser <- .acv2g_cli_parser("acv2g methods")
  options <- .acv2g_cli_parse(parser, args)
  if (inherits(options, "acv2g_cli_help")) return(0L)
  .acv2g_cli_json(list(
    status = "success",
    methods = data.frame(
      name = c("locus_pair_z", "gene_wise_z"),
      default = c(TRUE, FALSE),
      scope = c("fixed_locus", "explicit_gene_reference_universe"),
      stringsAsFactors = FALSE
    )
  ))
  0L
}

.acv2g_cli_validate <- function(args) {
  parser <- .acv2g_cli_config_parser("validate")
  options <- .acv2g_cli_parse(parser, args)
  if (inherits(options, "acv2g_cli_help")) return(0L)
  config_path <- .acv2g_cli_required(options, "config")
  resolved <- .acv2g_cli_read_config(config_path)
  config <- resolved$config
  .acv2g_cli_json(list(
    status = "valid",
    config = resolved$path,
    schema_version = config$schema_version,
    analysis_build = config$analysis_build,
    method = config$method$name,
    profiles = names(config$profiles)
  ))
  0L
}

.acv2g_cli_prediction <- function(args, command, expected_method = NULL) {
  parser <- .acv2g_cli_config_parser(command, include_output = TRUE)
  options <- .acv2g_cli_parse(parser, args)
  if (inherits(options, "acv2g_cli_help")) return(0L)
  config_path <- .acv2g_cli_required(options, "config")
  output <- .acv2g_cli_required(options, "output")
  resolved <- .acv2g_cli_read_config(config_path)
  config <- resolved$config
  if (!is.null(expected_method) && !identical(config$method$name, expected_method)) {
    .acv2g_cli_usage_error(paste0(
      "The ", command, " command requires method.name: ", expected_method,
      " in the configuration; observed ", config$method$name, "."
    ))
  }

  value <- .acv2g_cli_predict_api(
    config = config,
    output = output,
    overwrite = isTRUE(options$overwrite)
  )
  if (inherits(value, "v2g_result")) {
    if (!dir.exists(output)) {
      write_v2g_result(value, output, overwrite = isTRUE(options$overwrite))
    }
    bundle <- normalizePath(output, mustWork = TRUE)
    run_status <- value$run_status$status %||% "success"
  } else {
    bundle <- if (is.character(value) && length(value) == 1L && nzchar(value)) {
      normalizePath(value, mustWork = TRUE)
    } else {
      normalizePath(output, mustWork = TRUE)
    }
    run_status <- "success"
  }
  .acv2g_cli_json(list(
    status = run_status,
    command = command,
    method = config$method$name,
    config = resolved$path,
    output = bundle
  ))
  0L
}

.acv2g_cli_table_rows <- function(result) {
  counts <- result$manifest$table_rows
  if (is.null(counts)) counts <- list()
  counts <- as.list(counts)
  absent <- setdiff(names(result$tables), names(counts))
  for (name in absent) {
    counts[[name]] <- nrow(v2g_table(result, name))
  }
  counts
}

.acv2g_cli_summarize <- function(args) {
  parser <- .acv2g_cli_parser(
    "acv2g summarize [options]",
    options = list(
      .acv2g_cli_file_option(
        c("-r", "--result"), "result", "ACV2G result bundle directory"
      )
    )
  )
  options <- .acv2g_cli_parse(parser, args)
  if (inherits(options, "acv2g_cli_help")) return(0L)
  path <- .acv2g_cli_existing_path(
    .acv2g_cli_required(options, "result"), "Result bundle"
  )
  result <- read_v2g_result(path, lazy = TRUE, verify_checksums = TRUE)
  profiles <- v2g_table(result, "profiles")
  .acv2g_cli_json(list(
    status = "success",
    result = path,
    schema_version = result$schema_version,
    run_status = result$run_status$status,
    method = result$config$method$name %||% NULL,
    profile_ids = if ("profile_id" %in% names(profiles)) {
      unique(as.character(profiles$profile_id))
    } else {
      character()
    },
    table_rows = .acv2g_cli_table_rows(result)
  ))
  0L
}

.acv2g_cli_plot <- function(args) {
  parser <- .acv2g_cli_parser(
    "acv2g plot-locus-z [options]",
    options = list(
      .acv2g_cli_file_option(
        c("-r", "--result"), "result", "ACV2G result bundle directory"
      ),
      .acv2g_cli_file_option(
        c("-o", "--output"), "output", "Plot path (.pdf, .png, or .svg)"
      ),
      .acv2g_cli_option(
        "--profile", dest = "profile_id", type = "character",
        metavar = "PROFILE_ID",
        help = "Profile to plot (required when the result has multiple profiles)"
      ),
      .acv2g_cli_option(
        "--display-threshold", dest = "display_threshold", type = "double",
        metavar = "NUMBER", help = "Optional visual-only Z-score threshold"
      ),
      .acv2g_cli_option(
        "--title", dest = "title", type = "character", metavar = "TEXT",
        help = "Optional user-supplied plot title"
      ),
      .acv2g_cli_option(
        "--point-mode", dest = "point_mode", type = "character",
        metavar = "MODE", help = "Point display: all, sample, or none"
      ),
      .acv2g_cli_option(
        "--max-points", dest = "max_points", type = "double", default = Inf,
        metavar = "N", help = "Maximum jitter points per locus (default: all)"
      ),
      .acv2g_cli_option(
        "--seed", dest = "seed", type = "integer", default = 1L,
        metavar = "INTEGER", help = "Deterministic point-sampling seed"
      ),
      .acv2g_cli_option(
        "--width", dest = "width", type = "double", default = 14,
        metavar = "INCHES", help = "Plot width in inches"
      ),
      .acv2g_cli_option(
        "--height", dest = "height", type = "double", default = 4.5,
        metavar = "INCHES", help = "Plot height in inches"
      )
    )
  )
  options <- .acv2g_cli_parse(parser, args)
  if (inherits(options, "acv2g_cli_help")) return(0L)
  result_path <- .acv2g_cli_existing_path(
    .acv2g_cli_required(options, "result"), "Result bundle"
  )
  output <- .acv2g_cli_required(options, "output")
  point_mode <- options$point_mode %||%
    if (is.finite(options$max_points)) "sample" else "all"
  if (!point_mode %in% c("all", "sample", "none")) {
    .acv2g_cli_usage_error("--point-mode must be all, sample, or none.")
  }
  if (identical(point_mode, "sample") && !is.finite(options$max_points)) {
    .acv2g_cli_usage_error(
      "--point-mode sample requires a finite positive --max-points value."
    )
  }
  max_points_ok <- (is.infinite(options$max_points) && options$max_points > 0) ||
    (is.finite(options$max_points) && options$max_points >= 1 &&
     options$max_points == floor(options$max_points))
  if (!max_points_ok || is.na(options$seed) ||
      !is.finite(options$width) || options$width <= 0 ||
      !is.finite(options$height) || options$height <= 0) {
    .acv2g_cli_usage_error(
      "--max-points must be a positive integer or Inf, and plot dimensions must be positive."
    )
  }
  output_parent <- dirname(path.expand(output))
  if (!dir.exists(output_parent)) {
    dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
  }
  result <- read_v2g_result(result_path, lazy = TRUE, verify_checksums = TRUE)
  plot <- plot_locus_z(
    result,
    profile_id = options$profile_id,
    display_threshold = options$display_threshold,
    title = options$title,
    point_mode = point_mode,
    max_points_per_locus = options$max_points,
    seed = options$seed
  )
  save_locus_z_plot(plot, output, width = options$width, height = options$height)
  output <- normalizePath(output, mustWork = TRUE)
  .acv2g_cli_json(list(
    status = "success",
    command = "plot-locus-z",
    result = result_path,
    output = output,
    profile_id = options$profile_id,
    display_threshold = options$display_threshold,
    title = options$title,
    point_mode = point_mode,
    max_points_per_locus = options$max_points,
    seed = options$seed,
    width = options$width,
    height = options$height
  ))
  0L
}

.acv2g_cli_reference_common_options <- function() {
  list(
    .acv2g_cli_file_option(
      "--registry", "registry", "Optional pinned reference-registry YAML"
    ),
    .acv2g_cli_file_option(
      "--cache-dir", "cache_dir", "Optional platform-independent cache root"
    )
  )
}

.acv2g_cli_reference_list <- function(args) {
  parser <- .acv2g_cli_parser(
    "acv2g reference list [options]",
    options = .acv2g_cli_reference_common_options()[1L]
  )
  options <- .acv2g_cli_parse(parser, args)
  if (inherits(options, "acv2g_cli_help")) return(0L)
  references <- list_v2g_references(registry = options$registry)
  .acv2g_cli_json(list(status = "success", references = references))
  0L
}

.acv2g_cli_reference_fetch <- function(args) {
  parser <- .acv2g_cli_parser(
    "acv2g reference fetch --id ID [options]",
    options = c(
      list(
        .acv2g_cli_option(
          "--id", dest = "id", type = "character", metavar = "ID",
          help = "Exact pinned managed-reference ID"
        )
      ),
      .acv2g_cli_reference_common_options(),
      list(
        .acv2g_cli_option(
          "--force", dest = "force", action = "store_true", default = FALSE,
          type = "logical", help = "Replace an existing cache entry after verification"
        )
      )
    )
  )
  options <- .acv2g_cli_parse(parser, args)
  if (inherits(options, "acv2g_cli_help")) return(0L)
  managed_id <- .acv2g_cli_required(options, "id", "--id")
  path <- fetch_v2g_reference(
    managed_id = managed_id,
    registry = options$registry,
    cache_dir = options$cache_dir,
    force = isTRUE(options$force)
  )
  .acv2g_cli_json(list(
    status = "success",
    managed_id = managed_id,
    path = as.character(path)
  ))
  0L
}

.acv2g_cli_reference_status_table <- function(registry = NULL, cache_dir = NULL,
                                               managed_id = NULL) {
  references <- list_v2g_references(registry = registry)
  if (!is.null(managed_id)) {
    if (!managed_id %in% references$managed_id) {
      .acv2g_cli_usage_error(paste0(
        "Unknown managed reference ", managed_id, ". Available IDs: ",
        paste(references$managed_id, collapse = ", "), "."
      ))
    }
    references <- references[references$managed_id == managed_id, , drop = FALSE]
  }
  rows <- lapply(references$managed_id, function(id) {
    resolved <- tryCatch(
      resolve_v2g_reference(id, registry = registry, cache_dir = cache_dir),
      acv2g_validation_error = identity
    )
    if (inherits(resolved, "error")) {
      message <- conditionMessage(resolved)
      data.frame(
        managed_id = id,
        cache_status = if (grepl("not cached", message, fixed = TRUE)) {
          "not_cached"
        } else {
          "invalid"
        },
        path = NA_character_,
        message = message,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        managed_id = id,
        cache_status = "cached",
        path = as.character(resolved$path),
        message = "",
        stringsAsFactors = FALSE
      )
    }
  })
  if (!length(rows)) {
    return(data.frame(
      managed_id = character(), cache_status = character(),
      path = character(), message = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

.acv2g_cli_reference_status <- function(args) {
  parser <- .acv2g_cli_parser(
    "acv2g reference status [options]",
    options = c(
      list(
        .acv2g_cli_option(
          "--id", dest = "id", type = "character", metavar = "ID",
          help = "Inspect one exact pinned ID (default: all)"
        )
      ),
      .acv2g_cli_reference_common_options()
    )
  )
  options <- .acv2g_cli_parse(parser, args)
  if (inherits(options, "acv2g_cli_help")) return(0L)
  status <- .acv2g_cli_reference_status_table(
    registry = options$registry,
    cache_dir = options$cache_dir,
    managed_id = options$id
  )
  .acv2g_cli_json(list(status = "success", references = status))
  0L
}

.acv2g_cli_reference <- function(args) {
  if (!length(args) || args[[1L]] %in% c("-h", "--help")) {
    cat(.acv2g_cli_reference_help(), "\n", sep = "")
    return(0L)
  }
  subcommand <- args[[1L]]
  rest <- args[-1L]
  switch(
    subcommand,
    list = .acv2g_cli_reference_list(rest),
    fetch = .acv2g_cli_reference_fetch(rest),
    status = .acv2g_cli_reference_status(rest),
    .acv2g_cli_usage_error(paste0(
      "Unknown reference command '", subcommand,
      "'. Expected list, fetch, or status."
    ))
  )
}

.acv2g_cli_version <- function() {
  version <- tryCatch(
    as.character(utils::packageVersion("ACV2G")),
    error = function(e) "0.0.0.9000"
  )
  .acv2g_cli_json(list(
    package = "ACV2G",
    version = version,
    schema_version = acv2g_schema_version()
  ))
  0L
}

.acv2g_cli_dispatch <- function(args) {
  if (!length(args) || args[[1L]] %in% c("-h", "--help")) {
    cat(.acv2g_cli_top_help(), "\n", sep = "")
    return(0L)
  }
  if (identical(args[[1L]], "--version")) {
    if (length(args) != 1L) {
      .acv2g_cli_usage_error("--version does not accept additional arguments.")
    }
    return(.acv2g_cli_version())
  }
  command <- args[[1L]]
  rest <- args[-1L]
  switch(
    command,
    methods = .acv2g_cli_methods(rest),
    validate = .acv2g_cli_validate(rest),
    predict = .acv2g_cli_prediction(rest, "predict"),
    `gene-wise-z` = .acv2g_cli_prediction(
      rest, "gene-wise-z", expected_method = "gene_wise_z"
    ),
    summarize = .acv2g_cli_summarize(rest),
    `plot-locus-z` = .acv2g_cli_plot(rest),
    reference = .acv2g_cli_reference(rest),
    .acv2g_cli_usage_error(paste0(
      "Unknown command '", command, "'. Run 'acv2g --help' for usage."
    ))
  )
}

#' Run the ACV2G command-line interface
#'
#' This is a thin, portable dispatcher over the public ACV2G R API. Successful
#' commands write machine-readable JSON to standard output. Diagnostics are
#' written to standard error. The function returns an integer exit code rather
#' than terminating R, which makes the same interface usable from tests and
#' embedded R sessions.
#'
#' @param args Character vector of command-line arguments, excluding the
#'   executable name.
#'
#' @return Exit code `0L` on success, `2L` for invalid user input or validation
#'   failure, and `1L` for an unexpected failure, invisibly.
#' @export
acv2g_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  old_options <- options(cli.num_colors = 0)
  on.exit(options(old_options), add = TRUE)
  if (!is.character(args) || anyNA(args)) {
    return(invisible(.acv2g_cli_error(
      simpleError("args must be a character vector without missing values."),
      2L,
      "usage"
    )))
  }
  status <- tryCatch(
    .acv2g_cli_dispatch(args),
    acv2g_cli_usage_error = function(e) .acv2g_cli_error(e, 2L, "usage"),
    acv2g_validation_error = function(e) {
      .acv2g_cli_error(e, 2L, "validation")
    },
    error = function(e) .acv2g_cli_error(e, 1L, "unexpected")
  )
  invisible(as.integer(status))
}
