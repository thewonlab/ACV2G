.v2g_required_tables <- c(
  profiles = "profiles.tsv.gz",
  pairs = "pairs.tsv.gz",
  tss_bin_contacts = "tss_bin_contacts.tsv.gz",
  variant_record_anchor_map = "variant_record_anchor_map.tsv.gz",
  variant_locus_qc = "variant_locus_qc.tsv.gz",
  locus_status = "locus_status.tsv.gz",
  normalization_qc = "normalization_qc.tsv.gz",
  activity_qc = "activity_qc.tsv.gz",
  contact_qc = "contact_qc.tsv.gz",
  annotation_qc = "annotation_qc.tsv.gz",
  liftover_qc = "liftover_qc.tsv.gz",
  warnings = "warnings.tsv.gz"
)

.with_schema_version <- function(x) {
  x <- .as_plain_data_frame(x)
  if ("schema_version" %in% names(x)) {
    bad <- !is.na(x$schema_version) & x$schema_version != .acv2g_schema_version
    if (any(bad)) {
      .acv2g_validation_error(
        "A result table contains a schema version incompatible with ACV2G {.val {.acv2g_schema_version}}."
      )
    }
    x$schema_version[is.na(x$schema_version)] <- .acv2g_schema_version
  } else {
    x <- data.frame(
      schema_version = rep(.acv2g_schema_version, nrow(x)),
      x,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  x
}

.result_column_type <- function(x) {
  if (is.factor(x) || is.character(x)) return("character")
  if (is.integer(x)) return("integer")
  if (is.logical(x)) return("logical")
  if (is.numeric(x)) return("double")
  "character"
}

.result_table_schema <- function(x) {
  as.list(vapply(x, .result_column_type, character(1L)))
}

.fread_column_classes <- function(schema) {
  if (is.null(schema) || !length(schema)) return(NULL)
  types <- unlist(schema, use.names = TRUE)
  types[types == "double"] <- "numeric"
  types[!types %in% c("character", "integer", "numeric", "logical")] <-
    "character"
  types
}

.read_result_table <- function(path, schema = NULL) {
  classes <- .fread_column_classes(schema)
  expected_columns <- names(classes)
  fast <- tryCatch(
    suppressWarnings(data.table::fread(
      path,
      na.strings = "NA",
      colClasses = classes,
      check.names = FALSE
    )),
    error = function(e) NULL
  )
  fast_has_schema <- !is.null(fast) &&
    (!length(expected_columns) || all(expected_columns %in% names(fast)))
  if (fast_has_schema) {
    return(as.data.frame(
      fast,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)
  utils::read.delim(
    connection,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    na.strings = "NA",
    colClasses = if (length(classes)) unname(classes) else NA,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.result_default_schema <- function(name) {
  switch(
    name,
    profiles = c(
      profile_id = "character", activity_source = "character",
      contact_source = "character"
    ),
    pairs = c(
      profile_id = "character", pair_id = "character",
      variant_anchor_id = "character", gene_id_normalized = "character"
    ),
    tss_bin_contacts = c(
      profile_id = "character", tss_bin_id = "character",
      pair_id = "character", variant_anchor_id = "character",
      gene_id_normalized = "character"
    ),
    variant_record_anchor_map = c(
      input_role = "character", record_id = "character",
      variant_anchor_id = "character"
    ),
    variant_locus_qc = c(
      variant_anchor_id = "character", variant_key = "character",
      variant_locus_status = "character", locus_membership_count = "integer",
      locus_assignment_mode = "character"
    ),
    locus_status = c(
      profile_id = "character", locus_id = "character", status = "character"
    ),
    normalization_qc = c(
      profile_id = "character", method = "character",
      group_id = "character", status = "character"
    ),
    activity_qc = c(
      profile_id = "character", variant_anchor_id = "character",
      status = "character"
    ),
    contact_qc = c(profile_id = "character", status = "character"),
    annotation_qc = c(
      gene_id_raw = "character", gene_id_normalized = "character",
      status = "character"
    ),
    liftover_qc = c(
      input_role = "character", source_record_id = "character",
      status = "character"
    ),
    warnings = c(
      warning_id = "character", class = "character", message = "character",
      profile_id = "character", locus_id = "character"
    ),
    character()
  )
}

.empty_result_table <- function(schema) {
  columns <- lapply(unname(schema), function(type) {
    switch(
      type,
      character = character(),
      integer = integer(),
      double = double(),
      logical = logical(),
      character()
    )
  })
  names(columns) <- names(schema)
  as.data.frame(
    columns,
    stringsAsFactors = FALSE,
    optional = TRUE,
    check.names = FALSE
  )
}

.default_result_table <- function(name) {
  .with_schema_version(.empty_result_table(.result_default_schema(name)))
}

.result_table_contract <- function(name) {
  key <- switch(
    name,
    profiles = "profile_id",
    pairs = c("profile_id", "pair_id"),
    tss_bin_contacts = c("profile_id", "tss_bin_id"),
    variant_record_anchor_map = c("input_role", "record_id"),
    variant_locus_qc = "variant_anchor_id",
    locus_status = c("profile_id", "locus_id"),
    normalization_qc = c("profile_id", "method", "group_id"),
    character()
  )
  list(
    required = c("schema_version", names(.result_default_schema(name))),
    key = key
  )
}

.result_key_value <- function(x) {
  value <- as.character(x)
  paste0(nchar(value, type = "bytes"), ":", value)
}

.result_contract_issues <- function(x, name) {
  contract <- .result_table_contract(name)
  schema_issues <- character()
  key_issues <- character()
  x <- tryCatch(
    .as_plain_data_frame(x, name),
    error = function(e) e
  )
  if (inherits(x, "error")) {
    detail <- paste0("could not read table: ", conditionMessage(x))
    return(list(schema = detail, key = detail))
  }

  missing_columns <- setdiff(contract$required, names(x))
  if (length(missing_columns)) {
    schema_issues <- c(
      schema_issues,
      paste0("missing required column(s): ", paste(missing_columns, collapse = ", "))
    )
  }
  if ("schema_version" %in% names(x)) {
    version <- as.character(x$schema_version)
    invalid_version <- is.na(version) | !nzchar(version) |
      version != .acv2g_schema_version
    if (any(invalid_version)) {
      schema_issues <- c(
        schema_issues,
        paste0(
          sum(invalid_version),
          " row(s) have a missing or incompatible schema_version"
        )
      )
    }
  }

  missing_key <- setdiff(contract$key, names(x))
  if (length(missing_key)) {
    key_issues <- c(
      key_issues,
      paste0("missing key column(s): ", paste(missing_key, collapse = ", "))
    )
  } else if (length(contract$key) && nrow(x)) {
    invalid_by_column <- vapply(contract$key, function(column) {
      value <- x[[column]]
      if (!is.atomic(value) || is.list(value)) return(TRUE)
      invalid <- is.na(value)
      if (is.character(value) || is.factor(value)) {
        invalid <- invalid | !nzchar(trimws(as.character(value)))
      } else if (is.numeric(value)) {
        invalid <- invalid | !is.finite(value)
      }
      any(invalid)
    }, logical(1L))
    if (any(invalid_by_column)) {
      key_issues <- c(
        key_issues,
        paste0(
          "missing or invalid key value(s) in: ",
          paste(contract$key[invalid_by_column], collapse = ", ")
        )
      )
    } else {
      pieces <- lapply(x[contract$key], .result_key_value)
      composite <- do.call(paste, c(pieces, sep = "|"))
      duplicate_n <- sum(duplicated(composite))
      if (duplicate_n) {
        key_issues <- c(
          key_issues,
          paste0(
            duplicate_n, " duplicate composite key row(s) for ",
            paste(contract$key, collapse = " + ")
          )
        )
      }
    }
  }
  list(schema = schema_issues, key = key_issues)
}

.assert_result_table_contracts <- function(tables) {
  names_to_check <- intersect(names(.v2g_required_tables), names(tables))
  violations <- lapply(names_to_check, function(name) {
    issues <- .result_contract_issues(tables[[name]], name)
    detail <- c(
      if (length(issues$schema)) paste0("schema: ", issues$schema),
      if (length(issues$key)) paste0("key: ", issues$key)
    )
    if (!length(detail)) return(NULL)
    paste0(name, " [", paste(detail, collapse = "; "), "]")
  })
  violations <- unlist(violations, use.names = FALSE)
  if (length(violations)) {
    .acv2g_validation_error(
      paste0(
        "Result table contract validation failed: ",
        paste(violations, collapse = "; "), "."
      ),
      data = violations
    )
  }
  invisible(tables)
}

#' Construct a schema-versioned ACV2G result
#'
#' @param pairs Complete eligible pair table.
#' @param profiles Resolved scoring profiles.
#' @param tss_bin_contacts TSS-bin-level contact detail.
#' @param variant_record_anchor_map Submitted-record-to-anchor mapping.
#' @param variant_locus_qc Per-anchor fixed-locus assignment and exclusion QC.
#' @param locus_status,normalization_qc,activity_qc,contact_qc,annotation_qc,liftover_qc
#'   Status and quality-control tables.
#' @param warnings Structured warning table.
#' @param config Fully resolved configuration list.
#' @param manifest Provenance manifest list.
#' @param run_status Run-status list with `status` equal to `success`, `partial`,
#'   or `failed`.
#' @param called_pairs Optional convenience subset produced by an explicit call
#'   rule. It never replaces `pairs`.
#'
#' @return An object of class `v2g_result`.
#' @export
new_v2g_result <- function(
    pairs,
    profiles = .default_result_table("profiles"),
    tss_bin_contacts = .default_result_table("tss_bin_contacts"),
    variant_record_anchor_map = .default_result_table("variant_record_anchor_map"),
    variant_locus_qc = .default_result_table("variant_locus_qc"),
    locus_status = .default_result_table("locus_status"),
    normalization_qc = .default_result_table("normalization_qc"),
    activity_qc = .default_result_table("activity_qc"),
    contact_qc = .default_result_table("contact_qc"),
    annotation_qc = .default_result_table("annotation_qc"),
    liftover_qc = .default_result_table("liftover_qc"),
    warnings = .default_result_table("warnings"),
    config = list(),
    manifest = list(),
    run_status = list(status = "success"),
    called_pairs = NULL) {
  tables <- list(
    profiles = profiles,
    pairs = pairs,
    tss_bin_contacts = tss_bin_contacts,
    variant_record_anchor_map = variant_record_anchor_map,
    variant_locus_qc = variant_locus_qc,
    locus_status = locus_status,
    normalization_qc = normalization_qc,
    activity_qc = activity_qc,
    contact_qc = contact_qc,
    annotation_qc = annotation_qc,
    liftover_qc = liftover_qc,
    warnings = warnings
  )
  tables <- lapply(tables, .with_schema_version)
  if (!is.null(called_pairs)) {
    tables$called_pairs <- .with_schema_version(called_pairs)
  }
  .assert_result_table_contracts(tables)

  status <- run_status$status %||% NA_character_
  .assert_scalar_choice(status, c("success", "partial", "failed"),
                        "run_status$status")
  run_status$schema_version <- .acv2g_schema_version
  run_status$updated_at_utc <- run_status$updated_at_utc %||% .utc_now()

  manifest$schema_version <- .acv2g_schema_version
  manifest$package <- "ACV2G"
  manifest$package_version <- manifest$package_version %||%
    tryCatch(as.character(utils::packageVersion("ACV2G")),
             error = function(e) "0.0.0.9000")
  manifest$created_at_utc <- manifest$created_at_utc %||% .utc_now()
  manifest$table_rows <- lapply(tables, nrow)
  manifest$table_schema <- lapply(tables, .result_table_schema)

  structure(
    list(
      schema_version = .acv2g_schema_version,
      tables = tables,
      config = config,
      manifest = manifest,
      run_status = run_status,
      file_backed = FALSE,
      bundle_path = NULL
    ),
    class = "v2g_result"
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

#' Access a table in an ACV2G result
#'
#' @param x A `v2g_result`.
#' @param name Table name, such as `pairs` or `normalization_qc`.
#'
#' @return A data frame. File-backed tables are collected on demand.
#' @export
v2g_table <- function(x, name = "pairs") {
  if (!inherits(x, "v2g_result")) {
    .acv2g_validation_error("{.arg x} must be a {.cls v2g_result}.")
  }
  if (!is.character(name) || length(name) != 1L || !name %in% names(x$tables)) {
    .acv2g_validation_error(
      "Unknown result table {.val {name}}. Available tables: {.val {names(x$tables)}}."
    )
  }
  value <- x$tables[[name]]
  if (is.character(value) && length(value) == 1L && file.exists(value)) {
    return(.read_result_table(value, x$manifest$table_schema[[name]] %||% NULL))
  }
  .as_plain_data_frame(value)
}

#' @export
as.data.frame.v2g_result <- function(x, row.names = NULL, optional = FALSE, ...) {
  v2g_table(x, "pairs")
}

#' @export
print.v2g_result <- function(x, ...) {
  pair_n <- if (isTRUE(x$file_backed)) {
    x$manifest$table_rows$pairs %||% NA_integer_
  } else {
    nrow(x$tables$pairs)
  }
  cat("<v2g_result schema ", x$schema_version, ">\n", sep = "")
  cat("status: ", x$run_status$status, "\n", sep = "")
  cat("pairs: ", pair_n, "\n", sep = "")
  cat("storage: ", if (isTRUE(x$file_backed)) "file-backed" else "in-memory",
      "\n", sep = "")
  invisible(x)
}

#' Audit an ACV2G result object or bundle
#'
#' @param x A `v2g_result` or bundle directory.
#' @param verify_checksums Whether to validate `checksums.sha256` for bundles.
#'
#' @return A data frame of checks and statuses.
#' @export
audit_v2g_result <- function(x, verify_checksums = TRUE) {
  checks <- list()
  add_check <- function(check, status, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check, status = status, detail = detail,
      stringsAsFactors = FALSE
    )
  }
  audit_table <- function(name, table = NULL, read_error = NULL) {
    if (!is.null(read_error)) {
      add_check(paste0("table_schema:", name), "fail", read_error)
      add_check(paste0("table_key:", name), "fail", read_error)
      return(invisible(NULL))
    }
    issues <- .result_contract_issues(table, name)
    add_check(
      paste0("table_schema:", name),
      if (length(issues$schema)) "fail" else "pass",
      paste(issues$schema, collapse = "; ")
    )
    add_check(
      paste0("table_key:", name),
      if (length(issues$key)) "fail" else "pass",
      paste(issues$key, collapse = "; ")
    )
    invisible(NULL)
  }

  if (is.character(x) && length(x) == 1L) {
    path <- normalizePath(x, mustWork = TRUE)
    required_files <- c(
      "manifest.json", "config.resolved.yaml", "run_status.json",
      unname(.v2g_required_tables), "checksums.sha256"
    )
    missing_files <- required_files[!file.exists(file.path(path, required_files))]
    add_check("required_files", if (length(missing_files)) "fail" else "pass",
              paste(missing_files, collapse = ","))

    manifest <- NULL
    manifest_path <- file.path(path, "manifest.json")
    if (file.exists(manifest_path)) {
      manifest_read <- tryCatch(
        jsonlite::read_json(manifest_path, simplifyVector = FALSE),
        error = function(e) e
      )
      if (inherits(manifest_read, "error")) {
        add_check("manifest_metadata", "fail", conditionMessage(manifest_read))
      } else {
        manifest <- manifest_read
        manifest_issues <- character()
        if (!identical(manifest$schema_version, .acv2g_schema_version)) {
          manifest_issues <- c(manifest_issues, "incompatible schema_version")
        }
        if (!identical(manifest$package, "ACV2G")) {
          manifest_issues <- c(manifest_issues, "package must be ACV2G")
        }
        missing_schema <- setdiff(
          names(.v2g_required_tables), names(manifest$table_schema %||% list())
        )
        if (length(missing_schema)) {
          manifest_issues <- c(
            manifest_issues,
            paste0(
              "table_schema omits: ", paste(missing_schema, collapse = ", ")
            )
          )
        }
        missing_rows <- setdiff(
          names(.v2g_required_tables), names(manifest$table_rows %||% list())
        )
        if (length(missing_rows)) {
          manifest_issues <- c(
            manifest_issues,
            paste0("table_rows omits: ", paste(missing_rows, collapse = ", "))
          )
        }
        add_check(
          "manifest_metadata",
          if (length(manifest_issues)) "fail" else "pass",
          paste(manifest_issues, collapse = "; ")
        )
      }
    } else {
      add_check("manifest_metadata", "fail", "manifest.json is absent")
    }

    for (name in names(.v2g_required_tables)) {
      table_path <- file.path(path, .v2g_required_tables[[name]])
      if (!file.exists(table_path)) {
        audit_table(name, read_error = "required table file is absent")
        next
      }
      schema <- if (is.null(manifest)) NULL else
        manifest$table_schema[[name]] %||% NULL
      table <- tryCatch(
        .read_result_table(table_path, schema),
        error = function(e) e
      )
      if (inherits(table, "error")) {
        audit_table(
          name,
          read_error = paste0("could not read table: ", conditionMessage(table))
        )
      } else {
        audit_table(name, table)
      }
    }

    if (verify_checksums && !length(missing_files) &&
        file.exists(file.path(path, "checksums.sha256"))) {
      lines <- readLines(file.path(path, "checksums.sha256"), warn = FALSE)
      parsed <- strsplit(lines[nzchar(lines)], "  ", fixed = TRUE)
      structure_ok <- vapply(parsed, function(z) length(z) == 2L, logical(1L))
      listed <- vapply(parsed, function(z) {
        if (length(z) == 2L) z[[2L]] else NA_character_
      }, character(1L))
      expected <- list.files(
        path, recursive = TRUE, full.names = FALSE, all.files = FALSE
      )
      expected <- setdiff(expected, "checksums.sha256")
      coverage_ok <- all(structure_ok) && !anyNA(listed) &&
        !anyDuplicated(listed) && identical(sort(listed), sort(expected))
      ok <- vapply(parsed, function(z) {
        if (length(z) != 2L || !file.exists(file.path(path, z[[2L]]))) return(FALSE)
        identical(digest::digest(file = file.path(path, z[[2L]]),
                                 algo = "sha256", serialize = FALSE), z[[1L]])
      }, logical(1L))
      checksum_pass <- coverage_ok && length(ok) > 0L && all(ok)
      detail <- c(
        if (!coverage_ok) "checksum manifest does not cover each payload exactly once",
        if (length(ok) && any(!ok)) paste0(
          "invalid checksum row(s): ", paste(which(!ok), collapse = ",")
        ),
        if (!length(ok)) "checksum manifest is empty"
      )
      add_check(
        "checksums", if (checksum_pass) "pass" else "fail",
        paste(detail, collapse = "; ")
      )
    }
  } else if (inherits(x, "v2g_result")) {
    missing_tables <- setdiff(names(.v2g_required_tables), names(x$tables))
    add_check("required_tables", if (length(missing_tables)) "fail" else "pass",
              paste(missing_tables, collapse = ","))
    versions <- character()
    for (name in names(.v2g_required_tables)) {
      if (!name %in% names(x$tables)) {
        audit_table(name, read_error = "required result table is absent")
        next
      }
      table <- tryCatch(v2g_table(x, name), error = function(e) e)
      if (inherits(table, "error")) {
        audit_table(
          name,
          read_error = paste0("could not read table: ", conditionMessage(table))
        )
        next
      }
      audit_table(name, table)
      if ("schema_version" %in% names(table)) {
        versions <- c(versions, as.character(table$schema_version))
      }
    }
    versions <- unique(versions)
    versions <- versions[!is.na(versions)]
    add_check("schema_versions",
              if (length(versions) <= 1L &&
                  (!length(versions) || versions == .acv2g_schema_version)) "pass" else "fail",
              paste(versions, collapse = ","))
  } else {
    .acv2g_validation_error("{.arg x} must be a result object or bundle path.")
  }
  do.call(rbind, checks)
}

.validate_output_target <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    .acv2g_validation_error("{.arg path} must be one non-empty path.")
  }
  expanded <- path.expand(path)
  parent <- dirname(expanded)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  parent <- normalizePath(parent, mustWork = TRUE)
  target <- file.path(parent, basename(expanded))
  forbidden <- unique(c(normalizePath("/", mustWork = TRUE), parent))
  if (target %in% forbidden || basename(target) %in% c("", ".", "..")) {
    .acv2g_validation_error("Refusing unsafe output directory {.path {target}}.")
  }
  target
}

.write_tsv_gz <- function(x, path) {
  data.table::fwrite(x, file = path, sep = "\t", quote = TRUE,
                     na = "NA", compress = "gzip")
}

#' Write a transactional ACV2G result bundle
#'
#' @param x A `v2g_result`.
#' @param path Final output directory.
#' @param overwrite Explicitly replace an existing directory. Defaults to
#'   `FALSE`.
#'
#' @return The normalized final directory, invisibly.
#' @export
write_v2g_result <- function(x, path, overwrite = FALSE) {
  if (!inherits(x, "v2g_result")) {
    .acv2g_validation_error("{.arg x} must be a {.cls v2g_result}.")
  }
  object_audit <- audit_v2g_result(x, verify_checksums = FALSE)
  if (any(object_audit$status == "fail")) {
    .acv2g_validation_error(
      "Refusing to write a result object with schema or key contract failures.",
      data = object_audit[object_audit$status == "fail", , drop = FALSE]
    )
  }
  target <- .validate_output_target(path)
  if (file.exists(target) && !isTRUE(overwrite)) {
    .acv2g_validation_error(
      "Output {.path {target}} already exists; use {.arg overwrite = TRUE} explicitly to replace it."
    )
  }

  temp_dir <- tempfile(
    pattern = paste0(".", basename(target), ".tmp-"),
    tmpdir = dirname(target)
  )
  if (!dir.create(temp_dir, recursive = FALSE, showWarnings = FALSE)) {
    .acv2g_abort("Could not create temporary output directory {.path {temp_dir}}.")
  }
  published <- FALSE
  on.exit(if (!published && dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  for (name in names(.v2g_required_tables)) {
    .write_tsv_gz(v2g_table(x, name),
                  file.path(temp_dir, .v2g_required_tables[[name]]))
  }
  if ("called_pairs" %in% names(x$tables)) {
    .write_tsv_gz(v2g_table(x, "called_pairs"),
                  file.path(temp_dir, "called_pairs.tsv.gz"))
  }

  manifest <- x$manifest
  manifest$table_rows <- lapply(x$tables, function(z) {
    if (is.character(z) && length(z) == 1L) nrow(data.table::fread(z)) else nrow(z)
  })
  manifest$table_schema <- lapply(names(x$tables), function(name) {
    .result_table_schema(v2g_table(x, name))
  })
  names(manifest$table_schema) <- names(x$tables)
  yaml::write_yaml(x$config, file.path(temp_dir, "config.resolved.yaml"))
  jsonlite::write_json(x$run_status, file.path(temp_dir, "run_status.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null")

  payload_files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE,
                              all.files = FALSE)
  payload_relative <- substring(payload_files, nchar(temp_dir) + 2L)
  payload_hashes <- vapply(
    payload_files,
    digest::digest,
    character(1L),
    file = TRUE,
    algo = "sha256",
    serialize = FALSE
  )
  manifest$output_files <- lapply(seq_along(payload_files), function(i) {
    list(path = payload_relative[[i]], sha256 = payload_hashes[[i]])
  })
  jsonlite::write_json(manifest, file.path(temp_dir, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null")

  files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE,
                      all.files = FALSE)
  relative <- substring(files, nchar(temp_dir) + 2L)
  hashes <- vapply(files, digest::digest, character(1L), file = TRUE,
                   algo = "sha256", serialize = FALSE)
  writeLines(paste0(hashes, "  ", relative),
             file.path(temp_dir, "checksums.sha256"), useBytes = TRUE)

  audit <- audit_v2g_result(temp_dir, verify_checksums = TRUE)
  if (any(audit$status == "fail")) {
    failed_audit <- audit[audit$status == "fail", , drop = FALSE]
    failed_detail <- paste0(
      failed_audit$check,
      ifelse(nzchar(failed_audit$detail), paste0(": ", failed_audit$detail), ""),
      collapse = "; "
    )
    .acv2g_abort(
      paste0(
        "The staged result bundle failed its integrity audit: ",
        failed_detail, "."
      ),
      data = failed_audit
    )
  }

  if (file.exists(target)) {
    unlink(target, recursive = TRUE, force = TRUE)
  }
  if (!file.rename(temp_dir, target)) {
    .acv2g_abort("Could not publish result bundle to {.path {target}}.")
  }
  published <- TRUE
  invisible(normalizePath(target, mustWork = TRUE))
}

#' Read an ACV2G result bundle
#'
#' @param path Bundle directory.
#' @param lazy If `TRUE`, retain table paths and collect each table through
#'   [v2g_table()].
#' @param verify_checksums Validate bundle checksums before reading.
#'
#' @return A `v2g_result`.
#' @export
read_v2g_result <- function(path, lazy = FALSE, verify_checksums = TRUE) {
  path <- normalizePath(path, mustWork = TRUE)
  audit <- audit_v2g_result(path, verify_checksums = verify_checksums)
  if (any(audit$status == "fail")) {
    .acv2g_validation_error("Result bundle integrity validation failed.", data = audit)
  }
  table_files <- .v2g_required_tables
  if (file.exists(file.path(path, "called_pairs.tsv.gz"))) {
    table_files <- c(table_files, called_pairs = "called_pairs.tsv.gz")
  }
  manifest <- jsonlite::read_json(
    file.path(path, "manifest.json"), simplifyVector = FALSE
  )
  tables <- lapply(names(table_files), function(name) {
    filename <- table_files[[name]]
    full <- file.path(path, filename)
    if (isTRUE(lazy)) full else
      .read_result_table(full, manifest$table_schema[[name]] %||% NULL)
  })
  names(tables) <- names(table_files)
  config <- yaml::read_yaml(file.path(path, "config.resolved.yaml"))
  run_status <- jsonlite::read_json(file.path(path, "run_status.json"), simplifyVector = TRUE)
  structure(
    list(
      schema_version = manifest$schema_version,
      tables = tables,
      config = config,
      manifest = manifest,
      run_status = run_status,
      file_backed = isTRUE(lazy),
      bundle_path = path
    ),
    class = "v2g_result"
  )
}
