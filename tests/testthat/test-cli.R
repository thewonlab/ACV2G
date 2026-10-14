capture_acv2g_cli <- function(args) {
  stdout_lines <- character()
  status <- NA_integer_
  stderr_lines <- capture.output(
    stdout_lines <- capture.output(
      status <- acv2g_cli(args),
      type = "output"
    ),
    type = "message"
  )
  list(
    status = status,
    stdout = paste(stdout_lines, collapse = "\n"),
    stderr = paste(stderr_lines, collapse = "\n")
  )
}

read_cli_json <- function(x) {
  jsonlite::fromJSON(x, simplifyVector = TRUE)
}

write_cli_config <- function(directory, method = "locus_pair_z") {
  input_paths <- file.path(
    directory,
    c(
      "loci.tsv", "variants.tsv", "genes.tsv", "activity.bedgraph",
      "contact.cool", "chrom_sizes.tsv"
    )
  )
  invisible(lapply(input_paths, function(path) writeLines(character(), path)))
  config <- list(
    analysis_build = "GRCh38",
    inputs = list(
      loci = list(path = input_paths[[1L]], assembly = "GRCh38"),
      variants = list(path = input_paths[[2L]], assembly = "GRCh38"),
      genes = list(path = input_paths[[3L]], assembly = "GRCh38"),
      chrom_sizes = list(path = input_paths[[6L]], assembly = "GRCh38")
    ),
    profiles = list(
      test_profile = list(
        activity = list(path = input_paths[[4L]], assembly = "GRCh38"),
        contact = list(path = input_paths[[5L]], assembly = "GRCh38")
      )
    ),
    method = list(name = method)
  )
  if (identical(method, "gene_wise_z")) {
    config$method$reference_mode <- "input_variants"
  }
  path <- file.path(directory, paste0(method, ".yaml"))
  yaml::write_yaml(config, path)
  path
}

rewrite_cli_config <- function(path, mutate, name = "modified.yaml") {
  config <- yaml::read_yaml(path)
  config <- mutate(config)
  output <- file.path(dirname(path), name)
  yaml::write_yaml(config, output)
  output
}

cli_error_message <- function(result) {
  read_cli_json(result$stderr)$message
}

write_cli_registry <- function(directory) {
  path <- file.path(directory, "registry.yaml")
  yaml::write_yaml(
    list(references = list(
      test_annotation_v1 = list(
        url = "https://example.invalid/test.gtf.gz",
        sha256 = paste(rep("0", 64L), collapse = ""),
        filename = "test.gtf.gz",
        assembly = "GRCh38",
        format = "gtf"
      )
    )),
    path
  )
  path
}

test_that("methods emits JSON and user-input failures return exit code 2", {
  methods <- capture_acv2g_cli("methods")
  expect_identical(methods$status, 0L)
  expect_identical(methods$stderr, "")
  payload <- read_cli_json(methods$stdout)
  expect_identical(payload$status, "success")
  expect_identical(payload$methods$name, c("locus_pair_z", "gene_wise_z"))
  expect_identical(payload$methods$default, c(TRUE, FALSE))

  invalid <- capture_acv2g_cli("not-a-command")
  expect_identical(invalid$status, 2L)
  expect_identical(invalid$stdout, "")
  error <- read_cli_json(invalid$stderr)
  expect_identical(error$exit_code, 2L)
  expect_identical(error$error_type, "usage")
})

test_that("validate resolves a portable YAML configuration", {
  directory <- withr::local_tempdir()
  config <- write_cli_config(directory)

  result <- capture_acv2g_cli(c("validate", "--config", config))

  expect_identical(result$status, 0L)
  expect_identical(result$stderr, "")
  payload <- read_cli_json(result$stdout)
  expect_identical(payload$status, "valid")
  expect_identical(payload$analysis_build, "GRCh38")
  expect_identical(payload$method, "locus_pair_z")
  expect_identical(payload$profiles, "test_profile")

  missing <- capture_acv2g_cli(c(
    "validate", "--config", file.path(directory, "missing.yaml")
  ))
  expect_identical(missing$status, 2L)
  expect_identical(read_cli_json(missing$stderr)$error_type, "usage")

  malformed_path <- file.path(directory, "malformed.yaml")
  writeLines("profiles: [", malformed_path)
  malformed <- capture_acv2g_cli(c("validate", "--config", malformed_path))
  expect_identical(malformed$status, 2L)
})

test_that("validate requires profile source type when suffix inference is impossible", {
  directory <- withr::local_tempdir()
  base <- write_cli_config(directory)
  unknown_activity <- file.path(directory, "activity.signal")
  writeLines("placeholder", unknown_activity)
  activity_config <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$activity$path <- unknown_activity
    config
  }, "unknown-activity.yaml")

  activity <- capture_acv2g_cli(c("validate", "--config", activity_config))
  expect_identical(activity$status, 2L)
  expect_match(
    cli_error_message(activity),
    "profiles.test_profile.activity.type is required",
    fixed = TRUE
  )

  unknown_contact <- file.path(directory, "contact.matrix")
  writeLines("placeholder", unknown_contact)
  contact_config <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$contact$path <- unknown_contact
    config
  }, "unknown-contact.yaml")
  contact <- capture_acv2g_cli(c("validate", "--config", contact_config))
  expect_identical(contact$status, 2L)
  expect_match(
    cli_error_message(contact),
    "profiles.test_profile.contact.type is required",
    fixed = TRUE
  )

  invalid_type <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$activity$type <- "wiggle"
    config
  }, "invalid-type.yaml")
  invalid <- capture_acv2g_cli(c("validate", "--config", invalid_type))
  expect_identical(invalid$status, 2L)
  expect_match(cli_error_message(invalid), "activity.type")
})

test_that("validate rejects source-specific options that would be ignored", {
  directory <- withr::local_tempdir()
  base <- write_cli_config(directory)
  activity_option <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$activity$table_mode <- "variant"
    config
  }, "invalid-activity-option.yaml")
  activity <- capture_acv2g_cli(c("validate", "--config", activity_option))
  expect_identical(activity$status, 2L)
  expect_match(cli_error_message(activity), "table_mode is valid only")

  contact_option <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$contact$score_col <- "raw"
    config$profiles$test_profile$contact$score_type <- "raw"
    config
  }, "invalid-contact-option.yaml")
  contact <- capture_acv2g_cli(c("validate", "--config", contact_option))
  expect_identical(contact$status, 2L)
  expect_match(cli_error_message(contact), "valid only for quantitative contact tables")

  invalid_diagonal <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$contact$diagonal_policy <- "drop"
    config
  }, "invalid-diagonal.yaml")
  diagonal <- capture_acv2g_cli(c("validate", "--config", invalid_diagonal))
  expect_identical(diagonal$status, 2L)
  expect_match(cli_error_message(diagonal), "diagonal_policy")
})

test_that("validate inspects activity and contact table score mappings", {
  directory <- withr::local_tempdir()
  base <- write_cli_config(directory)
  activity_path <- file.path(directory, "activity.tsv")
  contact_path <- file.path(directory, "contact.tsv")
  data.table::fwrite(
    data.frame(
      variant_id = "v1", chrom = "chr1", pos1 = 1001L,
      activity_raw = 1
    ),
    activity_path,
    sep = "\t"
  )
  data.table::fwrite(
    data.frame(
      chrom1 = "chr1", start1_0 = 0L, end1_0 = 10000L,
      chrom2 = "chr1", start2_0 = 10000L, end2_0 = 20000L,
      raw = 2
    ),
    contact_path,
    sep = "\t"
  )
  table_config <- function(config) {
    config$profiles$test_profile$activity <- list(
      path = activity_path,
      assembly = "GRCh38",
      type = "table",
      table_mode = "variant",
      value_col = "activity_raw"
    )
    config$profiles$test_profile$contact <- list(
      path = contact_path,
      assembly = "GRCh38",
      type = "table",
      score_col = "raw",
      score_type = "raw",
      table_representation = "complete",
      contact_value = "raw"
    )
    config
  }

  valid_path <- rewrite_cli_config(base, table_config, "valid-tables.yaml")
  valid <- capture_acv2g_cli(c("validate", "--config", valid_path))
  expect_identical(valid$status, 0L)

  bad_activity_path <- rewrite_cli_config(base, function(config) {
    config <- table_config(config)
    config$profiles$test_profile$activity$value_col <- "missing_signal"
    config
  }, "bad-activity-mapping.yaml")
  bad_activity <- capture_acv2g_cli(c("validate", "--config", bad_activity_path))
  expect_identical(bad_activity$status, 2L)
  expect_match(cli_error_message(bad_activity), "unavailable column missing_signal")

  bad_contact_path <- rewrite_cli_config(base, function(config) {
    config <- table_config(config)
    config$profiles$test_profile$contact$score_col <- "missing_score"
    config
  }, "bad-contact-mapping.yaml")
  bad_contact <- capture_acv2g_cli(c("validate", "--config", bad_contact_path))
  expect_identical(bad_contact$status, 2L)
  expect_match(cli_error_message(bad_contact), "unavailable contact score.*missing_score")

  no_mapping_path <- rewrite_cli_config(base, function(config) {
    config <- table_config(config)
    config$profiles$test_profile$contact$score_col <- NULL
    config$profiles$test_profile$contact$score_type <- NULL
    config
  }, "missing-contact-mapping.yaml")
  no_mapping <- capture_acv2g_cli(c("validate", "--config", no_mapping_path))
  expect_identical(no_mapping$status, 2L)
  expect_match(cli_error_message(no_mapping), "requires.*score_col.*score_type")

  bad_semantics_path <- rewrite_cli_config(base, function(config) {
    config <- table_config(config)
    config$profiles$test_profile$contact$score_type <- "normalized_somehow"
    config
  }, "bad-contact-semantics.yaml")
  bad_semantics <- capture_acv2g_cli(c(
    "validate", "--config", bad_semantics_path
  ))
  expect_identical(bad_semantics$status, 2L)
  expect_match(cli_error_message(bad_semantics), "balanced, raw, or named_score")
})

test_that("validate checks every optional mapping path before prediction", {
  directory <- withr::local_tempdir()
  base <- write_cli_config(directory)
  for (mapping_name in c("variant_locus_map", "reference_gene_variant_map")) {
    path <- rewrite_cli_config(base, function(config) {
      config$inputs[[mapping_name]] <- file.path(
        directory, paste0("missing-", mapping_name, ".tsv")
      )
      config
    }, paste0("missing-", mapping_name, ".yaml"))
    result <- capture_acv2g_cli(c("validate", "--config", path))
    expect_identical(result$status, 2L, info = mapping_name)
    expect_match(cli_error_message(result), "does not exist", info = mapping_name)
    expect_match(cli_error_message(result), mapping_name, info = mapping_name)
  }
})

test_that("validate enforces cross-build activity chain direction and files", {
  directory <- withr::local_tempdir()
  base <- write_cli_config(directory)
  without_chain <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$activity$assembly <- "hg19"
    config
  }, "cross-build-no-chain.yaml")
  missing <- capture_acv2g_cli(c("validate", "--config", without_chain))
  expect_identical(missing$status, 2L)
  expect_match(cli_error_message(missing), "activity.chain is required")
  expect_match(
    cli_error_message(missing),
    "GRCh38 (analysis_build) -> hg19 (activity assembly)",
    fixed = TRUE
  )

  nonexistent <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$activity$assembly <- "hg19"
    config$profiles$test_profile$activity$chain <- file.path(
      directory, "missing-forward.chain.gz"
    )
    config
  }, "cross-build-missing-chain.yaml")
  absent <- capture_acv2g_cli(c("validate", "--config", nonexistent))
  expect_identical(absent$status, 2L)
  expect_match(cli_error_message(absent), "Activity query chain does not exist")

  forward <- file.path(directory, "analysis-to-activity.chain.gz")
  reverse <- file.path(directory, "activity-to-analysis.chain.gz")
  activity_sizes <- file.path(directory, "hg19.chrom.sizes")
  writeLines(character(), forward)
  writeLines(character(), reverse)
  writeLines(character(), activity_sizes)
  with_forward <- rewrite_cli_config(base, function(config) {
    config$profiles$test_profile$activity$assembly <- "hg19"
    config$profiles$test_profile$activity$chain <- forward
    config$profiles$test_profile$activity$chrom_sizes <- list(
      path = activity_sizes,
      assembly = "hg19"
    )
    config
  }, "cross-build-forward.yaml")
  valid <- capture_acv2g_cli(c("validate", "--config", with_forward))
  expect_identical(valid$status, 0L)

  reciprocal_missing <- rewrite_cli_config(with_forward, function(config) {
    config$liftover$reciprocal_liftover_policy <- "require"
    config
  }, "cross-build-reciprocal-missing.yaml")
  missing_reverse <- capture_acv2g_cli(c(
    "validate", "--config", reciprocal_missing
  ))
  expect_identical(missing_reverse$status, 2L)
  expect_match(cli_error_message(missing_reverse), "activity.reverse_chain is required")
  expect_match(
    cli_error_message(missing_reverse),
    "hg19 (activity assembly) -> GRCh38 (analysis_build)",
    fixed = TRUE
  )

  reciprocal_valid <- rewrite_cli_config(reciprocal_missing, function(config) {
    config$profiles$test_profile$activity$reverse_chain <- reverse
    config
  }, "cross-build-reciprocal-valid.yaml")
  expect_identical(
    capture_acv2g_cli(c("validate", "--config", reciprocal_valid))$status,
    0L
  )
})

test_that("predict is a thin wrapper and forwards overwrite explicitly", {
  directory <- withr::local_tempdir()
  config <- write_cli_config(directory)
  output <- file.path(directory, "result-bundle")
  observed <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .acv2g_cli_predict_api = function(config, output = NULL, overwrite = FALSE) {
      observed$method <- config$method$name
      observed$overwrite <- overwrite
      dir.create(output, recursive = TRUE)
      normalizePath(output, mustWork = TRUE)
    },
    .package = "ACV2G"
  )

  result <- capture_acv2g_cli(c(
    "predict", "--config", config, "--output", output, "--overwrite"
  ))

  expect_identical(result$status, 0L)
  expect_identical(observed$method, "locus_pair_z")
  expect_true(observed$overwrite)
  payload <- read_cli_json(result$stdout)
  expect_identical(payload$command, "predict")
  expect_identical(payload$output, normalizePath(output))
})

test_that("gene-wise-z requires that exact method in the configuration", {
  directory <- withr::local_tempdir()
  locus_config <- write_cli_config(directory, "locus_pair_z")
  rejected <- capture_acv2g_cli(c(
    "gene-wise-z", "--config", locus_config,
    "--output", file.path(directory, "rejected")
  ))
  expect_identical(rejected$status, 2L)
  expect_match(rejected$stderr, "requires method.name: gene_wise_z", fixed = TRUE)

  gene_config <- write_cli_config(directory, "gene_wise_z")
  output <- file.path(directory, "gene-result")
  observed <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .acv2g_cli_predict_api = function(config, output = NULL, overwrite = FALSE) {
      observed$method <- config$method$name
      observed$reference_mode <- config$method$reference_mode
      dir.create(output, recursive = TRUE)
      normalizePath(output, mustWork = TRUE)
    },
    .package = "ACV2G"
  )
  accepted <- capture_acv2g_cli(c(
    "gene-wise-z", "--config", gene_config, "--output", output
  ))
  expect_identical(accepted$status, 0L)
  expect_identical(observed$method, "gene_wise_z")
  expect_identical(observed$reference_mode, "input_variants")
})

test_that("summarize reads bundle metadata without reproducing scoring logic", {
  directory <- withr::local_tempdir()
  result_object <- structure(
    list(
      schema_version = "1.0.0",
      run_status = list(status = "success"),
      config = list(method = list(name = "locus_pair_z")),
      manifest = list(table_rows = list(profiles = 1L, pairs = 3L)),
      tables = list(
        profiles = data.frame(profile_id = "p1"),
        pairs = data.frame(pair_id = paste0("pair", 1:3))
      )
    ),
    class = "v2g_result"
  )
  local_mocked_bindings(
    read_v2g_result = function(path, lazy = FALSE, verify_checksums = TRUE) {
      expect_true(lazy)
      expect_true(verify_checksums)
      result_object
    },
    .package = "ACV2G"
  )

  result <- capture_acv2g_cli(c("summarize", "--result", directory))

  expect_identical(result$status, 0L)
  payload <- read_cli_json(result$stdout)
  expect_identical(payload$run_status, "success")
  expect_identical(payload$profile_ids, "p1")
  expect_identical(payload$table_rows$pairs, 3L)
})

test_that("plot-locus-z forwards only presentation arguments", {
  directory <- withr::local_tempdir()
  result_directory <- file.path(directory, "result")
  dir.create(result_directory)
  output <- file.path(directory, "plots", "locus.pdf")
  observed <- new.env(parent = emptyenv())
  local_mocked_bindings(
    read_v2g_result = function(path, lazy = FALSE, verify_checksums = TRUE) {
      structure(list(), class = "v2g_result")
    },
    plot_locus_z = function(x, display_threshold = NULL, title = NULL,
                            max_points_per_locus = Inf, seed = 1L, ...) {
      observed$threshold <- display_threshold
      observed$title <- title
      observed$max_points <- max_points_per_locus
      observed$seed <- seed
      structure(list(), class = "ggplot")
    },
    save_locus_z_plot = function(plot, path, width = 14, height = 4.5, ...) {
      observed$width <- width
      observed$height <- height
      writeLines("mock plot", path)
      invisible(path)
    },
    .package = "ACV2G"
  )

  result <- capture_acv2g_cli(c(
    "plot-locus-z", "--result", result_directory, "--output", output,
    "--display-threshold", "2", "--title", "My locus plot",
    "--max-points", "20", "--seed", "19", "--width", "8", "--height", "5"
  ))

  expect_identical(result$status, 0L)
  expect_equal(observed$threshold, 2)
  expect_identical(observed$title, "My locus plot")
  expect_equal(observed$max_points, 20)
  expect_identical(observed$seed, 19L)
  expect_equal(observed$width, 8)
  expect_equal(observed$height, 5)
  expect_true(file.exists(output))

  invalid <- capture_acv2g_cli(c(
    "plot-locus-z", "--result", result_directory,
    "--output", file.path(directory, "bad.pdf"), "--max-points", "-Inf"
  ))
  expect_identical(invalid$status, 2L)
})

test_that("reference list and status are offline and fetch is an explicit action", {
  directory <- withr::local_tempdir()
  registry <- write_cli_registry(directory)
  cache <- file.path(directory, "cache")

  listed <- capture_acv2g_cli(c(
    "reference", "list", "--registry", registry
  ))
  expect_identical(listed$status, 0L)
  expect_identical(
    read_cli_json(listed$stdout)$references$managed_id,
    "test_annotation_v1"
  )

  status <- capture_acv2g_cli(c(
    "reference", "status", "--id", "test_annotation_v1",
    "--registry", registry, "--cache-dir", cache
  ))
  expect_identical(status$status, 0L)
  expect_identical(
    read_cli_json(status$stdout)$references$cache_status,
    "not_cached"
  )

  local_mocked_bindings(
    fetch_v2g_reference = function(managed_id, registry = NULL,
                                   cache_dir = NULL, force = FALSE) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      path <- file.path(cache_dir, paste0(managed_id, ".gtf"))
      writeLines("mock reference", path)
      normalizePath(path)
    },
    .package = "ACV2G"
  )
  fetched <- capture_acv2g_cli(c(
    "reference", "fetch", "--id", "test_annotation_v1",
    "--registry", registry, "--cache-dir", cache
  ))
  expect_identical(fetched$status, 0L)
  expect_true(file.exists(read_cli_json(fetched$stdout)$path))
})
