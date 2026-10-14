config_validation_fixture <- function(directory, create_files = FALSE) {
  paths <- file.path(
    directory,
    c(
      "loci.tsv", "variants.tsv", "genes.tsv", "activity.bedgraph",
      "contact.cool", "chrom_sizes.tsv"
    )
  )
  if (isTRUE(create_files)) {
    invisible(lapply(paths, function(path) writeLines("placeholder", path)))
  }
  list(
    analysis_build = "GRCh38",
    inputs = list(
      loci = list(path = paths[[1L]], assembly = "GRCh38"),
      variants = list(path = paths[[2L]], assembly = "GRCh38"),
      genes = list(path = paths[[3L]], assembly = "GRCh38"),
      chrom_sizes = list(path = paths[[6L]], assembly = "GRCh38")
    ),
    profiles = list(
      test_profile = list(
        activity = list(path = paths[[4L]], assembly = "GRCh38"),
        contact = list(path = paths[[5L]], assembly = "GRCh38")
      )
    )
  )
}

test_that("direct validation runs source preflight without reading files", {
  directory <- withr::local_tempdir()
  config <- config_validation_fixture(directory)

  resolved <- validate_v2g_config(config, check_files = FALSE)
  expect_identical(resolved$profiles$test_profile$activity$type, "bedgraph")
  expect_identical(resolved$profiles$test_profile$contact$type, "cool")

  invalid_activity <- config
  invalid_activity$profiles$test_profile$activity$table_mode <- "variant"
  expect_error(
    validate_v2g_config(invalid_activity, check_files = FALSE),
    "table_mode is valid only",
    class = "acv2g_validation_error"
  )

  invalid_contact <- config
  invalid_contact$profiles$test_profile$contact$score_col <- "count"
  invalid_contact$profiles$test_profile$contact$score_type <- "raw"
  expect_error(
    validate_v2g_config(invalid_contact, check_files = FALSE),
    "valid only for quantitative contact tables",
    class = "acv2g_validation_error"
  )

  unknown <- config
  unknown$profiles$test_profile$activity$path <- file.path(directory, "signal.unknown")
  expect_error(
    validate_v2g_config(unknown, check_files = FALSE),
    "activity.type is required",
    class = "acv2g_validation_error"
  )
})

test_that("file sampling is skipped only when check_files is false", {
  directory <- withr::local_tempdir()
  config <- config_validation_fixture(directory)
  config$profiles$test_profile$activity <- list(
    path = file.path(directory, "missing-activity.tsv"),
    assembly = "GRCh38",
    type = "table",
    table_mode = "variant",
    value_col = "activity_raw"
  )
  config$profiles$test_profile$contact <- list(
    path = file.path(directory, "missing-contact.tsv"),
    assembly = "GRCh38",
    type = "table",
    score_col = "count",
    score_type = "raw",
    table_representation = "complete",
    contact_value = "raw"
  )

  expect_no_error(validate_v2g_config(config, check_files = FALSE))
  expect_error(
    validate_v2g_config(config, check_files = TRUE),
    "does not exist",
    class = "acv2g_validation_error"
  )

  config$profiles$test_profile$contact$table_representation <- NULL
  expect_error(
    validate_v2g_config(config, check_files = FALSE),
    "table_representation must be explicitly set to complete",
    class = "acv2g_validation_error"
  )
})

test_that("configured variant activity tables require coordinates", {
  directory <- withr::local_tempdir()
  config <- config_validation_fixture(directory, create_files = TRUE)
  activity <- file.path(directory, "activity.tsv")
  contact <- file.path(directory, "contact.tsv")
  data.table::fwrite(
    data.frame(variant_id = "submitted-v1", activity_raw = 2),
    activity,
    sep = "\t"
  )
  data.table::fwrite(
    data.frame(
      chrom1 = "chr1", start1_0 = 0L, end1_0 = 10000L,
      chrom2 = "chr1", start2_0 = 10000L, end2_0 = 20000L,
      count = 3
    ),
    contact,
    sep = "\t"
  )
  config$profiles$test_profile$activity <- list(
    path = activity,
    assembly = "GRCh38",
    type = "table",
    table_mode = "variant",
    value_col = "activity_raw"
  )
  config$profiles$test_profile$contact <- list(
    path = contact,
    assembly = "GRCh38",
    type = "table",
    score_col = "count",
    score_type = "raw",
    table_representation = "complete",
    contact_value = "raw"
  )

  expect_error(
    validate_v2g_config(config),
    "missing required columns?:.*chrom.*pos1",
    class = "acv2g_validation_error"
  )
  data.table::fwrite(
    data.frame(
      variant_id = "submitted-v1", chrom = "chr1", pos1 = 1001L,
      activity_raw = 2
    ),
    activity,
    sep = "\t"
  )
  expect_no_error(validate_v2g_config(config))
})

test_that("every cross-build structural input declares directed chains", {
  directory <- withr::local_tempdir()
  base <- config_validation_fixture(directory)
  roles <- c("loci", "variants", "genes", "background_variants")

  for (role in roles) {
    config <- base
    if (role == "background_variants") {
      config$inputs[[role]] <- list(
        path = file.path(directory, "background.tsv"), assembly = "GRCh37"
      )
    } else {
      config$inputs[[role]]$assembly <- "GRCh37"
    }
    expect_error(
      validate_v2g_config(config, check_files = FALSE),
      paste0(
        "inputs\\.", role, "\\.chain is required and must map ",
        "GRCh37 \\(inputs\\.", role,
        " assembly\\) -> GRCh38 \\(analysis_build\\)"
      ),
      class = "acv2g_validation_error",
      info = role
    )
  }

  config <- base
  config$inputs$variants$assembly <- "GRCh37"
  config$inputs$variants$chain <- file.path(directory, "forward.chain.gz")
  expect_no_error(validate_v2g_config(config, check_files = FALSE))

  config$liftover <- list(reciprocal_liftover_policy = "require")
  expect_error(
    validate_v2g_config(config, check_files = FALSE),
    paste0(
      "inputs\\.variants\\.reverse_chain is required and must map ",
      "GRCh38 \\(analysis_build\\) -> GRCh37 ",
      "\\(inputs\\.variants assembly\\)"
    ),
    class = "acv2g_validation_error"
  )

  fake_chain <- structure(list(), class = "Chain")
  config$inputs$variants$chain <- fake_chain
  config$inputs$variants$reverse_chain <- fake_chain
  resolved <- validate_v2g_config(config, check_files = FALSE)
  expect_s3_class(resolved$inputs$variants$chain, "Chain")
  expect_s3_class(resolved$inputs$variants$reverse_chain, "Chain")
})

test_that("structural chain existence is gated by check_files", {
  directory <- withr::local_tempdir()
  config <- config_validation_fixture(directory, create_files = TRUE)
  config$inputs$loci$assembly <- "GRCh37"
  config$inputs$loci$chain <- file.path(directory, "missing.chain.gz")

  expect_no_error(validate_v2g_config(config, check_files = FALSE))
  expect_error(
    validate_v2g_config(config, check_files = TRUE),
    "inputs.loci liftover chain does not exist",
    class = "acv2g_validation_error"
  )
})

test_that("analysis chromosome sizes are normalized and assembly locked", {
  directory <- withr::local_tempdir()
  base <- config_validation_fixture(directory)
  sizes <- file.path(directory, "GRCh38.chrom.sizes")

  config <- base
  config$inputs$chrom_sizes <- sizes
  expect_error(
    validate_v2g_config(config, check_files = FALSE),
    "resource mapping with explicit",
    class = "acv2g_validation_error"
  )

  config$inputs$chrom_sizes <- list(path = sizes)
  expect_error(
    validate_v2g_config(config, check_files = FALSE),
    "must explicitly equal analysis_build",
    class = "acv2g_validation_error"
  )

  config$inputs$chrom_sizes <- list(path = sizes, assembly = "GRCh37")
  expect_error(
    validate_v2g_config(config, check_files = FALSE),
    "assembly must equal analysis_build",
    class = "acv2g_validation_error"
  )

  config <- config_validation_fixture(directory, create_files = TRUE)
  config$inputs$chrom_sizes <- list(path = sizes, assembly = "GRCh38")
  expect_error(
    validate_v2g_config(config, check_files = TRUE),
    "does not exist",
    class = "acv2g_validation_error"
  )
})

test_that("interval activity requires chromosome lengths in its own assembly", {
  directory <- withr::local_tempdir()
  base <- config_validation_fixture(directory)

  missing_lengths <- base
  missing_lengths$inputs$chrom_sizes <- NULL
  expect_error(
    validate_v2g_config(missing_lengths, check_files = FALSE),
    "requires assembly-matched chromosome lengths",
    class = "acv2g_validation_error"
  )

  cross_build <- base
  cross_build$profiles$test_profile$activity$assembly <- "GRCh37"
  cross_build$profiles$test_profile$activity$chain <-
    file.path(directory, "GRCh38-to-GRCh37.chain.gz")
  expect_error(
    validate_v2g_config(cross_build, check_files = FALSE),
    "requires assembly-matched chromosome lengths",
    class = "acv2g_validation_error"
  )

  cross_build$profiles$test_profile$activity$chrom_sizes <- list(
    path = file.path(directory, "GRCh37.chrom.sizes"),
    assembly = "GRCh37"
  )
  expect_no_error(validate_v2g_config(cross_build, check_files = FALSE))

  cross_build$profiles$test_profile$activity$chrom_sizes$assembly <- "GRCh38"
  expect_error(
    validate_v2g_config(cross_build, check_files = FALSE),
    "activity.chrom_sizes assembly must equal the activity assembly",
    class = "acv2g_validation_error"
  )
})

test_that("coordinate mappings declare and match their coordinate assembly", {
  mapping <- data.frame(
    locus_id = "locus_1", chrom = "chr1", pos1 = 101L,
    stringsAsFactors = FALSE
  )
  expect_error(
    ACV2G:::.read_coordinate_mapping(
      mapping, expected_assembly = "GRCh38", name = "inputs.variant_locus_map"
    ),
    "must explicitly declare assembly",
    class = "acv2g_validation_error"
  )

  declared <- mapping
  declared$assembly <- "GRCh38"
  observed <- ACV2G:::.read_coordinate_mapping(
    declared, expected_assembly = "GRCh38", name = "inputs.variant_locus_map"
  )
  expect_identical(unique(observed$assembly), "GRCh38")

  declared$assembly <- "GRCh37"
  expect_error(
    ACV2G:::.read_coordinate_mapping(
      declared, expected_assembly = "GRCh38", name = "inputs.variant_locus_map"
    ),
    "coordinate assembly must equal GRCh38",
    class = "acv2g_validation_error"
  )

  id_only <- data.frame(
    locus_id = "locus_1", variant_anchor_id = "anchor_1",
    stringsAsFactors = FALSE
  )
  expect_identical(
    ACV2G:::.read_coordinate_mapping(
      id_only, expected_assembly = "GRCh38", name = "inputs.variant_locus_map"
    ),
    id_only
  )

  directory <- withr::local_tempdir()
  config <- config_validation_fixture(directory, create_files = TRUE)
  mapping_path <- file.path(directory, "variant-locus.tsv")
  data.table::fwrite(mapping, mapping_path, sep = "\t")
  config$inputs$variant_locus_map <- list(path = mapping_path)
  expect_error(
    validate_v2g_config(config, check_files = TRUE),
    "must explicitly declare assembly",
    class = "acv2g_validation_error"
  )
  config$inputs$variant_locus_map$assembly <- "GRCh38"
  expect_no_error(validate_v2g_config(config, check_files = TRUE))

  data.table::fwrite(id_only, mapping_path, sep = "\t")
  config$inputs$variant_locus_map <- list(path = mapping_path)
  expect_no_error(validate_v2g_config(config, check_files = TRUE))
})

test_that("schema versions and unknown configuration keys fail closed", {
  directory <- withr::local_tempdir()
  base <- config_validation_fixture(directory)

  incompatible <- base
  incompatible$schema_version <- "0.9.0"
  expect_error(
    validate_v2g_config(incompatible, check_files = FALSE),
    "schema_version must exactly equal 1.0.0",
    class = "acv2g_validation_error"
  )

  cases <- list(
    top = function(x) { x$analysis_buld <- "GRCh38"; x },
    method = function(x) { x$method <- list(gene_flank_bpp = 1L); x },
    input = function(x) { x$inputs$variants$bed_anchor_colum <- 4L; x },
    activity = function(x) {
      x$profiles$test_profile$activity$table_mod <- "variant"; x
    },
    contact = function(x) {
      x$profiles$test_profile$contact$resoluton_bp <- 10000L; x
    }
  )
  for (name in names(cases)) {
    expect_error(
      validate_v2g_config(cases[[name]](base), check_files = FALSE),
      "Unknown .*field.*did you mean",
      class = "acv2g_validation_error",
      info = name
    )
  }

  supported <- base
  supported$inputs$loci$format <- "tsv"
  supported$inputs$loci$column_map <- list(locus_id = "region")
  supported$inputs$variants$format <- "bed"
  supported$inputs$variants$bed_anchor_column <- 4L
  supported$inputs$genes$gene_biotype_field <- "gene_type"
  supported$inputs$genes$id_version_policy <- "preserve"
  expect_no_error(validate_v2g_config(supported, check_files = FALSE))
})

test_that("activity and contact transforms resolve independently", {
  directory <- withr::local_tempdir()
  base <- config_validation_fixture(directory)

  resolved <- validate_v2g_config(base, check_files = FALSE)
  expect_identical(resolved$method$activity_transform, "log1p")
  expect_identical(resolved$method$contact_transform, "log1p")
  expect_null(resolved$method$proxy_transform)

  shorthand <- base
  shorthand$method <- list(proxy_transform = "identity")
  resolved <- validate_v2g_config(shorthand, check_files = FALSE)
  expect_identical(resolved$method$activity_transform, "identity")
  expect_identical(resolved$method$contact_transform, "identity")
  expect_null(resolved$method$proxy_transform)

  independent <- base
  independent$method <- list(activity_transform = "identity")
  resolved <- validate_v2g_config(independent, check_files = FALSE)
  expect_identical(resolved$method$activity_transform, "identity")
  expect_identical(resolved$method$contact_transform, "log1p")

  conflicting <- base
  conflicting$method <- list(
    proxy_transform = "log1p", contact_transform = "identity"
  )
  expect_error(
    validate_v2g_config(conflicting, check_files = FALSE),
    "backward-compatible shorthand only",
    class = "acv2g_validation_error"
  )
})

test_that("YAML-relative cache paths resolve while UNC paths stay absolute", {
  directory <- withr::local_tempdir()
  yaml_path <- file.path(directory, "config.yaml")
  yaml::write_yaml(
    list(reference_cache_dir = "cache", reference_registry = "registry.yaml"),
    yaml_path
  )
  config <- read_v2g_config(yaml_path)
  expect_identical(config$reference_cache_dir, file.path(directory, "cache"))
  expect_identical(
    config$reference_registry, file.path(directory, "registry.yaml")
  )

  unc <- "\\\\server\\share\\cache"
  resolved <- ACV2G:::.resolve_config_paths(
    list(reference_cache_dir = unc), directory
  )
  expect_identical(resolved$reference_cache_dir, unc)
})

test_that("command help renders the executable name", {
  parser <- ACV2G:::.acv2g_cli_config_parser("validate")
  help <- paste(capture.output(optparse::print_help(parser)), collapse = "\n")
  expect_match(help, "Usage: acv2g validate", fixed = TRUE)
  expect_false(grepl("%prog", help, fixed = TRUE))
})
