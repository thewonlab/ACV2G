.example_config_path <- function() {
  source_path <- testthat::test_path(
    "..", "..", "inst", "extdata", "example-config.yaml"
  )
  if (file.exists(source_path)) {
    return(normalizePath(source_path, mustWork = TRUE))
  }
  system.file(
    "extdata", "example-config.yaml",
    package = "ACV2G",
    mustWork = TRUE
  )
}

test_that("example configuration runs the complete locus_pair_z workflow", {
  result <- predict_v2g(.example_config_path())
  expect_s3_class(result, "v2g_result")
  expect_equal(result$run_status$status, "success")
  pairs <- v2g_table(result, "pairs")
  expect_equal(nrow(pairs), 3L)
  expect_true(all(is.finite(pairs$locus_pair_z)))
  expect_equal(round(mean(pairs$locus_pair_z), 12), 0)
  expect_equal(round(stats::sd(pairs$locus_pair_z), 12), 1)
  expect_equal(pairs$method, rep("locus_pair_z", 3L))
  expect_equal(nrow(v2g_table(result, "variant_record_anchor_map")), 3L)
  expect_identical(
    unique(v2g_table(result, "variant_record_anchor_map")$input_role),
    "target_variants"
  )
  expect_equal(nrow(v2g_table(result, "variant_locus_qc")), 3L)
  expect_true(all(c("activity_qc", "contact_qc", "normalization_qc") %in%
                    names(result$tables)))
})

test_that("gene_wise_z supports both explicit reference source modes", {
  config <- read_v2g_config(.example_config_path())
  config$method$name <- "gene_wise_z"
  config$method$reference_mode <- "input_variants"
  input_result <- predict_v2g(config)
  input_pairs <- v2g_table(input_result, "pairs")
  expect_true(all(is.finite(input_pairs$gene_wise_z)))
  expect_equal(unique(input_pairs$reference_mode), "input_variants")

  config$method$reference_mode <- "background_variants"
  config$inputs$background_variants <- list(
    path = config$inputs$variants$path,
    format = "tsv",
    assembly = "hg38"
  )
  background_result <- predict_v2g(config)
  background_pairs <- v2g_table(background_result, "pairs")
  expect_equal(unique(background_pairs$reference_mode), "background_variants")
  expect_equal(background_pairs$gene_wise_z, input_pairs$gene_wise_z)
  record_map <- v2g_table(background_result, "variant_record_anchor_map")
  expect_setequal(
    unique(record_map$input_role),
    c("target_variants", "background_variants")
  )
  expect_equal(nrow(record_map), 6L)
  expect_true(all(
    record_map$reference_source_member[
      record_map$input_role == "background_variants"
    ]
  ))
  expect_false(any(
    record_map$reference_source_member[
      record_map$input_role == "target_variants"
    ]
  ))
})

test_that("explicit gene reference mapping completely overrides distance", {
  config <- read_v2g_config(.example_config_path())
  config$method$name <- "gene_wise_z"
  config$method$reference_mode <- "input_variants"
  mapping <- tempfile(fileext = ".tsv")
  on.exit(unlink(mapping), add = TRUE)
  data.table::fwrite(data.frame(
    gene_id_normalized = "ENSG_TEST001",
    chrom = "chr1",
    pos1 = c(100001L, 105001L)
  ), mapping, sep = "\t")
  config$inputs$reference_gene_variant_map <- list(
    path = mapping, assembly = "hg38"
  )
  result <- predict_v2g(config)
  qc <- v2g_table(result, "normalization_qc")
  expect_equal(qc$reference_universe_rule, "explicit_mapping")
  expect_equal(qc$reference_variant_count, 2L)
  expect_true(all(is.finite(v2g_table(result, "pairs")$gene_wise_z)))
})

test_that("explicit reference mappings are filtered to target genes", {
  config <- read_v2g_config(.example_config_path())
  genes <- data.table::fread(config$inputs$genes$path)
  genes <- rbind(
    genes,
    data.frame(
      gene_id = "ENSG_OUTSIDE",
      chrom = "chr1",
      tss_pos1 = 5000001L,
      strand = "+",
      gene_biotype = "protein_coding",
      gene_symbol = "OUTSIDE"
    )
  )
  gene_path <- withr::local_tempfile(fileext = ".tsv")
  data.table::fwrite(genes, gene_path, sep = "\t")
  mapping_path <- withr::local_tempfile(fileext = ".tsv")
  data.table::fwrite(data.frame(
    gene_id_normalized = c(
      "ENSG_TEST001", "ENSG_TEST001", "ENSG_OUTSIDE"
    ),
    chrom = "chr1",
    pos1 = c(100001L, 105001L, 110001L)
  ), mapping_path, sep = "\t")
  config$inputs$genes$path <- gene_path
  config$inputs$reference_gene_variant_map <- list(
    path = mapping_path, assembly = "hg38"
  )
  config$method$name <- "gene_wise_z"
  config$method$reference_mode <- "input_variants"

  result <- predict_v2g(config)

  expect_setequal(
    unique(v2g_table(result, "pairs")$gene_id_normalized),
    "ENSG_TEST001"
  )
  expect_equal(
    unique(v2g_table(result, "normalization_qc")$reference_variant_count),
    2L
  )
})

test_that("automatic reference identities are canonical and universe-sensitive", {
  base <- read_v2g_config(.example_config_path())
  base$method$name <- "gene_wise_z"
  base$method$reference_mode <- "input_variants"
  mapping_a <- withr::local_tempfile(fileext = ".tsv")
  mapping_b <- withr::local_tempfile(fileext = ".tsv")
  rows <- data.frame(
    gene_id_normalized = "ENSG_TEST001",
    chrom = "chr1",
    pos1 = c(100001L, 105001L, 110001L)
  )
  data.table::fwrite(rows, mapping_a, sep = "\t")
  data.table::fwrite(rows[3:1, ], mapping_b, sep = "\t")
  config_a <- base
  config_b <- base
  config_a$inputs$reference_gene_variant_map <- list(
    path = mapping_a, assembly = "hg38"
  )
  config_b$inputs$reference_gene_variant_map <- list(
    path = mapping_b, assembly = "hg38"
  )
  result_a <- predict_v2g(config_a)
  result_b <- predict_v2g(config_b)
  expect_identical(
    result_a$manifest$reference_set_id,
    result_b$manifest$reference_set_id
  )
  expect_identical(
    result_a$manifest$reference_gene_variant_map_id,
    result_b$manifest$reference_gene_variant_map_id
  )

  distance_30k <- base
  distance_25k <- base
  distance_30k$method$reference_gene_radius_bp <- 30000L
  distance_25k$method$reference_gene_radius_bp <- 25000L
  result_30k <- predict_v2g(distance_30k)
  result_25k <- predict_v2g(distance_25k)
  expect_false(identical(
    result_30k$manifest$reference_set_id,
    result_25k$manifest$reference_set_id
  ))
})

test_that("preparation warnings and excluded anchors are retained and rekeyed", {
  config <- read_v2g_config(.example_config_path())
  variants <- data.table::fread(config$inputs$variants$path)
  variants <- rbind(
    variants,
    data.frame(variant_id = "outside", chrom = "chr1", pos1 = 5000000L)
  )
  variant_path <- withr::local_tempfile(fileext = ".tsv")
  data.table::fwrite(variants, variant_path, sep = "\t")
  activity <- data.table::fread(config$profiles$example_profile$activity$path)
  activity <- activity[activity$pos1 != 110001L, ]
  activity_path <- withr::local_tempfile(fileext = ".tsv")
  data.table::fwrite(activity, activity_path, sep = "\t")
  config$inputs$variants$path <- variant_path
  config$profiles$example_profile$activity$path <- activity_path
  config$profiles$replicate_profile <- config$profiles$example_profile

  result <- suppressWarnings(predict_v2g(config))
  variant_qc <- v2g_table(result, "variant_locus_qc")
  warnings <- v2g_table(result, "warnings")

  expect_equal(nrow(variant_qc), 4L)
  expect_equal(
    variant_qc$variant_locus_status[variant_qc$pos1 == 5000000L],
    "outside_all_loci"
  )
  expect_true(any(grepl("outside all fixed loci", warnings$message)))
  expect_true(any(grepl("missing", warnings$message, ignore.case = TRUE)))
  expect_false(anyDuplicated(warnings$warning_id) > 0L)
  expect_identical(
    warnings$warning_id,
    paste0("warning_", sprintf("%06d", seq_len(nrow(warnings))))
  )
})

test_that("R and CLI prediction publish the same pair table", {
  config <- .example_config_path()
  expected <- predict_v2g(config)
  output <- tempfile("acv2g-cli-bundle-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  expect_output(
    status <- acv2g_cli(c(
      "predict", "--config", config, "--output", output
    )),
    '"status":"success"'
  )
  expect_equal(status, 0L)
  observed <- read_v2g_result(output)
  expect_equal(v2g_table(observed, "pairs"), v2g_table(expected, "pairs"),
               ignore_attr = TRUE)
})

test_that("example locus plot saves PDF, PNG, and SVG without a default title", {
  result <- predict_v2g(.example_config_path())
  plot <- plot_locus_z(result)
  expect_null(plot$labels$title)
  for (extension in c("pdf", "png", "svg")) {
    output <- tempfile(fileext = paste0(".", extension))
    on.exit(unlink(output), add = TRUE)
    save_locus_z_plot(plot, output, width = 6, height = 3, dpi = 72)
    expect_gt(file.info(output)$size, 0)
  }
})
