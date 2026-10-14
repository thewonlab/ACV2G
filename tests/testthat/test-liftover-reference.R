test_that("assembly contract keeps contact space authoritative", {
  plan <- validate_v2g_assemblies(
    "hg38",
    c(contact = "hg38", variants = "hg19"),
    chains = list(variants = "forward.chain.gz")
  )
  expect_equal(plan$operation, c("identity", "liftover"))
  expect_error(
    validate_v2g_assemblies("hg38", c(contact = "hg19")),
    class = "acv2g_validation_error"
  )
  expect_error(
    validate_v2g_assemblies("hg38", c(contact = "hg38", variants = "hg19")),
    class = "acv2g_validation_error"
  )
})

test_that("identity liftover preserves exact point and interval geometry", {
  points <- liftover_v2g_coordinates(
    data.frame(id = "v1", chrom = "chr1", pos1 = 100L),
    chain = NULL,
    source_build = "hg38",
    target_build = "hg38",
    id_col = "id"
  )
  expect_equal(points$data$pos1, 100L)
  expect_equal(points$qc$status, "identity")

  intervals <- liftover_v2g_coordinates(
    data.frame(id = "l1", chrom = "chr2", start0 = 10L, end0 = 20L),
    chain = NULL,
    source_build = "hg38",
    target_build = "hg38",
    record_type = "interval",
    id_col = "id"
  )
  expect_equal(intervals$data$start0, 10L)
  expect_equal(intervals$data$end0, 20L)
})

test_that("identity liftover is vectorized for large gene-TSS tables", {
  n <- 25000L
  input <- data.frame(
    id = paste0("tss_", seq_len(n)),
    chrom = rep(c("chr1", "chr2"), length.out = n),
    pos1 = seq_len(n) + 100L,
    strand = rep(c("+", "-"), length.out = n),
    stringsAsFactors = FALSE
  )
  observed <- liftover_v2g_coordinates(
    input,
    chain = NULL,
    source_build = "hg38",
    target_build = "hg38",
    id_col = "id",
    strand_col = "strand"
  )

  expect_equal(observed$data$chrom, input$chrom)
  expect_equal(observed$data$pos1, input$pos1)
  expect_equal(observed$data$strand, input$strand)
  expect_equal(observed$qc$mapping_count, rep(1L, n))
  expect_true(all(observed$qc$status == "identity"))
})

test_that("managed references are exact, pinned, and offline during resolve", {
  expect_error(
    resolve_v2g_reference("latest", registry = list(references = list())),
    class = "acv2g_validation_error"
  )
  expect_equal(nrow(list_v2g_references(list(references = list()))), 0L)

  source <- tempfile("annotation-", fileext = ".tsv")
  writeLines("gene_id\ttss", source)
  on.exit(unlink(source), add = TRUE)
  sha <- digest::digest(file = source, algo = "sha256", serialize = FALSE)
  registry <- list(references = list(
    test_annotation_v1 = list(
      url = paste0("file://", source),
      sha256 = sha,
      filename = "genes.tsv",
      assembly = "hg38",
      format = "gene_tss"
    )
  ))
  cache <- tempfile("acv2g-cache-")
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  dir.create(file.path(cache, "test_annotation_v1"), recursive = TRUE)
  file.copy(source, file.path(cache, "test_annotation_v1", "genes.tsv"))

  resolved <- resolve_v2g_reference(
    "test_annotation_v1", registry = registry, cache_dir = cache
  )
  expect_equal(resolved$assembly, "hg38")
  expect_true(file.exists(resolved$path))
})
