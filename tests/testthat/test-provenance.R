test_that("manifest hashes structural and activity liftover resources", {
  make_resource <- function(label) {
    path <- tempfile()
    writeLines(label, path)
    path
  }
  paths <- vapply(
    c("loci", "variants", "genes", "forward", "reverse", "activity",
      "activity_chain", "activity_reverse", "contact"),
    make_resource,
    character(1L)
  )
  activity_index <- paste0(paths[["activity"]], ".bai")
  writeLines("activity index", activity_index)
  on.exit(unlink(paths), add = TRUE)
  on.exit(unlink(activity_index), add = TRUE)

  config <- list(
    analysis_build = "GRCh38",
    inputs = list(
      loci = list(
        path = paths[["loci"]], assembly = "GRCh37",
        chain = paths[["forward"]], reverse_chain = paths[["reverse"]]
      ),
      variants = list(path = paths[["variants"]], assembly = "GRCh38"),
      genes = list(path = paths[["genes"]], assembly = "GRCh38")
    ),
    profiles = list(
      neuron = list(
        activity = list(
          path = paths[["activity"]], assembly = "GRCh37",
          type = "bam",
          chain = paths[["activity_chain"]],
          reverse_chain = paths[["activity_reverse"]]
        ),
        contact = list(path = paths[["contact"]], assembly = "GRCh38")
      )
    ),
    method = list(name = "locus_pair_z"),
    activity = list(), contact = list(), liftover = list()
  )
  manifest <- build_v2g_manifest(
    config,
    profiles = data.frame(profile_id = "neuron", stringsAsFactors = FALSE)
  )
  resources <- data.table::rbindlist(manifest$input_resources, fill = TRUE) |>
    as.data.frame(stringsAsFactors = FALSE)
  expect_true(all(c(
    "loci", "loci:chain", "loci:reverse_chain",
    "profile:neuron:activity", "profile:neuron:activity:chain",
    "profile:neuron:activity:reverse_chain", "profile:neuron:contact"
  ) %in% resources$role))
  expect_true(any(grepl(
    "^profile:neuron:activity_index:", resources$role
  )))
  expect_true(all(grepl("^[0-9a-f]{64}$", resources$sha256)))
})
