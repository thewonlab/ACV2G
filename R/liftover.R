.assert_build <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .acv2g_validation_error("{.arg {arg}} must explicitly name one genome assembly.")
  }
  x
}

.read_chain <- function(chain, arg = "chain") {
  if (inherits(chain, "Chain")) return(chain)
  if (!is.character(chain) || length(chain) != 1L || is.na(chain) ||
      !file.exists(chain)) {
    .acv2g_validation_error(
      "{.arg {arg}} must be an existing chain file or an imported Chain object."
    )
  }
  rtracklayer::import.chain(chain)
}

.chain_provenance <- function(chain) {
  if (is.null(chain)) {
    return(list(path = NA_character_, sha256 = NA_character_))
  }
  if (is.character(chain) && length(chain) == 1L && !is.na(chain) &&
      file.exists(chain)) {
    path <- normalizePath(chain, mustWork = TRUE)
    return(list(
      path = path,
      sha256 = digest::digest(
        file = path, algo = "sha256", serialize = FALSE
      )
    ))
  }
  if (inherits(chain, "Chain")) {
    return(list(
      path = "<Chain>",
      sha256 = digest::digest(chain, algo = "sha256", serialize = TRUE)
    ))
  }
  list(path = NA_character_, sha256 = NA_character_)
}

#' Validate the ACV2G genome-build contract
#'
#' @param analysis_build Assembly of the contact source and structural target
#'   coordinate system.
#' @param input_builds Named character vector or list giving an explicit build
#'   for each genomic input role. The `contact` role is required.
#' @param chains Named list of chain files for roles whose build differs from
#'   `analysis_build`. Structural chains map each input build to the analysis
#'   build. The activity chain instead maps analysis-build query anchors to the
#'   activity build.
#'
#' @return A data frame describing required identity or liftover operations.
#' @export
validate_v2g_assemblies <- function(analysis_build, input_builds, chains = list()) {
  analysis_build <- .assert_build(analysis_build, "analysis_build")
  if (is.list(input_builds)) input_builds <- unlist(input_builds, use.names = TRUE)
  if (!is.character(input_builds) || is.null(names(input_builds)) ||
      any(!nzchar(names(input_builds)))) {
    .acv2g_validation_error("{.arg input_builds} must be a named set of assemblies.")
  }
  if (!"contact" %in% names(input_builds)) {
    .acv2g_validation_error("{.arg input_builds} must include the contact build.")
  }
  if (anyNA(input_builds) || any(!nzchar(input_builds))) {
    .acv2g_validation_error("Every genomic input must have an explicit assembly.")
  }
  if (!identical(unname(input_builds[["contact"]]), analysis_build)) {
    .acv2g_validation_error(
      "Contact data define {.arg analysis_build} and are never lifted; contact build {.val {input_builds[['contact']]}} does not equal {.val {analysis_build}}."
    )
  }
  roles <- names(input_builds)
  needs_liftover <- unname(input_builds) != analysis_build
  missing_chain <- roles[needs_liftover & !roles %in% names(chains)]
  if (length(missing_chain)) {
    .acv2g_validation_error(
      "Input role{?s} {.val {missing_chain}} require explicit chain file{?s}."
    )
  }
  source_build <- unname(input_builds)
  target_build <- rep(analysis_build, length(roles))
  activity_query <- roles == "activity" & needs_liftover
  source_build[activity_query] <- analysis_build
  target_build[activity_query] <- unname(input_builds[activity_query])
  data.frame(
    input_role = roles,
    source_build = source_build,
    target_build = target_build,
    operation = ifelse(needs_liftover, "liftover", "identity"),
    chain = vapply(roles, function(role) {
      value <- chains[[role]]
      if (is.null(value)) NA_character_ else if (is.character(value)) value else "<Chain>"
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}

#' Lift ACV2G point or interval coordinates with complete QC
#'
#' Only unique one-to-one mappings are accepted. Structural failures stop the
#' run by default; evidence-query failures are retained with missing target
#' coordinates and a warning.
#'
#' @param x Input data frame.
#' @param chain Forward chain file or imported Chain object.
#' @param source_build,target_build Explicit source and target assemblies.
#' @param record_type Either `point` or `interval`.
#' @param role Either `structural` or `evidence_query`.
#' @param chrom_col,pos1_col,start0_col,end0_col,id_col Column mappings.
#' @param strand_col Optional strand column. When supplied, strand is carried
#'   through the chain and updated for mappings across inversions.
#' @param failure_policy Structural failure behavior: `error` or `fail_locus`.
#' @param reciprocal_policy `none` by default or `require` for exact round trip.
#' @param reverse_chain Required when `reciprocal_policy = "require"`.
#'
#' @return A list with mapped `data` and `qc` tables.
#' @export
liftover_v2g_coordinates <- function(
    x,
    chain,
    source_build,
    target_build,
    record_type = c("point", "interval"),
    role = c("structural", "evidence_query"),
    chrom_col = "chrom",
    pos1_col = "pos1",
    start0_col = "start0",
    end0_col = "end0",
    id_col = NULL,
    strand_col = NULL,
    failure_policy = c("error", "fail_locus"),
    reciprocal_policy = c("none", "require"),
    reverse_chain = NULL) {
  record_type <- match.arg(record_type)
  role <- match.arg(role)
  failure_policy <- match.arg(failure_policy)
  reciprocal_policy <- match.arg(reciprocal_policy)
  source_build <- .assert_build(source_build, "source_build")
  target_build <- .assert_build(target_build, "target_build")
  x <- .as_plain_data_frame(x, "x")
  required <- if (record_type == "point") c(chrom_col, pos1_col) else
    c(chrom_col, start0_col, end0_col)
  .require_columns(x, required, "x")
  source_strand <- rep("*", nrow(x))
  if (!is.null(strand_col)) {
    .require_columns(x, strand_col, "x")
    source_strand <- as.character(x[[strand_col]])
    if (anyNA(source_strand) || any(!source_strand %in% c("+", "-", "*"))) {
      .acv2g_validation_error(
        "Liftover strand values must be '+', '-', or '*'."
      )
    }
  }
  source_id <- if (is.null(id_col)) paste0("record_", seq_len(nrow(x))) else {
    .require_columns(x, id_col, "x")
    as.character(x[[id_col]])
  }
  if (anyNA(source_id) || any(!nzchar(source_id)) || anyDuplicated(source_id)) {
    .acv2g_validation_error("Liftover record IDs must be unique and non-missing.")
  }

  if (record_type == "point") {
    pos1 <- suppressWarnings(as.integer(x[[pos1_col]]))
    if (anyNA(pos1) || any(pos1 < 1L)) {
      .acv2g_validation_error("Point coordinates must be positive 1-based integers.")
    }
    source_start1 <- pos1
    source_end1 <- pos1
  } else {
    start0 <- suppressWarnings(as.integer(x[[start0_col]]))
    end0 <- suppressWarnings(as.integer(x[[end0_col]]))
    if (anyNA(start0) || anyNA(end0) || any(start0 < 0L) || any(end0 <= start0)) {
      .acv2g_validation_error("Intervals must be valid 0-based half-open coordinates.")
    }
    source_start1 <- start0 + 1L
    source_end1 <- end0
  }
  chain_provenance <- .chain_provenance(NULL)
  reverse_chain_provenance <- .chain_provenance(NULL)

  if (identical(source_build, target_build)) {
    status <- rep("identity", nrow(x))
    unique_map <- rep(TRUE, nrow(x))
    mapping_count <- rep(1L, nrow(x))
    target_chrom <- as.character(x[[chrom_col]])
    target_start1 <- as.integer(source_start1)
    target_end1 <- as.integer(source_end1)
    target_strand <- source_strand
  } else {
    source_gr <- GenomicRanges::GRanges(
      seqnames = as.character(x[[chrom_col]]),
      ranges = IRanges::IRanges(start = source_start1, end = source_end1),
      strand = source_strand,
      source_index = seq_len(nrow(x))
    )
    chain_provenance <- .chain_provenance(chain)
    chain <- .read_chain(chain)
    lifted <- rtracklayer::liftOver(source_gr, chain)
    status <- vapply(lifted, function(gr) {
      if (length(gr) == 0L) "unmapped" else if (length(gr) > 1L) "split_or_non_unique" else "mapped"
    }, character(1L))
    mapping_count <- lengths(lifted)
    target_chrom <- rep(NA_character_, nrow(x))
    target_start1 <- target_end1 <- rep(NA_integer_, nrow(x))
    target_strand <- rep(NA_character_, nrow(x))
    unique_map <- mapping_count == 1L
    for (i in which(unique_map)) {
      mapped_df <- as.data.frame(lifted[[i]])
      target_chrom[[i]] <- as.character(mapped_df$seqnames[[1L]])
      target_start1[[i]] <- as.integer(mapped_df$start[[1L]])
      target_end1[[i]] <- as.integer(mapped_df$end[[1L]])
      target_strand[[i]] <- as.character(mapped_df$strand[[1L]])
      if (record_type == "point" && target_start1[[i]] != target_end1[[i]]) {
        status[[i]] <- "non_unit_point"
        unique_map[[i]] <- FALSE
      }
    }
  }

  reciprocal_ok <- rep(NA, nrow(x))
  reciprocal_status <- if (identical(source_build, target_build)) {
    rep("not_applicable", nrow(x))
  } else {
    rep("not_checked", nrow(x))
  }
  if (reciprocal_policy == "require" &&
      identical(source_build, target_build)) {
    reciprocal_ok[] <- TRUE
    reciprocal_status[] <- "not_applicable"
  } else if (reciprocal_policy == "require") {
    if (is.null(reverse_chain)) {
      .acv2g_validation_error(
        "{.arg reverse_chain} is required when reciprocal liftover is required."
      )
    }
    reverse_chain_provenance <- .chain_provenance(reverse_chain)
    reverse_chain <- .read_chain(reverse_chain, "reverse_chain")
    reciprocal_status[!unique_map] <- "not_testable"
    for (i in which(unique_map)) {
      target_gr <- GenomicRanges::GRanges(
        seqnames = target_chrom[[i]],
        ranges = IRanges::IRanges(start = target_start1[[i]], end = target_end1[[i]]),
        strand = target_strand[[i]]
      )
      back <- rtracklayer::liftOver(target_gr, reverse_chain)[[1L]]
      back_df <- as.data.frame(back)
      reciprocal_ok[[i]] <- nrow(back_df) == 1L &&
        as.character(back_df$seqnames[[1L]]) == as.character(x[[chrom_col]][[i]]) &&
        as.integer(back_df$start[[1L]]) == source_start1[[i]] &&
        as.integer(back_df$end[[1L]]) == source_end1[[i]] &&
        (is.null(strand_col) ||
           as.character(back_df$strand[[1L]]) == source_strand[[i]])
      if (!isTRUE(reciprocal_ok[[i]])) {
        status[[i]] <- "reciprocal_failed"
        unique_map[[i]] <- FALSE
        reciprocal_status[[i]] <- "failed"
      } else {
        reciprocal_status[[i]] <- "passed"
      }
    }
  }

  mapped <- x
  mapped[[paste0(chrom_col, "_source")]] <- as.character(x[[chrom_col]])
  mapped[[chrom_col]] <- target_chrom
  if (record_type == "point") {
    mapped[[paste0(pos1_col, "_source")]] <- x[[pos1_col]]
    mapped[[pos1_col]] <- ifelse(unique_map, target_start1, NA_integer_)
  } else {
    mapped[[paste0(start0_col, "_source")]] <- x[[start0_col]]
    mapped[[paste0(end0_col, "_source")]] <- x[[end0_col]]
    mapped[[start0_col]] <- ifelse(unique_map, target_start1 - 1L, NA_integer_)
    mapped[[end0_col]] <- ifelse(unique_map, target_end1, NA_integer_)
  }
  if (!is.null(strand_col)) {
    mapped[[paste0(strand_col, "_source")]] <- source_strand
    mapped[[strand_col]] <- ifelse(unique_map, target_strand, NA_character_)
  }
  mapped$liftover_status <- status
  mapped$source_build <- source_build
  mapped$target_build <- target_build

  qc <- data.frame(
    schema_version = .acv2g_schema_version,
    input_role = role,
    source_record_id = source_id,
    source_build = source_build,
    target_build = target_build,
    record_type = record_type,
    source_chrom = as.character(x[[chrom_col]]),
    source_start1 = source_start1,
    source_end1 = source_end1,
    target_chrom = target_chrom,
    target_start1 = target_start1,
    target_end1 = target_end1,
    source_strand = source_strand,
    target_strand = target_strand,
    mapping_count = mapping_count,
    chain_path = chain_provenance$path,
    chain_sha256 = chain_provenance$sha256,
    chain_direction = paste0(source_build, "->", target_build),
    reciprocal_policy = reciprocal_policy,
    reciprocal_ok = reciprocal_ok,
    reciprocal_status = reciprocal_status,
    reverse_chain_path = reverse_chain_provenance$path,
    reverse_chain_sha256 = reverse_chain_provenance$sha256,
    reverse_chain_direction = if (reciprocal_policy == "require" &&
                                  !identical(source_build, target_build)) {
      paste0(target_build, "->", source_build)
    } else {
      NA_character_
    },
    status = status,
    stringsAsFactors = FALSE
  )
  failed <- !unique_map
  if (any(failed)) {
    detail <- paste(source_id[failed], collapse = ", ")
    if (role == "structural" && failure_policy == "error") {
      .acv2g_validation_error(
        "Structural liftover failed unique one-to-one validation for: {detail}.",
        data = qc[failed, , drop = FALSE]
      )
    }
    .acv2g_warn(
      "Liftover retained {sum(failed)} failed {.val {role}} record{?s} with unavailable target coordinates."
    )
  }
  structure(list(data = mapped, qc = qc), class = "v2g_liftover")
}
