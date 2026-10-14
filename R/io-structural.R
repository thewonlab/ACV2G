# Structural input readers -------------------------------------------------

.detect_structural_format <- function(path, allowed) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    .acv2g_validation_error("Input path must be one non-empty string.")
  }
  lower <- tolower(path)
  lower <- sub("\\.(gz|bgz)$", "", lower)
  extension <- tools::file_ext(lower)
  resolved <- switch(
    extension,
    bed = "bed",
    vcf = "vcf",
    gtf = "gtf",
    gff = "gff",
    gff3 = "gff",
    csv = "csv",
    tsv = "tsv",
    txt = "tsv",
    NA_character_
  )
  if (is.na(resolved) || !resolved %in% allowed) {
    .acv2g_validation_error(
      paste0(
        "Could not infer a supported format for ", path,
        ". Supply format explicitly; accepted values are ",
        paste(allowed, collapse = ", "), "."
      )
    )
  }
  resolved
}

.resolve_input_format <- function(x, format, allowed) {
  if (is.data.frame(x) || inherits(x, "DataFrame")) {
    if (!is.null(format) && !identical(format, "table")) {
      .acv2g_validation_error(
        "In-memory data frames require {.code format = 'table'} or no format."
      )
    }
    return("table")
  }
  if (is.null(format)) {
    return(.detect_structural_format(x, allowed))
  }
  .assert_scalar_choice(format, c(allowed, "table"), arg = "format")
}

.read_delimited_table <- function(path, format) {
  if (!file.exists(path)) {
    .acv2g_validation_error(paste0("Input file does not exist: ", path, "."))
  }
  separator <- if (identical(format, "csv")) "," else "\t"
  as.data.frame(
    data.table::fread(
      path,
      sep = separator,
      header = TRUE,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = FALSE
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.apply_column_map <- function(x, column_map, object_name) {
  if (is.null(column_map)) {
    return(x)
  }
  if (is.list(column_map)) {
    column_map <- unlist(column_map, use.names = TRUE)
  }
  if (!is.character(column_map) || is.null(names(column_map)) ||
      any(!nzchar(names(column_map))) || any(!nzchar(column_map))) {
    .acv2g_validation_error(
      paste0(
        "column_map for ", object_name,
        " must be a named character vector from canonical names to input names."
      )
    )
  }
  absent <- setdiff(unname(column_map), names(x))
  if (length(absent)) {
    .acv2g_validation_error(
      paste0(
        "Mapped source column(s) absent from ", object_name, ": ",
        paste(absent, collapse = ", "), "."
      )
    )
  }
  duplicated_targets <- intersect(names(column_map), names(x))
  duplicated_targets <- setdiff(duplicated_targets, unname(column_map))
  if (length(duplicated_targets)) {
    .acv2g_validation_error(
      paste0(
        "Column mapping would overwrite existing canonical column(s): ",
        paste(duplicated_targets, collapse = ", "), "."
      )
    )
  }
  names(x)[match(unname(column_map), names(x))] <- names(column_map)
  x
}

.validate_declared_assembly_column <- function(x, assembly, object_name) {
  if ("assembly" %in% names(x)) {
    observed <- unique(as.character(x$assembly[!is.na(x$assembly)]))
    observed <- observed[nzchar(observed)]
    if (length(observed) && !identical(observed, assembly)) {
      .acv2g_validation_error(
        paste0(
          object_name, " assembly values (", paste(observed, collapse = ", "),
          ") do not match declared assembly ", assembly, "."
        )
      )
    }
  }
  invisible(x)
}

.read_bed4 <- function(path, object_name) {
  if (!file.exists(path)) {
    .acv2g_validation_error(paste0("Input file does not exist: ", path, "."))
  }
  connection <- if (grepl("\\.(gz|bgz)$", tolower(path))) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)
  leading <- readLines(connection, n = 1000L, warn = FALSE)
  is_metadata <- !nzchar(trimws(leading)) |
    grepl("^(#|track([[:space:]]|$)|browser([[:space:]]|$))", trimws(leading))
  first_data <- which(!is_metadata)[1L]
  skip <- if (length(first_data) && !is.na(first_data)) first_data - 1L else 0L
  x <- as.data.frame(
    data.table::fread(
      path,
      sep = "\t",
      header = FALSE,
      skip = skip,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = FALSE
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(x) < 3L) {
    .acv2g_validation_error(
      paste0("BED ", object_name, " must contain at least three columns.")
    )
  }
  x
}

#' Read fixed user-supplied loci
#'
#' Locus boundaries are always interpreted and returned as 0-based half-open
#' intervals. ACV2G never constructs, merges, or expands these boundaries.
#'
#' @param x A data frame or path to BED, TSV, or CSV input.
#' @param assembly Explicit genome assembly for the locus coordinates.
#' @param format Optional format override: `"bed"`, `"tsv"`, or `"csv"`.
#' @param column_map Optional named character vector mapping canonical column
#'   names (`locus_id`, `chrom`, `start0`, `end0`) to input columns.
#' @return A standardized locus data frame.
#' @export
read_v2g_loci <- function(x, assembly, format = NULL, column_map = NULL) {
  assembly <- .validate_assembly(assembly, "assembly")
  format <- .resolve_input_format(x, format, c("bed", "tsv", "csv"))

  if (identical(format, "bed")) {
    bed <- .read_bed4(x, "loci")
    if (ncol(bed) < 4L) {
      .acv2g_validation_error(
        "Locus BED requires column 4 as a unique, non-missing locus_id."
      )
    }
    loci <- bed[c(4L, 1L, 2L, 3L)]
    names(loci) <- c("locus_id", "chrom", "start0", "end0")
    input_convention <- "bed_0_based_half_open"
  } else {
    loci <- if (identical(format, "table")) {
      .as_plain_data_frame(x, "x")
    } else {
      .read_delimited_table(x, format)
    }
    loci <- .apply_column_map(loci, column_map, "loci")
    .require_columns(loci, c("locus_id", "chrom", "start0", "end0"), "loci")
    .validate_declared_assembly_column(loci, assembly, "loci")
    input_convention <- "declared_0_based_half_open"
  }

  loci$locus_id <- as.character(loci$locus_id)
  invalid_id <- is.na(loci$locus_id) | !nzchar(trimws(loci$locus_id))
  if (any(invalid_id)) {
    .acv2g_validation_error(
      paste0(
        "locus_id contains ", sum(invalid_id), " missing or empty value(s)."
      )
    )
  }
  loci$locus_id <- trimws(loci$locus_id)
  if (anyDuplicated(loci$locus_id)) {
    duplicates <- unique(loci$locus_id[duplicated(loci$locus_id)])
    .acv2g_validation_error(
      paste0(
        "locus_id must be globally unique; duplicate(s): ",
        paste(head(duplicates, 10L), collapse = ", "), "."
      )
    )
  }
  loci$chrom <- .validate_chrom(loci$chrom)
  intervals <- .validate_half_open(loci$start0, loci$end0, "loci")

  out <- data.frame(
    locus_id = loci$locus_id,
    locus_assembly = assembly,
    locus_chrom = loci$chrom,
    locus_start0 = intervals$start0,
    locus_end0 = intervals$end0,
    locus_input_convention = input_convention,
    stringsAsFactors = FALSE
  )
  class(out) <- c("v2g_loci", class(out))
  out
}

.standardize_variant_records <- function(records, assembly, input_convention) {
  .require_columns(records, c("chrom", "pos1"), "variants")
  .validate_declared_assembly_column(records, assembly, "variants")
  records$chrom <- .validate_chrom(records$chrom)
  records$pos1 <- .as_integer_coordinate(records$pos1, "variants$pos1", 1)

  optional <- c("variant_id", "ref", "alt", "rsid", "submitted_record_id")
  for (column in optional) {
    if (!column %in% names(records)) {
      records[[column]] <- NA_character_
    }
    records[[column]] <- as.character(records[[column]])
    records[[column]][!is.na(records[[column]]) &
                        !nzchar(records[[column]])] <- NA_character_
  }

  n <- nrow(records)
  records$record_id <- sprintf("record_%08d", seq_len(n))
  missing_submitted <- is.na(records$submitted_record_id)
  records$submitted_record_id[missing_submitted] <- records$record_id[missing_submitted]
  records$assembly <- assembly
  records$anchor0 <- .pos1_to_anchor0(records$pos1)
  records$variant_key <- paste(assembly, records$chrom, records$pos1, sep = ":")
  records$variant_anchor_id <- .stable_row_id(
    records$variant_key,
    prefix = "anchor"
  )
  has_alleles <- !is.na(records$ref) & !is.na(records$alt)
  records$allele_key <- NA_character_
  records$allele_key[has_alleles] <- paste(
    records$variant_key[has_alleles],
    records$ref[has_alleles],
    records$alt[has_alleles],
    sep = ":"
  )
  records$variant_input_convention <- input_convention

  records <- records[c(
    "record_id", "submitted_record_id", "variant_id", "variant_anchor_id",
    "variant_key", "allele_key", "assembly", "chrom", "pos1", "anchor0",
    "ref", "alt", "rsid", "variant_input_convention"
  )]

  anchor_columns <- c(
    "variant_anchor_id", "variant_key", "assembly", "chrom", "pos1", "anchor0"
  )
  anchors <- unique(records[anchor_columns])
  anchors <- anchors[order(anchors$assembly, anchors$chrom, anchors$pos1), , drop = FALSE]
  rownames(anchors) <- NULL
  multiplicity <- table(records$variant_anchor_id)
  anchors$anchor_multiplicity <- as.integer(multiplicity[anchors$variant_anchor_id])

  record_map <- records[c(
    "record_id", "submitted_record_id", "variant_id", "variant_anchor_id",
    "variant_key", "allele_key", "assembly", "chrom", "pos1", "ref", "alt",
    "rsid"
  )]
  record_map$anchor_multiplicity <- as.integer(
    multiplicity[record_map$variant_anchor_id]
  )
  qc <- anchors[c(
    "variant_anchor_id", "variant_key", "assembly", "chrom", "pos1",
    "anchor_multiplicity"
  )]
  qc$variant_input_status <- "valid_anchor"

  structure(
    list(
      records = records,
      anchors = anchors,
      record_map = record_map,
      qc = qc,
      assembly = assembly
    ),
    class = "v2g_variants"
  )
}

.read_vcf_records <- function(path, assembly) {
  if (!file.exists(path)) {
    .acv2g_validation_error(paste0("Input file does not exist: ", path, "."))
  }
  vcf <- VariantAnnotation::readVcf(path, genome = assembly)
  rr <- GenomicRanges::granges(vcf)
  rr_table <- as.data.frame(rr, stringsAsFactors = FALSE)
  refs <- as.character(VariantAnnotation::ref(vcf))
  alts <- VariantAnnotation::alt(vcf)
  ids <- names(rr)
  if (is.null(ids)) {
    ids <- rep(NA_character_, length(rr))
  }
  ids[ids == "." | !nzchar(ids)] <- NA_character_

  pieces <- lapply(seq_along(rr), function(i) {
    alt_i <- as.character(alts[[i]])
    if (!length(alt_i)) {
      alt_i <- NA_character_
    }
    data.frame(
      submitted_record_id = paste0("vcf_row_", i, "_alt_", seq_along(alt_i)),
      variant_id = rep(ids[[i]], length(alt_i)),
      chrom = rep(as.character(rr_table$seqnames[[i]]), length(alt_i)),
      pos1 = rep(rr_table$start[[i]], length(alt_i)),
      ref = rep(refs[[i]], length(alt_i)),
      alt = alt_i,
      rsid = rep(if (!is.na(ids[[i]]) && grepl("^rs[0-9]+$", ids[[i]])) {
        ids[[i]]
      } else {
        NA_character_
      }, length(alt_i)),
      stringsAsFactors = FALSE
    )
  })
  if (!length(pieces)) {
    return(data.frame(
      submitted_record_id = character(), variant_id = character(),
      chrom = character(), pos1 = integer(), ref = character(),
      alt = character(), rsid = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, pieces)
}

#' Read variants and deduplicate coordinate-level scoring anchors
#'
#' The score is coordinate-level, so repeated IDs and multiple alleles at the
#' same assembly-qualified chromosome/position share one anchor. The returned
#' object retains a separate many-to-one submitted-record map.
#'
#' @param x A data frame or path to VCF, BED, TSV, or CSV input.
#' @param assembly Explicit genome assembly for the variant coordinates.
#' @param format Optional input format override.
#' @param column_map Optional named mapping of canonical table columns to input
#'   columns.
#' @param bed_anchor_column Optional one-based anchor-position column in a BED
#'   file. Without it, every BED interval must be exactly one base wide.
#' @return A `v2g_variants` object containing `records`, unique `anchors`,
#'   `record_map`, and `qc` tables.
#' @export
read_v2g_variants <- function(
    x,
    assembly,
    format = NULL,
    column_map = NULL,
    bed_anchor_column = NULL) {
  assembly <- .validate_assembly(assembly, "assembly")
  format <- .resolve_input_format(x, format, c("vcf", "bed", "tsv", "csv"))

  if (identical(format, "vcf")) {
    records <- .read_vcf_records(x, assembly)
    input_convention <- "vcf_1_based"
  } else if (identical(format, "bed")) {
    bed <- .read_bed4(x, "variants")
    intervals <- .validate_half_open(bed[[2L]], bed[[3L]], "variant BED")
    if (is.null(bed_anchor_column)) {
      if (any(intervals$end0 - intervals$start0 != 1)) {
        .acv2g_validation_error(
          paste0(
            "Variant BED intervals must be one base wide unless ",
            "bed_anchor_column explicitly supplies a 1-based anchor."
          )
        )
      }
      pos1 <- intervals$start0 + 1
    } else {
      if (is.character(bed_anchor_column)) {
        .acv2g_validation_error(
          paste0(
            "Headerless BED input requires bed_anchor_column as a ",
            "one-based column number."
          )
        )
      }
      bed_anchor_column <- .assert_scalar_number(
        bed_anchor_column,
        "bed_anchor_column",
        lower = 1,
        upper = ncol(bed),
        integer = TRUE
      )
      pos1 <- bed[[bed_anchor_column]]
      explicit_anchor0 <- .pos1_to_anchor0(pos1)
      outside_interval <- explicit_anchor0 < intervals$start0 |
        explicit_anchor0 >= intervals$end0
      if (any(outside_interval)) {
        .acv2g_validation_error(
          paste0(
            "bed_anchor_column places ", sum(outside_interval),
            " anchor(s) outside their submitted BED interval."
          )
        )
      }
    }
    records <- data.frame(
      chrom = bed[[1L]],
      pos1 = pos1,
      variant_id = if (ncol(bed) >= 4L) bed[[4L]] else NA_character_,
      stringsAsFactors = FALSE
    )
    input_convention <- "bed_0_based_half_open"
  } else {
    records <- if (identical(format, "table")) {
      .as_plain_data_frame(x, "x")
    } else {
      .read_delimited_table(x, format)
    }
    records <- .apply_column_map(records, column_map, "variants")
    input_convention <- "declared_1_based_position"
  }

  .standardize_variant_records(records, assembly, input_convention)
}

.as_variant_anchors <- function(x) {
  if (inherits(x, "v2g_variants")) {
    x <- x$anchors
  } else {
    x <- .as_plain_data_frame(x, "variants")
  }
  .require_columns(
    x,
    c("variant_anchor_id", "variant_key", "assembly", "chrom", "pos1"),
    "variants"
  )
  if (!"anchor0" %in% names(x)) {
    x$anchor0 <- .pos1_to_anchor0(x$pos1)
  }
  x$chrom <- .validate_chrom(x$chrom, "variants$chrom")
  x$pos1 <- .as_integer_coordinate(x$pos1, "variants$pos1", 1)
  expected_anchor0 <- .pos1_to_anchor0(x$pos1)
  supplied_anchor0 <- .as_integer_coordinate(x$anchor0, "variants$anchor0", 0)
  if (any(supplied_anchor0 != expected_anchor0)) {
    .acv2g_validation_error(
      "Standardized variant anchor0 values disagree with their 1-based positions."
    )
  }
  x$anchor0 <- supplied_anchor0
  invalid_assembly <- is.na(x$assembly) |
    !nzchar(trimws(as.character(x$assembly)))
  if (any(invalid_assembly)) {
    .acv2g_validation_error(
      "Standardized variants require an explicit assembly on every row."
    )
  }
  expected_key <- paste(x$assembly, x$chrom, x$pos1, sep = ":")
  key <- as.character(x$variant_key)
  if (any(is.na(key) | key != expected_key)) {
    .acv2g_validation_error(
      "Standardized variant_key values disagree with assembly/chrom/pos1."
    )
  }
  anchor_id <- as.character(x$variant_anchor_id)
  invalid_anchor_id <- is.na(anchor_id) | !nzchar(anchor_id)
  if (any(invalid_anchor_id) || anyDuplicated(anchor_id) || anyDuplicated(key)) {
    .acv2g_validation_error(
      "Standardized variant anchor IDs and keys must each be unique."
    )
  }
  if (!"anchor_multiplicity" %in% names(x)) x$anchor_multiplicity <- 1L
  x
}
