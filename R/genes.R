# Gene and TSS annotation --------------------------------------------------

.nonempty_character <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & nzchar(trimws(x))]
}

.collapse_values <- function(x) {
  x <- sort(unique(.nonempty_character(x)))
  if (length(x)) paste(x, collapse = ";") else NA_character_
}

.normalize_biotype_token <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[[:space:]-]+", "_", x)
}

.biotype_values <- function(x) {
  x <- .nonempty_character(x)
  if (!length(x)) return(character())
  unique(trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)))
}

.normalize_gene_id <- function(x, policy) {
  policy <- .assert_scalar_choice(
    policy,
    c("preserve", "strip_terminal_numeric_version"),
    arg = "id_version_policy"
  )
  raw <- as.character(x)
  bad <- is.na(raw) | !nzchar(trimws(raw))
  if (any(bad)) {
    .acv2g_validation_error(
      paste0("gene_id contains ", sum(bad), " missing or empty value(s).")
    )
  }
  raw <- trimws(raw)
  if (identical(policy, "strip_terminal_numeric_version")) {
    has_version <- grepl("\\.[0-9]+$", raw)
    version <- ifelse(has_version, sub("^.*\\.([0-9]+)$", "\\1", raw), NA_character_)
    normalized <- sub("\\.[0-9]+$", "", raw)
  } else {
    version <- rep(NA_character_, length(raw))
    normalized <- raw
  }
  data.frame(
    gene_id_raw = raw,
    gene_id_normalized = normalized,
    gene_id_version = version,
    id_version_policy = policy,
    stringsAsFactors = FALSE
  )
}

.parse_logical_field <- function(x, field = "is_protein_coding") {
  if (is.logical(x)) {
    return(x)
  }
  value <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[value %in% c("false", "f", "0", "no", "n")] <- FALSE
  invalid <- !is.na(value) & nzchar(value) & is.na(out)
  if (any(invalid)) {
    .acv2g_validation_error(
      paste0(
        field, " contains unrecognized logical value(s): ",
        paste(unique(value[invalid]), collapse = ", "), "."
      )
    )
  }
  out
}

.resolve_annotation_field <- function(x, explicit, candidates, role) {
  if (!is.null(explicit)) {
    if (!is.character(explicit) || length(explicit) != 1L ||
        is.na(explicit) || !nzchar(explicit)) {
      .acv2g_validation_error(paste0(role, " must name exactly one field."))
    }
    if (!explicit %in% names(x)) {
      .acv2g_validation_error(
        paste0("Configured ", role, " field ", explicit, " is absent.")
      )
    }
    return(explicit)
  }
  present <- candidates[candidates %in% names(x)]
  if (length(present)) present[[1L]] else NA_character_
}

.scalar_metadata <- function(x) {
  if (inherits(x, "List") || is.list(x)) {
    return(vapply(x, function(value) {
      value <- as.character(value)
      if (!length(value)) NA_character_ else paste(value, collapse = ",")
    }, character(1L)))
  }
  as.character(x)
}

.annotation_column <- function(metadata, field, n) {
  if (is.na(field) || !field %in% names(metadata)) {
    return(rep(NA_character_, n))
  }
  .scalar_metadata(metadata[[field]])
}

.import_gene_annotation <- function(
    path,
    format,
    gene_biotype_field,
    transcript_biotype_field) {
  if (!file.exists(path)) {
    .acv2g_validation_error(paste0("Annotation file does not exist: ", path, "."))
  }
  gr <- rtracklayer::import(path, format = if (format == "gff") "gff3" else "gtf")
  range_table <- as.data.frame(gr, stringsAsFactors = FALSE)
  metadata <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
  n <- length(gr)
  type <- if ("type" %in% names(metadata)) {
    as.character(metadata$type)
  } else {
    rep(NA_character_, n)
  }

  raw_gene_id <- if ("gene_id" %in% names(metadata)) {
    .scalar_metadata(metadata$gene_id)
  } else {
    rep(NA_character_, n)
  }
  gene_rows <- type == "gene"
  if ("ID" %in% names(metadata)) {
    missing_gene <- (is.na(raw_gene_id) | !nzchar(raw_gene_id)) & gene_rows
    raw_gene_id[missing_gene] <- .scalar_metadata(metadata$ID)[missing_gene]
  }
  if ("Parent" %in% names(metadata)) {
    missing_gene <- (is.na(raw_gene_id) | !nzchar(raw_gene_id)) & !gene_rows
    raw_gene_id[missing_gene] <- .scalar_metadata(metadata$Parent)[missing_gene]
  }

  transcript_id <- if ("transcript_id" %in% names(metadata)) {
    .scalar_metadata(metadata$transcript_id)
  } else if ("ID" %in% names(metadata)) {
    .scalar_metadata(metadata$ID)
  } else {
    rep(NA_character_, n)
  }

  gene_biotype_field <- .resolve_annotation_field(
    metadata,
    gene_biotype_field,
    c("gene_biotype", "gene_type"),
    "gene_biotype_field"
  )
  transcript_biotype_field <- .resolve_annotation_field(
    metadata,
    transcript_biotype_field,
    c("transcript_biotype", "transcript_type"),
    "transcript_biotype_field"
  )
  gene_symbol_field <- .resolve_annotation_field(
    metadata,
    NULL,
    c("gene_name", "gene_symbol", "Name"),
    "gene_symbol_field"
  )

  base <- data.frame(
    gene_id = raw_gene_id,
    transcript_id = transcript_id,
    chrom = as.character(range_table$seqnames),
    feature_start1 = range_table$start,
    feature_end1 = range_table$end,
    strand = as.character(range_table$strand),
    feature_type = type,
    gene_symbol = .annotation_column(metadata, gene_symbol_field, n),
    gene_biotype = .annotation_column(metadata, gene_biotype_field, n),
    transcript_biotype = .annotation_column(metadata, transcript_biotype_field, n),
    stringsAsFactors = FALSE
  )

  transcript_rows <- base$feature_type %in% c("transcript", "mRNA") &
    !is.na(base$gene_id) & nzchar(base$gene_id)
  if (any(transcript_rows)) {
    tx <- base[transcript_rows, , drop = FALSE]
  } else {
    exon_rows <- base$feature_type == "exon" &
      !is.na(base$gene_id) & nzchar(base$gene_id) &
      !is.na(base$transcript_id) & nzchar(base$transcript_id)
    if (any(exon_rows)) {
      exon <- base[exon_rows, , drop = FALSE]
      split_exon <- split(exon, paste(exon$gene_id, exon$transcript_id, sep = "\u001f"))
      tx <- do.call(rbind, lapply(split_exon, function(group) {
        strands <- unique(group$strand)
        chroms <- unique(group$chrom)
        if (length(strands) != 1L || length(chroms) != 1L) {
          .acv2g_validation_error(
            paste0(
              "Transcript ", group$transcript_id[[1L]],
              " spans multiple chromosomes or strands."
            )
          )
        }
        data.frame(
          gene_id = group$gene_id[[1L]],
          transcript_id = group$transcript_id[[1L]],
          chrom = chroms,
          feature_start1 = min(group$feature_start1),
          feature_end1 = max(group$feature_end1),
          strand = strands,
          feature_type = "transcript_from_exons",
          gene_symbol = .collapse_values(group$gene_symbol),
          gene_biotype = .collapse_values(group$gene_biotype),
          transcript_biotype = .collapse_values(group$transcript_biotype),
          stringsAsFactors = FALSE
        )
      }))
    } else {
      usable_gene <- gene_rows & !is.na(base$gene_id) & nzchar(base$gene_id)
      if (!any(usable_gene)) {
        .acv2g_validation_error(
          paste0(
            "GTF/GFF contains no usable transcript, exon, or gene features ",
            "with explicit gene identifiers."
          )
        )
      }
      tx <- base[usable_gene, , drop = FALSE]
      tx$transcript_id <- NA_character_
      tx$feature_type <- "gene_as_tss"
    }
  }

  # A gene-level field is semantically gene-level even when repeated on
  # transcript/exon records. Aggregate it over every record carrying gene_id.
  gene_metadata <- base[!is.na(base$gene_id) & nzchar(base$gene_id), , drop = FALSE]
  gene_meta_split <- split(gene_metadata, gene_metadata$gene_id)
  gene_meta <- do.call(rbind, lapply(gene_meta_split, function(group) {
    data.frame(
      gene_id = group$gene_id[[1L]],
      gene_symbol_meta = .collapse_values(group$gene_symbol),
      gene_biotype_meta = .collapse_values(group$gene_biotype),
      stringsAsFactors = FALSE
    )
  }))
  tx <- merge(tx, gene_meta, by = "gene_id", all.x = TRUE, sort = FALSE)
  tx$gene_symbol <- ifelse(
    is.na(tx$gene_symbol_meta), tx$gene_symbol, tx$gene_symbol_meta
  )
  tx$gene_biotype <- ifelse(
    is.na(tx$gene_biotype_meta), tx$gene_biotype, tx$gene_biotype_meta
  )
  tx$tss_pos1 <- ifelse(tx$strand == "-", tx$feature_end1, tx$feature_start1)
  tx$is_protein_coding <- NA
  tx$annotation_version <- NA_character_
  tx <- tx[c(
    "gene_id", "transcript_id", "chrom", "tss_pos1", "strand",
    "gene_symbol", "gene_biotype", "transcript_biotype",
    "is_protein_coding", "annotation_version"
  )]
  attr(tx, "gene_biotype_field") <- gene_biotype_field
  attr(tx, "transcript_biotype_field") <- transcript_biotype_field
  tx
}

.classify_raw_genes <- function(
    records,
    normalized_ids,
    protein_coding_values,
    protein_coding_gene_ids,
    gene_biotype_field,
    transcript_biotype_field) {
  pc_values <- unique(.normalize_biotype_token(protein_coding_values))
  if (!length(pc_values) || anyNA(pc_values) || any(!nzchar(pc_values))) {
    .acv2g_validation_error(
      "{.arg protein_coding_values} must contain explicit non-empty values."
    )
  }
  override_declared <- !is.null(protein_coding_gene_ids)
  override_ids <- if (override_declared) unique(protein_coding_gene_ids) else character()

  records$.gene_id_normalized <- normalized_ids
  groups <- split(records, records$gene_id)
  classified <- lapply(groups, function(group) {
    raw_id <- group$gene_id[[1L]]
    normalized_id <- unique(group$.gene_id_normalized)
    gene_values <- .biotype_values(group$gene_biotype)
    transcript_values <- .biotype_values(group$transcript_biotype)
    gene_classes <- unique(ifelse(
      .normalize_biotype_token(gene_values) %in% pc_values,
      "protein_coding", "non_protein_coding"
    ))
    transcript_classes <- unique(ifelse(
      .normalize_biotype_token(transcript_values) %in% pc_values,
      "protein_coding", "non_protein_coding"
    ))
    flags <- unique(group$is_protein_coding[!is.na(group$is_protein_coding)])
    if (length(gene_classes) > 1L) {
      .acv2g_validation_error(
        paste0("Gene ", raw_id, " has conflicting gene-level biotype values.")
      )
    }
    if (length(flags) > 1L) {
      .acv2g_validation_error(
        paste0("Gene ", raw_id, " has conflicting is_protein_coding values.")
      )
    }
    override_member <- override_declared && normalized_id %in% override_ids
    conflict_reasons <- character()

    if (length(gene_classes)) {
      normalized_class <- gene_classes[[1L]]
      classification_level <- "gene"
      biotype_source <- "gene_biotype"
      override_source <- "none"
      if (normalized_class == "non_protein_coding" &&
          "protein_coding" %in% transcript_classes) {
        conflict_reasons <- c(conflict_reasons, "gene_non_coding_transcript_protein_coding")
      }
      if (length(flags) && flags[[1L]] != (normalized_class == "protein_coding")) {
        conflict_reasons <- c(conflict_reasons, "gene_biotype_is_protein_coding_disagree")
      }
      if (override_member && normalized_class != "protein_coding") {
        conflict_reasons <- c(conflict_reasons, "gene_biotype_gene_id_override_disagree")
      }
    } else if (length(transcript_classes)) {
      normalized_class <- if ("protein_coding" %in% transcript_classes) {
        "protein_coding"
      } else {
        "non_protein_coding"
      }
      classification_level <- "transcript"
      biotype_source <- "transcript_fallback"
      override_source <- "none"
      if (length(flags) && flags[[1L]] != (normalized_class == "protein_coding")) {
        conflict_reasons <- c(conflict_reasons, "transcript_biotype_is_protein_coding_disagree")
      }
      if (override_member && normalized_class != "protein_coding") {
        conflict_reasons <- c(conflict_reasons, "transcript_biotype_gene_id_override_disagree")
      }
    } else if (length(flags)) {
      normalized_class <- if (flags[[1L]]) "protein_coding" else "non_protein_coding"
      classification_level <- "explicit_override"
      biotype_source <- "is_protein_coding"
      override_source <- "gene_tss_table:is_protein_coding"
      if (override_declared && flags[[1L]] != override_member) {
        conflict_reasons <- c(
          conflict_reasons,
          "is_protein_coding_gene_id_override_disagree"
        )
      }
    } else if (override_declared) {
      normalized_class <- if (override_member) "protein_coding" else "non_protein_coding"
      classification_level <- "explicit_override"
      biotype_source <- "protein_coding_gene_ids"
      override_source <- "protein_coding_gene_ids"
    } else {
      .acv2g_validation_error(
        paste0(
          "Gene ", raw_id,
          " has no usable gene/transcript biotype. Supply ",
          "protein_coding_gene_ids or is_protein_coding."
        )
      )
    }

    data.frame(
      gene_id_raw = raw_id,
      gene_id_normalized = normalized_id,
      gene_biotype_raw = .collapse_values(gene_values),
      transcript_biotype_raw = .collapse_values(transcript_values),
      biotype_normalized = normalized_class,
      classification_level = classification_level,
      biotype_source = biotype_source,
      gene_biotype_field = ifelse(is.na(gene_biotype_field), NA_character_, gene_biotype_field),
      transcript_biotype_field = ifelse(
        is.na(transcript_biotype_field), NA_character_, transcript_biotype_field
      ),
      protein_coding_values = paste(pc_values, collapse = ";"),
      biotype_conflict = length(conflict_reasons) > 0L,
      conflict_reason = if (length(conflict_reasons)) {
        paste(conflict_reasons, collapse = ";")
      } else {
        NA_character_
      },
      override_source = override_source,
      is_protein_coding_override_value = if (length(flags)) flags[[1L]] else NA,
      protein_coding_gene_ids_supplied = override_declared,
      protein_coding_gene_id_member = if (override_declared) override_member else NA,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, classified)
}

.validate_normalized_gene_collisions <- function(records, classification) {
  records <- merge(
    records,
    classification[c(
      "gene_id_raw", "gene_id_normalized", "biotype_normalized", "override_source"
    )],
    by.x = "gene_id",
    by.y = "gene_id_raw",
    all.x = TRUE,
    sort = FALSE
  )
  groups <- split(records, records$gene_id_normalized)
  conflict <- vapply(groups, function(group) {
    length(unique(group$chrom)) > 1L ||
      length(unique(group$strand)) > 1L ||
      length(unique(group$biotype_normalized)) > 1L ||
      length(unique(group$override_source)) > 1L
  }, logical(1L))
  if (any(conflict)) {
    bad <- names(conflict)[conflict]
    .acv2g_validation_error(
      paste0(
        "ID normalization produced incompatible gene collision(s): ",
        paste(head(bad, 10L), collapse = ", "),
        ". Chromosome, strand, biotype class, and override provenance must agree."
      )
    )
  }
  invisible(records)
}

#' Read and classify gene/TSS annotation
#'
#' Gene-level biotype is authoritative. Transcript biotype is used only when
#' gene-level classification is unavailable. Every unique supplied TSS for an
#' eligible stable gene is retained; symbols are optional annotations.
#'
#' @param x A standardized gene--TSS data frame or path to GTF, GFF/GFF3, TSV,
#'   or CSV input.
#' @param assembly Explicit genome assembly for the annotation.
#' @param format Optional format override.
#' @param column_map Optional canonical-to-input column mapping for tabular
#'   input.
#' @param gene_biotype_field,transcript_biotype_field Optional explicit field
#'   names. Defaults recognize `gene_biotype`/`gene_type` and
#'   `transcript_biotype`/`transcript_type`, respectively.
#' @param protein_coding_values Explicit values recognized as protein coding
#'   after whitespace/case/hyphen normalization.
#' @param protein_coding_gene_ids Optional complete protein-coding ID whitelist
#'   used only when gene and transcript classification are both absent.
#' @param id_version_policy Either preserve IDs or strip only terminal numeric
#'   versions such as `.12`.
#' @return A `v2g_genes` object containing eligible `tss`, all classified
#'   `genes`, and `annotation_qc`.
#' @export
read_v2g_genes <- function(
    x,
    assembly,
    format = NULL,
    column_map = NULL,
    gene_biotype_field = NULL,
    transcript_biotype_field = NULL,
    protein_coding_values = c("protein_coding", "protein-coding"),
    protein_coding_gene_ids = NULL,
    id_version_policy = c("preserve", "strip_terminal_numeric_version")) {
  assembly <- .validate_assembly(assembly, "assembly")
  id_version_policy <- match.arg(id_version_policy)
  format <- .resolve_input_format(x, format, c("gtf", "gff", "tsv", "csv"))

  if (format %in% c("gtf", "gff")) {
    records <- .import_gene_annotation(
      x, format, gene_biotype_field, transcript_biotype_field
    )
    resolved_gene_field <- attr(records, "gene_biotype_field")
    resolved_transcript_field <- attr(records, "transcript_biotype_field")
    input_convention <- paste0(format, "_1_based_closed")
  } else {
    records <- if (identical(format, "table")) {
      .as_plain_data_frame(x, "x")
    } else {
      .read_delimited_table(x, format)
    }
    records <- .apply_column_map(records, column_map, "genes")
    .require_columns(
      records,
      c("gene_id", "chrom", "tss_pos1", "strand"),
      "genes"
    )
    .validate_declared_assembly_column(records, assembly, "genes")
    resolved_gene_field <- .resolve_annotation_field(
      records,
      gene_biotype_field,
      c("gene_biotype", "gene_type"),
      "gene_biotype_field"
    )
    resolved_transcript_field <- .resolve_annotation_field(
      records,
      transcript_biotype_field,
      c("transcript_biotype", "transcript_type"),
      "transcript_biotype_field"
    )
    records$gene_biotype <- .annotation_column(
      records, resolved_gene_field, nrow(records)
    )
    records$transcript_biotype <- .annotation_column(
      records, resolved_transcript_field, nrow(records)
    )
    input_convention <- "declared_tss_1_based"
  }

  if (!nrow(records)) {
    .acv2g_validation_error("Gene/TSS annotation contains no records.")
  }

  for (column in c("transcript_id", "gene_symbol", "annotation_version")) {
    if (!column %in% names(records)) records[[column]] <- NA_character_
    records[[column]] <- as.character(records[[column]])
  }
  if (!"is_protein_coding" %in% names(records)) {
    records$is_protein_coding <- NA
  }
  records$is_protein_coding <- .parse_logical_field(records$is_protein_coding)
  records$gene_id <- as.character(records$gene_id)
  records$chrom <- .validate_chrom(records$chrom, "genes$chrom")
  records$tss_pos1 <- .as_integer_coordinate(records$tss_pos1, "genes$tss_pos1", 1)
  records$strand <- as.character(records$strand)
  invalid_strand <- is.na(records$strand) | !records$strand %in% c("+", "-")
  if (any(invalid_strand)) {
    .acv2g_validation_error(
      paste0(
        "Every TSS must have explicit '+' or '-' strand; ",
        sum(invalid_strand), " row(s) failed validation."
      )
    )
  }

  ids <- .normalize_gene_id(records$gene_id, id_version_policy)
  if (!is.null(protein_coding_gene_ids)) {
    override_ids <- .normalize_gene_id(protein_coding_gene_ids, id_version_policy)
    protein_coding_gene_ids <- unique(override_ids$gene_id_normalized)
  }
  classification <- .classify_raw_genes(
    records,
    ids$gene_id_normalized,
    protein_coding_values,
    protein_coding_gene_ids,
    resolved_gene_field,
    resolved_transcript_field
  )
  .validate_normalized_gene_collisions(records, classification)

  records$gene_id_raw <- ids$gene_id_raw
  records$gene_id_normalized <- ids$gene_id_normalized
  records$gene_id_version <- ids$gene_id_version
  records$assembly <- assembly
  records$tss_chrom <- records$chrom
  records$tss_pos1 <- as.numeric(records$tss_pos1)
  records$tss_anchor0 <- .pos1_to_anchor0(records$tss_pos1)
  records$tss_strand <- records$strand
  records$annotation_input_convention <- input_convention
  records <- merge(
    records,
    classification,
    by = c("gene_id_raw", "gene_id_normalized"),
    all.x = TRUE,
    sort = FALSE
  )

  normalized_groups <- split(records, records$gene_id_normalized)
  genes <- do.call(rbind, lapply(normalized_groups, function(group) {
    data.frame(
      gene_id = group$gene_id_normalized[[1L]],
      gene_id_normalized = group$gene_id_normalized[[1L]],
      gene_id_raw = sort(unique(group$gene_id_raw))[[1L]],
      gene_id_raw_values = .collapse_values(group$gene_id_raw),
      gene_id_version_values = .collapse_values(group$gene_id_version),
      id_version_policy = id_version_policy,
      gene_id_collision_count = length(unique(group$gene_id_raw)),
      assembly = assembly,
      gene_chrom = unique(group$tss_chrom)[[1L]],
      gene_strand = unique(group$tss_strand)[[1L]],
      gene_symbol = .collapse_values(group$gene_symbol),
      annotation_version = .collapse_values(group$annotation_version),
      gene_biotype_raw = .collapse_values(group$gene_biotype_raw),
      transcript_biotype_raw = .collapse_values(group$transcript_biotype_raw),
      biotype_normalized = unique(group$biotype_normalized)[[1L]],
      is_protein_coding = unique(group$biotype_normalized)[[1L]] == "protein_coding",
      classification_level = .collapse_values(group$classification_level),
      biotype_source = .collapse_values(group$biotype_source),
      biotype_conflict = any(group$biotype_conflict),
      conflict_reason = .collapse_values(group$conflict_reason),
      override_source = .collapse_values(group$override_source),
      stringsAsFactors = FALSE
    )
  }))
  rownames(genes) <- NULL

  eligible <- records$biotype_normalized == "protein_coding"
  eligible_records <- records[eligible, , drop = FALSE]
  if (nrow(eligible_records)) {
    tss_key <- paste(
      eligible_records$gene_id_normalized,
      eligible_records$tss_chrom,
      eligible_records$tss_pos1,
      eligible_records$tss_strand,
      sep = "\u001f"
    )
    tss_groups <- split(eligible_records, tss_key)
    tss <- do.call(rbind, lapply(tss_groups, function(group) {
      gene_id <- group$gene_id_normalized[[1L]]
      chrom <- group$tss_chrom[[1L]]
      pos1 <- group$tss_pos1[[1L]]
      strand <- group$tss_strand[[1L]]
      data.frame(
        tss_id = .stable_row_id(assembly, gene_id, chrom, pos1, strand, prefix = "tss"),
        gene_id = gene_id,
        gene_id_normalized = gene_id,
        gene_id_raw = sort(unique(group$gene_id_raw))[[1L]],
        gene_id_raw_values = .collapse_values(group$gene_id_raw),
        transcript_id = .first_or_na(group$transcript_id),
        transcript_ids = .collapse_values(group$transcript_id),
        assembly = assembly,
        tss_chrom = chrom,
        tss_pos1 = pos1,
        tss_anchor0 = pos1 - 1,
        tss_strand = strand,
        gene_symbol = .collapse_values(group$gene_symbol),
        annotation_version = .collapse_values(group$annotation_version),
        biotype_normalized = "protein_coding",
        annotation_input_convention = input_convention,
        exact_tss_record_count = nrow(group),
        stringsAsFactors = FALSE
      )
    }))
    rownames(tss) <- NULL
    tss <- tss[order(tss$gene_id_normalized, tss$tss_chrom, tss$tss_pos1), , drop = FALSE]
  } else {
    tss <- data.frame(
      tss_id = character(), gene_id = character(), gene_id_normalized = character(),
      gene_id_raw = character(), gene_id_raw_values = character(),
      transcript_id = character(), transcript_ids = character(), assembly = character(),
      tss_chrom = character(), tss_pos1 = numeric(), tss_anchor0 = numeric(),
      tss_strand = character(), gene_symbol = character(), annotation_version = character(),
      biotype_normalized = character(), annotation_input_convention = character(),
      exact_tss_record_count = integer(), stringsAsFactors = FALSE
    )
  }

  annotation_qc <- classification
  raw_id_map <- unique(ids[c("gene_id_raw", "gene_id_normalized")])
  collision_count <- vapply(
    split(raw_id_map$gene_id_raw, raw_id_map$gene_id_normalized),
    function(value) length(unique(value)),
    integer(1L)
  )
  annotation_qc$gene_id_collision_count <- unname(
    collision_count[annotation_qc$gene_id_normalized]
  )
  version_map <- split(ids$gene_id_version, ids$gene_id_raw)
  annotation_qc$gene_id_version <- vapply(
    annotation_qc$gene_id_raw,
    function(id) .collapse_values(version_map[[id]]),
    character(1L)
  )
  annotation_qc$id_version_policy <- id_version_policy
  annotation_qc$annotation_status <- ifelse(
    annotation_qc$biotype_conflict,
    "classified_with_conflict",
    ifelse(
      annotation_qc$biotype_normalized == "protein_coding",
      "eligible_protein_coding",
      "ineligible_non_protein_coding"
    )
  )

  structure(
    list(
      tss = tss,
      genes = genes,
      annotation_qc = annotation_qc,
      assembly = assembly,
      id_version_policy = id_version_policy
    ),
    class = "v2g_genes"
  )
}

.first_or_na <- function(x) {
  x <- sort(unique(.nonempty_character(x)))
  if (length(x)) x[[1L]] else NA_character_
}

.as_gene_tss <- function(x) {
  if (inherits(x, "v2g_genes")) {
    x <- x$tss
  } else {
    x <- .as_plain_data_frame(x, "genes")
  }
  .require_columns(
    x,
    c("gene_id_normalized", "assembly", "tss_chrom", "tss_pos1", "tss_strand"),
    "genes"
  )
  if (!"gene_id" %in% names(x)) x$gene_id <- x$gene_id_normalized
  if (!"tss_anchor0" %in% names(x)) x$tss_anchor0 <- .pos1_to_anchor0(x$tss_pos1)
  if (!"tss_id" %in% names(x)) {
    x$tss_id <- .stable_row_id(
      x$assembly, x$gene_id_normalized, x$tss_chrom, x$tss_pos1,
      x$tss_strand, prefix = "tss"
    )
  }
  if (!nrow(x)) return(x)
  x$tss_chrom <- .validate_chrom(x$tss_chrom, "genes$tss_chrom")
  x$tss_pos1 <- .as_integer_coordinate(x$tss_pos1, "genes$tss_pos1", 1)
  expected_anchor0 <- .pos1_to_anchor0(x$tss_pos1)
  supplied_anchor0 <- .as_integer_coordinate(x$tss_anchor0, "genes$tss_anchor0", 0)
  if (any(supplied_anchor0 != expected_anchor0)) {
    .acv2g_validation_error(
      "Standardized TSS anchor0 values disagree with their 1-based positions."
    )
  }
  x$tss_anchor0 <- supplied_anchor0
  invalid_gene <- is.na(x$gene_id_normalized) |
    !nzchar(trimws(as.character(x$gene_id_normalized)))
  invalid_assembly <- is.na(x$assembly) |
    !nzchar(trimws(as.character(x$assembly)))
  invalid_strand <- is.na(x$tss_strand) | !x$tss_strand %in% c("+", "-")
  invalid_tss_id <- is.na(x$tss_id) | !nzchar(as.character(x$tss_id))
  if (any(invalid_gene) || any(invalid_assembly) || any(invalid_strand) ||
      any(invalid_tss_id) || anyDuplicated(x$tss_id)) {
    .acv2g_validation_error(
      paste0(
        "Standardized TSS rows require non-missing gene/assembly/TSS IDs, ",
        "explicit strand, and unique exact-coordinate TSS IDs."
      )
    )
  }
  x
}
