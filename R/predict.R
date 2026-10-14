.read_quantitative_table <- function(path) {
  extension <- tolower(tools::file_ext(sub("\\.gz$", "", path, ignore.case = TRUE)))
  separator <- if (extension == "csv") "," else "\t"
  as.data.frame(
    data.table::fread(path, sep = separator, data.table = FALSE,
                      check.names = FALSE, na.strings = c("NA", "")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.read_optional_table <- function(spec) {
  if (is.null(spec)) return(NULL)
  if (is.data.frame(spec)) return(spec)
  path <- if (is.list(spec)) spec$path else spec
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    .acv2g_validation_error("Optional mapping/table input must be an existing file or data frame.")
  }
  .read_quantitative_table(path)
}

.read_coordinate_mapping <- function(spec, expected_assembly, name) {
  mapping <- .read_optional_table(spec)
  if (is.null(mapping) ||
      !all(c("chrom", "pos1") %in% names(mapping))) {
    return(mapping)
  }
  expected_assembly <- .assert_build(
    expected_assembly, paste0(name, " expected assembly")
  )
  declared <- if (is.list(spec) && !is.data.frame(spec)) {
    spec$assembly %||% NULL
  } else {
    NULL
  }
  if (!is.null(declared)) {
    declared <- .assert_build(declared, paste0(name, ".assembly"))
  }
  table_assembly <- NULL
  if ("assembly" %in% names(mapping)) {
    values <- unique(as.character(mapping$assembly))
    if (length(values) != 1L || is.na(values) || !nzchar(values)) {
      .acv2g_validation_error(
        paste0(name, " coordinate rows must declare one non-missing assembly.")
      )
    }
    table_assembly <- values[[1L]]
  }
  if (is.null(declared) && is.null(table_assembly)) {
    .acv2g_validation_error(
      paste0(
        name,
        " uses chrom + pos1 coordinates and must explicitly declare assembly ",
        "in its input specification or assembly column."
      )
    )
  }
  if (!is.null(declared) && !is.null(table_assembly) &&
      !identical(declared, table_assembly)) {
    .acv2g_validation_error(
      paste0(name, " input and row-level assembly declarations disagree.")
    )
  }
  resolved <- declared %||% table_assembly
  if (!identical(resolved, expected_assembly)) {
    .acv2g_validation_error(
      paste0(
        name, " coordinate assembly must equal ", expected_assembly,
        "; observed ", resolved, "."
      )
    )
  }
  mapping$assembly <- resolved
  mapping
}

.spec_call <- function(fun, primary, spec, aliases = list(), exclude = character()) {
  if (is.null(spec)) spec <- list()
  extras <- spec[setdiff(names(spec), exclude)]
  for (target in names(aliases)) {
    source <- aliases[[target]]
    if (!is.null(extras[[source]]) && is.null(extras[[target]])) {
      extras[[target]] <- extras[[source]]
    }
    extras[[source]] <- NULL
  }
  accepted <- names(formals(fun))
  extras <- extras[intersect(names(extras), accepted)]
  do.call(fun, c(primary, extras))
}

.resolve_gene_spec_path <- function(spec, config) {
  if (!is.null(spec$path)) return(spec)
  resolved <- resolve_v2g_reference(
    spec$managed_id,
    registry = config$reference_registry %||% NULL,
    cache_dir = config$reference_cache_dir %||% NULL
  )
  spec$path <- resolved$path
  spec$format <- spec$format %||% resolved$format
  if (!identical(spec$assembly, resolved$assembly)) {
    .acv2g_validation_error(
      "Managed annotation assembly does not match its explicit input declaration."
    )
  }
  spec$managed_reference_sha256 <- resolved$sha256
  spec
}

.read_chrom_lengths <- function(spec) {
  if (is.null(spec)) return(NULL)
  table <- .read_optional_table(spec)
  if (ncol(table) < 2L) {
    .acv2g_validation_error("Chromosome-length input requires at least two columns.")
  }
  if (!all(c("chrom", "length") %in% names(table))) {
    names(table)[1:2] <- c("chrom", "length")
  }
  table[c("chrom", "length")]
}

.liftover_loci <- function(loci, spec, analysis_build, liftover_config) {
  lifted <- liftover_v2g_coordinates(
    loci,
    chain = spec$chain %||% NULL,
    source_build = spec$assembly,
    target_build = analysis_build,
    record_type = "interval",
    role = "structural",
    chrom_col = "locus_chrom",
    start0_col = "locus_start0",
    end0_col = "locus_end0",
    id_col = "locus_id",
    failure_policy = liftover_config$structural_failure_policy,
    reciprocal_policy = liftover_config$reciprocal_liftover_policy,
    reverse_chain = spec$reverse_chain %||% NULL
  )
  keep <- lifted$data$liftover_status %in% c("identity", "mapped")
  lifted$qc$liftover_role <- lifted$qc$input_role
  lifted$qc$input_role <- "loci"
  mapped <- lifted$data[keep, , drop = FALSE]
  mapped$locus_assembly <- analysis_build
  class(mapped) <- c("v2g_loci", "data.frame")
  list(value = mapped, qc = lifted$qc,
       failed_loci = as.character(loci$locus_id[!keep]))
}

.liftover_variants <- function(
    variants,
    spec,
    analysis_build,
    liftover_config,
    qc_role = "target_variants") {
  lifted <- liftover_v2g_coordinates(
    variants$anchors,
    chain = spec$chain %||% NULL,
    source_build = spec$assembly,
    target_build = analysis_build,
    role = "structural",
    chrom_col = "chrom",
    pos1_col = "pos1",
    id_col = "variant_anchor_id",
    failure_policy = liftover_config$structural_failure_policy,
    reciprocal_policy = liftover_config$reciprocal_liftover_policy,
    reverse_chain = spec$reverse_chain %||% NULL
  )
  mapped <- lifted$data
  keep <- mapped$liftover_status %in% c("identity", "mapped")
  lifted$qc$liftover_role <- lifted$qc$input_role
  lifted$qc$input_role <- qc_role
  old_anchor <- as.character(mapped$variant_anchor_id)
  mapped <- mapped[keep, , drop = FALSE]
  mapped$variant_anchor_id_source <- as.character(mapped$variant_anchor_id)
  mapped$assembly <- analysis_build
  mapped$anchor0 <- mapped$pos1 - 1L
  mapped$variant_key <- paste(analysis_build, mapped$chrom, mapped$pos1, sep = ":")
  mapped$variant_anchor_id <- .stable_row_id(mapped$variant_key, prefix = "anchor")

  mapping <- data.frame(
    variant_anchor_id_source = old_anchor,
    variant_anchor_id = NA_character_,
    variant_key = NA_character_,
    assembly = analysis_build,
    chrom = NA_character_,
    pos1 = NA_real_,
    stringsAsFactors = FALSE
  )
  target_index <- match(mapping$variant_anchor_id_source,
                        mapped$variant_anchor_id_source)
  present <- !is.na(target_index)
  for (column in c("variant_anchor_id", "variant_key", "assembly", "chrom", "pos1")) {
    mapping[[column]][present] <- mapped[[column]][target_index[present]]
  }

  record_map <- variants$record_map
  record_map$variant_anchor_id_source <- record_map$variant_anchor_id
  target_index <- match(record_map$variant_anchor_id_source,
                        mapping$variant_anchor_id_source)
  record_map$variant_anchor_id <- mapping$variant_anchor_id[target_index]
  record_map$variant_key_source <- record_map$variant_key
  record_map$variant_key <- mapping$variant_key[target_index]
  record_map$assembly_source <- record_map$assembly
  record_map$assembly <- analysis_build
  record_map$chrom_source <- record_map$chrom
  record_map$chrom <- mapping$chrom[target_index]
  record_map$pos1_source <- record_map$pos1
  record_map$pos1 <- mapping$pos1[target_index]
  record_map$liftover_status <- lifted$qc$status[
    match(record_map$variant_anchor_id_source, lifted$qc$source_record_id)
  ]
  record_map$input_role <- qc_role

  mapped <- mapped[!duplicated(mapped$variant_anchor_id), , drop = FALSE]
  multiplicity <- table(record_map$variant_anchor_id[!is.na(record_map$variant_anchor_id)])
  mapped$anchor_multiplicity <- as.integer(multiplicity[mapped$variant_anchor_id])
  mapped$anchor_multiplicity[is.na(mapped$anchor_multiplicity)] <- 0L
  mapped <- mapped[c(
    "variant_anchor_id", "variant_key", "assembly", "chrom", "pos1",
    "anchor0", "anchor_multiplicity", "variant_anchor_id_source"
  )]
  structure_result <- structure(
    list(
      records = variants$records,
      anchors = mapped,
      record_map = record_map,
      qc = record_map,
      assembly = analysis_build
    ),
    class = "v2g_variants"
  )
  list(
    value = structure_result,
    qc = lifted$qc,
    failed_source_record_ids = lifted$qc$source_record_id[!keep],
    failed_loci = character()
  )
}

.liftover_genes <- function(genes, spec, analysis_build, liftover_config) {
  tss <- genes$tss
  lifted <- liftover_v2g_coordinates(
    tss,
    chain = spec$chain %||% NULL,
    source_build = spec$assembly,
    target_build = analysis_build,
    role = "structural",
    chrom_col = "tss_chrom",
    pos1_col = "tss_pos1",
    id_col = "tss_id",
    strand_col = "tss_strand",
    failure_policy = liftover_config$structural_failure_policy,
    reciprocal_policy = liftover_config$reciprocal_liftover_policy,
    reverse_chain = spec$reverse_chain %||% NULL
  )
  mapped <- lifted$data
  keep <- mapped$liftover_status %in% c("identity", "mapped")
  lifted$qc$liftover_role <- lifted$qc$input_role
  lifted$qc$input_role <- "genes"
  mapped <- mapped[keep, , drop = FALSE]
  mapped$tss_id_source <- mapped$tss_id
  mapped$assembly <- analysis_build
  mapped$tss_anchor0 <- mapped$tss_pos1 - 1L
  mapped$tss_id <- .stable_row_id(
    analysis_build, mapped$gene_id_normalized, mapped$tss_chrom,
    mapped$tss_pos1, mapped$tss_strand, prefix = "tss"
  )
  exact_key <- paste(mapped$gene_id_normalized, mapped$tss_chrom,
                     mapped$tss_pos1, mapped$tss_strand, sep = "\u001f")
  mapped <- mapped[!duplicated(exact_key), , drop = FALSE]

  gene_chrom <- split(mapped$tss_chrom, mapped$gene_id_normalized)
  genes$genes$gene_chrom <- vapply(genes$genes$gene_id_normalized, function(id) {
    value <- unique(gene_chrom[[id]])
    if (!length(value)) NA_character_ else if (length(value) == 1L) value else {
      .acv2g_validation_error(
        "Lifted TSSs for gene {.val {id}} map to more than one chromosome."
      )
    }
  }, character(1L))
  genes$genes$assembly <- analysis_build
  genes$tss <- mapped
  genes$assembly <- analysis_build
  list(
    value = genes,
    qc = lifted$qc,
    failed_source_record_ids = lifted$qc$source_record_id[!keep],
    failed_loci = character()
  )
}

.resolve_source_variant_locus_memberships <- function(
    loci,
    variants,
    variant_locus_map = NULL) {
  anchors <- variants$anchors
  if (is.null(variant_locus_map)) {
    locus_assemblies <- unique(as.character(loci$locus_assembly))
    variant_assemblies <- unique(as.character(anchors$assembly))
    if (length(locus_assemblies) != 1L ||
        length(variant_assemblies) != 1L ||
        !identical(locus_assemblies, variant_assemblies)) {
      .acv2g_validation_error(
        paste0(
          "Failed target variant anchors cannot be safely attributed to ",
          "fixed loci because their source assembly differs from the locus ",
          "source assembly. Supply an explicit inputs.variant_locus_map or ",
          "use liftover.structural_failure_policy = 'error'."
        )
      )
    }
    memberships <- .coordinate_memberships(loci, anchors)
    counts <- table(memberships$variant_anchor_id)
    ambiguous <- names(counts[counts > 1L])
    if (length(ambiguous)) {
      .acv2g_validation_error(
        paste0(
          "Failed target variant anchors cannot be safely attributed because ",
          "source-coordinate loci overlap. Supply an explicit ",
          "inputs.variant_locus_map; ambiguous anchor(s): ",
          paste(head(ambiguous, 10L), collapse = ", "), "."
        )
      )
    }
    return(unique(memberships[c("variant_anchor_id", "locus_id")]))
  }

  mapping <- .as_plain_data_frame(variant_locus_map, "variant_locus_map")
  .require_columns(mapping, "locus_id", "variant_locus_map")
  mapping$locus_id <- as.character(mapping$locus_id)
  mapping <- .resolve_explicit_membership_anchors(mapping, variants, anchors)
  unknown_locus <- is.na(mapping$locus_id) |
    !mapping$locus_id %in% loci$locus_id
  unknown_anchor <- is.na(mapping$variant_anchor_id) |
    !mapping$variant_anchor_id %in% anchors$variant_anchor_id
  if (any(unknown_locus) || any(unknown_anchor)) {
    .acv2g_validation_error(
      paste0(
        "Explicit variant-locus mapping contains ", sum(unknown_locus),
        " unknown locus reference(s) and ", sum(unknown_anchor),
        " variant reference(s) in the declared source variant universe."
      )
    )
  }
  unique(mapping[c("variant_anchor_id", "locus_id")])
}

.failed_tss_locus_memberships <- function(
    loci,
    genes,
    failed_tss_ids,
    gene_flank_bp) {
  failed_tss_ids <- unique(as.character(failed_tss_ids))
  if (!length(failed_tss_ids)) {
    return(data.frame(
      source_record_id = character(), locus_id = character(),
      stringsAsFactors = FALSE
    ))
  }
  locus_assemblies <- unique(as.character(loci$locus_assembly))
  gene_assemblies <- unique(as.character(genes$tss$assembly))
  if (length(locus_assemblies) != 1L || length(gene_assemblies) != 1L ||
      !identical(locus_assemblies, gene_assemblies)) {
    .acv2g_validation_error(
      paste0(
        "Failed gene TSSs cannot be safely attributed to fixed loci because ",
        "their source assembly differs from the locus source assembly. ",
        "Use liftover.structural_failure_policy = 'error' or harmonize the ",
        "structural input assemblies before running ACV2G."
      )
    )
  }
  failed_tss <- genes$tss[
    genes$tss$tss_id %in% failed_tss_ids, , drop = FALSE
  ]
  unresolved <- setdiff(failed_tss_ids, as.character(failed_tss$tss_id))
  if (length(unresolved)) {
    .acv2g_validation_error(
      paste0(
        "Internal structural lifecycle error: failed TSS ID(s) are absent ",
        "from the source gene universe: ",
        paste(head(unresolved, 10L), collapse = ", "), "."
      )
    )
  }
  gene_flank_bp <- .assert_scalar_number(
    gene_flank_bp, "gene_flank_bp", lower = 0, integer = TRUE
  )
  pieces <- lapply(seq_len(nrow(loci)), function(i) {
    start0 <- max(0, as.numeric(loci$locus_start0[[i]]) - gene_flank_bp)
    end0 <- as.numeric(loci$locus_end0[[i]]) + gene_flank_bp
    keep <- failed_tss$tss_chrom == loci$locus_chrom[[i]] &
      failed_tss$tss_anchor0 >= start0 & failed_tss$tss_anchor0 < end0
    if (!any(keep)) return(NULL)
    data.frame(
      source_record_id = as.character(failed_tss$tss_id[keep]),
      locus_id = as.character(loci$locus_id[[i]]),
      stringsAsFactors = FALSE
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) {
    return(data.frame(
      source_record_id = character(), locus_id = character(),
      stringsAsFactors = FALSE
    ))
  }
  unique(do.call(rbind, pieces))
}

.attribute_structural_failures <- function(
    loci,
    variants,
    genes,
    failed_locus_ids = character(),
    failed_variant_ids = character(),
    failed_tss_ids = character(),
    variant_locus_map = NULL,
    gene_flank_bp = 1e6) {
  failed_locus_ids <- unique(as.character(failed_locus_ids))
  failed_variant_ids <- unique(as.character(failed_variant_ids))
  failed_tss_ids <- unique(as.character(failed_tss_ids))
  unknown_loci <- setdiff(failed_locus_ids, as.character(loci$locus_id))
  if (length(unknown_loci)) {
    .acv2g_validation_error(
      "Internal structural lifecycle error: a failed locus ID is unknown."
    )
  }

  locus_rows <- data.frame(
    input_role = rep("loci", length(failed_locus_ids)),
    source_record_id = failed_locus_ids,
    locus_id = failed_locus_ids,
    stringsAsFactors = FALSE
  )
  variant_rows <- data.frame(
    input_role = character(), source_record_id = character(),
    locus_id = character(), stringsAsFactors = FALSE
  )
  if (length(failed_variant_ids)) {
    memberships <- .resolve_source_variant_locus_memberships(
      loci, variants, variant_locus_map
    )
    unresolved <- setdiff(
      failed_variant_ids, as.character(variants$anchors$variant_anchor_id)
    )
    if (length(unresolved)) {
      .acv2g_validation_error(
        "Internal structural lifecycle error: a failed variant anchor ID is unknown."
      )
    }
    memberships <- memberships[
      memberships$variant_anchor_id %in% failed_variant_ids, , drop = FALSE
    ]
    variant_rows <- data.frame(
      input_role = rep("target_variants", nrow(memberships)),
      source_record_id = as.character(memberships$variant_anchor_id),
      locus_id = as.character(memberships$locus_id),
      stringsAsFactors = FALSE
    )
  }
  gene_memberships <- .failed_tss_locus_memberships(
    loci, genes, failed_tss_ids, gene_flank_bp
  )
  gene_rows <- data.frame(
    input_role = rep("genes", nrow(gene_memberships)),
    source_record_id = gene_memberships$source_record_id,
    locus_id = gene_memberships$locus_id,
    stringsAsFactors = FALSE
  )
  rows <- data.table::rbindlist(
    list(locus_rows, variant_rows, gene_rows), fill = TRUE
  ) |> as.data.frame(stringsAsFactors = FALSE)
  rows <- unique(rows)
  list(
    failed_loci = unique(as.character(rows$locus_id)),
    record_locus_map = rows
  )
}

.translate_variant_locus_map <- function(
    variant_locus_map,
    source_loci,
    source_variants,
    mapped_variants) {
  if (is.null(variant_locus_map)) return(NULL)
  memberships <- .resolve_source_variant_locus_memberships(
    source_loci, source_variants, variant_locus_map
  )
  translation <- unique(mapped_variants$anchors[c(
    "variant_anchor_id_source", "variant_anchor_id"
  )])
  names(translation)[[1L]] <- "variant_anchor_id_source"
  names(memberships)[names(memberships) == "variant_anchor_id"] <-
    "variant_anchor_id_source"
  mapped <- merge(
    memberships, translation, by = "variant_anchor_id_source",
    all.x = TRUE, sort = FALSE
  )
  mapped <- mapped[!is.na(mapped$variant_anchor_id), , drop = FALSE]
  unique(mapped[c("variant_anchor_id", "locus_id")])
}

.annotate_liftover_locus_effects <- function(qc, record_locus_map) {
  qc$affected_locus_ids <- NA_character_
  if (!nrow(record_locus_map)) return(qc)
  keys <- paste(record_locus_map$input_role,
                record_locus_map$source_record_id, sep = "\u001f")
  collapsed <- vapply(split(record_locus_map$locus_id, keys), function(x) {
    paste(sort(unique(as.character(x))), collapse = ";")
  }, character(1L))
  qc_key <- paste(qc$input_role, qc$source_record_id, sep = "\u001f")
  qc$affected_locus_ids <- unname(collapsed[qc_key])
  qc
}

.prepare_structural_inputs <- function(config) {
  loci_spec <- config$inputs$loci
  variants_spec <- config$inputs$variants
  genes_spec <- .resolve_gene_spec_path(config$inputs$genes, config)
  protein_ids <- genes_spec$protein_coding_gene_ids %||% NULL
  if (is.character(protein_ids) && length(protein_ids) == 1L &&
      file.exists(protein_ids)) {
    ids_table <- .read_quantitative_table(protein_ids)
    protein_ids <- ids_table[[1L]]
  }
  loci <- .spec_call(
    read_v2g_loci,
    list(x = loci_spec$path, assembly = loci_spec$assembly),
    loci_spec,
    exclude = c("path", "assembly", "chain", "reverse_chain")
  )
  variants <- .spec_call(
    read_v2g_variants,
    list(x = variants_spec$path, assembly = variants_spec$assembly),
    variants_spec,
    exclude = c("path", "assembly", "chain", "reverse_chain")
  )
  genes_spec$protein_coding_gene_ids <- protein_ids
  genes <- .spec_call(
    read_v2g_genes,
    list(x = genes_spec$path, assembly = genes_spec$assembly),
    genes_spec,
    exclude = c(
      "path", "assembly", "chain", "reverse_chain", "managed_id",
      "managed_reference_sha256"
    )
  )
  variant_locus_map <- .read_coordinate_mapping(
    config$inputs$variant_locus_map %||% NULL,
    expected_assembly = variants_spec$assembly,
    name = "inputs.variant_locus_map"
  )

  liftover <- config$liftover
  mapped_loci <- .liftover_loci(loci, loci_spec, config$analysis_build, liftover)
  mapped_variants <- .liftover_variants(
    variants, variants_spec, config$analysis_build, liftover
  )
  mapped_genes <- .liftover_genes(genes, genes_spec, config$analysis_build, liftover)
  attribution <- .attribute_structural_failures(
    loci = loci,
    variants = variants,
    genes = genes,
    failed_locus_ids = mapped_loci$failed_loci,
    failed_variant_ids = mapped_variants$failed_source_record_ids,
    failed_tss_ids = mapped_genes$failed_source_record_ids,
    variant_locus_map = variant_locus_map,
    gene_flank_bp = config$method$gene_flank_bp
  )
  failed_loci <- attribution$failed_loci
  mapped_loci$value <- mapped_loci$value[
    !mapped_loci$value$locus_id %in% failed_loci, , drop = FALSE
  ]
  mapped_variant_locus_map <- .translate_variant_locus_map(
    variant_locus_map, loci, variants, mapped_variants$value
  )
  if (!is.null(mapped_variant_locus_map) && length(failed_loci)) {
    mapped_variant_locus_map <- mapped_variant_locus_map[
      !mapped_variant_locus_map$locus_id %in% failed_loci, , drop = FALSE
    ]
  }
  liftover_qc <- do.call(rbind, lapply(
    list(mapped_loci$qc, mapped_variants$qc, mapped_genes$qc),
    function(x) {
      x$schema_version <- .acv2g_schema_version
      x
    }
  ))
  liftover_qc <- .annotate_liftover_locus_effects(
    liftover_qc, attribution$record_locus_map
  )
  list(
    loci = mapped_loci$value,
    variants = mapped_variants$value,
    genes = mapped_genes$value,
    variant_locus_map = mapped_variant_locus_map,
    liftover_qc = liftover_qc,
    failed_loci = failed_loci,
    all_locus_ids = as.character(loci$locus_id)
  )
}

.key_string <- function(x, columns) {
  .require_columns(x, columns, "join table")
  do.call(paste, c(lapply(x[columns], as.character), sep = "\u001f"))
}

.left_join_unique <- function(x, y, by, y_name = "evidence") {
  x <- .as_plain_data_frame(x)
  y <- .as_plain_data_frame(y)
  x_key <- .key_string(x, names(by))
  y_columns <- unname(by)
  y_key <- .key_string(y, y_columns)
  if (anyDuplicated(y_key)) {
    .acv2g_validation_error("{y_name} contains duplicate join keys.")
  }
  index <- match(x_key, y_key)
  additions <- setdiff(names(y), y_columns)
  additions <- setdiff(additions, names(x))
  for (column in additions) x[[column]] <- y[[column]][index]
  x
}

.profile_source <- function(spec, type) {
  source_type <- spec$type %||% "auto"
  if (source_type == "table") .read_quantitative_table(spec$path) else spec$path
}

.activity_query_chrom_lengths <- function(
    spec,
    analysis_chrom_lengths,
    analysis_build) {
  supplied <- spec$chrom_sizes %||% NULL
  if (!is.null(supplied)) {
    if (is.numeric(supplied) && !is.null(names(supplied))) return(supplied)
    if (is.data.frame(supplied)) return(supplied)
    return(.read_chrom_lengths(supplied))
  }
  if (identical(spec$assembly, analysis_build)) analysis_chrom_lengths else NULL
}

.activity_extraction_call <- function(
    query_anchors,
    source,
    spec,
    config,
    chrom_lengths,
    profile_id) {
  extract_v2g_activity(
    variants = query_anchors,
    activity = source,
    source_type = spec$type %||% "auto",
    chrom_lengths = chrom_lengths,
    activity_window_bp = spec$activity_window_bp %||%
      config$activity$activity_window_bp,
    activity_aggregate = spec$aggregation %||% config$activity$aggregation,
    activity_uncovered_policy = spec$uncovered_policy %||% NULL,
    activity_transform = config$method$activity_transform %||%
      config$method$proxy_transform,
    activity_assay_id = spec$assay_id %||% profile_id,
    activity_assembly = spec$assembly,
    activity_units = spec$units %||% NA_character_,
    value_col = spec$value_col %||% NULL,
    table_mode = spec$table_mode %||% "auto",
    bam_mapq = spec$bam_mapq %||% config$activity$bam_mapq_min,
    bam_library_size = spec$bam_library_size %||% NULL
  )
}

.typed_missing <- function(x, n) {
  if (is.character(x) || is.factor(x)) return(rep(NA_character_, n))
  if (is.integer(x)) return(rep(NA_integer_, n))
  if (is.logical(x)) return(rep(NA, n))
  if (is.numeric(x)) return(rep(NA_real_, n))
  rep(NA_character_, n)
}

.complete_activity_query_result <- function(
    source_anchors,
    lifted,
    extracted,
    spec,
    profile_id,
    config,
    source_type) {
  n <- nrow(source_anchors)
  out <- data.frame(
    variant_anchor_id = as.character(source_anchors$variant_anchor_id),
    variant_id = as.character(source_anchors$variant_anchor_id),
    chrom = as.character(source_anchors$chrom),
    pos1 = as.numeric(source_anchors$pos1),
    stringsAsFactors = FALSE
  )
  if (nrow(extracted)) {
    names(extracted)[names(extracted) == "chrom"] <- "activity_query_chrom"
    names(extracted)[names(extracted) == "pos1"] <- "activity_query_pos1"
    extracted_index <- match(out$variant_anchor_id, extracted$variant_anchor_id)
    for (column in setdiff(
      names(extracted), c("variant_anchor_id", "variant_id")
    )) {
      out[[column]] <- .typed_missing(extracted[[column]], n)
      present <- !is.na(extracted_index)
      out[[column]][present] <- extracted[[column]][extracted_index[present]]
    }
  }
  required_numeric <- c(
    "activity_query_pos1", "activity_window_requested_start0",
    "activity_window_requested_end0", "activity_window_effective_start0",
    "activity_window_effective_end0", "activity_effective_window_bp",
    "activity_raw", "activity_transformed", "activity_covered_bp",
    "activity_coverage_fraction", "activity_bam_library_size"
  )
  required_character <- c(
    "activity_query_chrom", "activity_assay_id", "activity_assembly",
    "activity_units", "activity_source_type", "activity_aggregate",
    "activity_uncovered_policy", "activity_transform", "activity_status",
    "activity_query_failure_message", "activity_bam_mode",
    "activity_bam_filter_policy"
  )
  required_integer <- c("activity_window_bp", "activity_bam_mapq")
  for (column in setdiff(required_numeric, names(out))) out[[column]] <- NA_real_
  for (column in setdiff(required_character, names(out))) {
    out[[column]] <- NA_character_
  }
  for (column in setdiff(required_integer, names(out))) out[[column]] <- NA_integer_
  if (!"activity_window_clipped" %in% names(out)) {
    out$activity_window_clipped <- NA
  }
  if (!"activity_query_failed" %in% names(out)) {
    out$activity_query_failed <- NA
  }

  liftover_index <- match(
    out$variant_anchor_id, lifted$data$variant_anchor_id
  )
  out$activity_query_chrom <- as.character(
    lifted$data$chrom[liftover_index]
  )
  out$activity_query_pos1 <- as.numeric(lifted$data$pos1[liftover_index])
  out$activity_query_liftover_status <- as.character(
    lifted$data$liftover_status[liftover_index]
  )
  out$activity_query_source_build <- config$analysis_build
  out$activity_query_target_build <- spec$assembly
  out$activity_assay_id <- spec$assay_id %||% profile_id
  out$activity_assembly <- spec$assembly
  out$activity_source_type <- source_type
  if (identical(source_type, "bam")) {
    out$activity_units <- "fragment_coverage_cpm"
  } else {
    out$activity_units <- spec$units %||% NA_character_
  }
  out$activity_window_bp <- as.integer(
    spec$activity_window_bp %||% config$activity$activity_window_bp
  )
  out$activity_aggregate <- spec$aggregation %||% config$activity$aggregation
  out$activity_uncovered_policy <- spec$uncovered_policy %||%
    if (identical(out$activity_aggregate[[1L]], "window_mean")) "zero" else
      "exclude"
  out$activity_transform <- config$method$activity_transform %||%
    config$method$proxy_transform
  out$activity_bam_mapq <- if (identical(source_type, "bam")) {
    as.integer(spec$bam_mapq %||% config$activity$bam_mapq_min)
  } else {
    NA_integer_
  }
  out$activity_bam_mode <- if (identical(source_type, "bam")) {
    "fragment_coverage_cpm"
  } else {
    NA_character_
  }
  out$activity_bam_filter_policy <- if (identical(source_type, "bam")) {
    "mapped_primary_qcpass_nonduplicate_nonsupplementary_v1"
  } else {
    NA_character_
  }
  failed <- !out$activity_query_liftover_status %in% c("identity", "mapped")
  out$activity_status[failed] <- paste0(
    "liftover_", out$activity_query_liftover_status[failed]
  )
  for (column in c(
    "activity_raw", "activity_transformed", "activity_covered_bp",
    "activity_coverage_fraction", "activity_window_requested_start0",
    "activity_window_requested_end0", "activity_window_effective_start0",
    "activity_window_effective_end0", "activity_effective_window_bp",
    "activity_bam_library_size"
  )) {
    out[[column]][failed] <- NA_real_
  }
  out$activity_window_clipped[failed] <- NA
  out$activity_query_failed[failed] <- NA
  out$activity_query_failure_message[failed] <- NA_character_
  out
}

.activity_for_profile <- function(candidates, profile, config, chrom_lengths) {
  anchors <- unique(candidates$pairs[c("variant_anchor_id", "chrom", "pos1")])
  anchors$variant_id <- anchors$variant_anchor_id
  spec <- profile$activity
  source <- .profile_source(spec, "activity")
  lifted <- liftover_v2g_coordinates(
    anchors,
    chain = spec$chain %||% NULL,
    source_build = config$analysis_build,
    target_build = spec$assembly,
    role = "evidence_query",
    chrom_col = "chrom",
    pos1_col = "pos1",
    id_col = "variant_anchor_id",
    failure_policy = config$liftover$structural_failure_policy,
    reciprocal_policy = config$liftover$reciprocal_liftover_policy,
    reverse_chain = spec$reverse_chain %||% NULL
  )
  keep <- lifted$data$liftover_status %in% c("identity", "mapped")
  query_anchors <- lifted$data[keep, c(
    "variant_anchor_id", "chrom", "pos1"
  ), drop = FALSE]
  query_anchors$variant_id <- query_anchors$variant_anchor_id
  query_lengths <- .activity_query_chrom_lengths(
    spec, chrom_lengths, config$analysis_build
  )
  extracted <- if (nrow(query_anchors)) {
    .activity_extraction_call(
      query_anchors, source, spec, config, query_lengths,
      profile$profile_id
    )
  } else {
    data.frame()
  }
  source_type <- .activity_source_type(source, spec$type %||% "auto")
  activity <- .complete_activity_query_result(
    anchors, lifted, extracted, spec, profile$profile_id, config, source_type
  )
  lifted$qc$liftover_role <- lifted$qc$input_role
  lifted$qc$input_role <- "activity_query"
  lifted$qc$profile_id <- profile$profile_id
  list(value = activity, qc = lifted$qc)
}

.contact_for_profile <- function(candidates, profile, config) {
  pair_tss <- candidates$tss
  pair_tss$variant_id <- pair_tss$variant_anchor_id
  pair_tss$gene_id <- pair_tss$gene_id_normalized
  spec <- profile$contact
  source <- .profile_source(spec, "contact")
  extract_v2g_contacts(
    pairs = pair_tss,
    contacts = source,
    source_type = spec$type %||% "auto",
    contact_resolution_bp = spec$resolution_bp %||%
      config$contact$contact_resolution_bp,
    contact_value = spec$contact_value %||% config$contact$contact_value,
    contact_layer = spec[["layer"]] %||% NULL,
    contact_score_col = spec$score_col %||% NULL,
    contact_score_type = spec$score_type %||% NULL,
    contact_layers = spec$layers %||% NULL,
    table_representation = spec$table_representation %||% "complete",
    diagonal_policy = spec$diagonal_policy %||%
      config$contact$diagonal_policy,
    distance_policy = config$contact$distance_policy,
    contact_aggregate = "max",
    contact_transform = config$method$contact_transform %||%
      config$method$proxy_transform,
    contact_assay_id = spec$assay_id %||% profile$profile_id,
    contact_assembly = spec$assembly,
    contact_units = spec$units %||% NA_character_,
    contact_normalization = spec$normalization %||% NA_character_,
    contact_weights_source_provided = spec$weights_source_provided %||% NA
  )
}

.evidence_for_profile <- function(
    candidates,
    profile,
    config,
    chrom_lengths,
    evidence_scope = "target") {
  activity_query <- .activity_for_profile(
    candidates, profile, config, chrom_lengths
  )
  activity <- activity_query$value
  contacts <- .contact_for_profile(candidates, profile, config)
  pairs <- candidates$pairs
  pairs$profile_id <- profile$profile_id
  activity$profile_id <- profile$profile_id
  pairs <- .left_join_unique(
    pairs,
    activity,
    by = c(variant_anchor_id = "variant_anchor_id"),
    y_name = "activity evidence"
  )
  contact_pairs <- contacts$pairs
  contact_pairs$variant_anchor_id <- contact_pairs$variant_id
  contact_pairs$gene_id_normalized <- contact_pairs$gene_id
  contact_by <- c(
    variant_anchor_id = "variant_anchor_id",
    gene_id_normalized = "gene_id_normalized"
  )
  if ("locus_id" %in% names(pairs) && "locus_id" %in% names(contact_pairs)) {
    contact_by <- c(locus_id = "locus_id", contact_by)
  }
  pairs <- .left_join_unique(
    pairs,
    contact_pairs,
    by = contact_by,
    y_name = "contact evidence"
  )
  pairs <- compute_v2g_proxy(
    pairs,
    activity_col = "activity_raw",
    contact_col = "contact_raw",
    activity_transform = config$method$activity_transform %||%
      config$method$proxy_transform,
    contact_transform = config$method$contact_transform %||%
      config$method$proxy_transform,
    activity_output_col = "activity_transformed",
    contact_output_col = "contact_transformed",
    output_col = "proxy_score",
    row_id_col = "pair_id"
  )
  activity$status <- activity$activity_status
  contacts$qc$status <- "ok"
  contacts$qc$profile_id <- profile$profile_id
  contacts$tss_bin_contacts$profile_id <- profile$profile_id
  activity_query$qc$evidence_scope <- evidence_scope
  list(
    pairs = pairs,
    activity_qc = activity,
    contact_qc = contacts$qc,
    tss_bin_contacts = contacts$tss_bin_contacts,
    activity_liftover_qc = activity_query$qc,
    contact_result = contacts
  )
}

.score_profile_locus <- function(candidates, profile, config, chrom_lengths) {
  evidence <- .evidence_for_profile(candidates, profile, config, chrom_lengths)
  normalized <- compute_locus_pair_z(
    evidence$pairs,
    proxy_col = "proxy_score",
    group_cols = c("profile_id", "locus_id"),
    warn_missing_fraction = config$method$max_missing_fraction,
    z_col = "locus_pair_z"
  )
  scored <- call_v2g_links(
    normalized$pairs,
    call_rule = config$method$call_rule,
    z_col = "locus_pair_z",
    z_threshold = if (config$method$call_rule == "z_threshold")
      config$method$z_threshold else NULL
  )
  normalized$normalization_qc$profile_id <- profile$profile_id
  contacts <- evidence$contact_result
  list(
    pairs = scored,
    normalization_qc = normalized$normalization_qc,
    activity_qc = evidence$activity_qc,
    contact_qc = evidence$contact_qc,
    tss_bin_contacts = evidence$tss_bin_contacts,
    liftover_qc = evidence$activity_liftover_qc,
    profile = .profile_row_from_evidence(profile, evidence)
  )
}

.combine_tables <- function(results, name) {
  pieces <- lapply(results, `[[`, name)
  pieces <- pieces[!vapply(pieces, is.null, logical(1L))]
  if (!length(pieces)) return(data.frame())
  data.table::rbindlist(pieces, fill = TRUE, use.names = TRUE) |>
    as.data.frame(stringsAsFactors = FALSE)
}

.combine_warning_tables <- function(...) {
  pieces <- list(...)
  pieces <- pieces[!vapply(pieces, is.null, logical(1L))]
  pieces <- pieces[vapply(pieces, nrow, integer(1L)) > 0L]
  if (!length(pieces)) return(.default_result_table("warnings"))
  warnings <- data.table::rbindlist(
    pieces, fill = TRUE, use.names = TRUE
  ) |>
    as.data.frame(stringsAsFactors = FALSE)
  warnings$warning_id <- paste0(
    "warning_", sprintf("%06d", seq_len(nrow(warnings)))
  )
  warnings
}

.canonical_reference_rows <- function(x, columns) {
  columns <- intersect(columns, names(x))
  if (!length(columns)) return(data.frame())
  out <- unique(.as_plain_data_frame(x)[columns])
  if (nrow(out)) {
    ordering <- lapply(out, function(value) {
      value <- as.character(value)
      value[is.na(value)] <- "<NA>"
      value
    })
    out <- out[do.call(order, c(ordering, list(na.last = TRUE))), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.reference_mapping_id <- function(mapping) {
  if (is.null(mapping)) return(NA_character_)
  canonical <- .canonical_reference_rows(
    mapping, c("gene_id_normalized", "variant_anchor_id")
  )
  paste0(
    "mapping_",
    digest::digest(canonical, algo = "xxhash64", serialize = TRUE)
  )
}

.reference_identity <- function(reference_candidates, config, explicit_mapping) {
  mapping <- .canonical_reference_rows(
    reference_candidates$mapping,
    c("gene_id_normalized", "variant_anchor_id")
  )
  target_gene_ids <- unique(as.character(
    reference_candidates$target_gene_ids %||%
      reference_candidates$pairs$gene_id_normalized
  ))
  tss <- reference_candidates$genes$tss
  tss <- tss[tss$gene_id_normalized %in% target_gene_ids, , drop = FALSE]
  tss <- .canonical_reference_rows(
    tss,
    c(
      "gene_id_normalized", "tss_id", "assembly", "tss_chrom",
      "tss_pos1", "tss_strand"
    )
  )
  anchors <- .canonical_reference_rows(
    reference_candidates$variants$anchors,
    c("variant_anchor_id", "variant_key", "assembly", "chrom", "pos1")
  )
  rule <- if (is.null(explicit_mapping)) "distance" else "explicit_mapping"
  payload <- list(
    schema_version = .acv2g_schema_version,
    reference_mode = config$method$reference_mode,
    reference_universe_rule = rule,
    reference_gene_radius_bp = if (identical(rule, "distance")) {
      config$method$reference_gene_radius_bp
    } else {
      NA_real_
    },
    source_variant_anchors = anchors,
    gene_variant_associations = mapping,
    candidate_tss = tss
  )
  list(
    reference_set_id = paste0(
      config$method$reference_mode, "_",
      digest::digest(payload, algo = "xxhash64", serialize = TRUE)
    ),
    mapping_id = .reference_mapping_id(explicit_mapping),
    payload = payload
  )
}

.result_variant_record_map <- function(
    target_variants,
    candidates,
    reference_candidates = NULL,
    background_variants = NULL,
    reference_mode = NULL) {
  target_pair_anchors <- unique(as.character(candidates$pairs$variant_anchor_id))
  reference_pair_anchors <- if (is.null(reference_candidates)) {
    character()
  } else {
    unique(as.character(reference_candidates$mapping$variant_anchor_id))
  }

  annotate <- function(variants, input_role, is_reference_source) {
    if (is.null(variants)) return(NULL)
    out <- .as_plain_data_frame(variants$record_map)
    out$input_role <- input_role
    out$target_pair_member <- !is.na(out$variant_anchor_id) &
      out$variant_anchor_id %in% target_pair_anchors
    out$reference_source_member <- isTRUE(is_reference_source) &
      !is.na(out$variant_anchor_id)
    out$reference_universe_member <- if (is.null(reference_candidates)) {
      rep(NA, nrow(out))
    } else {
      isTRUE(is_reference_source) & !is.na(out$variant_anchor_id) &
        out$variant_anchor_id %in% reference_pair_anchors
    }
    out
  }

  target_is_reference <- identical(reference_mode, "input_variants")
  background_is_reference <- identical(reference_mode, "background_variants")
  maps <- list(
    annotate(target_variants, "target_variants", target_is_reference),
    annotate(background_variants, "background_variants", background_is_reference)
  )
  maps <- maps[!vapply(maps, is.null, logical(1L))]
  out <- data.table::rbindlist(maps, fill = TRUE, use.names = TRUE) |>
    as.data.frame(stringsAsFactors = FALSE)
  out <- out[order(out$input_role, out$record_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.profile_with_warnings <- function(expr, profile_id) {
  captured <- .capture_warnings(expr)
  if (nrow(captured$warnings)) {
    missing_profile <- is.na(captured$warnings$profile_id) |
      !nzchar(captured$warnings$profile_id)
    captured$warnings$profile_id[missing_profile] <- profile_id
  }
  captured
}

.candidate_locus_status <- function(
    candidates,
    profile_ids,
    failed_loci,
    all_locus_ids = NULL) {
  locus <- candidates$qc$loci
  if (is.null(locus)) return(.default_result_table("locus_status"))
  failed_loci <- unique(as.character(failed_loci))
  missing_failed <- setdiff(failed_loci, as.character(locus$locus_id))
  if (length(missing_failed)) {
    failed_rows <- data.frame(
      locus_id = missing_failed,
      candidate_gene_count = 0L,
      candidate_pair_count = 0L,
      candidate_status = "failed_structural_liftover",
      stringsAsFactors = FALSE
    )
    locus <- data.table::rbindlist(
      list(locus, failed_rows), fill = TRUE, use.names = TRUE
    ) |> as.data.frame(stringsAsFactors = FALSE)
  }
  if (!is.null(all_locus_ids)) {
    absent <- setdiff(as.character(all_locus_ids), as.character(locus$locus_id))
    if (length(absent)) {
      .acv2g_validation_error(
        "Internal structural lifecycle error: locus status omitted an input locus."
      )
    }
  }
  rows <- lapply(profile_ids, function(profile_id) {
    out <- locus
    out$profile_id <- profile_id
    out$status <- out$candidate_status
    out
  })
  out <- .combine_tables(lapply(rows, function(x) list(value = x)), "value")
  if (length(failed_loci)) {
    failed <- out$locus_id %in% failed_loci
    out$status[failed] <- "failed_structural_liftover"
  }
  out
}

.all_fixed_loci_failed <- function(structural) {
  length(structural$all_locus_ids) > 0L &&
    !nrow(structural$loci) &&
    setequal(structural$failed_loci, structural$all_locus_ids)
}

.failed_locus_profiles <- function(config) {
  rows <- lapply(config$profiles, function(profile) {
    data.frame(
      profile_id = profile$profile_id,
      activity_source = profile$activity$path,
      contact_source = profile$contact$path,
      activity_type = profile$activity$type %||% "auto",
      contact_type = profile$contact$type %||% "auto",
      status = "not_run_no_analyzable_loci",
      stringsAsFactors = FALSE
    )
  })
  data.table::rbindlist(rows, fill = TRUE, use.names = TRUE) |>
    as.data.frame(stringsAsFactors = FALSE)
}

.failed_locus_variant_qc <- function(structural) {
  anchors <- structural$variants$anchors
  if (!nrow(anchors)) return(.default_result_table("variant_locus_qc"))
  out <- anchors[c(
    "variant_anchor_id", "variant_key", "assembly", "chrom", "pos1"
  )]
  out$variant_locus_status <- "not_scored_failed_locus_scope"
  out$locus_membership_count <- 0L
  out$locus_assignment_mode <- "unavailable_structural_liftover"
  out
}

.all_failed_loci_result <- function(config, structural, warnings) {
  profiles <- .failed_locus_profiles(config)
  profile_locus <- expand.grid(
    profile_id = names(config$profiles),
    locus_id = structural$all_locus_ids,
    stringsAsFactors = FALSE
  )
  profile_locus$status <- "failed_structural_liftover"
  normalization_qc <- profile_locus
  normalization_qc$method <- config$method$name
  normalization_qc$group_id <- normalization_qc$locus_id
  normalization_qc <- normalization_qc[c(
    "profile_id", "method", "group_id", "status", "locus_id"
  )]
  annotation_qc <- structural$genes$annotation_qc
  annotation_qc$status <- annotation_qc$annotation_status
  empty_candidates <- list(pairs = data.frame(
    variant_anchor_id = character(), stringsAsFactors = FALSE
  ))
  new_v2g_result(
    pairs = .default_result_table("pairs"),
    profiles = profiles,
    variant_record_anchor_map = .result_variant_record_map(
      structural$variants, empty_candidates
    ),
    variant_locus_qc = .failed_locus_variant_qc(structural),
    locus_status = profile_locus,
    normalization_qc = normalization_qc,
    annotation_qc = annotation_qc,
    liftover_qc = structural$liftover_qc,
    warnings = warnings,
    config = config,
    manifest = build_v2g_manifest(config, profiles, extra = list(
      no_analyzable_loci = TRUE,
      failed_loci = structural$failed_loci
    )),
    run_status = list(
      status = "partial",
      failed_loci = structural$failed_loci,
      analyzable_locus_count = 0L,
      completed_at_utc = .utc_now()
    )
  )
}

.build_candidates_from_config <- function(structural, config) {
  build_v2g_candidates(
    structural$loci,
    structural$variants,
    structural$genes,
    gene_flank_bp = config$method$gene_flank_bp,
    variant_locus_map = structural$variant_locus_map
  )
}

.read_and_map_background_variants <- function(config) {
  spec <- config$inputs$background_variants
  variants <- .spec_call(
    read_v2g_variants,
    list(x = spec$path, assembly = spec$assembly),
    spec,
    exclude = c("path", "assembly", "chain", "reverse_chain")
  )
  mapped <- .liftover_variants(
    variants, spec, config$analysis_build, config$liftover,
    qc_role = "background_variants"
  )
  if (length(mapped$failed_source_record_ids)) {
    .acv2g_validation_error(
      paste0(
        "Background/reference variant liftover failed for ",
        length(mapped$failed_source_record_ids),
        " anchor(s). These failures change one or more gene reference ",
        "universes and cannot be attributed safely to fixed target loci; ",
        "use a complete analysis-build background or ",
        "liftover.structural_failure_policy = 'error'."
      ),
      data = mapped$qc[
        mapped$qc$source_record_id %in% mapped$failed_source_record_ids,
        , drop = FALSE
      ]
    )
  }
  mapped
}

.resolve_reference_mapping <- function(mapping, variants, genes) {
  mapping <- .as_plain_data_frame(mapping, "reference_gene_variant_map")
  gene_field <- intersect(c("gene_id_normalized", "gene_id"), names(mapping))
  if (!length(gene_field)) {
    .acv2g_validation_error(
      "Reference mapping requires gene_id_normalized or gene_id."
    )
  }
  mapping$gene_id_normalized <- as.character(mapping[[gene_field[[1L]]]])
  gene_table <- genes$genes
  raw_match <- match(mapping$gene_id_normalized, gene_table$gene_id_raw)
  normalized_match <- match(mapping$gene_id_normalized,
                            gene_table$gene_id_normalized)
  use_raw <- is.na(normalized_match) & !is.na(raw_match)
  mapping$gene_id_normalized[use_raw] <-
    gene_table$gene_id_normalized[raw_match[use_raw]]
  unknown <- !mapping$gene_id_normalized %in% gene_table$gene_id_normalized
  if (any(unknown)) {
    .acv2g_validation_error(
      "Reference mapping contains unknown gene ID{?s}: {paste(unique(mapping$gene_id_normalized[unknown]), collapse = ', ')}."
    )
  }
  mapping <- .resolve_explicit_membership_anchors(
    mapping, variants, variants$anchors
  )
  if (anyNA(mapping$variant_anchor_id)) {
    .acv2g_validation_error(
      "Reference mapping contains variant records absent from its declared reference variant universe."
    )
  }
  unique(mapping[c("gene_id_normalized", "variant_anchor_id")])
}

.distance_reference_pairs <- function(anchors, tss, radius) {
  if (!nrow(anchors) || !nrow(tss)) {
    return(data.frame(
      variant_anchor_id = character(), gene_id_normalized = character(),
      stringsAsFactors = FALSE
    ))
  }
  variant_gr <- GenomicRanges::GRanges(
    seqnames = anchors$chrom,
    ranges = IRanges::IRanges(start = anchors$pos1, width = 1L)
  )
  tss_start <- pmax(1, tss$tss_pos1 - radius)
  tss_end <- tss$tss_pos1 + radius
  reference_gr <- GenomicRanges::GRanges(
    seqnames = tss$tss_chrom,
    ranges = IRanges::IRanges(start = tss_start, end = tss_end)
  )
  hits <- GenomicRanges::findOverlaps(
    variant_gr, reference_gr, type = "within", ignore.strand = TRUE
  )
  if (!length(hits)) {
    return(data.frame(
      variant_anchor_id = character(), gene_id_normalized = character(),
      stringsAsFactors = FALSE
    ))
  }
  unique(data.frame(
    variant_anchor_id = anchors$variant_anchor_id[S4Vectors::queryHits(hits)],
    gene_id_normalized = tss$gene_id_normalized[S4Vectors::subjectHits(hits)],
    stringsAsFactors = FALSE
  ))
}

.reference_candidates <- function(
    variants,
    genes,
    target_gene_ids,
    radius,
    explicit_mapping = NULL) {
  target_gene_ids <- unique(as.character(target_gene_ids))
  tss <- genes$tss[genes$tss$gene_id_normalized %in% target_gene_ids, , drop = FALSE]
  anchors <- variants$anchors
  mapping <- if (is.null(explicit_mapping)) {
    .distance_reference_pairs(anchors, tss, radius)
  } else {
    .resolve_reference_mapping(explicit_mapping, variants, genes)
  }
  mapping <- mapping[mapping$gene_id_normalized %in% target_gene_ids, , drop = FALSE]
  pairs <- merge(mapping, anchors, by = "variant_anchor_id", all.x = TRUE,
                 sort = FALSE)
  gene_meta <- genes$genes[c(
    "gene_id_normalized", "gene_id_raw", "gene_symbol", "gene_strand"
  )]
  pairs <- merge(pairs, gene_meta, by = "gene_id_normalized", all.x = TRUE,
                 sort = FALSE)
  pairs$gene_id <- pairs$gene_id_normalized
  pairs$pair_id <- .stable_row_id(
    pairs$variant_anchor_id, pairs$gene_id_normalized, prefix = "reference_pair"
  )
  pairs$reference_pair <- TRUE
  pairs <- pairs[c(
    "pair_id", "variant_anchor_id", "variant_key", "assembly", "chrom",
    "pos1", "anchor0", "gene_id", "gene_id_normalized", "gene_id_raw",
    "gene_symbol", "gene_strand", "reference_pair"
  )]
  pair_tss <- merge(
    pairs[c(
      "pair_id", "variant_anchor_id", "variant_key", "chrom", "pos1",
      "anchor0", "gene_id", "gene_id_normalized"
    )],
    tss,
    by = c("gene_id", "gene_id_normalized"),
    all = FALSE,
    sort = FALSE
  )
  structure(
    list(pairs = pairs, tss = pair_tss, mapping = mapping,
         variants = variants, genes = genes,
         target_gene_ids = target_gene_ids),
    class = "v2g_reference_candidates"
  )
}

.profile_row_from_evidence <- function(profile, evidence) {
  qc <- evidence$contact_qc
  activity_liftover <- evidence$activity_liftover_qc
  activity_source_build <- unique(activity_liftover$source_build)
  activity_target_build <- unique(activity_liftover$target_build)
  activity_operation <- if (length(activity_source_build) &&
                            length(activity_target_build) &&
                            identical(activity_source_build[[1L]],
                                      activity_target_build[[1L]])) {
    "identity"
  } else {
    "liftover"
  }
  data.frame(
    profile_id = profile$profile_id,
    activity_source = profile$activity$path,
    contact_source = profile$contact$path,
    activity_type = profile$activity$type %||% "auto",
    activity_query_source_build = activity_source_build[[1L]],
    activity_query_target_build = activity_target_build[[1L]],
    activity_query_operation = activity_operation,
    activity_query_chain_direction = activity_liftover$chain_direction[[1L]],
    contact_type = profile$contact$type %||% "auto",
    contact_value_requested = qc$contact_value_requested[[1L]],
    contact_value_resolved = qc$contact_value_resolved[[1L]],
    contact_layer_resolved = qc$contact_layer_resolved[[1L]],
    contact_resolution_resolved = qc$contact_resolution_resolved[[1L]],
    stringsAsFactors = FALSE
  )
}

.score_profile_gene_wise <- function(
    candidates,
    reference_candidates,
    profile,
    config,
    chrom_lengths,
    explicit_mapping,
    reference_set_id) {
  target <- .evidence_for_profile(
    candidates, profile, config, chrom_lengths, evidence_scope = "target"
  )
  reference <- if (nrow(reference_candidates$pairs)) {
    .evidence_for_profile(
      reference_candidates, profile, config, chrom_lengths,
      evidence_scope = "reference"
    )
  } else {
    empty_pairs <- reference_candidates$pairs
    empty_pairs$profile_id <- character()
    empty_pairs$activity_raw <- numeric()
    empty_pairs$contact_raw <- numeric()
    empty_pairs$activity_transformed <- numeric()
    empty_pairs$contact_transformed <- numeric()
    empty_pairs$proxy_score <- numeric()
    list(
      pairs = empty_pairs,
      activity_qc = data.frame(),
      contact_qc = data.frame(
        profile_id = profile$profile_id,
        status = "no_reference_pairs",
        stringsAsFactors = FALSE
      ),
      tss_bin_contacts = data.frame(),
      activity_liftover_qc = data.frame()
    )
  }
  target$pairs$.target_row <- TRUE
  target$pairs$.target_order <- seq_len(nrow(target$pairs))
  reference$pairs$.target_row <- FALSE
  reference$pairs$.target_order <- NA_integer_
  reference_mode <- config$method$reference_mode
  gene_tss <- reference_candidates$genes$tss[c(
    "gene_id_normalized", "tss_chrom", "tss_pos1"
  )]

  if (reference_mode == "input_variants") {
    group_cols <- c("profile_id", "gene_id_normalized")
    if (nrow(reference$pairs)) {
      reference_normalized <- compute_gene_wise_z(
        reference$pairs,
        reference_mode = "input_variants",
        gene_tss = gene_tss,
        reference_gene_radius_bp = config$method$reference_gene_radius_bp,
        reference_gene_variant_map = explicit_mapping,
        reference_gene_variant_map_id = .reference_mapping_id(explicit_mapping),
        reference_set_id = reference_set_id,
        proxy_col = "proxy_score",
        group_cols = group_cols,
        gene_col = "gene_id_normalized",
        variant_col = "variant_anchor_id",
        variant_chrom_col = "chrom",
        variant_pos_col = "pos1",
        gene_chrom_col = "tss_chrom",
        gene_tss_col = "tss_pos1",
        warn_missing_fraction = config$method$max_missing_fraction,
        z_col = "gene_wise_z"
      )
      qc <- reference_normalized$normalization_qc
    } else {
      qc <- .acv2g_empty_normalization_qc(target$pairs, group_cols)
    }
    qc <- .acv2g_add_absent_target_groups(
      qc,
      target$pairs,
      group_cols,
      method = "gene_wise_z",
      warn_missing_fraction = config$method$max_missing_fraction
    )
    universe_rule <- if (is.null(explicit_mapping)) "distance" else
      "explicit_mapping"
    mapping_id <- .reference_mapping_id(explicit_mapping)
    radius_record <- if (is.null(explicit_mapping))
      config$method$reference_gene_radius_bp else NA_real_
    defaults <- list(
      reference_mode = "input_variants",
      reference_universe_rule = universe_rule,
      reference_set_id = reference_set_id,
      reference_gene_radius_bp = radius_record,
      reference_gene_variant_map_id = mapping_id,
      reference_variant_count = qc$normalization_n_total,
      reference_pair_count = qc$normalization_n_total
    )
    for (column in names(defaults)) {
      value <- defaults[[column]]
      if (!column %in% names(qc)) {
        qc[[column]] <- rep_len(value, nrow(qc))
      } else {
        missing <- is.na(qc[[column]])
        qc[[column]][missing] <- rep_len(value, nrow(qc))[missing]
      }
    }
    qc$normalization_universe_id <- paste0(
      reference_set_id, "::", qc$normalization_group
    )
    scored <- .acv2g_attach_group_statistics(
      target$pairs,
      qc,
      group_cols,
      proxy_col = "proxy_score",
      z_col = "gene_wise_z",
      method = "gene_wise_z"
    )
    qc_index <- .acv2g_match_groups(scored, qc, group_cols)
    for (column in c(
      "reference_mode", "reference_universe_rule", "reference_set_id",
      "reference_gene_radius_bp", "reference_gene_variant_map_id"
    )) {
      scored[[column]] <- qc[[column]][qc_index]
    }
    normalized <- .acv2g_normalization_result(scored, qc)
  } else {
    normalized <- compute_gene_wise_z(
      target$pairs,
      reference_mode = "background_variants",
      background_pairs = reference$pairs,
      gene_tss = gene_tss,
      reference_gene_radius_bp = config$method$reference_gene_radius_bp,
      reference_gene_variant_map = explicit_mapping,
      reference_gene_variant_map_id = .reference_mapping_id(explicit_mapping),
      reference_set_id = reference_set_id,
      proxy_col = "proxy_score",
      group_cols = c("profile_id", "gene_id_normalized"),
      gene_col = "gene_id_normalized",
      variant_col = "variant_anchor_id",
      variant_chrom_col = "chrom",
      variant_pos_col = "pos1",
      gene_chrom_col = "tss_chrom",
      gene_tss_col = "tss_pos1",
      warn_missing_fraction = config$method$max_missing_fraction,
      z_col = "gene_wise_z"
    )
    scored <- normalized$pairs
  }
  scored$.target_row <- NULL
  scored$.target_order <- NULL
  scored <- call_v2g_links(
    scored,
    call_rule = config$method$call_rule,
    z_col = "gene_wise_z",
    z_threshold = if (config$method$call_rule == "z_threshold")
      config$method$z_threshold else NULL
  )
  target$activity_qc$evidence_scope <- "target"
  reference$activity_qc$evidence_scope <- rep(
    "reference", nrow(reference$activity_qc)
  )
  target$contact_qc$evidence_scope <- "target"
  reference$contact_qc$evidence_scope <- rep(
    "reference", nrow(reference$contact_qc)
  )
  target$tss_bin_contacts$evidence_scope <- "target"
  reference$tss_bin_contacts$evidence_scope <- rep(
    "reference", nrow(reference$tss_bin_contacts)
  )
  list(
    pairs = scored,
    normalization_qc = normalized$normalization_qc,
    activity_qc = data.table::rbindlist(
      list(target$activity_qc, reference$activity_qc), fill = TRUE
    ) |> as.data.frame(),
    contact_qc = data.table::rbindlist(
      list(target$contact_qc, reference$contact_qc), fill = TRUE
    ) |> as.data.frame(),
    tss_bin_contacts = data.table::rbindlist(
      list(target$tss_bin_contacts, reference$tss_bin_contacts), fill = TRUE
    ) |> as.data.frame(),
    liftover_qc = data.table::rbindlist(
      list(target$activity_liftover_qc,
           reference$activity_liftover_qc), fill = TRUE
    ) |> as.data.frame(),
    profile = .profile_row_from_evidence(profile, target)
  )
}

.predict_v2g_gene_wise <- function(config) {
  preparation <- .capture_warnings({
    structural <- .prepare_structural_inputs(config)
    if (.all_fixed_loci_failed(structural)) {
      .acv2g_warn(
        paste0(
          "All fixed loci failed structural liftover; no pair scoring was ",
          "attempted and a zero-pair partial result is being returned."
        ),
        class = "acv2g_all_loci_failed_warning"
      )
      list(structural = structural, all_loci_failed = TRUE)
    } else {
      candidates <- .build_candidates_from_config(structural, config)
      if (!nrow(candidates$pairs)) {
        .acv2g_validation_error(
          "The validated inputs produce no target variant-gene pairs."
        )
      }
      background <- NULL
      if (config$method$reference_mode == "input_variants") {
        reference_variants <- structural$variants
        background_liftover_qc <- NULL
      } else {
        background <- .read_and_map_background_variants(config)
        reference_variants <- background$value
        background_liftover_qc <- background$qc
      }
    mapping_input <- .read_coordinate_mapping(
      config$inputs$reference_gene_variant_map %||% NULL,
      expected_assembly = config$analysis_build,
      name = "inputs.reference_gene_variant_map"
      )
      resolved_mapping <- if (is.null(mapping_input)) NULL else
        .resolve_reference_mapping(
          mapping_input, reference_variants, structural$genes
        )
      reference_candidates <- .reference_candidates(
        reference_variants,
        structural$genes,
        unique(candidates$pairs$gene_id_normalized),
        config$method$reference_gene_radius_bp,
        explicit_mapping = resolved_mapping
      )
      scoring_mapping <- if (is.null(mapping_input)) {
        NULL
      } else {
        reference_candidates$mapping
      }
      if (!nrow(reference_candidates$pairs)) {
        .acv2g_warn(
          paste0(
            "The declared gene reference universe contains no eligible ",
            "variant-gene pairs; gene-wise Z-scores will be unavailable."
          )
        )
      }
      identity <- .reference_identity(
        reference_candidates, config, scoring_mapping
      )
      reference_set_id <- config$method$reference_set_id %||%
        identity$reference_set_id
      list(
        structural = structural,
        all_loci_failed = FALSE,
        candidates = candidates,
        background = background,
        background_liftover_qc = background_liftover_qc,
        reference_candidates = reference_candidates,
        scoring_mapping = scoring_mapping,
        reference_identity = identity,
        reference_set_id = reference_set_id
      )
    }
  })
  structural <- preparation$value$structural
  if (isTRUE(preparation$value$all_loci_failed)) {
    return(.all_failed_loci_result(
      config, structural,
      .combine_warning_tables(preparation$warnings)
    ))
  }
  candidates <- preparation$value$candidates
  background <- preparation$value$background
  background_liftover_qc <- preparation$value$background_liftover_qc
  reference_candidates <- preparation$value$reference_candidates
  scoring_mapping <- preparation$value$scoring_mapping
  reference_identity <- preparation$value$reference_identity
  reference_set_id <- preparation$value$reference_set_id
  chrom_lengths <- .read_chrom_lengths(config$inputs$chrom_sizes %||% NULL)
  profile_runs <- lapply(config$profiles, function(profile) {
    .profile_with_warnings(
      .score_profile_gene_wise(
        candidates, reference_candidates, profile, config, chrom_lengths,
        scoring_mapping, reference_set_id
      ),
      profile$profile_id
    )
  })
  scored <- lapply(profile_runs, `[[`, "value")
  warnings <- .combine_warning_tables(
    preparation$warnings,
    .combine_tables(profile_runs, "warnings")
  )
  pairs <- .combine_tables(scored, "pairs")
  profiles <- .combine_tables(scored, "profile")
  annotation_qc <- structural$genes$annotation_qc
  annotation_qc$status <- annotation_qc$annotation_status
  liftover_qc <- structural$liftover_qc
  liftover_qc <- data.table::rbindlist(
    list(
      liftover_qc,
      background_liftover_qc,
      .combine_tables(scored, "liftover_qc")
    ),
    fill = TRUE
  ) |> as.data.frame()
  called_pairs <- if (config$method$call_rule == "z_threshold") {
    pairs[pairs$link_call %in% TRUE, , drop = FALSE]
  } else NULL
  record_map <- .result_variant_record_map(
    target_variants = structural$variants,
    candidates = candidates,
    reference_candidates = reference_candidates,
    background_variants = if (is.null(background)) NULL else background$value,
    reference_mode = config$method$reference_mode
  )
  new_v2g_result(
    pairs = pairs,
    profiles = profiles,
    tss_bin_contacts = .combine_tables(scored, "tss_bin_contacts"),
    variant_record_anchor_map = record_map,
    variant_locus_qc = candidates$qc$variant_locus,
    locus_status = .candidate_locus_status(
      candidates, names(config$profiles), structural$failed_loci,
      structural$all_locus_ids
    ),
    normalization_qc = .combine_tables(scored, "normalization_qc"),
    activity_qc = .combine_tables(scored, "activity_qc"),
    contact_qc = .combine_tables(scored, "contact_qc"),
    annotation_qc = annotation_qc,
    liftover_qc = liftover_qc,
    warnings = warnings,
    config = config,
    manifest = build_v2g_manifest(config, profiles, extra = list(
      reference_set_id = reference_set_id,
      reference_universe_identity = reference_identity$reference_set_id,
      reference_gene_variant_map_id = reference_identity$mapping_id,
      reference_universe_rule = if (is.null(scoring_mapping)) {
        "distance"
      } else {
        "explicit_mapping"
      },
      reference_variant_count = length(unique(
        reference_candidates$pairs$variant_anchor_id
      )),
      reference_pair_count = nrow(reference_candidates$pairs)
    )),
    run_status = list(
      status = if (length(structural$failed_loci)) "partial" else "success",
      failed_loci = structural$failed_loci,
      started_at_utc = .utc_now(),
      completed_at_utc = .utc_now()
    ),
    called_pairs = called_pairs
  )
}

#' Run the ACV2G workflow from a resolved configuration
#'
#' @param config Configuration list or YAML path.
#'
#' @return An in-memory `v2g_result`.
#' @export
predict_v2g <- function(config) {
  config <- validate_v2g_config(config)
  if (config$method$name != "locus_pair_z") {
    return(.predict_v2g_gene_wise(config))
  }
  preparation <- .capture_warnings({
    structural <- .prepare_structural_inputs(config)
    if (.all_fixed_loci_failed(structural)) {
      .acv2g_warn(
        paste0(
          "All fixed loci failed structural liftover; no pair scoring was ",
          "attempted and a zero-pair partial result is being returned."
        ),
        class = "acv2g_all_loci_failed_warning"
      )
      list(structural = structural, all_loci_failed = TRUE)
    } else {
      candidates <- .build_candidates_from_config(structural, config)
      if (!nrow(candidates$pairs)) {
        .acv2g_validation_error(
          "The validated inputs produce no eligible variant-gene pairs."
        )
      }
      list(
        structural = structural,
        candidates = candidates,
        all_loci_failed = FALSE
      )
    }
  })
  structural <- preparation$value$structural
  if (isTRUE(preparation$value$all_loci_failed)) {
    return(.all_failed_loci_result(
      config, structural,
      .combine_warning_tables(preparation$warnings)
    ))
  }
  candidates <- preparation$value$candidates
  chrom_lengths <- .read_chrom_lengths(config$inputs$chrom_sizes %||% NULL)
  profile_runs <- lapply(config$profiles, function(profile) {
    .profile_with_warnings(
      .score_profile_locus(candidates, profile, config, chrom_lengths),
      profile$profile_id
    )
  })
  scored <- lapply(profile_runs, `[[`, "value")
  warnings <- .combine_warning_tables(
    preparation$warnings,
    .combine_tables(profile_runs, "warnings")
  )
  pairs <- .combine_tables(scored, "pairs")
  profiles <- .combine_tables(scored, "profile")
  annotation_qc <- structural$genes$annotation_qc
  annotation_qc$status <- annotation_qc$annotation_status
  locus_status <- .candidate_locus_status(
    candidates, names(config$profiles), structural$failed_loci,
    structural$all_locus_ids
  )
  run_status <- list(
    status = if (length(structural$failed_loci)) "partial" else "success",
    failed_loci = structural$failed_loci,
    started_at_utc = .utc_now(),
    completed_at_utc = .utc_now()
  )
  manifest <- build_v2g_manifest(config, profiles)
  called_pairs <- if (config$method$call_rule == "z_threshold") {
    pairs[pairs$link_call %in% TRUE, , drop = FALSE]
  } else {
    NULL
  }
  new_v2g_result(
    pairs = pairs,
    profiles = profiles,
    tss_bin_contacts = .combine_tables(scored, "tss_bin_contacts"),
    variant_record_anchor_map = .result_variant_record_map(
      target_variants = structural$variants,
      candidates = candidates
    ),
    variant_locus_qc = candidates$qc$variant_locus,
    locus_status = locus_status,
    normalization_qc = .combine_tables(scored, "normalization_qc"),
    activity_qc = .combine_tables(scored, "activity_qc"),
    contact_qc = .combine_tables(scored, "contact_qc"),
    annotation_qc = annotation_qc,
    liftover_qc = data.table::rbindlist(
      list(structural$liftover_qc, .combine_tables(scored, "liftover_qc")),
      fill = TRUE
    ) |> as.data.frame(),
    warnings = warnings,
    config = config,
    manifest = manifest,
    run_status = run_status,
    called_pairs = called_pairs
  )
}

.write_failed_validation_bundle <- function(condition, config, path, overwrite) {
  target <- .validate_output_target(path)
  if (file.exists(target) && !isTRUE(overwrite)) return(invisible(NULL))
  temp_dir <- tempfile(paste0(".", basename(target), ".failed-"),
                       tmpdir = dirname(target))
  dir.create(temp_dir, recursive = FALSE)
  report <- data.frame(
    schema_version = .acv2g_schema_version,
    class = paste(class(condition), collapse = ";"),
    message = conditionMessage(condition),
    timestamp_utc = .utc_now(),
    stringsAsFactors = FALSE
  )
  .write_tsv_gz(report, file.path(temp_dir, "validation_report.tsv.gz"))
  jsonlite::write_json(
    list(schema_version = .acv2g_schema_version, status = "failed",
         error = conditionMessage(condition), updated_at_utc = .utc_now()),
    file.path(temp_dir, "run_status.json"), auto_unbox = TRUE, pretty = TRUE
  )
  yaml::write_yaml(config %||% list(), file.path(temp_dir, "config.resolved.yaml"))
  jsonlite::write_json(
    list(schema_version = .acv2g_schema_version, package = "ACV2G",
         created_at_utc = .utc_now(), status = "failed_validation"),
    file.path(temp_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE
  )
  files <- list.files(temp_dir, full.names = TRUE)
  relative <- basename(files)
  hashes <- vapply(files, digest::digest, character(1L), file = TRUE,
                   algo = "sha256", serialize = FALSE)
  writeLines(paste0(hashes, "  ", relative),
             file.path(temp_dir, "checksums.sha256"))
  if (file.exists(target)) unlink(target, recursive = TRUE, force = TRUE)
  file.rename(temp_dir, target)
  invisible(target)
}

#' Validate, run, and optionally write ACV2G from configuration
#'
#' @param config Configuration list or YAML path.
#' @param output Optional result-bundle directory.
#' @param overwrite Explicitly replace an existing output.
#'
#' @return A `v2g_result` when `output` is `NULL`; otherwise the normalized
#'   output path, invisibly.
#' @export
predict_v2g_from_config <- function(config, output = NULL, overwrite = FALSE) {
  raw_config <- config
  result <- tryCatch(
    predict_v2g(config),
    acv2g_validation_error = function(e) {
      if (!is.null(output)) {
        .write_failed_validation_bundle(e, raw_config, output, overwrite)
      }
      stop(e)
    }
  )
  if (is.null(output)) return(result)
  write_v2g_result(result, output, overwrite = overwrite)
}
