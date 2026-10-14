# Locus assignment and candidate universes --------------------------------

.as_v2g_loci <- function(x) {
  x <- .as_plain_data_frame(x, "loci")
  if (all(c(
    "locus_id", "locus_assembly", "locus_chrom", "locus_start0",
    "locus_end0"
  ) %in% names(x))) {
    invalid_id <- is.na(x$locus_id) |
      !nzchar(trimws(as.character(x$locus_id)))
    if (any(invalid_id) || anyDuplicated(x$locus_id)) {
      .acv2g_validation_error(
        "Standardized loci require unique, non-missing locus_id values."
      )
    }
    x$locus_chrom <- .validate_chrom(x$locus_chrom, "loci$locus_chrom")
    intervals <- .validate_half_open(
      x$locus_start0, x$locus_end0, "standardized loci"
    )
    x$locus_start0 <- intervals$start0
    x$locus_end0 <- intervals$end0
    invalid_assembly <- is.na(x$locus_assembly) |
      !nzchar(trimws(as.character(x$locus_assembly)))
    if (any(invalid_assembly)) {
      .acv2g_validation_error(
        "Standardized loci require an explicit locus_assembly on every row."
      )
    }
    return(x)
  }
  .acv2g_validation_error(
    "loci must be a standardized table returned by read_v2g_loci()."
  )
}

.resolve_explicit_membership_anchors <- function(mapping, variants, anchors) {
  if ("variant_anchor_id" %in% names(mapping)) {
    mapping$variant_anchor_id <- as.character(mapping$variant_anchor_id)
    return(mapping)
  }
  if ("variant_key" %in% names(mapping)) {
    key_to_anchor <- anchors[c("variant_key", "variant_anchor_id")]
    if (anyDuplicated(key_to_anchor$variant_key)) {
      .acv2g_validation_error("Internal variant keys are not unique.")
    }
    mapping <- merge(mapping, key_to_anchor, by = "variant_key", all.x = TRUE, sort = FALSE)
    return(mapping)
  }
  if (all(c("chrom", "pos1") %in% names(mapping))) {
    mapping$chrom <- .validate_chrom(mapping$chrom, "variant_locus_map$chrom")
    mapping$pos1 <- .as_integer_coordinate(
      mapping$pos1, "variant_locus_map$pos1", 1
    )
    assemblies <- unique(anchors$assembly)
    if ("assembly" %in% names(mapping)) {
      map_assembly <- as.character(mapping$assembly)
    } else if (length(assemblies) == 1L) {
      map_assembly <- rep(assemblies, nrow(mapping))
    } else {
      .acv2g_validation_error(
        paste0(
          "variant_locus_map must include assembly when variants span ",
          "multiple assemblies."
        )
      )
    }
    mapping$variant_key <- paste(
      map_assembly, mapping$chrom, mapping$pos1, sep = ":"
    )
    key_to_anchor <- anchors[c("variant_key", "variant_anchor_id")]
    mapping <- merge(mapping, key_to_anchor, by = "variant_key", all.x = TRUE, sort = FALSE)
    return(mapping)
  }

  record_field <- intersect(
    c("record_id", "submitted_record_id", "variant_id"), names(mapping)
  )
  if (length(record_field) && inherits(variants, "v2g_variants")) {
    field <- record_field[[1L]]
    invalid_identifier <- is.na(mapping[[field]]) |
      !nzchar(trimws(as.character(mapping[[field]])))
    if (any(invalid_identifier)) {
      .acv2g_validation_error(
        paste0(
          "Explicit variant-locus map contains ", sum(invalid_identifier),
          " missing ", field, " value(s)."
        )
      )
    }
    record_map <- variants$record_map[c(field, "variant_anchor_id")]
    record_map <- unique(record_map)
    ambiguous <- names(which(table(record_map[[field]]) > 1L))
    if (length(ambiguous)) {
      .acv2g_validation_error(
        paste0(
          "Explicit membership identifiers resolve to multiple anchors: ",
          paste(head(ambiguous, 10L), collapse = ", "),
          ". Use variant_anchor_id."
        )
      )
    }
    mapping <- merge(mapping, record_map, by = field, all.x = TRUE, sort = FALSE)
    return(mapping)
  }

  .acv2g_validation_error(
    paste0(
      "variant_locus_map must identify variants by variant_anchor_id, ",
      "variant_key, assembly/chrom/pos1, or an unambiguous submitted-record ",
      "identifier."
    )
  )
}

.coordinate_memberships <- function(loci, anchors) {
  locus_ranges <- GenomicRanges::GRanges(
    seqnames = loci$locus_chrom,
    ranges = IRanges::IRanges(
      start = as.numeric(loci$locus_start0) + 1,
      end = as.numeric(loci$locus_end0)
    )
  )
  anchor_ranges <- GenomicRanges::GRanges(
    seqnames = anchors$chrom,
    ranges = IRanges::IRanges(start = as.numeric(anchors$pos1), width = 1L)
  )
  hits <- GenomicRanges::findOverlaps(
    anchor_ranges, locus_ranges, type = "within", ignore.strand = TRUE
  )
  if (!length(hits)) {
    return(data.frame(
      variant_anchor_id = character(), locus_id = character(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    variant_anchor_id = anchors$variant_anchor_id[S4Vectors::queryHits(hits)],
    locus_id = loci$locus_id[S4Vectors::subjectHits(hits)],
    stringsAsFactors = FALSE
  )
}

.membership_qc <- function(anchors, memberships, mode) {
  membership_count <- table(memberships$variant_anchor_id)
  count <- as.integer(membership_count[anchors$variant_anchor_id])
  count[is.na(count)] <- 0L
  status <- if (identical(mode, "coordinate")) {
    ifelse(count == 0L, "outside_all_loci", "assigned_unique_locus")
  } else {
    ifelse(count == 0L, "not_in_explicit_mapping", "explicit_membership")
  }
  data.frame(
    anchors[c("variant_anchor_id", "variant_key", "assembly", "chrom", "pos1")],
    variant_locus_status = status,
    locus_membership_count = count,
    locus_assignment_mode = mode,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Assign unique coordinate anchors to fixed loci
#'
#' Without an explicit map, an anchor must overlap exactly one locus. Multiple
#' overlaps are an error; anchors outside all loci are retained in QC and
#' excluded with a warning. An explicit map may assign one anchor to several
#' loci but every membership must lie inside each unchanged locus boundary.
#'
#' @param loci Standardized loci from [read_v2g_loci()].
#' @param variants A `v2g_variants` object or standardized anchor table.
#' @param variant_locus_map Optional explicit membership data frame.
#' @return A list with validated `memberships` and per-anchor `qc`.
#' @export
assign_v2g_loci <- function(loci, variants, variant_locus_map = NULL) {
  loci <- .as_v2g_loci(loci)
  anchors <- .as_variant_anchors(variants)
  if (!nrow(anchors)) {
    .acv2g_validation_error("No variant anchors are available for locus assignment.")
  }

  locus_assemblies <- unique(loci$locus_assembly)
  variant_assemblies <- unique(anchors$assembly)
  if (length(locus_assemblies) != 1L || length(variant_assemblies) != 1L ||
      !identical(locus_assemblies, variant_assemblies)) {
    .acv2g_validation_error(
      paste0(
        "Loci and variants must be in the same single analysis assembly before ",
        "candidate construction; observed loci ",
        paste(locus_assemblies, collapse = ", "), " and variants ",
        paste(variant_assemblies, collapse = ", "), "."
      )
    )
  }

  if (is.null(variant_locus_map)) {
    memberships <- .coordinate_memberships(loci, anchors)
    counts <- table(memberships$variant_anchor_id)
    ambiguous_ids <- names(counts[counts > 1L])
    if (length(ambiguous_ids)) {
      details <- vapply(ambiguous_ids, function(anchor_id) {
        locus_ids <- memberships$locus_id[memberships$variant_anchor_id == anchor_id]
        key <- anchors$variant_key[match(anchor_id, anchors$variant_anchor_id)]
        paste0(key, " -> ", paste(sort(locus_ids), collapse = ","))
      }, character(1L))
      .acv2g_validation_error(
        paste0(
          "Variant anchor(s) overlap multiple fixed loci and require an ",
          "explicit map: ", paste(head(details, 10L), collapse = "; "), "."
        )
      )
    }
    mode <- "coordinate"
  } else {
    mapping <- .as_plain_data_frame(variant_locus_map, "variant_locus_map")
    .require_columns(mapping, "locus_id", "variant_locus_map")
    mapping$locus_id <- as.character(mapping$locus_id)
    mapping <- .resolve_explicit_membership_anchors(mapping, variants, anchors)
    unknown_locus <- is.na(mapping$locus_id) | !mapping$locus_id %in% loci$locus_id
    unknown_anchor <- is.na(mapping$variant_anchor_id) |
      !mapping$variant_anchor_id %in% anchors$variant_anchor_id
    if (any(unknown_locus) || any(unknown_anchor)) {
      .acv2g_validation_error(
        paste0(
          "Explicit variant-locus map contains ", sum(unknown_locus),
          " unknown locus reference(s) and ", sum(unknown_anchor),
          " unknown variant reference(s)."
        )
      )
    }
    memberships <- unique(mapping[c("variant_anchor_id", "locus_id")])
    membership_detail <- merge(
      merge(memberships, anchors, by = "variant_anchor_id", all.x = TRUE),
      loci,
      by = "locus_id",
      all.x = TRUE
    )
    inside <- membership_detail$chrom == membership_detail$locus_chrom &
      membership_detail$anchor0 >= membership_detail$locus_start0 &
      membership_detail$anchor0 < membership_detail$locus_end0
    if (any(!inside)) {
      bad <- paste0(
        membership_detail$variant_key[!inside],
        " -> ", membership_detail$locus_id[!inside]
      )
      .acv2g_validation_error(
        paste0(
          "Explicit membership places a variant outside its fixed locus: ",
          paste(head(bad, 10L), collapse = "; "), "."
        )
      )
    }
    mode <- "explicit_mapping"
  }

  memberships <- unique(memberships[c("variant_anchor_id", "locus_id")])
  memberships <- merge(
    merge(memberships, anchors, by = "variant_anchor_id", all.x = TRUE, sort = FALSE),
    loci,
    by = "locus_id",
    all.x = TRUE,
    sort = FALSE
  )
  memberships$variant_locus_status <- if (identical(mode, "coordinate")) {
    "assigned_unique_locus"
  } else {
    "explicit_membership"
  }
  memberships$locus_assignment_mode <- mode
  memberships <- memberships[order(
    memberships$locus_chrom, memberships$locus_start0,
    memberships$anchor0, memberships$locus_id
  ), , drop = FALSE]
  rownames(memberships) <- NULL

  qc <- .membership_qc(anchors, memberships, mode)
  excluded <- qc$locus_membership_count == 0L
  if (identical(mode, "coordinate") && any(excluded)) {
    .acv2g_warn(
      paste0(
        "Excluded ", sum(excluded),
        " variant anchor(s) outside all fixed loci; details are retained ",
        "in variant-locus QC."
      ),
      class = "acv2g_outside_loci_warning"
    )
  }
  if (!nrow(memberships)) {
    .acv2g_validation_error(
      "No variant anchors remain after fixed-locus assignment."
    )
  }
  list(memberships = memberships, qc = qc, mode = mode)
}

.minimal_gene_table <- function(tss) {
  split_tss <- split(tss, tss$gene_id_normalized)
  if (!length(split_tss)) {
    return(data.frame(
      gene_id = character(), gene_id_normalized = character(),
      gene_id_raw = character(), gene_symbol = character(),
      gene_strand = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(split_tss, function(group) {
    data.frame(
      gene_id = group$gene_id_normalized[[1L]],
      gene_id_normalized = group$gene_id_normalized[[1L]],
      gene_id_raw = if ("gene_id_raw" %in% names(group)) {
        sort(unique(group$gene_id_raw))[[1L]]
      } else {
        group$gene_id_normalized[[1L]]
      },
      gene_symbol = if ("gene_symbol" %in% names(group)) {
        .collapse_values(group$gene_symbol)
      } else {
        NA_character_
      },
      gene_strand = unique(group$tss_strand)[[1L]],
      stringsAsFactors = FALSE
    )
  }))
}

.empty_candidate_pairs <- function() {
  data.frame(
    schema_version = character(), pair_id = character(), locus_id = character(),
    locus_assembly = character(), locus_chrom = character(),
    locus_start0 = numeric(), locus_end0 = numeric(),
    variant_anchor_id = character(), variant_key = character(),
    assembly = character(), chrom = character(), pos1 = numeric(), anchor0 = numeric(),
    anchor_multiplicity = integer(), variant_locus_status = character(),
    locus_assignment_mode = character(), gene_id = character(),
    gene_id_normalized = character(), gene_id_raw = character(),
    gene_symbol = character(), gene_strand = character(),
    candidate_method = character(), gene_flank_bp = numeric(),
    candidate_region_start0 = numeric(), candidate_region_end0 = numeric(),
    eligible_tss_count = integer(), stringsAsFactors = FALSE
  )
}

.empty_candidate_tss <- function() {
  data.frame(
    pair_id = character(), locus_id = character(), variant_anchor_id = character(),
    variant_key = character(), chrom = character(), pos1 = numeric(), anchor0 = numeric(),
    gene_id = character(), gene_id_normalized = character(),
    gene_id_raw = character(), gene_id_raw_values = character(),
    transcript_id = character(), transcript_ids = character(),
    gene_symbol = character(), annotation_version = character(),
    assembly = character(), tss_id = character(),
    tss_chrom = character(), tss_pos1 = numeric(), tss_anchor0 = numeric(),
    tss_strand = character(), distance_to_tss_bp = numeric(),
    signed_distance_to_tss_bp = numeric(), gene_flank_bp = numeric(),
    candidate_region_start0 = numeric(), candidate_region_end0 = numeric(),
    exact_tss_record_count = integer(), annotation_input_convention = character(),
    stringsAsFactors = FALSE
  )
}

#' Build the complete locus anchor-by-gene candidate universe
#'
#' Variants remain restricted to original fixed loci. Candidate genes are the
#' protein-coding stable genes with at least one TSS in the locus expanded by
#' `gene_flank_bp`; only TSSs in that expanded interval enter that locus's
#' contact universe. No activity or contact values are inspected here.
#'
#' @param loci Standardized loci from [read_v2g_loci()].
#' @param variants Variants from [read_v2g_variants()] or a standardized anchor
#'   table.
#' @param genes Genes from [read_v2g_genes()] or a standardized eligible TSS
#'   table.
#' @param gene_flank_bp Non-negative integer flank applied only to the TSS
#'   search interval. It never changes locus boundaries or variant membership.
#' @param variant_locus_map Optional explicit many-to-many membership map.
#' @return A `v2g_candidates` object containing pair and pair-TSS tables,
#'   memberships, submitted-record mapping, and structural QC.
#' @export
build_v2g_candidates <- function(
    loci,
    variants,
    genes,
    gene_flank_bp = 1e6,
    variant_locus_map = NULL) {
  gene_flank_bp <- .assert_scalar_number(
    gene_flank_bp,
    arg = "gene_flank_bp",
    lower = 0,
    integer = TRUE
  )
  loci <- .as_v2g_loci(loci)
  anchors <- .as_variant_anchors(variants)
  tss <- .as_gene_tss(genes)
  gene_table <- if (inherits(genes, "v2g_genes")) {
    genes$genes[genes$genes$is_protein_coding, , drop = FALSE]
  } else {
    .minimal_gene_table(tss)
  }
  assignment <- assign_v2g_loci(loci, variants, variant_locus_map)
  memberships <- assignment$memberships

  assemblies <- list(
    loci = unique(loci$locus_assembly),
    variants = unique(anchors$assembly),
    genes = unique(tss$assembly)
  )
  assembly_values <- unique(unlist(assemblies, use.names = FALSE))
  if (length(assembly_values) != 1L) {
    .acv2g_validation_error(
      paste0(
        "Loci, variants, and genes must share one analysis assembly before ",
        "candidate construction; observed ",
        paste(assembly_values, collapse = ", "), "."
      )
    )
  }

  expanded <- loci
  expanded$candidate_region_start0 <- pmax(
    0, as.numeric(expanded$locus_start0) - gene_flank_bp
  )
  expanded$candidate_region_end0 <- as.numeric(expanded$locus_end0) + gene_flank_bp
  expanded$gene_flank_bp <- gene_flank_bp

  expanded_ranges <- GenomicRanges::GRanges(
    seqnames = expanded$locus_chrom,
    ranges = IRanges::IRanges(
      start = expanded$candidate_region_start0 + 1,
      end = expanded$candidate_region_end0
    )
  )
  tss_ranges <- GenomicRanges::GRanges(
    seqnames = tss$tss_chrom,
    ranges = IRanges::IRanges(start = tss$tss_pos1, width = 1L)
  )
  candidate_hits <- GenomicRanges::findOverlaps(
    tss_ranges, expanded_ranges, type = "within", ignore.strand = TRUE
  )

  if (!length(candidate_hits)) {
    pairs <- .empty_candidate_pairs()
    pair_tss <- .empty_candidate_tss()
  } else {
    tss_index <- S4Vectors::queryHits(candidate_hits)
    locus_index <- S4Vectors::subjectHits(candidate_hits)
    locus_tss <- tss[tss_index, , drop = FALSE]
    locus_tss$locus_id <- expanded$locus_id[locus_index]
    locus_tss$candidate_region_start0 <-
      expanded$candidate_region_start0[locus_index]
    locus_tss$candidate_region_end0 <-
      expanded$candidate_region_end0[locus_index]
    locus_tss$gene_flank_bp <- gene_flank_bp
    rownames(locus_tss) <- NULL
    candidate_genes <- unique(locus_tss[c(
      "locus_id", "gene_id_normalized", "candidate_region_start0",
      "candidate_region_end0", "gene_flank_bp"
    )])
    tss_count <- aggregate(
      locus_tss$tss_id,
      locus_tss[c("locus_id", "gene_id_normalized")],
      function(value) length(unique(value))
    )
    names(tss_count)[[3L]] <- "eligible_tss_count"
    candidate_genes <- merge(
      candidate_genes, tss_count,
      by = c("locus_id", "gene_id_normalized"),
      all.x = TRUE,
      sort = FALSE
    )
    candidate_genes <- merge(
      candidate_genes,
      gene_table,
      by = "gene_id_normalized",
      all.x = TRUE,
      sort = FALSE
    )

    pairs <- merge(
      memberships,
      candidate_genes,
      by = "locus_id",
      all = FALSE,
      sort = FALSE
    )
    if (!"assembly" %in% names(pairs)) {
      variant_assembly_column <- intersect(
        c("assembly.x", "assembly"), names(pairs)
      )
      if (length(variant_assembly_column)) {
        pairs$assembly <- pairs[[variant_assembly_column[[1L]]]]
      }
    }
    pairs$gene_id <- pairs$gene_id_normalized
    pairs$pair_id <- .stable_row_id(
      pairs$locus_id, pairs$variant_anchor_id, pairs$gene_id_normalized,
      prefix = "pair"
    )
    pairs$schema_version <- .acv2g_schema_version
    pairs$candidate_method <- "locus_tss_flank"
    required_pair_columns <- c(
      "schema_version", "pair_id", "locus_id", "locus_assembly",
      "locus_chrom", "locus_start0", "locus_end0", "variant_anchor_id",
      "variant_key", "assembly", "chrom", "pos1", "anchor0",
      "anchor_multiplicity", "variant_locus_status", "locus_assignment_mode",
      "gene_id", "gene_id_normalized", "gene_id_raw",
      "gene_symbol", "gene_strand", "candidate_method", "gene_flank_bp",
      "candidate_region_start0", "candidate_region_end0", "eligible_tss_count"
    )
    for (column in setdiff(required_pair_columns, names(pairs))) {
      pairs[[column]] <- NA
    }
    pairs <- pairs[required_pair_columns]
    pairs <- pairs[order(
      pairs$locus_chrom, pairs$locus_start0, pairs$anchor0,
      pairs$gene_id_normalized
    ), , drop = FALSE]
    rownames(pairs) <- NULL

    pair_keys <- pairs[c(
      "pair_id", "locus_id", "variant_anchor_id", "variant_key", "chrom",
      "pos1", "anchor0", "gene_id", "gene_id_normalized"
    )]
    pair_tss <- merge(
      pair_keys,
      locus_tss,
      by = c("locus_id", "gene_id_normalized"),
      all = FALSE,
      sort = FALSE
    )
    pair_tss$gene_id <- pair_tss$gene_id_normalized
    pair_tss$distance_to_tss_bp <- abs(pair_tss$pos1 - pair_tss$tss_pos1)
    pair_tss$signed_distance_to_tss_bp <- pair_tss$pos1 - pair_tss$tss_pos1
    required_tss_columns <- names(.empty_candidate_tss())
    for (column in setdiff(required_tss_columns, names(pair_tss))) {
      pair_tss[[column]] <- NA
    }
    pair_tss <- pair_tss[required_tss_columns]
    pair_tss <- pair_tss[order(
      pair_tss$locus_id, pair_tss$anchor0, pair_tss$gene_id_normalized,
      pair_tss$tss_anchor0
    ), , drop = FALSE]
    rownames(pair_tss) <- NULL
  }

  locus_candidate_count <- if (nrow(pairs)) {
    pair_gene <- unique(pairs[c("locus_id", "gene_id_normalized")])
    pair_count <- table(pairs$locus_id)
    gene_count <- table(pair_gene$locus_id)
    data.frame(
      locus_id = loci$locus_id,
      candidate_gene_count = as.integer(gene_count[loci$locus_id]),
      candidate_pair_count = as.integer(pair_count[loci$locus_id]),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      locus_id = loci$locus_id,
      candidate_gene_count = 0L,
      candidate_pair_count = 0L,
      stringsAsFactors = FALSE
    )
  }
  locus_candidate_count$candidate_gene_count[
    is.na(locus_candidate_count$candidate_gene_count)
  ] <- 0L
  locus_candidate_count$candidate_pair_count[
    is.na(locus_candidate_count$candidate_pair_count)
  ] <- 0L
  locus_candidate_count$candidate_status <- ifelse(
    locus_candidate_count$candidate_pair_count > 0L,
    "candidate_universe_built",
    "no_candidate_pairs"
  )

  structure(
    list(
      pairs = pairs,
      tss = pair_tss,
      memberships = memberships,
      record_map = if (inherits(variants, "v2g_variants")) {
        variants$record_map
      } else {
        NULL
      },
      qc = list(
        variant_locus = assignment$qc,
        loci = locus_candidate_count,
        annotation = if (inherits(genes, "v2g_genes")) genes$annotation_qc else NULL
      ),
      parameters = list(
        gene_flank_bp = gene_flank_bp,
        candidate_method = "locus_tss_flank",
        locus_assignment_mode = assignment$mode,
        analysis_assembly = assembly_values[[1L]]
      )
    ),
    class = "v2g_candidates"
  )
}
