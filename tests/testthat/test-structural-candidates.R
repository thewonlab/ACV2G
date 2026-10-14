test_that("coordinate primitives use exact 0-based half-open geometry", {
  expect_equal(ACV2G:::.pos1_to_anchor0(c(1, 10000, 10001)), c(0, 9999, 10000))
  expect_equal(ACV2G:::.anchor0_to_pos1(c(0, 9999)), c(1, 10000))

  bins <- ACV2G:::.fixed_contact_bins(c(10000, 10001), 10000)
  expect_equal(bins$bin_start0, c(0, 10000))
  expect_equal(bins$bin_end0, c(10000, 20000))

  even <- ACV2G:::.activity_window_geometry("chr1", 1001, 500)
  expect_equal(even$activity_window_requested_start0, 750)
  expect_equal(even$activity_window_requested_end0, 1250)
  expect_equal(even$activity_effective_window_bp, 500)

  odd <- ACV2G:::.activity_window_geometry("chr1", 1001, 501)
  expect_equal(odd$activity_window_requested_start0, 750)
  expect_equal(odd$activity_window_requested_end0, 1251)

  clipped <- ACV2G:::.activity_window_geometry(
    "chr1", 1, 500, chromosome_lengths = c(chr1 = 1000)
  )
  expect_equal(clipped$activity_window_requested_start0, -250)
  expect_equal(clipped$activity_window_effective_start0, 0)
  expect_equal(clipped$activity_window_effective_end0, 250)
  expect_true(clipped$activity_window_clipped)
  expect_error(
    ACV2G:::.activity_window_geometry(
      "chr1", 1001, 500, chromosome_lengths = c(chr1 = 1000)
    ),
    class = "acv2g_validation_error"
  )
})

test_that("loci are fixed, unique 0-based half-open intervals", {
  loci <- read_v2g_loci(
    data.frame(
      locus_id = c("L1", "L2"), chrom = c("chr1", "chr2"),
      start0 = c(0, 10), end0 = c(100, 20)
    ),
    assembly = "test1"
  )
  expect_s3_class(loci, "v2g_loci")
  expect_equal(loci$locus_start0, c(0, 10))
  expect_equal(loci$locus_end0, c(100, 20))
  expect_true(all(loci$locus_assembly == "test1"))

  expect_error(
    read_v2g_loci(
      data.frame(
        locus_id = c("same", "same"), chrom = "chr1",
        start0 = c(0, 20), end0 = c(10, 30)
      ),
      "test1"
    ),
    class = "acv2g_validation_error"
  )
  expect_error(
    read_v2g_loci(
      data.frame(locus_id = "L", chrom = "chr1", start0 = 10, end0 = 10),
      "test1"
    ),
    class = "acv2g_validation_error"
  )
})

test_that("BED loci require a fourth-column locus ID", {
  bed <- tempfile(fileext = ".bed")
  writeLines(c(
    "track name=synthetic",
    "chr1\t0\t10\tL1",
    "chr1\t20\t30\tL2"
  ), bed)
  loci <- read_v2g_loci(bed, "test1")
  expect_equal(loci$locus_id, c("L1", "L2"))

  invalid <- tempfile(fileext = ".bed")
  writeLines("chr1\t0\t10", invalid)
  expect_error(
    read_v2g_loci(invalid, "test1"),
    class = "acv2g_validation_error"
  )
})

test_that("variant records map many-to-one onto coordinate anchors", {
  variants <- read_v2g_variants(
    data.frame(
      submitted_record_id = c("r1", "r2", "r3"),
      chrom = c("chr1", "chr1", "chr1"),
      pos1 = c(101, 101, 202),
      ref = c("A", "A", "G"), alt = c("C", "G", "T"),
      variant_id = c("id1", "id2", "id3")
    ),
    assembly = "test1"
  )
  expect_s3_class(variants, "v2g_variants")
  expect_equal(nrow(variants$records), 3)
  expect_equal(nrow(variants$anchors), 2)
  expect_equal(sort(variants$anchors$anchor_multiplicity), c(1L, 2L))
  expect_equal(nrow(variants$record_map), 3)
  expect_equal(
    length(unique(variants$record_map$variant_anchor_id[1:2])),
    1
  )
  expect_match(variants$record_map$allele_key[[1]], "test1:chr1:101:A:C")
})

test_that("variant BED does not infer anchors for wide intervals", {
  bed <- tempfile(fileext = ".bed")
  writeLines("chr1\t10\t20\twide\t15", bed)
  expect_error(
    read_v2g_variants(bed, "test1"),
    class = "acv2g_validation_error"
  )
  anchored <- read_v2g_variants(bed, "test1", bed_anchor_column = 5)
  expect_equal(anchored$anchors$pos1, 15)
})

test_that("VCF alleles share one coordinate-level anchor", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t10\trs1\tA\tC,G\t.\tPASS\t."
  ), vcf)
  variants <- read_v2g_variants(vcf, "test1")
  expect_equal(nrow(variants$records), 2)
  expect_equal(nrow(variants$anchors), 1)
  expect_equal(variants$anchors$anchor_multiplicity, 2L)
  expect_equal(sort(variants$records$alt), c("C", "G"))
})

test_that("coordinate locus assignment rejects ambiguity and audits outsiders", {
  loci <- read_v2g_loci(
    data.frame(
      locus_id = c("L1", "L2"), chrom = "chr1",
      start0 = c(0, 40), end0 = c(100, 80)
    ),
    "test1"
  )
  variants <- read_v2g_variants(
    data.frame(chrom = "chr1", pos1 = c(51, 201)),
    "test1"
  )
  expect_error(
    assign_v2g_loci(loci, variants),
    class = "acv2g_validation_error"
  )

  anchor_id <- variants$anchors$variant_anchor_id[variants$anchors$pos1 == 51]
  mapped <- assign_v2g_loci(
    loci,
    variants,
    data.frame(
      variant_anchor_id = rep(anchor_id, 2),
      locus_id = c("L1", "L2")
    )
  )
  expect_equal(nrow(mapped$memberships), 2)
  expect_true(all(mapped$memberships$variant_locus_status == "explicit_membership"))

  nonoverlap_loci <- loci[loci$locus_id == "L1", , drop = FALSE]
  expect_warning(
    unique_assignment <- assign_v2g_loci(nonoverlap_loci, variants),
    class = "acv2g_outside_loci_warning"
  )
  outside <- unique_assignment$qc$pos1 == 201
  expect_equal(unique_assignment$qc$variant_locus_status[outside], "outside_all_loci")
})

test_that("explicit membership is validated against unchanged locus bounds", {
  loci <- read_v2g_loci(
    data.frame(locus_id = "L1", chrom = "chr1", start0 = 0, end0 = 100),
    "test1"
  )
  variants <- read_v2g_variants(
    data.frame(chrom = "chr1", pos1 = 101),
    "test1"
  )
  expect_error(
    assign_v2g_loci(
      loci,
      variants,
      data.frame(
        variant_anchor_id = variants$anchors$variant_anchor_id,
        locus_id = "L1"
      )
    ),
    class = "acv2g_validation_error"
  )
})

test_that("gene-level biotype is authoritative and transcript fallback is explicit", {
  genes <- read_v2g_genes(
    data.frame(
      gene_id = c("gene_nonpc", "gene_fallback"),
      transcript_id = c("tx1", "tx2"), chrom = "chr1",
      tss_pos1 = c(100, 200), strand = "+",
      gene_biotype = c("lncRNA", NA),
      transcript_biotype = c("protein_coding", "protein-coding")
    ),
    "test1"
  )
  nonpc <- genes$genes$gene_id == "gene_nonpc"
  fallback <- genes$genes$gene_id == "gene_fallback"
  expect_false(genes$genes$is_protein_coding[nonpc])
  expect_true(genes$genes$biotype_conflict[nonpc])
  expect_true(genes$genes$is_protein_coding[fallback])
  expect_equal(genes$genes$biotype_source[fallback], "transcript_fallback")
  expect_equal(genes$tss$gene_id_normalized, "gene_fallback")
})

test_that("unclassified genes require an explicit eligibility source", {
  annotation <- data.frame(
    gene_id = c("g1", "g2"), chrom = "chr1", tss_pos1 = c(10, 20),
    strand = "+"
  )
  expect_error(
    read_v2g_genes(annotation, "test1"),
    class = "acv2g_validation_error"
  )
  genes <- read_v2g_genes(
    annotation,
    "test1",
    protein_coding_gene_ids = "g1"
  )
  expect_equal(genes$genes$is_protein_coding, c(TRUE, FALSE))
  expect_equal(genes$tss$gene_id_normalized, "g1")

  annotation$is_protein_coding <- c(TRUE, FALSE)
  genes_by_flag <- read_v2g_genes(annotation, "test1")
  expect_equal(genes_by_flag$tss$gene_id_normalized, "g1")
})

test_that("terminal numeric ID versions strip only by explicit policy", {
  compatible <- data.frame(
    gene_id = c("ENSG1.1", "ENSG1.2"), transcript_id = c("t1", "t2"),
    chrom = "chr1", tss_pos1 = c(100, 200), strand = "+",
    gene_biotype = "protein_coding"
  )
  preserved <- read_v2g_genes(compatible, "test1")
  expect_equal(nrow(preserved$genes), 2)

  stripped <- read_v2g_genes(
    compatible,
    "test1",
    id_version_policy = "strip_terminal_numeric_version"
  )
  expect_equal(nrow(stripped$genes), 1)
  expect_equal(stripped$genes$gene_id_normalized, "ENSG1")
  expect_equal(stripped$genes$gene_id_collision_count, 2L)
  expect_equal(nrow(stripped$tss), 2)
  expect_equal(stripped$tss$gene_id_raw_values[[1]], "ENSG1.1")

  conflicting <- compatible
  conflicting$chrom[[2]] <- "chr2"
  expect_error(
    read_v2g_genes(
      conflicting,
      "test1",
      id_version_policy = "strip_terminal_numeric_version"
    ),
    class = "acv2g_validation_error"
  )
})

test_that("GTF reader retains transcript TSSs and strand-aware positions", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines(c(
    paste(
      "chr1", "synthetic", "gene", "100", "300", ".", "+", ".",
      'gene_id "g1"; gene_type "protein_coding"; gene_name "G1";',
      sep = "\t"
    ),
    paste(
      "chr1", "synthetic", "transcript", "100", "200", ".", "+", ".",
      'gene_id "g1"; transcript_id "tx1"; gene_type "protein_coding";',
      sep = "\t"
    ),
    paste(
      "chr1", "synthetic", "transcript", "150", "300", ".", "-", ".",
      'gene_id "g1"; transcript_id "tx2"; gene_type "protein_coding";',
      sep = "\t"
    )
  ), gtf)
  expect_error(
    read_v2g_genes(gtf, "test1"),
    class = "acv2g_validation_error"
  )

  # A stable gene cannot carry incompatible strands. Use separate genes to
  # verify the strand-aware TSS rule without permitting such a collision.
  writeLines(c(
    paste(
      "chr1", "synthetic", "transcript", "100", "200", ".", "+", ".",
      'gene_id "gplus"; transcript_id "tx1"; gene_type "protein_coding";',
      sep = "\t"
    ),
    paste(
      "chr1", "synthetic", "transcript", "150", "300", ".", "-", ".",
      'gene_id "gminus"; transcript_id "tx2"; gene_type "protein_coding";',
      sep = "\t"
    )
  ), gtf)
  genes <- read_v2g_genes(gtf, "test1")
  expect_equal(
    genes$tss$tss_pos1[match(c("gplus", "gminus"), genes$tss$gene_id)],
    c(100, 300)
  )
})

test_that("candidate construction is complete and locus centered", {
  loci <- read_v2g_loci(
    data.frame(locus_id = "L1", chrom = "chr1", start0 = 100, end0 = 200),
    "test1"
  )
  variants <- read_v2g_variants(
    data.frame(
      submitted_record_id = c("a1", "a2", "b", "outside"),
      chrom = "chr1", pos1 = c(101, 101, 200, 201)
    ),
    "test1"
  )
  genes <- read_v2g_genes(
    data.frame(
      gene_id = c("gA", "gA", "gA", "gB", "gC"),
      transcript_id = paste0("tx", 1:5), chrom = "chr1",
      # With a 50-bp flank, anchors 50 and 249 are included; 250 is excluded.
      tss_pos1 = c(51, 250, 251, 251, 150), strand = "+",
      gene_biotype = "protein_coding"
    ),
    "test1"
  )
  expect_warning(
    candidates <- build_v2g_candidates(
      loci, variants, genes, gene_flank_bp = 50
    ),
    class = "acv2g_outside_loci_warning"
  )
  expect_s3_class(candidates, "v2g_candidates")
  expect_equal(nrow(candidates$memberships), 2)
  expect_equal(nrow(candidates$pairs), 4)
  expect_equal(sort(unique(candidates$pairs$gene_id)), c("gA", "gC"))
  expect_equal(nrow(candidates$tss), 6)
  expect_equal(
    sort(unique(candidates$tss$tss_pos1[candidates$tss$gene_id == "gA"])),
    c(51, 250)
  )
  expect_false(251 %in% candidates$tss$tss_pos1)
  expect_true(all(candidates$pairs$assembly == "test1"))
  expect_equal(
    candidates$pairs$anchor_multiplicity[candidates$pairs$pos1 == 101],
    c(2L, 2L)
  )
  expect_equal(nrow(candidates$record_map), 4)

  expect_warning(
    narrow <- build_v2g_candidates(loci, variants, genes, gene_flank_bp = 0),
    class = "acv2g_outside_loci_warning"
  )
  expect_equal(
    narrow$memberships[c("variant_anchor_id", "locus_id")],
    candidates$memberships[c("variant_anchor_id", "locus_id")]
  )
  expect_equal(unique(narrow$pairs$gene_id), "gC")

  expect_warning(
    repeated <- build_v2g_candidates(loci, variants, genes, gene_flank_bp = 50),
    class = "acv2g_outside_loci_warning"
  )
  expect_identical(repeated$pairs$pair_id, candidates$pairs$pair_id)
})
