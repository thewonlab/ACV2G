# ACV2G 0.0.0.9000

- Prepared private GitHub onboarding documentation: installation and updates,
  runnable R/CLI examples, input/YAML guidance, and result interpretation.
  The scoring engine, APIs, and scientific defaults are unchanged.
- Began private/local implementation from the confirmed ACV2G specification.
- Added the default fixed-locus `locus_pair_z` workflow and independent
  `gene_wise_z` workflow with explicit input/background reference modes.
- Added activity readers for BigWig, bedGraph, indexed BAM, and quantitative
  tables; contact readers support complete tables, native `.cool`/`.mcool`
  through `HiCExperiment`, and bounded native `.hic` queries through `strawr`.
- Added fixed candidate-universe construction, protein-coding gene/TSS
  classification, coordinate/build validation, role-aware liftover, and
  zero-versus-missing evidence semantics.
- Added schema-versioned transactional result bundles, checksums, lazy reads,
  structured QC/warnings, locus-Z plotting, and portable R/CLI entry points.
- Added synthetic tests, real-format evidence smoke tests, an end-to-end
  example, R/CLI equivalence checks, and private-development documentation.
- Added bounded, exact-key native contact queries and vectorized table,
  candidate, activity, and normalization paths for practical locus-scale use.
- Vectorized same-build identity liftover so large analysis-build gene/TSS
  tables no longer construct and unpack one `GRanges` object per record.
- Preserved partial native contact and activity query failures as explicit
  missing evidence with structured warnings instead of discarding successful
  queries or coercing failures to zero.
- Enforced negative-activity rejection before interval aggregation under
  `log1p`, and recorded native normalization names, BAM/filter policy, and BAM
  index hashes in provenance/QC.
- Added generated real-format BigWig, indexed BAM, and Cooler tests plus a real
  `.hic` smoke test and an executable installed-data vignette example.
- Public licensing and release channels remain intentionally deferred until a
  separate release-readiness review.
