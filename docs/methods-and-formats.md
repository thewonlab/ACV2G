# ACV2G methods and file-format reference

Start with the [README](../README.md) for installation and a runnable example.

**ACV2G: Activity–Contact Variant-to-Gene Prioritization** is an R package and
command-line interface under authorized private/local development. It combines
variant-anchor regulatory activity with variant-to-promoter chromatin contact
to prioritize the complete variant–gene pair universe within each fixed,
user-defined locus.

ACV2G is inspired by the Activity-by-Contact idea, but it is not the canonical
enhancer-centric ABC model. It uses fixed, user-defined loci and variant
anchors rather than an independently defined enhancer universe.

> **Development and licensing status:** this source tree is for private group
> development, testing, and evaluation. It grants no public open-source
> license and is not an approved public release. Copyright, contributor
> authorization, license, repository, archive, container, and public CI
> decisions are deferred to a separate release-readiness review.

## Scientific contract

The default analysis is locus-centric:

1. users supply fixed locus boundaries;
2. only variant anchors inside each original locus enter that locus's target
   pair universe;
3. candidate protein-coding genes have an eligible TSS within the locus plus
   or minus `gene_flank_bp` (1 Mb by default), and that TSS may lie outside the
   locus boundary;
4. every eligible locus–variant × gene pair is created before contact lookup;
   and
5. the original locus, not the candidate-gene flank, remains the Z-score
   normalization group.

Both implemented methods begin with the visible pair-level proxy

\[
P(v,g)=T_A(A_v)\,T_C(C_{v,g}).
\]

`activity_transform` and `contact_transform` independently default to
`log1p`; either can instead be `identity`. A finite negative input is a
validation error whenever its selected transform is `log1p`.

- `locus_pair_z` is the default. It standardizes all finite variant–gene proxy
  values within each original user-supplied locus using the sample standard
  deviation.
- `gene_wise_z` is an independent optional method. It standardizes proxy
  values within an explicitly declared reference variant universe for each
  gene. It does not run inside or alter `locus_pair_z`.

Per-variant gene allocation, fractions, p-values, FDR values, and canonical
ABC scores are intentionally not implemented.

## Current development status

The independent R engine, YAML validation, standard evidence readers, both
normalization methods, schema-versioned result bundle, CLI dispatcher,
locus-Z plotting, and synthetic/format smoke tests are implemented. Current
work remains private hardening: group testing, Codex-assisted review,
performance profiling, API/CLI simplification, and documentation.

ACV2G has no required Slurm, environment-module, Longleaf, lab-path, or
container dependency. A standard R installation on Linux, macOS, or Windows
is the target runtime. Cluster module commands used by developers are not part
of ACV2G itself.

## Install from the local source tree

Install declared CRAN and Bioconductor dependencies through their normal
package managers, then install this directory. For example:

```r
install.packages(c(
  "cli", "data.table", "digest", "ggplot2", "jsonlite", "optparse",
  "strawr", "svglite", "yaml"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "GenomeInfoDb", "GenomicAlignments", "GenomicRanges", "HiCExperiment", "IRanges",
  "Rsamtools", "rtracklayer", "S4Vectors", "VariantAnnotation"
))
```

From a shell at the parent of this source directory:

```bash
R CMD INSTALL ACV2G
```

No public package archive is configured at this stage.

## Invoke the CLI portably

`R CMD INSTALL` installs the bundled launchers inside the R package; it does
**not** automatically add an `acv2g` executable to the user's `PATH`. The most
portable invocation therefore calls the exported dispatcher with `Rscript`:

```bash
Rscript -e 'quit(save="no", status=ACV2G::acv2g_cli(commandArgs(trailingOnly=TRUE)), runLast=FALSE)' methods
```

On Windows `cmd.exe`, use double quotes around the R expression:

```bat
Rscript -e "quit(save='no', status=ACV2G::acv2g_cli(commandArgs(trailingOnly=TRUE)), runLast=FALSE)" methods
```

The installed Unix/macOS and Windows launcher paths can also be discovered
without assuming an R library location:

```r
launcher <- system.file(
  "exec",
  if (.Platform$OS.type == "windows") "acv2g.cmd" else "acv2g",
  package = "ACV2G",
  mustWork = TRUE
)
if (.Platform$OS.type == "windows") {
  system2("cmd.exe", c("/c", shQuote(launcher), "methods"))
} else {
  system2(launcher, "methods")
}
```

Users may create their own shell wrapper or link named `acv2g` on `PATH`.
Examples below use `acv2g` only as shorthand for one of these portable
invocations; ACV2G does not modify `PATH` during installation.

Implemented commands are:

```text
acv2g methods
acv2g validate --config CONFIG.yaml
acv2g predict --config CONFIG.yaml --output RESULT_DIR [--overwrite]
acv2g gene-wise-z --config CONFIG.yaml --output RESULT_DIR [--overwrite]
acv2g summarize --result RESULT_DIR
acv2g plot-locus-z --result RESULT_DIR --output PLOT.pdf [plot options]
acv2g reference list|fetch|status [reference options]
```

The CLI is a thin wrapper around the same R functions and emits
machine-readable JSON. Input/usage/validation errors return status 2;
unexpected software errors return status 1.

## Minimal runnable example

Relative paths in YAML are resolved against the YAML file, never against the
caller's working directory.

```r
library(ACV2G)

config_path <- system.file(
  "extdata", "example-config.yaml", package = "ACV2G",
  mustWork = TRUE
)
config <- validate_v2g_config(read_v2g_config(config_path))
result <- predict_v2g(config)

v2g_table(result, "pairs")
summarize_v2g_result(result)

write_v2g_result(result, "acv2g-example-result")
plot <- plot_locus_z(result)
save_locus_z_plot(plot, "acv2g-example-locus-z.pdf")
```

Equivalent CLI shorthand:

```bash
acv2g validate --config analysis.yaml
acv2g predict --config analysis.yaml --output results/neuron
acv2g summarize --result results/neuron
acv2g plot-locus-z --result results/neuron --output results/neuron-locus-z.pdf
```

Keep plots and other derived files outside the checksummed result-bundle
directory so that `audit_v2g_result()` continues to verify the bundle exactly.

## Configuration schema

ACV2G validates unknown or misspelled keys rather than silently ignoring them.
A same-build `locus_pair_z` configuration can be written as follows:

```yaml
schema_version: "1.0.0"
analysis_build: GRCh38

inputs:
  loci:
    path: loci.bed
    format: bed
    assembly: GRCh38
  variants:
    path: variants.tsv
    format: tsv
    assembly: GRCh38
  genes:
    path: genes.gtf.gz
    format: gtf
    assembly: GRCh38
    id_version_policy: preserve
  chrom_sizes:
    path: GRCh38.chrom.sizes
    assembly: GRCh38

profiles:
  neuron:
    activity:
      path: neuron.bigWig
      type: bigwig
      assembly: GRCh38
      units: signal
    contact:
      path: neuron.mcool
      type: mcool
      assembly: GRCh38

method:
  name: locus_pair_z
  gene_flank_bp: 1000000
  activity_transform: log1p
  contact_transform: log1p
  max_missing_fraction: 0.20
  call_rule: none

activity:
  activity_window_bp: 500
  aggregation: window_mean
  bam_mapq_min: 30

contact:
  contact_value: balanced
  contact_resolution_bp: 10000
  diagonal_policy: neighbor_max
  distance_policy: none

liftover:
  structural_failure_policy: error
  reciprocal_liftover_policy: none
```

`contact_resolution_bp` may be omitted for a single-resolution native matrix,
which supplies its intrinsic resolution. A multi-resolution source requires
an exact requested resolution; ACV2G never substitutes a nearest resolution.

One named profile is exactly one activity source plus one contact source.
Multiple profiles are scored independently; ACV2G does not infer replicates,
pool assays, or form an activity-by-contact cross-product.

### Complete tabular evidence

Variant-level activity tables use `chrom`, `pos1`, and a declared numeric value
column. Interval tables use `chrom`, `start0`, and `end0`. For example:

```yaml
profiles:
  neuron:
    activity:
      path: activity.tsv
      type: table
      table_mode: variant
      value_col: activity_raw
      assembly: GRCh38
    contact:
      path: contacts.tsv
      type: table
      table_representation: complete
      assembly: GRCh38
      layers:
        balanced: KR
        raw: count
```

A tabular contact source must contain
`chrom1,start1_0,end1_0,chrom2,start2_0,end2_0` plus mapped numeric scores. It
must explicitly declare `table_representation: complete`; sparse called-loop
lists are not accepted as the sole contact source. Score semantics are never
guessed: declare either `score_col` plus `score_type` (`balanced`, `raw`, or
`named_score`) or a `layers` mapping as above.

### Chromosome lengths, builds, chains, and mappings

Every genomic input declares its assembly. `analysis_build` is the contact
build and contacts are never lifted.

- `inputs.chrom_sizes` must explicitly declare `assembly: analysis_build`.
  Table/bedGraph activity also needs chromosome lengths in the activity
  source's own assembly. When activity uses another build, place an
  assembly-matched `chrom_sizes: {path: ..., assembly: ...}` under that
  profile's `activity` specification.
- Structural `inputs.loci`, `inputs.variants`, `inputs.genes`, and optional
  `inputs.background_variants` use `chain` in the direction **input assembly →
  analysis_build** when their assembly differs.
- Native activity remains in its source assembly. Its query `chain` therefore
  runs in the opposite direction: **analysis_build → activity assembly**.
- With `reciprocal_liftover_policy: require`, every affected source also needs
  `reverse_chain` in the exact reverse direction. With the default `none`, no
  round-trip claim is made.
- `structural_failure_policy: error` stops the run. The explicit
  `fail_locus` alternative retains only safely unaffected loci and marks the
  result partial; it never silently changes a failed locus's normalization
  universe.

Optional coordinate-bearing mappings also declare their coordinate assembly:

```yaml
inputs:
  variant_locus_map:
    path: variant-locus-map.tsv
    assembly: GRCh37       # must match the source variants input
  reference_gene_variant_map:
    path: gene-variant-map.tsv
    assembly: GRCh38       # must match analysis_build
```

An `assembly` column in the mapping table can provide the same declaration,
but it must agree with the input specification. An ID-only mapping that has no
`chrom` + `pos1` coordinates does not need a coordinate assembly.

## Candidate genes and annotation

Loci are unique 0-based half-open intervals. Variant anchors require
chromosome plus 1-based position; optional IDs and alleles are annotations.
Repeated submitted records at one assembly-qualified coordinate map to one
scoring anchor and cannot reweight normalization.

Gene/TSS annotation is supplied as GTF, GFF/GFF3, TSV/CSV, or an exact pinned
managed reference fetched separately. Prediction never selects or downloads a
reference. The private-development build ships an empty managed registry; use
a local annotation file or an explicitly supplied, checksum-pinned custom
registry. Provider-specific adapters are not currently implemented. All unique
supplied TSSs for each eligible gene are retained.

Gene-level biotype is authoritative. Transcript-level protein-coding labels
are used only when gene-level classification is unavailable, and an explicit
gene-level non-protein-coding label is never overridden. If neither level can
classify a gene, provide `protein_coding_gene_ids` or an explicit
`is_protein_coding` field. Symbols and ID prefixes never determine
eligibility. Raw IDs, normalized IDs, classification source, conflicts, and
overrides remain auditable.

## Activity and contact evidence

Activity supports BigWig, bedGraph, indexed coordinate-sorted BAM, and
quantitative tables. The default 500-bp total window around `anchor0 = pos1 -
1` is

```text
[anchor0 - floor(width/2), anchor0 + ceiling(width/2))
```

It is clipped, never shifted, at valid chromosome boundaries. The default
`window_mean` treats valid uncovered bases as zero; `covered_mean` excludes
them. Missing chromosomes, incompatible builds, and failed queries remain
missing rather than becoming zero.

BAM uses one fixed `fragment_coverage_cpm` estimand. Users may set the minimum
MAPQ and optionally provide the already filtered library size; primary,
mapped, QC-passing, nonduplicate, nonsupplementary flag exclusions are a
versioned fixed policy. Output records the mode, filter policy, library size,
coverage, and BAM/index resource hashes.

Contacts support `.cool`, `.mcool`, `.hic`, and complete quantitative tables.
Cooler inputs use `HiCExperiment`; native `.hic` uses bounded direct `strawr`
queries and retains the selected normalization layer in provenance.
`contact_value: balanced` is the default. If no balanced layer exists, ACV2G
warns and uses the unique raw layer while recording the fallback. If raw is
also unavailable or ambiguous, validation fails and lists available layers.
An explicitly named missing layer and multiple ambiguous balanced layers fail
without fallback.

The initial implementation has only `distance_policy: none`. It directly uses
the resolved balanced value or warning-recorded raw fallback; there is no
observed/expected or expected-distance model. Same-bin contact defaults to
`neighbor_max`, with `observed` and `keep_na` alternatives. Contact across
multiple TSS bins is aggregated only by finite maximum, with detail and tie QC.

## Zero, missingness, and normalization

Observed contact zero is finite and remains zero. Missing/unavailable contact,
activity, balancing weights, pixels, or mappings remain `NA`; missing-to-zero
coercion is unsupported.

Native `.hic` requires one conservative distinction. If a bounded raw
(`NONE`) query succeeds but omits an otherwise valid requested pixel, that
sparse omission is recorded as `implicit_structural_zero` with value zero. If
a balanced `.hic` query omits the pixel or returns `NA`, ACV2G cannot inspect
the underlying per-bin normalization weights through this reader, so the value
remains `NA` with `missing_balance_weight` rather than being guessed to be
zero. A failed query remains missing under either layer.

Z-score statistics use finite proxies only. A normalization group needs at
least two finite values and a finite, nonzero sample standard deviation.
Otherwise its pair Z-scores remain `NA` with an explicit status. A missing
fraction strictly greater than `method.max_missing_fraction` (20% by default)
emits a warning and `high_missing_fraction` QC but does not suppress an
otherwise computable score.

## Independent `gene_wise_z`

The user must declare the per-gene reference source; there is no implicit
mode or fallback:

```yaml
method:
  name: gene_wise_z
  reference_mode: input_variants
  reference_gene_radius_bp: 1000000
  activity_transform: log1p
  contact_transform: log1p
  max_missing_fraction: 0.20
  call_rule: none
```

`input_variants` uses the complete current input variant set.
`background_variants` instead requires `inputs.background_variants`, including
its explicit assembly and any required structural chain. Within the selected
source, the default reference for a gene is variants within at least one of
its TSSs ± `reference_gene_radius_bp` (1 Mb). An explicit
`inputs.reference_gene_variant_map` completely replaces this distance rule;
it is not merged with it. Reference mode, rule, set identity, counts, radius,
and mapping identity are retained in QC/provenance.

## Calling and plotting

Continuous scores are primary. `call_rule: none` is the overall default. The
only optional call rule is `z_threshold`, which uses `z >= 1.5` unless
`method.z_threshold` is changed. These are dataset-relative prioritization
calls, not statistical significance, p-values, or FDR.

`plot_locus_z()` accepts only `locus_pair_z`. It never pools profiles: one
profile is selected automatically only when exactly one exists; otherwise
`profile_id` is required. A `z_threshold` result automatically displays its
resolved dashed threshold and highlights the stored `link_call` values. A
`call_rule: none` result has no implicit line or red calls. An explicit
`display_threshold` adds or overrides a labelled visual line only; it does not
change stored calls. The default title is blank and the axis label is generic.

```r
p <- plot_locus_z(result, profile_id = "neuron")
p_sample <- plot_locus_z(
  result,
  profile_id = "neuron",
  point_mode = "sample",
  max_points_per_locus = 2000,
  seed = 1
)
save_locus_z_plot(p, "locus-z.pdf")
```

All finite points are shown by default. `point_mode = "sample"` uses
deterministic per-locus sampling and offsets; `"none"` suppresses points.
Captions/provenance record missing/non-finite pairs, non-normalizable pairs and
loci, sampling settings, ordering, package version, and source-result hash.
The CLI additionally exposes `--width` and `--height` in inches for the saved
PDF, PNG, or SVG.

CLI shorthand for a multi-profile bundle:

```bash
acv2g plot-locus-z \
  --result results/all-profiles \
  --profile neuron \
  --point-mode sample \
  --max-points 2000 \
  --seed 1 \
  --output neuron-locus-z.svg
```

## Result bundle

`write_v2g_result()` publishes transactionally and refuses to overwrite an
existing bundle unless `overwrite = TRUE`. The bundle includes the complete
pair table, TSS-bin contact detail, record-to-anchor map, variant–locus QC,
locus status, normalization/activity/contact/annotation/liftover QC,
structured warnings, resolved YAML, run status, JSON manifest, and SHA-256
checksums. Tables carry schema versions. `read_v2g_result(..., lazy = TRUE)`
keeps large tables file-backed, and `audit_v2g_result()` verifies integrity.

## Implemented R API

The exported API is grouped around one computation engine:

- configuration and orchestration: `read_v2g_config()`,
  `validate_v2g_config()`, `validate_v2g_assemblies()`,
  `write_resolved_v2g_config()`, `acv2g_schema_version()`, `predict_v2g()`,
  and `predict_v2g_from_config()`;
- inputs and universes: `read_v2g_loci()`, `read_v2g_variants()`,
  `read_v2g_genes()`, `validate_v2g_profiles()`, `assign_v2g_loci()`, and
  `build_v2g_candidates()`, plus explicit `liftover_v2g_coordinates()`;
- evidence and scoring: `extract_v2g_activity()`, `extract_v2g_contacts()`,
  `compute_v2g_proxy()`, `compute_locus_pair_z()`, `compute_gene_wise_z()`, and
  `call_v2g_links()`;
- results, provenance, and plots: `new_v2g_result()`, `build_v2g_manifest()`,
  `hash_v2g_resource()`, `v2g_table()`, `summarize_v2g_result()`,
  `write_v2g_result()`, `read_v2g_result()`, `audit_v2g_result()`,
  `plot_locus_z()`, and `save_locus_z_plot()`; and
- managed references and interfaces: `list_v2g_references()`,
  `fetch_v2g_reference()`, `resolve_v2g_reference()`, `v2g_methods()`, and
  `acv2g_cli()`.

## Development validation and future release-readiness

The tree includes synthetic unit tests, generated real-format BigWig, indexed
BAM, and `.cool` smoke tests, a real `.hic` smoke test, both normalization
methods, coordinate boundaries, zero-versus-missing behavior, R/CLI
equivalence, transactional bundles, and PDF/PNG/SVG plotting. Historical checks are not a guarantee for a new installation. See the
[validation notes](validation.md) for checks performed when preparing this
repository. Cross-platform validation and biological benchmarking remain
development work.

Only after release-candidate status will a separate review decide copyright
holders and authorizations, an open-source license, a canonical public
repository, Bioconductor versus CRAN, possible Bioconda/BioContainer support,
tagged releases, public checksums/CI, and maintenance. MIT, Bioconductor,
Bioconda, and BioContainer are provisional candidates—not commitments.
