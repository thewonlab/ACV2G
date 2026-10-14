# Preparing inputs for your first analysis

This guide uses one profile named `NPC`. The name is only a label; replace
it with your own cell type and use matching activity and contact data.
For two profiles, add a second named entry under `profiles` with its own
activity/contact files. ACV2G scores profiles independently.

## Decide the analysis universe first

Supply fixed locus intervals and the variants you want to use for scoring.
For locus-wise normalization, every eligible pair in the supplied universe
contributes to the mean and SD. Filtering variants by a GWAS P-value changes
that universe; ACV2G neither selects a threshold nor filters by disorder.

If you are interested in emVars, score your intended reference universe
first, then join the output to your emVar list. Coordinate matching must use
the same genome build. Locus membership is determined by coordinates, not
by whether an rsID occurs in a GWAS credible-set file.

For independent disorder-specific backgrounds, use separate configurations
and result bundles. A profile denotes an activity/contact combination, not a
GWAS disorder filter.

## Minimal input tables

Examples below show TSV columns. Use tabs between fields and include headers.
For all genomic sources, supply an explicit assembly label. Use consistent
labels throughout (for example `GRCh38`) and matching chromosome naming
(`chr1` versus `1`). Changing an assembly label does not lift coordinates.

**Loci — `loci.tsv`:** 0-based, half-open intervals. The example contains
1-based positions 100001 through 120000, inclusive.

```text
locus_id	chrom	start0	end0
locus_1	chr1	100000	120000
```

**Variants — `variants.tsv`:** `chrom` and 1-based `pos1` are essential;
`variant_id` is a useful label. One genomic position is one scoring anchor,
even if multiple submitted records or rsIDs map to it.

```text
variant_id	chrom	pos1
variant_A	chr1	100001
variant_B	chr1	105001
variant_C	chr1	110001
```

**Gene annotation — `genes.tsv`:** stable gene IDs and 1-based TSSs.
Multiple TSS rows per gene are allowed. Symbols are optional annotations.

```text
gene_id	chrom	tss_pos1	strand	gene_biotype	gene_symbol
ENSG_TEST001	chr1	130001	+	protein_coding	TEST1
```

This TSS lies outside the locus, but inside its default ±1 Mb gene-search
flank. Gene-level classification controls protein-coding eligibility;
transcript labels are a fallback only when gene classification is missing.
Unclassifiable genes need explicit protein-coding eligibility information.
See [annotation details](methods-and-formats.md#candidate-genes-and-annotation).

**Chromosome sizes — `chrom_sizes.tsv`:** lengths in bases in the analysis build.

```text
chrom	length
chr1	1000000
```

This length is synthetic. Use the actual chromosome lengths from your own
assembly for biological analyses. The headered format above is recommended.

## A same-build BigWig + Cooler YAML

Save this as `analysis.yaml`. Replace the file names with your own files;
this template is not runnable until those files exist. The complete runnable
synthetic alternative ships in `inst/extdata/example-config.yaml`.

```yaml
schema_version: "1.0.0"
analysis_build: GRCh38

inputs:
  loci:
    path: loci.tsv
    format: tsv
    assembly: GRCh38
  variants:
    path: variants.tsv
    format: tsv
    assembly: GRCh38
  genes:
    path: genes.tsv
    format: tsv
    assembly: GRCh38
    id_version_policy: preserve
  chrom_sizes:
    path: chrom_sizes.tsv
    assembly: GRCh38

profiles:
  NPC:
    activity:
      path: NPC.bigWig
      type: bigwig
      assembly: GRCh38
      units: signal
    contact:
      path: NPC.mcool
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

contact:
  contact_value: balanced
  contact_resolution_bp: 10000
  diagonal_policy: neighbor_max
  distance_policy: none

liftover:
  structural_failure_policy: error
  reciprocal_liftover_policy: none
```

For GTF annotation, replace `genes.path` with your GTF and set `format: gtf`.
For a single-resolution `.cool`, set `type: cool`; you may omit the resolution
and use the matrix's native resolution. A `.mcool` requires an exact available
resolution: the package does not choose the nearest one.

## Settings most users should understand

| Setting | Default | Meaning |
| --- | --- | --- |
| `gene_flank_bp` | 1000000 | Candidate TSS search beyond each end of a fixed locus |
| `activity_window_bp` | 500 | Total window width around a variant anchor, not ±500 bp |
| `aggregation` | `window_mean` | Mean activity over valid window bases; uncovered bases contribute zero |
| `activity_transform`, `contact_transform` | `log1p` | Transform each evidence value before multiplying; negative finite values are rejected |
| `contact_value` | `balanced` | Balanced layer if available; warning-recorded raw fallback if none exists |
| `diagonal_policy` | `neighbor_max` | Special handling when the variant and TSS occupy the same Hi-C bin |
| `distance_policy` | `none` | No observed/expected or other distance correction |
| `max_missing_fraction` | 0.20 | Warning threshold, not a filter that discards the whole group |
| `call_rule` | `none` | Continuous scores only; use `z_threshold` for optional calls |
| `z_threshold` | 1.5 when calling is enabled | Calls use Z greater than or equal to the threshold |

For a `.hic` file with several balanced normalizations, specify the exact
`contact.contact_layer` reported by the source (for example `KR`, if available).
An ambiguous or explicitly missing layer fails validation; the package does
not silently substitute another named normalization.

## Using quantitative tables

Copy the package's `inst/extdata` directory to a writable location if you
want to edit the synthetic example. Keep its YAML and tables together;
do not edit the installed library's files.

Activity at variant anchors uses `chrom`, `pos1`, and a numeric value column:

```yaml
activity:
  path: activity.tsv
  type: table
  table_mode: variant
  value_col: activity_raw
  assembly: GRCh38
```

The block above replaces `profiles.NPC.activity`. Interval activity tables
instead use `chrom`, `start0`, `end0` and `table_mode: interval`.

Contact tables use the following coordinate columns:

```text
chrom1  start1_0  end1_0  chrom2  start2_0  end2_0  balanced  raw
```

Both bins use 0-based half-open coordinates. Replace `profiles.NPC.contact`
with an explicit score mapping:

```yaml
contact:
  path: contacts.tsv
  type: table
  table_representation: complete
  layers:
    balanced: balanced
    raw: raw
  assembly: GRCh38
```

Supply every required pixel, including valid zero values. Do not use a
thresholded significant-loop list and assume absent entries are zero.
Unmeasured/unavailable values must remain missing. `table_representation:
complete` is a declaration about your input, not a command to reconstruct
omitted measurements.

## Genome builds and liftover

The `analysis_build` is the contact matrix's build; the contact matrix itself
is never lifted. If structural inputs use a different build, provide chains
from their input build to the analysis build. Activity remains in its source
build: its query chain runs from the analysis build to the activity build.

For example, querying a GRCh37 BigWig for GRCh38 variants requires a
GRCh38 → GRCh37 activity chain. It is not valid to query the original GRCh38
positions directly in that BigWig. See the
[chain and chromosome-size configuration](methods-and-formats.md#chromosome-lengths-builds-chains-and-mappings)
and inspect the saved liftover QC.

## Gene-wise reference variants

For a gene-wise run, change the method block to include:

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

The entire input variant set is the reference source, then each gene's
reference is distance-limited to its TSSs ±1 Mb by default. To use a separate
reference, choose `reference_mode: background_variants` and add an
`inputs.background_variants` file specification with its format and assembly.
An explicit gene–variant mapping can replace the distance rule entirely.
Keep the loci and other structural inputs: they still define the target pairs.

## Before scaling up

Run the bundled synthetic example, then a small real-data configuration.
Inspect `variant_locus_qc`, `normalization_qc`, activity/contact QC, and
warnings before launching a larger job. Runtime and memory depend on the
complete pair universe and source matrices; a short list of emVars alone
does not determine the workload. On a cluster, follow its compute-node policy.
