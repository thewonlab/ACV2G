# ACV2G

**Activity–Contact Variant-to-Gene Prioritization**

ACV2G helps researchers prioritize genes that may be regulated by variants of
interest. It combines local ATAC-seq activity with Hi-C contact between each
variant and candidate gene transcription start sites (TSSs). You supply the
genomic regions, variants, annotation, and cell-type-specific profiles; ACV2G
returns scored variant–gene pairs, quality-control tables, and analysis provenance.

The same analysis is available as an **R package** and a **command-line
interface (CLI)**. RStudio, R scripts, and cluster jobs can use the same
configuration.

## Start here

- [Install](#installation) and run the [small example](#quick-start-in-r).
- [Prepare your own inputs and YAML](docs/input-guide.md).
- [Understand the scores](#how-the-scores-work) and [read the results](#understanding-your-results).
- Use the [CLI](#command-line-use) or the [step-by-step R workflow](vignettes/acv2g-workflow.Rmd).
- Consult [installation troubleshooting](INSTALL.md) and the
  [detailed methods and formats](docs/methods-and-formats.md).

## Installation

You need **R 4.4.0 or newer**.
The first installation also installs CRAN and Bioconductor dependencies;
it can take considerably longer than running the small example.

Run this in R or RStudio after configuring your GitHub credentials as
described in [INSTALL.md](INSTALL.md#access-to-the-private-repository):

```r
install.packages(c("BiocManager", "remotes"))

remotes::install_github(
  "thewonlab/ACV2G",
  repos = BiocManager::repositories(),
  dependencies = NA,
  upgrade = "never",
  build_vignettes = FALSE
)

library(ACV2G)
packageVersion("ACV2G")
```

This installs into your R library. After installation, a normal session
starts with `library(ACV2G)` (capitalization matters).

**Updates are not automatic.** When the team pushes an update, restart R
and rerun `remotes::install_github(...)` above. Development commits may share
a version number; record `packageDescription("ACV2G")$RemoteSha` and
`sessionInfo()` with your analysis.

For local source installation, user-library permissions, dependency details,
and cluster use, see [INSTALL.md](INSTALL.md).

## Quick start in R

This example uses tiny synthetic tables bundled with the package: one locus,
three variant anchors, one gene, and one activity/contact profile.
No GWAS, ATAC-seq, or Hi-C download is needed.

```r
library(ACV2G)

config_file <- system.file(
  "extdata", "example-config.yaml",
  package = "ACV2G", mustWork = TRUE
)

# Check configuration keys, declared builds, and required files.
config <- validate_v2g_config(config_file)

# Run and keep the complete result in R memory.
result <- predict_v2g_from_config(config)

summarize_v2g_result(result)
pairs <- v2g_table(result, "pairs")
head(pairs[c("locus_id", "variant_id", "gene_id", "locus_pair_z")])
```

Expected summary: `run_status = "success"`, 1 profile, 1 locus, 3 variants,
1 gene, 3 pairs, and 3 finite scores. The scores have mean zero and sample
standard deviation one within this example locus.

Validation does not replace execution: unreadable native matrices or
unavailable evidence layers may still be detected during a run.

Save the results and a plot:

```r
# Choose a new directory name; existing bundles are protected by default.
write_v2g_result(result, "acv2g-example-result")

p <- plot_locus_z(result)
save_locus_z_plot(p, "acv2g-example-locus-z.pdf")

# Read the saved bundle in a later session.
saved <- read_v2g_result("acv2g-example-result")
audit_v2g_result("acv2g-example-result")
```

Store plots and exported subsets **outside** the checksummed result-bundle
directory. To intentionally replace a bundle, use `overwrite = TRUE`.

For a one-step run that saves directly to disk:

```r
predict_v2g_from_config(
  config_file,
  output = "acv2g-example-result-second-run"
)
```

With `output = ...`, the return value is the saved path, invisibly.
Omit `output` when you want the result object in memory.

## Running your own data

Prepare fixed loci, variants, protein-coding gene/TSS annotation, chromosome
lengths, and one activity/contact pair per profile. Start with the
[YAML template and input examples](docs/input-guide.md).

```r
# Replace analysis.yaml with your own configuration.
# config <- validate_v2g_config("analysis.yaml")
# result <- predict_v2g_from_config(config)
# write_v2g_result(result, "results/my-analysis")
```

All relative input paths are interpreted relative to the YAML file.
A profile such as `NPC` or `neuron` is a user-chosen name for one activity
source plus one contact source. Profiles are scored independently.

| Input | Supported file formats |
| --- | --- |
| Loci | BED, TSV, CSV |
| Variants | TSV, CSV, BED, VCF/VCF.gz |
| Gene/TSS annotation | GTF, GFF/GFF3, TSV, CSV |
| Activity | BigWig, bedGraph, indexed coordinate-sorted BAM, quantitative tables |
| Contact | Cooler `.cool`/`.mcool`, Juicer `.hic`, complete quantitative contact tables |

The [input guide](docs/input-guide.md) explains columns and coordinates.
A list of significant Hi-C loops alone is not a complete contact matrix.

### Choosing variants and loci

A locus is a **fixed, user-defined genomic interval**, for example a GWAS
region defined by a study's lead signals and LD analysis. ACV2G does not call
GWAS loci or automatically apply a GWAS P-value filter.

Choose locus boundaries and the reference variant set before scoring.
For locus-wise analysis, all eligible pairs from supplied variants inside a
locus contribute to normalization. To interpret emVars against a larger GWAS
universe, include that universe in the run and subset the scored table to
emVars afterward. Running only emVars changes the reference and the question.

## How the scores work

```text
Fixed loci + variants + gene annotation
    → complete candidate variant–gene pairs
    → activity and contact extraction
    → activity × contact proxy
    → normalization
    → scored pairs, optional calls, QC, and plots
```

The default proxy and normalization are:

```text
P(v, g) = log(1 + activity(v)) × log(1 + contact(v, g))
Z       = (P − mean of reference proxies) / sample SD of reference proxies
```

Each transform can instead be `identity`. ACV2G is inspired by
Activity-by-Contact, but these variant-based Z-scores are not canonical
enhancer-based ABC scores.

| Method | Reference for mean and SD | Interpretation |
| --- | --- | --- |
| `locus_pair_z` (default) | All finite candidate-pair proxies in the original locus, per profile | Relative strength of a pair within a genomic region |
| `gene_wise_z` (optional) | Finite proxies for one gene across its explicitly declared reference variants, per profile | Relative strength of a variant for that gene |

Variants must lie inside the original locus. Candidate genes can have a TSS
within the **locus boundaries ±1 Mb** (`gene_flank_bp`), including outside
the locus. This gene-search flank does not expand the normalization locus.

Gene-wise scoring is independent. Try it with the synthetic inputs:

```r
gene_config <- read_v2g_config(config_file)
gene_config$method$name <- "gene_wise_z"
gene_config$method$reference_mode <- "input_variants"
gene_config$method$reference_gene_radius_bp <- 1000000
gene_result <- predict_v2g_from_config(gene_config)
gene_pairs <- v2g_table(gene_result, "pairs")
head(gene_pairs[c("chrom", "pos1", "gene_id", "gene_wise_z")])
```

`input_variants` uses the complete input variant set as the reference source.
Alternatively, declare `background_variants` and provide a separate file.
By default, references are restricted to gene TSS ±1 Mb. Explicit gene–variant
mapping can override this distance rule. See the
[reference-universe details](docs/methods-and-formats.md#independent-gene_wise_z).

## Understanding your results

`v2g_table(result, "pairs")` is the main table. It retains missing scores
and uncalled pairs to distinguish weak evidence from unavailable evidence.

| Column | Meaning |
| --- | --- |
| `profile_id`, `locus_id` | Cell profile and locus |
| `variant_id`, `variant_anchor_id` | Scoring identifiers; the pipeline's `variant_id` may be an internal anchor ID rather than your submitted rsID |
| `gene_id`, `gene_id_normalized` | Gene identifiers; use stable IDs when joining datasets |
| `activity_raw`, `contact_raw` | Evidence before transforms; `contact_raw` can hold balanced contact values, not necessarily raw counts |
| `proxy_score` | Product of transformed activity and contact |
| `locus_pair_z` or `gene_wise_z` | Score relative to the method's reference group |
| `normalization_status` | Whether the group statistics could be calculated |
| `link_call` | Optional descriptive threshold call; no calls are assigned by default |

Repeated input rows at one genomic coordinate share a scoring anchor.
Use `variant_record_anchor_map` to recover every submitted record/ID.

```r
record_map <- subset(
  v2g_table(result, "variant_record_anchor_map"),
  input_role == "target_variants"
)
labels <- unique(record_map[c("variant_anchor_id", "variant_id")])
names(labels)[2] <- "input_variant_id"
display_pairs <- merge(pairs, labels, by = "variant_anchor_id", sort = FALSE)
head(display_pairs[c("input_variant_id", "gene_id", "chrom", "pos1",
                     "proxy_score", "locus_pair_z")])
```

For the bundled example, rounded values are:

| Input variant | Gene | Chromosome | Position (1-based) | Proxy | Locus Z |
| --- | --- | --- | ---: | ---: | ---: |
| variant_A | ENSG_TEST001 | chr1 | 100001 | 0.761500 | -0.794936 |
| variant_B | ENSG_TEST001 | chr1 | 105001 | 1.206949 | -0.327830 |
| variant_C | ENSG_TEST001 | chr1 | 110001 | 2.590290 | 1.122766 |

This display join can expand rows when a coordinate has several submitted
IDs; use the original pair table to count distinct scored pairs.

```r
head(v2g_table(result, "normalization_qc"))
head(v2g_table(result, "activity_qc"))
head(v2g_table(result, "contact_qc"))
v2g_table(result, "warnings")
```

Observed zero remains zero; missing activity/contact remains `NA`.
Normalization needs at least two finite proxies and nonzero variation.
A missing fraction above 20% warns but does not suppress computable scores.
Review QC before interpreting high scores.

### Ranking and optional calls

Higher Z-scores indicate stronger proxies relative to a group's reference.
Scores from different reference universes are not calibrated probabilities.

```r
one_locus <- subset(
  pairs, locus_id == "locus_1" & profile_id == "example_profile"
)
ranked <- one_locus[order(one_locus$locus_pair_z, decreasing = TRUE,
                          na.last = TRUE), ]
head(ranked[c("chrom", "pos1", "gene_id", "locus_pair_z")])
```

To enable calls, add `call_rule: z_threshold` and `z_threshold: 1.5` to
the YAML's existing `method` block, or edit the configuration in R:

```r
called_config <- read_v2g_config(config_file)
called_config$method$call_rule <- "z_threshold"
called_config$method$z_threshold <- 1.5
called_result <- predict_v2g_from_config(called_config)
called_pairs <- subset(v2g_table(called_result, "pairs"), link_call %in% TRUE)
write.csv(called_pairs, "called-pairs.csv", row.names = FALSE)
```

Calls use **Z ≥ 1.5** by default under this rule. This is a prioritization
cutoff, **not a statistical significance test, P-value, or FDR**.
Missing scores keep a missing call. The three-pair example can yield
no calls at this threshold; that is expected.

For multi-profile plots, select a profile explicitly:

```r
# For your own result with a profile named NPC:
# p <- plot_locus_z(result, profile_id = "NPC")
# save_locus_z_plot(p, "NPC-locus-z.pdf")
```

Locus plots support `locus_pair_z` only. Titles are blank by default;
add one with `p + ggplot2::labs(title = "Your title")`.

## Command-line use

The CLI uses the same installed R package and YAML as the R API.
From the root of a cloned repository on Linux/macOS:

```bash
Rscript exec/acv2g methods
Rscript exec/acv2g validate --config inst/extdata/example-config.yaml
Rscript exec/acv2g predict --config inst/extdata/example-config.yaml --output cli-example-result
Rscript exec/acv2g summarize --result cli-example-result
Rscript exec/acv2g plot-locus-z --result cli-example-result --output cli-example-locus-z.pdf
```

On Windows Command Prompt, use `exec/acv2g.cmd` with Windows path separators
instead of `Rscript exec/acv2g`. If you installed from R without cloning,
invoke the dispatcher from any directory (Linux/macOS shell):

```bash
Rscript -e 'quit(save="no", status=ACV2G::acv2g_cli(commandArgs(trailingOnly=TRUE)), runLast=FALSE)' methods
```

Installation does not automatically add `acv2g` to your shell's PATH.
[INSTALL.md](INSTALL.md#using-the-cli-after-installation) explains installed
launchers and Windows quoting.

`predict` runs the method in the YAML. `gene-wise-z` requires a
configuration declaring `method.name: gene_wise_z`.
Use `--help` for options. CLI summaries/errors use JSON; exit codes are
0 for success, 2 for input/usage errors, and 1 for unexpected software errors.

## Help and development

Function documentation is installed with the package:

```r
?predict_v2g_from_config
?compute_gene_wise_z
?plot_locus_z
```

Open an [issue](https://github.com/thewonlab/ACV2G/issues) with your R version,
operating system, package version/commit, error, and a small reproducible
configuration. Remove credentials and restricted data.

The [workflow vignette source](vignettes/acv2g-workflow.Rmd) shows the
step-by-step API. [Validation notes](docs/validation.md) describe this
snapshot's checks. [DESCRIPTION](DESCRIPTION) is the authoritative dependency
list; a Python requirements file is not needed. Public distribution,
license, release channels, and formal citation details remain future decisions.
