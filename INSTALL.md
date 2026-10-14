# Installing and updating ACV2G

Start in the R or RStudio session where you intend to use ACV2G.
You need R >= 4.4.0 and internet access for the first dependency installation.
The package has been exercised on Linux with R 4.4.0. macOS and Windows are
intended platforms, but have not been independently validated for this snapshot.

## Access to the private repository

Ask the repository maintainers to grant your GitHub account access to
`thewonlab/ACV2G`. Confirm that you can open the repository in your browser.
Browser access and R/Git authentication are separate: a private repository
can return HTTP 404 when the installing tool is unauthenticated.

For `remotes::install_github()`, configure an authorized GitHub token in your
local `GITHUB_PAT` environment variable using your institution's credential
practice. A fine-grained token, where permitted, needs this repository's
Contents read permission; organization approval or SSO authorization may be
required. Never put a token in a committed script, YAML, notebook, or issue.

If Git authentication is already configured, you can instead clone the
repository and use the local installation below. This avoids needing an
additional token in R. You can also use GitHub's **Code → Download ZIP** in
your authenticated browser, extract it, and install that source directory.

## Install from GitHub

Run in R/RStudio:

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
```

`BiocManager::repositories()` supplies CRAN plus the Bioconductor repositories
appropriate for your R installation. `dependencies = NA` installs required
dependencies; it does not install every optional development tool.
`upgrade = "never"` avoids updating unrelated installed packages during this
operation. It does not suppress installation of missing required packages.
Vignettes are not built during this installation; R help pages are installed,
and the workflow vignette source remains available on GitHub.

The package version is currently `0.0.0.9000`. It is a development version,
not a stable release. A GitHub install records its commit as `RemoteSha`:

```r
packageVersion("ACV2G")
packageDescription("ACV2G")$RemoteSha
sessionInfo()
```

## Install from a local source directory

Clone in a terminal, or extract the downloaded ZIP:

```bash
git clone https://github.com/thewonlab/ACV2G.git
```

Then set R's working directory to the directory containing `DESCRIPTION`
(the repository root) and run:

```r
install.packages(c("BiocManager", "remotes"))
remotes::install_local(
  ".",
  repos = BiocManager::repositories(),
  dependencies = NA,
  upgrade = "never",
  build_vignettes = FALSE
)
library(ACV2G)
```

With dependencies already installed, `R CMD INSTALL /path/to/ACV2G` is
also sufficient. A local install does not normally record `RemoteSha`;
record `git rev-parse HEAD` from the source checkout for reproducibility.

## Dependencies and system requirements

[DESCRIPTION](DESCRIPTION) is the authoritative dependency manifest.
The package currently imports the following non-base R packages:

| Source | Packages |
| --- | --- |
| CRAN | cli, data.table, digest, ggplot2, jsonlite, optparse, strawr, svglite, yaml |
| Bioconductor | GenomeInfoDb, GenomicAlignments, GenomicRanges, HiCExperiment, IRanges, Rsamtools, rtracklayer, S4Vectors, VariantAnnotation |
| Optional development/documentation | knitr, rmarkdown, testthat (>= 3.0.0), withr |

These are declared package imports, including native-format readers: they are
required for installation even if your first analysis uses only TSV files.
There is no separate Python, conda, Cooler command-line, or Slurm dependency.

Some dependencies contain compiled code. Use compatible binary packages when
available. Source installation can require a C/C++/Fortran toolchain and
development libraries used by the R dependencies (for example libcurl,
OpenSSL, XML, compression, fonts, and text rendering libraries). On Windows,
use the Rtools version matching your R version if compilation is needed;
on macOS, the relevant command-line compiler tools; on Linux or a cluster,
ask the system administrator to provide libraries named in the installation
error. ACV2G itself has no compiled source directory.

Do not combine arbitrary Bioconductor releases. `BiocManager::valid()` helps
identify a mixed installation. See the official
[Bioconductor installation guide](https://bioconductor.org/install/) and
[remotes installation options](https://remotes.r-lib.org/reference/install_github.html).

## Your R library and RStudio

Normally, accept R's offer to install into a personal library if the system
library is not writable. Check your active search paths with `.libPaths()`.
You do not need administrator privileges or a copy of the developers' library.

If you need to create the default user library manually:

```r
user_lib <- Sys.getenv("R_LIBS_USER")
stopifnot(nzchar(user_lib))
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
```

This creates a directory for packages and puts it on this session's R library
search path. For the standard `R_LIBS_USER` location, a new R session includes
it once the directory exists. Install and run with the same R major/minor
version; RStudio and a cluster job can otherwise see different libraries.

On Longleaf, loading an R module is part of starting an R environment in a
shell/job, not an ACV2G command. Once RStudio is running the appropriate R
version, installation and use happen inside R as above. Do not reinstall the
package in every analysis job.

## Updates and reproducible analyses

Restart R before updating a loaded package. For GitHub installs, rerun the
installation command; for a clone, run `git pull --ff-only` and reinstall
the local source. Pulling files alone does not update the installed R library.

For an exact version, replace `HEAD` with an actual recorded commit SHA:

```r
# remotes::install_github(
#   "thewonlab/ACV2G", ref = "<full-commit-sha>",
#   repos = BiocManager::repositories(),
#   dependencies = NA, upgrade = "never", build_vignettes = FALSE
# )
```

Keep the source commit, `sessionInfo()`, input identifiers, YAML, and complete
result bundle together in your analysis records.

## Using the CLI after installation

An R installation does not automatically add `acv2g` to your shell's PATH.
You can locate its bundled launcher from R:

```r
launcher <- system.file(
  "exec", if (.Platform$OS.type == "windows") "acv2g.cmd" else "acv2g",
  package = "ACV2G", mustWork = TRUE
)
launcher
```

Use `Rscript /printed/path/to/exec/acv2g methods` on Linux/macOS, or invoke
the printed `.cmd` path on Windows. You can also call the dispatcher directly:

```bash
Rscript -e 'quit(save="no", status=ACV2G::acv2g_cli(commandArgs(trailingOnly=TRUE)), runLast=FALSE)' methods
```

Windows Command Prompt:

```bat
Rscript -e "quit(save='no', status=ACV2G::acv2g_cli(commandArgs(trailingOnly=TRUE)), runLast=FALSE)" methods
```

Replace `methods` with a command and its arguments, for example
`predict --config analysis.yaml --output result-bundle`.

## Common problems

| Symptom | What to check |
| --- | --- |
| GitHub returns 404/401/403 | Repository access, R/Git credentials, organization approval/SSO; browser login alone is not enough |
| `there is no package called 'ACV2G'` | Correct capitalization, R version, `.libPaths()`, and whether installation completed successfully |
| A Bioconductor dependency is unavailable | Use `BiocManager::repositories()` and a compatible R/Bioconductor release |
| A dependency fails to compile | Read the first compiler/system-library error; later missing-package errors are often consequences |
| `acv2g: command not found` | Use `Rscript` with the launcher or the dispatcher above |
| Input file not found | Relative input paths are relative to the YAML file |
| All/some scores are `NA` | Inspect activity/contact and normalization QC; check builds and reference-group variation |
| No pairs are called | Default `call_rule` is `none`; explicit threshold calls also need finite scores meeting the cutoff |
| Output directory already exists | Choose a new path or explicitly request overwrite |

Use the bundled quick start before debugging a large biological dataset.
