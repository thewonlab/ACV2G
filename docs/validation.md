# Validation of the private repository snapshot

Date: 2026-10-14. Package version: `0.0.0.9000`.

These checks exercise the current source and onboarding examples; they are
not a biological benchmark or a public release certification.

## Environment and scope

- Linux x86_64, Red Hat Enterprise Linux 9.8, R 4.4.0.
- Existing CRAN/Bioconductor dependencies were available locally.
- The package was installed from source into a separate test library.
- Only synthetic and bundled format-test data were used. No research
  analysis was rerun, no cluster job was submitted, and no data download
  was needed for the examples or tests.

## Checks performed

| Check | Result |
| --- | --- |
| Source installation with `R CMD INSTALL` | Passed |
| Existing automated test suite | 121 cases, 654 passing expectations, 0 failures, 0 errors, 0 skips, 0 unexpected warnings |
| README R examples | Locus-wise and gene-wise scoring, input-ID display join, ranking, optional calling, result write/read/audit, and PDF generation passed |
| Synthetic example summary | 1 profile, 1 locus, 3 variant anchors, 1 gene, 3 pairs, 3 finite locus scores |
| CLI validate / predict / summarize / plot-locus-z | Passed with the bundled YAML and installed package |
| Source archive build | Passed; tar emitted a host UID-mapping warning |
| `R CMD check --no-manual --ignore-vignettes` | `Status: OK`; installed-package tests also reported 654 passes |

The check could not reach external package indexes from the restricted
checking environment. It completed using installed dependencies, so it does
not verify downloading dependencies into a fresh machine. The PDF manual
and vignette rebuild were deliberately excluded from this lightweight check.

The optional 1.5 threshold produces zero calls on the tiny synthetic dataset.
Its locus Z-scores are approximately -0.794936, -0.327830, and 1.122766.

## Repeating development checks

From a source checkout with dependencies installed, use:

```r
install.packages(c("testthat", "withr", "knitr", "rmarkdown"))
testthat::test_local(".")
```

For a lightweight package check without rebuilding the vignette or PDF manual,
run from the parent of the checkout:

```bash
R CMD build --no-build-vignettes ACV2G
R CMD check --no-manual --ignore-vignettes ACV2G_0.0.0.9000.tar.gz
```

These checks do not establish a fresh-machine dependency installation,
macOS/Windows compatibility, large-dataset performance, or scientific
calibration. Source installation was tested locally; installation directly
from the private GitHub repository requires working access credentials.
Dependency versions are not locked. Record your own `sessionInfo()` and
source commit when reporting results.
