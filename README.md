# ASRgwas

`ASRgwas` is a free add-on R package developed by [VSN International](https://vsni.co.uk/asrgwas/) that provides a comprehensive set of tools for performing Genome-Wide Association Studies (GWAS) using the modelling flexibility available in `ASReml-R`. This library assists with preparing and auditing phenotypic and genomic data; fitting GWAS models and identifying significant markers; and evaluating and exploring output.

---

## Prerequisites

Before installing `ASRgwas`, make sure the following requirements:

1. **R** (version 4.0 or later) — [https://cran.r-project.org/](https://cran.r-project.org/)
2. **Rtools** (Windows only) — required to build packages from source. Available at [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/).

---

## Step 1 — Install `ASReml-R`

`ASRgwas` depends on `ASReml-R`, which is **not** available on CRAN and requires a valid licence from VSN International.

Check the instructions for installing `ASReml-R` [here](https://asreml.kb.vsni.co.uk/knowledge-base/asreml-r-installation-guide/)

## Step 2 — Install `ASRtools` from GitHub

```r
install.packages("remotes")
remotes::install_github("VSNinternational/asrtools")
```

or

```r
install.packages("devtools")
devtools::install_github("VSNinternational/asrtools")
```

## Step 3 — Install `ASRgwas` from GitHub

```r
install.packages("remotes")
remotes::install_github("VSNinternational/asrgwas")
```

or

```r
install.packages("devtools")
devtools::install_github("VSNinternational/asrgwas")
```

> If you are behind a corporate proxy, you may need to set the `HTTP_PROXY` and `HTTPS_PROXY` environment variables before running the commands above.

---

## Step 4 — Verify the installation

After installation, load the package to confirm everything worked:

```r
library(asrgwas)

# Check the installed version
packageVersion("asrgwas")
```

If no errors appear, the package is ready to use.

---
 
