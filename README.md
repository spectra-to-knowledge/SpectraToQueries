
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpectraToQueries <img src='https://raw.githubusercontent.com/spectra-to-knowledge/SpectraToQueries/main/man/figures/logo.svg' align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/spectra-to-knowledge/SpectraToQueries/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/spectra-to-knowledge/SpectraToQueries/actions/workflows/R-CMD-check.yaml)
[![R-Universe](https://spectra-to-knowledge.r-universe.dev/badges/SpectraToQueries)](https://spectra-to-knowledge.r-universe.dev/SpectraToQueries)
[![Codecov test
coverage](https://codecov.io/gh/spectra-to-knowledge/SpectraToQueries/graph/badge.svg)](https://app.codecov.io/gh/spectra-to-knowledge/SpectraToQueries)
<!-- badges: end -->

Repository to translate spectra to queries.

## Requirements

Here is what you *minimally* need:

- A file containing spectra grouped by classes to extract class-specific
  queries (the one of the MIADB is attached as demo).

## Installation

As the package is not (yet) available on CRAN, you will need to install
with:

``` r
install.packages(
  "SpectraToQueries",
  repos = c(
    "https://spectra-to-knowledge.r-universe.dev",
    "https://bioc.r-universe.dev",
    "https://cloud.r-project.org"
  )
)
```

## Use

To generate diagnostic ions queries from your spectra:

``` r
SpectraToQueries::spectra_to_queries()
```

## Main Citations

TODO

### Others

- The *RforMassSpectrometry* packages suite:
  <https://doi.org/10.3390/metabo12020173>
