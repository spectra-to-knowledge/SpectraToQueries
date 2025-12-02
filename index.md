

<!-- README.md is generated from README.qmd. Please edit that file -->

# SpectraToQueries <img src='https://raw.githubusercontent.com/spectra-to-knowledge/SpectraToQueries/main/man/figures/logo.svg' align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/spectra-to-knowledge/SpectraToQueries/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/spectra-to-knowledge/SpectraToQueries/actions/workflows/R-CMD-check.yaml)
[![r-universe
badge](https://spectra-to-knowledge.r-universe.dev/SpectraToQueries/badges/version?&color=blue&style=classic.png)](https://spectra-to-knowledge.r-universe.dev/SpectraToQueries)
[![Codecov test
coverage](https://codecov.io/gh/spectra-to-knowledge/SpectraToQueries/graph/badge.svg)](https://app.codecov.io/gh/spectra-to-knowledge/SpectraToQueries)
<!-- badges: end -->

Repository to translate spectra to queries.

## Requirements

Here is what you *minimally* need:

- **A file containing MS/MS spectra with associated skeleton
  information** (or any other relevant chemical classification)
  **provided as metadata**. This structure information, stored in the
  metadata field “skeleton”, allows the generation of queries specific
  to a given skeleton by extracting repetitive skeleton-specific
  fragmentation patterns. The MIADB file is provided as an example.

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

To reproduce the example that uses the Monoterpene Indole Alkaloids
Database (.mgf) file by default, which includes the annotation of
spectral skeletons:

``` r
SpectraToQueries::spectra_to_queries()
```

To reproduce the “grouped” example that uses the MIADB file, which
includes an expert-based annotation of spectral “super skeletons”
(combination of skeletons exhibiting a high structural similarity):

``` r
SpectraToQueries::spectra_to_queries(
  spectra = system.file(
    "extdata",
    "spectra_grouped.rds",
    package = "SpectraToQueries"
  ),
  export = "data/interim/queries-grouped.tsv"
)
```

To generate diagnostic ions queries from your spectra:

``` r
SpectraToQueries::spectra_to_queries(
  spectra = "yourAwesomeSpectra.mgf",
  export = "path/yourEvenBetterResults.tsv"
)
```

Showing all parameters:

``` r
SpectraToQueries::spectra_to_queries(
  spectra = NULL,
  export = "data/interim/queries.tsv",
  beta_1 = 1.0,
  beta_2 = 0.5,
  dalton = 0.01,
  decimals = 4L,
  intensity_min = 0.0,
  ions_max = 10L,
  n_skel_min = 5L,
  n_spec_min = 3L,
  ppm = 30.0,
  fscore_min = 0.0,
  precision_min = 0.0,
  recall_min = 0.0,
  zero_val = 0.0
)
```

## Main Citations

Translating community-wide spectral library into actionable chemical
knowledge: a proof of concept with monoterpene indole alkaloids:
<https://doi.org/10.1186/s13321-025-01009-0>

## Additional software credits

| Package | Version | Citation |
|:---|:---|:---|
| base | 4.5.2 | R Core Team (2025) |
| BiocManager | 1.30.27 | Morgan and Ramos (2025) |
| BiocParallel | 1.44.0 | Wang et al. (2025) |
| BiocVersion | 3.22.0 | Morgan (2025) |
| knitr | 1.50 | Xie (2014); Xie (2015); Xie (2025) |
| MsBackendMgf | 1.18.0 | Gatto, Rainer, and Gibb (2025) |
| pkgload | 1.4.1 | Wickham et al. (2025) |
| progress | 1.2.3 | Csárdi and FitzJohn (2023) |
| rmarkdown | 2.30 | Xie, Allaire, and Grolemund (2018); Xie, Dervieux, and Riederer (2020); Allaire et al. (2025) |
| Spectra | 1.19.11 | Rainer et al. (2022) |
| testthat | 3.3.1 | Wickham (2011) |
| tidytable | 0.11.2 | Fairbanks (2024) |
| tidyverse | 2.0.0 | Wickham et al. (2019) |

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-rmarkdown2025" class="csl-entry">

Allaire, JJ, Yihui Xie, Christophe Dervieux, Jonathan McPherson, Javier
Luraschi, Kevin Ushey, Aron Atkins, et al. 2025.
*<span class="nocase">rmarkdown</span>: Dynamic Documents for r*.
<https://github.com/rstudio/rmarkdown>.

</div>

<div id="ref-progress" class="csl-entry">

Csárdi, Gábor, and Rich FitzJohn. 2023.
*<span class="nocase">progress</span>: Terminal Progress Bars*.
<https://github.com/r-lib/progress#readme>.

</div>

<div id="ref-tidytable" class="csl-entry">

Fairbanks, Mark. 2024. *<span class="nocase">tidytable</span>: Tidy
Interface to “<span class="nocase">data.table</span>”*.
<https://markfairbanks.github.io/tidytable/>.

</div>

<div id="ref-MsBackendMgf" class="csl-entry">

Gatto, Laurent, Johannes Rainer, and Sebastian Gibb. 2025.
*MsBackendMgf: Mass Spectrometry Data Backend for Mascot Generic Format
(Mgf) Files*. <https://doi.org/10.18129/B9.bioc.MsBackendMgf>.

</div>

<div id="ref-BiocVersion" class="csl-entry">

Morgan, Martin. 2025. *BiocVersion: Set the Appropriate Version of
Bioconductor Packages*. <https://doi.org/10.18129/B9.bioc.BiocVersion>.

</div>

<div id="ref-BiocManager" class="csl-entry">

Morgan, Martin, and Marcel Ramos. 2025. *BiocManager: Access the
Bioconductor Project Package Repository*.
<https://doi.org/10.32614/CRAN.package.BiocManager>.

</div>

<div id="ref-base" class="csl-entry">

R Core Team. 2025. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

<div id="ref-Spectra" class="csl-entry">

Rainer, Johannes, Andrea Vicini, Liesa Salzer, Jan Stanstrup, Josep M.
Badia, Steffen Neumann, Michael A. Stravs, et al. 2022. “A Modular and
Expandable Ecosystem for Metabolomics Data Annotation in r.”
*Metabolites* 12: 173. <https://doi.org/10.3390/metabo12020173>.

</div>

<div id="ref-BiocParallel" class="csl-entry">

Wang, Jiefei, Martin Morgan, Valerie Obenchain, Michel Lang, Ryan
Thompson, and Nitesh Turaga. 2025. *BiocParallel: Bioconductor
Facilities for Parallel Evaluation*.
<https://doi.org/10.18129/B9.bioc.BiocParallel>.

</div>

<div id="ref-testthat" class="csl-entry">

Wickham, Hadley. 2011. “<span class="nocase">testthat</span>: Get
Started with Testing.” *The R Journal* 3: 5–10.
<https://journal.r-project.org/articles/RJ-2011-002/>.

</div>

<div id="ref-tidyverse" class="csl-entry">

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019.
“Welcome to the <span class="nocase">tidyverse</span>.” *Journal of Open
Source Software* 4 (43): 1686. <https://doi.org/10.21105/joss.01686>.

</div>

<div id="ref-pkgload" class="csl-entry">

Wickham, Hadley, Winston Chang, Jim Hester, and Lionel Henry. 2025.
*<span class="nocase">pkgload</span>: Simulate Package Installation and
Attach*. <https://github.com/r-lib/pkgload>.

</div>

<div id="ref-knitr2014" class="csl-entry">

Xie, Yihui. 2014. “<span class="nocase">knitr</span>: A Comprehensive
Tool for Reproducible Research in R.” In *Implementing Reproducible
Computational Research*, edited by Victoria Stodden, Friedrich Leisch,
and Roger D. Peng. Chapman; Hall/CRC.

</div>

<div id="ref-knitr2015" class="csl-entry">

———. 2015. *Dynamic Documents with R and Knitr*. 2nd ed. Boca Raton,
Florida: Chapman; Hall/CRC. <https://yihui.org/knitr/>.

</div>

<div id="ref-knitr2025" class="csl-entry">

———. 2025. *<span class="nocase">knitr</span>: A General-Purpose Package
for Dynamic Report Generation in R*. <https://yihui.org/knitr/>.

</div>

<div id="ref-rmarkdown2018" class="csl-entry">

Xie, Yihui, J. J. Allaire, and Garrett Grolemund. 2018. *R Markdown: The
Definitive Guide*. Boca Raton, Florida: Chapman; Hall/CRC.
<https://bookdown.org/yihui/rmarkdown>.

</div>

<div id="ref-rmarkdown2020" class="csl-entry">

Xie, Yihui, Christophe Dervieux, and Emily Riederer. 2020. *R Markdown
Cookbook*. Boca Raton, Florida: Chapman; Hall/CRC.
<https://bookdown.org/yihui/rmarkdown-cookbook>.

</div>

</div>
