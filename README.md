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

`SpectraToQueries` discovers diagnostic MS/MS fragment (`_frag`) and
neutral-loss (`_nl`) ion combinations from annotated tandem mass spectrometry
libraries and translates them into
[MassQL](https://mwang87.github.io/MassQueryLanguage_Documentation/) queries.
Candidate queries for each chemical class are ranked by Matthews Correlation
Coefficient (MCC), selecting the smallest ion combination tied at the highest
score.

## Requirements

- **R (\>= 4.5.0)** with `MsBackendMgf` and `Spectra`.
- **MS2 spectra with class metadata**: an `.mgf` file, a `Spectra` object, or a
  compatible `data.frame` containing MS2 scans, `precursorMz` values, and a
  `SKELETON` metadata column (lowercase `skeleton` is accepted as a fallback).
  Skeleton classes with fewer than `n_skel_min` non-empty MS2 spectra are
  skipped.

## Installation

`SpectraToQueries` is distributed via r-universe:

```r
install.packages(
  "SpectraToQueries",
  repos = c(
    "https://spectra-to-knowledge.r-universe.dev",
    "https://cran.r-universe.dev",
    "https://bioc.r-universe.dev"
  )
)
```

## Usage

Run on the bundled raw Monoterpene Indole Alkaloids Database (`mia_spectra_df`):

```r
SpectraToQueries::spectra_to_queries(export = NULL)
```

Run on the bundled grouped skeleton dataset (`mia_spectra_grouped_df`):

```r
SpectraToQueries::spectra_to_queries(spectra = "grouped", export = NULL)
```

Generate queries from your own `.mgf` library and save the TSV output:

```r
SpectraToQueries::spectra_to_queries(
  spectra = "path/to/library.mgf",
  export = "results/queries.tsv"
)
```

Pass `export = NULL` to return the `data.frame` without writing a file to disk.

### Parameters

```r
SpectraToQueries::spectra_to_queries(
  spectra = NULL,
  export = "data/interim/queries.tsv",
  dalton = 0.01,
  decimals = 4L,
  intensity_min = 0.0,
  ions_max = 10L,
  n_skel_min = 5L,
  n_spec_min = 3L,
  ppm = 30.0,
  mcc_min = 0.0,
  precision_min = 0.0,
  recall_min = 0.0,
  zero_val = 0.0
)
```

### Output, MassQL translation, and known behavior

`spectra_to_queries()` returns (and optionally exports as TSV) one row per
qualifying skeleton class:

- `skeleton`: chemical class identifier.
- `mcc`: Matthews Correlation Coefficient on the input MS2 library.
- `combined_ions`: minimal conjunctive/disjunctive combination of fragments
  (`_frag`) and neutral losses (`_nl`).
- `massql`: executable
  [MassQL](https://mwang87.github.io/MassQueryLanguage_Documentation/) query
  (`QUERY scaninfo(MS2DATA) WHERE ...`).

You can also convert ion expressions or vectors directly with `to_massql()`:

```r
SpectraToQueries::to_massql(
  "108.0804_frag & (241.095_nl | 144.0807_frag)",
  dalton = 0.01,
  ppm = 30.0
)
```

Notes on query semantics:

- MassQL ignores `TOLERANCEMZ` when `TOLERANCEPPM` is present on the same
  condition. When `dalton > 0`, `to_massql()` emits a single `:TOLERANCEMZ=`
  modifier (`round(dalton + max(mz) * ppm / 1e6, 6)` on `MS2PROD` and static
  `dalton` on `MS2NL`). When `dalton == 0` and `ppm > 0`, it emits
  `:TOLERANCEPPM=`.
- Disjunctive clauses in `combined_ions` are expanded to disjunctive normal form
  and separated by `" ; "` in `massql`.
- Peak matching evaluates binary peak presence after relative base-peak
  normalization (`0` to `100`) and `intensity_min` filtering; relative intensity
  ratios between fragments are not encoded in the queries.
- `+` characters in `SKELETON` labels are replaced with `.` for internal
  grouping (a warning is emitted if distinct labels collide).

### Example data

The package ships two MIADB datasets as portable `data.frame` objects (with `mz`
and `intensity` list columns):

- `mia_spectra_df`: 321 MS/MS spectra across 81 skeletons (24,255 peaks total).
- `mia_spectra_grouped_df`: 321 MS/MS spectra across 53 grouped skeleton
  classes.

```r
utils::data(mia_spectra_df, package = "SpectraToQueries")
utils::data(mia_spectra_grouped_df, package = "SpectraToQueries")
```

## Contributing

Bug reports, tests, and pull requests are welcome. See
[CONTRIBUTING.md](CONTRIBUTING.md) for local setup and check commands, and
[TESTING.md](TESTING.md) for testing conventions and error-handling contracts.

## Main Citation

Translating community-wide spectral library into actionable chemical knowledge:
a proof of concept with monoterpene indole alkaloids:
<https://doi.org/10.1186/s13321-025-01009-0>

## Additional software credits

  | Package          | Version    | Citation                                                                 |
  | :---             | :---       | :---                                                                     |
  | base             | 4.6.2      | R Core Team (2026)                                                       |
  | BiocManager      | 1.30.27    | Morgan and Ramos (2025)                                                  |
  | BiocVersion      | 3.23.1     | Morgan (2025)                                                            |
  | grateful         | 0.3.1      | Rodriguez-Sanchez and Jackson (2025)                                     |
  | MsBackendMgf     | 1.20.1     | Gatto et al. (2026)                                                      |
  | Spectra          | 1.22.2     | Rainer et al. (2022)                                                     |
  | SpectraToQueries | 0.0.0.9002 | Rutz and Szwarc (2025); <span class="nocase">Szwarc et al.</span> (2025) |

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-MsBackendMgf" class="csl-entry">

Gatto, Laurent, Johannes Rainer, and Sebastian Gibb. 2026.
*MsBackendMgf: Mass Spectrometry Data Backend for Mascot Generic Format
(Mgf) Files*. <https://github.com/RforMassSpectrometry/MsBackendMgf>.

</div>

<div id="ref-BiocVersion" class="csl-entry">

Morgan, Martin. 2025. *BiocVersion: Set the Appropriate Version of
Bioconductor Packages*.

</div>

<div id="ref-BiocManager" class="csl-entry">

Morgan, Martin, and Marcel Ramos. 2025. *BiocManager: Access the
Bioconductor Project Package Repository*.
<https://doi.org/10.32614/CRAN.package.BiocManager>.

</div>

<div id="ref-base" class="csl-entry">

R Core Team. 2026. *R: A Language and Environment for Statistical
Computing*. R Foundation for Statistical Computing.
<https://doi.org/10.32614/R.manuals>.

</div>

<div id="ref-Spectra" class="csl-entry">

Rainer, Johannes, Andrea Vicini, Liesa Salzer, et al. 2022. “A Modular
and Expandable Ecosystem for Metabolomics Data Annotation in r.”
*Metabolites* 12: 173. <https://doi.org/10.3390/metabo12020173>.

</div>

<div id="ref-grateful" class="csl-entry">

Rodriguez-Sanchez, Francisco, and Connor P. Jackson. 2025.
*<span class="nocase">grateful</span>: Facilitate Citation of R
Packages*. <https://pakillo.github.io/grateful/>.

</div>

<div id="ref-SpectraToQueries2025a" class="csl-entry">

Rutz, Adriano, and Sarah Szwarc. 2025. *SpectraToQueries: Spectra to
Queries*.

</div>

<div id="ref-SpectraToQueries2025b" class="csl-entry">

<span class="nocase">Szwarc, Sarah, Adriano Rutz, Kyungha Lee, et
al.</span> 2025. “Translating Community-Wide Spectral Library into
Actionable Chemical Knowledge: A Proof of Concept with Monoterpene
Indole Alkaloids.” *J. Cheminform.* 17.
<https://doi.org/10.1186/s13321-025-01009-0>.

</div>

</div>
