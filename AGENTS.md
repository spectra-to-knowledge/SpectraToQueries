# AGENTS.md

Instructions for coding agents working in `SpectraToQueries` (`R >= 4.5.0`).
Human contributor and testing guides live in [CONTRIBUTING.md](CONTRIBUTING.md)
and [TESTING.md](TESTING.md). Keep human-facing files independent: never
reference `AGENTS.md` or `.agents/` from human docs, `R/`, `man/`, or
`vignettes/`.

## Verification Commands (Cheapest First)

1. Static formatting and linting (CLI, fast):
   - `air format --check .`
   - `jarl check .`
   - `tombi lint .`
   - `panache --no-cache format --check . && panache lint .`
2. R package lint and spelling (load workspace first so `lintr` sees internal
   helpers):
   - `Rscript --vanilla -e 'pkgload::load_all(".", quiet = TRUE); stopifnot(length(lintr::lint_package()) == 0L, nrow(spelling::spell_check_package()) == 0L)'`
3. Targeted single test file:
   - `Rscript --vanilla -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-architecture.R", stop_on_failure = TRUE)'`
4. Full test suite and pre-commit gate:
   - `Rscript --vanilla -e 'testthat::test_local(stop_on_failure = TRUE)'`
   - `prek run --all-files`
5. Coverage and tarball check (before finishing substantial work):
   - `Rscript --vanilla -e 'cov <- covr::package_coverage(); stopifnot(covr::percent_coverage(cov) >= 95)'`
   - `R CMD build . && R CMD check --no-manual SpectraToQueries_*.tar.gz`

## Repository Map

- `R/`: Package source (`spectra_to_queries.R`, `combine_ions_minimal.R`,
  `harmonize_mzs.R`, `fix_binned_mzs.R`, `create_matrix.R`, `filter_matrix.R`,
  `normalize_peaks.R`, `perform_list_of_queries.R`, `create_dir.R`,
  `spectra_utils.R`, `data.R`, `SpectraToQueries-package.R`).
- `tests/testthat/`: `test-architecture.R` (contracts and characterization) and
  `test-functions.R` (unit, boundary, backend, and pipeline tests).
- `data/` & `data-raw/`: `mia_spectra_df.rda`, `mia_spectra_grouped_df.rda`, and
  `data-raw/generate_spectra_data.R` reading `data/source/spectra/*.mgf`.
- `data/interim/`: Checked-in reference outputs `queries.tsv` and
  `queries-grouped.tsv`.
- `man/`: Public roxygen `.Rd` pages only (5 `.Rd` files) plus
  `man/figures/logo.svg`.
- `vignettes/`: `spectra-to-queries.qmd` (excluded from tarball via
  `.Rbuildignore`, built by `altdoc` on CI).
- `.github/workflows/`: CI for `R-CMD-check`, `format-check`, `test-coverage`
  (`>= 95%`), `external-use`, `altdoc`, and `update-citation-cff`.

## Non-Negotiable Rules

- Never add unrequested runtime dependencies (`Imports` is locked to
  `MsBackendMgf`, `Spectra`, `stats`, `utils`).
- Never leave `browser()`, `debug()`, `debugonce()`, or `print()` calls in `R/`
  (enforced by `prek.toml`).
- Always follow [architecture invariants](.agents/architecture-invariants.md)
  and [playbooks](.agents/playbooks.md).
- Never commit `.agent-work/`, build tarballs, `.Rcheck/`, or temporary
  artifacts. Do not push unless explicitly asked.
