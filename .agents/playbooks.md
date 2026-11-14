# Agent Playbooks

## 1. Fixing a Bug or Changing Pipeline Logic

1. Add a failing test in `tests/testthat/test-functions.R` covering exact
   boundary values and error classes (`stq_input_error` /
   `stq_file_not_found_error`).
2. Run only the targeted test file:
   `Rscript --vanilla -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-functions.R")'`
3. Implement the smallest clean fix in `R/`, reusing `abort_input()`,
   `check_num()`, `calc_mcc_vec()`, `empty_query_table()`, and
   `write_query_table()`.
4. If the change intentionally alters query output formatting or values,
   regenerate checked-in reference TSVs and inspect the diff carefully:
   - `Rscript --vanilla -e 'pkgload::load_all(".", quiet = TRUE); spectra_to_queries(export = "data/interim/queries.tsv"); spectra_to_queries(spectra = "grouped", export = "data/interim/queries-grouped.tsv")'`
5. Format with `air format .`, check `jarl check .`, load-check `lintr`, run
   `testthat::test_local(stop_on_failure = TRUE)`, and run
   `prek run --all-files`.

## 2. Adding, Renaming, or Removing `R/*.R` Files or Roxygen Tags

1. Update `Collate` in `DESCRIPTION` and `Imports`/`Suggests` if needed.
2. Mark all non-exported functions `@keywords internal` and `@noRd`.
3. Run `Rscript --vanilla -e 'roxygen2::roxygenise()'` and confirm zero warnings
   and no new `.Rd` files beyond the 5 allowed in `man/`.
4. Run
   `Rscript --vanilla -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-architecture.R", stop_on_failure = TRUE)'`.

## 3. Regenerating Bundled MIADB Datasets (`data/*.rda`)

1. Run `Rscript --vanilla data-raw/generate_spectra_data.R`.
2. Verify `data/mia_spectra_df.rda` and `data/mia_spectra_grouped_df.rda` load
   cleanly and contain zero `/` or `\` in `dataOrigin` by running
   `testthat::test_file("tests/testthat/test-architecture.R")`.

## 4. Rendering `README.qmd`, Vignettes, and `CITATION.cff`

The host R library may contain an older installation of `SpectraToQueries`. Do
not modify the host system library; install into a temporary directory before
running Quarto:

```bash
TMP_LIB="$(mktemp -d)"
R CMD INSTALL --no-multiarch -l "$TMP_LIB" .
R_LIBS_USER="$TMP_LIB" R_LIBS="$TMP_LIB" quarto render README.qmd --to gfm
TMP_VIG="$(mktemp -d)"
R_LIBS_USER="$TMP_LIB" R_LIBS="$TMP_LIB" quarto render vignettes/spectra-to-queries.qmd --output-dir "$TMP_VIG"
rm -rf "$TMP_LIB" "$TMP_VIG" README.html README_files vignettes/.gitignore
```

When `DESCRIPTION` or `inst/CITATION` changes, regenerate `CITATION.cff` with
the explicit r-universe repository key so `cffr` does not default to
`bioconductor.org`:

```bash
Rscript --vanilla -e 'cffr::cff_write(keys = list(repository = "https://spectra-to-knowledge.r-universe.dev"))'
```

Then run `prek run --files CITATION.cff` (twice if hooks trim trailing spaces).

## 5. Definition of Done Before Committing

- [ ] `air format --check .` passes
- [ ] `jarl check .` passes
- [ ] `tombi lint .` and `panache --no-cache format --check . && panache lint .`
  pass
- [ ] `pkgload::load_all(".", quiet = TRUE); lintr::lint_package()` returns 0
  lints
- [ ] `spelling::spell_check_package()` returns 0 rows
- [ ] `testthat::test_local(stop_on_failure = TRUE)` passes with 0
  fails/warnings/skips
- [ ] `covr::percent_coverage(covr::package_coverage()) >= 95`
- [ ] `prek run --all-files` passes on clean second run
- [ ] Clean archive check passes (`git archive HEAD` into temp dir ->
  `R CMD build .` -> `R CMD check --no-manual` with 0 errors and 0 warnings)
- [ ] Agent-facing files (`AGENTS.md`, `.agents/`) are committed separately from
  code and human docs, with semantic commit subject (< 72 chars) and zero AI
  attribution trailers
