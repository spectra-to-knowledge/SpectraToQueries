# Testing Guide for `SpectraToQueries`

This document describes the test structure, error-handling contracts, and
commands used to verify changes in `SpectraToQueries`.

## Test Suite Layout

Tests use `testthat` edition 3 and live under `tests/testthat/`:

- `tests/testthat/test-architecture.R`: Structural and compatibility checks
  covering namespace exports (`spectra_to_queries`, `to_massql`), MIADB subset
  characterization results, absence of internal `@backend@` slot access, public
  `.Rd` inventory in `man/`, clean `basename()` values in bundled `dataOrigin`
  columns, and `DESCRIPTION` `Collate`/`Imports`/`Suggests` synchronization.
- `tests/testthat/test-functions.R`: Unit, boundary, and end-to-end tests across
  preprocessing (`normalize_peaks`, `harmonize_mzs`, `fix_binned_mzs`,
  `create_matrix`, `filter_matrix`), candidate selection and Branch-and-Bound
  search, Boolean minimization (`combine_ions_minimal`), MassQL translation
  (`to_massql`), input validation, non-`MsBackendMgf` backends
  (`MsBackendDataFrame`), and full `spectra_to_queries()` execution.

## Testing Conventions

1. **Assertion of Independently Expected Results**
   - Assert expected values, types, column names, and ordering explicitly with
     `testthat::expect_identical()` or
     `testthat::expect_equal(..., tolerance = ...)`.
   - Do not compute expected test values by calling the same function under
     test.

2. **Tolerance Boundary Checks**
   - Test numeric thresholds (`dalton`, `ppm`, `intensity_min`, `mcc_min`,
     `precision_min`, `recall_min`) at both sides of the boundary and at exact
     equality where applicable.
   - Prefer exact binary floating-point values (such as `0.25` and `0.5`) when
     testing floating-point boundaries around `dalton + mz * ppm / 1e6`.

3. **Classed Error Conditions**
   - Input validation errors are raised via the internal `abort_input()` helper
     in `R/spectra_utils.R` with condition class `"stq_input_error"`.
   - Missing input file paths in `load_input_spectra()` raise class
     `"stq_file_not_found_error"`.
   - Assert both condition class and message text:

   ```r
   testthat::expect_error(
     SpectraToQueries::to_massql("100_frag", dalton = -1),
     regexp = "`dalton` must be a single finite non-negative number",
     class = "stq_input_error"
   )
   ```

4. **Temporary Files and Cleanup**
   - Write temporary `.mgf` or `.tsv` files under `tempdir()` using `tempfile()`
     and clean them up with `on.exit(unlink(...), add = TRUE)`.

5. **Regression Tests Before Fixes**
   - Add a failing unit or boundary test reproducing the defect before applying
     a code change, then verify the fix makes the test pass without breaking
     characterization checks.

## Running Tests, Coverage, and Static Checks Locally

```bash
# Run the full testthat suite against the local source tree
Rscript --vanilla -e 'testthat::test_local(stop_on_failure = TRUE)'

# Compute line coverage (CI enforces >= 95%)
Rscript --vanilla -e '
cov <- covr::package_coverage(quiet = FALSE)
print(cov)
covr::zero_coverage(cov)
pct <- covr::percent_coverage(cov)
if (pct < 95) stop(sprintf("Coverage %.2f%% is below 95%%", pct))
'

# Lint and spell check (load workspace first so lintr resolves internal helpers)
Rscript --vanilla -e '
pkgload::load_all(".", quiet = TRUE)
lints <- lintr::lint_package()
print(lints)
if (length(lints) > 0L) stop("Lint errors found.")
sp <- spelling::spell_check_package()
print(sp)
if (nrow(sp) > 0L) stop("Spelling errors found.")
'

# Formatting, static analysis, and pre-commit hooks
air format --check .
jarl check .
tombi lint .
panache lint .
prek run --all-files

# Clean tarball build and R CMD check
R CMD build .
R CMD check --no-manual SpectraToQueries_*.tar.gz
```
