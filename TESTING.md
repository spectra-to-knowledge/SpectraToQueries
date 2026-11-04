# Testing Guide for `SpectraToQueries`

This document describes the testing conventions, error-handling contracts, and local commands for running tests, coverage, and mutation checks in `SpectraToQueries`.

## Testing Conventions

1. **Mutation-Resistant Assertions Over Coverage-Only Execution**
   - Every test must assert an independently computed expected outcome (`expect_identical()` or `expect_equal(..., tolerance = ...)`).
   - Check types, dimensions, row/column names, and ordering in addition to numeric values.
   - Never compute expected values by re-running the function under test.

2. **Exact Boundary Testing**
   - Every relational operator (`<`, `<=`, `>`, `>=`, `==`, `!=`) must be tested at `boundary - eps`, `boundary`, and `boundary + eps`.
   - When testing floating-point tolerance boundaries (`dalton`, `ppm`), prefer exact binary floating-point values (e.g., `0.25`, `0.5`) so boundary equality (`diff == tol`) is exact in IEEE 754 arithmetic.

3. **Classed Error Conditions**
   - Input validation errors use `cli::cli_abort(..., class = "stq_input_error")` via the internal `abort_bad_arg()` helper in `R/spectra_utils.R`.
   - Missing file errors use `class = "stq_file_not_found_error"`.
   - Always assert both the condition class and message pattern:
     ```r
     testthat::expect_error(
       normalize_peaks(c(1, 2, 3)),
       regexp = "`x` must be a 2-column numeric matrix",
       class = "stq_input_error"
     )
     ```

4. **Filesystem Hygiene**
   - Always write temporary files/directories under `tempdir()` (`tempfile()`) and register cleanup via `on.exit(unlink(...), add = TRUE)`.

5. **Separate Bug-Fix Commits**
   - When fixing a bug, commit a failing test (`test: ...`) before committing the fix (`fix: ...`).

## Running Tests & Quality Checks Locally

```bash
# Run the testthat suite
Rscript -e "devtools::test()"

# Compute line and function coverage (target: >= 90% floor; currently 100%)
Rscript -e "cov <- covr::package_coverage(); print(cov); covr::zero_coverage(cov)"

# Run static analysis and formatting checks
Rscript -e "lintr::lint_package()"
air format --check .
jarl check .

# Run full R CMD check
Rscript -e "devtools::check()"
```

## Running Mutation Checks Locally

To verify that tests catch small semantic regressions (operator flips, constant changes, branch removals), swap modified function bodies into both the package namespace and attached package environment and verify that `testthat::test_file()` fails:

```r
pkgload::load_all(".", quiet = TRUE)

check_mutant <- function(fn_name, pattern, replacement) {
  ns <- asNamespace("SpectraToQueries")
  pkg_env <- as.environment("package:SpectraToQueries")
  orig_fn <- get(fn_name, envir = ns)
  on.exit({
    unlockBinding(fn_name, ns)
    assign(fn_name, orig_fn, envir = ns)
    lockBinding(fn_name, ns)
    if (exists(fn_name, envir = pkg_env, inherits = FALSE)) {
      unlockBinding(fn_name, pkg_env)
      assign(fn_name, orig_fn, envir = pkg_env)
      lockBinding(fn_name, pkg_env)
    }
  }, add = TRUE)

  src <- paste(deparse(orig_fn, width.cutoff = 500L), collapse = "\n")
  mut_fn <- eval(parse(text = sub(pattern, replacement, src, fixed = TRUE)), envir = ns)

  unlockBinding(fn_name, ns)
  assign(fn_name, mut_fn, envir = ns)
  lockBinding(fn_name, ns)
  if (exists(fn_name, envir = pkg_env, inherits = FALSE)) {
    unlockBinding(fn_name, pkg_env)
    assign(fn_name, mut_fn, envir = pkg_env)
    lockBinding(fn_name, pkg_env)
  }

  out <- as.data.frame(testthat::test_file(
    "tests/testthat/test-functions.R",
    reporter = "silent",
    env = new.env(parent = ns),
    load_package = "none"
  ))
  sum(out$failed) > 0L || any(out$error)
}
```
