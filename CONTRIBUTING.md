# Contributing to `SpectraToQueries`

Thank you for contributing to `SpectraToQueries`. This guide describes the local
development tools, repository checks, and pull-request expectations.

## Prerequisites

- **R (>= 4.5.0)** with runtime and development packages installed:

```r
install.packages(
  c(
    "MsBackendMgf",
    "Spectra",
    "cffr",
    "covr",
    "grateful",
    "knitr",
    "lintr",
    "pkgload",
    "roxygen2",
    "spelling",
    "testthat"
  ),
  repos = c(
    "https://cran.r-universe.dev",
    "https://bioc.r-universe.dev"
  )
)
```

- **Command-line tools** used by pre-commit hooks and documentation builds:
  - [air](https://posit-dev.github.io/air/) (R formatting)
  - [jarl](https://jarl.etiennebacher.com/) (R static analysis)
  - [tombi](https://tombi-toml.github.io/tombi/) (TOML formatting and linting)
  - [panache](https://github.com/jolars/panache) (Markdown formatting and
    linting)
  - [prek](https://prek.j178.dev/) (pre-commit hook runner)
  - [Quarto](https://quarto.org/) (rendering `README.qmd` and vignettes)

## Setting up git hooks

Install the hooks defined in `prek.toml` once after cloning:

```bash
prek install
prek run --all-files
```

## Development workflow

Run formatting, static analysis, spell checking, tests, and package checks
before opening a pull request:

```bash
# Format code, TOML, and Markdown
air format .
tombi format prek.toml
panache format .

# Verify formatting and static lints
air format --check .
jarl check .
tombi lint .
panache --no-cache format --check .
panache lint .
prek run --all-files

# Load the workspace namespace so lintr sees development helpers
Rscript --vanilla -e '
pkgload::load_all(".", quiet = TRUE)
stopifnot(length(lintr::lint_package()) == 0L)
stopifnot(nrow(spelling::spell_check_package()) == 0L)
'

# Run the test suite
Rscript --vanilla -e 'testthat::test_local(stop_on_failure = TRUE)'

# Build and run R CMD check
R CMD build .
R CMD check --no-manual SpectraToQueries_*.tar.gz
```

See [TESTING.md](TESTING.md) for testing conventions, error condition classes,
and coverage verification (`>= 95%` line coverage enforced in CI).

### Updating roxygen documentation, README, and citations

- Exported functions (`spectra_to_queries`, `to_massql`) and datasets
  (`mia_spectra_df`, `mia_spectra_grouped_df`) generate `.Rd` files in `man/`.
  Mark all internal helpers `@keywords internal` and `@noRd` so internal pages
  are not emitted. Regenerate `man/` and `NAMESPACE` with:

```r
roxygen2::roxygenise()
```

- If you edit `README.qmd`, install the current workspace into a temporary
  library before rendering so Quarto executes against your local changes rather
  than an older system installation:

```bash
TMP_LIB="$(mktemp -d)"
R CMD INSTALL --no-multiarch -l "$TMP_LIB" .
R_LIBS_USER="$TMP_LIB" R_LIBS="$TMP_LIB" quarto render README.qmd --to gfm
rm -rf "$TMP_LIB" README.html README_files
```

- When `DESCRIPTION` or `inst/CITATION` changes, refresh `CITATION.cff`:

```bash
Rscript --vanilla -e '
cffr::cff_write(keys = list(repository = "https://spectra-to-knowledge.r-universe.dev"))
'
```

- Genuine technical or chemical terms flagged by
  `spelling::spell_check_package()` belong in `inst/WORDLIST` (one word per
  line, sorted).

## Codebase conventions

- Keep runtime `Imports` minimal (`MsBackendMgf`, `Spectra`, `stats`, `utils`).
- Keep `Collate` and `Imports`/`Suggests` in `DESCRIPTION` synchronized with
  files in `R/` and `pkg::` calls (enforced by
  `tests/testthat/test-architecture.R`).
- Interact with `Spectra` objects strictly through public `Spectra` methods,
  never internal `@backend@` slots, so alternative backends (`MsBackendMemory`,
  `MsBackendDataFrame`) work without errors.
- Avoid leaving `browser()`, `debug()`, `debugonce()`, or bare `print()` calls
  in `R/`. Use `message()` or `warning()` when user-facing diagnostics are
  needed.
- Use concise semantic commit messages in the imperative mood (for example,
  `fix: ...`, `test: ...`, `docs: ...`) under 72 characters.
