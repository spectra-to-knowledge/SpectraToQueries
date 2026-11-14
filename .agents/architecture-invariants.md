# Architecture Invariants

Every invariant below is tested in `tests/testthat/test-architecture.R` or
`tests/testthat/test-functions.R`. Do not weaken tests to pass a change.

## 1. Public Surface and Roxygen Hygiene

- `NAMESPACE` exports only `spectra_to_queries` and `to_massql`.
- `man/` must contain exactly 5 `.Rd` files: `SpectraToQueries-package.Rd`,
  `mia_spectra_df.Rd`, `mia_spectra_grouped_df.Rd`, `spectra_to_queries.Rd`,
  `to_massql.Rd`.
- Mark every non-exported helper `@keywords internal` and `@noRd`. Do not emit
  internal `.Rd` pages.

## 2. DESCRIPTION Sync

- `Collate` in `DESCRIPTION` must match the exact set of `*.R` files in `R/`
  (`test-architecture.R`).
- Every `pkg::fn` reference in `R/` must appear in `Imports`; every `pkg::fn`
  reference in `tests/testthat/` must appear in `Imports` or `Suggests`.
- Keep `.Rbuildignore` sorted in C locale (`prek` `file-contents-sorter` hook).

## 3. Spectra Backend Independence & S4 Error Propagation

- Never access `@backend@` S4 slots anywhere in `R/`. Use public `Spectra`
  getters/methods (`spectraData`, `peaksData`, `precursorMz`, `msLevel`,
  `filterMsLevel`, `filterEmptySpectra`, `addProcessing`, `applyProcessing`,
  `neutralLoss`, `PrecursorMzParam`).
- Always evaluate input/file checks (raising `stq_input_error` or
  `stq_file_not_found_error`) into a variable **before** passing the `Spectra`
  object to an S4 generic such as `Spectra::filterMsLevel()`, because S4 generic
  dispatch wraps conditions thrown from lazy argument evaluation into
  `simpleError`.
- Keep `MsBackendDataFrame` and `MsBackendMemory` tests passing alongside `.mgf`
  inputs.
- Always check `"COL" %in% Spectra::spectraVariables(x)` before accessing
  optional metadata columns with `x$COL` (such as `TITLE`, `PRECURSOR_MZ`, or
  `skeleton`), because `Spectra` `$` throws an error on missing variables
  instead of returning `NULL`.

## 4. Label Handling and Intensity Filtering

- Flatten named character lists with `rep(names(x), lengths(x))`, never
  `names(unlist(x))` (which corrupts skeleton names ending in digits by
  appending occurrence suffixes).
- Normalize `+` to `.` in `SKELETON` labels only with an explicit collision
  warning (`Distinct SKELETON labels collide after '+' to '.' normalization`).
- `intensity_min` is validated in `[0, 100]` and enforced relative to base peak
  (`0..100`) during `prepare_query_data()`.
- Filter to `msLevel == 2L` and non-empty `SKELETON` metadata before computing
  `total_spectra` or MCC denominators. Accept lowercase `skeleton` as fallback
  in `load_input_spectra()`.

## 5. Boolean Precedence and MassQL Semantics (`to_massql`, `combine_ions_minimal`)

- Match `_frag` and `_nl` suffixes strictly at token end (`endsWith()`), never
  by unanchored substring search.
- Parenthesize multi-term OR clauses produced by `combine_ions_minimal()` so
  `A & (B | C)` never collapses into `A & B | C`.
- Parse expressions in `to_massql()` using `parse_massql_dnf()` (explicit env
  cursor state, no `<<-` superassignment), distribute to DNF, apply superset
  absorption, and join independent OR clauses with `" ; "`.
- Upstream MassQL ignores `TOLERANCEMZ` whenever `TOLERANCEPPM` is present on
  the same condition. Therefore emit a single effective tolerance modifier per
  condition:
  - When `dalton > 0`: `:TOLERANCEMZ=round(dalton + max(mz) * ppm / 1e6, 6L)` on
    `MS2PROD` and static `:TOLERANCEMZ=dalton` on `MS2NL`.
  - When `dalton == 0 && ppm > 0`: `:TOLERANCEPPM=ppm`.
  - Omit tolerance modifier when both `dalton == 0` and `ppm == 0`.
- Format numeric literals via `fmt_num()` (`scientific = FALSE`, `trim = TRUE`,
  `digits = 15L`) so values like `1e-04` never enter MassQL queries.

## 6. Dataset Provenance Hygiene

- `spectra_to_df()` in `R/spectra_utils.R` must strip directory prefixes from
  non-`NA` `dataOrigin` values with `basename()`.
- Bundled `mia_spectra_df` and `mia_spectra_grouped_df` must never contain `/`
  or `\` in `dataOrigin` (`test-architecture.R`).
