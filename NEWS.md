# SpectraToQueries

# SpectraToQueries 0.0.9002

## Data & Documentation

- **Added clean `data/` objects** for standard spectra:
  - `spectra_grouped` (352KB): Example dataset, 321 grouped Spectra from MIADB,
    53 unique skeletons
  - `spectra` (354KB): Ungrouped variant for comparison
  - Both fully documented with `@format`, `@source`, and usage examples
  - Load via `data(spectra_grouped)` or `data(spectra)` instead of `readRDS()` +
    `system.file()`

- **`spectra_to_queries()` API improvements**:
  - `spectra = NULL` now loads `spectra_grouped` automatically (cleaner default)
  - `spectra = "grouped"` alias for explicit grouped-data loading
  - Accepts file paths and Spectra objects as before
  - Updated roxygen with concrete usage examples

- **Cleaned `R/data.R`**:
  - Removed all `inst/extdata/` system paths from examples
  - Added concrete example showing `data(spectra_grouped)` usage
  - Marked examples as `\dontrun` to avoid long test times in check

## Dependencies

- All fastverse toolkit already declared:
  - `fastmatch` (set membership) in `combine_ions_minimal.R`
  - `stringi` (fixed-string operations) across multiple functions
  - `cheapr` (NA handling) in `perform_list_of_queries.R`
  - `tidytable` (grouping/summarization) now used in `fix_binned_mzs.R`
- Removed direct `data.table` dependency (tidytable already in imports)

# SpectraToQueries 0.0.9001

- Updated minimal R version to `4.4.0` (and related Bioconductor dependencies)

# SpectraToQueries 0.0.9000

- Initial version.
