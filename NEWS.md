# SpectraToQueries

# SpectraToQueries 0.0.0.9002

## Performance & Optimization

- **~30x end-to-end speedup** in `spectra_to_queries()` with 100% identical
  outputs:
  - Precomputed single-ion logical hit vectors in
    `perform_list_of_queries_progress()` and converted S4 `SimpleNumericList`
    m/z lists to plain R vectors once in `prepare_query_data()`.
  - Replaced nested loop linear scans in `harmonize_mzs()` with vectorized
    `findInterval()` nearest-neighbor lookup.
  - Switched `Spectra::bin()` to sparse mode (`zero.rm = TRUE`) with
    sparse-aware `create_matrix()` and column-major distance calculation in
    `fix_binned_mzs()`.
  - Avoided redundant `Spectra::combineSpectra()` and duplicate lazy queue
    evaluation before `Spectra::applyProcessing()`.

## Dependencies

- Removed `tidytable` from `Imports`, leaving only `MsBackendMgf`, `Spectra`,
  and Base R `utils`.
- Removed `fastmatch`, `stringi`, `BiocParallel`, `data.table`, `progress`,
  `purrr`, and `tibble`.

## Data, Documentation & Metadata

- **Portable `data/` datasets** regenerated from source MGF files:
  - `mia_spectra_df`: 321 MS/MS spectra from MIADB across 81 chemical skeletons
    (24,255 total peaks, median 66 peaks per spectrum).
  - `mia_spectra_grouped_df`: 321 MS/MS spectra from MIADB across 53 grouped
    skeleton classes.
- Regenerated `data/interim/queries.tsv` and `data/interim/queries-grouped.tsv`
  using Matthews Correlation Coefficient (MCC) ranking.
- Added verified EDAM ontology terms (`X-EDAM-Topic`, `X-EDAM-Operation`,
  `X-EDAM-Data`, `X-EDAM-Format`, `X-schema.org-keywords`, `biocViews`) to
  `DESCRIPTION`.

# SpectraToQueries 0.0.9001

- Updated minimal R version to `4.4.0` (and related Bioconductor dependencies)

# SpectraToQueries 0.0.9000

- Initial version.
