# SpectraToQueries

# SpectraToQueries 0.0.0.9002

## Algorithms & MassQL Support

- Added `to_massql()` and a `massql` column in `spectra_to_queries()` output to
  translate diagnostic fragment (`_frag`) and neutral-loss (`_nl`) expressions
  into directly executable MassQL queries (`QUERY scaninfo(MS2DATA) WHERE ...`).
- Improved peak alignment and candidate query discovery:
  - Intensity-weighted consensus centroiding in `harmonize_mzs()` and
    `fix_binned_mzs()` within `dalton + mz * ppm / 1e6`.
  - Enforced physical neutral-loss bounds (`nl > 0` and `precursor_mz - nl > 0`)
    during preprocessing and query evaluation.
  - Replaced top-N single-ion ranking with greedy hitting-set pruning on
    100%-recall ions and forward conjunctive MCC selection penalized by
    negative-class correlation (mRMR).
  - Added Branch-and-Bound pruning of ion combinations using the zero-FP MCC
    upper bound, Boolean superset absorption in `combine_ions_minimal()`, and
    selection of the single most parsimonious query per skeleton.

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
