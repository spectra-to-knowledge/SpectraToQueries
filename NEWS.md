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

## Bug Fixes & Compatibility

- Fixed `to_massql()` and `combine_ions_minimal()` Boolean precedence and MassQL
  qualifier handling: parenthesized multi-term `OR` groups, parsed disjunctive
  normal form clauses joined by `" ; "`, and emitted a single effective
  tolerance modifier (`:TOLERANCEMZ=` when `dalton > 0`, or `:TOLERANCEPPM=`
  when `dalton == 0` and `ppm > 0`).
- Enforced relative base-peak intensity filtering `[0, 100]` via `intensity_min`
  during `prepare_query_data()`.
- Replaced internal `@backend@` slot access with public `Spectra` methods so
  in-memory backends (`MsBackendMemory`, `MsBackendDataFrame`) work reliably,
  and added classed errors (`stq_input_error`, `stq_file_not_found_error`) plus
  lowercase `skeleton` metadata fallback.
- Fixed trailing-digit corruption on skeleton labels ending in digits by
  replacing `names(unlist(...))` with `rep(names(...), lengths(...))`, and
  warned if `+` to `.` normalization causes skeleton label collisions.
- Fixed `create_dir()` directory creation, duplicate m/z accumulation in
  `create_matrix()`, and non-finite intensity filtering in `fix_binned_mzs()`.
- Evaluated helper validations before calling S4 `Spectra` generics so classed
  conditions propagate directly without wrapping.

## Performance & Optimization

- Simplified the pipeline to single-pass candidate evaluation in
  `prune_combinations_bnb()` and `summarize_best_queries()` using precomputed
  ion hit vectors (`precompute_ion_hits()`).
- Removed dead helper code (`generate_combinations*()`, `contains_all_mz*()`,
  `perform_query()`, `perform_list_of_queries*()`, `merge_diagnostic_ions()`,
  `show_progress()`, and unused `inst/scripts/spectra_to_queries.R`), keeping
  internal helpers hidden with `@noRd`.
- Replaced nested linear scans in `harmonize_mzs()` with vectorized
  `findInterval()` lookup and avoided redundant `Spectra` processing passes.

## Dependencies

- Narrowed `Imports` to `MsBackendMgf`, `Spectra`, `stats`, and `utils`
  (removing `tidytable`, `fastmatch`, `stringi`, `BiocParallel`, `data.table`,
  `progress`, `purrr`, and `tibble`).

## Data, Documentation & Metadata

- Regenerated `data/mia_spectra_df.rda` (321 MS/MS spectra across 81 skeletons;
  24,255 peaks) and `data/mia_spectra_grouped_df.rda` (321 MS/MS spectra across
  53 grouped classes) with `basename()` applied to `dataOrigin` so no local
  paths are embedded.
- Regenerated `data/interim/queries.tsv` and `data/interim/queries-grouped.tsv`.
- Added executable `@examples` for `spectra_to_queries()` and updated vignette,
  `README`, `CONTRIBUTING.md`, `TESTING.md`, and `CITATION.cff`.
- Added verified EDAM ontology terms (`X-EDAM-Topic`, `X-EDAM-Operation`,
  `X-EDAM-Data`, `X-EDAM-Format`, `X-schema.org-keywords`, `biocViews`) to
  `DESCRIPTION`.

# SpectraToQueries 0.0.9001

- Updated minimal R version to `4.4.0` (and related Bioconductor dependencies)

# SpectraToQueries 0.0.9000

- Initial version.
