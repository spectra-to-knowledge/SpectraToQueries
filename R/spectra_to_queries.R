#' @title Normalize and intensity-filter spectra
#' @keywords internal
#' @noRd
normalize_and_filter_spectra <- function(sp, intensity_min) {
  sp |>
    Spectra::addProcessing(normalize_peaks()) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing()
}

#' @title Preprocess spectra (reduce, combine by title, deisotope, normalize)
#' @keywords internal
#' @noRd
preprocess_spectra <- function(spectra_obj, dalton, ppm, intensity_min) {
  sp <- Spectra::reduceSpectra(spectra_obj, tolerance = dalton, ppm = ppm)
  if ("TITLE" %in% Spectra::spectraVariables(sp)) {
    titles <- sp$TITLE
    if (anyDuplicated(titles, incomparables = NA) > 0L) {
      # `combineSpectra()` drops spectra whose grouping key is NA
      no_title <- is.na(titles)
      titles[no_title] <- sprintf("\r%d", which(no_title))
      sp <- Spectra::combineSpectra(sp, f = titles)
    }
  }
  sp |>
    Spectra::deisotopeSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::filterPrecursorPeaks(
      tolerance = dalton,
      ppm = ppm,
      mz = ">="
    ) |>
    normalize_and_filter_spectra(intensity_min) |>
    Spectra::filterEmptySpectra()
}

#' @title Build feature matrix from harmonized spectra and fix m/z centroids
#' @keywords internal
#' @noRd
build_feature_matrix <- function(
  spectra_obj,
  dalton,
  ppm,
  decimals,
  n_spec_min
) {
  mat <- create_matrix(spectra_obj, name = spectra_obj$SKELETON)
  refine_feature_mzs(
    mat,
    spectra = spectra_obj,
    decimals = decimals,
    dalton = dalton,
    ppm = ppm
  ) |>
    filter_matrix(n = n_spec_min)
}

#' @title Build feature matrix with suffixed ion column names
#' @keywords internal
#' @noRd
extract_feature_matrix <- function(
  spectra_obj,
  suffix,
  dalton,
  ppm,
  decimals,
  n_spec_min
) {
  mat <- build_feature_matrix(
    spectra_obj,
    dalton = dalton,
    ppm = ppm,
    decimals = decimals,
    n_spec_min = n_spec_min
  )
  colnames(mat) <- sprintf(
    "%s_%s",
    round(as.numeric(colnames(mat)), decimals),
    suffix
  )
  mat
}

#' @title Spectra to queries
#'
#' @description Converts tandem mass spectra into diagnostic fragment (`_frag`)
#'   and neutral-loss (`_nl`) ion queries ranked by Matthews Correlation
#'   Coefficient (MCC), together with ready-to-run MassQL queries.
#'
#' @details
#' For each skeleton class with at least `n_skel_min` spectra, the pipeline:
#' 1. Normalizes MS2 intensities, removes sub-threshold and non-positive
#'    neutral-loss peaks, and aligns m/z values across spectra using
#'    intensity-weighted consensus centroids within `dalton + mz * ppm / 1e6`.
#' 2. Evaluates candidate fragment and neutral-loss centroids against the input
#'    spectra after `intensity_min` filtering only. A neutral loss `nl` matches
#'    a peak within `dalton` (or `nl * ppm / 1e6` when `dalton == 0`) of
#'    `precursor_mz - nl > 0`, the window of the emitted MassQL `MS2NL`
#'    condition.
#' 3. Selects up to `ions_max` candidate ions per skeleton using greedy
#'    hitting-set pruning on 100%-recall ions followed by forward conjunctive
#'    MCC selection with negative-class correlation penalty (mRMR).
#' 4. Prunes the power set of ion combinations via Branch-and-Bound when the
#'    optimistic zero-FP upper bound cannot reach the current best MCC.
#' 5. Returns, per skeleton, the conjunction with the highest MCC (rounded to
#'    `decimals`); ties go to fewer ions, then higher unrounded MCC, higher
#'    precision, and more fragments. It is translated to MassQL via
#'    [to_massql()].
#'
#' MCC is computed over all MS2 spectra with a non-empty `SKELETON`, with the
#' skeleton as the positive class and all other skeletons as negatives:
#' \deqn{
#'   MCC = \frac{TP \cdot TN - FP \cdot FN}
#'   {\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}}
#' }{
#'   MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
#' }
#' and is set to 0 when the denominator is 0 (Matthews 1975,
#' \doi{10.1016/0005-2795(75)90109-9}; Chicco and Jurman 2020,
#' \doi{10.1186/s12864-019-6413-7}).
#'
#' @param spectra A [Spectra::Spectra] object, a portable `data.frame` (such as
#'   [mia_spectra_df] or [mia_spectra_grouped_df]), or a file path to an MGF
#'   file. If `NULL`, loads the example [mia_spectra_df] dataset. If the string
#'   `"grouped"`, loads the grouped example [mia_spectra_grouped_df] dataset.
#' @param export Export file path for the resulting TSV table (or `NULL` to
#'   skip writing to disk). Default is `"data/interim/queries.tsv"`.
#' @param dalton Absolute m/z tolerance in Da, added to the `ppm` term: two m/z
#'   values match when `|a - b| <= dalton + b * ppm / 1e6`. Default is `0.01`.
#' @param decimals Decimal places for reported m/z values and MCC; MCC ties are
#'   resolved at this precision. Default is `4L`.
#' @param intensity_min Relative intensity threshold in percent of the base peak
#'   (0 to 100). Candidate discovery keeps peaks `>= intensity_min`; MCC is
#'   computed with peaks `> min(intensity_min, 99)`, as MassQL applies
#'   `INTENSITYPERCENT`. Default is `0.0` (all peaks).
#' @param ions_max Maximal number of ions in a query combination. Default is
#'   `10L`.
#' @param n_skel_min Minimal number of MS2 spectra per skeleton class
#'   (spectra without peaks count). Default is `5L`.
#' @param n_spec_min Minimal number of spectra (all classes) in which a fragment
#'   or neutral loss must be found to be a candidate. Default is `3L`.
#' @param ppm Relative m/z tolerance in parts per million, added to `dalton`.
#'   Default is `30.0`.
#' @param mcc_min Minimal single-ion Matthews Correlation Coefficient, in
#'   `[-1, 1]`. Default is `0.0`.
#' @param precision_min Minimal single-ion precision (pre-filter before MCC).
#'   Default is `0.0`.
#' @param recall_min Minimal single-ion recall (pre-filter before MCC). Default
#'   is `0.0`.
#' @param zero_val Feature filter: fragment and neutral-loss features whose
#'   relative intensity (0 to 100 per spectrum), summed over all spectra, is
#'   `<= zero_val` are not evaluated. Default is `0.0` (drop all-zero features).
#'
#' @return Invisibly returns a `data.frame` with columns `skeleton`, `mcc`,
#'   `combined_ions`, and `massql`, and writes the TSV table to `export` when
#'   `export` is non-`NULL`.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   msLevel = c(2L, 2L, 2L, 2L),
#'   SKELETON = c("SkelA", "SkelA", "SkelB", "SkelB"),
#'   precursorMz = c(300.0, 300.0, 300.0, 300.0)
#' )
#' df$mz <- list(
#'   c(100.0, 150.0),
#'   c(100.0, 160.0),
#'   c(200.0, 250.0),
#'   c(200.0, 260.0)
#' )
#' df$intensity <- list(c(80, 100), c(70, 100), c(90, 100), c(85, 100))
#' spectra_to_queries(
#'   spectra = df,
#'   export = NULL,
#'   n_skel_min = 2L,
#'   n_spec_min = 2L
#' )
spectra_to_queries <- function(
  spectra = NULL,
  export = "data/interim/queries.tsv",
  dalton = 0.01,
  decimals = 4L,
  intensity_min = 0.0,
  ions_max = 10L,
  n_skel_min = 5L,
  n_spec_min = 3L,
  ppm = 30.0,
  mcc_min = 0.0,
  precision_min = 0.0,
  recall_min = 0.0,
  zero_val = 0.0
) {
  validate_stq_inputs(
    spectra = spectra,
    export = export,
    dalton = dalton,
    decimals = decimals,
    intensity_min = intensity_min,
    ions_max = ions_max,
    n_skel_min = n_skel_min,
    n_spec_min = n_spec_min,
    ppm = ppm,
    mcc_min = mcc_min,
    precision_min = precision_min,
    recall_min = recall_min,
    zero_val = zero_val
  )

  mia_spectra <- load_input_spectra(spectra)
  mia_spectra <- Spectra::filterMsLevel(mia_spectra, 2L)
  skel_raw <- trimws(as.character(mia_spectra$SKELETON))
  valid_skel <- !is.na(skel_raw) & nzchar(skel_raw)
  if (!all(valid_skel)) {
    mia_spectra <- mia_spectra[valid_skel]
    skel_raw <- skel_raw[valid_skel]
  }

  total_spectra <- length(mia_spectra)
  if (total_spectra == 0L) {
    return(write_query_table(empty_query_table(), export))
  }
  mia_spectra$SKELETON <- skel_raw

  skel_norm <- normalize_skeleton_labels(skel_raw)
  if (length(unique(skel_norm)) < length(unique(skel_raw))) {
    warning(
      "Distinct `SKELETON` labels collide after '+' to '.' normalization.",
      call. = FALSE
    )
  }

  message(
    "Preprocessing ",
    total_spectra,
    " MS2 spectra (relative intensity >= ",
    intensity_min,
    "%)."
  )
  mia_spectra_1 <- preprocess_spectra(
    mia_spectra,
    dalton = dalton,
    ppm = ppm,
    intensity_min = intensity_min
  )

  message(
    "Harmonizing fragment and neutral-loss m/z values (tolerance ",
    dalton,
    " Da + ",
    ppm,
    " ppm)."
  )
  mia_spectra_w <- harmonize_mzs(mia_spectra_1, dalton = dalton, ppm = ppm)

  mia_spectra_nl <- mia_spectra_1 |>
    Spectra::neutralLoss(
      Spectra::PrecursorMzParam(
        filterPeaks = c("abovePrecursor"),
        msLevel = 2L,
        tolerance = dalton,
        ppm = ppm
      )
    ) |>
    Spectra::reduceSpectra(tolerance = dalton, ppm = ppm) |>
    normalize_and_filter_spectra(intensity_min)

  mia_spectra_w_nl <- harmonize_mzs(mia_spectra_nl, dalton = dalton, ppm = ppm)

  spectra_mat <- extract_feature_matrix(
    mia_spectra_w,
    suffix = "frag",
    dalton = dalton,
    ppm = ppm,
    decimals = decimals,
    n_spec_min = n_spec_min
  )
  spectra_nl_mat <- extract_feature_matrix(
    mia_spectra_w_nl,
    suffix = "nl",
    dalton = dalton,
    ppm = ppm,
    decimals = decimals,
    n_spec_min = n_spec_min
  )

  merged_mat <- cbind(spectra_mat, spectra_nl_mat)
  rm(spectra_mat, spectra_nl_mat)

  zeros <- colSums(merged_mat) <= zero_val
  if (all(zeros)) {
    return(write_query_table(empty_query_table(), export))
  }

  cand_cols <- unique(colnames(merged_mat)[!zeros])
  qd <- prepare_query_data(mia_spectra, intensity_min = intensity_min)
  hits <- precompute_ion_hits(
    unique_ions = cand_cols,
    all_mz = qd$mz,
    precursor_mz = qd$precursor_mz,
    valid_prec_and_spec = is.finite(qd$precursor_mz) & (lengths(qd$mz) > 0L),
    dalton = dalton,
    ppm = ppm
  )
  merged_mat <- matrix(
    vapply(
      cand_cols,
      function(u) hits[[u]],
      logical(length(mia_spectra))
    ),
    nrow = length(mia_spectra),
    dimnames = list(mia_spectra$SKELETON, cand_cols)
  )

  message("Selecting candidate ions per skeleton.")
  ions_list <- extract_candidate_ions(
    merged_mat = merged_mat,
    total_spectra = total_spectra,
    n_skel_min = n_skel_min,
    n_spec_min = n_spec_min,
    ions_max = ions_max,
    precision_min = precision_min,
    recall_min = recall_min,
    mcc_min = mcc_min
  )
  rm(merged_mat)

  if (length(ions_list) == 0L) {
    return(write_query_table(empty_query_table(), export))
  }

  message("Searching ion combinations.")
  candidates <- prune_combinations_bnb(
    ions_list = ions_list,
    skeleton = qd$skeleton,
    hits_lookup = hits,
    max_ions = ions_max,
    decimals = decimals
  )
  rm(ions_list, hits)

  skeleton_counts <- table(qd$skeleton)

  final_queries <- summarize_best_queries(
    candidates = candidates,
    skeleton_counts = skeleton_counts,
    total_spectra = total_spectra,
    decimals = decimals,
    dalton = dalton,
    ppm = ppm,
    intensity_min = intensity_min
  )

  write_query_table(final_queries, export)
}
