#' @title Spectra to queries
#'
#' @description This function converts spectra to queries.
#'
#' @param spectra Spectra path
#' @param export Export path
#' @param dalton Tolerance in Dalton. Default to 0.01
#' @param decimals Number of decimals for rounding. Default to 4L
#' @param intensity_min Minimal intensity. Default to 0.0
#' @param ions_max Maximal number of ions in the query. Default to 10L
#' @param n_skel_min Minimal number of individuals per skeleton. Default to 5L
#' @param n_spec_min Minimal number of individuals where a signal has to be found. Default to 3L
#' @param ppm Tolerance in parts per million Default to 30.0
#' @param mcc_min Minimal single-ion Matthews Correlation Coefficient. Default to 0.0
#' @param precision_min Minimal single ion precision (pre-filter before MCC). Default to 0.0
#' @param recall_min Minimal single ion recall (pre-filter before MCC). Default to 0.0
#' @param zero_val Zero value for intensity. Default to 0.0
#'
#' @return A file with diagnostic query ions
#'
#' @export
#'
#' @examples NULL
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
  if (is.null(spectra)) {
    message("No spectra given, loading example spectra.")
    mia_spectra <- readRDS(system.file(
      "extdata",
      "spectra.rds",
      package = "SpectraToQueries"
    ))
  } else if (
    spectra ==
      system.file(
        "extdata",
        "spectra_grouped.rds",
        package = "SpectraToQueries"
      )
  ) {
    message("Loading grouped spectra.")
    mia_spectra <- readRDS(spectra)
  } else {
    message("Loading spectra.")
    mia_spectra <- spectra |>
      MsBackendMgf::readMgf() |>
      Spectra::Spectra(BPPARAM = BiocParallel::SerialParam()) |>
      Spectra::setBackend(
        backend = Spectra::MsBackendMemory(),
        BPPARAM = BiocParallel::SerialParam()
      )
  }
  mia_spectra@backend@spectraData$precursorMz <-
    mia_spectra@backend@spectraData$PRECURSOR_MZ |>
    as.numeric()

  # Total number of spectra — required as the denominator for TN in MCC:
  #   TN = total_spectra - TP - FP - FN
  total_spectra <- length(mia_spectra)

  message(
    "Cut the fragments lower than ",
    intensity_min,
    " and rough preprocessing"
  )
  mia_spectra_1 <- mia_spectra |>
    Spectra::filterMsLevel(2L) |>
    Spectra::reduceSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::combineSpectra(
      f = mia_spectra$TITLE,
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    Spectra::deisotopeSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::filterPrecursorPeaks(
      tolerance = dalton,
      ppm = ppm,
      mz = ">="
    ) |>
    Spectra::filterEmptySpectra() |>
    Spectra::addProcessing(
      normalize_peaks(),
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing()

  message(
    "Harmonize m/z values across spectra, given ",
    dalton,
    " Da or ",
    ppm,
    " ppm tolerance."
  )
  mia_spectra_w <- mia_spectra_1 |>
    harmonize_mzs(dalton = dalton, ppm = ppm)

  message(
    "Calculate neutral losses, given ",
    dalton,
    " Da or ",
    ppm,
    " ppm tolerance."
  )
  message("Remove the ones above the precursor.")
  mia_spectra_nl <- mia_spectra_w |>
    Spectra::neutralLoss(
      Spectra::PrecursorMzParam(
        filterPeaks = c("abovePrecursor"),
        msLevel = 2L,
        tolerance = dalton,
        ppm = ppm
      )
    ) |>
    Spectra::reduceSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::combineSpectra(
      f = mia_spectra$TITLE,
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    Spectra::deisotopeSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::filterPrecursorPeaks(
      tolerance = dalton,
      ppm = ppm,
      mz = ">="
    ) |>
    Spectra::filterEmptySpectra() |>
    Spectra::addProcessing(
      normalize_peaks(),
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing()
  mia_spectra_w_nl <- mia_spectra_nl |>
    harmonize_mzs(dalton = dalton, ppm = ppm)

  message("Bin spectra to get a matrix.")
  message("The window is ", dalton, " Dalton")
  mia_spectra_binned <- mia_spectra_w |>
    Spectra::bin(
      binSize = dalton,
      zero.rm = FALSE
    ) |>
    Spectra::applyProcessing()

  message("Create fragments and neutral losses matrices.")
  message(
    "Remove features not appearing in at least ",
    n_spec_min,
    " spectra."
  )
  spectra_mat <- mia_spectra_binned |>
    create_matrix(name = mia_spectra_binned$SKELETON) |>
    filter_matrix(n = n_spec_min)
  message("Fixing mzs (fragments)")
  spectra_mat <- spectra_mat |>
    fix_binned_mzs(
      original_mzs = mia_spectra_w,
      decimals = decimals,
      dalton = dalton,
      ppm = ppm
    )
  rm(mia_spectra_binned)
  gc(verbose = FALSE)

  mia_spectra_binned_nl <- mia_spectra_w_nl |>
    Spectra::reset() |>
    Spectra::bin(
      binSize = dalton,
      zero.rm = FALSE
    ) |>
    Spectra::applyProcessing()
  spectra_nl_mat <- mia_spectra_binned_nl |>
    create_matrix(name = mia_spectra_binned_nl$SKELETON) |>
    filter_matrix(n = n_spec_min)
  message("Fixing mzs (losses)")
  spectra_nl_mat <- spectra_nl_mat |>
    fix_binned_mzs(
      original_mzs = mia_spectra_w_nl,
      decimals = decimals,
      dalton = dalton,
      ppm = ppm
    )
  rm(mia_spectra_binned_nl)
  gc(verbose = FALSE)

  message("Create a matrix containing fragments and neutral losses.")
  message("Round the values to ", decimals, ".")

  # Pre-calculate column names to avoid redundant operations
  frag_names <- paste0(
    round(as.numeric(colnames(spectra_mat)), decimals),
    "_frag"
  )
  nl_names <- paste0(
    round(as.numeric(colnames(spectra_nl_mat)), decimals),
    "_nl"
  )

  # Merge matrices
  merged_mat <- cbind(spectra_mat, spectra_nl_mat)

  # Filter zero columns
  zeros <- colSums(merged_mat) <= zero_val
  if (!all(zeros)) {
    merged_mat <- merged_mat[, !zeros, drop = FALSE]
    # Update column names for non-zero columns
    all_names <- c(frag_names, nl_names)
    colnames(merged_mat) <- all_names[!zeros]
  } else {
    # Handle edge case where all columns are zeros
    merged_mat <- matrix(0, nrow = nrow(merged_mat), ncol = 0)
    colnames(merged_mat) <- character(0)
  }

  message("Count the number of members per skeleton and pivot the matrix.")
  ions_table <- merged_mat |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "group") |>
    tidytable::mutate(
      group = gsub(
        pattern = "\\.[0-9]{1,3}",
        replacement = "",
        x = group
      )
    ) |>
    tidytable::group_by(group) |>
    tidytable::add_count(name = "group_count") |>
    tidytable::ungroup() |>
    tidytable::filter(group_count >= n_skel_min) |>
    tidytable::pivot_longer(
      cols = !tidytable::starts_with("group"),
      names_to = "ion"
    ) |>
    tidytable::filter(value != 0)

  # Clean up large matrices to free memory
  rm(merged_mat, spectra_mat, spectra_nl_mat)
  gc(verbose = FALSE)

  message("Extract the best ions to perform a query.")
  ions_table_filtered_1 <- ions_table |>
    tidytable::group_by(ion) |>
    tidytable::add_count(name = "count_per_ion") |>
    tidytable::filter(count_per_ion >= n_spec_min) |>
    tidytable::ungroup() |>
    tidytable::group_by(group, ion) |>
    tidytable::add_count(name = "count_per_ion_per_group") |>
    tidytable::ungroup() |>
    tidytable::select(-value) |>
    tidytable::distinct() |>
    # Optional pre-filters on precision and recall before MCC computation.
    # These reduce the candidate set cheaply and remain meaningful on their own.
    tidytable::mutate(
      precision = count_per_ion_per_group / count_per_ion,
      recall = count_per_ion_per_group / group_count
    ) |>
    tidytable::filter(precision >= precision_min) |>
    tidytable::filter(recall >= recall_min) |>
    # -----------------------------------------------------------------------
    # Matthews Correlation Coefficient (MCC) for a single ion / group pair.
    #
    #   TP = spectra in this group that contain the ion
    #   FP = spectra in other groups that contain the ion
    #   FN = spectra in this group that do NOT contain the ion
    #   TN = all remaining spectra (neither in this group nor containing the ion)
    #
    #   MCC = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    #
    # When the denominator is zero (degenerate partition), MCC is set to 0.
    # -----------------------------------------------------------------------
    tidytable::mutate(
      tp = count_per_ion_per_group,
      fp = count_per_ion - count_per_ion_per_group,
      fn_ = group_count - count_per_ion_per_group, # fn_ avoids clash with base fn()
      tn = total_spectra - tp - fp - fn_,
      mcc = tidytable::if_else(
        condition = (tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_) > 0,
        true = (tp * tn - fp * fn_) /
          sqrt((tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_)),
        false = 0.0 # degenerate case: ion present in every or no spectrum
      )
    ) |>
    tidytable::select(-tp, -fp, -fn_, -tn) |> # drop intermediate columns
    tidytable::arrange(tidytable::desc(mcc)) |>
    # Keep ions that pass the MCC threshold OR are perfectly diagnostic (recall == 1)
    tidytable::filter(mcc >= mcc_min | recall == 1) |>
    tidytable::mutate(value = 1)

  # Clean up ions_table to free memory
  rm(ions_table)
  gc(verbose = FALSE)

  ions_table_diagnostic <- ions_table_filtered_1 |>
    tidytable::filter(recall == 1) |>
    tidytable::distinct(group, ion, value) |>
    tidytable::pivot_wider(
      names_from = ion,
      values_from = value,
      values_fn = mean
    ) |>
    tibble::column_to_rownames("group")

  # Pivot back again.
  ions_table_final <- ions_table_filtered_1 |>
    tidytable::group_by(group) |>
    # Retain only non-perfectly-diagnostic ions that clear the MCC threshold,
    # keeping the top `ions_max` per skeleton ranked by MCC.
    tidytable::filter(mcc >= mcc_min & recall != 1) |>
    tidytable::arrange(tidytable::desc(mcc)) |>
    tidytable::slice_head(n = ions_max) |>
    tidytable::ungroup() |>
    tidytable::distinct(group, ion, value) |>
    tidytable::pivot_wider(
      names_from = ion,
      values_from = value,
      values_fn = mean
    ) |>
    tibble::column_to_rownames("group")

  # Extract the matching ions per skeleton
  ions_list_diagnostic <- split(
    colnames(ions_table_diagnostic)[col(ions_table_diagnostic)[
      ions_table_diagnostic > 0
    ]],
    row(ions_table_diagnostic)[ions_table_diagnostic > 0]
  )
  names(ions_list_diagnostic) <- rownames(ions_table_diagnostic)

  ions_list <- split(
    colnames(ions_table_final)[col(ions_table_final)[ions_table_final > 0]],
    row(ions_table_final)[ions_table_final > 0]
  )
  names(ions_list) <- rownames(ions_table_final)

  message("Generate all combinations of queries.")
  combinations <- names(ions_list) |>
    generate_combinations_progress(ions_list = ions_list, max_ions = ions_max)
  names(combinations) <- names(ions_list)

  # Combination merging
  message("Merging diagnostic ions with combinations.")

  # Pre-process diagnostic ions
  has_diagnostic <- names(ions_list) %in% names(ions_list_diagnostic)

  if (any(has_diagnostic)) {
    # Only process items that have diagnostic ions
    for (name in names(combinations)[has_diagnostic]) {
      x <- combinations[[name]]
      diagnostic_ions <- ions_list_diagnostic[[name]]

      if (length(diagnostic_ions) > 0) {
        # In-place update more efficient than creating new list
        combinations[[name]] <- lapply(x, function(sublist) {
          unique(c(unlist(sublist), diagnostic_ions))
        })
      }
    }
  }

  # Free memory from intermediate objects
  rm(ions_list, ions_list_diagnostic)
  gc(verbose = FALSE)

  all_combinations <- unlist(combinations, recursive = FALSE)
  names(all_combinations) <- gsub("\\d+$", "", names(all_combinations))

  # Clean up
  rm(combinations)
  gc(verbose = FALSE)

  message("Test the queries. (This is the longest step)")
  queries_results <- all_combinations |>
    perform_list_of_queries_progress(
      spectra = mia_spectra,
      dalton = dalton,
      ppm = ppm
    )
  names(queries_results) <- names(all_combinations)

  message("Evaluate the performance of each query using MCC.")
  results_stats <- queries_results |>
    seq_along() |>
    purrr::map(
      .f = function(result) {
        # -----------------------------------------------------------------------
        # Per-query confusion matrix entries:
        #   TP = query hits that belong to the target skeleton
        #   FP = query hits that belong to a different skeleton
        #   FN = members of the target skeleton missed by the query
        #   TN = all other spectra (not the target skeleton, not hit by query)
        # -----------------------------------------------------------------------
        tp <- nrow(
          queries_results[[result]] |>
            tidytable::filter(target == value)
        )
        fp <- nrow(
          queries_results[[result]] |>
            tidytable::filter(target != value)
        )
        fn_ <-
          length(mia_spectra$SKELETON[
            mia_spectra$SKELETON |>
              gsub(
                pattern = "+",
                replacement = ".",
                fixed = TRUE
              ) ==
              names(queries_results)[result]
          ]) -
          tp
        tn <- total_spectra - tp - fp - fn_

        # MCC — returns 0 for degenerate partitions (zero denominator)
        denom <- sqrt((tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_))
        mcc <- if (denom > 0) (tp * tn - fp * fn_) / denom else 0.0

        return(round(mcc, decimals))
      }
    )

  best_queries <- data.frame(
    skeleton = names(all_combinations),
    mcc = as.numeric(unlist(results_stats)),
    ions = I(all_combinations)
  ) |>
    tidytable::filter(!is.nan(mcc)) |>
    tidytable::arrange(tidytable::desc(mcc)) |>
    tidytable::group_by(skeleton) |>
    # For each skeleton keep only the query with the highest MCC
    tidytable::mutate(id = match(mcc, unique(mcc))) |>
    tidytable::filter(id == 1) |>
    tidytable::select(-id)

  final_queries <- best_queries |>
    tidytable::group_by(skeleton, mcc) |>
    tidytable::summarize(
      combined_ions = combine_ions_minimal(ions),
      .groups = "drop"
    )

  message("Export the best queries for further use.")
  create_dir(export)
  final_queries |>
    tidytable::arrange(tidytable::desc(mcc)) |>
    tidytable::fwrite(file = export, sep = "\t")
}
