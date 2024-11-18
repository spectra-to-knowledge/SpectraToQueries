#' @title Spectra to queries
#'
#' @description This function converts spectra to queries.
#'
#' @param spectra Spectra path
#' @param export Export path
#' @param beta Beta parameter of the F-score calculation
#' @param dalton Tolerance in Dalton. Default to 0.01
#' @param decimals Number of decimals for rounding. Default to 4
#' @param intensity_min Minimal intensity. Default to 0
#' @param ions_max Maximal number of ions in the query. Default to 10
#' @param n_skel_min Minimal number of individuals per skeleton. Default to 5
#' @param n_spec_min Minimal number of individuals where a signal has to be found. Default to 3
#' @param ppm Tolerance in parts per million Default to 25
#' @param senspe_min Minimal product of inner and outer ratios. Default to 0.1
#' @param sensitivity_min Minimal sensitivity. Default to 0.3
#' @param specificity_min Minimal specificity. Default to 0
#' @param zero_val Zero value for intensity. Default to 0
#'
#' @return A file with diagnostic query ions
#'
#' @export
#'
#' @examples NULL
spectra_to_queries <- function(spectra = NULL,
                               export = "data/interim/queries.tsv",
                               beta = 0.5,
                               dalton = 0.01,
                               decimals = 4L,
                               intensity_min = 0L,
                               ions_max = 10L,
                               n_skel_min = 5L,
                               n_spec_min = 3L,
                               ppm = 25L,
                               senspe_min = 0.1,
                               sensitivity_min = 0.3,
                               specificity_min = 0L,
                               zero_val = 0L) {
  if (is.null(spectra)) {
    message("No spectra given, loading example spectra.")
    mia_spectra <- readRDS(system.file("extdata", "spectra.rds", package = "SpectraToQueries"))
  } else {
    message("Loading spectra.")
    mia_spectra <- spectra |>
      MsBackendMgf::readMgf() |>
      Spectra::Spectra()
  }

  message(
    "Cut the fragments lower than ",
    intensity_min,
    " and rough preprocessing"
  )
  mia_spectra_1 <- mia_spectra |>
    Spectra::filterMsLevel(2L) |>
    Spectra::reduceSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::combineSpectra(f = mia_spectra$TITLE) |>
    Spectra::deisotopeSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::filterPrecursorPeaks(
      tolerance = dalton,
      ppm = ppm,
      mz = ">="
    ) |>
    Spectra::filterEmptySpectra() |>
    Spectra::addProcessing(normalize_peaks()) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing()

  message(
    "Harmonize m/z values across spectra, given ",
    dalton,
    " Da or ",
    ppm,
    " ppm tolerance."
  )
  mia_spectra@backend@spectraData$precursorMz <-
    mia_spectra@backend@spectraData$PRECURSOR_MZ |>
    as.numeric()
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
    Spectra::combineSpectra(f = mia_spectra$TITLE) |>
    Spectra::deisotopeSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::filterPrecursorPeaks(
      tolerance = dalton,
      ppm = ppm,
      mz = ">="
    ) |>
    Spectra::filterEmptySpectra() |>
    Spectra::addProcessing(normalize_peaks()) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing()
  mia_spectra_w_nl <- mia_spectra_nl |>
    harmonize_mzs(dalton = dalton, ppm = ppm)

  message("Bin spectra to get a matrix.")
  message("The window is ", dalton, " Da or ", ppm, " ppm tolerance.")
  mia_spectra_binned <- mia_spectra_w |>
    Spectra::bin(binSize = dalton, zero.rm = FALSE) |>
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
  message("Fixing mzs")
  spectra_mat <- spectra_mat |>
    fix_binned_mzs(
      original_mzs = mia_spectra_w,
      decimals = decimals,
      dalton = dalton
    )
  rm(mia_spectra_binned)

  mia_spectra_binned_nl <- mia_spectra_w_nl |>
    Spectra::reset() |>
    Spectra::bin(binSize = dalton, zero.rm = FALSE) |>
    Spectra::applyProcessing()
  spectra_nl_mat <- mia_spectra_binned_nl |>
    create_matrix(name = mia_spectra_binned_nl$SKELETON) |>
    filter_matrix(n = n_spec_min)
  message("Fixing mzs")
  spectra_nl_mat <- spectra_nl_mat |>
    fix_binned_mzs(
      original_mzs = mia_spectra_w_nl,
      decimals = decimals,
      dalton = dalton
    )
  rm(mia_spectra_binned_nl)

  message("Create a matrix containing fragments and neutral losses.")
  message("Round the values to ", decimals, ".")
  merged_mat <- cbind(spectra_mat, spectra_nl_mat) |>
    as.matrix()
  colnames(merged_mat) <- NULL
  zeros <- colSums(merged_mat) <= zero_val
  merged_mat <- merged_mat[, !zeros]
  tmp <- merged_mat |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "group")
  rownames(merged_mat) <- tmp$group
  colnames(merged_mat) <-
    c(paste0(round(
      colnames(spectra_mat) |>
        as.numeric(), decimals
    ), "_frag"), paste0(round(
      colnames(spectra_nl_mat) |>
        as.numeric(), decimals
    ), "_nl"))

  message("Count the number of members per skeleton and pivot the matrix.")
  ions_table <- merged_mat |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "group") |>
    tidytable::mutate(group = gsub(
      pattern = "\\.[0-9]{1,3}",
      replacement = "",
      x = group
    )) |>
    tidytable::group_by(group) |>
    tidytable::add_count(name = "group_count") |>
    tidytable::pivot_longer(cols = !starts_with("group"), names_to = "ion") |>
    tidytable::filter(value != 0)

  message("Extract the best ions to perform a query.")
  # Calculate the sensitivity (and specificity) of the features.
  # Filter only features which sensitivity is at least `sensitivity_min`.
  # Filter only features which specificity is at least `specificity_min`.
  # Filter only ions occurring in at least `n_spec_min` spectra.
  # Filter only groups with at least `n_skel_min` spectra.
  ions_table_filtered_1 <- ions_table |>
    tidytable::group_by(ion) |>
    tidytable::add_count(name = "count_per_ion") |>
    tidytable::ungroup() |>
    tidytable::group_by(group, ion) |>
    tidytable::add_count(name = "count_per_ion_per_group") |>
    tidytable::ungroup() |>
    tidytable::select(-value) |>
    tidytable::distinct() |>
    tidytable::filter(count_per_ion >= n_spec_min) |>
    tidytable::filter(group_count >= n_skel_min) |>
    tidytable::mutate(
      ratio_inter = count_per_ion_per_group / count_per_ion,
      ratio_intra = count_per_ion_per_group / group_count
    ) |>
    tidytable::filter(ratio_inter >= specificity_min) |>
    tidytable::filter(ratio_intra >= sensitivity_min) |>
    tidytable::mutate(senspe = ratio_intra * ratio_inter) |>
    tidytable::arrange(tidytable::desc(senspe)) |>
    tidytable::filter(senspe >= senspe_min | ratio_intra == 1) |>
    tidytable::mutate(value = 1)

  ions_table_diagnostic <- ions_table_filtered_1 |>
    tidytable::filter(ratio_intra == 1) |>
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
    tidytable::filter(senspe >= senspe_min & ratio_intra != 1) |>
    tidytable::slice_head(n = ions_max) |>
    tidytable::ungroup() |>
    tidytable::distinct(group, ion, value) |>
    tidytable::pivot_wider(
      names_from = ion,
      values_from = value,
      values_fn = mean
    ) |>
    tibble::column_to_rownames("group")

  # Extract the matching ions per skeleton.
  ions_list_diagnostic <-
    apply(
      X = ions_table_diagnostic[, seq_len(ncol(ions_table_diagnostic))],
      MARGIN = 1,
      FUN = function(x) {
        names(which(x > 0))
      }
    )
  ions_list <-
    apply(
      X = ions_table_final[, seq_len(ncol(ions_table_final))],
      MARGIN = 1,
      FUN = function(x) {
        names(which(x > 0))
      }
    )

  message("Generate all combinations of queries.")
  combinations <- names(ions_list) |>
    generate_combinations_progress(ions_list = ions_list, max_ions = ions_max) |>
    progressr::with_progress(enable = TRUE)
  names(combinations) <- names(ions_list)

  all_combinations <- combinations |>
    unlist(recursive = FALSE)
  names(all_combinations) <- names(all_combinations) |>
    gsub(pattern = "\\d", replacement = "")

  message("Test the queries.")
  queries_results <- seq_along(all_combinations) |>
    perform_list_of_queries_progress(ions_list = all_combinations, spectra = mia_spectra) |>
    progressr::with_progress(enable = TRUE)
  names(queries_results) <- names(all_combinations)

  message("Evaluate the performance of the query based on F-score.")
  results_stats <- seq_along(queries_results) |>
    furrr::future_map(
      .f = function(result, beta) {
        tp <- nrow(queries_results[[result]] |>
          tidytable::filter(target == value))
        fp <- nrow(queries_results[[result]] |>
          tidytable::filter(target != value))
        fn <-
          length(mia_spectra$SKELETON[mia_spectra$SKELETON |> gsub(
            pattern = "+",
            replacement = ".",
            fixed = TRUE
          ) == names(queries_results)[result]]) - tp
        recall <- tp / (tp + fn)
        precision <- tp / (tp + fp)
        f_beta <-
          (1 + beta^2) * (precision * recall) / ((precision * beta^
            2) + recall)
        return(round(f_beta, decimals))
      },
      beta = beta
    )

  best_queries <- data.frame(
    skeleton = names(all_combinations),
    fscore = c(results_stats |> as.character()),
    ions = I(all_combinations)
  ) |>
    tidytable::filter(fscore != NaN) |>
    tidytable::arrange(tidytable::desc(fscore)) |>
    tidytable::group_by(skeleton) |>
    tidytable::mutate(id = match(fscore, unique(fscore))) |>
    tidytable::filter(id == 1) |>
    tidytable::select(-id)

  message("Export the best queries for further use.")
  create_dir(export)
  best_queries |>
    tidytable::arrange(tidytable::desc(fscore)) |>
    tidytable::fwrite(file = export, sep = "\t")
}
