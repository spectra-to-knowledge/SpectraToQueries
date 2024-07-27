start <- Sys.time()
library(dplyr)
require(package = "spectra2queries", quietly = TRUE)
# weird namespace issue
source("R/fix_binned_mzs.R")
strat <- ifelse(test = .Platform$OS.type == "unix",
  yes = "multicore",
  no = "multisession"
)
future::plan(strategy = strat, workers = future::nbrOfWorkers())
progressr::handlers(progressr::handler_txtprogressbar(
  enable = TRUE,
  char = cli::col_yellow(cli::symbol$star)
))

paths <- "inst/paths.yaml" |>
  parse_yaml()
params <- "inst/params.yaml" |>
  parse_yaml()

message("Loading parameters")
mc.cores <- parallel::detectCores()
BETA <- params$misc$beta
DALTON <- params$ms$tolerances$mass$dalton
DECIMALS <- params$misc$decimals
INTENSITY_MIN <- params$ms$thresholds$ms2$intensity
IONS_MAX <- params$misc$max_ions
N_SKEL_MIN <- params$misc$min_n_skeleton
N_SPEC_MIN <- params$misc$min_n_spectra
PPM <- params$ms$tolerances$mass$ppm
SENSITIVITY_MIN <- params$misc$min_sensitivity
SPECIFICITY_MIN <- params$misc$min_specificity
ZERO_VAL <- params$misc$min_total_intensity

message("Import the spectra and ensure they are normalized.")
message("Cut the fragments lower than ", INTENSITY_MIN, " off.")
mia_spectra <- paths$data$source$spectra$mia |>
  MsBackendMgf::readMgf() |>
  Spectra::Spectra()
mia_spectra_1 <- mia_spectra |>
  Spectra::filterMsLevel(2L) |>
  Spectra::reduceSpectra(tolerance = DALTON, ppm = PPM) |>
  Spectra::combineSpectra(f = mia_spectra$TITLE) |>
  Spectra::deisotopeSpectra(tolerance = DALTON, ppm = PPM) |>
  Spectra::filterPrecursorPeaks(
    tolerance = DALTON,
    ppm = PPM,
    mz = ">="
  ) |>
  Spectra::filterEmptySpectra() |>
  Spectra::addProcessing(normalize_peaks()) |>
  Spectra::filterIntensity(intensity = c(INTENSITY_MIN, Inf)) |>
  Spectra::applyProcessing()

message(
  "Harmonize m/z values across spectra, given ",
  DALTON,
  " Da or ",
  PPM,
  " ppm tolerance."
)
mia_spectra@backend@spectraData$precursorMz <-
  mia_spectra@backend@spectraData$PRECURSOR_MZ |>
  as.numeric()
mia_spectra_w <- mia_spectra |>
  harmonize_mzs(dalton = DALTON, ppm = PPM)

message(
  "Calculate neutral losses, given ",
  DALTON,
  " Da or ",
  PPM,
  " ppm tolerance."
)
message("Remove the ones above the precursor.")
mia_spectra_nl <- mia_spectra_w |>
  Spectra::neutralLoss(Spectra::PrecursorMzParam(
    filterPeaks = c("abovePrecursor"),
    msLevel = 2L,
    tolerance = DALTON,
    ppm = PPM
  )) |>
  Spectra::reduceSpectra(tolerance = DALTON, ppm = PPM) |>
  Spectra::combineSpectra(f = mia_spectra$TITLE) |>
  Spectra::deisotopeSpectra(tolerance = DALTON, ppm = PPM) |>
  Spectra::filterPrecursorPeaks(
    tolerance = DALTON,
    ppm = PPM,
    mz = ">="
  ) |>
  Spectra::filterEmptySpectra() |>
  Spectra::addProcessing(normalize_peaks()) |>
  Spectra::filterIntensity(intensity = c(INTENSITY_MIN, Inf)) |>
  Spectra::applyProcessing()
mia_spectra_w_nl <- mia_spectra_nl |>
  harmonize_mzs(dalton = DALTON, ppm = PPM)

message("Bin spectra to get a matrix.")
message("The window is ", DALTON, " Da or ", PPM, " ppm tolerance.")
mia_spectra_binned <- mia_spectra_w |>
  Spectra::bin(binSize = DALTON, zero.rm = FALSE) |>
  Spectra::applyProcessing()

message("Create fragments and neutral losses matrices.")
message(
  "Remove features not appearing in at least ",
  N_SPEC_MIN,
  " spectra."
)
spectra_mat <- mia_spectra_binned |>
  create_matrix(name = mia_spectra_binned$SKELETON) |>
  filter_matrix(n = N_SPEC_MIN)
message("Fixing mzs")
spectra_mat <- spectra_mat |>
  fix_binned_mzs(
    original_mzs = mia_spectra_w,
    decimals = DECIMALS,
    dalton = DALTON
  )
rm(mia_spectra_binned)

mia_spectra_binned_nl <- mia_spectra_w_nl |>
  Spectra::reset() |>
  Spectra::bin(binSize = DALTON, zero.rm = FALSE) |>
  Spectra::applyProcessing()
spectra_nl_mat <- mia_spectra_binned_nl |>
  create_matrix(name = mia_spectra_binned_nl$SKELETON) |>
  filter_matrix(n = N_SPEC_MIN)
message("Fixing mzs")
spectra_nl_mat <- spectra_nl_mat |>
  fix_binned_mzs(
    original_mzs = mia_spectra_w_nl,
    decimals = DECIMALS,
    dalton = DALTON
  )
rm(mia_spectra_binned_nl)

message("Create a matrix containing fragments and neutral losses.")
message("Round the values to ", DECIMALS, ".")
merged_mat <- cbind(spectra_mat, spectra_nl_mat) |>
  as.matrix()
colnames(merged_mat) <- NULL
zeros <- colSums(merged_mat) <= ZERO_VAL
merged_mat <- merged_mat[, !zeros]
tmp <- merged_mat |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "group")
rownames(merged_mat) <- tmp$group
colnames(merged_mat) <-
  c(paste0(round(
    colnames(spectra_mat) |>
      as.numeric(), DECIMALS
  ), "_frag"), paste0(round(
    colnames(spectra_nl_mat) |>
      as.numeric(), DECIMALS
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
# Filter only features which sensitivity is at least `SENSITIVITY_MIN`.
# Filter only features which specificity is at least `SPECIFICITY_MIN`.
# Filter only ions occurring in at least `N_SPEC_MIN` spectra.
# Filter only groups with at least `N_SKEL_MIN` spectra.
ions_table_filtered <- ions_table |>
  tidytable::group_by(ion) |>
  tidytable::add_count(name = "count_per_ion") |>
  tidytable::ungroup() |>
  tidytable::group_by(group, ion) |>
  tidytable::add_count(name = "count_per_ion_per_group") |>
  tidytable::ungroup() |>
  tidytable::select(-value) |>
  tidytable::distinct() |>
  tidytable::filter(count_per_ion >= N_SPEC_MIN) |>
  tidytable::filter(group_count >= N_SKEL_MIN) |>
  tidytable::mutate(
    ratio_inter = count_per_ion_per_group / count_per_ion,
    ratio_intra = count_per_ion_per_group / group_count
  ) |>
  tidytable::filter(ratio_inter >= SPECIFICITY_MIN) |>
  tidytable::filter(ratio_intra >= SENSITIVITY_MIN) |>
  # tidytable::arrange(tidytable::desc(ratio_intra)) |>
  # tidytable::arrange(tidytable::desc(ratio_inter)) |>
  # tidytable::group_by(group) |>
  # tidytable::slice_head(n = IONS_MAX) |>
  # tidytable::ungroup() |>
  tidytable::mutate(value = 1)

# Pivot back again.
ions_table_final <- ions_table_filtered |>
  tidytable::distinct(group, ion, value) |>
  tidytable::pivot_wider(
    names_from = ion,
    values_from = value,
    values_fn = mean
  ) |>
  tibble::column_to_rownames("group")

# Extract the matching ions per skeleton.
ions_list <-
  apply(
    X = ions_table_final[, 1:ncol(ions_table_final)],
    MARGIN = 1,
    FUN = function(x) {
      names(which(x > 0))
    }
  )

best_queries <- seq_along(ions_list) |>
  lapply(
    FUN = function(x) {
      message("Generate all combinations of queries.")
      combinations <-
        generate_combinations(x = ions_list[[x]], max_ions = IONS_MAX)

      names(combinations) <-
        rep(names(ions_list)[x], length(combinations))

      message("Test the queries.")
      queries_results <- seq_along(1:length(combinations)) |>
        perform_list_of_queries_progress(ions_list = combinations, spectra = mia_spectra) |>
        progressr::with_progress()
      names(queries_results) <-
        rep(names(ions_list)[x], length(combinations))

      message("Evaluate the performance of the query based on F-score.")
      results_stats <- seq_along(1:length(queries_results)) |>
        lapply(
          FUN = function(result) {
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
            tpfn <- tp + fn
            tpfp <- tp + fp
            recall <- tp / tpfn
            precision <- tp / tpfp
            rxp <- recall * precision
            rpp <- recall + precision
            f_beta <-
              (1 + BETA^2) * (precision * recall) / ((precision * BETA^
                2) + recall)
            return(round(f_beta, DECIMALS))
          }
        )

      message("Finished evaluating the query for ", names(ions_list)[x])
      message("Extract the best query (with ties) based on its F-score.")
      return(
        data.frame(
          skeleton = names(combinations),
          fscore = c(results_stats |> as.character()),
          ions = I(combinations)
        ) |>
          tidytable::bind_rows() |>
          tidytable::filter(fscore != NaN) |>
          tidytable::arrange(tidytable::desc(fscore)) |>
          tidytable::mutate(id = match(fscore, unique(fscore))) |>
          tidytable::filter(id == 1) |>
          tidytable::select(-id)
      )
    }
  )

message("Export the best queries for further use.")
create_dir(paths$data$interim$queries)
best_queries |>
  tidytable::bind_rows() |>
  tidytable::arrange(tidytable::desc(fscore)) |>
  tidytable::fwrite(file = paths$data$interim$queries, sep = "\t")

end <- Sys.time()

message("Script finished in ", format(end - start))
