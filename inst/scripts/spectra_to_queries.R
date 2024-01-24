start <- Sys.time()

require(
  package = "spectra2queries",
  quietly = TRUE
)

paths <- "inst/paths.yaml" |>
  parse_yaml()
params <- "inst/params.yaml" |>
  parse_yaml()

message("Loading parameters")
if (Sys.info()["sysname"] == "Windows") {
  mc.cores <- 1
} else {
  mc.cores <- parallel::detectCores()
}
BETA <- params$misc$beta
BIN_WINDOWS <- params$misc$bin_windows
DALTON <- params$ms$tolerances$mass$dalton$ms2
DECIMALS <- params$misc$decimals
INTENSITY_MIN <- params$ms$thresholds$ms2$intensity
IONS_MAX <- params$misc$max_ions
N_SPEC_MIN <- params$misc$min_n_spectra
SENSITIVITY_MIN <- params$misc$min_sensitivity
ZERO_VAL <- params$misc$min_total_intensity

message("Import the spectra and ensure they are normalized.")
message("Cut the fragments lower than ", INTENSITY_MIN, " off.")
mia_spectra <- paths$data$source$spectra$mia |>
  MsBackendMgf::readMgf() |>
  Spectra::Spectra() |>
  Spectra::filterMsLevel(2L) |>
  Spectra::addProcessing(normalize_peaks()) |>
  Spectra::filterIntensity(intensity = c(INTENSITY_MIN, Inf)) |>
  Spectra::applyProcessing()

message(
  "Harmonize m/z values across spectra, given ",
  DALTON,
  " tolerance."
)
mia_spectra <- mia_spectra |>
  harmonize_mzs(dalton = DALTON)

message("Calculate neutral losses, given ", DALTON, " tolerance.")
message("Remove the ones above the precursor.")
mia_spectra_nl <- mia_spectra |>
  Spectra::neutralLoss(Spectra::PrecursorMzParam(
    filterPeaks = c("abovePrecursor"),
    msLevel = 2L,
    tolerance = DALTON
  )) |>
  Spectra::applyProcessing()

message("Bin spectra to get a matrix.")
message("The window is ", DALTON, " divided by ", BIN_WINDOWS, ".")
mia_spectra_binned <- mia_spectra |>
  Spectra::bin(binSize = DALTON / BIN_WINDOWS)
mia_spectra_binned_nl <- mia_spectra_nl |>
  Spectra::bin(binSize = DALTON / BIN_WINDOWS)

message("Create fragments and neutral losses matrices.")
message(
  "Remove features not appearing in at least ",
  N_SPEC_MIN,
  " spectra."
)
spectra_mat <- mia_spectra_binned |>
  create_matrix(name = mia_spectra_binned$SKELETON) |>
  filter_matrix(n = N_SPEC_MIN)
spectra_nl_mat <- mia_spectra_binned_nl |>
  create_matrix(name = mia_spectra_binned_nl$SKELETON) |>
  filter_matrix(n = N_SPEC_MIN)

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
  c(
    paste0(round(
      colnames(spectra_mat) |>
        as.numeric(), DECIMALS
    ), "_frag"),
    paste0(round(
      colnames(spectra_nl_mat) |>
        as.numeric(), DECIMALS
    ), "_nl")
  )

message("Count the number of members per skeleton and pivot the matrix.")
ions_table <- merged_mat |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "group") |>
  tidytable::mutate(group = gsub(
    pattern = "\\.",
    replacement = " ",
    x = group
  )) |>
  tidytable::mutate(group = gsub(
    pattern = " .*",
    replacement = "",
    x = group
  )) |>
  tidytable::pivot_longer(cols = !starts_with("group"), names_to = "ion") |>
  tidytable::group_by(group, ion) |>
  tidytable::add_count(name = "count_members_per_group") |>
  tidytable::filter(value != 0)

message("Extract the top ", IONS_MAX, " ions to perform a query.")
# Calculate the sensitivity (and specificity) of the features.
# Filter only features which sensitivity is at least `SENSITIVITY`.
# Filter only ions occurring in at least `N_SPEC_MIN` spectra.
# Filter only groups with at least `N_SPEC_MIN` spectra.
ions_table_filtered <- ions_table |>
  tidytable::group_by(ion) |>
  tidytable::add_count(name = "count_ion") |>
  tidytable::group_by(group, ion) |>
  tidytable::add_count(name = "count_ion_per_group") |>
  tidytable::ungroup() |>
  tidytable::select(-value) |>
  tidytable::distinct() |>
  tidytable::filter(count_ion >= N_SPEC_MIN) |>
  tidytable::mutate(
    ratio_inter = count_ion_per_group / count_ion,
    ratio_intra = count_ion_per_group / count_members_per_group
  ) |>
  tidytable::mutate(ratio = ratio_intra / ratio_inter) |>
  tidytable::filter(ratio_intra >= SENSITIVITY_MIN) |>
  tidytable::arrange(tidytable::desc(ratio_inter)) |>
  tidytable::slice_head(
    n = IONS_MAX,
    weight_by = ratio,
    by = c(group)
  ) |>
  tidytable::filter(count_members_per_group >= N_SPEC_MIN) |>
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
        generate_combinations(x = ions_list[[x]])
      names(combinations) <-
        rep(names(ions_list)[x], length(combinations))

      # Test the query.
      queries_results <- seq_along(1:length(combinations)) |>
        pbmcapply::pbmclapply(
          FUN = perform_list_of_queries,
          ions_list = combinations,
          spectra = mia_spectra,
          mc.cores = mc.cores
        )
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
              length(mia_spectra$SKELETON[mia_spectra$SKELETON == names(queries_results)[result]]) - tp
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
