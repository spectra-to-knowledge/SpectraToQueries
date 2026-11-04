#' @title Validate inputs for `spectra_to_queries`
#' @keywords internal
#' @noRd
validate_stq_inputs <- function(
  spectra,
  export,
  dalton,
  decimals,
  intensity_min,
  ions_max,
  n_skel_min,
  n_spec_min,
  ppm,
  mcc_min,
  precision_min,
  recall_min,
  zero_val
) {
  abort_input <- function(msg) {
    stop(errorCondition(msg, class = "stq_input_error"))
  }
  is_scalar_num <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x)
  }
  is_valid_str <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }

  valid_spectra <- is.null(spectra) ||
    is.data.frame(spectra) ||
    inherits(spectra, "Spectra") ||
    is_valid_str(spectra)
  if (!valid_spectra) {
    abort_input(
      paste(
        "`spectra` must be NULL, a non-empty string,",
        "a data.frame, or a Spectra object."
      )
    )
  }
  if (!is.null(export) && !is_valid_str(export)) {
    abort_input("`export` must be NULL or a non-empty character string.")
  }
  if (!is_scalar_num(dalton) || dalton < 0) {
    abort_input("`dalton` must be a single finite non-negative number.")
  }
  if (!is_scalar_num(ppm) || ppm < 0) {
    abort_input("`ppm` must be a single finite non-negative number.")
  }
  if (!is_scalar_num(decimals) || decimals < 0) {
    abort_input("`decimals` must be a single finite non-negative number.")
  }
  if (!is_scalar_num(ions_max) || ions_max < 1) {
    abort_input("`ions_max` must be a single finite number >= 1.")
  }
  if (!is_scalar_num(n_skel_min) || n_skel_min < 1) {
    abort_input("`n_skel_min` must be a single finite number >= 1.")
  }
  if (!is_scalar_num(n_spec_min) || n_spec_min < 1) {
    abort_input("`n_spec_min` must be a single finite number >= 1.")
  }
  if (!is_scalar_num(precision_min) || precision_min < 0 || precision_min > 1) {
    abort_input("`precision_min` must be a single finite number in [0, 1].")
  }
  if (!is_scalar_num(recall_min) || recall_min < 0 || recall_min > 1) {
    abort_input("`recall_min` must be a single finite number in [0, 1].")
  }
  if (!is_scalar_num(mcc_min) || mcc_min < -1 || mcc_min > 1) {
    abort_input("`mcc_min` must be a single finite number in [-1, 1].")
  }
  if (!is_scalar_num(intensity_min) || !is_scalar_num(zero_val)) {
    abort_input("`intensity_min` and `zero_val` must be single finite numbers.")
  }
  invisible(TRUE)
}

#' @title Spectra to queries
#'
#' @description This function converts tandem mass spectra into diagnostic
#'   fragment and neutral-loss ion queries ranked by Matthews Correlation
#'   Coefficient (MCC).
#'
#' @param spectra A [Spectra::Spectra] object, a portable `data.frame` (such as
#'   [mia_spectra_df] or [mia_spectra_grouped_df]), or a file path to an MGF
#'   file. If `NULL`, loads the example [mia_spectra_df] dataset. If the string
#'   `"grouped"`, loads the grouped example [mia_spectra_grouped_df] dataset.
#' @param export Export file path for the resulting TSV table (or `NULL` to
#'   skip writing to disk). Default is `"data/interim/queries.tsv"`.
#' @param dalton Tolerance in Dalton. Default is `0.01`.
#' @param decimals Number of decimals for rounding. Default is `4L`.
#' @param intensity_min Minimal relative intensity threshold (0 to 100). Default
#'   is `0.0`.
#' @param ions_max Maximal number of ions in a query combination. Default is
#'   `10L`.
#' @param n_skel_min Minimal number of spectra per skeleton class. Default is
#'   `5L`.
#' @param n_spec_min Minimal number of spectra where a signal has to be found.
#'   Default is `3L`.
#' @param ppm Tolerance in parts per million. Default is `30.0`.
#' @param mcc_min Minimal single-ion Matthews Correlation Coefficient. Default
#'   is `0.0`.
#' @param precision_min Minimal single-ion precision (pre-filter before MCC).
#'   Default is `0.0`.
#' @param recall_min Minimal single-ion recall (pre-filter before MCC). Default
#'   is `0.0`.
#' @param zero_val Zero value threshold for intensity. Default is `0.0`.
#'
#' @return Invisibly returns a `data.frame` with columns `skeleton`, `mcc`, and
#'   `combined_ions`, and writes the TSV table to `export` when `export` is
#'   non-`NULL`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Run with default example spectra (mia_spectra_df)
#'   spectra_to_queries(export = tempfile(fileext = ".tsv"))
#'
#'   # Run with grouped example spectra (mia_spectra_grouped_df)
#'   spectra_to_queries(
#'     spectra = "grouped",
#'     export = tempfile(fileext = ".tsv")
#'   )
#' }
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

  if (is.null(spectra)) {
    message("No spectra given, loading example spectra.")
    data_env <- new.env(parent = emptyenv())
    utils::data(
      "mia_spectra_df",
      package = "SpectraToQueries",
      envir = data_env
    )
    mia_spectra <- df_to_spectra(data_env$mia_spectra_df)
  } else if (
    is.character(spectra) && length(spectra) == 1L && spectra == "grouped"
  ) {
    message("Loading grouped example spectra.")
    data_env <- new.env(parent = emptyenv())
    utils::data(
      "mia_spectra_grouped_df",
      package = "SpectraToQueries",
      envir = data_env
    )
    mia_spectra <- df_to_spectra(data_env$mia_spectra_grouped_df)
  } else if (is.character(spectra)) {
    message("Loading spectra from file.")
    mia_spectra <- spectra |>
      MsBackendMgf::readMgf() |>
      Spectra::Spectra()
  } else if (is.data.frame(spectra)) {
    mia_spectra <- df_to_spectra(spectra)
  } else {
    mia_spectra <- spectra
  }

  if (!is.null(mia_spectra@backend@spectraData$PRECURSOR_MZ)) {
    mia_spectra@backend@spectraData$precursorMz <-
      as.numeric(mia_spectra@backend@spectraData$PRECURSOR_MZ)
  }

  # Total number of spectra (as numeric to prevent integer overflow in MCC
  # when computing true negatives from total_spectra minus TP, FP, and FN)
  total_spectra <- as.numeric(length(mia_spectra))

  combine_by_title <- function(x) {
    titles <- x$TITLE
    if (!is.null(titles) && (anyDuplicated(titles) || is.unsorted(titles))) {
      Spectra::combineSpectra(x, f = titles)
    } else {
      x
    }
  }

  message(
    "Cut the fragments lower than ",
    intensity_min,
    " and rough preprocessing"
  )
  mia_spectra_1 <- mia_spectra |>
    Spectra::filterMsLevel(2L) |>
    Spectra::reduceSpectra(tolerance = dalton, ppm = ppm) |>
    combine_by_title() |>
    Spectra::deisotopeSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::filterPrecursorPeaks(
      tolerance = dalton,
      ppm = ppm,
      mz = ">="
    ) |>
    Spectra::addProcessing(normalize_peaks()) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing() |>
    Spectra::filterEmptySpectra()

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
    combine_by_title() |>
    Spectra::deisotopeSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::filterPrecursorPeaks(
      tolerance = dalton,
      ppm = ppm,
      mz = ">="
    ) |>
    Spectra::addProcessing(normalize_peaks()) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing() |>
    Spectra::filterEmptySpectra()

  mia_spectra_w_nl <- mia_spectra_nl |>
    harmonize_mzs(dalton = dalton, ppm = ppm)

  message("Bin spectra to get a matrix.")
  message("The window is ", dalton, " Dalton")
  message("Create fragments and neutral losses matrices.")
  message(
    "Remove features not appearing in at least ",
    n_spec_min,
    " spectra."
  )
  if (sum(lengths(mia_spectra_w)) > 0L) {
    mia_spectra_binned <- mia_spectra_w |>
      Spectra::bin(
        binSize = dalton,
        zero.rm = TRUE
      ) |>
      Spectra::applyProcessing()

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
  } else {
    message("Fixing mzs (fragments)")
    spectra_mat <- matrix(
      0,
      nrow = length(mia_spectra_w),
      ncol = 0L,
      dimnames = list(mia_spectra_w$SKELETON, character(0))
    )
  }

  if (sum(lengths(mia_spectra_w_nl)) > 0L) {
    mia_spectra_binned_nl <- mia_spectra_w_nl |>
      Spectra::reset() |>
      Spectra::bin(
        binSize = dalton,
        zero.rm = TRUE
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
  } else {
    message("Fixing mzs (losses)")
    spectra_nl_mat <- matrix(
      0,
      nrow = length(mia_spectra_w_nl),
      ncol = 0L,
      dimnames = list(mia_spectra_w_nl$SKELETON, character(0))
    )
  }

  message("Create a matrix containing fragments and neutral losses.")
  message("Round the values to ", decimals, ".")

  frag_names <- if (ncol(spectra_mat) > 0L) {
    paste0(round(as.numeric(colnames(spectra_mat)), decimals), "_frag")
  } else {
    character(0)
  }
  nl_names <- if (ncol(spectra_nl_mat) > 0L) {
    paste0(round(as.numeric(colnames(spectra_nl_mat)), decimals), "_nl")
  } else {
    character(0)
  }

  merged_mat <- cbind(spectra_mat, spectra_nl_mat)
  rm(spectra_mat, spectra_nl_mat)

  zeros <- colSums(merged_mat) <= zero_val
  if (!all(zeros)) {
    merged_mat <- merged_mat[, !zeros, drop = FALSE]
    all_names <- c(frag_names, nl_names)
    colnames(merged_mat) <- all_names[!zeros]
  } else {
    merged_mat <- matrix(0, nrow = nrow(merged_mat), ncol = 0L)
    colnames(merged_mat) <- character(0)
  }

  message("Count the number of members per skeleton and pivot the matrix.")
  groups_vec <- gsub(
    pattern = "+",
    replacement = ".",
    x = as.character(rownames(merged_mat)),
    fixed = TRUE
  )
  group_counts_tab <- table(groups_vec)
  group_count_vec <- as.integer(group_counts_tab[groups_vec])

  keep_rows <- group_count_vec >= n_skel_min
  merged_sub <- merged_mat[keep_rows, , drop = FALSE]
  groups_sub <- groups_vec[keep_rows]
  group_count_sub <- group_count_vec[keep_rows]
  rm(merged_mat)

  # Pivot non-zero entries to long format (column-major order)
  nz_mat <- merged_sub != 0
  col_counts <- colSums(nz_mat)
  keep_ions <- col_counts >= n_spec_min
  nz_mat <- nz_mat[, keep_ions, drop = FALSE]
  col_counts <- col_counts[keep_ions]
  rm(merged_sub)

  nz_idx <- which(nz_mat, arr.ind = TRUE)
  row_idx <- nz_idx[, 1L]
  col_idx <- nz_idx[, 2L]

  ions_df <- data.frame(
    group = groups_sub[row_idx],
    group_count = group_count_sub[row_idx],
    ion = as.character(colnames(nz_mat))[col_idx],
    count_per_ion = as.integer(col_counts[col_idx]),
    stringsAsFactors = FALSE
  )
  rm(nz_mat, nz_idx)

  message("Extract the best ions to perform a query.")
  grp_ion_key <- paste(ions_df$group, ions_df$ion, sep = "\r")
  u_grp_ion <- unique(grp_ion_key)
  grp_ion_tab <- tabulate(
    match(grp_ion_key, u_grp_ion),
    nbins = length(u_grp_ion)
  )
  is_first <- !duplicated(grp_ion_key)
  ions_df <- ions_df[is_first, , drop = FALSE]
  ions_df$count_per_ion_per_group <- grp_ion_tab

  ions_df$precision <- ions_df$count_per_ion_per_group / ions_df$count_per_ion
  ions_df$recall <- ions_df$count_per_ion_per_group / ions_df$group_count

  ions_df <- ions_df[
    ions_df$precision >= precision_min & ions_df$recall >= recall_min,
    ,
    drop = FALSE
  ]

  # -------------------------------------------------------------------------
  # Matthews Correlation Coefficient (MCC) for a single ion and group pair
  # -------------------------------------------------------------------------
  tp <- as.numeric(ions_df$count_per_ion_per_group)
  fp <- as.numeric(ions_df$count_per_ion) - tp
  fn_ <- as.numeric(ions_df$group_count) - tp
  tn <- total_spectra - tp - fp - fn_
  denom_sq <- (tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_)
  ions_df$mcc <- ifelse(
    denom_sq > 0,
    (tp * tn - fp * fn_) / sqrt(denom_sq),
    0.0
  )

  ions_df <- ions_df[order(ions_df$mcc, decreasing = TRUE), , drop = FALSE]
  ions_df <- ions_df[
    ions_df$mcc >= mcc_min | ions_df$recall == 1,
    ,
    drop = FALSE
  ]

  extract_ions_by_group <- function(df_sub) {
    if (nrow(df_sub) == 0L) {
      return(list())
    }
    groups <- sort(unique(df_sub$group), method = "radix")
    u_ions <- unique(df_sub$ion)
    df_sub <- df_sub[order(match(df_sub$ion, u_ions)), , drop = FALSE]
    res <- lapply(groups, function(g) {
      unique(df_sub$ion[df_sub$group == g])
    })
    names(res) <- groups
    res
  }

  diag_df <- ions_df[ions_df$recall == 1, c("group", "ion"), drop = FALSE]
  ions_list_diagnostic <- extract_ions_by_group(diag_df)

  non_diag_df <- ions_df[
    ions_df$mcc >= mcc_min & ions_df$recall != 1,
    ,
    drop = FALSE
  ]
  non_diag_df <- non_diag_df[
    order(non_diag_df$mcc, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  if (nrow(non_diag_df) > 0L) {
    u_grp <- unique(non_diag_df$group)
    sliced_list <- lapply(u_grp, function(g) {
      utils::head(
        non_diag_df[non_diag_df$group == g, c("group", "ion"), drop = FALSE],
        ions_max
      )
    })
    non_diag_df <- do.call(rbind, sliced_list)
  }
  ions_list <- extract_ions_by_group(non_diag_df)
  rm(ions_df, diag_df, non_diag_df)

  if (length(ions_list) == 0L) {
    final_queries <- data.frame(
      skeleton = character(0),
      mcc = numeric(0),
      combined_ions = character(0),
      stringsAsFactors = FALSE
    )
    if (!is.null(export)) {
      create_dir(export)
      utils::write.table(
        final_queries,
        file = export,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
    return(invisible(final_queries))
  }

  message("Generate all combinations of queries.")
  combinations <- names(ions_list) |>
    generate_combinations_progress(
      ions_list = ions_list,
      max_ions = ions_max,
      show_progress = TRUE
    )
  names(combinations) <- names(ions_list)

  message("Merging diagnostic ions with combinations.")
  has_diagnostic <- match(
    names(ions_list),
    names(ions_list_diagnostic),
    nomatch = 0L
  ) >
    0L

  if (any(has_diagnostic)) {
    for (name in names(combinations)[has_diagnostic]) {
      x <- combinations[[name]]
      diagnostic_ions <- ions_list_diagnostic[[name]]

      if (length(diagnostic_ions) > 0L) {
        combinations[[name]] <- lapply(x, function(sublist) {
          unique(c(unlist(sublist), diagnostic_ions))
        })
      }
    }
  }

  rm(ions_list, ions_list_diagnostic)

  all_combinations <- unlist(combinations, recursive = FALSE)
  names(all_combinations) <- gsub("\\d+$", "", names(all_combinations))
  rm(combinations)

  message("Test the queries. (This is the longest step)")
  queries_results <- all_combinations |>
    perform_list_of_queries_progress(
      spectra = mia_spectra,
      dalton = dalton,
      ppm = ppm,
      show_progress = TRUE
    )
  names(queries_results) <- names(all_combinations)

  message("Evaluate the performance of each query using MCC.")
  normalized_skeleton <- gsub(
    pattern = "+",
    replacement = ".",
    x = mia_spectra$SKELETON,
    fixed = TRUE
  )
  skeleton_counts <- table(normalized_skeleton)

  results_stats <- vapply(
    seq_along(queries_results),
    function(result) {
      target_name <- names(queries_results)[result]
      vals <- queries_results[[result]]$value
      tp <- as.numeric(sum(vals == target_name))
      fp <- as.numeric(length(vals)) - tp
      group_total <- as.numeric(skeleton_counts[[target_name]])
      fn_ <- group_total - tp
      tn <- total_spectra - tp - fp - fn_

      denom <- sqrt((tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_))
      mcc <- if (denom > 0) (tp * tn - fp * fn_) / denom else 0.0
      round(mcc, decimals)
    },
    numeric(1L)
  )

  valid_mcc <- !is.nan(results_stats)
  skel_vec <- names(all_combinations)[valid_mcc]
  mcc_vec <- results_stats[valid_mcc]
  ions_vec <- all_combinations[valid_mcc]

  ord <- order(mcc_vec, decreasing = TRUE)
  skel_vec <- skel_vec[ord]
  mcc_vec <- mcc_vec[ord]
  ions_vec <- ions_vec[ord]

  max_mcc_by_skel <- tapply(mcc_vec, skel_vec, max)
  keep_best <- mcc_vec == as.numeric(max_mcc_by_skel[skel_vec])

  skel_best <- skel_vec[keep_best]
  mcc_best <- mcc_vec[keep_best]
  ions_best <- ions_vec[keep_best]

  u_skel <- unique(skel_best)
  final_queries <- data.frame(
    skeleton = u_skel,
    mcc = vapply(
      u_skel,
      function(s) mcc_best[skel_best == s][1L],
      numeric(1L),
      USE.NAMES = FALSE
    ),
    combined_ions = vapply(
      u_skel,
      function(s) combine_ions_minimal(ions_best[skel_best == s]),
      character(1L),
      USE.NAMES = FALSE
    ),
    stringsAsFactors = FALSE
  )
  final_queries <- final_queries[
    order(final_queries$mcc, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  rownames(final_queries) <- NULL

  message("Export the best queries for further use.")
  if (!is.null(export)) {
    create_dir(export)
    utils::write.table(
      final_queries,
      file = export,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  invisible(final_queries)
}
