#' @title Fix binned mzs
#'
#' @param binned_m Binned matrix
#' @param original_mzs Original m/z values from spectra
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#' @param decimals Number of decimal places for rounding
#'
#' @return Matrix with corrected m/z column names
#'
#' @examples NULL
fix_binned_mzs <- function(binned_m, original_mzs, dalton, ppm, decimals) {
  all_mzs <- original_mzs |>
    Spectra::peaksData(
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    Spectra::combinePeaksData(
      tolerance = dalton,
      ppm = ppm,
      peaks = "union",
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    data.frame() |>
    tidytable::pull("mz")

  old_colnames <- as.numeric(colnames(binned_m))

  # Pre-compute ppm tolerance component
  ppm_component <- all_mzs * ppm / 1e6
  tolerance_vec <- dalton + ppm_component

  # Vectorized: compute distance matrix for all old_colnames at once
  diffs_matrix <- abs(outer(old_colnames, all_mzs, "-"))

  # For each old column, find matching m/z values within tolerance and average
  new_colnames <- numeric(length(old_colnames))
  for (i in seq_along(old_colnames)) {
    within_tolerance <- diffs_matrix[i, ] <= tolerance_vec
    if (any(within_tolerance)) {
      new_colnames[i] <- mean(all_mzs[within_tolerance])
    } else {
      new_colnames[i] <- old_colnames[i]
    }
  }

  # Round new column names
  new_colnames_rounded <- round(new_colnames, decimals)

  # Use tidytable for efficient grouping of columns with same rounded m/z
  col_mapping <- tidytable::tidytable(
    old_idx = seq_along(new_colnames_rounded),
    rounded = new_colnames_rounded
  ) |>
    tidytable::group_by(rounded) |>
    tidytable::summarise(
      indices = list(old_idx),
      .by = rounded
    ) |>
    tidytable::ungroup()

  # Create result matrix with merged columns
  unique_cols <- col_mapping$rounded
  result_matrix <- matrix(
    0,
    nrow = nrow(binned_m),
    ncol = length(unique_cols)
  )
  rownames(result_matrix) <- rownames(binned_m)

  # Merge columns: sum values for columns with same rounded m/z
  for (i in seq_len(nrow(col_mapping))) {
    matching_cols <- col_mapping$indices[[i]]
    if (length(matching_cols) == 1L) {
      result_matrix[, i] <- binned_m[, matching_cols[[1L]]]
    } else {
      result_matrix[, i] <- rowSums(binned_m[, matching_cols, drop = FALSE])
    }
  }

  colnames(result_matrix) <- as.character(unique_cols)

  return(result_matrix)
}
