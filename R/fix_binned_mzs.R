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

  # Vectorized calculation of new column names
  old_colnames <- as.numeric(colnames(binned_m))

  # Pre-allocate result vector
  new_colnames <- numeric(length(old_colnames))

  # Vectorized approach: for each old column, find all matching mzs and average
  for (i in seq_along(old_colnames)) {
    # Find matching m/z values within tolerance
    within_tolerance <- abs(all_mzs - old_colnames[i]) <= dalton
    if (any(within_tolerance)) {
      new_colnames[i] <- mean(all_mzs[within_tolerance])
    } else {
      new_colnames[i] <- old_colnames[i]
    }
  }

  # Round new column names
  new_colnames_rounded <- round(new_colnames, decimals)

  # Group columns with same rounded m/z more efficiently
  # Use aggregate instead of multiple pivot operations
  unique_cols <- unique(new_colnames_rounded)

  if (length(unique_cols) < ncol(binned_m)) {
    # Need to merge some columns
    result_matrix <- matrix(
      0,
      nrow = nrow(binned_m),
      ncol = length(unique_cols)
    )
    rownames(result_matrix) <- rownames(binned_m)

    # Sum values for columns with same rounded m/z
    for (i in seq_along(unique_cols)) {
      matching_cols <- which(new_colnames_rounded == unique_cols[i])
      if (length(matching_cols) == 1) {
        result_matrix[, i] <- binned_m[, matching_cols]
      } else {
        result_matrix[, i] <- rowSums(binned_m[, matching_cols, drop = FALSE])
      }
    }

    colnames(result_matrix) <- as.character(unique_cols)
  } else {
    # No merging needed, just update column names
    result_matrix <- binned_m
    colnames(result_matrix) <- as.character(new_colnames_rounded)
  }

  return(result_matrix)
}
