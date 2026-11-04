#' @title Fix binned mzs
#'
#' @param binned_m Binned matrix
#' @param original_mzs Original `Spectra` object with unbinned m/z values
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#' @param decimals Number of decimal places for rounding
#'
#' @return Matrix with corrected m/z column names
#'
#' @examples NULL
fix_binned_mzs <- function(binned_m, original_mzs, dalton, ppm, decimals) {
  if (ncol(binned_m) == 0L) {
    return(binned_m)
  }

  all_mzs <- if (length(original_mzs) > 0L) {
    combined_peaks <- original_mzs |>
      Spectra::peaksData() |>
      Spectra::combinePeaksData(
        tolerance = dalton,
        ppm = ppm,
        peaks = "union"
      )
    mz_vals <- combined_peaks[, "mz"]
    mz_vals[is.finite(mz_vals)]
  } else {
    numeric(0)
  }

  old_colnames <- as.numeric(colnames(binned_m))
  new_colnames <- old_colnames

  if (length(all_mzs) > 0L) {
    tolerance_vec <- dalton + all_mzs * (ppm / 1e6)
    max_tol <- max(tolerance_vec)
    left_bounds <- pmax(1L, findInterval(old_colnames - max_tol, all_mzs))
    right_bounds <- findInterval(old_colnames + max_tol, all_mzs)

    for (i in seq_along(old_colnames)) {
      lo <- left_bounds[i]
      hi <- right_bounds[i]
      if (hi >= lo) {
        idx <- lo:hi
        within <- abs(all_mzs[idx] - old_colnames[i]) <= tolerance_vec[idx]
        if (any(within)) {
          new_colnames[i] <- mean(all_mzs[idx[within]])
        }
      }
    }
  }

  new_colnames_rounded <- round(new_colnames, decimals)
  unique_cols <- unique(new_colnames_rounded)

  # Fast path when all rounded column names are already unique
  if (length(unique_cols) == length(new_colnames_rounded)) {
    result_matrix <- binned_m
    colnames(result_matrix) <- as.character(unique_cols)
    return(result_matrix)
  }

  col_groups <- split(
    seq_along(new_colnames_rounded),
    factor(new_colnames_rounded, levels = unique_cols)
  )

  result_matrix <- matrix(
    0,
    nrow = nrow(binned_m),
    ncol = length(unique_cols),
    dimnames = list(rownames(binned_m), as.character(unique_cols))
  )

  for (i in seq_along(col_groups)) {
    matching_cols <- col_groups[[i]]
    if (length(matching_cols) == 1L) {
      result_matrix[, i] <- binned_m[, matching_cols]
    } else {
      result_matrix[, i] <- rowSums(binned_m[, matching_cols, drop = FALSE])
    }
  }

  result_matrix
}
