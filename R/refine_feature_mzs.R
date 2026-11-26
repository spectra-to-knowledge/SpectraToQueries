#' @title Refine feature matrix m/z values
#'
#' @description Replaces each column m/z `c` of a spectra x m/z feature matrix
#'   by the weighted mean of the consensus peaks of `spectra` (from
#'   [Spectra::combinePeaksData()]) within `dalton + mz * ppm / 1e6` of `c`,
#'   with weights `max(0, intensity) + 1`; the `+ 1` keeps zero-intensity peaks
#'   in the mean. Columns whose refined m/z values round to the same value at
#'   `decimals` are summed.
#'
#' @param feature_mat Spectra x m/z intensity matrix
#' @param spectra `Spectra` object the matrix was built from
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#' @param decimals Number of decimal places for rounding
#'
#' @return Matrix with corrected m/z column names
#'
#' @keywords internal
#' @noRd
refine_feature_mzs <- function(feature_mat, spectra, dalton, ppm, decimals) {
  if (ncol(feature_mat) == 0L) {
    return(feature_mat)
  }

  if (length(spectra) > 0L) {
    combined_peaks <- spectra |>
      Spectra::peaksData() |>
      Spectra::combinePeaksData(
        tolerance = dalton,
        ppm = ppm,
        peaks = "union",
        weighted = TRUE
      )
    keep_fin <- is.finite(combined_peaks[, "mz"]) &
      is.finite(combined_peaks[, "intensity"])
    all_mzs <- combined_peaks[keep_fin, "mz"]
    all_ints <- combined_peaks[keep_fin, "intensity"]
  } else {
    all_mzs <- numeric(0)
    all_ints <- numeric(0)
  }

  old_colnames <- as.numeric(colnames(feature_mat))
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
          hit_idx <- idx[within]
          w <- pmax(0, all_ints[hit_idx]) + 1
          new_colnames[i] <- sum(all_mzs[hit_idx] * w) / sum(w)
        }
      }
    }
  }

  new_colnames_rounded <- round(new_colnames, decimals)

  # Fast path when all rounded column names are already unique
  if (!anyDuplicated(new_colnames_rounded)) {
    colnames(feature_mat) <- as.character(new_colnames_rounded)
    return(feature_mat)
  }

  t(rowsum(t(feature_mat), new_colnames_rounded, reorder = FALSE))
}
