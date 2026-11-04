#' @title Harmonize mzs
#'
#' @param spectra Spectra object
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Spectra with harmonized m/z values
#'
#' @examples NULL
harmonize_mzs <- function(spectra, dalton, ppm) {
  spectra_new <- spectra

  combined_peaks <- spectra_new |>
    Spectra::peaksData() |>
    Spectra::combinePeaksData(
      tolerance = dalton,
      ppm = ppm,
      peaks = "union"
    )
  averaged_intensities <- combined_peaks[, "mz"]
  n_avg <- length(averaged_intensities)

  if (n_avg == 0L) {
    return(spectra_new)
  }

  pks <- spectra_new@backend@peaksData
  for (i in seq_along(pks)) {
    mzs <- pks[[i]][, 1L]
    if (length(mzs) > 0L) {
      idx_left <- findInterval(mzs, averaged_intensities)
      idx_left_clamped <- pmax(1L, idx_left)
      idx_right <- pmin(n_avg, idx_left + 1L)

      d_left <- abs(mzs - averaged_intensities[idx_left_clamped])
      d_right <- abs(mzs - averaged_intensities[idx_right])

      use_right <- d_right < d_left
      min_diff <- ifelse(use_right, d_right, d_left)
      best_idx <- ifelse(use_right, idx_right, idx_left_clamped)

      within <- min_diff <= dalton
      if (any(within)) {
        mzs[within] <- averaged_intensities[best_idx[within]]
        pks[[i]][, 1L] <- mzs
      }
    }
  }
  spectra_new@backend@peaksData <- pks

  spectra_new
}
