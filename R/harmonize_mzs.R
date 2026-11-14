#' @title Harmonize mzs
#'
#' @description Aligns peak m/z values across spectra to intensity-weighted
#'   consensus centroids within a combined `dalton + mz * ppm / 1e6` window.
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
  if (length(spectra_new) == 0L) {
    return(spectra_new)
  }

  combined_peaks <- spectra_new |>
    Spectra::peaksData() |>
    Spectra::combinePeaksData(
      tolerance = dalton,
      ppm = ppm,
      peaks = "union",
      weighted = TRUE
    )
  consensus_mzs <- combined_peaks[, "mz"]
  consensus_mzs <- consensus_mzs[is.finite(consensus_mzs)]
  n_avg <- length(consensus_mzs)

  if (n_avg == 0L) {
    return(spectra_new)
  }

  align_mzs <- function(x, ...) {
    if (nrow(x) > 0L) {
      mzs <- x[, 1L]
      idx_left <- findInterval(mzs, consensus_mzs)
      idx_left_clamped <- pmax(1L, idx_left)
      idx_right <- pmin(n_avg, idx_left + 1L)

      d_left <- abs(mzs - consensus_mzs[idx_left_clamped])
      d_right <- abs(mzs - consensus_mzs[idx_right])

      use_right <- d_right < d_left
      min_diff <- ifelse(use_right, d_right, d_left)
      best_idx <- ifelse(use_right, idx_right, idx_left_clamped)

      target_mzs <- consensus_mzs[best_idx]
      within <- min_diff <= (dalton + target_mzs * (ppm / 1e6))
      if (any(within)) {
        mzs[within] <- target_mzs[within]
        x[, 1L] <- mzs
      }
    }
    x
  }

  sp_aligned <- Spectra::addProcessing(spectra_new, align_mzs)
  tryCatch(
    Spectra::applyProcessing(sp_aligned),
    error = function(e) {
      Spectra::setBackend(sp_aligned, Spectra::MsBackendMemory()) |>
        Spectra::applyProcessing()
    }
  )
}
