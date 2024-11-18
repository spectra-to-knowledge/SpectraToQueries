#' @title Harmonize mzs
#'
#' @param spectra Spectra
#' @param dalton Dalton
#' @param ppm PPM
#'
#' @return NULL
#'
#' @examples NULL
harmonize_mzs <- function(spectra, dalton, ppm) {
  spectra_new <- spectra
  averaged_intensities <- spectra_new |>
    Spectra::peaksData() |>
    Spectra::combinePeaksData(tolerance = dalton, ppm = ppm, peaks = "union") |>
    data.frame() |>
    tidytable::pull("mz")
  for (i in seq_along(seq_along(spectra_new))) {
    mzs <- spectra_new@backend@peaksData[[i]] |>
      data.frame() |>
      tidytable::pull("mz")

    for (j in seq_along(mzs)) {
      matching_index <- which.min(abs(mzs[j] - averaged_intensities))
      if (abs(mzs[[j]] - averaged_intensities[matching_index]) <= dalton) {
        mzs[j] <- averaged_intensities[matching_index]
      }
    }
    spectra_new@backend@peaksData[[i]][, 1] <- mzs
  }
  return(spectra_new)
}
