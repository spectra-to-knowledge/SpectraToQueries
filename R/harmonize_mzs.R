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

  # Averaged m/z values from combined spectra
  averaged_intensities <- spectra_new |>
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
