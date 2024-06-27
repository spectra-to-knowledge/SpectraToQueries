#' Title
#'
#' @param spectra Spectra
#' @param dalton Dalton
#'
#' @return NULL
#' @export
#'
#' @examples NULL
harmonize_mzs <- function(spectra, dalton) {
  spectra_new <- spectra
  averaged_intensities <- spectra_new |>
    Spectra::peaksData() |>
    Spectra::combinePeaksData(tolerance = dalton, peaks = "union") |>
    data.frame() |>
    tidytable::pull("mz")
  for (i in seq_along(1:length(spectra_new))) {
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
