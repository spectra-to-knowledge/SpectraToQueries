#' Title
#'
#' @param spectra
#' @param dalton
#'
#' @return
#' @export
#'
#' @examples
harmonize_mzs <- function(spectra, dalton) {
  averaged_intensities <- spectra |>
    Spectra::peaksData() |>
    Spectra::combinePeaksData(tolerance = dalton, peaks = "union") |>
    data.frame() |>
    tidytable::pull("mz")

  for (i in seq_along(1:length(spectra))) {
    mzs <- spectra@backend@peaksData[[i]] |>
      data.frame() |>
      tidytable::pull("mz")

    for (j in seq_along(mzs)) {
      matching_index <- which.min(abs(mzs[j] - averaged_intensities))
      if (abs(mzs[[j]] - averaged_intensities[matching_index]) <= dalton) {
        mzs[j] <- averaged_intensities[matching_index]
      }
    }
    spectra@backend@peaksData[[i]][, 1] <- mzs
  }
  return(spectra)
}
