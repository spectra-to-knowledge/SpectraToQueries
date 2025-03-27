#' @title Perform query
#'
#' @param spectra Spectra
#' @param frags Fragments
#' @param nls Neutral losses
#' @param dalton Dalton
#' @param ppm PPM
#'
#' @return NULL
#'
#' @examples NULL
perform_query <- function(spectra, frags, nls, dalton, ppm) {
  if (length(frags) > 0) {
    spectra <- spectra[Spectra::containsMz(
      spectra,
      mz = frags,
      tolerance = dalton,
      ppm = ppm,
      which = "all",
      BPPARAM = BiocParallel::SerialParam()
    )]
  }

  # Return early if spectra is empty after filtering
  if (length(spectra) == 0) {
    return(spectra)
  }

  if (length(nls) > 0) {
    mzs <- (Spectra::precursorMz(spectra) - nls) |>
      unique() |>
      sort()
    spectra <- spectra[Spectra::containsMz(
      spectra,
      mz = mzs,
      tolerance = dalton,
      ppm = ppm,
      which = "all",
      BPPARAM = BiocParallel::SerialParam()
    )]
  }
  return(spectra)
}
