#' @title Perform query
#'
#' @param spectra Spectra object
#' @param frags Fragment masses to search for
#' @param nls Neutral losses to search for
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Filtered spectra object
#'
#' @examples NULL
perform_query <- function(spectra, frags, nls, dalton, ppm) {
  # Early return if no spectra
  if (length(spectra) == 0) {
    return(spectra)
  }

  normalize_idx <- function(idx, n) {
    # Spectra::containsMz can return logical vector or matrix
    if (is.matrix(idx)) {
      idx <- as.vector(idx)
    }
    # Ensure length matches spectra
    if (length(idx) != n) {
      idx <- rep(FALSE, n)
    }
    idx
  }

  # Filter by fragments if any are specified
  if (length(frags) > 0) {
    idx <- Spectra::containsMz(
      spectra,
      mz = frags,
      tolerance = dalton,
      ppm = ppm,
      which = "all",
      BPPARAM = BiocParallel::SerialParam()
    )

    # Early return if no spectra remain after fragment filtering
    spectra <- spectra[normalize_idx(idx, length(spectra))]
    if (length(spectra) == 0) {
      return(spectra)
    }
  }

  # Filter by neutral losses if any are specified
  if (length(nls) > 0) {
    mzs <- (Spectra::precursorMz(spectra) - nls) |>
      unique() |>
      sort()
    idx <- Spectra::containsMz(
      spectra,
      mz = mzs,
      tolerance = dalton,
      ppm = ppm,
      which = "all",
      BPPARAM = BiocParallel::SerialParam()
    )
    spectra <- spectra[normalize_idx(idx, length(spectra))]
  }

  return(spectra)
}
