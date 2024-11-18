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
perform_query <- function(spectra, frags, nls, dalton = 0.01, ppm = 25) {
  if (length(frags) != 0) {
    spectra <- spectra[spectra |>
      Spectra::containsMz(
        mz = frags,
        tolerance = dalton,
        ppm = ppm,
        which = "all"
      )]
  }
  if (length(nls) != 0) {
    if (length(spectra) != 0) {
      # TODO report issue not taking "all" argument for now
      for (nl in nls) {
        if (length(spectra) != 0) {
          spectra <- spectra[spectra |>
            Spectra::containsNeutralLoss(
              neutralLoss = nl,
              tolerance = dalton,
              ppm = ppm
            )]
        }
      }
    }
  }
  # message("Query done.")
  return(spectra)
}
