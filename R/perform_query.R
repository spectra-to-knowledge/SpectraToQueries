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
perform_query <- function(spectra,
                          frags,
                          nls,
                          dalton = 0.01,
                          ppm = 25) {
  if (length(frags) > 0) {
    spectra <- spectra[Spectra::containsMz(
      spectra,
      mz = frags,
      tolerance = dalton,
      ppm = ppm,
      which = "all"
    )]
  }

  # Return early if spectra is empty after filtering
  if (length(spectra) == 0) {
    # message("Remaining spectra empty, very early return")
    return(spectra)
  }

  if (length(nls) > 0) {
    # TODO report issue not taking "all" argument for now
    for (nl in nls) {
      spectra <- spectra[Spectra::containsNeutralLoss(
        spectra,
        neutralLoss = nl,
        tolerance = dalton,
        ppm = ppm
      )]
      # Return early if spectra becomes empty during the loop
      if (length(spectra) == 0) {
        # message("Remaining spectra empty, early return")
        return(spectra)
      }
    }
  }
  return(spectra)
}
