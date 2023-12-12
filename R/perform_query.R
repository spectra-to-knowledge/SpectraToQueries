#' Title
#'
#' @param spectra
#' @param frags
#' @param nls
#'
#' @return
#' @export
#'
#' @examples
perform_query <- function(spectra, frags, nls) {
  if (length(frags) != 0) {
    spectra <- spectra[spectra |>
      Spectra::containsMz(
        mz = frags,
        tolerance = DALTON,
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
              tolerance = DALTON
            )]
        }
      }
    }
  }
  message("Query done.")
  return(spectra)
}
