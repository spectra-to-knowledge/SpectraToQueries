#' Convert Spectra object to portable data.frame
#'
#' Converts a [Spectra::Spectra] object to a portable `data.frame` format where
#' metadata columns are preserved and peaks (`mz`/`intensity`) are added as list
#' columns. [Spectra::Spectra()] converts the result back.
#'
#' @param spectra A [Spectra::Spectra] object
#'
#' @return A `data.frame` with metadata columns plus `mz` and `intensity`
#'   as list columns (one entry per spectrum)
#'
#' @keywords internal
#' @noRd
spectra_to_df <- function(spectra) {
  meta_df <- Spectra::spectraData(spectra) |>
    as.data.frame()
  if ("dataOrigin" %in% names(meta_df) && is.character(meta_df$dataOrigin)) {
    non_na <- !is.na(meta_df$dataOrigin)
    meta_df$dataOrigin[non_na] <- basename(meta_df$dataOrigin[non_na])
  }
  peaks_data <- Spectra::peaksData(spectra)

  meta_df$mz <- lapply(peaks_data, \(x) x[, "mz"])
  meta_df$intensity <- lapply(peaks_data, \(x) x[, "intensity"])

  rownames(meta_df) <- NULL
  meta_df
}
