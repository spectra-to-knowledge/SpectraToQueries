#' @title Abort with stq_input_error condition
#' @keywords internal
#' @noRd
abort_input <- function(msg) {
  stop(errorCondition(msg, class = "stq_input_error"))
}

#' @title Validate a scalar finite numeric argument
#' @keywords internal
#' @noRd
check_num <- function(x, msg, lower = -Inf, upper = Inf) {
  ok <- is.numeric(x) &&
    length(x) == 1L &&
    is.finite(x) &&
    x >= lower &&
    x <= upper
  if (!ok) {
    abort_input(msg)
  }
}

#' Convert Spectra object to portable data.frame
#'
#' Converts a [Spectra::Spectra] object to a portable `data.frame` format where
#' metadata columns are preserved and peaks (`mz`/`intensity`) are added as list
#' columns. This format is ideal for export to portable formats and can be
#' converted back to a [Spectra::Spectra] object via `df_to_spectra()`.
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
  peaks_data <- Spectra::peaksData(spectra)

  meta_df$mz <- lapply(peaks_data, \(x) x[, "mz"])
  meta_df$intensity <- lapply(peaks_data, \(x) x[, "intensity"])

  rownames(meta_df) <- NULL
  meta_df
}

#' Convert data.frame back to Spectra object
#'
#' Reconstructs a [Spectra::Spectra] object from a portable `data.frame` format
#' created by `spectra_to_df()`. The `data.frame` must have `mz` and `intensity`
#' columns (as lists) along with metadata columns.
#'
#' @param df A `data.frame` with metadata columns plus `mz` and `intensity`
#'   as list columns
#'
#' @return A [Spectra::Spectra] object
#'
#' @keywords internal
#' @noRd
df_to_spectra <- function(df) {
  Spectra::Spectra(df)
}
