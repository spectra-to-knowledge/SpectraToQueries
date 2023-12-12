#' Title
#'
#' @param spectra
#' @param zero_val
#' @param name
#'
#' @return
#' @export
#'
#' @examples
create_matrix <- function(spectra, zero_val = ZERO_VAL, name) {
  spectra_mat <- do.call(
    BiocGenerics::rbind,
    Spectra::intensity(spectra) |>
      as.list()
  )
  zeros <- colSums(spectra_mat) <= zero_val
  spectra_mat <- spectra_mat[, !zeros]
  colnames(spectra_mat) <-
    Spectra::mz(spectra)[[1L]][!zeros]
  rownames(spectra_mat) <- name
  return(spectra_mat)
}
