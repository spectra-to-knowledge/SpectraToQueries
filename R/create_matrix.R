#' @title Create matrix
#'
#' @param spectra Spectra
#' @param zero_val Zero value
#' @param name Name
#'
#' @return A matrix with spectra intensity data
#'
#' @examples NULL
create_matrix <- function(spectra, zero_val = 0, name) {
  # Get intensity data
  intensity_list <- Spectra::intensity(
    spectra,
    BPPARAM = BiocParallel::SerialParam()
  )

  # Pre-allocate matrix
  n_spectra <- length(intensity_list)
  n_peaks <- length(intensity_list[[1L]])

  # First pass: identify non-zero columns to avoid storing unnecessary data
  col_sums <- numeric(n_peaks)
  for (i in seq_len(n_spectra)) {
    col_sums <- col_sums + intensity_list[[i]]
  }

  keep_cols <- col_sums > zero_val
  n_keep <- sum(keep_cols)

  if (n_keep == 0L) {
    # Handle edge case where all columns are zeros
    spectra_mat <- matrix(0, nrow = n_spectra, ncol = 0)
    rownames(spectra_mat) <- name
    return(spectra_mat)
  }

  # Create smaller matrix with only non-zero columns
  spectra_mat <- matrix(0, nrow = n_spectra, ncol = n_keep)

  # Fill matrix with only kept columns
  for (i in seq_len(n_spectra)) {
    spectra_mat[i, ] <- intensity_list[[i]][keep_cols]
  }

  # Get mz values only once and filter
  mz_vals <- Spectra::mz(
    spectra,
    BPPARAM = BiocParallel::SerialParam()
  )[[1L]][keep_cols]

  colnames(spectra_mat) <- mz_vals
  rownames(spectra_mat) <- name

  return(spectra_mat)
}
