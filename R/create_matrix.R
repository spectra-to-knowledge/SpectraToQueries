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
    spectra
  )

  # Pre-allocate matrix
  n_spectra <- length(intensity_list)
  n_peaks <- length(intensity_list[[1L]])

  # Use colSums approach: create matrix directly and filter
  # Initialize matrix with all values
  spectra_mat <- matrix(0, nrow = n_spectra, ncol = n_peaks)

  # Fill matrix row by row (more cache-friendly than column-wise)
  for (i in seq_len(n_spectra)) {
    spectra_mat[i, ] <- intensity_list[[i]]
  }

  # Vectorized: identify non-zero columns using colSums
  col_sums <- colSums(spectra_mat)
  keep_cols <- col_sums > zero_val

  if (!any(keep_cols)) {
    # Handle edge case where all columns are zeros
    spectra_mat <- matrix(0, nrow = n_spectra, ncol = 0)
    rownames(spectra_mat) <- name
    return(spectra_mat)
  }

  # Subset matrix using logical indexing (vectorized)
  spectra_mat <- spectra_mat[, keep_cols, drop = FALSE]

  # Get mz values only once and filter
  mz_vals <- Spectra::mz(
    spectra
  )[[1L]][keep_cols]

  colnames(spectra_mat) <- mz_vals
  rownames(spectra_mat) <- name

  return(spectra_mat)
}
