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
  spectra_mat <- matrix(0, nrow = n_spectra, ncol = n_peaks)

  # Fill matrix row by row (more memory efficient than rbind)
  for (i in seq_len(n_spectra)) {
    spectra_mat[i, ] <- intensity_list[[i]]
  }

  # Filter out zero columns
  zeros <- colSums(spectra_mat) <= zero_val
  if (any(!zeros)) {
    spectra_mat <- spectra_mat[, !zeros, drop = FALSE]
    colnames(spectra_mat) <- Spectra::mz(
      spectra,
      BPPARAM = BiocParallel::SerialParam()
    )[[1L]][!zeros]
  } else {
    # Handle edge case where all columns are zeros
    spectra_mat <- matrix(0, nrow = n_spectra, ncol = 0)
  }

  rownames(spectra_mat) <- name
  return(spectra_mat)
}
