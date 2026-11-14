#' @title Create matrix
#'
#' @param spectra Spectra object (binned)
#' @param zero_val Zero value threshold for filtering columns
#' @param name Row names for the resulting matrix
#'
#' @return A matrix with spectra intensity data and m/z column names
#'
#' @examples NULL
#' @keywords internal
#' @noRd
create_matrix <- function(spectra, zero_val = 0, name) {
  pks <- Spectra::peaksData(spectra)
  n_spectra <- length(pks)

  lens <- vapply(pks, nrow, integer(1L))
  if (all(lens == 0L)) {
    spectra_mat <- matrix(0, nrow = n_spectra, ncol = 0L)
    rownames(spectra_mat) <- name
    return(spectra_mat)
  }

  # Fast column-major path when all spectra share the exact same sorted m/z grid
  ref_mz <- pks[[1L]][, 1L]
  if (
    length(unique(lens)) == 1L &&
      lens[1L] > 0L &&
      !is.unsorted(ref_mz, strictly = TRUE) &&
      all(vapply(
        pks[-1L],
        function(x) identical(x[, 1L], ref_mz),
        logical(1L)
      ))
  ) {
    all_int <- unlist(lapply(pks, function(x) x[, 2L]), use.names = FALSE)
    dim(all_int) <- c(lens[1L], n_spectra)
    keep_cols <- rowSums(all_int) > zero_val
    if (!any(keep_cols)) {
      spectra_mat <- matrix(0, nrow = n_spectra, ncol = 0L)
      rownames(spectra_mat) <- name
      return(spectra_mat)
    }
    spectra_mat <- t(all_int[keep_cols, , drop = FALSE])
    colnames(spectra_mat) <- pks[[1L]][keep_cols, 1L]
    rownames(spectra_mat) <- name
    return(spectra_mat)
  }

  # Sparse path (e.g. binned with zero.rm = TRUE)
  all_mz <- unlist(lapply(pks, function(x) x[, 1L]), use.names = FALSE)
  all_int <- unlist(lapply(pks, function(x) x[, 2L]), use.names = FALSE)
  u_mz <- sort(unique(all_mz))
  col_idx <- match(all_mz, u_mz)
  row_idx <- rep.int(seq_len(n_spectra), lens)

  spectra_mat <- matrix(0, nrow = n_spectra, ncol = length(u_mz))
  has_dups <- any(vapply(
    pks,
    function(x) anyDuplicated(x[, 1L]) > 0L,
    logical(1L)
  ))
  if (has_dups) {
    lin_idx <- row_idx + (col_idx - 1L) * n_spectra
    dup_sums <- rowsum(all_int, lin_idx, reorder = FALSE)
    spectra_mat[as.integer(rownames(dup_sums))] <- dup_sums[, 1L]
  } else {
    spectra_mat[cbind(row_idx, col_idx)] <- all_int
  }

  keep_cols <- colSums(spectra_mat) > zero_val
  if (!any(keep_cols)) {
    spectra_mat <- matrix(0, nrow = n_spectra, ncol = 0L)
    rownames(spectra_mat) <- name
    return(spectra_mat)
  }

  spectra_mat <- spectra_mat[, keep_cols, drop = FALSE]
  colnames(spectra_mat) <- u_mz[keep_cols]
  rownames(spectra_mat) <- name

  spectra_mat
}
