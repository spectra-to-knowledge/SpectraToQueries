#' @title Filter matrix
#'
#' @param matrix Matrix to filter
#' @param n Minimum number of non-zero values per column
#'
#' @return Filtered matrix with columns having at least `n` non-zero values
#'
#' @examples NULL
filter_matrix <- function(matrix, n) {
  keep_cols <- colSums(matrix > 0) >= n
  matrix[, keep_cols, drop = FALSE]
}
