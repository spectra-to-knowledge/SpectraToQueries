#' @title Filter matrix
#'
#' @param matrix Matrix to filter
#' @param n Minimum number of non-zero values per column
#'
#' @return Filtered matrix with columns having at least n non-zero values
#'
#' @examples NULL
filter_matrix <- function(matrix, n) {
  non_null_count <- colSums(matrix > 0)
  keep_cols <- non_null_count >= n
  filtered_matrix <- matrix[, keep_cols, drop = FALSE]
  return(filtered_matrix)
}
