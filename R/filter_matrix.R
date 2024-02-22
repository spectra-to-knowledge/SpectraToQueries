#' Title
#'
#' @param matrix Matrix
#' @param n N
#'
#' @return NULL
#' @export
#'
#' @examples NULL
filter_matrix <- function(matrix, n) {
  non_null_count <-
    apply(
      X = matrix,
      MARGIN = 2,
      FUN = function(column) {
        sum(column > 0)
      }
    )
  filtered_matrix <- matrix[, non_null_count >= n, drop = FALSE]

  return(filtered_matrix)
}
