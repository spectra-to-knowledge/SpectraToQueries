#' Title
#'
#' @param x x
#'
#' @return NULL
#' @export
#'
#' @examples NULL
generate_combinations <- function(x, max_ions) {
  1:min(length(x), max_ions) |>
    lapply(
      FUN = function(k) {
        combn(x, k, simplify = FALSE)
      }
    ) |>
    purrr::flatten()
}
