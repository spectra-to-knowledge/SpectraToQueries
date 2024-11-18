#' @title Generate combinations
#'
#' @param x x
#' @param max_ions MAx ions
#'
#' @return NULL
#'
#' @examples NULL
generate_combinations <- function(x, max_ions) {
  1:min(length(x), max_ions) |>
    furrr::future_map(
      .f = function(k) {
        combn(x, k, simplify = FALSE)
      }
    ) |>
    unlist(recursive = FALSE)
}
