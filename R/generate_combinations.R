#' Title
#'
#' @param x x
#'
#' @return NULL
#' @export
#'
#' @examples NULL
generate_combinations <- function(x) {
  1:length(x) |>
    lapply(
      FUN = function(k) {
        combn(x, k, simplify = FALSE)
      }
    ) |>
    purrr::flatten()
}
