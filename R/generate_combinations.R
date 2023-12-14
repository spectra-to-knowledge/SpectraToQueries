#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
generate_combinations <- function(x) {
  1:length(x) |>
    lapply(
      FUN = function(k) {
        combn(x, k, simplify = FALSE)
      }
    ) |>
    purrr::flatten()
}
