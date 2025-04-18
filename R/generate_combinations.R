#' @title Generate combinations
#'
#' @param x x
#' @param max_ions Max ions
#'
#' @return NULL
#'
#' @examples NULL
generate_combinations <- function(x, max_ions) {
  1:min(length(x), max_ions) |>
    purrr::map(
      .f = function(k) {
        utils::combn(x, k, simplify = FALSE)
      }
    ) |>
    unlist(recursive = FALSE)
}

#' @title Perform list of queries (progress)
#'
#' @param indices Indices
#' @param ions_list Ions list
#' @param max_ions Max ions
#'
#' @return NULL
#'
#' @examples NULL
generate_combinations_progress <- function(indices, ions_list, max_ions) {
  purrr::map(
    .progress = TRUE,
    .x = indices,
    .f = function(index, ions_list, max_ions) {
      x <- ions_list[[index]]
      generate_combinations(x = x, max_ions = max_ions)
    },
    ions_list = ions_list,
    max_ions = max_ions
  )
}
