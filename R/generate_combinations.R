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
    furrr::future_map(
      .f = function(k) {
        combn(x, k, simplify = FALSE)
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
  p <- progressr::progressor(along = indices)
  furrr::future_map(
    .x = indices,
    ions_list = ions_list,
    max_ions = max_ions,
    .f = function(x, ions_list, max_ions) {
      p()
      generate_combinations(x = ions_list[[x]], max_ions = max_ions)
    }
  )
}
