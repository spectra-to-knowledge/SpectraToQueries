#' @title Generate combinations
#'
#' @param x Vector of elements to combine
#' @param max_ions Maximum number of ions in combinations
#'
#' @return List of all combinations
#'
#' @examples NULL
generate_combinations <- function(x, max_ions) {
  n <- length(x)
  if (n == 0) {
    return(list())
  }
  
  max_k <- min(n, max_ions)
  
  # Pre-allocate result list with estimated size
  estimated_size <- sum(choose(n, seq_len(max_k)))
  result <- vector("list", estimated_size)
  idx <- 1L
  
  # Generate combinations more efficiently
  for (k in seq_len(max_k)) {
    combs <- utils::combn(x, k, simplify = FALSE)
    n_combs <- length(combs)
    result[idx:(idx + n_combs - 1L)] <- combs
    idx <- idx + n_combs
  }
  
  # Trim to actual size
  result[seq_len(idx - 1L)]
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
