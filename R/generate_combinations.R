#' @title Generate combinations
#'
#' @param x Vector of elements to combine
#' @param max_ions Maximum number of ions in combinations
#'
#' @return List of all combinations
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


#' @title Generate combinations with progress
#'
#' @param indices Vector of indices to process
#' @param ions_list List of ions for each index
#' @param max_ions Maximum number of ions in combinations
#'
#' @return List of all combinations for each index
#'
#' @examples NULL
generate_combinations_progress <- function(indices, ions_list, max_ions) {
  n_indices <- length(indices)
  use_parallel <- n_indices > 10L # Use parallel for larger sets

  if (use_parallel) {
    # Use parallel processing if beneficial
    message(
      "Using parallel processing for combination generation with ",
      n_indices,
      " groups."
    )
    BiocParallel::bplapply(
      X = indices,
      FUN = function(index) {
        x <- ions_list[[index]]
        generate_combinations(x = x, max_ions = max_ions)
      },
      BPPARAM = BiocParallel::bpparam(),
      BPOPTIONS = BiocParallel::bpoptions(progressbar = TRUE)
    )
  } else {
    # Fall back to sequential processing with progress
    purrr::map(
      .progress = TRUE,
      .x = indices,
      .f = function(index) {
        x <- ions_list[[index]]
        generate_combinations(x = x, max_ions = max_ions)
      }
    )
  }
}
