#' @title Generate combinations
#'
#' @param x Vector of elements to combine
#' @param max_ions Maximum number of ions in combinations
#'
#' @return List of all combinations
#'
#' @examples NULL
generate_combinations <- function(x, max_ions) {
  lapply(
    1:min(length(x), max_ions),
    function(k) {
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
#' @param show_progress Logical: show progress bar (default TRUE)
#'
#' @return List of all combinations for each index
#'
#' @examples NULL
generate_combinations_progress <- function(
  indices,
  ions_list,
  max_ions,
  show_progress = TRUE
) {
  n <- length(indices)

  lapply(
    seq_along(indices),
    function(idx) {
      index <- indices[idx]
      result <- generate_combinations(
        x = ions_list[[index]],
        max_ions = max_ions
      )

      if (show_progress) {
        show_progress(idx, n)
      }

      result
    }
  )
}
