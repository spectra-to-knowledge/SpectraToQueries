#' @title Normalize peaks
#'
#' @return A function that normalizes a spectrum peak matrix to a maximum intensity of 100
#'
#' @examples NULL
normalize_peaks <- function() {
  function(x, ...) {
    if (nrow(x) > 0L) {
      max_int <- max(x[, 2L])
      if (max_int > 0) {
        x[, 2L] <- 100 * x[, 2L] / max_int
      }
    }
    x
  }
}
