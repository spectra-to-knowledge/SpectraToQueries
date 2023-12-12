#' Title
#'
#' @return
#' @export
#'
#' @examples
normalize_peaks <- function() {
  function(x, ...) {
    x[, 2] <- 100 * x[, 2] / max(x[, 2])
    return(x)
  }
}
