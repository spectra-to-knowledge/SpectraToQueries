#' @title Simple progress indicator (base R only)
#'
#' @param current Current iteration number
#' @param total Total iterations
#'
#' @return `NULL` (invisibly). Prints progress message at 10% milestones.
#'
#' @keywords internal
#' @examples NULL
show_progress <- function(current, total) {
  if (total <= 0L) {
    return(invisible(NULL))
  }

  if (current == total) {
    message(sprintf("Progress: %d/%d (100%%) - Done!", current, total))
    return(invisible(NULL))
  }

  decile <- (10L * current) %/% total
  prev_decile <- (10L * (current - 1L)) %/% total

  if (decile > prev_decile && decile < 10L) {
    pct <- round(100 * current / total)
    message(sprintf("Progress: %d/%d (%d%%)", current, total, pct))
  }

  invisible(NULL)
}
