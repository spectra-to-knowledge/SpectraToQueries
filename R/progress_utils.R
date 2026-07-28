#' @title Simple progress indicator (base R only)
#'
#' @param current Current iteration number
#' @param total Total iterations
#'
#' @return NULL (invisibly). Prints progress message at milestones.
#'
#' @keywords internal
#' @examples NULL
show_progress <- function(current, total) {
  # Show progress at 10% intervals
  pct <- round(100 * current / total)
  prev_pct <- round(100 * (current - 1) / total)

  # Print message when crossing 10% threshold
  if (pct > prev_pct && pct %% 10 == 0) {
    message(sprintf("Progress: %d/%d (%d%%)", current, total, pct))
  }

  # Print final message
  if (current == total) {
    message(sprintf("Progress: %d/%d (100%%) - Done!", current, total))
  }

  invisible(NULL)
}
