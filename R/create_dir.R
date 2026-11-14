#' @title Create directory
#'
#' @description Ensures the directory at `path` exists, creating it recursively
#'   if needed.
#'
#' @param path Path to the directory to create
#'
#' @return Path to the directory (invisibly)
#'
# nolint start: line_length_linter.
#' @source https://github.com/taxonomicallyinformedannotation/tima/blob/main/R/create_dir.R
# nolint end
#' @keywords internal
#' @noRd
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Directory ", path, " created.")
  }

  invisible(path)
}
