#' @title Create directory
#'
#' @description This function creates a directory at the specified path
#'    if it does not already exist.
#'
#' @param export Path to the directory to be created
#'
#' @return Message indicating the status of directory creation
#'
#' @source https://github.com/taxonomicallyinformedannotation/tima/blob/main/R/create_dir.R
create_dir <- function(export) {
  ## Check if the export path includes a file name
  if (
    grepl(
      pattern = ".",
      x = export,
      fixed = TRUE
    )
  ) {
    dirname_path <- dirname(export)
  } else {
    dirname_path <- export
  }

  ## Create the directory at the specified path if it does not exist
  if (!dir.exists(dirname_path)) {
    dir.create(dirname_path, recursive = TRUE)
    message("Directory ", dirname_path, " created.")
  }
}
