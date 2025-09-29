#' @title Perform list of queries
#'
#' @include perform_query.R
#'
#' @param index Index of the ion list to process
#' @param ions_list List of ions for queries
#' @param spectra Spectra object to search
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Data frame with target and value columns
#'
#' @examples NULL
perform_list_of_queries <- function(index, ions_list, spectra, dalton, ppm) {
  target <- names(ions_list)[index]
  ions <- ions_list[[index]]
  
  # Vectorized extraction of fragments and neutral losses
  is_frag <- grepl("_frag$", ions)
  is_nl <- grepl("_nl$", ions)
  
  frags <- as.numeric(sub("_frag$", "", ions[is_frag]))
  nls <- as.numeric(sub("_nl$", "", ions[is_nl]))
  
  # Perform the query
  result_spectra <- spectra |>
    perform_query(frags = frags, nls = nls, dalton = dalton, ppm = ppm)
  
  # Extract skeleton information more efficiently
  if (length(result_spectra) > 0) {
    value <- Spectra::spectraData(result_spectra)$SKELETON
    value <- gsub("+", ".", value, fixed = TRUE)
  } else {
    value <- character(0)
  }
  
  return(tidytable::tidytable(target = target, value = value))
}

#' @title Perform list of queries (progress)
#'
#' @param ions_list Ions list
#' @param spectra Spectra
#' @param dalton Dalton
#' @param ppm PPM
#'
#' @return NULL
#'
#' @examples NULL
perform_list_of_queries_progress <- function(ions_list, spectra, dalton, ppm) {
  purrr::map(
    .progress = TRUE,
    .x = seq_along(ions_list),
    .f = function(index, ions_list, spectra, dalton, ppm) {
      perform_list_of_queries(index, ions_list, spectra, dalton, ppm)
    },
    ions_list = ions_list,
    spectra = spectra,
    dalton = dalton,
    ppm = ppm
  )
}
