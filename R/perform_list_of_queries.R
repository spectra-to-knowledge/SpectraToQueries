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

  # Extract skeleton information
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
#' @param ions_list List of ion combinations for queries
#' @param spectra Spectra object to search
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return List of query results
#'
#' @examples NULL
perform_list_of_queries_progress <- function(ions_list, spectra, dalton, ppm) {
  # Determine if parallel processing would be beneficial
  n_queries <- length(ions_list)
  use_parallel <- n_queries > 1000L # Use parallel for more than 10 queries

  if (use_parallel) {
    # Use parallel processing if beneficial
    message("Using parallel processing for ", n_queries, " queries.")
    BiocParallel::bplapply(
      X = seq_along(ions_list),
      FUN = function(index) {
        perform_list_of_queries(index, ions_list, spectra, dalton, ppm)
      },
      BPPARAM = BiocParallel::bpparam(),
      BPOPTIONS = BiocParallel::bpoptions(progressbar = TRUE)
    )
  } else {
    # Fall back to sequential processing with progress
    purrr::map(
      .progress = TRUE,
      .x = seq_along(ions_list),
      .f = function(index) {
        perform_list_of_queries(index, ions_list, spectra, dalton, ppm)
      }
    )
  }
}
