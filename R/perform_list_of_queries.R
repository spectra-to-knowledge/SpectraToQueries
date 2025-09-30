#' @title Perform query
#'
#' @param spectra Spectra object
#' @param frags Fragment masses to search for
#' @param nls Neutral losses to search for
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Filtered spectra object
#'
#' @examples NULL
perform_query <- function(spectra, frags, nls, dalton, ppm) {
  # Filter by fragments if specified
  if (length(frags) > 0L) {
    idx <- Spectra::containsMz(
      spectra,
      mz = frags,
      tolerance = dalton,
      ppm = ppm,
      which = "all",
      BPPARAM = BiocParallel::SerialParam()
    )

    if (is.matrix(idx)) {
      idx <- as.vector(idx)
    }

    # Subset if valid index
    if (length(idx) == length(spectra)) {
      spectra <- spectra[idx]
    } else {
      return(spectra[FALSE]) # empty
    }

    if (length(spectra) == 0L) {
      return(spectra)
    }
  }

  # Filter by neutral losses if specified
  if (length(nls) > 0L) {
    mzs <- sort(unique(Spectra::precursorMz(spectra) - nls))

    idx <- Spectra::containsMz(
      spectra,
      mz = mzs,
      tolerance = dalton,
      ppm = ppm,
      which = "all",
      BPPARAM = BiocParallel::SerialParam()
    )

    if (is.matrix(idx)) {
      idx <- as.vector(idx)
    }

    if (length(idx) == length(spectra)) {
      spectra <- spectra[idx]
    } else {
      return(spectra[FALSE])
    }
  }

  spectra
}

#' @title Perform list of queries
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
  if (is.null(target) || is.na(target) || target == "") {
    target <- paste0("query_", index)
  }

  ions <- ions_list[[index]]

  # Extract frags and nls
  is_frag <- grepl("_frag$", ions)
  is_nl <- grepl("_nl$", ions)

  frags <- as.numeric(sub("_frag$", "", ions[is_frag]))
  nls <- as.numeric(sub("_nl$", "", ions[is_nl]))

  # Perform query
  result_spectra <- perform_query(
    spectra = spectra,
    frags = frags,
    nls = nls,
    dalton = dalton,
    ppm = ppm
  )

  # Early return if no spectra remain
  if (length(result_spectra) == 0L) {
    return(tidytable::tidytable(target = target, value = character(0)))
  }

  spectra_data <- result_spectra |>
    Spectra::spectraData()

  value <- character(0)
  if (nrow(spectra_data) > 0L && "SKELETON" %in% names(spectra_data)) {
    value <- spectra_data$SKELETON
    if (length(value) > 0L) {
      value <- gsub("+", ".", value, fixed = TRUE)
    }
  }

  tidytable::tidytable(target = target, value = value)
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
  n_queries <- length(ions_list)
  message("Processing ", n_queries, " queries...")

  results <- vector("list", n_queries)

  for (i in seq_len(n_queries)) {
    if (
      i == 1L || i %% max(100L, floor(n_queries / 10)) == 0L || i == n_queries
    ) {
      message(sprintf(
        "Progress: %d/%d (%.1f%%)",
        i,
        n_queries,
        100 * i / n_queries
      ))
    }

    results[[i]] <- perform_list_of_queries(i, ions_list, spectra, dalton, ppm)
  }

  message("Processing complete!")
  results
}
