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
  # Early return if no spectra to search
  if (length(spectra) == 0L) {
    return(spectra)
  }

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

    # Normalize to vector
    if (is.matrix(idx)) {
      idx <- as.vector(idx)
    }

    # Validate and subset
    if (length(idx) != length(spectra)) {
      return(spectra[FALSE])
    }

    spectra <- spectra[idx]

    # Early return if empty
    if (length(spectra) == 0L) {
      return(spectra)
    }
  }

  # Filter by neutral losses if specified
  if (length(nls) > 0L) {
    # Calculate neutral loss m/z values
    precursor_mz <- Spectra::precursorMz(spectra)
    mzs <- sort(unique(precursor_mz - nls))

    idx <- Spectra::containsMz(
      spectra,
      mz = mzs,
      tolerance = dalton,
      ppm = ppm,
      which = "all",
      BPPARAM = BiocParallel::SerialParam()
    )

    # Normalize to vector
    if (is.matrix(idx)) {
      idx <- as.vector(idx)
    }

    # Validate and subset
    if (length(idx) != length(spectra)) {
      return(spectra[FALSE])
    }

    spectra <- spectra[idx]
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
  # Get and validate target name
  target <- names(ions_list)[index]
  if (is.null(target) || is.na(target) || target == "") {
    target <- paste0("query_", index)
  }

  ions <- ions_list[[index]]

  # Extract fragments and neutral losses
  is_frag <- grepl("_frag$", ions)
  is_nl <- grepl("_nl$", ions)

  frags <- if (any(is_frag)) {
    as.numeric(sub("_frag$", "", ions[is_frag]))
  } else {
    numeric(0)
  }
  nls <- if (any(is_nl)) {
    as.numeric(sub("_nl$", "", ions[is_nl]))
  } else {
    numeric(0)
  }

  # Perform query
  result_spectra <- perform_query(
    spectra = spectra,
    frags = frags,
    nls = nls,
    dalton = dalton,
    ppm = ppm
  )

  # Early return if no results
  if (length(result_spectra) == 0L) {
    return(tidytable::tidytable(target = target, value = character(0)))
  }

  # Extract skeleton values
  spectra_data <- result_spectra |>
    Spectra::spectraData()

  value <- if (nrow(spectra_data) > 0L && "SKELETON" %in% names(spectra_data)) {
    skeleton <- spectra_data$SKELETON
    if (length(skeleton) > 0L) {
      gsub("+", ".", skeleton, fixed = TRUE)
    } else {
      character(0)
    }
  } else {
    character(0)
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
  purrr::map(
    .progress = interactive(),
    .x = seq_along(ions_list),
    .f = function(index) {
      perform_list_of_queries(index, ions_list, spectra, dalton, ppm)
    }
  )
}
