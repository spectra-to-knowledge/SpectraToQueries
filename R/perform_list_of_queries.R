#' @title Check if spectrum contains all target m/z values
#'
#' @param spec_mz Numeric vector of spectrum m/z values
#' @param target_mz Numeric vector of target m/z values to find
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Logical indicating if all targets are found
#'
#' @examples NULL
contains_all_mz <- function(spec_mz, target_mz, dalton, ppm) {
  spec_mz <- spec_mz[is.finite(spec_mz)]
  target_mz <- target_mz[is.finite(target_mz)]

  if (length(spec_mz) == 0L || length(target_mz) == 0L) {
    return(FALSE)
  }

  tol <- dalton + target_mz * ppm / 1e6

  all(vapply(
    X = seq_along(target_mz),
    FUN = function(i) {
      any(abs(spec_mz - target_mz[i]) <= tol[i])
    },
    FUN.VALUE = logical(1)
  ))
}

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
  if (length(spectra) == 0L) {
    return(spectra)
  }

  all_mz <- Spectra::mz(spectra)
  keep <- rep(TRUE, length(spectra))

  # --- Filter by fragments ---
  if (length(frags) > 0L) {
    frags <- frags[is.finite(frags)]
    if (length(frags) > 0L) {
      keep <- keep &
        vapply(
          X = all_mz,
          FUN = contains_all_mz,
          FUN.VALUE = logical(1),
          target_mz = frags,
          dalton = dalton,
          ppm = ppm
        )
    }
  }

  # --- Filter by neutral losses (per spectrum) ---
  if (length(nls) > 0L) {
    precursor_mz <- Spectra::precursorMz(spectra)
    valid <- is.finite(precursor_mz)

    keep <- keep &
      vapply(
        X = seq_along(spectra),
        FUN = function(i) {
          if (!valid[i]) {
            return(FALSE)
          }
          spec <- all_mz[[i]]
          if (length(spec) == 0L || all(is.na(spec))) {
            return(FALSE)
          }
          target_mz <- precursor_mz[i] - nls
          target_mz <- target_mz[is.finite(target_mz) & target_mz > 0]
          if (length(target_mz) == 0L) {
            return(FALSE)
          }
          contains_all_mz(spec, target_mz, dalton, ppm)
        },
        FUN.VALUE = logical(1)
      )
  }

  spectra[keep]
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
  if (length(ions) == 0L || all(is.na(ions))) {
    return(tidytable::tidytable(target = target, value = character(0)))
  }

  # Extract fragments and neutral losses
  is_frag <- grepl("_frag$", ions)
  is_nl <- grepl("_nl$", ions)

  frags <- if (any(is_frag)) {
    vals <- suppressWarnings(as.numeric(sub("_frag$", "", ions[is_frag])))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  nls <- if (any(is_nl)) {
    vals <- suppressWarnings(as.numeric(sub("_nl$", "", ions[is_nl])))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  # Early return if no valid ions
  if (length(frags) == 0L && length(nls) == 0L) {
    return(tidytable::tidytable(target = target, value = character(0)))
  }

  # Perform query with error handling
  result_spectra <- tryCatch(
    {
      perform_query(
        spectra = spectra,
        frags = frags,
        nls = nls,
        dalton = dalton,
        ppm = ppm
      )
    },
    error = function(e) {
      warning("Error in perform_query for ", target, ": ", e$message)
      return(spectra[FALSE])
    }
  )

  if (length(result_spectra) == 0L) {
    return(tidytable::tidytable(target = target, value = character(0)))
  }

  # Extract skeleton values with error handling
  value <- tryCatch(
    {
      spectra_data <- Spectra::spectraData(result_spectra)
      if (nrow(spectra_data) > 0L && "SKELETON" %in% names(spectra_data)) {
        skeleton <- spectra_data$SKELETON
        if (length(skeleton) > 0L && !all(is.na(skeleton))) {
          gsub("+", ".", as.character(skeleton), fixed = TRUE)
        } else {
          character(0)
        }
      } else {
        character(0)
      }
    },
    error = function(e) {
      warning("Error extracting SKELETON for ", target, ": ", e$message)
      character(0)
    }
  )

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
  invisible(gc(verbose = FALSE, full = TRUE))
  n <- length(ions_list)
  results <- vector("list", n)

  pb <- if (requireNamespace("progress", quietly = TRUE)) {
    progress::progress_bar$new(
      total = n,
      format = "[:bar] :current/:total (:percent) eta: :eta"
    )
  } else {
    NULL
  }

  for (i in seq_len(n)) {
    results[[i]] <- tryCatch(
      {
        perform_list_of_queries(i, ions_list, spectra, dalton, ppm)
      },
      error = function(e) {
        warning("Error in query ", i, ": ", e$message)
        tidytable::tidytable(
          target = if (!is.null(names(ions_list)[i])) {
            names(ions_list)[i]
          } else {
            paste0("query_", i)
          },
          value = character(0)
        )
      }
    )
    if (!is.null(pb)) pb$tick()
  }
  results
}
