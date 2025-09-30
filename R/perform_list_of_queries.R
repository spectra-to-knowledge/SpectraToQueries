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
  if (length(spec_mz) == 0L || length(target_mz) == 0L) {
    return(FALSE)
  }
  abs_tol <- dalton + (target_mz * ppm / 1e6)
  # Vectorized outer subtraction
  hits <- outer(spec_mz, target_mz, function(x, y) abs(x - y) <= abs_tol)
  all(colSums(hits) > 0)
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
  # Early return if no spectra to search
  if (length(spectra) == 0L) {
    return(spectra)
  }

  # --- Filter by fragments ---
  if (length(frags) > 0L) {
    # Remove any NA or infinite values from frags
    frags <- frags[is.finite(frags)]

    if (length(frags) > 0L) {
      # Extract all m/z values at once
      idx <- tryCatch(
        {
          all_mz <- Spectra::mz(spectra)

          # Ensure we got the right number of m/z lists
          if (length(all_mz) != length(spectra)) {
            warning("m/z list length mismatch")
            return(rep(FALSE, length(spectra)))
          }

          # Vectorized check across all spectra
          vapply(
            all_mz,
            function(spec_mz) {
              contains_all_mz(spec_mz, frags, dalton, ppm)
            },
            logical(1),
            USE.NAMES = FALSE
          )
        },
        error = function(e) {
          warning("Error in fragment filtering: ", e$message)
          return(rep(FALSE, length(spectra)))
        }
      )

      # Ensure logical vector and correct length
      if (length(idx) != length(spectra)) {
        warning("Index length mismatch after fragment filtering")
        return(spectra[FALSE])
      }

      idx[is.na(idx)] <- FALSE

      # Subset safely
      if (any(idx)) {
        spectra <- spectra[idx]
      } else {
        return(spectra[FALSE])
      }

      # Early return if empty
      if (length(spectra) == 0L) {
        return(spectra)
      }
    }
  }

  # --- Filter by neutral losses ---
  if (length(nls) > 0L) {
    # Remove any NA or infinite values from nls
    nls <- nls[is.finite(nls)]

    if (length(nls) > 0L) {
      precursor_mz <- tryCatch(
        {
          Spectra::precursorMz(spectra)
        },
        error = function(e) {
          warning("Error getting precursor m/z: ", e$message)
          return(rep(NA_real_, length(spectra)))
        }
      )

      # Ensure correct length
      if (length(precursor_mz) != length(spectra)) {
        warning("Precursor m/z length mismatch")
        return(spectra[FALSE])
      }

      # Handle NA precursor m/z
      valid_idx <- !is.na(precursor_mz) & is.finite(precursor_mz)

      if (!any(valid_idx)) {
        return(spectra[FALSE])
      }

      # Extract all m/z values at once
      idx <- tryCatch(
        {
          all_mz <- Spectra::mz(spectra)

          # Ensure we got the right number of m/z lists
          if (length(all_mz) != length(spectra)) {
            warning("m/z list length mismatch in neutral loss filtering")
            return(rep(FALSE, length(spectra)))
          }

          # Vectorized check across all spectra
          result <- vapply(
            seq_along(spectra),
            function(i) {
              if (!valid_idx[i]) {
                return(FALSE)
              }

              spec_mz <- all_mz[[i]]
              if (length(spec_mz) == 0L) {
                return(FALSE)
              }

              # Calculate target m/z values for this spectrum
              target_mz <- precursor_mz[i] - nls
              target_mz <- target_mz[is.finite(target_mz) & target_mz > 0]

              if (length(target_mz) == 0L) {
                return(FALSE)
              }

              contains_all_mz(spec_mz, target_mz, dalton, ppm)
            },
            logical(1),
            USE.NAMES = FALSE
          )

          result
        },
        error = function(e) {
          warning("Error in neutral loss filtering: ", e$message)
          return(rep(FALSE, length(spectra)))
        }
      )

      # Ensure logical vector and correct length
      if (length(idx) != length(spectra)) {
        warning("Index length mismatch after neutral loss filtering")
        return(spectra[FALSE])
      }

      idx[is.na(idx)] <- FALSE

      # Subset safely
      if (any(idx)) {
        spectra <- spectra[idx]
      } else {
        return(spectra[FALSE])
      }
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
  # Get and validate target name
  target <- names(ions_list)[index]
  if (is.null(target) || is.na(target) || target == "") {
    target <- paste0("query_", index)
  }

  ions <- ions_list[[index]]

  # Validate ions
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

  # Early return if no results
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
  # Force garbage collection before starting
  invisible(gc(verbose = FALSE, full = TRUE))

  results <- vector("list", length(ions_list))

  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      total = length(ions_list),
      format = "[:bar] :current/:total (:percent) eta: :eta"
    )
  } else {
    pb <- NULL
  }

  for (i in seq_along(ions_list)) {
    results[[i]] <- tryCatch(
      {
        perform_list_of_queries(i, ions_list, spectra, dalton, ppm)
      },
      error = function(e) {
        warning("Critical error in query ", i, ": ", e$message)
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
