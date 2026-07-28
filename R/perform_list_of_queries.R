#' @title Check if spectrum contains all target m/z values (internal variant with pre-computed tolerances)
#'
#' @param spec_mz Numeric vector of spectrum m/z values
#' @param target_mz Numeric vector of target m/z values to find
#' @param tol Numeric vector of pre-computed tolerances (length must match target_mz)
#'
#' @return Logical indicating if all targets are found
#'
#' @examples NULL
contains_all_mz_tol <- function(spec_mz, target_mz, tol) {
  spec_mz <- spec_mz[is.finite(spec_mz)]
  target_mz <- target_mz[is.finite(target_mz)]

  if (length(spec_mz) == 0L || length(target_mz) == 0L) {
    return(FALSE)
  }

  all(vapply(
    X = seq_along(target_mz),
    FUN = function(i) {
      any(abs(spec_mz - target_mz[i]) <= tol[i])
    },
    FUN.VALUE = logical(1)
  ))
}

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
  target_mz <- target_mz[cheapr::which_not_na(target_mz)]

  if (length(target_mz) == 0L) {
    return(FALSE)
  }

  tol <- dalton + target_mz * ppm / 1e6
  contains_all_mz_tol(spec_mz, target_mz, tol)
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
    frags <- frags[cheapr::which_not_na(frags)]
    if (length(frags) > 0L) {
      # Pre-compute tolerances for fragments to avoid recalculating for each spectrum
      frags_tol <- dalton + frags * ppm / 1e6

      keep <- keep &
        vapply(
          X = all_mz,
          FUN = contains_all_mz_tol,
          FUN.VALUE = logical(1),
          target_mz = frags,
          tol = frags_tol
        )
    }
  }

  # --- Filter by neutral losses (per spectrum) ---
  if (length(nls) > 0L) {
    precursor_mz <- Spectra::precursorMz(spectra)
    valid <- cheapr::which_not_na(precursor_mz)

    keep <- keep &
      vapply(
        X = seq_along(spectra),
        FUN = function(i) {
          if (!(i %in% valid)) {
            return(FALSE)
          }
          spec <- all_mz[[i]]
          if (length(spec) == 0L || cheapr::all_na(spec)) {
            return(FALSE)
          }
          target_mz <- precursor_mz[i] - nls
          # Filter for finite and positive values
          valid_idx <- cheapr::which_not_na(target_mz)
          if (length(valid_idx) == 0L) {
            return(FALSE)
          }
          target_mz <- target_mz[valid_idx]
          target_mz <- target_mz[target_mz > 0]
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

  # Extract fragments and neutral losses using stringi (fixed-pattern detection)
  is_frag <- stringi::stri_detect_fixed(ions, "_frag")
  is_nl <- stringi::stri_detect_fixed(ions, "_nl")

  frags <- if (any(is_frag)) {
    vals <- suppressWarnings(as.numeric(stringi::stri_replace_all_fixed(
      ions[is_frag],
      "_frag",
      ""
    )))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  nls <- if (any(is_nl)) {
    vals <- suppressWarnings(as.numeric(stringi::stri_replace_all_fixed(
      ions[is_nl],
      "_nl",
      ""
    )))
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
      if (
        nrow(spectra_data) > 0L &&
          fastmatch::fmatch("SKELETON", names(spectra_data), nomatch = 0L) > 0L
      ) {
        skeleton <- spectra_data$SKELETON
        if (length(skeleton) > 0L && !cheapr::all_na(skeleton)) {
          stringi::stri_replace_all_fixed(as.character(skeleton), "+", ".")
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

  # Process in batches to allow periodic garbage collection
  batch_size <- 100L
  n_batches <- ceiling(n / batch_size)

  for (batch in seq_len(n_batches)) {
    start_idx <- (batch - 1L) * batch_size + 1L
    end_idx <- min(batch * batch_size, n)

    for (i in start_idx:end_idx) {
      results[[i]] <- tryCatch(
        {
          perform_list_of_queries(i, ions_list, spectra, dalton, ppm)
        },
        error = function(e) {
          warning("Error in query ", i, ": ", e$message)
          tidytable::tidytable(
            target = names(ions_list)[i] %||% paste0("query_", i),
            value = character(0)
          )
        }
      )
      if (!is.null(pb)) pb$tick()
    }

    # Garbage collect after each batch
    if (batch < n_batches) {
      invisible(gc(verbose = FALSE))
    }
  }

  results
}
