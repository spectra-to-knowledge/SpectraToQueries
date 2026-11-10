#' @title Prepare spectra data for querying
#' @param spectra Spectra object
#' @description Extracts mz, precursor m/z, and skeleton labels once as plain R
#'   vectors/lists so `perform_query()` never has to touch S4 objects again.
#' @examples NULL
prepare_query_data <- function(spectra) {
  raw_mz <- as.list(Spectra::mz(spectra))
  clean_mz <- lapply(raw_mz, function(x) {
    if (all(is.finite(x))) x else x[is.finite(x)]
  })
  list(
    mz = clean_mz,
    precursor_mz = as.numeric(Spectra::precursorMz(spectra)),
    skeleton = gsub(
      "+",
      ".",
      as.character(Spectra::spectraData(spectra)$SKELETON),
      fixed = TRUE
    )
  )
}

#' @title Check if spectrum contains all target m/z values (pre-computed
#'   tolerances)
#'
#' @param spec_mz Numeric vector of spectrum m/z values
#' @param target_mz Numeric vector of target m/z values to find
#' @param tol Numeric vector of pre-computed tolerances (same length as
#'   `target_mz`)
#'
#' @return Logical indicating if all targets are found
#'
#' @examples NULL
contains_all_mz_tol <- function(spec_mz, target_mz, tol) {
  if (!all(is.finite(spec_mz))) {
    spec_mz <- spec_mz[is.finite(spec_mz)]
  }

  if (length(spec_mz) == 0L || length(target_mz) == 0L) {
    return(FALSE)
  }

  for (i in seq_along(target_mz)) {
    if (!any(abs(spec_mz - target_mz[i]) <= tol[i])) {
      return(FALSE)
    }
  }
  TRUE
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
  target_mz <- target_mz[is.finite(target_mz)]

  if (length(target_mz) == 0L) {
    return(FALSE)
  }

  tol <- dalton + target_mz * (ppm / 1e6)
  contains_all_mz_tol(spec_mz, target_mz, tol)
}

#' @title Perform query
#'
#' @param all_mz List of numeric vectors of m/z values for each spectrum
#' @param precursor_mz Numeric vector of precursor m/z values for each spectrum
#' @param frags Fragment masses to search for
#' @param nls Neutral losses to search for
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Logical vector indicating matching spectra
#'
#' @examples NULL
perform_query <- function(all_mz, precursor_mz, frags, nls, dalton, ppm) {
  n <- length(all_mz)
  if (n == 0L) {
    return(logical(0))
  }

  if (isS4(all_mz)) {
    all_mz <- as.list(all_mz)
  }

  keep <- rep(TRUE, n)

  if (length(frags) > 0L) {
    frags <- frags[is.finite(frags)]
    if (length(frags) > 0L) {
      frags_tol <- dalton + frags * (ppm / 1e6)
      keep <- vapply(
        all_mz,
        contains_all_mz_tol,
        FUN.VALUE = logical(1L),
        target_mz = frags,
        tol = frags_tol
      )
    }
  }

  if (length(nls) > 0L && any(keep)) {
    nls <- nls[is.finite(nls)]
    if (length(nls) > 0L) {
      valid_mask <- is.finite(precursor_mz)
      idx <- which(keep & valid_mask)
      for (i in idx) {
        spec <- all_mz[[i]]
        if (length(spec) == 0L || !any(is.finite(spec))) {
          keep[i] <- FALSE
          next
        }
        target_mz <- precursor_mz[i] - nls
        keep[i] <- all(target_mz > 0) &&
          contains_all_mz(spec, target_mz, dalton, ppm)
      }
      keep[!valid_mask] <- FALSE
    }
  }

  keep
}

#' @title Parse fragment and neutral-loss ions from query character vector
#' @keywords internal
#' @noRd
parse_query_ions <- function(ions) {
  is_frag <- grepl("_frag", ions, fixed = TRUE)
  is_nl <- grepl("_nl", ions, fixed = TRUE)

  frags <- if (any(is_frag)) {
    vals <- suppressWarnings(as.numeric(
      sub("_frag", "", ions[is_frag], fixed = TRUE)
    ))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  nls <- if (any(is_nl)) {
    vals <- suppressWarnings(as.numeric(
      sub("_nl", "", ions[is_nl], fixed = TRUE)
    ))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  list(frags = frags, nls = nls)
}

#' @title Perform list of queries
#'
#' @param index Index of the ion list to process
#' @param ions_list List of ions for queries
#' @param query_data List with elements `mz`, `precursor_mz`, and `skeleton`
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Data frame with `target` and `value` columns
#'
#' @examples NULL
perform_list_of_queries <- function(index, ions_list, query_data, dalton, ppm) {
  target <- names(ions_list)[index]
  if (is.null(target) || is.na(target) || target == "") {
    target <- paste0("query_", index)
  }

  empty_df <- data.frame(
    target = character(0),
    value = character(0),
    stringsAsFactors = FALSE
  )

  ions <- ions_list[[index]]
  if (length(ions) == 0L || all(is.na(ions))) {
    return(empty_df)
  }

  parsed <- parse_query_ions(ions)
  frags <- parsed$frags
  nls <- parsed$nls

  if (length(frags) == 0L && length(nls) == 0L) {
    return(empty_df)
  }

  keep <- tryCatch(
    perform_query(
      all_mz = query_data$mz,
      precursor_mz = query_data$precursor_mz,
      frags = frags,
      nls = nls,
      dalton = dalton,
      ppm = ppm
    ),
    error = function(e) {
      warning("Error in perform_query for ", target, ": ", e$message)
      rep(FALSE, length(query_data$mz))
    }
  )

  if (!any(keep)) {
    return(empty_df)
  }

  value <- query_data$skeleton[keep]
  value <- value[!is.na(value)]

  if (length(value) == 0L) {
    return(empty_df)
  }

  data.frame(
    target = rep.int(target, length(value)),
    value = value,
    stringsAsFactors = FALSE
  )
}

#' @title Precompute single-ion logical hit vectors across all spectra
#' @keywords internal
#' @noRd
precompute_ion_hits <- function(
  unique_ions,
  all_mz,
  precursor_mz,
  valid_prec_and_spec,
  dalton,
  ppm
) {
  n_spec <- length(all_mz)
  valid_idx <- which(valid_prec_and_spec)
  hits <- list()

  for (u in unique_ions) {
    if (grepl("_frag", u, fixed = TRUE)) {
      val <- suppressWarnings(as.numeric(sub("_frag", "", u, fixed = TRUE)))
      if (is.finite(val)) {
        tol <- dalton + val * (ppm / 1e6)
        hits[[u]] <- vapply(
          all_mz,
          function(spec_mz) any(abs(spec_mz - val) <= tol),
          logical(1L)
        )
      }
    } else if (grepl("_nl", u, fixed = TRUE)) {
      val <- suppressWarnings(as.numeric(sub("_nl", "", u, fixed = TRUE)))
      if (is.finite(val)) {
        vec <- logical(n_spec)
        for (idx in valid_idx) {
          tmz <- precursor_mz[idx] - val
          if (tmz > 0) {
            tol <- dalton + tmz * (ppm / 1e6)
            vec[idx] <- any(abs(all_mz[[idx]] - tmz) <= tol)
          }
        }
        hits[[u]] <- vec
      }
    }
  }

  list(
    names = names(hits),
    hits = hits
  )
}

#' @title Perform list of queries with progress
#'
#' @param ions_list List of ion combinations for queries
#' @param spectra Spectra object to search
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#' @param show_progress Logical: show progress bar (default `TRUE`)
#'
#' @return List of query result data frames
#'
#' @examples NULL
# nolint start: object_length_linter.
perform_list_of_queries_progress <- function(
  # nolint end
  ions_list,
  spectra,
  dalton,
  ppm,
  show_progress = TRUE
) {
  n <- length(ions_list)
  results <- vector("list", n)
  if (n == 0L) {
    return(results)
  }

  query_data <- prepare_query_data(spectra)
  all_mz <- query_data$mz
  precursor_mz <- query_data$precursor_mz
  skeleton <- query_data$skeleton

  # Precompute single-ion logical hit vectors across all spectra once
  unique_ions <- unique(unlist(ions_list, use.names = FALSE))
  unique_ions <- unique_ions[!is.na(unique_ions) & nzchar(unique_ions)]

  valid_prec_and_spec <- is.finite(precursor_mz) & (lengths(all_mz) > 0L)
  precomp <- precompute_ion_hits(
    unique_ions = unique_ions,
    all_mz = all_mz,
    precursor_mz = precursor_mz,
    valid_prec_and_spec = valid_prec_and_spec,
    dalton = dalton,
    ppm = ppm
  )
  all_valid_names <- precomp$names
  hits_lookup <- precomp$hits

  empty_df <- data.frame(
    target = character(0),
    value = character(0),
    stringsAsFactors = FALSE
  )
  target_names <- names(ions_list)

  for (i in seq_len(n)) {
    target <- if (!is.null(target_names)) target_names[i] else NULL
    if (is.null(target) || is.na(target) || target == "") {
      target <- paste0("query_", i)
    }

    idx_vec <- match(ions_list[[i]], all_valid_names, nomatch = 0L)
    idx_vec <- idx_vec[idx_vec > 0L]

    if (length(idx_vec) == 0L) {
      results[[i]] <- empty_df
    } else {
      keep <- hits_lookup[[idx_vec[1L]]]
      for (code in idx_vec[-1L]) {
        keep <- keep & hits_lookup[[code]]
        if (!any(keep)) {
          break
        }
      }

      if (!any(keep)) {
        results[[i]] <- empty_df
      } else {
        val <- skeleton[keep]
        val <- val[!is.na(val)]
        n_val <- length(val)
        if (n_val == 0L) {
          results[[i]] <- empty_df
        } else {
          results[[i]] <- structure(
            list(
              target = rep.int(target, n_val),
              value = val
            ),
            class = "data.frame",
            row.names = c(NA, -n_val)
          )
        }
      }
    }

    if (show_progress) {
      show_progress(i, n)
    }
  }

  results
}
