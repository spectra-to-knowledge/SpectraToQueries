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

#' @title Check if spectrum contains all target m/z values (internal variant with pre-computed tolerances)
#'
#' @param spec_mz Numeric vector of spectrum m/z values
#' @param target_mz Numeric vector of target m/z values to find
#' @param tol Numeric vector of pre-computed tolerances (length must match `target_mz`)
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
        target_mz <- target_mz[is.finite(target_mz) & target_mz > 0]
        keep[i] <- length(target_mz) > 0L &&
          contains_all_mz(spec, target_mz, dalton, ppm)
      }
      keep[!valid_mask] <- FALSE
    }
  }

  keep
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
perform_list_of_queries_progress <- function(
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

  n_spec <- length(all_mz)
  valid_prec_and_spec <- is.finite(precursor_mz) & (lengths(all_mz) > 0L)

  frag_hits <- list()
  nl_hits <- list()
  nl_vals <- list()

  for (u in unique_ions) {
    if (grepl("_frag", u, fixed = TRUE)) {
      val <- suppressWarnings(as.numeric(sub("_frag", "", u, fixed = TRUE)))
      if (is.finite(val)) {
        frag_hits[[u]] <- perform_query(
          all_mz = all_mz,
          precursor_mz = precursor_mz,
          frags = val,
          nls = numeric(0),
          dalton = dalton,
          ppm = ppm
        )
      }
    } else if (grepl("_nl", u, fixed = TRUE)) {
      val <- suppressWarnings(as.numeric(sub("_nl", "", u, fixed = TRUE)))
      if (is.finite(val)) {
        nl_vals[[u]] <- val
        hit_vec <- logical(n_spec)
        for (idx in which(valid_prec_and_spec)) {
          tmz <- precursor_mz[idx] - val
          if (tmz <= 0) {
            # In multi-NL queries, target_mz filters to > 0; non-emptiness of
            # target_mz is enforced below via precursor_mz > min(nls).
            hit_vec[idx] <- TRUE
          } else {
            hit_vec[idx] <- contains_all_mz(
              all_mz[[idx]],
              tmz,
              dalton = dalton,
              ppm = ppm
            )
          }
        }
        nl_hits[[u]] <- hit_vec
      }
    }
  }

  empty_df <- data.frame(
    target = character(0),
    value = character(0),
    stringsAsFactors = FALSE
  )
  target_names <- names(ions_list)
  frag_names_set <- names(frag_hits)
  nl_names_set <- names(nl_hits)

  for (i in seq_len(n)) {
    target <- if (!is.null(target_names)) target_names[i] else NULL
    if (is.null(target) || is.na(target) || target == "") {
      target <- paste0("query_", i)
    }

    ions <- ions_list[[i]]
    q_frags <- ions[ions %in% frag_names_set]
    q_nls <- ions[ions %in% nl_names_set]

    if (length(q_frags) == 0L && length(q_nls) == 0L) {
      results[[i]] <- empty_df
    } else {
      keep <- NULL
      for (f_ion in q_frags) {
        keep <- if (is.null(keep)) {
          frag_hits[[f_ion]]
        } else {
          keep & frag_hits[[f_ion]]
        }
      }
      if (length(q_nls) > 0L) {
        min_nl <- min(unlist(nl_vals[q_nls], use.names = FALSE))
        nl_nonempty <- valid_prec_and_spec & (precursor_mz > min_nl)
        keep <- if (is.null(keep)) nl_nonempty else keep & nl_nonempty
        for (n_ion in q_nls) {
          keep <- keep & nl_hits[[n_ion]]
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
