#' @title Prepare spectra data for querying
#' @param spectra Spectra object
#' @description Extracts mz, precursor m/z, and skeleton labels once, so
#'   perform_query() never has to touch the Spectra backend again.
#' @examples NULL
prepare_query_data <- function(spectra) {
  list(
    mz = Spectra::mz(spectra),
    precursor_mz = Spectra::precursorMz(spectra),
    skeleton = stringi::stri_replace_all_fixed(
      as.character(Spectra::spectraData(spectra)$SKELETON),
      "+",
      "."
    )
  )
}

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
  target_mz <- target_mz[!is.na(target_mz)]

  if (length(target_mz) == 0L) {
    return(FALSE)
  }

  tol <- dalton + target_mz * ppm / 1e6
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
#' @return Filtered spectra object
#'
#' @examples NULL
perform_query <- function(all_mz, precursor_mz, frags, nls, dalton, ppm) {
  n <- length(all_mz)
  if (n == 0L) {
    return(logical(0))
  }

  keep <- rep(TRUE, n)

  if (length(frags) > 0L) {
    frags <- frags[!is.na(frags)]
    if (length(frags) > 0L) {
      frags_tol <- dalton + frags * ppm / 1e6
      keep <- vapply(
        all_mz,
        contains_all_mz_tol,
        FUN.VALUE = logical(1),
        target_mz = frags,
        tol = frags_tol
      )
    }
  }

  if (length(nls) > 0L && any(keep)) {
    valid_mask <- !is.na(precursor_mz)
    # Only test the spectra that are still candidates — skip anything
    # already ruled out by the fragment filter.
    idx <- which(keep & valid_mask)
    for (i in idx) {
      spec <- all_mz[[i]]
      if (length(spec) == 0L || all(is.na(spec))) {
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

  keep
}

#' @title Perform list of queries
#'
#' @param index Index of the ion list to process
#' @param ions_list List of ions for queries
#' @param query_data Data frame with columns `mz`, `precursor_mz`, and `skeleton`
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Data frame with target and value columns
#'
#' @examples NULL
perform_list_of_queries <- function(index, ions_list, query_data, dalton, ppm) {
  target <- names(ions_list)[index]
  if (is.null(target) || is.na(target) || target == "") {
    target <- paste0("query_", index)
  }

  ions <- ions_list[[index]]
  if (length(ions) == 0L || all(is.na(ions))) {
    return(tidytable::tidytable(target = target, value = character(0)))
  }

  is_frag <- stringi::stri_detect_fixed(ions, "_frag")
  is_nl <- stringi::stri_detect_fixed(ions, "_nl")

  frags <- if (any(is_frag)) {
    vals <- suppressWarnings(as.numeric(
      stringi::stri_replace_all_fixed(ions[is_frag], "_frag", "")
    ))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  nls <- if (any(is_nl)) {
    vals <- suppressWarnings(as.numeric(
      stringi::stri_replace_all_fixed(ions[is_nl], "_nl", "")
    ))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  if (length(frags) == 0L && length(nls) == 0L) {
    return(tidytable::tidytable(target = target, value = character(0)))
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
    return(tidytable::tidytable(target = target, value = character(0)))
  }

  value <- query_data$skeleton[keep]
  value <- value[!is.na(value)]

  tidytable::tidytable(target = target, value = value)
}

#' @title Perform list of queries
#'
#' @param ions_list List of ion combinations for queries
#' @param spectra Spectra object to search
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#' @param show_progress Logical: show progress bar (default TRUE)
#'
#' @return List of query results
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

  # Extracted once — this used to happen inside perform_query() on every
  # single query.
  query_data <- prepare_query_data(spectra)

  batch_size <- max(50L, min(200L, n %/% 10L + 1L))
  n_batches <- (n + batch_size - 1L) %/% batch_size

  for (batch in seq_len(n_batches)) {
    start_idx <- (batch - 1L) * batch_size + 1L
    end_idx <- min(batch * batch_size, n)

    for (i in start_idx:end_idx) {
      results[[i]] <- tryCatch(
        perform_list_of_queries(i, ions_list, query_data, dalton, ppm),
        error = function(e) {
          warning("Error in query ", i, ": ", e$message)
          tidytable::tidytable(
            target = names(ions_list)[i] %||% paste0("query_", i),
            value = character(0)
          )
        }
      )
      if (show_progress) show_progress(i, n)
    }
    # No more per-batch gc() — we're no longer creating/discarding large
    # Spectra subsets each iteration, so it was pure overhead here.
  }

  results
}
