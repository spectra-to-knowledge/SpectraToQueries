#' @title Prepare spectra data for querying
#' @param spectra Spectra object
#' @param intensity_min Minimal relative intensity threshold in `[0, 100]`.
#' @description Extracts m/z, precursor m/z, and normalized skeleton labels once
#'   as plain R vectors/lists so query precomputation avoids repeated S4
#'   dispatch.
#' @keywords internal
#' @noRd
prepare_query_data <- function(spectra, intensity_min = 0.0) {
  if (is.finite(intensity_min) && intensity_min > 0) {
    pks <- Spectra::peaksData(spectra)
    clean_mz <- lapply(pks, function(x) {
      if (nrow(x) == 0L) {
        return(numeric(0))
      }
      mz_v <- x[, 1L]
      int_v <- x[, 2L]
      fin <- is.finite(mz_v) & is.finite(int_v)
      if (!any(fin)) {
        return(numeric(0))
      }
      max_int <- max(int_v[fin])
      if (max_int <= 0) {
        return(numeric(0))
      }
      mz_v[fin & ((100 * int_v / max_int) >= intensity_min)]
    })
  } else {
    raw_mz <- as.list(Spectra::mz(spectra))
    clean_mz <- lapply(raw_mz, function(x) {
      if (all(is.finite(x))) x else x[is.finite(x)]
    })
  }
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

#' @title Parse fragment and neutral-loss ions from query character vector
#' @keywords internal
#' @noRd
parse_query_ions <- function(ions) {
  ions <- trimws(as.character(ions))
  ions <- ions[!is.na(ions) & nzchar(ions)]
  is_frag <- endsWith(ions, "_frag")
  is_nl <- endsWith(ions, "_nl")

  frags <- if (any(is_frag)) {
    raw_f <- substr(ions[is_frag], 1L, nchar(ions[is_frag]) - 5L)
    vals <- suppressWarnings(as.numeric(raw_f))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  nls <- if (any(is_nl)) {
    raw_n <- substr(ions[is_nl], 1L, nchar(ions[is_nl]) - 3L)
    vals <- suppressWarnings(as.numeric(raw_n))
    vals[is.finite(vals)]
  } else {
    numeric(0)
  }

  list(frags = frags, nls = nls)
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
    if (endsWith(u, "_frag")) {
      val <- suppressWarnings(as.numeric(substr(u, 1L, nchar(u) - 5L)))
      if (is.finite(val)) {
        tol <- dalton + val * (ppm / 1e6)
        hits[[u]] <- vapply(
          all_mz,
          function(spec_mz) any(abs(spec_mz - val) <= tol),
          logical(1L)
        )
      }
    } else if (endsWith(u, "_nl")) {
      val <- suppressWarnings(as.numeric(substr(u, 1L, nchar(u) - 3L)))
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
