#' @title Prepare spectra data for querying
#' @param spectra Spectra object
#' @param intensity_min Relative intensity threshold in `[0, 100]` (percent of
#'   the spectrum's base peak).
#' @description Extracts m/z, precursor m/z, and normalized skeleton labels once
#'   as plain R vectors/lists so query precomputation avoids repeated S4
#'   dispatch. When `intensity_min > 0`, keeps peaks with relative intensity
#'   `> min(intensity_min, 99)`, the rule MassQL applies to
#'   `INTENSITYPERCENT` (strict comparison, threshold capped at 99%).
#' @keywords internal
#' @noRd
prepare_query_data <- function(spectra, intensity_min = 0.0) {
  if (is.finite(intensity_min) && intensity_min > 0) {
    threshold <- min(intensity_min, 99)
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
      mz_v[fin & ((100 * int_v / max_int) > threshold)]
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

#' @title Parse fragment and neutral-loss ion tokens
#' @description Splits tokens of the form `<m/z>_frag` (fragment m/z) and
#'   `<mass>_nl` (neutral loss, Da) into numeric vectors. Empty and `NA` tokens
#'   are skipped; any other token without a positive finite value and an exact
#'   `_frag` or `_nl` suffix raises `stq_input_error`.
#' @keywords internal
#' @noRd
parse_query_ions <- function(ions) {
  ions <- trimws(as.character(ions))
  ions <- ions[!is.na(ions) & nzchar(ions)]
  is_frag <- endsWith(ions, "_frag")
  is_nl <- endsWith(ions, "_nl")
  vals <- suppressWarnings(as.numeric(
    substr(ions, 1L, nchar(ions) - ifelse(is_frag, 5L, 3L))
  ))
  bad <- !(is_frag | is_nl) | !is.finite(vals) | vals <= 0
  if (any(bad)) {
    abort_input(sprintf(
      paste(
        "Malformed ion token(s): %s.",
        "Expected `<positive number>_frag` or `<positive number>_nl`."
      ),
      paste0("'", ions[bad], "'", collapse = ", ")
    ))
  }
  list(frags = vals[is_frag], nls = vals[is_nl])
}

#' @title Precompute single-ion hit vectors across all spectra
#' @description For each ion token returns a logical vector over spectra. A
#'   fragment `f_frag` hits a spectrum with a peak within
#'   `dalton + f * ppm / 1e6` of `f`. A neutral loss `n_nl` hits a spectrum
#'   with finite precursor m/z `p` when `p - n > 0` and a peak lies within
#'   `tol` of `p - n`, with `tol = dalton` if `dalton > 0`, else
#'   `n * ppm / 1e6`. This is the window of the `MS2NL` condition emitted by
#'   [to_massql()] (MassQL applies `TOLERANCEPPM` to the loss value), so the
#'   reported MCC is reproducible by the emitted query. Intervals are closed.
#' @return Named list of logical vectors, one per element of `unique_ions`.
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
  hit_any <- function(spec_mz, target, tol) {
    any(abs(spec_mz - target) <= tol)
  }

  hits <- lapply(unique_ions, function(u) {
    parsed <- parse_query_ions(u)
    if (length(parsed$frags) == 1L) {
      f <- parsed$frags
      return(vapply(all_mz, hit_any, logical(1L), f, dalton + f * ppm / 1e6))
    }
    nl_tol <- if (dalton > 0) dalton else parsed$nls * ppm / 1e6
    vec <- logical(n_spec)
    for (idx in valid_idx) {
      tmz <- precursor_mz[idx] - parsed$nls
      if (tmz > 0) {
        vec[idx] <- hit_any(all_mz[[idx]], tmz, nl_tol)
      }
    }
    vec
  })
  names(hits) <- unique_ions
  hits
}
