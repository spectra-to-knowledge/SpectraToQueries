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
    skeleton = normalize_skeleton_labels(Spectra::spectraData(spectra)$SKELETON)
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
  is_frag <- is_frag_ion(ions)
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

#' @title Match fragment ions against a globally sorted peak m/z index
#' @keywords internal
#' @noRd
match_fragment_hits <- function(f_vals, all_mz, dalton, ppm, n_spec) {
  out <- rep(list(logical(n_spec)), length(f_vals))
  lens <- lengths(all_mz)
  if (!any(lens > 0L)) {
    return(out)
  }

  flat_mz <- unlist(all_mz, use.names = FALSE)
  flat_sp <- rep.int(seq_len(n_spec), lens)
  fin <- is.finite(flat_mz)
  if (!any(fin)) {
    return(out)
  }

  flat_mz <- flat_mz[fin]
  flat_sp <- flat_sp[fin]
  ord <- order(flat_mz)
  sorted_mz <- flat_mz[ord]
  sorted_sp <- flat_sp[ord]

  f_tols <- mz_tolerance(f_vals, dalton, ppm)
  pad <- pmax(f_tols, abs(f_vals), 1) * 1e-9 + 1e-9
  lo_vec <- findInterval(f_vals - f_tols - pad, sorted_mz, left.open = TRUE) +
    1L
  hi_vec <- findInterval(f_vals + f_tols + pad, sorted_mz)

  for (j in seq_along(f_vals)) {
    lo <- lo_vec[j]
    hi <- hi_vec[j]
    if (hi >= lo) {
      idx <- lo:hi
      matched <- idx[abs(sorted_mz[idx] - f_vals[j]) <= f_tols[j]]
      if (length(matched) > 0L) {
        out[[j]][sorted_sp[matched]] <- TRUE
      }
    }
  }
  out
}

#' @title Match neutral-loss ions against a globally sorted loss index
#' @keywords internal
#' @noRd
match_nl_hits <- function(
  n_vals,
  all_mz,
  precursor_mz,
  valid_idx,
  dalton,
  ppm,
  n_spec
) {
  out <- rep(list(logical(n_spec)), length(n_vals))
  nl_lens <- lengths(all_mz[valid_idx])
  if (!any(nl_lens > 0L)) {
    return(out)
  }

  nl_mz <- unlist(all_mz[valid_idx], use.names = FALSE)
  nl_sp <- rep.int(valid_idx, nl_lens)
  nl_prec <- rep.int(precursor_mz[valid_idx], nl_lens)
  fin_nl <- is.finite(nl_mz)
  if (!any(fin_nl)) {
    return(out)
  }

  nl_mz <- nl_mz[fin_nl]
  nl_sp <- nl_sp[fin_nl]
  nl_prec <- nl_prec[fin_nl]
  obs_nl <- nl_prec - nl_mz
  ord_nl <- order(obs_nl)
  sorted_obs_nl <- obs_nl[ord_nl]
  sorted_nl_mz <- nl_mz[ord_nl]
  sorted_nl_prec <- nl_prec[ord_nl]
  sorted_nl_sp <- nl_sp[ord_nl]

  n_tols <- nl_tolerance(n_vals, dalton, ppm)
  max_prec <- max(abs(precursor_mz[valid_idx]))
  pad_nl <- pmax(n_tols, abs(n_vals), max_prec, 1) * 1e-9 + 1e-9
  lo_nl <- findInterval(
    n_vals - n_tols - pad_nl,
    sorted_obs_nl,
    left.open = TRUE
  ) +
    1L
  hi_nl <- findInterval(n_vals + n_tols + pad_nl, sorted_obs_nl)

  for (j in seq_along(n_vals)) {
    lo <- lo_nl[j]
    hi <- hi_nl[j]
    if (hi >= lo) {
      idx <- lo:hi
      tmz <- sorted_nl_prec[idx] - n_vals[j]
      matched <- idx[tmz > 0 & abs(sorted_nl_mz[idx] - tmz) <= n_tols[j]]
      if (length(matched) > 0L) {
        out[[j]][sorted_nl_sp[matched]] <- TRUE
      }
    }
  }
  out
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
  n_ions <- length(unique_ions)
  if (n_ions == 0L) {
    return(stats::setNames(list(), character(0)))
  }

  ions_clean <- trimws(as.character(unique_ions))
  parsed <- parse_query_ions(ions_clean)
  non_empty <- !is.na(ions_clean) & nzchar(ions_clean)
  is_frag <- is_frag_ion(ions_clean)
  frag_pos <- which(non_empty & is_frag)
  nl_pos <- which(non_empty & !is_frag)

  n_spec <- length(all_mz)
  hits <- stats::setNames(rep(list(logical(n_spec)), n_ions), unique_ions)

  if (length(frag_pos) > 0L) {
    hits[frag_pos] <- match_fragment_hits(
      parsed$frags,
      all_mz,
      dalton,
      ppm,
      n_spec
    )
  }

  valid_idx <- which(valid_prec_and_spec)
  if (length(nl_pos) > 0L && length(valid_idx) > 0L) {
    hits[nl_pos] <- match_nl_hits(
      parsed$nls,
      all_mz,
      precursor_mz,
      valid_idx,
      dalton,
      ppm,
      n_spec
    )
  }

  hits
}
