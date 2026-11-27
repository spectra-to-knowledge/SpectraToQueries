#' @title Abort with stq_input_error condition
#' @keywords internal
#' @noRd
abort_input <- function(msg) {
  stop(errorCondition(msg, class = "stq_input_error"))
}

#' @title Format numeric bound description for validation errors
#' @keywords internal
#' @noRd
format_num_constraint <- function(lower, upper, whole) {
  if (whole) {
    return(sprintf("whole number >= %s", lower))
  }
  if (is.finite(lower) && is.finite(upper)) {
    return(sprintf("number in [%s, %s]", lower, upper))
  }
  if (lower == 0) {
    return("non-negative number")
  }
  if (is.finite(lower)) {
    return(sprintf("number >= %s", lower))
  }
  "number"
}

#' @title Validate a scalar finite numeric argument
#' @keywords internal
#' @noRd
check_num <- function(x, name, lower = -Inf, upper = Inf, whole = FALSE) {
  ok <- is.numeric(x) &&
    length(x) == 1L &&
    is.finite(x) &&
    x >= lower &&
    x <= upper &&
    (!whole || x == round(x))
  if (!ok) {
    abort_input(sprintf(
      "`%s` must be a single finite %s.",
      name,
      format_num_constraint(lower, upper, whole)
    ))
  }
}

#' @title Validate tolerance and relative intensity arguments
#' @keywords internal
#' @noRd
check_tolerance_args <- function(dalton, ppm, intensity_min) {
  check_num(dalton, "dalton", lower = 0)
  check_num(ppm, "ppm", lower = 0)
  check_num(intensity_min, "intensity_min", lower = 0, upper = 100)
}

#' @title Validate inputs for `spectra_to_queries`
#' @keywords internal
#' @noRd
validate_stq_inputs <- function(
  spectra,
  export,
  dalton,
  decimals,
  intensity_min,
  ions_max,
  n_skel_min,
  n_spec_min,
  ppm,
  mcc_min,
  precision_min,
  recall_min,
  zero_val
) {
  is_valid_str <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }

  valid_spectra <- is.null(spectra) ||
    is.data.frame(spectra) ||
    inherits(spectra, "Spectra") ||
    is_valid_str(spectra)
  if (!valid_spectra) {
    abort_input(
      paste(
        "`spectra` must be NULL, a non-empty string,",
        "a data.frame, or a Spectra object."
      )
    )
  }
  if (!is.null(export) && !is_valid_str(export)) {
    abort_input("`export` must be NULL or a non-empty character string.")
  }

  check_tolerance_args(dalton, ppm, intensity_min)
  check_num(decimals, "decimals", lower = 0, whole = TRUE)
  check_num(ions_max, "ions_max", lower = 1, whole = TRUE)
  check_num(n_skel_min, "n_skel_min", lower = 1, whole = TRUE)
  check_num(n_spec_min, "n_spec_min", lower = 1, whole = TRUE)
  check_num(precision_min, "precision_min", lower = 0, upper = 1)
  check_num(recall_min, "recall_min", lower = 0, upper = 1)
  check_num(mcc_min, "mcc_min", lower = -1, upper = 1)
  check_num(zero_val, "zero_val")
  invisible(TRUE)
}

#' @title Compute fragment m/z tolerance window in Dalton
#' @keywords internal
#' @noRd
mz_tolerance <- function(mz, dalton, ppm) {
  dalton + mz * (ppm / 1e6)
}

#' @title Compute neutral-loss tolerance window in Dalton
#' @keywords internal
#' @noRd
nl_tolerance <- function(nl, dalton, ppm) {
  if (dalton > 0) dalton else nl * (ppm / 1e6)
}

#' @title Normalize skeleton labels by replacing '+' with '.'
#' @keywords internal
#' @noRd
normalize_skeleton_labels <- function(x) {
  gsub("+", ".", as.character(x), fixed = TRUE)
}

#' @title Identify fragment ion tokens by exact '_frag' suffix
#' @keywords internal
#' @noRd
is_frag_ion <- function(ions) {
  endsWith(ions, "_frag")
}
