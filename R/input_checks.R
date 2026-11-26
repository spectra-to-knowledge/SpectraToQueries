#' @title Abort with stq_input_error condition
#' @keywords internal
#' @noRd
abort_input <- function(msg) {
  stop(errorCondition(msg, class = "stq_input_error"))
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
    desc <- if (whole) {
      sprintf("whole number >= %s", lower)
    } else if (is.finite(lower) && is.finite(upper)) {
      sprintf("number in [%s, %s]", lower, upper)
    } else if (lower == 0) {
      "non-negative number"
    } else if (is.finite(lower)) {
      sprintf("number >= %s", lower)
    } else {
      "number"
    }
    abort_input(sprintf("`%s` must be a single finite %s.", name, desc))
  }
}
