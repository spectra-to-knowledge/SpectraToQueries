#' @title Compute Matthews Correlation Coefficient (vectorized)
#' @keywords internal
#' @noRd
calc_mcc_vec <- function(tp, fp, n_pos, total_spectra) {
  tp <- as.numeric(tp)
  fp <- as.numeric(fp)
  n_pos <- as.numeric(n_pos)
  total_spectra <- as.numeric(total_spectra)
  fn_ <- n_pos - tp
  tn <- total_spectra - tp - fp - fn_
  denom <- sqrt((tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_))
  ifelse(denom > 0, (tp * tn - fp * fn_) / denom, 0.0)
}
