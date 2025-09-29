#' @title Fix binned mzs
#'
#' @param binned_m Binned matrix
#' @param original_mzs Original m/z values from spectra
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#' @param decimals Number of decimal places for rounding
#'
#' @return Matrix with corrected m/z column names
#'
#' @examples NULL
fix_binned_mzs <- function(binned_m, original_mzs, dalton, ppm, decimals) {
  # Get all m/z values more efficiently
  all_mzs <- original_mzs |>
    Spectra::peaksData() |>
    Spectra::combinePeaksData(
      tolerance = dalton,
      ppm = ppm,
      peaks = "union"
    ) |>
    data.frame() |>
    tidytable::pull("mz")

  # Vectorized calculation of new column names
  old_colnames <- as.numeric(colnames(binned_m))
  
  # Pre-allocate result vector
  new_colnames <- numeric(length(old_colnames))
  
  for (i in seq_along(old_colnames)) {
    # Find matching m/z values within tolerance
    within_tolerance <- abs(all_mzs - old_colnames[i]) <= dalton
    if (any(within_tolerance)) {
      new_colnames[i] <- mean(all_mzs[within_tolerance])
    } else {
      new_colnames[i] <- old_colnames[i]
    }
  }

  # Update column names
  new_m <- binned_m
  colnames(new_m) <- new_colnames

  # More efficient pivoting and grouping
  new_m_df <- as.data.frame(new_m)
  new_m_df$rowname <- rownames(new_m)
  
  # Convert to long format
  long_df <- tidytable::pivot_longer(
    new_m_df,
    cols = -rowname,
    names_to = "name",
    values_to = "value",
    names_repair = "minimal"
  )
  
  # Clean and round names
  long_df$name <- round(as.numeric(long_df$name), decimals)
  
  # Group and summarize
  result_df <- long_df |>
    tidytable::filter(!is.na(name)) |>
    tidytable::group_by(rowname, name) |>
    tidytable::summarize(value = sum(value), .groups = "drop") |>
    tidytable::pivot_wider(names_from = name, values_from = value, values_fill = 0)
  
  # Convert back to matrix
  result_matrix <- as.matrix(result_df[, -1, drop = FALSE])
  rownames(result_matrix) <- result_df$rowname

  return(result_matrix)
}
