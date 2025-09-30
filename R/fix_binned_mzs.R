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
  all_mzs <- original_mzs |>
    Spectra::peaksData(
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    Spectra::combinePeaksData(
      tolerance = dalton,
      ppm = ppm,
      peaks = "union",
      BPPARAM = BiocParallel::SerialParam()
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

  result_matrix <- new_m |>
    data.frame() |>
    tibble::rownames_to_column() |>
    tidytable::pivot_longer(
      cols = -rowname,
      names_to = "name",
      values_to = "value",
      names_repair = "minimal"
    ) |>
    tidytable::mutate(
      name = name |>
        gsub(pattern = "\\.[0-9]{1,2}$", replacement = "")
    ) |>
    # tidytable::filter(!name |>
    #   grepl(pattern = "V", fixed = TRUE)) |>
    tidytable::mutate(
      name = name |>
        gsub(
          pattern = "X",
          replacement = "",
          fixed = TRUE
        ) |>
        as.numeric() |>
        round(digits = decimals) |>
        as.character()
    ) |>
    tidytable::group_by(rowname, name) |>
    tidytable::summarize(value = sum(value)) |>
    tidytable::filter(!is.na(name)) |>
    tidytable::pivot_wider(names_from = name, values_from = value) |>
    tibble::column_to_rownames() |>
    as.matrix()

  return(result_matrix)
}
