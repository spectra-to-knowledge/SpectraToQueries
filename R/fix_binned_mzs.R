#' @title Fix binned mzs
#'
#' @param binned_m binned_matrix
#' @param original_mzs original mzs
#' @param dalton dalton
#' @param ppm PPM
#' @param decimals decimals
#'
#' @return NULL
#'
#' @examples NULL
fix_binned_mzs <- function(binned_m, original_mzs, dalton, ppm, decimals) {
  all_mzs <- original_mzs |>
    Spectra::peaksData() |>
    Spectra::combinePeaksData(tolerance = dalton, ppm = ppm, peaks = "union") |>
    data.frame() |>
    tidytable::pull("mz")

  new_m <- binned_m

  new_colnames <- sapply(colnames(new_m), function(i) {
    val <- mean(all_mzs[abs(all_mzs - as.numeric(i)) <= dalton])
    if (!is.nan(val)) val else i
  })

  colnames(new_m) <- new_colnames

  new_m_new <- new_m |>
    data.frame() |>
    tibble::rownames_to_column() |>
    tidytable::pivot_longer(
      cols = -rowname,
      names_to = "name",
      values_to = "value"
    ) |>
    tidytable::mutate(name = gsub(
      pattern = ".[0-9]{1,2}$",
      replacement = "",
      x = name
    )) |>
    tidytable::mutate(name = round(
      x = gsub(
        pattern = "X",
        replacement = "",
        x = name,
        fixed = TRUE
      ) |> as.numeric(),
      digits = decimals
    ) |> as.character()) |>
    tidytable::group_by(rowname, name) |>
    tidytable::summarize(value = sum(value)) |>
    tidytable::pivot_wider(names_from = name, values_from = value) |>
    tibble::column_to_rownames() |>
    as.matrix()

  return(new_m_new)
}
