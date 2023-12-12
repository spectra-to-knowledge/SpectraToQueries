#' Title
#'
#' @param index
#' @param ions_list
#' @param spectra
#'
#' @return
#' @export
#'
#' @examples
perform_list_of_queries <- function(index, ions_list, spectra) {
  message("Extract the fragments of the list.")
  frags <-
    ions_list[[index]][ions_list[[index]] |> grepl(pattern = "frag")] |>
    gsub(pattern = "_frag", replacement = "") |>
    as.numeric()

  message("Extract the neutral losses of the list.")
  nls <-
    ions_list[[index]][ions_list[[index]] |> grepl(pattern = "nl")] |>
    gsub(pattern = "_nl", replacement = "") |>
    as.numeric()

  message("Perform the query.")
  matching_spectra <- spectra |>
    perform_query(frags = frags, nls = nls)

  message("Return the corresponding tidytable.")
  return(tidytable::tidytable(
    target = names(ions_list)[index],
    value = c(
      matching_spectra |>
        Spectra::spectraData() |>
        tidytable::pull(SKELETON)
    )
  ))
}
