#' Title
#'
#' @param index Index
#' @param ions_list Ions list
#' @param spectra Spectra
#'
#' @return NULL
#' @export
#'
#' @examples NULL
perform_list_of_queries <- function(index, ions_list, spectra) {
  return(tidytable::tidytable(
    target = names(ions_list)[index],
    value = c(
      spectra |>
        perform_query(
          frags = ions_list[[index]][ions_list[[index]] |>
            grepl(pattern = "frag")] |>
            gsub(pattern = "_frag", replacement = "") |>
            as.numeric(),
          nls = ions_list[[index]][ions_list[[index]] |>
            grepl(pattern = "nl")] |>
            gsub(pattern = "_nl", replacement = "") |>
            as.numeric()
        ) |>
        Spectra::spectraData() |>
        tidytable::pull(SKELETON)
    )
  ))
}
