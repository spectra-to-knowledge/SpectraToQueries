#' @title Perform list of queries
#'
#' @include perform_query.R
#'
#' @param index Index
#' @param ions_list Ions list
#' @param spectra Spectra
#' @param dalton Dalton
#' @param ppm PPM
#'
#' @return NULL
#'
#' @examples NULL
perform_list_of_queries <- function(index, ions_list, spectra, dalton, ppm) {
  target <- names(ions_list)[index]
  frags <- ions_list[[index]][
    ions_list[[index]] |>
      grepl(pattern = "frag")
  ] |>
    gsub(pattern = "_frag", replacement = "") |>
    as.numeric()
  nls <- ions_list[[index]][
    ions_list[[index]] |>
      grepl(pattern = "nl")
  ] |>
    gsub(pattern = "_nl", replacement = "") |>
    as.numeric()
  value <- spectra |>
    perform_query(frags = frags, nls = nls, dalton, ppm) |>
    Spectra::spectraData() |>
    tidytable::pull(SKELETON) |>
    gsub(
      pattern = "+",
      replacement = ".",
      fixed = TRUE
    )
  return(tidytable::tidytable(target = target, value = value))
}

#' @title Perform list of queries (progress)
#'
#' @param ions_list Ions list
#' @param spectra Spectra
#' @param dalton Dalton
#' @param ppm PPM
#'
#' @return NULL
#'
#' @examples NULL
perform_list_of_queries_progress <- function(ions_list, spectra, dalton, ppm) {
  purrr::map(
    .progress = TRUE,
    .x = seq_along(ions_list),
    .f = function(index, ions_list, spectra, dalton, ppm) {
      perform_list_of_queries(index, ions_list, spectra, dalton, ppm)
    },
    ions_list = ions_list,
    spectra = spectra,
    dalton = dalton,
    ppm = ppm
  )
}
