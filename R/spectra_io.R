#' @title Return empty query table
#' @keywords internal
#' @noRd
empty_query_table <- function() {
  data.frame(
    skeleton = character(0),
    mcc = numeric(0),
    combined_ions = character(0),
    massql = character(0),
    stringsAsFactors = FALSE
  )
}

#' @title Write query table to TSV when export path is given
#' @keywords internal
#' @noRd
write_query_table <- function(queries, export) {
  if (!is.null(export)) {
    create_dir(dirname(export))
    utils::write.table(
      queries,
      file = export,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  invisible(queries)
}

#' @title Load or coerce spectra input
#' @keywords internal
#' @noRd
load_input_spectra <- function(spectra) {
  load_dataset <- function(name) {
    data_env <- new.env(parent = emptyenv())
    utils::data(
      list = name,
      package = "SpectraToQueries",
      envir = data_env
    )
    Spectra::Spectra(data_env[[name]])
  }

  if (is.null(spectra)) {
    message("No spectra given, loading example spectra.")
    mia_spectra <- load_dataset("mia_spectra_df")
  } else if (identical(spectra, "grouped")) {
    message("Loading grouped example spectra.")
    mia_spectra <- load_dataset("mia_spectra_grouped_df")
  } else if (is.character(spectra)) {
    if (!file.exists(spectra)) {
      stop(errorCondition(
        sprintf("Spectra file not found: '%s'.", spectra),
        class = "stq_file_not_found_error"
      ))
    }
    message("Loading spectra from file.")
    mia_spectra <- spectra |>
      MsBackendMgf::readMgf() |>
      Spectra::Spectra()
  } else if (is.data.frame(spectra)) {
    mia_spectra <- Spectra::Spectra(spectra)
  } else {
    mia_spectra <- spectra
  }

  vars <- Spectra::spectraVariables(mia_spectra)
  if (!("SKELETON" %in% vars)) {
    if ("skeleton" %in% vars) {
      mia_spectra$SKELETON <- mia_spectra$skeleton
    } else {
      abort_input("`spectra` must contain a `SKELETON` metadata variable.")
    }
  }
  if ("PRECURSOR_MZ" %in% vars) {
    mia_spectra$precursorMz <- as.numeric(mia_spectra$PRECURSOR_MZ)
  }
  mia_spectra
}
