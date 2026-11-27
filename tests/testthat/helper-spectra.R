make_spectra_df <- function(
  skeleton,
  mz,
  intensity = lapply(mz, function(x) rep(100.0, length(x))),
  ms_level = rep(2L, length(skeleton)),
  ...
) {
  df <- data.frame(
    msLevel = ms_level,
    SKELETON = skeleton,
    ...,
    stringsAsFactors = FALSE
  )
  df$mz <- mz
  df$intensity <- intensity
  df
}

load_miadb_df <- function(grouped = FALSE) {
  ds_name <- if (grouped) "mia_spectra_grouped_df" else "mia_spectra_df"
  env <- new.env(parent = emptyenv())
  utils::data(list = ds_name, package = "SpectraToQueries", envir = env)
  env[[ds_name]]
}

make_miadb_subset <- function(
  skeletons = c("Corynanthean", "Iboga", "Aspidosperma")
) {
  df <- load_miadb_df()
  df[df$SKELETON %in% skeletons, , drop = FALSE]
}
