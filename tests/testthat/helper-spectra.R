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
