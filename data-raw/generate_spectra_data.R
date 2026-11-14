# Run from the package root directory: Rscript data-raw/generate_spectra_data.R

pkg_root <- normalizePath(".", mustWork = TRUE)
source(file.path(pkg_root, "R", "spectra_utils.R"))

mgf_raw <- file.path(
  pkg_root,
  "data",
  "source",
  "spectra",
  "export_traitement_python_matchms-squelettes-tries.mgf"
)
mgf_grouped <- file.path(
  pkg_root,
  "data",
  "source",
  "spectra",
  "export_traitement_python_matchms-squelettes-tries-grouped.mgf"
)

message("Loading raw spectra from MGF...")
spectra_obj <- mgf_raw |>
  MsBackendMgf::readMgf() |>
  Spectra::Spectra()
message("Loaded ", length(spectra_obj), " raw spectra.")

mia_spectra_df <- spectra_to_df(spectra_obj)
message("Created mia_spectra_df with ", nrow(mia_spectra_df), " rows.")

message("Loading grouped spectra from MGF...")
spectra_grouped_obj <- mgf_grouped |>
  MsBackendMgf::readMgf() |>
  Spectra::Spectra()
message("Loaded ", length(spectra_grouped_obj), " grouped spectra.")

mia_spectra_grouped_df <- spectra_to_df(spectra_grouped_obj)
message(
  "Created mia_spectra_grouped_df with ",
  nrow(mia_spectra_grouped_df),
  " rows."
)

save(
  mia_spectra_df,
  file = file.path(pkg_root, "data", "mia_spectra_df.rda"),
  compress = "bzip2"
)
save(
  mia_spectra_grouped_df,
  file = file.path(pkg_root, "data", "mia_spectra_grouped_df.rda"),
  compress = "bzip2"
)
message("Saved sanitized datasets to data/.")
