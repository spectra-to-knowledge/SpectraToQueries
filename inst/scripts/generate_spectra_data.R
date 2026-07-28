# Generate package data from raw MGF spectra files
# This script demonstrates how to create the mia_spectra_df and mia_spectra_grouped_df
# data objects from the raw MGF source files in data/source/spectra/
#
# To reproduce:
#   Rscript inst/scripts/generate_spectra_data.R

library(Spectra)
library(MsBackendMgf)

# Paths (assumes script is run from package root)
pkg_root <- here::here()
mgf_raw <- file.path(
  pkg_root,
  "data/source/spectra/export_traitement_python_matchms-squelettes-tries.mgf"
)
mgf_grouped <- file.path(
  pkg_root,
  "data/source/spectra/export_traitement_python_matchms-squelettes-tries-grouped.mgf"
)

cat("========================================\n")
cat("Generating spectra data from MGF files\n")
cat("========================================\n\n")

# 1. Load raw spectra from MGF file
cat("Step 1: Loading raw spectra from MGF...\n")
spectra_obj <- mgf_raw |>
  MsBackendMgf::readMgf() |>
  Spectra::Spectra()
cat("  ✓ Loaded", length(spectra_obj), "raw spectra\n\n")

# 2. Convert raw Spectra to portable data.frame
cat("Step 2: Converting spectra to portable data.frame...\n")
mia_spectra_df <- spectra_obj |>
  spectra_to_df()
cat("  ✓ Created data.frame with", nrow(mia_spectra_df), "rows\n")
cat(
  "  ✓ Columns:",
  paste(head(colnames(mia_spectra_df), 5), collapse = ", "),
  "...\n\n"
)

# 3. Load grouped spectra from MGF file
cat("Step 3: Loading grouped spectra from MGF...\n")
spectra_grouped_obj <- mgf_grouped |>
  MsBackendMgf::readMgf() |>
  Spectra::Spectra()
cat("  ✓ Loaded", length(spectra_grouped_obj), "grouped spectra\n\n")

# 4. Convert grouped Spectra to portable data.frame
cat("Step 4: Converting grouped spectra to portable data.frame...\n")
mia_spectra_grouped_df <- spectra_grouped_obj |>
  spectra_to_df()
cat("  ✓ Created data.frame with", nrow(mia_spectra_grouped_df), "rows\n\n")

# 5. Save as RDA files in data/ directory
cat("Step 5: Saving as RDA files...\n")
usethis::use_data(mia_spectra_df, overwrite = TRUE)
usethis::use_data(mia_spectra_grouped_df, overwrite = TRUE)
cat("  ✓ Saved data/mia_spectra_df.rda\n")
cat("  ✓ Saved data/mia_spectra_grouped_df.rda\n\n")

cat("========================================\n")
cat("✅ Data generation complete!\n")
cat("========================================\n")
cat("\nGenerated files:\n")
cat("  - data/mia_spectra_df.rda (raw spectra as portable data.frame)\n")
cat(
  "  - data/mia_spectra_grouped_df.rda (grouped spectra as portable data.frame)\n"
)
cat("\nThese can be loaded via:\n")
cat("  utils::data(mia_spectra_df, package = 'SpectraToQueries')\n")
cat("  utils::data(mia_spectra_grouped_df, package = 'SpectraToQueries')\n")
