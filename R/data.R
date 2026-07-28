#' Example mass spectrometry spectra from MIADB (as portable data.frame)
#'
#' A `data.frame` containing 321 spectra of monoterpene
#' indole alkaloids (MIAs) from the MIADB (Monoterpene Indole Alkaloid
#' Database). This is stored as a portable data.frame with m/z and intensity
#' peaks as list columns, and is automatically converted to a `Spectra` object
#' when used in `spectra_to_queries()`.
#'
#' @format A `data.frame` with 321 rows and 37 columns:
#'   - **Metadata columns**: `msLevel`, `rtime`, `acquisitionNum`, `TITLE`,
#'     `PRECURSOR_MZ`, `SKELETON`, and other MS/MS spectrum attributes
#'   - **`mz`** (list): Fragment m/z values for each spectrum
#'   - **`intensity`** (list): Fragment intensities for each spectrum
#'
#'   Total of approximately 48,510 peaks across all spectra.
#'
#' @source Data extracted from MIADB, a public spectral library
#'   of monoterpene indole alkaloids. Original source:
#'   https://miadb.imsb.qcb.usp.br/
#'
#'   Accompanying publication:
#'   Szwarc, S., Rutz, A., Lee, K. et al. (2025). "Translating community-wide
#'   spectral library into actionable chemical knowledge: a proof of concept
#'   with monoterpene indole alkaloids." *Journal of Cheminformatics*, 17:62.
#'   https://doi.org/10.1186/s13321-025-01009-0
#'
#' @details
#' The `mia_spectra_df` object represents the raw spectral data with
#' minimal preprocessing. Use `mia_spectra_grouped_df` for the pre-grouped variant.
#'
#' This portable data.frame format enables export to CSV, Parquet, JSON, and
#' other formats for cross-platform distribution via r-universe. It is
#' automatically converted to a `Spectra` object by `spectra_to_queries()`.
#'
#' @examples
#' \dontrun{
#'   # Load the raw spectra
#'   utils::data(mia_spectra_df)
#'   nrow(mia_spectra_df)  # 321 spectra
#'
#'   # Convert to Spectra object (for advanced users)
#'   # spectra_obj <- SpectraToQueries:::df_to_spectra(mia_spectra_df)
#' }
"mia_spectra_df"

#' Grouped mass spectrometry spectra from MIADB (as portable data.frame)
#'
#' A `data.frame` containing 321 grouped mass spectra from the MIADB
#' (Monoterpene Indole Alkaloid Database), combined by spectrum title.
#' This is the primary example dataset used in the `spectra_to_queries()`
#' function. Stored as a portable data.frame with m/z and intensity peaks
#' as list columns.
#'
#' @format A `data.frame` with 321 rows and 36 columns:
#'   - **Metadata columns**: `msLevel`, `rtime`, `acquisitionNum`, `TITLE`,
#'     `PRECURSOR_MZ`, `SKELETON`, and other MS/MS spectrum attributes
#'   - **`mz`** (list): Fragment m/z values for each spectrum
#'   - **`intensity`** (list): Fragment intensities for each spectrum
#'
#'   Total of approximately 48,510 peaks across all spectra (median ~132 peaks per spectrum).
#'
#' @source Data extracted and pre-grouped from MIADB, a public spectral
#'   library of monoterpene indole alkaloids. Original source:
#'   https://miadb.imsb.qcb.usp.br/
#'
#'   Accompanying publication:
#'   Szwarc, S., Rutz, A., Lee, K. et al. (2025). "Translating community-wide
#'   spectral library into actionable chemical knowledge: a proof of concept
#'   with monoterpene indole alkaloids." *Journal of Cheminformatics*, 17:62.
#'   https://doi.org/10.1186/s13321-025-01009-0
#'
#' @details
#' This dataset is grouped by spectrum title (via `Spectra::combineSpectra()`)
#' to simplify downstream analysis. It is the recommended example for
#' demonstrating the `spectra_to_queries()` workflow.
#'
#' This portable data.frame format enables export to CSV, Parquet, JSON, and
#' other formats for cross-platform distribution via r-universe. It is
#' automatically converted to a `Spectra` object by `spectra_to_queries()`.
#'
#' @examples
#' \dontrun{
#'   # Load the grouped spectra
#'   utils::data(mia_spectra_grouped_df)
#'   nrow(mia_spectra_grouped_df)  # 321 spectra
#'   table(mia_spectra_grouped_df$SKELETON)  # Skeleton distribution
#'
#'   # Convert to Spectra object (for advanced users)
#'   # spectra_obj <- SpectraToQueries:::df_to_spectra(mia_spectra_grouped_df)
#' }
"mia_spectra_grouped_df"
