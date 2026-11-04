#' Example mass spectrometry spectra from MIADB (as portable data.frame)
#'
#' A `data.frame` containing 321 MS/MS spectra of monoterpene
#' indole alkaloids (MIAs) from the MIADB (Monoterpene Indole Alkaloid
#' Database) across 81 fine-grained chemical skeletons (`SKELETON`). This is
#' stored as a portable `data.frame` with `mz` and `intensity` peaks as list
#' columns, and is automatically converted to a [Spectra::Spectra] object
#' when used in [spectra_to_queries()].
#'
#' @format A `data.frame` with 321 rows and 37 columns:
#'   - **Metadata columns**: `msLevel`, `rtime`, `acquisitionNum`, `TITLE`,
#'     `PRECURSOR_MZ`, `SKELETON` (81 unique skeletons), and other MS/MS
#'     spectrum attributes
#'   - **`mz`** (list): Fragment m/z values for each spectrum
#'   - **`intensity`** (list): Fragment intensities for each spectrum
#'
#'   Total of 24,255 peaks across all 321 spectra (median 66 peaks per spectrum).
#'
#' @source Data extracted from MIADB, a public spectral library
#'   of monoterpene indole alkaloids.
#'
#'   Accompanying publication:
#'   Szwarc, S., Rutz, A., Lee, K. et al. (2025). "Translating community-wide
#'   spectral library into actionable chemical knowledge: a proof of concept
#'   with monoterpene indole alkaloids." *Journal of Cheminformatics*, 17:62.
#'   \doi{10.1186/s13321-025-01009-0}
#'
#' @details
#' The `mia_spectra_df` object represents the spectral data with fine-grained
#' skeleton annotations. Use `mia_spectra_grouped_df` for the variant with
#' broader grouped skeleton annotations.
#'
#' This portable `data.frame` format enables cross-platform distribution and
#' is automatically converted to a [Spectra::Spectra] object by
#' [spectra_to_queries()].
#'
#' @examples
#' utils::data(mia_spectra_df, package = "SpectraToQueries")
#' nrow(mia_spectra_df) # 321 spectra
"mia_spectra_df"

#' Grouped mass spectrometry spectra from MIADB (as portable data.frame)
#'
#' A `data.frame` containing the 321 MS/MS spectra from the MIADB
#' (Monoterpene Indole Alkaloid Database) with `SKELETON` annotations grouped
#' into 53 broader chemical skeleton classes. Stored as a portable `data.frame`
#' with `mz` and `intensity` peaks as list columns.
#'
#' @format A `data.frame` with 321 rows and 36 columns:
#'   - **Metadata columns**: `msLevel`, `rtime`, `acquisitionNum`, `TITLE`,
#'     `PRECURSOR_MZ`, `SKELETON` (53 grouped skeleton classes), and other
#'     MS/MS spectrum attributes
#'   - **`mz`** (list): Fragment m/z values for each spectrum
#'   - **`intensity`** (list): Fragment intensities for each spectrum
#'
#'   Total of 24,255 peaks across all 321 spectra (median 66 peaks per spectrum).
#'
#' @source Data extracted from MIADB, a public spectral library
#'   of monoterpene indole alkaloids.
#'
#'   Accompanying publication:
#'   Szwarc, S., Rutz, A., Lee, K. et al. (2025). "Translating community-wide
#'   spectral library into actionable chemical knowledge: a proof of concept
#'   with monoterpene indole alkaloids." *Journal of Cheminformatics*, 17:62.
#'   \doi{10.1186/s13321-025-01009-0}
#'
#' @details
#' In this dataset, closely related sub-skeletons are merged into broader
#' skeleton classes (e.g., `corynanthean_1+corynanthean_2+...`) to evaluate
#' diagnostic ion queries at a higher structural hierarchy level.
#'
#' This portable `data.frame` format enables cross-platform distribution and
#' is automatically converted to a [Spectra::Spectra] object by
#' [spectra_to_queries()].
#'
#' @examples
#' utils::data(mia_spectra_grouped_df, package = "SpectraToQueries")
#' nrow(mia_spectra_grouped_df) # 321 spectra
#' length(unique(mia_spectra_grouped_df$SKELETON)) # 53 grouped classes
"mia_spectra_grouped_df"
