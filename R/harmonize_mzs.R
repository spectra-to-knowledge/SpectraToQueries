#' @title Harmonize mzs
#'
#' @param spectra Spectra object
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Spectra with harmonized m/z values
#'
#' @examples NULL
harmonize_mzs <- function(spectra, dalton, ppm) {
  spectra_new <- spectra

  # Averaged m/z values from combined spectra
  averaged_intensities <- spectra_new |>
    Spectra::peaksData(
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    Spectra::combinePeaksData(
      tolerance = dalton,
      ppm = ppm,
      peaks = "union",
      BPPARAM = BiocParallel::SerialParam()
    ) |>
    data.frame() |>
    tidytable::pull("mz")

  # Process all spectra more efficiently
  for (i in seq_along(spectra_new)) {
    # Direct access to mz values without intermediate data frame
    mzs <- spectra_new@backend@peaksData[[i]][, 1]

    # Vectorized matching for all m/z values at once
    if (length(mzs) > 0) {
      # For each mz, find the closest averaged value
      for (j in seq_along(mzs)) {
        # Calculate distances to all averaged values
        diffs <- abs(mzs[j] - averaged_intensities)
        min_idx <- which.min(diffs)

        # Only update if within tolerance
        if (diffs[min_idx] <= dalton) {
          mzs[j] <- averaged_intensities[min_idx]
        }
      }

      # Update the m/z values in place
      spectra_new@backend@peaksData[[i]][, 1] <- mzs
    }
  }

  return(spectra_new)
}
