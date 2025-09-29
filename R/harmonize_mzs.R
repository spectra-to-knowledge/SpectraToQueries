#' @title Harmonize mzs
#'
#' @param spectra Spectra object
#' @param dalton Tolerance in Dalton
#' @param ppm Tolerance in parts per million
#'
#' @return Spectra object with harmonized m/z values
#'
#' @examples NULL
harmonize_mzs <- function(spectra, dalton, ppm) {
  spectra_new <- spectra
  
  # Get averaged intensities more efficiently
  averaged_intensities <- spectra_new |>
    Spectra::peaksData() |>
    Spectra::combinePeaksData(tolerance = dalton, ppm = ppm, peaks = "union") |>
    data.frame() |>
    tidytable::pull("mz")
  
  # Pre-calculate for efficiency
  n_spectra <- length(spectra_new)
  
  # Vectorized approach for finding closest matches
  for (i in seq_len(n_spectra)) {
    peaks_data <- spectra_new@backend@peaksData[[i]]
    mzs <- peaks_data[, 1]
    
    if (length(mzs) > 0 && length(averaged_intensities) > 0) {
      # Vectorized distance calculation using outer()
      distances <- abs(outer(mzs, averaged_intensities, "-"))
      matching_indices <- apply(distances, 1, which.min)
      min_distances <- apply(distances, 1, min)
      
      # Update m/z values where distance is within tolerance
      within_tolerance <- min_distances <= dalton
      mzs[within_tolerance] <- averaged_intensities[matching_indices[within_tolerance]]
      
      spectra_new@backend@peaksData[[i]][, 1] <- mzs
    }
  }
  
  return(spectra_new)
}
