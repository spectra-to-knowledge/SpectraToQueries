library(testthat)

testthat::test_that(desc = "Spectra to queries", code = {
  message("\n")
  spectra_to_queries()
  succeed()
})
