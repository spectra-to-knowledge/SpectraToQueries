testthat::test_that("namespace exports only public API", {
  testthat::expect_identical(
    sort(getNamespaceExports("SpectraToQueries")),
    c("spectra_to_queries", "to_massql")
  )
})

testthat::test_that("spectra_to_queries matches MIADB subset baseline", {
  utils::data(
    "mia_spectra_df",
    package = "SpectraToQueries",
    envir = environment()
  )
  keep <- mia_spectra_df$SKELETON %in%
    c("Corynanthean", "Iboga", "Aspidosperma")
  res <- suppressMessages(spectra_to_queries(
    spectra = mia_spectra_df[keep, ],
    export = NULL,
    ions_max = 4L
  ))

  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_identical(
    colnames(res),
    c("skeleton", "mcc", "combined_ions", "massql")
  )
  testthat::expect_identical(
    res$skeleton,
    c("Aspidosperma", "Iboga", "Corynanthean")
  )
  testthat::expect_equal(
    res$mcc,
    c(0.6947, 0.6896, 0.5021),
    tolerance = 1e-12
  )
  testthat::expect_identical(
    res$combined_ions,
    c(
      "134.0951_frag & 199.1217_nl",
      "121.091_frag & 224.1054_frag",
      "116.0562_frag & 143.0735_nl"
    )
  )
})
