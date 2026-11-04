testthat::test_that("combine_ions_minimal simplifies and sorts ion expressions", {
  testthat::expect_identical(combine_ions_minimal(list()), "")
  testthat::expect_identical(combine_ions_minimal(list(c("", NA))), "")

  # Sorting: fragments ascending then neutral losses ascending
  testthat::expect_identical(
    combine_ions_minimal(list(c("200_frag", "50_nl", "100_frag", "20_nl"))),
    "100_frag & 200_frag & 20_nl & 50_nl"
  )

  # Common ions factored out with OR branches
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "200_frag"),
      c("100_frag", "50_nl", "60_nl")
    )),
    "100_frag & 200_frag | (50_nl & 60_nl)"
  )

  # Subset query absorbs superset queries into common ions
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag"),
      c("100_frag", "200_frag")
    )),
    "100_frag"
  )

  # Duplicate ions within a query and order-variant remaining sets are deduplicated
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "100_frag", "200_frag"),
      c("100_frag", "300_frag")
    )),
    "100_frag & 200_frag | 300_frag"
  )
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "200_frag"),
      c("300_frag", "400_frag"),
      c("400_frag", "300_frag")
    )),
    "(100_frag & 200_frag) | (300_frag & 400_frag)"
  )
})

testthat::test_that("create_dir creates parent or target directories", {
  tmp_base <- file.path(tempdir(), paste0("stq_test_", Sys.getpid()))
  on.exit(unlink(tmp_base, recursive = TRUE), add = TRUE)

  file_target <- file.path(tmp_base, "sub1", "queries.tsv")
  testthat::expect_message(create_dir(file_target), "created")
  testthat::expect_true(dir.exists(file.path(tmp_base, "sub1")))

  dir_target <- file.path(tmp_base, "sub2")
  testthat::expect_message(create_dir(dir_target), "created")
  testthat::expect_true(dir.exists(dir_target))
})

testthat::test_that("filter_matrix and normalize_peaks handle standard and edge cases", {
  mat <- matrix(
    c(1, 0, 0, 1, 2, 0, 1, 2, 3),
    nrow = 3,
    dimnames = list(c("A", "B", "C"), c("100", "200", "300"))
  )
  filtered <- filter_matrix(mat, n = 2L)
  testthat::expect_identical(colnames(filtered), c("200", "300"))
  testthat::expect_identical(nrow(filtered), 3L)

  norm_fn <- normalize_peaks()
  pks <- cbind(mz = c(100, 200), intensity = c(25, 50))
  res <- norm_fn(pks)
  testthat::expect_equal(res[, 2L], c(50, 100))

  empty_pks <- matrix(numeric(0), nrow = 0L, ncol = 2L)
  testthat::expect_identical(norm_fn(empty_pks), empty_pks)

  zero_pks <- cbind(mz = c(100, 200), intensity = c(0, 0))
  testthat::expect_identical(norm_fn(zero_pks), zero_pks)
})

testthat::test_that("generate_combinations and show_progress work as expected", {
  combs <- generate_combinations(c("a", "b", "c"), max_ions = 2L)
  testthat::expect_length(combs, 6L)
  testthat::expect_identical(generate_combinations(character(0), 2L), list())
  testthat::expect_identical(generate_combinations(c("a", "b"), 0L), list())

  ions_list <- list(s1 = c("a", "b"), s2 = c("c"))
  res <- generate_combinations_progress(
    c("s1", "s2"),
    ions_list = ions_list,
    max_ions = 2L,
    show_progress = FALSE
  )
  testthat::expect_length(res, 2L)
  testthat::expect_length(res[[1L]], 3L)
  testthat::expect_length(res[[2L]], 1L)

  testthat::expect_null(show_progress(0L, 0L))
  testthat::expect_message(show_progress(10L, 10L), "Done!")
})

testthat::test_that("spectra_to_df, df_to_spectra, create_matrix, harmonize_mzs, fix_binned_mzs work", {
  utils::data(
    mia_spectra_df,
    package = "SpectraToQueries",
    envir = environment()
  )
  sub_df <- mia_spectra_df[1:10, ]
  sp <- df_to_spectra(sub_df)
  testthat::expect_s4_class(sp, "Spectra")
  testthat::expect_identical(length(sp), 10L)

  roundtrip_df <- spectra_to_df(sp)
  testthat::expect_identical(nrow(roundtrip_df), 10L)
  testthat::expect_equal(roundtrip_df$mz, sub_df$mz)

  sp_harm <- harmonize_mzs(sp, dalton = 0.01, ppm = 30)
  testthat::expect_s4_class(sp_harm, "Spectra")

  sp_bin_sparse <- Spectra::bin(sp_harm, binSize = 0.01, zero.rm = TRUE) |>
    Spectra::applyProcessing()
  sp_bin_dense <- Spectra::bin(sp_harm, binSize = 0.01, zero.rm = FALSE) |>
    Spectra::applyProcessing()

  m_sparse <- create_matrix(sp_bin_sparse, name = sp_bin_sparse$SKELETON)
  m_dense <- create_matrix(sp_bin_dense, name = sp_bin_dense$SKELETON)
  testthat::expect_equal(m_sparse, m_dense)

  m_filt <- filter_matrix(m_sparse, n = 2L)
  m_fixed <- fix_binned_mzs(
    m_filt,
    original_mzs = sp_harm,
    dalton = 0.01,
    ppm = 30,
    decimals = 4L
  )
  testthat::expect_identical(nrow(m_fixed), 10L)
  testthat::expect_true(ncol(m_fixed) <= ncol(m_filt))
})

testthat::test_that("create_matrix detects middle-spectrum grid differences and unsorted m/z", {
  df_mid_diff <- data.frame(
    msLevel = c(2L, 2L, 2L),
    SKELETON = c("S1", "S2", "S3")
  )
  df_mid_diff$mz <- list(c(100, 200), c(100, 300), c(100, 200))
  df_mid_diff$intensity <- list(c(10, 20), c(30, 40), c(50, 60))
  sp_mid_diff <- df_to_spectra(df_mid_diff)

  mat_mid <- create_matrix(
    sp_mid_diff,
    zero_val = 0,
    name = df_mid_diff$SKELETON
  )
  expected_mid <- matrix(
    c(
      10,
      30,
      50,
      20,
      0,
      60,
      0,
      40,
      0
    ),
    nrow = 3L,
    ncol = 3L,
    dimnames = list(c("S1", "S2", "S3"), c("100", "200", "300"))
  )
  testthat::expect_identical(mat_mid, expected_mid)
})

testthat::test_that("query matching helpers work accurately", {
  spec_mz <- c(100.001, 150.0, 200.0)
  testthat::expect_true(contains_all_mz(
    spec_mz,
    c(100.0, 150.0),
    dalton = 0.01,
    ppm = 10
  ))
  testthat::expect_false(contains_all_mz(
    spec_mz,
    c(100.0, 300.0),
    dalton = 0.01,
    ppm = 10
  ))
  testthat::expect_false(contains_all_mz(
    numeric(0),
    100.0,
    dalton = 0.01,
    ppm = 10
  ))

  all_mz <- list(c(100.0, 150.0), c(100.0, 200.0), c(250.0))
  prec_mz <- c(200.0, 250.0, 300.0)
  # Neutral loss of 50 Da matches spectrum 1 (200 - 50 = 150) and spectrum 2 (250 - 50 = 200)
  hits <- perform_query(
    all_mz,
    prec_mz,
    frags = 100.0,
    nls = 50.0,
    dalton = 0.01,
    ppm = 10
  )
  testthat::expect_identical(hits, c(TRUE, TRUE, FALSE))

  qd <- list(
    mz = all_mz,
    precursor_mz = prec_mz,
    skeleton = c("S1", "S2", "S1")
  )
  ions_list <- list(S1 = c("100_frag", "50_nl"), S2 = c("999_frag"))
  r1 <- perform_list_of_queries(1L, ions_list, qd, dalton = 0.01, ppm = 10)
  testthat::expect_identical(r1$value, c("S1", "S2"))
  r2 <- perform_list_of_queries(2L, ions_list, qd, dalton = 0.01, ppm = 10)
  testthat::expect_identical(nrow(r2), 0L)
})

testthat::test_that("spectra_to_queries reproduces exact raw and grouped results", {
  tmp_raw <- tempfile(fileext = ".tsv")
  tmp_grp <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_raw, tmp_grp)), add = TRUE)

  res_raw <- spectra_to_queries(export = tmp_raw)
  testthat::expect_true(file.exists(tmp_raw))
  testthat::expect_s3_class(res_raw, "data.frame")
  testthat::expect_identical(nrow(res_raw), 21L)
  testthat::expect_identical(
    colnames(res_raw),
    c("skeleton", "mcc", "combined_ions")
  )
  testthat::expect_identical(res_raw$skeleton[1L], "Lapidilectine")
  testthat::expect_equal(res_raw$mcc[1L], 1.0)

  res_grp <- spectra_to_queries(spectra = "grouped", export = tmp_grp)
  testthat::expect_true(file.exists(tmp_grp))
  testthat::expect_s3_class(res_grp, "data.frame")
  testthat::expect_identical(nrow(res_grp), 17L)
  testthat::expect_identical(res_grp$skeleton[1L], "Lapidilectine.Lunduran")
  testthat::expect_equal(res_grp$mcc[1L], 1.0)
})
