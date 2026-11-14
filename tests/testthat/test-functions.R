testthat::test_that("combine_ions_minimal simplifies and sorts expressions", {
  testthat::expect_identical(combine_ions_minimal(list()), "")
  testthat::expect_identical(combine_ions_minimal(list(character(0))), "")
  testthat::expect_identical(combine_ions_minimal(list(c("", NA))), "")
  testthat::expect_identical(combine_ions_minimal(list("not_an_ion")), "")

  # Single ion and numeric (not lexicographic) sorting of fragments then NLs
  testthat::expect_identical(combine_ions_minimal(list("100_frag")), "100_frag")
  testthat::expect_identical(combine_ions_minimal(list("50_nl")), "50_nl")
  testthat::expect_identical(
    combine_ions_minimal(list(c(
      "200_frag",
      "50_nl",
      "100_frag",
      "5_frag",
      "20_nl",
      "5_nl"
    ))),
    "5_frag & 100_frag & 200_frag & 5_nl & 20_nl & 50_nl"
  )

  # Strict suffix matching ignores malformed tokens without warnings
  testthat::expect_no_warning(
    testthat::expect_identical(
      combine_ions_minimal(list(
        "100_fragment",
        "abc_frag",
        "Inf_nl",
        "100_frag"
      )),
      "100_frag"
    )
  )

  # Common ions factored out with OR branches (both single and multi-ion)
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "200_frag"),
      c("100_frag", "50_nl", "60_nl")
    )),
    "100_frag & (200_frag | (50_nl & 60_nl))"
  )

  # Disjoint sets (no common ions)
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("200_frag", "100_frag"),
      c("50_nl")
    )),
    "(100_frag & 200_frag) | 50_nl"
  )

  # Subset query absorbs superset queries into common ions (kills any -> all)
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "200_frag"),
      c("100_frag"),
      c("100_frag", "300_frag")
    )),
    "100_frag"
  )
  # Multi-ion superset absorbed by smaller subset when other clauses remain
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "200_frag"),
      c("100_frag", "300_frag"),
      c("100_frag", "200_frag", "300_frag")
    )),
    "100_frag & (200_frag | 300_frag)"
  )

  # Duplicate ions within a query and order-variant remaining sets
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "100_frag", "200_frag"),
      c("100_frag", "300_frag")
    )),
    "100_frag & (200_frag | 300_frag)"
  )
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "200_frag"),
      c("300_frag", "400_frag"),
      c("400_frag", "300_frag")
    )),
    "(100_frag & 200_frag) | (300_frag & 400_frag)"
  )

  # to_massql translations, precedence, and tolerance qualifier rules
  testthat::expect_identical(to_massql(character(0)), character(0))
  testthat::expect_identical(
    to_massql(c("", NA_character_, "bad_ion", "100_fragment")),
    c("", "", "", "")
  )
  testthat::expect_identical(
    to_massql("100_frag & 50_nl", dalton = 0.01, ppm = 30, intensity_min = 5),
    paste0(
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=100:TOLERANCEMZ=0.013:INTENSITYPERCENT=5 AND ",
      "MS2NL=50:TOLERANCEMZ=0.01:INTENSITYPERCENT=5"
    )
  )
  testthat::expect_identical(
    to_massql("100_frag", dalton = 0.01, ppm = 0),
    "QUERY scaninfo(MS2DATA) WHERE MS2PROD=100:TOLERANCEMZ=0.01"
  )
  testthat::expect_identical(
    to_massql("100_frag & 50_nl", dalton = 0, ppm = 30),
    paste0(
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=100:TOLERANCEPPM=30 AND ",
      "MS2NL=50:TOLERANCEPPM=30"
    )
  )
  testthat::expect_identical(
    to_massql("100_frag", dalton = 0, ppm = 0),
    "QUERY scaninfo(MS2DATA) WHERE MS2PROD=100:TOLERANCEMZ=0"
  )
  testthat::expect_identical(
    to_massql(list(c("200_frag", "100_frag")), dalton = 0.02, ppm = 10),
    paste0(
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=100:TOLERANCEMZ=0.021 AND ",
      "MS2PROD=200:TOLERANCEMZ=0.022"
    )
  )
  testthat::expect_identical(
    to_massql(list("100_frag & 50_nl"), dalton = 0.01, ppm = 30),
    to_massql(list(c("100_frag", "50_nl")), dalton = 0.01, ppm = 30)
  )
  expected_or_prod <- paste0(
    "QUERY scaninfo(MS2DATA) WHERE ",
    "MS2PROD=100:TOLERANCEMZ=0.013 AND ",
    "MS2PROD=(200 OR 300):TOLERANCEMZ=0.019"
  )
  testthat::expect_identical(
    to_massql("100_frag & (200_frag | 300_frag)", dalton = 0.01, ppm = 30),
    expected_or_prod
  )
  testthat::expect_identical(
    to_massql("(200_frag | 300_frag) & 100_frag", dalton = 0.01, ppm = 30),
    expected_or_prod
  )
  testthat::expect_identical(
    to_massql("100_frag & (50_nl | 60_nl)", dalton = 0.01, ppm = 30),
    paste0(
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=100:TOLERANCEMZ=0.013 AND ",
      "MS2NL=(50 OR 60):TOLERANCEMZ=0.01"
    )
  )
  testthat::expect_identical(
    to_massql("100_frag & 200_frag | 300_frag", dalton = 0.01, ppm = 30),
    paste0(
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=100:TOLERANCEMZ=0.013 AND ",
      "MS2PROD=200:TOLERANCEMZ=0.016 ; ",
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=300:TOLERANCEMZ=0.019"
    )
  )
  testthat::expect_identical(
    to_massql(
      "(100_frag | 200_frag) & (300_frag | 400_frag)",
      dalton = 0.01,
      ppm = 0
    ),
    paste(
      "QUERY scaninfo(MS2DATA) WHERE",
      c(
        "MS2PROD=100:TOLERANCEMZ=0.01 AND MS2PROD=300:TOLERANCEMZ=0.01",
        "MS2PROD=100:TOLERANCEMZ=0.01 AND MS2PROD=400:TOLERANCEMZ=0.01",
        "MS2PROD=200:TOLERANCEMZ=0.01 AND MS2PROD=300:TOLERANCEMZ=0.01",
        "MS2PROD=200:TOLERANCEMZ=0.01 AND MS2PROD=400:TOLERANCEMZ=0.01"
      ),
      collapse = " ; "
    )
  )
  testthat::expect_identical(
    to_massql(
      "100_frag & ((200_frag & 50_nl) | 300_frag)",
      dalton = 0.01,
      ppm = 30
    ),
    paste0(
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=100:TOLERANCEMZ=0.013 AND ",
      "MS2PROD=200:TOLERANCEMZ=0.016 AND ",
      "MS2NL=50:TOLERANCEMZ=0.01 ; ",
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=100:TOLERANCEMZ=0.013 AND ",
      "MS2PROD=300:TOLERANCEMZ=0.019"
    )
  )
  testthat::expect_identical(
    to_massql("(100_frag & 200_frag) | 50_nl", dalton = 0.01, ppm = 30),
    paste0(
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2PROD=100:TOLERANCEMZ=0.013 AND ",
      "MS2PROD=200:TOLERANCEMZ=0.016 ; ",
      "QUERY scaninfo(MS2DATA) WHERE ",
      "MS2NL=50:TOLERANCEMZ=0.01"
    )
  )
  testthat::expect_error(
    to_massql("100_frag", dalton = -0.01),
    class = "stq_input_error"
  )
  testthat::expect_error(
    to_massql("100_frag", ppm = -1),
    class = "stq_input_error"
  )
  testthat::expect_error(
    to_massql("100_frag", intensity_min = -1),
    class = "stq_input_error"
  )
  testthat::expect_error(
    to_massql("100_frag", intensity_min = 101),
    class = "stq_input_error"
  )
  testthat::expect_error(to_massql("(100_frag"), class = "stq_input_error")
})

testthat::test_that("create_dir creates target directories", {
  tmp_base <- file.path(tempdir(), paste0("stq_test_", Sys.getpid()))
  on.exit(unlink(tmp_base, recursive = TRUE), add = TRUE)

  dir_target <- file.path(tmp_base, "a.b", "sub")
  testthat::expect_message(
    vis_dir <- withVisible(create_dir(dir_target)),
    "Directory .* created\\."
  )
  testthat::expect_false(vis_dir$visible)
  testthat::expect_identical(vis_dir$value, dir_target)
  testthat::expect_true(dir.exists(dir_target))

  # Second call on existing directory emits no message
  testthat::expect_no_message(
    testthat::expect_no_warning(create_dir(dir_target))
  )
})

testthat::test_that("filter_matrix and normalize_peaks handle boundaries", {
  mat <- matrix(
    c(1, 0, -5, 1, 2, 0, 1, 2, 3),
    nrow = 3L,
    dimnames = list(c("A", "B", "C"), c("100", "200", "300"))
  )
  # Column "100" has only 1 strictly positive value (0 and -5 are not > 0)
  filtered_2 <- filter_matrix(mat, n = 2L)
  testthat::expect_identical(colnames(filtered_2), c("200", "300"))
  testthat::expect_identical(dim(filtered_2), c(3L, 2L))

  filtered_3 <- filter_matrix(mat, n = 3L)
  testthat::expect_identical(colnames(filtered_3), "300")
  testthat::expect_identical(dim(filtered_3), c(3L, 1L))

  filtered_4 <- filter_matrix(mat, n = 4L)
  testthat::expect_identical(dim(filtered_4), c(3L, 0L))
  testthat::expect_identical(rownames(filtered_4), c("A", "B", "C"))

  norm_fn <- normalize_peaks()
  pks <- cbind(mz = c(100, 200), intensity = c(25, 50))
  res <- testthat::expect_no_warning(norm_fn(pks))
  testthat::expect_equal(res[, 1L], c(100, 200), tolerance = 1e-12)
  testthat::expect_equal(res[, 2L], c(50, 100), tolerance = 1e-12)
  # Idempotence: normalizing an already-normalized peak matrix is a no-op
  testthat::expect_identical(norm_fn(res), res)

  empty_pks <- matrix(numeric(0), nrow = 0L, ncol = 2L)
  testthat::expect_no_warning(
    testthat::expect_identical(norm_fn(empty_pks), empty_pks)
  )

  zero_pks <- cbind(mz = c(100, 200), intensity = c(0, 0))
  testthat::expect_identical(norm_fn(zero_pks), zero_pks)

  neg_pks <- cbind(mz = c(100, 200), intensity = c(-10, -5))
  testthat::expect_identical(norm_fn(neg_pks), neg_pks)
})

testthat::test_that("spectra conversion and matrix helpers round-trip", {
  utils::data(
    "mia_spectra_df",
    package = "SpectraToQueries",
    envir = environment()
  )
  sub_df <- mia_spectra_df[3:12, ]
  sp <- df_to_spectra(sub_df)
  rownames(sp@backend@spectraData) <- paste0("custom_", seq_along(sp))
  testthat::expect_s4_class(sp, "Spectra")
  testthat::expect_identical(length(sp), 10L)

  roundtrip_df <- spectra_to_df(sp)
  testthat::expect_identical(nrow(roundtrip_df), 10L)
  testthat::expect_identical(.row_names_info(roundtrip_df, type = 1L), -10L)
  testthat::expect_equal(roundtrip_df$mz, sub_df$mz)
  testthat::expect_equal(roundtrip_df$intensity, sub_df$intensity)

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

testthat::test_that("create_matrix handles empty, fast, and sparse paths", {
  # 0 spectra
  mat_0 <- create_matrix(Spectra::Spectra(), zero_val = 0, name = character(0))
  testthat::expect_identical(dim(mat_0), c(0L, 0L))

  # Spectra with 0 peaks
  df_no_pks <- data.frame(msLevel = c(2L, 2L), SKELETON = c("A", "B"))
  df_no_pks$mz <- list(numeric(0), numeric(0))
  df_no_pks$intensity <- list(numeric(0), numeric(0))
  mat_no_pks <- create_matrix(
    df_to_spectra(df_no_pks),
    zero_val = 0,
    name = c("A", "B")
  )
  testthat::expect_identical(dim(mat_no_pks), c(2L, 0L))
  testthat::expect_identical(rownames(mat_no_pks), c("A", "B"))

  # Mix of non-empty and empty spectra (all(lens == 0L) vs any(lens == 0L))
  df_one_empty <- data.frame(msLevel = c(2L, 2L), SKELETON = c("A", "B"))
  df_one_empty$mz <- list(c(100, 200), numeric(0))
  df_one_empty$intensity <- list(c(10, 20), numeric(0))
  mat_one_empty <- create_matrix(
    df_to_spectra(df_one_empty),
    zero_val = 0,
    name = c("A", "B")
  )
  testthat::expect_identical(
    mat_one_empty,
    matrix(
      c(10, 0, 20, 0),
      nrow = 2L,
      ncol = 2L,
      dimnames = list(c("A", "B"), c("100", "200"))
    )
  )

  # Fast path (shared sorted m/z grid) zero_val boundary
  df_same <- data.frame(msLevel = c(2L, 2L), SKELETON = c("S1", "S2"))
  df_same$mz <- list(c(100, 200), c(100, 200))
  df_same$intensity <- list(c(4, 10), c(6, 20))
  sp_same <- df_to_spectra(df_same)

  # Column "100" has sum 10; with zero_val = 10 it must be dropped (> not >=)
  mat_same_z10 <- create_matrix(sp_same, zero_val = 10, name = c("S1", "S2"))
  testthat::expect_identical(
    mat_same_z10,
    matrix(
      c(10, 20),
      nrow = 2L,
      ncol = 1L,
      dimnames = list(c("S1", "S2"), "200")
    )
  )
  # All columns <= zero_val in fast path
  mat_same_z50 <- create_matrix(sp_same, zero_val = 50, name = c("S1", "S2"))
  testthat::expect_identical(dim(mat_same_z50), c(2L, 0L))
  testthat::expect_identical(rownames(mat_same_z50), c("S1", "S2"))

  # Middle-spectrum grid difference (sparse path) + zero_val boundaries
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

  # Sparse path boundary: colSums are 90, 80, 40; zero_val = 40 drops "300"
  mat_mid_z40 <- create_matrix(
    sp_mid_diff,
    zero_val = 40,
    name = df_mid_diff$SKELETON
  )
  testthat::expect_identical(colnames(mat_mid_z40), c("100", "200"))

  # Sparse path all columns <= zero_val
  mat_mid_z100 <- create_matrix(
    sp_mid_diff,
    zero_val = 100,
    name = df_mid_diff$SKELETON
  )
  testthat::expect_identical(dim(mat_mid_z100), c(3L, 0L))
  testthat::expect_identical(rownames(mat_mid_z100), c("S1", "S2", "S3"))
})

testthat::test_that("harmonize_mzs and fix_binned_mzs respect tolerances", {
  # Spectra with 0 peaks across all spectra (n_avg == 0L)
  df_empty_pks <- data.frame(msLevel = c(2L, 2L), SKELETON = c("A", "B"))
  df_empty_pks$mz <- list(numeric(0), numeric(0))
  df_empty_pks$intensity <- list(numeric(0), numeric(0))
  sp_empty_pks <- df_to_spectra(df_empty_pks)
  testthat::expect_identical(
    Spectra::mz(harmonize_mzs(sp_empty_pks, dalton = 0.01, ppm = 30)),
    Spectra::mz(sp_empty_pks)
  )

  # Left vs right neighbor selection and exact binary dalton + ppm boundary
  df_h <- data.frame(msLevel = c(2L, 2L, 2L), SKELETON = c("A", "B", "C"))
  df_h$mz <- list(
    c(50.0, 99.875, 1000.000, 2000.000),
    c(50.0, 100.125, 1000.500, 2000.750),
    numeric(0)
  )
  df_h$intensity <- list(
    c(10, 20, 30, 40),
    c(10, 20, 30, 40),
    numeric(0)
  )
  sp_h <- harmonize_mzs(df_to_spectra(df_h), dalton = 0.25, ppm = 500)
  mzs_h <- as.list(Spectra::mz(sp_h))
  # With dalton = 0.25, ppm = 500:
  # - 99.875 and 100.125 average to 100.0 (diff 0.125 <= 0.30)
  # - 1000.000 and 1000.500 average to 1000.250 (diff 0.25 <= 0.750125)
  # - 2000.000 and 2000.750 average to 2000.375 (diff 0.375 <= 1.2501875)
  testthat::expect_equal(
    mzs_h[[1L]],
    c(50.0, 100.0, 1000.25, 2000.375),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    mzs_h[[2L]],
    c(50.0, 100.0, 1000.25, 2000.375),
    tolerance = 1e-10
  )
  testthat::expect_identical(mzs_h[[3L]], numeric(0))

  # fix_binned_mzs: 0-column input
  mat_0col <- matrix(
    0,
    nrow = 2L,
    ncol = 0L,
    dimnames = list(c("A", "B"), NULL)
  )
  testthat::expect_identical(
    fix_binned_mzs(mat_0col, sp_h, dalton = 0.01, ppm = 30, decimals = 2L),
    mat_0col
  )

  # fix_binned_mzs: fast path + exact binary tolerance boundary
  df_orig <- data.frame(msLevel = 2L, SKELETON = "A")
  df_orig$mz <- list(c(100.5, 200.125))
  df_orig$intensity <- list(c(50, 100))
  sp_orig <- df_to_spectra(df_orig)

  binned_unique <- matrix(
    c(10, 20, 30, 40),
    nrow = 2L,
    ncol = 2L,
    dimnames = list(c("S1", "S2"), c("100.00", "200.00"))
  )
  fixed_unique <- fix_binned_mzs(
    binned_unique,
    original_mzs = sp_orig,
    dalton = 0.5,
    ppm = 0,
    decimals = 3L
  )
  testthat::expect_identical(colnames(fixed_unique), c("100.5", "200.125"))
  testthat::expect_identical(as.vector(fixed_unique), c(10, 20, 30, 40))

  # fix_binned_mzs: column-merging path when two columns round to the same value
  binned_dup <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2L,
    ncol = 3L,
    dimnames = list(c("S1", "S2"), c("100.001", "100.004", "300.5"))
  )
  fixed_dup <- fix_binned_mzs(
    binned_dup,
    original_mzs = sp_orig,
    dalton = 0.01,
    ppm = 0,
    decimals = 2L
  )
  testthat::expect_identical(colnames(fixed_dup), c("100", "300.5"))
  testthat::expect_identical(
    fixed_dup,
    matrix(
      c(4, 6, 5, 6),
      nrow = 2L,
      ncol = 2L,
      dimnames = list(c("S1", "S2"), c("100", "300.5"))
    )
  )
})

testthat::test_that("query matching helpers work across all branches", {
  # prepare_query_data with non-finite m/z values and "+" in SKELETON
  df_qd <- data.frame(
    msLevel = c(2L, 2L),
    SKELETON = c("SkelA+SkelB", "SkelC"),
    precursorMz = c(200.0, 300.0)
  )
  df_qd$mz <- list(c(100.0, 150.0), c(Inf, 250.0))
  df_qd$intensity <- list(c(50.0, 100.0), c(50.0, 100.0))
  sp_qd <- df_to_spectra(df_qd)
  prep <- prepare_query_data(sp_qd)
  testthat::expect_identical(prep$mz, list(c(100.0, 150.0), 250.0))
  testthat::expect_identical(prep$precursor_mz, c(200.0, 300.0))
  testthat::expect_identical(prep$skeleton, c("SkelA.SkelB", "SkelC"))

  # prepare_query_data with intensity_min filtering and edge cases
  df_int <- data.frame(
    msLevel = c(2L, 2L, 2L, 2L),
    SKELETON = c("A", "B", "C", "D"),
    precursorMz = c(300.0, 300.0, 300.0, 300.0)
  )
  df_int$mz <- list(
    numeric(0),
    c(100.0, Inf),
    c(100.0, 200.0),
    c(100.0, 200.0)
  )
  df_int$intensity <- list(
    numeric(0),
    c(NaN, 100.0),
    c(0.0, 0.0),
    c(10.0, 100.0)
  )
  prep_int <- prepare_query_data(df_to_spectra(df_int), intensity_min = 50)
  testthat::expect_identical(
    prep_int$mz,
    list(numeric(0), numeric(0), numeric(0), 200.0)
  )

  # precompute_ion_hits: non-finite precursor, NL bounds, and malformed tokens
  all_mz <- list(c(100.0, 150.0), c(0.0, 100.0, 200.0), c(250.0), numeric(0))
  prec_mz <- c(Inf, 250.0, 300.0, 250.0)
  valid_ps <- is.finite(prec_mz) & (lengths(all_mz) > 0L)
  pre <- precompute_ion_hits(
    unique_ions = c(
      "100_frag",
      "50_nl",
      "300_nl",
      "250_nl",
      "bad_ion",
      "Inf_frag",
      "100_fragment"
    ),
    all_mz = all_mz,
    precursor_mz = prec_mz,
    valid_prec_and_spec = valid_ps,
    dalton = 0.01,
    ppm = 10
  )
  testthat::expect_identical(
    pre$names,
    c("100_frag", "50_nl", "300_nl", "250_nl")
  )
  testthat::expect_identical(
    pre$hits[["100_frag"]],
    c(TRUE, TRUE, FALSE, FALSE)
  )
  # 50_nl on spec 2 (250-50=200) & spec 3 (300-50=250); Inf prec skipped
  testthat::expect_identical(
    pre$hits[["50_nl"]],
    c(FALSE, TRUE, TRUE, FALSE)
  )
  # 300_nl > precursor_mz (250) and 250_nl == precursor_mz (0.0 peak) rejected
  testthat::expect_identical(pre$hits[["300_nl"]], rep(FALSE, 4L))
  testthat::expect_identical(pre$hits[["250_nl"]], rep(FALSE, 4L))
})

testthat::test_that("spectra_to_queries validates inputs and synthetic cases", {
  testthat::expect_error(
    spectra_to_queries(spectra = 123),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(spectra = ""),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(export = ""),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(dalton = -0.01),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(ppm = -1),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(decimals = -1L),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(ions_max = 0L),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(n_skel_min = 0L),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(n_spec_min = 0L),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(precision_min = 1.5),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(recall_min = -0.1),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(mcc_min = 2.0),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(intensity_min = Inf),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(intensity_min = -1),
    class = "stq_input_error"
  )
  testthat::expect_error(
    spectra_to_queries(intensity_min = 101),
    class = "stq_input_error"
  )

  # Valid boundary values in validate_stq_inputs must succeed without error
  testthat::expect_true(validate_stq_inputs(
    spectra = NULL,
    export = NULL,
    dalton = 0,
    decimals = 0L,
    intensity_min = 0,
    ions_max = 1L,
    n_skel_min = 1L,
    n_spec_min = 1L,
    ppm = 0,
    mcc_min = -1,
    precision_min = 0,
    recall_min = 0,
    zero_val = 0
  ))
  testthat::expect_true(validate_stq_inputs(
    spectra = "grouped",
    export = "out.tsv",
    dalton = 0.01,
    decimals = 4L,
    intensity_min = 100,
    ions_max = 2L,
    n_skel_min = 2L,
    n_spec_min = 2L,
    ppm = 30,
    mcc_min = 1,
    precision_min = 1,
    recall_min = 1,
    zero_val = 0
  ))

  # Startup message coverage
  testthat::expect_message(
    .onAttach(NULL, "SpectraToQueries"),
    "Welcome to SpectraToQueries"
  )

  # Synthetic 6-spectrum dataset with distinct PRECURSOR_MZ and TITLEs:
  # SkelA: 100_frag (3/3, FP=0 -> MCC = 1.0)
  # SkelB: 150_frag (TP=2, FP=1, FN=1, TN=2) -> MCC = (4-1)/9 = 0.3333
  synth_df <- data.frame(
    msLevel = rep(2L, 6L),
    rtime = seq_len(6L) * 1.0,
    acquisitionNum = seq_len(6L),
    TITLE = paste0("Spec_", seq_len(6L)),
    PRECURSOR_MZ = c(300.0, 400.0, 500.0, 600.0, 700.0, 800.0),
    SKELETON = c("SkelA", "SkelA", "SkelA", "SkelB", "SkelB", "SkelB"),
    stringsAsFactors = FALSE
  )
  synth_df$mz <- list(
    c(100.0, 120.0),
    c(100.0, 120.0),
    c(100.0, 150.0),
    c(150.0, 160.0),
    c(150.0, 170.0),
    c(180.0, 190.0)
  )
  synth_df$intensity <- rep(list(c(100.0, 80.0)), 6L)

  synth_sp <- df_to_spectra(synth_df)
  res_sp <- spectra_to_queries(
    spectra = synth_sp,
    export = NULL,
    n_skel_min = 2L,
    n_spec_min = 2L
  )
  testthat::expect_identical(nrow(res_sp), 2L)
  testthat::expect_equal(
    res_sp$mcc,
    c(1.0, round(1 / 3, 4L)),
    tolerance = 1e-12
  )

  # Write synthetic dataset to temporary MGF file and load via file path branch
  tmp_mgf <- tempfile(fileext = ".mgf")
  tmp_empty_tsv <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_mgf, tmp_empty_tsv)), add = TRUE)
  mgf_lines <- unlist(lapply(seq_len(nrow(synth_df)), function(i) {
    c(
      "BEGIN IONS",
      paste0("TITLE=", synth_df$TITLE[i]),
      paste0("PEPMASS=", synth_df$PRECURSOR_MZ[i]),
      paste0("PRECURSOR_MZ=", synth_df$PRECURSOR_MZ[i]),
      "MSLEVEL=2",
      paste0("SKELETON=", synth_df$SKELETON[i]),
      paste(synth_df$mz[[i]], synth_df$intensity[[i]]),
      "END IONS",
      ""
    )
  }))
  writeLines(mgf_lines, tmp_mgf)
  res_mgf <- spectra_to_queries(
    spectra = tmp_mgf,
    export = NULL,
    n_skel_min = 2L,
    n_spec_min = 2L
  )
  testthat::expect_identical(res_mgf, res_sp)

  # All columns filtered by high zero_val + empty export file written
  res_zero_filt <- spectra_to_queries(
    spectra = synth_df,
    export = tmp_empty_tsv,
    n_skel_min = 2L,
    n_spec_min = 2L,
    zero_val = 1e6
  )
  testthat::expect_identical(nrow(res_zero_filt), 0L)
  testthat::expect_true(file.exists(tmp_empty_tsv))
})

testthat::test_that("pipeline functions handle empty and filtered-out inputs", {
  sp_empty <- Spectra::Spectra()
  testthat::expect_s4_class(
    harmonize_mzs(sp_empty, dalton = 0.01, ppm = 30),
    "Spectra"
  )
  mat_one <- matrix(
    1,
    nrow = 1L,
    ncol = 1L,
    dimnames = list("S1", "100.12345")
  )
  fixed_empty <- fix_binned_mzs(
    mat_one,
    original_mzs = sp_empty,
    dalton = 0.01,
    ppm = 30,
    decimals = 2L
  )
  testthat::expect_identical(colnames(fixed_empty), "100.12")

  utils::data(
    "mia_spectra_df",
    package = "SpectraToQueries",
    envir = environment()
  )
  sub_df <- mia_spectra_df[1:10, ]
  res_no_skel <- spectra_to_queries(
    sub_df,
    export = NULL,
    n_skel_min = 100L
  )
  testthat::expect_identical(
    res_no_skel,
    data.frame(
      skeleton = character(0),
      mcc = numeric(0),
      combined_ions = character(0),
      massql = character(0),
      stringsAsFactors = FALSE
    )
  )

  df_zero_int <- data.frame(
    msLevel = c(2L, 2L),
    TITLE = c("Z1", "Z2"),
    PRECURSOR_MZ = c(300.0, 400.0),
    SKELETON = c("A", "A"),
    stringsAsFactors = FALSE
  )
  df_zero_int$mz <- list(c(100.0, 120.0), c(100.0, 120.0))
  df_zero_int$intensity <- list(c(0.0, 0.0), c(0.0, 0.0))
  res_no_peaks <- spectra_to_queries(
    df_zero_int,
    export = NULL,
    intensity_min = 50,
    n_skel_min = 2L,
    n_spec_min = 2L
  )
  testthat::expect_identical(
    res_no_peaks,
    data.frame(
      skeleton = character(0),
      mcc = numeric(0),
      combined_ions = character(0),
      massql = character(0),
      stringsAsFactors = FALSE
    )
  )
})

testthat::test_that("extracted pipeline helpers handle boundaries and ties", {
  # parse_query_ions
  parsed <- parse_query_ions(c("100.5_frag", "bad_frag", "50.25_nl", "Inf_nl"))
  testthat::expect_identical(parsed$frags, 100.5)
  testthat::expect_identical(parsed$nls, 50.25)
  parsed_empty <- parse_query_ions(character(0))
  testthat::expect_identical(parsed_empty$frags, numeric(0))
  testthat::expect_identical(parsed_empty$nls, numeric(0))

  # extract_candidate_ions:
  # 4 spectra total: 2 in G+1 ("G.1"), 2 in G2
  # - "100_frag": in all 4 spectra -> denom_sq == 0 (tn + fn_ == 0) -> mcc = 0
  #   for G.1: precision = 0.5, recall = 1.0 (diagnostic!)
  # - "200_frag": in 1 G.1 spectrum and 1 G2 spectrum -> precision = 0.5,
  #   recall = 0.5, tp=1, fp=1, fn=1, tn=1 -> mcc = 0.0
  # - "300_frag": in 1 G.1 spectrum and 1 G2 spectrum -> mcc = 0.0
  mat_cand <- matrix(
    c(
      10,
      10,
      10,
      10,
      10,
      0,
      10,
      0,
      0,
      10,
      0,
      10
    ),
    nrow = 4L,
    ncol = 3L,
    dimnames = list(
      c("G+1", "G+1", "G2", "G2"),
      c("100_frag", "200_frag", "300_frag")
    )
  )
  cand <- suppressMessages(extract_candidate_ions(
    merged_mat = mat_cand,
    total_spectra = 4,
    n_skel_min = 2L,
    n_spec_min = 2L,
    ions_max = 1L,
    precision_min = 0.5,
    recall_min = 0.5,
    mcc_min = 0.0
  ))
  # ions_max = 1L slices candidate ions per group to 1 ion
  testthat::expect_identical(lengths(cand), c(G.1 = 1L, G2 = 1L))

  # Strict precision_min / recall_min filtering (> 0.5 drops 0.5)
  cand_filtered <- suppressMessages(extract_candidate_ions(
    merged_mat = mat_cand,
    total_spectra = 4,
    n_skel_min = 2L,
    n_spec_min = 2L,
    ions_max = 5L,
    precision_min = 0.6,
    recall_min = 0.5,
    mcc_min = 0.0
  ))
  testthat::expect_identical(cand_filtered, list())

  # prune_combinations_bnb: empty list, zero spectra, and no FP drop
  bnb_empty <- prune_combinations_bnb(
    ions_list = list(),
    skeleton = c("G1", "G2"),
    hits_lookup = list("100_frag" = c(TRUE, FALSE)),
    max_ions = 2L,
    decimals = 4L
  )
  testthat::expect_identical(bnb_empty$ions, list())

  bnb_res <- prune_combinations_bnb(
    ions_list = list(Empty = character(0), G1 = c("100_frag", "200_frag")),
    skeleton = c("G1", "G1", "G2"),
    hits_lookup = list(
      "100_frag" = c(TRUE, TRUE, TRUE),
      "200_frag" = c(TRUE, TRUE, TRUE)
    ),
    max_ions = 2L,
    decimals = 4L
  )
  testthat::expect_identical(bnb_res$skeleton, c("G1", "G1"))
  testthat::expect_identical(bnb_res$ions, list("100_frag", "200_frag"))
  testthat::expect_equal(bnb_res$tp, c(2, 2))
  testthat::expect_equal(bnb_res$fp, c(1, 1))

  # summarize_best_queries: empty input + minimal best query selection + massql
  empty_sum <- summarize_best_queries(
    candidates = bnb_empty,
    skeleton_counts = table(c("G1", "G2")),
    total_spectra = 2,
    decimals = 4L
  )
  testthat::expect_identical(nrow(empty_sum), 0L)

  cands_in <- list(
    skeleton = c("SkelA", "SkelA", "SkelB"),
    ions = list(
      c("100_frag", "200_frag"),
      c("100_frag", "300_frag"),
      "400_frag"
    ),
    tp = c(2, 2, 0),
    fp = c(0, 0, 0)
  )
  skel_counts <- table(c("SkelA", "SkelA", "SkelB", "SkelB"))
  summary_df <- summarize_best_queries(
    candidates = cands_in,
    skeleton_counts = skel_counts,
    total_spectra = 4,
    decimals = 4L
  )
  testthat::expect_identical(summary_df$skeleton, c("SkelA", "SkelB"))
  testthat::expect_equal(summary_df$mcc, c(1.0, 0.0))
  testthat::expect_identical(
    summary_df$combined_ions,
    c("100_frag & 200_frag", "400_frag")
  )
  testthat::expect_identical(
    summary_df$massql,
    to_massql(c("100_frag & 200_frag", "400_frag"))
  )
})

testthat::test_that("spectra_to_queries reproduces raw and grouped results", {
  tmp_raw <- tempfile(fileext = ".tsv")
  tmp_grp <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_raw, tmp_grp)), add = TRUE)

  res_raw <- spectra_to_queries(export = tmp_raw)
  testthat::expect_true(file.exists(tmp_raw))
  testthat::expect_s3_class(res_raw, "data.frame")
  testthat::expect_identical(nrow(res_raw), 21L)
  testthat::expect_identical(
    colnames(res_raw),
    c("skeleton", "mcc", "combined_ions", "massql")
  )
  testthat::expect_identical(res_raw$skeleton[1L], "Lapidilectine")
  testthat::expect_equal(res_raw$mcc[1L], 1.0)

  res_grp <- spectra_to_queries(spectra = "grouped", export = tmp_grp)
  testthat::expect_true(file.exists(tmp_grp))
  testthat::expect_s3_class(res_grp, "data.frame")
  testthat::expect_identical(nrow(res_grp), 17L)
  testthat::expect_identical(res_grp$skeleton[1L], "Lapidilectine.Lunduran")
  testthat::expect_equal(res_grp$mcc[1L], 1.0)

  tsv_raw_path <- testthat::test_path(
    "..",
    "..",
    "data",
    "interim",
    "queries.tsv"
  )
  if (file.exists(tsv_raw_path)) {
    tsv_raw <- utils::read.delim(tsv_raw_path, stringsAsFactors = FALSE)
    testthat::expect_equal(res_raw, tsv_raw, tolerance = 1e-12)
  }
  tsv_grp_path <- testthat::test_path(
    "..",
    "..",
    "data",
    "interim",
    "queries-grouped.tsv"
  )
  if (file.exists(tsv_grp_path)) {
    tsv_grp <- utils::read.delim(tsv_grp_path, stringsAsFactors = FALSE)
    testthat::expect_equal(res_grp, tsv_grp, tolerance = 1e-12)
  }
})

testthat::test_that("pipeline handles numeric skeletons and backends", {
  # 1. Skeletons ending in digits must not be stripped by unlist regex
  df_num <- data.frame(
    msLevel = rep(2L, 6L),
    TITLE = paste0("T", 1:6),
    PRECURSOR_MZ = rep(300.0, 6L),
    SKELETON = c("Skel_1", "Skel_1", "Skel_1", "Skel_2", "Skel_2", "Skel_2"),
    stringsAsFactors = FALSE
  )
  df_num$mz <- list(
    c(100.0, 120.0),
    c(100.0, 120.0),
    c(100.0, 150.0),
    c(200.0, 210.0),
    c(200.0, 220.0),
    c(200.0, 230.0)
  )
  df_num$intensity <- rep(list(c(100.0, 80.0)), 6L)

  res_num <- suppressMessages(spectra_to_queries(
    df_num,
    export = NULL,
    n_skel_min = 2L,
    n_spec_min = 2L
  ))
  testthat::expect_identical(sort(res_num$skeleton), c("Skel_1", "Skel_2"))
  testthat::expect_equal(res_num$mcc, c(1.0, 1.0), tolerance = 1e-12)

  # 2. Lowercase skeleton accepted; missing SKELETON/skeleton throws error
  df_low <- df_num
  names(df_low)[names(df_low) == "SKELETON"] <- "skeleton"
  res_low <- suppressMessages(spectra_to_queries(
    df_low,
    export = NULL,
    n_skel_min = 2L,
    n_spec_min = 2L
  ))
  testthat::expect_identical(res_low, res_num)

  df_no_skel_col <- df_num[, setdiff(names(df_num), "SKELETON"), drop = FALSE]
  testthat::expect_error(
    spectra_to_queries(df_no_skel_col, export = NULL),
    class = "stq_input_error"
  )

  # 3. Missing MGF path raises stq_file_not_found_error
  testthat::expect_error(
    spectra_to_queries(
      spectra = file.path(tempdir(), "does_not_exist_12345.mgf"),
      export = NULL
    ),
    class = "stq_file_not_found_error"
  )

  # 4. MsBackendDataFrame works in harmonize_mzs() and spectra_to_queries()
  sp_df_backend <- Spectra::Spectra(
    Spectra::backendInitialize(
      Spectra::MsBackendDataFrame(),
      data = df_num
    )
  )
  sp_df_harm <- harmonize_mzs(sp_df_backend, dalton = 0.01, ppm = 30)
  testthat::expect_s4_class(sp_df_harm, "Spectra")
  res_df_backend <- suppressMessages(spectra_to_queries(
    sp_df_backend,
    export = NULL,
    n_skel_min = 2L,
    n_spec_min = 2L
  ))
  testthat::expect_identical(res_df_backend, res_num)

  # 5. Extensionless export path with dotted parent directory writes file
  tmp_dir <- file.path(tempdir(), paste0("stq_noext_", Sys.getpid()), "a.b")
  on.exit(unlink(dirname(tmp_dir), recursive = TRUE), add = TRUE)
  out_noext <- file.path(tmp_dir, "queries_out")
  suppressMessages(spectra_to_queries(
    df_num,
    export = out_noext,
    n_skel_min = 2L,
    n_spec_min = 2L
  ))
  testthat::expect_true(file.exists(out_noext))
  testthat::expect_false(dir.exists(out_noext))

  # 6. intensity_min filters sub-threshold peaks during query evaluation
  df_int <- data.frame(
    msLevel = rep(2L, 4L),
    TITLE = paste0("I", 1:4),
    PRECURSOR_MZ = rep(300.0, 4L),
    SKELETON = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  df_int$mz <- list(
    c(100.0, 150.0),
    c(100.0, 160.0),
    c(100.0, 200.0),
    c(100.0, 200.0)
  )
  df_int$intensity <- list(
    c(100.0, 50.0),
    c(100.0, 50.0),
    c(1.0, 100.0),
    c(1.0, 100.0)
  )
  res_int <- suppressMessages(spectra_to_queries(
    df_int,
    export = NULL,
    intensity_min = 50,
    n_skel_min = 2L,
    n_spec_min = 2L
  ))
  testthat::expect_identical(sort(res_int$skeleton), c("A", "B"))
  testthat::expect_equal(res_int$mcc, c(1.0, 1.0), tolerance = 1e-12)

  # 7. msLevel != 2 and NA/empty SKELETON excluded before total_spectra
  df_mixed <- data.frame(
    msLevel = c(1L, 2L, 2L, 2L, 2L, 2L, 2L),
    TITLE = paste0("M", 1:7),
    PRECURSOR_MZ = rep(300.0, 7L),
    SKELETON = c("A", NA_character_, "", "A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  df_mixed$mz <- list(
    c(100.0, 200.0),
    c(100.0, 200.0),
    c(100.0, 200.0),
    c(100.0, 150.0),
    c(100.0, 160.0),
    c(200.0, 250.0),
    c(200.0, 260.0)
  )
  df_mixed$intensity <- rep(list(c(100.0, 80.0)), 7L)
  res_mixed <- suppressMessages(spectra_to_queries(
    df_mixed,
    export = NULL,
    n_skel_min = 2L,
    n_spec_min = 2L
  ))
  testthat::expect_equal(res_mixed$mcc, c(1.0, 1.0), tolerance = 1e-12)

  # All spectra filtered out by msLevel/SKELETON returns empty result cleanly
  df_all_bad <- df_mixed[1:3, ]
  out_empty_ms1 <- file.path(tmp_dir, "empty_ms1.tsv")
  res_all_bad <- suppressMessages(spectra_to_queries(
    df_all_bad,
    export = out_empty_ms1,
    n_skel_min = 2L,
    n_spec_min = 2L
  ))
  testthat::expect_identical(nrow(res_all_bad), 0L)
  testthat::expect_true(file.exists(out_empty_ms1))

  # Collision warning when '+' -> '.' maps distinct skeletons together
  df_col <- df_num
  df_col$SKELETON <- c("A+B", "A+B", "A+B", "A.B", "A.B", "A.B")
  testthat::expect_warning(
    suppressMessages(spectra_to_queries(
      df_col,
      export = NULL,
      n_skel_min = 2L,
      n_spec_min = 2L
    )),
    "collide"
  )

  # 8. create_matrix accumulates duplicate m/z intensities in sparse path
  df_dup_mz <- data.frame(
    msLevel = c(2L, 2L),
    SKELETON = c("A", "B"),
    stringsAsFactors = FALSE
  )
  df_dup_mz$mz <- list(c(100.0, 100.0, 200.0), c(100.0, 300.0))
  df_dup_mz$intensity <- list(c(10.0, 25.0, 50.0), c(5.0, 10.0))
  mat_dup_mz <- create_matrix(
    df_to_spectra(df_dup_mz),
    zero_val = 0,
    name = c("A", "B")
  )
  testthat::expect_identical(mat_dup_mz["A", "100"], 35.0)
})

testthat::test_that("spectra_to_queries works when TITLE is absent", {
  df_no_title <- data.frame(
    msLevel = c(2L, 2L, 2L, 2L),
    SKELETON = c("SkelA", "SkelA", "SkelB", "SkelB"),
    precursorMz = c(300.0, 300.0, 300.0, 300.0),
    stringsAsFactors = FALSE
  )
  df_no_title$mz <- list(
    c(100.0, 150.0),
    c(100.0, 160.0),
    c(200.0, 250.0),
    c(200.0, 260.0)
  )
  df_no_title$intensity <- list(
    c(80, 100),
    c(70, 100),
    c(90, 100),
    c(85, 100)
  )
  res <- suppressMessages(spectra_to_queries(
    spectra = df_no_title,
    export = NULL,
    n_skel_min = 2L,
    n_spec_min = 2L
  ))
  testthat::expect_identical(res$skeleton, c("SkelA", "SkelB"))
  testthat::expect_equal(res$mcc, c(1.0, 1.0), tolerance = 1e-12)
  testthat::expect_identical(res$combined_ions, c("100_frag", "200_frag"))
})
