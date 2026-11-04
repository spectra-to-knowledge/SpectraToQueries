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

  # Common ions factored out with OR branches (both single and multi-ion)
  testthat::expect_identical(
    combine_ions_minimal(list(
      c("100_frag", "200_frag"),
      c("100_frag", "50_nl", "60_nl")
    )),
    "100_frag & 200_frag | (50_nl & 60_nl)"
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

  # Duplicate ions within a query and order-variant remaining sets
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
  testthat::expect_message(
    vis_file <- withVisible(create_dir(file_target)),
    "Directory .* created\\."
  )
  testthat::expect_false(vis_file$visible)
  testthat::expect_identical(vis_file$value, dirname(file_target))
  testthat::expect_true(dir.exists(dirname(file_target)))

  # Second call on existing directory emits no message
  testthat::expect_no_message(
    testthat::expect_no_warning(create_dir(file_target))
  )

  dir_target <- file.path(tmp_base, "sub2")
  testthat::expect_message(
    vis_dir <- withVisible(create_dir(dir_target)),
    "Directory .* created\\."
  )
  testthat::expect_false(vis_dir$visible)
  testthat::expect_identical(vis_dir$value, dir_target)
  testthat::expect_true(dir.exists(dir_target))
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

testthat::test_that("generate_combinations and show_progress work", {
  combs <- generate_combinations(c("a", "b", "c"), max_ions = 2L)
  testthat::expect_identical(
    combs,
    list("a", "b", "c", c("a", "b"), c("a", "c"), c("b", "c"))
  )
  testthat::expect_identical(
    generate_combinations(c("a", "b"), max_ions = 5L),
    list("a", "b", c("a", "b"))
  )
  testthat::expect_identical(
    generate_combinations(c("a", "b"), max_ions = 1L),
    list("a", "b")
  )
  testthat::expect_identical(generate_combinations(character(0), 2L), list())
  testthat::expect_identical(generate_combinations(c("a", "b"), 0L), list())
  testthat::expect_identical(generate_combinations(c("a", "b"), -2L), list())

  ions_list <- list(s1 = c("a", "b"), s2 = c("c"))
  # Reordered indices verify `indices[idx]` lookup rather than positional `idx`
  res <- testthat::expect_no_message(
    generate_combinations_progress(
      c("s2", "s1"),
      ions_list = ions_list,
      max_ions = 2L,
      show_progress = FALSE
    )
  )
  testthat::expect_identical(
    res,
    list(list("c"), list("a", "b", c("a", "b")))
  )

  testthat::expect_message(
    generate_combinations_progress(
      "s2",
      ions_list = ions_list,
      max_ions = 1L,
      show_progress = TRUE
    ),
    "Progress: 1/1 \\(100%\\) - Done!"
  )

  testthat::expect_no_message(
    vis0 <- withVisible(show_progress(0L, 0L))
  )
  testthat::expect_false(vis0$visible)
  testthat::expect_null(vis0$value)
  testthat::expect_no_message(show_progress(1L, -5L))
  testthat::expect_no_message(show_progress(9L, 100L))
  testthat::expect_message(
    show_progress(10L, 100L),
    "Progress: 10/100 \\(10%\\)"
  )
  testthat::expect_no_message(show_progress(11L, 100L))
  testthat::expect_no_message(show_progress(10L, 9L))
  testthat::expect_message(
    vis_done <- withVisible(show_progress(10L, 10L)),
    "Progress: 10/10 \\(100%\\) - Done!"
  )
  testthat::expect_false(vis_done$visible)
  testthat::expect_null(vis_done$value)
})

testthat::test_that("spectra conversion and matrix helpers round-trip", {
  utils::data(
    "mia_spectra_df",
    package = "SpectraToQueries",
    envir = environment()
  )
  sub_df <- mia_spectra_df[3:12, ]
  sp <- df_to_spectra(sub_df)
  testthat::expect_s4_class(sp, "Spectra")
  testthat::expect_identical(length(sp), 10L)

  roundtrip_df <- spectra_to_df(sp)
  testthat::expect_identical(nrow(roundtrip_df), 10L)
  testthat::expect_identical(rownames(roundtrip_df), as.character(1:10))
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

  # Left vs right neighbor selection and exact dalton boundary in harmonize_mzs
  df_h <- data.frame(msLevel = c(2L, 2L, 2L), SKELETON = c("A", "B", "C"))
  df_h$mz <- list(
    c(50.0, 99.996, 1000.000, 2000.000),
    c(50.0, 100.004, 1000.020, 2000.030),
    numeric(0)
  )
  df_h$intensity <- list(
    c(10, 20, 30, 40),
    c(10, 20, 30, 40),
    numeric(0)
  )
  sp_h <- harmonize_mzs(df_to_spectra(df_h), dalton = 0.01, ppm = 30)
  mzs_h <- as.list(Spectra::mz(sp_h))
  # 99.996 and 100.004 average to 100.0 (diff 0.004 <= 0.01)
  # 1000.000 and 1000.020 average to 1000.010 (diff 0.01 == 0.01, within!)
  # 2000.000 and 2000.030 average to 2000.015 (diff 0.015 > 0.01, not within!)
  testthat::expect_equal(
    mzs_h[[1L]],
    c(50.0, 100.0, 1000.01, 2000.000),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    mzs_h[[2L]],
    c(50.0, 100.0, 1000.01, 2000.030),
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

  # fix_binned_mzs: fast path (all unique rounded names)
  df_orig <- data.frame(msLevel = 2L, SKELETON = "A")
  df_orig$mz <- list(c(100.002, 200.004))
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
    dalton = 0.01,
    ppm = 0,
    decimals = 3L
  )
  testthat::expect_identical(colnames(fixed_unique), c("100.002", "200.004"))
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
    PRECURSOR_MZ = c(200.0, 300.0)
  )
  df_qd$mz <- list(c(100.0, 150.0), c(100.0, 250.0))
  df_qd$intensity <- list(c(50.0, 100.0), c(50.0, 100.0))
  sp_qd <- df_to_spectra(df_qd)
  sp_qd@backend@spectraData$precursorMz <- df_qd$PRECURSOR_MZ
  sp_qd@backend@peaksData[[2L]][1L, 1L] <- Inf
  prep <- prepare_query_data(sp_qd)
  testthat::expect_identical(prep$mz, list(c(100.0, 150.0), 250.0))
  testthat::expect_identical(prep$precursor_mz, c(200.0, 300.0))
  testthat::expect_identical(prep$skeleton, c("SkelA.SkelB", "SkelC"))

  # contains_all_mz_tol & contains_all_mz exact binary boundary and non-finite
  testthat::expect_true(contains_all_mz_tol(c(100.0, Inf), 100.5, tol = 0.5))
  testthat::expect_false(contains_all_mz_tol(c(100.0, Inf), 100.6, tol = 0.5))
  testthat::expect_false(contains_all_mz_tol(numeric(0), 100.0, tol = 0.5))
  testthat::expect_false(
    contains_all_mz_tol(100.0, numeric(0), tol = numeric(0))
  )

  # At target_mz = 1000, dalton = 0.01, ppm = 20 -> tol = 0.03
  testthat::expect_true(
    contains_all_mz(1000.025, 1000.0, dalton = 0.01, ppm = 20)
  )
  testthat::expect_false(
    contains_all_mz(1000.025, 1000.0, dalton = 0.01, ppm = 10)
  )
  testthat::expect_false(
    contains_all_mz(100.0, c(NA_real_, Inf), dalton = 0.01, ppm = 10)
  )
  testthat::expect_false(contains_all_mz(
    numeric(0),
    100.0,
    dalton = 0.01,
    ppm = 10
  ))

  # perform_query edge cases: 0 spectra, S4 NumericList, empty spec in NL loop
  testthat::expect_identical(
    perform_query(list(), numeric(0), 100, numeric(0), 0.01, 10),
    logical(0)
  )
  testthat::expect_identical(
    perform_query(
      Spectra::mz(sp_qd),
      c(200.0, 300.0),
      150.0,
      numeric(0),
      0.01,
      10
    ),
    c(TRUE, FALSE)
  )

  all_mz <- list(c(100.0, 150.0), c(100.0, 200.0), c(250.0), numeric(0))
  prec_mz <- c(200.0, 250.0, 300.0, 250.0)
  # NL of 50 Da + 100_frag matches spectrum 1 (200-50=150) and 2 (250-50=200)
  hits <- perform_query(
    all_mz,
    prec_mz,
    frags = 100.0,
    nls = 50.0,
    dalton = 0.01,
    ppm = 10
  )
  testthat::expect_identical(hits, c(TRUE, TRUE, FALSE, FALSE))

  # NL-only query where spec 4 is empty and spec 1 has precursor == nl
  hits_nl_only <- perform_query(
    all_mz,
    prec_mz,
    frags = numeric(0),
    nls = 200.0,
    dalton = 0.01,
    ppm = 10
  )
  testthat::expect_identical(hits_nl_only, c(FALSE, FALSE, FALSE, FALSE))

  qd <- list(
    mz = all_mz,
    precursor_mz = prec_mz,
    skeleton = c("S1", "S2", "S1", NA_character_)
  )
  ions_list <- list(
    S1 = c("100_frag", "50_nl"),
    S2 = c("999_frag"),
    c("50_nl"),
    Empty = character(0),
    Invalid = c("bad_frag", "bad_nl"),
    OnlyNA = c("250_frag")
  )
  r1 <- perform_list_of_queries(1L, ions_list, qd, dalton = 0.01, ppm = 10)
  testthat::expect_identical(r1$target, c("S1", "S1"))
  testthat::expect_identical(r1$value, c("S1", "S2"))
  r2 <- perform_list_of_queries(2L, ions_list, qd, dalton = 0.01, ppm = 10)
  testthat::expect_identical(nrow(r2), 0L)
  # Unnamed entry defaults to "query_3" and matches specs 1, 2, 3 via 50_nl
  r3 <- perform_list_of_queries(3L, ions_list, qd, dalton = 0.01, ppm = 10)
  testthat::expect_identical(r3$target, c("query_3", "query_3", "query_3"))
  testthat::expect_identical(r3$value, c("S1", "S2", "S1"))
  # Empty and invalid ion entries return 0-row data.frame
  testthat::expect_identical(
    nrow(perform_list_of_queries(4L, ions_list, qd, dalton = 0.01, ppm = 10)),
    0L
  )
  testthat::expect_identical(
    nrow(perform_list_of_queries(5L, ions_list, qd, dalton = 0.01, ppm = 10)),
    0L
  )
  # Entry matching only a spectrum whose skeleton is NA returns 0-row data.frame
  qd_na <- qd
  qd_na$skeleton <- rep(NA_character_, 4L)
  testthat::expect_identical(
    nrow(perform_list_of_queries(
      6L,
      ions_list,
      qd_na,
      dalton = 0.01,
      ppm = 10
    )),
    0L
  )
  # Error inside perform_query is caught and warns
  testthat::expect_warning(
    r_err <- perform_list_of_queries(
      1L,
      ions_list,
      qd,
      dalton = "invalid_dalton",
      ppm = 10
    ),
    "Error in perform_query for S1"
  )
  testthat::expect_identical(nrow(r_err), 0L)

  # Non-finite target_mz or precursor_mz must not produce false-positive matches
  testthat::expect_false(contains_all_mz(
    c(100.0, 200.0),
    Inf,
    dalton = 0.01,
    ppm = 10
  ))
  df_inf <- data.frame(
    msLevel = c(2L, 2L, 2L),
    SKELETON = c("S1", "S2", NA_character_),
    PRECURSOR_MZ = c(Inf, 250.0, 300.0)
  )
  df_inf$mz <- list(c(100.0, 150.0), c(100.0, 200.0), c(250.0))
  df_inf$intensity <- list(c(50.0, 100.0), c(50.0, 100.0), c(100.0))
  sp_inf <- df_to_spectra(df_inf)
  sp_inf@backend@spectraData$precursorMz <- df_inf$PRECURSOR_MZ

  testthat::expect_identical(
    perform_list_of_queries_progress(
      list(),
      spectra = sp_inf,
      dalton = 0.01,
      ppm = 10,
      show_progress = FALSE
    ),
    list()
  )

  res_inf <- perform_list_of_queries_progress(
    list(
      Q1 = "50_nl",
      c("50_nl", "300_nl"),
      Invalid = "bad_ion",
      OnlyNA = "250_frag"
    ),
    spectra = sp_inf,
    dalton = 0.01,
    ppm = 10,
    show_progress = FALSE
  )
  testthat::expect_identical(res_inf[[1L]]$value, "S2")
  # Multi-NL query: 300_nl > precursor_mz (250) ignored while 50_nl matches S2
  testthat::expect_identical(res_inf[[2L]]$target, "query_2")
  testthat::expect_identical(res_inf[[2L]]$value, "S2")
  testthat::expect_identical(nrow(res_inf[[3L]]), 0L)
  testthat::expect_identical(nrow(res_inf[[4L]]), 0L)
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

  # Startup message coverage
  testthat::expect_message(
    SpectraToQueries:::.onAttach(NULL, "SpectraToQueries"),
    "Welcome to SpectraToQueries"
  )

  # Synthetic 6-spectrum dataset with sorted unique TITLEs:
  # 3 in SkelA (all have 100_frag, 2 have 120_frag; precursor 300)
  # 3 in SkelB (all have 150_frag, 2 have 160_frag; precursor 300)
  # Combined query (100_frag and 120_frag and 180_nl) has TP 2, FP 0, FN 1, TN 3
  # giving exact MCC of 6 divided by sqrt(72), which is 1 / sqrt(2) (0.7071)
  synth_df <- data.frame(
    msLevel = rep(2L, 6L),
    rtime = seq_len(6L) * 1.0,
    acquisitionNum = seq_len(6L),
    TITLE = paste0("Spec_", seq_len(6L)),
    PRECURSOR_MZ = rep(300.0, 6L),
    SKELETON = c("SkelA", "SkelA", "SkelA", "SkelB", "SkelB", "SkelB"),
    stringsAsFactors = FALSE
  )
  synth_df$mz <- list(
    c(100.0, 120.0),
    c(100.0, 120.0),
    c(100.0, 130.0),
    c(150.0, 160.0),
    c(150.0, 160.0),
    c(150.0, 170.0)
  )
  synth_df$intensity <- rep(list(c(100.0, 80.0)), 6L)

  synth_sp <- df_to_spectra(synth_df)
  res_sp <- spectra_to_queries(
    spectra = synth_sp,
    export = NULL,
    n_skel_min = 2L,
    n_spec_min = 2L
  )
  expected_mcc <- round(1 / sqrt(2), 4L)
  testthat::expect_identical(nrow(res_sp), 2L)
  testthat::expect_equal(
    res_sp$mcc,
    c(expected_mcc, expected_mcc),
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
      stringsAsFactors = FALSE
    )
  )

  res_no_peaks <- spectra_to_queries(
    sub_df,
    export = NULL,
    intensity_min = 101
  )
  testthat::expect_identical(
    res_no_peaks,
    data.frame(
      skeleton = character(0),
      mcc = numeric(0),
      combined_ions = character(0),
      stringsAsFactors = FALSE
    )
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
