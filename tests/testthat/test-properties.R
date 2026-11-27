testthat::test_that("golden MIADB raw and grouped tables match in tarball", {
  raw_res <- suppressMessages(spectra_to_queries(export = NULL))
  expected_raw <- data.frame(
    skeleton = c(
      "Lapidilectine",
      "Cinchona",
      "Ajmalicinespiro",
      "Corynantheanspiro",
      "Chanofruticosinate",
      "Kopsane",
      "Sarpagine",
      "Vobasine",
      "Macroline",
      "Akuammiline",
      "Eburna",
      "Strychnidine",
      "Akuammicine",
      "Iboga",
      "Aspidosperma",
      "Aspidofractane",
      "Ajmalicine",
      "Ervatamia",
      "Ajmaline",
      "Yohimbinoid",
      "Corynanthean"
    ),
    mcc = c(
      1,
      1,
      0.9559,
      0.9114,
      0.8139,
      0.7722,
      0.7415,
      0.7059,
      0.6733,
      0.6678,
      0.6463,
      0.6264,
      0.6026,
      0.5963,
      0.5681,
      0.5681,
      0.5681,
      0.5303,
      0.5156,
      0.5144,
      0.3872
    ),
    combined_ions = c(
      "223.0805_nl & 253.0909_nl",
      "134.0957_frag & 137.1209_nl & 166.1226_nl",
      "124.0426_frag & 227.115_nl",
      "253.1657_nl & 270.1587_nl",
      "210.0999_frag & 88.0173_nl",
      "135.0664_frag & 320.1035_frag",
      "134.0597_frag & 151.0527_frag & 106.0712_nl & 134.088_nl",
      "152.0603_frag & 188.1291_nl",
      "131.085_frag & 144.1059_frag & 169.0841_frag",
      "191.0968_frag & 195.0907_nl",
      "77.0842_nl & 107.0593_nl",
      "294.1048_frag & 144.1014_nl",
      "196.1118_frag & 17.0067_nl",
      "150.1262_frag & 101.1184_nl",
      "168.1938_frag",
      "277.1246_frag & 367.1656_frag",
      "173.0824_frag & 230.1189_frag & 253.133_nl",
      "140.108_frag & 226.1185_nl",
      "394.1874_nl",
      "212.1303_frag & 227.1477_nl",
      "143.065_nl & 169.1433_nl"
    ),
    massql = c(
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2NL=223.0805:TOLERANCEMZ=0.01 AND MS2NL=253.0909:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=134.0957:TOLERANCEMZ=0.014023 AND",
        "MS2NL=137.1209:TOLERANCEMZ=0.01 AND MS2NL=166.1226:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=124.0426:TOLERANCEMZ=0.013722 AND",
        "MS2NL=227.115:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2NL=253.1657:TOLERANCEMZ=0.01 AND MS2NL=270.1587:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=210.0999:TOLERANCEMZ=0.016303 AND",
        "MS2NL=88.0173:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=135.0664:TOLERANCEMZ=0.014052 AND",
        "MS2PROD=320.1035:TOLERANCEMZ=0.019604"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=134.0597:TOLERANCEMZ=0.014022 AND",
        "MS2PROD=151.0527:TOLERANCEMZ=0.014532 AND",
        "MS2NL=106.0712:TOLERANCEMZ=0.01 AND MS2NL=134.088:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=152.0603:TOLERANCEMZ=0.014562 AND",
        "MS2NL=188.1291:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=131.085:TOLERANCEMZ=0.013933 AND",
        "MS2PROD=144.1059:TOLERANCEMZ=0.014324 AND",
        "MS2PROD=169.0841:TOLERANCEMZ=0.015073"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=191.0968:TOLERANCEMZ=0.015733 AND",
        "MS2NL=195.0907:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2NL=77.0842:TOLERANCEMZ=0.01 AND MS2NL=107.0593:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=294.1048:TOLERANCEMZ=0.018824 AND",
        "MS2NL=144.1014:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=196.1118:TOLERANCEMZ=0.015884 AND",
        "MS2NL=17.0067:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=150.1262:TOLERANCEMZ=0.014504 AND",
        "MS2NL=101.1184:TOLERANCEMZ=0.01"
      ),
      "QUERY scaninfo(MS2DATA) WHERE MS2PROD=168.1938:TOLERANCEMZ=0.015046",
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=277.1246:TOLERANCEMZ=0.018314 AND",
        "MS2PROD=367.1656:TOLERANCEMZ=0.021015"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=173.0824:TOLERANCEMZ=0.015193 AND",
        "MS2PROD=230.1189:TOLERANCEMZ=0.016904 AND",
        "MS2NL=253.133:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=140.108:TOLERANCEMZ=0.014204 AND",
        "MS2NL=226.1185:TOLERANCEMZ=0.01"
      ),
      "QUERY scaninfo(MS2DATA) WHERE MS2NL=394.1874:TOLERANCEMZ=0.01",
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=212.1303:TOLERANCEMZ=0.016364 AND",
        "MS2NL=227.1477:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2NL=143.065:TOLERANCEMZ=0.01 AND MS2NL=169.1433:TOLERANCEMZ=0.01"
      )
    ),
    stringsAsFactors = FALSE
  )
  testthat::expect_equal(raw_res, expected_raw, tolerance = 1e-12)

  grp_res <- suppressMessages(spectra_to_queries(
    spectra = "grouped",
    export = NULL
  ))
  expected_grp <- data.frame(
    skeleton = c(
      "Cinchona",
      "Lapidilectine.Lunduran",
      "Ajmalicinespiro.Corynantheanspiro",
      "Chanofruticosinate",
      "Kopsane",
      "Akuammiline.Vincorinan",
      "Eburna",
      "Sarpagine.Vobasine",
      "Strychnidine",
      "Akagerane.Macroline",
      "Akuammicine.Condylocarpan",
      "Arboricininan.Antirhinan.Vallesiachotaman",
      "Quebrachamine.Andranginan.Iboga",
      "Ervatamia.Indolocinchona",
      "Ajmaline.Apogeissoschizane",
      "Leuconoxane.Secoaspidosperma.Aspidofractane.Aspidosperma",
      "Yohimbinoid.Arboricinan.Corynanthean.Ajmalicine.Arbophyllininan"
    ),
    mcc = c(
      1,
      1,
      0.8234,
      0.8139,
      0.7722,
      0.6992,
      0.6463,
      0.6367,
      0.6264,
      0.6207,
      0.5502,
      0.5303,
      0.526,
      0.4953,
      0.4929,
      0.4338,
      0.3866
    ),
    combined_ions = c(
      "134.0957_frag & 137.1209_nl & 166.1226_nl",
      "240.0973_frag & 199.0833_nl & 223.0805_nl",
      "132.0446_frag & 142.0647_frag & 253.133_nl",
      "210.0999_frag & 88.0173_nl",
      "135.0664_frag & 320.1035_frag",
      "191.0968_frag & 220.076_frag & 133.0745_nl",
      "77.0842_nl & 107.0593_nl",
      "164.068_frag & 217.1126_nl",
      "294.1048_frag & 144.1014_nl",
      "131.085_frag & 144.1059_frag & 169.0841_frag",
      "196.1118_frag & 17.0067_nl",
      "107.049_frag & 338.1574_frag",
      "116.0538_frag & 138.1283_nl",
      "140.108_frag & 226.1185_nl",
      "394.1874_nl",
      "217.0881_frag & 253.1163_frag & 309.1728_frag",
      "145.0882_frag & 143.065_nl"
    ),
    massql = c(
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=134.0957:TOLERANCEMZ=0.014023 AND",
        "MS2NL=137.1209:TOLERANCEMZ=0.01 AND MS2NL=166.1226:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=240.0973:TOLERANCEMZ=0.017203 AND",
        "MS2NL=199.0833:TOLERANCEMZ=0.01 AND MS2NL=223.0805:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=132.0446:TOLERANCEMZ=0.013962 AND",
        "MS2PROD=142.0647:TOLERANCEMZ=0.014262 AND",
        "MS2NL=253.133:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=210.0999:TOLERANCEMZ=0.016303 AND",
        "MS2NL=88.0173:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=135.0664:TOLERANCEMZ=0.014052 AND",
        "MS2PROD=320.1035:TOLERANCEMZ=0.019604"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=191.0968:TOLERANCEMZ=0.015733 AND",
        "MS2PROD=220.076:TOLERANCEMZ=0.016603 AND",
        "MS2NL=133.0745:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2NL=77.0842:TOLERANCEMZ=0.01 AND MS2NL=107.0593:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=164.068:TOLERANCEMZ=0.014923 AND",
        "MS2NL=217.1126:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=294.1048:TOLERANCEMZ=0.018824 AND",
        "MS2NL=144.1014:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=131.085:TOLERANCEMZ=0.013933 AND",
        "MS2PROD=144.1059:TOLERANCEMZ=0.014324 AND",
        "MS2PROD=169.0841:TOLERANCEMZ=0.015073"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=196.1118:TOLERANCEMZ=0.015884 AND",
        "MS2NL=17.0067:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=107.049:TOLERANCEMZ=0.013212 AND",
        "MS2PROD=338.1574:TOLERANCEMZ=0.020145"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=116.0538:TOLERANCEMZ=0.013482 AND",
        "MS2NL=138.1283:TOLERANCEMZ=0.01"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=140.108:TOLERANCEMZ=0.014204 AND",
        "MS2NL=226.1185:TOLERANCEMZ=0.01"
      ),
      "QUERY scaninfo(MS2DATA) WHERE MS2NL=394.1874:TOLERANCEMZ=0.01",
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=217.0881:TOLERANCEMZ=0.016513 AND",
        "MS2PROD=253.1163:TOLERANCEMZ=0.017594 AND",
        "MS2PROD=309.1728:TOLERANCEMZ=0.019276"
      ),
      paste(
        "QUERY scaninfo(MS2DATA) WHERE",
        "MS2PROD=145.0882:TOLERANCEMZ=0.014353 AND",
        "MS2NL=143.065:TOLERANCEMZ=0.01"
      )
    ),
    stringsAsFactors = FALSE
  )
  testthat::expect_equal(grp_res, expected_grp, tolerance = 1e-12)
})

testthat::test_that("golden subset pins ppm-only, da-only, intensity_min", {
  sub <- make_miadb_subset()

  res_ppm <- suppressMessages(spectra_to_queries(
    sub,
    export = NULL,
    ions_max = 4L,
    dalton = 0,
    ppm = 30
  ))
  testthat::expect_equal(
    res_ppm,
    data.frame(
      skeleton = c("Iboga", "Aspidosperma", "Corynanthean"),
      mcc = c(0.6896, 0.6576, 0.4588),
      combined_ions = c(
        "198.0893_frag & 224.1045_frag",
        "167.096_nl",
        "145.0832_frag & 143.0735_nl"
      ),
      massql = c(
        paste(
          "QUERY scaninfo(MS2DATA) WHERE",
          "MS2PROD=198.0893:TOLERANCEPPM=30 AND",
          "MS2PROD=224.1045:TOLERANCEPPM=30"
        ),
        "QUERY scaninfo(MS2DATA) WHERE MS2NL=167.096:TOLERANCEPPM=30",
        paste(
          "QUERY scaninfo(MS2DATA) WHERE",
          "MS2PROD=145.0832:TOLERANCEPPM=30 AND MS2NL=143.0735:TOLERANCEPPM=30"
        )
      ),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-12
  )

  res_da <- suppressMessages(spectra_to_queries(
    sub,
    export = NULL,
    ions_max = 4L,
    dalton = 0.01,
    ppm = 0
  ))
  testthat::expect_equal(
    res_da,
    data.frame(
      skeleton = c("Aspidosperma", "Iboga", "Corynanthean"),
      mcc = c(0.6947, 0.6896, 0.5021),
      combined_ions = c(
        "134.0951_frag & 199.1217_nl",
        "121.0922_frag & 224.1051_frag",
        "116.0562_frag & 143.0735_nl"
      ),
      massql = c(
        paste(
          "QUERY scaninfo(MS2DATA) WHERE",
          "MS2PROD=134.0951:TOLERANCEMZ=0.01 AND",
          "MS2NL=199.1217:TOLERANCEMZ=0.01"
        ),
        paste(
          "QUERY scaninfo(MS2DATA) WHERE",
          "MS2PROD=121.0922:TOLERANCEMZ=0.01 AND",
          "MS2PROD=224.1051:TOLERANCEMZ=0.01"
        ),
        paste(
          "QUERY scaninfo(MS2DATA) WHERE",
          "MS2PROD=116.0562:TOLERANCEMZ=0.01 AND",
          "MS2NL=143.0735:TOLERANCEMZ=0.01"
        )
      ),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-12
  )

  res_int <- suppressMessages(spectra_to_queries(
    sub,
    export = NULL,
    ions_max = 4L,
    intensity_min = 10
  ))
  testthat::expect_equal(
    res_int,
    data.frame(
      skeleton = c("Iboga", "Corynanthean", "Aspidosperma"),
      mcc = c(0.7347, 0.6796, 0.6425),
      combined_ions = c("122.0964_frag", "115.0536_frag", "156.0792_frag"),
      massql = c(
        paste0(
          "QUERY scaninfo(MS2DATA) WHERE ",
          "MS2PROD=122.0964:TOLERANCEMZ=0.013663:INTENSITYPERCENT=10"
        ),
        paste0(
          "QUERY scaninfo(MS2DATA) WHERE ",
          "MS2PROD=115.0536:TOLERANCEMZ=0.013452:INTENSITYPERCENT=10"
        ),
        paste0(
          "QUERY scaninfo(MS2DATA) WHERE ",
          "MS2PROD=156.0792:TOLERANCEMZ=0.014683:INTENSITYPERCENT=10"
        )
      ),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-12
  )
})

testthat::test_that("precompute_ion_hits pins closed tolerance boundaries", {
  # Binary-exact values: 1953.125 / 1e6 == 2^-9, so 128 * 2^-9 == 0.25 exact
  ppm_exact <- 1953.125
  prec <- rep(384.0, 4L)
  ok <- rep(TRUE, 4L)

  # dalton = 0.25, ppm = 1953.125 -> frag tol = 0.50; nl tol = 0.25
  all_mz_both <- list(128.5, 128.5 + 1e-9, 256.25, 256.25 + 1e-9)
  hits_both <- precompute_ion_hits(
    c("128_frag", "128_nl"),
    all_mz_both,
    prec,
    ok,
    dalton = 0.25,
    ppm = ppm_exact
  )
  testthat::expect_identical(
    hits_both[["128_frag"]],
    c(TRUE, FALSE, FALSE, FALSE)
  )
  testthat::expect_identical(
    hits_both[["128_nl"]],
    c(FALSE, FALSE, TRUE, FALSE)
  )

  # dalton = 0, ppm = 1953.125 -> frag tol = 0.25; nl tol = 0.25
  all_mz_quarter <- list(128.25, 128.25 + 1e-9, 256.25, 256.25 + 1e-9)
  hits_ppm <- precompute_ion_hits(
    c("128_frag", "128_nl"),
    all_mz_quarter,
    prec,
    ok,
    dalton = 0,
    ppm = ppm_exact
  )
  testthat::expect_identical(
    hits_ppm[["128_frag"]],
    c(TRUE, FALSE, FALSE, FALSE)
  )
  testthat::expect_identical(hits_ppm[["128_nl"]], c(FALSE, FALSE, TRUE, FALSE))

  # dalton = 0.25, ppm = 0 -> tol = 0.25 for both
  hits_da <- precompute_ion_hits(
    c("128_frag", "128_nl"),
    all_mz_quarter,
    prec,
    ok,
    dalton = 0.25,
    ppm = 0
  )
  testthat::expect_identical(
    hits_da[["128_frag"]],
    c(TRUE, FALSE, FALSE, FALSE)
  )
  testthat::expect_identical(hits_da[["128_nl"]], c(FALSE, FALSE, TRUE, FALSE))

  # dalton = 0, ppm = 0 -> exact match only
  all_mz_zero <- list(128.0, 128.0 + 1e-9, 256.0, 256.0 + 1e-9)
  hits_zero <- precompute_ion_hits(
    c("128_frag", "128_nl"),
    all_mz_zero,
    prec,
    ok,
    dalton = 0,
    ppm = 0
  )
  testthat::expect_identical(
    hits_zero[["128_frag"]],
    c(TRUE, FALSE, FALSE, FALSE)
  )
  testthat::expect_identical(
    hits_zero[["128_nl"]],
    c(FALSE, FALSE, TRUE, FALSE)
  )
})

testthat::test_that("round-trip, idempotence, permutation, monotonicity hold", {
  # Round-trip spectra_to_df strips dataOrigin directory prefixes
  df_rt <- make_spectra_df(
    c("A", "B"),
    list(c(100.5, 200.25), c(150.0, 250.5)),
    list(c(40, 100), c(80, 100)),
    dataOrigin = c("/tmp/dir/a.mgf", NA_character_)
  )
  rt <- spectra_to_df(Spectra::Spectra(df_rt))
  testthat::expect_identical(rt$dataOrigin, c("a.mgf", NA_character_))
  testthat::expect_equal(rt$mz, df_rt$mz)
  testthat::expect_equal(rt$intensity, df_rt$intensity)

  # Idempotence of canonical_ion_conjunction, absorb_dnf_clauses, to_massql
  ions <- c("200_frag", "50_nl", "100_frag", "25_nl")
  c1 <- canonical_ion_conjunction(ions)
  toks <- strsplit(c1, " & ", fixed = TRUE)[[1L]]
  testthat::expect_identical(canonical_ion_conjunction(toks), c1)

  clauses <- list(c("200_frag", "100_frag"), c("100_frag", "200_frag", "50_nl"))
  testthat::expect_identical(
    absorb_dnf_clauses(absorb_dnf_clauses(clauses)),
    absorb_dnf_clauses(clauses)
  )
  expr <- "100_frag & (200_frag | 300_frag)"
  testthat::expect_identical(
    to_massql(paste0("(", expr, ") | (", expr, ")")),
    to_massql(expr)
  )

  # Permutation invariance
  set.seed(11L)
  for (i in seq_len(10L)) {
    testthat::expect_identical(canonical_ion_conjunction(sample(ions)), c1)
  }

  # Permutation invariance of precompute_ion_hits & extract_candidate_ions
  sub_df <- make_miadb_subset()
  perm <- sample(nrow(sub_df))
  qd_orig <- prepare_query_data(Spectra::Spectra(sub_df), 0)
  qd_perm <- prepare_query_data(Spectra::Spectra(sub_df[perm, ]), 0)
  test_ions <- c("122.0964_frag", "144.0807_frag", "143.065_nl")
  h_orig <- precompute_ion_hits(
    test_ions,
    qd_orig$mz,
    qd_orig$precursor_mz,
    is.finite(qd_orig$precursor_mz) & lengths(qd_orig$mz) > 0L,
    0.01,
    30
  )
  h_perm <- precompute_ion_hits(
    test_ions,
    qd_perm$mz,
    qd_perm$precursor_mz,
    is.finite(qd_perm$precursor_mz) & lengths(qd_perm$mz) > 0L,
    0.01,
    30
  )
  for (u in test_ions) {
    testthat::expect_identical(h_perm[[u]], h_orig[[u]][perm])
  }

  hm_orig <- vapply(
    test_ions,
    function(u) as.numeric(h_orig[[u]]),
    numeric(nrow(sub_df))
  )
  rownames(hm_orig) <- qd_orig$skeleton
  hm_perm <- hm_orig[perm, , drop = FALSE]
  ec_orig <- extract_candidate_ions(hm_orig, nrow(sub_df), 5L, 3L, 10L, 0, 0, 0)
  ec_perm <- extract_candidate_ions(hm_perm, nrow(sub_df), 5L, 3L, 10L, 0, 0, 0)
  testthat::expect_identical(ec_perm, ec_orig)

  # Monotonicity: wider tolerance keeps >= hits; higher intensity_min <= peaks
  h_narrow <- precompute_ion_hits(
    test_ions,
    qd_orig$mz,
    qd_orig$precursor_mz,
    is.finite(qd_orig$precursor_mz) & lengths(qd_orig$mz) > 0L,
    0.005,
    10
  )
  h_wide <- precompute_ion_hits(
    test_ions,
    qd_orig$mz,
    qd_orig$precursor_mz,
    is.finite(qd_orig$precursor_mz) & lengths(qd_orig$mz) > 0L,
    0.02,
    50
  )
  for (u in test_ions) {
    testthat::expect_true(all(h_wide[[u]] >= h_narrow[[u]]))
  }
  sp_sub <- Spectra::Spectra(sub_df)
  qd_0 <- prepare_query_data(sp_sub, 0)
  qd_20 <- prepare_query_data(sp_sub, 20)
  testthat::expect_true(all(lengths(qd_20$mz) <= lengths(qd_0$mz)))

  # Monotonicity: stricter n_skel_min or mcc_min produces subset of skeletons
  r_loose <- suppressMessages(spectra_to_queries(
    sub_df,
    export = NULL,
    ions_max = 3L,
    n_skel_min = 5L,
    mcc_min = 0.0
  ))
  r_strict <- suppressMessages(spectra_to_queries(
    sub_df,
    export = NULL,
    ions_max = 3L,
    n_skel_min = 20L,
    mcc_min = 0.2
  ))
  testthat::expect_true(all(r_strict$skeleton %in% r_loose$skeleton))
})

testthat::test_that("error paths cover check_num and MassQL parser", {
  testthat::expect_error(
    check_num(-5, "x", lower = 1),
    regexp = "`x` must be a single finite number >= 1\\.",
    class = "stq_input_error"
  )
  for (bad_q in c(
    "100_frag &",
    "& 100_frag",
    ")100_frag",
    "100_frag 200_frag",
    "100_frag)"
  )) {
    testthat::expect_error(to_massql(bad_q), class = "stq_input_error")
  }
})
