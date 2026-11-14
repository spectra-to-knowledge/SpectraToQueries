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

testthat::test_that("package functions do not access Spectra @backend@ slots", {
  ns <- asNamespace("SpectraToQueries")
  obj_names <- ls(ns, all.names = TRUE)
  has_slot <- vapply(
    obj_names,
    function(nm) {
      fn <- get(nm, envir = ns, inherits = FALSE)
      if (!is.function(fn)) {
        return(FALSE)
      }
      txt <- paste(deparse(body(fn)), collapse = "\n")
      grepl("@backend@", txt, fixed = TRUE)
    },
    logical(1L)
  )
  testthat::expect_false(any(has_slot))

  r_dir <- testthat::test_path("..", "..", "R")
  if (dir.exists(r_dir)) {
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
    file_hits <- vapply(
      r_files,
      function(f) {
        any(grepl("@backend@", readLines(f, warn = FALSE), fixed = TRUE))
      },
      logical(1L)
    )
    testthat::expect_false(any(file_hits))
  }
})

testthat::test_that("only public Rd help pages are shipped", {
  expected_rd <- sort(
    c(
      "SpectraToQueries-package.Rd",
      "mia_spectra_df.Rd",
      "mia_spectra_grouped_df.Rd",
      "spectra_to_queries.Rd",
      "to_massql.Rd"
    ),
    method = "radix"
  )
  man_dir <- testthat::test_path("..", "..", "man")
  actual_rd <- if (dir.exists(man_dir)) {
    sort(list.files(man_dir, pattern = "\\.Rd$"), method = "radix")
  } else {
    sort(names(tools::Rd_db("SpectraToQueries")), method = "radix")
  }
  testthat::expect_identical(actual_rd, expected_rd)
})

testthat::test_that("bundled datasets have no directory paths in dataOrigin", {
  for (ds_name in c("mia_spectra_df", "mia_spectra_grouped_df")) {
    env <- new.env(parent = emptyenv())
    utils::data(list = ds_name, package = "SpectraToQueries", envir = env)
    df <- env[[ds_name]]
    if ("dataOrigin" %in% names(df)) {
      origins <- df$dataOrigin[!is.na(df$dataOrigin)]
      testthat::expect_false(any(grepl("[/\\\\]", origins)))
    }
  }
})

testthat::test_that("DESCRIPTION Collate and namespace imports stay in sync", {
  desc_path <- testthat::test_path("..", "..", "DESCRIPTION")
  r_dir <- testthat::test_path("..", "..", "R")
  tests_dir <- testthat::test_path("..", "..", "tests", "testthat")

  if (file.exists(desc_path) && dir.exists(r_dir)) {
    dcf <- read.dcf(desc_path)

    collate_files <- scan(
      text = dcf[1L, "Collate"],
      what = character(),
      quiet = TRUE
    )
    actual_r_files <- sort(
      list.files(r_dir, pattern = "\\.R$"),
      method = "radix"
    )
    testthat::expect_identical(
      sort(collate_files, method = "radix"),
      actual_r_files
    )

    parse_pkg_field <- function(val) {
      if (is.na(val) || !nzchar(trimws(val))) {
        return(character())
      }
      parts <- trimws(unlist(strsplit(val, ",", fixed = TRUE)))
      parts <- sub("\\s*\\(.*$", "", parts)
      parts[nzchar(parts)]
    }

    extract_pkg_calls <- function(paths) {
      found <- character()
      for (p in paths) {
        lines <- readLines(p, warn = FALSE)
        lines <- lines[!grepl("^\\s*#", lines)]
        m <- gregexpr("[A-Za-z][A-Za-z0-9.]*::", lines, perl = TRUE)
        hits <- unlist(regmatches(lines, m))
        found <- c(found, sub("::$", "", hits))
      }
      unique(found)
    }

    imports_pkgs <- parse_pkg_field(dcf[1L, "Imports"])
    suggests_pkgs <- parse_pkg_field(dcf[1L, "Suggests"])

    r_paths <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
    r_pkgs <- extract_pkg_calls(r_paths)
    testthat::expect_true(
      all(r_pkgs %in% c(imports_pkgs, "SpectraToQueries", "base"))
    )

    if (dir.exists(tests_dir)) {
      test_paths <- list.files(tests_dir, pattern = "\\.R$", full.names = TRUE)
      test_pkgs <- extract_pkg_calls(test_paths)
      testthat::expect_true(
        all(
          test_pkgs %in%
            c(imports_pkgs, suggests_pkgs, "SpectraToQueries", "base")
        )
      )
    }
  }
})
