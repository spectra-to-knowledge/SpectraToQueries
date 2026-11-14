#' @title Format numbers without scientific notation
#' @keywords internal
#' @noRd
fmt_num <- function(x) {
  vapply(
    x,
    function(v) format(v, scientific = FALSE, trim = TRUE, digits = 15L),
    character(1L)
  )
}

#' @title Canonicalize and absorb DNF ion clauses
#' @keywords internal
#' @noRd
absorb_dnf_clauses <- function(ion_lists) {
  sort_ions <- function(ions) {
    parsed <- parse_query_ions(unique(ions))
    frags_sorted <- if (length(parsed$frags) > 0L) {
      paste0(fmt_num(sort(parsed$frags)), "_frag")
    } else {
      character(0)
    }
    nls_sorted <- if (length(parsed$nls) > 0L) {
      paste0(fmt_num(sort(parsed$nls)), "_nl")
    } else {
      character(0)
    }
    c(frags_sorted, nls_sorted)
  }

  clauses <- unique(Filter(length, lapply(ion_lists, sort_ions)))
  n_clauses <- length(clauses)
  if (n_clauses <= 1L) {
    return(clauses)
  }

  # Boolean absorption: A | (A & B) = A (drop any strict superset of a clause)
  lens <- lengths(clauses)
  keep_min <- rep(TRUE, n_clauses)
  for (i in seq_len(n_clauses)) {
    ci <- clauses[[i]]
    li <- lens[i]
    for (j in seq_len(n_clauses)) {
      if (
        keep_min[j] &&
          li < lens[j] &&
          all(match(ci, clauses[[j]], nomatch = 0L) > 0L)
      ) {
        keep_min[j] <- FALSE
      }
    }
  }
  clauses[keep_min]
}

#' @title Combine ions into minimal AND/OR expression with sorting
#'
#' @description Simplifies a list of conjunctive ion clauses into a minimal
#'   Boolean expression by canonical sorting (`_frag` then `_nl`), Boolean
#'   absorption (`A | (A & B) = A`), and common-factor extraction.
#'
#' @param ion_lists List of character vectors representing ion combinations
#'
#' @return Single character string representing nested AND/OR expression
#'
#' @keywords internal
#' @noRd
combine_ions_minimal <- function(ion_lists) {
  ion_lists <- absorb_dnf_clauses(ion_lists)
  if (length(ion_lists) == 0L) {
    return("")
  }
  if (length(ion_lists) == 1L) {
    return(paste(ion_lists[[1L]], collapse = " & "))
  }

  common_sorted <- Reduce(
    function(x, y) x[match(x, y, nomatch = 0L) > 0L],
    ion_lists
  )

  remaining <- lapply(ion_lists, function(x) {
    x[match(x, common_sorted, nomatch = 0L) == 0L]
  })

  build_and <- function(sorted) {
    if (length(sorted) == 1L) {
      return(sorted)
    }
    paste0("(", paste(sorted, collapse = " & "), ")")
  }

  or_parts <- vapply(remaining, build_and, character(1L))
  or_expr <- paste(or_parts, collapse = " | ")

  if (length(common_sorted) > 0L) {
    common_expr <- paste(common_sorted, collapse = " & ")
    paste0(common_expr, " & (", or_expr, ")")
  } else {
    or_expr
  }
}

#' @title Parse Boolean ion query string into absorbed DNF clauses
#' @keywords internal
#' @noRd
parse_massql_dnf <- function(q) {
  tokens <- regmatches(
    q,
    gregexpr("\\(|\\)|&|\\||[^()&|[:space:]]+", q)
  )[[1L]]
  n_tok <- length(tokens)
  if (n_tok == 0L) {
    return(list())
  }

  state <- new.env(parent = emptyenv())
  state$pos <- 1L
  peek <- function() if (state$pos <= n_tok) tokens[state$pos] else ""

  dnf_and <- function(lhs, rhs) {
    if (length(lhs) == 0L || length(rhs) == 0L) {
      return(list())
    }
    out <- vector("list", length(lhs) * length(rhs))
    idx <- 0L
    for (a in lhs) {
      for (b in rhs) {
        idx <- idx + 1L
        out[[idx]] <- c(a, b)
      }
    }
    absorb_dnf_clauses(out)
  }

  parse_or <- function() {
    res <- parse_and()
    while (peek() == "|") {
      state$pos <- state$pos + 1L
      res <- absorb_dnf_clauses(c(res, parse_and()))
    }
    res
  }

  parse_and <- function() {
    res <- parse_atom()
    while (peek() == "&") {
      state$pos <- state$pos + 1L
      res <- dnf_and(res, parse_atom())
    }
    res
  }

  parse_atom <- function() {
    if (state$pos > n_tok) {
      abort_input("Unexpected end of MassQL query expression.")
    }
    tok <- peek()
    if (tok == "(") {
      state$pos <- state$pos + 1L
      inner <- parse_or()
      if (peek() != ")") {
        abort_input("Unmatched '(' in MassQL query expression.")
      }
      state$pos <- state$pos + 1L
      return(inner)
    }
    if (tok %in% c(")", "&", "|")) {
      abort_input(sprintf(
        "Unexpected token '%s' in MassQL query expression.",
        tok
      ))
    }
    state$pos <- state$pos + 1L
    absorb_dnf_clauses(list(tok))
  }

  clauses <- parse_or()
  if (state$pos <= n_tok) {
    abort_input(sprintf(
      "Unexpected token '%s' in MassQL query expression.",
      peek()
    ))
  }
  clauses
}

#' @title Convert ion queries to MassQL syntax
#'
#' @description Translates one or more ion query expressions (such as
#'   `"108.0804_frag & 241.095_nl"`, `"100_frag & (200_frag | 300_frag)"`, or a
#'   list of ion vectors/expressions) into MassQL queries
#'   (`QUERY scaninfo(MS2DATA) WHERE ...`).
#'
#' @details MassQL ignores `TOLERANCEMZ` when `TOLERANCEPPM` is present on the
#'   same condition. When `dalton > 0`, `to_massql()` emits a single effective
#'   `TOLERANCEMZ` modifier (`round(dalton + max(mz) * ppm / 1e6, 6)` for
#'   `MS2PROD`, and static `dalton` for `MS2NL`). When `dalton == 0` and
#'   `ppm > 0`, it emits `TOLERANCEPPM`. Disjunctive normal form clauses that
#'   cannot share a single `OR` group are joined with `" ; "`.
#'
#' @param query Character vector of query expressions, or a `list` of character
#'   vectors of ions or expression strings.
#' @param dalton Tolerance in Dalton (`TOLERANCEMZ`). Default is `0.01`.
#' @param ppm Tolerance in parts per million (`TOLERANCEPPM`). Default is
#'   `30.0`.
#' @param intensity_min Minimal relative intensity threshold in percent
#'   (`INTENSITYPERCENT`). Omitted when `<= 0`. Default is `0.0`.
#'
#' @return Character vector of MassQL query strings.
#'
#' @export
#'
#' @examples
#' to_massql("108.0804_frag & 241.095_nl", dalton = 0.01, ppm = 30)
to_massql <- function(
  query,
  dalton = 0.01,
  ppm = 30.0,
  intensity_min = 0.0
) {
  nonneg_msg <- "must be a single finite non-negative number."
  check_num(dalton, paste("`dalton`", nonneg_msg), lower = 0)
  check_num(ppm, paste("`ppm`", nonneg_msg), lower = 0)
  check_num(
    intensity_min,
    "`intensity_min` must be a single finite number in [0, 100].",
    lower = 0,
    upper = 100
  )

  if (is.null(query) || length(query) == 0L) {
    return(character(0))
  }
  if (is.list(query)) {
    query <- vapply(
      query,
      function(x) {
        if (is.null(x) || length(x) == 0L) {
          return("")
        }
        if (length(x) == 1L) {
          return(as.character(x))
        }
        combine_ions_minimal(list(as.character(x)))
      },
      character(1L)
    )
  } else {
    query <- as.character(query)
  }

  int_mod <- if (intensity_min > 0) {
    paste0(":INTENSITYPERCENT=", fmt_num(intensity_min))
  } else {
    ""
  }

  tol_mod <- function(mzs, is_nl = FALSE) {
    if (dalton == 0 && ppm > 0) {
      return(paste0(":TOLERANCEPPM=", fmt_num(ppm), int_mod))
    }
    if (dalton == 0) {
      return(paste0(":TOLERANCEMZ=0", int_mod))
    }
    eff_da <- if (is_nl) {
      dalton
    } else {
      round(dalton + max(mzs) * (ppm / 1e6), 6L)
    }
    paste0(":TOLERANCEMZ=", fmt_num(eff_da), int_mod)
  }

  format_conjunction <- function(ions_vec, extra_clause = NULL) {
    parsed <- parse_query_ions(ions_vec)
    frag_clauses <- if (length(parsed$frags) > 0L) {
      vapply(
        parsed$frags,
        function(mz) paste0("MS2PROD=", fmt_num(mz), tol_mod(mz, FALSE)),
        character(1L)
      )
    } else {
      character(0)
    }
    nl_clauses <- if (length(parsed$nls) > 0L) {
      paste0("MS2NL=", fmt_num(parsed$nls), tol_mod(parsed$nls, TRUE))
    } else {
      character(0)
    }
    clauses <- c(frag_clauses, nl_clauses, extra_clause)
    if (length(clauses) == 0L) {
      return("")
    }
    paste0(
      "QUERY scaninfo(MS2DATA) WHERE ",
      paste(clauses, collapse = " AND ")
    )
  }

  vapply(
    query,
    function(q) {
      if (is.na(q) || !nzchar(trimws(q))) {
        return("")
      }
      clauses <- parse_massql_dnf(q)
      if (length(clauses) == 0L) {
        return("")
      }
      if (length(clauses) == 1L) {
        return(format_conjunction(clauses[[1L]]))
      }

      common_ions <- Reduce(
        function(x, y) x[match(x, y, nomatch = 0L) > 0L],
        clauses
      )
      remaining <- lapply(clauses, function(x) {
        x[match(x, common_ions, nomatch = 0L) == 0L]
      })

      if (all(lengths(remaining) == 1L)) {
        alt_ions <- unlist(remaining, use.names = FALSE)
        alt_parsed <- parse_query_ions(alt_ions)
        if (length(alt_parsed$frags) == length(alt_ions)) {
          or_clause <- paste0(
            "MS2PROD=(",
            paste(fmt_num(alt_parsed$frags), collapse = " OR "),
            ")",
            tol_mod(alt_parsed$frags, FALSE)
          )
          return(format_conjunction(common_ions, or_clause))
        }
        if (length(alt_parsed$nls) == length(alt_ions)) {
          or_clause <- paste0(
            "MS2NL=(",
            paste(fmt_num(alt_parsed$nls), collapse = " OR "),
            ")",
            tol_mod(alt_parsed$nls, TRUE)
          )
          return(format_conjunction(common_ions, or_clause))
        }
      }

      sub_queries <- vapply(
        clauses,
        format_conjunction,
        character(1L)
      )
      sub_queries <- sub_queries[nzchar(sub_queries)]
      paste(sub_queries, collapse = " ; ")
    },
    character(1L),
    USE.NAMES = FALSE
  )
}
