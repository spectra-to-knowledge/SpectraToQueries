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
    c(
      sprintf("%s_frag", fmt_num(sort(parsed$frags))),
      sprintf("%s_nl", fmt_num(sort(parsed$nls)))
    )
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

#' @title Format one ion conjunction in canonical order
#' @description Deduplicates ion tokens and joins them with `" & "`, fragments
#'   first then neutral losses, each in ascending numeric order.
#' @param ions Character vector of `_frag` / `_nl` tokens.
#' @return Single string; `""` when `ions` has no non-empty tokens.
#' @keywords internal
#' @noRd
canonical_ion_conjunction <- function(ions) {
  clause <- absorb_dnf_clauses(list(ions))
  if (length(clause) == 0L) "" else paste(clause[[1L]], collapse = " & ")
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

  state <- new.env(parent = emptyenv())
  state$pos <- 1L
  peek <- function() if (state$pos <= n_tok) tokens[state$pos] else ""
  unexpected_tok <- function(tok) {
    abort_input(sprintf(
      "Unexpected token '%s' in MassQL query expression.",
      tok
    ))
  }

  dnf_and <- function(lhs, rhs) {
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
      unexpected_tok(tok)
    }
    state$pos <- state$pos + 1L
    absorb_dnf_clauses(list(tok))
  }

  clauses <- parse_or()
  if (state$pos <= n_tok) {
    unexpected_tok(peek())
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
#' @details Ion tokens are `<m/z>_frag` (fragment m/z) and `<mass>_nl`
#'   (neutral loss in Da), combined with `&`, `|` and parentheses (`&` binds
#'   tighter). Malformed tokens raise an error of class `stq_input_error`.
#'
#'   MassQL uses `TOLERANCEPPM` and ignores `TOLERANCEMZ` when both are given,
#'   and defaults to 0.1 Da when neither is given, so one tolerance modifier is
#'   emitted per condition:
#'   - `dalton > 0`: `TOLERANCEMZ` = `dalton + max(mz) * ppm / 1e6` rounded up
#'     to 1e-6 Da for `MS2PROD`, and `dalton` for `MS2NL`;
#'   - `dalton == 0`, `ppm > 0`: `TOLERANCEPPM` = `ppm`;
#'   - `dalton == 0`, `ppm == 0`: `TOLERANCEMZ=0`.
#'
#'   Differences from the matching used by [spectra_to_queries()] to compute
#'   MCC:
#'   - `MS2NL`: the pipeline matches a peak within
#'     `dalton + (precursor - nl) * ppm / 1e6` of `precursor - nl`; the emitted
#'     `TOLERANCEMZ` is `dalton` alone, and MassQL applies `TOLERANCEPPM` to the
#'     neutral-loss value. A query with neutral losses can match fewer spectra
#'     than its reported MCC implies.
#'   - `INTENSITYPERCENT`: MassQL keeps peaks strictly above the threshold; the
#'     pipeline keeps peaks at or above `intensity_min`.
#'
#'   Expressions are expanded to disjunctive normal form and absorbed
#'   (`A | (A & B) = A`). When all clauses share a common conjunction and differ
#'   by one ion of the same type, one `MS2PROD=(a OR b)` or `MS2NL=(a OR b)`
#'   condition is emitted. Otherwise one MassQL query per clause is returned,
#'   joined by `" ; "`; split on `" ; "` and run each query separately, as the
#'   MassQL grammar has no query separator.
#'
#' @references Damiani et al. (2025) A universal language for finding mass
#'   spectrometry data patterns. \doi{10.1038/s41592-025-02660-z};
#'   <https://mwang87.github.io/MassQueryLanguage_Documentation/>
#'
#' @param query Character vector of query expressions, or a `list` of character
#'   vectors of ions or expression strings.
#' @param dalton Absolute m/z tolerance in Da. Default is `0.01`.
#' @param ppm Relative m/z tolerance in parts per million, added to `dalton`.
#'   Default is `30.0`.
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
  check_num(dalton, "dalton", lower = 0)
  check_num(ppm, "ppm", lower = 0)
  check_num(intensity_min, "intensity_min", lower = 0, upper = 100)

  if (length(query) == 0L) {
    return(character(0))
  }
  if (is.list(query)) {
    query <- vapply(
      query,
      function(x) {
        if (length(x) == 1L) {
          as.character(x)
        } else {
          canonical_ion_conjunction(as.character(x))
        }
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
      # Round up to 1e-6 Da so the emitted window contains the match window;
      # the inner round() absorbs floating-point noise before ceiling().
      ceiling(round((dalton + max(mzs) * (ppm / 1e6)) * 1e6, 6L)) / 1e6
    }
    paste0(":TOLERANCEMZ=", fmt_num(eff_da), int_mod)
  }

  format_conjunction <- function(ions_vec, extra_clause = NULL) {
    parsed <- parse_query_ions(ions_vec)
    frag_clauses <- vapply(
      parsed$frags,
      function(mz) paste0("MS2PROD=", fmt_num(mz), tol_mod(mz, FALSE)),
      character(1L)
    )
    nl_clauses <- sprintf(
      "MS2NL=%s%s",
      fmt_num(parsed$nls),
      tol_mod(parsed$nls, TRUE)
    )
    clauses <- c(frag_clauses, nl_clauses, extra_clause)
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

      common_ions <- Reduce(intersect, clauses)
      remaining <- lapply(clauses, setdiff, common_ions)

      if (all(lengths(remaining) == 1L)) {
        alt_ions <- unlist(remaining, use.names = FALSE)
        alt_parsed <- parse_query_ions(alt_ions)
        is_all_frag <- length(alt_parsed$frags) == length(alt_ions)
        if (is_all_frag || length(alt_parsed$nls) == length(alt_ions)) {
          vals <- if (is_all_frag) alt_parsed$frags else alt_parsed$nls
          or_clause <- sprintf(
            "%s=(%s)%s",
            if (is_all_frag) "MS2PROD" else "MS2NL",
            paste(fmt_num(vals), collapse = " OR "),
            tol_mod(vals, !is_all_frag)
          )
          return(format_conjunction(common_ions, or_clause))
        }
      }

      paste(
        vapply(clauses, format_conjunction, character(1L)),
        collapse = " ; "
      )
    },
    character(1L),
    USE.NAMES = FALSE
  )
}
