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

  common_sorted <- Reduce(intersect, ion_lists)
  remaining <- lapply(ion_lists, setdiff, common_sorted)

  or_parts <- vapply(
    remaining,
    function(s) {
      if (length(s) == 1L) s else sprintf("(%s)", paste(s, collapse = " & "))
    },
    character(1L)
  )
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
          combine_ions_minimal(list(as.character(x)))
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
