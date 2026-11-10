#' @title Combine ions into minimal AND/OR expression with sorting
#'
#' @description Simplifies a list of conjunctive ion clauses into a minimal
#'   Boolean expression by canonical sorting (`_frag` then `_nl`), Boolean
#'   absorption (`A | (A & B) = A`), and common-factor extraction.
#'
#' @param ion_lists List of character vectors representing ion combinations
#'
#' @return Single character string representing nested AND/OR expression
combine_ions_minimal <- function(ion_lists) {
  # Helper to split and sort fragments/nls, return only non-empty unique ions
  sort_ions <- function(ions) {
    ions <- unique(ions[!is.na(ions) & nzchar(ions)])
    frags <- ions[grepl("_frag", ions, fixed = TRUE)]
    nls <- ions[grepl("_nl", ions, fixed = TRUE)]
    frags_sorted <- if (length(frags) > 0L) {
      paste0(
        sort(as.numeric(sub("_frag", "", frags, fixed = TRUE))),
        "_frag"
      )
    } else {
      character(0)
    }
    nls_sorted <- if (length(nls) > 0L) {
      paste0(
        sort(as.numeric(sub("_nl", "", nls, fixed = TRUE))),
        "_nl"
      )
    } else {
      character(0)
    }
    c(frags_sorted, nls_sorted)
  }

  # Canonicalize and remove empty entries
  ion_lists <- unique(Filter(length, lapply(ion_lists, sort_ions)))
  if (length(ion_lists) == 0L) {
    return("")
  }

  # Boolean absorption: A | (A & B) = A (drop any strict superset of a clause)
  n_clauses <- length(ion_lists)
  if (n_clauses > 1L) {
    lens <- lengths(ion_lists)
    keep_min <- rep(TRUE, n_clauses)
    for (i in seq_len(n_clauses)) {
      ci <- ion_lists[[i]]
      li <- lens[i]
      for (j in seq_len(n_clauses)) {
        if (
          keep_min[j] &&
            li < lens[j] &&
            all(match(ci, ion_lists[[j]], nomatch = 0L) > 0L)
        ) {
          keep_min[j] <- FALSE
        }
      }
    }
    ion_lists <- ion_lists[keep_min]
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

  # Build AND expression
  build_and <- function(sorted) {
    if (length(sorted) == 1L) {
      return(sorted)
    }
    paste0("(", paste(sorted, collapse = " & "), ")")
  }

  # Combine remaining sets with OR
  or_parts <- vapply(remaining, build_and, character(1L))

  # Combine with common ions
  if (length(common_sorted) > 0L) {
    common_expr <- paste(common_sorted, collapse = " & ")
    paste(common_expr, paste(or_parts, collapse = " | "), sep = " & ")
  } else {
    paste(or_parts, collapse = " | ")
  }
}

#' @title Convert ion queries to MassQL syntax
#'
#' @description Translates one or more ion query expressions (such as
#'   `"108.0804_frag & 241.095_nl"` or a list of ion vectors) into directly
#'   executable MassQL queries (`QUERY scaninfo(MS2DATA) WHERE ...`).
#'
#' @param query Character vector of query expressions, or a `list` of character
#'   vectors of ions.
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
  if (is.list(query)) {
    query <- vapply(
      query,
      function(x) combine_ions_minimal(list(as.character(x))),
      character(1L)
    )
  }
  query <- as.character(query)
  if (length(query) == 0L) {
    return(character(0))
  }

  mods <- paste0(
    ":TOLERANCEMZ=",
    dalton,
    ":TOLERANCEPPM=",
    ppm,
    if (is.finite(intensity_min) && intensity_min > 0) {
      paste0(":INTENSITYPERCENT=", intensity_min)
    } else {
      ""
    }
  )

  format_conjunction <- function(ions_vec, extra_clause = NULL) {
    parsed <- parse_query_ions(ions_vec)
    clauses <- c(
      if (length(parsed$frags) > 0L) paste0("MS2PROD=", parsed$frags, mods),
      if (length(parsed$nls) > 0L) paste0("MS2NL=", parsed$nls, mods),
      extra_clause
    )
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
      or_parts <- trimws(strsplit(q, "|", fixed = TRUE)[[1L]])
      if (length(or_parts) == 1L) {
        tokens <- trimws(gsub(
          "[()]",
          "",
          strsplit(or_parts, "&", fixed = TRUE)[[1L]]
        ))
        return(format_conjunction(tokens))
      }

      # Extract common prefix before first OR branch if unparenthesized
      first_part <- or_parts[1L]
      if (startsWith(first_part, "(") && endsWith(first_part, ")")) {
        common_ions <- character(0)
      } else {
        first_tokens <- trimws(strsplit(first_part, "&", fixed = TRUE)[[1L]])
        paren_open <- which(startsWith(first_tokens, "("))[1L]
        if (!is.na(paren_open)) {
          common_ions <- first_tokens[seq_len(paren_open - 1L)]
          or_parts[1L] <- paste(
            first_tokens[paren_open:length(first_tokens)],
            collapse = " & "
          )
        } else {
          n_ft <- length(first_tokens)
          common_ions <- if (n_ft > 1L) first_tokens[-n_ft] else character(0)
          or_parts[1L] <- first_tokens[n_ft]
        }
      }

      alt_lists <- lapply(or_parts, function(p) {
        trimws(gsub("[()]", "", strsplit(p, "&", fixed = TRUE)[[1L]]))
      })

      # Single-ion OR alternatives of same type map to MS2PROD/MS2NL=(a OR b)
      if (all(lengths(alt_lists) == 1L)) {
        alt_ions <- unlist(alt_lists, use.names = FALSE)
        alt_parsed <- parse_query_ions(alt_ions)
        if (
          length(alt_parsed$frags) == length(alt_ions) &&
            length(alt_ions) > 0L
        ) {
          or_clause <- paste0(
            "MS2PROD=(",
            paste(alt_parsed$frags, collapse = " OR "),
            ")",
            mods
          )
          return(format_conjunction(common_ions, or_clause))
        }
        if (
          length(alt_parsed$nls) == length(alt_ions) &&
            length(alt_ions) > 0L
        ) {
          or_clause <- paste0(
            "MS2NL=(",
            paste(alt_parsed$nls, collapse = " OR "),
            ")",
            mods
          )
          return(format_conjunction(common_ions, or_clause))
        }
      }

      sub_queries <- vapply(
        alt_lists,
        function(alt) format_conjunction(c(common_ions, alt)),
        character(1L)
      )
      sub_queries <- sub_queries[nzchar(sub_queries)]
      paste(sub_queries, collapse = " ; ")
    },
    character(1L),
    USE.NAMES = FALSE
  )
}
