#' @title Combine ions into minimal AND/OR expression with sorting
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
  ion_lists <- lapply(ion_lists, sort_ions)
  ion_lists <- Filter(length, ion_lists)
  if (length(ion_lists) == 0L) {
    return("")
  }

  common_sorted <- Reduce(
    function(x, y) x[match(x, y, nomatch = 0L) > 0L],
    ion_lists
  )

  remaining <- lapply(ion_lists, function(x) {
    x[match(x, common_sorted, nomatch = 0L) == 0L]
  })

  # If all or any query is reduced to the common set, the OR branch simplifies
  # to TRUE, so the minimal expression is just the common ions.
  if (any(lengths(remaining) == 0L)) {
    return(paste(common_sorted, collapse = " & "))
  }

  # Generate unique remaining sets (already canonically sorted)
  unique_sets <- unique(remaining)

  # Build AND expression
  build_and <- function(sorted) {
    if (length(sorted) == 1L) {
      return(sorted)
    }
    paste0("(", paste(sorted, collapse = " & "), ")")
  }

  # Combine unique sets with OR
  or_parts <- vapply(unique_sets, build_and, character(1L))

  # Combine with common ions
  if (length(common_sorted) > 0L) {
    common_expr <- paste(common_sorted, collapse = " & ")
    paste(common_expr, paste(or_parts, collapse = " | "), sep = " & ")
  } else {
    paste(or_parts, collapse = " | ")
  }
}
