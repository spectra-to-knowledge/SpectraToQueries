#' @title Combine ions into minimal AND/OR expression with sorting
#'
#' @param ion_lists List of character vectors representing ion combinations
#'
#' @return Single character string representing nested AND/OR expression
combine_ions_minimal <- function(ion_lists) {
  # Remove empty/NA entries
  ion_lists <- lapply(ion_lists, function(x) x[!is.na(x) & nzchar(x)])
  ion_lists <- Filter(length, ion_lists)
  if (length(ion_lists) == 0L) {
    return("")
  }

  # Helper to split and sort fragments/nls, return only non-empty
  sort_ions <- function(ions) {
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

  common <- Reduce(
    function(x, y) x[match(x, y, nomatch = 0L) > 0L],
    ion_lists
  )
  common_sorted <- sort_ions(common)

  remaining <- lapply(ion_lists, function(x) {
    x[match(x, common, nomatch = 0L) == 0L]
  })

  # If all or any query is reduced to the common set, the OR branch simplifies
  # to TRUE, so the minimal expression is just the common ions.
  if (any(lengths(remaining) == 0L)) {
    return(paste(common_sorted, collapse = " & "))
  }

  # Generate unique remaining sets
  unique_sets <- unique(remaining)

  # Build AND expression
  build_and <- function(vec) {
    sorted <- sort_ions(vec)
    if (length(sorted) == 0L) {
      return(NA_character_)
    }
    if (length(sorted) == 1L) {
      return(sorted)
    }
    paste0("(", paste(sorted, collapse = " & "), ")")
  }

  # Combine unique sets with OR
  or_parts <- vapply(unique_sets, build_and, character(1L))
  or_parts <- or_parts[!is.na(or_parts) & nzchar(or_parts)]

  # Combine with common ions
  if (length(common_sorted) > 0L) {
    common_expr <- paste(common_sorted, collapse = " & ")
    if (length(or_parts) == 0L) {
      return(common_expr)
    }
    paste(common_expr, paste(or_parts, collapse = " | "), sep = " & ")
  } else {
    paste(or_parts, collapse = " | ")
  }
}
