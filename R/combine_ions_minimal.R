#' @title Combine ions into minimal AND/OR expression with sorting
#'
#' @param ion_lists List of character vectors representing ion combinations
#'
#' @return Single character string representing nested AND/OR expression
combine_ions_minimal <- function(ion_lists) {
  # Remove empty/NA
  ion_lists <- lapply(ion_lists, function(x) x[x != "" & !is.na(x)])
  if (length(ion_lists) == 0) {
    return("")
  }

  # Helper to split and sort fragments/nls, return only non-empty
  sort_ions <- function(ions) {
    frags <- grep("_frag$", ions, value = TRUE)
    nls <- grep("_nl$", ions, value = TRUE)
    frags_sorted <- if (length(frags)) {
      paste0(sort(as.numeric(sub("_frag$", "", frags))), "_frag")
    } else {
      character(0)
    }
    nls_sorted <- if (length(nls)) {
      paste0(sort(as.numeric(sub("_nl$", "", nls))), "_nl")
    } else {
      character(0)
    }
    c(frags_sorted, nls_sorted)
  }

  # Find common ions
  common <- Reduce(intersect, ion_lists)
  common_sorted <- sort_ions(common)

  # Remove common from each list
  remaining <- lapply(ion_lists, function(x) setdiff(x, common))

  if (all(lengths(remaining) == 0)) {
    return(paste(common_sorted, collapse = " & "))
  }

  # Generate unique remaining sets
  unique_sets <- unique(remaining)

  # Build AND expression
  build_and <- function(vec) {
    sorted <- sort_ions(vec)
    if (length(sorted) == 1) {
      return(sorted)
    }
    paste0("(", paste(sorted, collapse = " & "), ")")
  }

  # Combine unique sets with OR
  or_parts <- sapply(unique_sets, build_and)

  # Combine with common ions
  if (length(common_sorted) > 0) {
    common_expr <- paste(common_sorted, collapse = " & ")
    paste(common_expr, paste(or_parts, collapse = " | "), sep = " & ")
  } else {
    paste(or_parts, collapse = " | ")
  }
}
