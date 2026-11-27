#' @title Search ion combinations for one skeleton via Branch-and-Bound
#' @keywords internal
#' @noRd
prune_skeleton_bnb <- function(
  s,
  cands,
  skeleton,
  hits_lookup,
  max_ions,
  decimals,
  total_spectra
) {
  m <- length(cands)
  if (m == 0L) {
    return(NULL)
  }

  is_pos <- !is.na(skeleton) & (skeleton == s)
  is_neg <- !is.na(skeleton) & (skeleton != s)
  n_pos <- sum(is_pos)

  pos_hits <- lapply(cands, function(u) hits_lookup[[u]][is_pos])
  neg_hits <- lapply(cands, function(u) hits_lookup[[u]][is_neg])

  kept_combs <- as.list(cands)
  kept_tp <- as.numeric(vapply(pos_hits, sum, integer(1L)))
  kept_fp <- as.numeric(vapply(neg_hits, sum, integer(1L)))
  active <- lapply(seq_len(m), function(j) {
    list(
      idx = j,
      last = j,
      ph = pos_hits[[j]],
      nh = neg_hits[[j]],
      tp = kept_tp[j],
      fp = kept_fp[j]
    )
  })
  best_mcc <- max(round(
    calc_mcc_vec(kept_tp, kept_fp, n_pos, total_spectra),
    decimals
  ))
  max_k <- min(m, as.integer(max_ions))

  for (k in seq_len(max(0L, max_k - 1L))) {
    # Prune nodes that cannot strictly beat best_mcc at larger cardinality
    expandable <- Filter(
      function(node) {
        node$last < m &&
          node$tp > 0 &&
          node$fp > 0 &&
          round(
            calc_mcc_vec(node$tp, 0.0, n_pos, total_spectra),
            decimals
          ) >
            best_mcc
      },
      active
    )
    if (length(expandable) == 0L) {
      break
    }

    next_active <- list()
    for (node in expandable) {
      for (j in (node$last + 1L):m) {
        nh_new <- node$nh & neg_hits[[j]]
        fp_new <- sum(nh_new)
        # Prune if adding ion j fails to eliminate any false positives
        if (fp_new < node$fp) {
          ph_new <- node$ph & pos_hits[[j]]
          tp_new <- sum(ph_new)
          if (tp_new > 0) {
            mcc_new <- round(
              calc_mcc_vec(tp_new, fp_new, n_pos, total_spectra),
              decimals
            )
            if (mcc_new > best_mcc) {
              best_mcc <- mcc_new
            }
            idx_new <- c(node$idx, j)
            kept_idx <- length(kept_combs) + 1L
            kept_combs[[kept_idx]] <- cands[idx_new]
            kept_tp[kept_idx] <- tp_new
            kept_fp[kept_idx] <- fp_new
            next_active[[length(next_active) + 1L]] <- list(
              idx = idx_new,
              last = j,
              ph = ph_new,
              nh = nh_new,
              tp = tp_new,
              fp = fp_new
            )
          }
        }
      }
    }
    if (length(next_active) == 0L) {
      break
    }
    active <- next_active
  }

  list(
    skeleton = rep(s, length(kept_combs)),
    ions = kept_combs,
    tp = kept_tp,
    fp = kept_fp
  )
}

#' @title Prune redundant superset combinations via exact Branch-and-Bound
#' @keywords internal
#' @noRd
prune_combinations_bnb <- function(
  ions_list,
  skeleton,
  hits_lookup,
  max_ions,
  decimals
) {
  total_spectra <- sum(!is.na(skeleton))
  valid_names <- names(hits_lookup)
  skel_names <- names(ions_list)
  empty_eval <- list(
    skeleton = character(0),
    ions = list(),
    tp = numeric(0),
    fp = numeric(0)
  )
  if (length(skel_names) == 0L || total_spectra == 0) {
    return(empty_eval)
  }

  parts <- Filter(
    Negate(is.null),
    lapply(skel_names, function(s) {
      prune_skeleton_bnb(
        s = s,
        cands = intersect(ions_list[[s]], valid_names),
        skeleton = skeleton,
        hits_lookup = hits_lookup,
        max_ions = max_ions,
        decimals = decimals,
        total_spectra = total_spectra
      )
    })
  )

  list(
    skeleton = unlist(
      c(list(character(0)), lapply(parts, `[[`, "skeleton")),
      use.names = FALSE
    ),
    ions = do.call(c, c(list(list()), lapply(parts, `[[`, "ions"))),
    tp = unlist(
      c(list(numeric(0)), lapply(parts, `[[`, "tp")),
      use.names = FALSE
    ),
    fp = unlist(
      c(list(numeric(0)), lapply(parts, `[[`, "fp")),
      use.names = FALSE
    )
  )
}

#' @title Evaluate query MCC and summarize minimal best queries per skeleton
#' @keywords internal
#' @noRd
summarize_best_queries <- function(
  candidates,
  skeleton_counts,
  total_spectra,
  decimals,
  dalton,
  ppm,
  intensity_min
) {
  if (length(candidates$ions) == 0L) {
    return(empty_query_table())
  }

  skel_vec <- candidates$skeleton
  ions_vec <- candidates$ions
  tp <- candidates$tp
  fp <- candidates$fp
  n_pos <- skeleton_counts[skel_vec]

  mcc_raw_vec <- calc_mcc_vec(tp, fp, n_pos, total_spectra)
  mcc_vec <- round(mcc_raw_vec, decimals)
  prec_vec <- ifelse(tp + fp > 0, tp / (tp + fp), 0.0)

  n_ions_vec <- lengths(ions_vec)
  n_frags_vec <- vapply(
    ions_vec,
    function(x) sum(is_frag_ion(x)),
    integer(1L)
  )

  # Rank by rounded MCC, parsimony, unrounded MCC, precision, and fragments
  ord <- order(
    mcc_vec,
    -n_ions_vec,
    mcc_raw_vec,
    prec_vec,
    n_frags_vec,
    decreasing = TRUE
  )
  skel_vec <- skel_vec[ord]
  mcc_vec <- mcc_vec[ord]
  ions_vec <- ions_vec[ord]

  first <- !duplicated(skel_vec)
  combined_vec <- vapply(
    ions_vec[first],
    canonical_ion_conjunction,
    character(1L),
    USE.NAMES = FALSE
  )
  data.frame(
    skeleton = skel_vec[first],
    mcc = mcc_vec[first],
    combined_ions = combined_vec,
    massql = to_massql(
      combined_vec,
      dalton = dalton,
      ppm = ppm,
      intensity_min = intensity_min
    ),
    stringsAsFactors = FALSE
  )
}
