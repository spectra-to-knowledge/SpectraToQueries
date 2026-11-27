#' Weight of the redundancy penalty in `select_mrmr_ions()`: score =
#' best conjunctive MCC - weight * mean |phi| with already selected ions, phi
#' computed on negative-class spectra only. Empirical choice, not derived.
#' @keywords internal
#' @noRd
mrmr_redundancy_weight <- 0.3

#' @title Select greedy hitting-set cover among 100%-recall diagnostic ions
#' @keywords internal
#' @noRd
select_hitting_set_ions <- function(cand_recall, neg_mat, ions_max) {
  selected <- 1L
  diag_pool <- which(cand_recall == 1)
  if (length(diag_pool) == 0L) {
    return(selected)
  }

  d_first <- diag_pool[1L]
  if (!(d_first %in% selected) && length(selected) < ions_max) {
    selected <- c(selected, d_first)
  }
  diag_fp_vec <- neg_mat[, d_first]
  rem_diag <- setdiff(diag_pool, selected)

  while (
    length(rem_diag) > 0L &&
      any(diag_fp_vec) &&
      length(selected) < ions_max
  ) {
    fp_after <- colSums(diag_fp_vec & neg_mat[, rem_diag, drop = FALSE])
    min_pos <- which.min(fp_after)
    if (fp_after[min_pos] < sum(diag_fp_vec)) {
      chosen_d <- rem_diag[min_pos]
      selected <- c(selected, chosen_d)
      diag_fp_vec <- diag_fp_vec & neg_mat[, chosen_d]
      rem_diag <- rem_diag[-min_pos]
    } else {
      break
    }
  }
  selected
}

#' @title Fill remaining ion slots via forward conjunctive MCC and mRMR
#' @keywords internal
#' @noRd
select_forward_mrmr_ions <- function(
  selected,
  cand_mcc,
  sub_mat,
  neg_mat,
  phi_mat,
  is_group,
  ions_max
) {
  n_pos <- sum(is_group)
  total_spectra <- nrow(sub_mat)
  neg_mask <- !is_group
  pos_mat <- sub_mat[is_group, , drop = FALSE]

  n_unique <- ncol(sub_mat)
  remaining <- setdiff(seq_len(n_unique), selected)
  active_conjs <- lapply(selected, function(idx) sub_mat[, idx])
  active_mccs <- cand_mcc[selected]

  # Track raw which.max(mccs_k) across active conjunctions incrementally
  cache <- new.env(parent = emptyenv())
  cache$mcc <- rep(-Inf, n_unique)
  cache$k <- integer(n_unique)

  update_conj_cache <- function(ac_vec, k_id) {
    if (length(remaining) == 0L) {
      return(invisible(NULL))
    }
    ac_pos <- which(ac_vec[is_group])
    ac_neg <- which(ac_vec[neg_mask])
    tp_v <- if (length(ac_pos) > 0L) {
      colSums(pos_mat[ac_pos, remaining, drop = FALSE])
    } else {
      numeric(length(remaining))
    }
    fp_v <- if (length(ac_neg) > 0L) {
      colSums(neg_mat[ac_neg, remaining, drop = FALSE])
    } else {
      numeric(length(remaining))
    }
    mcc_v <- calc_mcc_vec(tp_v, fp_v, n_pos, total_spectra)
    better <- mcc_v > cache$mcc[remaining]
    if (any(better)) {
      idx_b <- remaining[better]
      cache$mcc[idx_b] <- mcc_v[better]
      cache$k[idx_b] <- k_id
    }
  }

  for (k_id in seq_along(active_conjs)) {
    update_conj_cache(active_conjs[[k_id]], k_id)
  }

  while (length(selected) < ions_max && length(remaining) > 0L) {
    redundancy <- rowMeans(phi_mat[remaining, selected, drop = FALSE])
    rk <- cache$k[remaining]
    rmcc <- cache$mcc[remaining]
    improved <- rk > 0L & rmcc > active_mccs[pmax(1L, rk)]
    best_conj_mcc <- cand_mcc[remaining]
    best_conj_mcc[improved] <- rmcc[improved]
    best_conj_k <- integer(length(remaining))
    best_conj_k[improved] <- rk[improved]

    scores <- best_conj_mcc - mrmr_redundancy_weight * redundancy
    best_pos <- which.max(scores)
    chosen <- remaining[best_pos]
    selected <- c(selected, chosen)
    remaining <- remaining[-best_pos]

    chosen_vec <- sub_mat[, chosen]
    active_conjs[[length(active_conjs) + 1L]] <- chosen_vec
    active_mccs <- c(active_mccs, cand_mcc[chosen])
    update_conj_cache(chosen_vec, length(active_conjs))

    if (best_conj_mcc[best_pos] > cand_mcc[chosen]) {
      conj_vec <- active_conjs[[best_conj_k[best_pos]]] & chosen_vec
      active_conjs[[length(active_conjs) + 1L]] <- conj_vec
      active_mccs <- c(active_mccs, best_conj_mcc[best_pos])
      update_conj_cache(conj_vec, length(active_conjs))
    }
  }
  selected
}

#' @title Select candidate ions via forward conjunction gain and mRMR
#' @keywords internal
#' @noRd
select_mrmr_ions <- function(
  cand_ions,
  cand_mcc,
  cand_recall,
  nz_mat,
  is_group,
  ions_max
) {
  n_cand <- length(cand_ions)
  if (n_cand <= 1L) {
    return(cand_ions)
  }

  is_frag <- is_frag_ion(cand_ions)
  ord <- order(cand_mcc, cand_recall, is_frag, decreasing = TRUE)
  cand_ions <- cand_ions[ord]
  cand_mcc <- cand_mcc[ord]
  cand_recall <- cand_recall[ord]

  # Deduplicate ions with identical binary presence/absence vectors
  sub_mat <- nz_mat[, cand_ions, drop = FALSE]
  keep_unique <- !duplicated(sub_mat, MARGIN = 2L)
  cand_ions <- cand_ions[keep_unique]
  cand_mcc <- cand_mcc[keep_unique]
  cand_recall <- cand_recall[keep_unique]
  sub_mat <- sub_mat[, keep_unique, drop = FALSE]

  if (length(cand_ions) <= ions_max) {
    return(cand_ions)
  }

  # False-positive (negative-class) binary phi correlation for redundancy
  neg_mat <- sub_mat[!is_group, , drop = FALSE]
  phi_mat <- suppressWarnings(abs(stats::cor(neg_mat)))
  phi_mat[!is.finite(phi_mat)] <- 0.0

  selected <- select_hitting_set_ions(cand_recall, neg_mat, ions_max)
  selected <- select_forward_mrmr_ions(
    selected = selected,
    cand_mcc = cand_mcc,
    sub_mat = sub_mat,
    neg_mat = neg_mat,
    phi_mat = phi_mat,
    is_group = is_group,
    ions_max = ions_max
  )

  cand_ions[selected]
}

#' @title Extract candidate ions by group with mRMR redundancy filtering
#' @keywords internal
#' @noRd
extract_candidate_ions <- function(
  merged_mat,
  total_spectra,
  n_skel_min,
  n_spec_min,
  ions_max,
  precision_min,
  recall_min,
  mcc_min
) {
  groups_vec <- normalize_skeleton_labels(rownames(merged_mat))
  group_counts_tab <- table(groups_vec)
  group_count_vec <- as.integer(group_counts_tab[groups_vec])

  # Count ion presence across all spectra so small-group FPs are included
  nz_mat <- if (is.logical(merged_mat)) merged_mat else merged_mat != 0
  col_counts <- colSums(nz_mat)
  keep_ions <- col_counts >= n_spec_min
  nz_mat <- nz_mat[, keep_ions, drop = FALSE]
  col_counts <- col_counts[keep_ions]

  keep_rows <- group_count_vec >= n_skel_min
  nz_sub <- nz_mat[keep_rows, , drop = FALSE]
  groups_sub <- groups_vec[keep_rows]
  if (nrow(nz_sub) == 0L || ncol(nz_sub) == 0L) {
    return(list())
  }

  counts_mat <- rowsum(+nz_sub, groups_sub, reorder = FALSE)
  nz_idx <- which(counts_mat > 0L, arr.ind = TRUE)
  if (nrow(nz_idx) == 0L) {
    return(list())
  }

  row_idx <- nz_idx[, 1L]
  col_idx <- nz_idx[, 2L]
  tp <- as.integer(counts_mat[nz_idx])
  group_count <- as.integer(group_counts_tab[rownames(counts_mat)][row_idx])
  count_per_ion <- as.integer(col_counts[col_idx])

  ions_df <- data.frame(
    group = rownames(counts_mat)[row_idx],
    group_count = group_count,
    ion = as.character(colnames(nz_mat))[col_idx],
    count_per_ion = count_per_ion,
    count_per_ion_per_group = tp,
    precision = tp / count_per_ion,
    recall = tp / group_count,
    mcc = calc_mcc_vec(
      tp = tp,
      fp = count_per_ion - tp,
      n_pos = group_count,
      total_spectra = total_spectra
    ),
    stringsAsFactors = FALSE
  )

  ions_df <- ions_df[
    ions_df$precision >= precision_min &
      ions_df$recall >= recall_min &
      ions_df$mcc >= mcc_min,
    ,
    drop = FALSE
  ]
  if (nrow(ions_df) == 0L) {
    return(list())
  }

  groups <- sort(unique(ions_df$group), method = "radix")
  ions_list <- lapply(groups, function(g) {
    sub_df <- ions_df[ions_df$group == g, , drop = FALSE]
    select_mrmr_ions(
      cand_ions = sub_df$ion,
      cand_mcc = sub_df$mcc,
      cand_recall = sub_df$recall,
      nz_mat = nz_mat,
      is_group = groups_vec == g,
      ions_max = ions_max
    )
  })
  names(ions_list) <- groups
  ions_list
}
