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
  calc_mcc <- function(vec) {
    calc_mcc_vec(
      tp = sum(vec[is_group]),
      fp = sum(vec[neg_mask]),
      n_pos = n_pos,
      total_spectra = total_spectra
    )
  }

  remaining <- setdiff(seq_len(ncol(sub_mat)), selected)
  active_conjs <- lapply(selected, function(idx) sub_mat[, idx])
  active_mccs <- cand_mcc[selected]

  while (length(selected) < ions_max && length(remaining) > 0L) {
    redundancy <- rowMeans(phi_mat[remaining, selected, drop = FALSE])
    best_conj_mcc <- cand_mcc[remaining]
    best_conj_k <- integer(length(remaining))
    for (r_i in seq_along(remaining)) {
      col_v <- sub_mat[, remaining[r_i]]
      mccs_k <- vapply(
        active_conjs,
        function(ac) calc_mcc(ac & col_v),
        numeric(1L)
      )
      bk <- which.max(mccs_k)
      if (mccs_k[bk] > active_mccs[bk]) {
        best_conj_mcc[r_i] <- mccs_k[bk]
        best_conj_k[r_i] <- bk
      }
    }
    scores <- best_conj_mcc - mrmr_redundancy_weight * redundancy
    best_pos <- which.max(scores)
    chosen <- remaining[best_pos]
    selected <- c(selected, chosen)
    chosen_vec <- sub_mat[, chosen]
    active_conjs[[length(active_conjs) + 1L]] <- chosen_vec
    active_mccs <- c(active_mccs, cand_mcc[chosen])
    if (best_conj_mcc[best_pos] > cand_mcc[chosen]) {
      conj_vec <- active_conjs[[best_conj_k[best_pos]]] & chosen_vec
      active_conjs[[length(active_conjs) + 1L]] <- conj_vec
      active_mccs <- c(active_mccs, best_conj_mcc[best_pos])
    }
    remaining <- remaining[-best_pos]
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
  nz_mat <- merged_mat != 0
  col_counts <- colSums(nz_mat)
  keep_ions <- col_counts >= n_spec_min
  nz_mat <- nz_mat[, keep_ions, drop = FALSE]
  col_counts <- col_counts[keep_ions]

  keep_rows <- group_count_vec >= n_skel_min
  nz_sub <- nz_mat[keep_rows, , drop = FALSE]
  groups_sub <- groups_vec[keep_rows]
  group_count_sub <- group_count_vec[keep_rows]

  # Pivot non-zero entries of qualifying skeleton groups to long format
  nz_idx <- which(nz_sub, arr.ind = TRUE)
  row_idx <- nz_idx[, 1L]
  col_idx <- nz_idx[, 2L]

  ions_df <- data.frame(
    group = groups_sub[row_idx],
    group_count = group_count_sub[row_idx],
    ion = as.character(colnames(nz_mat))[col_idx],
    count_per_ion = as.integer(col_counts[col_idx]),
    stringsAsFactors = FALSE
  )

  grp_ion_key <- paste(ions_df$group, ions_df$ion, sep = "\r")
  u_grp_ion <- unique(grp_ion_key)
  grp_ion_tab <- tabulate(
    match(grp_ion_key, u_grp_ion),
    nbins = length(u_grp_ion)
  )
  ions_df <- ions_df[!duplicated(grp_ion_key), , drop = FALSE]
  ions_df$count_per_ion_per_group <- grp_ion_tab

  ions_df$precision <- ions_df$count_per_ion_per_group / ions_df$count_per_ion
  ions_df$recall <- ions_df$count_per_ion_per_group / ions_df$group_count

  ions_df <- ions_df[
    ions_df$precision >= precision_min & ions_df$recall >= recall_min,
    ,
    drop = FALSE
  ]

  # Matthews Correlation Coefficient (MCC) for a single ion and group pair
  tp <- ions_df$count_per_ion_per_group
  fp <- ions_df$count_per_ion - tp
  ions_df$mcc <- calc_mcc_vec(
    tp = tp,
    fp = fp,
    n_pos = ions_df$group_count,
    total_spectra = total_spectra
  )

  ions_df <- ions_df[ions_df$mcc >= mcc_min, , drop = FALSE]
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
