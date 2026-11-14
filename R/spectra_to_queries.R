#' @title Validate inputs for `spectra_to_queries`
#' @keywords internal
#' @noRd
validate_stq_inputs <- function(
  spectra,
  export,
  dalton,
  decimals,
  intensity_min,
  ions_max,
  n_skel_min,
  n_spec_min,
  ppm,
  mcc_min,
  precision_min,
  recall_min,
  zero_val
) {
  is_valid_str <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }

  valid_spectra <- is.null(spectra) ||
    is.data.frame(spectra) ||
    inherits(spectra, "Spectra") ||
    is_valid_str(spectra)
  if (!valid_spectra) {
    abort_input(
      paste(
        "`spectra` must be NULL, a non-empty string,",
        "a data.frame, or a Spectra object."
      )
    )
  }
  if (!is.null(export) && !is_valid_str(export)) {
    abort_input("`export` must be NULL or a non-empty character string.")
  }

  nonneg_msg <- "must be a single finite non-negative number."
  check_num(dalton, paste("`dalton`", nonneg_msg), lower = 0)
  check_num(ppm, paste("`ppm`", nonneg_msg), lower = 0)
  check_num(decimals, paste("`decimals`", nonneg_msg), lower = 0)

  ge1_msg <- "must be a single finite number >= 1."
  check_num(ions_max, paste("`ions_max`", ge1_msg), lower = 1)
  check_num(n_skel_min, paste("`n_skel_min`", ge1_msg), lower = 1)
  check_num(n_spec_min, paste("`n_spec_min`", ge1_msg), lower = 1)

  unit_msg <- "must be a single finite number in [0, 1]."
  check_num(precision_min, paste("`precision_min`", unit_msg), 0, 1)
  check_num(recall_min, paste("`recall_min`", unit_msg), 0, 1)
  check_num(
    mcc_min,
    "`mcc_min` must be a single finite number in [-1, 1].",
    lower = -1,
    upper = 1
  )

  check_num(
    intensity_min,
    "`intensity_min` must be a single finite number in [0, 100].",
    lower = 0,
    upper = 100
  )
  check_num(zero_val, "`zero_val` must be a single finite number.")
  invisible(TRUE)
}

#' @title Load or coerce spectra input
#' @keywords internal
#' @noRd
load_input_spectra <- function(spectra) {
  if (is.null(spectra)) {
    message("No spectra given, loading example spectra.")
    data_env <- new.env(parent = emptyenv())
    utils::data(
      "mia_spectra_df",
      package = "SpectraToQueries",
      envir = data_env
    )
    mia_spectra <- df_to_spectra(data_env$mia_spectra_df)
  } else if (
    is.character(spectra) && length(spectra) == 1L && spectra == "grouped"
  ) {
    message("Loading grouped example spectra.")
    data_env <- new.env(parent = emptyenv())
    utils::data(
      "mia_spectra_grouped_df",
      package = "SpectraToQueries",
      envir = data_env
    )
    mia_spectra <- df_to_spectra(data_env$mia_spectra_grouped_df)
  } else if (is.character(spectra)) {
    if (!file.exists(spectra)) {
      stop(errorCondition(
        sprintf("Spectra file not found: '%s'.", spectra),
        class = "stq_file_not_found_error"
      ))
    }
    message("Loading spectra from file.")
    mia_spectra <- spectra |>
      MsBackendMgf::readMgf() |>
      Spectra::Spectra()
  } else if (is.data.frame(spectra)) {
    mia_spectra <- df_to_spectra(spectra)
  } else {
    mia_spectra <- spectra
  }

  vars <- Spectra::spectraVariables(mia_spectra)
  if (!("SKELETON" %in% vars)) {
    if ("skeleton" %in% vars) {
      mia_spectra$SKELETON <- mia_spectra$skeleton
    } else {
      stop(errorCondition(
        "`spectra` must contain a `SKELETON` metadata variable.",
        class = "stq_input_error"
      ))
    }
  }
  if ("PRECURSOR_MZ" %in% vars) {
    mia_spectra$precursorMz <- as.numeric(mia_spectra$PRECURSOR_MZ)
  }
  mia_spectra
}

#' @title Preprocess spectra (reduce, combine by title, deisotope, normalize)
#' @keywords internal
#' @noRd
preprocess_spectra <- function(spectra_obj, dalton, ppm, intensity_min) {
  combine_by_title <- function(x) {
    titles <- x$TITLE
    if (!is.null(titles) && (anyDuplicated(titles) || is.unsorted(titles))) {
      Spectra::combineSpectra(x, f = titles)
    } else {
      x
    }
  }

  spectra_obj |>
    Spectra::reduceSpectra(tolerance = dalton, ppm = ppm) |>
    combine_by_title() |>
    Spectra::deisotopeSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::filterPrecursorPeaks(
      tolerance = dalton,
      ppm = ppm,
      mz = ">="
    ) |>
    Spectra::addProcessing(normalize_peaks()) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing() |>
    Spectra::filterEmptySpectra()
}

#' @title Preprocess neutral-loss spectra (reduce, normalize, filter intensity)
#' @keywords internal
#' @noRd
preprocess_nl_spectra <- function(spectra_obj, dalton, ppm, intensity_min) {
  spectra_obj |>
    Spectra::reduceSpectra(tolerance = dalton, ppm = ppm) |>
    Spectra::addProcessing(normalize_peaks()) |>
    Spectra::filterIntensity(intensity = c(intensity_min, Inf)) |>
    Spectra::applyProcessing()
}

#' @title Build feature matrix from harmonized spectra and fix m/z centroids
#' @keywords internal
#' @noRd
bin_and_fix_spectra <- function(
  spectra_obj,
  label,
  dalton,
  ppm,
  decimals,
  n_spec_min,
  reset_first = FALSE
) {
  if (sum(lengths(spectra_obj)) == 0L) {
    message("Fixing mzs (", label, ")")
    return(matrix(
      0,
      nrow = length(spectra_obj),
      ncol = 0L,
      dimnames = list(spectra_obj$SKELETON, character(0))
    ))
  }

  to_mat <- if (reset_first) Spectra::reset(spectra_obj) else spectra_obj
  mat <- create_matrix(to_mat, name = to_mat$SKELETON)
  message("Fixing mzs (", label, ")")
  fix_binned_mzs(
    mat,
    original_mzs = spectra_obj,
    decimals = decimals,
    dalton = dalton,
    ppm = ppm
  ) |>
    filter_matrix(n = n_spec_min)
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

  is_frag <- grepl("_frag", cand_ions, fixed = TRUE)
  ord <- order(cand_mcc, cand_recall, is_frag, decreasing = TRUE)
  cand_ions <- cand_ions[ord]
  cand_mcc <- cand_mcc[ord]
  cand_recall <- cand_recall[ord]

  # Deduplicate ions with identical binary presence/absence vectors
  sub_mat <- nz_mat[, cand_ions, drop = FALSE]
  col_keys <- apply(sub_mat, 2L, paste, collapse = "")
  keep_unique <- !duplicated(col_keys)
  cand_ions <- cand_ions[keep_unique]
  cand_mcc <- cand_mcc[keep_unique]
  cand_recall <- cand_recall[keep_unique]
  sub_mat <- sub_mat[, keep_unique, drop = FALSE]

  n_unique <- length(cand_ions)
  if (n_unique <= ions_max) {
    return(cand_ions)
  }

  # False-positive (negative-class) binary phi correlation for redundancy
  neg_mask <- !is_group
  neg_mat <- sub_mat[neg_mask, , drop = FALSE]
  storage.mode(neg_mat) <- "double"
  phi_mat <- suppressWarnings(abs(stats::cor(neg_mat)))
  phi_mat[!is.finite(phi_mat)] <- 0.0

  n_pos <- as.numeric(sum(is_group))
  n_neg <- as.numeric(sum(neg_mask))
  calc_mcc <- function(vec) {
    tp <- as.numeric(sum(vec[is_group]))
    fp <- as.numeric(sum(vec[neg_mask]))
    fn_ <- n_pos - tp
    tn <- n_neg - fp
    denom <- sqrt((tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_))
    if (denom > 0) (tp * tn - fp * fn_) / denom else 0.0
  }

  selected <- 1L
  diag_pool <- which(cand_recall == 1)
  if (length(diag_pool) > 0L) {
    d_first <- diag_pool[1L]
    if (!(d_first %in% selected) && length(selected) < ions_max) {
      selected <- c(selected, d_first)
    }
    diag_fp_vec <- sub_mat[neg_mask, d_first]
    rem_diag <- setdiff(diag_pool, selected)
    while (
      length(rem_diag) > 0L &&
        any(diag_fp_vec) &&
        length(selected) < ions_max
    ) {
      fp_after <- vapply(
        rem_diag,
        function(idx) sum(diag_fp_vec & sub_mat[neg_mask, idx]),
        integer(1L)
      )
      min_pos <- which.min(fp_after)
      if (fp_after[min_pos] < sum(diag_fp_vec)) {
        chosen_d <- rem_diag[min_pos]
        selected <- c(selected, chosen_d)
        diag_fp_vec <- diag_fp_vec & sub_mat[neg_mask, chosen_d]
        rem_diag <- rem_diag[-min_pos]
      } else {
        break
      }
    }
  }

  remaining <- setdiff(seq_len(n_unique), selected)
  active_conjs <- lapply(selected, function(idx) sub_mat[, idx])
  active_mccs <- cand_mcc[selected]

  while (length(selected) < ions_max && length(remaining) > 0L) {
    redundancy <- rowMeans(phi_mat[remaining, selected, drop = FALSE])
    best_conj_mcc <- numeric(length(remaining))
    best_conj_vec <- vector("list", length(remaining))
    for (r_i in seq_along(remaining)) {
      idx <- remaining[r_i]
      col_v <- sub_mat[, idx]
      mccs_k <- vapply(
        active_conjs,
        function(ac) calc_mcc(ac & col_v),
        numeric(1L)
      )
      bk <- which.max(mccs_k)
      if (mccs_k[bk] > active_mccs[bk]) {
        best_conj_mcc[r_i] <- mccs_k[bk]
        best_conj_vec[[r_i]] <- active_conjs[[bk]] & col_v
      } else {
        best_conj_mcc[r_i] <- cand_mcc[idx]
        best_conj_vec[[r_i]] <- col_v
      }
    }
    scores <- best_conj_mcc - 0.3 * redundancy
    best_pos <- which.max(scores)
    chosen <- remaining[best_pos]
    selected <- c(selected, chosen)
    active_conjs[[length(active_conjs) + 1L]] <- sub_mat[, chosen]
    active_mccs <- c(active_mccs, cand_mcc[chosen])
    if (best_conj_mcc[best_pos] > cand_mcc[chosen]) {
      active_conjs[[length(active_conjs) + 1L]] <- best_conj_vec[[best_pos]]
      active_mccs <- c(active_mccs, best_conj_mcc[best_pos])
    }
    remaining <- remaining[-best_pos]
  }

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
  groups_vec <- gsub(
    pattern = "+",
    replacement = ".",
    x = as.character(rownames(merged_mat)),
    fixed = TRUE
  )
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

  message("Extract the best ions to perform a query.")
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
  tp <- as.numeric(ions_df$count_per_ion_per_group)
  fp <- as.numeric(ions_df$count_per_ion) - tp
  fn_ <- as.numeric(ions_df$group_count) - tp
  tn <- total_spectra - tp - fp - fn_
  denom_sq <- (tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_)
  ions_df$mcc <- ifelse(
    denom_sq > 0,
    (tp * tn - fp * fn_) / sqrt(denom_sq),
    0.0
  )

  ions_df <- ions_df[
    ions_df$mcc >= mcc_min,
    ,
    drop = FALSE
  ]
  if (nrow(ions_df) == 0L) {
    return(list(
      ions_list = list(),
      ions_list_diagnostic = list()
    ))
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

  diag_df <- ions_df[ions_df$recall == 1, c("group", "ion"), drop = FALSE]
  d_groups <- sort(unique(diag_df$group), method = "radix")
  ions_list_diagnostic <- stats::setNames(
    lapply(d_groups, function(g) unique(diag_df$ion[diag_df$group == g])),
    d_groups
  )

  list(
    ions_list = ions_list,
    ions_list_diagnostic = ions_list_diagnostic
  )
}

#' @title Merge diagnostic ions into candidate combinations
#' @keywords internal
#' @noRd
merge_diagnostic_ions <- function(combinations, ions_list_diagnostic) {
  has_diagnostic <- match(
    names(combinations),
    names(ions_list_diagnostic),
    nomatch = 0L
  ) >
    0L

  if (any(has_diagnostic)) {
    for (name in names(combinations)[has_diagnostic]) {
      diagnostic_ions <- ions_list_diagnostic[[name]]
      if (length(diagnostic_ions) > 0L) {
        combinations[[name]] <- lapply(combinations[[name]], function(sublist) {
          unique(c(unlist(sublist), diagnostic_ions))
        })
      }
    }
  }
  combinations
}

#' @title Prune redundant superset combinations via exact Branch-and-Bound
#' @keywords internal
#' @noRd
prune_combinations_bnb <- function(
  ions_list,
  spectra,
  max_ions,
  dalton,
  ppm,
  decimals,
  intensity_min = 0.0
) {
  query_data <- prepare_query_data(spectra, intensity_min = intensity_min)
  all_mz <- query_data$mz
  precursor_mz <- query_data$precursor_mz
  skeleton <- query_data$skeleton
  total_spectra <- as.numeric(length(all_mz))

  unique_ions <- unique(unlist(ions_list, use.names = FALSE))
  unique_ions <- unique_ions[!is.na(unique_ions) & nzchar(unique_ions)]
  valid_prec_and_spec <- is.finite(precursor_mz) & (lengths(all_mz) > 0L)

  precomp <- precompute_ion_hits(
    unique_ions = unique_ions,
    all_mz = all_mz,
    precursor_mz = precursor_mz,
    valid_prec_and_spec = valid_prec_and_spec,
    dalton = dalton,
    ppm = ppm
  )
  hits_lookup <- precomp$hits

  skel_names <- names(ions_list)
  out <- stats::setNames(vector("list", length(skel_names)), skel_names)

  for (s in skel_names) {
    cands <- intersect(ions_list[[s]], names(hits_lookup))
    m <- length(cands)
    if (m == 0L) {
      out[[s]] <- list()
      next
    }

    is_pos <- !is.na(skeleton) & (skeleton == s)
    is_neg <- !is.na(skeleton) & (skeleton != s)
    n_pos <- as.numeric(sum(is_pos))
    n_neg <- total_spectra - n_pos

    pos_hits <- lapply(cands, function(u) hits_lookup[[u]][is_pos])
    neg_hits <- lapply(cands, function(u) hits_lookup[[u]][is_neg])

    eval_mcc <- function(tp, fp) {
      fn_ <- n_pos - tp
      tn <- n_neg - fp
      denom <- sqrt((tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_))
      if (denom > 0) (tp * tn - fp * fn_) / denom else 0.0
    }

    # Level 1: all single-ion queries
    kept_combs <- as.list(cands)
    active <- lapply(seq_len(m), function(j) {
      ph <- pos_hits[[j]]
      nh <- neg_hits[[j]]
      tp <- as.numeric(sum(ph))
      fp <- as.numeric(sum(nh))
      list(
        idx = j,
        last = j,
        ph = ph,
        nh = nh,
        tp = tp,
        fp = fp,
        mcc = round(eval_mcc(tp, fp), decimals)
      )
    })

    best_mcc <- max(vapply(active, `[[`, numeric(1L), "mcc"))
    max_k <- min(m, as.integer(max_ions))

    if (max_k >= 2L) {
      for (k in 2L:max_k) {
        # Prune nodes that cannot strictly beat best_mcc at larger cardinality
        expandable <- Filter(
          function(node) {
            node$last < m &&
              node$tp > 0 &&
              node$fp > 0 &&
              round(eval_mcc(node$tp, 0.0), decimals) > best_mcc
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
            fp_new <- as.numeric(sum(nh_new))
            # Prune if adding ion j fails to eliminate any false positives
            if (fp_new < node$fp) {
              ph_new <- node$ph & pos_hits[[j]]
              tp_new <- as.numeric(sum(ph_new))
              if (tp_new > 0) {
                mcc_new <- round(eval_mcc(tp_new, fp_new), decimals)
                idx_new <- c(node$idx, j)
                kept_combs[[length(kept_combs) + 1L]] <- cands[idx_new]
                next_active[[length(next_active) + 1L]] <- list(
                  idx = idx_new,
                  last = j,
                  ph = ph_new,
                  nh = nh_new,
                  tp = tp_new,
                  fp = fp_new,
                  mcc = mcc_new
                )
              }
            }
          }
        }
        if (length(next_active) == 0L) {
          break
        }
        lvl_max <- max(vapply(next_active, `[[`, numeric(1L), "mcc"))
        if (lvl_max > best_mcc) {
          best_mcc <- lvl_max
        }
        active <- next_active
      }
    }

    out[[s]] <- kept_combs
  }

  out
}

#' @title Evaluate query MCC and summarize minimal best queries per skeleton
#' @keywords internal
#' @noRd
summarize_best_queries <- function(
  queries_results,
  all_combinations,
  skeleton_counts,
  total_spectra,
  decimals,
  dalton = 0.01,
  ppm = 30.0,
  intensity_min = 0.0
) {
  stats_mat <- vapply(
    seq_along(queries_results),
    function(result) {
      target_name <- names(queries_results)[result]
      vals <- queries_results[[result]]$value
      tp <- as.numeric(sum(vals == target_name))
      fp <- as.numeric(length(vals)) - tp
      group_total <- as.numeric(skeleton_counts[[target_name]])
      fn_ <- group_total - tp
      tn <- total_spectra - tp - fp - fn_

      denom <- sqrt((tp + fp) * (tp + fn_) * (tn + fp) * (tn + fn_))
      mcc_raw <- if (denom > 0) (tp * tn - fp * fn_) / denom else 0.0
      prec <- if (tp + fp > 0) tp / (tp + fp) else 0.0
      c(round(mcc_raw, decimals), mcc_raw, prec)
    },
    numeric(3L)
  )

  results_stats <- stats_mat[1L, ]
  valid_mcc <- !is.nan(results_stats)
  skel_vec <- names(all_combinations)[valid_mcc]
  mcc_vec <- results_stats[valid_mcc]
  mcc_raw_vec <- stats_mat[2L, valid_mcc]
  prec_vec <- stats_mat[3L, valid_mcc]
  ions_vec <- all_combinations[valid_mcc]
  n_ions_vec <- lengths(ions_vec)
  n_frags_vec <- vapply(
    ions_vec,
    function(x) sum(grepl("_frag", x, fixed = TRUE)),
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

  max_mcc_by_skel <- tapply(mcc_vec, skel_vec, max)
  keep_best <- mcc_vec == as.numeric(max_mcc_by_skel[skel_vec])

  skel_best <- skel_vec[keep_best]
  mcc_best <- mcc_vec[keep_best]
  ions_best <- ions_vec[keep_best]

  u_skel <- unique(skel_best)
  combined_vec <- vapply(
    u_skel,
    function(s) combine_ions_minimal(ions_best[skel_best == s][1L]),
    character(1L),
    USE.NAMES = FALSE
  )
  final_queries <- data.frame(
    skeleton = u_skel,
    mcc = vapply(
      u_skel,
      function(s) mcc_best[skel_best == s][1L],
      numeric(1L),
      USE.NAMES = FALSE
    ),
    combined_ions = combined_vec,
    massql = to_massql(
      combined_vec,
      dalton = dalton,
      ppm = ppm,
      intensity_min = intensity_min
    ),
    stringsAsFactors = FALSE
  )
  final_queries <- final_queries[
    order(final_queries$mcc, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  rownames(final_queries) <- NULL
  final_queries
}

#' @title Spectra to queries
#'
#' @description Converts tandem mass spectra into diagnostic fragment (`_frag`)
#'   and neutral-loss (`_nl`) ion queries ranked by Matthews Correlation
#'   Coefficient (MCC), together with ready-to-run MassQL queries.
#'
#' @details
#' For each skeleton class with at least `n_skel_min` spectra, the pipeline:
#' 1. Normalizes MS2 intensities, removes sub-threshold and non-positive
#'    neutral-loss peaks, and aligns m/z values across spectra using
#'    intensity-weighted consensus centroids within `dalton + mz * ppm / 1e6`.
#' 2. Evaluates candidate fragment and neutral-loss centroids against the raw
#'    unbinned spectra (requiring `precursor_mz - nl > 0`).
#' 3. Selects up to `ions_max` candidate ions per skeleton using greedy
#'    hitting-set pruning on 100%-recall ions followed by forward conjunctive
#'    MCC selection with negative-class correlation penalty (mRMR).
#' 4. Prunes the power set of ion combinations via Branch-and-Bound when the
#'    optimistic zero-FP upper bound cannot reach the current best MCC.
#' 5. Returns the single most parsimonious query achieving the maximum MCC per
#'    skeleton and translates it to MassQL via [to_massql()].
#'
#' @param spectra A [Spectra::Spectra] object, a portable `data.frame` (such as
#'   [mia_spectra_df] or [mia_spectra_grouped_df]), or a file path to an MGF
#'   file. If `NULL`, loads the example [mia_spectra_df] dataset. If the string
#'   `"grouped"`, loads the grouped example [mia_spectra_grouped_df] dataset.
#' @param export Export file path for the resulting TSV table (or `NULL` to
#'   skip writing to disk). Default is `"data/interim/queries.tsv"`.
#' @param dalton Tolerance in Dalton. Default is `0.01`.
#' @param decimals Number of decimals for rounding. Default is `4L`.
#' @param intensity_min Minimal relative intensity threshold (0 to 100). Default
#'   is `0.0`.
#' @param ions_max Maximal number of ions in a query combination. Default is
#'   `10L`.
#' @param n_skel_min Minimal number of spectra per skeleton class. Default is
#'   `5L`.
#' @param n_spec_min Minimal number of spectra where a signal has to be found.
#'   Default is `3L`.
#' @param ppm Tolerance in parts per million. Default is `30.0`.
#' @param mcc_min Minimal single-ion Matthews Correlation Coefficient. Default
#'   is `0.0`.
#' @param precision_min Minimal single-ion precision (pre-filter before MCC).
#'   Default is `0.0`.
#' @param recall_min Minimal single-ion recall (pre-filter before MCC). Default
#'   is `0.0`.
#' @param zero_val Zero value threshold for intensity. Default is `0.0`.
#'
#' @return Invisibly returns a `data.frame` with columns `skeleton`, `mcc`,
#'   `combined_ions`, and `massql`, and writes the TSV table to `export` when
#'   `export` is non-`NULL`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Run with default example spectra (mia_spectra_df)
#'   spectra_to_queries(export = tempfile(fileext = ".tsv"))
#'
#'   # Run with grouped example spectra (mia_spectra_grouped_df)
#'   spectra_to_queries(
#'     spectra = "grouped",
#'     export = tempfile(fileext = ".tsv")
#'   )
#' }
spectra_to_queries <- function(
  spectra = NULL,
  export = "data/interim/queries.tsv",
  dalton = 0.01,
  decimals = 4L,
  intensity_min = 0.0,
  ions_max = 10L,
  n_skel_min = 5L,
  n_spec_min = 3L,
  ppm = 30.0,
  mcc_min = 0.0,
  precision_min = 0.0,
  recall_min = 0.0,
  zero_val = 0.0
) {
  validate_stq_inputs(
    spectra = spectra,
    export = export,
    dalton = dalton,
    decimals = decimals,
    intensity_min = intensity_min,
    ions_max = ions_max,
    n_skel_min = n_skel_min,
    n_spec_min = n_spec_min,
    ppm = ppm,
    mcc_min = mcc_min,
    precision_min = precision_min,
    recall_min = recall_min,
    zero_val = zero_val
  )

  mia_spectra <- load_input_spectra(spectra)
  mia_spectra <- Spectra::filterMsLevel(mia_spectra, 2L)
  skel_raw <- as.character(mia_spectra$SKELETON)
  valid_skel <- !is.na(skel_raw) & nzchar(trimws(skel_raw))
  if (!all(valid_skel)) {
    mia_spectra <- mia_spectra[valid_skel]
    skel_raw <- skel_raw[valid_skel]
  }

  total_spectra <- as.numeric(length(mia_spectra))
  if (total_spectra == 0) {
    empty_queries <- data.frame(
      skeleton = character(0),
      mcc = numeric(0),
      combined_ions = character(0),
      massql = character(0),
      stringsAsFactors = FALSE
    )
    if (!is.null(export)) {
      create_dir(dirname(export))
      utils::write.table(
        empty_queries,
        file = export,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
    return(invisible(empty_queries))
  }

  skel_norm <- gsub("+", ".", skel_raw, fixed = TRUE)
  if (length(unique(skel_norm)) < length(unique(skel_raw))) {
    warning(
      "Distinct `SKELETON` labels collide after '+' to '.' normalization.",
      call. = FALSE
    )
  }

  message(
    "Cut the fragments lower than ",
    intensity_min,
    " and rough preprocessing"
  )
  mia_spectra_1 <- preprocess_spectra(
    mia_spectra,
    dalton = dalton,
    ppm = ppm,
    intensity_min = intensity_min
  )

  message(
    "Harmonize m/z values across spectra, given ",
    dalton,
    " Da or ",
    ppm,
    " ppm tolerance."
  )
  mia_spectra_w <- harmonize_mzs(mia_spectra_1, dalton = dalton, ppm = ppm)

  message(
    "Calculate neutral losses, given ",
    dalton,
    " Da or ",
    ppm,
    " ppm tolerance."
  )
  message("Remove the ones above the precursor.")
  mia_spectra_nl <- mia_spectra_1 |>
    Spectra::neutralLoss(
      Spectra::PrecursorMzParam(
        filterPeaks = c("abovePrecursor"),
        msLevel = 2L,
        tolerance = dalton,
        ppm = ppm
      )
    ) |>
    preprocess_nl_spectra(
      dalton = dalton,
      ppm = ppm,
      intensity_min = intensity_min
    )

  mia_spectra_w_nl <- harmonize_mzs(mia_spectra_nl, dalton = dalton, ppm = ppm)

  message("Bin spectra to get a matrix.")
  message("The window is ", dalton, " Dalton")
  message("Create fragments and neutral losses matrices.")
  message("Remove features not appearing in at least ", n_spec_min, " spectra.")

  spectra_mat <- bin_and_fix_spectra(
    mia_spectra_w,
    label = "fragments",
    dalton = dalton,
    ppm = ppm,
    decimals = decimals,
    n_spec_min = n_spec_min,
    reset_first = FALSE
  )
  spectra_nl_mat <- bin_and_fix_spectra(
    mia_spectra_w_nl,
    label = "losses",
    dalton = dalton,
    ppm = ppm,
    decimals = decimals,
    n_spec_min = n_spec_min,
    reset_first = FALSE
  )

  message("Create a matrix containing fragments and neutral losses.")
  message("Round the values to ", decimals, ".")
  frag_names <- if (ncol(spectra_mat) > 0L) {
    paste0(round(as.numeric(colnames(spectra_mat)), decimals), "_frag")
  } else {
    character(0)
  }
  nl_names <- if (ncol(spectra_nl_mat) > 0L) {
    paste0(round(as.numeric(colnames(spectra_nl_mat)), decimals), "_nl")
  } else {
    character(0)
  }

  merged_mat <- cbind(spectra_mat, spectra_nl_mat)
  rm(spectra_mat, spectra_nl_mat)

  zeros <- colSums(merged_mat) <= zero_val
  if (!all(zeros)) {
    cand_cols <- unique(c(frag_names, nl_names)[!zeros])
    qd <- prepare_query_data(mia_spectra, intensity_min = intensity_min)
    pre <- precompute_ion_hits(
      unique_ions = cand_cols,
      all_mz = qd$mz,
      precursor_mz = qd$precursor_mz,
      valid_prec_and_spec = is.finite(qd$precursor_mz) & (lengths(qd$mz) > 0L),
      dalton = dalton,
      ppm = ppm
    )
    merged_mat <- matrix(
      vapply(
        cand_cols,
        function(u) as.numeric(pre$hits[[u]]),
        numeric(length(mia_spectra))
      ),
      nrow = length(mia_spectra),
      dimnames = list(mia_spectra$SKELETON, cand_cols)
    )
  } else {
    merged_mat <- matrix(
      0,
      nrow = nrow(merged_mat),
      ncol = 0L,
      dimnames = list(rownames(merged_mat), character(0))
    )
  }

  message("Count the number of members per skeleton and pivot the matrix.")
  candidates <- extract_candidate_ions(
    merged_mat = merged_mat,
    total_spectra = total_spectra,
    n_skel_min = n_skel_min,
    n_spec_min = n_spec_min,
    ions_max = ions_max,
    precision_min = precision_min,
    recall_min = recall_min,
    mcc_min = mcc_min
  )
  rm(merged_mat)

  ions_list <- candidates$ions_list
  if (length(ions_list) == 0L) {
    final_queries <- data.frame(
      skeleton = character(0),
      mcc = numeric(0),
      combined_ions = character(0),
      massql = character(0),
      stringsAsFactors = FALSE
    )
    if (!is.null(export)) {
      create_dir(dirname(export))
      utils::write.table(
        final_queries,
        file = export,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
    return(invisible(final_queries))
  }

  message("Generate all combinations of queries.")
  combinations <- prune_combinations_bnb(
    ions_list = ions_list,
    spectra = mia_spectra,
    max_ions = ions_max,
    dalton = dalton,
    ppm = ppm,
    decimals = decimals,
    intensity_min = intensity_min
  )
  rm(candidates, ions_list)

  all_combinations <- unlist(combinations, recursive = FALSE)
  names(all_combinations) <- rep(names(combinations), lengths(combinations))
  rm(combinations)

  message("Test the queries. (This is the longest step)")
  queries_results <- all_combinations |>
    perform_list_of_queries_progress(
      spectra = mia_spectra,
      dalton = dalton,
      ppm = ppm,
      show_progress = TRUE,
      intensity_min = intensity_min
    )
  names(queries_results) <- names(all_combinations)

  message("Evaluate the performance of each query using MCC.")
  normalized_skeleton <- gsub(
    pattern = "+",
    replacement = ".",
    x = mia_spectra$SKELETON,
    fixed = TRUE
  )
  skeleton_counts <- table(normalized_skeleton)

  final_queries <- summarize_best_queries(
    queries_results = queries_results,
    all_combinations = all_combinations,
    skeleton_counts = skeleton_counts,
    total_spectra = total_spectra,
    decimals = decimals,
    dalton = dalton,
    ppm = ppm,
    intensity_min = intensity_min
  )

  message("Export the best queries for further use.")
  if (!is.null(export)) {
    create_dir(dirname(export))
    utils::write.table(
      final_queries,
      file = export,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  invisible(final_queries)
}
