# combine_fs() -- combine several get_fs() results into one block-diagonal
# multi-factor result (PLAN: combining stage-1 results across models / data
# subsets). Each input is one block; the combined fsL/fsT/psi are block
# diagonal across the inputs (every cross-block entry is exactly zero). The
# per-row machinery reuses resolve_fs_per_row() + fs_row_cols() + block_diag(),
# so the derived _se / _by_ / ev_ / ecov_ columns are the same values get_fs()
# itself emits.

#' Combine `get_fs()` results into a block-diagonal multi-factor result
#'
#' @description
#' `combine_fs()` combines several [get_fs()] results -- possibly from
#' different measurement models and different subsets of the data (for
#' example a `mirt` fit for some items together with a `lavaan` fit for the
#' rest) -- into a single [get_fs()]-shaped result. The combined `fsL`,
#' `fsT`, and `psi` attributes are **block-diagonal across the inputs**: each
#' input contributes its own latent block (a multi-latent input contributes a
#' `k x k` block with its internal structure intact), and every entry that
#' crosses from one input to another is exactly zero. The result is
#' downstream-transparent: it feeds [tspa()] directly (no explicit
#' `fsL`/`fsT` needed).
#'
#' @details
#' Each element of `fs_list` is one **block**, so the block structure is
#' exactly the list structure. A two-factor `mirt` result is one `2 x 2`
#' block (its two factors stay correlated, as estimated by `mirt`); two
#' separate single-factor fits are two `1 x 1` blocks. Both ways of supplying
#' the `mirt` factors therefore work unchanged.
#'
#' **Row alignment.** With `id = NULL` (the default) every input must carry
#' the same number of rows in the same order -- the "different items, same
#' observations" case. With `id = "<column>"` the rows are aligned on that
#' column and the output rows are the **union** of the inputs' `id` values in
#' first-appearance order; a row absent from input *k* gets an `NA` score and
#' an all-`NA` block for that input's latents (the same NA-row convention
#' `mirt` uses for unscorable observations).
#'
#' **Derived columns.** The `fs_<name>` scores and `<name>_se` columns are
#' joined on `id`. The `<name_j>_by_fs_<name_i>` implied-loadings and the
#' `ev_`/`ecov_` error columns are recomputed from the combined block-diagonal
#' `fsL`/`fsT`/`fsb`, so every cross-block loading and cross-block error
#' covariance is exactly zero (within-block values are the inputs' own).
#'
#' **`scoring_matrix`.** The combined `scoring_matrix` is a per-row `q x p`
#' matrix, populated from the item-level scoring matrices of the inputs that
#' carry one (a complete-data `lavaan` result). Inputs without an item-level
#' matrix -- `mirt` (no item model in the result) and `lme4` (a per-cluster
#' random-effects matrix, not an item matrix) -- leave their rows as `NA`.
#' The attribute is always present (with `p = 0` when no input contributes
#' item columns).
#'
#' **Why not base `cbind()`.** `cbind()` on a data frame drops the `fsL` /
#' `fsT` / `psi` attributes, so it cannot build the block-diagonal
#' measurement inputs; `combine_fs()` is a distinct operation from the
#' single-factor column-join that the [tspa()] `se_fs` derivation is built
#' around.
#'
#' **Scope (v1).** Single-group inputs only: a multi-group result (a list of
#' per-group data frames, or a unified frame with a `group_col`) is rejected.
#' `lme4` inputs are one row per cluster, so align them on the cluster id.
#'
#' @param fs_list A list of at least two [get_fs()] results (lavaan, mirt, or
#'        lme4; single- or multi-latent). Each element is one block (see
#'        Details). Alternatively, the results may be passed as separate
#'        arguments (cbind-style): `combine_fs(a, b, c)` is equivalent to
#'        `combine_fs(list(a, b, c))`.
#' @param ... Additional [get_fs()] results (cbind-style). Use only with a
#'        non-list first argument: `combine_fs(a, b, c)`, not
#'        `combine_fs(list(a), b)`.
#' @param id A character string naming a column present in every input, used
#'        to align the rows (see Details). Default `NULL` (same rows, same
#'        order required).
#' @param latent_names Optional character vector. When supplied it must equal
#'        the auto-derived global latent order (each input's latent names, in
#'        `fs_list` order); a mismatch is an error. When `NULL` (default) the
#'        names are derived from each input.
#' @return A single-group data frame with one row per observation (marked
#'        `per_obs = TRUE`). Columns, in order: the combined `fs_<name>`
#'        scores, the `<name>_se` standard errors, the
#'        `<name_j>_by_fs_<name_i>` implied loadings, and the `ev_`/`ecov_`
#'        error terms (all cross-block terms exactly zero). Attributes:
#'        per-row `fsL`/`fsT`/`fsb`/`scoring_matrix` (block diagonal),
#'        `fs_pattern`, the block-diagonal `psi`, the concatenated `alpha`,
#'        and `per_obs = TRUE`. The result feeds [tspa()] directly.
#'
#' @examples
#' library(lavaan)
#' fit1 <- cfa("vis =~ x1 + x2 + x3", data = HolzingerSwineford1939)
#' fit2 <- cfa("qua =~ x4 + x5 + x6", data = HolzingerSwineford1939)
#' # list form and the cbind-style form are equivalent:
#' comb <- combine_fs(list(get_fs(fit1), get_fs(fit2)))
#' stopifnot(identical(comb, combine_fs(get_fs(fit1), get_fs(fit2))))
#' head(comb, 3)
#' # The psi attribute is block-diagonal (a vis block and a qua block with a
#' # zero cross block); the cross-block fsL / fsT entries are zero too:
#' attr(comb, "psi")
#'
#' @export
combine_fs <- function(fs_list, ..., id = NULL, latent_names = NULL) {
  # Accept either one list of results or cbind-style separate arguments:
  # combine_fs(a, b, c) == combine_fs(list(a, b, c)). A data frame is itself
  # a list in R, so "a results list" means a list that is NOT a data frame
  # (a multi-group result also looks like this; it is rejected below on its
  # group_col attribute, before the per-input loop).
  if (is.list(fs_list) && !is.data.frame(fs_list)) {
    if (!is.null(attr(fs_list, "group_col"))) {
      stop(
        "combine_fs() (v1) supports single-group inputs only; the first ",
        "argument is a multi-group (list-of-groups) get_fs() result.",
        call. = FALSE
      )
    }
    if (...length() > 0L) {
      stop(
        "combine_fs(): pass either one list of get_fs() results or the ",
        "results as separate arguments, not both.",
        call. = FALSE
      )
    }
  } else {
    fs_list <- c(list(fs_list), list(...))
  }
  if (length(fs_list) < 2L) {
    stop(
      "combine_fs() needs at least two get_fs() results.",
      call. = FALSE
    )
  }
  K <- length(fs_list)

  # -- Resolve each input to per-row form and collect its latent names. ---
  resolved <- vector("list", K)
  lv_all <- character(0L)
  for (k in seq_len(K)) {
    x <- fs_list[[k]]
    if (is.list(x) && !is.data.frame(x)) {
      stop(
        "combine_fs() (v1) supports single-group inputs only; input ", k,
        " is a multi-group (list-of-groups) result.",
        call. = FALSE
      )
    }
    if (!is.data.frame(x)) {
      stop(
        "combine_fs() input ", k,
        " is not a data frame as returned by get_fs().",
        call. = FALSE
      )
    }
    if (!is.null(attr(x, "group_col"))) {
      stop(
        "combine_fs() (v1) supports single-group inputs only; input ", k,
        " is a multi-group result.",
        call. = FALSE
      )
    }
    resolved[[k]] <- tryCatch(
      resolve_fs_per_row(x),
      error = function(e) {
        stop(
          "combine_fs() input ", k,
          " is not a resolvable get_fs() result: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
    lv_all <- c(lv_all, colnames(resolved[[k]]$blocks[[1L]]$fsL))
  }
  dup <- lv_all[duplicated(lv_all)]
  if (length(dup) > 0L) {
    stop(
      "combine_fs(): duplicate latent name(s) across inputs: ",
      paste(unique(dup), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.null(latent_names) && !identical(as.character(latent_names), lv_all)) {
    stop(
      "combine_fs(): 'latent_names' does not match the auto-derived latent ",
      "order. Expected: ", paste(lv_all, collapse = ", "), ".",
      call. = FALSE
    )
  }
  q <- length(lv_all)
  fs_names <- paste0("fs_", lv_all)

  # -- Row alignment. ------------------------------------------------------
  if (is.null(id)) {
    n0 <- resolved[[1L]]$n
    if (any(vapply(resolved, function(r) r$n, integer(1L)) != n0)) {
      stop(
        "combine_fs() with id = NULL requires every input to have the same ",
        "number of rows (got ",
        paste(vapply(resolved, function(r) r$n, integer(1L)),
              collapse = ", "),
        "); pass 'id' to align rows with different observations.",
        call. = FALSE
      )
    }
    n <- n0
    pos <- lapply(seq_len(K), function(k) seq_len(n))
  } else {
    if (!is.character(id) || length(id) != 1L || !nzchar(id)) {
      stop("'id' must be a character string naming a column.", call. = FALSE)
    }
    idv <- vector("list", K)
    for (k in seq_len(K)) {
      if (!id %in% names(fs_list[[k]])) {
        stop("combine_fs(): input ", k, " has no column named '", id, "'.",
             call. = FALSE)
      }
      idv[[k]] <- as.character(fs_list[[k]][[id]])
    }
    all_ids <- unique(do.call(c, idv))
    n <- length(all_ids)
    pos <- lapply(idv, function(v) match(all_ids, v))
  }

  # -- Per-input latent counts and global row offsets. --------------------
  sk <- vapply(resolved, function(r) ncol(r$blocks[[1L]]$fsL), integer(1L))
  off <- c(0L, cumsum(sk))

  # -- Column names via a reference block-diagonal pair. ------------------
  refL <- matrix(0, q, q, dimnames = list(fs_names, lv_all))
  refT <- matrix(0, q, q, dimnames = list(fs_names, fs_names))
  nm <- fs_row_colnames(refL, refT)
  se_nm <- nm$se
  ld_nm <- nm$ld
  ev_nm <- nm$ev
  k_ev <- q * (q + 1L) / 2L

  # -- Per-row combined blocks and derived columns. -----------------------
  S <- matrix(NA_real_, n, q, dimnames = list(NULL, fs_names))
  L_list <- vector("list", n)
  T_list <- vector("list", n)
  B_list <- vector("list", n)
  se_mat <- matrix(NA_real_, n, q)
  ld_mat <- matrix(NA_real_, n, q * q)
  ev_mat <- matrix(NA_real_, n, k_ev)

  dummy <- data.frame(x = 0L)
  for (r in seq_len(n)) {
    Lk <- vector("list", K)
    Tk <- vector("list", K)
    Bk <- vector("list", K)
    for (k in seq_len(K)) {
      rk <- resolved[[k]]
      s_k <- sk[k]
      lv_k <- colnames(rk$blocks[[1L]]$fsL)
      pk <- pos[[k]][r]
      if (is.na(pk)) {
        Lk[[k]] <- matrix(
          NA_real_, s_k, s_k,
          dimnames = list(paste0("fs_", lv_k), lv_k)
        )
        Tk[[k]] <- matrix(
          NA_real_, s_k, s_k,
          dimnames = list(paste0("fs_", lv_k), paste0("fs_", lv_k))
        )
        Bk[[k]] <- setNames(rep(NA_real_, s_k), paste0("fs_", lv_k))
        next
      }
      S[r, (off[k] + 1L):(off[k] + s_k)] <-
        as.numeric(rk$scores[pk, seq_len(s_k), drop = FALSE])
      b <- rk$blocks[[rk$pattern_idx[pk]]]
      Lk[[k]] <- b$fsL
      Tk[[k]] <- b$fsT
      Bk[[k]] <- if (is.null(b$fsb)) {
        setNames(rep(NA_real_, s_k), paste0("fs_", lv_k))
      } else {
        b$fsb
      }
    }
    Lr <- block_diag(Lk)
    Tr <- block_diag(Tk)
    L_list[[r]] <- Lr
    T_list[[r]] <- Tr
    B_list[[r]] <- do.call(c, Bk)
    vals <- fs_row_cols(dummy, Lr, Tr, NULL)
    se_mat[r, ] <- vals[1L, seq_len(q), drop = FALSE]
    ld_mat[r, ] <- vals[1L, (q + 1L):(q + q * q), drop = FALSE]
    ev_mat[r, ] <- vals[1L, (q + q * q + 1L):(q + q * q + k_ev), drop = FALSE]
  }

  # -- scoring_matrix (item-level where the input carries one). ----------
  sm_k <- lapply(fs_list, combine_fs_sm)
  p_k <- vapply(sm_k, function(m) if (is.null(m)) 0L else ncol(m),
                integer(1L))
  p_total <- sum(p_k)
  col_off <- c(0L, cumsum(p_k))
  sm_list <- vector("list", n)
  for (r in seq_len(n)) {
    M <- matrix(NA_real_, q, p_total, dimnames = list(fs_names, NULL))
    for (k in seq_len(K)) {
      mk <- sm_k[[k]]
      if (is.null(mk) || nrow(mk) != sk[k]) next
      pk <- pos[[k]][r]
      if (is.na(pk)) next
      M[
        (off[k] + 1L):(off[k] + sk[k]), (col_off[k] + 1L):(col_off[k] + p_k[k])
      ] <- mk
    }
    sm_list[[r]] <- M
  }

  # -- Model-level psi (block diagonal) and alpha (concatenated). --------
  psi_k <- lapply(seq_len(K), function(k) {
    p <- combine_fs_psi(fs_list[[k]])
    s_k <- sk[k]
    lv_k <- colnames(resolved[[k]]$blocks[[1L]]$fsL)
    if (is.null(p) || nrow(p) != s_k) {
      matrix(NA_real_, s_k, s_k, dimnames = list(lv_k, lv_k))
    } else {
      as.matrix(p)
    }
  })
  alpha_all <- do.call(
    c,
    lapply(seq_len(K), function(k) {
      a <- combine_fs_alpha(fs_list[[k]])
      s_k <- sk[k]
      if (is.null(a) || length(a) != s_k) rep(NA_real_, s_k) else as.numeric(a)
    })
  )

  # -- Assemble. ----------------------------------------------------------
  big <- cbind(S, se_mat, ld_mat, ev_mat)
  colnames(big) <- c(fs_names, se_nm, ld_nm, ev_nm)
  out <- as.data.frame(big, check.names = FALSE)
  rownames(out) <- NULL
  attr(out, "fsL") <- L_list
  attr(out, "fsT") <- T_list
  attr(out, "fsb") <- B_list
  attr(out, "scoring_matrix") <- sm_list
  attr(out, "fs_pattern") <- list(label = seq_len(n), pat = NULL)
  attr(out, "psi") <- block_diag(psi_k)
  attr(out, "alpha") <- setNames(alpha_all, lv_all)
  attr(out, "per_obs") <- TRUE
  out
}

# Item-level scoring matrix for an input, if it carries a single (complete-
# data) one: a plain matrix or a length-1 list holding one. Per-pattern
# (missing data) and per-cluster (lme4) matrices are a different shape and
# yield NULL (the input's combined scoring rows stay NA).
combine_fs_sm <- function(x) {
  sm <- attr(x, "scoring_matrix")
  if (is.null(sm)) return(NULL)
  if (is.matrix(sm)) return(as.matrix(sm))
  if (is.list(sm) && length(sm) == 1L && is.matrix(sm[[1L]])) {
    return(as.matrix(sm[[1L]]))
  }
  NULL
}

# Normalise the model-level psi attribute to a matrix (a unified single-group
# result wraps it in a length-1 list; mirt carries a bare matrix).
combine_fs_psi <- function(x) {
  p <- attr(x, "psi")
  if (is.null(p)) return(NULL)
  if (is.list(p) && !is.matrix(p)) p <- p[[1L]]
  as.matrix(p)
}

# Normalise the model-level alpha attribute to a numeric vector.
combine_fs_alpha <- function(x) {
  a <- attr(x, "alpha")
  if (is.null(a)) return(NULL)
  if (is.list(a) && !is.numeric(a)) a <- a[[1L]]
  as.numeric(a)
}
