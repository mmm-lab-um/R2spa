# S3 methods for get_fs()
#
# Future methods (e.g. get_fs.mirt()) should be added to this file.

normalize_fs_method <- function(method) {
  method <- match.arg(
    method,
    c("regression", "Bartlett", "ML", "EB", "mean")
  )
  switch(method, ML = "Bartlett", EB = "regression", method)
}

validate_fs_priors <- function(prior_mean, prior_cov, lv_names) {
  q <- length(lv_names)

  if (!is.null(prior_mean)) {
    prior_nms <- names(prior_mean)
    prior_mean <- as.numeric(prior_mean)
    if (!is.null(prior_nms)) {
      if (!setequal(prior_nms, lv_names)) {
        stop(
          "'prior_mean' names must match the latent variable names (",
          paste(lv_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      prior_mean <- prior_mean[match(lv_names, prior_nms)]
    } else {
      if (length(prior_mean) != q) {
        stop(
          "'prior_mean' must have length ", q,
          " (the number of latent variables).",
          call. = FALSE
        )
      }
    }
    if (any(!is.finite(prior_mean))) {
      stop("'prior_mean' must contain only finite values.", call. = FALSE)
    }
    names(prior_mean) <- lv_names
  }

  if (!is.null(prior_cov)) {
    if (!is.matrix(prior_cov)) {
      prior_cov <- as.matrix(prior_cov)
    }
    if (any(!is.finite(prior_cov))) {
      stop("'prior_cov' must contain only finite values.", call. = FALSE)
    }
    if (nrow(prior_cov) != ncol(prior_cov)) {
      stop("'prior_cov' must be a square matrix.", call. = FALSE)
    }
    if (nrow(prior_cov) != q) {
      stop(
        "'prior_cov' must be a ", q, " x ", q, " matrix (one row and ",
        "column per latent variable).",
        call. = FALSE
      )
    }
    pcv <- prior_cov
    dimnames(pcv) <- NULL
    if (!isTRUE(all.equal(pcv, t(pcv)))) {
      stop("'prior_cov' must be symmetric.", call. = FALSE)
    }
    eig_vals <- eigen(pcv, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig_vals <= 0)) {
      stop("'prior_cov' must be positive definite.", call. = FALSE)
    }
    rn <- rownames(prior_cov)
    cn <- colnames(prior_cov)
    if (q == 1) {
      if (is.null(rn) && !is.null(cn)) {
        rn <- cn
      }
      if (is.null(cn) && !is.null(rn)) {
        cn <- rn
      }
    }
    if (!is.null(rn) || !is.null(cn)) {
      if (!setequal(rn, lv_names) || !setequal(cn, lv_names)) {
        stop(
          "'prior_cov' names must match the latent variable names (",
          paste(lv_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      prior_cov <- prior_cov[match(lv_names, rn),
                             match(lv_names, cn),
                             drop = FALSE]
    }
    rownames(prior_cov) <- colnames(prior_cov) <- lv_names
  }

  list(mean = prior_mean, cov = prior_cov)
}

#' @rdname get_fs
#' @export
get_fs.data.frame <- function(
  object,
  model = NULL,
  group = NULL,
  method = c("regression", "Bartlett", "ML", "EB", "mean"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  format = c("unified", "list"),
  prior_mean = NULL,
  prior_cov = NULL,
  sum_items = NULL,
  ...
) {
  if (is.null(model)) {
    ind_names <- colnames(object)
    if (!is.null(group)) {
      ind_names <- setdiff(ind_names, group)
    }
    model <- paste("f1 =~", paste(ind_names, collapse = " + "))
  }
  fit <- cfa(model, data = object, group = group, ...)
  get_fs(
    fit,
    method = method,
    corrected_fsT = corrected_fsT,
    vfsLT = vfsLT,
    reliability = reliability,
    format = format,
    prior_mean = prior_mean,
    prior_cov = prior_cov,
    sum_items = sum_items
  )
}

get_fs_blocks.lavaan <- function(
  object,
  method,
  add_to_evfs,
  prior_mean = NULL,
  prior_cov = NULL,
  sum_items = NULL,
  ...
) {
  method <- match.arg(method, c("regression", "Bartlett", "mean"))
  est <- lavInspect(object, what = "est")
  y <- lavInspect(object, what = "data")
  miss_pat <- object@Data@Mp

  prepare_fs <- function(y, est, add, mp, method) {
    psi_use <- if (is.null(prior_cov)) est$psi else prior_cov
    alpha_use <- if (is.null(prior_mean)) est$alpha else prior_mean
    if (is.null(mp)) {
      fscore <-
        compute_fscore(
          y,
          lambda = est$lambda,
          theta = est$theta,
          psi = psi_use,
          nu = est$nu,
          alpha = alpha_use,
          method = method,
          fs_matrices = TRUE,
          sum_items = sum_items
        )
      list(
        case_idx = seq_len(nrow(y)),
        fs = fscore,
        fsT = attr(fscore, "fsT") + add,
        fsL = attr(fscore, "fsL"),
        fsb = attr(fscore, "fsb"),
        scoring_matrix = attr(fscore, "scoring_matrix"),
        pat_label = paste0(colnames(y), collapse = "+"),
        pat = setNames(rep_len(TRUE, ncol(y)), colnames(y))
      )
    } else {
      npat <- mp$npatterns
      pats <- mp$pat
      mis_idx <- mp$case.idx
      blocks <- vector("list", npat)
      for (m in seq_len(npat)) {
        idx_m <- mis_idx[[m]]
        pat_m <- pats[m, ]
        fs_m <-
          compute_fscore(
            y[idx_m, pat_m, drop = FALSE],
            lambda = est$lambda[pat_m, , drop = FALSE],
            theta = est$theta[pat_m, pat_m, drop = FALSE],
            psi = psi_use,
            nu = est$nu[pat_m, , drop = FALSE],
            alpha = alpha_use,
            method = method,
            fs_matrices = TRUE,
            sum_items = sum_items
          )
        blocks[[m]] <- list(
          case_idx = idx_m,
          fs = fs_m,
          fsT = attr(fs_m, "fsT") + add,
          fsL = attr(fs_m, "fsL"),
          fsb = attr(fs_m, "fsb"),
          scoring_matrix = attr(fs_m, "scoring_matrix"),
          pat_label = paste0(colnames(y)[pat_m], collapse = "+"),
          pat = setNames(pat_m, colnames(y))
        )
      }
      blocks
    }
  }

  # Direct slot access avoids lavInspect()'s per-call version-compatibility
  # check (lav_object_check_version(), which re-reads lavaan's DESCRIPTION
  # file via read.dcf() every time) -- a major source of overhead when
  # get_fs() is called repeatedly. @Data@ngroups is already relied upon
  # elsewhere in this file (e.g. @Data@Mp, @Data@group.label).
  ngroups <- object@Data@ngroups

  if (ngroups == 1) {
    blocks <- prepare_fs(y, est, add_to_evfs[[1]], miss_pat[[1]], method)
    if (is.null(miss_pat[[1]])) {
      blocks <- list(blocks)
    }
    return(setNames(list(blocks), ""))
  }

  group_labels <- object@Data@group.label
  blocks_by_group <- setNames(vector("list", ngroups), group_labels)

  for (g in seq_len(ngroups)) {
    grp_est <- est[[g]]
    grp_y <- y[[g]]
    grp_add <- add_to_evfs[[g]]
    grp_mp <- miss_pat[[g]]
    blocks_by_group[[g]] <- prepare_fs(grp_y, grp_est, grp_add, grp_mp, method)
    if (is.null(grp_mp)) {
      blocks_by_group[[g]] <- list(blocks_by_group[[g]])
    }
  }

  blocks_by_group
}

#' @rdname get_fs
#' @param format Output format: `"unified"` returns a single data frame with
#'        a `group` column; `"list"` returns a list of data frames per group.
#' @export
get_fs.default <- function(object, ...) {
  if (is.matrix(object)) {
    object <- as.data.frame(object)
    return(get_fs(object, ...))
  }
  stop(
    "get_fs() is not implemented for objects of class '",
    paste(class(object), collapse = "', '"),
    "'. Currently supported: 'data.frame', 'lavaan', ",
    "and 'lmerMod'. Support for 'mirt' models is planned.",
    call. = FALSE
  )
}

#' @rdname get_fs
#' @param format Output format: `"unified"` returns a single data frame with
#'        a `group` column; `"list"` returns a list of data frames per group.
#' @export
get_fs.lavaan <- function(
  object,
  method = c("regression", "Bartlett", "ML", "EB", "mean"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  format = c("unified", "list"),
  prior_mean = NULL,
  prior_cov = NULL,
  sum_items = NULL,
  ...
) {
  method <- normalize_fs_method(method)
  if (!inherits(object, "lavaan")) {
    stop("`object` must be a `lavaan` model object.")
  }
  if (method == "mean") {
    # "mean" scores are raw item means: they bypass the corrected-FS-T /
    # vfsLT / reliability / prior machinery, which only supports
    # regression and Bartlett scores.
    # A FIML-style fit carries @Data@Mp even for complete data (a single
    # all-TRUE pattern), so missingness is read from the patterns, not from
    # the slot's presence.
    if (!all(vapply(object@Data@Mp, function(mp) {
      is.null(mp) || all(mp$pat)
    }, logical(1)))) {
      stop(
        "method = 'mean' does not support models fitted with missing ",
        "data; use method = 'regression' or 'Bartlett'.",
        call. = FALSE
      )
    }
    incompatible <- c(
      if (corrected_fsT) "corrected_fsT",
      if (vfsLT) "vfsLT",
      if (reliability) "reliability",
      if (!is.null(prior_mean)) "prior_mean",
      if (!is.null(prior_cov)) "prior_cov"
    )
    if (length(incompatible) > 0) {
      stop(
        "method = 'mean' is not supported together with: ",
        paste(incompatible, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  format <- match.arg(format)

  # The sampling-error SE paths (corrected_fsT/reliability/vfsLT) require
  # the full observed-indicator set, so they cannot be applied to lavaan's
  # per-missing-data-pattern blocks; fail here instead of deep inside
  # compute_fspars()/correct_evfs().
  has_miss_patterns <- any(vapply(
    object@Data@Mp,
    function(m) !is.null(m),
    logical(1)
  ))
  if (has_miss_patterns && (corrected_fsT || reliability || vfsLT)) {
    stop(
      "'corrected_fsT', 'reliability', and 'vfsLT' are not supported when ",
      "the data contain missing values: lavaan scores the cases on their ",
      "observed-indicator patterns, and these SE paths require the full ",
      "indicator set. Fit on complete data or re-run with ",
      "corrected_fsT = FALSE, reliability = FALSE, vfsLT = FALSE.",
      call. = FALSE
    )
  }

  has_priors <- !is.null(prior_mean) || !is.null(prior_cov)
  if (has_priors) {
    if (method == "Bartlett") {
      stop(
        "'prior_mean'/'prior_cov' are only supported for 'regression' ",
        "(EB) scoring.",
        call. = FALSE
      )
    }
    if (reliability) {
      stop(
        "'reliability = TRUE' is not supported with user-supplied ",
        "'prior_mean'/'prior_cov'.",
        call. = FALSE
      )
    }
    est_first <- lavInspect(object, what = "est")
    if (object@Data@ngroups > 1) {
      est_first <- est_first[[1]]
    }
    lv_names <- colnames(est_first$lambda)
    priors <- validate_fs_priors(prior_mean, prior_cov, lv_names)
  } else {
    priors <- list(mean = NULL, cov = NULL)
  }

  if (reliability) {
    corrected_fsT <- TRUE
  }
  if (corrected_fsT) {
    add_to_evfs <- correct_evfs(
      object,
      method = method,
      psi_override = priors$cov
    )
  } else {
    # Direct slot access; see note in get_fs_blocks.lavaan() -- avoids
    # lavInspect()'s expensive per-call version check.
    add_to_evfs <- rep(0, object@Data@ngroups)
  }

  blocks_by_group <- get_fs_blocks.lavaan(
    object,
    method = method,
    add_to_evfs = add_to_evfs,
    prior_mean = priors$mean,
    prior_cov = priors$cov,
    sum_items = sum_items
  )

  group_var <- object@Data@group
  group_col <- if (length(group_var) > 0) group_var else NULL

  out <- assemble_fs_blocks(
    blocks_by_group,
    format = format,
    group_col = group_col
  )

  # Group-level latent moments (effective / prior-adjusted), attached
  # post-assemble beside the vfsLT/reliability attachments below. The
  # values are the same psi_use/alpha_use the scoring above used
  # (see get_fs_blocks.lavaan(): prior_cov/prior_mean when supplied,
  # otherwise the per-group estimates; a named zero vector when the model
  # has no estimated mean structure). The shape mirrors fsT: unified ->
  # named list by group label; list -> direct attribute on each group data
  # frame plus a list-valued attribute on the outer list (multigroup).
  est_all <- lavInspect(object, what = "est", drop.list.single.group = FALSE)
  ngroups <- object@Data@ngroups
  group_labels <- object@Data@group.label
  if (length(group_labels) != ngroups) {
    # Single-group fits carry @Data@group.label as character(0); use the
    # same "" wrapper convention assemble_fs_blocks() applies to the
    # fsT/fsL attributes (psi/alpha must mirror the fsT shape).
    group_labels <- rep("", ngroups)
  }
  psi_g <- vector("list", ngroups)
  alpha_g <- vector("list", ngroups)
  names(psi_g) <- names(alpha_g) <- group_labels
  for (g in seq_len(ngroups)) {
    est_g <- est_all[[g]]
    psi_g[[g]] <- if (is.null(priors$cov)) est_g$psi else priors$cov
    if (!is.null(priors$mean)) {
      alpha_g[[g]] <- priors$mean
    } else if (is.null(est_g$alpha)) {
      # No (estimated) mean structure: the compute_fscore() zero-alpha
      # convention, named to match the score columns.
      alpha_g[[g]] <- setNames(
        rep(0, ncol(est_g$lambda)),
        colnames(est_g$lambda)
      )
    } else {
      a_est <- est_g$alpha
      alpha_g[[g]] <- setNames(
        as.numeric(a_est),
        if (!is.null(rownames(a_est))) {
          rownames(a_est)
        } else {
          colnames(est_g$lambda)
        }
      )
    }
  }
  if (format == "unified") {
    attr(out, "psi") <- psi_g
    attr(out, "alpha") <- alpha_g
  } else if (ngroups == 1) {
    attr(out, "psi") <- psi_g[[1L]]
    attr(out, "alpha") <- alpha_g[[1L]]
  } else {
    for (g in seq_len(ngroups)) {
      attr(out[[g]], "psi") <- psi_g[[g]]
      attr(out[[g]], "alpha") <- alpha_g[[g]]
    }
    attr(out, "psi") <- psi_g
    attr(out, "alpha") <- alpha_g
  }

  if (vfsLT) {
    attr(out, "vfsLT") <- vcov_ld_evfs(
      object,
      method = method,
      psi_override = priors$cov
    )
  }
  if (reliability) {
    est <- lavInspect(object, what = "est")
    group_labels <- object@Data@group.label
    ngroups <- object@Data@ngroups
    # Dimensionality from the model estimates, not from the output
    # attributes: format = "list" carries the SG fsb attribute as a bare
    # vector while "unified" wraps it in a per-group list, so the attr
    # shape cannot be tested uniformly. psi is q x q for the q latent
    # variables (SG: top-level est element; MG: per-group list, and all
    # groups share the latent structure, so group 1 suffices).
    n_latent <- if (ngroups > 1) nrow(est[[1]]$psi) else nrow(est$psi)
    multifactor <- n_latent > 1
    if (multifactor) {
      warning(
        "Computation of reliability for a multi-factor model is not ",
        "currently supported. "
      )
    } else {
      if (ngroups == 1) {
        is_std.lv <- all(est$psi == 1)
        attr(out, "reliability") <-
          compute_fsrel(object, method = method)[[1]]
      } else {
        is_std.lv <- all(unlist(lapply(est, function(x) x$psi)) == 1)
        rels <- compute_fsrel(object, method = method)
        # The norig slot is a per-group list; unlist to the vector form
        # that lavInspect(what = "norig") returns.
        group_n <- unlist(object@Data@norig)
        rels[ngroups + 1] <- sum(unlist(rels) * group_n / sum(group_n))
        attr(out, "reliability") <-
          setNames(rels, c(group_labels, "overall"))
      }
      if (!is_std.lv) {
        warning(
          "Computation of reliability may not be accurate when ",
          "the latent variables are not standardized. "
        )
      }
    }
  }
  out
}

get_fs_blocks.merMod <- function(
  object,
  method = c("EB", "ML"),
  legacy_names = FALSE,
  ...
) {
  method <- match.arg(method)
  num_re <- length(object@cnms[[1]])
  base_names <- paste0("u", seq_len(num_re) - 1)
  re_names <- if (legacy_names) paste0(base_names, "_eb") else base_names
  fs_names <- paste0("fs_", re_names)

  # Grouping factor for the first RE term. lme4 canonicalizes its levels
  # (sorted for atomic vectors, user-specified order for pre-factorized
  # factors) and EVERY structure read below follows that level order --
  # not the row order of the data: split() on a factor, the Z columns,
  # the b/u random-effect vector, and ranef() rows are all indexed by
  # level position. Names must therefore come from the levels, never
  # from first-appearance order in the data.
  f1 <- as.factor(object@flist[[1]])
  n_clus <- nlevels(f1)

  # Row indices per cluster, in level order (split() follows factor
  # levels, so block j corresponds to level j regardless of row order).
  case_idx <- split(seq_len(length(f1)), f1)

  # Random-effects design Z (sparse n x sum_p matrix from lme4, column-
  # major). First-term fold invariant -- Zden folds the first RE term's
  # columns onto the num_re coefficient columns of each row's own cluster:
  # (a) term-major layout: the first term's columns are the FIRST
  #     num_re * n_clus columns of getME("Z"); level j of its grouping
  #     factor occupies columns (j - 1) * num_re + seq_len(num_re) --
  #     the same layout the per-cluster slicing below relies on;
  # (b) first-term nonzeros have disjoint support per level (every
  #     nonzero lies in a row of its own level's cluster), so row i of
  #     Zden holds the first term's coefficients for i's own cluster;
  # (c) multi-term models make Z wider; only the first term's block is
  #     folded, matching the `[[1]]` convention used for b/cnms/flist
  #     below. The random design -- not the fixed design -- determines
  #     Kz: with Z != X the fixed-design code produces non-conformable
  #     products. Zden is n x num_re (tiny) instead of the
  #     n x (num_re * n_clus) dense matrix as.matrix() used to
  #     allocate (multi-GB on large fits).
  Zsp <- lme4::getME(object, "Z")
  if (!inherits(Zsp, "CsparseMatrix")) {
    stop("`lme4::getME(object, 'Z')` must return a column-compressed ",
         "sparse matrix (got: ", class(Zsp)[1], ").", call. = FALSE)
  }
  stopifnot(ncol(Zsp) >= num_re * n_clus)
  n1 <- Zsp@p[(num_re * n_clus) + 1L]   # nnz in the first p * G columns
  cc0 <- rep(seq_len(num_re * n_clus),
            diff(Zsp@p[seq_len(num_re * n_clus + 1L)]))
  # Fold safety net: the scatter below (m[cbind(i, j)] <- v) silently
  # keeps the LAST value for duplicate (i, j) pairs, so a future change
  # to lme4's Z column layout would mis-score silently instead of
  # erroring. Under the invariants above, a row has at most one nonzero
  # per folded column (its own level's block is the only block with
  # nnz in that row, and each folded column maps from exactly one raw
  # column within the level), so this fires only if the invariant
  # breaks; all-zero rows/levels simply yield no nnz and cannot trip
  # it. Cost: one vectorized O(nnz) pass -- negligible.
  if (any(duplicated(cbind(Zsp@i[seq_len(n1)], (cc0 - 1L) %% num_re)))) {
    stop(
      "internal error: first-term Z fold produced duplicate (row, ",
      "coefficient) scatter pairs; the lme4 Z column layout assumed by ",
      "get_fs_blocks.merMod() (first term's ", num_re * n_clus,
      " columns, level-major) no longer holds.",
      call. = FALSE
    )
  }
  Zden <- matrix(0, nrow = nrow(Zsp), ncol = num_re)
  Zden[cbind(Zsp@i[seq_len(n1)] + 1L, (cc0 - 1L) %% num_re + 1L)] <-
    Zsp@x[seq_len(n1)]
  s <- stats::sigma(object)

  if (method == "EB") {
    # Scaled random-effects covariance of the first RE term; the lme4-2.x
    # theta convention and the first-term-only contract are documented at
    # the definition of get_D().
    D <- get_D(object)
    # EB scores for the first term: getME("b") = crossprod(Lambdat, u),
    # bit-identical to ranef(object)[[1]] values but computed level-major
    # without ranef()'s per-term work (cheaper for multi-term models).
    b <- lme4::getME(object, "b")
    stopifnot(length(b) >= num_re * n_clus)
    u_b <- matrix(
      b[seq_len(num_re * n_clus)],
      nrow = n_clus,
      ncol = num_re,
      byrow = TRUE
    )
    rownames(u_b) <- levels(f1)
  } else {
    # ML (prior-free) scores: per cluster, the MLE of u_j with u treated as
    # fixed is the OLS fit of the cluster's fixed-effects-adjusted residuals
    # on its random-effects design block. Residuals are reconstructed as
    # y_j - X_j %*% beta (no offset()/weights support, same as the EB
    # scoring identity).
    y <- stats::model.response(stats::model.frame(object))
    X <- as.matrix(lme4::getME(object, "X"))
    beta <- lme4::fixef(object)
  }

  blocks <- vector("list", n_clus)

  for (j in seq_len(n_clus)) {
    idx <- case_idx[[j]]
    zj <- Zden[idx, seq_len(num_re), drop = FALSE]
    Kz <- crossprod(zj)

    if (method == "ML") {
      # u_hat_j = (Z'Z)^+ Z' r_j: Bartlett-analog (estimator uses no prior D);
      # fsL = I (score = u_j + (Z'Z)^+ Z' e_j), fsT = sigma^2 (Z'Z)^+ (Penrose:
      # (Z'Z)^+ Z'Z (Z'Z)^+ = (Z'Z)^+). ginv handles rank-deficient Z blocks
      # (e.g. a random slope on a within-cluster constant predictor) with the
      # minimum-norm solution.
      rj <- y[idx] - as.numeric(X[idx, , drop = FALSE] %*% beta)
      Gz <- MASS::ginv(Kz)
      fs_row <- t(Gz %*% crossprod(zj, rj))
      fsL_j <- diag(num_re)
      fsT_j <- s^2 * Gz
      scoring_matrix_j <- Gz %*% t(zj)
    } else {
      DKz <- D %*% Kz
      inv_W <- solve(DKz + diag(nrow(Kz)))
      fsL_j <- DKz - DKz %*% inv_W %*% DKz
      fsT_j <- s^2 * inv_W %*% DKz %*% D %*% t(inv_W)
      fs_row <- u_b[j, , drop = FALSE]
      scoring_matrix_j <- inv_W %*% D %*% t(zj)
    }

    colnames(fs_row) <- re_names

    # colnames = indicator/lv names, rownames = fs names (augment_fs convention)
    colnames(fsL_j) <- re_names
    rownames(fsL_j) <- fs_names
    attr(fs_row, "fsL") <- fsL_j

    rownames(fsT_j) <- colnames(fsT_j) <- fs_names

    # Scoring matrix: S_j %*% (y_j - X_j %*% beta) reproduces the scores
    # (EB / ranef for method "EB", per-cluster OLS for "ML"), where y_j/X_j
    # are the cluster's rows of the model frame and fixed-effects design.
    # See vignettes/scoring-matrices.Rmd.
    rownames(scoring_matrix_j) <- fs_names
    colnames(scoring_matrix_j) <- as.character(seq_len(nrow(zj)))

    blocks[[j]] <- list(
      case_idx = idx,
      fs = fs_row,
      fsL = fsL_j,
      fsT = fsT_j,
      fsb = NULL,
      scoring_matrix = scoring_matrix_j
    )
  }

  # Names in canonical level order (matches ranef() row names / @cnms).
  setNames(blocks, levels(f1))
}
# Reconstruct the SCALED random-effects covariance of the FIRST random-
# effects term of a merMod fit from its @theta parameters.
#
# lme4 >= 2.0 convention (read from the lme4 2.0.6 source; same scale as
# 1.x, which also stored the Cholesky of the scaled covariance):
#   - x@theta packs one block per RE term, in cnms() (formula) order; the
#     block for a term with p coefficients holds p * (p + 1) / 2 entries
#     and is the packed lower-triangular (column-major filled) Cholesky
#     factor L of the SCALED covariance D / sigma^2, where D =
#     VarCorr(x)[[term]] is the unscaled term covariance;
#   - equivalently VarCorr(x)[[term]] == sigma(x)^2 * tcrossprod(L). That
#     is exactly what lme4 2.x implements: mkVarCorr() splits @theta into
#     per-term blocks of length p * (p + 1) / 2 (the same idiom used
#     below), builds one "Covariance.us" object per block (whose
#     getLambda() fills the packed entries into lower.tri(column-major)),
#     and returns sc^2 * tcrossprod(L) for non-GLMM fits.
#   - lme4 2.x no longer attaches the "clen" (block lengths) attribute to
#     @theta -- it is missing even on single-term fits. A fallback parser
#     such as lme4::vec2mlist() therefore re-parses the whole mixed theta
#     as a SINGLE block of (sqrt(8L + 1) - 1) / 2 rows: exact only while
#     L is a triangular number (single term), a replacement-length warning
#     for a 2+1 term split, a mixed 3x3 block for 2+2. The split is thus
#     done explicitly from @cnms here, independent of that attribute.
#
# Returns the first term's p x p scaled covariance (p = the first term's
# number of coefficients), which equals VarCorr(x)[[1]] / sigma(x)^2 for
# LMMs. It must be the first term's, not a block-diagonal combination of
# all terms: get_fs_blocks.merMod() scores only the first term -- its
# clusters (flist[[1]]), Kz block (first-term Z fold) and scores
# (getME(x, "b"), term-major) all follow the flist[[1]] / cnms[[1]]
# convention, and the EB formulas consuming this D carry the explicit
# sigma^2 scale themselves (fsT_j = s^2 * (D Kz + I)^-1 D Kz D ...).
get_D <- function(object) {
  n1 <- length(object@cnms[[1L]])
  blk <- object@theta[seq_len(n1 * (n1 + 1L) / 2L)]
  L_D <- diag(nrow = n1)
  L_D[lower.tri(L_D, diag = TRUE)] <- blk
  tcrossprod(L_D)
}

#' @rdname get_fs
#' @param fsm Currently not used.
#' @param format Currently not used for `merMod` objects: the output is
#'        always a single data frame with one row per cluster (no `group`
#'        column).
#' @param legacy_names Logical. Random-effect score naming convention for
#'        `merMod` objects. `FALSE` (default) uses `fs_u0`-style names
#'        (`fs_u0`/`fs_u1`/..., with loadings `u0_by_fs_u0` and error terms
#'        `ev_fs_u0`, `ecov_fs_u1_fs_u0`). `TRUE` reproduces the pre-refactor
#'        `u0_eb`-style *column names* (`u0_eb`, `u0_by_u0_eb`, `ev_u0_eb`,
#'        `ecov_u0_eb_u1_eb`) in the legacy column order. Note the legacy
#'        output is name-compatible, not byte-identical, with the
#'        pre-refactor `get_fs_lmer()` result: it additionally carries
#'        score-error columns (`u0_eb_se`, ...), per-cluster `fsL`/`fsT`
#'        array attributes, a per-cluster `scoring_matrix` list attribute
#'        (see [get_fs()]), and has NULL row names (the pre-refactor
#'        output had none of these and used the ranef subject IDs as row
#'        names).
#' @export
get_fs.merMod <- function(
  object,
  method = c("EB", "ML"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  fsm = FALSE,
  format = c("unified", "list"),
  legacy_names = FALSE,
  ...
) {
  method <- match.arg(method)
  if (!inherits(object, "merMod")) {
    stop("`object` must be an `lmerMod` model object.", call. = FALSE)
  }
  prior_dots <- list(...)[c("prior_mean", "prior_cov")]
  if (!is.null(prior_dots$prior_mean) || !is.null(prior_dots$prior_cov)) {
    stop(
      "'prior_mean'/'prior_cov' are not supported for `lmerMod` objects.",
      call. = FALSE
    )
  }

  blocks <- get_fs_blocks.merMod(
    object,
    method = method,
    legacy_names = legacy_names
  )

  # NOTE: merMod does NOT route through assemble_fs_blocks(). The shared
  # assembler assumes one data row per individual case (nrow = n_cases),
  # filling a template DataFrame by case_idx. merMod's semantics produce
  # one row per cluster (nrow = n_clusters), because each cluster has a
  # single EB estimate shared across its cases. Forcing it through the
  # assembler would produce n_cases rows with repeated values, breaking
  # backward compatibility with ranef()-aligned output. This is an
  # intentional architectural exception.
  aug_list <- lapply(blocks, function(b) {
    augment_fs(b$fs, b$fsT)
  })
  out <- do.call(rbind, aug_list)
  rownames(out) <- NULL

  # Legacy `u<k>_eb`-style names are not produced by augment_fs() (which
  # unconditionally prefixes scores with "fs_" and emits `ecov` columns in
  # row/col iteration order), so translate the column names when requested.
  if (legacy_names) {
    colnames(out) <- rename_legacy_fs_cols(colnames(out))
  }

  # Per-cluster fsL/fsT as array attributes (merMod-specific: each cluster
  # has its own covariance matrix, unlike lavaan groups where attributes are
  # shared within a group). Names come from the first block's fsL.
  n_clus <- length(blocks)
  fsL_1 <- blocks[[1]]$fsL
  re_names <- colnames(fsL_1)
  fs_names <- rownames(fsL_1)
  fsL_arr <- array(
    0,
    dim = c(nrow(fsL_1), ncol(fsL_1), n_clus),
    dimnames = list(fs_names, re_names, names(blocks))
  )
  fsT_arr <- array(
    0,
    dim = c(nrow(fsL_1), ncol(fsL_1), n_clus),
    dimnames = list(fs_names, fs_names, names(blocks))
  )
  for (j in seq_len(n_clus)) {
    fsL_arr[,, j] <- blocks[[j]]$fsL
    fsT_arr[,, j] <- blocks[[j]]$fsT
  }
  attr(out, "fsL") <- fsL_arr
  attr(out, "fsT") <- fsT_arr

  # Group-level latent moments: the (shared) prior covariance of the first
  # random-effects term, with dimnames renamed to the score names
  # (re_names) so they align with the fsL column names; random effects are
  # mean zero, so alpha is the named zero vector.
  psi_re <- as.matrix(lme4::VarCorr(object)[[1L]])
  rownames(psi_re) <- colnames(psi_re) <- re_names
  attr(out, "psi") <- psi_re
  attr(out, "alpha") <- setNames(rep(0, length(re_names)), re_names)

  # Per-cluster scoring matrices as a named list (one p x n_j matrix per
  # cluster; list, not array, because cluster sizes can differ).
  attr(out, "scoring_matrix") <- setNames(
    lapply(blocks, function(b) b$scoring_matrix),
    names(blocks)
  )

  out
}

# Translate augment_fs() column names into the pre-refactor `u<k>_eb`-style
# names. augment_fs() unconditionally prefixes scores with "fs_" and emits
# `ecov` columns in matrix row/col iteration order; the legacy convention
# needs neither that prefix (scores/SEs/ev), nor a bare (no "_eb") indicator
# side on loadings, nor the ascending `ecov` order — hence the per-category
# rules below.
rename_legacy_fs_cols <- function(x) {
  vapply(
    x,
    function(nm) {
      if (grepl("_by_", nm, fixed = TRUE)) {
        # <ind>_eb_by_fs_<fs>_eb -> <ind>_by_<fs>_eb
        sub("^(u[0-9]+)_eb_by_fs_(u[0-9]+)_eb$", "\\1_by_\\2_eb", nm)
      } else if (grepl("^ecov_", nm)) {
        # ecov_fs_uA_eb_fs_uB_eb -> ecov_uA_eb_uB_eb, ascending in u<k>
        reorder_ecov_col(
          sub("^ecov_fs_(u[0-9]+)_eb_fs_(u[0-9]+)_eb$",
              "ecov_\\1_eb_\\2_eb",
              nm
          )
        )
      } else if (grepl("^ev_", nm)) {
        # ev_fs_u0_eb -> ev_u0_eb
        sub("^ev_fs_", "ev_", nm)
      } else {
        # scores (fs_u0_eb) and SEs (fs_u0_eb_se): drop the leading "fs_"
        sub("^fs_", "", nm)
      }
    },
    character(1),
    USE.NAMES = FALSE
  )
}

# Reorder `ecov_uA_eb_uB_eb` so the u<k> indices are ascending (A < B).
# Comparing the extracted integers (not the strings) keeps this correct for
# 10+ random effects as well.
reorder_ecov_col <- function(nm) {
  m <- regmatches(nm, regexec("^ecov_u([0-9]+)_eb_u([0-9]+)_eb$", nm))[[1]]
  if (length(m) != 3L) return(nm)
  a <- as.integer(m[2])
  b <- as.integer(m[3])
  if (a <= b) nm else paste0("ecov_u", b, "_eb_u", a, "_eb")
}

# ===========================================================================
# mirt (Item Response Theory) support
#
# get_fs() methods for mirt's S4 item-response fits. A fitted
# SingleGroupClass has mirt's DEFAULT unit-variance / zero-mean factor prior
# (psi = I, alpha = 0). The score for each observation is its EAP posterior
# mean; its EAP posterior covariance (Vpost) feeds the shared regression-form
# matrix engine compute_lav_fs_matrices() (R/get_fscore_math.R), giving the
# per-observation implied loading / error-covariance:
#   fsL_i = I - Vpost_i        (univariate: 1 - SE^2)
#   fsT_i = fsL_i %*% Vpost_i  (univariate: (1 - SE^2) * SE^2)
#   fsb   = 0                  (alpha = 0)
# Because Vpost_i varies per observation, fsL/fsT are attached as PER-ROW
# (list) attributes plus the `mirt_per_obs` marker so that fs_indiv()
# (resolve_per_obs(), R/fs_indiv.R) mints one block per row.
# ===========================================================================

# mirt is a Suggests-only dependency; guard every mirt:: call at entry with a
# clear, actionable message (never library()/require() a function body).
require_mirt <- function() {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop(
      "'mirt' is required to extract factor scores from a mirt model. ",
      "Install it with install.packages('mirt') and retry.",
      call. = FALSE
    )
  }
}

#' @rdname get_fs
#' @param format Currently not used for `mirt` objects: the output is always a
#'        single data frame with one row per observation (no `group` column).
#' @export
get_fs.SingleGroupClass <- function(object, format = c("unified", "list"), ...) {
  require_mirt()
  if (!inherits(object, "SingleGroupClass")) {
    stop("`object` must be a mirt `SingleGroupClass` model object.", call. = FALSE)
  }
  format <- match.arg(format)  # accepted but unused: single-group mirt -> one df

  q <- mirt::extract.mirt(object, "nfact")
  fn <- mirt::extract.mirt(object, "factorNames")
  if (length(fn) != q) {
    stop(
      "internal error: nfact (", q, ") does not match the number of factor ",
      "names (", length(fn), ").",
      call. = FALSE
    )
  }
  fs_names <- paste0("fs_", fn)

  # EAP posterior means (+ per-observation SEs). mirt re-adds completely-
  # missing rows here, so `full` has one row per observation with NA scores
  # for the rows it could not score.
  full <- mirt::fscores(object, full.scores = TRUE, full.scores.SE = TRUE)
  full <- as.data.frame(full)

  # Per-observation EAP posterior covariance (acov). This early-returns BEFORE
  # mirt re-adds completely-missing rows, so the list has one q x q matrix per
  # SCORABLE observation only (named by scorable-row index).
  acov <- mirt::fscores(object, full.scores = TRUE, return.acov = TRUE)
  n <- nrow(full)
  # Plain numeric score matrix (one row per observation); column subsetting on
  # a matrix (not a data frame) keeps the `drop` argument honoured.
  score_mx <- as.matrix(full[, fn])
  if (n < 1 || length(acov) < 1) {
    stop("mirt returned no factor scores; the fit has no scorable observations.")
  }

  # Unit-variance / zero-mean factor prior (mirt's default). NOTE:
  # object@Model$Theta is the quadrature NODE grid, not the factor covariance,
  # so psi/alpha are not read from the fit.
  psi <- diag(q)
  rownames(psi) <- colnames(psi) <- fn
  alpha <- setNames(rep(0, q), fn)

  # Reconcile row alignment. `full` (score/SE call) includes the
  # completely-missing rows (NA scores); `acov` skips them. extract.mirt(
  # "completely_missing") names those original-data positions, so the scorable
  # rows are their complement -- in which order acov[[k]] is the k-th scorable
  # row (verified: diag(Vpost_k) == SE_row^2 for every scorable row).
  cm <- mirt::extract.mirt(object, "completely_missing")
  if (is.null(cm)) cm <- integer(0)
  keep <- !seq_len(n) %in% cm
  scorsc <- which(keep)
  if (length(acov) != length(scorsc)) {
    stop(
      "internal error: mirt posterior covariances (", length(acov),
      ") do not match the number of scorable rows (", length(scorsc), ").",
      call. = FALSE
    )
  }

  # Per-row regression-form matrices (one source of truth with the lavaan /
  # merMod paths: compute_lav_fs_matrices() with psi = I, alpha = 0).
  fsL_list <- vector("list", n)
  fsT_list <- vector("list", n)
  naL <- matrix(NA_real_, q, q, dimnames = list(fs_names, fn))
  naT <- matrix(NA_real_, q, q, dimnames = list(fs_names, fs_names))
  for (k in seq_len(length(scorsc))) {
    i <- scorsc[k]
    Vpost_i <- as.matrix(acov[[k]])
    m_i <- compute_lav_fs_matrices(Vpost_i, psi, alpha, method = "regression")
    L_i <- m_i$fsL
    T_i <- m_i$fsT
    rownames(L_i) <- fs_names
    colnames(L_i) <- fn
    rownames(T_i) <- colnames(T_i) <- fs_names
    fsL_list[[i]] <- L_i
    fsT_list[[i]] <- T_i
  }
  # Completely-missing rows: R2spa's NA-row convention (all-NA per-row block,
  # keeping reference dimnames for the column-name resolver).
  for (i in setdiff(seq_len(n), scorsc)) {
    fsL_list[[i]] <- naL
    fsT_list[[i]] <- naT
  }
  fsb <- setNames(rep(0, q), fs_names)  # alpha = 0 => zero intercept (constant)

  # Per-row se / loadings / error terms -- the shared value-only engine
  # fs_row_cols() (R/fs_indiv.R), applied row by row. fsb is passed as NULL so
  # get_fs() itself emits no intercept columns (fs_indiv() may emit them via
  # include_intercept = TRUE using the attached fsb).
  # fs_row_cols() layout (no intercept here): se (q) + loadings (q^2) +
  # error terms (q*(q+1)/2). The score columns are carried separately in
  # scores_df, so this width excludes them.
  K <- q + q * q + q * (q + 1L) / 2L
  vals <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(length(scorsc))) {
    i <- scorsc[k]
    vals[i, ] <- fs_row_cols(
      as.data.frame(score_mx[i, , drop = FALSE]),
      fsL_list[[i]], fsT_list[[i]], NULL
    )[1L, , drop = FALSE]
  }

  # Assemble the canonical data frame. Column set + order is identical to what
  # fs_indiv(get_fs(m)) emits: fs_<fn> | fs_<fn>_se | <fn_j>_by_fs_<fn>
  # (q^2, column-major per latent) | ev_fs_<fn>/ecov_<a>_<b> (lower-tri
  # row-major, i-outer j<=i). No group/id columns.
  scores_df <- as.data.frame(score_mx)
  colnames(scores_df) <- fs_names
  se_nm <- paste0(fs_names, "_se")
  ld_nm <- c(create_fsL_names(fn, fs_names))
  ev_nm <- character(q * (q + 1L) / 2L)
  count <- 1L
  for (i in seq_len(q)) {
    for (j in seq_len(i)) {
      ev_nm[count] <- if (i == j) {
        paste0("ev_", fs_names[i])
      } else {
        paste0("ecov_", fs_names[i], "_", fs_names[j])
      }
      count <- count + 1L
    }
  }
  out <- as.data.frame(
    cbind(as.matrix(scores_df), vals),
    check.names = FALSE
  )
  colnames(out) <- c(fs_names, se_nm, ld_nm, ev_nm)

  # Per-row attributes + group-level latent moments + the per-obs marker that
  # fs_indiv()'s resolve_per_obs() dispatches on.
  attr(out, "fsT") <- fsT_list
  attr(out, "fsL") <- fsL_list
  attr(out, "fsb") <- fsb
  attr(out, "fs_pattern") <- list(label = seq_len(n), pat = NULL)
  attr(out, "psi") <- psi
  attr(out, "alpha") <- alpha
  attr(out, "mirt_per_obs") <- TRUE
  out
}

#' @rdname get_fs
#' @export
get_fs.MultipleGroupClass <- function(object, ...) {
  stop(
    "Multi-group mirt models are not supported by get_fs(). Fit or extract a ",
    "single group first (e.g. with mirt::extract.group()) and call get_fs() ",
    "on the resulting SingleGroupClass object.",
    call. = FALSE
  )
}
