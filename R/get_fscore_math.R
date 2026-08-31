# Math / statistics engine for factor-score computation
#
# Pure computational helpers: delta-method corrections, scoring-matrix code,
# and the internal machinery that backs get_fs.lavaan() and augment_lav_predict().

sqrt_or_na <- function(x) {
  sqrt(replace(x, x < 0, NA_real_))
}

# Legacy per-row layout (augment_lav_predict() contract): scores followed by
# se, loadings and error terms, with the error terms in UPPER-triangular
# order and a trailing intercept block. Values come from the shared
# value-only engine fs_row_cols() (R/fs_indiv.R); relative to that engine
# only the error-term layout differs (upper-triangular re-slice of the same
# fsT, which for a symmetric matrix carries the same values).
augment_fs2 <- function(fs, fsL, fsT, fsb = NULL) {
  vals <- fs_row_cols(fs, fsL, fsT, fsb)
  q <- ncol(fsT)
  k_ld_end <- q + q * q
  k_ev_end <- k_ld_end + q * (q + 1L) / 2
  fs_vec <- c(
    vals[1L, seq_len(q)],
    vals[1L, (q + 1L):k_ld_end],
    fsT[upper.tri(fsT, diag = TRUE)]
  )
  if (!is.null(fsb)) {
    fs_vec <- c(fs_vec, vals[1L, (k_ev_end + 1L):ncol(vals)])
  }
  cbind(as.data.frame(fs), matrix(fs_vec, nrow = 1))
}

compute_lav_fs_matrices <- function(
  acov,
  psi = NULL,
  alpha = NULL,
  method = c("regression", "Bartlett")
) {
  method <- match.arg(method)
  if (method == "regression") {
    # Solve psi X = t(acov) instead of forming solve(psi): identical for
    # symmetric psi (all call sites pass a symmetric latent covariance)
    # without an explicit inverse.
    fsL <- diag(nrow(acov)) - t(solve(psi, t(acov)))
    fsT <- fsL %*% acov
    if (is.null(alpha)) {
      fsb <- NULL
    } else {
      fsb <- alpha - fsL %*% alpha
    }
  } else if (method == "Bartlett") {
    fsL <- diag(nrow(acov))
    fsT <- acov
    if (is.null(alpha)) {
      fsb <- NULL
    } else {
      fsb <- rep(0, nrow(acov))
    }
  }
  list(fsL = fsL, fsT = fsT, fsb = fsb)
}

create_fsT_names <- function(fs_names) {
  out <- outer(fs_names, Y = fs_names, FUN = paste, sep = "_")
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  out[] <- paste0("ecov_", out)
  diag(out) <- paste0("ev_", fs_names)
  out
}

create_fsL_names <- function(lv_names, fs_names) {
  out <- outer(lv_names, Y = fs_names, FUN = paste, sep = "_by_")
  t(out)
}

get_fs_mat_names <- function(lv_names, int = TRUE) {
  # Initialize data frame
  fs_names <- paste0("fs_", lv_names)
  se_names <- paste0(fs_names, "_se")
  ev_names <- create_fsT_names(fs_names)
  dimnames(ev_names) <- rep(list(fs_names), 2)
  ld_names <- create_fsL_names(lv_names, fs_names = fs_names)
  dimnames(ld_names) <- list(fs_names, lv_names)
  out <- list(
    fs = fs_names,
    se = se_names,
    ld = ld_names,
    ev = ev_names
  )
  if (int) {
    c(out, int = paste0("int_", fs_names))
  } else {
    out
  }
}

#' Obtain factor scores and related definition variables from
#' a `lavaan` object for 2S-PA analyses.
#'
#' The score columns are obtained by calling the
#' [lavaan::lavPredict()] function; the per-row definitions of those
#' scores (standard errors, loading/cross-loading matrix, error
#' variance-covariance matrix, intercepts) come from the canonical
#' `get_fs()` pattern blocks (the same `compute_fscore()` engine
#' [get_fs()] itself uses), so there is one source of truth for
#' `fsL`/`fsT`/SE shared with [get_fs()] and [fs_indiv()].
#'
#' @param lavobj A fitted [`lavaan::lavaan-class`] object
#' @param method A character string indicating the scoring method to use.
#'               Must be either `"regression"` or `"Bartlett"`. The score
#'               columns are taken from \code{\link[lavaan]{lavPredict}}()
#'               (which reproduces lavaan's `NA` convention for Bartlett
#'               factors with no observed indicator); this function
#'               intentionally does not support the `"mean"` (sum-score)
#'               method of [get_fs()].
#' @param drop_list_single logical. Should the results be unlisted
#'                         for single-group models?
#' @param ... Additional arguments passed to [lavaan::lavPredict()]
#' @return A `data.frame` containing the factor scores, the corresponding
#'         standard errors, the loadings and cross-loadings of the factor
#'         scores as indicators of the latent variables, the
#'         error variance-covariance matrix of the factor scores,
#'         and the measurement intercepts (legacy column layout:
#'         `fs_*`, `fs_*_se`, `<indicator>_by_fs_*`, `ev_*`/`ecov_*` in
#'         upper-triangular order, and `int_*`). The per-row values are
#'         identical to those returned by `fs_indiv(get_fs(lavobj, method =
#'         method))`; the standard-error columns share the canonical
#'         `fs_*_se` naming of `get_fs()`/`fs_indiv()`, and the only
#'         difference from them is the ordering of the `ev_*`/`ecov_*`
#'         columns (upper-triangular here, lower-triangular in `get_fs()`
#'         and `fs_indiv()`), together with the three character-matrix
#'         attributes below.
#'         In addition, three character matrices are added as attributes
#'         that can be used as input to `tspa_mx_model()`:
#' * `ld`: cross-loading matrix
#' * `ev`: error variance-covariance matrix
#' * `int`: measurement intercepts
#' @export
#' @examples
#' library(lavaan)
#' hs_model <- ' visual  =~ x1 + x2 + x3 '
#' fit <- cfa(hs_model,
#'            data = HolzingerSwineford1939,
#'            group = "school")
#' augment_lav_predict(fit)
augment_lav_predict <- function(
  lavobj,
  method = c("regression", "Bartlett"),
  drop_list_single = TRUE,
  ...
) {
  method <- match.arg(method)
  # Per-row fsL/fsT/fsb come from the canonical get_fs() pattern blocks
  # (one source of truth with get_fs()/fs_indiv()); the score columns are
  # lavaan::lavPredict() output, which preserves lavaan's NA convention for
  # Bartlett factors with no observed indicator. `...` is forwarded to
  # lavPredict() only.
  ngroups <- lavaan::lavInspect(lavobj, what = "ngroups")
  blocks_by_group <- get_fs_blocks.lavaan(
    lavobj,
    method = method,
    add_to_evfs = rep(0, ngroups)
  )
  fs_lst <- lavaan::lavPredict(
    lavobj,
    type = "lv",
    method = method,
    ...
  )
  if (ngroups == 1) {
    fs_lst <- list(fs_lst)
  }
  has_means <- lavaan::lavInspect(lavobj, what = "meanstructure")
  out <- vector("list", ngroups)
  if (ngroups > 1) names(out) <- names(fs_lst)
  # The column layout depends only on the latent names (identical across
  # groups), so it is computed once, outside the group loop.
  fs_matnames <- get_fs_mat_names(
    colnames(as.matrix(fs_lst[[1L]])), int = has_means
  )
  fs_matnames_flat <- fs_matnames
  fs_matnames_flat$ld <- c(fs_matnames_flat$ld)
  fs_matnames_flat$ev <- fs_matnames_flat$ev[upper.tri(fs_matnames_flat$ev,
                                                        diag = TRUE)]
  fs_colnames <- unlist(fs_matnames_flat)
  for (g in seq_len(ngroups)) {
    fs_g <- as.matrix(fs_lst[[g]])
    blocks <- blocks_by_group[[g]]
    n_cases <- max(unlist(lapply(blocks, function(b) max(b$case_idx))))
    fs_dat <- data.frame(
      matrix(
        NA,
        nrow = n_cases,
        ncol = length(fs_colnames),
        dimnames = list(NULL, fs_colnames)
      )
    )
    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      case_idx <- blk$case_idx
      # Positional assignment: fs_dat keeps the legacy column names, while
      # augment_fs2() fills the (name-free) values in the legacy layout.
      # The int block exists in fs_dat only when the fit has a (estimated)
      # mean structure; compute_fscore() always attaches an fsb (a zero
      # vector without mean structure), so pass NULL to suppress it otherwise.
      fs_dat[case_idx, ] <- augment_fs2(
        fs_g[case_idx, , drop = FALSE],
        fsL = blk$fsL,
        fsT = blk$fsT,
        fsb = if (has_means) blk$fsb else NULL
      )
    }
    out[[g]] <- fs_dat
  }
  if (drop_list_single && ngroups == 1) {
    out <- out[[1]]
  }
  attr(out, "ld") <- fs_matnames$ld
  attr(out, "ev") <- fs_matnames$ev
  attr(out, "int") <- fs_matnames$int
  out
}

#' Compute factor scores
#'
#' @param y An N x p matrix where each row is a response vector. If there
#'          is only one observation, it should be a matrix of one row.
#' @param lambda A p x q matrix of factor loadings.
#' @param theta A p x p matrix of unique variance-covariances.
#' @param psi A q x q matrix of latent factor variance-covariances.
#' @param nu A vector of length p of measurement intercepts.
#' @param alpha A vector of length q of latent means.
#' @param method A character string indicating the method for computing
#'               factor scores: `"regression"` (default), `"Bartlett"`, or
#'               `"mean"` (sum scores, i.e. the plain uncentered item means;
#'               see [get_fs()] for the full description).
#' @param center_y Logical indicating whether \code{y} should be
#'                 mean-centered. Default to \code{TRUE}. Ignored for
#'                 \code{method = "mean"}, whose scores are raw item means.
#' @param fs_matrices Logical indicating whether covariances of the error
#'                    portion of factor scores (\code{fsT}), factor score
#'                    loading matrix (\eqn{L}; \code{fsL}) and intercept vector
#'                    (\eqn{b}; \code{fsb}) should be returned.
#'                    The loading and intercept matrices are the implied
#'                    loadings and intercepts by the model when using the
#'                    factor scores as indicators of the latent variables.
#'                    If \code{TRUE}, these matrices will be added as
#'                    attributes.
#' @param acov Logical indicating whether the asymptotic covariance matrix
#'             of factor scores should be returned as an attribute.
#' @param sum_items For `method = "mean"` only: a named list mapping factor
#'                  names to the items included in each factor's sum score.
#'                  `NULL` (default) auto-derives the assignment from the
#'                  loadings. See [get_fs()] for the full description.
#'
#' @return An N x p matrix of factor scores.
#' @export
#'
#' @examples
#' library(lavaan)
#' fit <- cfa(" ind60 =~ x1 + x2 + x3
#'              dem60 =~ y1 + y2 + y3 + y4 ",
#'            data = PoliticalDemocracy)
#' fs_lavaan <- lavPredict(fit, method = "Bartlett")
#' # Using R2spa::compute_fscore()
#' est <- lavInspect(fit, what = "est")
#' fs_hand <- compute_fscore(lavInspect(fit, what = "data"),
#'                           lambda = est$lambda,
#'                           theta = est$theta,
#'                           psi = est$psi,
#'                           method = "Bartlett")
#' fs_hand - fs_lavaan  # same scores
compute_fscore <- function(
  y,
  lambda,
  theta,
  psi = NULL,
  nu = NULL,
  alpha = NULL,
  method = c("regression", "Bartlett", "mean"),
  center_y = TRUE,
  acov = FALSE,
  fs_matrices = FALSE,
  sum_items = NULL
) {
  method <- match.arg(method)
  # "mean" scores are raw (uncentered) item means by definition; center_y is
  # ignored for that method.
  center_y <- if (method == "mean") FALSE else center_y
  if (is.null(nu)) {
    nu <- colMeans(y)
  }
  if (is.null(alpha)) {
    alpha <- matrix(0, nrow = ncol(as.matrix(lambda)))
  }
  y1c <- t(as.matrix(y))
  if (center_y) {
    meany <- lambda %*% alpha + nu
    y1c <- y1c - as.vector(meany)
  }
  a_mat <- compute_a_from_mat(
    method,
    lambda = lambda,
    psi = psi,
    theta = theta,
    sum_items = sum_items
  )
  fs <- if (method == "mean") {
    t(a_mat %*% y1c)
  } else {
    t(a_mat %*% y1c + as.vector(alpha))
  }
  # Bartlett scores are undefined for a factor with no observed indicator
  # among `lambda`'s rows (its a-matrix row is the zero row); regression
  # scores stay defined for such factors through the cross-factor
  # covariances. NA mirrors lavaan::lavPredict()'s convention so
  # missing-data output stays comparable to it.
  if (method == "Bartlett" && ncol(fs) > 0) {
    no_item <- which(colSums(abs(as.matrix(lambda))) == 0)
    if (length(no_item) > 0) {
      fs[, no_item] <- NA
    }
  }
  if (acov) {
    if (method == "regression") {
      covy <- lambda %*% psi %*% t(lambda) + theta
      attr(fs, "acov") <-
        unclass(psi - a_mat %*% covy %*% t(a_mat))
    } else if (method %in% c("Bartlett", "mean")) {
      attr(fs, "acov") <-
        unclass(a_mat %*% theta %*% t(a_mat))
    }
  }
  if (fs_matrices) {
    attr(fs, "scoring_matrix") <- a_mat
    fsL <- unclass(a_mat %*% lambda)
    fs_names <- paste0("fs_", colnames(fsL))
    rownames(fsL) <- fs_names
    attr(fs, "fsL") <- fsL
    # Intercept of the score regressed on the uncentered latent:
    # fsb = E[fs] - fsL %*% alpha, consistent with all other methods.
    # For raw mean scores this is M * nu (the mean of the factor's item
    # intercepts); it equals E[fs] when there is no mean structure (alpha = 0).
    fsb <- if (method == "mean") {
      as.numeric(a_mat %*% nu)
    } else {
      as.numeric(alpha - fsL %*% alpha)
    }
    names(fsb) <- fs_names
    attr(fs, "fsb") <- fsb
    fsT <- a_mat %*% theta %*% t(a_mat)
    rownames(fsT) <- colnames(fsT) <- fs_names
    attr(fs, "fsT") <- fsT
  }
  fs
}

compute_fspars <- function(
  par,
  lavobj,
  method = c("regression", "Bartlett", "mean"),
  what = c("a", "evfs", "ldfs"),
  psi_override = NULL,
  frees = NULL,
  mats = NULL,
  sum_items = NULL
) {
  method <- match.arg(method)
  what <- match.arg(what)
  # Direct slot access; avoids lavInspect()'s per-call version check.
  ngrp <- lavobj@Data@ngroups
  # frees/mats may be pre-fetched by a caller that evaluates this repeatedly
  # over a perturbed `par` (correct_evfs()): the free mask and the base est
  # matrices are identical on every evaluation -- only `par`'s free values
  # move -- so the two lavInspect() file-reads are hoisted out of that loop.
  if (is.null(frees)) {
    frees <- lavInspect(lavobj, what = "free")
  }
  if (is.null(mats)) {
    mats <- lavInspect(lavobj, what = "est")
  }
  if (ngrp == 1) {
    frees <- list(frees)
    mats <- list(mats)
  }
  out <- vector("list", ngrp)
  mp <- lavobj@Data@Mp
  for (g in seq_len(ngrp)) {
    free <- frees[[g]]
    mat <- mats[[g]]
    # free[[l]] maps each est-matrix cell to its parameter index in `par`
    # (0 when fixed); the re-injection below is vectorized -- one
    # assignment per matrix instead of one scalar assignment (plus one
    # mask scan) per free cell on every Jacobian evaluation.
    for (l in seq_along(free)) {
      pos <- which(free[[l]] > 0)
      mat[[l]][pos] <- par[free[[l]][pos]]
    }
    if (!is.null(psi_override)) {
      mat$psi <- psi_override
    }
    pat <- mp[[g]]$pat
    if (is.null(pat)) {
      pat <- matrix(TRUE, nrow = 1, ncol = ncol(mat$theta))
    }
    num_mp <- nrow(pat)
    out[[g]] <- vector("list", num_mp)
    for (m in seq_len(num_mp)) {
      idx <- which(pat[m, ])
      # do.call() routes method positionally and lambda/psi/theta/idx by
      # their element names below -- do not drop the names from the list.
      a <- do.call(
        compute_a_from_mat,
        args = c(
          method,
          mat[c("lambda", "psi", "theta")],
          idx = list(idx),
          sum_items = list(sum_items)
        )
      )
      if (what == "a") {
        out[[g]][[m]] <- a
      } else if (what == "evfs") {
        out[[g]][[m]] <- a %*% mat$theta[idx, idx, drop = FALSE] %*% t(a)
      } else if (what == "ldfs") {
        out[[g]][[m]] <- a %*% mat$lambda[idx, , drop = FALSE]
      }
    }
    if (num_mp == 1) {
      out[[g]] <- out[[g]][[1]]
    }
  }
  out
}

compute_a <- function(
  par,
  lavobj,
  method = c("regression", "Bartlett"),
  psi_override = NULL,
  frees = NULL,
  mats = NULL
) {
  # Resolve to a single value before forwarding: compute_fspars() runs
  # match.arg() against its own (longer) method choices, which only accepts a
  # multi-element default when it equals that default verbatim -- this keeps
  # compute_a()'s shorter default working when callers omit `method`.
  method <- match.arg(method)
  compute_fspars(
    par,
    lavobj = lavobj,
    method = method,
    what = "a",
    psi_override = psi_override,
    frees = frees,
    mats = mats
  )
}

compute_a_from_mat <- function(
  method = c("regression", "Bartlett", "mean"),
  lambda,
  theta,
  psi = NULL,
  idx = NULL,
  sum_items = NULL
) {
  if (!is.null(idx)) {
    lambda <- lambda[idx, , drop = FALSE]
    theta <- theta[idx, idx, drop = FALSE]
  }
  method <- match.arg(method)
  if (method == "regression") {
    if (is.null(psi)) {
      stop("input of psi (latent covariance) is needed for regression scores")
    }
    compute_a_reg(lambda, theta = theta, psi = psi)
  } else if (method == "Bartlett") {
    compute_a_bartlett(lambda, theta = theta, psi = psi)
  } else if (method == "mean") {
    compute_a_mean(lambda, sum_items = sum_items)
  }
}

compute_a_reg <- function(lambda, theta, psi) {
  covy <- lambda %*% psi %*% t(lambda) + theta
  ginvcovy <- MASS::ginv(covy)
  tlam_invcov <- crossprod(lambda, ginvcovy)
  psi %*% tlam_invcov
}

compute_a_bartlett <- function(lambda, theta, psi = NULL) {
  ginvth <- MASS::ginv(theta)
  tlam_invth <- crossprod(lambda, ginvth)
  A <- tlam_invth %*% lambda
  if (qr(A)$rank < nrow(A)) {
    # A is singular when some factor has no observed indicator among the
    # rows of `lambda` (e.g. a missing-data pattern in which a factor's
    # indicators are all NA). Fall back to the Moore-Penrose
    # (minimum-norm) solution: it leaves the remaining factors' Bartlett
    # weights unchanged and gives the unscoreable factor(s) identically
    # zero weights. compute_fscore() turns those zero rows into NA scores,
    # mirroring lavaan::lavPredict()'s NA convention.
    a <- MASS::ginv(A) %*% tlam_invth
  } else {
    a <- solve(A, tlam_invth)
  }
  # ginv(A) %*% drops the dimnames that solve() preserves; restore the
  # latent-variable names so the score columns downstream stay named.
  rownames(a) <- rownames(tlam_invth)
  a
}

# Mean (sum-score) scoring matrix: q x p (factor x item); row k holds
# 1/|I_k| on the items assigned to factor k, so the scores are the raw
# (uncentered) item means, fs = M y. When sum_items is NULL the
# item -> sum assignment is auto-derived from the estimated loadings (an
# indicator must load on exactly one factor).
compute_a_mean <- function(lambda, sum_items = NULL) {
  lambda <- as.matrix(lambda)
  p <- nrow(lambda)
  q <- ncol(lambda)
  lv_names <- colnames(lambda)
  ind_names <- rownames(lambda)

  if (is.null(sum_items)) {
    nload <- rowSums(lambda != 0)
    if (any(nload > 1)) {
      stop(
        "The following indicator(s) load on more than one factor: ",
        paste(deparse(ind_names[nload > 1]), collapse = ", "),
        ". Specify which sum each belongs to via 'sum_items'.",
        call. = FALSE
      )
    }
    sum_items <- vector("list", q)
    names(sum_items) <- lv_names
    for (k in seq_len(q)) {
      items <- ind_names[lambda[, k] != 0]
      if (length(items) == 0) {
        stop(
          "Factor '", lv_names[k], "' has no items. ",
          "Specify the item-to-sum assignment via 'sum_items'.",
          call. = FALSE
        )
      }
      sum_items[[k]] <- items
    }
  } else {
    if (!is.list(sum_items) || is.null(names(sum_items))) {
      stop(
        "'sum_items' must be a named list mapping factor names to ",
        "item names.",
        call. = FALSE
      )
    }
    unknown_lv <- setdiff(names(sum_items), lv_names)
    if (length(unknown_lv) > 0) {
      stop(
        "Unknown factor name(s) in 'sum_items': ",
        paste(deparse(unknown_lv), collapse = ", "),
        ". Model factors are: ",
        paste(deparse(lv_names), collapse = ", "), ".",
        call. = FALSE
      )
    }
    missing_lv <- setdiff(lv_names, names(sum_items))
    if (length(missing_lv) > 0) {
      stop(
        "'sum_items' must cover all model factors; no items given for: ",
        paste(deparse(missing_lv), collapse = ", "), ".",
        call. = FALSE
      )
    }
    all_items <- unname(unlist(sum_items))
    unknown_item <- setdiff(all_items, ind_names)
    if (length(unknown_item) > 0) {
      stop(
        "Unknown item name(s) in 'sum_items': ",
        paste(deparse(unknown_item), collapse = ", "),
        ". Model indicators are: ",
        paste(deparse(ind_names), collapse = ", "), ".",
        call. = FALSE
      )
    }
    dup_items <- all_items[duplicated(all_items)]
    if (length(dup_items) > 0) {
      stop(
        "The following item(s) are assigned to more than one sum: ",
        paste(deparse(dup_items), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    # The user's list order is arbitrary; reorder it to the model's factor
    # order so the positional indexing below (zero-item check and M build)
    # always refers to model factor k. (The auto-derive branch is already
    # built in model order.)
    sum_items <- sum_items[match(lv_names, names(sum_items))]
    for (k in seq_len(q)) {
      if (length(sum_items[[k]]) == 0) {
        stop(
          "Factor '", lv_names[k], "' has no items in 'sum_items'.",
          call. = FALSE
        )
      }
    }
  }

  M <- matrix(
    0,
    nrow = q,
    ncol = p,
    dimnames = list(lv_names, ind_names)
  )
  for (k in seq_len(q)) {
    M[k, match(sum_items[[k]], ind_names)] <- 1 / length(sum_items[[k]])
  }
  M
}

correct_evfs <- function(
  fit,
  method = c("regression", "Bartlett"),
  psi_override = NULL
) {
  method <- match.arg(method)
  # Direct slot access; avoids lavInspect()'s per-call version check.
  ngrp <- fit@Data@ngroups
  est_raw <- lavInspect(fit, what = "est")
  frees_raw <- lavInspect(fit, what = "free")
  est_fits <- est_raw
  if (ngrp == 1) {
    est_fits <- list(est_fits)
  }
  outs <- vector("list", ngrp)
  # vcov(fit) is group-independent; hoist it out of the per-group loop.
  vc_fit <- vcov(fit)
  for (g in seq_len(ngrp)) {
    est_fit <- est_fits[[g]]
    p <- nrow(est_fit$psi)
    th <- est_fit$theta
    # One complex-step Jacobian over the FULL a matrix (p x c_col) instead of
    # one call per row: at p x fewer a-matrix evaluations and no lavInspect()
    # inside the per-evaluation compute_a() (est/free pre-fetched above).
    # lavaan flattens f's matrix output column-major (`dx[, p] <- Im(...)`), so
    # row i's Jacobian is the slice J[i + p*(0:(c_col-1)), ] -- entry-identical
    # to a per-row lav_func_jacobian_complex() (complex steps depend only on the
    # perturbed parameter and are linear in the output entries, so slicing the
    # stacked result reproduces each row Jacobian exactly).
    J <- lavaan::lav_func_jacobian_complex(
      function(x, fit, method, psi_override, frees, mats) {
        compute_a(x, lavobj = fit, method = method, psi_override = psi_override,
                  frees = frees, mats = mats)[[g]]
      },
      coef(fit),
      fit = fit,
      method = method,
      psi_override = psi_override,
      frees = frees_raw,
      mats = est_raw
    )
    c_col <- nrow(J) %/% p
    # tr(th Ji vc Jj') = sum((th Ji) * (Jj vc')) by tr(AB) = sum(A * t(B))
    # with vc symmetric: each (i, j) pair pays only the element-wise sum,
    # and the per-i / per-j products are hoisted out of the pair loop.
    Jrows <- function(k) J[k + p * (0:(c_col - 1)), , drop = FALSE]
    thJ <- lapply(seq_len(p), function(i) th %*% Jrows(i))
    Jvc <- lapply(seq_len(p), function(j) Jrows(j) %*% vc_fit)
    out <- matrix(nrow = p, ncol = p)
    for (j in seq_len(p)) {
      for (i in j:p) {
        out[i, j] <- sum(thJ[[i]] * Jvc[[j]])
        if (i > j) {
          out[j, i] <- out[i, j]
        }
      }
    }
    outs[[g]] <- out
  }
  outs
}

compute_evfs <- function(
  par,
  lavobj,
  method = c("regression", "Bartlett", "mean"),
  psi_override = NULL,
  sum_items = NULL
) {
  method <- match.arg(method)
  compute_fspars(
    par,
    lavobj = lavobj,
    method = method,
    what = "evfs",
    psi_override = psi_override,
    sum_items = sum_items
  )
}

compute_ldfs <- function(
  par,
  lavobj,
  method = c("regression", "Bartlett", "mean"),
  psi_override = NULL,
  sum_items = NULL
) {
  method <- match.arg(method)
  compute_fspars(
    par,
    lavobj = lavobj,
    method = method,
    what = "ldfs",
    psi_override = psi_override,
    sum_items = sum_items
  )
}

compute_grad_ld_evfs <- function(
  fit,
  method = c("regression", "Bartlett", "mean"),
  psi_override = NULL,
  sum_items = NULL
) {
  method <- match.arg(method)
  lavaan::lav_func_jacobian_complex(
    function(x, fit, method, psi_override, sum_items) {
      evfs <- compute_evfs(x, lavobj = fit, method = method,
                           psi_override = psi_override,
                           sum_items = sum_items)
      evfs_lower <- lapply(evfs, function(x) {
        x[lower.tri(x, diag = TRUE)]
      })
      c(
        unlist(compute_ldfs(x, lavobj = fit, method = method,
                            psi_override = psi_override,
                            sum_items = sum_items)),
        unlist(evfs_lower)
      )
    },
    coef(fit),
    fit = fit,
    method = method,
    psi_override = psi_override,
    sum_items = sum_items
  )
}

vcov_ld_evfs <- function(
  fit,
  method = c("regression", "Bartlett", "mean"),
  psi_override = NULL,
  sum_items = NULL
) {
  method <- match.arg(method)
  jac <- compute_grad_ld_evfs(fit, method = method,
                              psi_override = psi_override,
                              sum_items = sum_items)
  jac %*% vcov(fit) %*% t(jac)
}

compute_fsrel <- function(fit, method = c("regression", "Bartlett")) {
  method <- match.arg(method)
  # Direct slot access; avoids lavInspect()'s per-call version check.
  ngrp <- fit@Data@ngroups
  est_raw <- lavInspect(fit, what = "est")
  frees_raw <- lavInspect(fit, what = "free")
  sigmas <- lavInspect(fit, "implied")
  est_fits <- if (ngrp == 1) list(est_raw) else est_raw
  sigmas <- if (ngrp == 1) list(sigmas) else sigmas
  # vcov(fit) is group-independent; hoist it out of the per-group loop.
  vc_fit <- vcov(fit)
  # est/free hoisted into compute_a(): the Jacobian below re-evaluates it
  # once per free parameter, and lavInspect() re-reads lavaan's DESCRIPTION
  # on every call (same idiom as correct_evfs()).
  a <- compute_a(coef(fit), lavobj = fit, method = method,
                 frees = frees_raw, mats = est_raw)
  outs <- vector("list", ngrp)
  for (g in seq_len(ngrp)) {
    est_fit <- est_fits[[g]]
    lam <- est_fit$lambda
    psi <- est_fit$psi
    if (ncol(lam) > 1) {
      stop("reliability is only supported for unidimensional models.")
    }
    jac_a <- lavaan::lav_func_jacobian_complex(
      function(x, fit, method, frees, mats) {
        compute_a(x, lavobj = fit, method = method,
                  frees = frees, mats = mats)[[g]]
      },
      coef(fit),
      fit = fit,
      method = method,
      frees = frees_raw,
      mats = est_raw
    )
    va <- jac_a %*% vc_fit %*% t(jac_a)
    aa <- crossprod(a[[g]]) + va
    # tr(M aa) = sum(M * aa): aa is symmetric by construction, the implied
    # covariance is symmetric, and lam %*% psi %*% t(lam) is an outer
    # product (lam is n x 1, psi is 1 x 1).
    lam_psi_lamT <- tcrossprod(lam %*% psi, lam)
    outs[[g]] <- sum(lam_psi_lamT * aa) / sum(sigmas[[g]]$cov * aa)
  }
  outs
}
