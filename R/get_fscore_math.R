# Math / statistics engine for factor-score computation
#
# Pure computational helpers: delta-method corrections, scoring-matrix code,
# and the internal machinery that backs get_fs.lavaan() and augment_lav_predict().

sqrt_or_na <- function(x) {
  sqrt(ifelse(x >= 0, x, NA))
}

augment_fs2 <- function(fs, fsL, fsT, fsb = NULL) {
  fs_se <- sqrt_or_na(diag(fsT))
  fs_lds <- c(fsL)
  fs_evs <- fsT[upper.tri(fsT, diag = TRUE)]
  fs_vec <- c(fs_se, fs_lds, fs_evs)
  if (!is.null(fsb)) {
    fs_vec <- c(fs_vec, fsb)
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
    fsL <- diag(nrow(acov)) - acov %*% solve(psi)
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
  se_names <- paste0("se_", fs_names)
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
#' This function obtained the factor scores, standard errors,
#' loading matrix, and variance covariance matrix by calling
#' the [lavaan::lavPredict()] function.
#'
#' @param lavobj A fitted [`lavaan::lavaan-class`] object
#' @param method A character string indicating the scoring method to use.
#'               Must be either `"regression"` or `"Bartlett"`.
#' @param drop_list_single logical. Should the results be unlisted
#'                         for single-group models?
#' @param ... Additional arguments passed to [lavaan::lavPredict()]
#' @return A `data.frame` containing the factor scores, the corresponding
#'         standard errors, the loadings and cross-loadings of the factor
#'         scores as indicators of the latent variables, the
#'         error variance-covariance matrix of the factor scores,
#'         and the measurement intercepts.
#'         In addition, three character matrices are added as attributes
#'         that can be used as input to [tspa_mx_model()]:
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
  mp_lst <- lavobj@Data@Mp
  fs_lst <- lavaan::lavPredict(
    lavobj,
    type = "lv",
    method = method,
    acov = TRUE,
    ...
  )
  if (lavInspect(lavobj, what = "ngroups") == 1) {
    fs_lst <- list(fs_lst)
    attr(fs_lst, "acov") <- attr(fs_lst[[1]], "acov")
  }
  pars <- lavInspect(lavobj, what = "est", drop.list.single.group = FALSE)
  out <- vector("list", length = length(fs_lst))
  names(out) <- names(fs_lst)
  has_means <- lavInspect(lavobj, what = "meanstructure")
  for (g in seq_along(fs_lst)) {
    mp <- mp_lst[[g]]
    fs <- fs_lst[[g]]
    if (is.null(mp)) {
      case_idx <- list(seq_len(nrow(fs)))
      acov_g <- list(attr(fs_lst, "acov")[[g]])
      acov_rank <- 1
    } else {
      case_idx <- mp$case.idx
      # Somehow lavaan sort the `acov` output by the missing data pattern and
      # does not match the order of the missing pattern
      # So need to find the order first
      acov_g <- attr(fs_lst, "acov")[[g]]
      acov_rank <- rank(mp$id)
    }
    # Initialize empty data frame
    fs_matnames <- get_fs_mat_names(colnames(fs), int = has_means)
    fs_matnames_flat <- fs_matnames
    fs_matnames_flat$ld <- c(fs_matnames_flat$ld)
    fs_matnames_flat$ev <- fs_matnames_flat$ev[upper.tri(fs_matnames_flat$ev,
                                                         diag = TRUE)]
    fs_colnames <- unlist(fs_matnames_flat)
    fs_dat <- data.frame(
      matrix(
        NA,
        nrow = nrow(fs),
        ncol = length(fs_colnames),
        dimnames = list(NULL, fs_colnames)
      )
    )
    psi <- pars[[g]]$psi
    alpha <- pars[[g]]$alpha
    for (i in seq_along(case_idx)) {
      mat_idx <- acov_rank[i]
      fs_matrices <- compute_lav_fs_matrices(
        acov = acov_g[[mat_idx]],
        psi = psi,
        alpha = alpha,
        method = method
      )
      fs_dat[case_idx[[i]], ] <- augment_fs2(
        fs[case_idx[[i]], , drop = FALSE],
        fsL = fs_matrices$fsL,
        fsT = fs_matrices$fsT,
        fsb = fs_matrices$fsb
      )
    }
    out[[g]] <- fs_dat
  }
  if (drop_list_single && length(out) == 1) {
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
#' @param method A character string indicating the method for computing factor
#'               scores. Currently, only "regression" is supported.
#' @param center_y Logical indicating whether \code{y} should be mean-centered.
#'                 Default to \code{TRUE}.
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
  method = c("regression", "Bartlett"),
  center_y = TRUE,
  acov = FALSE,
  fs_matrices = FALSE
) {
  method <- match.arg(method)
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
  a_mat <- compute_a_from_mat(method, lambda = lambda, psi = psi, theta = theta)
  fs <- t(a_mat %*% y1c + as.vector(alpha))
  if (acov) {
    if (method == "regression") {
      covy <- lambda %*% psi %*% t(lambda) + theta
      attr(fs, "acov") <-
        unclass(psi - a_mat %*% covy %*% t(a_mat))
    } else if (method == "Bartlett") {
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
    fsb <- as.numeric(alpha - fsL %*% alpha)
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
  method = c("regression", "Bartlett"),
  what = c("a", "evfs", "ldfs")
) {
  method <- match.arg(method)
  what <- match.arg(what)
  # Direct slot access; avoids lavInspect()'s per-call version check.
  ngrp <- lavobj@Data@ngroups
  frees <- lavInspect(lavobj, what = "free")
  mats <- lavInspect(lavobj, what = "est")
  if (ngrp == 1) {
    frees <- list(frees)
    mats <- list(mats)
  }
  out <- vector("list", ngrp)
  mp <- lavobj@Data@Mp
  for (g in seq_len(ngrp)) {
    free <- frees[[g]]
    mat <- mats[[g]]
    free_list <- lapply(free, FUN = \(x) x[which(x > 0)])
    for (l in seq_along(free_list)) {
      for (i in free_list[[l]]) {
        mat[[l]][which(free[[l]] == i)] <- par[i]
      }
    }
    pat <- mp[[g]]$pat
    if (is.null(pat)) {
      pat <- matrix(TRUE, nrow = 1, ncol = ncol(mat$theta))
    }
    num_mp <- nrow(pat)
    out[[g]] <- vector("list", num_mp)
    for (m in seq_len(num_mp)) {
      idx <- which(pat[m, ])
      a <- do.call(
        compute_a_from_mat,
        args = c(method, mat[c("lambda", "psi", "theta")], idx = list(idx))
      )
      if (what == "a") {
        out[[g]][[m]] <- a
      } else if (what == "evfs") {
        out[[g]][[m]] <- a %*% mat$theta[idx, idx, drop = FALSE] %*% t(a)
      } else if (what == "ldfs") {
        out[[g]][[m]] <- a %*% mat$lambda[idx, , drop = FALSE]
      }
      if (num_mp == 1) {
        out[[g]] <- out[[g]][[1]]
      }
    }
  }
  out
}

compute_a <- function(par, lavobj, method = c("regression", "Bartlett")) {
  compute_fspars(par, lavobj = lavobj, method = method, what = "a")
}

compute_a_from_mat <- function(
  method = c("regression", "Bartlett"),
  lambda,
  theta,
  psi = NULL,
  idx = NULL
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
  solve(tlam_invth %*% lambda, tlam_invth)
}

correct_evfs <- function(fit, method = c("regression", "Bartlett")) {
  method <- match.arg(method)
  # Direct slot access; avoids lavInspect()'s per-call version check.
  ngrp <- fit@Data@ngroups
  est_fits <- lavInspect(fit, what = "est")
  if (ngrp == 1) {
    est_fits <- list(est_fits)
  }
  outs <- vector("list", ngrp)
  for (g in seq_len(ngrp)) {
    est_fit <- est_fits[[g]]
    p <- nrow(est_fit$psi)
    jac_a <- vector("list", length = p)
    for (i in seq_len(p)) {
      jac_a[[i]] <- lavaan::lav_func_jacobian_complex(
        function(x, fit, method) {
          compute_a(x, lavobj = fit, method = method)[[g]][i, ]
        },
        coef(fit),
        fit = fit,
        method = method
      )
    }
    out <- matrix(nrow = p, ncol = p)
    th <- est_fit$theta
    vc_fit <- vcov(fit)
    for (j in seq_len(p)) {
      for (i in j:p) {
        out[i, j] <- sum(diag(th %*% jac_a[[i]] %*% vc_fit %*% t(jac_a[[j]])))
        if (i > j) {
          out[j, i] <- out[i, j]
        }
      }
    }
    outs[[g]] <- out
  }
  outs
}

compute_evfs <- function(par, lavobj, method = c("regression", "Bartlett")) {
  compute_fspars(par, lavobj = lavobj, method = method, what = "evfs")
}

compute_ldfs <- function(par, lavobj, method = c("regression", "Bartlett")) {
  compute_fspars(par, lavobj = lavobj, method = method, what = "ldfs")
}

compute_grad_ld_evfs <- function(fit, method = c("regression", "Bartlett")) {
  method <- match.arg(method)
  lavaan::lav_func_jacobian_complex(
    function(x, fit, method) {
      evfs <- compute_evfs(x, lavobj = fit, method = method)
      evfs_lower <- lapply(evfs, function(x) {
        x[lower.tri(x, diag = TRUE)]
      })
      c(
        unlist(compute_ldfs(x, lavobj = fit, method = method)),
        unlist(evfs_lower)
      )
    },
    coef(fit),
    fit = fit,
    method = method
  )
}

vcov_ld_evfs <- function(fit, method = c("regression", "Bartlett")) {
  method <- match.arg(method)
  jac <- compute_grad_ld_evfs(fit, method = method)
  jac %*% lavaan::vcov(fit) %*% t(jac)
}

compute_fsrel <- function(fit, method = c("regression", "Bartlett")) {
  method <- match.arg(method)
  # Direct slot access; avoids lavInspect()'s per-call version check.
  ngrp <- fit@Data@ngroups
  est_fits <- lavInspect(fit, what = "est")
  sigmas <- lavInspect(fit, "implied")
  if (ngrp == 1) {
    est_fits <- list(est_fits)
    sigmas <- list(sigmas)
  }
  vc_fit <- vcov(fit)
  a <- compute_a(coef(fit), lavobj = fit, method = method)
  outs <- vector("list", ngrp)
  for (g in seq_len(ngrp)) {
    est_fit <- est_fits[[g]]
    lam <- est_fit$lambda
    psi <- est_fit$psi
    if (ncol(lam) > 1) {
      stop("reliability is only supported for unidimensional models.")
    }
    jac_a <- lavaan::lav_func_jacobian_complex(
      function(x, fit, method) {
        compute_a(x, lavobj = fit, method = method)[[g]]
      },
      coef(fit),
      fit = fit,
      method = method
    )
    va <- jac_a %*% vc_fit %*% t(jac_a)
    aa <- crossprod(a[[g]]) + va
    outs[[g]] <- sum(diag(lam %*% psi %*% t(lam) %*% aa)) /
      sum(diag(sigmas[[g]]$cov %*% aa))
  }
  outs
}
