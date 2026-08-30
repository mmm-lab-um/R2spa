#' Grand Standardized Solution
#'
#' Grand standardized solution of a two-stage path analysis model.
#'
#' @details
#' For a **multigroup** fit the grand standardized solution is computed
#' from the grand (pooled) means and variance-covariances and is **not**
#' equivalent to [`lavaan::standardizedSolution()`], which matches it only
#' in the single-group case (the function emits that message when it
#' applies). When the input fit is an SE-corrected `tspa()` fit (produced
#' with `tspa(corrected_se = TRUE)`, attribute `tspa_corrected = TRUE`) —
#' or when `acov_par = vcov(corrected_fit)` is supplied — the reported
#' standard errors are the first-order corrected grand-standardized
#' standard errors, while the point estimates (`est.std`) are unchanged.
#'
#' Structural slopes may be free or user-fixed. A user-fixed slope is
#' reported alongside the free ones: its `est.std` is the user value rescaled
#' by the (grand) SD ratio, and its `se` is the first-order delta
#' approximation, both matching `lavaan::standardizedSolution()` (which also
#' reports a delta SE for fixed slopes). A structural regression whose outcome
#' or predictor is not part of the beta matrix (e.g. observed-on-observed) is
#' still rejected with an error. The delta-method SEs are propagated over
#' the free `beta`/`psi` (and `alpha`, single group) elements as marked by
#' the lavaan free-position matrices; a free parameter whose block matrix
#' carries no mark contributes no first-order term (in lavaan 0.7-x the
#' free-position matrices do not mark, e.g., some observed-variable
#' intercept blocks, whose omission only matters when their sampling
#' uncertainty affects the group means pooled into the grand covariance).
#'
#' @param object An object of class lavaan.
#' @param model_list A list of string variable describing the structural path
#'                   model, in \code{lavaan} syntax.
#' @param se A Boolean variable. If TRUE, standard errors for the grand
#'                   standardized parameters will be computed.
#' @param acov_par An asymptotic variance-covariance matrix for a fitted
#'                 model object; defaults to `vcov(object)`. Supplying the
#'                 covariance of an SE-corrected `tspa()` fit (produced
#'                 with `tspa(corrected_se = TRUE)`) carries the first-order
#'                 correction through to the grand-standardized standard
#'                 errors reported by this function.
#' @param free_list A list of model matrices that indicate the position of
#'                  the free parameters in the parameter vector.
#' @param level The confidence level required.
#' @return A data frame (class `lavaan.data.frame`) with one row per structural
#' (`~`) parameter in partable order: `lhs`, `op`, `rhs`, `exo`, `group`,
#' `block`, `label`, the grand standardized estimate `est.std`, and, when
#' `se = TRUE`, `se`, `z`, `pvalue`, `ci.lower`, and `ci.upper`.
#'
#' @importFrom stats pnorm qnorm
#' @importFrom utils tail
#' @importFrom lavaan vcov lav_func_jacobian_complex
#'
#' @seealso
#' - `vignette("Grand Standardized Coefficients", package = "R2spa")` for the grand-standardization workflow.
#'
#' @export
#'
#' @examples
#' library(lavaan)
#'
#' ## A single-group, two-factor example
#' mod1 <- '
#'    # latent variables
#'      ind60 =~ x1 + x2 + x3
#'      dem60 =~ y1 + y2 + y3 + y4
#'    # regressions
#'      dem60 ~ ind60
#' '
#' fit1 <- sem(model = mod1,
#'           data  = PoliticalDemocracy)
#' grand_standardized_solution(fit1)
#'
#' ## A single-group, three-factor example
#' mod2 <- '
#'     # latent variables
#'       ind60 =~ x1 + x2 + x3
#'       dem60 =~ y1 + y2 + y3 + y4
#'       dem65 =~ y5 + y6 + y7 + y8
#'     # regressions
#'       dem60 ~ ind60
#'       dem65 ~ ind60 + dem60
#' '
#' fit2 <- sem(model = mod2,
#'             data  = PoliticalDemocracy)
#' grand_standardized_solution(fit2)
#'
#' ## A multigroup, two-factor example
#' mod3 <- '
#'   # latent variable definitions
#'     visual =~ x1 + x2 + x3
#'     speed =~ x7 + x8 + x9
#'   # regressions
#'     visual ~ c(b1, b1) * speed
#' '
#' fit3 <- sem(mod3, data = HolzingerSwineford1939,
#'             group = "school",
#'             group.equal = c("loadings", "intercepts"))
#' grand_standardized_solution(fit3)
#'
#' ## A multigroup, three-factor example
#' mod4 <- '
#'   # latent variable definitions
#'     visual =~ x1 + x2 + x3
#'     textual =~ x4 + x5 + x6
#'     speed =~ x7 + x8 + x9
#'
#'   # regressions
#'     visual ~ c(b1, b1) * textual + c(b2, b2) * speed
#' '
#' fit4 <- sem(mod4, data = HolzingerSwineford1939,
#'             group = "school",
#'             group.equal = c("loadings", "intercepts"))
#' grand_standardized_solution(fit4)


grand_standardized_solution <- function(object, model_list = NULL,
                                      se = TRUE, acov_par = NULL,
                                      free_list = NULL, level = .95) {
  if (is.null(model_list)) model_list <- tsp_model_matrices(object)
  ns <- tsp_nobs(object)
  if (length(ns) == 1) ns <- NULL
  if (is.null(ns)) {
    message(
      "The grand standardized solution is equivalent to ",
      "lavaan::standardizedSolution() for a model with a single group."
    )
  }
  if (is.null(acov_par)) acov_par <- vcov(object)
  if (is.null(free_list)) free_list <- tsp_free_matrices(object)

  partable <- tsp_partable_read(object)
  positions <- tsp_partable_positions(object)
  i_struct <- which(partable$op == "~")
  out <- partable[i_struct, c("lhs", "op", "rhs", "exo", "group",
                              "block", "label")]
  out_positions <- positions[i_struct]
  if (!length(i_struct)) {
    stop("grand_standardized_solution(): the model contains no structural (",
         "'~') parameters, so there is nothing to standardize.",
         call. = FALSE)
  }

  # Get standardized betas
  if (is.null(ns)) {
    tmp_std_beta <- unlist(std_beta_est(model_list))
  } else {
    tmp_std_beta <- unlist(grand_std_beta_est(model_list, ns))
  }
  # The standardized estimates live in the per-group beta matrices in
  # column-major order, which does not in general follow the partable row
  # order (model-statement order). Free rows are matched to their matrix
  # position through the global free position (the partable `free` column and
  # the free-position matrix entries are the same bijection). A user-FIXED
  # slope has no free position, so it is anchored by (lhs, rhs) variable
  # identity instead: the beta matrix is the full lhs x rhs grid, so the cell
  # exists whether or not it is free. The column-major position of (lhs = row
  # r, rhs = col c) in an nrow x ncol matrix is (c - 1) * nrow + r. The two
  # anchors agree on every free cell, so this only changes fixed rows.
  beta_free <- free_list[which(names(free_list) == "beta")]
  size_beta <- nrow(beta_free[[1]]) * ncol(beta_free[[1]])
  bnames <- tsp_beta_names(object)
  out_idx <- vapply(seq_len(nrow(out)), function(i) {
    g <- partable$group[i_struct[i]]
    off <- (g - 1L) * size_beta
    fp <- which(beta_free[[g]] == out_positions[i])
    if (length(fp) == 1L) {
      return(off + fp)
    }
    rn <- bnames[[g]]$rnm
    cm <- bnames[[g]]$clm
    r <- match(partable$lhs[i_struct[i]], rn)
    c <- match(partable$rhs[i_struct[i]], cm)
    if (is.na(r) || is.na(c)) {
      return(0L)
    }
    off + (c - 1L) * nrow(beta_free[[g]]) + r
  }, integer(1))
  if (any(out_idx == 0L)) {
    stop("grand_standardized_solution(): a structural ('~') parameter is not ",
         "a beta-matrix slope (its outcome or predictor is not part of the ",
         "structural regression block); only free or fixed latent-variable ",
         "slopes are supported.")
  }
  out$est.std <- tmp_std_beta[out_idx]

  # Get SEs for the standardized betas
  if (se) {
    if (is.null(ns)) {
      free_beta_psi <- free_list[c("beta", "psi")]
      est <- .combine_est(model_list[c("beta", "psi")],
                          free = free_beta_psi)
      jac <- lav_func_jacobian_complex(
        function(x) std_beta_est(model_list, free_list = free_list, est = x),
        x = est
      )
      pos_par <- .combine_est(free_beta_psi, free = free_beta_psi)
    } else {
      free_beta_psi_alpha <- free_list[which(names(model_list) %in%
                                               c("beta", "psi", "alpha"))]
      est <- .combine_est(model_list[which(names(model_list) %in%
                                             c("beta", "psi", "alpha"))],
                          free = free_beta_psi_alpha)
      jac <- lav_func_jacobian_complex(
        function(x) {
          unlist(grand_std_beta_est(model_list,
            ns = ns,
            free_list = free_list, est = x
          ))
        },
        x = est
      )
      pos_par <- .combine_est(free_beta_psi_alpha,
                              free = free_beta_psi_alpha)
    }
    acov_par <- acov_par[pos_par, pos_par]
    # Only the diagonal of the structural-path block of jac %*% acov_par %*%
    # t(jac) is needed. Subset the Jacobian first (out_idx rows of the full
    # q^2 x p Jacobian), then use rowSums(A * B) == diag(A %*% t(B)) to avoid
    # forming the full q^2 x q^2 covariance matrix.
    jac_sub <- jac[out_idx, , drop = FALSE]
    out$se <- sqrt(rowSums((jac_sub %*% acov_par) * jac_sub))
    out$z <- out$est.std / out$se
    out$pvalue <- 2 * (1 - pnorm(abs(out$z)))
    ci <- out$est.std +
      out$se %o% qnorm(c((1 - level) / 2, 1 - (1 - level) / 2))
    out$ci.lower <- ci[, 1]
    out$ci.upper <- ci[, 2]
  }
  class(out) <- c("lavaan.data.frame", "data.frame")
  out
}


# Latent variances
veta <- function(beta, psi, gamma = NULL, cov_x = NULL) {
  inv_Imb <- solve(diag(nrow = nrow(beta)) - beta)
  if (!is.null(gamma) && !is.null(cov_x)) {
    psi_plus_vgammax <- psi + gamma %*% cov_x %*% t(gamma)
  } else {
    psi_plus_vgammax <- psi
  }
  inv_Imb %*% psi_plus_vgammax %*% t(inv_Imb)
}

.fill_matrix_list <- function(mod, free, est) {
  start_idx <- 0
  for (m in seq_along(mod)) {
    len_m <- sum(free[[m]] != 0)
    if (len_m == 0) next
    m_idx <- start_idx + seq_len(len_m)
    mod[[m]][free[[m]] != 0] <- est[m_idx]
    start_idx <- tail(m_idx, n = 1)
  }
  mod
}

std_beta_est <- function(model_list, free_list = NULL, est = NULL) {
  # The `est` argument is used to evaluate how changes in parameters
  # affect the standardized estimates, and is used to obtain the
  # derivatives for the delta method.
  if (!is.null(est) && !is.null(free_list)) {
    model_list <- .fill_matrix_list(model_list[c("beta", "psi")],
                                    free = free_list[c("beta", "psi")],
                                    est = est)
  }
  beta <- model_list$beta
  psi <- model_list$psi
  v_eta <- veta(beta, psi = psi)
  s_eta <- sqrt(diag(v_eta))
  inv_s_eta <- 1 / s_eta
  # diag(inv_s_eta) %*% beta %*% diag(s_eta) as O(q^2) recycling: a length-nrow
  # vector recycles down the columns (row scaling), rep(v, each = nrow) scales
  # the columns. Also immune to diag()'s length-1-numeric-as-size trap.
  (inv_s_eta * beta) * rep(s_eta, each = nrow(beta))
}

# Function for combining free estimates into a vector
.combine_est <- function(mod, free) {
  out <- vector("list", length(free))
  for (m in seq_along(free)) {
    out[[m]] <- mod[[m]][free[[m]] != 0]
  }
  unlist(out)
}

eeta <- function(beta, alpha, gamma = NULL, mean_x = NULL) {
  inv_Imb <- solve(diag(nrow = nrow(beta)) - beta)
  if (!is.null(gamma) && !is.null(mean_x)) {
    alpha_plus_gammax <- alpha + gamma %*% mean_x
  } else {
    alpha_plus_gammax <- alpha
  }
  inv_Imb %*% alpha_plus_gammax
}

veta_grand <- function(ns, beta_list, psi_list, alpha_list,
                       gamma_list = vector("list", length(beta_list)),
                       cov_x_list = vector("list", length(beta_list)),
                       mean_x_list = vector("list", length(beta_list))) {
  # Within-group variance-covariances
  vetas <- mapply(veta, beta = beta_list, psi = psi_list,
                  gamma = gamma_list, cov_x = cov_x_list,
                  SIMPLIFY = FALSE)
  # Group means
  eetas <- mapply(eeta, beta = beta_list, alpha = alpha_list,
                  gamma = gamma_list, mean_x = mean_x_list,
                  SIMPLIFY = FALSE)
  # Grand mean
  eeta_grand <- do.call(cbind, eetas) %*% ns / sum(ns)
  Reduce(
    `+`,
    mapply(function(v, m, n) n * (v + tcrossprod(m - eeta_grand)),
           v = vetas, m = eetas, n = ns, SIMPLIFY = FALSE)
  ) / sum(ns)
}

grand_std_beta_est <- function(model_list, ns, free_list = NULL, est = NULL) {
  if (!is.null(est) && !is.null(free_list)) {
    mat_idx <- which(names(model_list) %in% c("beta", "psi", "alpha"))
    model_list <- .fill_matrix_list(model_list[mat_idx],
                                    free = free_list[mat_idx],
                                    est = est)
  }
  beta_list <- model_list[which(names(model_list) == "beta")]
  psi_list <- model_list[which(names(model_list) == "psi")]
  alpha_list <- model_list[which(names(model_list) == "alpha")]
  v_eta <- veta_grand(ns,
                      beta_list,
                      psi_list = psi_list,
                      alpha_list = alpha_list)
  s_eta <- sqrt(diag(v_eta))
  inv_s_eta <- 1 / s_eta
  lapply(beta_list, function(x) {
    (inv_s_eta * x) * rep(s_eta, each = nrow(x))
  })
}

#' @rdname grand_standardized_solution
#' @export
grandStandardizedSolution <- grand_standardized_solution
