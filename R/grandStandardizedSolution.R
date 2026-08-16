#' Grand Standardized Solution
#'
#' Grand standardized solution of a two-stage path analysis model.
#'
#' @param object An object of class lavaan.
#' @param model_list A list of string variable describing the structural path
#'                   model, in \code{lavaan} syntax.
#' @param se A Boolean variable. If TRUE, standard errors for the grand
#'                   standardized parameters will be computed.
#' @param acov_par An asymptotic variance-covariance matrix for a fitted
#'                 model object.
#' @param free_list A list of model matrices that indicate the position of
#'                  the free parameters in the parameter vector.
#' @param level The confidence level required.
#' @return A matrix of the standardized model parameters and standard errors.
#'
#' @importFrom stats pnorm qnorm
#' @importFrom utils tail
#' @importFrom lavaan vcov lavInspect lav_func_jacobian_complex
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
  if (is.null(model_list)) model_list <- lavTech(object, what = "est")
  # The nobs slot is a per-group list; unlist to the vector form that
  # lavInspect(what = "nobs") returns.
  ns <- unlist(object@Data@nobs)
  if (length(ns) == 1) ns <- NULL
  if (is.null(ns)) {
    message(
      "The grand standardized solution is equivalent to ",
      "lavaan::standardizedSolution() for a model with a single group."
    )
  }
  if (is.null(acov_par)) acov_par <- vcov(object)
  if (is.null(free_list)) free_list <- lavTech(object, what = "free")

  partable <- subset(lavInspect(object, what = "list"), op == "~")
  out <- partable[, c("lhs", "op", "rhs", "exo", "group",
                      "block", "label")]
  partable_beta <- lavTech(object, what = "partable", list.by.group = TRUE)

  # Get standardized betas
  if (is.null(ns)) {
    tmp_std_beta <- std_beta_est(model_list)
    all_beta_pos <- partable_beta[[1]]$beta
  } else {
    tmp_std_beta <- unlist(grand_std_beta_est(model_list, ns))
    group_names <- names(partable_beta)
    all_beta_pos <- sapply(group_names, function(x) {
      partable_beta[[x]]$beta
    })
  }
  beta_pos <- which(all_beta_pos != 0)
  out$est.std <- tmp_std_beta[beta_pos]

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
    tmp_acov_std_beta <- jac %*% acov_par %*% t(jac)
    out$se <- sqrt(diag(as.matrix(tmp_acov_std_beta[beta_pos, beta_pos])))
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
  diag(inv_s_eta) %*% beta %*% diag(s_eta)
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
    diag(inv_s_eta) %*% x %*% diag(s_eta)
  })
}

#' @rdname grand_standardized_solution
#' @export
grandStandardizedSolution <- grand_standardized_solution
