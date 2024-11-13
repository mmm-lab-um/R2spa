#' Get Factor Scores and the Corresponding Standard Error of Measurement
#' @param data A data frame containing indicators.
#' @param model An optional string specifying the measurement model
#'              in \code{lavaan} syntax.
#'              See \code{\link[lavaan]{model.syntax}} for more information.
#' @param group Character. Name of the grouping variable for multiple group
#'              analysis, which is passed to \code{\link[lavaan]{cfa}}.
#' @param method Character. Method for computing factor scores (options are
#'               "regression" or "Bartlett"). Currently, the default is
#'               "regression" to be consistent with
#'               \code{\link[lavaan]{lavPredict}}, but the Bartlett scores have
#'               more desirable properties and may be preferred for 2S-PA.
#' @param corrected_fsT Logical. Whether to correct for the sampling
#'                      error in the factor score weights when computing
#'                      the error variance estimates of factor scores.
#' @param vfsLT Logical. Whether to return the covariance matrix of `fsT`
#'              and `fsL`, which can be used as input for [vcov_corrected()]
#'              to obtain corrected covariances and standard errors for
#'              [tspa()] results. This is currently ignored.
#' @param reliability Logical. Whether to return the reliability of factor
#'                    scores.
#' @param ... additional arguments passed to \code{\link[lavaan]{cfa}}. See
#'            \code{\link[lavaan]{lavOptions}} for a complete list.
#' @return A data frame containing the factor scores (with prefix `"fs_"`),
#'         the standard errors (with suffix `"_se"`), the implied loadings
#'         of factor `"_by_"` factor scores, and the error variance-covariance
#'         of the factor scores (with prefix `"evfs_"`). The following are
#'         also returned as attributes:
#'         * `fsT`: error covariance of factor scores
#'         * `fsL`: loading matrix of factor scores
#'         * `fsb`: intercepts of factor scores
#'         * `scoring_matrix`: weights for computing factor scores from items
#'
#' @importFrom lavaan cfa sem lavInspect lavTech coef
#' @importFrom stats setNames
#'
#' @export
#'
#' @examples
#' library(lavaan)
#' get_fs(PoliticalDemocracy[c("x1", "x2", "x3")])
#'
#' # Multiple factors
#' get_fs(PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")],
#'        model = " ind60 =~ x1 + x2 + x3
#'                  dem60 =~ y1 + y2 + y3 + y4 ")
#'
#' # Multiple-group
#' hs_model <- ' visual  =~ x1 + x2 + x3 '
#' fit <- cfa(hs_model,
#'            data = HolzingerSwineford1939,
#'            group = "school")
#' get_fs(HolzingerSwineford1939, hs_model, group = "school")
#' # Or without the model
#' get_fs(HolzingerSwineford1939[c("school", "x4", "x5", "x6")],
#'        group = "school")

get_fs <- function(data, model = NULL, group = NULL,
                   method = c("regression", "Bartlett"),
                   corrected_fsT = FALSE,
                   vfsLT = FALSE,
                   reliability = FALSE,
                   ...) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (is.null(model)) {
    ind_names <- colnames(data)
    if (!is.null(group)) {
      ind_names <- setdiff(ind_names, group)
    }
    model <- paste("f1 =~",
                   paste(ind_names, collapse = " + "))
  }
  fit <- cfa(model, data = data, group = group, ...)
  get_fs_lavaan(lavobj = fit, method = method,
                corrected_fsT = corrected_fsT,
                vfsLT = vfsLT,
                reliability = reliability)
}

#' @inherit get_fs
#' @param lavobj A lavaan model object when using [get_fs_lavaan()].
#' @export
get_fs_lavaan <- function(lavobj,
                          method = c("regression", "Bartlett"),
                          corrected_fsT = FALSE,
                          vfsLT = FALSE,
                          reliability = FALSE) {
  if (!inherits(lavobj, "lavaan")) {
    stop("`lavobj` must be a `lavaan` model object.")
  }
  est <- lavInspect(lavobj, what = "est")
  y <- lavInspect(lavobj, what = "data")
  if (reliability) corrected_fsT <- TRUE
  if (corrected_fsT) {
    add_to_evfs <- correct_evfs(lavobj, method = method)
  } else {
    add_to_evfs <- rep(0, lavInspect(lavobj, what = "ngroups"))
  }
  miss_pat <- lavobj@Data@Mp
  prepare_fs_dat <- function(y, est, add, case_idx, mp) {
    if (is.null(mp)) {
      fscore <-
        compute_fscore(y,
                       lambda = est$lambda,
                       theta = est$theta,
                       psi = est$psi,
                       nu = est$nu,
                       alpha = est$alpha,
                       method = method,
                       fs_matrices = TRUE
        )
      augment_fs(fscore, attr(fscore, which = "fsT") + add)
    } else {
      fscore <- matrix(NA, nrow = nrow(y), ncol = ncol(est$psi))
      npat <- mp$npatterns
      pats <- mp$pat
      mis_idx <- mp$case.idx
      for (m in seq_len(npat)) {
        idx_m <- mis_idx[[m]]
        pat_m <- pats[m, ]
        fs_m <-
          compute_fscore(y[idx_m, pat_m, drop = FALSE],
                         lambda = est$lambda[pat_m, , drop = FALSE],
                         theta = est$theta[pat_m, pat_m, drop = FALSE],
                         psi = est$psi,
                         nu = est$nu[pat_m, , drop = FALSE],
                         alpha = est$alpha,
                         method = method,
                         fs_matrices = TRUE
          )
        fs_dat <- augment_fs(fs_m, attr(fs_m, which = "fsT") + add)
        if (m == 1) {
          fscore <- as.data.frame(
            matrix(NA, nrow = nrow(y), ncol = ncol(fs_dat)))
          attributes(fscore) <-
            c(attributes(fs_dat)[1:2],
              list(row.names = rownames(fscore)),
              attributes(fs_dat)[-(1:3)])
        }
        fscore[idx_m, ] <- fs_dat
      }
      fscore
    }
  }
  group <- lavInspect(lavobj, what = "group")
  if (length(group) == 0) {
    out <- prepare_fs_dat(y, est = est, add = add_to_evfs[[1]],
                          mp = miss_pat[[1]])
  } else {
    fs_lst <- setNames(
      vector("list", length = length(est)),
      lavobj@Data@group.label
    )
    for (g in seq_along(fs_lst)) {
      fs_lst[[g]] <- prepare_fs_dat(
        y[[g]], est = est[[g]], add = add_to_evfs[[g]],
        mp = miss_pat[[g]]
      )
      fs_lst[[g]][[group]] <- names(est[g])
    }
    attr_names <- setdiff(names(attributes(fs_lst[[1]])),
                          c("names", "class", "row.names", "col.names"))
    attr_lst <- rep(
      list(
        setNames(vector("list", length = length(est)), lavobj@Data@group.label)
      ),
      length(attr_names)
    )
    for (j in seq_along(attr_lst)) {
      for (g in seq_along(fs_lst)) {
        attr_lst[[j]][[g]] <- attr(fs_lst[[g]], which = attr_names[j])
      }
      attr(fs_lst, which = attr_names[j]) <- attr_lst[[j]]
    }
    out <- fs_lst
  }
  if (vfsLT) {
    attr(out, "vfsLT") <- vcov_ld_evfs(lavobj, method = method)
  }
  if (reliability) {
    multifactor <- ifelse(inherits(attr(out, "fsb"), "list"),
                          length(attr(out, "fsb")[[1]]) > 1,
                          length(attr(out, "fsb")) > 1)
    if (multifactor) {
      warning("Compution of reliability for a multi-factor model is not ",
              "currently supported. ")
    } else {
      if (length(group) == 0) {
        is_std.lv <- est$psi == 1
        attr(out, "reliability") <-
          compute_fsrel(lavobj, method = method)[[1]]
      } else {
        is_std.lv = all(unlist(lapply(est, function(x) x$psi)) == 1)
        rels <- compute_fsrel(lavobj, method = method)
        group_n <- lavInspect(lavobj, what = "norig")
        rels[g + 1] <- sum(unlist(rels) * group_n / sum(group_n))
        attr(out, "reliability") <-
          setNames(rels, c(lavobj@Data@group.label, "overall"))
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

augment_fs <- function(fs, fs_ev) {
  fs_se <- t(as.matrix(sqrt(diag(fs_ev))))
  # fs_se[is.nan(fs_se)] <- 0
  colnames(fs) <- paste0("fs_", colnames(fs))
  colnames(fs_se) <- paste0(colnames(fs_se), "_se")
  num_lvs <- ncol(fs_ev)
  fs_evs <- rep(NA, num_lvs * (num_lvs + 1) / 2)
  count <- 1
  for (i in seq_len(num_lvs)) {
    for (j in seq_len(i)) {
      fs_evs[count] <- fs_ev[i, j]
      if (i == j) {
        names(fs_evs)[count] <- paste0("ev_", rownames(fs_ev)[i])
      } else {
        names(fs_evs)[count] <- paste0("ecov_",
                                       rownames(fs_ev)[i], "_",
                                       colnames(fs_ev)[j])
      }
      count <- count + 1
    }
  }
  fsL <- attr(fs, "fsL")
  fs_names <- paste0("fs_", colnames(fsL))
  fs_lds <- lapply(seq_len(ncol(fsL)), function(i) {
    setNames(fsL[, i],
             paste(colnames(fsL)[i], fs_names, sep = "_by_"))
  })
  fs_lds <- unlist(fs_lds)
  fs_dat <- cbind(as.data.frame(fs), fs_se, t(as.matrix(fs_lds)),
                  t(as.matrix(fs_evs)))
  attr(fs_dat, "fsT") <- fs_ev
  attr(fs_dat, "fsL") <- fsL
  attr(fs_dat, "fsb") <- attr(fs, "fsb")
  attr(fs_dat, "scoring_matrix") <- attr(fs, "scoring_matrix")
  fs_dat
}

get_fsLT_lmer <- function(object, tidy = FALSE) {
  # Predictor matrices for each cluster
  Xlist <- split(data.frame(object@pp$X), object@flist[[1]])
  D <- get_D(object@theta)
  lapply(Xlist, FUN = function(xj) {
    get_fsLT(D, Kz = crossprod(as.matrix(xj)), s = stats::sigma(object),
             tidy = tidy)
  })
}
get_D <- function(theta) {
    L_D <- lme4::vec2mlist(theta, symm = FALSE)[[1]]
    tcrossprod(L_D)
}
get_fsLT <- function(D, Kz, s, tidy) {
    DKz <- D %*% Kz
    inv_W <- solve(DKz + diag(nrow(Kz)))
    fsL <- DKz - DKz %*% inv_W %*% DKz
    fsT <- s^2 * inv_W %*% DKz %*% D %*% t(inv_W)
    if (tidy) {
        return(c(fsL, fsT[lower.tri(fsT, diag = TRUE)]))
    }
    list(fsL = fsL, fsT = fsT)
}

#' Get Factor Scores and the Corresponding Scoring Matrices for
#' Mixed-Effect Models
#' @param object A fiited model object of class [lme4::lmerMod-class].
#' @param method Currently only `"EB"` for empirical Bayes.
#' @param corrected_fsT Currently not used.
#' @param vfsLT Currently not used.
#' @param fsm Currently not used.
#'
#' @importFrom lme4 getME ranef
#' @importFrom Matrix crossprod solve t tcrossprod
#' 
#' @export
get_fs_lmer <- function(object,
                        method = c("EB"),
                        corrected_fsT = FALSE,
                        vfsLT = FALSE,
                        fsm = FALSE) {
  if (!inherits(object, "lmerMod")) {
    stop("`object` must be a `lmerMod` object.")
  }
  # TODO: check number of columns
  # tilde_e <- object@resp$y - as.vector(object@pp$X %*% object@beta)
  # Zt <- getME(object, "Zt")
  # Sigma <- Matrix::crossprod(getME(object, "Lambdat"))
  # V <- Matrix::crossprod(getME(object, "A")) + diag(length(tilde_e))
  # score_mat_eb <- Sigma %*% Zt %*% Matrix::solve(V)
  # fsL_eb <- score_mat_eb %*% Matrix::t(Zt)
  # fsT_eb <- Matrix::tcrossprod(score_mat_eb) * stats::sigma(object)^2
  num_re <- length(object@cnms[[1]]) # Number of random effect
  # num_clus <- nlevels(object@flist[[1]])
  # fsL_arr_eb <- fsT_arr_eb <-
  #   array(dim = c(num_re, num_re, num_clus))
  # for (i in seq_len(num_clus)) {
  #   idx <- seq_len(num_re) + num_re * (i - 1)
  #   fsL_arr_eb[, , i] <- as.matrix(fsL_eb[idx, idx])
  #   fsT_arr_eb[, , i] <- as.matrix(fsT_eb[idx, idx])
  # }
  fsLT <- get_fsLT_lmer(object, tidy = TRUE)
  fs_dat <- cbind(lme4::ranef(object)[[1]], do.call(rbind, fsLT))
  fs_names <- paste0("u", seq_len(num_re) - 1, "_eb")
  lv_names <- paste0("u", seq_len(num_re) - 1)
  fsT_names <- create_fsT_names(fs_names)
  colnames(fs_dat) <- c(
    fs_names,
    create_fsL_names(lv_names, fs_names = fs_names),
    fsT_names[lower.tri(fsT_names, diag = TRUE)]
  )
  # fsL_cols <- t(apply(fsL_arr_eb, MARGIN = 3, FUN = c))
  # fsT_cols <- t(apply(fsT_arr_eb,
  #   MARGIN = 3,
  #   FUN = \(x) x[lower.tri(x, diag = TRUE)]
  # ))
  # lcount <- tcount <- 1
  # fsL_names <- rep("", ncol(fsL_cols))
  # fsT_names <- rep("", ncol(fsT_cols))
  # for (j in seq_along(lv_names)) {
  #   for (i in seq_along(fs_names)) {
  #     fsL_names[lcount] <- paste(lv_names[j], fs_names[i], sep = "_by_")
  #     lcount <- lcount + 1
  #     if (i <= j) {
  #       if (i == j) {
  #         fsT_names[tcount] <- paste0("ev_", fs_names[i])
  #       } else {
  #         fsT_names[tcount] <- paste0(
  #           "ecov_",
  #           fs_names[i], "_",
  #           fs_names[j]
  #         )
  #       }
  #       tcount <- tcount + 1
  #     }
  #   }
  # }
  # colnames(fsL_cols) <- fsL_names
  # colnames(fsT_cols) <- fsT_names
  # fs_dat <- cbind(
  #   t(matrix(score_mat_eb %*% tilde_e, nrow = num_re)),
  #   fsL_cols, fsT_cols
  # )
  # colnames(fs_dat)[seq_len(num_re)] <- fs_names
  # if (fsm) {
  #   attr(fs_dat, "fsT") <- fsT_eb
  #   attr(fs_dat, "fsL") <- fsL_eb
  #   attr(fs_dat, "scoring_matrix") <- score_mat_eb
  # }
  fs_dat
}

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
  acov, psi = NULL, alpha = NULL,
  method = c("regression", "Bartlett")) {
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
  return(list(fsL = fsL, fsT = fsT, fsb = fsb))
}

create_fsT_names <- function(fs_names) {
  out <- outer(fs_names,
    Y = fs_names,
    FUN = paste, sep = "_"
  )
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  out[] <- paste0("ecov_", out)
  diag(out) <- paste0("ev_", fs_names)
  out
}

create_fsL_names <- function(lv_names, fs_names) {
  out <- outer(lv_names,
    Y = fs_names, FUN = paste,
    sep = "_by_"
  )
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
    fs = fs_names, se = se_names, ld = ld_names, ev = ev_names
  )
  if (int) {
    return(c(out, int = paste0("int_", fs_names)))
  } else {
    return(out)
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
    lavobj, method = c("regression", "Bartlett"),
    drop_list_single = TRUE, ...) {
  method <- match.arg(method)
  mp_lst <- lavobj@Data@Mp
  fs_lst <- lavaan::lavPredict(
    lavobj,
    type = "lv", method = method,
    acov = TRUE,
    ...
  )
  if (lavInspect(lavobj, what = "ngroups") == 1) {
    fs_lst <- list(fs_lst)
    attr(fs_lst, "acov") <- attr(fs_lst[[1]], "acov")
  }
  pars <- lavInspect(lavobj, what = "est",
                     drop.list.single.group = FALSE)
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
    fs_matnames <- get_fs_mat_names(colnames(fs),
                                    int = has_means)
    fs_colnames <- unlist(
      within(fs_matnames, expr = {
        ld <- c(ld)
        ev <- ev[upper.tri(ev, diag = TRUE)]
      })
    )
    fs_dat <- data.frame(
      matrix(NA,
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
compute_fscore <- function(y, lambda, theta, psi = NULL,
                           nu = NULL, alpha = NULL,
                           method = c("regression", "Bartlett"),
                           center_y = TRUE,
                           acov = FALSE,
                           fs_matrices = FALSE) {
  method <- match.arg(method)
  if (is.null(nu)) nu <- colMeans(y)
  if (is.null(alpha)) alpha <- matrix(0, nrow = ncol(as.matrix(lambda)))
  y1c <- t(as.matrix(y))
  if (center_y) {
    meany <- lambda %*% alpha + nu
    y1c <- y1c - as.vector(meany)
  }
  a_mat <- compute_a_from_mat(method,
                              lambda = lambda, psi = psi, theta = theta)
  fs <- t(a_mat %*% y1c + as.vector(alpha))
  if (acov) {
    # if (is.null(psi)) {
    #   stop("input of psi (latent covariance) is needed for acov")
    # }
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
    # tv <- fsL %*% psi %*% t(fsL)
    # fsv <- a_mat %*% covy %*% t(a_mat)
    # attr(fs, "fsT") <- fsv - tv
    fsT <- a_mat %*% theta %*% t(a_mat)
    rownames(fsT) <- colnames(fsT) <- fs_names
    attr(fs, "fsT") <- fsT
  }
  fs
}

compute_fspars <- function(par, lavobj, method = c("regression", "Bartlett"),
                           what = c("a", "evfs", "ldfs")) {
  method <- match.arg(method)
  what <- match.arg(what)
  ngrp <- lavInspect(lavobj, what = "ngroups")
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
      a <- do.call(compute_a_from_mat,
                   args = c(method, mat[c("lambda", "psi", "theta")],
                            idx = list(idx)))
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

compute_a_from_mat <- function(method = c("regression", "Bartlett"),
                               lambda, theta, psi = NULL, idx = NULL) {
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
  ngrp <- lavInspect(fit, what = "ngroups")
  est_fits <- lavInspect(fit, what = "est")
  if (ngrp == 1) est_fits <- list(est_fits)
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
      c(unlist(compute_ldfs(x, lavobj = fit, method = method)),
        unlist(evfs_lower))
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
  ngrp <- lavInspect(fit, what = "ngroups")
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
