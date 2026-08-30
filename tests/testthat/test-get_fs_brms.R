# get_fs() for brms brmsfit fits (Gaussian mixed models).
#
# get_fs() on a brmsfit returns one row per level of the (single) random-effects
# term, mirroring get_fs.merMod(): the score columns (fs_u0, ..., fs_u<p-1>) are
# the posterior-mean random effects; fsL / fsT / scoring_matrix come from the
# same EB / ML formulas get_fs_blocks.merMod() uses, generalised to a p x p
# random-effects covariance D reconstructed from the posterior of the term's
# sd_ / cor_ hyperparameters. brms is Suggests-only (guarded by require_brms());
# the fixed/random design and cluster structure are built directly from the
# model formula + data (reformulas, bit-identical to lme4::getME), so no
# throwaway fit is required.
#
# These tests cover: S3 dispatch + the output contract, the EB score == ranef
# identity, the EB scoring identity (S_j %*% (y_j - X_j beta) ~ b_j), the D
# reconstruction (PSD + exact brms_re_cov against synthetic draws, incl. p = 3
# for the Cholesky extension), the ML method (fsL = I, per-cluster OLS score),
# the brms-vs-lme4 cross-check, and the guards (Gaussian-only, single random
# term, prior_* rejection).
#
# The fixtures are committed .rds (cached) so no MCMC runs at test time. brms is
# Suggests-only; the whole file skips when brms is not installed. No library().

skip_if_not_installed("brms")

fix <- function(name) readRDS(test_path(name))

# Number of coefficients in the single random-effects term.
n_coefs <- function(fit) {
  rt <- reformulas::mkReTrms(reformulas::findbars(fit$formula$formula),
                             fit$data, calc.lambdat = FALSE)
  length(rt$cnms[[1L]])
}
# Posterior-mean random effects (n_clus x p) in the term's canonical level order
# (per-level coefficients minus fixed effects, the same identity the method uses).
re_mat <- function(fit) {
  rt <- reformulas::mkReTrms(reformulas::findbars(fit$formula$formula),
                             fit$data, calc.lambdat = FALSE)
  group <- names(rt$flist)[[1L]]
  co3d <- stats::coef(fit)[[group]]
  all_terms <- dimnames(co3d)[[3L]]
  m <- as.matrix(co3d[, "Estimate", , drop = FALSE][, 1L, ])
  colnames(m) <- all_terms
  m <- m - matrix(as.numeric(brms::fixef(fit)[all_terms, "Estimate"]),
                  nrow = nrow(m), ncol = length(all_terms), byrow = TRUE)
  m <- m[, match(gsub("[()]", "", rt$cnms[[1L]]), all_terms), drop = FALSE]
  m[match(levels(rt$flist[[1L]]), dimnames(co3d)[[1L]]), , drop = FALSE]
}
# Dense random design (n x (p * n_clus)), level-major, from reformulas (== the
# lme4 Z fold invariant): row i holds i's own cluster's p coefficients.
z_dense <- function(fit) {
  rt <- reformulas::mkReTrms(reformulas::findbars(fit$formula$formula),
                             fit$data, calc.lambdat = FALSE)
  Zt <- rt$Zt
  Z <- matrix(0, nrow = ncol(Zt), ncol = nrow(Zt))
  Z[cbind(rep(seq_len(ncol(Zt)), diff(Zt@p)), Zt@i + 1L)] <- Zt@x
  Z
}
score_cols <- function(fs, p) as.matrix(fs[, paste0("fs_u", seq_len(p) - 1L), drop = FALSE])

core_fixtures <- c("brms_sleep_p1.rds", "brms_sleep_p2.rds")

# ============================================================================
# 1. S3 dispatch + output contract
# ============================================================================

test_that("get_fs.brmsfit: S3 dispatch + output contract (p=1)", {
  fit <- fix("brms_sleep_p1.rds")
  expect_true(inherits(fit, "brmsfit"))
  fs <- get_fs(fit)
  expect_s3_class(fs, "data.frame")
  expect_equal(nrow(fs), 18)
  expect_true(all(c("fs_u0", "fs_u0_se", "u0_by_fs_u0", "ev_fs_u0") %in% colnames(fs)))
  expect_equal(dim(attr(fs, "fsL")), c(1L, 1L, 18L))
  expect_equal(dim(attr(fs, "fsT")), c(1L, 1L, 18L))
  expect_equal(dim(attr(fs, "psi")), c(1L, 1L))
  expect_equal(length(attr(fs, "alpha")), 1L)
  expect_equal(length(attr(fs, "scoring_matrix")), 18L)
  # one row per level of the RE term, named by level (level order)
  expect_equal(names(attr(fs, "scoring_matrix")),
               levels(as.factor(fit$data$Subject)))
})

test_that("get_fs.brmsfit: output contract (p=2, random slopes)", {
  fit <- fix("brms_sleep_p2.rds")
  expect_equal(n_coefs(fit), 2L)
  fs <- get_fs(fit)
  expect_equal(nrow(fs), 18)
  expect_true(all(c("fs_u0", "fs_u1", "ev_fs_u0", "ev_fs_u1") %in% colnames(fs)))
  expect_equal(dim(attr(fs, "fsL")), c(2L, 2L, 18L))
  expect_equal(dim(attr(fs, "fsT")), c(2L, 2L, 18L))
  expect_equal(dim(attr(fs, "psi")), c(2L, 2L))
  expect_equal(dim(attr(fs, "scoring_matrix")[[1L]]), c(2L, 10L))
})

# ============================================================================
# 2. EB scores == posterior-mean random effects
# ============================================================================

test_that("get_fs.brmsfit: EB scores equal the posterior-mean random effects (p=1, p=2)", {
  for (tag in core_fixtures) {
    fit <- fix(tag)
    p <- n_coefs(fit)
    expect_equal(score_cols(get_fs(fit), p), re_mat(fit), tolerance = 1e-8,
                 ignore_attr = TRUE)
  }
})

# ============================================================================
# 3. EB scoring identity (loose: posterior hyperparameter + MCMC noise)
# ============================================================================

test_that("get_fs.brmsfit: EB scoring identity S_j %*% (y_j - X_j beta) ~ b_j", {
  for (tag in core_fixtures) {
    fit <- fix(tag)
    f <- fit$formula$formula
    data <- fit$data
    rt <- reformulas::mkReTrms(reformulas::findbars(f), data, calc.lambdat = FALSE)
    f1 <- rt$flist[[1L]]
    p <- length(rt$cnms[[1L]])
    Z <- z_dense(fit)
    ci <- split(seq_len(nrow(data)), f1)
    fx <- reformulas::nobars(f)
    X <- as.matrix(stats::model.matrix(fx, data))
    y <- stats::model.response(stats::model.frame(fx, data))
    beta <- as.numeric(brms::fixef(fit)[, "Estimate"])
    b <- re_mat(fit)
    sm <- attr(get_fs(fit), "scoring_matrix")
    maxdev <- 0
    for (j in seq_len(nlevels(f1))) {
      idx <- ci[[j]]
      zj <- Z[idx, (j - 1L) * p + seq_len(p), drop = FALSE]
      rj <- y[idx] - as.numeric(X[idx, , drop = FALSE] %*% beta)
      maxdev <- max(maxdev, max(abs(as.numeric(sm[[j]] %*% rj) - as.numeric(b[j, , drop = FALSE]))))
    }
    expect_lt(maxdev, 10)
  }
})

# ============================================================================
# 4. D reconstruction: PSD + exact p=1 + symmetric
# ============================================================================

test_that("get_fs.brmsfit: D is a symmetric PSD random-effects covariance", {
  for (tag in core_fixtures) {
    fit <- fix(tag)
    p <- n_coefs(fit)
    D <- attr(get_fs(fit), "psi")
    expect_equal(dim(D), c(p, p))
    expect_equal(D, t(D), tolerance = 1e-10)
    expect_true(all(eigen(D, symmetric = TRUE)$values > 0))
  }
  # p=1: D[1, 1] == (posterior-mean sd)^2 exactly
  fit <- fix("brms_sleep_p1.rds")
  draws <- posterior::as_draws_df(fit)
  sd1 <- mean(draws[["sd_Subject__Intercept"]])
  expect_equal(attr(get_fs(fit), "psi")[1L, 1L], sd1^2, tolerance = 1e-8)
})

# ============================================================================
# 5. brms_re_cov: exact D from sd_/cor_ draws (p = 1, 2, 3)
# ============================================================================

test_that("brms_re_cov: reconstructs D from the sd_/cor_ posterior (p=1,2,3)", {
  # p = 1: D = sd^2
  d1 <- data.frame(sd_G__a = 2)
  expect_equal(R2spa:::brms_re_cov(d1, "G", "(a)"), matrix(4, 1, 1),
               tolerance = 1e-10, ignore_attr = TRUE)

  # p = 2: cor_G__a__b is the Cholesky off-diagonal L[2, 1]
  s2 <- c(2, 3)
  L21 <- 0.5
  d2 <- data.frame(sd_G__a = s2[1], sd_G__b = s2[2], cor_G__a__b = L21)
  L <- diag(2); L[2, 1] <- L21; L[2, 2] <- sqrt(1 - L21^2)
  S <- diag(2); diag(S) <- s2
  expect_equal(R2spa:::brms_re_cov(d2, "G", c("(a)", "(b)")),
               S %*% (L %*% t(L)) %*% S, tolerance = 1e-10, ignore_attr = TRUE)

  # p = 3: three Cholesky off-diagonals (L[2,1], L[3,1], L[3,2])
  s3 <- c(1, 2, 3)
  cL <- c(L21 = 0.3, L31 = -0.2, L32 = 0.4)
  d3 <- data.frame(sd_G__a = s3[1], sd_G__b = s3[2], sd_G__c = s3[3],
                   cor_G__a__b = cL["L21"], cor_G__a__c = cL["L31"], cor_G__b__c = cL["L32"])
  L3 <- matrix(0, 3, 3)
  L3[2, 1] <- cL["L21"]; L3[3, 1] <- cL["L31"]; L3[3, 2] <- cL["L32"]
  L3[1, 1] <- 1; L3[2, 2] <- sqrt(1 - cL["L21"]^2); L3[3, 3] <- sqrt(1 - cL["L31"]^2 - cL["L32"]^2)
  S3 <- diag(3); diag(S3) <- s3
  expect_equal(R2spa:::brms_re_cov(d3, "G", c("(a)", "(b)", "(c)")),
               S3 %*% (L3 %*% t(L3)) %*% S3, tolerance = 1e-10, ignore_attr = TRUE)
})

# ============================================================================
# 6. ML method: fsL = I, per-cluster OLS score
# ============================================================================

test_that("get_fs.brmsfit: ML method (fsL = I, per-cluster OLS score)", {
  for (tag in core_fixtures) {
    fit <- fix(tag)
    p <- n_coefs(fit)
    fs <- get_fs(fit, method = "ML")
    for (j in seq_len(nrow(fs))) {
      expect_equal(as.numeric(attr(fs, "fsL")[, , j]), as.numeric(diag(p)),
                   tolerance = 1e-10)
    }
    f <- fit$formula$formula
    data <- fit$data
    rt <- reformulas::mkReTrms(reformulas::findbars(f), data, calc.lambdat = FALSE)
    f1 <- rt$flist[[1L]]
    Z <- z_dense(fit)
    ci <- split(seq_len(nrow(data)), f1)
    fx <- reformulas::nobars(f)
    X <- as.matrix(stats::model.matrix(fx, data))
    y <- stats::model.response(stats::model.frame(fx, data))
    beta <- as.numeric(brms::fixef(fit)[, "Estimate"])
    sc <- score_cols(fs, p)
    for (j in seq_len(nlevels(f1))) {
      idx <- ci[[j]]
      zj <- Z[idx, (j - 1L) * p + seq_len(p), drop = FALSE]
      rj <- y[idx] - as.numeric(X[idx, , drop = FALSE] %*% beta)
      Gz <- MASS::ginv(crossprod(zj))
      expect_equal(as.numeric(sc[j, , drop = FALSE]),
                   as.numeric(t(Gz %*% crossprod(zj, rj))), tolerance = 1e-6)
    }
  }
})

# ============================================================================
# 7. brms-vs-lme4 cross-check (loose)
# ============================================================================

test_that("get_fs.brmsfit: EB scores are close to the lme4 reference (loose)", {
  for (tag in core_fixtures) {
    fit <- fix(tag)
    p <- n_coefs(fit)
    f <- fit$formula$formula
    maxdev <- max(abs(score_cols(get_fs(fit), p) -
                       score_cols(get_fs(lme4::lmer(f, fit$data)), p)))
    expect_lt(maxdev, 3)
  }
})

# ============================================================================
# 8. Guards: Gaussian-only, single random term, prior_* rejection
# ============================================================================

test_that("get_fs.brmsfit: guards (Gaussian-only, single term, prior_*)", {
  expect_error(get_fs(fix("brms_guard_binom.rds")), "Gaussian")
  expect_error(get_fs(fix("brms_guard_multibar.rds")), "single random-effects term")
  expect_error(get_fs(fix("brms_sleep_p1.rds"), prior_mean = 0), "prior_mean")
  expect_error(get_fs(fix("brms_sleep_p1.rds"), prior_cov = diag(1)), "prior_cov")
})
