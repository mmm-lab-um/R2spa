# get_fs() for brms location-scale fits (random sigma coefficient), incl.
# random slopes.
#
# get_fs() on a Gaussian brmsfit with a RANDOM SIGMA coefficient (e.g.
# bf(Reaction ~ Days + (Days|p|Subject), sigma ~ (1|p|Subject))) routes to the
# posterior-based (EAP-analog) path get_fs_blocks.brmsfit_ls(): one row per
# level of the single RE grouping factor, one factor per RE coefficient in
# the model's coefficient order (fs_u0, fs_u1, ...), the scores are the
# posterior means of the level's raw r_ draws, psi is the full correlated RE
# covariance (posterior-mean sd_ / cor_ of the joint mu + sigma term), and
# the per-level fsL / fsT come from the per-level posterior covariance of
# those r_ draws fed through the shared regression-form engine
# (fsL = I - Vpost %*% solve(psi), fsT = fsL %*% Vpost). No design matrix is
# built, so scoring_matrix is NULL and the method is EB-only ("ML",
# corrected_fsT / vfsLT and prior_* are rejected). brms is Suggests-only
# (guarded by require_brms()); the whole file skips when brms is not
# installed. No library().
#
# The fixtures are committed .rds (cached) so no MCMC runs at test time:
#   brms_ls_sleep_p1.rds - 3-coef model {Intercept, Days, sigma_Intercept},
#                          4 chains x 1000 draws, brms seed = 1;
#   brms_ls_sleep_p4.rds - 4-coef model {Intercept, Days, sigma_Intercept,
#                          sigma_Days} (sigma random slope), 4 chains x 1000
#                          draws, brms seed = 2.
# Both are fitted on lme4::sleepstudy with p = factor(rep(c("B1","B2","B3"),
# length.out = nlevels(Subject)))[as.integer(Subject)] added to the data.
# NOTE: brms 2.23 re-encodes the RE grouping factor's levels to
# numeric-string ids in fit$data / the r_ draw names, so the level values
# below are always read from the fit's own coef() array, never hard-coded.

skip_if_not_installed("brms")

fix <- function(name) readRDS(test_path(name))

score_cols <- function(fs, p) as.matrix(fs[, paste0("fs_u", seq_len(p) - 1L), drop = FALSE])

# The RE structure of a location-scale brmsfit: the single grouping factor
# (the only array in stats::coef()), its canonical coefficient order (incl.
# sigma_* coefs) and its level order (the get_fs() row order).
ls_structure <- function(fit) {
  co <- stats::coef(fit)
  group <- names(co)[[1L]]
  co3d <- co[[group]]
  list(group = group, cnms = dimnames(co3d)[[3L]], lvls = dimnames(co3d)[[1L]])
}

# The S x p matrix of one level's raw RE draws, columns in the cnms order:
# mu-submodel draws are r_<group>[<level>,<coef>], sigma-submodel draws
# r_<group>__sigma[<level>,<coef>] (the __sigma infix -> canonical
# sigma_<coef>).
ls_re_draws <- function(draws, group, lvl, cnms) {
  idx <- vapply(cnms, function(cnm) {
    if (startsWith(cnm, "sigma_")) {
      sprintf("r_%s__sigma[%s,%s]", group, lvl, sub("^sigma_", "", cnm))
    } else {
      sprintf("r_%s[%s,%s]", group, lvl, cnm)
    }
  }, character(1L))
  as.matrix(as.data.frame(draws)[, idx, drop = FALSE])
}

# Per-level posterior-mean RE deviations: the per-level coefficients
# (stats::coef, dispatched to brms's method) minus the fixed effects
# (brms::fixef), selecting the term's coefficients in the cnms order.
ls_re_means <- function(fit, group, cnms) {
  co3d <- stats::coef(fit)[[group]]
  b <- as.matrix(co3d[, "Estimate", , drop = FALSE][, 1L, ])
  colnames(b) <- cnms
  b - matrix(as.numeric(brms::fixef(fit)[cnms, "Estimate"]),
             nrow = nrow(b), ncol = length(cnms), byrow = TRUE)
}

# The (unscaled) joint mu + sigma RE covariance rebuilt from the posterior
# MEANS of the sd_ / cor_ hyperparameters: brms's "generic" parameterisation
# carries the lower-triangular elements L[i, j] (i > j) of the correlation
# matrix's Cholesky factor in cor_<group>__<co_j>__<co_i>; D = diag(sd) R
# diag(sd) with R = tcrossprod(L).
ls_re_cov_manual <- function(draws, group, cnms) {
  p <- length(cnms)
  s <- vapply(cnms, function(cnm) mean(draws[[paste0("sd_", group, "__", cnm)]]),
              numeric(1L))
  L <- matrix(0, p, p)
  for (i in 2L:p) {
    for (j in seq_len(i - 1L)) {
      L[i, j] <- mean(draws[[paste0("cor_", group, "__", cnms[j], "__", cnms[i])]])
    }
  }
  for (i in seq_len(p)) {
    L[i, i] <- sqrt(max(0, 1 - sum(L[i, seq_len(i - 1L)]^2)))
  }
  S <- diag(p)
  diag(S) <- s
  list(D = S %*% tcrossprod(L) %*% S, R = tcrossprod(L), s = s)
}

# ============================================================================
# 1. S3 dispatch + output contract (3-coef model)
# ============================================================================

test_that("get_fs.brmsfit_ls: S3 dispatch + output contract (3 coefs)", {
  fit <- fix("brms_ls_sleep_p1.rds")
  expect_true(inherits(fit, "brmsfit"))
  expect_equal(fit$family$family, "gaussian")
  st <- ls_structure(fit)
  expect_equal(st$cnms, c("Intercept", "Days", "sigma_Intercept"))
  fs <- get_fs(fit)
  expect_s3_class(fs, "data.frame")
  expect_equal(nrow(fs), 18L)
  expect_true(all(c("fs_u0", "fs_u1", "fs_u2", "fs_u0_se", "fs_u1_se",
                    "fs_u2_se", "u0_by_fs_u0", "ev_fs_u0",
                    "ecov_fs_u1_fs_u0", "ecov_fs_u2_fs_u1") %in% colnames(fs)))
  # no linear scoring map exists for location-scale models
  expect_null(attr(fs, "scoring_matrix"))
  # one row per level of the RE grouping factor, in the coef() level order
  expect_equal(dim(attr(fs, "fsL")), c(3L, 3L, 18L))
  expect_equal(dim(attr(fs, "fsT")), c(3L, 3L, 18L))
  expect_equal(dimnames(attr(fs, "fsT"))[[3L]], st$lvls)
  expect_equal(dim(attr(fs, "psi")), c(3L, 3L))
  # alpha: the named zero vector over the factor names
  expect_equal(attr(fs, "alpha"), setNames(c(0, 0, 0), c("u0", "u1", "u2")))
})

# ============================================================================
# 2. Scores == posterior means of the RE draws == coef() - fixef() deviations
# ============================================================================

test_that("get_fs.brmsfit_ls: scores equal the posterior-mean random effects", {
  fit <- fix("brms_ls_sleep_p1.rds")
  st <- ls_structure(fit)
  draws <- posterior::as_draws_df(fit)
  fs <- get_fs(fit)
  sc <- score_cols(fs, 3L)
  for (j in seq_along(st$lvls)) {
    M <- ls_re_draws(draws, st$group, st$lvls[j], st$cnms)
    expect_equal(sc[j, , drop = FALSE], colMeans(M), tolerance = 1e-8,
                 ignore_attr = TRUE)
  }
  # the same identity the implementation uses: per-level coefficients minus
  # the fixed effects
  expect_equal(sc, ls_re_means(fit, st$group, st$cnms), tolerance = 1e-8,
               ignore_attr = TRUE)
})

# ============================================================================
# 3. psi: the full correlated RE covariance (PSD, off-diagonals from sd_*cor_)
# ============================================================================

test_that("get_fs.brmsfit_ls: psi is the PSD correlated RE covariance", {
  fit <- fix("brms_ls_sleep_p1.rds")
  st <- ls_structure(fit)
  draws <- posterior::as_draws_df(fit)
  psi <- attr(get_fs(fit), "psi")
  expect_equal(dim(psi), c(3L, 3L))
  expect_equal(psi, t(psi), tolerance = 1e-10)
  expect_true(all(eigen(psi, symmetric = TRUE)$values > 0))
  man <- ls_re_cov_manual(draws, st$group, st$cnms)
  expect_equal(psi, unname(man$D), tolerance = 1e-8, ignore_attr = TRUE)
  # each off-diagonal psi[i, j] == sd_i * sd_j * cor_ij (cor from the
  # Cholesky-rebuilt correlation matrix); the off-diagonals are nonzero
  for (i in 1:3) {
    for (j in setdiff(1:3, i)) {
      expect_equal(psi[i, j], unname(man$s[i]) * unname(man$s[j]) * man$R[i, j],
                   tolerance = 1e-8)
      expect_gt(abs(psi[i, j]), 0)
    }
  }
})

# ============================================================================
# 4. Per level: Vpost = cov(r_ draws) is PSD; fsL / fsT regression identities
# ============================================================================

test_that("get_fs.brmsfit_ls: per-level fsL / fsT from the posterior covariance", {
  fit <- fix("brms_ls_sleep_p1.rds")
  st <- ls_structure(fit)
  draws <- posterior::as_draws_df(fit)
  fs <- get_fs(fit)
  psi <- attr(fs, "psi")
  fsL_arr <- attr(fs, "fsL")
  fsT_arr <- attr(fs, "fsT")
  for (j in seq_along(st$lvls)) {
    M <- ls_re_draws(draws, st$group, st$lvls[j], st$cnms)
    Vpost <- cov(M)
    expect_true(all(eigen(Vpost, symmetric = TRUE)$values > 0))
    fsL_exp <- diag(3L) - Vpost %*% solve(psi)
    fsT_exp <- fsL_exp %*% Vpost
    expect_equal(fsL_arr[, , j], fsL_exp, tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(fsT_arr[, , j], fsT_exp, tolerance = 1e-8, ignore_attr = TRUE)
  }
})

# ============================================================================
# 5. The 3-coef model is fully finite (_se / ev_ / ecov_ columns)
# ============================================================================

test_that("get_fs.brmsfit_ls: all _se / ev_ / ecov_ columns are finite (3 coefs)", {
  fit <- fix("brms_ls_sleep_p1.rds")
  fs <- get_fs(fit)
  sev_cols <- colnames(fs)[grepl("_se$|^ev_|^ecov_", colnames(fs))]
  expect_true(length(sev_cols) > 0)
  expect_true(all(is.finite(as.matrix(fs[, sev_cols]))))
})

# ============================================================================
# 6. 4-coef model: sigma random slope (sigma_Days parsing)
# ============================================================================

test_that("get_fs.brmsfit_ls: 4-coef model with a sigma random slope", {
  fit <- fix("brms_ls_sleep_p4.rds")
  st <- ls_structure(fit)
  expect_equal(length(st$cnms), 4L)
  expect_true("sigma_Days" %in% st$cnms)
  # the sigma-submodel slope shows up in the posterior hyperparameters and
  # in the raw RE draws
  draws <- posterior::as_draws_df(fit)
  expect_true(paste0("sd_", st$group, "__sigma_Days") %in% colnames(draws))
  expect_true(any(grepl(paste0("^r_", st$group, "__sigma\\[.*,Days\\]$"),
                        colnames(draws))))
  fs <- get_fs(fit)
  expect_equal(nrow(fs), 18L)
  expect_true(all(paste0("fs_u", 0:3) %in% colnames(fs)))
  expect_equal(dim(attr(fs, "fsL")), c(4L, 4L, 18L))
  expect_equal(dim(attr(fs, "fsT")), c(4L, 4L, 18L))
  psi <- attr(fs, "psi")
  expect_equal(dim(psi), c(4L, 4L))
  expect_equal(psi, t(psi), tolerance = 1e-10)
  expect_true(all(eigen(psi, symmetric = TRUE)$values > 0))
  # NOTE: some _se / ev_ cells may be NA here (a negative regression-form
  # fsT diagonal goes through the engine's sqrt_or_na); the finiteness of
  # every cell is therefore NOT asserted for this model.
})

# ============================================================================
# 7. Guards: ML / corrected_fsT / vfsLT / prior_* rejected
# ============================================================================

test_that("get_fs.brmsfit_ls: guards (ML, corrected_fsT, vfsLT, prior_*)", {
  fit <- fix("brms_ls_sleep_p1.rds")
  expect_error(get_fs(fit, method = "ML"), "location-scale")
  expect_error(get_fs(fit, corrected_fsT = TRUE), "corrected_fsT")
  expect_error(get_fs(fit, vfsLT = TRUE), "vfsLT")
  expect_error(get_fs(fit, prior_mean = c(u0 = 0)), "prior_mean")
  expect_error(get_fs(fit, prior_cov = diag(3)), "prior_cov")
})

# ============================================================================
# 8. Guard: non-Gaussian family rejected (shared binomial guard fixture)
# ============================================================================

test_that("get_fs.brmsfit_ls: Gaussian-only family guard", {
  expect_error(get_fs(fix("brms_guard_binom.rds")), "Gaussian")
})

# ============================================================================
# 9. No regression: a plain Gaussian fit (no random sigma) keeps the existing
#    path (scoring matrix + ML)
# ============================================================================

test_that("get_fs.brmsfit: plain Gaussian fit (no random sigma) is unchanged", {
  fit <- fix("brms_sleep_p1.rds")
  fs <- get_fs(fit)
  expect_false(is.null(attr(fs, "scoring_matrix")))
  expect_equal(length(attr(fs, "scoring_matrix")), 18L)
  expect_no_error(get_fs(fit, method = "EB"))
  expect_no_error(get_fs(fit, method = "ML"))
})

# ============================================================================
# 10. Downstream smoke: tspa() on the V1 location factors
# ============================================================================

test_that("tspa() runs on a brms location-scale get_fs() result", {
  fit <- fix("brms_ls_sleep_p1.rds")
  fs <- get_fs(fit)
  # multi-factor path: the per-level fsL / fsT attributes are pooled
  # (mirrors the merMod smoke test tspa("u1 ~ u0", data = fs_mer))
  fit_mf <- suppressWarnings(tspa("u1 ~ u0", data = fs))
  expect_s4_class(fit_mf, "lavaan")
  # single-factor path: one location factor (fs_u0) with its SE column and a
  # small synthetic dependent variable
  set.seed(2026)
  df1 <- data.frame(y = fs$fs_u0 + rnorm(nrow(fs)), fs_u0 = fs$fs_u0)
  fit_sf <- suppressWarnings(tspa("y ~ u0", data = df1,
                                  se_fs = list(u0 = fs$fs_u0_se[1L])))
  expect_s4_class(fit_sf, "lavaan")
})
