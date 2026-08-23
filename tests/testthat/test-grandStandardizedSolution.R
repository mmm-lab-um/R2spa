library(lavaan)

## Maximum expected error of estimates
err <- .0001

# Single-group, two-factor -----------------------------------------------------

mod1 <- '
   # latent variables
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
   # regressions
     dem60 ~ ind60
'
fit1 <- sem(model = mod1,
            data  = PoliticalDemocracy)
s2_std_beta <- grandStandardizedSolution(fit1)

### lavaan::standardizedSolution()
s2_std_beta_lav <- subset(standardizedSolution(fit1), op == "~")

test_that("Test for single group warning message", {
  expect_message(grandStandardizedSolution(fit1))
})

test_that("Standardized beta in a model with single group, two factors",
          { expect_equal(s2_std_beta$est.std, s2_std_beta_lav$est.std) })
test_that("SE of standardized beta in a model with single group, two-factors",
          { expect_equal(s2_std_beta$se, s2_std_beta_lav$se, tolerance = 1e-7) })

# Single-group, three-factor ---------------------------------------------------

mod2 <- '
    # latent variables
      ind60 =~ x1 + x2 + x3
      dem60 =~ y1 + y2 + y3 + y4
      dem65 =~ y5 + y6 + y7 + y8
    # regressions
      dem60 ~ ind60
      dem65 ~ ind60 + dem60
'
fit2 <- sem(model = mod2,
            data  = PoliticalDemocracy)
s3_std_beta <- grandStandardizedSolution(fit2)

### lavaan::standardizedSolution()
s3_std_beta_lav <- subset(standardizedSolution(fit2), op == "~")

test_that("Standardized beta in a model with single group, three factors", {
  expect_equal(s3_std_beta_lav$est.std, s3_std_beta$est.std)
})
test_that("SE of standardized beta in a model with single group, three factors", {
  expect_equal(s3_std_beta_lav$se, s3_std_beta$se, tolerance = err / 100)
})

# Multigroup, two-factor -------------------------------------------------------

mod3 <- '
  # latent variable definitions
    visual =~ x1 + x2 + x3
    speed =~ x7 + x8 + x9
  # regressions
    visual ~ c(b1, b1) * speed
'
fit3 <- sem(mod3, data = HolzingerSwineford1939,
            group = "school",
            group.equal = c("loadings", "intercepts"))

m2_std_beta <- grandStandardizedSolution(fit3)

test_that("grand std est in the range of std estimates", code = {
  gp_std_est <- subset(standardizedSolution(fit3), subset = op == "~")$est
  gr_std_est <- m2_std_beta$est
  expect_gte(mean(gr_std_est), min(gp_std_est))
  expect_lte(mean(gr_std_est), max(gp_std_est))
})

test_that("grand std SE similar to average std SE", code = {
  gp_std_se <- subset(standardizedSolution(fit3), subset = op == "~")$se
  gr_std_se <- m2_std_beta$se
  expect_equal(mean(gp_std_se), mean(gr_std_se), tolerance = 0.1)
})

## Hand calculation
### std.est
model_list <- lavTech(fit3, what = "est")
ns <- lavInspect(fit3, what = "nobs")
beta_list <- model_list[which(names(model_list) == "beta")]
psi_list <- model_list[which(names(model_list) == "psi")]
alpha_list <- model_list[which(names(model_list) == "alpha")]
v_eta <- veta_grand(ns,
                    beta_list,
                    psi_list = psi_list,
                    alpha_list = alpha_list)
s_eta <- sqrt(diag(v_eta))
inv_s_eta <- 1 / s_eta
std_betas <- lapply(beta_list, function(x) {
  diag(inv_s_eta) %*% x %*% diag(s_eta)
})
m2_std_betas_h <- unlist(std_betas)[unlist(std_betas) != 0]
### se
free_list <- lavTech(fit3, what = "free")
acov_par <- vcov(fit3)
free_beta_psi_alpha <- free_list[which(names(model_list) %in%
                                         c("beta", "psi", "alpha"))]
est <- .combine_est(model_list[which(names(model_list) %in%
                                       c("beta", "psi", "alpha"))],
                    free = free_beta_psi_alpha)
jac <- lav_func_jacobian_complex(function(x)
  unlist(grand_std_beta_est(model_list, ns = ns, free_list = free_list, est = x)),
  x = est)
pos_beta_psi_alpha <- .combine_est(free_beta_psi_alpha,
                                   free = free_beta_psi_alpha)
acov_beta_psi_alpha <- acov_par[pos_beta_psi_alpha, pos_beta_psi_alpha]
acov <- jac %*% acov_beta_psi_alpha %*% t(jac)
m2_std_se_h <- sqrt(acov[3, c(3, 7)])

test_that("veta_grand() gives expected results", code = {
  lvcov_fit3 <- lavInspect(fit3, what = "cov.lv")
  lvmean_fit3 <- lavInspect(fit3, what = "mean.lv")
  ns_fit3 <- lavInspect(fit3, what = "nobs")
  v_eta_hand <- Reduce(
    `+`,
    mapply(function(v, m, n) n *
             (v + tcrossprod(m - do.call(cbind, lvmean_fit3) %*%
                               ns_fit3 / sum(ns_fit3))),
           v = lvcov_fit3, m = lvmean_fit3, n = ns_fit3, SIMPLIFY = FALSE)
  ) / sum(ns_fit3)
  expect_true(all(v_eta > Reduce(`+`, x = lvcov_fit3) / 2))
  expect_equal(v_eta, v_eta_hand, ignore_attr = TRUE)
})

test_that("Standardized beta in a model with multiple groups, two factors", {
  expect_equal(m2_std_betas_h, m2_std_beta$est.std,
                ignore_attr = TRUE)
})
test_that("SE of standardized beta in a model with multiple groups, two factors", {
  expect_equal(m2_std_se_h, m2_std_beta$se)
})


# Multigroup, three-factor -----------------------------------------------------

mod4 <- '
  # latent variable definitions
    visual =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed =~ x7 + x8 + x9

  # regressions
    visual ~ c(b1, b1) * textual + c(b2, b2) * speed
'
fit4 <- sem(mod4, data = HolzingerSwineford1939,
            group = "school",
            group.equal = c("loadings", "intercepts"))
m3_std_beta <- grandStandardizedSolution(fit4)

## Hand calculation
### std.est
model_list <- lavTech(fit4, what = "est")
ns <- lavInspect(fit4, what = "nobs")
beta_list <- model_list[which(names(model_list) == "beta")]
psi_list <- model_list[which(names(model_list) == "psi")]
alpha_list <- model_list[which(names(model_list) == "alpha")]
v_eta <- veta_grand(ns,
                    beta_list,
                    psi_list = psi_list,
                    alpha_list = alpha_list)
s_eta <- sqrt(diag(v_eta))
inv_s_eta <- 1 / s_eta
std_betas <- lapply(beta_list, function(x) {
  diag(inv_s_eta) %*% x %*% diag(s_eta)
})
m3_std_betas_h <- unlist(std_betas)[unlist(std_betas) != 0]
### se
free_list <- lavTech(fit4, what = "free")
acov_par <- vcov(fit4)
free_beta_psi_alpha <- free_list[which(names(model_list) %in%
                                         c("beta", "psi", "alpha"))]
est <- .combine_est(model_list[which(names(model_list) %in%
                                       c("beta", "psi", "alpha"))],
                    free = free_beta_psi_alpha)
jac <- lav_func_jacobian_complex(function(x)
  unlist(grand_std_beta_est(model_list, ns = ns, free_list = free_list, est = x)),
  x = est)
pos_beta_psi_alpha <- .combine_est(free_beta_psi_alpha,
                                   free = free_beta_psi_alpha)
acov_beta_psi_alpha <- acov_par[pos_beta_psi_alpha, pos_beta_psi_alpha]
std_se <- jac %*% acov_beta_psi_alpha %*% t(jac)
m3_std_se_h <- sqrt(diag(std_se[c(4, 7, 13, 16), c(4, 7, 13, 16)]))

test_that("Standardized beta in a model with multiple groups, three factors", {
  expect_equal(m3_std_betas_h, m3_std_beta$est.std, ignore_attr = TRUE)
})

test_that("SE of standardized beta in a model with multiple groups, three factors", {
  expect_equal(m3_std_se_h, m3_std_beta$se, ignore_attr = TRUE)
})

# Quarantined with R/grandStandardizedSolution.R (see _PLAN_QUARANTINE.md).
# Extracted from tests/testthat/test-lavaan_compat.R: the wrapper A/B test for
# grandStandardizedSolution (lines 107-121 as of 2026-08-17) plus its setup
# (lines 9-18: canon_mod/canon_fit, mg_mod/mg_fit).

canon_mod <- "ind60 =~ x1 + x2 + x3
              dem60 =~ y1 + y2 + y3 + y4
              dem60 ~ ind60"
canon_fit <- sem(canon_mod, data = PoliticalDemocracy, std.lv = TRUE)

mg_mod <- "visual =~ x1 + x2 + x3
           speed =~ x7 + x8 + x9
           visual ~ speed"
mg_fit <- sem(mg_mod, data = HolzingerSwineford1939, group = "school",
              std.lv = TRUE)

test_that("grand_standardized_solution output is unchanged by the wrapper", {
  # Single-group: must match lavaan::standardizedSolution() as before.
  got <- suppressMessages(grandStandardizedSolution(canon_fit))
  lav <- subset(standardizedSolution(canon_fit), op == "~")
  expect_equal(got$est.std, lav$est.std)
  expect_equal(got$se, lav$se, tolerance = 1e-7)
  # The returned frame keeps its historical columns (incl. `exo`).
  expect_true(all(c("lhs", "op", "rhs", "exo", "group", "block", "label",
                    "est.std", "se") %in% names(got)))
  # Multigroup path runs through the wrappers and stays finite; the MG
  # hand-calculation A/B lives in test-grandStandardizedSolution.R.
  mg <- grandStandardizedSolution(mg_fit)
  expect_true(all(is.finite(mg$est.std)))
  expect_true(all(is.finite(mg$se)))
})

## MG 2S-PA: corrected grand-standardized SEs (SE-only correction) ------
## A distinct model from fit3 above (`visual ~ speed` free in both
## `school` groups, no group.equal constraints): fit by two-stage path
## analysis on factor scores so a corrected stage-2 covariance exists to
## thread through grandStandardizedSolution() -- the existing MG fixtures
## in this file are plain sem() fits over raw data, which cannot be
## SE-corrected. The corrected fit is built once at file scope for the
## same cost reason as in test-tspa_corrected_se.R (~28 stage-2 refits).
mod_tspa_mg <- "visual =~ x1 + x2 + x3\nspeed =~ x7 + x8 + x9"
fs_mg_gs <- get_fs(HolzingerSwineford1939, model = mod_tspa_mg, std.lv = TRUE,
                   group = "school", vfsLT = TRUE, format = "list")
tspa_mg_gs_plain <- tspa("visual ~ speed", data = do.call(rbind, fs_mg_gs),
                         fsT = attr(fs_mg_gs, "fsT"), fsL = attr(fs_mg_gs, "fsL"),
                         group = "school")
tspa_mg_gs_corr <- tspa("visual ~ speed", data = do.call(rbind, fs_mg_gs),
                        fsT = attr(fs_mg_gs, "fsT"), fsL = attr(fs_mg_gs, "fsL"),
                        vfsLT = attr(fs_mg_gs, "vfsLT"), corrected_se = TRUE,
                        group = "school")
gr_gs_plain <- grandStandardizedSolution(tspa_mg_gs_plain)
gr_gs_corr  <- grandStandardizedSolution(tspa_mg_gs_corr)
gr_gs_inj   <- grandStandardizedSolution(tspa_mg_gs_plain,
                                         acov_par = vcov(tspa_mg_gs_corr))

test_that("grand standardization threads the corrected covariance (MG, 2S-PA)", {
  keep <- gr_gs_corr$op == "~"
  # (a) SE-only: corrected grand-std SEs >= uncorrected
  expect_true(all(gr_gs_corr$se[keep] >= gr_gs_plain$se[keep] - 1e-8))
  # (b) A/B invariant: injecting the corrected covariance into the plain
  #     fit equals calling grandStandardizedSolution on the corrected fit
  expect_equal(gr_gs_inj$se[keep], gr_gs_corr$se[keep], tolerance = 1e-8)
  expect_equal(gr_gs_inj$est.std[keep], gr_gs_corr$est.std[keep],
               tolerance = 1e-8)
  # (c) point estimates unchanged by the correction
  expect_equal(gr_gs_corr$est.std[keep], gr_gs_plain$est.std[keep],
               tolerance = 1e-8)
})
