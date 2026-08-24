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

## Row assignment: partable row order, not column-major beta order --------
##
## The standardized estimates live in the per-group beta matrices in
## column-major order; partable rows follow model-statement order. They
## coincide only by accident (single-predictor / single-endogenous models
## -- all fixtures above), so a model with >= 2 endogenous variables and
## interleaved predictors is needed to pin the row <-> value mapping.

mod5 <- '
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # regressions (statement order deliberately differs from the beta
  # matrix column-major free order: dem60 ~ ind60 sorts first in
  # column-major but last here)
    dem65 ~ dem60 + ind60
    dem60 ~ ind60
'
fit5 <- sem(mod5, data = PoliticalDemocracy)
s5_std_beta <- suppressMessages(grandStandardizedSolution(fit5))
s5_std_beta_lav <- subset(standardizedSolution(fit5), op == "~")

test_that("SG row assignment follows partable order (column-major swap guard)", {
  expect_identical(
    paste(s5_std_beta$lhs, s5_std_beta$rhs),
    paste(s5_std_beta_lav$lhs, s5_std_beta_lav$rhs)
  )
  expect_equal(s5_std_beta$est.std, s5_std_beta_lav$est.std)
  expect_equal(s5_std_beta$se, s5_std_beta_lav$se, tolerance = 1e-7)
})

## MG 2S-PA, two endogenous scores per group (vignette shape): the
## independent per-row check of the same mapping on a stage-2 fit. The
## hand calculation uses only partable estimates + nobs (no veta_grand):
## build the per-group implied moments from the structural equations and
## pool them into the grand covariance, then standardize each slope by
## the matching grand SDs.
mod6 <- "visual =~ x1 + x2 + x3
         textual =~ x4 + x5 + x6
         speed =~ x7 + x8 + x9"
cfa6 <- cfa(mod6, data = HolzingerSwineford1939, std.lv = TRUE,
            group = "school",
            group.equal = c("loadings", "intercepts"),
            group.partial = c("visual=~x2", "x7~1"))
fs6 <- get_fs_lavaan(cfa6)
fit6 <- tspa("textual ~ visual + speed
              visual ~ speed", data = fs6, group = "school",
             fsL = attr(fs6, "fsL"), fsT = attr(fs6, "fsT"))
g6_std_beta <- suppressMessages(grandStandardizedSolution(fit6))

pt6 <- lavaan::partable(fit6)
ns6 <- lavInspect(fit6, what = "nobs")
e6 <- function(lhs, op, rhs = NULL, g) {
  i <- pt6$lhs == lhs & pt6$op == op & pt6$group == g
  if (!is.null(rhs)) i <- i & pt6$rhs == rhs
  as.numeric(pt6$est[i])
}
mom6 <- function(g) {
  # (t = textual, v = visual, s = speed); s is exogenous, t/v endogenous
  btv <- e6("textual", "~", "visual", g)
  cts <- e6("textual", "~", "speed", g)
  cvs <- e6("visual",  "~", "speed", g)
  pv <- e6("visual",  "~~", NULL, g)
  pt_ <- e6("textual", "~~", NULL, g)
  ps <- e6("speed",   "~~", NULL, g)
  at <- e6("textual", "~1", NULL, g)
  av <- e6("visual",  "~1", NULL, g)
  as_ <- e6("speed",   "~1", NULL, g)
  var_s <- ps
  var_v <- pv + cvs^2 * var_s
  var_t <- pt_ + cts^2 * var_s + btv^2 * var_v +
    2 * btv * cts * cvs * var_s
  cov_tv <- btv * var_v + cts * cvs * var_s
  cov_ts <- (cts + btv * cvs) * var_s
  cov_vs <- cvs * var_s
  cov <- rbind(
    c(var_t, cov_tv, cov_ts),
    c(cov_tv, var_v, cov_vs),
    c(cov_ts, cov_vs, var_s)
  )
  mu <- c(at + cts * as_ + btv * (av + cvs * as_),
          av + cvs * as_,
          as_)
  list(cov = cov, mu = mu)
}
m6 <- lapply(seq_along(ns6), mom6)
mu_grand <- Reduce(`+`, mapply(function(m, n) m * n,
                                lapply(m6, `[[`, "mu"), ns6)) / sum(ns6)
grand_cov6 <- Reduce(
  `+`,
  mapply(function(v, m, n) n * (v + tcrossprod(m - mu_grand)),
         v = lapply(m6, `[[`, "cov"), m = lapply(m6, `[[`, "mu"),
         n = ns6, SIMPLIFY = FALSE)
) / sum(ns6)
sd_grand6 <- sqrt(diag(grand_cov6))
var_ord6 <- c("textual", "visual", "speed")
std6_hand <- vapply(seq_len(nrow(g6_std_beta)), function(i) {
  raw <- e6(g6_std_beta$lhs[i], "~", g6_std_beta$rhs[i],
            g6_std_beta$group[i])
  # standardized coefficient = raw * sd(predictor) / sd(outcome)
  raw * sd_grand6[match(g6_std_beta$rhs[i], var_ord6)] /
    sd_grand6[match(g6_std_beta$lhs[i], var_ord6)]
}, numeric(1))

test_that("MG 2S-PA est.std rows match the independent grand-std hand calc", {
  # row identities first (the table must stay in partable order)
  expect_identical(g6_std_beta[["lhs"]], pt6$lhs[pt6$op == "~"])
  expect_identical(g6_std_beta[["rhs"]], pt6$rhs[pt6$op == "~"])
  expect_identical(g6_std_beta[["group"]], pt6$group[pt6$op == "~"])
  expect_equal(g6_std_beta$est.std, std6_hand, tolerance = 1e-10)
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

## Fixed structural slopes ----------------------------------------------------
## A user-fixed slope (`~ k*var`, k a literal) has no free position, so the
## global-free-position anchor (free cell <-> matrix bijection) is undefined
## for it. grand_standardized_solution() now anchors fixed slopes by (lhs, rhs)
## variable identity in the beta matrix and reports the user value rescaled by
## the (grand) SD ratio plus the first-order delta SE, matching
## lavaan::standardizedSolution() (which also reports a delta SE for fixed
## slopes). All-free slopes still use the free-position anchor.

modF <- '
   ind60 =~ x1 + x2 + x3
   dem60 =~ y1 + y2 + y3 + y4
   dem60 ~ 1.5*ind60
'
fitF <- sem(modF, data = PoliticalDemocracy)
sF <- suppressMessages(grandStandardizedSolution(fitF))
sF_lav <- subset(standardizedSolution(fitF), op == "~")

test_that("SG: a fixed structural slope is reported (est.std + delta se)", {
  expect_identical(paste(sF$lhs, sF$rhs), paste(sF_lav$lhs, sF_lav$rhs))
  expect_equal(sF$est.std, sF_lav$est.std)
  expect_equal(sF$se, sF_lav$se, tolerance = 1e-7)
})

# Fixed + free slopes in one model; also pins the column-major row order with a
# fixed cell interleaved among free ones. dem65 ~ 2*dem60 is the fixed slope;
# dem60 ~ ind60 and dem65 ~ ind60 are free.
modFM <- '
   ind60 =~ x1 + x2 + x3
   dem60 =~ y1 + y2 + y3 + y4
   dem65 =~ y5 + y6 + y7 + y8
   dem60 ~ ind60
   dem65 ~ 2*dem60 + ind60
'
fitFM <- sem(modFM, data = PoliticalDemocracy)
sFM <- suppressMessages(grandStandardizedSolution(fitFM))
sFM_lav <- standardizedSolution(fitFM)
sFM_lav <- sFM_lav[sFM_lav$op == "~", ]
ordF <- order(paste(sFM$lhs, sFM$rhs))
ordL <- order(paste(sFM_lav$lhs, sFM_lav$rhs))
test_that("SG: mixed free + fixed slopes all match standardizedSolution", {
  expect_identical(paste(sFM$lhs[ordF], sFM$rhs[ordF]),
                   paste(sFM_lav$lhs[ordL], sFM_lav$rhs[ordL]))
  expect_equal(sFM$est.std[ordF], sFM_lav$est.std[ordL])
  expect_equal(sFM$se[ordF], sFM_lav$se[ordL], tolerance = 1e-6)
})
test_that("SG: a free structural slope is unchanged (regression guard)", {
  free_r <- sFM$lhs == "dem60" & sFM$rhs == "ind60"
  expect_true(any(free_r))
  lav_free <- sFM_lav[sFM_lav$lhs == "dem60" & sFM_lav$rhs == "ind60", ]
  expect_equal(sFM$est.std[free_r], lav_free$est.std)
  expect_equal(sFM$se[free_r], lav_free$se, tolerance = 1e-7)
})

# MG: a fixed slope (equal raw value in both groups) reports the SAME grand
# est.std and grand SE in every group (grand SD is pooled), and the point
# estimate matches an independent hand calc from the pooled grand latent SDs.
modMF <- '
   visual =~ x1 + x2 + x3
   speed  =~ x7 + x8 + x9
   visual ~ 0.5*speed
'
fitMF <- sem(modMF, data = HolzingerSwineford1939, group = "school")
sMF <- suppressMessages(grandStandardizedSolution(fitMF))
kMF <- sMF$op == "~"
lvcov <- lavInspect(fitMF, what = "cov.lv")
lvmean <- lavInspect(fitMF, what = "mean.lv")
nsg <- lavInspect(fitMF, what = "nobs")
mu_g <- Reduce("+", mapply(function(m, n) n * m, m = lvmean, n = nsg,
                           SIMPLIFY = FALSE)) / sum(nsg)
grand_cov <- Reduce("+", mapply(function(v, m, n) n * (v + tcrossprod(m - mu_g)),
                                v = lvcov, m = lvmean, n = nsg,
                                SIMPLIFY = FALSE)) / sum(nsg)
sd_g <- sqrt(diag(grand_cov))
vrn <- rownames(lvcov[[1]])
hand_stdMF <- 0.5 * sd_g[match("speed", vrn)] / sd_g[match("visual", vrn)]

test_that("MG: a fixed slope is identical across groups and matches the hand calc", {
  e <- sMF$est.std[kMF]
  s_ <- sMF$se[kMF]
  expect_length(e, 2)
  expect_equal(e[1], e[2], tolerance = 1e-12)
  expect_equal(s_[1], s_[2], tolerance = 1e-12)
  expect_equal(unname(e[1]), unname(hand_stdMF), tolerance = 1e-8)
})

# 2S-PA corrected SG: the fixed slope's grand SE is threaded from the corrected
# stage-2 covariance; the SG grand SE equals a plain per-group standardized
# solution, so it must match standardizedSolution() on the corrected fit.
mod_tspa_f <- "ind60 =~ x1 + x2 + x3
               dem60 =~ y1 + y2 + y3 + y4"
fs_f <- get_fs(PoliticalDemocracy, model = mod_tspa_f, std.lv = TRUE,
               vfsLT = TRUE)
fitF_corr <- tspa("dem60 ~ 1.5*ind60", data = fs_f,
                  fsT = attr(fs_f, "fsT"), fsL = attr(fs_f, "fsL"),
                  vfsLT = attr(fs_f, "vfsLT"), corrected_se = TRUE)
sF_corr <- suppressMessages(grandStandardizedSolution(fitF_corr))
sF_corr_lav <- standardizedSolution(fitF_corr)
sF_corr_lav <- sF_corr_lav[sF_corr_lav$op == "~", ]
test_that("2S-PA corrected SG: a fixed slope SE matches standardizedSolution", {
  expect_equal(sF_corr$est.std, sF_corr_lav$est.std, tolerance = 1e-10)
  expect_equal(sF_corr$se, sF_corr_lav$se, tolerance = 1e-6)
})

# Defensive guard: a structural row that cannot be resolved to a beta-matrix
# cell (its lhs or rhs is not among the beta variable names) is rejected with
# an error rather than silently producing NA rows. In lavaan 0.7-x every '~'
# row lands in the beta block (observed-on-observed included; latent/observed
# means are op '~1', not '~'), so the guard is a drift safeguard against a
# future layout where a '~' row is not a beta cell. Exercised by mocking the
# dimname anchor to drop the predictor.
test_that("defensive guard: an unresolvable structural slope errors, not NAs", {
  local_mocked_bindings(
    tsp_beta_names = function(fit) list(
      list(rnm = c("ind60", "dem60"), clm = "dem60")
    ),
    .package = "R2spa"
  )
  fitGuard <- sem("ind60 =~ x1 + x2 + x3
                     dem60 =~ y1 + y2 + y3 + y4
                     dem60 ~ 1.5*ind60",
                  data = PoliticalDemocracy)
  expect_error(grandStandardizedSolution(fitGuard), "beta-matrix slope")
})

# tsp_beta_names(): the fixed-slope dimname anchor (canary for 0.7-x format).
test_that("tsp_beta_names: SG + MG beta dimnames resolve the fixed cell", {
  bn <- tsp_beta_names(fitF)
  expect_length(bn, 1)
  expect_setequal(bn[[1]]$rnm, c("ind60", "dem60"))
  expect_setequal(bn[[1]]$clm, c("ind60", "dem60"))
  expect_false(is.na(match("dem60", bn[[1]]$rnm)))
  expect_false(is.na(match("ind60", bn[[1]]$clm)))
  bnM <- tsp_beta_names(fitMF)
  expect_length(bnM, 2)
  expect_setequal(bnM[[1]]$rnm, c("visual", "speed"))
  expect_setequal(bnM[[1]]$clm, c("visual", "speed"))
  expect_setequal(bnM[[2]]$rnm, bnM[[1]]$rnm)
  expect_setequal(bnM[[2]]$clm, bnM[[1]]$clm)
})
