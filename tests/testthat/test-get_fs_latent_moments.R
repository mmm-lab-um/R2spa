# get_fs() latent-moment attributes: effective (prior-adjusted) psi / alpha.
#
# get_fs() now attaches the effective latent covariance (psi) and latent means
# (alpha) as attributes, group-level, with shape mirroring fsT. These tests
# cover: single-group no-mean-structure, prior-adjusted (regression) values,
# multi-group (per-group estimates and shared priors, shape vs fsT, unified
# vs list nesting), the merMod case, and the fs_to_group_list() round-trip.
library(lavaan)
library(lme4)

# ---- fixtures -----------------------------------------------------------
# SG, no (estimated) mean structure (q = 1)
hs_model_1f <- "visual =~ x1 + x2 + x3"
fit_1f <- cfa(hs_model_1f, data = HolzingerSwineford1939)
fs_1f <- get_fs(fit_1f)

# Priors (regression only) on PoliticalDemocracy, 2 factors, no mean structure
pri_model <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
'
pri_fit <- cfa(pri_model, data = PoliticalDemocracy)
pri_lv <- c("ind60", "dem60")
pm2 <- c(ind60 = 0.3, dem60 = -0.4)
pc2 <- matrix(c(1.2, 0.25, 0.25, 0.8), 2, 2, dimnames = list(pri_lv, pri_lv))

# MG on HS, 2 factors, group = school
mg_model <- 'ind =~ x1 + x2 + x3
             dem =~ x4 + x5 + x6 + x7'
fit_mg <- cfa(mg_model, data = HolzingerSwineford1939, group = "school")

# merMod (2 random effects per cluster)
lmod <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fs_mer <- get_fs(lmod)

# ---- helpers ------------------------------------------------------------
# Single-group unified (or single-group list) attribute value.
sg_attr <- function(fs, ak) attr(fs, ak)[[1L]]

# ============================================================================
# 1. Single-group, no mean structure
# ============================================================================

test_that("get_fs(): SG no-mean-structure -- alpha is the named zero vector, psi == est$psi", {
  est <- lavaan::lavInspect(fit_1f, what = "est")
  q <- ncol(est$lambda)       # number of latent factors
  expect_null(est$alpha)   # no estimated mean structure
  expect_equal(length(attr(fs_1f, "psi")), 1L)   # unified single-group wrapper
  psi <- sg_attr(fs_1f, "psi")
  alpha <- sg_attr(fs_1f, "alpha")
  expect_equal(as.matrix(psi), as.matrix(est$psi), tolerance = 1e-8, ignore_attr = TRUE)
  expect_length(alpha, q)
  expect_true(all(unname(alpha) == 0))
  expect_named(alpha, colnames(est$lambda))  # the latent (factor) columns
})

# ============================================================================
# 2. Priors (regression)
# ============================================================================

test_that("get_fs(): prior_mean + prior_cov (regression) -- psi == prior_cov, alpha == prior_mean", {
  fs_pp <- get_fs(pri_fit, prior_mean = pm2, prior_cov = pc2, format = "unified")
  expect_equal(as.matrix(sg_attr(fs_pp, "psi")), as.matrix(pc2),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_named(sg_attr(fs_pp, "alpha"), pri_lv)
  expect_equal(unname(sg_attr(fs_pp, "alpha")), unname(pm2), tolerance = 1e-10)
})

test_that("get_fs(): prior_cov only -- psi == prior_cov, alpha == named zeros (no mean structure)", {
  est <- lavaan::lavInspect(pri_fit, what = "est")
  expect_null(est$alpha)
  fs_pc <- get_fs(pri_fit, prior_cov = pc2, format = "unified")
  expect_equal(as.matrix(sg_attr(fs_pc, "psi")), as.matrix(pc2),
               tolerance = 1e-8, ignore_attr = TRUE)
  alpha <- sg_attr(fs_pc, "alpha")  # named zero vector when no mean structure
  expect_length(alpha, length(pri_lv))
  expect_true(all(unname(alpha) == 0))
  expect_named(alpha, pri_lv)
})

test_that("get_fs(): prior_mean only -- psi == est$psi, alpha == prior_mean", {
  est <- lavaan::lavInspect(pri_fit, what = "est")
  fs_pm <- get_fs(pri_fit, prior_mean = pm2, format = "unified")
  expect_equal(as.matrix(sg_attr(fs_pm, "psi")), as.matrix(est$psi),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_named(sg_attr(fs_pm, "alpha"), pri_lv)
  expect_equal(unname(sg_attr(fs_pm, "alpha")), unname(pm2), tolerance = 1e-10)
})

# ============================================================================
# 3. Multi-group
# ============================================================================

test_that("get_fs(): MG -- per-group psi/alpha equal the per-group estimates; shape mirrors fsT", {
  est_g <- lavaan::lavInspect(fit_mg, what = "est", drop.list.single.group = FALSE)
  fs_u <- get_fs(fit_mg, format = "unified")
  grp_labels <- names(attr(fs_u, "fsT"))
  expect_equal(length(grp_labels), 2L)
  # shape mirrors fsT: named list by group label, same labels as fsT
  expect_identical(names(attr(fs_u, "psi")), grp_labels)
  expect_identical(names(attr(fs_u, "alpha")), grp_labels)
  expect_equal(length(attr(fs_u, "psi")), 2L)
  for (g in grp_labels) {
    expect_equal(as.matrix(attr(fs_u, "psi")[[g]]), as.matrix(est_g[[g]]$psi),
                 tolerance = 1e-8, ignore_attr = TRUE)
    # no mean structure -> each group's alpha is the named zero vector
    ag <- attr(fs_u, "alpha")[[g]]
    expect_true(all(unname(ag) == 0))
    expect_named(ag, colnames(est_g[[g]]$lambda))
  }
})

test_that("get_fs(): MG with priors -- every group's psi/alpha equals the shared prior", {
  est_g <- lavaan::lavInspect(fit_mg, what = "est", drop.list.single.group = FALSE)
  pm_g <- c(ind = 0.5, dem = -0.4)
  pc_g <- matrix(c(1.5, 0.3, 0.3, 0.9), 2, 2, dimnames = list(c("ind", "dem"), c("ind", "dem")))
  fs_u <- get_fs(fit_mg, prior_mean = pm_g, prior_cov = pc_g, format = "unified")
  for (g in names(attr(fs_u, "psi"))) {
    expect_equal(as.matrix(attr(fs_u, "psi")[[g]]), as.matrix(pc_g),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(unname(attr(fs_u, "alpha")[[g]]), unname(pm_g), tolerance = 1e-10)
    # and the prior is distinct from the per-group estimate (so this is meaningful)
    expect_false(isTRUE(all.equal(as.matrix(attr(fs_u, "psi")[[g]]),
                                  as.matrix(est_g[[g]]$psi), tolerance = 1e-6)))
  }
})

test_that("get_fs(): MG list format -- psi/alpha on each element AND as a list attribute on the outer list", {
  fs_l <- get_fs(fit_mg, format = "list")
  # list-valued attribute on the outer list, named by group
  expect_equal(length(attr(fs_l, "psi")), 2L)
  expect_equal(length(attr(fs_l, "alpha")), 2L)
  expect_identical(names(attr(fs_l, "psi")), names(fs_l))
  expect_identical(names(attr(fs_l, "alpha")), names(fs_l))
  # direct attribute on each group data frame, and it equals the outer entry
  for (g in names(fs_l)) {
    expect_true(is.matrix(attr(fs_l[[g]], "psi")))
    expect_true(is.numeric(attr(fs_l[[g]], "alpha")))
    expect_equal(as.matrix(attr(fs_l, "psi")[[g]]), as.matrix(attr(fs_l[[g]], "psi")),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(unname(attr(fs_l, "alpha")[[g]]), unname(attr(fs_l[[g]], "alpha")),
                 tolerance = 1e-10)
  }
})

# ============================================================================
# 4. merMod
# ============================================================================

test_that("get_fs(): merMod -- psi == VarCorr()[[1]] renamed to the RE/score names; alpha named zeros", {
  pc1 <- as.matrix(lme4::VarCorr(lmod)[[1L]])
  psi <- attr(fs_mer, "psi")
  # values equal the first random-effects term's covariance
  expect_equal(unname(psi), unname(pc1), tolerance = 1e-8)
  # dimnames renamed to the re/score names so they align with the fsL columns
  expect_identical(rownames(psi), c("u0", "u1"))
  expect_identical(colnames(psi), c("u0", "u1"))
  # random effects are mean zero
  alpha <- attr(fs_mer, "alpha")
  expect_length(alpha, 2L)
  expect_true(all(unname(alpha) == 0))
  expect_named(alpha, c("u0", "u1"))
})

# ============================================================================
# 5. fs_to_group_list() round-trip
# ============================================================================

test_that("fs_to_group_list(): psi/alpha survive the unified <-> list conversion", {
  # multi-group: unified -> list, per-group values preserved
  uni <- get_fs(fit_mg, format = "unified")
  lst <- fs_to_group_list(uni)
  expect_true(is.list(lst) && !is.data.frame(lst))
  expect_identical(names(attr(lst, "psi")), names(lst))
  for (g in names(lst)) {
    expect_equal(as.matrix(attr(lst[[g]], "psi")), as.matrix(attr(uni, "psi")[[g]]),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(unname(attr(lst[[g]], "alpha")), unname(attr(uni, "alpha")[[g]]),
                 tolerance = 1e-10)
  }
  # and back: list -> unified, values preserved
  uni_back <- fs_to_group_list(lst)
  expect_true(is.data.frame(uni_back))
  for (g in names(lst)) {
    expect_equal(as.matrix(attr(uni_back, "psi")[[g]]), as.matrix(attr(uni, "psi")[[g]]),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(unname(attr(uni_back, "alpha")[[g]]), unname(attr(uni, "alpha")[[g]]),
                 tolerance = 1e-10)
  }

  # single-group: unified is a data frame; the converter unwraps to a plain value
  uni1 <- get_fs(fit_1f, format = "unified")
  lst1 <- fs_to_group_list(uni1)
  expect_true(is.data.frame(lst1))
  psi_before <- attr(uni1, "psi")[[1L]]
  psi_after <- attr(lst1, "psi")
  if (is.list(psi_after)) psi_after <- psi_after[[1L]]
  expect_equal(as.matrix(psi_before), as.matrix(psi_after), tolerance = 1e-8, ignore_attr = TRUE)
  alpha_before <- attr(uni1, "alpha")[[1L]]
  alpha_after <- attr(lst1, "alpha")
  if (is.list(alpha_after)) alpha_after <- alpha_after[[1L]]
  expect_equal(unname(alpha_before), unname(alpha_after), tolerance = 1e-10)
})
