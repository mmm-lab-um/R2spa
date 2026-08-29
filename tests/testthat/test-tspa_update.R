# =====================================================================
# lavaan::update() on a tspa() fit.
#
# Regression tests for the self-contained @call fix: tspa() now inlines
# the internally-built tspaModel / data / se into the fit's @call, so
# lavaan::update(fit, ...) re-evaluates the stored call without
# referencing tspa()-local symbols (before the fix: update(fit,
# meanstructure = TRUE) failed with `object 'se' not found`). An
# updated fit must be bit-identical to calling tspa(...,
# meanstructure = TRUE) from the start.
# =====================================================================

library(lavaan)

# --- Case 1: base 3-factor (multi-factor path) --------------------------
f1 <- get_fs(lavaan::PoliticalDemocracy, model = "ind60 =~ x1 + x2 + x3")
f2 <- get_fs(lavaan::PoliticalDemocracy, model = "dem60 =~ y1 + y2 + y3 + y4")
f3 <- get_fs(lavaan::PoliticalDemocracy, model = "dem65 =~ y5 + y6 + y7 + y8")
fs_dat <- cbind(f1, f2, f3)
sefs <- list(ind60 = 0.1213615, dem60 = 0.6756472, dem65 = 0.5724405)
mod <- "dem60 ~ ind60\ndem65 ~ ind60 + dem60"
fit <- tspa(mod, fs_dat, se_fs = sefs)

test_that("U1: base 3-factor @call is self-contained; update() works and matches a fresh meanstructure fit", {
  # The fix: model / data / se are inlined literals in the stored call,
  # not tspa()-local symbols (out of scope when update() re-evaluates).
  expect_false(is.symbol(fit@call[["model"]]))
  expect_false(is.symbol(fit@call[["data"]]))
  expect_false(is.symbol(fit@call[["se"]]))
  # The natural update now works (previously: `object 'se' not found`).
  upd <- update(fit, meanstructure = TRUE)
  expect_s4_class(upd, "lavaan")
  # The issue's exact form: re-feed the inspected data.
  upd2 <- update(fit, data = lavInspect(fit, "data"), meanstructure = TRUE)
  expect_s4_class(upd2, "lavaan")
  # Equivalence: update() must match a fresh meanstructure fit.
  ref <- tspa(mod, fs_dat, se_fs = sefs, meanstructure = TRUE)
  common <- intersect(names(coef(upd)), names(coef(ref)))
  expect_gt(length(common), 0)
  expect_equal(unname(coef(upd)[common]), unname(coef(ref)[common]),
               tolerance = 1e-6)
})

# --- Case 2: growth model (the vignette example) ------------------------
test_that("U2: growth-model @call is self-contained; update() works and matches a fresh meanstructure fit", {
  p <- test_path("../../vignettes/eclsk.rds")
  skip_if(!file.exists(p), "eclsk.rds not shipped")
  eclsk <- readRDS(p)

  # Stage 1 (strict invariance) — vignettes/tspa-growth-vignette.Rmd.
  strict_mod <- "
eta1 =~ 15.1749088 * s_g3 + l2 * r_g3 + l3 * m_g3
eta2 =~ 15.1749088 * s_g5 + l2 * r_g5 + l3 * m_g5
eta3 =~ 15.1749088 * s_g8 + l2 * r_g8 + l3 * m_g8
eta1 ~~ 1 * eta1 + eta2 + eta3
eta2 ~~ eta2 + eta3
eta3 ~~ eta3
s_g3 ~~ u1 * s_g3 + s_g5 + s_g8
s_g5 ~~ u1 * s_g5 + s_g8
s_g8 ~~ u1 * s_g8
r_g3 ~~ u2 * r_g3 + r_g5 + r_g8
r_g5 ~~ u2 * r_g5 + r_g8
r_g8 ~~ u2 * r_g8
m_g3 ~~ u3 * m_g3 + m_g5 + m_g8
m_g5 ~~ u3 * m_g5 + m_g8
m_g8 ~~ u3 * m_g8
eta1 ~ 0 * 1
eta2 ~ 1
eta3 ~ 1
s_g3 ~ i1 * 1
s_g5 ~ i1 * 1
s_g8 ~ i1 * 1
r_g3 ~ i2 * 1
r_g5 ~ i2 * 1
r_g8 ~ i2 * 1
m_g3 ~ i3 * 1
m_g5 ~ i3 * 1
m_g8 ~ i3 * 1
"
  # Stage 2 (growth) — vignettes/tspa-growth-vignette.Rmd.
  growth_mod <- "
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + start(.5) * eta2 + 1 * eta3
eta1 ~~ psi * eta1
eta2 ~~ psi * eta2
eta3 ~~ psi * eta3
i ~~ start(.8) * i
s ~~ start(.5) * s
i ~~ start(0) * s
i ~ 1
s ~ 1
"
  fs_dat <- get_fs(eclsk, model = strict_mod)
  fit <- tspa(growth_mod, fs_dat,
              fsT = attr(fs_dat, "fsT"), fsL = attr(fs_dat, "fsL"),
              fsb = attr(fs_dat, "fsb"), estimator = "ML")

  # The fix: model / data / se are inlined literals in the stored call.
  expect_false(is.symbol(fit@call[["model"]]))
  expect_false(is.symbol(fit@call[["data"]]))
  expect_false(is.symbol(fit@call[["se"]]))
  # The natural update now works (previously: `object 'se' not found`).
  upd <- update(fit, meanstructure = TRUE)
  expect_s4_class(upd, "lavaan")
  # Equivalence: update() must match a fresh meanstructure fit.
  ref <- tspa(growth_mod, fs_dat,
              fsT = attr(fs_dat, "fsT"), fsL = attr(fs_dat, "fsL"),
              fsb = attr(fs_dat, "fsb"), estimator = "ML",
              meanstructure = TRUE)
  common <- intersect(names(coef(upd)), names(coef(ref)))
  expect_gt(length(common), 0)
  expect_equal(unname(coef(upd)[common]), unname(coef(ref)[common]),
               tolerance = 1e-6)
})
