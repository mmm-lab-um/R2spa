# =====================================================================
# tspa_mx_model(): OpenMx two-stage path analysis (no umx).
#
# A/B gates: the OpenMx fit is validated against lavaan::tspa() (the identical
# 2S-PA model in lavaan) for the fixed / per-row-standard-error form, and
# against the fixed-value equivalent of the per-row definition-variable form
# fed by fs_indiv(). Tolerance is tight (FIML raw likelihood in both) except
# where a dataset's corrected latent variance is negative (then OpenMx reports
# a non-convex Hessian but still reaches the same point estimate).
# =====================================================================

library(lavaan)
library(OpenMx)

# --- coefficient extraction by (from -> to), robust to RAM var ordering ----
mx_path_val <- function(m, from, to, model = "m1") {
  v <- c(m$manifestVars, m$latentVars)
  unname(coef(m)[sprintf("%s.A[%d,%d]", model, match(to, v), match(from, v))])
}
mx_var_val <- function(m, x, model = "m1") {
  v <- c(m$manifestVars, m$latentVars)
  unname(coef(m)[sprintf("%s.S[%d,%d]", model, match(x, v), match(x, v))])
}

# --- Setup: 2-variable and 3-variable PoliticalDemocracy CFAs ---------------
cfa_ind <- "ind60 =~ x1 + x2 + x3"
cfa_dem <- "dem60 =~ y1 + y2 + y3 + y4"
cfa_dem65 <- "dem65 =~ y5 + y6 + y7 + y8"
se2 <- c(ind60 = 0.1213615, dem60 = 0.6756472)
se3 <- c(ind60 = 0.1213615, dem60 = 0.6756472, dem65 = 0.5724405)

fsd2 <- get_fs(lavaan::PoliticalDemocracy,
               paste(cfa_ind, cfa_dem, sep = "\n"), format = "unified")
fsd3 <- get_fs(lavaan::PoliticalDemocracy,
               paste(cfa_ind, cfa_dem, cfa_dem65, sep = "\n"), format = "unified")
model2 <- "dem60 ~ ind60"
model3 <- "dem60 ~ ind60\ndem65 ~ ind60 + dem60"

test_that("T1a: 2-var se_fs model matches lavaan::tspa() (A/B)", {
  m <- suppressWarnings(tspa_mx_model(model2, data = fsd2, se_fs = se2))
  ref <- tspa(model = model2, data = fsd2, se_fs = se2)
  expect_equal(mx_path_val(m, "ind60", "dem60"),
               as.numeric(coef(ref)["dem60~ind60"]), tolerance = 1e-4)
  expect_equal(mx_var_val(m, "ind60"),
               as.numeric(coef(ref)["ind60~~ind60"]), tolerance = 1e-4)
  expect_equal(mx_var_val(m, "dem60"),
               as.numeric(coef(ref)["dem60~~dem60"]), tolerance = 1e-4)
})

test_that("T1b: 3-var se_fs model matches lavaan::tspa() (A/B)", {
  m <- suppressWarnings(tspa_mx_model(model3, data = fsd3, se_fs = se3))
  ref <- tspa(model = model3, data = fsd3, se_fs = se3)
  expect_equal(mx_path_val(m, "ind60", "dem60"),
               as.numeric(coef(ref)["dem60~ind60"]), tolerance = 1e-4)
  expect_equal(mx_path_val(m, "ind60", "dem65"),
               as.numeric(coef(ref)["dem65~ind60"]), tolerance = 1e-4)
  expect_equal(mx_path_val(m, "dem60", "dem65"),
               as.numeric(coef(ref)["dem65~dem60"]), tolerance = 1e-4)
  # dem65's corrected variance is negative in this data; the point estimate
  # still matches the reference to a loose tolerance.
  expect_equal(mx_var_val(m, "dem65"),
               as.numeric(coef(ref)["dem65~~dem65"]), tolerance = 1e-3)
})

# --- fixed fsL/fsT (identity loadings + diag se^2) is the same model --------
test_that("T2: fixed fsL/fsT (identity + diag se^2) matches the se_fs form", {
  Lm <- `dimnames<-`(diag(2),
                     list(c("fs_ind60", "fs_dem60"),
                          c("ind60", "dem60")))
  Tm <- `dimnames<-`(diag(unname(se2)^2),
                     rep(list(c("fs_ind60", "fs_dem60")), 2))
  m_fix <- suppressWarnings(tspa_mx_model(model2, data = fsd2,
                                          fsL = Lm, fsT = Tm))
  m_se  <- suppressWarnings(tspa_mx_model(model2, data = fsd2, se_fs = se2))
  for (p in list(c("ind60", "dem60"))) {
    expect_equal(mx_path_val(m_fix, p[1], p[2]),
                 mx_path_val(m_se, p[1], p[2]), tolerance = 1e-4)
  }
  for (v in c("ind60", "dem60")) {
    expect_equal(mx_var_val(m_fix, v), mx_var_val(m_se, v), tolerance = 1e-4)
  }
})

# --- per-row definition-variable form fed by fs_indiv() (data flow) ---------
cfa2fit <- cfa(paste(cfa_ind, cfa_dem, sep = "\n"),
               data = lavaan::PoliticalDemocracy)
fs2 <- get_fs(cfa2fit)
dat2 <- fs_indiv(fs2, include_intercept = TRUE)
L2 <- `dimnames<-`(diag(2),
                   list(c("fs_ind60", "fs_dem60"),
                        c("ind60", "dem60")))
T2dv <- matrix(c("ev_fs_ind60", NA, NA, "ev_fs_dem60"),
               nrow = 2, dimnames = list(c("fs_ind60", "fs_dem60"),
                                         c("fs_ind60", "fs_dem60")))
b2dv <- c(fs_ind60 = "int_fs_ind60", fs_dem60 = "int_fs_dem60")

test_that("T9: fs_indiv() output columns are the data-flow for def vars", {
  # the definition-variable columns named in fsT/fsb all exist in the table
  dv <- c(unname(unlist(T2dv)), unname(b2dv))
  dv <- dv[!is.na(dv)]
  expect_true(all(dv %in% colnames(dat2)))
  # and they are NA-free (OpenMx definition variables must be complete)
  expect_true(all(!is.na(dat2[, dv, drop = FALSE])))
})

test_that("T4: per-row def-var fit equals its fixed-value equivalent", {
  m_dv <- suppressWarnings(tspa_mx_model(model2, data = dat2,
                                         fsL = L2, fsT = T2dv, fsb = b2dv))
  T2fix <- matrix(c(unname(mean(dat2$ev_fs_ind60)), 0,
                    0, unname(mean(dat2$ev_fs_dem60))),
                  nrow = 2, dimnames = dimnames(T2dv))
  b2fix <- c(fs_ind60 = unname(mean(dat2$int_fs_ind60)),
             fs_dem60 = unname(mean(dat2$int_fs_dem60)))
  m_fx <- suppressWarnings(tspa_mx_model(model2, data = dat2[, 1:2],
                                         fsL = L2, fsT = T2fix, fsb = b2fix))
  expect_equal(mx_path_val(m_dv, "ind60", "dem60"),
               mx_path_val(m_fx, "ind60", "dem60"), tolerance = 1e-6)
  expect_equal(mx_var_val(m_dv, "ind60"),
               mx_var_val(m_fx, "ind60"), tolerance = 1e-6)
  expect_equal(mx_var_val(m_dv, "dem60"),
               mx_var_val(m_fx, "dem60"), tolerance = 1e-6)
})

# --- missing-data FIML: recover the lavaan FIML latent structure ------------
hs <- HolzingerSwineford1939
set.seed(1334)
hs[!rbinom(301, size = 1, prob = 0.7), 7] <- NA
hs[!rbinom(301, size = 1, prob = 0.7), 8] <- NA
hs[!rbinom(301, size = 1, prob = 0.7), 9] <- NA
mod_hs <- "visual =~ x1 + x2 + x3\nspeed  =~ x7 + x8 + x9"
cfa_hs <- suppressWarnings(cfa(mod_hs, data = hs, missing = "fiml"))
# Per-row (FIML) form goes through the fs_indiv() definition-variable data
# flow; the per-pattern fsL/fsT list is resolved to one row per observation.
fs_hs <- get_fs(cfa_hs)
dat_hs <- fs_indiv(fs_hs, include_intercept = TRUE)

test_that("T8: missing-data per-row def-var fit runs via fs_indiv (FIML over scores)", {
  Lh <- `dimnames<-`(diag(2),
                     list(c("fs_visual", "fs_speed"),
                          c("visual", "speed")))
  Th <- matrix(c("ev_fs_visual", NA, NA, "ev_fs_speed"), nrow = 2,
               dimnames = list(c("fs_visual", "fs_speed"),
                               c("fs_visual", "fs_speed")))
  bh <- c(fs_visual = "int_fs_visual", fs_speed = "int_fs_speed")
  expect_true(all(c("ev_fs_visual", "ev_fs_speed", "int_fs_visual",
                    "int_fs_speed") %in% colnames(dat_hs)))
  expect_no_error(
    m <- suppressWarnings(tspa_mx_model("visual ~ speed", data = dat_hs,
                                        fsL = Lh, fsT = Th, fsb = bh))
  )
  expect_true(all(is.finite(coef(m))))
  expect_gt(mx_var_val(m, "visual"), 0)
  expect_gt(mx_var_val(m, "speed"), 0)
})

# --- input guards ------------------------------------------------------------
test_that("guards: model and measurement inputs are validated", {
  expect_error(tspa_mx_model(123, data = fsd2, se_fs = se2),
               "must be a lavaan syntax string")
  expect_error(tspa_mx_model(model2, data = 1:3, se_fs = se2),
               "must be a data frame")
  # a structural latent must be a corrected latent
  expect_error(tspa_mx_model("dem60 ~ x1", data = fsd2, se_fs = se2),
               "must be a corrected latent")
  # score column missing from the data
  expect_error(
    tspa_mx_model("dem60 ~ ind60", data = fsd2[, colnames(fsd2) != "fs_ind60"],
                  se_fs = se2),
    "missing factor-score column"
  )
  # definition-variable column missing
  expect_error(
    tspa_mx_model(model2, data = dat2[, 1:2],
                  fsL = L2, fsT = T2dv, fsb = b2dv),
    "missing definition-variable column"
  )
  # definition-variable column with NA
  da <- dat2; da$ev_fs_ind60[1] <- NA
  expect_error(
    suppressWarnings(tspa_mx_model(model2, data = da,
                                   fsL = L2, fsT = T2dv, fsb = b2dv)),
    "contain NA"
  )
})

test_that("guards: multigroup fsL/fsT are rejected (Phase 1 single-group)", {
  lm <- list(`p1` = L2, `p2` = L2)
  tm <- list(`p1` = T2dv, `p2` = T2dv)
  expect_error(
    tspa_mx_model(model2, data = dat2[, 1:2], fsL = lm, fsT = tm),
    "not supported yet"
  )
})

test_that("guards: incomplete or mistyped measurement inputs are rejected", {
  # NA in se_fs (an NA error variance would otherwise yield a RED-status fit)
  expect_error(
    tspa_mx_model(model2, data = fsd2, se_fs = c(ind60 = NA, dem60 = 0.6756472)),
    "must not contain NA"
  )
  t_ok <- matrix(c(0.25, 0, 0, 0.4), nrow = 2,
                 dimnames = list(c("fs_ind60", "fs_dem60"),
                                 c("fs_ind60", "fs_dem60")))
  # NA on the fsT diagonal (a score without a known error variance)
  t_na <- t_ok; t_na[2, 2] <- NA
  expect_error(
    tspa_mx_model(model2, data = dat2[, 1:2], fsL = L2, fsT = t_na),
    "error variance"
  )
  # fsT without column names (previously: opaque 'subscript out of bounds')
  t_nc <- t_ok; colnames(t_nc) <- NULL
  expect_error(
    tspa_mx_model(model2, data = dat2[, 1:2], fsL = L2, fsT = t_nc),
    "columns must be named"
  )
  # non-numeric definition-variable column
  dn <- dat2; dn$ev_fs_ind60 <- as.character(dn$ev_fs_ind60)
  expect_error(
    suppressWarnings(tspa_mx_model(model2, data = dn,
                                   fsL = L2, fsT = T2dv, fsb = b2dv)),
    "must be numeric"
  )
  # a score loading on no latent (previously dropped silently)
  l_na <- L2; l_na[1, ] <- NA
  expect_error(
    tspa_mx_model(model2, data = dat2[, 1:2], fsL = l_na, fsT = t_ok),
    "load on no latent"
  )
  # zero-row data
  expect_error(
    tspa_mx_model(model2, data = fsd2[0, ], se_fs = se2),
    "at least one row"
  )
})

test_that("auto-seed release is user-gated (user == 0): a user-fixed zero mean stays fixed", {
  spec <- R2spa:::tspa_mx_spec(se2, NULL, NULL, NULL)
  vcol <- function(m, x) {
    v <- c(m$manifestVars, m$latentVars)
    match(x, v)
  }
  # user-declared 'ind60 ~ 0*1' (fixed zero mean): must not be re-freed
  ms <- R2spa:::tspa_mx_model_string("dem60 ~ ind60; ind60 ~ 0*1", spec)
  m_fix <- suppressWarnings(mxRun(R2spa:::lav_to_mx_ram(ms, spec, fsd2),
                                  silent = TRUE))
  i <- vcol(m_fix, "ind60")
  expect_false(sprintf("m1.M[1,%d]", i) %in% names(coef(m_fix)))
  # the auto-seeded latent variance is still released (free)
  expect_true(sprintf("m1.S[%d,%d]", i, i) %in% names(coef(m_fix)))
  # control: no user mean structure -> the latent mean is free (need_mean)
  m_free <- suppressWarnings(mxRun(
    R2spa:::lav_to_mx_ram(R2spa:::tspa_mx_model_string("dem60 ~ ind60", spec),
                          spec, fsd2), silent = TRUE))
  j <- vcol(m_free, "ind60")
  expect_true(sprintf("m1.M[1,%d]", j) %in% names(coef(m_free)))
})

# --- mean-structure convention (documented): same model, different split ----
# Both routes fit the same model, so the covariance/structural quantities
# (paths, latent variances) agree. The unidentifiable mean split differs:
# tspa() (lavaan) fixes the exogenous latent mean at 0 and estimates the
# factor-score mean, whereas tspa_mx_model() (OpenMx) fixes the score
# residual mean at 0 and estimates the latent mean. The same data value ends
# up in the latent mean (OpenMx) vs. the factor-score mean (lavaan).
test_that("mean-structure convention: same model, latent vs score mean split", {
  set.seed(11)
  n <- 2000
  fa <- rnorm(n, 1.5, 1)
  a  <- fa + rnorm(n, 0.5)
  b  <- 2*a + rnorm(n, 0.8)
  dat <- data.frame(fs_a = fa, fs_b = b)
  se  <- c(a = 0.5, b = 0.8)
  mod <- "b ~ a; b ~ 1"
  # only the mean split is overparameterized (the total is identifiable), so
  # lavaan warns on vcov(); the point estimates are deterministic and tested.
  la <- suppressWarnings(tspa(model = mod, data = dat, se_fs = se))
  m  <- suppressWarnings(tspa_mx_model(mod, data = dat, se_fs = se))
  v  <- c(m$manifestVars, m$latentVars)
  # covariance/structural quantities agree to optimizer tolerance:
  expect_equal(mx_path_val(m, "a", "b"),
               as.numeric(coef(la)["b~a"]), tolerance = 1e-4)
  expect_equal(mx_var_val(m, "a"),
               as.numeric(coef(la)["a~~a"]), tolerance = 1e-4)
  expect_equal(mx_var_val(m, "b"),
               as.numeric(coef(la)["b~~b"]), tolerance = 1e-4)
  # the mean split: the data value (~1.5) is the LATENT mean in OpenMx and
  # the factor-SCORE mean in lavaan:
  expect_equal(unname(m$M$values[1, match("a", v)]),
               as.numeric(coef(la)["fs_a~1"]), tolerance = 1e-4)
  # ...and OpenMx pins the score's own (residual) mean at zero:
  expect_equal(unname(m$M$values[1, match("fs_a", v)]), 0)
})
