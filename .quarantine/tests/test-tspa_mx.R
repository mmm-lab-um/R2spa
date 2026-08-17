# Quarantined with R/tspa_mx.R (see _PLAN_QUARANTINE.md).
# Extracted from tests/testthat/test-tspa.R, "# Compare to Mx" block
# (lines 159-239 as of 2026-08-17) plus its setup (lines 116-157).
# The setup is copied verbatim and is also retained in test-tspa.R; the
# copy below keeps this file self-contained for re-integration.

# Loading packages and functions
library(lavaan)
library(OpenMx)
if (requireNamespace("umx", quietly = TRUE)) {
  library(umx)
}

# Example 2: Single group with three variables

# CFA model
cfa_3var1 <- '
                            # latent variables
                            ind60 =~ x1 + x2 + x3
                           '
cfa_3var2 <- '
                            # latent variables
                            dem60 =~ y1 + y2 + y3 + y4
                           '
cfa_3var3 <- '
                            # latent variables
                            dem65 =~ y5 + y6 + y7 + y8
                           '

# get factor scores
fs_3var1 <- get_fs(PoliticalDemocracy, cfa_3var1, format = "list")
fs_3var2 <- get_fs(PoliticalDemocracy, cfa_3var2, format = "list")
fs_3var3 <- get_fs(PoliticalDemocracy, cfa_3var3, format = "list")
fs_dat_3var <- cbind(fs_3var1, fs_3var2, fs_3var3)

sem_model_3var <- '
                           # latent variables (indicated by factor scores)
                             ind60 =~ x1 + x2 + x3
                             dem60 =~ y1 + y2 + y3 + y4
                             dem65 =~ y5 + y6 + y7 + y8
                           # regressions
                             dem60 ~ ind60
                             dem65 ~ ind60 + dem60
                        '

sem_3var <- sem(model = sem_model_3var, data  = PoliticalDemocracy)

# tspa model
tspa_3var <- tspa(
  model = "dem60 ~ ind60
               dem65 ~ ind60 + dem60",
  data = fs_dat_3var,
  se_fs = c(
    ind60 = 0.1213615,
    dem60 = 0.6756472,
    dem65 = 0.5724405
  )
)

# Compare to Mx
if (requireNamespace("umx", quietly = TRUE)) {
  model_umx <- umxLav2RAM("
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  dem65 + dem60 + ind60 ~ 1
  ", printTab = FALSE)
# Loading
matL <- mxMatrix(
  type = "Iden", nrow = 3,
  free = FALSE,
  name = "L"
)
# Error
matE <- mxMatrix(
  type = "Diag", nrow = 3, ncol = 3,
  free = FALSE,
  values = c(0.6756472, 0.5724405, 0.1213615)^2,
  name = "E"
)
tspa_mx <- tspa_mx_model(model_umx, data = fs_dat_3var,
                         mat_ld = matL, mat_ev = matE,
                         fs_lv_names = c(ind60 = "fs_ind60",
                                         dem60 = "fs_dem60",
                                         dem65 = "fs_dem65"))
tspa_mx_fit <- mxRun(tspa_mx)
# Check same coefficients and standard errors
test_that("test same regression coefficients with Mx", {
  expect_equal(
    coef(tspa_mx_fit)[c(2, 3, 1, 6, 4, 5)],
    expected = coef(tspa_3var),
    tolerance = 1e-5,
    ignore_attr = TRUE
  )
})
test_that("test same standard errors with Mx", {
  vc_mx <- diag(vcov(tspa_mx_fit))
  vc_lavaan <- diag(vcov(tspa_3var))
  expect_equal(
    vc_mx[c(2, 3, 1, 6, 4, 5)],
    expected = vc_lavaan,
    tolerance = 1e-4,
    ignore_attr = TRUE
  )
})
# Use numeric matrices
tspa_mx2 <- tspa_mx_model(
  model_umx,
  data = fs_dat_3var,
  mat_ld = diag(3) |>
    `dimnames<-`(list(
      c("fs_ind60", "fs_dem60", "fs_dem65"),
      c("ind60", "dem60", "dem65")
    )),
  mat_ev = diag(c(0.1213615, 0.6756472, 0.5724405)^2) |>
    `dimnames<-`(rep(list(c("fs_ind60", "fs_dem60", "fs_dem65")), 2))
)
tspa_mx_fit2 <- mxRun(tspa_mx2)
# Use column names for VC
err_cov <- matrix(c("ev_fs_ind60", NA, NA,
                    NA, "ev_fs_dem60", NA,
                    NA, NA, "ev_fs_dem65"), nrow = 3) |>
  `dimnames<-`(rep(list(c("fs_ind60", "fs_dem60", "fs_dem65")), 2))
tspa_mx3 <- tspa_mx_model(model_umx, data = fs_dat_3var,
                          mat_ld = matL, mat_ev = err_cov,
                          fs_lv_names = c(ind60 = "fs_ind60",
                                          dem60 = "fs_dem60",
                                          dem65 = "fs_dem65"))
tspa_mx_fit3 <- mxRun(tspa_mx3)
test_that("Same results with different Mx matrices input", {
  expect_equal(
    coef(tspa_mx_fit2),
    expected = coef(tspa_mx_fit)
  )
  expect_equal(
    coef(tspa_mx_fit3),
    expected = coef(tspa_mx_fit),
    tolerance = 1e-5
  )
})
}

# ---------------------------------------------------------------------------
# Quarantined with R/tspa_mx.R (see _PLAN_QUARANTINE.md).
# Extracted from tests/testthat/test-get_fscore.R, umx/OpenMx missing-data
# block (lines 675-700 as of 2026-08-17) plus minimal setup: hs_model_2
# (line 240), NA-injected hs (lines 637-642), cfa_fiml + a2 (lines 650-654).
# ---------------------------------------------------------------------------

# Prepare for test objects
hs_model_2 <- ' visual =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed =~ x7 + x8 + x9 '

hs <- HolzingerSwineford1939
# introduce missing data
set.seed(1334)
hs[!rbinom(301, size = 1, prob = 0.7), 7] <- NA
hs[!rbinom(301, size = 1, prob = 0.7), 8] <- NA
hs[!rbinom(301, size = 1, prob = 0.7), 9] <- NA

cfa_fiml <- cfa(
  model = hs_model_2, data = hs, missing = "fiml",
  estimator = "MLR"
)
a2 <- augment_lav_predict(cfa_fiml, method = "Bartlett")

if (requireNamespace("umx", quietly = TRUE)) {
  lcov_umx <- umxLav2RAM(
  "
    visual ~~ textual + speed
    textual ~~ speed
    visual + textual + speed ~ 1
  ",
  printTab = FALSE
)
tspab_mx <- tspa_mx_model(lcov_umx,
  data = a2,
  mat_ld = attr(a2, which = "ld"),
  mat_ev = attr(a2, which = "ev")
)
# Run OpenMx
tspab_mx_fit <- mxRun(tspab_mx)

test_that("tspa_mx() gives similar results as lavaan with missing data",
  code = {
    expect_equal(tspab_mx_fit$m1$S$values,
                 expected = lavInspect(cfa_fiml, what = "est")$psi,
                 tolerance = 1e-5,
                 ignore_attr = TRUE)
  }
)
}
