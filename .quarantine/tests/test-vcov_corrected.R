# Quarantined with R/tspa_corrected_se.R (see _PLAN_QUARANTINE.md).
# Two extractions, each self-contained with its setup copied verbatim:
#  - tests/testthat/test-tspa_render.R: vcov_corrected() multigroup test
#    (lines 322-344 as of 2026-08-17) plus mod2g (line 17).
#  - tests/testthat/test-get_fs_priors.R: vcov_corrected() prior test
#    (lines 175-205) plus single-group setup (lines 1-18).

library(lavaan)
library(lme4)

## From tests/testthat/test-tspa_render.R -------------------------------------

mod2g <- "visual =~ x1 + x2 + x3\nspeed =~ x7 + x8 + x9"

test_that("tspa_call stays re-callable and vcov_corrected() runs on an MG fit", {
  fs_v <- get_fs(HolzingerSwineford1939, model = mod2g, std.lv = TRUE,
                 group = "school", vfsLT = TRUE, format = "list")
  # vcov_corrected() re-evaluates the recorded tspa() call, so the referenced
  # objects must live in globalenv during the test (mirroring vignette-style
  # top-level usage).
  gobjs <- list(fs_dat_v = do.call(rbind, fs_v),
                fsL_v = attr(fs_v, "fsL"),
                fsT_v = attr(fs_v, "fsT"))
  for (nm in names(gobjs)) {
    assign(nm, gobjs[[nm]], envir = globalenv())
  }
  on.exit(for (nm in names(gobjs)) rm(list = nm, envir = globalenv()),
          add = TRUE)
  fit <- tspa("visual ~ speed", data = fs_dat_v, fsT = fsT_v, fsL = fsL_v,
              group = "school")
  vc <- vcov_corrected(fit, vfsLT = attr(fs_v, "vfsLT"))
  expect_s3_class(vc, "matrix")
  expect_equal(dim(vc), dim(vcov(fit)))
  expect_true(all(is.finite(vc)))
  expect_equal(vc, t(vc), tolerance = 1e-10)
  expect_true(all(diag(vc) > 0))
})

## From tests/testthat/test-get_fs_priors.R ------------------------------------

########## Single-group example ##########

prior_model <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
'

prior_fit <- cfa(prior_model, data = PoliticalDemocracy)
prior_est <- lavInspect(prior_fit, what = "est")
prior_data <- lavInspect(prior_fit, what = "data")
prior_lv_names <- colnames(prior_est$lambda)

pm <- c(ind60 = 0.3, dem60 = -0.4)
pc <- matrix(c(1.2, 0.25, 0.25, 0.8), 2, 2,
             dimnames = list(prior_lv_names, prior_lv_names))

test_that("vcov_corrected() works with prior-adjusted factor scores", {
  fs_dat <- get_fs(prior_fit, prior_mean = pm, prior_cov = pc,
                   corrected_fsT = TRUE, vfsLT = TRUE,
                   format = "list")
  # vcov_corrected() re-evaluates the recorded tspa() call from the package
  # namespace, so the referenced objects must live in globalenv during the
  # test (mirroring vignette-style top-level usage).
  gobjs <- list(
    fs_dat_p = fs_dat,
    fsL_p = attr(fs_dat, "fsL"),
    fsT_p = attr(fs_dat, "fsT"),
    fsb_p = attr(fs_dat, "fsb")
  )
  for (nm in names(gobjs)) {
    assign(nm, gobjs[[nm]], envir = globalenv())
  }
  on.exit(
    for (nm in names(gobjs)) rm(list = nm, envir = globalenv()),
    add = TRUE
  )
  tspa_fit <- tspa("dem60 ~ ind60",
                   data = fs_dat_p,
                   fsT = fsT_p,
                   fsL = fsL_p,
                   fsb = fsb_p)
  vc <- vcov_corrected(tspa_fit,
                       vfsLT = attr(fs_dat, "vfsLT"))
  expect_s3_class(vc, "matrix")
  expect_equal(dim(vc), c(3, 3))
  expect_true(all(diag(vc) > 0))
})
