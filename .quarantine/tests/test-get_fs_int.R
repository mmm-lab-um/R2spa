########################### Test get_fs_int function ##############################

test_that("get_fs_int() works with toy data set", code = {
  # Artificial data set
  set.seed(1230)
  dd <- data.frame(
    "a" = rnorm(5),
    "b" = runif(5),
    "w" = 1:5,
    "se_a" = rep(.15, 5),
    "se_b" = rgamma(5, 1, 1),
    "se_w" = rep(.2, 5),
    "ld_a" = 1,
    "ld_b" = c(1, 1, 1, .8, .9),
    "ld_w" = 0.7
  )

  fsint1 <- get_fs_int(dd,
    fs_name = c("a", "b", "w"), se = c("se_a", "se_b", "se_w"),
    loading = c("ld_a", "ld_b", "ld_w"), model = "a:w + a:b"
  )
  expect_equal(cor(fsint1[["a:b"]], dd[["a"]] * dd[["b"]]), 1)
  expect_in(c("a:w", "a:b", "a:b_ld", "a:w_se"), names(fsint1))
})

# Simulate a dat
sample_dat <- data.frame(
  fs1 = rnorm(100),
  fs1_se = rnorm(100),
  fs1_ld = rnorm(100),
  fs2 = rnorm(100),
  fs2_se = rnorm(100),
  fs2_ld = rnorm(100),
  fs3 = rnorm(100),
  fs3_se = rnorm(100),
  fs3_ld = rnorm(100)
)

# Test that data input must be a data frame
test_that("Data input must be a data frame", {
  expect_error(get_fs_int(dat = list(),
                          fs_name = c("fs1", "fs2", "fs3"),
                          se_fs = c("fs1_se", "fs2_se", "fs3_se"),
                          loading_fs = c("fs1_ld", "fs2_ld", "fs3_ld")))
})

# Test for correct input checks
test_that("Inputs are checked for correct types and existence", {
  expect_error(get_fs_int(dat = sample_dat,
                          fs_name = c("fs1", "abc", "fs3"),
                          se_fs = c("fs1_se", "fs2_se", "fs3_se"),
                          loading_fs = c("fs1_ld", "fs2_ld", "fs3_ld")))
})

# Test for handling numeric lat_var and its length
test_that("lat_var must be numeric and match the length of fs_name", {
  expect_error(get_fs_int(dat = sample_dat,
                          fs_name = c("fs1", "fs2", "fs3"),
                          se_fs = c("fs1_se", "fs2_se", "fs3_se"),
                          loading_fs = c("fs1_ld", "fs2_ld", "fs3_ld"),
                          lat_var = "1"))
  expect_error(get_fs_int(dat = sample_dat,
                          fs_name = c("fs1", "fs2"),
                          se_fs = c("fs1_se", "fs2_se"),
                          loading_fs = c("fs1_ld", "fs2_ld"),
                          lat_var = c(1, 2, 3)))
})

# Test function with no `lat_var` or `model`
test_that("Product indicators are correctly calculated", {
  result <- get_fs_int(dat = sample_dat,
                       fs_name = c("fs1", "fs2", "fs3"),
                       se_fs = c("fs1_se", "fs2_se", "fs3_se"),
                       loading_fs = c("fs1_ld", "fs2_ld", "fs3_ld"))
  expect_true("fs1:fs2" %in% names(result))
  expect_true("fs1:fs2_se" %in% names(result))
  expect_true("fs1:fs2_ld" %in% names(result))
  expect_true("fs1:fs3" %in% names(result))
  expect_true("fs1:fs3_se" %in% names(result))
  expect_true("fs1:fs3_ld" %in% names(result))
  expect_true("fs2:fs3" %in% names(result))
  expect_true("fs2:fs3_se" %in% names(result))
  expect_true("fs2:fs3_ld" %in% names(result))
})

# Quarantined with R/get_fs_int.R (see _PLAN_QUARANTINE.md).
# Extracted from tests/testthat/test-tspa_render.R: the product-score
# auto-alias section (lines 346-405 as of 2026-08-17: section header,
# int_setup(), "Product-score columns are auto-aliased" test) and the
# tspa_sf_alias no-op test (lines 419-427).

## Interaction (product-score) auto-alias ---------------------------------------

int_setup <- function(n = 500) {
  set.seed(2116)
  cov_xmz_ey <- matrix(c(1, 0.1, 0.15, 0,
                         0.1, 1, 0.12, 0,
                         0.15, 0.12, 1, 0,
                         0, 0, 0, 0.481351), nrow = 4)
  eta <- as.data.frame(MASS::mvrnorm(n, mu = rep(0, 4),
                                     Sigma = cov_xmz_ey))
  names(eta) <- c("x", "m", "z", "ey")
  eta <- transform(eta, xm = x * m, xz = x * z, mz = m * z)
  etay <- 0.3 * eta$x + 0.4 * eta$m + 0.2 * eta$z +
    0.1 * eta$xm + 0.15 * eta$xz + 0.12 * eta$mz + eta$ey
  lk <- list(x = c(0.9, 0.8, 0.7), m = c(0.85, 0.75, 0.65),
             z = c(0.8, 0.7, 0.6), y = c(0.75, 0.7, 0.65))
  obs <- setNames(lapply(c("x", "m", "z"), function(v0) {
    eta[[v0]] %*% t(lk[[v0]]) + rnorm(n * 3)
  }), c("x", "m", "z"))
  obs$y <- etay %*% t(lk$y) + rnorm(n * 3)
  df <- as.data.frame(do.call(cbind, obs[c("x", "m", "z", "y")]))
  names(df) <- c(paste0("x", 1:3), paste0("m", 1:3), paste0("z", 1:3),
                 paste0("y", 1:3))
  fs_dat <- get_fs(df, model = "x =~ x1 + x2 + x3
                                  m =~ m1 + m2 + m3
                                  z =~ z1 + z2 + z3
                                  y =~ y1 + y2 + y3",
                   std.lv = TRUE, method = "Bartlett")
  ind <- get_fs_int(dat = fs_dat,
                    fs_name = c("fs_x", "fs_m", "fs_z"),
                    se_fs = c("fs_x_se", "fs_m_se", "fs_z_se"),
                    loading_fs = c("x_by_fs_x", "m_by_fs_m", "z_by_fs_z"))
  list(ind = ind,
       se = c(y = ind[1, "fs_y_se"], x = ind[1, "fs_x_se"],
              m = ind[1, "fs_m_se"], z = ind[1, "fs_z_se"],
              xm = ind[1, "fs_x:fs_m_se"], xz = ind[1, "fs_x:fs_z_se"],
              mz = ind[1, "fs_m:fs_z_se"]))
}

test_that("Product-score columns are auto-aliased (no manual rename needed)", {
  setup <- int_setup()
  m <- "y ~ x + m + z + xm + xz + mz"
  fit_new <- tspa(m, data = setup$ind, se_fs = setup$se)
  # the manual-rename workaround must give a bit-identical result
  ind_old <- setup$ind
  ind_old$fs_xm <- ind_old[["fs_x:fs_m"]]
  ind_old$fs_xz <- ind_old[["fs_x:fs_z"]]
  ind_old$fs_mz <- ind_old[["fs_m:fs_z"]]
  fit_old <- tspa(m, data = ind_old, se_fs = setup$se)
  expect_identical(attr(fit_new, "tspaModel"), attr(fit_old, "tspaModel"))
  expect_identical(coef(fit_new), coef(fit_old))
  expect_identical(vcov(fit_new), vcov(fit_old))
  # the generated model names (not the `:`-separated data columns) are used
  # in the rendered model
  m0 <- attr(fit_new, "tspaModel")
  expect_match(m0, "fs_xm~~", fixed = TRUE)
  expect_match(m0, "fs_xz~~", fixed = TRUE)
  expect_match(m0, "fs_mz~~", fixed = TRUE)
  expect_no_match(m0, "fs_x:fs_m", fixed = TRUE)
})

test_that("tspa_sf_alias is a no-op when the score column already exists", {
  setup <- int_setup()
  ind_renamed <- setup$ind
  ind_renamed$fs_xm <- ind_renamed[["fs_x:fs_m"]]
  out <- tspa_sf_alias(ind_renamed,
                       data.frame(xm = 0.1, x = 0.2, m = 0.3))
  expect_identical(out$data, ind_renamed)
  expect_length(out$aliases, 0)
})

