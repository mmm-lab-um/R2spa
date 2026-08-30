library(lavaan)

test_that("compute_fscore() gives same scores as lavaan::lavPredict()", {
  fit <- cfa(" ind60 =~ x1 + x2 + x3 ",
             data = PoliticalDemocracy, std.lv = TRUE)
  fs_lavaan1 <- lavPredict(fit, method = "regression")
  fs_lavaan2 <- lavPredict(fit, method = "Bartlett")
  est <- lavInspect(fit, what = "est")
  fs_hand1 <- compute_fscore(lavInspect(fit, what = "data"),
                             lambda = est$lambda,
                             theta = est$theta,
                             psi = est$psi)
  fs_hand2 <- compute_fscore(lavInspect(fit, what = "data"),
                             lambda = est$lambda,
                             theta = est$theta,
                             psi = est$psi,
                             method = "Bartlett")
  expect_equal(unclass(fs_lavaan1), unclass(fs_hand1))
  expect_equal(unclass(fs_lavaan2), unclass(fs_hand2))
})

test_that("compute_fscore() works for multiple factors", {
  fit <- cfa(
    " ind60 =~ x1 + x2 + x3
      dem60 =~ y1 + y2 + y3 + y4 ",
    data = PoliticalDemocracy)
  fs_lavaan1 <- lavPredict(fit, method = "regression",
                           acov = "standard")
  est <- lavInspect(fit, what = "est")
  fs_hand1 <- compute_fscore(lavInspect(fit, what = "data"),
                             lambda = est$lambda,
                             theta = est$theta,
                             psi = est$psi,
                             acov = TRUE)
  expect_equal(fs_lavaan1, fs_hand1, ignore_attr = TRUE)
  expect_equal(attr(fs_lavaan1, "acov")[[1]],
               attr(fs_hand1, "acov"))
})

test_that("ACOV = Var(e.fs) for Bartlett scores", {
  fit <- cfa(
    " ind60 =~ x1 + x2 + x3
      dem60 =~ y1 + y2 + y3 + y4
      dem65 =~ y5 + y6 + y7 + y8 ",
    data = PoliticalDemocracy)
  fs_lavaan2 <- lavPredict(fit, method = "Bartlett",
                           acov = "standard")
  est <- lavInspect(fit, what = "est")
  fs_hand2 <- compute_fscore(lavInspect(fit, what = "data"),
                             lambda = est$lambda,
                             theta = est$theta,
                             psi = est$psi,
                             acov = TRUE,
                             method = "Bartlett",
                             fs_matrices = TRUE)
  expect_equal(fs_lavaan2, fs_hand2, ignore_attr = TRUE)
  expect_equal(attr(fs_hand2, "acov"), attr(fs_lavaan2, "acov")[[1]])
  expect_equal(attr(fs_hand2, "fsT"), attr(fs_lavaan2, "acov")[[1]],
               ignore_attr = TRUE)
  # From Issue 56
  expect_equal(rownames(attr(fs_hand2, "fsT")),
               paste0("fs_", rownames(attr(fs_lavaan2, "acov")[[1]])))
})

test_that("Same factor scores with same priors", {
  hs_model <- "
  visual  =~ x1 + x2 + x3
  "
  fit <- cfa(hs_model,
             data = HolzingerSwineford1939,
             group = "school",
             group.equal = c("loadings", "intercepts", "residuals"))
  two_cases <-
    c(which(subset(HolzingerSwineford1939, subset = school == "Pasteur",
                   select = id) == 47),
      which(subset(HolzingerSwineford1939, subset = school == "Grant-White",
                   select = id) == 275))
  fs_lavaan <- lavPredict(fit, method = "regression",
                          acov = "standard")
  # Different factor scores in lavaan
  expect_false(fs_lavaan[[1]][two_cases[1], ] ==
                 fs_lavaan[[2]][two_cases[2], ])
  est <- lavInspect(fit, what = "est")
  y <- lavInspect(fit, what = "data")
  fs1_hand <- compute_fscore(y[[1]],
                             lambda = est[[1]]$lambda,
                             theta = est[[1]]$theta,
                             psi = est[[1]]$psi,
                             nu = est[[1]]$nu,
                             alpha = est[[1]]$alpha,
                             acov = TRUE)
  fs2_hand <- compute_fscore(y[[2]],
                             lambda = est[[2]]$lambda,
                             theta = est[[2]]$theta,
                             psi = est[[1]]$psi,
                             nu = est[[2]]$nu,
                             alpha = est[[1]]$alpha,
                             acov = TRUE)
  expect_equal(fs1_hand[two_cases[1], ], fs2_hand[two_cases[2], ])
})

# Moved from `test-get_fscore.R`
# Prepare for test objects
fscore_model <- " ind60 =~ x1 + x2 + x3
                    dem60 =~ y1 + y2 + y3 + y4 "
fit <- cfa(fscore_model, data = PoliticalDemocracy)
fs_lavaan <- lavPredict(fit, method = "regression")
est <- lavInspect(fit, what = "est")
fscore_data <- lavInspect(fit, what = "data")
test_object_fscore <- compute_fscore(fscore_data, lambda = est$lambda,
                                     theta = est$theta, psi = est$psi)

########## Testing section ############

test_that("Test the length of output is equal", {
  expect_equal(nrow(test_object_fscore), nrow(fs_lavaan))
})

test_that("Test the output is the same for fscore and lavaan funtion", {
  expect_equal(test_object_fscore, fs_lavaan, ignore_attr = TRUE)
})

# Check scoring matrices
hs_model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
  "
fit <- cfa(hs_model,
           data = HolzingerSwineford1939,
           group = "school")
est <- lavInspect(fit, what = "est")
y <- lavInspect(fit, what = "data")

test_that("Correct scoring matrix for regression scores", {
  fs1_hand <- compute_fscore(y[[1]],
                             lambda = est[[1]]$lambda,
                             theta = est[[1]]$theta,
                             psi = est[[1]]$psi,
                             nu = est[[1]]$nu,
                             alpha = est[[1]]$alpha,
                             fs_matrices = TRUE)
  a_mat1 <- est[[1]]$psi %*%
    crossprod(est[[1]]$lambda, MASS::ginv(fit@implied$cov[[1]]))
  fsL1 <- attr(fs1_hand, "fsL")
  expect_equal(a_mat1 %*% est[[1]]$lambda, fsL1, ignore_attr = TRUE)
  # Issue 56
  expect_equal(colnames(fsL1), rownames(a_mat1))
  expect_equal(rownames(fsL1), paste0("fs_", colnames(fsL1)))
  expect_equal(a_mat1 %*% cov(y[[1]]) %*% t(a_mat1),
               expected = cov(fs1_hand))
  implied_covfs1 <- a_mat1 %*% fit@implied$cov[[1]] %*% t(a_mat1)
  expect_equal(attr(fs1_hand, "fsT"),
               expected = implied_covfs1 - fsL1 %*% est[[1]]$psi %*% t(fsL1),
               ignore_attr = TRUE)
  fs2_hand <- compute_fscore(y[[2]],
                             lambda = est[[2]]$lambda,
                             theta = est[[2]]$theta,
                             psi = est[[1]]$psi,
                             nu = est[[2]]$nu,
                             alpha = est[[1]]$alpha,
                             fs_matrices = TRUE)
  implied_cov2 <- est[[2]]$lambda %*% est[[1]]$psi %*% t(est[[2]]$lambda) +
    est[[2]]$theta
  a_mat2 <- est[[1]]$psi %*%
    crossprod(est[[2]]$lambda, MASS::ginv(implied_cov2))
  fsL2 <- attr(fs2_hand, "fsL")
  expect_equal(a_mat2 %*% est[[2]]$lambda, fsL2, ignore_attr = TRUE)
  expect_equal(a_mat2 %*% cov(y[[2]]) %*% t(a_mat2),
               expected = cov(fs2_hand))
  implied_covfs2 <- a_mat2 %*% fit@implied$cov[[2]] %*% t(a_mat2)
  expect_equal(attr(fs2_hand, "fsT"),
               expected = implied_covfs2 - fsL2 %*% est[[2]]$psi %*% t(fsL2),
               ignore_attr = TRUE)
})

test_that("Correct scoring matrix for Bartlett scores", {
  fs1_hand <- compute_fscore(y[[1]],
                             lambda = est[[1]]$lambda,
                             theta = est[[1]]$theta,
                             psi = est[[1]]$psi,
                             nu = est[[1]]$nu,
                             alpha = est[[1]]$alpha,
                             method = "Bartlett",
                             fs_matrices = TRUE)
  tlam_invth1 <- crossprod(est[[1]]$lambda,
                           MASS::ginv(est[[1]]$theta))
  a_mat1 <- solve(tlam_invth1 %*% est[[1]]$lambda, tlam_invth1)
  fsL1 <- attr(fs1_hand, "fsL")
  expect_equal(fsL1, diag(3), ignore_attr = TRUE)
  expect_equal(a_mat1 %*% cov(y[[1]]) %*% t(a_mat1),
               expected = cov(fs1_hand),
               ignore_attr = TRUE)
  implied_covfs1 <- a_mat1 %*% fit@implied$cov[[1]] %*% t(a_mat1)
  expect_equal(attr(fs1_hand, "fsT"),
               expected = implied_covfs1 - fsL1 %*% est[[1]]$psi %*% t(fsL1),
               ignore_attr = TRUE)
  fs2_hand <- compute_fscore(y[[2]],
                             lambda = est[[2]]$lambda,
                             theta = est[[2]]$theta,
                             psi = est[[1]]$psi,
                             nu = est[[2]]$nu,
                             alpha = est[[1]]$alpha,
                             method = "Bartlett",
                             fs_matrices = TRUE)
  implied_cov2 <- est[[2]]$lambda %*% est[[1]]$psi %*% t(est[[2]]$lambda) +
    est[[2]]$theta
  tlam_invth2 <- crossprod(est[[2]]$lambda,
                           MASS::ginv(est[[2]]$theta))
  a_mat2 <- solve(tlam_invth2 %*% est[[2]]$lambda, tlam_invth2)
  fsL2 <- attr(fs2_hand, "fsL")
  expect_equal(fsL2, diag(3), ignore_attr = TRUE)
  expect_equal(a_mat2 %*% cov(y[[2]]) %*% t(a_mat2),
               expected = cov(fs2_hand), ignore_attr = TRUE)
  implied_covfs2 <- a_mat2 %*% fit@implied$cov[[2]] %*% t(a_mat2)
  expect_equal(attr(fs2_hand, "fsT"),
               expected = implied_covfs2 - fsL2 %*% est[[2]]$psi %*% t(fsL2),
               ignore_attr = TRUE)
})

test_that("Correction factor shrinks to zero in large sample", {
  cov1 <- lavInspect(fit, what = "implied")[[1]]$cov
  fit_small <- cfa(" visual  =~ x1 + x2 + x3 ",
                   sample.cov = cov1, sample.nobs = 50)
  c1 <- correct_evfs(fit_small, method = "regression")[[1]]
  fit_medium <- cfa(" visual  =~ x1 + x2 + x3 ",
                     sample.cov = cov1, sample.nobs = 60)
  c2 <- correct_evfs(fit_medium, method = "regression")[[1]]
  fit_large <- cfa(" visual  =~ x1 + x2 + x3 ",
                   sample.cov = cov1, sample.nobs = 1e6)
  c3 <- correct_evfs(fit_large, method = "regression")[[1]]
  expect_lt(c2, c1)
  expect_lt(c3, 1e-4)
  c1b <- correct_evfs(fit_small, method = "Bartlett")[[1]]
  expect_gt(c1b, c1)
})

########## method = "mean" (sum scores) ############

# The default cfa() (in this lavaan version) has no mean structure, so
# est$nu / est$alpha are NULL there and the mean-structure variants below
# refit with meanstructure = TRUE.
ms_fit_1f <- cfa("ind60 =~ x1 + x2 + x3",
                 data = PoliticalDemocracy, meanstructure = TRUE)
ms_est_1f <- lavInspect(ms_fit_1f, what = "est")
ms_y_1f <- as.matrix(lavInspect(ms_fit_1f, what = "data"))
ms_fit_1f_nm <- cfa("ind60 =~ x1 + x2 + x3", data = PoliticalDemocracy)
ms_est_1f_nm <- lavInspect(ms_fit_1f_nm, what = "est")
ms_y_1f_nm <- as.matrix(lavInspect(ms_fit_1f_nm, what = "data"))
ms_fit_2f <- cfa("ind60 =~ x1 + x2 + x3
                  dem60 =~ y1 + y2 + y3 + y4",
                 data = PoliticalDemocracy, meanstructure = TRUE)
ms_est_2f <- lavInspect(ms_fit_2f, what = "est")
ms_y_2f <- as.matrix(lavInspect(ms_fit_2f, what = "data"))

test_that("compute_fscore(method = 'mean') gives raw item means (single factor)", {
  fs_hand <- compute_fscore(ms_y_1f,
                            lambda = ms_est_1f$lambda,
                            theta = ms_est_1f$theta,
                            nu = ms_est_1f$nu,
                            alpha = ms_est_1f$alpha,
                            method = "mean",
                            fs_matrices = TRUE)
  # M: one row of 1/3 on the factor's three items
  M1 <- matrix(1 / 3, nrow = 1, ncol = 3,
               dimnames = list("ind60", c("x1", "x2", "x3")))
  # Raw (uncentered) scores: the plain item means
  expect_equal(unname(as.numeric(fs_hand)), unname(rowMeans(ms_y_1f)),
               tolerance = 1e-12)
  lam1 <- unname(as.matrix(ms_est_1f$lambda))
  th1 <- as.matrix(ms_est_1f$theta)
  nu1 <- as.numeric(unname(ms_est_1f$nu))
  alpha1 <- as.numeric(unname(ms_est_1f$alpha))
  # fsL == M %*% lambda == mean of the factor's loadings
  expect_equal(unname(attr(fs_hand, "fsL")), unname(M1 %*% lam1),
               tolerance = 1e-10)
  expect_equal(unname(as.numeric(attr(fs_hand, "fsL"))), mean(lam1),
               tolerance = 1e-10)
  # fsT == M %*% theta %*% t(M) (sum of the residual block, / 9)
  expect_equal(unname(attr(fs_hand, "fsT")),
               unname(M1 %*% th1 %*% t(M1)),
               tolerance = 1e-10)
  # fsb == M %*% nu: the mean of the item intercepts, i.e. the measurement
  # intercept of the score regressed on the uncentered latent
  # (fsb = E[fs] - fsL %*% alpha; here alpha = 0, so fsb = E[fs])
  expect_equal(unname(attr(fs_hand, "fsb")),
               as.numeric(unname(M1 %*% nu1)),
               tolerance = 1e-10)
  expect_equal(unname(attr(fs_hand, "fsb")), mean(nu1), tolerance = 1e-10)
  expect_equal(rownames(attr(fs_hand, "fsL")), "fs_ind60")
  expect_equal(names(attr(fs_hand, "fsb")), "fs_ind60")
})

test_that(
  "compute_fscore(method = 'mean') with no mean structure: fsb falls back to the score column mean",
  {
    expect_null(ms_est_1f_nm$nu)
    expect_null(ms_est_1f_nm$alpha)
    fs_hand <- compute_fscore(ms_y_1f_nm,
                              lambda = ms_est_1f_nm$lambda,
                              theta = ms_est_1f_nm$theta,
                              method = "mean",
                              fs_matrices = TRUE)
    expect_equal(unname(as.numeric(fs_hand)), unname(rowMeans(ms_y_1f_nm)),
                 tolerance = 1e-12)
    # nu falls back to the sample item means, so fsb = E[fs] is exactly
    # the column mean of the (raw) score
    expect_equal(unname(attr(fs_hand, "fsb")),
                 mean(unname(as.numeric(fs_hand))),
                 tolerance = 1e-10)
    expect_equal(unname(attr(fs_hand, "fsb")),
                 mean(colMeans(ms_y_1f_nm)),
                 tolerance = 1e-10)
  }
)

# Hand scoring matrix for the two-factor model (q x p weights)
ms_M2 <- matrix(
  c(rep(1 / 3, 3), rep(0, 4), rep(0, 3), rep(1 / 4, 4)),
  nrow = 2, ncol = 7, byrow = TRUE
)
dimnames(ms_M2) <- list(c("ind60", "dem60"),
                        c("x1", "x2", "x3", "y1", "y2", "y3", "y4"))

test_that("compute_fscore(method = 'mean') gives raw item means (two factors)", {
  fs_hand <- compute_fscore(ms_y_2f,
                            lambda = ms_est_2f$lambda,
                            theta = ms_est_2f$theta,
                            nu = ms_est_2f$nu,
                            alpha = ms_est_2f$alpha,
                            method = "mean",
                            fs_matrices = TRUE)
  # Math-level score columns are the bare factor names
  expect_equal(colnames(fs_hand), c("ind60", "dem60"))
  # One score column per factor: per-factor row means of its items
  expect_equal(unname(fs_hand[, "ind60"]),
               unname(rowMeans(ms_y_2f[, c("x1", "x2", "x3")])),
               tolerance = 1e-12)
  expect_equal(unname(fs_hand[, "dem60"]),
               unname(rowMeans(ms_y_2f[, c("y1", "y2", "y3", "y4")])),
               tolerance = 1e-12)
  lam2 <- unname(as.matrix(ms_est_2f$lambda))
  th2 <- as.matrix(ms_est_2f$theta)
  # fsL == M %*% lambda (full 2 x 2 check)
  expect_equal(unname(attr(fs_hand, "fsL")), unname(ms_M2 %*% lam2),
               tolerance = 1e-10)
  expect_equal(dimnames(attr(fs_hand, "fsL"))[[1]],
               c("fs_ind60", "fs_dem60"))
  # fsT == M %*% theta %*% t(M); theta is diagonal for a plain CFA and the
  # item sets are disjoint, so the off-diagonal is exactly zero
  expect_equal(unname(attr(fs_hand, "fsT")),
               unname(ms_M2 %*% th2 %*% t(ms_M2)),
               tolerance = 1e-10)
  expect_equal(attr(fs_hand, "fsT")[1, 2], 0, tolerance = 0)
  expect_equal(attr(fs_hand, "fsT")[2, 1], 0, tolerance = 0)
  nu2 <- as.numeric(unname(ms_est_2f$nu))
  alpha2 <- as.numeric(unname(ms_est_2f$alpha))
  # fsb == M %*% nu: the mean of each factor's item intercepts, i.e. the
  # measurement intercept of the scores regressed on the uncentered latents
  # (fsb = E[fs] - fsL %*% alpha; here alpha = 0, so fsb = E[fs])
  expect_equal(as.numeric(unname(attr(fs_hand, "fsb"))),
               as.numeric(unname(ms_M2 %*% matrix(nu2, ncol = 1))),
               tolerance = 1e-10)
  expect_equal(names(attr(fs_hand, "fsb")), c("fs_ind60", "fs_dem60"))
})

test_that("compute_fscore(method = 'mean', acov = TRUE) carries fsT as acov", {
  fs_ac <- compute_fscore(ms_y_2f,
                          lambda = ms_est_2f$lambda,
                          theta = ms_est_2f$theta,
                          method = "mean",
                          acov = TRUE,
                          fs_matrices = TRUE)
  # Same values as fsT; acov keeps the bare factor names, fsT the "fs_" names
  expect_equal(attr(fs_ac, "acov"), attr(fs_ac, "fsT"), ignore_attr = TRUE)
})

test_that("compute_a_mean() auto-derives the scoring weights from loadings", {
  lam2 <- as.matrix(lavInspect(ms_fit_2f, what = "est")$lambda)
  M_auto <- R2spa:::compute_a_mean(lam2)
  expect_equal(dim(M_auto), c(2L, 7L))
  expect_setequal(rownames(M_auto), c("ind60", "dem60"))
  expect_setequal(colnames(M_auto),
                  c("x1", "x2", "x3", "y1", "y2", "y3", "y4"))
  expect_equal(unname(M_auto["ind60", c("x1", "x2", "x3")]), rep(1 / 3, 3))
  expect_equal(unname(M_auto["ind60", c("y1", "y2", "y3", "y4")]), rep(0, 4))
  expect_equal(unname(M_auto["dem60", c("y1", "y2", "y3", "y4")]), rep(1 / 4, 4))
  expect_equal(unname(M_auto["dem60", c("x1", "x2", "x3")]), rep(0, 3))
  expect_equal(rowSums(unname(M_auto)), c(1, 1))
  # An explicit sum_items reproducing the loading structure gives the same M
  M_explicit <- R2spa:::compute_a_mean(
    lam2,
    sum_items = list(ind60 = c("x1", "x2", "x3"),
                     dem60 = c("y1", "y2", "y3", "y4"))
  )
  expect_equal(M_explicit, M_auto)
})

test_that("compute_a_mean() auto-derivation errors", {
  # Indicators loading on more than one factor
  lam_x <- matrix(
    0, nrow = 3, ncol = 2,
    dimnames = list(c("x1", "x2", "x3"), c("f1", "f2"))
  )
  lam_x[1, 1] <- 1.0
  lam_x[2, 1:2] <- c(0.5, 0.3)
  lam_x[3, 1] <- 0.6
  expect_error(R2spa:::compute_a_mean(lam_x),
               "load on more than one factor.*x2")
  lam_x[3, 2] <- 0.4
  expect_error(R2spa:::compute_a_mean(lam_x),
               "load on more than one factor.*x2.*x3")
  # A factor with no indicators
  lam_z <- matrix(
    0, nrow = 3, ncol = 2,
    dimnames = list(c("x1", "x2", "x3"), c("f1", "f2"))
  )
  lam_z[1:2, 1] <- c(0.8, 0.9)
  expect_error(R2spa:::compute_a_mean(lam_z),
               "Factor 'f2' has no items")
})

test_that("compute_a_mean() sum_items validation errors", {
  lam2 <- as.matrix(lavInspect(ms_fit_2f, what = "est")$lambda)
  # Unknown factor name
  expect_error(
    R2spa:::compute_a_mean(lam2,
                           sum_items = list(ind60 = c("x1", "x2", "x3"),
                                            nope = c("y1"))),
    "Unknown factor name.*nope"
  )
  # Not all model factors covered
  expect_error(
    R2spa:::compute_a_mean(lam2,
                           sum_items = list(ind60 = c("x1", "x2", "x3"))),
    "must cover all model factors.*dem60"
  )
  # Unknown item name
  expect_error(
    R2spa:::compute_a_mean(lam2,
                           sum_items = list(ind60 = c("x1", "zzz"),
                                            dem60 = c("y1", "y2", "y3", "y4"))),
    "Unknown item name.*zzz"
  )
  # One item assigned to two sums
  expect_error(
    R2spa:::compute_a_mean(lam2,
                           sum_items = list(ind60 = c("x1", "x2", "x3", "y2"),
                                            dem60 = c("y1", "y2"))),
    "assigned to more than one sum.*y2"
  )
  # Unnamed list
  expect_error(
    R2spa:::compute_a_mean(
      lam2,
      sum_items = list(c("x1", "x2", "x3"), c("y1", "y2", "y3", "y4"))
    ),
    "must be a named list mapping factor names to item names"
  )
  # Non-list input
  expect_error(
    R2spa:::compute_a_mean(lam2, sum_items = "ind60"),
    "must be a named list mapping factor names to item names"
  )
  # A factor with zero items
  expect_error(
    R2spa:::compute_a_mean(lam2,
                           sum_items = list(ind60 = character(0),
                                            dem60 = c("y1", "y2", "y3", "y4"))),
    "Factor 'ind60' has no items in 'sum_items'"
  )
})

test_that("compute_fscore(method = 'mean') ignores center_y", {
  fs_raw <- compute_fscore(ms_y_1f,
                           lambda = ms_est_1f$lambda,
                           theta = ms_est_1f$theta,
                           nu = ms_est_1f$nu,
                           alpha = ms_est_1f$alpha,
                           method = "mean",
                           fs_matrices = TRUE)
  fs_cy <- compute_fscore(ms_y_1f,
                          lambda = ms_est_1f$lambda,
                          theta = ms_est_1f$theta,
                          nu = ms_est_1f$nu,
                          alpha = ms_est_1f$alpha,
                          method = "mean",
                          center_y = TRUE,
                          fs_matrices = TRUE)
   expect_identical(fs_raw, fs_cy)
   # The scores are raw (uncentered): their column mean is the item level,
   # not zero
   expect_gt(mean(unname(as.numeric(fs_raw))), 1)
})

test_that(
  "compute_fscore(method = 'mean') fsb is the mean of the item intercepts (M %*% nu), not E[fs]",
  {
    # Hand-built (non-simulated) inputs with a NONZERO latent mean, so the
    # two candidate conventions E[fs] = M(nu + lambda %*% alpha) and the
    # measurement intercept M %*% nu differ by exactly fsL %*% alpha.
    # This avoids the CFA mean-identifiability issue (a plain CFA's latent
    # mean is not identified and comes back 0, masking the distinction).
    set.seed(1)
    n <- 200
    lam <- matrix(c(1.0, 2.0, 1.8), ncol = 1,
                  dimnames = list(c("x1", "x2", "x3"), "f1"))
    nuc <- c(0.5, -1.0, 2.0)
    alphac <- 3.0
    thc <- diag(c(0.4, 0.6, 0.5))
    mus <- as.numeric(lam) * alphac + nuc
    y <- sweep(matrix(rnorm(n * 3), ncol = 3), 2, mus, "+")
    colnames(y) <- c("x1", "x2", "x3")

    fs_hand <- compute_fscore(y, lambda = lam, theta = thc,
                              nu = nuc, alpha = alphac,
                              method = "mean", fs_matrices = TRUE)
    M1 <- matrix(1 / 3, nrow = 1, ncol = 3,
                 dimnames = list("f1", c("x1", "x2", "x3")))
    fsL <- unname(attr(fs_hand, "fsL"))
    fsb <- unname(attr(fs_hand, "fsb"))
    # Scores are the raw item means (unaffected by nu/alpha)
    expect_equal(unname(as.numeric(fs_hand)), unname(rowMeans(y)),
                 tolerance = 1e-12)
    # fsL == M %*% lambda (mean of the loadings)
    expect_equal(as.numeric(fsL), as.numeric(M1 %*% lam), tolerance = 1e-12)
    # THE CONVENTION: fsb is the mean of the item intercepts, M %*% nu
    expect_equal(fsb, mean(nuc), tolerance = 1e-12)
    expect_equal(fsb, as.numeric(M1 %*% nuc), tolerance = 1e-12)
    # ... NOT E[fs] = M %*% (nu + lambda %*% alpha); they differ by fsL*alpha
    efs <- as.numeric(M1 %*% (nuc + as.numeric(lam) * alphac))
    expect_equal(fsb, efs - as.numeric(fsL) * alphac,
                 tolerance = 1e-12)
    expect_false(isTRUE(all.equal(fsb, efs, tolerance = 1e-8)))
  }
)
