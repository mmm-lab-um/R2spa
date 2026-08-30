# Product-indicator error covariances (shared factor scores)
#
# Provenance: covers the stage-2 fix of the measurement-error covariances
# between DMC product factor-score indicators (the `fs_a:fs_b` columns from
# compute_fs_prod()/get_fs(product = )). When two or more product latents
# share a factor score (e.g. `x:m` and `x:z` share `x`), their measurement
# errors are correlated; tspa() computes the covariance of every unordered
# pair from the stage-1 fsL/fsT/psi attributes (fs_prod_ecov(),
# tspa_prod_ecov()) and emits fixed `fs_v1 ~~ <val> * fs_v2` rows into the
# rendered model -- on both the single-factor (se_fs) and multi-factor
# (fsT/fsL) paths, and both the product = TRUE auto-compute flow and the
# manual flow (product columns pre-computed and listed in se_fs; detected
# via tspa_sf_alias()'s prod_map). Pairs whose value is ~0 are dropped
# (a single product -> no rows); data lacking the stage-1 attributes with
# two or more product latents is rejected with an informative error.
library(lavaan)

## Shared unit matrices (hand-picked, fixed) ----------------------------------
## Rows of L are factor scores, columns are latent variables; T is the
## score-error (co)variance, P the latent (co)variance (the 'psi' attribute).
## All are positive definite (Gershgorin), so they also serve as the
## joint-normal DGP of the MC moment test below.
ecov_L3 <- matrix(c(1, 0.2, 0,
                    0, 1, 0.1,
                    0.1, 0, 1), nrow = 3, byrow = TRUE)
rownames(ecov_L3) <- c("fs_x", "fs_m", "fs_z")
colnames(ecov_L3) <- c("x", "m", "z")
ecov_T3 <- matrix(c(0.20, 0.05, 0,
                    0.05, 0.30, 0.02,
                    0, 0.02, 0.40), nrow = 3, byrow = TRUE)
ecov_P3 <- matrix(c(1, 0.3, 0.2,
                    0.3, 1, 0.4,
                    0.2, 0.4, 1), nrow = 3)
ecov_L4 <- matrix(c(1, 0.2, 0, 0.1,
                    0, 1, 0.1, 0,
                    0.1, 0, 1, 0.2,
                    0, 0.2, 0, 1), nrow = 4, byrow = TRUE)
rownames(ecov_L4) <- paste0("fs_", c("x", "m", "z", "w"))
colnames(ecov_L4) <- c("x", "m", "z", "w")
ecov_T4 <- matrix(c(0.20, 0.05, 0, 0.03,
                    0.05, 0.30, 0.02, 0,
                    0, 0.02, 0.40, 0.04,
                    0.03, 0, 0.04, 0.50), nrow = 4, byrow = TRUE)
ecov_P4 <- matrix(c(1, 0.3, 0.2, 0.1,
                    0.3, 1, 0.4, 0.2,
                    0.2, 0.4, 1, 0.3,
                    0.1, 0.2, 0.3, 1), nrow = 4)
ecov_T3d <- diag(c(0.1, 0.2, 0.3))
ecov_T4d <- diag(c(0.2, 0.3, 0.4, 0.5))

## The rendered error rows with lhs != rhs (the cross-covariance rows, in
## contrast to the diagonal `fs_v ~~ <val> * fs_v` rows).
ecov_cross_lines <- function(m0) {
  lines <- strsplit(m0, "\n")[[1L]]
  err <- lines[grepl("^fs_[A-Za-z0-9_.]+ ~~ [^ ]+ \\* fs_[A-Za-z0-9_.]+$",
                     lines)]
  keep <- vapply(
    err, function(ln) sub(" ~~ .*$", "", ln) != sub("^.* \\* ", "", ln),
    logical(1L)
  )
  err[keep]
}

# ---------------------------------------------------------------------------
# Group 1: fs_prod_ecov() units
# ---------------------------------------------------------------------------

test_that("the diagonal (i, j) = (k, l) reproduces fs_prod_se2() exactly", {
  # (a) q = 3, T diagonal, L identity-ish, full psi
  La <- diag(c(1, 0.9, 0.8))
  for (ij in list(c(1L, 2L), c(1L, 3L), c(2L, 3L))) {
    expect_equal(R2spa:::fs_prod_ecov(La, ecov_T3d, ecov_P3,
                                      ij[1L], ij[2L], ij[1L], ij[2L]),
                 R2spa:::fs_prod_se2(La, ecov_T3d, ecov_P3, ij[1L], ij[2L]),
                 tolerance = 1e-12,
                 info = paste0("diagonal T, ij = ", paste(ij, collapse = ",")))
  }
  # (b) q = 3, full L and full (non-diagonal) T
  for (ij in list(c(1L, 2L), c(1L, 3L), c(2L, 3L))) {
    expect_equal(R2spa:::fs_prod_ecov(ecov_L3, ecov_T3, ecov_P3,
                                      ij[1L], ij[2L], ij[1L], ij[2L]),
                 R2spa:::fs_prod_se2(ecov_L3, ecov_T3, ecov_P3,
                                     ij[1L], ij[2L]),
                 tolerance = 1e-12,
                 info = paste0("full 3 x 3, ij = ", paste(ij, collapse = ",")))
  }
  # (c) q = 4, full L, T, and psi
  for (ij in list(c(1L, 2L), c(1L, 4L), c(3L, 4L))) {
    expect_equal(R2spa:::fs_prod_ecov(ecov_L4, ecov_T4, ecov_P4,
                                      ij[1L], ij[2L], ij[1L], ij[2L]),
                 R2spa:::fs_prod_se2(ecov_L4, ecov_T4, ecov_P4,
                                     ij[1L], ij[2L]),
                 tolerance = 1e-12,
                 info = paste0("full 4 x 4, ij = ", paste(ij, collapse = ",")))
  }
})

test_that("symmetric in the two product pairs (i, j) <-> (k, l)", {
  expect_equal(R2spa:::fs_prod_ecov(ecov_L3, ecov_T3, ecov_P3, 1L, 2L, 1L, 3L),
               R2spa:::fs_prod_ecov(ecov_L3, ecov_T3, ecov_P3, 1L, 3L, 1L, 2L),
               tolerance = 1e-12)
  expect_equal(R2spa:::fs_prod_ecov(ecov_L4, ecov_T4, ecov_P4, 1L, 2L, 3L, 4L),
               R2spa:::fs_prod_ecov(ecov_L4, ecov_T4, ecov_P4, 3L, 4L, 1L, 2L),
               tolerance = 1e-12)
})

test_that("diagonal T: a shared-score pair reduces to tau_other * s_i^2", {
  # The pairs (i, j) and (i, k) share the score i; with T diagonal every
  # term of the 6-term formula that carries an off-diagonal T entry is
  # zero, leaving ecov = L_j psi L_k' * T[i, i] -- the other pair's tau
  # times the shared score's error variance.
  for (ijk in list(c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
                   c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L))) {
    i <- ijk[1L]; j <- ijk[2L]; k <- ijk[3L]
    rhs <- as.numeric(ecov_L3[j, , drop = FALSE] %*% ecov_P3 %*%
                      t(ecov_L3[k, , drop = FALSE])) * ecov_T3d[i, i]
    expect_equal(R2spa:::fs_prod_ecov(ecov_L3, ecov_T3d, ecov_P3, i, j, i, k),
                 rhs, tolerance = 1e-12,
                 info = paste0("ijk = ", paste(ijk, collapse = ",")))
  }
})

test_that("diagonal T: disjoint pairs (four distinct scores) give exactly 0", {
  for (ijkl in list(c(1L, 2L, 3L, 4L), c(2L, 1L, 4L, 3L),
                    c(1L, 4L, 2L, 3L))) {
    expect_equal(
      R2spa:::fs_prod_ecov(ecov_L4, ecov_T4d, ecov_P4,
                           ijkl[1L], ijkl[2L], ijkl[3L], ijkl[4L]),
      0,
      info = paste0("ijkl = ", paste(ijkl, collapse = ","))
    )
  }
})

test_that("a hand-expanded 6-term value is pinned exactly", {
  L5 <- matrix(c(1, 0, 0,
                 1, 1, 0,
                 0, 0, 2), nrow = 3, byrow = TRUE)
  T5 <- diag(c(0.1, 0.2, 0.3))
  P5 <- matrix(c(1, 0.5, 0,
                 0.5, 1, 0.25,
                 0, 0.25, 1), nrow = 3)
  # Pairs (1, 2) and (1, 3); T5 diagonal, so the 6 terms
  #   tau_ik c_jl + tau_il c_jk + tau_jk c_il + tau_jl c_ik
  #   + c_ik c_jl + c_il c_jk
  # reduce to the single nonzero product tau_jl * c_ik:
  #   tau_11 * T5[2, 3] = 1   * 0   = 0
  #   tau_13 * T5[2, 1] = 0   * 0   = 0
  #   tau_21 * T5[1, 3] = 1.5 * 0   = 0
  #   tau_23 * T5[1, 1] = 0.5 * 0.1 = 0.05
  #       (tau_23 = (1, 1, 0) P5 (0, 0, 2)' = 2 * P5[2, 3] = 0.5)
  #   T5[1, 1] * T5[2, 3] = 0.1 * 0 = 0
  #   T5[1, 3] * T5[2, 1] = 0   * 0   = 0
  expect_equal(R2spa:::fs_prod_ecov(L5, T5, P5, 1L, 2L, 1L, 3L), 0.05,
               tolerance = 1e-14)
})

test_that("the population moment matches the formula (MC validation)", {
  # (xi, e) independent joint-normal with zero means -- the
  # scoring-orthogonality condition Cov(e, xi) = 0. For the DMC products
  # P_12, P_13 of the scores of the pairs (1, 2) and (1, 3),
  #   cov(P_12, P_13) - cov(Q_12, Q_13) -> Cov(u_12, u_13)
  # where Q_ij = (L_i xi)(L_j xi) - tau_ij is the population-centered
  # (alpha = 0) "signal" part and u_ij the measurement error, so the
  # difference isolates fs_prod_ecov(L, T, psi, 1, 2, 1, 3) including the
  # double-mean-centering constants. Observed with seed 2118, n = 200000:
  # |deviation| = 0.00315; SE of the moment difference
  # sd(p12 * p13 - q12 * q13) / sqrt(n) = 0.00714 (the proxy
  # sd((p12 - q12) * (p13 - q13)) / sqrt(n) of the difference series is
  # 0.00263); tolerance 0.03 ~= 4.2 SE.
  set.seed(2118)
  n <- 200000L
  xi <- MASS::mvrnorm(n, mu = rep(0, 3), Sigma = ecov_P3)
  e <- MASS::mvrnorm(n, mu = rep(0, 3), Sigma = ecov_T3)
  fs <- xi %*% t(ecov_L3) + e
  f1 <- fs[, 1L]; f2 <- fs[, 2L]; f3 <- fs[, 3L]
  p12 <- (f1 - mean(f1)) * (f2 - mean(f2))
  p12 <- p12 - mean(p12)
  p13 <- (f1 - mean(f1)) * (f3 - mean(f3))
  p13 <- p13 - mean(p13)
  x1 <- as.numeric(xi %*% ecov_L3[1L, ])
  x2 <- as.numeric(xi %*% ecov_L3[2L, ])
  x3 <- as.numeric(xi %*% ecov_L3[3L, ])
  q12 <- x1 * x2 - as.numeric(ecov_L3[1L, , drop = FALSE] %*% ecov_P3 %*%
                              t(ecov_L3[2L, , drop = FALSE]))
  q13 <- x1 * x3 - as.numeric(ecov_L3[1L, , drop = FALSE] %*% ecov_P3 %*%
                              t(ecov_L3[3L, , drop = FALSE]))
  expect_equal(cov(p12, p13) - cov(q12, q13),
               R2spa:::fs_prod_ecov(ecov_L3, ecov_T3, ecov_P3, 1L, 2L, 1L,
                                    3L),
               tolerance = 0.03)
})

# ---------------------------------------------------------------------------
# Group 2: tspa_prod_ecov() + schema rendering
# ---------------------------------------------------------------------------

test_that("tspa_prod_ecov(): shared-score pairs give tau_other * s_i^2, and zero pairs are dropped", {
  prods3 <- data.frame(v = c("xm", "xz", "mz"), a = c("x", "x", "m"),
                       b = c("m", "z", "z"), stringsAsFactors = FALSE)
  Tp <- diag(c(0.2, 0.3, 0.4))
  pe <- R2spa:::tspa_prod_ecov(prods3, ecov_L3, Tp, ecov_P3)
  # (r, p) loop order over the prods rows
  expect_equal(pe$v1, c("xm", "xm", "xz"))
  expect_equal(pe$v2, c("xz", "mz", "mz"))
  # Hand-computed with diagonal T (each = the other pair's tau * s_i^2):
  #   tau_mz = (0, 1, 0.1) P3 (0.1, 0, 1)'   = 0.532 -> xm:xz 0.532 * 0.2
  #   tau_xz = (1, 0.2, 0) P3 (0.1, 0, 1)'   = 0.386 -> xm:mz 0.386 * 0.3
  #   tau_xm = (1, 0.2, 0) P3 (0, 1, 0.1)'   = 0.528 -> xz:mz 0.528 * 0.4
  expect_equal(pe$ecov, c(0.1064, 0.1158, 0.2112), tolerance = 1e-12)
  # The same values from the shared-factor reduction, on the same matrices
  expect_equal(
    pe$ecov[1L],
    as.numeric(ecov_L3[2L, , drop = FALSE] %*% ecov_P3 %*%
               t(ecov_L3[3L, , drop = FALSE])) * Tp[1L, 1L],
    tolerance = 1e-12
  )
  expect_equal(
    pe$ecov[2L],
    as.numeric(ecov_L3[1L, , drop = FALSE] %*% ecov_P3 %*%
               t(ecov_L3[3L, , drop = FALSE])) * Tp[2L, 2L],
    tolerance = 1e-12
  )
  expect_equal(
    pe$ecov[3L],
    as.numeric(ecov_L3[1L, , drop = FALSE] %*% ecov_P3 %*%
               t(ecov_L3[2L, , drop = FALSE])) * Tp[3L, 3L],
    tolerance = 1e-12
  )
  # a single product: nothing to pair -> NULL
  expect_null(R2spa:::tspa_prod_ecov(prods3[1L, , drop = FALSE],
                                     ecov_L3, Tp, ecov_P3))
  # two products over four latents sharing no score: all zero -> NULL
  prods2 <- data.frame(v = c("xm", "zw"), a = c("x", "z"),
                       b = c("m", "w"), stringsAsFactors = FALSE)
  L4d <- diag(4)
  rownames(L4d) <- paste0("fs_", c("x", "m", "z", "w"))
  colnames(L4d) <- c("x", "m", "z", "w")
  expect_null(R2spa:::tspa_prod_ecov(prods2, L4d, ecov_T4d, ecov_P4))
})

test_that("sf schema: the fixed cross rows render, and the no-ecov call is unchanged", {
  m8 <- "y ~ x + m + z + xm + xz + mz"
  se6 <- data.frame(x = 0.1, m = 0.2, z = 0.3, xm = 0.5, xz = 0.6, mz = 0.7)
  pe3 <- data.frame(v1 = c("xm", "xm", "xz"), v2 = c("xz", "mz", "mz"),
                    ecov = c(0.06, 0.09, 0.07), stringsAsFactors = FALSE)
  exp3 <- c("fs_xm ~~ 0.06 * fs_xz", "fs_xm ~~ 0.09 * fs_mz",
            "fs_xz ~~ 0.07 * fs_mz")
  got <- R2spa:::tspa_sf(m8, data.frame(x = 1), se6, pe3)
  got0 <- R2spa:::tspa_sf(m8, data.frame(x = 1), se6)
  expect_setequal(ecov_cross_lines(got), exp3)
  # the pre-feature call form (prod_ecov = NULL) has no cross rows, and
  # the ecov render is exactly the plain render with the three cross
  # lines appended at the end of the error block (set-difference would
  # collapse the duplicated blank section lines)
  expect_length(ecov_cross_lines(got0), 0L)
  lg0 <- strsplit(got0, "\n")[[1L]]
  i0 <- which(grepl("^# structural model$", lg0))
  expect_equal(strsplit(got, "\n")[[1L]],
               append(lg0, exp3, after = i0 - 2L))
})

test_that("mf schema: the fixed cross rows render alongside the product rows", {
  m8 <- "y ~ x + m + z + xm + xz + mz"
  prods_mf <- data.frame(v = c("xm", "xz", "mz"), ld = c(1, 1, 1),
                         se2 = c(0.5, 0.6, 0.7), stringsAsFactors = FALSE)
  pe3 <- data.frame(v1 = c("xm", "xm", "xz"), v2 = c("xz", "mz", "mz"),
                    ecov = c(0.06, 0.09, 0.07), stringsAsFactors = FALSE)
  exp3 <- c("fs_xm ~~ 0.06 * fs_xz", "fs_xm ~~ 0.09 * fs_mz",
            "fs_xz ~~ 0.07 * fs_mz")
  got <- R2spa:::tspa_mf(m8, NULL, ecov_T3d, ecov_L3, NULL, prods_mf, pe3)
  ln <- strsplit(got, "\n")[[1L]]
  expect_true(all(exp3 %in% ln))
  got0 <- R2spa:::tspa_mf(m8, NULL, ecov_T3d, ecov_L3, NULL, prods_mf)
  expect_false(any(exp3 %in% strsplit(got0, "\n")[[1L]]))
  # the product loading / error-variance rows are still there
  expect_match(got, "xm =~ 1 * fs_xm", fixed = TRUE)
  expect_match(got, "fs_xm ~~ 0.5 * fs_xm", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# Group 3: end-to-end (the vignette DGP)
# ---------------------------------------------------------------------------
#
# Four correlated unit-variance latents x, m, z and the outcome latent y,
# driven by the main effects plus the three pairwise products. The residual
# variance 0.58890039 makes the population variance of the y latent exactly
# 1, so the std.lv = TRUE score scale is the simulation scale and the true
# stage-2 coefficients are exactly the simulated ones.

ecov_dgp <- function(n = 50000L) {
  set.seed(2119)
  cov_xmz_ey <- matrix(c(1, 0.1, 0.15, 0,
                         0.1, 1, 0.12, 0,
                         0.15, 0.12, 1, 0,
                         0, 0, 0, 0.58890039), nrow = 4)
  eta <- as.data.frame(MASS::mvrnorm(n, mu = rep(0, 4),
                                     Sigma = cov_xmz_ey))
  names(eta) <- c("x", "m", "z", "ey")
  lk <- list(x = c(0.9, 0.8, 0.7), m = c(0.85, 0.75, 0.65),
             z = c(0.8, 0.7, 0.6), y = c(0.75, 0.7, 0.65))
  obs <- setNames(lapply(c("x", "m", "z"), function(v0) {
    eta[[v0]] %*% t(lk[[v0]]) + rnorm(n * 3)
  }), c("x", "m", "z"))
  etay <- 0.3 * eta$x + 0.4 * eta$m + 0.2 * eta$z +
    0.1 * eta$x * eta$m + 0.15 * eta$x * eta$z +
    0.12 * eta$m * eta$z + eta$ey
  obs$y <- etay %*% t(lk$y) + rnorm(n * 3)
  df <- as.data.frame(do.call(cbind, obs[c("x", "m", "z", "y")]))
  names(df) <- c(paste0("x", 1:3), paste0("m", 1:3), paste0("z", 1:3),
                 paste0("y", 1:3))
  df
}

ecov_stage1_model <- "x =~ x1 + x2 + x3
                m =~ m1 + m2 + m3
                z =~ z1 + z2 + z3
                y =~ y1 + y2 + y3"

ecov_prods3 <- data.frame(v = c("xm", "xz", "mz"), a = c("x", "x", "m"),
                          b = c("m", "z", "z"), stringsAsFactors = FALSE)

# Cached end-to-end fixtures: the DGP + stage-1 fits + stage-2 fits are
# built at most once per session (the n = 50000 CFA and sem fits dominate
# the runtime of this file).
ecov_fx_env <- new.env()
ecov_fx <- function() {
  fx <- ecov_fx_env$fx
  if (is.null(fx)) {
    df <- ecov_dgp()
    fs <- get_fs(df, model = ecov_stage1_model, std.lv = TRUE,
                 method = "Bartlett")
    fs_prod <- get_fs(df, model = ecov_stage1_model, std.lv = TRUE,
                      method = "Bartlett", product = "x:m + x:z + m:z")
    se_reg <- c(y = fs[1L, "fs_y_se"], x = fs[1L, "fs_x_se"],
                m = fs[1L, "fs_m_se"], z = fs[1L, "fs_z_se"])
    se_man <- c(se_reg, xm = fs_prod[1L, "fs_x:fs_m_se"],
                xz = fs_prod[1L, "fs_x:fs_z_se"],
                mz = fs_prod[1L, "fs_m:fs_z_se"])
    fx <- list(
      df = df, fs = fs, fs_prod = fs_prod, se_reg = se_reg, se_man = se_man,
      fit_auto = suppressWarnings(
        tspa("y ~ x + m + z + x:m + x:z + m:z", data = fs,
             se_fs = se_reg, product = TRUE)
      ),
      fit_man = suppressWarnings(
        tspa("y ~ x + m + z + xm + xz + mz", data = fs_prod, se_fs = se_man)
      ),
      fit_one = suppressWarnings(
        tspa("y ~ x + m + z + x:m", data = fs, se_fs = se_reg,
             product = TRUE)
      )
    )
    ecov_fx_env$fx <- fx
  }
  fx
}

# The no-std.lv (marker-scale) multi-factor fixtures.
ecov_fx_ns_env <- new.env()
ecov_fx_ns <- function() {
  fx <- ecov_fx_ns_env$fx
  if (is.null(fx)) {
    df <- ecov_fx()$df
    fs_ns <- get_fs(df, model = ecov_stage1_model, method = "Bartlett")
    L1 <- attr(fs_ns, "fsL"); if (is.list(L1)) L1 <- L1[[1L]]
    T1 <- attr(fs_ns, "fsT"); if (is.list(T1)) T1 <- T1[[1L]]
    P1 <- R2spa:::fs_psi_matrix(attr(fs_ns, "psi"))
    fx <- list(
      fs_ns = fs_ns,
      pe = R2spa:::tspa_prod_ecov(ecov_prods3, L1, T1, P1),
      fit_mf = suppressWarnings(
        tspa("y ~ x + m + z + x:m + x:z + m:z", data = fs_ns,
             product = TRUE)
      )
    )
    ecov_fx_ns_env$fx <- fx
  }
  fx
}

test_that("sf auto-compute: the shared-score ecov rows are fixed and the simulated coefficients are recovered", {
  fx <- ecov_fx()
  fit <- fx$fit_auto
  m0 <- attr(fit, "tspaModel")
  # (a) the three cross lines equal fs_prod_ecov() recomputed from the
  #     stage-1 attributes of the very data the fit used (the replay
  #     record carries them)
  fdat <- attr(fit, "tspa_args")$data
  L1 <- attr(fdat, "fsL"); if (is.list(L1)) L1 <- L1[[1L]]
  T1 <- attr(fdat, "fsT"); if (is.list(T1)) T1 <- T1[[1L]]
  P1 <- R2spa:::fs_psi_matrix(attr(fdat, "psi"))
  pe <- R2spa:::tspa_prod_ecov(ecov_prods3, L1, T1, P1)
  expect_equal(nrow(pe), 3L)
  expect_equal(pe$v1, c("xm", "xm", "xz"))
  expect_equal(pe$v2, c("xz", "mz", "mz"))
  for (k in seq_len(3L)) {
    expect_match(
      m0,
      paste0("fs_", pe$v1[k], " ~~ ", as.character(pe$ecov[k]),
             " * fs_", pe$v2[k]),
      fixed = TRUE, info = paste("ecov line", k)
    )
  }
  # (b) the six structural coefficients recover the simulated values.
  # Observed |diff| at n = 50000 (seed 2119): y~x 0.0104, y~m 0.0048,
  # y~z 0.0059, y~xm 0.0025, y~xz 0.0096, y~mz 0.0105 (max 0.0105, a
  # 3x margin under the 0.03 tolerance; pure MC noise at this n).
  cf <- coef(fit)
  truth <- c("y~x" = 0.3, "y~m" = 0.4, "y~z" = 0.2,
             "y~xm" = 0.1, "y~xz" = 0.15, "y~mz" = 0.12)
  for (nm in names(truth)) {
    expect_true(abs(as.numeric(cf[nm]) - truth[nm]) <= 0.03, info = nm)
  }
})

test_that("sf manual flow: the prod_map detection renders the identical model", {
  fx <- ecov_fx()
  expect_identical(attr(fx$fit_man, "tspaModel"),
                   attr(fx$fit_auto, "tspaModel"))
  expect_equal(coef(fx$fit_man), coef(fx$fit_auto))
  expect_equal(vcov(fx$fit_man), vcov(fx$fit_auto), ignore_attr = TRUE)
})

test_that("mf path (no std.lv): the same three cross rows are fixed at the marker scale", {
  fxn <- ecov_fx_ns()
  fx <- ecov_fx()
  m0 <- attr(fxn$fit_mf, "tspaModel")
  pe <- fxn$pe
  expect_equal(nrow(pe), 3L)
  for (k in seq_len(3L)) {
    expect_match(
      m0,
      paste0("fs_", pe$v1[k], " ~~ ", as.character(pe$ecov[k]),
             " * fs_", pe$v2[k]),
      fixed = TRUE, info = paste("ecov line", k)
    )
  }
  # the marker-scale values differ from the std.lv (score-scale) ones
  L1 <- attr(fx$fs, "fsL"); if (is.list(L1)) L1 <- L1[[1L]]
  T1 <- attr(fx$fs, "fsT"); if (is.list(T1)) T1 <- T1[[1L]]
  P1 <- R2spa:::fs_psi_matrix(attr(fx$fs, "psi"))
  pe_std <- R2spa:::tspa_prod_ecov(ecov_prods3, L1, T1, P1)
  expect_false(isTRUE(all.equal(unname(pe$ecov), unname(pe_std$ecov))))
  expect_true(all(is.finite(as.numeric(coef(fxn$fit_mf)))))
})

test_that("data lacking the stage-1 attributes is rejected with two or more product latents", {
  fx <- ecov_fx()
  # a cbind()ed frame carries only the columns (cbind drops the
  # stage-1 attributes) but keeps the pre-computed product columns
  cb <- cbind(
    fx$fs_prod[, c("fs_y", "fs_y_se", "fs_x", "fs_x_se", "fs_m", "fs_m_se",
                   "fs_z", "fs_z_se", "fs_x:fs_m", "fs_x:fs_m_se",
                   "fs_x:fs_z", "fs_x:fs_z_se")]
  )
  expect_null(attr(cb, "fsL"))
  # manual flow: two products in the model, the product SEs in se_fs
  se_man2 <- c(fx$se_reg, xm = fx$fs_prod[1L, "fs_x:fs_m_se"],
               xz = fx$fs_prod[1L, "fs_x:fs_z_se"])
  expect_error(
    tspa("y ~ x + m + z + xm + xz", data = cb, se_fs = se_man2),
    "stage-1 attributes"
  )
  # product = TRUE: the columns pre-exist so auto-compute is skipped, but
  # the error covariances still need the attributes
  expect_error(
    tspa("y ~ x + m + z + x:m + x:z", data = cb, se_fs = fx$se_reg,
         product = TRUE),
    "stage-1 attributes"
  )
})

test_that("a single product latent emits no error-covariance rows", {
  fx <- ecov_fx()
  m0 <- attr(fx$fit_one, "tspaModel")
  expect_length(ecov_cross_lines(m0), 0L)
  # the product's own (diagonal) error-variance row is still present
  expect_match(m0, "fs_xm ~~ [0-9.]+ \\* fs_xm")
  expect_true(all(is.finite(as.numeric(coef(fx$fit_one)))))
})

test_that("a fit with ecov rows replays identically via tspa_args", {
  # The replay data already carries the aliased product columns, so the
  # product pairs are recovered from tspa_sf_alias()'s prod_map on replay.
  fx <- ecov_fx()
  fit_re <- do.call(tspa, attr(fx$fit_auto, "tspa_args"))
  expect_identical(coef(fit_re), coef(fx$fit_auto))
  expect_identical(vcov(fit_re), vcov(fx$fit_auto))
  expect_identical(attr(fit_re, "tspaModel"), attr(fx$fit_auto, "tspaModel"))
})
