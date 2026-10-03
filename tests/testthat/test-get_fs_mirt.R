# Single-group mirt fixtures; shared per-row identities live in helper-mirt_fs.R.

# ---- fixtures -----------------------------------------------------------
set.seed(2024)
n <- 120
# 1-factor, 5 items
d1 <- data.frame(lapply(1:5, function(k) rbinom(n, 1, sample(0.3:0.7))),
                 check.names = FALSE)
m1 <- suppressWarnings(mirt::mirt(d1, 1))
# 2-factor, 8 items
d2 <- data.frame(lapply(1:8, function(k) rbinom(n, 1, sample(0.3:0.7))),
                 check.names = FALSE)
m2 <- suppressWarnings(mirt::mirt(d2, 2))
# 1-factor with a few completely-missing (all-NA) rows
na_rows <- c(3L, 4L, 5L, 100L)
d_na <- d1
d_na[na_rows, ] <- NA
m_na <- suppressWarnings(mirt::mirt(d_na, 1))

fs <- get_fs(m1)
fs2 <- get_fs(m2)
fs_na <- get_fs(m_na)

test_that("get_fs(): S3 dispatch -- mirt S4 objects route to the mirt methods", {
  expect_true(inherits(m1, "SingleGroupClass"))
  # SingleGroupClass -> a data.frame carrying the per-observation marker
  expect_true(is.data.frame(fs))
  expect_true(isTRUE(attr(fs, "mirt_per_obs")))
  # Both mirt methods are registered (multi-group behaviour is exercised in
  # test-get_fs_mirt_multigroup.R)
  expect_true(exists("get_fs.SingleGroupClass", where = asNamespace("R2spa")))
  expect_true(exists("get_fs.MultipleGroupClass", where = asNamespace("R2spa")))
  # Non-mirt input still routes to get_fs.default (unchanged behaviour)
  expect_error(get_fs(42L), regexp = "not implemented for objects of class")
})

test_that("get_fs(): per-row identities (1-factor + 2-factor)", {
  assert_mirt_row_identities(fs, m1, n_rows = nrow(d1))
  assert_mirt_row_identities(fs2, m2, n_rows = nrow(d2))
})

test_that("get_fs(): group-level psi/alpha/fsb + fs_pattern (per-row shapes in helper)", {
  q1 <- mirt::extract.mirt(m1, "nfact")
  fn1 <- mirt::extract.mirt(m1, "factorNames")
  n1 <- nrow(fs)
  # fs_pattern: label = 1..n, pat = NULL
  fp <- attr(fs, "fs_pattern")
  expect_equal(fp$label, seq_len(n1))
  expect_null(fp$pat)
  # psi == diag(q), named by the factor names
  psi <- attr(fs, "psi")
  expect_equal(as.matrix(psi), diag(q1), tolerance = 1e-10, ignore_attr = TRUE)
  expect_identical(rownames(psi), fn1)
  expect_identical(colnames(psi), fn1)
  # alpha: a named zero vector
  alpha <- attr(fs, "alpha")
  expect_length(alpha, q1)
  expect_named(alpha, fn1)
  expect_true(all(unname(alpha) == 0))
  # fsb: a per-row list of named q-vectors, all zero when alpha = 0
  fsb <- attr(fs, "fsb")
  expect_true(is.list(fsb))
  expect_length(fsb, n1)
  expect_true(all(vapply(fsb, function(x) {
    length(x) == q1 && identical(names(x), paste0("fs_", fn1)) &&
      all(unname(x) == 0)
  }, logical(1L))))
})

test_that("get_fs(): 2-factor off-diagonals non-identity; column order (engine in helper)", {
  q2 <- mirt::extract.mirt(m2, "nfact")
  fn2 <- mirt::extract.mirt(m2, "factorNames")
  fsn2 <- paste0("fs_", fn2)
  expect_equal(q2, 2L)
  expect_setequal(fn2, c("F1", "F2"))
  # off-diagonal implied loadings / error covariance are not all zero
  expect_true(any(abs(fs2[["F1_by_fs_F2"]]) > 1e-8))
  expect_true(any(abs(fs2[["F2_by_fs_F1"]]) > 1e-8))
  expect_true(any(abs(fs2[["ecov_fs_F2_fs_F1"]]) > 1e-8))
  # implied-loadings columns are column-major per latent (as in test-fs_indiv.R)
  ld_cols <- grep("_by_fs_", names(fs2), value = TRUE)
  expect_identical(
    ld_cols,
    unlist(lapply(seq_len(q2), function(j) paste(fn2[j], fsn2, sep = "_by_")))
  )
  # ev/ecov columns are the lower triangle in i-outer / j<=i-inner order
  ev_cols <- grep("^ev_|^ecov_", names(fs2), value = TRUE)
  exp_ev <- character(as.integer(q2 * (q2 + 1L) / 2L))
  cnt <- 1L
  for (i in seq_len(q2)) {
    for (j in seq_len(i)) {
      exp_ev[cnt] <- if (i == j) {
        paste0("ev_", fsn2[i])
      } else {
        paste0("ecov_", fsn2[i], "_", fsn2[j])
      }
      cnt <- cnt + 1L
    }
  }
  expect_identical(ev_cols, exp_ev)
})

test_that("get_fs(): completely-missing rows -> NA score/SE/ev, rows preserved", {
  expect_equal(nrow(fs_na), nrow(d_na))
  # the all-NA rows are NA in score / SE / ev / implied loading
  expect_true(all(is.na(fs_na[na_rows, "fs_F1"])))
  expect_true(all(is.na(fs_na[na_rows, "fs_F1_se"])))
  expect_true(all(is.na(fs_na[na_rows, "ev_fs_F1"])))
  expect_true(all(is.na(fs_na[na_rows, "F1_by_fs_F1"])))
  # but the scorable rows are not NA
  scorable <- setdiff(seq_len(nrow(d_na)), na_rows)
  expect_true(!any(is.na(fs_na[scorable, "fs_F1"])))
  # the per-row fsT for a missing row is an all-NA block
  for (i in na_rows) {
    expect_true(all(is.na(attr(fs_na, "fsT")[[i]])))
  }
  # fs_indiv() preserves the rows with NA per-row values
  ind_na <- fs_indiv(fs_na)
  expect_equal(nrow(ind_na), nrow(fs_na))
  expect_true(all(is.na(ind_na[na_rows, "fs_F1"])))
  expect_true(all(is.na(ind_na[na_rows, "fs_F1_se"])))
  expect_true(all(is.na(ind_na[na_rows, "ev_fs_F1"])))
})

# ============================================================================
# 9. prior_mean: non-zero factor prior mean -> per-row non-zero fsb
# ============================================================================

test_that("get_fs(): prior_mean -> per-row fsb == Vpost_i * alpha (1-factor)", {
  al <- 2.0
  n1 <- nrow(fs)
  fs0 <- get_fs(m1)                          # default (alpha = 0)
  fsp <- get_fs(m1, prior_mean = c(F1 = al))
  fb0 <- attr(fs0, "fsb")
  fb <- attr(fsp, "fsb")
  # per-row list of named (by score name) q-vectors
  expect_true(is.list(fb))
  expect_length(fb, n1)
  expect_true(all(vapply(fb, function(x) identical(names(x), "fs_F1"),
                        logical(1L))))
  # default (alpha = 0): per-row fsb all zero
  expect_true(all(vapply(fb0, function(x) all(unname(x) == 0), logical(1L))))
  # non-zero: per-row fsb_i == Vpost_i * alpha, cross-checked against mirt's
  # posterior covariance (independent of fsL)
  acov <- mirt::fscores(m1, full.scores = TRUE, return.acov = TRUE,
                        mean = c(F1 = al))
  scorable <- mirt_scorable(m1, n1)
  for (i in c(1L, 40L, n1)) {
    k <- mirt_acov_k(scorable, i)
    Vpost_i <- as.numeric(as.matrix(acov[[k]])[1L, 1L])
    expect_identical(length(fb[[i]]), 1L)
    expect_equal(unname(as.numeric(fb[[i]])), Vpost_i * al, tolerance = 1e-8)
  }
  # the group-level alpha moment reflects the supplied prior
  expect_equal(unname(as.numeric(attr(fsp, "alpha"))), al)
  # the EAP scores are extracted under the prior, so they differ from default
  expect_false(isTRUE(all.equal(unname(fsp$fs_F1), unname(fs0$fs_F1))))
})

test_that("get_fs(): prior_mean -> fs_indiv(include_intercept) carries per-row fsb", {
  al <- 2.0
  fsp <- get_fs(m1, prior_mean = c(F1 = al))
  fb <- unlist(attr(fsp, "fsb"))
  iv <- fs_indiv(fsp, include_intercept = TRUE)
  expect_true("int_fs_F1" %in% names(iv))
  expect_equal(unname(iv$int_fs_F1), unname(fb), tolerance = 1e-8)
  # without the flag, no intercept column
  expect_false("int_fs_F1" %in% names(fs_indiv(fsp)))
})

test_that("get_fs(): prior_mean validation (mirt)", {
  expect_error(get_fs(m2, prior_mean = c(1.0, 0.5, 0.3)),
               regexp = "must have length 2")
  expect_error(get_fs(m1, prior_mean = c(F9 = 1.0)),
               regexp = "names must match")
  expect_error(get_fs(m1, prior_mean = c(F1 = NA_real_)),
               regexp = "finite")
})

test_that("get_fs(): unsupported options passed to a mirt fit are rejected (not silently ignored)", {
  # prior_cov, product, and method are not mirt options; each must error and
  # name the offending option rather than being a silent no-op.
  expect_error(get_fs(m1, prior_cov = matrix(2)), "prior_cov")
  expect_error(get_fs(m2, product = "F1:F2"), "product")
  expect_error(get_fs(m1, method = "ML"), "method")
  # multiple unsupported options are all named
  expect_error(get_fs(m1, prior_cov = matrix(2), method = "ML"),
               "prior_cov.*method|method.*prior_cov")
  # the supported arguments still work
  expect_no_error(get_fs(m1, prior_mean = c(F1 = 0.5)))
  expect_no_error(get_fs(m1, format = "unified"))
})

# ============================================================================
# 10. 2-factor with CORRELATED factors: psi must be the full mirt covariance
# ============================================================================

test_that("get_fs(): 2-factor correlated -> psi is the full mirt covariance, not I", {
  # correlated latent factors so the estimated factor covariance is non-diagonal
  set.seed(1235)
  nn <- 1000
  eta <- MASS::mvrnorm(nn, c(0, 0), diag(c(1, 0.75)), empirical = TRUE)
  th1 <- eta[, 1]
  th2 <- -1 + 0.5 * th1 + eta[, 2]
  dat1 <- mirt::simdata(a = matrix(1, 10), d = matrix(rnorm(10)), N = nn,
                        itemtype = "2PL", Theta = th1)
  dat2 <- mirt::simdata(a = matrix(runif(10, 0.5, 1.5)), d = matrix(rnorm(10)),
                        N = nn, itemtype = "2PL", Theta = th2)
  datc <- cbind(dat1, dat2)
  colnames(datc) <- paste0("Item_", 1:20)
  mc <- suppressWarnings(mirt::mirt(
    datc, "F1=1-10\nF2=11-20\nCOV=F1*F2\nCONSTRAIN=(1-10,a1)",
    itemtype = "2PL", verbose = FALSE
  ))
  fsc <- get_fs(mc)

  # independent reconstruction of the mirt factor covariance from GroupPars
  gp <- coef(mc)$GroupPars[1L, , drop = TRUE]
  cv <- as.numeric(unname(gp[grepl("^COV_", names(gp))]))  # COV_11, COV_21, COV_22
  Vp <- matrix(c(cv[1], cv[2], cv[2], cv[3]), 2)
  dimnames(Vp) <- list(c("F1", "F2"), c("F1", "F2"))
  expect_gt(abs(Vp[1L, 2L]), 0.2)   # the factors really do correlate in this data

  # (a) get_fs psi == the full mirt covariance, and NOT the identity
  expect_equal(unname(as.matrix(attr(fsc, "psi"))), unname(Vp), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(as.matrix(attr(fsc, "psi")), diag(2), tolerance = 1e-6)))

  # (b) per-row fsL == I - Vpost_i %*% solve(full cov) (uses the covariances)
  acov <- mirt::fscores(mc, full.scores = TRUE, return.acov = TRUE)
  i <- 10L
  k <- mirt_acov_k(mirt_scorable(mc, nrow(fsc)), i)
  Lexp <- diag(2) - as.matrix(acov[[k]]) %*% solve(Vp)
  expect_equal(unname(c(attr(fsc, "fsL")[[i]])), unname(c(Lexp)), tolerance = 1e-8)

  # (c) independent identity: cov(EAP scores) == full_cov - mean(Vpost), and
  #     NOT the identity-prior version (diag - mean Vpost) -- the old bug
  sc <- as.data.frame(mirt::fscores(mc, full.scores = TRUE))
  mvp <- Reduce("+", lapply(acov, as.matrix)) / length(acov)
  expect_equal(unname(cov(sc)), unname(Vp - mvp), tolerance = 0.02)
  expect_gt(max(abs(cov(sc) - (diag(2) - mvp))), 0.2)
})
