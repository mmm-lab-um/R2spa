# get_fs() for mirt SingleGroupClass fits (per-observation factor scores).
#
# get_fs() on a fitted mirt SingleGroupClass returns a data.frame whose
# per-row columns are identical (in set AND order) to the output of fs_indiv().
# Each row carries its own implied loadings / error-covariance, computed
# per-observation from mirt's EAP posterior covariance (Vpost) through the
# shared regression-form matrix engine compute_lav_fs_matrices() (psi = I,
# alpha = 0):
#   fsL_i = I - Vpost_i        (1-factor: F1_by_fs_F1 == 1 - SE_i^2)
#   fsT_i = fsL_i %*% Vpost_i  (1-factor: (1 - SE_i^2) * SE_i^2)
#   fsb   = 0                  (alpha = 0)
# where diag(Vpost_i) == SE_i^2 (the SE column of mirt::fscores(full.scores.SE
# = TRUE)). In addition to the per-row list-valued attributes (fsL/fsT) and the
# `mirt_per_obs` marker, psi/alpha are attached as group-level moments.
#
# These tests cover: S3 dispatch + the MultipleGroupClass guard, the
# column-set / row-count contract, the 1-factor regression identities, the
# attribute shapes, the 2-factor off-diagonal identity against
# compute_lav_fs_matrices(), the column naming order, and the
# completely-missing-row handling (R2spa NA-row convention).
#
# mirt is Suggests-only; every call is namespaced (mirt::). No library().

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
ind <- fs_indiv(fs)

# mirt's per-observation SEs for the 1-factor fit (SE_F1 = sqrt(diag(Vpost))).
se1 <- mirt::fscores(m1, full.scores = TRUE, full.scores.SE = TRUE)[, "SE_F1"]

# ---- helpers ------------------------------------------------------------
# Map a full (n-row) result row index i to the 1-based position of its entry in
# a mirt acov list. mirt::fscores(return.acov = TRUE) skips completely-missing
# rows, so the scorable rows are the complement of `completely_missing` and
# acov[[k]] is the k-th scorable row.
acov_index_for_row <- function(acov, cm, i, n) {
  keep <- !seq_len(n) %in% (if (is.null(cm)) integer(0) else cm)
  scorsc <- which(keep)
  k <- which(scorsc == i)
  if (length(k) != 1L) stop("row ", i, " has no scorable posterior covariance")
  k
}

# ============================================================================
# 1. S3 dispatch + type + MultipleGroupClass guard
# ============================================================================

test_that("get_fs(): S3 dispatch -- mirt S4 objects route to the mirt methods", {
  expect_true(inherits(m1, "SingleGroupClass"))
  # SingleGroupClass -> a data.frame carrying the per-observation marker
  expect_true(is.data.frame(fs))
  expect_true(isTRUE(attr(fs, "mirt_per_obs")))
  # The MultipleGroupClass guard is registered and fires with a clear message
  expect_true(exists("get_fs.MultipleGroupClass", where = asNamespace("R2spa")))
  stub <- structure(list(), class = "MultipleGroupClass")
  expect_error(
    get_fs(stub),
    regexp = "Multi-group mirt models are not supported by get_fs\\(\\)"
  )
  # Non-mirt input still routes to get_fs.default (unchanged behaviour)
  expect_error(get_fs(42L), regexp = "not implemented for objects of class")
})

# ============================================================================
# 2. Column-set / order identity with fs_indiv()
# ============================================================================

test_that("get_fs(): mirt per-row column set and order match fs_indiv()", {
  expect_identical(sort(names(fs)), sort(names(fs_indiv(fs))))
  expect_identical(sort(names(fs2)), sort(names(fs_indiv(fs2))))
  # and the exact column ORDER is preserved, not just the set
  expect_identical(names(fs), names(ind))
  expect_identical(names(fs2), names(fs_indiv(fs2)))
})

# ============================================================================
# 3. Row count
# ============================================================================

test_that("get_fs(): one row per observation, preserved through fs_indiv()", {
  expect_equal(nrow(fs), nrow(d1))
  expect_equal(nrow(ind), nrow(d1))
  expect_equal(nrow(fs2), nrow(d2))
  expect_equal(nrow(fs_indiv(fs2)), nrow(d2))
})

# ============================================================================
# 4. 1-factor: implied loading == 1 - SE_i^2
# ============================================================================

test_that("get_fs(): 1-factor implied loading F1_by_fs_F1 == 1 - SE_i^2 (all rows)", {
  expect_equal(unname(fs[["F1_by_fs_F1"]]), 1 - se1^2, tolerance = 1e-8)
  expect_true(all(fs[["F1_by_fs_F1"]] < 1))
  # the per-row fsL diagonal is also 1 - SE_i^2, and strictly < 1
  Ldiag <- vapply(attr(fs, "fsL"), function(Lm) Lm[1L, 1L], numeric(1L))
  expect_true(all(Ldiag < 1))
  expect_equal(unname(Ldiag), 1 - se1^2, tolerance = 1e-8)
})

# ============================================================================
# 5. 1-factor: ev and score SE == (1 - SE_i^2) * SE_i^2 (and its square root)
# ============================================================================

test_that("get_fs(): 1-factor ev and score SE == (1 - SE_i^2) * SE_i^2 (and sqrt)", {
  ev_exp <- (1 - se1^2) * se1^2
  expect_equal(unname(fs[["ev_fs_F1"]]), ev_exp, tolerance = 1e-8)
  expect_equal(unname(fs[["fs_F1_se"]]), sqrt(ev_exp), tolerance = 1e-8)
  # the per-row fsT diagonal equals the per-row ev value
  Tdiag <- vapply(attr(fs, "fsT"), function(Tm) Tm[1L, 1L], numeric(1L))
  expect_equal(unname(Tdiag), unname(fs[["ev_fs_F1"]]), tolerance = 1e-8)
  # the score column is mirt's EAP posterior mean
  eap <- mirt::fscores(m1, full.scores = TRUE, full.scores.SE = TRUE)[, "F1"]
  expect_equal(unname(fs[["fs_F1"]]), unname(eap), tolerance = 1e-8)
})

# ============================================================================
# 6. Attribute shapes
# ============================================================================

test_that("get_fs(): per-observation list attributes + group-level psi/alpha/fsb", {
  q1 <- mirt::extract.mirt(m1, "nfact")
  fn1 <- mirt::extract.mirt(m1, "factorNames")
  n1 <- nrow(fs)
  Tl <- attr(fs, "fsT")
  Ll <- attr(fs, "fsL")
  # fsT / fsL are per-row lists, each element a q x q matrix
  expect_true(is.list(Tl) && is.list(Ll))
  expect_length(Tl, n1)
  expect_length(Ll, n1)
  sq <- paste(c(q1, q1), collapse = "x")
  expect_true(all(vapply(Tl, function(x) paste(dim(x), collapse = "x"),
                         character(1L)) == sq))
  expect_true(all(vapply(Ll, function(x) paste(dim(x), collapse = "x"),
                         character(1L)) == sq))
  # the fsT diagonal equals the per-row ev value
  Tdiag <- vapply(Tl, function(Tm) Tm[1L, 1L], numeric(1L))
  expect_equal(unname(Tdiag), unname(fs[["ev_fs_F1"]]), tolerance = 1e-8)
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
  # fsb: a named zero vector
  fsb <- attr(fs, "fsb")
  expect_length(fsb, q1)
  expect_named(fsb, paste0("fs_", fn1))
  expect_true(all(unname(fsb) == 0))
})

# ============================================================================
# 7. 2-factor: off-diagonals are NOT identity; per-row fsL/fsT == engine
# ============================================================================

test_that("get_fs(): 2-factor off-diagonals non-identity; per-row fsL/fsT == engine", {
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
  # per-row full-matrix identity against the shared regression engine for a
  # few sampled (scorable) rows, reconciling the acov row order
  acov2 <- mirt::fscores(m2, full.scores = TRUE, return.acov = TRUE)
  cm <- mirt::extract.mirt(m2, "completely_missing")
  if (is.null(cm)) cm <- integer(0)
  psi2 <- diag(q2)
  rownames(psi2) <- colnames(psi2) <- fn2
  alpha2 <- setNames(rep(0, q2), fn2)
  for (i in c(1L, 5L, 60L, nrow(fs2))) {
    k <- acov_index_for_row(acov2, cm, i, nrow(fs2))
    m_i <- R2spa:::compute_lav_fs_matrices(
      as.matrix(acov2[[k]]),
      psi = psi2, alpha = alpha2,
      method = "regression"
    )
    L_act <- attr(fs2, "fsL")[[i]]
    T_act <- attr(fs2, "fsT")[[i]]
    expect_equal(unname(c(L_act)), unname(c(m_i$fsL)), tolerance = 1e-8)
    expect_equal(unname(c(T_act)), unname(c(m_i$fsT)), tolerance = 1e-8)
    # the per-row fsL is not the identity (off-diagonal cross-loadings live here)
    expect_false(isTRUE(all.equal(unname(c(L_act)), c(diag(q2)), tolerance = 1e-8)))
  }
})

# ============================================================================
# 8. Completely-missing rows
# ============================================================================

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
