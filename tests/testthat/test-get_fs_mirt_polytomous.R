# get_fs() for mirt SingleGroupClass fits with POLYTOMOUS (graded-response)
# items (itemtype = "graded").
#
# Binary mirt is covered in test-get_fs_mirt.R and
# test-get_fs_mirt_multigroup.R; graded items had no coverage. The extraction
# path is itemtype-agnostic: get_fs() scores a fitted SingleGroupClass by its
# EAP posterior means and feeds the per-observation EAP posterior covariance
# (Vpost_i, from mirt::fscores(return.acov = TRUE)) into the shared
# regression-form engine compute_lav_fs_matrices() -- exactly as for binary.
# The regression-form closed forms therefore hold for graded items too, and
# for a single factor, with SE_i^2 = diag(Vpost_i):
#   fs_F1_by_fs_F1_i == 1 - SE_i^2
#   ev_fs_F1_i       == (1 - SE_i^2) * SE_i^2
#   fs_F1_se_i       == sqrt((1 - SE_i^2) * SE_i^2)
#   fs_F1_i          == mirt's EAP posterior mean
#
# Graded data are generated manually (item-by-item uniform-threshold sampling
# from the GRM category CDFs) rather than via mirt::simdata. Note a single
# 4-category graded item is under-identified in a 1-factor model (4 free
# parameters vs 3 df), so each factor gets 5 items, mirroring the binary
# fixture size.
#
# These tests cover: S3 dispatch + the per-observation contract, the exact
# 1-factor and 2-factor column sets (in order), the 1-factor regression
# identities against mirt's posterior covariance, the per-row attribute
# shapes, the 2-factor off-diagonals against the shared engine (and against
# the independent identity fsL_i == I - Vpost_i %*% solve(psi) with psi the
# FULL estimated factor covariance), and column/score equivalence with
# fs_indiv().
#
# mirt is Suggests-only; every call is namespaced (mirt::). No library().

skip_if_not_installed("mirt")

# ---- fixtures -------------------------------------------------------------
# Returns the integer category codes (0..K-1) for each person: draw one
# uniform per person and count how many cumulative thresholds it falls below
# (GRM sampling: the shared uniform makes the category monotone in theta).
gen_grm <- function(theta, a, d) {
  cums <- 1 / (1 + exp(-outer(theta, d, function(th, dk) a * th - dk)))
  u <- runif(nrow(cums))
  rowSums(u <= cums)
}
# K 4-category graded items (3 thresholds each) on one latent vector.
gen_graded_items <- function(theta, K, a0, sd_a) {
  do.call(cbind, lapply(seq_len(K), function(i) {
    gen_grm(theta, a0 + rnorm(1, 0, sd_a), sort(rnorm(3, 0, 1.2)))
  }))
}

set.seed(2026)
n <- 200

# 1-factor, 5 graded items (4 categories each)
theta1 <- rnorm(n)
dat1f <- as.data.frame(gen_graded_items(theta1, 5L, 1.5, 0.3))
colnames(dat1f) <- paste0("x", 1:5)
m1f <- suppressWarnings(
  mirt::mirt(dat1f, 1, itemtype = "graded", verbose = FALSE)
)

# 2-factor, 5 graded items per factor, correlated latent vectors
set.seed(2026)
th1 <- rnorm(n)
th2 <- 0.5 * th1 + rnorm(n)
dat2f <- data.frame(gen_graded_items(th1, 5L, 1.0, 0.3),
                    gen_graded_items(th2, 5L, 1.2, 0.3))
colnames(dat2f) <- paste0("x", 1:10)
m2f <- suppressWarnings(
  mirt::mirt(dat2f, "F1 = 1-5\nF2 = 6-10\nCOV = F1*F2",
             itemtype = "graded", verbose = FALSE)
)

# The fixtures are genuinely polytomous: every item uses 4 categories.
expect_true(all(apply(as.matrix(dat1f), 2L, function(x) length(unique(x)) == 4L)))
expect_true(all(apply(as.matrix(dat2f), 2L, function(x) length(unique(x)) == 4L)))

fs1 <- get_fs(m1f)
fs2 <- get_fs(m2f)
ind1 <- fs_indiv(fs1)

# mirt's per-observation posterior covariances (no missing data, so one
# q x q matrix per observation, in data-row order).
acov1 <- mirt::fscores(m1f, full.scores = TRUE, return.acov = TRUE)
acov2 <- mirt::fscores(m2f, full.scores = TRUE, return.acov = TRUE)
se1sq <- vapply(acov1, function(M) as.numeric(as.matrix(M)[1L, 1L]),
                numeric(1L), USE.NAMES = FALSE)

# ============================================================================
# 1. S3 dispatch + per-observation contract
# ============================================================================

test_that("get_fs(): graded mirt SingleGroupClass dispatches to the mirt method", {
  expect_true(inherits(m1f, "SingleGroupClass"))
  # the fit really is graded (per-item itemtype)
  expect_true(all(mirt::extract.mirt(m1f, "itemtype") == "graded"))
  # per-observation data frame, one row per person, marker set
  expect_true(is.data.frame(fs1))
  expect_true(isTRUE(attr(fs1, "mirt_per_obs")))
  expect_equal(nrow(fs1), nrow(dat1f))
  expect_equal(nrow(fs2), nrow(dat2f))
  # non-mirt input still routes to get_fs.default (unchanged behaviour)
  expect_error(get_fs(42L), regexp = "not implemented for objects of class")
})

# ============================================================================
# 2. Exact column sets (set AND order)
# ============================================================================

test_that("get_fs(): 1-factor graded column set and order are exact", {
  expect_identical(
    colnames(fs1),
    c("fs_F1", "fs_F1_se", "F1_by_fs_F1", "ev_fs_F1")
  )
})

test_that("get_fs(): 2-factor graded column set and order are exact", {
  expect_identical(
    colnames(fs2),
    c("fs_F1", "fs_F2", "fs_F1_se", "fs_F2_se",
      "F1_by_fs_F1", "F1_by_fs_F2", "F2_by_fs_F1", "F2_by_fs_F2",
      "ev_fs_F1", "ecov_fs_F2_fs_F1", "ev_fs_F2")
  )
})

# ============================================================================
# 3. 1-factor regression identities (IRF-independent)
# ============================================================================

test_that("get_fs(): 1-factor graded implied loading == 1 - SE_i^2 (all rows)", {
  expect_equal(unname(fs1[["F1_by_fs_F1"]]), 1 - se1sq, tolerance = 1e-8)
  expect_true(all(fs1[["F1_by_fs_F1"]] < 1))
  # the per-row fsL diagonal is also 1 - SE_i^2, and strictly < 1
  Ldiag <- vapply(attr(fs1, "fsL"), function(Lm) Lm[1L, 1L], numeric(1L))
  expect_true(all(Ldiag < 1))
  expect_equal(unname(Ldiag), 1 - se1sq, tolerance = 1e-8)
})

test_that("get_fs(): 1-factor graded ev/SE == (1 - SE_i^2) * SE_i^2 (and sqrt)", {
  ev_exp <- (1 - se1sq) * se1sq
  expect_equal(unname(fs1[["ev_fs_F1"]]), ev_exp, tolerance = 1e-8)
  expect_equal(unname(fs1[["fs_F1_se"]]), sqrt(ev_exp), tolerance = 1e-8)
  # the per-row fsT diagonal equals the per-row ev value
  Tdiag <- vapply(attr(fs1, "fsT"), function(Tm) Tm[1L, 1L], numeric(1L))
  expect_equal(unname(Tdiag), unname(fs1[["ev_fs_F1"]]), tolerance = 1e-8)
  # the score column is mirt's EAP posterior mean
  eap <- mirt::fscores(m1f, full.scores = TRUE, full.scores.SE = TRUE)[, "F1"]
  expect_equal(unname(fs1[["fs_F1"]]), unname(eap), tolerance = 1e-8)
})

# ============================================================================
# 4. Per-row attribute shapes + group-level moments
# ============================================================================

test_that("get_fs(): graded per-row list attributes + group-level psi/alpha/fsb", {
  n1 <- nrow(fs1)
  Tl <- attr(fs1, "fsT")
  Ll <- attr(fs1, "fsL")
  fsb <- attr(fs1, "fsb")
  # fsT / fsL are per-row lists, each element a 1 x 1 matrix
  expect_true(is.list(Tl) && is.list(Ll))
  expect_length(Tl, n1)
  expect_length(Ll, n1)
  expect_true(all(vapply(Tl, function(x) paste(dim(x), collapse = "x"),
                         character(1L)) == "1x1"))
  expect_true(all(vapply(Ll, function(x) paste(dim(x), collapse = "x"),
                         character(1L)) == "1x1"))
  # fsb: a per-row list of named 1-vectors, all zero (alpha = 0)
  expect_true(is.list(fsb))
  expect_length(fsb, n1)
  expect_true(all(vapply(fsb, function(x) {
    length(x) == 1L && identical(names(x), "fs_F1") && all(unname(x) == 0)
  }, logical(1L))))
  # fs_pattern: label = 1..n, pat = NULL (no missing data)
  fp <- attr(fs1, "fs_pattern")
  expect_equal(fp$label, seq_len(n1))
  expect_null(fp$pat)
  # psi == diag(1), named F1 (a 1-factor mirt model fixes COV_11 = 1)
  psi <- attr(fs1, "psi")
  expect_equal(as.matrix(psi), diag(1), tolerance = 1e-10, ignore_attr = TRUE)
  expect_identical(rownames(psi), "F1")
  expect_identical(colnames(psi), "F1")
  # alpha: a named zero vector
  alpha <- attr(fs1, "alpha")
  expect_length(alpha, 1L)
  expect_named(alpha, "F1")
  expect_true(all(unname(alpha) == 0))
})

# ============================================================================
# 5. 2-factor: off-diagonals non-identity; per-row fsL/fsT == shared engine
# ============================================================================

test_that("get_fs(): 2-factor graded off-diagonals non-identity; per-row fsL/fsT == engine", {
  fn2 <- mirt::extract.mirt(m2f, "factorNames")
  fsn2 <- paste0("fs_", fn2)
  expect_equal(mirt::extract.mirt(m2f, "nfact"), 2L)
  expect_setequal(fn2, c("F1", "F2"))
  expect_true(isTRUE(attr(fs2, "mirt_per_obs")))
  # off-diagonal implied loadings / error covariance are not all zero
  expect_true(any(abs(fs2[["F1_by_fs_F2"]]) > 1e-8))
  expect_true(any(abs(fs2[["F2_by_fs_F1"]]) > 1e-8))
  expect_true(any(abs(fs2[["ecov_fs_F2_fs_F1"]]) > 1e-8))
  # implied-loadings columns are column-major per latent
  ld_cols <- grep("_by_fs_", names(fs2), value = TRUE)
  expect_identical(
    ld_cols,
    unlist(lapply(seq_len(2L), function(j) paste(fn2[j], fsn2, sep = "_by_")))
  )
  # ev/ecov columns are the lower triangle in i-outer / j<=i-inner order
  ev_cols <- grep("^ev_|^ecov_", names(fs2), value = TRUE)
  expect_identical(ev_cols, c("ev_fs_F1", "ecov_fs_F2_fs_F1", "ev_fs_F2"))
  # the factors really do correlate in this data (from GroupPars)
  gp <- coef(m2f)$GroupPars[1L, , drop = TRUE]
  cv <- as.numeric(unname(gp[grepl("^COV_", names(gp))]))
  Vp <- matrix(c(cv[1], cv[2], cv[2], cv[3]), 2)
  dimnames(Vp) <- list(c("F1", "F2"), c("F1", "F2"))
  expect_gt(abs(Vp[1L, 2L]), 0.2)
  # get_fs psi == the full estimated factor covariance
  expect_equal(unname(as.matrix(attr(fs2, "psi"))), unname(Vp), tolerance = 1e-8)
  # per-row full-matrix identity against the shared regression engine for a
  # few sampled rows (no missing data: acov[[i]] is row i)
  alpha2 <- setNames(rep(0, 2L), c("F1", "F2"))
  for (i in c(1L, 5L, 60L, nrow(fs2))) {
    m_i <- R2spa:::compute_lav_fs_matrices(
      as.matrix(acov2[[i]]),
      psi = Vp, alpha = alpha2,
      method = "regression"
    )
    L_act <- attr(fs2, "fsL")[[i]]
    T_act <- attr(fs2, "fsT")[[i]]
    expect_equal(unname(c(L_act)), unname(c(m_i$fsL)), tolerance = 1e-8)
    expect_equal(unname(c(T_act)), unname(c(m_i$fsT)), tolerance = 1e-8)
    # independent identity: fsL_i == I - Vpost_i %*% solve(psi), using the
    # FULL (non-diagonal) factor covariance
    Lexp <- diag(2L) - as.matrix(acov2[[i]]) %*% solve(Vp)
    expect_equal(unname(c(L_act)), unname(c(Lexp)), tolerance = 1e-8)
    # the per-row fsL is not the identity (cross-factor cross-loadings live here)
    expect_false(isTRUE(all.equal(unname(c(L_act)), c(diag(2L)), tolerance = 1e-8)))
  }
  # the score columns are mirt's EAP posterior means (both factors)
  eap2 <- mirt::fscores(m2f, full.scores = TRUE)
  expect_equal(unname(fs2[["fs_F1"]]), unname(as.data.frame(eap2)[["F1"]]),
               tolerance = 1e-8)
  expect_equal(unname(fs2[["fs_F2"]]), unname(as.data.frame(eap2)[["F2"]]),
               tolerance = 1e-8)
})

# ============================================================================
# 6. Equivalence with fs_indiv() (set, order, and values)
# ============================================================================

test_that("get_fs(): graded per-row column set and order match fs_indiv()", {
  expect_identical(sort(names(fs1)), sort(names(ind1)))
  expect_identical(names(fs1), names(ind1))
  expect_identical(names(fs2), names(fs_indiv(fs2)))
  # scores and SEs agree row by row
  expect_equal(unname(fs1[["fs_F1"]]), unname(ind1[["fs_F1"]]), tolerance = 1e-8)
  expect_equal(unname(fs1[["fs_F1_se"]]), unname(ind1[["fs_F1_se"]]),
               tolerance = 1e-8)
  ind2 <- fs_indiv(fs2)
  expect_equal(unname(fs2[["fs_F1"]]), unname(ind2[["fs_F1"]]), tolerance = 1e-8)
  expect_equal(unname(fs2[["fs_F2_se"]]), unname(ind2[["fs_F2_se"]]),
               tolerance = 1e-8)
  expect_equal(nrow(ind2), nrow(dat2f))
})
