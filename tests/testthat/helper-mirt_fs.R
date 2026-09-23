# Shared binary/graded and single-/multigroup mirt identities. Per-file tests
# independently pin psi/alpha values; mirt's posterior covariances omit
# completely-missing rows, so map scorable positions back to full-row indices.

# Full-row indices of the scorable rows (the complement of completely-missing).
mirt_scorable <- function(fit, n) {
  cm <- mirt::extract.mirt(fit, "completely_missing")
  which(!seq_len(n) %in% (if (is.null(cm)) integer(0) else cm))
}

# Position of full-row i within the scorable rows (i.e. the acov list index);
# NA when i is completely missing.
mirt_acov_k <- function(scorable, i) match(i, scorable)

# Per-row psi: a plain matrix (SG) or, for MG, the entry of the per-group named
# list matching the row's `group` value (group_i ignored when psi is a matrix).
mirt_row_psi <- function(psi, group_i) {
  if (is.list(psi) && !is.matrix(psi)) psi[[as.character(group_i)]] else psi
}

# `samples`, when supplied, are full-row indices rather than acov positions.
assert_mirt_row_identities <- function(fs, fit, tol = 1e-8, samples = NULL,
                                       n_rows = NULL) {
  n <- if (is.null(n_rows)) nrow(fs) else n_rows
  q <- mirt::extract.mirt(fit, "nfact")
  fn <- mirt::extract.mirt(fit, "factorNames")
  psi <- attr(fs, "psi")
  alpha <- attr(fs, "alpha")
  scorable <- mirt_scorable(fit, n)
  Tl <- attr(fs, "fsT")
  Ll <- attr(fs, "fsL")

  # row count: get_fs and fs_indiv both keep every observation (and n_rows, the
  # data's row count, when supplied)
  expect_equal(nrow(fs), n)
  expect_equal(nrow(fs_indiv(fs)), n)

  # column set + order vs fs_indiv (fs_indiv drops the trailing `group` in MG),
  # and the score/SE columns agree row by row (NA rows compare as NA)
  ind <- fs_indiv(fs)
  keep_cols <- if ("group" %in% names(fs)) setdiff(names(fs), "group") else names(fs)
  expect_identical(keep_cols, names(ind))
  for (v in paste0("fs_", fn)) {
    expect_equal(unname(fs[[v]]), unname(ind[[v]]), tolerance = tol)
    expect_equal(unname(fs[[paste0(v, "_se")]]), unname(ind[[paste0(v, "_se")]]),
                 tolerance = tol)
  }

  # per-row fsT / fsL: n blocks, each q x q
  expect_true(is.list(Tl) && is.list(Ll))
  expect_length(Tl, n)
  expect_length(Ll, n)
  sq <- paste(c(q, q), collapse = "x")
  expect_true(all(vapply(Tl, function(x) paste(dim(x), collapse = "x"),
                         character(1L)) == sq))
  expect_true(all(vapply(Ll, function(x) paste(dim(x), collapse = "x"),
                         character(1L)) == sq))

  acov <- mirt::fscores(fit, full.scores = TRUE, return.acov = TRUE)

  # 1-factor regression closed forms (NA rows compare as NA)
  if (q == 1L) {
    full <- mirt::fscores(fit, full.scores = TRUE, full.scores.SE = TRUE)
    se <- full[, "SE_F1"]
    ev_exp <- (1 - se^2) * se^2
    expect_equal(unname(fs[["F1_by_fs_F1"]]), 1 - se^2, tolerance = tol)
    expect_equal(unname(fs[["ev_fs_F1"]]), ev_exp, tolerance = tol)
    expect_equal(unname(fs[["fs_F1_se"]]), sqrt(ev_exp), tolerance = tol)
    Tdiag <- vapply(Tl, function(Tm) Tm[1L, 1L], numeric(1L))
    expect_equal(unname(Tdiag), unname(fs[["ev_fs_F1"]]), tolerance = tol)
    expect_equal(unname(fs[["fs_F1"]]), unname(full[, "F1"]), tolerance = tol)
  }

  # per-row full-matrix identity against the shared regression engine
  idx <- unique(c(1L, 2L, ceiling(length(scorable) / 2), length(scorable)))
  rows <- if (is.null(samples)) scorable[idx] else samples
  for (i in rows) {
    k <- mirt_acov_k(scorable, i)
    expect_false(is.na(k), label = paste("row", i, "is scorable"))
    grp <- if ("group" %in% names(fs)) fs$group[i] else NULL
    m_i <- R2spa:::compute_lav_fs_matrices(
      as.matrix(acov[[k]]),
      psi = mirt_row_psi(psi, grp), alpha = alpha, method = "regression"
    )
    expect_equal(unname(c(attr(fs, "fsL")[[i]])), unname(c(m_i$fsL)), tolerance = tol)
    expect_equal(unname(c(attr(fs, "fsT")[[i]])), unname(c(m_i$fsT)), tolerance = tol)
  }
  invisible(NULL)
}
