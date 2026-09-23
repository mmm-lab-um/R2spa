# Multigroup mirt fixtures. Missing rows retain NA-valued score blocks and
# group; the shared per-row identities live in helper-mirt_fs.R.

# ---- fixtures -----------------------------------------------------------
set.seed(2024)
n <- 120

# 1-factor, 5 binary items, 2 groups ("A" / "B"); one completely-missing row.
sim2 <- function(a, d, th, itemtype = "2PL") {
  as.data.frame(mirt::simdata(a = a, d = d, N = length(th),
                              itemtype = itemtype, Theta = th))
}
a1 <- matrix(runif(5, 0.5, 1.5), 5, 1L)
d1 <- matrix(rnorm(5), 5, 1L)
dat1 <- rbind(sim2(a1, d1, rnorm(n)), sim2(a1, d1, rnorm(n)))
grp1 <- factor(rep(c("A", "B"), each = n))
# a few cell-wise NAs + one completely-missing row (row 1)
mtmp <- as.matrix(dat1)
mtmp[sample(length(mtmp), 8L)] <- NA
dat1 <- as.data.frame(mtmp)
dat1[1L, ] <- NA

mg1 <- suppressWarnings(mirt::multipleGroup(dat1, 1L, group = grp1,
                                            invariance = "slopes",
                                            verbose = FALSE))
# 2-factor, 8 binary items, 2 groups ("1" / "2").
a2 <- matrix(0, 8L, 2L)
a2[1:4, 1L] <- runif(4L, 0.5, 1.5)
a2[5:8, 2L] <- runif(4L, 0.5, 1.5)
d2 <- matrix(rnorm(8L), 8L, 1L)
dat2 <- rbind(sim2(a2, d2, cbind(rnorm(n), rnorm(n))),
              sim2(a2, d2, cbind(rnorm(n), rnorm(n))))
grp2 <- factor(rep(1:2, each = n))
mg2 <- suppressWarnings(mirt::multipleGroup(dat2, 2L, group = grp2,
                                            invariance = "slopes",
                                            verbose = FALSE))

fs1 <- get_fs(mg1)
fs2 <- get_fs(mg2)
ind1 <- fs_indiv(fs1)

# scorable rows of the 1-factor fit (complement of completely-missing)
cm1 <- mirt::extract.mirt(mg1, "completely_missing")
sc1 <- which(!seq_len(nrow(fs1)) %in% cm1)

test_that("get_fs(): multi-group mirt routes to the MultipleGroupClass method", {
  expect_true(inherits(mg1, "MultipleGroupClass"))
  expect_true(is.data.frame(fs1))
  expect_true(isTRUE(attr(fs1, "mirt_per_obs")))
  # one row per observation (full, includes the completely-missing row)
  expect_equal(nrow(fs1), 2L * n)
  expect_equal(nrow(ind1), 2L * n)
})

test_that("get_fs(): multi-group columns = single-group columns + trailing group", {
  expect_identical(names(fs1),
                   c("fs_F1", "fs_F1_se", "F1_by_fs_F1", "ev_fs_F1", "group"))
  expect_identical(names(fs2),
                   c("fs_F1", "fs_F2", "fs_F1_se", "fs_F2_se",
                     "F1_by_fs_F1", "F1_by_fs_F2", "F2_by_fs_F1", "F2_by_fs_F2",
                     "ev_fs_F1", "ecov_fs_F2_fs_F1", "ev_fs_F2", "group"))
  expect_identical(names(fs1)[ncol(fs1)], "group")
  expect_identical(names(fs2)[ncol(fs2)], "group")
  # the per-row block (everything but `group`) matches fs_indiv() exactly
  expect_identical(setdiff(names(fs1), "group"), names(ind1))
})

test_that("get_fs(): multi-group `group` column carries the model groups", {
  expect_true(is.factor(fs1$group))
  expect_identical(levels(fs1$group), c("A", "B"))
  # scorable rows carry their own group; the completely-missing row is NA
  grp_exp <- as.character(grp1)
  grp_exp[cm1] <- NA
  expect_equal(unname(as.character(fs1$group)), grp_exp)
  expect_identical(levels(fs2$group), c("1", "2"))
  expect_true(all(as.character(fs2$group) %in% c("1", "2")))
})

test_that("get_fs(): multi-group per-row identities (1-factor + 2-factor)", {
  assert_mirt_row_identities(fs1, mg1, n_rows = 2L * n)
  assert_mirt_row_identities(fs2, mg2, n_rows = 2L * n)
})

test_that("get_fs(): multi-group per-group psi + alpha/fsb/fs_pattern (per-row shapes in helper)", {
  n1 <- nrow(fs1)
  # psi is a NAMED LIST (one 1x1 per group), keyed by the model group names --
  # the one structural difference from the single-group (single matrix) result.
  psig <- attr(fs1, "psi")
  expect_true(is.list(psig) && !is.matrix(psig))
  expect_length(psig, 2L)
  expect_named(psig, c("A", "B"))
  # metric invariance fixes each group's 1-factor covariance to 1
  expect_equal(unname(vapply(psig, function(V) V[1L, 1L], numeric(1L))),
               c(1, 1), tolerance = 1e-10)

  # alpha: named zero vector; fs_pattern: label 1..n, pat NULL
  alpha <- attr(fs1, "alpha")
  expect_length(alpha, 1L); expect_named(alpha, "F1")
  expect_true(all(unname(alpha) == 0))
  fp <- attr(fs1, "fs_pattern")
  expect_equal(fp$label, seq_len(n1))
  expect_null(fp$pat)

  # fsb: a per-row list, all zeros (alpha = 0)
  fsb <- attr(fs1, "fsb")
  expect_length(fsb, n1)
  # scorable rows are all-zero (alpha = 0); the missing row's NA fsb is
  # asserted in the completely-missing test below
  expect_true(all(vapply(fsb[sc1], function(b) all(unname(b) == 0),
                    logical(1L))))
})

test_that("get_fs(): multi-group 2-factor per-group psi (engine in helper)", {
  fn2 <- c("F1", "F2")
  gname <- mirt::extract.mirt(mg2, "groupNames")
  psig <- attr(fs2, "psi")
  expect_named(psig, gname)
  # per-group psi == mirt_full_cov of that group; dimnames are the factor names
  for (k in seq_along(gname)) {
    expect_equal(unname(as.matrix(psig[[k]])),
                 unname(as.matrix(R2spa:::mirt_full_cov(mirt::extract.group(mg2, k)))),
                 tolerance = 1e-8)
    expect_identical(rownames(psig[[k]]), fn2)
    expect_identical(colnames(psig[[k]]), fn2)
  }
})

test_that("get_fs(): multi-group result resolves through fs_indiv()", {
  expect_true(is.data.frame(ind1))
  expect_false("group" %in% names(ind1))
  expect_equal(nrow(ind1), nrow(fs1))
  expect_equal(ncol(ind1), ncol(fs1) - 1L)
  # per-row values agree with the get_fs() data frame
  expect_equal(unname(ind1[["F1_by_fs_F1"]]), unname(fs1[["F1_by_fs_F1"]]),
               tolerance = 1e-8)
  expect_equal(unname(ind1[["ev_fs_F1"]]), unname(fs1[["ev_fs_F1"]]),
               tolerance = 1e-8)
})

test_that("get_fs(): multi-group completely-missing row is all-NA + NA group", {
  na_rows <- cm1  # extract.mirt("completely_missing") gives original row indices
  expect_true(length(na_rows) >= 1L)
  for (i in na_rows) {
    numeric_cols <- names(fs1)[names(fs1) != "group"]
    expect_true(all(is.na(as.numeric(fs1[i, numeric_cols, drop = FALSE]))))
    expect_true(is.na(fs1$group[i]))
    # the per-row block is all-NA too
    expect_true(all(is.na(attr(fs1, "fsT")[[i]])))
    expect_true(all(is.na(attr(fs1, "fsL")[[i]])))
  }
})
