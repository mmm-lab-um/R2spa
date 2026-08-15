# Step 7: Numeric equivalence + timing vs lavPredict(acov=TRUE)
#
# These tests validate that the refactored get_fs() pipeline produces results
# numerically equivalent to lavaan::lavPredict(..., acov = TRUE), which is the
# canonical reference for factor score computation in lavaan.

library(lavaan)

# ---- Shared test fixtures ----
single_model_1f <- 'ind60 =~ x1 + x2 + x3'
single_model_2f <- 'ind60 =~ x1 + x2 + x3
                     dem60 =~ y1 + y2 + y3 + y4'

hs_model_1f <- 'visual  =~ x1 + x2 + x3'
hs_model_3f <- 'visual  =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed   =~ x7 + x8 + x9'

# ---- Helper: extract fs matrix from get_fs() output ----
get_fs_scores <- function(fs_out) {
  fs_names <- grep("^fs_[^_]+$", colnames(fs_out), value = TRUE)
  as.matrix(fs_out[, fs_names, drop = FALSE])
}

# ---- Helper: extract fsT/fsL from a single-group or first element ----
get_single_fsT <- function(fs_out) {
  a <- attr(fs_out, "fsT")
  if (is.list(a)) a[[1]] else a
}
get_single_fsL <- function(fs_out) {
  a <- attr(fs_out, "fsL")
  if (is.list(a)) a[[1]] else a
}

# ---- Helper: compute reference fsT/fsL from lavPredict acov ----
lav_group_matrices <- function(fit, g, method) {
  acov_lst <- attr(
    lavPredict(fit, type = "lv", method = method, acov = TRUE), "acov")
  pars <- lavInspect(fit, "est", drop.list.single.group = TRUE)
  psi <- if (lavInspect(fit, "ngroups") == 1) pars$psi else pars[[g]]$psi
  alpha <- if (lavInspect(fit, "ngroups") == 1) pars$alpha else pars[[g]]$alpha
  R2spa:::compute_lav_fs_matrices(acov_lst[[g]], psi = psi, alpha = alpha,
                                   method = method)
}

# ---- Helper: get lavPredict factor scores for a specific group ----
lav_group_fs <- function(fit, g, method) {
  lp <- lavPredict(fit, type = "lv", method = method)
  if (lavInspect(fit, "ngroups") == 1) {
    as.matrix(lp)
  } else {
    as.matrix(lp[[g]])
  }
}

# ============================================================================
# 1. Single-group, complete data
# ============================================================================

test_that("get_fs() ~ lavPredict: SG 1-factor complete regression", {
  fit <- cfa(single_model_1f, data = PoliticalDemocracy)
  fs_out <- get_fs(fit, method = "regression")
  ref <- lav_group_matrices(fit, 1, "regression")

  expect_equal(get_fs_scores(fs_out), lav_group_fs(fit, 1, "regression"),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(get_single_fsT(fs_out), ref$fsT, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(get_single_fsL(fs_out), ref$fsL, tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("get_fs() ~ lavPredict: SG 1-factor complete Bartlett", {
  fit <- cfa(single_model_1f, data = PoliticalDemocracy)
  fs_out <- get_fs(fit, method = "Bartlett")
  ref <- lav_group_matrices(fit, 1, "Bartlett")

  expect_equal(get_fs_scores(fs_out), lav_group_fs(fit, 1, "Bartlett"),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(get_single_fsT(fs_out), ref$fsT, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(get_single_fsL(fs_out), ref$fsL, tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("get_fs() ~ lavPredict: SG 2-factor complete regression", {
  fit <- cfa(single_model_2f, data = PoliticalDemocracy)
  fs_out <- get_fs(fit, method = "regression")
  ref <- lav_group_matrices(fit, 1, "regression")

  expect_equal(get_fs_scores(fs_out), lav_group_fs(fit, 1, "regression"),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(get_single_fsT(fs_out), ref$fsT, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(get_single_fsL(fs_out), ref$fsL, tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("get_fs() ~ lavPredict: SG 2-factor complete Bartlett", {
  fit <- cfa(single_model_2f, data = PoliticalDemocracy)
  fs_out <- get_fs(fit, method = "Bartlett")
  ref <- lav_group_matrices(fit, 1, "Bartlett")

  expect_equal(get_fs_scores(fs_out), lav_group_fs(fit, 1, "Bartlett"),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(get_single_fsT(fs_out), ref$fsT, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(get_single_fsL(fs_out), ref$fsL, tolerance = 1e-8,
               ignore_attr = TRUE)
})

# ============================================================================
# 2. Single-group, missing data
# ============================================================================

set.seed(1334)
hs_miss <- HolzingerSwineford1939
hs_miss[!rbinom(301, size = 1, prob = 0.7), 7] <- NA
hs_miss[!rbinom(301, size = 1, prob = 0.7), 8] <- NA
hs_miss[!rbinom(301, size = 1, prob = 0.7), 9] <- NA

test_that("get_fs blocks ~ lavPredict: SG missing data regression", {
  fit <- cfa(hs_model_3f, data = hs_miss, missing = "fiml")
  mp <- fit@Data@Mp

  if (is.null(mp)) skip("No missing-data patterns detected")

  case_idx <- mp$case.idx
  acov_rank <- rank(mp$id)
  acov_lst <- attr(
    lavPredict(fit, type = "lv", method = "regression", acov = TRUE), "acov")
  pars <- lavInspect(fit, "est")

  our_blocks <- R2spa:::get_fs_blocks.lavaan(fit, method = "regression",
                                              add_to_evfs = rep(0, 1))[[""]]

  for (m in seq_along(case_idx)) {
    mat_idx <- acov_rank[m]
    mats <- R2spa:::compute_lav_fs_matrices(
      acov = acov_lst[[mat_idx]],
      psi = pars$psi,
      alpha = pars$alpha,
      method = "regression"
    )
    expect_equal(our_blocks[[m]]$fsT, mats$fsT, tolerance = 1e-5,
                 ignore_attr = TRUE)
    expect_equal(our_blocks[[m]]$fsL, mats$fsL, tolerance = 1e-5,
                 ignore_attr = TRUE)
  }

  # Row-level factor scores
  fs_out <- get_fs(fit, method = "regression")
  ref_fs <- lav_group_fs(fit, 1, "regression")
  valid_rows <- rowSums(!is.na(ref_fs)) > 0
  if (any(valid_rows)) {
    expect_equal(get_fs_scores(fs_out)[valid_rows, , drop = FALSE],
                 ref_fs[valid_rows, , drop = FALSE],
                 tolerance = 1e-5, ignore_attr = TRUE)
  }
})

# ============================================================================
# 3. Multigroup, complete data
# ============================================================================

test_that("get_fs() ~ lavPredict: MG 1-factor complete regression", {
  fit <- cfa(hs_model_1f, data = HolzingerSwineford1939, group = "school")
  fs_list <- fs_to_group_list(get_fs(fit, method = "regression"))

  ngroups <- lavInspect(fit, "ngroups")
  for (g in seq_len(ngroups)) {
    lbl <- fit@Data@group.label[g]
    ref_mat <- lav_group_matrices(fit, g, "regression")

    expect_equal(get_fs_scores(fs_list[[lbl]]), lav_group_fs(fit, g, "regression"),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(get_single_fsT(fs_list[[lbl]]), ref_mat$fsT,
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(get_single_fsL(fs_list[[lbl]]), ref_mat$fsL,
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("get_fs() ~ lavPredict: MG 1-factor complete Bartlett", {
  fit <- cfa(hs_model_1f, data = HolzingerSwineford1939, group = "school")
  fs_list <- fs_to_group_list(get_fs(fit, method = "Bartlett"))

  ngroups <- lavInspect(fit, "ngroups")
  for (g in seq_len(ngroups)) {
    lbl <- fit@Data@group.label[g]
    ref_mat <- lav_group_matrices(fit, g, "Bartlett")

    expect_equal(get_fs_scores(fs_list[[lbl]]), lav_group_fs(fit, g, "Bartlett"),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(get_single_fsT(fs_list[[lbl]]), ref_mat$fsT,
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(get_single_fsL(fs_list[[lbl]]), ref_mat$fsL,
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("get_fs() ~ lavPredict: MG 3-factor complete regression", {
  fit <- cfa(hs_model_3f, data = HolzingerSwineford1939, group = "school")
  fs_list <- fs_to_group_list(get_fs(fit, method = "regression"))

  ngroups <- lavInspect(fit, "ngroups")
  for (g in seq_len(ngroups)) {
    lbl <- fit@Data@group.label[g]
    ref_mat <- lav_group_matrices(fit, g, "regression")

    expect_equal(get_fs_scores(fs_list[[lbl]]), lav_group_fs(fit, g, "regression"),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(get_single_fsT(fs_list[[lbl]]), ref_mat$fsT,
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(get_single_fsL(fs_list[[lbl]]), ref_mat$fsL,
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("get_fs() ~ lavPredict: MG 3-factor complete Bartlett", {
  fit <- cfa(hs_model_3f, data = HolzingerSwineford1939, group = "school")
  fs_list <- fs_to_group_list(get_fs(fit, method = "Bartlett"))

  ngroups <- lavInspect(fit, "ngroups")
  for (g in seq_len(ngroups)) {
    lbl <- fit@Data@group.label[g]
    ref_mat <- lav_group_matrices(fit, g, "Bartlett")

    expect_equal(get_fs_scores(fs_list[[lbl]]), lav_group_fs(fit, g, "Bartlett"),
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(get_single_fsT(fs_list[[lbl]]), ref_mat$fsT,
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(get_single_fsL(fs_list[[lbl]]), ref_mat$fsL,
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})

# ============================================================================
# 4. Multigroup, missing data
# ============================================================================

set.seed(42)
hs_miss_grp <- HolzingerSwineford1939
hs_miss_grp[!rbinom(301, size = 1, prob = 0.7), 7] <- NA
hs_miss_grp[!rbinom(301, size = 1, prob = 0.7), 8] <- NA
hs_miss_grp[!rbinom(301, size = 1, prob = 0.7), 9] <- NA

test_that("get_fs blocks ~ lavPredict: MG missing data regression", {
  fit <- cfa(hs_model_3f, data = hs_miss_grp, group = "school", missing = "fiml")
  ngroups <- lavInspect(fit, "ngroups")
  mp_lst <- fit@Data@Mp
  acov_all <- attr(
    lavPredict(fit, type = "lv", method = "regression", acov = TRUE), "acov")
  pars <- lavInspect(fit, "est", drop.list.single.group = FALSE)

  our_all_blocks <- R2spa:::get_fs_blocks.lavaan(
    fit, method = "regression", add_to_evfs = rep(0, ngroups))

  for (g in seq_len(ngroups)) {
    lbl <- fit@Data@group.label[g]
    mp <- mp_lst[[g]]
    our_blocks <- our_all_blocks[[lbl]]

    if (is.null(mp)) {
      ref_mat <- lav_group_matrices(fit, g, "regression")
      expect_equal(our_blocks[[1]]$fsT, ref_mat$fsT,
                   tolerance = 1e-5, ignore_attr = TRUE)
      expect_equal(our_blocks[[1]]$fsL, ref_mat$fsL,
                   tolerance = 1e-5, ignore_attr = TRUE)
    } else {
      case_idx <- mp$case.idx
      acov_rank <- rank(mp$id)

      for (m in seq_along(case_idx)) {
        mat_idx <- acov_rank[m]
        mats <- R2spa:::compute_lav_fs_matrices(
          acov = acov_all[[g]][[mat_idx]],
          psi = pars[[g]]$psi,
          alpha = pars[[g]]$alpha,
          method = "regression"
        )
        expect_equal(our_blocks[[m]]$fsT, mats$fsT,
                     tolerance = 1e-5, ignore_attr = TRUE)
        expect_equal(our_blocks[[m]]$fsL, mats$fsL,
                     tolerance = 1e-5, ignore_attr = TRUE)
      }
    }
  }
})

# ============================================================================
# 5. Timing: get_fs() vs lavPredict(acov=TRUE) on pre-fitted models
# ============================================================================

test_that("get_fs() comparable speed to lavPredict: SG 2-factor", {
  fit_fitted <- cfa(single_model_2f, data = PoliticalDemocracy)
  t1 <- system.time(for (i in 1:20) get_fs(fit_fitted))["elapsed"]
  t2 <- system.time(
    for (i in 1:20) lavPredict(fit_fitted, type = "lv", acov = TRUE))["elapsed"]
  message("SG get_fs: ", round(t1, 3), "s, lavPredict: ", round(t2, 3),
          "s, ratio: ", round(t1 / t2, 2))
  expect_lt(t1 / t2, 5)
})

test_that("get_fs() comparable speed to lavPredict: MG 3-factor", {
  fit_multi <- cfa(hs_model_3f, data = HolzingerSwineford1939,
                   group = "school")
  t1 <- system.time(for (i in 1:20) get_fs(fit_multi))["elapsed"]
  t2 <- system.time(
    for (i in 1:20) lavPredict(fit_multi, type = "lv", acov = TRUE))["elapsed"]
  message("MG get_fs: ", round(t1, 3), "s, lavPredict: ", round(t2, 3),
          "s, ratio: ", round(t1 / t2, 2))
  expect_lt(t1 / t2, 5)
})
