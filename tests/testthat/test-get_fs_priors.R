library(lavaan)
library(lme4)

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

fs_prior <- get_fs(prior_fit, prior_mean = pm, prior_cov = pc,
                   format = "list")
fs_none <- get_fs(prior_fit, format = "list")

hand_scores <- function(psi, alpha, est = prior_est, y = prior_data) {
  compute_fscore(y,
                 lambda = est$lambda,
                 theta = est$theta,
                 psi = psi,
                 nu = est$nu,
                 alpha = alpha,
                 method = "regression",
                 fs_matrices = TRUE)
}

########## Testing section ############

test_that("prior_mean + prior_cov match manual compute_fscore()", {
  hand <- hand_scores(psi = pc, alpha = pm)
  expect_equal(fs_prior[, c("fs_ind60", "fs_dem60")],
               as.data.frame(hand),
               ignore_attr = TRUE)
  expect_equal(attr(fs_prior, "fsT"), attr(hand, "fsT"))
  expect_equal(attr(fs_prior, "fsL"), attr(hand, "fsL"))
  expect_equal(attr(fs_prior, "fsb"), attr(hand, "fsb"))
  expect_equal(attr(fs_prior, "scoring_matrix"),
               attr(hand, "scoring_matrix"))
})

test_that("prior_cov only keeps model-estimated alpha", {
  fs_pc <- get_fs(prior_fit, prior_cov = pc, format = "list")
  hand <- hand_scores(psi = pc, alpha = prior_est$alpha)
  expect_equal(fs_pc[, c("fs_ind60", "fs_dem60")],
               as.data.frame(hand),
               ignore_attr = TRUE)
  expect_equal(attr(fs_pc, "fsb"), attr(hand, "fsb"))
})

test_that("prior_mean only keeps model-estimated psi", {
  fs_pm <- get_fs(prior_fit, prior_mean = pm, format = "list")
  hand <- hand_scores(psi = prior_est$psi, alpha = pm)
  expect_equal(fs_pm[, c("fs_ind60", "fs_dem60")],
               as.data.frame(hand),
               ignore_attr = TRUE)
  expect_equal(attr(fs_pm, "fsT"), attr(hand, "fsT"))
})

test_that("named priors are reordered to model order", {
  fs_rev <- get_fs(prior_fit, prior_mean = rev(pm),
                   prior_cov = pc[2:1, 2:1],
                   format = "list")
  expect_equal(fs_rev, fs_prior)
})

test_that("unnamed priors are used in model order", {
  fs_unnamed <- get_fs(prior_fit, prior_mean = c(0.3, -0.4),
                       prior_cov = pc,
                       format = "list")
  expect_equal(fs_unnamed, fs_prior)
})

test_that("prior input validation errors", {
  expect_error(
    get_fs(prior_fit, prior_mean = c(0.1, 0.2, 0.3)),
    "length"
  )
  expect_error(
    get_fs(prior_fit, prior_mean = c(foo = 0.1, bar = 0.2)),
    "names must match"
  )
  expect_error(
    get_fs(prior_fit, prior_mean = c(ind60 = 1)),
    "names must match"
  )
  expect_error(
    get_fs(prior_fit, prior_mean = c(ind60 = NA, dem60 = 0)),
    "finite"
  )
  expect_error(
    get_fs(prior_fit, prior_cov = matrix(1:6, 2, 3)),
    "square"
  )
  expect_error(
    get_fs(prior_fit, prior_cov = matrix(1.2, 1, 1)),
    "2 x 2"
  )
  pc_asym <- matrix(c(1, 0.3, 0.1, 1), 2, 2)
  expect_error(
    get_fs(prior_fit, prior_cov = pc_asym),
    "symmetric"
  )
  pc_npd <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(
    get_fs(prior_fit, prior_cov = pc_npd),
    "positive definite"
  )
  expect_error(
    get_fs(prior_fit, prior_cov = matrix(c(1, NA, NA, 1), 2, 2)),
    "finite"
  )
  expect_error(
    get_fs(prior_fit,
           prior_cov = matrix(c(1, 0.1, 0.1, 1), 2, 2,
                              dimnames = list(c("a", "b"), c("a", "b")))),
    "names must match"
  )
  # an unnamed q x q (q > 1) matrix is ambiguous (its row/column order might
  # not be the model's latent order), so it is rejected rather than fitted
  # silently in whatever order it happens to carry
  expect_error(
    get_fs(prior_fit, prior_cov = matrix(c(1.2, 0.25, 0.25, 0.8), 2, 2)),
    "must be a named matrix"
  )
})

test_that("priors are rejected for Bartlett scoring", {
  expect_error(
    get_fs(prior_fit, method = "Bartlett", prior_mean = pm),
    "only supported"
  )
  expect_error(
    get_fs(prior_fit, method = "ML", prior_cov = pc),
    "only supported"
  )
})

test_that("priors are rejected with reliability = TRUE", {
  expect_error(
    get_fs(prior_fit, reliability = TRUE, prior_mean = pm),
    "not supported"
  )
  expect_error(
    get_fs(prior_fit, reliability = TRUE, prior_cov = pc),
    "not supported"
  )
})

test_that("corrected_fsT = TRUE uses prior_cov", {
  fs_c <- get_fs(prior_fit, prior_cov = pc, corrected_fsT = TRUE,
                 format = "list")
  add <- correct_evfs(prior_fit, method = "regression",
                      psi_override = pc)[[1]]
  expect_equal(attr(fs_c, "fsT"), attr(fs_prior, "fsT") + add)
})

test_that("corrected_fsT = TRUE with prior_mean only uses model psi", {
  fs_cm <- get_fs(prior_fit, prior_mean = pm, corrected_fsT = TRUE,
                  format = "list")
  fs_pm <- get_fs(prior_fit, prior_mean = pm, format = "list")
  add <- correct_evfs(prior_fit, method = "regression")[[1]]
  expect_equal(attr(fs_cm, "fsT"), attr(fs_pm, "fsT") + add)
})

test_that("vfsLT = TRUE uses prior_cov", {
  fs_v <- get_fs(prior_fit, prior_cov = pc, vfsLT = TRUE,
                 format = "list")
  expect_equal(attr(fs_v, "vfsLT"),
               vcov_ld_evfs(prior_fit, method = "regression",
                            psi_override = pc),
               ignore_attr = TRUE)
})

test_that("data.frame entry point matches lavaan method entry point", {
  fs_df <- get_fs(PoliticalDemocracy, prior_model,
                  prior_mean = pm, prior_cov = pc,
                  format = "list")
  expect_equal(fs_df, fs_prior)
})

test_that("matrix entry point matches lavaan method entry point", {
  ind_cols <- c("x1", "x2", "x3", "y1", "y2", "y3", "y4")
  fs_mat <- get_fs(as.matrix(PoliticalDemocracy[, ind_cols]),
                   prior_model,
                   prior_mean = pm, prior_cov = pc,
                   format = "list")
  expect_equal(fs_mat, fs_prior)
})

test_that("get_fs_lavaan() forwards priors", {
  fs_leg <- get_fs_lavaan(prior_fit, prior_mean = pm, prior_cov = pc)
  expect_equal(fs_leg, fs_prior)
})

test_that("NULL priors preserve existing behavior", {
  hand <- hand_scores(psi = prior_est$psi, alpha = prior_est$alpha)
  expect_equal(fs_none[, c("fs_ind60", "fs_dem60")],
               as.data.frame(hand),
               ignore_attr = TRUE)
})

########## Single-factor (q = 1) example ##########

hs_model1 <- "visual =~ x1 + x2 + x3"
hs_fit1 <- cfa(hs_model1, data = HolzingerSwineford1939)

test_that("q = 1 accepts scalar, 1 x 1, and named-scalar prior_cov", {
  fs_a <- get_fs(hs_fit1, prior_cov = 1.2, format = "list")
  fs_b <- get_fs(hs_fit1, prior_cov = matrix(1.2, 1, 1),
                 format = "list")
  fs_c <- get_fs(hs_fit1, prior_cov = c(visual = 1.2),
                 format = "list")
  expect_equal(fs_a, fs_b)
  expect_equal(fs_a, fs_c)
  expect_equal(attr(fs_a, "fsT"), attr(fs_b, "fsT"))
  expect_equal(attr(fs_a, "fsL"), attr(fs_b, "fsL"))
})

########## Multi-group example ##########

hs_model3 <- "ind =~ x1 + x2 + x3
              dem =~ x4 + x5 + x6 + x7"
hs_data <- lavaan::HolzingerSwineford1939

pm_3 <- c(ind = 0.5, dem = -0.4)
pc_3 <- matrix(c(1.5, 0.3, 0.3, 0.9), 2, 2,
               dimnames = list(c("ind", "dem"), c("ind", "dem")))

hs_fit_g <- cfa(hs_model3, data = hs_data, group = "school")
fs_g_prior <- get_fs(hs_fit_g, prior_mean = pm_3, prior_cov = pc_3,
                     format = "list")
fs_g_unified <- get_fs(hs_fit_g, prior_mean = pm_3, prior_cov = pc_3,
                       format = "unified")

test_that("multi-group priors match manual compute_fscore() per group", {
  est_g <- lavInspect(hs_fit_g, what = "est")
  data_g <- lavInspect(hs_fit_g, what = "data")
  for (g in names(fs_g_prior)) {
    hand <- compute_fscore(data_g[[g]],
                           lambda = est_g[[g]]$lambda,
                           theta = est_g[[g]]$theta,
                           psi = pc_3,
                           nu = est_g[[g]]$nu,
                           alpha = pm_3,
                           method = "regression",
                           fs_matrices = TRUE)
    expect_equal(fs_g_prior[[g]][, c("fs_ind", "fs_dem")],
                 as.data.frame(hand),
                 ignore_attr = TRUE)
    expect_equal(attr(fs_g_prior[[g]], "fsT"), attr(hand, "fsT"))
    expect_equal(attr(fs_g_prior[[g]], "fsL"), attr(hand, "fsL"))
    expect_equal(attr(fs_g_prior[[g]], "fsb"), attr(hand, "fsb"))
  }
})

test_that("multi-group list and unified formats hold the same scores", {
  group_col <- names(fs_g_unified)[ncol(fs_g_unified)]
  for (g in names(fs_g_prior)) {
    expect_equal(
      unlist(fs_g_prior[[g]][c("fs_ind", "fs_dem")], use.names = FALSE),
      unlist(fs_g_unified[fs_g_unified[[group_col]] == g,
                          c("fs_ind", "fs_dem")], use.names = FALSE)
    )
  }
})

test_that("identical group data give identical prior-adjusted scores", {
  base <- hs_data[1:145, ]
  cols <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7")
  d2 <- rbind(base[cols], base[cols])
  d2$school <- factor(rep(c("GA", "GB"), each = 145))
  fit2 <- cfa(hs_model3, data = d2, group = "school")
  fs2 <- get_fs(fit2, prior_mean = pm_3, prior_cov = pc_3,
                format = "list")
  expect_equal(fs2[["GA"]][c("fs_ind", "fs_dem")],
               fs2[["GB"]][c("fs_ind", "fs_dem")],
               tolerance = 1e-8)
})

test_that("shared prior applies to all groups with missing data", {
  d_miss <- hs_data
  d_miss[1, "x4"] <- NA
  fit_miss <- cfa(hs_model3, data = d_miss, group = "school")
  est_miss <- lavInspect(fit_miss, what = "est")
  data_miss <- lavInspect(fit_miss, what = "data")
  blocks <- get_fs_blocks.lavaan(
    fit_miss,
    method = "regression",
    add_to_evfs = list(0, 0),
    prior_mean = pm_3,
    prior_cov = pc_3
  )
  for (g in names(blocks)) {
    est <- est_miss[[g]]
    grp_miss <- fit_miss@Data@Mp[[match(g, fit_miss@Data@group.label)]]
    pat <- if (is.null(grp_miss)) {
      setNames(rep(TRUE, nrow(est$lambda)), rownames(est$lambda))
    } else {
      grp_miss$pat
    }
    for (i in seq_along(blocks[[g]])) {
      blk <- blocks[[g]][[i]]
      idx <- blk$case_idx
      pat_i <- if (is.null(grp_miss)) pat else pat[i, ]
      y_i <- data_miss[[g]][idx, pat_i, drop = FALSE]
      hand <- compute_fscore(
        y_i,
        lambda = est$lambda[pat_i, , drop = FALSE],
        theta = est$theta[pat_i, pat_i, drop = FALSE],
        psi = pc_3,
        nu = est$nu[pat_i, , drop = FALSE],
        alpha = pm_3,
        method = "regression",
        fs_matrices = TRUE
      )
      expect_equal(blk$fs, hand, ignore_attr = TRUE)
      expect_equal(blk$fsT, attr(hand, "fsT"))
      expect_equal(blk$fsL, attr(hand, "fsL"))
      expect_equal(blk$fsb, attr(hand, "fsb"))
    }
  }
})

########## Mermod rejection ##########

lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

test_that("merMod objects reject prior_mean / prior_cov", {
  expect_error(
    get_fs(lme1, prior_mean = c(Reaction = 1)),
    "not supported"
  )
  expect_error(
    get_fs(lme1, prior_cov = 1),
    "not supported"
  )
  expect_error(
    get_fs_lmer(lme1, prior_mean = c(Reaction = 1)),
    "not supported"
  )
  expect_error(
    get_fs_lmer(lme1, prior_cov = 1),
    "not supported"
  )
})
