# Loading packages and functions
library(lavaan)
library(lme4)
if (requireNamespace("umx", quietly = TRUE)) {
  library(umx)
}

########## Single-group example ##########

# Prepare test objects
single_model <- '
                 # latent variables
                   ind60 =~ x1 + x2 + x3
                   dem60 =~ y1 + y2 + y3 + y4
                '

# model = cfa_3var
test_object_fs <- get_fs(PoliticalDemocracy, single_model, format = "list")

########## Testing section ############

# Class of input

test_that("test the model input", {
  expect_type(single_model, "character")
})

test_that("test the data input", {
  expect_s3_class(PoliticalDemocracy, "data.frame")
})

test_that("matrix input is converted and re-dispatched identically", {
  fs_mat <- get_fs(as.matrix(PoliticalDemocracy), single_model,
                  format = "list")
  expect_identical(fs_mat, test_object_fs)
})

# Class of output

test_that("Gives an output of data frame", {
  expect_s3_class(test_object_fs, "data.frame")
})

ncol_cfa <- function(x) {
  nfac <- nrow(lavInspect(x, what = "cor.lv"))
  # fs + se + loadings + ev
  nfac * 2 + nfac ^ 2 + nfac * (nfac + 1) / 2
}

test_that("Test the number of factors is equal", {
  expected_cols <- ncol_cfa(
    cfa(model = single_model, data = PoliticalDemocracy)
  )
  expect_equal(length(test_object_fs), expected_cols)
})

test_that("Number of rows is the same as the original data", {
  expect_identical(nrow(test_object_fs), nrow(PoliticalDemocracy))
})

# Test standard error

test_that("Standard errors for each observation are the same", {
  fs_names <- colnames(test_object_fs)
  names_se <- grep("_se", fs_names, value = TRUE)
  for (j in names_se) {
    expect_identical(var(test_object_fs[[j]]), 0)
  }
})

test_that(
  "Standard errors for each observation are positive numbers and within 1",
  {
    fs_names <- colnames(test_object_fs)
    names_se <- grep("_se", fs_names, value = TRUE)
    for (j in names_se) {
      expect_gt(min(test_object_fs[[j]]), 0)
      expect_lt(max(test_object_fs[[j]]), 1)
    }
  }
)

# Bartlett scores
test_object_fs_bar <- get_fs(
  PoliticalDemocracy, single_model,
  method = "Bartlett",
  format = "list"
)

test_that("Standard errors for each observation are the same", {
  fs_names <- colnames(test_object_fs_bar)
  names_se <- grep("_se", fs_names, value = TRUE)
  for (j in names_se) {
    expect_identical(var(test_object_fs_bar[[j]]), 0)
  }
})

test_that(
  "Standard errors for each observation are positive numbers and within 1",
  {
    fs_names <- colnames(test_object_fs)
    names_se <- grep("_se", fs_names, value = TRUE)
    for (j in names_se) {
      expect_gt(min(test_object_fs_bar[[j]]), 0)
      expect_lt(max(test_object_fs_bar[[j]]), 1)
    }
  }
)

# Method aliases
test_that("method = 'EB' is an alias for 'regression'", {
  fs_eb <- get_fs(PoliticalDemocracy, single_model,
    method = "EB",
    format = "list"
  )
  expect_equal(fs_eb, test_object_fs, ignore_attr = TRUE)
})

test_that("method = 'ML' is an alias for 'Bartlett'", {
  fs_ml <- get_fs(PoliticalDemocracy, single_model,
    method = "ML",
    format = "list"
  )
  expect_equal(fs_ml, test_object_fs_bar, ignore_attr = TRUE)
})

test_that("aliases also work on a fitted lavaan object", {
  fit <- cfa(single_model, data = PoliticalDemocracy)
  fs_ml <- get_fs(fit, method = "ML")
  fs_bart <- get_fs(fit, method = "Bartlett")
  expect_equal(fs_ml, fs_bart, ignore_attr = TRUE)
  fs_eb <- get_fs(fit, method = "EB")
  fs_reg <- get_fs_lavaan(fit)
  expect_equal(as.data.frame(fs_eb), fs_reg, ignore_attr = TRUE)
})

test_that("invalid method errors with the full choice set", {
  expect_error(
    get_fs(PoliticalDemocracy, single_model, method = "MLR"),
    "should be one of"
  )
})


########## multi-group examples ##########

###### One-factor example #####

# Prepare for test objects
hs_model <- 'visual  =~ x1 + x2 + x3'
# multi_fit <- cfa(hs_model,
#                  data = HolzingerSwineford1939,
#                  group = "school")
#get_fs(HolzingerSwineford1939, hs_model, group = "school")
test_object_fs_multi <- get_fs(
  HolzingerSwineford1939[c("school", "x1", "x2", "x3")],
  hs_model,
  group = "school",
  format = "list"
)

########## Testing section ############

# Class of output

test_that("Number of factors is equal", {
  expected_cols <- ncol_cfa(
    cfa(model = hs_model, data = PoliticalDemocracy)
  )
  # Add 1 for group
  expect_equal(length(test_object_fs_multi[[1]]), expected_cols + 1)
})

test_that("Number of rows is the same as the original data", {
  expect_identical(nrow(test_object_fs_multi[[1]]) +
                     nrow(test_object_fs_multi[[2]]),
                   nrow(HolzingerSwineford1939))
})

# Test standard errors

test_that(
  "Standard errors for each observation are the same within groups",
  {
    test_se <- vapply(test_object_fs_multi,
      FUN = function(x) var(x$fs_visual_se),
      FUN.VALUE = numeric(1)
    )
    for (i in length(test_se)) {
      expect_identical(unname(test_se[i]), 0)
    }
  }
)

test_that(
  "Standard errors for each observation are positive numbers and within 1",
  {
    expect_gt(min(test_object_fs_multi[[1]]$fs_visual_se), 0)
    expect_lt(max(test_object_fs_multi[[1]]$fs_visual_se), 1)
    expect_gt(min(test_object_fs_multi[[2]]$fs_visual_se), 0)
    expect_lt(max(test_object_fs_multi[[2]]$fs_visual_se), 1)
  }
)

# Bartlett scores
test_object_fs_multi_bar <- get_fs(
  HolzingerSwineford1939[c("school", "x1", "x2", "x3")],
  hs_model,
  group = "school",
  method = "Bartlett",
  format = "list"
)

test_that(
  "Standard errors for each observation are the same within groups",
  {
    test_se <- vapply(test_object_fs_multi_bar,
      FUN = function(x) var(x$fs_visual_se),
      FUN.VALUE = numeric(1)
    )
    for (i in length(test_se)) {
      expect_identical(unname(test_se[i]), 0)
    }
  }
)

test_that(
  "Standard errors for each observation are positive numbers and within 1",
  {
    expect_gt(min(test_object_fs_multi_bar[[1]]$fs_visual_se), 0)
    expect_lt(max(test_object_fs_multi_bar[[1]]$fs_visual_se), 1)
    expect_gt(min(test_object_fs_multi_bar[[2]]$fs_visual_se), 0)
    expect_lt(max(test_object_fs_multi_bar[[2]]$fs_visual_se), 1)
  }
)

###### Multiple factors example #####

# Prepare for test objects
hs_model_2 <- ' visual =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed =~ x7 + x8 + x9 '
test_object_fs_multi_2 <- get_fs(HolzingerSwineford1939,
                                  hs_model_2,
                                  group = "school",
                                  format = "list")

########## Testing section ############

# Class of output

test_that("Number of factors is equal", {
  expected_cols <- ncol_cfa(
    cfa(model = hs_model_2, data = HolzingerSwineford1939)
  )
  expect_equal(length(test_object_fs_multi_2[[1]]), expected_cols + 1)
})

test_that("Number of rows is the same as the original data", {
  expect_identical(nrow(test_object_fs_multi_2[[1]]) +
                     nrow(test_object_fs_multi_2[[2]]),
                   nrow(HolzingerSwineford1939))
})

# Test standard errors

test_that(
  "Standard errors for each observation are the same within groups",
  {
    fs_names <- colnames(test_object_fs_multi_2[[1]])
    names_se <- grep("_se", fs_names, value = TRUE)
    for (i in names_se) {
      test_se <- vapply(test_object_fs_multi_2,
        FUN = function(x) var(x[[i]]),
        FUN.VALUE = numeric(1)
      )
      expect_identical(max(test_se), 0)
    }
  }
)

test_that(
  "Standard errors for each observation are positive numbers and within 1",
  {
    rb_test_object <- do.call(rbind, test_object_fs_multi_2)
    fs_names <- colnames(rb_test_object)
    names_se <- grep("_se", fs_names, value = TRUE)
    for (i in names_se) {
      expect_gt(min(rb_test_object[, i]), 0)
      expect_lt(max(rb_test_object[, i]), 1)
    }
  }
)

# Bartlett scores
test_object_fs_multi_2_bar <- get_fs(HolzingerSwineford1939,
                                      hs_model_2,
                                      group = "school",
                                      format = "list")

test_that(
  "Standard errors for each observation are the same within groups",
  {
    rb_test_object <- do.call(rbind, test_object_fs_multi_2_bar)
    fs_names <- colnames(rb_test_object)
    names_se <- grep("_se", fs_names, value = TRUE)
    for (i in names_se) {
      test_se <- tapply(
        rb_test_object[[i]],
        rb_test_object$school, var
      )
      expect_identical(max(test_se), 0)
    }
  }
)

test_that(
  "Standard errors for each observation are positive numbers and within 1",
  {
    rb_test_object <- do.call(rbind, test_object_fs_multi_2_bar)
    fs_names <- colnames(rb_test_object)
    names_se <- grep("_se", fs_names, value = TRUE)
    for (i in names_se) {
      expect_gt(min(rb_test_object[, i]), 0)
      expect_lt(max(rb_test_object[, i]), 1)
    }
  }
)

fs_config <- get_fs(HolzingerSwineford1939,
                    hs_model_2,
                    group = "school",
                    corrected_fsT = TRUE,
                    format = "list"
)

fs_metric <- get_fs(HolzingerSwineford1939,
                    hs_model_2,
                    group = "school",
                    group.equal = "loadings",
                    corrected_fsT = TRUE,
                    format = "list"
)

fs_single <- get_fs(
  HolzingerSwineford1939 |>
    subset(school == "Grant-White"),
  hs_model_2,
  corrected_fsT = TRUE,
  format = "list"
)

test_that("Correction factor is similar with single or multiple groups", {
  fst1 <- attr(fs_config, "fsT")
  fst2 <- attr(fs_metric, "fsT")
  fst3 <- attr(fs_single, "fsT")
  expect_equal(fst1[[2]], fst3, tolerance = 1e-5)
  d1 <- fst1[[1]] - fst1[[2]]
  d2 <- fst2[[1]] - fst2[[2]]
  expect_lt(mean(abs(d2)), mean(abs(d1)))
})

test_that("Same factor scores as `lme4::ranef()`", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  expect_equal(as.data.frame(get_fs_lmer(lme1)[, 1:2]),
               ranef(lme1)[[1]],
               ignore_attr = TRUE)
})

test_that("get_fs() S3 dispatch on merMod object", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  fs_direct <- get_fs(lme1)
  fs_wrapper <- get_fs_lmer(lme1)
  expect_equal(as.data.frame(fs_direct[, 1:2]),
               as.data.frame(fs_wrapper[, 1:2]),
               ignore_attr = TRUE)
})

test_that("get_fs_blocks.merMod() produces correct block structure", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  blocks <- R2spa:::get_fs_blocks.merMod(lme1)
  n_clusters <- nlevels(sleepstudy$Subject)
  expect_equal(length(blocks), n_clusters)
  for (b in blocks) {
    expect_true(is.list(b))
    expect_true(all(c("case_idx", "fs", "fsL", "fsT") %in% names(b)))
    expect_true(length(b$case_idx) > 0)
    expect_equal(nrow(b$fs), 1)
  }
})

test_that("get_fs.merMod() default uses fs_u0-style names (no _eb)", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  fs <- get_fs(lme1)
  expect_equal(
    colnames(fs),
    c("fs_u0", "fs_u1",
      "fs_u0_se", "fs_u1_se",
      "u0_by_fs_u0", "u0_by_fs_u1", "u1_by_fs_u0", "u1_by_fs_u1",
      "ev_fs_u0", "ecov_fs_u1_fs_u0", "ev_fs_u1")
  )
})

test_that("get_fs.merMod(legacy_names = TRUE) reproduces u0_eb-style names", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  fs <- get_fs(lme1, legacy_names = TRUE)
  expect_equal(
    colnames(fs),
    c("u0_eb", "u1_eb",
      "u0_eb_se", "u1_eb_se",
      "u0_by_u0_eb", "u0_by_u1_eb", "u1_by_u0_eb", "u1_by_u1_eb",
      "ev_u0_eb", "ecov_u0_eb_u1_eb", "ev_u1_eb")
  )
})

test_that("get_fs_lmer() defaults to legacy u0_eb-style names", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  fs <- get_fs_lmer(lme1)
  expect_equal(
    colnames(fs),
    c("u0_eb", "u1_eb",
      "u0_eb_se", "u1_eb_se",
      "u0_by_u0_eb", "u0_by_u1_eb", "u1_by_u0_eb", "u1_by_u1_eb",
      "ev_u0_eb", "ecov_u0_eb_u1_eb", "ev_u1_eb")
  )
})

test_that("merMod scoring_matrix reproduces EB scores (score identity)", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  sm <- attr(get_fs(lme1), "scoring_matrix")
  mf <- model.frame(lme1)
  y <- model.response(mf)
  Xfull <- as.matrix(lme1@pp$X)
  beta <- fixef(lme1)
  flist <- lme1@flist[[1]]
  case_idx <- split(seq_len(nrow(mf)), flist)
  for (j in seq_along(case_idx)) {
    lv <- levels(flist)[j]
    idx <- case_idx[[j]]
    rec <- t(sm[[lv]] %*% (y[idx] - as.numeric(Xfull[idx, , drop = FALSE] %*% beta)))
    expect_equal(
      as.numeric(rec),
      as.numeric(ranef(lme1)[[1]][j, ]),
      tolerance = 1e-10
    )
  }
})

test_that("merMod scoring_matrix has correct structure (per-cluster list)", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  sm <- attr(get_fs(lme1), "scoring_matrix")
  expect_named(sm, as.character(levels(sleepstudy$Subject)))
  expect_length(sm, 18)
  for (m in sm) {
    expect_equal(dim(m), c(2L, 10L))
    expect_equal(rownames(m), c("fs_u0", "fs_u1"))
  }
})

test_that("get_fs.merMod() works when random and fixed designs differ (Z != X)", {
  lme2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  fs <- get_fs(lme2)
  sm <- attr(fs, "scoring_matrix")
  mf <- model.frame(lme2)
  y <- model.response(mf)
  Xfull <- as.matrix(lme2@pp$X)
  beta <- fixef(lme2)
  flist <- lme2@flist[[1]]
  case_idx <- split(seq_len(nrow(mf)), flist)
  expect_length(sm, nlevels(flist))
  for (j in seq_along(case_idx)) {
    lv <- levels(flist)[j]
    idx <- case_idx[[j]]
    expect_equal(dim(sm[[lv]]), c(1L, 10L))
    rec <- t(sm[[lv]] %*% (y[idx] - as.numeric(Xfull[idx, , drop = FALSE] %*% beta)))
    expect_equal(
      as.numeric(rec),
      as.numeric(ranef(lme2)[[1]][j, 1]),
      tolerance = 1e-10
    )
  }
})

test_that("merMod scoring_matrix works with unbalanced clusters", {
  # drop one row from each of two different subjects: 9 vs 10 obs/cluster
  ss <- sleepstudy[-c(1, 111), ]
  lmu <- lmer(Reaction ~ Days + (Days | Subject), ss)
  sm <- attr(get_fs(lmu), "scoring_matrix")
  mf <- model.frame(lmu)
  y <- model.response(mf)
  Xfull <- as.matrix(lmu@pp$X)
  beta <- fixef(lmu)
  flist <- lmu@flist[[1]]
  case_idx <- split(seq_len(nrow(mf)), flist)
  expect_length(sm, nlevels(flist))
  sizes <- vapply(sm, ncol, integer(1))
  expect_setequal(unique(sizes), c(9L, 10L))
  for (j in seq_along(case_idx)) {
    lv <- levels(flist)[j]
    idx <- case_idx[[j]]
    expect_equal(ncol(sm[[lv]]), length(idx))
    rec <- t(sm[[lv]] %*% (y[idx] - as.numeric(Xfull[idx, , drop = FALSE] %*% beta)))
    expect_equal(
      as.numeric(rec),
      as.numeric(ranef(lmu)[[1]][j, ]),
      tolerance = 1e-10
    )
  }
})

test_that("get_fs_lmer() carries the scoring_matrix attribute", {
  lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  sm_new <- attr(get_fs(lme1), "scoring_matrix")
  sm_leg <- attr(get_fs_lmer(lme1), "scoring_matrix")
  expect_false(is.null(sm_leg))
  expect_named(sm_leg, names(sm_new))
  for (j in seq_along(sm_new)) {
    expect_equal(sm_leg[[j]], sm_new[[j]], ignore_attr = TRUE)
  }
})

test_that("get_fs.merMod() is robust to row order (shuffled data)", {
  set.seed(42)
  sh <- sleepstudy[sample(nrow(sleepstudy)), ]
  lms <- lmer(Reaction ~ Days + (1 | Subject), sh)
  fs <- get_fs(lms)
  sm <- attr(fs, "scoring_matrix")
  # Cluster names follow the factor's level (canonical) order, not the
  # first-appearance order induced by the shuffled row order.
  expect_named(sm, as.character(sort(unique(sh$Subject))))
  expect_equal(names(sm), dimnames(attr(fs, "fsT"))[[3]])
  expect_equal(names(sm), dimnames(attr(fs, "fsL"))[[3]])
  # Scores (one row per cluster, in that same order) align with ranef.
  rn <- as.data.frame(ranef(lms)[[1]])
  score_col <- grep("^fs_u0$", colnames(fs))[1]
  expect_equal(as.numeric(fs[[score_col]]), as.numeric(rn[[1]]),
               tolerance = 1e-10)
})

test_that("get_fs.merMod() follows user-specified factor level order", {
  d <- sleepstudy
  d$SubF <- factor(d$Subject, levels = rev(sort(unique(d$Subject))))
  fm <- lmer(Reaction ~ Days + (1 | SubF), d)
  fs <- get_fs(fm)
  sm <- attr(fs, "scoring_matrix")
  expect_named(sm, levels(d$SubF))
  rn <- as.data.frame(ranef(fm)[[1]])
  expect_equal(rownames(rn), levels(d$SubF))
  score_col <- grep("^fs_u0$", colnames(fs))[1]
  expect_equal(as.numeric(fs[[score_col]]), as.numeric(rn[[1]]),
               tolerance = 1e-10)
})

test_that("get_fs.merMod() is robust to non-monotonic numeric cluster ids", {
  d <- sleepstudy
  sids <- unique(d$Subject)
  idmap <- c(9L, 4L, 7L, 2L, 8L, 1L, 11L, 3L, 19L, 5L, 17L, 6L,
             21L, 10L, 13L, 15L, 16L, 20L)
  d$nid <- idmap[match(d$Subject, sids)]
  fn <- lmer(Reaction ~ Days + (1 | nid), d)
  fs <- get_fs(fn)
  sm <- attr(fs, "scoring_matrix")
  # Names are in lme4's canonical (sorted) numeric level order, not in
  # data-appearance order (9, 4, 7, ...).
  expect_equal(names(sm), rownames(as.data.frame(ranef(fn)[[1]])))
  # canonical order starts at "1"; appearance order would start at "9"
  expect_equal(names(sm)[1], "1")
  score_col <- grep("^fs_u0$", colnames(fs))[1]
  rn <- as.data.frame(ranef(fn)[[1]])
  expect_equal(as.numeric(fs[[score_col]]), as.numeric(rn[[1]]),
               tolerance = 1e-10)
})
########## Computing reliability ##########

test_that("Reliability of regression factor scores", {
  fs <- get_fs(PoliticalDemocracy[c("x1", "x2", "x3")],
               corrected_fsT = TRUE, reliability = TRUE, std.lv = TRUE)
  expect_equal(attr(fs, "reliability"), .9607411,
               tolerance = 1e-7)
})

test_that("Reliability of Bartlett factor scores", {
  fs <- get_fs(PoliticalDemocracy[c("x1", "x2", "x3")],
               corrected_fsT = TRUE, reliability = TRUE, std.lv = TRUE,
               method = "Bartlett")
  expect_equal(attr(fs, "reliability"), .9607457,
               tolerance = 1e-7)
})

test_that("Reliability with non-diagonal theta", {
  mod1 <- "visual =~ x1 + x2 + x3 + x9"
  fit1 <- cfa(model = mod1, data = HolzingerSwineford1939)
  fit2 <- cfa(model = c(mod1, "x1 ~~ x9"), data = HolzingerSwineford1939)
  expect_gt(compute_fsrel(fit1)[[1]], compute_fsrel(fit2)[[1]])
})

# test_that("Reliability of regression fs > reliability of Bartlett fs", {
#   rel_reg <- get_fs(PoliticalDemocracy[c("x1", "x2", "x3")],
#                     corrected_fsT = TRUE, reliability = TRUE, std.lv = TRUE)
#   rel_bart <- get_fs(PoliticalDemocracy[c("x1", "x2", "x3")],
#                      corrected_fsT = TRUE, reliability = TRUE, std.lv = TRUE,
#                      method = "Bartlett")
#   expect_gt(attr(rel_reg, "reliability"), attr(rel_bart, "reliability"))
# })



test_that("augment_lav_predict() works for complete data",
  code = {
    cfa_single <- cfa(model = hs_model_2,
                      data = HolzingerSwineford1939)
    a0 <- augment_lav_predict(cfa_single, drop_list_single = FALSE)
    expect_type(a0, type = "list")
    a1 <- augment_lav_predict(cfa_single)
    expect_equal(as.matrix(a1[, 1:3]),
                 as.matrix(lavPredict(cfa_single)),
                 ignore_attr = TRUE)
    a2 <- augment_lav_predict(cfa_single, method = "Bartlett")
    expect_equal(a2,
                 get_fs_lavaan(cfa_single, method = "Bartlett"),
                 ignore_attr = TRUE)
    # Mean structure
    cfa_multi <- cfa(model = hs_model_2,
                     data = HolzingerSwineford1939[sample(301), ],
                     group = "school")
    a3 <- augment_lav_predict(cfa_multi)
    fs_multi <- get_fs_lavaan(cfa_multi)
    expect_equal(a3[[1]][, 1:21],
                 fs_multi[[1]][, 1:21],
                 ignore_attr = TRUE)
    expect_equal(a3[[2]][, 1:21],
                 fs_multi[[2]][, 1:21],
                 ignore_attr = TRUE)
  }
)

hs <- HolzingerSwineford1939
# introduce missing data
set.seed(1334)
hs[!rbinom(301, size = 1, prob = 0.7), 7] <- NA
hs[!rbinom(301, size = 1, prob = 0.7), 8] <- NA
hs[!rbinom(301, size = 1, prob = 0.7), 9] <- NA
# Mean structure
cfa_multi <- cfa(model = hs_model_2,
                 data = hs[sample(301), ],
                 group = "school",
                 missing = "fiml")
a3 <- augment_lav_predict(cfa_multi)

cfa_fiml <- cfa(
  model = hs_model_2, data = hs, missing = "fiml",
  estimator = "MLR"
)
a2 <- augment_lav_predict(cfa_fiml, method = "Bartlett")

test_that("augment_lav_predict() works for missing data",
  code = {
    cfa_lw <- cfa(model = hs_model_2, data = hs, bounds = TRUE,
                  meanstructure = TRUE)
    expect_no_error(augment_lav_predict(cfa_lw))
    # NA for cases with x1, x2, x3 all missing
    expect_true(all(is.na(
      a2$fs_visual[which(rowSums(!is.na(hs[, 7:9])) == 0)]))
    )
    fs_multi <- get_fs_lavaan(cfa_multi)
    expect_equal(a3[[1]][, 1:21],
                 fs_multi[[1]][, 1:21],
                 ignore_attr = TRUE)
    expect_equal(a3[[2]][, 1:21],
                 fs_multi[[2]][, 1:21],
                 ignore_attr = TRUE)
  }
)

test_that("SE paths error explicitly on multi-pattern missing data", {
  # cfa_fiml (defined above) contains missing values on x1:x3, so lavaan
  # scores cases on multiple observed-indicator patterns.
  expect_error(get_fs(cfa_fiml, corrected_fsT = TRUE),
               "not supported when the data contain missing values")
  expect_error(get_fs(cfa_fiml, reliability = TRUE),
               "not supported when the data contain missing values")
  expect_error(get_fs(cfa_fiml, vfsLT = TRUE),
               "not supported when the data contain missing values")
})
