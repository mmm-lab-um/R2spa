# Loading packages and functions
library(lavaan)
library(lme4)
library(umx)

########## Single-group example ##########

# Prepare test objects
single_model <- '
                 # latent variables
                   ind60 =~ x1 + x2 + x3
                   dem60 =~ y1 + y2 + y3 + y4
                '

# model = cfa_3var
test_object_fs <- get_fs(PoliticalDemocracy, single_model)

########## Testing section ############

# Class of input

test_that("test the model input", {
  expect_type(single_model, "character")
})

test_that("test the data input", {
  expect_s3_class(PoliticalDemocracy, "data.frame")
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
  method = "Bartlett"
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
  group = "school"
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
  method = "Bartlett"
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
                                 group = "school")

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
                                     group = "school")

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
                    corrected_fsT = TRUE
)

fs_metric <- get_fs(HolzingerSwineford1939,
                    hs_model_2,
                    group = "school",
                    group.equal = "loadings",
                    corrected_fsT = TRUE
)

fs_single <- get_fs(
  HolzingerSwineford1939 |>
    subset(school == "Grant-White"),
  hs_model_2,
  corrected_fsT = TRUE
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

lcov_umx <- umxLav2RAM(
  "
    visual ~~ textual + speed
    textual ~~ speed
    visual + textual + speed ~ 1
  ",
  printTab = FALSE
)
tspab_mx <- tspa_mx_model(lcov_umx,
  data = a2,
  mat_ld = attr(a2, which = "ld"),
  mat_ev = attr(a2, which = "ev")
)
# Run OpenMx
tspab_mx_fit <- mxRun(tspab_mx)

test_that("tspa_mx() gives similar results as lavaan with missing data", 
  code = {
    expect_equal(tspab_mx_fit$m1$S$values,
                 expected = lavInspect(cfa_fiml, what = "est")$psi,
                 tolerance = 1e-5,
                 ignore_attr = TRUE)
  }
)
