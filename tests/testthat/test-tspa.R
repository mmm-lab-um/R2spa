########################### Test 2S-PA function ##############################

# Loading packages and functions
library(lavaan)
if (requireNamespace("umx", quietly = TRUE)) {
  library(umx)
}

########## Single-group example ##########

# Prepare test objects

# Example 1: Single-group with two variables

# CFA model
cfa_single1 <- "
# latent variables
ind60 =~ x1 + x2 + x3
"
cfa_single2 <- "
# latent variables
dem60 =~ y1 + y2 + y3 + y4
"

# get factor scores
fs_single1 <- get_fs(PoliticalDemocracy, cfa_single1, format = "list")
fs_single2 <- get_fs(PoliticalDemocracy, cfa_single2, format = "list")
fs_dat_single <- cbind(fs_single1, fs_single2)

cfa_model_single <- '
                                 # latent variables (indicated by factor scores)
                                   ind60 =~ 1 * fs_ind60
                                   dem60 =~ 1 * fs_dem60
                                 # constrain the errors
                                   fs_ind60 ~~ 0.1213615^2 * fs_ind60
                                   fs_dem60 ~~ 0.6756472^2 * fs_dem60
                                 # latent variances
                                   ind60 ~~ v1 * ind60
                                   dem60 ~~ v2 * dem60
                                 # regressions
                                   dem60 ~ ind60
                              '
cfa_single <-
  sem(model = cfa_model_single, data  = fs_dat_single)

# tspa model
tspa_single <-
  tspa(
    model = "dem60 ~ ind60",
    data = fs_dat_single,
    se_fs = c(ind60 = 0.1213615, dem60 = 0.6756472)
  )

########## Testing section ############

# Class of input
var_len <- 2
se <- c(ind60 = 0.1213615, dem60 = 0.6756472)

# The tspa data should be composed of two parts: variable, and se
test_that(
  "Number of columns in tspa data are multiples of the variable length",
  {
    expect_gt(ncol(fs_dat_single), 1)
    expect_equal(ncol(fs_dat_single) %% var_len, 0)
  }
)

test_that("Test the data variable names should contain prefix (fs_)", {
  fs_names <- colnames(fs_dat_single)
  expect_true(all(grepl("fs_", fs_names)))
})

# Class of output

# Parameter estimates

test_that(
  "Regression coefficients of factors are the same for two methods",
  {
    expect_equal(
      coef(cfa_single)["dem60~ind60"],
      coef(tspa_single)["dem60~ind60"]
    )
  }
)

test_that(
  "se of regression coefficients are the same for two methods",
  {
    expect_equal(
      vcov(cfa_single)[
        c("dem60~ind60", "v1", "v2"),
        c("dem60~ind60", "v1", "v2")
      ],
      vcov(tspa_single)[
        c("dem60~ind60", "ind60~~ind60", "dem60~~dem60"),
        c("dem60~ind60", "ind60~~ind60", "dem60~~dem60")
      ],
      ignore_attr = TRUE
    )
  }
)

# Fit measures

test_that("test if fit indices are the same for two methods", {
  expect_equal(fitmeasures(cfa_single), fitmeasures(tspa_single))
  # We can add more comparisons of fitting measures by changing the name
})

# Example 2: Single group with three variables

# CFA model
cfa_3var1 <- '
                            # latent variables
                            ind60 =~ x1 + x2 + x3
                           '
cfa_3var2 <- '
                            # latent variables
                            dem60 =~ y1 + y2 + y3 + y4
                           '
cfa_3var3 <- '
                            # latent variables
                            dem65 =~ y5 + y6 + y7 + y8
                           '

# get factor scores
fs_3var1 <- get_fs(PoliticalDemocracy, cfa_3var1, format = "list")
fs_3var2 <- get_fs(PoliticalDemocracy, cfa_3var2, format = "list")
fs_3var3 <- get_fs(PoliticalDemocracy, cfa_3var3, format = "list")
fs_dat_3var <- cbind(fs_3var1, fs_3var2, fs_3var3)

sem_model_3var <- '
                           # latent variables (indicated by factor scores)
                             ind60 =~ x1 + x2 + x3
                             dem60 =~ y1 + y2 + y3 + y4
                             dem65 =~ y5 + y6 + y7 + y8
                           # regressions
                             dem60 ~ ind60
                             dem65 ~ ind60 + dem60
                      '

sem_3var <- sem(model = sem_model_3var, data  = PoliticalDemocracy)

# tspa model
tspa_3var <- tspa(
  model = "dem60 ~ ind60
               dem65 ~ ind60 + dem60",
  data = fs_dat_3var,
  se_fs = c(
    ind60 = 0.1213615,
    dem60 = 0.6756472,
    dem65 = 0.5724405
  )
)

########## Testing section #############

# Standardized parameter estimates
sem_path_3var <- subset(standardizedSolution(sem_3var),
                        subset = op == "~")
tspa_path_3var <- subset(standardizedSolution(tspa_3var),
                         subset = op == "~")

test_that(
  "Regression coefficients of factors are similar for two methods",
  {
    expect_lt(
      max(abs(sem_path_3var$est.std - tspa_path_3var$est.std)),
      expected = .05
    )
  }
)

test_that(
  "se of regression coefficients are similar for two methods",
  {
    expect_lt(
      max(abs(sem_path_3var$se - tspa_path_3var$se)),
      expected = .01
    )
  }
)

# Variance of factors
sem_var_3var <- subset(standardizedSolution(sem_3var),
                       subset = op == "~~" &
                         lhs %in% c("ind60", "dem60", "dem65"))
tspa_var_3var <- subset(standardizedSolution(tspa_3var),
                        subset = op == "~~" &
                          lhs %in% c("ind60", "dem60", "dem65"))

test_that("test if the variance of factor is similar for two methods", {
  expect_lt(
    max(abs(sem_var_3var$est.std - tspa_var_3var$est.std)),
    expected = .05
  )
})

test_that("test if the se of variance is similar for two methods", {
  expect_lt(
    max(abs(sem_var_3var$se - tspa_var_3var$se)),
    expected = .01
  )
})

########## Multi-group example ##########

# get factor scores
fs_dat_visual <- get_fs(HolzingerSwineford1939,
                        model = "visual =~ x1 + x2 + x3",
                        group = "school",
                        format = "list")
fs_dat_speed <- get_fs(HolzingerSwineford1939,
                       model = "speed =~ x7 + x8 + x9",
                       group = "school",
                       format = "list")
fs_dat_multi <- cbind(
  do.call(rbind, fs_dat_visual),
  do.call(rbind, fs_dat_speed)
)

# SEM model
sem_model_multi <- '
 # latent variables (indicated by factor scores)
   visual =~ c(1, 1) * fs_visual
   speed =~ c(1, 1) * fs_speed
 # constrain the errors
   fs_visual ~~ c(0.11501092038276, 0.097236701584) * fs_visual
   fs_speed ~~ c(0.07766672265625, 0.07510378617049) * fs_speed
 # latent variances
   visual ~~ c(v11, v12) * visual
   speed ~~ c(v21, v22) * speed
 # regressions
   visual ~ speed
'

sem_multi <-
  sem(model = sem_model_multi,
      data  = fs_dat_multi,
      group = "school")

# tspa model
tspa_multi <- tspa(
  model = "visual ~ speed",
  data = fs_dat_multi,
  se_fs = data.frame(
    visual = c(0.3391326, 0.3118280),
    speed = c(0.2786875, 0.2740507)
  ),
  group = "school"
  # group.equal = "regressions"
)

# the same single-factor fit with the data passed as a list of per-group
# data frames (the get_fs(format = "list") results cbind()ed per group);
# tspa() coerces the list to a data frame before the stage-2 fit
fs_dat_multi_list <- Map(
  function(a, b) cbind(a, b), fs_dat_visual, fs_dat_speed
)
tspa_multi_list <- tspa(
  model = "visual ~ speed",
  data = fs_dat_multi_list,
  se_fs = data.frame(
    visual = c(0.3391326, 0.3118280),
    speed = c(0.2786875, 0.2740507)
  ),
  group = "school"
)

########## Testing section #############

# Standardized parameter estimates
sem_path_multi <- subset(standardizedSolution(sem_multi),
                         subset = op == "~")
tspa_path_multi <- subset(standardizedSolution(tspa_multi),
                          subset = op == "~")

test_that(
  "Regression coefficients of factors are similar for two methods",
  {
    expect_equal(
      sem_path_multi$est.std,
      tspa_path_multi$est.std
    )
  }
)

test_that("se of regression coefficients are similar for two methods", {
  expect_equal(
    sem_path_multi$se,
    tspa_path_multi$se
  )
})

# Variance of factors

sem_var_multi <- subset(standardizedSolution(sem_multi),
                        subset = op == "~~" &
                          lhs %in% c("ind60", "dem60", "dem65"))
tspa_var_multi <- subset(standardizedSolution(tspa_multi),
                         subset = op == "~~" &
                           lhs %in% c("ind60", "dem60", "dem65"))

test_that("test if the variance of factor is similar for two methods", {
  expect_equal(
    sem_var_multi$est.std,
    tspa_var_multi$est.std
  )
})

test_that("test if the se of variance is similar for two methods", {
  expect_equal(
    sem_var_multi$se,
    tspa_var_multi$se
  )
})

test_that("tspa(): single-factor list-of-frames data matches the rbind'd-frame fit", {
  expect_equal(
    lavInspect(tspa_multi_list, "est"),
    lavInspect(tspa_multi, "est")
  )
})

# Test tspa_mf()
mod4 <- "
  # latent variables
    visual =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed =~ x7 + x8 + x9

"
fs_dat4 <- get_fs(HolzingerSwineford1939, model = mod4, std.lv = TRUE,
                  group = "school", format = "list")
tspa_mod_m <- tspa_mf(
  model = "visual ~ speed
           textual ~ visual + speed",
  data = fs_dat4,
  fsT = attr(fs_dat4, "fsT"),
  fsL = attr(fs_dat4, "fsL"),
  fsb = NULL
)

factors_order_m <- subset(lavaan::lavaanify(tspa_mod_m, ngroup = 2),
                          op == "~")
loadings_order_m <- subset(lavaan::lavaanify(tspa_mod_m, ngroup = 2),
                           op == "=~")

test_that("The order of factors in the model from tspa_mf()", {
  expect_equal(rep(c("visual", "textual", "textual"), 2),
               factors_order_m$lhs)
  expect_equal(rep(c("speed", "visual", "speed"), 2),
               factors_order_m$rhs)
})
test_that("The order of loadings in the model from tspa_mf()", {
  expect_equal(rep(rep(c("visual", "textual", "speed"), each = 3), 2),
               loadings_order_m$lhs)
  expect_equal(rep(c("fs_visual", "fs_textual", "fs_speed"), 6),
               loadings_order_m$rhs)
})

# Compare results to using Bartlett's scores
tspa_fit_m <- tspa(
  model = "visual ~ speed
           textual ~ visual + speed",
  data = fs_dat4,
  group = "school",
  fsT = attr(fs_dat4, "fsT"),
  fsL = attr(fs_dat4, "fsL")
)
fs_dat4b <- get_fs(HolzingerSwineford1939, model = mod4,
                    group = "school", method = "Bartlett", format = "list")
sem_fit_m <- sem(
  model = "visual =~ fs_visual
           speed =~ fs_speed
           textual =~ fs_textual
           fs_visual ~~ c(0.2633962, 0.2827317) * fs_visual
           fs_textual ~~ c(0.1239827, 0.1282725) * fs_textual
           fs_speed ~~ c(0.2020107, 0.1332701) * fs_speed
           visual ~ speed
           textual ~ visual + speed",
  data = do.call(rbind, fs_dat4b),
  group = "school"
)

test_that("Multiple-group multiple-factor example", code = {
  sct <- standardizedSolution(tspa_fit_m)
  scs <- standardizedSolution(sem_fit_m)
  expect_equal(sct$est[sct$op == "~"], expected = scs$est[scs$op == "~"],
               tolerance = 0.0001)
  expect_equal(sct$se[sct$op == "~"], expected = scs$se[scs$op == "~"],
               tolerance = 0.0001)
})

# An example from Chapter 14 of Grimm et al. (2016)
# https://quantdev.ssri.psu.edu/tutorials/growth-modeling-chapter-14-modeling-change-latent-variables-measured-continuous

mean_vec <- c(50.99, 65.25, 84.89, 127.66, 151.09, 172.05,
              99.72, 124.35, 142.47)
cov_mat <- matrix(c(
  232.71, 207.92, 188.09, 319.68, 285.26, 277.85, 260.75, 249.28, 217.96,
  207.92, 254.88, 212.14, 331.88, 313.8, 314.91, 274.99, 281.29, 243.6,
  188.09, 212.14, 270.46, 325.97, 308.84, 346.36, 284.9, 291.28, 281.55,
  319.68, 331.88, 325.97, 797.86, 617.02, 581.17, 511.8, 470.36, 420.6,
  285.26, 313.8, 308.84, 617.02, 662.41, 555.9, 448.81, 449.25, 394.63,
  277.85, 314.91, 346.36, 581.17, 555.9, 736.45, 440.78, 439.33, 443.67,
  260.75, 274.99, 284.9, 511.8, 448.81, 440.78, 618.23, 528.01, 437.92,
  249.28, 281.29, 291.28, 470.36, 449.25, 439.33, 528.01, 583.24, 448.64,
  217.96, 243.6, 281.55, 420.6, 394.63, 443.67, 437.92, 448.64, 480.57
), nrow = 9, ncol = 9, byrow = TRUE)
set.seed(123)
sim_dat <- MASS::mvrnorm(n = 2000, mu = mean_vec, Sigma = cov_mat,
                         empirical = TRUE)
colnames(sim_dat) <- c("s_g3", "s_g5", "s_g8", "r_g3", "r_g5", "r_g8",
                       "m_g3", "m_g5", "m_g8")

strict_mod <- "
# factor loadings
eta1 =~ 15.1749088 * s_g3 + l2 * r_g3 + l3 * m_g3
eta2 =~ 15.1749088 * s_g5 + l2 * r_g5 + l3 * m_g5
eta3 =~ 15.1749088 * s_g8 + l2 * r_g8 + l3 * m_g8

# factor variances/covariances
eta1 ~~ 1 * eta1 + eta2 + eta3
eta2 ~~ eta2 + eta3
eta3 ~~ eta3

# unique variances/covariances
s_g3 ~~ u1 * s_g3 + s_g5 + s_g8
s_g5 ~~ u1 * s_g5 + s_g8
s_g8 ~~ u1 * s_g8
r_g3 ~~ u2 * r_g3 + r_g5 + r_g8
r_g5 ~~ u2 * r_g5 + r_g8
r_g8 ~~ u2 * r_g8
m_g3 ~~ u3 * m_g3 + m_g5 + m_g8
m_g5 ~~ u3 * m_g5 + m_g8
m_g8 ~~ u3 * m_g8

# latent variable intercepts
eta1 ~ 0 * 1
eta2 ~ 1
eta3 ~ 1

# observed variable intercepts
s_g3 ~ i1 * 1
s_g5 ~ i1 * 1
s_g8 ~ i1 * 1
r_g3 ~ i2 * 1
r_g5 ~ i2 * 1
r_g8 ~ i2 * 1
m_g3 ~ i3 * 1
m_g5 ~ i3 * 1
m_g8 ~ i3 * 1
"
fs_growth_dat <- get_fs(sim_dat, model = strict_mod, format = "list")

growth_mod <- "
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + start(.5) * eta2 + 1 * eta3

# factor variances
eta1 ~~ psi * eta1
eta2 ~~ psi * eta2
eta3 ~~ psi * eta3

i ~~ start(.8) * i
s ~~ start(.5) * s
i ~~ start(0) * s

i ~ 1
s ~ 1
"
growth_fit <- tspa(growth_mod, fs_growth_dat,
                   fsT = attr(fs_growth_dat, "fsT"),
                   fsL = attr(fs_growth_dat, "fsL"),
                   fsb = attr(fs_growth_dat, "fsb"))

########## Error messages ##########

test_that("Empty path model", {
  expect_error(
    tspa(model = 123,
         data = fs_growth_dat,
         fsT = attr(fs_growth_dat, "fsT")),
    "The structural path model provided is not a string."
  )
})

test_that("Need to provide none or both fsT and fsL", {
  expect_error(
    tspa(model = growth_mod,
         data = fs_growth_dat,
         fsT = attr(fs_growth_dat, "fsT")),
    "Please provide both or none of fsT and fsL"
  )
  expect_error(
    tspa(model = growth_mod,
         data = fs_growth_dat,
         fsL = attr(fs_growth_dat, "fsL")),
    "Please provide both or none of fsT and fsL"
  )
  expect_no_error(
    tspa(growth_mod,
         data = fs_growth_dat,
         fsT = attr(fs_growth_dat, "fsT"),
         fsL = attr(fs_growth_dat, "fsL"),
         fsb = attr(fs_growth_dat, "fsb"))
  )
})

test_that("tspa(): non-numeric or non-finite se_fs is a clear error", {
  # non-numeric: a cryptic 'non-numeric argument' from the schema otherwise
  expect_error(
    tspa("dem60 ~ ind60", data = fs_dat_single,
         se_fs = c(ind60 = "0.12", dem60 = "0.68")),
    "numeric standard errors"
  )
  # non-finite: a NaN fixed value in the model string is a parse failure
  # in lavaan, not an actionable error
  expect_error(
    tspa("dem60 ~ ind60", data = fs_dat_single,
         se_fs = c(ind60 = 0.12, dem60 = NA_real_)),
    "not all finite"
  )
})

test_that("tspa(): non-finite or unnamed fsT/fsL values are clear errors", {
  data("PoliticalDemocracy", package = "lavaan")
  fit <- cfa("ind60 =~ x1 + x2 + x3
              dem60 =~ y1 + y2 + y3 + y4",
             data = PoliticalDemocracy, std.lv = TRUE)
  fs <- get_fs(fit)
  Tm <- attr(fs, "fsT")[[1L]]
  Lm <- attr(fs, "fsL")[[1L]]
  # an NA in fsT renders as "NA" in the model string, which lavaan parses
  # as a parameter label: the fixed value silently becomes a free estimate
  T_bad <- Tm
  T_bad[1L, 1L] <- NA
  expect_error(
    tspa("dem60 ~ ind60", data = fs, fsT = T_bad, fsL = Lm),
    "non-finite"
  )
  # unnamed fsL would render NA indicator names into the model string
  # (a cryptic rbind failure downstream)
  L_bare <- Lm
  dimnames(L_bare) <- list(NULL, colnames(Lm))
  expect_error(
    tspa("dem60 ~ ind60", data = fs, fsT = Tm, fsL = L_bare),
    "row and column names"
  )
})

test_that(
  "Names of factor score variables need to match those in the input data",
  {
    data("PoliticalDemocracy", package = "lavaan")
    mod2 <- "
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  "
    fs_dat2 <- get_fs(PoliticalDemocracy, model = mod2, std.lv = TRUE, format = "list")
    ecov_fs <- attr(fs_dat2, "fsT")
    dimnames(ecov_fs) <- lapply(dimnames(ecov_fs),
      FUN = function(x) paste0("bs_", x)
    )
    expect_error(
      tspa(
        model = "dem60 ~ ind60
              dem65 ~ ind60 + dem60",
        data = fs_dat2,
        fsT = ecov_fs,
        fsL = attr(fs_dat2, "fsL")
      ),
      "Names of factor score variables do not match those in the input data."
    )
    expect_no_error(
      tspa(
        model = "dem60 ~ ind60
              dem65 ~ ind60 + dem60",
        data = fs_dat2,
        fsT = attr(fs_dat2, "fsT"),
        fsL = attr(fs_dat2, "fsL")
      )
    )
  }
)

# Single-group unified get_fs() output carries length-1 list attributes;
# tspa() must accept them in any combination with plain matrices (e.g. an
# identity `fsL` for Bartlett scores).
mod2sg <- "
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
"
fs_dat_uni <- get_fs(PoliticalDemocracy, model = mod2sg, std.lv = TRUE)

test_that("Single-group list-valued attributes mixed with plain matrices", {
  fsT_u <- attr(fs_dat_uni, "fsT")
  fsL_u <- attr(fs_dat_uni, "fsL")
  expect_true(is.list(fsT_u) && length(fsT_u) == 1)
  expect_true(is.list(fsL_u) && length(fsL_u) == 1)
  identity_ld <- `dimnames<-`(diag(2),
                              list(c("fs_ind60", "fs_dem60"),
                                   c("ind60", "dem60")))
  # list fsT + matrix fsL (the Bartlett identity-loadings case)
  fit_tl <- tspa(model = "dem60 ~ ind60", data = fs_dat_uni,
                 fsT = fsT_u, fsL = identity_ld)
  fit_mm <- tspa(model = "dem60 ~ ind60", data = fs_dat_uni,
                 fsT = fsT_u[[1]], fsL = identity_ld)
  expect_equal(
    parameterestimates(fit_tl)[c("est", "se")],
    parameterestimates(fit_mm)[c("est", "se")],
    tolerance = 1e-10
  )
  # matrix fsT + list fsL (the reverse combination)
  fit_lt <- tspa(model = "dem60 ~ ind60", data = fs_dat_uni,
                 fsT = fsT_u[[1]], fsL = fsL_u)
  fit_mm2 <- tspa(model = "dem60 ~ ind60", data = fs_dat_uni,
                  fsT = fsT_u[[1]], fsL = fsL_u[[1]])
  expect_equal(
    parameterestimates(fit_lt)[c("est", "se")],
    parameterestimates(fit_mm2)[c("est", "se")],
    tolerance = 1e-10
  )
})

test_that("Multigroup fsT/fsL shape mismatch is a clear error", {
  expect_error(
    tspa(model = "visual ~ speed
                 textual ~ visual + speed",
         data = fs_dat4,
         group = "school",
         fsT = attr(fs_dat4, "fsT"),
         fsL = attr(fs_dat4, "fsL")[[1]]),
    "must be a list of the same length"
  )
  expect_error(
    tspa(model = "visual ~ speed
                 textual ~ visual + speed",
         data = fs_dat4,
         group = "school",
         fsT = attr(fs_dat4, "fsT")[[1]],
         fsL = attr(fs_dat4, "fsL")),
    "must be a list of the same length"
  )
})

# Nested (per-pattern) fsT/fsL are what get_fs() emits for a group with
# k >= 2 observed-indicator patterns (missing data). tspa() pools them
# (PLAN 09), but only when the data frame is a resolvable get_fs() result;
# with a plain (non-get_fs) data frame the pooler must error informatively
# instead of silently mis-pooling a single-group k-pattern list as k groups.
nested_fs_matrices <- function(vT, vL, pats) {
  tlist <- lapply(pats, function(p) {
    matrix(vT, 2, 2, dimnames = list(c("fs_a", "fs_b"), c("fs_a", "fs_b")))
  })
  llist <- lapply(pats, function(p) {
    matrix(vL, 2, 2, dimnames = list(c("fs_a", "fs_b"), c("a", "b")))
  })
  list(fsT = setNames(tlist, pats), fsL = setNames(llist, pats))
}

test_that("tspa(): per-unit fsT/fsL with a non-get_fs data frame errors informatively (SG shape)", {
  pats <- c("a+b", "a")
  nm <- nested_fs_matrices(0.1, 1, pats)
  dat <- data.frame(fs_a = 0:3, fs_b = 0:3)
  # attr() shape emitted by a single-group k=2 missing-data get_fs(): a
  # length-1 list whose element is a named list of per-pattern matrices
  # (a plain data frame carries no fs_pattern attribute, so the pooler
  # cannot resolve the rows and stops instead of mis-pooling)
  expect_error(
    tspa("b ~ a", data = dat,
         fsT = list(nm$fsT), fsL = list(nm$fsL)),
    "get_fs\\(\\) result"
  )
})

test_that("tspa(): per-unit fsT/fsL with a non-get_fs data frame errors informatively (MG shape); nested fsb alone keeps the backstop", {
  pats <- c("a+b", "a")
  nm <- nested_fs_matrices(0.1, 1, pats)
  dat <- data.frame(fs_a = c(0, 1, 2, 3), fs_b = c(0, 1, 2, 3),
                    school = c("V", "V", "G", "G"))
  expect_error(
    tspa("b ~ a", data = dat,
         fsT = setNames(rep(list(nm$fsT), 2), c("V", "G")),
         fsL = setNames(rep(list(nm$fsL), 2), c("V", "G")),
         group = "school"),
    "get_fs\\(\\) result"
  )
  # nested fsb (per-pattern intercept vectors, the single-group attr()
  # shape) is rejected even with plain matrix fsT/fsL
  fsT_m <- matrix(0.1, 2, 2, dimnames = list(c("fs_a", "fs_b"), c("fs_a", "fs_b")))
  fsL_m <- matrix(1, 2, 2, dimnames = list(c("fs_a", "fs_b"), c("a", "b")))
  fsb_ng <- setNames(list(c(a = 0, b = 0), c(a = 0, b = 0)), c("a+b", "a"))
  expect_error(
    tspa("b ~ a", data = dat,
         fsT = list(fsT_m), fsL = list(fsL_m),
         fsb = list(fsb_ng)),
    "does not yet support groups with multiple missing-data patterns"
  )
})

test_that("tspa() accepts complete-data attributes unchanged", {
  pats <- c("a+b")
  nm <- nested_fs_matrices(0.1, 1, pats)
  # fs_a/fs_b must not be perfectly correlated (singular sample covariance)
  dat <- data.frame(fs_a = c(0, 1, 2, 3), fs_b = c(3, 1, 2, 0))
  # k = 1 single-group output: plain matrix and length-1 list both accepted
  expect_no_error(
    fit1 <- suppressWarnings(tspa("b ~ a", data = dat,
                                  fsT = nm$fsT[[1]], fsL = nm$fsL[[1]]))
  )
  expect_no_error(
    fit2 <- suppressWarnings(tspa("b ~ a", data = dat,
                                  fsT = list(nm$fsT[[1]]),
                                  fsL = list(nm$fsL[[1]])))
  )
  expect_equal(
    parameterestimates(fit1)["est"],
    parameterestimates(fit2)["est"],
    tolerance = 1e-10
  )
})

# Per deep-value checks of the pooled per-group matrices (test-tspa_pooled.R),
# this only asserts that the real FIML get_fs() output now FITS and that the
# pooled (not the nested) values are what get attached to the fit.
test_that("tspa() fits real get_fs() missing-data output", {
  hs <- HolzingerSwineford1939
  set.seed(1334)
  hs[!rbinom(301, size = 1, prob = 0.7), 7] <- NA
  hs[!rbinom(301, size = 1, prob = 0.7), 8] <- NA
  hs[!rbinom(301, size = 1, prob = 0.7), 9] <- NA
  mod2 <- "
  # latent variables
    visual =~ x1 + x2 + x3
    speed  =~ x7 + x8 + x9
  "
  # multigroup: both groups carry nested per-pattern matrices, pooled to one
  # plain matrix per group
  fit_mg <- suppressWarnings(
    cfa(mod2, data = hs, group = "school", missing = "fiml")
  )
  fs_mg <- get_fs(fit_mg)
  expect_no_error(
    fit_mg2 <- suppressWarnings(
      tspa("visual ~ speed", data = fs_mg, group = "school",
           fsT = attr(fs_mg, "fsT"), fsL = attr(fs_mg, "fsL"))
    )
  )
  T_mg <- attr(fit_mg2, "fsT")
  expect_true(is.list(T_mg) && !is.null(names(T_mg)))
  expect_true(all(vapply(T_mg, is.matrix, logical(1))))   # no nesting left
  expect_equal(attr(fit_mg2, "pooled_fs"), "mean")
  # single group: the length-1 attribute list wrapping a nested pattern
  # list (the "misread as k groups" trap) pools to a plain matrix
  fit_sg <- suppressWarnings(cfa(mod2, data = hs, missing = "fiml"))
  fs_sg <- get_fs(fit_sg)
  expect_no_error(
    fit_sg2 <- suppressWarnings(
      tspa("visual ~ speed", data = fs_sg,
           fsT = attr(fs_sg, "fsT"), fsL = attr(fs_sg, "fsL"))
    )
  )
  expect_true(is.matrix(attr(fit_sg2, "fsT")))
  expect_equal(attr(fit_sg2, "pooled_fs"), "mean")
})

test_that("Test indicator names not starting with 'fs_'", {
  data("PoliticalDemocracy", package = "lavaan")
  mod2 <- "
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  "
  fs_dat2 <- get_fs(PoliticalDemocracy, model = mod2, std.lv = TRUE, format = "list")
  names(fs_dat2) <- gsub("fs_", "bs_", names(fs_dat2))
  ecov_fs <- attr(fs_dat2, "fsT")
  dimnames(ecov_fs) <- lapply(dimnames(ecov_fs),
                              FUN = function(x) gsub("fs_", "bs_", x))
  mat_ld <- attr(fs_dat2, "fsL")
  rownames(mat_ld) <- gsub("fs_", "bs_", rownames(mat_ld))
  expect_no_error(
    bs_fit <- tspa(model = "dem60 ~ ind60
                            dem65 ~ ind60 + dem60",
                   data = fs_dat2,
                   fsT = ecov_fs,
                   fsL = mat_ld)
  )
  fs_fit <- tspa(model = "dem60 ~ ind60
                          dem65 ~ ind60 + dem60",
                  data = get_fs(PoliticalDemocracy, model = mod2, std.lv = TRUE, format = "list"),
                 fsT = attr(fs_dat2, "fsT"),
                 fsL = attr(fs_dat2, "fsL"))
  expect_identical(
    parameterestimates(bs_fit)["est"], parameterestimates(fs_fit)["est"]
  )
})

test_that("Missing group argument for a multigroup model", {
  expect_error(
    tspa(
      model = "visual ~ speed
               textual ~ visual + speed",
      data = fs_dat4,
      fsT = attr(fs_dat4, "fsT"),
      fsL = attr(fs_dat4, "fsL")
    ),
    "Please specify 'group = ' to fit a multigroup model in lavaan"
  )
})

########## tspa_args self-contained replay ##########

test_that("tspa() records a self-contained replayable argument list", {
  # Single-group multi-factor fit (the vcov_corrected() SG use case).
  mod2s <- "ind60 =~ x1 + x2 + x3
            dem60 =~ y1 + y2 + y3 + y4"
  fs_sg_mf <- get_fs(cfa(mod2s, data = PoliticalDemocracy))
  fit_sg_mf <- tspa(
    model = "dem60 ~ ind60",
    data = fs_sg_mf,
    fsT = attr(fs_sg_mf, "fsT"),
    fsL = attr(fs_sg_mf, "fsL")
  )
  for (fit in list(tspa_single, tspa_multi, tspa_fit_m, fit_sg_mf)) {
    args <- attr(fit, "tspa_args")
    expect_true(is.list(args))
    # Evaluated values only: no calls/symbols (a match.call() style record
    # would trip this).
    expect_false(any(vapply(args, is.language, logical(1))))
    refit <- do.call(tspa, args)
    expect_equal(coef(refit), coef(fit), tolerance = 1e-10)
    expect_equal(vcov(refit), vcov(fit), ignore_attr = TRUE, tolerance = 1e-10)
  }
})
