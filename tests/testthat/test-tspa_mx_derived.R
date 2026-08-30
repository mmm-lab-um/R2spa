# =====================================================================
# PLAN 15: tspa_mx_model() auto-derives the measurement inputs (fsL/fsT/fsb)
# from a get_fs() result when the caller omits se_fs/fsL/fsT/fsb entirely.
# Explicit arguments always win (D1); constant quantities become fixed
# numeric cells (D2); per-row (per_obs/mirt_per_obs) and per-pattern
# (SG FIML) quantities become definition-variable matrices (D3); a
# provenance gate rejects hand-rolled attributes (D4); the D5 fail-fast
# replaces the old misleading "'fsL' rows must be named..." error; the
# existing data-contract guards (incl. the NA-free defvar guard, D6) run on
# the possibly int_fs_-augmented frame; multigroup (group_col) input is
# refused (D7).
#
# A/B convention (mirrors test-tspa_derived.R): every derived fit is
# compared against the explicit-argument control fit on the same data, via
# coef()/vcov(). Derived and control calls build the identical spec and
# feed the identical column values to OpenMx, so the pairs below agree
# bit-exactly (expect_identical on unname(coef()) + tight-tolerance
# coefficient extractions, as in test-tspa_mx.R).
#
# q >= 2 off-diagonal defvar models were once pinned at the string level
# only (implementation finding V3c, plan section 5 item 3): they aborted
# with "implied covariance not positive definite" on the derived route.
# The abort was NOT an OpenMx limitation -- it was a '~~' defvar-lookup
# orientation bug: lavaanify() may present a covariance row with (lhs, rhs)
# reversed relative to the score order, so a lower-triangle-only fsT (the
# documented and the derived convention) was not found and the c(1) defvar
# sentinel leaked into the model as a fixed unit covariance between scores.
# Fixed in tspa_mx_defvar_col(); the q >= 2 cases (SG FIML per-pattern and
# mirt 2-factor per-row) are pinned numerically end-to-end below.
# =====================================================================

library(lavaan)
library(lme4)
skip_if_not_installed("OpenMx")

# --- coefficient extraction by (from -> to), robust to RAM var ordering ----
# (same pattern as test-tspa_mx.R)
mx_path_val <- function(m, from, to, model = "m1") {
  v <- c(m$manifestVars, m$latentVars)
  unname(coef(m)[sprintf("%s.A[%d,%d]", model, match(to, v), match(from, v))])
}
mx_var_val <- function(m, x, model = "m1") {
  v <- c(m$manifestVars, m$latentVars)
  unname(coef(m)[sprintf("%s.S[%d,%d]", model, match(x, v), match(x, v))])
}

# A copy of a get_fs() frame with the measurement attributes stripped: the
# "explicit arguments only" control frame (derivation cannot fire on it, so a
# fit on it isolates the explicit-args path). Same helper as
# test-tspa_derived.R.
strip_fs_attrs <- function(d) {
  at <- attributes(d)
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix", "fs_pattern",
               "group_col", "mirt_per_obs", "psi", "alpha")) {
    at[[ak]] <- NULL
  }
  attributes(d) <- at
  d
}

# The full measurement model string that derivation would feed to OpenMx for
# a get_fs() result (the user model plus the derived measurement block).
derived_model_string <- function(data, model) {
  der <- R2spa:::tspa_mx_derive_measurement(data)
  R2spa:::tspa_mx_model_string(
    model, R2spa:::tspa_mx_spec(NULL, der$fsL, der$fsT, der$fsb)
  )
}

########## Fixtures and A/B fits (file scope) ##########

# --- 1. SG complete, 2-factor (unified; constant attributes) --------------
mod2 <- "ind60 =~ x1 + x2 + x3\ndem60 =~ y1 + y2 + y3 + y4"
model2 <- "dem60 ~ ind60"
fs1 <- get_fs(PoliticalDemocracy, mod2)
fit_d1 <- suppressWarnings(tspa_mx_model(model2, data = fs1))
fit_e1 <- suppressWarnings(tspa_mx_model(model2, data = fs1,
                                         fsL = attr(fs1, "fsL"),
                                         fsT = attr(fs1, "fsT"),
                                         fsb = attr(fs1, "fsb")))

# --- 2. local = TRUE compact, 3-factor (plain matrices -> D2) --------------
mod3 <- paste("ind60 =~ x1 + x2 + x3", "dem60 =~ y1 + y2 + y3 + y4",
              "dem65 =~ y5 + y6 + y7 + y8", sep = "\n")
model3 <- "dem60 ~ ind60\ndem65 ~ ind60 + dem60"
fs_local3 <- get_fs(PoliticalDemocracy, mod3, local = TRUE)
fit_d2 <- suppressWarnings(tspa_mx_model(model3, data = fs_local3))
fit_e2 <- suppressWarnings(tspa_mx_model(model3, data = fs_local3,
                                         fsL = attr(fs_local3, "fsL"),
                                         fsT = attr(fs_local3, "fsT"),
                                         fsb = attr(fs_local3, "fsb")))

# --- 3. SG FIML joint, 2-factor, with (fixed, nonzero) factor means --------
# String/int_fs_ pins only: a full q = 2 per-pattern defvar fit does not
# optimize in OpenMx 2.22.11 (V3c, see header).
mod_hs <- "visual =~ x1 + x2 + x3\nspeed =~ x7 + x8 + x9"
hs_fiml2 <- HolzingerSwineford1939
set.seed(1334)
hs_fiml2$x2[!rbinom(nrow(hs_fiml2), 1L, 0.4)] <- NA
hs_fiml2$x8[!rbinom(nrow(hs_fiml2), 1L, 0.4)] <- NA
fit_fiml2 <- suppressWarnings(
  sem(paste(mod_hs, "visual ~ 1.2*1", "speed ~ 2*1", sep = "\n"),
      data = hs_fiml2, missing = "fiml")
)
fs_fiml2 <- suppressWarnings(get_fs(fit_fiml2))
# The manual definition-variable matrices the derivation must reproduce
# (the documented roxygen-example route, pattern-constant columns).
Lm_fiml <- matrix(c("visual_by_fs_visual", "visual_by_fs_speed",
                    "speed_by_fs_visual", "speed_by_fs_speed"),
                  nrow = 2L,
                  dimnames = list(c("fs_visual", "fs_speed"),
                                  c("visual", "speed")))
Tm_fiml <- matrix(c("ev_fs_visual", "ecov_fs_speed_fs_visual", NA,
                    "ev_fs_speed"),
                  nrow = 2L,
                  dimnames = list(c("fs_visual", "fs_speed"),
                                  c("fs_visual", "fs_speed")))
bm_fiml <- c(fs_visual = "int_fs_visual", fs_speed = "int_fs_speed")
# 3b: the same q = 2 case fitted end-to-end. The lower-triangle-only
# Tm_fiml is the sharper pin: it exercises the '~~' lookup fallback for the
# orientation in which lavaanify() presents the (lhs, rhs) pair.
fit_d3b <- suppressWarnings(tspa_mx_model("visual ~ speed", data = fs_fiml2))
fit_e3b <- suppressWarnings(tspa_mx_model(
  "visual ~ speed", data = fs_indiv(fs_fiml2, include_intercept = TRUE),
  fsL = Lm_fiml, fsT = Tm_fiml, fsb = bm_fiml
))

# --- 3c/4a. q = 1 fixtures (single latent -> diagonal-only defvars fit) ----
mod_hs1 <- "visual =~ x1 + x2 + x3"
hs_fiml1 <- HolzingerSwineford1939
set.seed(42)
hs_fiml1$x1[!rbinom(nrow(hs_fiml1), 1L, 0.4)] <- NA
hs_fiml1$x2[!rbinom(nrow(hs_fiml1), 1L, 0.4)] <- NA
fit_fiml1 <- suppressWarnings(
  sem(paste(mod_hs1, "visual ~ 2*1", sep = "\n"),
      data = hs_fiml1, missing = "fiml")
)
fs_fiml1 <- suppressWarnings(get_fs(fit_fiml1))
# local FIML: per-row (per_obs) result, same data
fs_locfiml <- suppressWarnings(
  get_fs(hs_fiml1, mod_hs1, local = TRUE, missing = "fiml")
)
Lq <- matrix("visual_by_fs_visual", dimnames = list("fs_visual", "visual"))
Tq <- matrix("ev_fs_visual", dimnames = list("fs_visual", "fs_visual"))
bq <- c(fs_visual = "int_fs_visual")
# 3c: SG FIML q=1 derived vs manual def-var A/B
fit_d3c <- suppressWarnings(tspa_mx_model("visual ~ 1", data = fs_fiml1))
fit_e3c <- suppressWarnings(tspa_mx_model(
  "visual ~ 1", data = fs_indiv(fs_fiml1, include_intercept = TRUE),
  fsL = Lq, fsT = Tq, fsb = bq
))
# 4a: local FIML q=1 derived vs manual def-var A/B (the per-row result
# carries an all-zero per-row fsb, so the derived spec also models the score
# mean through int_fs_visual; the manual control mirrors that)
fit_d4a <- suppressWarnings(tspa_mx_model("visual ~ 1", data = fs_locfiml))
fit_e4a <- suppressWarnings(tspa_mx_model(
  "visual ~ 1", data = fs_indiv(fs_locfiml, include_intercept = TRUE),
  fsL = Lq, fsT = Tq, fsb = bq
))

# --- 5a. Mean-structure CFA: nonzero constant fsb -> D2 numeric path -------
# (lavaan fixed factor means, `f ~ N*1`, so fsb = alpha - fsL %*% alpha is
# nonzero without relying on mean-structure estimation)
fit_ms <- suppressWarnings(
  sem(paste(mod2, "ind60 ~ 3*1", "dem60 ~ 4*1", sep = "\n"),
      data = PoliticalDemocracy)
)
fs_ms <- get_fs(fit_ms)
fit_d5 <- suppressWarnings(tspa_mx_model(model2, data = fs_ms))
fit_e5 <- suppressWarnings(tspa_mx_model(model2, data = fs_ms,
                                         fsL = attr(fs_ms, "fsL"),
                                         fsT = attr(fs_ms, "fsT"),
                                         fsb = attr(fs_ms, "fsb")))

# --- 5c. merMod (per-cluster 3-D fsL/fsT arrays; no fsb attached) -----------
fs_mer <- get_fs(lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
Lm_mer <- matrix(c("u0_by_fs_u0", "u0_by_fs_u1",
                   "u1_by_fs_u0", "u1_by_fs_u1"),
                 nrow = 2L,
                 dimnames = list(c("fs_u0", "fs_u1"), c("u0", "u1")))
Tm_mer <- matrix(c("ev_fs_u0", "ecov_fs_u1_fs_u0", NA, "ev_fs_u1"),
                 nrow = 2L,
                 dimnames = list(c("fs_u0", "fs_u1"), c("fs_u0", "fs_u1")))
fit_d5c <- suppressWarnings(tspa_mx_model("u0 ~ u1\nu0 + u1 ~ 1",
                                          data = fs_mer))
fit_e5c <- suppressWarnings(tspa_mx_model("u0 ~ u1\nu0 + u1 ~ 1",
                                          data = fs_mer,
                                          fsL = Lm_mer, fsT = Tm_mer))

# --- 6. Explicit-wins (D1) -------------------------------------------------
se2 <- c(ind60 = 0.1213615, dem60 = 0.6756472)
fit_6a_d <- suppressWarnings(tspa_mx_model(model2, data = fs1, se_fs = se2))
fit_6a_e <- suppressWarnings(tspa_mx_model(model2, data = strip_fs_attrs(fs1),
                                           se_fs = se2))
# fully explicit on the attribute-less frame: the explicit path must be
# unaffected by the attributes' presence (fit_e1 is the same explicit call
# on the attribute-bearing frame from item 1)
fit_6c_s <- suppressWarnings(tspa_mx_model(model2, data = strip_fs_attrs(fs1),
                                           fsL = attr(fs1, "fsL"),
                                           fsT = attr(fs1, "fsT"),
                                           fsb = attr(fs1, "fsb")))

# --- 7. Fail-fast (D4/D5/D7) fixtures --------------------------------------
set.seed(42)
d_hand <- data.frame(fs_ind60 = rnorm(100L), fs_dem60 = rnorm(100L))
attr(d_hand, "fsT") <- matrix(c(0.2, 0.1, 0.1, 0.3), 2L, 2L,
                              dimnames = list(c("fs_ind60", "fs_dem60"),
                                              c("fs_ind60", "fs_dem60")))
attr(d_hand, "fsL") <- `dimnames<-`(diag(2L),
                                    list(c("fs_ind60", "fs_dem60"),
                                         c("ind60", "dem60")))
d_plain <- data.frame(fs_ind60 = rnorm(50L), fs_dem60 = rnorm(50L))
fsmg <- get_fs(HolzingerSwineford1939, mod_hs, std.lv = TRUE,
               group = "school")

# --- 9. Product-score result (product columns inert to derivation) ---------
fs_prod <- get_fs(PoliticalDemocracy, mod2, product = "ind60:dem60")
# OpenMx mxData() rejects the ':' in the product column names (a pre-existing
# OpenMx limitation, hit by the explicit route as well); drop the product
# columns -- keeping every measurement attribute -- for the numerical A/B.
fs_prod_stripped <- fs_prod[, !grepl(":", names(fs_prod)), drop = FALSE]
for (ak in setdiff(names(attributes(fs_prod)),
                   c("names", "row.names", "class"))) {
  attr(fs_prod_stripped, ak) <- attr(fs_prod, ak)
}
fit_d9 <- suppressWarnings(tspa_mx_model(model2, data = fs_prod_stripped))
fit_e9 <- suppressWarnings(tspa_mx_model(model2, data = fs_prod_stripped,
                                         fsL = attr(fs_prod, "fsL"),
                                         fsT = attr(fs_prod, "fsT"),
                                         fsb = attr(fs_prod, "fsb")))

########## Tests ##########

test_that("derived SG 2-factor equals the explicit full-triple fit (D2 fixed numeric)", {
  expect_identical(unname(coef(fit_d1)), unname(coef(fit_e1)))
  expect_equal(mx_path_val(fit_d1, "ind60", "dem60"),
               mx_path_val(fit_e1, "ind60", "dem60"), tolerance = 1e-10)
  expect_equal(mx_var_val(fit_d1, "ind60"),
               mx_var_val(fit_e1, "ind60"), tolerance = 1e-10)
  expect_equal(mx_var_val(fit_d1, "dem60"),
               mx_var_val(fit_e1, "dem60"), tolerance = 1e-10)
  expect_equal(vcov(fit_d1), vcov(fit_e1), tolerance = 1e-10,
               ignore_attr = TRUE)
  # D2 dispatch: the unified plain-matrix attributes derive to fixed numeric
  # matrices (not definition variables), and the frame is returned untouched
  der1 <- R2spa:::tspa_mx_derive_measurement(fs1)
  expect_null(der1$prov_err)
  expect_true(is.matrix(der1$fsL) && is.matrix(der1$fsT))
  expect_identical(der1$data, fs1)
  s1 <- R2spa:::tspa_mx_model_string(
    model2, R2spa:::tspa_mx_spec(NULL, der1$fsL, der1$fsT, der1$fsb)
  )
  expect_false(any(grepl("data\\.", strsplit(s1, "\n")[[1L]])))
})

test_that("derived local = TRUE 3-factor equals the explicit full-triple fit", {
  expect_identical(unname(coef(fit_d2)), unname(coef(fit_e2)))
  expect_equal(mx_path_val(fit_d2, "ind60", "dem60"),
               mx_path_val(fit_e2, "ind60", "dem60"), tolerance = 1e-10)
  expect_equal(mx_path_val(fit_d2, "ind60", "dem65"),
               mx_path_val(fit_e2, "ind60", "dem65"), tolerance = 1e-10)
  expect_equal(mx_path_val(fit_d2, "dem60", "dem65"),
               mx_path_val(fit_e2, "dem60", "dem65"), tolerance = 1e-10)
  expect_equal(vcov(fit_d2), vcov(fit_e2), tolerance = 1e-10,
               ignore_attr = TRUE)
  # D2: the local (compact) attributes are plain matrices with exact-zero
  # cross terms -- derived as fixed numeric cells, no definition variables
  der2 <- R2spa:::tspa_mx_derive_measurement(fs_local3)
  expect_true(is.matrix(der2$fsL) && is.matrix(der2$fsT))
  off <- der2$fsL
  diag(off) <- NA
  expect_true(all(is.na(off) | off == 0))
  expect_false(any(grepl("data\\.",
                         strsplit(derived_model_string(fs_local3, model3),
                                  "\n")[[1L]])))
})

test_that("derived SG FIML (per-pattern) reproduces the manual def-var spec and int_fs_* columns (V3c: string-level pin)", {
  der3 <- R2spa:::tspa_mx_derive_measurement(fs_fiml2)
  expect_null(der3$prov_err)
  # the character matrices are exactly the manual roxygen-example route
  expect_identical(der3$fsL, Lm_fiml)
  expect_identical(der3$fsT, Tm_fiml)
  expect_identical(der3$fsb, bm_fiml)
  # ...and so is the full model string
  expect_identical(
    derived_model_string(fs_fiml2, "visual ~ speed"),
    R2spa:::tspa_mx_model_string(
      "visual ~ speed",
      R2spa:::tspa_mx_spec(NULL, Lm_fiml, Tm_fiml, bm_fiml)
    )
  )
  # every referenced definition-variable column exists in the derived frame
  dv <- c(unlist(der3$fsL), unlist(der3$fsT), unlist(der3$fsb))
  dv <- dv[!is.na(dv)]
  expect_true(all(dv %in% names(der3$data)))
  # the appended int_fs_* columns are the fs_indiv(include_intercept = TRUE)
  # values, row-wise (mapped through fs_pattern$label), and nonzero here
  dat3 <- fs_indiv(fs_fiml2, include_intercept = TRUE)
  expect_identical(der3$data[["int_fs_visual"]], dat3[["int_fs_visual"]])
  expect_identical(der3$data[["int_fs_speed"]], dat3[["int_fs_speed"]])
  expect_true(any(der3$data[["int_fs_visual"]] != 0))
})

test_that("derived SG FIML q = 2 (per-pattern) end-to-end equals the manual def-var fit", {
  # Regression for the '~~' defvar-lookup orientation bug: with the
  # lower-triangle-only Tm_fiml, the derived fit must equal the manual
  # character-matrix fit (both previously aborted with a leaked fixed unit
  # covariance whenever lavaanify() presented the pair reversed).
  expect_identical(unname(coef(fit_d3b)), unname(coef(fit_e3b)))
  expect_equal(mx_path_val(fit_d3b, "speed", "visual"),
               mx_path_val(fit_e3b, "speed", "visual"), tolerance = 1e-10)
  expect_equal(vcov(fit_d3b), vcov(fit_e3b), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("derived SG FIML q = 1 end-to-end equals the manual def-var fit", {
  # the q = 1 (diagonal-only) case (no off-diagonal defvars to look up)
  expect_identical(unname(coef(fit_d3c)), unname(coef(fit_e3c)))
  expect_equal(mx_var_val(fit_d3c, "visual"),
               mx_var_val(fit_e3c, "visual"), tolerance = 1e-10)
})

test_that("the '~~' defvar lookup is orientation-agnostic for lower-triangle-only fsT", {
  S2 <- c("fs_F1", "fs_F2")
  Lm <- `dimnames<-`(matrix(c("F1_by_fs_F1", NA, NA, "F2_by_fs_F2"), 2L),
                     list(S2, c("F1", "F2")))
  Tm <- matrix(NA_character_, 2L, 2L, dimnames = list(S2, S2))
  Tm[1L, 1L] <- "ev_fs_F1"
  Tm[2L, 1L] <- "ecov_fs_F2_fs_F1" # lower triangle only
  Tm[2L, 2L] <- "ev_fs_F2"
  spec <- R2spa:::tspa_mx_spec(NULL, Lm, Tm, NULL)
  expect_identical(R2spa:::tspa_mx_defvar_col(spec, "fs_F2", "~~", "fs_F1"),
                   "ecov_fs_F2_fs_F1")
  # the orientation in which lavaanify() actually presents the row
  expect_identical(R2spa:::tspa_mx_defvar_col(spec, "fs_F1", "~~", "fs_F2"),
                   "ecov_fs_F2_fs_F1")
  expect_identical(R2spa:::tspa_mx_defvar_col(spec, "fs_F1", "~~", "fs_F1"),
                   "ev_fs_F1")
})

test_that("derived local FIML (per_obs, q = 1) equals the manual def-var fit", {
  expect_true(isTRUE(attr(fs_locfiml, "per_obs")))
  expect_identical(unname(coef(fit_d4a)), unname(coef(fit_e4a)))
  expect_equal(mx_var_val(fit_d4a, "visual"),
               mx_var_val(fit_e4a, "visual"), tolerance = 1e-10)
  # the per-row (here all-zero) fsb list still appends int_fs_visual, equal
  # to the fs_indiv(include_intercept = TRUE) column
  der4a <- R2spa:::tspa_mx_derive_measurement(fs_locfiml)
  expect_identical(der4a$data[["int_fs_visual"]],
                   fs_indiv(fs_locfiml, include_intercept = TRUE)[[
                     "int_fs_visual"
                   ]])
})

test_that("derived mean-structure CFA (nonzero constant fsb) equals the explicit fsb fit (D2)", {
  b_ms <- attr(fs_ms, "fsb")[[1L]]
  expect_true(any(abs(b_ms) > 0))
  expect_identical(unname(coef(fit_d5)), unname(coef(fit_e5)))
  expect_equal(mx_path_val(fit_d5, "ind60", "dem60"),
               mx_path_val(fit_e5, "ind60", "dem60"), tolerance = 1e-10)
  expect_equal(vcov(fit_d5), vcov(fit_e5), tolerance = 1e-10,
               ignore_attr = TRUE)
  # the constant fsb attribute derives to a plain (fixed numeric) vector
  der5 <- R2spa:::tspa_mx_derive_measurement(fs_ms)
  expect_true(is.numeric(der5$fsb))
  expect_equal(unname(der5$fsb), unname(b_ms), tolerance = 0)
})

test_that("derived merMod (per-cluster 3-D arrays, no fsb) equals the manual def-var fit", {
  der5c <- R2spa:::tspa_mx_derive_measurement(fs_mer)
  expect_null(der5c$prov_err)
  # character defvar matrices over the frame's own per-cluster columns;
  # merMod attaches no fsb -> fixed-zero score intercepts, frame untouched
  expect_identical(der5c$fsL, Lm_mer)
  expect_identical(der5c$fsT, Tm_mer)
  expect_null(der5c$fsb)
  expect_identical(der5c$data, fs_mer)
  dv <- c(unlist(der5c$fsL), unlist(der5c$fsT))
  dv <- dv[!is.na(dv)]
  expect_true(all(dv %in% names(fs_mer)))
  expect_true(all(!is.na(fs_mer[, dv, drop = FALSE])))
  expect_identical(unname(coef(fit_d5c)), unname(coef(fit_e5c)))
  expect_equal(mx_path_val(fit_d5c, "u0", "u1"),
               mx_path_val(fit_e5c, "u0", "u1"), tolerance = 1e-10)
  expect_equal(mx_var_val(fit_d5c, "u0"), mx_var_val(fit_e5c, "u0"),
               tolerance = 1e-10)
  expect_equal(mx_var_val(fit_d5c, "u1"), mx_var_val(fit_e5c, "u1"),
               tolerance = 1e-10)
  expect_equal(vcov(fit_d5c), vcov(fit_e5c), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("explicit arguments always win over derivation (D1)", {
  # (a) se_fs given on a get_fs() frame: the fit equals the same call on the
  # attribute-less frame, proving no derivation from the attributes happened
  expect_identical(unname(coef(fit_6a_d)), unname(coef(fit_6a_e)))
  expect_equal(mx_path_val(fit_6a_d, "ind60", "dem60"),
               mx_path_val(fit_6a_e, "ind60", "dem60"), tolerance = 1e-10)
  # (b) exactly one of fsL/fsT: the pre-existing xor error, before derivation
  expect_error(
    tspa_mx_model(model2, data = fs1, fsL = attr(fs1, "fsL")),
    "Provide both 'fsL' and 'fsT', or use 'se_fs'.", fixed = TRUE
  )
  # (c) the fully explicit call is unchanged by the feature: identical on the
  # attribute-bearing and attribute-less frames
  expect_identical(unname(coef(fit_e1)), unname(coef(fit_6c_s)))
})

test_that("fail-fast errors (D4/D5/D7) name the accepted input forms", {
  # (a) hand-rolled frame: plain matrix fsT/fsL attributes but no get_fs()
  # provenance -- the D5 fail-fast carries the gate's message
  err <- expect_error(tspa_mx_model(model2, data = d_hand))
  msg <- conditionMessage(err)
  expect_match(msg, "No measurement inputs found", fixed = TRUE)
  expect_match(msg, "does not look like a get_fs() result", fixed = TRUE)
  expect_match(msg, "has no 'fs_pattern' attribute", fixed = TRUE)
  # (b) plain non-get_fs frame (no attributes): the base D5 message, with no
  # gate addendum
  err2 <- expect_error(tspa_mx_model(model2, data = d_plain))
  msg2 <- conditionMessage(err2)
  expect_match(msg2, "No measurement inputs found", fixed = TRUE)
  expect_false(grepl("does not look like a get_fs() result", msg2,
                     fixed = TRUE))
  # (c) unified multigroup result (group_col attribute): the D7 refusal
  expect_error(
    tspa_mx_model("visual ~ speed", data = fsmg),
    "Multigroup 'fsL'/'fsT' are not supported yet (Phase 1 is single-group).",
    fixed = TRUE
  )
})

test_that("product-score result: product columns are inert to derivation", {
  # the raw product frame cannot reach mxData at all: OpenMx rejects the ':'
  # in the product column names -- and the explicit route fails identically
  # (a pre-existing OpenMx naming limitation, not a derivation defect)
  err_d <- expect_error(
    tspa_mx_model(model2, data = fs_prod),
    "is illegal because it contains the ':' character", fixed = TRUE
  )
  err_e <- expect_error(
    tspa_mx_model(model2, data = fs_prod, fsL = attr(fs_prod, "fsL"),
                  fsT = attr(fs_prod, "fsT"), fsb = attr(fs_prod, "fsb")),
    "is illegal because it contains the ':' character", fixed = TRUE
  )
  expect_identical(conditionMessage(err_d), conditionMessage(err_e))
  # with the product columns dropped (attributes kept), the derived fit
  # equals the explicit-attribute control ...
  expect_identical(unname(coef(fit_d9)), unname(coef(fit_e9)))
  expect_equal(mx_path_val(fit_d9, "ind60", "dem60"),
               mx_path_val(fit_e9, "ind60", "dem60"), tolerance = 1e-10)
  # ... and the derived spec is identical with or without the product columns
  expect_identical(derived_model_string(fs_prod, model2),
                   derived_model_string(fs_prod_stripped, model2))
})

# ---- mirt (per-observation fsT/fsL; mirt is Suggests-only) -----------------

skip_if_not_installed("mirt")

set.seed(2026)
mirt_sim1f <- function(N) {
  as.data.frame(mirt::simdata(
    a = matrix(runif(5L, 0.5, 1.5), 5L, 1L),
    d = rnorm(5L), N = N, itemtype = "2PL", Theta = rnorm(N)
  ))
}

# --- 4b: mirt SG, single factor (mirt_per_obs) ------------------------------
mirt_sg <- suppressWarnings(mirt::mirt(mirt_sim1f(200L), 1L, verbose = FALSE))
fs_mirt <- get_fs(mirt_sg)
# the manual character-matrix route (q = 1 analogue of the roxygen example)
Lm1 <- matrix("F1_by_fs_F1", dimnames = list("fs_F1", "F1"))
Tm1 <- matrix("ev_fs_F1", dimnames = list("fs_F1", "fs_F1"))
bm1 <- c(fs_F1 = "int_fs_F1")
fit_m_d <- suppressWarnings(tspa_mx_model("F1 ~ 1", data = fs_mirt))
fit_m_e <- suppressWarnings(tspa_mx_model(
  "F1 ~ 1", data = fs_indiv(fs_mirt, include_intercept = TRUE),
  fsL = Lm1, fsT = Tm1, fsb = bm1
))

# --- 5b: mirt SG with prior_mean (nonzero per-row fsb) ----------------------
fs_mirt_pm <- get_fs(mirt_sg, prior_mean = c(F1 = 2))
fit_pm_d <- suppressWarnings(tspa_mx_model("F1 ~ 1", data = fs_mirt_pm))
fit_pm_e <- suppressWarnings(tspa_mx_model(
  "F1 ~ 1", data = fs_indiv(fs_mirt_pm, include_intercept = TRUE),
  fsL = Lm1, fsT = Tm1, fsb = bm1
))

# --- 4c: mirt MG (MultipleGroupClass: group column, no group_col) -----------
n_g <- 100L
set.seed(2026)
dat_mg <- rbind(mirt_sim1f(n_g), mirt_sim1f(n_g))
grp_mg <- factor(rep(c("A", "B"), each = n_g))
mirt_mg <- suppressWarnings(
  mirt::multipleGroup(dat_mg, 1L, group = grp_mg, invariance = "slopes",
                      verbose = FALSE)
)
fs_mirt_mg <- get_fs(mirt_mg)
fit_mg_d <- suppressWarnings(tspa_mx_model("F1 ~ 1", data = fs_mirt_mg))
fit_mg_e <- suppressWarnings(tspa_mx_model(
  "F1 ~ 1", data = fs_indiv(fs_mirt_mg, include_intercept = TRUE),
  fsL = Lm1, fsT = Tm1, fsb = bm1
))

# --- 4d: mirt SG, two factors (off-diagonal defvars) ------------------------
# The vignette case (Lai & Hsiao 2022): the q >= 2 per-row fit that aborted
# before the '~~' lookup fix. The manual control uses the documented
# lower-triangle-only fsT (the vignette's full-triangle matrix works too).
set.seed(1235)
n2d <- 300L
eta2d <- MASS::mvrnorm(n2d, mu = c(0, 0), Sigma = diag(c(1, 1 - 0.5^2)),
                       empirical = TRUE)
th2d_1 <- eta2d[, 1]
th2d_2 <- -1 + 0.5 * th2d_1 + eta2d[, 2]
dat2d_a <- as.data.frame(mirt::simdata(a = matrix(1, 5L), d = rnorm(5L),
                                       N = n2d, itemtype = "2PL",
                                       Theta = th2d_1))
dat2d_b <- as.data.frame(mirt::simdata(a = matrix(runif(5L, 0.5, 1.5), 5L),
                                       d = rnorm(5L), N = n2d,
                                       itemtype = "2PL", Theta = th2d_2))
dat2d <- cbind(dat2d_a, dat2d_b)
colnames(dat2d) <- paste0("Item_", 1:10)
mirt_2f <- suppressWarnings(
  mirt::mirt(dat2d, model = "F1 = 1-5\nF2 = 6-10\nCOV = F1*F2",
             itemtype = "2PL", verbose = FALSE)
)
fs_mirt2f <- get_fs(mirt_2f)
S2f <- c("fs_F1", "fs_F2")
Lm2f <- matrix(c("F1_by_fs_F1", "F1_by_fs_F2", "F2_by_fs_F1",
                 "F2_by_fs_F2"),
               nrow = 2L, dimnames = list(S2f, c("F1", "F2")))
# lower triangle only (column-major fill: [1,1], [2,1], [1,2], [2,2])
Tm2f <- `dimnames<-`(
  matrix(c("ev_fs_F1", "ecov_fs_F2_fs_F1", NA, "ev_fs_F2"), 2L),
  rep(list(S2f), 2)
)
bm2f <- c(fs_F1 = "int_fs_F1", fs_F2 = "int_fs_F2")
fit_2f_d <- suppressWarnings(tspa_mx_model("F2 ~ F1\nF2 + F1 ~ 1",
                                           data = fs_mirt2f))
fit_2f_e <- suppressWarnings(tspa_mx_model(
  "F2 ~ F1\nF2 + F1 ~ 1",
  data = fs_indiv(fs_mirt2f, include_intercept = TRUE),
  fsL = Lm2f, fsT = Tm2f, fsb = bm2f
))

# --- 8: mirt SG with a completely-missing row (D6) --------------------------
d_na <- mirt_sim1f(50L)
d_na[1L, ] <- NA
mirt_na <- suppressWarnings(mirt::mirt(d_na, 1L, verbose = FALSE))
fs_na <- get_fs(mirt_na)

test_that("derived mirt SG (mirt_per_obs) equals the manual def-var fit", {
  expect_identical(unname(coef(fit_m_d)), unname(coef(fit_m_e)))
  expect_equal(mx_var_val(fit_m_d, "F1"), mx_var_val(fit_m_e, "F1"),
               tolerance = 1e-10)
  # the per-row fsb list appends int_fs_F1, equal to the fs_indiv column
  der_m <- R2spa:::tspa_mx_derive_measurement(fs_mirt)
  expect_identical(der_m$data[["int_fs_F1"]],
                   fs_indiv(fs_mirt, include_intercept = TRUE)[["int_fs_F1"]])
})

test_that("mirt MG (group column, no group_col attribute) derives as a pooled per-row-corrected fit, equal to the manual def-var fit", {
  # Documented behavior: a MultipleGroupClass result carries a `group`
  # column but no `group_col` attribute, so the D7 refusal does NOT fire.
  # The result derives as a single pooled per-row-corrected fit: the exact
  # per-row measurement quantities, no per-group structural parameters, and
  # the `group` column is inert.
  expect_true(isTRUE(attr(fs_mirt_mg, "mirt_per_obs")))
  expect_null(attr(fs_mirt_mg, "group_col"))
  expect_true("group" %in% names(fs_mirt_mg))
  expect_identical(unname(coef(fit_mg_d)), unname(coef(fit_mg_e)))
  expect_equal(mx_var_val(fit_mg_d, "F1"), mx_var_val(fit_mg_e, "F1"),
               tolerance = 1e-10)
})

test_that("derived mirt 2-factor (per-row, off-diagonal defvars) equals the manual def-var fit", {
  # The vignette's multidimensional case end-to-end (q = 2 per-row). Before
  # the '~~' lookup fix the derived route aborted ("implied covariance not
  # positive definite") with a leaked fixed unit score covariance.
  expect_identical(unname(coef(fit_2f_d)), unname(coef(fit_2f_e)))
  expect_equal(mx_path_val(fit_2f_d, "F1", "F2"),
               mx_path_val(fit_2f_e, "F1", "F2"), tolerance = 1e-10)
  expect_equal(mx_var_val(fit_2f_d, "F1"), mx_var_val(fit_2f_e, "F1"),
               tolerance = 1e-10)
  expect_equal(mx_var_val(fit_2f_d, "F2"), mx_var_val(fit_2f_e, "F2"),
               tolerance = 1e-10)
  expect_equal(vcov(fit_2f_d), vcov(fit_2f_e), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("derived mirt prior_mean (nonzero per-row fsb) equals the explicit int_fs_ fit", {
  dat_pm <- fs_indiv(fs_mirt_pm, include_intercept = TRUE)
  expect_true(any(dat_pm[["int_fs_F1"]] != 0))
  expect_identical(unname(coef(fit_pm_d)), unname(coef(fit_pm_e)))
  expect_equal(mx_var_val(fit_pm_d, "F1"), mx_var_val(fit_pm_e, "F1"),
               tolerance = 1e-10)
})

test_that("a completely-missing row fires the D6 NA guard on the defvar columns", {
  err <- expect_error(
    tspa_mx_model("F1 ~ 1", data = fs_na),
    "Definition-variable column(s) contain NA", fixed = TRUE
  )
  expect_match(conditionMessage(err), "ev_fs_F1", fixed = TRUE)
  expect_match(conditionMessage(err), "int_fs_F1", fixed = TRUE)
})
