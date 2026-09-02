# PLAN 17: when se_fs is auto-derived from a cbind()'d get_fs() result
# (no fsL/fsT attributes left after cbind(), and no explicit se_fs), the
# per-latent loading of the stage-2 model is recovered from the
# <v>_by_fs_<v> implied-loading columns (which survive cbind()) instead
# of the pre-PLAN 17 hard-coded unit loading. The contracts that do not
# change: an explicit se_fs still means a unit loading, and a non-cbind'd
# get_fs() result (attributes intact) still takes the multi-factor
# (fsL/fsT) path. The derivation cases below are A/B-gated against the
# explicit fsL/fsT form (not just pinned).
library(lavaan)
library(lme4)

########## Fixtures (file scope) ##########

# Single random-intercept lme4 DGP (balanced: every cluster has n_i
# observations, so the per-cluster fsL/fsT blocks are identical and the
# multi-factor pooling, reduce = "mean", reduces to the constant block).
# seed 6000 / 200 clusters: the derived structural coefficient lands
# 0.01 SE from the single-stage full-SEM reference (house-pinned in T2).
set.seed(6000L)
n_j <- 200L
n_i <- 5L
u_j <- rnorm(n_j, 0, 1.5)
y_mat <- sweep(matrix(rnorm(n_j * n_i), nrow = n_j, ncol = n_i), 1L, u_j, "+")
d_j <- 0.8 * u_j + rnorm(n_j, 0, 0.5)
wide <- as.data.frame(y_mat)
names(wide) <- paste0("y", seq_len(n_i))
wide$d <- d_j
long <- data.frame(j = rep(seq_len(n_j), each = n_i), y = as.vector(t(y_mat)))
lmer_fit <- lmer(y ~ 1 + (1 | j), long)
fs_lmer <- get_fs(lmer_fit)
fsd_lmer <- cbind(fs_lmer, d = d_j)
psi_lmer <- unname(as.numeric(VarCorr(lmer_fit)[[1L]]))

# T1 reference: the explicit fsL/fsT form on the cbind()'d frame. The
# balanced DGP's per-cluster blocks are identical, so the mf-path pooling
# (reduce = "mean") is the constant block itself, passed as a plain 1x1
# (single-group) matrix. The raw 3-D per-cluster attributes cannot be
# passed on the cbind()'d frame: the per-unit pooling needs the get_fs()
# row provenance (fs_pattern) that cbind() drops.
L1 <- `dimnames<-`(matrix(attr(fs_lmer, "fsL")[1L, 1L, 1L], 1L, 1L),
                   list("fs_u0", "u0"))
T1 <- `dimnames<-`(matrix(attr(fs_lmer, "fsT")[1L, 1L, 1L], 1L, 1L),
                   list("fs_u0", "fs_u0"))
fit_lmer_der  <- tspa("d ~ u0", data = fsd_lmer)
fit_lmer_exp  <- tspa("d ~ u0", data = fsd_lmer, fsL = L1, fsT = T1)
fit_lmer_unit <- tspa("d ~ u0", data = fsd_lmer,
                      se_fs = list(u0 = mean(fs_lmer$fs_u0_se)))
# T2 single-stage full-SEM reference on the original indicators
sem_lmer <- sem(paste0("u0 =~ ",
                       paste(paste0("y", seq_len(n_i)), collapse = " + "),
                       "\nd ~ u0"),
                data = wide)
# T3 latent-variance recovery (derived vs the unit-loading contract)
fit_lmer_var      <- tspa("u0 ~~ u0", data = fsd_lmer)
fit_lmer_var_unit <- tspa("u0 ~~ u0", data = fsd_lmer,
                          se_fs = list(u0 = mean(fs_lmer$fs_u0_se)))
# T6(b) non-cbind'd get_fs() result: attributes intact, multi-factor path
fit_lmer_mf     <- tspa("u0 ~~ u0", data = fs_lmer)
fit_lmer_mf_exp <- tspa("u0 ~~ u0", data = fs_lmer,
                        fsL = attr(fs_lmer, "fsL"), fsT = attr(fs_lmer, "fsT"))

# T4/T5 hand-rolled single-factor frames (no attributes at all)
set.seed(42L)
n_hand <- 120L
fs_vv <- rnorm(n_hand)
d_vv <- 0.5 * fs_vv + rnorm(n_hand)
hand_unit_ld <- data.frame(fs_v = fs_vv, fs_v_se = 0.25, v_by_fs_v = 1,
                           d = d_vv)
hand_no_ld <- data.frame(fs_v = fs_vv, fs_v_se = 0.25, d = d_vv)
fit_hand_u_d <- tspa("d ~ v", data = hand_unit_ld)
fit_hand_u_e <- tspa("d ~ v", data = hand_unit_ld, se_fs = list(v = 0.25))
fit_hand_n_d <- tspa("d ~ v", data = hand_no_ld)
# T8: a non-numeric (character / factor) <v>_by_fs_<v> column is treated the
# same as a missing one (unit loading, no error in the derivation path)
hand_chr_ld <- data.frame(fs_v = fs_vv, fs_v_se = 0.25,
                          v_by_fs_v = "not-numeric", d = d_vv)
hand_fct_ld <- data.frame(fs_v = fs_vv, fs_v_se = 0.25,
                          v_by_fs_v = factor("x"), d = d_vv)
fit_hand_c_d <- tspa("d ~ v", data = hand_chr_ld)
fit_hand_f_d <- tspa("d ~ v", data = hand_fct_ld)
# T9: an all-NA numeric <v>_by_fs_<v> column -> NaN/NA pool -> unit loading
hand_nan_ld <- data.frame(fs_v = fs_vv, fs_v_se = 0.25,
                          v_by_fs_v = NA_real_, d = d_vv)
fit_hand_nan_d <- tspa("d ~ v", data = hand_nan_ld)

# T7 multigroup cbind()'d get_fs() result with deliberately reversed
# factor levels (first appearance: Pasteur, Grant-White; level order:
# Grant-White, Pasteur)
fs_v_mg <- get_fs(HolzingerSwineford1939, model = "visual =~ x1 + x2 + x3",
                  group = "school")
fs_s_mg <- get_fs(HolzingerSwineford1939, model = "speed =~ x7 + x8 + x9",
                  group = "school")
fs_hs_mg <- cbind(fs_v_mg, fs_s_mg)
fs_hs_rev <- fs_hs_mg
fs_hs_rev[["school"]] <- factor(fs_hs_rev[["school"]],
                                levels = c("Grant-White", "Pasteur"))
fit_mg_ord <- suppressWarnings(tspa("visual ~ speed", data = fs_hs_rev,
                                    group = "school"))

########## Tests ##########

test_that("T1: derived structural coefficient equals the explicit fsL/fsT form (A/B)", {
  # precondition: the balanced DGP's per-cluster blocks are identical, so
  # the mf-path pooling is the constant block L1/T1 above
  expect_equal(length(unique(as.numeric(attr(fs_lmer, "fsT")))), 1L)
  expect_equal(length(unique(as.numeric(attr(fs_lmer, "fsL")))), 1L)
  # A/B gate: both sides reach the same MLE
  expect_equal(coef(fit_lmer_der)[["d~u0"]], coef(fit_lmer_exp)[["d~u0"]],
               tolerance = 1e-8)
  expect_equal(coef(fit_lmer_der), coef(fit_lmer_exp), tolerance = 1e-8)
  expect_equal(vcov(fit_lmer_der), vcov(fit_lmer_exp), tolerance = 1e-8,
               ignore_attr = TRUE)
  # the derived model carries the recovered (non-unit) loading: the
  # u0_by_fs_u0 column value (the fsL value), not 1
  mdl <- strsplit(attr(fit_lmer_der, "tspaModel"), "\n")[[1L]]
  expect_false(any(grepl("^u0 =~ 1 \\* fs_u0$", mdl)))
  ld_line <- grep("^u0 =~ ", mdl, value = TRUE)[1L]
  ld_val <- as.numeric(sub("^u0 =~ (\\S+) \\* fs_u0$", "\\1", ld_line))
  expect_equal(ld_val, attr(fs_lmer, "fsL")[1L, 1L, 1L], tolerance = 1e-12)
})

test_that("T2: derived structural coefficient is within 0.15 SE of the single-stage full-SEM reference", {
  b_der <- coef(fit_lmer_der)[["d~u0"]]
  b_sem <- coef(sem_lmer)[["d~u0"]]
  se_der <- sqrt(vcov(fit_lmer_der)[["d~u0", "d~u0"]])
  # the +21-28% unit-loading bias of the pre-PLAN 17 hard-coded unit
  # loading is gone: the derived (recovered-loading) fit agrees with the
  # single-stage full-SEM reference to a fraction of a SE (0.01 SE for
  # this seed)
  expect_lt(abs(b_der - b_sem), 0.15 * se_der)
  # house-pinned (rounded) corrected value
  expect_equal(unname(b_der), 0.7879429, tolerance = 1e-6)
  # the unit-loading footgun (explicit se_fs) is still measurably biased
  # on the same data: 2.40 SE from the single-stage reference
  b_unit <- coef(fit_lmer_unit)[["d~u0"]]
  expect_gt(abs(b_unit - b_sem), 0.15 * se_der)
})

test_that("T3: derived latent-variance model recovers the lme4 random-effect variance", {
  psi2 <- coef(fit_lmer_var)[["u0~~u0"]]
  se2 <- sqrt(vcov(fit_lmer_var)[["u0~~u0", "u0~~u0"]])
  # statistical recovery: the stage-2 latent variance is within 1 SE of
  # the lme4 (REML) random-effect variance. The two are different sample
  # statistics (stage 2 re-estimates from the BLUP second moment), so
  # exact equality is not expected at finite n; the seeded gap here is
  # 0.05 SE.
  expect_lt(abs(psi2 - psi_lmer), se2)
  # house-pinned (rounded) recovered value
  expect_equal(unname(psi2), 1.785505, tolerance = 1e-6)
  # the unit-loading contract (explicit se_fs) does NOT recover it (2.03
  # SE from the lme4 variance): the recovery is what the recovered
  # loading does
  psi2_unit <- coef(fit_lmer_var_unit)[["u0~~u0"]]
  expect_gt(abs(psi2_unit - psi_lmer), abs(psi2 - psi_lmer))
})

test_that("T4: a unit implied-loading column is a no-op (derived == explicit se_fs)", {
  expect_identical(attr(fit_hand_u_d, "tspaModel"),
                   attr(fit_hand_u_e, "tspaModel"))
  expect_equal(coef(fit_hand_u_d), coef(fit_hand_u_e), tolerance = 1e-10)
  mdl <- strsplit(attr(fit_hand_u_d, "tspaModel"), "\n")[[1L]]
  expect_true(any(grepl("^v =~ 1 \\* fs_v$", mdl)))
})

test_that("T5: no <v>_by_fs_<v> column falls back to a unit loading", {
  expect_identical(attr(fit_hand_n_d, "tspaModel"),
                   attr(fit_hand_u_e, "tspaModel"))
  expect_equal(coef(fit_hand_n_d), coef(fit_hand_u_e), tolerance = 1e-10)
  mdl <- strsplit(attr(fit_hand_n_d, "tspaModel"), "\n")[[1L]]
  expect_true(any(grepl("^v =~ 1 \\* fs_v$", mdl)))
})

test_that("T6: an explicit se_fs keeps the unit-loading contract; a non-cbind'd result uses fsL", {
  # (a) explicit se_fs on a cbind()'d get_fs() result: unit loading, not
  # the recovered one
  mdl_u <- strsplit(attr(fit_lmer_unit, "tspaModel"), "\n")[[1L]]
  expect_true(any(grepl("^u0 =~ 1 \\* fs_u0$", mdl_u)))
  expect_false(identical(attr(fit_lmer_unit, "tspaModel"),
                         attr(fit_lmer_der, "tspaModel")))
  expect_false(isTRUE(all.equal(coef(fit_lmer_unit), coef(fit_lmer_der),
                                tolerance = 1e-6)))
  # (b) non-cbind'd get_fs() result (attributes intact, no se_fs): the
  # multi-factor path uses the fsL attribute, unchanged
  expect_identical(attr(fit_lmer_mf, "tspaModel"),
                   attr(fit_lmer_mf_exp, "tspaModel"))
  expect_equal(coef(fit_lmer_mf), coef(fit_lmer_mf_exp), tolerance = 1e-10)
})

test_that("T7: derived per-group loadings follow first-appearance group order", {
  lines <- strsplit(attr(fit_mg_ord, "tspaModel"), "\n")[[1L]]
  vis <- grep("^visual =~ ", lines, value = TRUE)[1L]
  # first appearance in the data: Pasteur, Grant-White (the reversed
  # factor's level order is Grant-White, Pasteur); the per-group loading
  # values must be in first-appearance order
  ld_p <- attr(fs_v_mg, "fsL")[["Pasteur"]][1L, 1L]
  ld_g <- attr(fs_v_mg, "fsL")[["Grant-White"]][1L, 1L]
  expect_equal(vis,
               paste0("visual =~ c(", as.character(ld_p), ", ",
                      as.character(ld_g), ") * fs_visual"))
  expect_false(identical(vis,
                         paste0("visual =~ c(", as.character(ld_g), ", ",
                                as.character(ld_p), ") * fs_visual")))
  # the loading rows share the derived per-group SE order
  se_der <- attr(fit_mg_ord, "tspa_args")$se_fs
  expect_equal(rownames(se_der), c("Pasteur", "Grant-White"))
})

test_that("T8: a non-numeric <v>_by_fs_<v> column falls back to a unit loading (no error)", {
  # a mistyped (character or factor) implied-loading column is treated the
  # same as a missing one: the derived fit uses a unit loading and equals the
  # explicit-se_fs unit-loading reference, rather than erroring during the
  # auto-derivation path
  expect_identical(attr(fit_hand_c_d, "tspaModel"),
                   attr(fit_hand_u_e, "tspaModel"))
  expect_equal(coef(fit_hand_c_d), coef(fit_hand_u_e), tolerance = 1e-10)
  expect_identical(attr(fit_hand_f_d, "tspaModel"),
                   attr(fit_hand_u_e, "tspaModel"))
  expect_equal(coef(fit_hand_f_d), coef(fit_hand_u_e), tolerance = 1e-10)
  mdl <- strsplit(attr(fit_hand_c_d, "tspaModel"), "\n")[[1L]]
  expect_true(any(grepl("^v =~ 1 \\* fs_v$", mdl)))
})

test_that("T9: an all-NA numeric <v>_by_fs_<v> column falls back to a unit loading (no NaN)", {
  # a numeric implied-loading column that is entirely NA pools to NaN/NA; the
  # non-finite value is coerced back to the unit loading, not rendered into the
  # lavaan syntax
  expect_identical(attr(fit_hand_nan_d, "tspaModel"),
                   attr(fit_hand_u_e, "tspaModel"))
  expect_equal(coef(fit_hand_nan_d), coef(fit_hand_u_e), tolerance = 1e-10)
  mdl <- strsplit(attr(fit_hand_nan_d, "tspaModel"), "\n")[[1L]]
  expect_true(any(grepl("^v =~ 1 \\* fs_v$", mdl)))
})
