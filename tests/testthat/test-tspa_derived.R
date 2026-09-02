# PLAN 13: tspa() auto-derives the measurement inputs (fsT/fsL/fsb or se_fs)
# from a get_fs() result when the caller omits them. Explicit arguments always
# win; a supplied se_fs suppresses the multi-factor derivation; the provenance
# gate (resolve_fs_per_row) rejects hand-rolled attributes; every call that
# reaches the new fail-fast error was erroring before the change.
#
# Convention: every "derived" fit is compared against the explicit-argument
# control fit on the same data. Multi-factor controls pass the FULL explicit
# triple (fsT + fsL + fsb) so both sides carry an identical free-parameter
# set (with group =, lavaan auto-enables a mean structure and a no-fsb form
# estimates the fs_* intercepts freely). The user model string is carried
# verbatim into tspaModel, so derived and control calls always use the same
# model object below.
library(lavaan)
library(lme4)

########## Shared fixtures (file scope) ##########

# --- Single-group single-factor: cbind()'d get_fs() results (no attributes) ---
cfa_ind60 <- "ind60 =~ x1 + x2 + x3"
cfa_dem60 <- "dem60 =~ y1 + y2 + y3 + y4"
fs_ind <- get_fs(PoliticalDemocracy, model = cfa_ind60)
fs_dem <- get_fs(PoliticalDemocracy, model = cfa_dem60)
fs_cb <- cbind(fs_ind, fs_dem)   # cbind() drops the get_fs() attributes

# --- Single-group multi-factor: unified and list format ---
mod3 <- "
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
"
mod_s3 <- "dem60 ~ ind60
           dem65 ~ ind60 + dem60"
fs3 <- get_fs(PoliticalDemocracy, model = mod3, std.lv = TRUE)
fs3l <- get_fs(PoliticalDemocracy, model = mod3, std.lv = TRUE, format = "list")

# A copy of a get_fs() frame with the measurement attributes stripped: the
# "explicit arguments only" control frame (derivation cannot fire on it, so a
# fit on it isolates the explicit-args path).
strip_fs_attrs <- function(d) {
  at <- attributes(d)
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix", "fs_pattern",
               "group_col", "mirt_per_obs", "psi", "alpha")) {
    at[[ak]] <- NULL
  }
  attributes(d) <- at
  d
}
fs3_stripped <- strip_fs_attrs(fs3)

# --- Multigroup multi-factor: unified and list format ---
mod2g <- "
  visual =~ x1 + x2 + x3
  speed  =~ x7 + x8 + x9
"
mod_m <- "visual ~ speed"
fsmg <- get_fs(HolzingerSwineford1939, model = mod2g, std.lv = TRUE,
               group = "school")
fsmgl <- get_fs(HolzingerSwineford1939, model = mod2g, std.lv = TRUE,
                group = "school", format = "list")

# --- Multigroup single-factor: cbind()'d MG frames (no attributes) ---
fs_v <- get_fs(HolzingerSwineford1939, model = "visual =~ x1 + x2 + x3",
               group = "school")
fs_s <- get_fs(HolzingerSwineford1939, model = "speed =~ x7 + x8 + x9",
               group = "school")
fs_hs <- cbind(fs_v, fs_s)

# --- FIML (per-pattern fsT/fsL; the derivation path pools them) ---
hs_fiml <- HolzingerSwineford1939
set.seed(1334)
hs_fiml$x2[!rbinom(nrow(hs_fiml), 1L, 0.4)] <- NA
hs_fiml$x8[!rbinom(nrow(hs_fiml), 1L, 0.4)] <- NA
fit_fiml_sg <- suppressWarnings(
  cfa(mod2g, data = hs_fiml, missing = "fiml")
)
fs_fin <- suppressWarnings(get_fs(fit_fiml_sg))
fit_fiml_mg <- suppressWarnings(
  cfa(mod2g, data = hs_fiml, group = "school", missing = "fiml")
)
fs_fin_mg <- suppressWarnings(get_fs(fit_fiml_mg))

# --- merMod (per-cluster fsT/fsL; no fsb) ---
fs_mer <- get_fs(lmer(Reaction ~ Days + (Days | Subject), sleepstudy))

# --- Corrected-SE stage-1 fits (q = 2 keeps the Jacobian refits cheap) ---
mod2c <- "ind60 =~ x1 + x2 + x3
          dem60 =~ y1 + y2 + y3 + y4"
fs2c <- get_fs(PoliticalDemocracy, model = mod2c, vfsLT = TRUE)
fs2g <- get_fs(HolzingerSwineford1939, model = mod2g, std.lv = TRUE,
               group = "school", vfsLT = TRUE)

########## Derived / control fits (file scope) ##########

# Single-group multi-factor (unified)
fit_sg_d <- tspa(mod_s3, data = fs3)
fit_sg_e <- tspa(mod_s3, data = fs3,
                 fsT = attr(fs3, "fsT"), fsL = attr(fs3, "fsL"),
                 fsb = attr(fs3, "fsb"))
# Single-group multi-factor (list format: a plain data frame with direct
# attributes for a single group)
fit_sgdl_d <- tspa(mod_s3, data = fs3l)
fit_sgdl_e <- tspa(mod_s3, data = fs3l,
                   fsT = attr(fs3l, "fsT"), fsL = attr(fs3l, "fsL"),
                   fsb = attr(fs3l, "fsb"))
# Multigroup multi-factor (unified and list format)
fit_mg_d <- tspa(mod_m, data = fsmg, group = "school")
fit_mg_e <- tspa(mod_m, data = fsmg, group = "school",
                 fsT = attr(fsmg, "fsT"), fsL = attr(fsmg, "fsL"),
                 fsb = attr(fsmg, "fsb"))
fit_mgdl_d <- tspa(mod_m, data = fsmgl, group = "school")

# Single-group single-factor from the cbind()'d frame
se_exp_sg <- list(ind60 = fs_ind$fs_ind60_se[1L],
                  dem60 = fs_dem$fs_dem60_se[1L])
fit_sg_sfd <- tspa("dem60 ~ ind60", data = fs_cb)
fit_sg_sfe <- tspa("dem60 ~ ind60", data = fs_cb, se_fs = se_exp_sg)
# Multigroup single-factor from the cbind()'d MG frames
se_exp_mg <- data.frame(visual = unique(fs_hs$fs_visual_se),
                        speed = unique(fs_hs$fs_speed_se))
fit_mg_sfd <- tspa(mod_m, data = fs_hs, group = "school")
fit_mg_sfe <- tspa(mod_m, data = fs_hs, se_fs = se_exp_mg, group = "school")
# PLAN 17 re-baseline: the derived single-factor fits now carry the
# recovered per-latent implied loading (the <v>_by_fs_<v> columns that
# survive cbind()) instead of a unit loading, so their reference is the
# explicit fsL/fsT (multi-factor) form, not the explicit-se_fs (unit)
# control. The per-construct 1x1 blocks are combined into the 2x2
# (block-diagonal) form the mf schema consumes; the zero off-diagonal
# cells render as extra `+ 0 *` terms in the explicit model (same
# implied covariance, a different optimizer path), so the A/B in the
# tests below is on coef/vcov, not on the model string.
bd2 <- function(a, b) {
  z_r <- `dimnames<-`(matrix(0, 1L, 1L), list(rownames(b), colnames(a)))
  z_c <- `dimnames<-`(matrix(0, 1L, 1L), list(rownames(a), colnames(b)))
  rbind(cbind(a, z_c), cbind(z_r, b))
}
fsL_cb <- bd2(attr(fs_ind, "fsL")[[1L]], attr(fs_dem, "fsL")[[1L]])
fsT_cb <- bd2(attr(fs_ind, "fsT")[[1L]], attr(fs_dem, "fsT")[[1L]])
fit_sg_mf <- tspa("dem60 ~ ind60", data = fs_cb, fsL = fsL_cb, fsT = fsT_cb)
# per-group explicit blocks in first-appearance group order (Pasteur,
# Grant-White), mirroring the derived per-group se/ld rows
grp_ord <- c("Pasteur", "Grant-White")
fsL_hs <- setNames(lapply(grp_ord, function(g) bd2(attr(fs_v, "fsL")[[g]],
                                                   attr(fs_s, "fsL")[[g]])),
                   grp_ord)
fsT_hs <- setNames(lapply(grp_ord, function(g) bd2(attr(fs_v, "fsT")[[g]],
                                                   attr(fs_s, "fsT")[[g]])),
                   grp_ord)
fit_mg_mf <- suppressWarnings(tspa(mod_m, data = fs_hs, group = "school",
                                   fsL = fsL_hs, fsT = fsT_hs))

# FIML (SG and MG); the derived and explicit-triple forms both pool
fit_fiml_sg_d <- suppressWarnings(tspa(mod_m, data = fs_fin))
fit_fiml_sg_e <- suppressWarnings(
  tspa(mod_m, data = fs_fin,
       fsT = attr(fs_fin, "fsT"), fsL = attr(fs_fin, "fsL"),
       fsb = attr(fs_fin, "fsb"))
)
fit_fiml_mg_d <- suppressWarnings(tspa(mod_m, data = fs_fin_mg, group = "school"))
fit_fiml_mg_e <- suppressWarnings(
  tspa(mod_m, data = fs_fin_mg, group = "school",
       fsT = attr(fs_fin_mg, "fsT"), fsL = attr(fs_fin_mg, "fsL"),
       fsb = attr(fs_fin_mg, "fsb"))
)

# merMod (no fsb on either side)
fit_mer_d <- suppressWarnings(tspa("u1 ~ u0", data = fs_mer))
fit_mer_e <- suppressWarnings(
  tspa("u1 ~ u0", data = fs_mer,
       fsT = attr(fs_mer, "fsT"), fsL = attr(fs_mer, "fsL"))
)

# Precedence: deliberately different explicit fsT (scaled) must bypass the
# derivation; the control fits the same explicit args on the attribute-less
# frame.
T_scaled <- lapply(attr(fs3, "fsT"), `*`, 1.5)
fit_prec <- suppressWarnings(tspa(mod_s3, data = fs3,
                                  fsT = T_scaled, fsL = attr(fs3, "fsL")))
fit_prec_ctrl <- suppressWarnings(tspa(mod_s3, data = fs3_stripped,
                                       fsT = T_scaled, fsL = attr(fs3, "fsL")))

# D3: an explicit se_fs on a multi-factor frame keeps the single-factor path
fit_d3 <- suppressWarnings(
  tspa(mod_s3, data = fs3,
       se_fs = list(ind60 = fs3[1L, "fs_ind60_se"],
                    dem60 = fs3[1L, "fs_dem60_se"],
                    dem65 = fs3[1L, "fs_dem65_se"]))
)
fit_d3_ctrl <- suppressWarnings(
  tspa(mod_s3, data = fs3_stripped,
       se_fs = list(ind60 = fs3[1L, "fs_ind60_se"],
                    dem60 = fs3[1L, "fs_dem60_se"],
                    dem65 = fs3[1L, "fs_dem65_se"]))
)

# corrected_se with derived vs explicit measurement inputs (each corrected
# build refits stage 2 once per free fsL/fsT element; built once at file
# scope to bound cost)
fit_corr_sg_d <- tspa("dem60 ~ ind60", data = fs2c, corrected_se = TRUE,
                      vfsLT = attr(fs2c, "vfsLT"))
fit_corr_sg_e <- tspa("dem60 ~ ind60", data = fs2c, corrected_se = TRUE,
                      vfsLT = attr(fs2c, "vfsLT"),
                      fsT = attr(fs2c, "fsT"), fsL = attr(fs2c, "fsL"),
                      fsb = attr(fs2c, "fsb"))
fit_corr_mg_d <- tspa(mod_m, data = fs2g, group = "school",
                      corrected_se = TRUE, vfsLT = attr(fs2g, "vfsLT"))
fit_corr_mg_e <- tspa(mod_m, data = fs2g, group = "school",
                      corrected_se = TRUE, vfsLT = attr(fs2g, "vfsLT"),
                      fsT = attr(fs2g, "fsT"), fsL = attr(fs2g, "fsL"),
                      fsb = attr(fs2g, "fsb"))
# Plain (uncorrected) derived fits: the correction must be non-trivial, and
# the tspa_args record must carry the resolved measurement inputs.
fit_2c_d_plain <- tspa("dem60 ~ ind60", data = fs2c)

# Empty structural model (the efa-score pattern): no model-string parsing
fit_empty_d <- tspa("", data = fs3)
fit_empty_e <- tspa("", data = fs3_stripped,
                    fsT = attr(fs3, "fsT"), fsL = attr(fs3, "fsL"),
                    fsb = attr(fs3, "fsb"))

# Provenance-gate fixture: a plain data frame with score columns and plain
# matrix fsT/fsL attributes, but no get_fs() provenance (no fs_pattern).
d_hand <- data.frame(fs_a = c(0, 1, 2, 3, 4, 5), fs_b = c(3, 1, 2, 0, 5, 4))
attr(d_hand, "fsT") <- matrix(c(0.2, 0.1, 0.1, 0.3), 2L, 2L,
                              dimnames = list(c("fs_a", "fs_b"),
                                              c("fs_a", "fs_b")))
attr(d_hand, "fsL") <- `dimnames<-`(diag(2L),
                                    list(c("fs_a", "fs_b"), c("a", "b")))

########## Tests ##########

test_that("derived SG multi-factor (unified) equals the explicit full-triple fit", {
  expect_identical(attr(fit_sg_d, "tspaModel"), attr(fit_sg_e, "tspaModel"))
  expect_equal(coef(fit_sg_d), coef(fit_sg_e), tolerance = 1e-10)
  expect_equal(vcov(fit_sg_d), vcov(fit_sg_e), tolerance = 1e-10,
               ignore_attr = TRUE)
  # The fit carries the (post-pooling, here identity) measurement inputs
  expect_equal(attr(fit_sg_d, "fsT"), attr(fs3, "fsT"), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(attr(fit_sg_d, "fsL"), attr(fs3, "fsL"), tolerance = 1e-10,
               ignore_attr = TRUE)
  # Complete data: no pooling happened on the derived path
  expect_null(attr(fit_sg_d, "pooled_fs"))
})

test_that("derived SG multi-factor (list format) equals the explicit full-triple fit", {
  # SG list format is a plain data frame with direct (unwrapped) attributes
  expect_identical(attr(fit_sgdl_d, "tspaModel"), attr(fit_sgdl_e, "tspaModel"))
  expect_equal(coef(fit_sgdl_d), coef(fit_sgdl_e), tolerance = 1e-10)
  expect_equal(vcov(fit_sgdl_d), vcov(fit_sgdl_e), tolerance = 1e-10,
               ignore_attr = TRUE)
  # list and unified formats derive to the same stage-2 model
  expect_identical(attr(fit_sgdl_d, "tspaModel"), attr(fit_sg_d, "tspaModel"))
})

test_that("derived MG multi-factor (unified and list) equals the explicit full-triple fit", {
  expect_identical(attr(fit_mg_d, "tspaModel"), attr(fit_mg_e, "tspaModel"))
  expect_equal(coef(fit_mg_d), coef(fit_mg_e), tolerance = 1e-10)
  expect_equal(vcov(fit_mg_d), vcov(fit_mg_e), tolerance = 1e-10,
               ignore_attr = TRUE)
  # stage-2 group order: first appearance of the group column
  expect_equal(names(attr(fit_mg_d, "fsT")), c("Pasteur", "Grant-White"))
  expect_identical(attr(fit_mgdl_d, "tspaModel"), attr(fit_mg_e, "tspaModel"))
  expect_equal(coef(fit_mgdl_d), coef(fit_mg_e), tolerance = 1e-10)
})

test_that("derived SG single-factor se_fs from cbind()'d frames equals the explicit fsL/fsT fit", {
  se_der <- attr(fit_sg_sfd, "tspa_args")$se_fs
  # the derived values are the (row-constant) fs_<v>_se columns, in
  # first-appearance (cbind) order of the score columns
  expect_equal(se_der, as.data.frame(as.list(se_exp_sg)))
  # house-pinned (rounded) values from the tspa() examples
  expect_equal(unname(se_der$ind60), 0.1213615, tolerance = 1e-6)
  expect_equal(unname(se_der$dem60), 0.6756472, tolerance = 1e-6)
  # PLAN 17: the derived fit carries the recovered per-latent implied
  # loading (the <v>_by_fs_<v> columns, 1 - Vpost/psi; here the fsL
  # diagonals), not a unit loading, so the reference is the explicit
  # fsL/fsT (multi-factor) form, not the unit-loading explicit-se_fs
  # control. Cross-form A/B: the explicit mf rendering carries the zero
  # off-diagonal score terms (same implied covariance, a different
  # optimizer path), so the MLEs agree to ~1e-8 relative, not
  # bit-for-bit.
  expect_equal(coef(fit_sg_sfd), coef(fit_sg_mf), tolerance = 1e-6)
  expect_equal(vcov(fit_sg_sfd), vcov(fit_sg_mf), tolerance = 1e-6,
               ignore_attr = TRUE)
  # ...and it differs from the pre-PLAN 17 unit-loading fit
  expect_false(isTRUE(all.equal(coef(fit_sg_sfd), coef(fit_sg_sfe),
                                tolerance = 1e-6)))
  # house-pinned (rounded) structural coefficient
  expect_equal(unname(coef(fit_sg_sfd)["dem60~ind60"]), 1.4478669,
               tolerance = 1e-6)
  # golden: the derived model string carries the recovered (non-unit)
  # per-latent loadings
  expect_identical(attr(fit_sg_sfd, "tspaModel"), paste0(
    "# latent variables (indicated by factor scores)\n",
    "ind60 =~ 0.965767270434244 * fs_ind60\n",
    "dem60 =~ 0.886804906625876 * fs_dem60\n",
    "\n",
    "# constrain the errors\n",
    "fs_ind60 ~~ 0.0147286194470876 * fs_ind60\n",
    "fs_dem60 ~~ 0.456499146148539 * fs_dem60\n",
    "\n",
    "# structural model\n",
    "dem60 ~ ind60"))
})

test_that("derived MG single-factor se_fs from cbind()'d MG frames equals the explicit per-group fsL/fsT fit", {
  se_der <- attr(fit_mg_sfd, "tspa_args")$se_fs
  # one row per group in first-appearance order: Pasteur then Grant-White
  expect_equal(rownames(se_der), c("Pasteur", "Grant-White"))
  expect_equal(unname(se_der$visual), unname(unique(fs_hs$fs_visual_se)))
  expect_equal(unname(se_der$speed), unname(unique(fs_hs$fs_speed_se)))
  # house-pinned (rounded) values from the tspa() examples, in group order
  expect_equal(unname(se_der$visual), c(0.3391326, 0.3118280), tolerance = 1e-6)
  expect_equal(unname(se_der$speed), c(0.2786875, 0.2740507), tolerance = 1e-6)
  # PLAN 17: the recovered per-group loadings make the derived fit equal
  # the explicit per-group fsL/fsT reference (cross-form A/B, tolerance as
  # in the SG test), not the unit-loading explicit-se_fs control
  expect_equal(coef(fit_mg_sfd), coef(fit_mg_mf), tolerance = 1e-6)
  expect_equal(vcov(fit_mg_sfd), vcov(fit_mg_mf), tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(coef(fit_mg_sfd), coef(fit_mg_sfe),
                                tolerance = 1e-6)))
  # house-pinned (rounded) structural coefficient
  expect_equal(unname(coef(fit_mg_sfd)["visual~speed"]), 0.3398422,
               tolerance = 1e-6)
  # golden: the per-group loadings are in first-appearance group order
  expect_identical(attr(fit_mg_sfd, "tspaModel"), paste0(
    "# latent variables (indicated by factor scores)\n",
    "visual =~ c(0.673482617533149, 0.699050944243108) * fs_visual\n",
    "speed =~ c(0.662316098427474, 0.812759419528942) * fs_speed\n",
    "\n",
    "# constrain the errors\n",
    "fs_visual ~~ c(0.115010888044477, 0.0972366713691954) * fs_visual\n",
    "fs_speed ~~ c(0.0776667057991782, 0.0751037642308851) * fs_speed\n",
    "\n",
    "# structural model\n",
    "visual ~ speed"))
})

test_that("derived FIML (SG and MG) equals the explicit full-triple fit; replay is self-contained", {
  # SG: the per-pattern attributes are pooled (reduce = "mean") on the
  # derived path exactly as on the explicit path
  expect_identical(attr(fit_fiml_sg_d, "tspaModel"),
                   attr(fit_fiml_sg_e, "tspaModel"))
  expect_equal(coef(fit_fiml_sg_d), coef(fit_fiml_sg_e), tolerance = 1e-10)
  expect_identical(attr(fit_fiml_sg_d, "pooled_fs"), "mean")
  # MG: one pooled plain matrix per group
  expect_identical(attr(fit_fiml_mg_d, "tspaModel"),
                   attr(fit_fiml_mg_e, "tspaModel"))
  expect_equal(coef(fit_fiml_mg_d), coef(fit_fiml_mg_e), tolerance = 1e-10)
  expect_identical(attr(fit_fiml_mg_d, "pooled_fs"), "mean")
  # tspa_args carries the RESOLVED (pooled) values: replay re-passes them
  # explicitly, so no re-derivation and no double-pooling (the pooled_fs
  # marker is not re-attached on replay)
  replay <- suppressWarnings(do.call(tspa, attr(fit_fiml_sg_d, "tspa_args")))
  expect_equal(coef(replay), coef(fit_fiml_sg_d), tolerance = 1e-10)
  expect_equal(vcov(replay), vcov(fit_fiml_sg_d), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_null(attr(replay, "pooled_fs"))
})

test_that("derived merMod equals the explicit fsT/fsL fit (no fsb)", {
  expect_identical(attr(fit_mer_d, "tspaModel"), attr(fit_mer_e, "tspaModel"))
  expect_equal(coef(fit_mer_d), coef(fit_mer_e), tolerance = 1e-10)
  expect_identical(attr(fit_mer_d, "pooled_fs"), "mean")
  # merMod carries no fsb; the derived path must not invent one
  expect_null(attr(fit_mer_d, "tspa_args")$fsb)
})

test_that("explicit fsT/fsL take precedence over the data's attributes", {
  # same fit as the explicit-args-only control (attributes irrelevant)
  expect_identical(attr(fit_prec, "tspaModel"),
                   attr(fit_prec_ctrl, "tspaModel"))
  expect_equal(coef(fit_prec), coef(fit_prec_ctrl), tolerance = 1e-10)
  # the scaled explicit value (not the attribute) is what got fitted
  expect_equal(attr(fit_prec, "fsT")[[1L]], T_scaled[[1L]],
               tolerance = 1e-10, ignore_attr = TRUE)
  # ...so the scaled fit differs from the plain derived (unscaled) fit
  expect_false(isTRUE(all.equal(coef(fit_prec), coef(fit_sg_d),
                                tolerance = 1e-8)))
})

test_that("a supplied se_fs suppresses the multi-factor derivation (D3)", {
  # single-factor path: no latent-variable line carries a second score term
  lines3 <- strsplit(attr(fit_d3, "tspaModel"), "\n")[[1L]]
  expect_false(any(grepl("=~", lines3) & grepl("\\+", lines3)))
  # identical to today's explicit behavior on an attribute-less frame
  expect_identical(attr(fit_d3, "tspaModel"), attr(fit_d3_ctrl, "tspaModel"))
  expect_equal(coef(fit_d3), coef(fit_d3_ctrl), tolerance = 1e-10)
})

test_that("provenance gate and fail-fast errors name the accepted input forms", {
  # hand-rolled fsT/fsL attributes without get_fs() provenance are not
  # derived; the fail-fast error carries the gate's message
  err <- expect_error(tspa("b ~ a", data = d_hand))
  msg <- conditionMessage(err)
  expect_match(msg, "No measurement inputs found", fixed = TRUE)
  expect_match(msg, "does not look like a get_fs() result", fixed = TRUE)
  expect_match(msg, "has no 'fs_pattern' attribute", fixed = TRUE)
  # nothing available at all: the base fail-fast message (no gate addendum)
  err2 <- expect_error(tspa("y ~ x", data = data.frame(x = 1:5, y = 1:5)))
  msg2 <- conditionMessage(err2)
  expect_match(msg2, "No measurement inputs found", fixed = TRUE)
  expect_false(grepl("does not look like a get_fs() result", msg2,
                     fixed = TRUE))
  # an empty-but-non-NULL se_fs suppresses the column derivation and hits
  # the same fail-fast error (not lavaan's "model is NULL")
  err3 <- expect_error(
    tspa("y ~ x", data = data.frame(x = 1:5, y = 1:5), se_fs = list()))
  msg3 <- conditionMessage(err3)
  expect_match(msg3, "No measurement inputs found", fixed = TRUE)
  expect_false(grepl("model is NULL", msg3, fixed = TRUE))
  # the pre-existing both-or-neither check fires before any derivation
  expect_error(tspa(mod_s3, data = fs3, fsT = attr(fs3, "fsT")),
               "Please provide both or none of fsT and fsL.")
})

test_that("corrected_se with derived fsT/fsL equals the explicit-triple corrected fit", {
  # SG
  expect_equal(vcov(fit_corr_sg_d), vcov(fit_corr_sg_e), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_true(isTRUE(attr(fit_corr_sg_d, "tspa_corrected")) &&
                isTRUE(attr(fit_corr_sg_e, "tspa_corrected")))
  # the correction is non-trivial (a no-op correction would pass the
  # equality above silently)
  expect_false(isTRUE(all.equal(vcov(fit_corr_sg_d), vcov(fit_2c_d_plain),
                                tolerance = 1e-8)))
  # MG
  expect_equal(vcov(fit_corr_mg_d), vcov(fit_corr_mg_e), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_true(isTRUE(attr(fit_corr_mg_d, "tspa_corrected")) &&
                isTRUE(attr(fit_corr_mg_e, "tspa_corrected")))
  # derived tspa_args are fully resolved (self-contained replay record)
  args_mf <- attr(fit_2c_d_plain, "tspa_args")
  expect_false(is.null(args_mf$fsT))
  expect_false(is.null(args_mf$fsL))
  args_sf <- attr(fit_sg_sfd, "tspa_args")
  expect_false(is.null(args_sf$se_fs))
  expect_gte(nrow(args_sf$se_fs), 1L)
})

test_that("model = '' (efa-score pattern) with derived attrs equals the explicit form", {
  expect_identical(attr(fit_empty_d, "tspaModel"),
                   attr(fit_empty_e, "tspaModel"))
  expect_equal(coef(fit_empty_d), coef(fit_empty_e), tolerance = 1e-10)
  expect_equal(vcov(fit_empty_d), vcov(fit_empty_e), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("derived per-group se_fs order follows first appearance, not factor level order", {
  # deliberately reversed factor levels (first appearance: Pasteur,
  # Grant-White; level order: Grant-White, Pasteur)
  d_rev <- fs_hs
  d_rev[["school"]] <- factor(d_rev[["school"]],
                              levels = c("Grant-White", "Pasteur"))
  se_der <- R2spa:::derive_sf_se_fs(d_rev, "school", "mean")
  # PLAN 17: derive_sf_se_fs() now returns list(se, ld); the se block's
  # group order is unchanged
  expect_equal(rownames(se_der$se), c("Pasteur", "Grant-White"))
  # the derived fit's per-group rows (se and the PLAN 17 recovered
  # loadings) are in first-appearance order: identical to the derived
  # fit on the unreversed frame ...
  fit_rd <- suppressWarnings(tspa(mod_m, data = d_rev, group = "school"))
  expect_identical(attr(fit_rd, "tspaModel"),
                   attr(fit_mg_sfd, "tspaModel"))
  # ... and not the level-order one
  se_lvl <- data.frame(
    visual = c(unique(fs_hs$fs_visual_se)[2L], unique(fs_hs$fs_visual_se)[1L]),
    speed = c(unique(fs_hs$fs_speed_se)[2L], unique(fs_hs$fs_speed_se)[1L]))
  fit_rl <- suppressWarnings(tspa(mod_m, data = d_rev, se_fs = se_lvl,
                                  group = "school"))
  expect_false(identical(attr(fit_rd, "tspaModel"),
                         attr(fit_rl, "tspaModel")))
})

# ---- mirt (per-observation fsT/fsL; the derivation path pools them) ----

skip_if_not_installed("mirt")

set.seed(2025)
NMF <- 120L
mrt_sim2f <- function(N) as.data.frame(mirt::simdata(
  a = matrix(c(runif(4L, 0.5, 1.5), runif(4L, 0.5, 1.5)), 8L, 2L),
  d = rnorm(8L), N = N, itemtype = "2PL",
  Theta = cbind(rnorm(N), rnorm(N))))

dat_mf_sg <- mrt_sim2f(NMF)
dat_mf_sg[1L, ] <- NA
mf_sg <- suppressWarnings(mirt::mirt(dat_mf_sg, 2L, invariance = "slopes",
                                     verbose = FALSE))
fs_mf_sg <- get_fs(mf_sg)

dat_mf_mg <- rbind(mrt_sim2f(NMF), mrt_sim2f(NMF))
dat_mf_mg[1L, ] <- NA
grp_mf <- factor(rep(c("A", "B"), each = NMF))
mf_mg <- suppressWarnings(mirt::multipleGroup(dat_mf_mg, 2L, group = grp_mf,
                                              invariance = "slopes",
                                              verbose = FALSE))
fs_mf_mg <- get_fs(mf_mg)

# The mirt result carries per-row fsb (pooled to a constant), so the
# explicit control must pass the full triple to match the derived fit's
# parameter set.
fit_mirt_sg_d <- suppressWarnings(tspa("F2 ~ F1", data = fs_mf_sg))
fit_mirt_sg_e <- suppressWarnings(
  tspa("F2 ~ F1", data = fs_mf_sg,
       fsT = attr(fs_mf_sg, "fsT"), fsL = attr(fs_mf_sg, "fsL"),
       fsb = attr(fs_mf_sg, "fsb"))
)
fit_mirt_mg_d <- suppressWarnings(tspa("F2 ~ F1", data = fs_mf_mg,
                                       group = "group"))
fit_mirt_mg_e <- suppressWarnings(
  tspa("F2 ~ F1", data = fs_mf_mg, group = "group",
       fsT = attr(fs_mf_mg, "fsT"), fsL = attr(fs_mf_mg, "fsL"),
       fsb = attr(fs_mf_mg, "fsb"))
)

test_that("derived mirt (SG and MG) equals the explicit full-triple fit", {
  expect_identical(attr(fit_mirt_sg_d, "tspaModel"),
                   attr(fit_mirt_sg_e, "tspaModel"))
  expect_equal(coef(fit_mirt_sg_d), coef(fit_mirt_sg_e), tolerance = 1e-10)
  expect_equal(vcov(fit_mirt_sg_d), vcov(fit_mirt_sg_e), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_identical(attr(fit_mirt_mg_d, "tspaModel"),
                   attr(fit_mirt_mg_e, "tspaModel"))
  expect_equal(coef(fit_mirt_mg_d), coef(fit_mirt_mg_e), tolerance = 1e-10)
  expect_equal(names(attr(fit_mirt_mg_d, "fsT")), c("A", "B"))
})
