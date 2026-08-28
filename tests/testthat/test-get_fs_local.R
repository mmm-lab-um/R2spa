# PLAN 14: get_fs(..., local = TRUE) per-construct ("local") stage-1 scoring.
#
# Each latent is scored from its own local single-factor measurement model
# (the canonical per-construct 2S-PA stage 1) instead of one joint
# multi-factor model. The merged output reproduces the joint layout
# (columns, attribute shapes) with exactly-zero cross terms (block-diagonal
# fsT/fsL/psi, zero off-diagonal _by_ and ecov_ columns) and is
# downstream-transparent: it feeds tspa() directly (no explicit fsT/fsL),
# works through fs_indiv() and fs_to_group_list(), and accepts the merged
# fsL in tspa_mx_model().
#
# Convention: every derived/local result is compared against an independent
# control built from the per-local single-factor get_fs() calls (the cbind
# vignette workflow). Multi-factor controls pass the full explicit triple
# (fsT + fsL + fsb) so both sides carry an identical free-parameter set.
# The A/B (local vs joint) bridge uses tolerance 1e-4 (the uncorrelated
# factorization agrees to optimizer tolerance, ~1e-5, not 1e-6).
library(lavaan)

# ---------------------------------------------------------------------------
# Shared fixtures (file scope)
# ---------------------------------------------------------------------------

pd <- PoliticalDemocracy
cols3 <- c("x1", "x2", "x3", "y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8")
cols2 <- c("x1", "x2", "x3", "y1", "y2", "y3", "y4")
model3 <- "ind60 =~ x1 + x2 + x3
           dem60 =~ y1 + y2 + y3 + y4
           dem65 =~ y5 + y6 + y7 + y8"
# the same three factors constrained uncorrelated: the likelihood
# factorizes, so each factor's MLE coincides with its local fit
model3_unc <- "ind60 =~ x1 + x2 + x3
               dem60 =~ y1 + y2 + y3 + y4
               dem65 =~ y5 + y6 + y7 + y8
               ind60 ~~ 0 * dem60
               ind60 ~~ 0 * dem65
               dem60 ~~ 0 * dem65"
mod2 <- "ind60 =~ x1 + x2 + x3
         dem60 =~ y1 + y2 + y3 + y4"

# --- Single-group: local, joint, uncorrelated-joint, and the per-local cbind ---
fs_local3 <- get_fs(pd[cols3], model3, local = TRUE)
fs_joint3 <- get_fs(pd[cols3], model3)
fs_unc3   <- get_fs(pd[cols3], model3_unc)
# the canonical per-construct pattern (the vignette workflow)
fs_p1 <- get_fs(pd[c("x1", "x2", "x3")], "ind60 =~ x1 + x2 + x3")
fs_p2 <- get_fs(pd[c("y1", "y2", "y3", "y4")], "dem60 =~ y1 + y2 + y3 + y4")
fs_p3 <- get_fs(pd[c("y5", "y6", "y7", "y8")], "dem65 =~ y5 + y6 + y7 + y8")
cb3 <- cbind(fs_p1, fs_p2, fs_p3)
score3 <- paste0("fs_", c("ind60", "dem60", "dem65"))
se3    <- paste0(score3, "_se")

# --- Two-factor SG (keeps corrected_fsT / std.lv fits cheap) ---
fit_pc1 <- cfa("ind60 =~ x1 + x2 + x3", data = pd[c("x1", "x2", "x3")])
fit_pc2 <- cfa("dem60 =~ y1 + y2 + y3 + y4", data = pd[c("y1", "y2", "y3", "y4")])
fs_local2  <- get_fs(pd[cols2], mod2, local = TRUE)
fs_pc1     <- get_fs(fit_pc1)
fs_pc2     <- get_fs(fit_pc2)
fs_local2c <- get_fs(pd[cols2], mod2, local = TRUE, corrected_fsT = TRUE)
fs_pc1c    <- get_fs(fit_pc1, corrected_fsT = TRUE)
fs_pc2c    <- get_fs(fit_pc2, corrected_fsT = TRUE)
fs_local2s <- get_fs(pd[cols2], mod2, local = TRUE, std.lv = TRUE)
fs_pc1s    <- get_fs(cfa("ind60 =~ x1 + x2 + x3", data = pd, std.lv = TRUE))
fs_pc2s    <- get_fs(cfa("dem60 =~ y1 + y2 + y3 + y4", data = pd, std.lv = TRUE))
fs_local2m <- get_fs(pd[cols2], mod2, local = TRUE, method = "mean")

# --- Multigroup (HS, school) ---
hs <- HolzingerSwineford1939
mod2g <- "visual =~ x1 + x2 + x3
          speed  =~ x7 + x8 + x9"
fs_mg <- get_fs(hs, mod2g, group = "school", local = TRUE)
fs_v_mg <- get_fs(hs, "visual =~ x1 + x2 + x3", group = "school", format = "list")
fs_s_mg <- get_fs(hs, "speed =~ x7 + x8 + x9", group = "school", format = "list")

# --- FIML (HS with NA injection + forced all-NA visual rows; fixed seed) ---
mk_fiml_hs <- function() {
  d <- HolzingerSwineford1939
  set.seed(1334)
  d$x2[!rbinom(nrow(d), 1L, 0.4)] <- NA
  d$x8[!rbinom(nrow(d), 1L, 0.4)] <- NA
  d[c("x1", "x2", "x3")][c(1:3), ] <- NA   # rows with all of visual's items missing
  d
}
fiml_d <- mk_fiml_hs()
fs_fiml   <- suppressWarnings(get_fs(fiml_d, mod2g, local = TRUE, missing = "fiml"))
fs_v_fiml <- suppressWarnings(get_fs(
  suppressWarnings(cfa("visual =~ x1 + x2 + x3", data = fiml_d, missing = "fiml"))))
fs_s_fiml <- suppressWarnings(get_fs(
  suppressWarnings(cfa("speed =~ x7 + x8 + x9", data = fiml_d, missing = "fiml"))))

# Per-row block-diagonal assembly of the two per-local FIML results: the
# independent control for the merged per-row fsT/fsL/fsb (and for the
# cbind-of-locals downstream fit, which carries these on the data frame).
rv_fiml <- R2spa:::resolve_fs_per_row(fs_v_fiml)
rs_fiml <- R2spa:::resolve_fs_per_row(fs_s_fiml)
n_fiml <- nrow(fs_fiml)
Texp_fiml <- lapply(seq_len(n_fiml), function(r) {
  pv <- rv_fiml$pattern_idx[r]; ps <- rs_fiml$pattern_idx[r]
  b <- rbind(c(rv_fiml$blocks[[pv]]$fsT[1, 1], 0),
             c(0, rs_fiml$blocks[[ps]]$fsT[1, 1]))
  dimnames(b) <- list(c("fs_visual", "fs_speed"), c("fs_visual", "fs_speed")); b
})
Lexp_fiml <- lapply(seq_len(n_fiml), function(r) {
  pv <- rv_fiml$pattern_idx[r]; ps <- rs_fiml$pattern_idx[r]
  b <- rbind(c(rv_fiml$blocks[[pv]]$fsL[1, 1], 0),
             c(0, rs_fiml$blocks[[ps]]$fsL[1, 1]))
  dimnames(b) <- list(c("fs_visual", "fs_speed"), c("visual", "speed")); b
})
Bexp_fiml <- lapply(seq_len(n_fiml), function(r) {
  pv <- rv_fiml$pattern_idx[r]; ps <- rs_fiml$pattern_idx[r]
  c(fs_visual = rv_fiml$blocks[[pv]]$fsb[1], fs_speed = rs_fiml$blocks[[ps]]$fsb[1])
})
cb_fiml <- cbind(fs_v_fiml, fs_s_fiml)
attr(cb_fiml, "fsT") <- Texp_fiml
attr(cb_fiml, "fsL") <- Lexp_fiml
attr(cb_fiml, "fsb") <- Bexp_fiml
attr(cb_fiml, "per_obs") <- TRUE

# Derived vs explicit downstream fits (built once at file scope to bound cost)
fit_mg_d <- tspa("visual ~ speed", data = fs_mg, group = "school")
fit_fiml_d <- suppressWarnings(tspa("visual ~ speed", data = fs_fiml))
fit_fiml_c <- suppressWarnings(tspa("visual ~ speed", data = cb_fiml))

# ---------------------------------------------------------------------------
# 1. Parser
# ---------------------------------------------------------------------------

test_that("parser: valid splits preserve latent order (comments, ';', trailing-+) {", {
  v <- R2spa:::split_local_models(
    "# a comment
    ind60 =~ x1 + x2 +
      x3; dem60 =~ y1 + y2
    dem65 =~ y5 + y6 + y7 + y8")
  expect_named(v, c("ind60", "dem60", "dem65"))
  expect_equal(v[["ind60"]], "ind60 =~ x1 + x2 + x3")
  expect_equal(v[["dem60"]], "dem60 =~ y1 + y2")
  expect_equal(v[["dem65"]], "dem65 =~ y5 + y6 + y7 + y8")
  # 3+ latents, statement order preserved
  v3 <- R2spa:::split_local_models(model3)
  expect_named(v3, c("ind60", "dem60", "dem65"))
})

# The full set of unsupported-syntax classes the strict grammar rejects;
# every one must die with the stable "unsupported model syntax on line"
# prefix (the previously-raw-error classes -- '~~', structural '~',
# multi-latent LHS, non-bare LHS, non-identifier item, empty RHS -- now
# reach the same message, r-architect fix 2026-08).
parser_bad <- list(
  tilde_lv   = "a =~ x1 + x2\na ~~ b",
  tilde_res  = "a =~ x1 + x2\nx1 ~~ x2",
  struct     = "a =~ x1 + x2\nb ~ a",
  ordered    = "a |~ x1 + x2",
  threshold  = "a $ x1 + x2",
  no_meas    = "a b",
  two_meas   = "a =~ b =~ c",
  multi_lhs  = "a b =~ x1 + x2",
  bare_lhs   = "1a =~ x1 + x2",
  empty_rhs  = "a =~",
  rhs_tilde  = "a =~ x1 + y1 ~ x2",
  rhs_star   = "a =~ x1 * 1 + x2",
  call_rhs   = "a =~ c(x1, x2)",
  bad_item   = "a =~ x1 + 1x2",
  # dangling '+' (strsplit would drop the trailing empty field and the typo
  # would be fitted silently); the ';' form is the non-EOF case -- a dangling
  # '+' at the end of the input is already caught by the continuation check
  dangling_plus = "a =~ x1 + x2 +; b =~ y1"
)

test_that("parser: every unsupported-syntax class is rejected with the syntax prefix", {
  for (nm in names(parser_bad)) {
    expect_error(
      R2spa:::split_local_models(parser_bad[[nm]]),
      "unsupported model syntax on line", fixed = TRUE, label = nm
    )
  }
})

test_that("parser: the class-specific detail fragments name the offending syntax", {
  # fragments pinned exactly as local_model_syntax_error() emits them (the
  # multi-fragment details are reassembled verbatim by the `...` collapse)
  for (nm in c("tilde_lv", "tilde_res")) {
    expect_error(
      R2spa:::split_local_models(parser_bad[[nm]]),
      "is not supported in the local-mode string form", fixed = TRUE, label = nm
    )
  }
  expect_error(
    R2spa:::split_local_models(parser_bad[["struct"]]),
    "a structural path", fixed = TRUE
  )
  expect_error(
    R2spa:::split_local_models(parser_bad[["multi_lhs"]]),
    "multiple latent names", fixed = TRUE
  )
  expect_error(
    R2spa:::split_local_models(parser_bad[["bare_lhs"]]),
    "must be a single bare identifier", fixed = TRUE
  )
  expect_error(
    R2spa:::split_local_models(parser_bad[["bad_item"]]),
    "must be bare identifiers", fixed = TRUE
  )
  expect_error(
    R2spa:::split_local_models(parser_bad[["empty_rhs"]]),
    "an empty right-hand side", fixed = TRUE
  )
  expect_error(
    R2spa:::split_local_models(parser_bad[["dangling_plus"]]),
    "a dangling '+' operator", fixed = TRUE
  )
})

test_that("parser: duplicate item and duplicate latent are named", {
  expect_error(
    R2spa:::split_local_models("a =~ x1 + x2\nb =~ x2 + x3"),
    "are used in more than one statement", fixed = TRUE
  )
  expect_error(
    R2spa:::split_local_models("a =~ x1 + x2\na =~ x3 + x4"),
    "is defined more than once", fixed = TRUE
  )
  # trailing '+' with nothing following
  expect_error(
    R2spa:::split_local_models("a =~ x1 +"),
    "unsupported model syntax on line", fixed = TRUE
  )
})

test_that("parser: vector/list form -- verbatim fit, validation, order", {
  # within-factor residual covariance the strict grammar rejects: fits verbatim
  # (the y1 ~~ y4 fit is near-singular on this data, so muffle the theta warning)
  fs_vv <- suppressWarnings(get_fs(
    pd[cols2], model = c("ind60 =~ x1 + x2 + x3",
                         "dem60 =~ y1 + y2 + y3 + y4\ny1 ~~ y4"),
    local = TRUE))
  f2v <- suppressWarnings(
    cfa("dem60 =~ y1 + y2 + y3 + y4\ny1 ~~ y4",
        data = pd[c("y1", "y2", "y3", "y4")]))
  expect_equal(
    attr(fs_vv, "fsT")[[1L]][2, 2],
    attr(get_fs(f2v), "fsT")[[1L]][1, 1],
    tolerance = 1e-12
  )
  # an element defining two latents (or none) is rejected (the covariance-only
  # element is a degenerate fit, so muffle its vcov warning)
  expect_error(
    suppressWarnings(get_fs(pd[cols2], model = c(mod2, "x1 ~~ x2"), local = TRUE)),
    "every element must define exactly one latent", fixed = TRUE
  )
  # the same latent name in two elements (muffles a vcov warning from the
  # second, near-degenerate fit)
  expect_error(
    suppressWarnings(get_fs(
      pd[cols2], model = c("ind60 =~ x1 + x2 + x3", "ind60 =~ y1 + y2"),
      local = TRUE)),
    "latent names must be unique", fixed = TRUE
  )
  # a named element whose name does not match its fitted latent
  expect_error(
    get_fs(pd[cols2], model = list(wrongname = "ind60 =~ x1 + x2 + x3",
                                   dem60 = "dem60 =~ y1 + y2 + y3 + y4"),
           local = TRUE),
    "are named", fixed = TRUE
  )
  # latent order = element order (not name order) for a named list
  fs_nl <- get_fs(pd[cols2], model = list(
    dem60 = "dem60 =~ y1 + y2 + y3 + y4",
    ind60 = "ind60 =~ x1 + x2 + x3"), local = TRUE)
  expect_equal(head(colnames(fs_nl), 2L), c("fs_dem60", "fs_ind60"))
})

# ---------------------------------------------------------------------------
# 2-3. Single-group: equivalence to the canonical pattern + layout pin
# ---------------------------------------------------------------------------

test_that("SG: local score/SE columns identical to the cbind per-construct pattern", {
  expect_identical(fs_local3[score3], cb3[score3])
  expect_identical(fs_local3[se3], cb3[se3])
})

test_that("SG: merged fsT/fsL/psi are block-diagonals of the per-local attributes", {
  Tm <- attr(fs_local3, "fsT")[[1L]]
  Lm <- attr(fs_local3, "fsL")[[1L]]
  Pm <- attr(fs_local3, "psi")[[1L]]
  T1 <- attr(fs_p1, "fsT")[[1L]]; T2 <- attr(fs_p2, "fsT")[[1L]]; T3 <- attr(fs_p3, "fsT")[[1L]]
  L1 <- attr(fs_p1, "fsL")[[1L]]; L2 <- attr(fs_p2, "fsL")[[1L]]; L3 <- attr(fs_p3, "fsL")[[1L]]
  P1 <- attr(fs_p1, "psi")[[1L]]; P2 <- attr(fs_p2, "psi")[[1L]]; P3 <- attr(fs_p3, "psi")[[1L]]
  expect_equal(Tm, block_diag(T1, T2, T3), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(Lm, block_diag(L1, L2, L3), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(Pm, block_diag(P1, P2, P3), tolerance = 1e-12, ignore_attr = TRUE)
  # correct dimnames
  expect_equal(rownames(Tm), score3)
  expect_equal(rownames(Lm), score3)
  expect_equal(colnames(Lm), c("ind60", "dem60", "dem65"))
  expect_equal(rownames(Pm), c("ind60", "dem60", "dem65"))
})

test_that("SG: layout equals the joint layout; off-diagonal _by_/ecov_ are exactly 0", {
  # same column set and order as the joint get_fs() for the same latents
  expect_equal(colnames(fs_local3), colnames(fs_joint3))
  # off-diagonal _by_ columns: latent k loading on score i, i != k
  lv <- c("ind60", "dem60", "dem65")
  off_by <- unlist(lapply(seq_along(lv), function(k) {
    paste0(lv[k], "_by_fs_", setdiff(lv, lv[k]))
  }))
  expect_true(all(fs_local3[off_by] == 0))
  expect_true(all(fs_local3[grep("^ecov_", colnames(fs_local3), value = TRUE)] == 0))
  # single-group attribute convention: length-1 list named ""
  expect_length(attr(fs_local3, "fsT"), 1L)
  expect_named(attr(fs_local3, "fsT"), "")
  expect_length(attr(fs_local3, "scoring_matrix"), 1L)
})

test_that("SG: scoring_matrix reproduces the scores exactly (regression / Bartlett / mean)", {
  y <- as.matrix(pd[cols3])
  S  <- attr(fs_local3, "scoring_matrix")[[1L]]
  # regression: centered (no mean structure -> centered on the item means)
  rec <- sweep(y, 2L, colMeans(y), "-") %*% t(S)
  expect_equal(as.numeric(rec), as.numeric(as.matrix(fs_local3[score3])), tolerance = 1e-12)
  # Bartlett: the same centered identity
  fs_l3b <- get_fs(pd[cols3], model3, local = TRUE, method = "Bartlett")
  Sb <- attr(fs_l3b, "scoring_matrix")[[1L]]
  recb <- sweep(y, 2L, colMeans(y), "-") %*% t(Sb)
  expect_equal(as.numeric(recb), as.numeric(as.matrix(fs_l3b[score3])), tolerance = 1e-12)
  # mean: raw (uncentered) item means
  Sm <- attr(fs_local2m, "scoring_matrix")[[1L]]
  y2 <- as.matrix(pd[cols2])
  recm <- y2 %*% t(Sm)
  expect_equal(as.numeric(recm),
               as.numeric(as.matrix(fs_local2m[c("fs_ind60", "fs_dem60")])),
               tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 4. Joint A/B
# ---------------------------------------------------------------------------

test_that("A/B: uncorrelated local scores/attrs == joint (1e-4); free local != joint", {
  # uncorrelated bridge: the likelihood factorizes, so local == joint to
  # optimizer tolerance (observed ~1e-5; use 1e-4 to be optimizer-robust)
  expect_equal(as.numeric(as.matrix(fs_local3[score3])),
               as.numeric(as.matrix(fs_unc3[score3])), tolerance = 1e-4)
  expect_equal(attr(fs_local3, "fsT")[[1L]], attr(fs_unc3, "fsT")[[1L]],
               tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(attr(fs_local3, "fsL")[[1L]], attr(fs_unc3, "fsL")[[1L]],
               tolerance = 1e-4, ignore_attr = TRUE)
  # free correlation: local scores are pure per-construct and differ from the
  # joint (which leaks every factor's items into every score)
  d_free <- max(abs(as.numeric(as.matrix(fs_joint3[score3])) -
                    as.numeric(as.matrix(fs_local3[score3]))))
  expect_gt(d_free, 0.1)
})

# ---------------------------------------------------------------------------
# 5. Multigroup
# ---------------------------------------------------------------------------

test_that("MG: group column, per-group block-diagonal attrs, row counts", {
  expect_equal(attr(fs_mg, "group_col"), "school")
  expect_equal(sort(unique(fs_mg$school)), c("Grant-White", "Pasteur"))
  # per-group row counts match the data
  expect_equal(as.integer(sum(fs_mg$school == "Pasteur")),
               as.integer(sum(hs$school == "Pasteur")))
  expect_equal(as.integer(sum(fs_mg$school == "Grant-White")),
               as.integer(sum(hs$school == "Grant-White")))
  # per-group block-diagonal attributes equal the per-local per-group values
  for (g in c("Pasteur", "Grant-White")) {
    Tm <- attr(fs_mg, "fsT")[[g]]
    Lm <- attr(fs_mg, "fsL")[[g]]
    expect_equal(Tm[1, 1], as.numeric(attr(fs_v_mg[[g]], "fsT")), label = paste("T v", g))
    expect_equal(Tm[2, 2], as.numeric(attr(fs_s_mg[[g]], "fsT")), label = paste("T s", g))
    expect_equal(Lm[1, 1], attr(fs_v_mg[[g]], "fsL")[1, 1], label = paste("L v", g))
    expect_equal(Lm[2, 2], attr(fs_s_mg[[g]], "fsL")[1, 1], label = paste("L s", g))
    expect_equal(Tm[1, 2], 0, label = paste("offdiag", g))
  }
})

test_that("MG: derived tspa == explicit full-triple control; group.equal is forwarded", {
  g <- c("Pasteur", "Grant-White")
  T_exp <- setNames(lapply(g, function(k) {
    b <- rbind(c(as.numeric(attr(fs_v_mg[[k]], "fsT")), 0),
               c(0, as.numeric(attr(fs_s_mg[[k]], "fsT"))))
    dimnames(b) <- list(c("fs_visual", "fs_speed"), c("fs_visual", "fs_speed")); b
  }), g)
  L_exp <- setNames(lapply(g, function(k) {
    b <- rbind(c(attr(fs_v_mg[[k]], "fsL")[1, 1], 0),
               c(0, attr(fs_s_mg[[k]], "fsL")[1, 1]))
    dimnames(b) <- list(c("fs_visual", "fs_speed"), c("visual", "speed")); b
  }), g)
  B_exp <- setNames(lapply(g, function(k) {
    c(fs_visual = as.numeric(attr(fs_v_mg[[k]], "fsb")),
      fs_speed = as.numeric(attr(fs_s_mg[[k]], "fsb")))
  }), g)
  fit_mg_e <- tspa("visual ~ speed", data = fs_mg, group = "school",
                   fsT = T_exp, fsL = L_exp, fsb = B_exp)
  expect_identical(attr(fit_mg_d, "tspaModel"), attr(fit_mg_e, "tspaModel"))
  expect_equal(coef(fit_mg_d), coef(fit_mg_e), tolerance = 1e-10)
  # group.equal is forwarded to every local cfa()
  fs_ge <- get_fs(hs, mod2g, group = "school", local = TRUE,
                  group.equal = c("loadings", "intercepts"))
  expect_equal(nrow(fs_ge), nrow(hs))
})

# ---------------------------------------------------------------------------
# 6. FIML (per-row merge)
# ---------------------------------------------------------------------------

test_that("FIML: per-row attribute lists, per_obs marker, per-row blocks == per-local", {
  T_f <- attr(fs_fiml, "fsT"); L_f <- attr(fs_fiml, "fsL"); B_f <- attr(fs_fiml, "fsb")
  S_f <- attr(fs_fiml, "scoring_matrix")
  expect_true(is.list(T_f) && !is.matrix(T_f))
  expect_length(T_f, nrow(fs_fiml))
  expect_length(L_f, nrow(fs_fiml))
  expect_length(B_f, nrow(fs_fiml))
  expect_length(S_f, nrow(fs_fiml))
  expect_true(isTRUE(attr(fs_fiml, "per_obs")))
  # every row's block is the block-diagonal of the two per-local blocks
  ok <- vapply(seq_len(n_fiml), function(r) {
    pv <- rv_fiml$pattern_idx[r]; ps <- rs_fiml$pattern_idx[r]
    Tv <- rv_fiml$blocks[[pv]]$fsT; Ts <- rs_fiml$blocks[[ps]]$fsT
    Tr <- T_f[[r]]
    isTRUE(all.equal(unname(Tr[1, 1]), unname(Tv[1, 1]), tolerance = 1e-14)) &&
      isTRUE(all.equal(unname(Tr[2, 2]), unname(Ts[1, 1]), tolerance = 1e-14)) &&
      isTRUE(all.equal(unname(Tr[1, 2]), 0))
  }, logical(1L))
  expect_true(all(ok))
  # row-aligned through fs_indiv() on the merged vs each per-local output
  ind_m <- fs_indiv(fs_fiml)
  ind_v <- fs_indiv(fs_v_fiml)
  ind_s <- fs_indiv(fs_s_fiml)
  expect_equal(max(abs(ind_m$fs_visual_se - ind_v$fs_visual_se), na.rm = TRUE), 0)
  expect_equal(max(abs(ind_m$ev_fs_speed - ind_s$ev_fs_speed), na.rm = TRUE), 0)
})

test_that("FIML: all-NA rows get the all-NA block for that latent; listwise default errors", {
  # rows 1:3 have all of visual's items missing -> all-NA 1x1 visual block
  r <- 1L
  expect_true(all(is.na(attr(fs_fiml, "fsT")[[r]][1, 1])))
  expect_true(is.na(fs_fiml$fs_visual[r]))
  expect_true(is.na(fs_fiml$fs_visual_se[r]))
  # ... while the speed block for the same row is still scored
  expect_true(is.finite(attr(fs_fiml, "fsT")[[r]][2, 2]))
  expect_true(is.finite(fs_fiml$fs_speed[r]))
  # listwise deletion (the cfa() default) drops different rows per local fit
  expect_error(
    suppressWarnings(get_fs(fiml_d, mod2g, local = TRUE)),
    "the per-latent fits have different row counts", fixed = TRUE
  )
})

test_that("FIML: derived tspa (pooled) == cbind-of-locals control; fs_indiv dispatches", {
  expect_identical(attr(fit_fiml_d, "pooled_fs"), "mean")
  expect_identical(attr(fit_fiml_d, "tspaModel"), attr(fit_fiml_c, "tspaModel"))
  expect_equal(coef(fit_fiml_d), coef(fit_fiml_c), tolerance = 1e-10)
  # fs_indiv() works on the per-obs merged result (marker dispatch)
  ind <- fs_indiv(fs_fiml)
  expect_equal(nrow(ind), nrow(fs_fiml))
  expect_true(all(c("fs_visual", "fs_speed", "fs_visual_se", "ev_fs_speed") %in%
                    colnames(ind)))
})

# ---------------------------------------------------------------------------
# 7. Downstream (complete data)
# ---------------------------------------------------------------------------

test_that("downstream: derived mf fit carries zero cross-loadings; se_fs-only; mx", {
  fit_l2 <- tspa("dem60 ~ ind60", data = fs_local2)
  m <- attr(fit_l2, "tspaModel")
  # off-diagonal score loadings and the error covariance are pinned at 0
  expect_true(grepl("0 * fs_dem60", m, fixed = TRUE))
  expect_true(grepl("0 * fs_ind60", m, fixed = TRUE))
  # the se_fs-only path (explicit se_fs suppresses the mf derivation) still fits
  fit_sf <- suppressWarnings(tspa("dem60 ~ ind60", data = fs_local2,
                                  se_fs = c(ind60 = fs_local2$fs_ind60_se[1L],
                                            dem60 = fs_local2$fs_dem60_se[1L])))
  lines_sf <- strsplit(attr(fit_sf, "tspaModel"), "\n")[[1L]]
  expect_false(any(grepl("=~", lines_sf) & grepl("\\+", lines_sf)))
  skip_if_not_installed("OpenMx")
  # tspa_mx_model() accepts the merged (block-diagonal) fsL
  Lm2 <- attr(fs_local2, "fsL")[[1L]]
  Tm2 <- attr(fs_local2, "fsT")[[1L]]
  dat <- fs_indiv(fs_local2, include_intercept = TRUE)
  fmx <- suppressWarnings(tspa_mx_model(
    "dem60 ~ ind60", data = dat, fsL = Lm2, fsT = Tm2,
    fsb = c(fs_ind60 = "int_fs_ind60", fs_dem60 = "int_fs_dem60")))
  expect_true(inherits(fmx, "MxModel"))
})

# ---------------------------------------------------------------------------
# 8. Guards
# ---------------------------------------------------------------------------

test_that("guards: vfsLT / prior_cov / reliability / fitted-object are rejected", {
  expect_error(
    get_fs(pd[cols2], mod2, local = TRUE, vfsLT = TRUE),
    "'vfsLT = TRUE' is not supported with 'local = TRUE'", fixed = TRUE
  )
  expect_error(
    get_fs(pd[cols2], mod2, local = TRUE, prior_cov = 0.33),
    "'prior_cov' is not supported with 'local = TRUE' (v1)", fixed = TRUE
  )
  expect_error(
    get_fs(pd[cols2], mod2, local = TRUE, reliability = TRUE),
    "'reliability = TRUE' is not supported with 'local = TRUE' (v1)", fixed = TRUE
  )
  fit_j <- cfa(mod2, data = pd[cols2])
  expect_error(
    get_fs(fit_j, local = TRUE),
    "'local = TRUE' is only supported when 'object' is a data frame", fixed = TRUE
  )
})

test_that("guards: model = NULL is a no-op (the auto single-factor model)", {
  a <- get_fs(pd[c("x1", "x2", "x3")], local = TRUE)
  b <- get_fs(pd[c("x1", "x2", "x3")])
  expect_identical(a, b)
})

test_that("guard: model = NULL with no indicator columns errors clearly", {
  g <- data.frame(g = rep(c("a", "b"), each = 3L))
  expect_error(
    get_fs(g, group = "g"),
    "no indicator columns"
  )
})

# ---------------------------------------------------------------------------
# 9-10. corrected_fsT and std.lv
# ---------------------------------------------------------------------------

test_that("corrected_fsT = TRUE: merged blocks equal the per-local corrected values", {
  Tm <- attr(fs_local2c, "fsT")[[1L]]
  expect_equal(Tm[1, 1], as.numeric(attr(fs_pc1c, "fsT")[[1L]]), tolerance = 1e-10)
  expect_equal(Tm[2, 2], as.numeric(attr(fs_pc2c, "fsT")[[1L]]), tolerance = 1e-10)
  expect_equal(Tm[1, 2], 0)
})

test_that("std.lv = TRUE: psi diagonal is 1 per latent; scores equal the per-local std.lv", {
  Pm <- attr(fs_local2s, "psi")[[1L]]
  expect_equal(unname(diag(Pm)), c(1, 1))
  expect_equal(Pm[1, 2], 0)
  expect_equal(as.numeric(fs_local2s$fs_ind60), as.numeric(fs_pc1s$fs_ind60),
               tolerance = 1e-12)
  expect_equal(as.numeric(fs_local2s$fs_dem60), as.numeric(fs_pc2s$fs_dem60),
               tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Edge: sum_items and prior_mean slicing
# ---------------------------------------------------------------------------

test_that("edge: user-supplied sum_items and prior_mean are sliced per latent", {
  # sum_items (method = "mean"): each factor's items define its sum score
  fs_si <- get_fs(pd[cols2], mod2, local = TRUE, method = "mean",
                  sum_items = list(ind60 = c("x1", "x2", "x3"),
                                   dem60 = c("y1", "y2", "y3", "y4")))
  expect_equal(as.numeric(fs_si$fs_ind60),
               as.numeric(rowMeans(pd[cols2][c("x1", "x2", "x3")])), tolerance = 1e-12)
  # an unknown factor in sum_items is rejected
  expect_error(
    get_fs(pd[cols2], mod2, local = TRUE, method = "mean",
           sum_items = list(ind60 = c("x1", "x2", "x3"), bogus = c("y1"))),
    "Unknown factor name", fixed = TRUE
  )
  # prior_mean (length-q named vector) is validated once, then sliced
  fs_pm <- get_fs(pd[cols2], mod2, local = TRUE,
                  prior_mean = c(ind60 = 0.5, dem60 = -0.5))
  fs_pm1 <- get_fs(cfa("ind60 =~ x1 + x2 + x3", data = pd), prior_mean = 0.5)
  expect_equal(as.numeric(fs_pm$fs_ind60), as.numeric(fs_pm1$fs_ind60), tolerance = 1e-12)
  # wrong names are rejected against the (combined) latent names
  expect_error(
    get_fs(pd[cols2], mod2, local = TRUE, prior_mean = c(foo = 0.5, bar = -0.5)),
    "'prior_mean' names must match the latent variable names", fixed = TRUE
  )
})
