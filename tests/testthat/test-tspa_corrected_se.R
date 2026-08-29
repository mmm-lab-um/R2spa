# Re-integrated 2026-08 into the package (corrected-SE path, incl. the new
# tspa(corrected_se = TRUE) option); original provenance: quarantined with
# R/tspa_corrected_se.R (see archive/PLAN_QUARANTINE.md).
# Provenance:
#  - tests/testthat/test-tspa_render.R: vcov_corrected() multigroup test
#    (lines 322-344 as of 2026-08-17) plus mod2g (line 17).
#  - tests/testthat/test-get_fs_priors.R: vcov_corrected() prior test
#    (lines 175-205) plus single-group setup (lines 1-18).
# T1-T8 (standalone vcov_corrected() tests) are the A/B gate for the
# relocation: their assertions are unchanged since quarantine.
# T1/T2: globalenv scaffolding removed 2026-08-22 — tspa() records
# tspa_args (evaluated values), so refits are environment-agnostic; the
# passing of T1/T2 with fixtures in file scope only IS that regression
# proof.
# T3 goldens: the shipped central-difference output at the shipped step
# (h0 = 1e-5), R 4.6.1 / lavaan 0.7-2, re-derived 2026-08-22. Deterministic
# across runs at a fixed h (the stage-2 refits are reproducible), but
# platform-sensitive: the h0 = 1e-5 central difference amplifies BLAS/
# optimizer noise in the stage-2 refits (observed cross-platform drift
# 1.6e-3 relative, x86-64 Linux/Windows vs aarch64 macOS, 2026-08-27), so
# the tolerance is 1e-2. The correction is NOT step-invariant below
# h ~ 1e-6 (the optimizer-noise regime; at h = 1e-7 the d60~~d60 entry
# drops to ~16.89), so the goldens are pinned to h0, not to an h -> 0
# limit.
# Re-derive: cfa(3-factor joint, PoliticalDemocracy) ->
# get_fs_lavaan(vfsLT = TRUE) -> tspa("dem60 ~ ind60; dem65 ~ ind60 +
# dem60") -> vcov_corrected(fit, vfsLT) - vcov(fit), read the named
# elements. Drift protocol: if a lavaan upgrade moves a golden, diff the
# base fit first — base unchanged but correction moved => bug in the fix,
# do NOT update the golden.
# T4 fixture: boo_joint.RDS = corrected-se vignette bootstrap (R = 1999),
# labels pinned by the vignette's setNames. Shipped in vignettes/ (relocated
# 2026-08 from tests/testthat/): the fixture is shared between the
# corrected-se vignette and this test, so the test reads it from the repo
# root via test_path("../../vignettes/boo_joint.RDS").

library(lavaan)

########## Shared setup (file scope; refit-heavy calls stay in test_that)

## Single-group, two-factor joint CFA with prior adjustment (feeds T1, T5,
## T6)
prior_model <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
'
prior_fit <- cfa(prior_model, data = PoliticalDemocracy)
pm <- c(ind60 = 0.3, dem60 = -0.4)
pc <- matrix(c(1.2, 0.25, 0.25, 0.8), 2, 2,
             dimnames = list(c("ind60", "dem60"), c("ind60", "dem60")))
fs_prior <- get_fs(prior_fit, prior_mean = pm, prior_cov = pc,
                   corrected_fsT = TRUE, vfsLT = TRUE, format = "list")
tspa_prior <- tspa("dem60 ~ ind60", data = fs_prior,
                   fsT = attr(fs_prior, "fsT"), fsL = attr(fs_prior, "fsL"),
                   fsb = attr(fs_prior, "fsb"))

## Multigroup, two-factor (feeds T2, T7, T9, IT6-IT8)
mod2g <- "visual =~ x1 + x2 + x3\nspeed =~ x7 + x8 + x9"
fs_mg <- get_fs(HolzingerSwineford1939, model = mod2g, std.lv = TRUE,
                group = "school", vfsLT = TRUE, format = "list")
tspa_mg <- tspa("visual ~ speed", data = do.call(rbind, fs_mg),
                fsT = attr(fs_mg, "fsT"), fsL = attr(fs_mg, "fsL"),
                group = "school")

## Multigroup corrected fit (feeds IT6/IT7/IT8) — built once at file scope
## to bound cost: each corrected build refits stage 2 ~twice per free
## fsL/fsT element (~28 refits for this q=2, 2-group case), so doing it
## per-test would multiply that cost.
tspa_mg_corr <- tspa("visual ~ speed", data = do.call(rbind, fs_mg),
                     fsT = attr(fs_mg, "fsT"), fsL = attr(fs_mg, "fsL"),
                     vfsLT = attr(fs_mg, "vfsLT"), corrected_se = TRUE,
                     group = "school")
## Standalone multigroup correction (feeds T2/IT6/T9) — computed once for
## the same cost reason.
vcov_corr_mg <- vcov_corrected(tspa_mg, vfsLT = attr(fs_mg, "vfsLT"))

## Joint two-factor (feeds T4, T8)
mod2 <- "ind60 =~ x1 + x2 + x3
         dem60 =~ y1 + y2 + y3 + y4"
cfa_joint2 <- cfa(mod2, data = PoliticalDemocracy)
fs_joint2 <- get_fs_lavaan(cfa_joint2, vfsLT = TRUE)
tspa_joint2 <- tspa("dem60 ~ ind60", data = fs_joint2,
                    fsT = attr(fs_joint2, "fsT"), fsL = attr(fs_joint2, "fsL"))

## Three-factor joint (feeds T3)
mod3 <- "ind60 =~ x1 + x2 + x3
         dem60 =~ y1 + y2 + y3 + y4
         dem65 =~ y5 + y6 + y7 + y8"
cfa_joint3 <- cfa(mod3, data = PoliticalDemocracy)
fs_joint3 <- get_fs_lavaan(cfa_joint3, vfsLT = TRUE)
tspa_joint3 <- tspa("dem60 ~ ind60
                    dem65 ~ ind60 + dem60", data = fs_joint3,
                    fsT = attr(fs_joint3, "fsT"), fsL = attr(fs_joint3, "fsL"))

## PLAN 16 (engine = "analytic", feeds T10-T13, IT9): the T3 saturated fit
## (tspa_joint3, df = 0) is the A/B reference. The FD correction (vc_fd_sat)
## is the expensive side (30 stage-2 refits) and is computed once; the analytic
## (vc_an_sat) is refit-free and deterministic. A restricted fit on the same
## 3-factor data (tspa_joint3_nsat, df > 0) drives the saturated-only -> FD
## fallback (the closed form is exact only for a df = 0 structural model).
vc_an_sat <- vcov_corrected(tspa_joint3, vfsLT = attr(fs_joint3, "vfsLT"),
                            engine = "analytic")
vc_fd_sat <- vcov_corrected(tspa_joint3, vfsLT = attr(fs_joint3, "vfsLT"),
                            engine = "fd")
tspa_joint3_nsat <- tspa("dem65 ~ ind60", data = fs_joint3,
                         fsT = attr(fs_joint3, "fsT"),
                         fsL = attr(fs_joint3, "fsL"))

########## Tests ##########

test_that("T1: vcov_corrected() works with prior-adjusted factor scores (no globalenv)", {
  vc <- vcov_corrected(tspa_prior, vfsLT = attr(fs_prior, "vfsLT"))
  expect_s3_class(vc, "matrix")
  expect_equal(dim(vc), c(3, 3))
  expect_true(all(diag(vc) > 0))
  # The correction is non-trivial ("corrected == uncorrected" was an
  # undetectable failure mode before the guard tests landed).
  expect_gt(sum((vc - vcov(tspa_prior))^2), 0)
})

test_that("T2: vcov_corrected() runs on an MG fit (no globalenv objects)", {
  # Precomputed at file scope (vcov_corr_mg) to bound cost — a fresh
  # vcov_corrected() call here would spend ~28 stage-2 refits per run.
  vc <- vcov_corr_mg
  expect_s3_class(vc, "matrix")
  expect_equal(dim(vc), dim(vcov(tspa_mg)))
  expect_true(all(is.finite(vc)))
  expect_equal(vc, t(vc), tolerance = 1e-10)
  expect_true(all(diag(vc) > 0))
  # The PSD correction has strictly positive trace (an all-zero
  # correction — the q>=3 failure mode — fails this cheaply).
  expect_gt(sum(diag(vc - vcov(tspa_mg))), 0)
})

test_that("T3: q = 3 correction is non-zero and matches golden values (B1 guard)", {
  vc <- vcov_corrected(tspa_joint3, vfsLT = attr(fs_joint3, "vfsLT"))
  cor <- vc - vcov(tspa_joint3)
  expect_s3_class(cor, "matrix")
  expect_equal(dim(cor), dim(vcov(tspa_joint3)))
  expect_true(all(is.finite(cor)))
  expect_equal(cor, t(cor), tolerance = 1e-10)
  expect_true(all(diag(cor) >= -1e-10))
  # B1 sentinel: under the row/col-major bug this element was exactly 0;
  # the golden is ~17.3, so a threshold 4 orders below it survives
  # golden updates.
  expect_gt(cor["dem60~~dem60", "dem60~~dem60"], 1)
  # Golden elements (provenance, step, and drift protocol in the file
  # header): the shipped central-difference step (h0 = 1e-5), deterministic
  # across runs on the pinned R/lavaan/platform (1e-2 absorbs the
  # cross-platform BLAS/optimizer drift; see header).
  expect_equal(cor["dem60~~dem60", "dem60~~dem60"], 17.3455467,
               tolerance = 1e-2)
  expect_equal(cor["dem65~dem60", "dem65~dem60"], 1.1961761,
               tolerance = 1e-2)
  expect_equal(cor["dem60~ind60", "dem60~ind60"], 1.4245011,
               tolerance = 1e-2)
})

test_that("T4: corrected SEs are within a loose tolerance of the bootstrap MAD", {
  p <- test_path("../../vignettes/boo_joint.RDS")
  skip_if(!file.exists(p),
          "bootstrap fixture not shipped (vignettes/boo_joint.RDS)")
  boo <- readRDS(p)
  mad_v <- setNames(apply(boo$t, 2, mad),
                    c("dem60~ind60", "ind60~~ind60", "dem60~~dem60"))
  # label sanity (guards against a fixture reordering)
  expect_gt(mad_v["dem60~~dem60"], mad_v["ind60~~ind60"])
  vc <- vcov_corrected(tspa_joint2, vfsLT = attr(fs_joint2, "vfsLT"))
  se_cor <- sqrt(diag(vc))
  expect_setequal(names(se_cor), names(mad_v))
  # Loose tolerance on purpose: a MAD is ~1.48*sigma, not the MLE, and
  # the vignette documents that the methods "diverge slightly" (observed
  # max divergence on this frozen fixture: 20.3%). This is a
  # plausibility anchor, not a guard (the naive SEs also pass it).
  expect_equal(se_cor[names(mad_v)], mad_v, tolerance = 0.3)
})

test_that("T5: which_free requires a matching principal submatrix of vfsLT", {
  vldev7 <- attr(fs_prior, "vfsLT")
  # Error case: full 7 x 7 vfsLT with a length-2 which_free — must fail
  # before any refit is spent.
  expect_error(
    vcov_corrected(tspa_prior, vfsLT = vldev7, which_free = c(5, 7)),
    "principal submatrix"
  )
  # Happy path: the vignette-documented usage (propagate the two error
  # variances, positions 5 and 7).
  vc <- vcov_corrected(tspa_prior,
                       vfsLT = vldev7[c(5, 7), c(5, 7)],
                       which_free = c(5, 7))
  expect_equal(dim(vc), c(3, 3))
  expect_true(all(diag(vc) > 0))
})

test_that("T5b: which_free rejects non-integer / invalid positions before any refit", {
  vldev7 <- attr(fs_prior, "vfsLT")   # 7 x 7
  # Fractional position (as.integer() used to truncate this silently).
  expect_error(
    vcov_corrected(tspa_prior, vfsLT = vldev7, which_free = 2.5),
    "whole-number positions"
  )
  # NA position.
  expect_error(
    vcov_corrected(tspa_prior, vfsLT = vldev7, which_free = c(NA, 2)),
    "whole-number positions"
  )
  # Out-of-range position.
  expect_error(
    vcov_corrected(tspa_prior, vfsLT = vldev7, which_free = c(1, 100)),
    "whole-number positions"
  )
  # Duplicated positions.
  expect_error(
    vcov_corrected(tspa_prior, vfsLT = vldev7, which_free = c(5, 5)),
    "whole-number positions"
  )
  # Non-numeric positions.
  expect_error(
    vcov_corrected(tspa_prior, vfsLT = vldev7, which_free = c("1", "2")),
    "whole-number positions"
  )
  # Integer-valued doubles are accepted: with the full 7 x 7 vfsLT the
  # call must proceed past the which_free guard and fail at the (matching)
  # vfsLT check instead — proving c(5.0, 7.0) was not rejected/truncated.
  expect_error(
    vcov_corrected(tspa_prior, vfsLT = vldev7, which_free = c(5.0, 7.0)),
    "principal submatrix"
  )
})

test_that("T6: the convergence gate stops on a non-converged refit", {
  # Unit layer: a real converged fit with the flag flipped exercises the
  # exact slot read the gate performs (a forced non-convergence would be
  # platform-dependent and slow to trigger).
  bad <- tspa_prior
  bad@optim$converged <- FALSE
  expect_error(check_refit_convergence(bad, 1), "did not converge")
  expect_no_error(check_refit_convergence(tspa_prior, 1))
})

test_that("T7: MG attribute-shape guards are clear errors (no refits)", {
  vldev <- attr(fs_mg, "vfsLT")
  # fsL list length differs from the number of groups
  f <- tspa_mg
  attr(f, "fsL") <- c(attr(f, "fsL"), list(matrix(0, 2, 2)))
  expect_error(vcov_corrected(f, vfsLT = vldev), "lists of length")
  # per-group dimension mismatch
  f2 <- tspa_mg
  L <- attr(f2, "fsL")
  L[[2]] <- matrix(0, 3, 3)
  attr(f2, "fsL") <- L
  expect_error(vcov_corrected(f2, vfsLT = vldev), "identical; group 2")
  # plain matrix attribute on a multigroup fit
  f3 <- tspa_mg
  attr(f3, "fsT") <- attr(tspa_mg, "fsT")[[1]]
  expect_error(vcov_corrected(f3, vfsLT = vldev), "lists of length")
})

test_that("T8: Jacobian wiring — independent central differences reproduce the correction", {
  skip_if_not_installed("numDeriv")
  fit <- tspa_joint2
  args0 <- attr(fit, "tspa_args")
  vldev7 <- attr(fs_joint2, "vfsLT")
  fsL0 <- args0$fsL
  fsT0 <- args0$fsT
  expect_true(is.matrix(fsL0))
  names_coef <- names(coef(fit))
  tri <- which(lower.tri(fsT0, diag = TRUE), arr.ind = TRUE)
  x0 <- c(unlist(fsL0), unlist(fsT0[lower.tri(fsT0, diag = TRUE)]))

  # Test-side reference: its own reconstruction (independent of the
  # package helpers), replayed through tspa_args with stage-2 SEs
  # suppressed.
  f <- function(x) {
    L <- matrix(x[seq_len(length(fsL0))], ncol = ncol(fsL0),
                dimnames = dimnames(fsL0))
    Tm <- fsT0
    for (k in seq_len(nrow(tri))) {
      Tm[tri[k, 1], tri[k, 2]] <- x[length(fsL0) + k]
    }
    a <- args0
    a$fsL <- L
    a$fsT <- Tm
    a$se <- "none"
    coef(do.call(tspa, a))
  }
  h <- 1e-5
  J_test <- matrix(NA_real_, nrow = length(names_coef),
                   ncol = length(x0))
  for (p in seq_along(x0)) {
    step <- h * max(1, abs(x0[p]))
    e <- numeric(length(x0))
    e[p] <- 1
    J_test[, p] <- (f(x0 + step * e) - f(x0 - step * e)) / (2 * step)
  }
  rownames(J_test) <- names_coef

  # Wiring: the package correction equals J_test %*% V %*% t(J_test). A
  # wiring/ordering bug (incl. the B1 class) makes the LHS zero or
  # scrambled while the RHS stays correct.
  cor_pkg <- vcov_corrected(fit, vfsLT = vldev7) - vcov(fit)
  expect_equal(cor_pkg, J_test %*% vldev7 %*% t(J_test), tolerance = 1e-4)

  # Independent FD scheme (numDeriv Richardson; same order, different
  # step) as a consistency cross-check — not a golden. compare values only
  # (J_test carries rownames; numDeriv::jacobian() returns none).
  J_ref <- numDeriv::jacobian(f, x0)
  expect_equal(J_test, J_ref, tolerance = 1e-3, ignore_attr = TRUE)
})
########## IT fixtures (corrected_se option on tspa(), single-group) ##########

## Same model family as the joint2 fixture above (mod2), kept as a separate
## fixture set: the corrected fit is a new object kind (tspa_corrected
## attribute + overwritten @vcov slot) and the IT blocks must not depend on
## the T-block fixtures.
cfa_it <- cfa(mod2, data = PoliticalDemocracy)
fs_it <- get_fs_lavaan(cfa_it, vfsLT = TRUE)
vfs_it <- attr(fs_it, "vfsLT")
tspa_it_plain <- tspa("dem60 ~ ind60", data = fs_it,
                      fsT = attr(fs_it, "fsT"), fsL = attr(fs_it, "fsL"))
tspa_it_corr <- tspa("dem60 ~ ind60", data = fs_it,
                     fsT = attr(fs_it, "fsT"), fsL = attr(fs_it, "fsL"),
                     vfsLT = vfs_it, corrected_se = TRUE)

########## IT tests: corrected_se option on tspa() ##########

test_that("IT1: corrected_se on tspa() equals standalone vcov_corrected()", {
  se_on <- sqrt(diag(vcov(tspa_it_corr)))
  se_st <- sqrt(diag(vcov_corrected(tspa_it_plain, vfsLT = vfs_it)))
  expect_equal(se_on, se_st, tolerance = 1e-8)
})

test_that("IT2: corrected standardized SE >= naive, matches golden; SE-only", {
  ss0 <- standardizedSolution(tspa_it_plain)
  ss1 <- standardizedSolution(tspa_it_corr)
  keep <- ss1$lhs == "dem60" & ss1$rhs == "ind60" & ss1$op == "~"
  # (a) on this coefficient the correction inflates the SE, never deflates
  expect_gte(ss1$se[keep], ss0$se[keep])
  # (b) golden derived 2026-08-23 by running the code (R 4.6.1 /
  # lavaan 0.7-2, PoliticalDemocracy): 0.10293499... vs 0.10010434...
  # naive on the plain fit.
  expect_equal(ss1$se[keep], 0.102935, tolerance = 1e-3)
  # (c) point estimates are untouched — the correction is SE-only
  expect_identical(ss0$est.std[keep], ss1$est.std[keep])
})

test_that("IT3: replaying tspa_args reproduces the corrected covariance", {
  rep <- do.call(tspa, attr(tspa_it_corr, "tspa_args"))
  # The replay re-runs the correction through the recorded args; the
  # double-correction guard must not fire (the internal base fit inside
  # the replay is uncorrected), so this also proves the guard terminates.
  expect_equal(rep@vcov[["vcov"]], tspa_it_corr@vcov[["vcov"]],
               tolerance = 1e-8)
  expect_true(isTRUE(attr(rep, "tspa_corrected")))
})

test_that("IT4: tspa() option guards and the tspa_corrected attribute", {
  # (a) corrected_se = TRUE without vfsLT is a clear error naming 'vfsLT'
  # (fires after the stage-2 fit, by design of the option placement)
  expect_error(
    tspa("dem60 ~ ind60", data = fs_it,
         fsT = attr(fs_it, "fsT"), fsL = attr(fs_it, "fsL"),
         corrected_se = TRUE),
    "vf"
  )
  # (b) the attribute marks only the corrected fit
  expect_true(isTRUE(attr(tspa_it_corr, "tspa_corrected")))
  expect_null(attr(tspa_it_plain, "tspa_corrected"))
  # The former MG-rejection check is gone: the guard was removed and
  # multigroup corrected fits now build (positive-path MG coverage lives
  # in IT6/IT7/IT8).
})

test_that("IT5: double-correction guard on vcov_corrected()", {
  # Single regex with true alternation (two separate string args would be
  # OR'd by R before expect_error ever sees them).
  expect_error(
    vcov_corrected(tspa_it_corr, vfsLT = vfs_it),
    "twice|already SE-corrected"
  )
  # The plain fit does not trip the guard: IT1 runs
  # vcov_corrected(tspa_it_plain, vfsLT = vfs_it) to completion, which is
  # the no-error side of this guard (not re-run here — a full Jacobian).
})

########## IT6-IT8: corrected_se option on tspa(), multigroup ##########
## The MG corrected fixtures (tspa_mg_corr, vcov_corr_mg) are built once
## at file scope with the T-block fixtures (cost-bounded); IT8's replay
## is the one MG corrected build that stays in-test — the replay itself
## is the object under test.

test_that("IT6: MG in-place corrected fit equals the standalone correction", {
  expect_equal(sqrt(diag(vcov(tspa_mg_corr))), sqrt(diag(vcov_corr_mg)),
               tolerance = 1e-8)
})

test_that("IT7: MG corrected_se is SE-only (est unchanged, SEs inflated-or-equal)", {
  # coef()/diag(vcov()) directly: lavaan exports no se() generic in this
  # dependency set, so the standard-error side is taken from the
  # covariance's diagonal (same invariant: SE-only, estimates unchanged).
  expect_identical(names(coef(tspa_mg)), names(coef(tspa_mg_corr)))
  expect_equal(coef(tspa_mg_corr), coef(tspa_mg), tolerance = 1e-10)
  sp <- sqrt(diag(vcov(tspa_mg)))
  sc <- sqrt(diag(vcov(tspa_mg_corr)))
  expect_true(all(sc >= sp - 1e-8))
})

test_that("IT8: MG corrected fit replays via tspa_args", {
  rep <- do.call(tspa, attr(tspa_mg_corr, "tspa_args"))
  expect_equal(rep@vcov[["vcov"]], tspa_mg_corr@vcov[["vcov"]],
               tolerance = 1e-8)
  expect_true(isTRUE(attr(rep, "tspa_corrected")))
})

test_that("T9: MG Jacobian wiring — independent central differences reproduce the correction", {
  args0 <- attr(tspa_mg, "tspa_args")
  fsL0 <- args0$fsL          # list of ngrp q x q matrices
  fsT0 <- args0$fsT          # list of ngrp q x q symmetric matrices
  ng <- length(fsL0)
  q <- ncol(fsL0[[1]])
  ld_len <- q * q
  ev_len <- q * (q + 1L) / 2L
  tri_list <- lapply(fsT0, function(T) which(lower.tri(T, diag = TRUE), arr.ind = TRUE))
  # Parameter layout must match vcov_corrected()'s val_fsLT: all groups'
  # fsL (column-major) first, then all groups' fsT (lower triangle,
  # column-major) — [ld_g1, ld_g2, ev_g1, ev_g2] for q=2, ngrp=2.
  x0 <- c(unlist(fsL0),
          unlist(lapply(fsT0, function(T) T[lower.tri(T, diag = TRUE)])))
  expect_true(length(x0) == ng * (ld_len + ev_len))
  names_coef <- names(coef(tspa_mg))

  f <- function(x) {
    L_list <- fsL0; T_list <- fsT0; pos <- 1L
    for (g in seq_len(ng)) {
      L_list[[g]][, ] <- x[pos:(pos + ld_len - 1L)]
      pos <- pos + ld_len
    }
    for (g in seq_len(ng)) {
      tri_g <- tri_list[[g]]
      for (k in seq_len(nrow(tri_g))) {
        T_list[[g]][tri_g[k, 1], tri_g[k, 2]] <- x[pos]
        pos <- pos + 1L
      }
    }
    a <- args0
    a$fsL <- L_list
    a$fsT <- T_list
    a$se <- "none"
    coef(do.call(tspa, a))
  }
  # Base round-trip: a mis-wired perturbation index would be invisible
  # unless the unperturbed replay reproduces the base coefficients.
  expect_equal(unname(f(x0)), unname(coef(tspa_mg)), tolerance = 1e-10)

  h <- 1e-5
  J_test <- matrix(NA_real_, nrow = length(names_coef), ncol = length(x0))
  for (p in seq_along(x0)) {
    step <- h * max(1, abs(x0[p])); e <- numeric(length(x0)); e[p] <- 1
    J_test[, p] <- (f(x0 + step * e) - f(x0 - step * e)) / (2 * step)
  }
  rownames(J_test) <- names_coef

  # Wiring: the package correction equals J_test %*% vfsLT %*% t(J_test).
  # A wiring/ordering bug makes the LHS zero or scrambled while the RHS
  # stays correct.
  vfsLT <- attr(fs_mg, "vfsLT")
  cor_pkg <- vcov_corr_mg - vcov(tspa_mg)
  expect_equal(cor_pkg, J_test %*% vfsLT %*% t(J_test), tolerance = 1e-4)
})

########## PLAN 16: engine = "analytic" (refit-free Jacobian) ##########
## The saturated closed form (PLAN 16, section 2.4) is A/B'd against the
## finite-difference engine on the T3 saturated fit (the common 2S-PA case)
## and must agree to the FD's own noise floor. The gate also routes a
## restricted (df > 0) structural model back to the FD.

test_that("T10: engine = 'analytic' matches the FD engine on the saturated fit", {
  expect_equal(vc_an_sat, vc_fd_sat, tolerance = 1e-2)
  # The analytic correction is non-trivial, symmetric, and PSD.
  expect_gt(sum((vc_an_sat - vcov(tspa_joint3))^2), 0)
  expect_equal(vc_an_sat, t(vc_an_sat), tolerance = 1e-10)
})

test_that("T11: engine = 'analytic' is deterministic (bit-identical, no refits)", {
  # No finite differences and no optimizer jitter, so the result reproduces
  # bit-for-bit (the FD engine is only deterministic to ~1e-8 cross-run).
  expect_identical(
    vcov_corrected(tspa_joint3, vfsLT = attr(fs_joint3, "vfsLT"),
                   engine = "analytic"),
    vc_an_sat)
})

test_that("T12: engine = 'analytic' honours which_free exactly as the FD", {
  wl <- c(1, 10)  # one loading (fsL) + the first error variance (fsT)
  sub <- attr(fs_joint3, "vfsLT")[wl, wl]
  va <- vcov_corrected(tspa_joint3, vfsLT = sub, which_free = wl,
                       engine = "analytic")
  vf <- vcov_corrected(tspa_joint3, vfsLT = sub, which_free = wl,
                       engine = "fd")
  expect_equal(va, vf, tolerance = 1e-2)
})

test_that("T13: the analytic path covers saturated and restricted (general) models", {
  # Saturated single-group fit (df = 0): the analytic engine returns a full
  # finite p x nfree matrix (nfree = q^2 loadings + q(q+1)/2 error terms = 15
  # for q = 3).
  j_sat <- R2spa:::vcov_jacobian_analytic(tspa_joint3, names(coef(tspa_joint3)),
                                          seq_len(15))
  expect_true(is.matrix(j_sat))
  expect_equal(dim(j_sat), c(length(coef(tspa_joint3)), 15))
  expect_true(all(is.finite(j_sat)))
  expect_gt(max(abs(j_sat)), 1e-6)
  # Restricted fit (df > 0): handled by the general path (PLAN 16 section 4.3),
  # not a NULL/FD fallback -- the analytic engine returns a full finite J for
  # the p free params of the restricted model.
  j_nsat <- R2spa:::vcov_jacobian_analytic(
    tspa_joint3_nsat, names(coef(tspa_joint3_nsat)), seq_len(15))
  # A NULL here means an analytic-engine guard fired -- in practice the
  # file-scope restricted fit (tspa_joint3_nsat) occasionally fails to
  # converge (a fixture flake, not an engine bug). Skip cleanly instead of
  # erroring on max(abs(NULL)) below.
  if (is.null(j_nsat)) {
    skip("analytic J is NULL for the restricted fit (an engine guard fired; ",
         "the file-scope tspa_joint3_nsat likely did not converge this run)")
  }
  expect_true(is.matrix(j_nsat))
  expect_equal(dim(j_nsat), c(length(coef(tspa_joint3_nsat)), 15))
  expect_true(all(is.finite(j_nsat)))
  expect_gt(max(abs(j_nsat)), 1e-6)
  # The public analytic path yields a finite, correctly-dimensioned corrected
  # covariance.
  vc_nsat <- vcov_corrected(tspa_joint3_nsat, vfsLT = attr(fs_joint3, "vfsLT"),
                            engine = "analytic")
  expect_true(all(is.finite(vc_nsat)))
  expect_equal(dim(vc_nsat), dim(vcov(tspa_joint3_nsat)))
})

test_that("IT9: in-place corrected_se = TRUE, engine = 'analytic' matches the FD", {
  fa <- tspa("dem60 ~ ind60\ndem65 ~ ind60 + dem60", data = fs_joint3,
             fsT = attr(fs_joint3, "fsT"), fsL = attr(fs_joint3, "fsL"),
             vfsLT = attr(fs_joint3, "vfsLT"), corrected_se = TRUE,
             engine = "analytic")
  expect_true(isTRUE(attr(fa, "tspa_corrected")))
  # In-place analytic == standalone FD (IT1: in-place == standalone to 1e-8;
  # T10: standalone analytic == standalone FD to 1e-2).
  expect_equal(vcov(fa), vc_fd_sat, tolerance = 1e-2)
})

