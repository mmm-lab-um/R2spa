# Quarantined with R/tspa_corrected_se.R (see archive/PLAN_QUARANTINE.md).
# Provenance:
#  - tests/testthat/test-tspa_render.R: vcov_corrected() multigroup test
#    (lines 322-344 as of 2026-08-17) plus mod2g (line 17).
#  - tests/testthat/test-get_fs_priors.R: vcov_corrected() prior test
#    (lines 175-205) plus single-group setup (lines 1-18).
# T1/T2: globalenv scaffolding removed 2026-08-22 — tspa() records
# tspa_args (evaluated values), so refits are environment-agnostic; the
# passing of T1/T2 with fixtures in file scope only IS that regression
# proof.
# T3 goldens: the shipped central-difference output at the shipped step
# (h0 = 1e-5), R 4.6.1 / lavaan 0.7-2, re-derived 2026-08-22. Deterministic
# across runs at a fixed h (the stage-2 refits are reproducible); the
# correction is NOT step-invariant below h ~ 1e-6 (the optimizer-noise
# regime; at h = 1e-7 the d60~~d60 entry drops to ~16.89), so the goldens
# are pinned to h0, not to an h -> 0 limit.
# Re-derive: cfa(3-factor joint, PoliticalDemocracy) ->
# get_fs_lavaan(vfsLT = TRUE) -> tspa("dem60 ~ ind60; dem65 ~ ind60 +
# dem60") -> vcov_corrected(fit, vfsLT) - vcov(fit), read the named
# elements. Drift protocol: if a lavaan upgrade moves a golden, diff the
# base fit first — base unchanged but correction moved => bug in the fix,
# do NOT update the golden.
# T4 fixture: boo_joint.RDS = corrected-se vignette bootstrap (R = 1999),
# labels pinned by the vignette's setNames.

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

## Multigroup, two-factor (feeds T2, T7)
mod2g <- "visual =~ x1 + x2 + x3\nspeed =~ x7 + x8 + x9"
fs_mg <- get_fs(HolzingerSwineford1939, model = mod2g, std.lv = TRUE,
                group = "school", vfsLT = TRUE, format = "list")
tspa_mg <- tspa("visual ~ speed", data = do.call(rbind, fs_mg),
                fsT = attr(fs_mg, "fsT"), fsL = attr(fs_mg, "fsL"),
                group = "school")

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
  vc <- vcov_corrected(tspa_mg, vfsLT = attr(fs_mg, "vfsLT"))
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
  # across runs on the pinned R/lavaan.
  expect_equal(cor["dem60~~dem60", "dem60~~dem60"], 17.3455467,
               tolerance = 1e-3)
  expect_equal(cor["dem65~dem60", "dem65~dem60"], 1.1961761,
               tolerance = 1e-3)
  expect_equal(cor["dem60~ind60", "dem60~ind60"], 1.4245011,
               tolerance = 1e-3)
})

test_that("T4: corrected SEs are within a loose tolerance of the bootstrap MAD", {
  p <- test_path("..", "vignettes", "boo_joint.RDS")
  skip_if(!file.exists(p), "bootstrap fixture not shipped")
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
