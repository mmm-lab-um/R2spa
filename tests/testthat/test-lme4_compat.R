library(lme4)

## lme4-2.x theta-packing canary (drift detector for get_D()) ---------------
##
## get_D() (R/get_fs_methods.R) parses a merMod fit's @theta directly: one
## block per RE term, in cnms() (formula) order, p * (p + 1) / 2 entries
## per term, packed into the lower-triangular (column-major) Cholesky
## factor L of the SCALED covariance D / sigma^2. If a future lme4 release
## changes the packing, the term order, or the scale, exactly these tests
## fail with a matrix diff; the fix is confined to get_D(). Without this
## canary, such a change would silently corrupt the EB fsT/fsL/
## scoring_matrix attributes (D enters only there) while the scores stay
## green: scores come from getME("b") and never touch D. The pinned
## identities were measured live against lme4 2.0.x (2.0.6).

fit_ss <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

## 2+2 multi-term fixture, mirroring test-get_fscore.R exactly (set.seed
## 43): its mixed theta (3 + 3 parameters) is the historical mis-parse
## case -- a single-block fallback parse lands on one mixed 3 x 3 block
## instead of two 2 x 2 blocks.
set.seed(43)
d22 <- data.frame(
  g1 = gl(40, 10),
  g2 = factor(rep_len(1:6, 400)),
  x1 = rnorm(400),
  x2 = rnorm(400)
)
d22$y <- d22$x1 + d22$x2 +
  rep(rnorm(40), each = 10) +
  rep(rnorm(40, sd = 0.5), each = 10) * d22$x1 +
  rnorm(6)[as.integer(d22$g2)] +
  rnorm(6, sd = 0.5)[as.integer(d22$g2)] * d22$x2 +
  rnorm(400)
fit22 <- suppressMessages(
  lmer(y ~ x1 + x2 + (x1 | g1) + (x2 | g2), d22)
)

test_that(
  "lme4 2.x theta packing: sigma()^2 * tcrossprod(L) reproduces VarCorr() per term",
  {
    # Raw lme4-side identity, independent of any R2spa code: the block for
    # term i is the packed lower-triangular (column-major) Cholesky factor
    # L of the scaled covariance, so VarCorr(m)[[i]] == sigma(m)^2 *
    # tcrossprod(L). ignore_attr: VarCorr() matrices carry lme4's
    # `stddev`/`correlation` metadata, which is incidental here.
    for (m in list(fit_ss, fit22)) {
      s2 <- stats::sigma(m)^2
      p <- lengths(m@cnms, use.names = FALSE)
      expect_equal(sum(p * (p + 1L) / 2L), length(m@theta))
      off <- 0L
      for (i in seq_along(p)) {
        blk <- m@theta[(off + 1L):(off + p[i] * (p[i] + 1L) / 2L)]
        off <- off + p[i] * (p[i] + 1L) / 2L
        L_i <- diag(nrow = p[i])
        L_i[lower.tri(L_i, diag = TRUE)] <- blk
        expect_equal(
          unname(s2 * tcrossprod(L_i)),
          as.matrix(VarCorr(m)[[i]]),
          tolerance = 1e-15,
          ignore_attr = TRUE
        )
      }
    }
  }
)

test_that("lme4 2.x theta parsing: get_D() equals VarCorr()[[1]] / sigma()^2", {
  # Pins R2spa's own parse to the same identity, on the first term -- the
  # only term get_D() (and hence the EB attributes) consumes. ignore_attr
  # as above: VarCorr() metadata is incidental.
  for (m in list(fit_ss, fit22)) {
    expect_equal(
      unname(R2spa:::get_D(m)),
      as.matrix(VarCorr(m)[[1L]]) / stats::sigma(m)^2,
      tolerance = 1e-15,
      ignore_attr = TRUE
    )
  }
})

test_that("lme4 2.x multi-term theta parse is warning-free", {
  # The 2.x parse guard: a fallback single-block parser of the mixed
  # multi-term theta would warn (replacement length) or mis-parse; the
  # explicit @cnms split in get_D() (R/get_fs_methods.R) must not.
  for (m in list(fit_ss, fit22)) {
    nw <- 0L
    withCallingHandlers(
      R2spa:::get_D(m),
      warning = function(wn) {
        nw <<- nw + 1L
        invokeRestart("muffleWarning")
      }
    )
    expect_identical(nw, 0L)
  }
})
