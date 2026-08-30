# Tests for get_fs() with method = "mean" (sum-score factor scores)
library(lavaan)

########## Shared fits ############

mean_mod_1f <- "ind60 =~ x1 + x2 + x3"
mean_mod_2f <- "ind60 =~ x1 + x2 + x3
                dem60 =~ y1 + y2 + y3 + y4"
mean_mod_hs <- "visual =~ x1 + x2 + x3"

mean_fit_1f <- cfa(mean_mod_1f, data = PoliticalDemocracy)
mean_fit_2f <- cfa(mean_mod_2f, data = PoliticalDemocracy)
mean_fit_hs <- cfa(mean_mod_hs, data = HolzingerSwineford1939,
                   group = "school")

########## Testing section ############

test_that("method = 'mean' single-factor column layout and unified attributes", {
  fs <- get_fs(mean_fit_1f, method = "mean")
  expect_s3_class(fs, "data.frame")
  expect_equal(
    colnames(fs),
    c("fs_ind60", "fs_ind60_se", "ind60_by_fs_ind60", "ev_fs_ind60")
  )
  expect_identical(nrow(fs), nrow(PoliticalDemocracy))
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix")) {
    expect_false(is.null(attr(fs, ak)), info = ak)
  }
  # Unified single-group: length-1 list attributes named ""
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix")) {
    a <- attr(fs, ak)
    expect_true(is.list(a) && length(a) == 1, info = ak)
    expect_identical(names(a), "")
  }
  expect_true(is.matrix(attr(fs, "fsT")[[1]]))
  expect_true(is.matrix(attr(fs, "fsL")[[1]]))
  expect_true(is.numeric(attr(fs, "fsb")[[1]]))
  expect_true(is.matrix(attr(fs, "scoring_matrix")[[1]]))
  # Raw item means
  y1 <- as.matrix(lavInspect(mean_fit_1f, what = "data"))
  expect_equal(unname(as.numeric(fs$fs_ind60)), unname(rowMeans(y1)),
               tolerance = 1e-12)
  # Constant _se == sqrt(diag(fsT))
  expect_identical(var(fs$fs_ind60_se), 0)
  expect_equal(unname(as.numeric(fs$fs_ind60_se)),
               rep(sqrt(attr(fs, "fsT")[[1]][1, 1]), nrow(fs)),
               tolerance = 1e-14)
  # The scoring_matrix rows are the item weights of the sum (q x p)
  sm1 <- attr(fs, "scoring_matrix")[[1]]
  expect_equal(dim(sm1), c(1L, 3L))
  expect_equal(rownames(sm1), "ind60")
  expect_setequal(colnames(sm1), c("x1", "x2", "x3"))
})

test_that("method = 'mean' two-factor column layout, with exact zero cross terms", {
  fs <- get_fs(mean_fit_2f, method = "mean")
  expect_equal(
    colnames(fs),
    c("fs_ind60", "fs_dem60", "fs_ind60_se", "fs_dem60_se",
      "ind60_by_fs_ind60", "ind60_by_fs_dem60",
      "dem60_by_fs_ind60", "dem60_by_fs_dem60",
      "ev_fs_ind60", "ecov_fs_dem60_fs_ind60", "ev_fs_dem60")
  )
  fsT <- attr(fs, "fsT")[[1]]
  fsL <- attr(fs, "fsL")[[1]]
  expect_equal(dim(fsT), c(2L, 2L))
  expect_equal(dim(fsL), c(2L, 2L))
  # Plain CFA: disjoint item sets + diagonal theta -> exact zeros
  expect_equal(fsT[1, 2], 0, tolerance = 0)
  expect_equal(fsT[2, 1], 0, tolerance = 0)
  expect_equal(fsL[1, 2], 0, tolerance = 0)
  expect_equal(fsL[2, 1], 0, tolerance = 0)
  # Constant _se per factor == sqrt of the corresponding fsT diagonal
  expect_identical(var(fs$fs_ind60_se), 0)
  expect_identical(var(fs$fs_dem60_se), 0)
  expect_equal(unname(as.numeric(fs$fs_ind60_se)),
               rep(sqrt(fsT[1, 1]), nrow(fs)), tolerance = 1e-14)
  expect_equal(unname(as.numeric(fs$fs_dem60_se)),
               rep(sqrt(fsT[2, 2]), nrow(fs)), tolerance = 1e-14)
})

test_that("method = 'mean' score column means equal the fsb attribute", {
  fs1 <- get_fs(mean_fit_1f, method = "mean")
  expect_equal(mean(unname(as.numeric(fs1$fs_ind60))),
               unname(attr(fs1, "fsb")[[1]]))
  fs2 <- get_fs(mean_fit_2f, method = "mean")
  expect_equal(
    c(mean(unname(as.numeric(fs2$fs_ind60))),
      mean(unname(as.numeric(fs2$fs_dem60)))),
    unname(attr(fs2, "fsb")[[1]])
  )
})

test_that("method = 'mean' scoring_matrix %*% t(y) reproduces the scores", {
  y1 <- as.matrix(lavInspect(mean_fit_1f, what = "data"))
  fs1 <- get_fs(mean_fit_1f, method = "mean")
  sm1 <- attr(fs1, "scoring_matrix")[[1]]
  expect_equal(unname(as.numeric(t(sm1 %*% t(y1)))),
               unname(as.numeric(fs1$fs_ind60)),
               tolerance = 1e-12)

  y2 <- as.matrix(lavInspect(mean_fit_2f, what = "data"))
  fs2 <- get_fs(mean_fit_2f, method = "mean")
  sm2 <- attr(fs2, "scoring_matrix")[[1]]
  expect_equal(dim(sm2), c(2L, 7L))
  expect_equal(unname(as.matrix(t(sm2 %*% t(y2)))),
               unname(as.matrix(fs2[, c("fs_ind60", "fs_dem60")])),
               tolerance = 1e-12)
})

test_that("method = 'mean' multigroup unified output", {
  fsu <- get_fs(mean_fit_hs, method = "mean")
  expect_equal(
    colnames(fsu),
    c("fs_visual", "fs_visual_se", "visual_by_fs_visual",
      "ev_fs_visual", "school")
  )
  expect_setequal(unique(fsu$school), c("Pasteur", "Grant-White"))
  # List attributes keyed by group label
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix")) {
    a <- attr(fsu, ak)
    expect_true(is.list(a), info = ak)
    expect_setequal(names(a), c("Pasteur", "Grant-White"))
    if (ak == "fsb") {
      expect_true(is.numeric(a[[1]]), info = ak)
    } else {
      expect_true(is.matrix(a[[1]]), info = ak)
    }
  }
  estg <- lavInspect(mean_fit_hs, what = "est")
  ydata <- lavInspect(mean_fit_hs, what = "data")
  for (g in names(estg)) {
    ym <- as.matrix(ydata[[g]])
    sub <- fsu[fsu$school == g, , drop = FALSE]
    # Raw scores: per-group item means
    expect_equal(unname(as.numeric(sub$fs_visual)),
                 unname(rowMeans(ym)), tolerance = 1e-12)
    # Per-group fsb == M_g colMeans(y_g)
    expect_equal(unname(attr(fsu, "fsb")[[g]]),
                 as.numeric(matrix(1 / 3, 1, ncol(ym)) %*% colMeans(ym)),
                 tolerance = 1e-12)
    # Per-group _se == sqrt(diag(per-group fsT))
    expect_equal(unname(as.numeric(sub$fs_visual_se)),
                 rep(sqrt(attr(fsu, "fsT")[[g]][1, 1]), nrow(sub)),
                 tolerance = 1e-14)
    # Per-group scoring identity
    sm_g <- attr(fsu, "scoring_matrix")[[g]]
    expect_equal(unname(as.numeric(t(sm_g %*% t(ym)))),
                 unname(as.numeric(sub$fs_visual)),
                 tolerance = 1e-12)
  }
})

test_that("method = 'mean' multigroup list format and fs_to_group_list round trip", {
  fsu <- get_fs(mean_fit_hs, method = "mean")
  fsl <- get_fs(mean_fit_hs, method = "mean", format = "list")
  expect_setequal(names(fsl), c("Pasteur", "Grant-White"))
  for (g in names(fsl)) {
    expect_s3_class(fsl[[g]], "data.frame")
    # The group label is kept as a column in list format
    expect_true("school" %in% colnames(fsl[[g]]))
    # Plain (non-list) per-group attributes in list format
    expect_true(is.matrix(attr(fsl[[g]], "fsT")))
    expect_true(is.matrix(attr(fsl[[g]], "fsL")))
    expect_true(is.numeric(attr(fsl[[g]], "fsb")))
    expect_true(is.matrix(attr(fsl[[g]], "scoring_matrix")))
  }
  # The outer list carries the per-group attributes as well
  expect_setequal(names(attr(fsl, "fsT")), c("Pasteur", "Grant-White"))

  # unified -> list-of-df: group column dropped, per-group matrix attributes
  fsl2 <- fs_to_group_list(fsu)
  expect_setequal(names(fsl2), c("Pasteur", "Grant-White"))
  for (g in names(fsl2)) {
    expect_false("school" %in% colnames(fsl2[[g]]))
    base_cols <- setdiff(colnames(fsl[[g]]), "school")
    expect_equal(unname(as.matrix(fsl2[[g]])),
                 unname(as.matrix(fsl[[g]][base_cols])),
                 ignore_attr = TRUE)
    for (ak in c("fsT", "fsL", "fsb", "scoring_matrix")) {
      expect_equal(unname(attr(fsl2[[g]], ak)), unname(attr(fsl[[g]], ak)))
    }
  }

  # list-of-df (with no group column) -> unified: values preserved, the
  # group column is added as "group"
  fsu2 <- fs_to_group_list(fsl2)
  expect_s3_class(fsu2, "data.frame")
  expect_true("group" %in% colnames(fsu2))
  expect_false("school" %in% colnames(fsu2))
  # Row order follows the group labels, as in the original unified df
  expect_equal(unname(fsu2$group), unname(fsu$school))
  expect_equal(
    unname(as.matrix(fsu2[setdiff(colnames(fsu2), "group")])),
    unname(as.matrix(fsu[setdiff(colnames(fsu), "school")])),
    tolerance = 0
  )
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix")) {
    expect_equal(unname(attr(fsu2, ak)), unname(attr(fsu, ak)),
                 ignore_attr = TRUE)
  }
})

test_that("get_fs.data.frame(method = 'mean') end-to-end, single factor auto", {
  d <- PoliticalDemocracy[c("x1", "x2", "x3")]
  fs <- get_fs(d, method = "mean")
  # Auto-derived model: f1 =~ x1 + x2 + x3
  expect_equal(
    colnames(fs),
    c("fs_f1", "fs_f1_se", "f1_by_fs_f1", "ev_fs_f1")
  )
  y <- as.matrix(d)
  expect_equal(unname(as.numeric(fs$fs_f1)), unname(rowMeans(y)),
               tolerance = 1e-12)
  # Same scores as fitting the model explicitly
  fit_auto <- cfa("f1 =~ x1 + x2 + x3", data = d)
  fs_fit <- get_fs(fit_auto, method = "mean")
  expect_equal(unname(as.numeric(fs$fs_f1)),
               unname(as.numeric(fs_fit$fs_f1)),
               tolerance = 0)
})

test_that("get_fs.data.frame(method = 'mean') with explicit model and sum_items", {
  # x3 cross-loads onto f2, so the auto derivation would error; an
  # explicit non-duplicate assignment resolves it (x3 -> f2 only)
  d <- PoliticalDemocracy[, c("x1", "x2", "x3", "y1", "y2")]
  mod <- "f1 =~ x1 + x2 + x3\nf2 =~ y1 + y2 + x3"
  fs <- get_fs(d,
               model = mod,
               method = "mean",
               sum_items = list(f1 = c("x1", "x2"),
                                f2 = c("y1", "y2", "x3")))
  expect_equal(unname(as.numeric(fs$fs_f1)),
               unname(rowMeans(as.matrix(d[, c("x1", "x2")]))),
               tolerance = 1e-12)
   expect_equal(unname(as.numeric(fs$fs_f2)),
                unname(rowMeans(as.matrix(d[, c("y1", "y2", "x3")]))),
                tolerance = 1e-12)
})

test_that("explicit sum_items is matched by factor name, not list position", {
  # Regression: an explicitly supplied sum_items list may name the factors in
  # any order; the assignment must be matched by factor name so a reordered
  # list cannot silently swap the per-factor scoring (and fsL/fsT).
  fs_auto <- get_fs(mean_fit_2f, method = "mean")
  fwd <- list(ind60 = c("x1", "x2", "x3"), dem60 = c("y1", "y2", "y3", "y4"))
  rev <- list(dem60 = c("y1", "y2", "y3", "y4"), ind60 = c("x1", "x2", "x3"))
  fs_f <- get_fs(mean_fit_2f, method = "mean", sum_items = fwd)
  fs_r <- get_fs(mean_fit_2f, method = "mean", sum_items = rev)
  expect_equal(unname(as.matrix(fs_r)), unname(as.matrix(fs_f)),
               tolerance = 0)
  expect_equal(unname(as.matrix(fs_r)), unname(as.matrix(fs_auto)),
               tolerance = 0)
  expect_equal(attr(fs_r, "fsL")[[1]], attr(fs_f, "fsL")[[1]], tolerance = 0)
  expect_equal(attr(fs_r, "fsT")[[1]], attr(fs_f, "fsT")[[1]], tolerance = 0)
})

test_that("get_fs_lavaan(method = 'mean') wrapper", {
  fs_u <- get_fs(mean_fit_1f, method = "mean")
  fs_l <- get_fs_lavaan(mean_fit_1f, method = "mean")
  # Single-group list format: a plain data frame, not a list
  expect_s3_class(fs_l, "data.frame")
  # Plain matrix/vector attributes (no length-1 list wrapper)
  expect_true(is.matrix(attr(fs_l, "fsT")))
  expect_true(is.matrix(attr(fs_l, "fsL")))
  expect_true(is.numeric(attr(fs_l, "fsb")))
  expect_true(is.matrix(attr(fs_l, "scoring_matrix")))
  expect_equal(unname(as.matrix(fs_l)), unname(as.matrix(fs_u)),
               tolerance = 0)
  # sum_items is forwarded
  fs_l2 <- get_fs_lavaan(
    mean_fit_2f, method = "mean",
    sum_items = list(ind60 = c("x1", "x2", "x3"),
                     dem60 = c("y1", "y2", "y3", "y4"))
  )
  fs_u2 <- get_fs(mean_fit_2f, method = "mean")
  expect_equal(unname(as.matrix(fs_l2)), unname(as.matrix(fs_u2)),
               tolerance = 0)
})

test_that(
  "method = 'mean' rejects FIML fits with missing data but accepts listwise",
  {
    set.seed(42)
    d_na <- PoliticalDemocracy
    d_na[51:nrow(d_na), "x1"] <- NA
    expect_true(anyNA(d_na[["x1"]]))
    n_complete <- sum(complete.cases(d_na[, c("x1", "x2", "x3")]))

    # FIML keeps the missing-data machinery -> the gate must fire
    fit_fiml <- cfa(mean_mod_1f, data = d_na, missing = "FIML")
    expect_error(get_fs(fit_fiml, method = "mean"),
                 "does not support models fitted with missing data")

    # Listwise (the default) deletes the NAs, so no missing-data pattern
    # is recorded and "mean" legitimately passes the gate
    fit_lw <- cfa(mean_mod_1f, data = d_na)
    fs_lw <- get_fs(fit_lw, method = "mean")
    expect_identical(nrow(fs_lw), n_complete)
    y_lw <- as.matrix(lavInspect(fit_lw, what = "data"))
    expect_equal(unname(as.numeric(fs_lw$fs_ind60)), unname(rowMeans(y_lw)),
                 tolerance = 1e-12)

    # A FIML fit with COMPLETE data carries a single all-TRUE missing-data
    # pattern, so the gate must not fire
    fit_fc <- cfa(mean_mod_1f, data = PoliticalDemocracy, missing = "FIML")
    fs_fc <- get_fs(fit_fc, method = "mean")
    expect_identical(nrow(fs_fc), nrow(PoliticalDemocracy))
    y_fc <- as.matrix(lavInspect(fit_fc, what = "data"))
    est_fc <- lavInspect(fit_fc, what = "est")
    M_fc <- matrix(1 / 3, nrow = 1, ncol = 3)
    expect_equal(unname(as.numeric(fs_fc$fs_ind60)), unname(rowMeans(y_fc)),
                 tolerance = 1e-12)
    # FIML estimates item intercepts (~ sample means, to optimizer
    # tolerance), so fsb = M nu (the mean of the intercepts) holds exactly;
    # this fit has alpha = 0, so the score column mean also matches up to
    # that tolerance.
    expect_equal(
      unname(attr(fs_fc, "fsb")[[1]]),
      unname(as.numeric(M_fc %*% est_fc$nu)),
      tolerance = 1e-12
    )
    expect_lt(
      abs(unname(attr(fs_fc, "fsb")[[1]]) -
            unname(mean(as.numeric(fs_fc$fs_ind60)))),
      1e-6
    )
  }
)

test_that("method = 'mean' is rejected together with the corrected-FS options", {
  expect_error(
    get_fs(mean_fit_1f, method = "mean", corrected_fsT = TRUE),
    "not supported together with: corrected_fsT"
  )
  expect_error(
    get_fs(mean_fit_1f, method = "mean", vfsLT = TRUE),
    "not supported together with: vfsLT"
  )
  expect_error(
    get_fs(mean_fit_1f, method = "mean", reliability = TRUE),
    "not supported together with: reliability"
  )
  expect_error(
    get_fs(mean_fit_1f, method = "mean", prior_mean = 0),
    "not supported together with: prior_mean"
  )
  expect_error(
    get_fs(mean_fit_1f, method = "mean", prior_cov = 1),
    "not supported together with: prior_cov"
  )
  # Two at once: the message lists both
  expect_error(
    get_fs(mean_fit_1f, method = "mean", corrected_fsT = TRUE,
           vfsLT = TRUE),
    "not supported together with: corrected_fsT, vfsLT"
  )
})

test_that("explicit sum_items equals auto derivation on a clean CFA", {
  fs_auto <- get_fs(mean_fit_2f, method = "mean")
  fs_exp <- get_fs(
    mean_fit_2f, method = "mean",
    sum_items = list(ind60 = c("x1", "x2", "x3"),
                     dem60 = c("y1", "y2", "y3", "y4"))
  )
  expect_equal(unname(as.matrix(fs_exp)), unname(as.matrix(fs_auto)),
               tolerance = 0)
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix")) {
    expect_equal(unname(attr(fs_exp, ak)[[1]]),
                 unname(attr(fs_auto, ak)[[1]]),
                 tolerance = 0)
  }
})

test_that(
  "cross-loading model: auto derivation errors naming the item; sum_items resolves it",
  {
    d <- PoliticalDemocracy[, c("x1", "x2", "x3", "y1", "y2")]
    mod <- "f1 =~ x1 + x2 + x3\nf2 =~ y1 + y2 + x3"
    fit <- cfa(mod, data = d)
    # x3 loads on both factors -> auto derivation must name it
    expect_error(get_fs(fit, method = "mean"),
                 "load on more than one factor.*x3")
    # Assigning x3 to BOTH sums is a duplicate by design and so is rejected
    expect_error(
      get_fs(fit, method = "mean",
             sum_items = list(f1 = c("x1", "x2", "x3"),
                              f2 = c("y1", "y2", "x3"))),
      "assigned to more than one sum.*x3"
    )
    # Assigning x3 to exactly one sum resolves the ambiguity
    fs_ok <- get_fs(
      fit, method = "mean",
      sum_items = list(f1 = c("x1", "x2"),
                       f2 = c("y1", "y2", "x3"))
    )
    expect_equal(unname(as.numeric(fs_ok$fs_f1)),
                 unname(rowMeans(as.matrix(d[, c("x1", "x2")]))),
                 tolerance = 1e-12)
    expect_equal(unname(as.numeric(fs_ok$fs_f2)),
                 unname(rowMeans(as.matrix(d[, c("y1", "y2", "x3")]))),
                 tolerance = 1e-12)
  }
)

test_that("stage-2 tspa() runs with method = 'mean' scores", {
  fs2 <- get_fs(mean_fit_2f, method = "mean")
  fit2s <- tspa(
    model = "dem60 ~ ind60",
    data = fs2,
    fsT = attr(fs2, "fsT"),
    fsL = attr(fs2, "fsL"),
    fsb = attr(fs2, "fsb")
  )
  expect_s4_class(fit2s, "lavaan")
  expect_true(is.finite(coef(fit2s)["dem60~ind60"]))
})

test_that("'mean' is a distinct method, not an alias of the others", {
  # normalize_fs_method() identity: "mean" stays "mean"
  expect_identical(R2spa:::normalize_fs_method("mean"), "mean")
  expect_identical(R2spa:::normalize_fs_method("ML"), "Bartlett")
  expect_identical(R2spa:::normalize_fs_method("EB"), "regression")
  # "ML" still means Bartlett scores (bit-exact vs lavPredict)
  fs_ml <- get_fs(mean_fit_1f, method = "ML")
  bar <- as.vector(lavPredict(mean_fit_1f, method = "Bartlett"))
  expect_equal(unname(as.numeric(fs_ml$fs_ind60)), unname(bar),
               tolerance = 1e-12)
  # "mean" scores are the raw item means, distinct from both
  fs_mean <- get_fs(mean_fit_1f, method = "mean")
  reg <- as.vector(lavPredict(mean_fit_1f, method = "regression"))
  expect_false(isTRUE(all.equal(unname(as.numeric(fs_mean$fs_ind60)),
                                unname(bar), tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(unname(as.numeric(fs_mean$fs_ind60)),
                                unname(reg), tolerance = 1e-6)))
})
