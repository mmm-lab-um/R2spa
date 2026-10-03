# tspa() per-unit factor-score pooling (PLAN 09).
#
# When the per-unit fsL/fsT/fsb attributes are heterogeneous -- per-pattern
# values from a group fitted with missing data (missing = "fiml") or
# per-cluster values from a merMod fit -- tspa() resolves them to
# individual-specific long-form values and reduces them to a single
# representative per-group set (reduce = "mean" | "median") that feeds the
# stage-2 model; the pooled values (not the nested/per-cluster ones) are what
# get attached to the returned fit. These tests validate the pooled values
# against hand computations, the no-op (identity) behavior on
# complete-data inputs, the se_fs single-factor FIML path, the PSD guard for
# reduce = "median", the all-NA-row exclusion, and the mirt per-obs
# out-of-scope guard.
library(lavaan)
library(lme4)

# ---- shared fixtures ------------------------------------------------------
# Two-factor CFA on HolzingerSwineford1939 with NA injection (fixed seed),
# single-group and multigroup FIML fits.
hs_model_2f <- "visual =~ x1 + x2 + x3
                speed  =~ x7 + x8 + x9"

mk_fiml_2f_hs <- function() {
  d <- HolzingerSwineford1939
  set.seed(1334)
  d$x2[!rbinom(nrow(d), 1L, 0.4)] <- NA
  d$x8[!rbinom(nrow(d), 1L, 0.4)] <- NA
  d
}
fiml_d_2f <- mk_fiml_2f_hs()

fit_fiml_2f_sg <- suppressWarnings(
  cfa(hs_model_2f, data = fiml_d_2f, missing = "fiml")
)
fs_fiml_2f_sg <- get_fs(fit_fiml_2f_sg)

fit_fiml_2f_mg <- suppressWarnings(
  cfa(hs_model_2f, data = fiml_d_2f, group = "school", missing = "fiml")
)
fs_fiml_2f_mg <- get_fs(fit_fiml_2f_mg)

# merMod: two random effects per cluster (sleepstudy -- one row per cluster).
lmod_pool <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fs_mer_pool <- get_fs(lmod_pool)

# Two-factor CFA with NA injected on EVERY indicator (fixed seed): produces
# all-missing (NA-pattern) rows, which get_fs() keeps and pooling must
# exclude from the reduction.
mk_fiml_all_hs <- function() {
  d <- HolzingerSwineford1939
  set.seed(1334)
  for (x in c("x1", "x2", "x3", "x7", "x8", "x9")) {
    d[[x]][!rbinom(nrow(d), 1L, 0.35)] <- NA
  }
  d
}
fit_fiml_all2f <- suppressWarnings(
  cfa(hs_model_2f, data = mk_fiml_all_hs(), missing = "fiml")
)
fs_fiml_all2f <- get_fs(fit_fiml_all2f)

# ---- helpers --------------------------------------------------------------
# Recover the full symmetric q x q matrix from its lower-triangle entries in
# the r2spa ordering the package uses (i-outer / j<=i-inner, row-major).
symm_from_lower <- function(vals, n) {
  M <- matrix(0, n, n); cte <- 1L
  for (i in seq_len(n)) for (j in seq_len(i)) {
    M[i, j] <- vals[cte]; M[j, i] <- vals[cte]; cte <- cte + 1L
  }
  M
}

# Case-count-weighted (per-row, na.rm) mean of a group's per-pattern
# fsT/fsL: each pattern matrix is weighted by the number of rows carrying
# that pattern label; rows with an NA label (all-missing rows) carry no
# matrix and are simply not counted. A single-pattern group (plain matrix
# attribute) is its own weighted mean.
pat_weighted_mean <- function(T_g, L_g, label) {
  if (is.matrix(T_g)) {
    return(list(fsT = T_g, fsL = L_g))
  }
  cnt <- c(table(factor(label, levels = names(T_g))))
  mean_of <- function(M_g) {
    Reduce("+", lapply(names(M_g), function(p) M_g[[p]] * cnt[p])) / sum(cnt)
  }
  list(fsT = mean_of(T_g), fsL = mean_of(L_g))
}

# Evaluate `expr`, collecting (and muffling) every warning message.
run_collecting_warnings <- function(expr) {
  warns <- character()
  out <- withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(fit = out, warnings = warns)
}

# ============================================================================
# 1. FIML missing data (per-pattern) -- multi-factor
# ============================================================================

test_that("tspa(): FIML missing data (SG) -- pooled fsT/fsL are the case-count-weighted per-row pattern mean (PSD)", {
  Tpat <- attr(fs_fiml_2f_sg, "fsT")[[1L]]
  Lpat <- attr(fs_fiml_2f_sg, "fsL")[[1L]]
  fp <- attr(fs_fiml_2f_sg, "fs_pattern")[[1L]]$label
  expect_gt(length(Tpat), 1L)   # genuinely multiple patterns
  expect_equal(sum(!is.na(fp)), nrow(fs_fiml_2f_sg))  # no all-NA rows here

  fit <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_2f_sg,
         fsT = attr(fs_fiml_2f_sg, "fsT"), fsL = attr(fs_fiml_2f_sg, "fsL"))
  )
  T_pooled <- attr(fit, "fsT")
  L_pooled <- attr(fit, "fsL")
  # a plain symmetric 2-D matrix is attached (not a nested pattern list)
  expect_true(is.matrix(T_pooled))
  expect_equal(dim(T_pooled), c(2L, 2L))
  expect_equal(T_pooled, t(T_pooled), tolerance = 1e-12)
  expect_true(is.matrix(L_pooled))
  expect_identical(attr(fit, "pooled_fs"), "mean")
  # PSD: the mean of PSD matrices stays PSD
  expect_gt(min(eigen(T_pooled, symmetric = TRUE, only.values = TRUE)$values),
            0)
  # equals the per-row (na.rm) mean of the per-pattern fsT / fsL
  ref <- pat_weighted_mean(Tpat, Lpat, fp)
  expect_equal(T_pooled, ref$fsT, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(L_pooled, ref$fsL, tolerance = 1e-8, ignore_attr = TRUE)
  # cross-check via the lower-triangle reconstruction (r2spa ev ordering)
  q <- nrow(T_pooled)
  lo <- numeric(q * (q + 1L) / 2L)
  cte <- 1L
  for (i in seq_len(q)) {
    for (j in seq_len(i)) {
      lo[cte] <- sum(table(factor(fp, levels = names(Tpat))) *
                       unname(vapply(names(Tpat), function(p) Tpat[[p]][i, j],
                                     numeric(1L)))) /
        sum(!is.na(fp))
      cte <- cte + 1L
    }
  }
  expect_equal(T_pooled, symm_from_lower(lo, q), tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("tspa(): FIML missing data (MG) -- pooled per group; one plain matrix per group", {
  Tpat <- attr(fs_fiml_2f_mg, "fsT")
  Lpat <- attr(fs_fiml_2f_mg, "fsL")
  pats_per_group <- vapply(Tpat, function(x) if (is.list(x)) length(x) else 1L,
                           integer(1L))
  expect_true(all(pats_per_group > 1L),
              label = "every group carries multiple patterns here")

  fit <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_2f_mg, group = "school",
         fsT = Tpat, fsL = Lpat)
  )
  T_pooled <- attr(fit, "fsT")
  expect_true(is.list(T_pooled))
  expect_equal(names(T_pooled), names(Tpat))   # one entry per group
  expect_equal(length(T_pooled), 2L)
  # flat: one plain matrix per group, no residual nesting
  expect_true(all(vapply(T_pooled, is.matrix, logical(1))))
  expect_identical(attr(fit, "pooled_fs"), "mean")
  for (g in names(T_pooled)) {
    lab_g <- attr(fs_fiml_2f_mg, "fs_pattern")[[g]]$label
    ref_g <- pat_weighted_mean(Tpat[[g]], Lpat[[g]], lab_g)
    expect_equal(T_pooled[[g]], ref_g$fsT, tolerance = 1e-8,
                 ignore_attr = TRUE, label = paste("pooled fsT for", g))
    expect_equal(attr(fit, "fsL")[[g]], ref_g$fsL, tolerance = 1e-8,
                 ignore_attr = TRUE, label = paste("pooled fsL for", g))
    expect_gt(min(eigen(T_pooled[[g]], symmetric = TRUE, only.values = TRUE)$values),
              0, label = paste("PSD for", g))
  }
})

# ============================================================================
# 2. merMod (per-cluster)
# ============================================================================

test_that("tspa(): merMod per-cluster fsT/fsL -- pooled to the unweighted cluster mean (2-D)", {
  T3d <- attr(fs_mer_pool, "fsT")
  L3d <- attr(fs_mer_pool, "fsL")
  expect_equal(length(dim(T3d)), 3L)
  expect_true(R2spa:::is_per_unit_fs(T3d, L3d))

  fit <- suppressWarnings(
    tspa("u1 ~ u0", data = fs_mer_pool, fsT = T3d, fsL = L3d)
  )
  T_pooled <- attr(fit, "fsT")
  # plain 2-D matrix attached (the 3-D input previously broke upper.tri())
  expect_true(is.matrix(T_pooled))
  expect_equal(dim(T_pooled), c(2L, 2L))
  expect_identical(attr(fit, "pooled_fs"), "mean")
  # unweighted (one weight per cluster) mean of the per-cluster slices
  T_ref <- apply(T3d, c(1L, 2L), mean)
  L_ref <- apply(L3d, c(1L, 2L), mean)
  expect_equal(T_pooled, T_ref, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(attr(fit, "fsL"), L_ref, tolerance = 1e-8, ignore_attr = TRUE)
})

# ============================================================================
# 3. Single-factor FIML: se_fs pooling of the per-row fs_<v>_se columns
# ============================================================================

test_that("tspa(): single-factor FIML -- an explicit se_fs is used as-is (not replaced by the pooled per-row SEs)", {
  d <- HolzingerSwineford1939
  set.seed(4242)
  d$x2[!rbinom(nrow(d), 1L, 0.5)] <- NA
  fit_sf <- suppressWarnings(
    cfa("visual =~ x1 + x2 + x3", data = d, group = "school", missing = "fiml")
  )
  fs_sf <- get_fs(fit_sf)
  expect_equal(attr(fs_sf, "group_col"), "school")
  # the per-row SE genuinely varies within every group (the FIML signal
  # that would have triggered the derived-path pooling)
  se_col <- fs_sf[["fs_visual_se"]]
  varied <- vapply(split(se_col, fs_sf$school),
                   function(x) length(unique(x)) > 1L, logical(1))
  expect_true(all(varied))

  # An explicit se_fs is used as-is: the 1-row value applies to every
  # group (the complete-data single-row convention), so the rendered model
  # is the plain 0.35 render, not the per-group pooled values.
  fit <- tspa("", data = fs_sf, se_fs = c(visual = 0.35), group = "school")
  expect_equal(attr(fit, "tspaModel"),
               R2spa:::tspa_render(R2spa:::tspa_schema_sf("",
                                   data.frame(visual = 0.35)),
                                   style = "sf"))
})

test_that("tspa(): single-factor FIML -- se_fs omitted keeps the derived (pooled) behavior", {
  d <- HolzingerSwineford1939
  set.seed(4242)
  d$x2[!rbinom(nrow(d), 1L, 0.5)] <- NA
  fit_sf <- suppressWarnings(
    cfa("visual =~ x1 + x2 + x3", data = d, group = "school", missing = "fiml")
  )
  fs_sf <- get_fs(fit_sf)
  # With se_fs omitted the multi-factor derivation fires (the result carries
  # fsT/fsL attributes): the per-pattern fsT/fsL are pooled to a single
  # representative matrix per group (case-count-weighted per-row mean), and
  # the rendered model is the mf render of those pooled values.
  Tpat <- attr(fs_sf, "fsT")
  Lpat <- attr(fs_sf, "fsL")
  Bpat <- attr(fs_sf, "fsb")
  fit <- tspa("", data = fs_sf, group = "school")
  gorder <- as.character(unique(fs_sf$school))
  expect_equal(names(attr(fit, "fsT")), gorder)
  # Per-group case-count-weighted mean of the per-pattern intercepts.
  b_pooled_g <- function(g) {
    lab_g <- attr(fs_sf, "fs_pattern")[[g]]$label
    b_g <- Bpat[[g]]
    setNames(vapply(seq_along(b_g), function(i) {
      sum(table(factor(lab_g, levels = names(b_g))) * b_g[[i]]) /
        sum(!is.na(lab_g))
    }, numeric(1L)), names(b_g[[1L]]))
  }
  B_pooled <- setNames(lapply(gorder, b_pooled_g), gorder)
  for (g in gorder) {
    lab_g <- attr(fs_sf, "fs_pattern")[[g]]$label
    ref_g <- pat_weighted_mean(Tpat[[g]], Lpat[[g]], lab_g)
    expect_equal(attr(fit, "fsT")[[g]], ref_g$fsT, tolerance = 1e-8,
                 ignore_attr = TRUE, label = paste("pooled fsT for", g))
    expect_equal(attr(fit, "fsL")[[g]], ref_g$fsL, tolerance = 1e-8,
                 ignore_attr = TRUE, label = paste("pooled fsL for", g))
  }
  # byte-identical to the hand render of the pooled per-group values
  expect_equal(attr(fit, "tspaModel"),
               R2spa:::tspa_render(
                 R2spa:::tspa_schema_mf("", attr(fit, "fsT"),
                                        attr(fit, "fsL"), B_pooled),
                 style = "mf"))
})

test_that("tspa(): single-factor FIML from a list-of-frames result -- the group= argument finds the group column", {
  d <- HolzingerSwineford1939
  set.seed(4242)
  d$x2[!rbinom(nrow(d), 1L, 0.5)] <- NA
  fit_sf <- suppressWarnings(
    cfa("visual =~ x1 + x2 + x3", data = d, group = "school", missing = "fiml")
  )
  fs_list <- get_fs(fit_sf, format = "list")
  # a multigroup list result carries no group_col attribute, so the group
  # column is resolved from the group= argument (the derivation convention)
  expect_null(attr(fs_list, "group_col"))
  expect_true(is.list(fs_list) && !is.data.frame(fs_list))
  fit_list <- tspa("", data = fs_list, se_fs = c(visual = 0.35),
                   group = "school")
  # reference: the rbind'd frame with the group_col attribute; the list fit
  # must pool the per-row SEs to the same per-group values
  fs_rb <- do.call(rbind, fs_list)
  attr(fs_rb, "group_col") <- "school"
  fit_rb <- tspa("", data = fs_rb, se_fs = c(visual = 0.35),
                 group = "school")
  expect_equal(attr(fit_list, "tspaModel"), attr(fit_rb, "tspaModel"))
})

test_that("tspa(): single-factor complete data -- se_fs untouched, byte-identical to the legacy path", {
  fit_sf <- suppressWarnings(
    cfa("visual =~ x1 + x2 + x3", data = HolzingerSwineford1939,
        group = "school")
  )
  fs_c <- get_fs(fit_sf)
  # constant within-group SE gives no pooling signal
  expect_true(all(vapply(split(fs_c[["fs_visual_se"]], fs_c$school),
                         function(x) length(unique(x)) == 1L, logical(1))))
  se_in <- c(visual = 0.34)
  fit <- tspa("", data = fs_c, se_fs = se_in, group = "school")
  expect_null(attr(fit, "pooled_fs"))
  expect_equal(attr(fit, "tspaModel"),
               R2spa:::tspa_render(
                 R2spa:::tspa_schema_sf("", data.frame(visual = 0.34)),
                 style = "sf"))
})

# ============================================================================
# 4. Identity / no-op guards (complete data)
# ============================================================================

test_that("tspa(): complete SG and MG fsT/fsL -- no pooling: no marker, unchanged model string", {
  m3f <- "visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9"
  mod3 <- "visual ~ speed
           textual ~ visual + speed"
  # single group: the length-1 attribute list is NOT per-unit
  fs_csg <- get_fs(cfa(m3f, data = HolzingerSwineford1939))
  Tsg <- attr(fs_csg, "fsT")
  Lsg <- attr(fs_csg, "fsL")
  expect_false(R2spa:::is_per_unit_fs(Tsg, Lsg))
  f1 <- suppressWarnings(tspa(mod3, data = fs_csg, fsT = Tsg, fsL = Lsg))
  f2 <- suppressWarnings(tspa(mod3, data = fs_csg,
                              fsT = Tsg[[1L]], fsL = Lsg[[1L]]))
  expect_null(attr(f1, "pooled_fs"))
  expect_null(attr(f2, "pooled_fs"))
  # identical stage-2 estimates to the legacy (unwrapped) call
  expect_equal(parameterestimates(f1)[c("est", "se")],
               parameterestimates(f2)[c("est", "se")], tolerance = 1e-10)
  # the model string is exactly the pre-pooling render of the raw attrs
  expect_equal(attr(f1, "tspaModel"),
               R2spa:::tspa_render(R2spa:::tspa_schema_mf(mod3, Tsg, Lsg, NULL),
                                   style = "mf"))
  # multigroup: per-group plain matrices are NOT per-unit either
  fs_cmg <- get_fs(cfa(m3f, data = HolzingerSwineford1939, group = "school"))
  Tmg <- attr(fs_cmg, "fsT")
  Lmg <- attr(fs_cmg, "fsL")
  expect_false(R2spa:::is_per_unit_fs(Tmg, Lmg))
  fit_mg <- suppressWarnings(tspa(mod3, data = fs_cmg, group = "school",
                                  fsT = Tmg, fsL = Lmg))
  expect_null(attr(fit_mg, "pooled_fs"))
  expect_equal(attr(fit_mg, "tspaModel"),
               R2spa:::tspa_render(R2spa:::tspa_schema_mf(mod3, Tmg, Lmg, NULL),
                                   style = "mf"))
})

# ============================================================================
# 5. reduce = "median"
# ============================================================================

test_that("tspa(): reduce = 'median' -- differs from 'mean' when per-unit values vary", {
  mean_res <- run_collecting_warnings(
    tspa("visual ~ speed", data = fs_fiml_2f_sg,
         fsT = attr(fs_fiml_2f_sg, "fsT"), fsL = attr(fs_fiml_2f_sg, "fsL"))
  )
  med_res <- run_collecting_warnings(
    tspa("visual ~ speed", data = fs_fiml_2f_sg,
         fsT = attr(fs_fiml_2f_sg, "fsT"), fsL = attr(fs_fiml_2f_sg, "fsL"),
         reduce = "median")
  )
  expect_identical(attr(med_res$fit, "pooled_fs"), "median")
  T_mean <- attr(mean_res$fit, "fsT")
  T_med <- attr(med_res$fit, "fsT")
  # the per-unit values vary here, so the two reductions must not agree
  expect_false(isTRUE(all.equal(unname(T_mean), unname(T_med), tolerance = 1e-8)))
})

# A hand-built per-unit (nested per-pattern) get_fs()-like frame: attributes
# in the unified single-group shape (length-1 outer list named "" wrapping a
# named per-pattern list).
mkT2 <- function(a, b2, c2) {
  matrix(c(a, b2, b2, c2), 2L, 2L,
         dimnames = list(c("fs_a", "fs_b"), c("fs_a", "fs_b")))
}
# A hand-built per-unit (nested per-pattern) get_fs()-like single-group data
# frame: unified single-group attribute shape (length-1 outer list named ""
# wrapping a named per-pattern list), continuous score values so the
# (exactly-identified) stage-2 model converges without warnings. Pattern
# names must match the entry labels.
mk_hand_fs <- function(Ts, labels, fa, fb) {
  L2 <- as.matrix(diag(2L))
  dimnames(L2) <- list(c("fs_a", "fs_b"), c("a", "b"))
  dat <- data.frame(fs_a = fa, fs_b = fb)
  structure(
    dat,
    fsT = structure(list(Ts), names = ""),
    fsL = structure(list(setNames(lapply(seq_along(Ts), function(p) L2),
                          names(Ts))), names = ""),
    fs_pattern = structure(list(list(label = labels, pat = NULL)), names = "")
  )
}

test_that("tspa(): reduce = 'median' -- PSD-poolable contrived input: exact hand pool, no PSD warning", {
  # 15 rows of T1 and 25 rows of T2 (both PSD): the per-row (na.rm) median is
  # the hand-computed [[2, 0.5], [0.5, 3]] (PSD); the mean is
  # [[13/8, 0.5], [0.5, 21/8]] -- the two must differ.
  T1 <- mkT2(1, 0.5, 2)
  T2 <- mkT2(2, 0.5, 3)
  Ts <- setNames(list(T1, T2), c("x1+x2", "x1"))
  set.seed(2024)
  dat <- mk_hand_fs(Ts, c(rep("x1+x2", 15L), rep("x1", 25L)),
                    rnorm(40L, 0, 2.5), rnorm(40L, 0, 3.5))
  res_med <- run_collecting_warnings(
    tspa("b ~ a", data = dat, fsT = attr(dat, "fsT"), fsL = attr(dat, "fsL"),
         reduce = "median")
  )
  res_mean <- run_collecting_warnings(
    tspa("b ~ a", data = dat, fsT = attr(dat, "fsT"), fsL = attr(dat, "fsL"))
  )
  expect_length(res_med$warnings, 0L)   # no (PSD or other) warning
  expect_length(res_mean$warnings, 0L)
  T_med <- attr(res_med$fit, "fsT")
  T_mean <- attr(res_mean$fit, "fsT")
  expect_equal(T_med,
               matrix(c(2, 0.5, 0.5, 3), 2L, 2L, dimnames = dimnames(T_med)),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(T_mean,
               matrix(c(13 / 8, 0.5, 0.5, 21 / 8), 2L, 2L,
                      dimnames = dimnames(T_mean)),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_false(
    isTRUE(all.equal(unname(T_med), unname(T_mean), tolerance = 1e-8)),
    label = "median and mean must differ when the per-unit values vary"
  )
})

test_that("tspa(): reduce = 'median' -- non-PSD pooled fsT (contrived) warns; reduce = 'mean' stays PSD", {
  # 30 rows of each PSD extreme plus 10 rows of the identity-ish PSD matrix:
  # the element-wise (per-row) median is [[1, 2], [2, 1]], whose smallest
  # eigenvalue is -1, while the (convex) mean [[16/7, 13/7], [13/7, 16/7]]
  # is PSD (smallest eigenvalue 3/7).
  Ta <- mkT2(1, 2, 4)
  Tb <- mkT2(4, 2, 1)
  Tc <- mkT2(1, 1, 1)
  Ts <- setNames(list(Ta, Tb, Tc), c("x1+x2", "x1", "x2"))
  set.seed(3030)
  dat <- mk_hand_fs(Ts, c(rep("x1+x2", 30L), rep("x1", 30L), rep("x2", 10L)),
                    rnorm(70L, 0, 3), rnorm(70L, 0, 3))
  med_res <- tryCatch(
    run_collecting_warnings(
      tspa("b ~ a", data = dat, fsT = attr(dat, "fsT"), fsL = attr(dat, "fsL"),
           reduce = "median")
    ),
    error = function(e) NULL
  )
  expect_true(any(grepl("positive semi-definite", med_res$warnings)),
              label = "the PSD guard must fire for the non-PSD median pool")
  mean_res <- run_collecting_warnings(
    tspa("b ~ a", data = dat, fsT = attr(dat, "fsT"), fsL = attr(dat, "fsL"))
  )
  expect_false(any(grepl("positive semi-definite", mean_res$warnings)),
               label = "the mean of PSD matrices must not trip the guard")
  T_mean <- attr(mean_res$fit, "fsT")
  expect_equal(T_mean,
               matrix(c(16 / 7, 13 / 7, 13 / 7, 16 / 7), 2L, 2L,
                      dimnames = dimnames(T_mean)),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_gt(min(eigen(T_mean, symmetric = TRUE, only.values = TRUE)$values),
            0)
})

# ============================================================================
# 6. All-missing (NA-pattern) rows
# ============================================================================

test_that("tspa(): all-missing (NA-pattern) rows are excluded from the pool, not dropped", {
  fp <- attr(fs_fiml_all2f, "fs_pattern")[[1L]]$label
  expect_gt(sum(is.na(fp)), 0L)   # the fixture genuinely has all-NA rows
  # get_fs() keeps those rows (they are present, not dropped)
  expect_equal(nrow(fs_fiml_all2f), nrow(HolzingerSwineford1939))
  Tpat <- attr(fs_fiml_all2f, "fsT")[[1L]]
  cnt <- sum(table(factor(fp, levels = names(Tpat))))
  expect_lt(cnt, nrow(fs_fiml_all2f))  # NA-labeled rows are not counted

  fit <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_all2f,
         fsT = attr(fs_fiml_all2f, "fsT"), fsL = attr(fs_fiml_all2f, "fsL"))
  )
  # the pooled value equals the hand mean over the NON-NA (labeled) rows only
  ref <- pat_weighted_mean(Tpat, attr(fs_fiml_all2f, "fsL")[[1L]], fp)
  expect_equal(attr(fit, "fsT"), ref$fsT, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(attr(fit, "fsL"), ref$fsL, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("tspa(): pooling with no scorable rows errors informatively (no misleading PSD warning)", {
  fp <- attr(fs_fiml_all2f, "fs_pattern")[[1L]]$label
  # keep only the all-missing rows (NA pattern label) and shorten the carried
  # label vector to match the rows (the get_fs() attribute/row consistency
  # invariant): nothing left to pool
  fs_na <- fs_fiml_all2f[is.na(fp), , drop = FALSE]
  fp_attr <- attr(fs_fiml_all2f, "fs_pattern")
  fp_attr[[1L]]$label <- rep(NA_character_, nrow(fs_na))
  attr(fs_na, "fs_pattern") <- fp_attr
  expect_error(
    R2spa:::pool_per_unit(fs_na, "mean", have_int = FALSE),
    "no scorable rows"
  )
})

# ============================================================================
# 7. mirt per-obs is out of scope; data must be a get_fs() result
# ============================================================================

test_that("is_per_unit_fs: mirt per-obs (flat list of bare matrices) is NOT treated as per-unit", {
  # mirt's shape: one bare matrix per observation (not nested per group)
  fsT_flat <- lapply(1:3L, function(i)
    matrix(0.5, 2, 2,
           dimnames = list(c("fs_a", "fs_b"), c("fs_a", "fs_b"))))
  fsL_flat <- lapply(1:3L, function(i)
    matrix(1, 2, 2, dimnames = list(c("fs_a", "fs_b"), c("a", "b"))))
  # a bare flat list is NOT per-unit on its own (marker-gated, PLAN 11)
  expect_false(R2spa:::is_per_unit_fs(fsT_flat, fsL_flat))
  # ...but IS per-unit when the authoritative mirt_per_obs marker is set
  expect_true(R2spa:::is_per_unit_fs(fsT_flat, fsL_flat, mirt_per_obs = TRUE))
  # positive controls: the poolable per-unit shapes
  T_nested <- structure(
    list(list(mkT2(0.5, 0.1, 0.6), mkT2(0.7, 0.1, 0.8))), names = "")
  expect_true(R2spa:::is_per_unit_fs(T_nested, fsL_flat))
  expect_true(R2spa:::is_per_unit_fs(array(0.5, c(2L, 2L, 3L)), NULL))
  expect_false(R2spa:::is_per_unit_fs(mkT2(0.5, 0.1, 0.6), fsL_flat))
})

test_that("tspa(): per-unit fsT/fsL with a non-get_fs data frame errors informatively", {
  Tn1 <- mkT2(0.1, 0.05, 0.1)
  Tn2 <- mkT2(0.15, 0.05, 0.2)
  L2 <- as.matrix(diag(2L))
  dimnames(L2) <- list(c("fs_a", "fs_b"), c("a", "b"))
  pats <- c("a+b", "a")
  # single-group attr() shape: a length-1 list wrapping a per-pattern list
  expect_error(
    tspa("b ~ a", data = data.frame(fs_a = 0:3, fs_b = 0:3),
         fsT = structure(list(setNames(list(Tn1, Tn2), pats)), names = ""),
         fsL = structure(list(setNames(list(L2, L2), pats)), names = "")),
    "get_fs\\(\\) result"
  )
  # multigroup attr() shape: one nested pattern list per group
  dat_mg <- data.frame(fs_a = 0:3, fs_b = 0:3,
                       school = c("V", "V", "G", "G"))
  expect_error(
    tspa("b ~ a", data = dat_mg, group = "school",
         fsT = setNames(rep(list(setNames(list(Tn1, Tn2), pats)), 2L),
                        c("V", "G")),
         fsL = setNames(rep(list(setNames(list(Tn1, Tn2), pats)), 2L),
                        c("V", "G"))),
    "get_fs\\(\\) result"
  )
})

# ============================================================================
# 8. Explicit measurement inputs always win (PLAN 18)
# ============================================================================

test_that("tspa(): an explicit fsb wins over the data's fsb attribute", {
  m2f <- "ind60 =~ x1 + x2 + x3
          dem60 =~ y1 + y2 + y3 + y4"
  fs2f <- get_fs(PoliticalDemocracy, model = m2f, std.lv = TRUE)
  # the data's own fsb attribute is zero (std.lv regression scores)
  expect_equal(unname(unlist(attr(fs2f, "fsb"))), c(0, 0))
  fit <- suppressWarnings(
    tspa("dem60 ~ ind60", data = fs2f,
         fsT = attr(fs2f, "fsT"), fsL = attr(fs2f, "fsL"),
         fsb = c(fs_ind60 = 2, fs_dem60 = 3))
  )
  # the supplied intercepts (not the attribute's zeros) are recorded ...
  expect_equal(attr(fit, "tspa_args")$fsb, c(fs_ind60 = 2, fs_dem60 = 3))
  # ... and carried into the rendered stage-2 model
  model_txt <- attr(fit, "tspaModel")
  expect_match(model_txt, "fs_ind60 ~ 2 \\* 1", fixed = FALSE)
  expect_match(model_txt, "fs_dem60 ~ 3 \\* 1", fixed = FALSE)
  # ...and the explicit-fsb fit differs from the derived (attribute-fsb) fit
  fit_d <- suppressWarnings(tspa("dem60 ~ ind60", data = fs2f))
  expect_false(identical(attr(fit, "tspaModel"), attr(fit_d, "tspaModel")))
})

test_that("tspa(): a constant explicit fsb is broadcast over per-unit (FIML) data", {
  # A constant (flat) fsb on per-unit (multi-pattern) data must be broadcast to
  # every pattern/group, not misread as a per-pattern intercept list (which
  # corrupted the values at two patterns and errored "subscript out of bounds"
  # with more; and was silently discarded for a multigroup list). Auto-derived
  # fsT/fsL are used here (no explicit fsT/fsL) so the original
  # auto-derivation path is exercised.
  # -- single-group FIML (multiple observed patterns) --
  n_pat <- length(attr(fs_fiml_2f_sg, "fsT")[[1L]])
  expect_true(n_pat > 1L, label = "fixture must carry >1 observed pattern")
  fit_sg <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_2f_sg,
         fsb = c(fs_visual = 2, fs_speed = 3))
  )
  expect_equal(attr(fit_sg, "tspa_args")$fsb, c(fs_visual = 2, fs_speed = 3))
  expect_match(attr(fit_sg, "tspaModel"), "fs_visual ~ 2 \\* 1", fixed = FALSE)
  expect_match(attr(fit_sg, "tspaModel"), "fs_speed ~ 3 \\* 1", fixed = FALSE)

  # -- multigroup FIML (unified data frame, group column) --
  fit_mg <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_2f_mg, group = "school",
         fsb = c(fs_visual = 2, fs_speed = 3))
  )
  fsb_mg <- attr(fit_mg, "tspa_args")$fsb
  expect_true(is.list(fsb_mg) && length(fsb_mg) == 2L)
  expect_true(all(vapply(fsb_mg,
                         function(g) isTRUE(all.equal(unname(g), c(2, 3))),
                         logical(1))))

  # -- control: the derived (no explicit fsb) fit keeps the attribute zeros --
  fit0 <- suppressWarnings(tspa("visual ~ speed", data = fs_fiml_2f_sg))
  expect_equal(unname(unlist(attr(fit0, "tspa_args")$fsb)), c(0, 0))
})

test_that("tspa(): a per-group constant fsb is broadcast within each group (FIML)", {
  # A per-group constant fsb (one flat intercept vector per group) must be
  # broadcast over each group's own patterns, not misread as a per-pattern list
  # (which corrupted the values at two patterns and errored "subscript out of
  # bounds" with more) -- in both the unified and list MG formats.
  glabs <- names(attr(fs_fiml_2f_mg, "fsb"))
  expect_length(glabs, 2L)
  pg_fsb <- setNames(
    list(c(fs_visual = 2, fs_speed = 3), c(fs_visual = 4, fs_speed = 5)),
    glabs
  )
  # unified MG data frame
  fit_u <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_2f_mg, group = "school",
         fsb = pg_fsb)
  )
  got_u <- attr(fit_u, "tspa_args")$fsb
  expect_equal(unname(got_u[[glabs[1L]]]), c(2, 3))
  expect_equal(unname(got_u[[glabs[2L]]]), c(4, 5))
  # list format
  fs_l <- get_fs(fit_fiml_2f_mg, format = "list")
  fit_l <- suppressWarnings(
    tspa("visual ~ speed", data = fs_l, group = "school", fsb = pg_fsb)
  )
  got_l <- attr(fit_l, "tspa_args")$fsb
  expect_equal(unname(got_l[[glabs[1L]]]), c(2, 3))
  expect_equal(unname(got_l[[glabs[2L]]]), c(4, 5))

  # a fully named per-group list is matched by label, so a reordered named
  # list must not be silently swapped (positional fallback would swap it)
  rev_fsb <- pg_fsb[c(2L, 1L)]
  expect_false(identical(names(rev_fsb), glabs))
  for (data in list(fs_fiml_2f_mg, fs_l)) {
    got <- suppressWarnings(
      attr(tspa("visual ~ speed", data = data, group = "school", fsb = rev_fsb),
           "tspa_args")
    )$fsb
    expect_equal(unname(got[[glabs[1L]]]), c(2, 3))
    expect_equal(unname(got[[glabs[2L]]]), c(4, 5))
  }
  # an unnamed list is matched positionally (unchanged behavior)
  got_un <- suppressWarnings(
    attr(tspa("visual ~ speed", data = fs_fiml_2f_mg, group = "school",
              fsb = unname(pg_fsb)), "tspa_args")
  )$fsb
  expect_equal(unname(got_un[[glabs[1L]]]), c(2, 3))
  expect_equal(unname(got_un[[glabs[2L]]]), c(4, 5))
  # duplicate, mismatched, and partially named group lists are rejected
  dup_fsb <- list(c(2, 3), c(4, 5))
  names(dup_fsb) <- c(glabs[1L], glabs[1L])
  expect_error(
    tspa("visual ~ speed", data = fs_fiml_2f_mg, group = "school",
         fsb = dup_fsb),
    "duplicate group names"
  )
  expect_error(
    tspa("visual ~ speed", data = fs_fiml_2f_mg, group = "school",
         fsb = setNames(list(c(2, 3), c(4, 5)), c("Other", glabs[2L]))),
    "do not match the data's groups"
  )
  part_fsb <- list(c(2, 3), c(4, 5))
  names(part_fsb) <- c(glabs[1L], "")
  expect_error(
    tspa("visual ~ speed", data = fs_fiml_2f_mg, group = "school",
         fsb = part_fsb),
    "only partially named"
  )
  # a per-group list with the wrong number of entries is rejected up front
  # (it must not be broadcast to every group / applied to the wrong intercepts)
  expect_error(
    tspa("visual ~ speed", data = fs_fiml_2f_mg, group = "school",
         fsb = list(c(2, 3))),
    "has length 1 but the model has 2 groups"
  )
})

test_that("tspa(): explicit fsb for a complete-data multigroup result (schema path)", {
  # Complete-data MG has plain per-group matrix lists (no per-unit pooling), so
  # the fsb broadcast/resolve happens in the stage-2 schema, not in
  # pool_per_unit(). A flat constant is broadcast to every group; a per-group
  # list is resolved by group label (order-independent) and a wrong count is
  # rejected. (PLAN 18; mirrors the per-unit path.)
  fit_mg_cd <- suppressWarnings(
    cfa(hs_model_2f, data = HolzingerSwineford1939, group = "school")
  )
  fs_cd <- get_fs(fit_mg_cd)
  glabs <- names(attr(fs_cd, "fsT"))
  expect_length(glabs, 2L)
  m <- "visual ~ speed"

  # (a) flat constant fsb -> broadcast to both groups (one value per group)
  fit_a <- suppressWarnings(
    tspa(m, data = fs_cd, group = "school", fsb = c(fs_visual = 2, fs_speed = 3))
  )
  ma <- attr(fit_a, "tspaModel")
  expect_match(ma, "fs_visual ~ c\\(2, 2\\) \\* 1")
  expect_match(ma, "fs_speed ~ c\\(3, 3\\) \\* 1")

  # (b) reordered NAMED per-group list -> resolved by label, not position
  pg <- setNames(list(c(fs_visual = 2, fs_speed = 3),
                      c(fs_visual = 4, fs_speed = 5)), glabs)
  exp_vis <- c(pg[[glabs[1L]]][1L], pg[[glabs[2L]]][1L])
  exp_spe <- c(pg[[glabs[1L]]][2L], pg[[glabs[2L]]][2L])
  fit_b <- suppressWarnings(
    tspa(m, data = fs_cd, group = "school", fsb = pg[c(2L, 1L)])
  )
  mb <- attr(fit_b, "tspaModel")
  expect_match(mb, paste0("fs_visual ~ c\\(", exp_vis[1L], ", ", exp_vis[2L], "\\) \\* 1"))
  expect_match(mb, paste0("fs_speed ~ c\\(", exp_spe[1L], ", ", exp_spe[2L], "\\) \\* 1"))

  # (c) mismatched-count per-group list -> rejected
  expect_error(
    tspa(m, data = fs_cd, group = "school",
         fsb = list(c(fs_visual = 2, fs_speed = 3))),
    "has length 1 but the model has 2 groups"
  )
})

test_that("tspa(): explicit per-unit fsT/fsL are pooled, not replaced by the attributes (FIML and merMod)", {
  # FIML per-pattern: double every supplied per-pattern fsT/fsL -> the pooled
  # result is exactly double the no-argument (derived) pooled result.
  dbl <- function(x) {
    lapply(x, function(g) if (is.list(g) && !is.matrix(g)) {
      lapply(g, function(m) 2 * m)
    } else {
      2 * g
    })
  }
  fit0 <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_2f_sg)
  )
  fit2 <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_2f_sg,
         fsT = dbl(attr(fs_fiml_2f_sg, "fsT")),
         fsL = dbl(attr(fs_fiml_2f_sg, "fsL")))
  )
  expect_equal(attr(fit2, "fsT"), 2 * attr(fit0, "fsT"), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(attr(fit2, "fsL"), 2 * attr(fit0, "fsL"), tolerance = 1e-10,
               ignore_attr = TRUE)
  # the explicit fit is NOT the derived (attribute-pooled) fit
  expect_false(isTRUE(all.equal(unname(attr(fit2, "fsT")),
                                unname(attr(fit0, "fsT")), tolerance = 1e-8)))

  # merMod per-cluster: same check on the 3-D array form
  T3d <- attr(fs_mer_pool, "fsT")
  L3d <- attr(fs_mer_pool, "fsL")
  fit0m <- suppressWarnings(
    tspa("u1 ~ u0", data = fs_mer_pool)
  )
  fit2m <- suppressWarnings(
    tspa("u1 ~ u0", data = fs_mer_pool, fsT = 2 * T3d, fsL = 2 * L3d)
  )
  expect_equal(attr(fit2m, "fsT"), 2 * attr(fit0m, "fsT"), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(attr(fit2m, "fsL"), 2 * attr(fit0m, "fsL"), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("tspa(): explicit per-unit fsT/fsL with an explicit fsb pool the supplied triple (FIML)", {
  # Doubling fsT/fsL and shifting fsb by a constant: the pooled fit must
  # reflect the supplied values (fsb pooled from the shifted per-pattern
  # vectors), not the data's own attributes.
  dbl <- function(x) {
    lapply(x, function(g) if (is.list(g) && !is.matrix(g)) {
      lapply(g, function(m) 2 * m)
    } else {
      2 * g
    })
  }
  Tpat <- attr(fs_fiml_2f_sg, "fsT")[[1L]]
  Lpat <- attr(fs_fiml_2f_sg, "fsL")[[1L]]
  Bpat <- attr(fs_fiml_2f_sg, "fsb")[[1L]]
  B_shift <- lapply(Bpat, function(b) b + 5)
  fp <- attr(fs_fiml_2f_sg, "fs_pattern")[[1L]]$label
  cnt <- c(table(factor(fp, levels = names(Tpat))))
  # per-latent case-count-weighted mean of the shifted per-pattern intercepts
  b_ref <- setNames(vapply(seq_len(length(B_shift[[1L]])), function(j) {
    sum(cnt * vapply(B_shift, function(b) unname(b)[j], numeric(1L))) /
      sum(cnt)
  }, numeric(1L)), names(B_shift[[1L]]))
  fit <- suppressWarnings(
    tspa("visual ~ speed", data = fs_fiml_2f_sg,
         fsT = structure(list(Tpat), names = ""),
         fsL = structure(list(Lpat), names = ""),
         fsb = structure(list(B_shift), names = ""))
  )
  # pooled fsb == the case-count-weighted mean of the shifted intercepts
  expect_equal(attr(fit, "tspa_args")$fsb, b_ref, tolerance = 1e-10)
  # and the fsT pooled from the (here unmodified) patterns is the plain mean
  ref_T <- Reduce("+", lapply(names(Tpat), function(p) Tpat[[p]] * cnt[p])) /
    sum(cnt)
  expect_equal(attr(fit, "fsT"), ref_T, tolerance = 1e-10, ignore_attr = TRUE)
})

# ============================================================================
# 9. mirt per-obs MULTI-FACTOR -> tspa() is pooled (SG + MG) (PLAN 11)
# ============================================================================

skip_if_not_installed("mirt")

set.seed(2025)
NMF <- 120L
mrt_sim2f <- function(N) as.data.frame(mirt::simdata(
  a = matrix(c(runif(4L, 0.5, 1.5), runif(4L, 0.5, 1.5)), 8L, 2L),
  d = rnorm(8L), N = N, itemtype = "2PL",
  Theta = cbind(rnorm(N), rnorm(N))))
isnaT_mf <- function(fs) vapply(attr(fs, "fsT"), function(m) all(is.na(m)),
                                logical(1))

dat_mf_sg <- mrt_sim2f(NMF); dat_mf_sg[1L, ] <- NA
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

test_that("tspa(): mirt multi-factor (SG) -- per-obs fsT/fsL == scorable row mean (PLAN 11)", {
  sc <- which(!isnaT_mf(fs_mf_sg))
  hand_T <- Reduce(`+`, attr(fs_mf_sg, "fsT")[sc]) / length(sc)
  hand_L <- Reduce(`+`, attr(fs_mf_sg, "fsL")[sc]) / length(sc)
  fit <- suppressWarnings(
    tspa("F2 ~ F1", data = fs_mf_sg,
         fsT = attr(fs_mf_sg, "fsT"), fsL = attr(fs_mf_sg, "fsL")))
  expect_identical(dim(attr(fit, "fsT")), c(2L, 2L))
  expect_equal(as.numeric(attr(fit, "fsT")), as.numeric(hand_T),
               tolerance = 1e-10)
  expect_equal(as.numeric(attr(fit, "fsL")), as.numeric(hand_L),
               tolerance = 1e-10)
})

test_that("tspa(): mirt multi-factor (MG) -- one matrix per group, first-appearance order, scorable means (PLAN 11)", {
  fit <- suppressWarnings(
    tspa("F2 ~ F1", data = fs_mf_mg, group = "group",
         fsT = attr(fs_mf_mg, "fsT"), fsL = attr(fs_mf_mg, "fsL")))
  ft <- attr(fit, "fsT")
  expect_type(ft, "list")
  expect_length(ft, 2L)
  # group order == first appearance of the group column (the lavaan
  # stage-2 order; for this blocked fixture also the level order)
  expect_identical(names(ft),
                   as.character(unique(fs_mf_mg$group[!is.na(fs_mf_mg$group)])))
  for (g in levels(fs_mf_mg$group)) {
    rows_g <- which(as.character(fs_mf_mg$group) == g)
    sc_g <- rows_g[!isnaT_mf(fs_mf_mg)[rows_g]]
    hand_T <- Reduce(`+`, attr(fs_mf_mg, "fsT")[sc_g]) / length(sc_g)
    expect_equal(as.numeric(ft[[g]]), as.numeric(hand_T), tolerance = 1e-10)
  }
})

test_that("tspa(): mirt MG with block order != level order -- pooled list follows stage-2 (first-appearance) order", {
  set.seed(2026)
  dat_i <- rbind(mrt_sim2f(NMF), mrt_sim2f(NMF))
  # mirt's data must be contiguous by group; here the first block is B and
  # the second A, while the factor's level order is A, B. mirt's stage-1
  # attribute/psi lists are in level (groupNames) order (A, B), but the
  # stage-2 lavaan fit groups in first-appearance order (B, A): the pooled
  # list must follow the stage-2 order, or the groups' matrices are
  # silently swapped.
  grp_i <- factor(rep(c("B", "A"), each = NMF), levels = c("A", "B"))
  expect_identical(as.character(levels(grp_i)), c("A", "B"))
  expect_identical(as.character(unique(grp_i)), c("B", "A"))
  mf_i <- suppressWarnings(mirt::multipleGroup(dat_i, 2L, group = grp_i,
                                               invariance = "slopes",
                                               verbose = FALSE))
  fs_i <- get_fs(mf_i)
  # the stage-1 per-group psi list is in mirt's level order (A, B)
  expect_identical(names(attr(fs_i, "psi")), c("A", "B"))
  fit <- suppressWarnings(
    tspa("F2 ~ F1", data = fs_i, group = "group",
         fsT = attr(fs_i, "fsT"), fsL = attr(fs_i, "fsL")))
  ft <- attr(fit, "fsT")
  # the pooled list must be in the stage-2 order (B, A), not level order:
  # element k feeds lavaan group k
  expect_identical(names(ft), c("B", "A"))
  # and each group's matrix is that group's own scorable-row mean
  rows_b <- which(as.character(fs_i$group) == "B")
  hand_T <- Reduce(`+`, attr(fs_i, "fsT")[rows_b]) / length(rows_b)
  expect_equal(as.numeric(ft[["B"]]), as.numeric(hand_T), tolerance = 1e-10)
  rows_a <- which(as.character(fs_i$group) == "A")
  hand_T <- Reduce(`+`, attr(fs_i, "fsT")[rows_a]) / length(rows_a)
  expect_equal(as.numeric(ft[["A"]]), as.numeric(hand_T), tolerance = 1e-10)
})

test_that("tspa(): mirt multi-factor -- completely-missing row (group NA) is excluded (PLAN 11)", {
  # the fully-NA row is reconciled to group NA + an all-NA per-row fsT, so it
  # cannot contribute to any group's finite pool
  expect_true(is.na(as.character(fs_mf_mg$group)[1L]))
  expect_true(all(is.na(unname(attr(fs_mf_mg, "fsT")[[1L]]))))
  # group A carries all of its NMF rows except that one
  nA <- sum(!is.na(as.character(fs_mf_mg$group)) &
            as.character(fs_mf_mg$group) == "A")
  expect_identical(nA, NMF - 1L)
  fit <- suppressWarnings(
    tspa("F2 ~ F1", data = fs_mf_mg, group = "group",
         fsT = attr(fs_mf_mg, "fsT"), fsL = attr(fs_mf_mg, "fsL")))
  expect_true(all(is.finite(unname(attr(fit, "fsT")[["A"]]))))
})

test_that("tspa(): mirt multi-factor (SG) -- reduce = 'median' runs (PLAN 11)", {
  fit <- suppressWarnings(suppressMessages(
    tspa("F2 ~ F1", data = fs_mf_sg,
         fsT = attr(fs_mf_sg, "fsT"), fsL = attr(fs_mf_sg, "fsL"),
         reduce = "median")))
  expect_identical(dim(attr(fit, "fsT")), c(2L, 2L))
})
