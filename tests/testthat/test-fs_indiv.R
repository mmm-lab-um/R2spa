# fs_indiv() -- individual-specific (per-row) factor-score definition quantities.
#
# fs_indiv() expands the per-block fsL/fsT/fsb attributes of a get_fs()
# result into one long data frame (one row per input row), reusing
# get_fs() column naming. These tests validate its row resolution and
# per-row values across: single-group complete, single-group FIML
# missing-data, multi-group complete, priors/mean (include_intercept),
# merMod (per-cluster), the legacy_names merMod variant, unified vs list
# input, and numeric equivalence with augment_lav_predict().
library(lavaan)
library(lme4)

# ---- shared fixtures ----------------------------------------------------
hs_model_3f <- 'visual  =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed   =~ x7 + x8 + x9'
lv_names <- c("visual", "textual", "speed")

# Single-group, complete data (q = 3, so the q >= 3 ev-ordering guard is
# exercised).
fit_sg <- cfa(hs_model_3f, data = HolzingerSwineford1939)
fs_sg <- get_fs(fit_sg)
ind_sg <- fs_indiv(fs_sg)

# Prior (regression) priors, shared by the include_intercept / flow-through
# tests.
pm_p <- c(visual = 0.3, textual = -0.4, speed = 0.1)
pc_p <- matrix(c(1.5, 0.4, 0.2, 0.4, 1.3, 0.3, 0.2, 0.3, 1.1), 3, 3,
               dimnames = list(lv_names, lv_names))
fs_prio <- get_fs(fit_sg, prior_mean = pm_p, prior_cov = pc_p)

# method = "mean" (fsb is non-null, so include_intercept can emit int).
fs_mean <- get_fs(fit_sg, method = "mean")

# merMod: two random effects per cluster (sleepstudy).
lmod <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fs_mer <- get_fs(lmod)

# ---- helpers ------------------------------------------------------------
# A (possibly multi-column) row slice of a data frame as a plain, unnamed
# numeric vector, in column order.
rowcols <- function(df, i, cols) unname(as.numeric(as.matrix(df[i, cols, drop = FALSE])))
# Single-group unified (or single-group list) attribute value.
sg_attr <- function(fs, ak) attr(fs, ak)[[1L]]
# Recover the full symmetric q x q error matrix from its per-row block values
# in the two triangular orderings the package uses (fs_indiv: i-outer/j<=i
# lower-tri; augment_lav_predict: column-major upper-tri incl. diag).
symm_from_lower <- function(vals, n) {
  M <- matrix(0, n, n); cte <- 1L
  for (i in seq_len(n)) for (j in seq_len(i)) {
    M[i, j] <- vals[cte]; M[j, i] <- vals[cte]; cte <- cte + 1L
  }
  M
}
symm_from_upper <- function(vals, n) {
  M <- matrix(0, n, n)
  M[upper.tri(M, diag = TRUE)] <- vals
  M + t(M) - diag(diag(M))
}

# ============================================================================
# 1. Single-group, complete data
# ============================================================================

test_that("fs_indiv(): SG complete data -- one row per observation, scores preserved", {
  expect_equal(nrow(ind_sg), nrow(HolzingerSwineford1939))
  expect_equal(nrow(ind_sg), nrow(fs_sg))
  score_cols <- paste0("fs_", lv_names)
  expect_true(all(score_cols %in% names(ind_sg)))
  # the input score columns are carried through unmodified
  expect_equal(round(ind_sg[score_cols], 8), round(fs_sg[score_cols], 8),
               ignore_attr = TRUE)
  # a complete single-group result carries no trailing group / id columns
  expect_false("group" %in% names(ind_sg))
  expect_false("id" %in% names(ind_sg))
})

test_that("fs_indiv(): SG complete data -- per-row SE == sqrt(diag(fsT)) for every row", {
  expect_equal(nrow(sg_attr(fs_sg, "fsT")), 3L)
  ref_se <- unname(sqrt(diag(sg_attr(fs_sg, "fsT"))))
  expect_true(all(is.finite(ref_se)))
  se_cols <- paste0("fs_", lv_names, "_se")
  for (i in seq_len(nrow(ind_sg))) {
    expect_equal(rowcols(ind_sg, i, se_cols), ref_se, tolerance = 1e-8)
  }
})

test_that("fs_indiv(): SG complete data -- ev/ecov equal fsT entries in lower-tri row-major order (q >= 3)", {
  T1 <- sg_attr(fs_sg, "fsT")
  q <- nrow(T1)
  rnm <- rownames(T1); cnm <- colnames(T1)
  # i-outer / j <= i-inner (row-major lower triangle) -- the r2spa ordering
  for (i in seq_len(q)) {
    for (j in seq_len(i)) {
      cn <- if (i == j) paste0("ev_", rnm[i]) else paste0("ecov_", rnm[i], "_", cnm[j])
      expect_equal(ind_sg[[cn]][1L], T1[i, j], tolerance = 1e-8,
                   label = sprintf("fsT[%d,%d] (%s)", i, j, cn))
    }
  }
  # the guard this ordering exists for: a lower.tri() column-major read gives
  # a DIFFERENT vector for q >= 3, so it must not be what we emit.
  colmaj <- unname(T1[lower.tri(T1, diag = TRUE)])
  rowmaj <- unname(ind_sg[1L, grep("^ev_|^ecov_", names(ind_sg))])
  expect_false(
    isTRUE(all.equal(rowmaj, colmaj, tolerance = 1e-8)),
    label = "ev/ecov block should be lower-tri row-major, not the column-major triangle"
  )
})

test_that("fs_indiv(): SG complete data -- implied loadings == c(as.matrix(fsL)) column-major", {
  L1 <- sg_attr(fs_sg, "fsL")
  q <- ncol(L1)
  # values match the flat column-major read of the loading matrix
  expect_equal(
    rowcols(ind_sg, 1L, grep("_by_fs_", names(ind_sg), value = TRUE)),
    unname(c(as.matrix(L1))),
    tolerance = 1e-8
  )
  # column-major: latent j occupies columns ((j-1)*q+1 : j*q)
  expect_equal(
    grep("_by_fs_", names(ind_sg), value = TRUE),
    unlist(lapply(seq_len(q), function(j) paste(colnames(L1)[j], rownames(L1), sep = "_by_")))
  )
})

test_that("fs_indiv(): SG complete data -- matches get_fs()'s own per-row columns (same engine)", {
  # identical column sets
  expect_identical(sort(names(fs_sg)), sort(names(ind_sg)))
  # and identical values on every shared per-row column
  expect_equal(
    round(as.matrix(fs_sg[, names(fs_sg)]), 8),
    round(as.matrix(ind_sg[, names(ind_sg)]), 8),
    ignore_attr = TRUE
  )
})

# ============================================================================
# 2. Multi-group, complete data
# ============================================================================

fit_mg <- cfa(hs_model_3f, data = HolzingerSwineford1939, group = "school")
fs_mg <- get_fs(fit_mg)
ind_mg <- fs_indiv(fs_mg)

test_that("fs_indiv(): MG complete data -- per-group resolution with an aligned group column", {
  # the input's own group-column name is used (unified MG)
  expect_equal(attr(fs_mg, "group_col"), "school")
  expect_true("school" %in% names(ind_mg))
  gcol <- attr(fs_mg, "group_col")
  # every group label present in fsT appears in the group column, and the
  # per-group row counts line up
  grp_labels <- names(attr(fs_mg, "fsT"))
  expect_setequal(unique(ind_mg[[gcol]]), grp_labels)
  for (g in grp_labels) {
    Tg <- attr(fs_mg, "fsT")[[g]]
    rows_g <- which(ind_mg[[gcol]] == g)
    expect_gt(length(rows_g), 0L)
    # every row of the group uses that group's own fsT for its SEs
    ref_se <- unname(sqrt(diag(Tg)))
    for (i in rows_g) {
      expect_equal(rowcols(ind_mg, i, paste0("fs_", lv_names, "_se")), ref_se,
                   tolerance = 1e-8)
    }
    # and the (i, j) ev mapping for the group
    rnm <- rownames(Tg); cnm <- colnames(Tg); q <- nrow(Tg)
    for (i in seq_len(q)) for (j in seq_len(i)) {
      cn <- if (i == j) paste0("ev_", rnm[i]) else paste0("ecov_", rnm[i], "_", cnm[j])
      expect_equal(ind_mg[[cn]][rows_g[1L]], Tg[i, j], tolerance = 1e-8)
    }
  }
  # the two groups actually carry different latent variances
  se_sets <- unlist(lapply(grp_labels, function(g) {
    round(unname(sqrt(diag(attr(fs_mg, "fsT")[[g]]))), 6)
  }))
  expect_gt(length(unique(se_sets)), 1L)
})

# ============================================================================
# 3. FIML missing data (single pattern attribute named by pattern label)
# ============================================================================

# Multiple observed patterns but NO all-missing rows (x1 always observed, so
# the per-pattern attributes are well defined for every row).
mk_fiml_hs <- function() {
  d <- HolzingerSwineford1939
  set.seed(1334)
  d$x2[!rbinom(nrow(d), 1L, 0.6)] <- NA
  d$x3[!rbinom(nrow(d), 1L, 0.6)] <- NA
  d
}
fit_fiml <- suppressWarnings(
  cfa("visual =~ x1 + x2 + x3", data = mk_fiml_hs(), missing = "fiml")
)
fs_fiml <- get_fs(fit_fiml)

test_that("fs_indiv(): FIML missing data -- per-row SE/ev match the row's own pattern", {
  fp <- attr(fs_fiml, "fs_pattern")[[1L]]
  Tpat <- attr(fs_fiml, "fsT")[[1L]]   # named list by pattern label
  ind_f <- fs_indiv(fs_fiml)
  expect_equal(nrow(ind_f), nrow(fs_fiml))
  expect_gt(length(names(Tpat)), 1L)   # genuinely multiple patterns
  # every row's single-factor SE / ev equals its OBSERVED pattern's fsT
  ref_se <- unname(vapply(fp$label, function(l) sqrt(diag(Tpat[[l]]))[1L], numeric(1L)))
  ref_ev <- unname(vapply(fp$label, function(l) Tpat[[l]][1L, 1L], numeric(1L)))
  expect_equal(unname(ind_f[["fs_visual_se"]]), ref_se, tolerance = 1e-6)
  expect_equal(unname(ind_f[["ev_fs_visual"]]), ref_ev, tolerance = 1e-6)
  # different observed patterns give different per-row SEs (the individualized
  # per-pattern behavior)
  expect_gt(length(unique(round(ref_se, 6L))), 1L)
})

# All-missing (NA-pattern) rows: get_fs() keeps these rows with NA
# scores/SE/ev (its per-row matrix starts all-NA and only fills block rows).
# fs_indiv() must do the same -- not error and not drop the rows.
mk_fiml_all <- function() {
  d <- HolzingerSwineford1939
  set.seed(1334)
  d[!rbinom(nrow(d), 1L, 0.7), 7L] <- NA
  d[!rbinom(nrow(d), 1L, 0.7), 8L] <- NA
  d[!rbinom(nrow(d), 1L, 0.7), 9L] <- NA
  d
}
fit_fiml_all <- suppressWarnings(
  cfa("visual =~ x1 + x2 + x3", data = mk_fiml_all(), missing = "fiml")
)
fs_fiml_all <- get_fs(fit_fiml_all)

test_that("fs_indiv(): FIML missing data -- all-missing (NA-pattern) rows are not an error", {
  fp <- attr(fs_fiml_all, "fs_pattern")[[1L]]
  na_rows <- which(is.na(fp$label))
  expect_gt(length(na_rows), 0L)
  # get_fs() itself keeps these rows with NA per-row values (the reference)
  expect_true(all(is.na(fs_fiml_all[na_rows, "fs_visual_se"])))
  # fs_indiv() must agree: no error, rows preserved, NA per-row values
  ind_fall <- tryCatch(fs_indiv(fs_fiml_all), error = function(e) e)
  ok <- is.data.frame(ind_fall)
  expect_true(
    ok,
    label = if (ok) {
      "n/a"
    } else {
      paste0("fs_indiv() errored on all-missing (NA-pattern) rows, which ",
             "get_fs() returns with NA se/ev (no error): ",
             conditionMessage(ind_fall))
    }
  )
  if (ok) {
    expect_equal(nrow(ind_fall), nrow(fs_fiml_all))
    expect_true(all(is.na(ind_fall[na_rows, "fs_visual_se"])))
    expect_true(all(is.na(ind_fall[na_rows, "ev_fs_visual"])))
  }
})

# ============================================================================
# 4. method = "mean" + priors + include_intercept
# ============================================================================

test_that("fs_indiv(): priors (regression) -- include_intercept emits int_fs_* == fsb; values flow through the attrs", {
  b1 <- sg_attr(fs_prio, "fsb")
  T1p <- sg_attr(fs_prio, "fsT")
  T1v <- sg_attr(fs_sg, "fsT")
  ind_pi <- fs_indiv(fs_prio, include_intercept = TRUE)
  int_cols <- paste0("int_fs_", lv_names)
  # int columns present and equal the (regressed) intercept vector fsb
  expect_true(all(int_cols %in% names(ind_pi)))
  expect_equal(rowcols(ind_pi, 1L, int_cols), unname(as.numeric(b1)), tolerance = 1e-8)
  # se / ev reflect the prior-adjusted fsT (not the vanilla regression fsT)
  expect_equal(
    rowcols(ind_pi, 1L, paste0("fs_", lv_names, "_se")),
    unname(sqrt(diag(T1p))), tolerance = 1e-8
  )
  expect_false(
    isTRUE(all.equal(unname(sqrt(diag(T1p))), unname(sqrt(diag(T1v))), tolerance = 1e-6)),
    label = "prior-adjusted fsT should differ from the vanilla regression fsT"
  )
  # include_intercept = FALSE omits the int columns
  ind_pi0 <- fs_indiv(fs_prio, include_intercept = FALSE)
  expect_false(any(grepl("^int_", names(ind_pi0))))
})

test_that("fs_indiv(): method = 'mean' -- include_intercept emits int from the non-null fsb", {
  b1 <- sg_attr(fs_mean, "fsb")
  expect_false(is.null(b1))
  ind_mn <- fs_indiv(fs_mean, include_intercept = TRUE)
  int_cols <- paste0("int_fs_", lv_names)
  expect_true(all(int_cols %in% names(ind_mn)))
  expect_equal(rowcols(ind_mn, 1L, int_cols), unname(as.numeric(b1)), tolerance = 1e-8)
})

test_that("fs_indiv(): merMod -- include_intercept = TRUE adds NO int columns (no fsb)", {
  ind_mit <- fs_indiv(fs_mer, include_intercept = TRUE)
  expect_false(any(grepl("^int_", names(ind_mit))))
})

# ============================================================================
# 5. merMod (one row per cluster)
# ============================================================================

test_that("fs_indiv(): merMod -- one row per cluster, SE matches per-cluster fsT, id column", {
  fsTm <- attr(fs_mer, "fsT")
  n_clus <- dim(fsTm)[3L]
  clus <- dimnames(fsTm)[[3L]]
  ind_m <- fs_indiv(fs_mer)
  expect_equal(nrow(ind_m), n_clus)
  expect_equal(nrow(ind_m), nlevels(sleepstudy$Subject))
  expect_true("id" %in% names(ind_m))
  expect_false("group" %in% names(ind_m))
  # the id column is the cluster/subject level, per cluster
  expect_equal(ind_m[["id"]], clus)
  se_cols <- paste0("fs_u", seq_len(2L) - 1L, "_se")
  for (j in seq_len(n_clus)) {
    expect_equal(rowcols(ind_m, j, se_cols), unname(sqrt(diag(fsTm[, , j]))),
                 tolerance = 1e-8)
  }
})

test_that("fs_indiv(): merMod legacy_names -- u<k>_eb column naming, values unchanged", {
  fs_mleg <- get_fs(lmod, legacy_names = TRUE)
  ind_lg <- fs_indiv(fs_mleg)
  # the legacy name set is present (score / SE / ev / ecov / loading)
  expect_true(all(c("u0_eb", "u1_eb", "u0_eb_se", "u1_eb_se",
                    "ev_u0_eb", "ecov_u0_eb_u1_eb",
                    "u0_by_u0_eb", "u0_by_u1_eb") %in% names(ind_lg)))
  # values are unchanged by the naming translation
  ind_std <- fs_indiv(fs_mer)
  expect_equal(rowcols(ind_lg, 1L, "u0_eb_se"), rowcols(ind_std, 1L, "fs_u0_se"),
               tolerance = 1e-8)
  expect_equal(rowcols(ind_lg, 1L, "u0_by_u1_eb"), rowcols(ind_std, 1L, "u0_by_fs_u1"),
               tolerance = 1e-8)
  expect_equal(rowcols(ind_lg, 1L, "ecov_u0_eb_u1_eb"),
               rowcols(ind_std, 1L, "ecov_fs_u1_fs_u0"),
               tolerance = 1e-8)
})

# ============================================================================
# 6. unified vs list input
# ============================================================================

test_that("fs_indiv(): unified vs list input give identical per-row values", {
  fs_l <- get_fs(fit_sg, format = "list")
  ind_l <- fs_indiv(fs_l)
  ind_u <- ind_sg
  # the shared per-row numeric definition columns compare equal
  shared <- intersect(
    grep("_se$|_by_fs_|^ev_|^ecov_", names(ind_u), value = TRUE),
    names(ind_l)
  )
  expect_gt(length(shared), 0L)
  expect_equal(
    round(as.matrix(ind_u[, shared]), 8),
    round(as.matrix(ind_l[, shared]), 8),
    ignore_attr = TRUE
  )
})

# ============================================================================
# 7. Equivalence with augment_lav_predict()
# ============================================================================

test_that("fs_indiv(): per-row values == augment_lav_predict() after reconciling naming/order", {
  for (m in c("regression", "Bartlett")) {
    ind_e <- fs_indiv(get_fs(fit_sg, method = m))
    alp <- augment_lav_predict(fit_sg, method = m)

    # SE: canonical fs_<f>_se columns, now the same names in both
    se_ind <- grep("_se$", names(ind_e), value = TRUE)
    se_alp <- grep("_se$", names(alp), value = TRUE)
    expect_setequal(se_ind, se_alp)
    expect_equal(
      unname(as.matrix(ind_e[, se_ind])), unname(as.matrix(alp[, se_ind])),
      tolerance = 1e-8
    )

    # loadings: identical <latent>_by_fs_<score> names (column-major)
    ld_ind <- grep("_by_fs_", names(ind_e), value = TRUE)
    ld_alp <- grep("_by_fs_", names(alp), value = TRUE)
    expect_setequal(ld_ind, ld_alp)
    ord <- order(ld_alp)
    expect_equal(
      unname(as.matrix(ind_e[, ld_alp[ord]])),
      unname(as.matrix(alp[, ld_alp[ord]])),
      tolerance = 1e-8
    )

    # error terms: fs_indiv lower-tri, augment_lav_predict upper-tri -- the
    # same symmetric values in different orders (reconcile into a full
    # symmetric matrix before comparing). n is solved from evq = n(n+1)/2.
    evq <- length(grep("^ev_|^ecov_", names(ind_e), value = TRUE))
    n <- as.integer(round((sqrt(1L + 8L * evq) - 1L) / 2L))
    symm_ind <- symm_from_lower(unname(as.matrix(ind_e[1L, grep("^ev_|^ecov_", names(ind_e))])), n)
    symm_alp <- symm_from_upper(unname(as.matrix(alp[1L, grep("^ev_|^ecov_", names(alp))])), n)
    expect_equal(symm_ind, symm_alp, tolerance = 1e-8, ignore_attr = TRUE)
    # and both recover the shared per-pattern fsT
    exp_T <- sg_attr(get_fs(fit_sg, method = m), "fsT")
    expect_equal(symm_ind, exp_T, tolerance = 1e-8, ignore_attr = TRUE)
  }
})
