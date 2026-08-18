# get_fs() per-pattern attributes with missing data (PLAN 06)
#
# lavaan partitions each group's cases into one block per distinct
# observed-indicator pattern (@Data@Mp); get_fs() must preserve one
# fsT/fsL/fsb/scoring_matrix entry per pattern (named list, keyed by the
# observed indicators joined with "+") plus a fs_pattern attribute holding
# the per-case pattern labels and the pattern x indicator matrix.
library(lavaan)

hs <- HolzingerSwineford1939
# introduce missing data (existing test recipe; in the current lavaan
# dataset columns 7:9 are x1, x2, x3)
set.seed(1334)
hs[!rbinom(301, size = 1, prob = 0.7), 7] <- NA
hs[!rbinom(301, size = 1, prob = 0.7), 8] <- NA
hs[!rbinom(301, size = 1, prob = 0.7), 9] <- NA
hs_ind <- c("x1", "x2", "x3")
visual_model <- "visual =~ x1 + x2 + x3"

# Pattern labels taken from lavaan's internal pattern structure, used only
# to cross-check get_fs()'s own pattern names (never index-matched:
# lavaan's pattern order is case-count-descending and internal).
pattern_labels <- function(fit, gi = 1) {
  mp <- fit@Data@Mp[[gi]]
  y <- lavInspect(fit, "data")
  # single-group fits return a data frame; multigroup fits a list of them
  if (is.list(y) && !is.data.frame(y)) {
    y <- y[[gi]]
  }
  vapply(seq_len(mp$npatterns), function(m) {
    paste0(colnames(y)[mp$pat[m, ]], collapse = "+")
  }, character(1))
}

label_for_row <- function(row) {
  if (rowSums(is.na(row[, hs_ind])) == length(hs_ind)) {
    return(NA_character_)
  }
  paste0(hs_ind[!is.na(row[, hs_ind])], collapse = "+")
}

test_that("single-group missing data: per-pattern named-list attributes", {
  fit <- suppressWarnings(cfa(visual_model, data = hs, missing = "fiml"))
  expect_no_message(fs <- get_fs(fit))
  fst <- attr(fs, "fsT")
  expect_type(fst, "list")
  expect_equal(length(fst), 1L)
  for (ak in c("fsL", "fsb", "scoring_matrix")) {
    expect_type(attr(fs, ak), "list")
    expect_equal(length(attr(fs, ak)), 1L)
  }
  pats <- fst[[1]]
  expect_type(pats, "list")
  expect_gt(length(pats), 1L)
  expect_true(all(vapply(pats, is.matrix, logical(1))))
  fsl_pats <- attr(fs, "fsL")[[1]]
  lav_labels <- pattern_labels(fit)
  expect_setequal(names(pats), lav_labels)
  # Each per-pattern matrix must equal the package's canonical reference
  # computed from lavPredict's acov for the SAME pattern, matched by pattern
  # label. lavPredict orders its acov list by pattern id (not case-count),
  # so index-matching is not allowed.
  acov_list <- attr(lavPredict(fit, acov = TRUE), "acov")[[1]]
  pars <- lavInspect(fit, "est")
  # acov_list is ordered by ascending pattern id; scatter labels into that
  # order (rank() is a permutation of pattern positions)
  acov_labels <- lav_labels[order(rank(fit@Data@Mp[[1]]$id))]
  expect_equal(length(acov_list), length(pats))
  for (nm in names(pats)) {
    m <- match(nm, acov_labels)
    expect_false(is.na(m), label = sprintf("pattern %s missing from acov", nm))
    ref <- R2spa:::compute_lav_fs_matrices(
      acov_list[[m]],
      psi = pars$psi,
      alpha = pars$alpha
    )
    expect_equal(pats[[nm]], ref$fsT, tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fsl_pats[[nm]], ref$fsL, tolerance = 1e-5, ignore_attr = TRUE)
  }
})

test_that("single-group fs_pattern labels agree with the raw NA pattern", {
  fit <- suppressWarnings(cfa(visual_model, data = hs, missing = "fiml"))
  fs <- get_fs(fit)
  fp <- attr(fs, "fs_pattern")
  expect_type(fp, "list")
  expect_equal(length(fp), 1L)
  fp <- fp[[1]]
  expect_named(fp, c("label", "pat"))
  expect_equal(length(fp$label), nrow(fs))
  exp_lab <- vapply(seq_len(nrow(hs)), function(r) {
    label_for_row(hs[r, , drop = FALSE])
  }, character(1))
  expect_equal(fp$label, exp_lab)
  expect_false(anyNA(fp$label[!is.na(exp_lab)]))
  expect_equal(rownames(fp$pat), hs_ind)
  expect_setequal(colnames(fp$pat), unique(exp_lab[!is.na(exp_lab)]))
  # each pattern column marks exactly its observed indicators
  for (nm in colnames(fp$pat)) {
    r <- which(!is.na(exp_lab) & exp_lab == nm)[1]
    expect_equal(fp$pat[, nm], !is.na(hs[r, hs_ind]), ignore_attr = TRUE)
  }
})

test_that("multigroup missing data: per-group nested attributes and fs_pattern", {
  fit <- suppressWarnings(cfa(visual_model, data = hs, group = "school", missing = "fiml"))
  fs <- get_fs(fit)
  fst <- attr(fs, "fsT")
  expect_type(fst, "list")
  expect_equal(length(fst), 2L)
  grp_lbls <- names(fst)
  for (g in grp_lbls) {
    expect_type(fst[[g]], "list")
    expect_gt(length(fst[[g]]), 1L)
  }
  expect_true("school" %in% names(fs))
  fp <- attr(fs, "fs_pattern")
  expect_equal(names(fp), grp_lbls)
  fsl_all <- attr(fs, "fsL")
  gi <- match(grp_lbls, fit@Data@group.label)
  pars <- lavInspect(fit, "est", drop.list.single.group = FALSE)
  acov_all <- attr(lavPredict(fit, acov = TRUE), "acov")
  for (k in seq_along(grp_lbls)) {
    g <- grp_lbls[k]
    gidx <- which(hs$school == g)
    exp_lab <- vapply(gidx, function(r) {
      label_for_row(hs[r, , drop = FALSE])
    }, character(1))
    expect_equal(sum(fs$school == g), length(exp_lab))
    expect_equal(fp[[g]]$label, exp_lab)
    # per-pattern matrices equal the canonical reference for the same
    # pattern (matched by label; lavPredict orders acov by pattern id)
    mp <- fit@Data@Mp[[gi[k]]]
    lav_labels <- pattern_labels(fit, gi[k])
    expect_setequal(names(fst[[g]]), lav_labels)
    acov_labels <- lav_labels[order(rank(mp$id))]
    for (nm in names(fst[[g]])) {
      ref <- R2spa:::compute_lav_fs_matrices(
        acov_all[[g]][[match(nm, acov_labels)]],
        psi = pars[[gi[k]]]$psi,
        alpha = pars[[gi[k]]]$alpha
      )
      expect_equal(fst[[g]][[nm]], ref$fsT, tolerance = 1e-5,
                   ignore_attr = TRUE)
      expect_equal(fsl_all[[g]][[nm]], ref$fsL, tolerance = 1e-5,
                   ignore_attr = TRUE)
    }
  }
})

test_that("missing data format=list: attributes sit directly on the data frame", {
  fit <- suppressWarnings(cfa(visual_model, data = hs, missing = "fiml"))
  fsl <- get_fs(fit, format = "list")
  expect_s3_class(fsl, "data.frame")
  pats <- attr(fsl, "fsT")
  expect_type(pats, "list")
  expect_gt(length(pats), 1L)
  fp <- attr(fsl, "fs_pattern")
  expect_type(fp, "list")
  expect_named(fp, c("label", "pat"))
  expect_equal(length(fp$label), nrow(fsl))
})

test_that("complete data keeps plain-matrix attributes and one-pattern fs_pattern", {
  fit <- cfa(visual_model, data = HolzingerSwineford1939)
  fs <- get_fs(fit)
  fst <- attr(fs, "fsT")
  expect_equal(length(fst), 1L)
  expect_true(is.matrix(fst[[1]]))
  expect_equal(nrow(fst[[1]]), 1L)
  expect_true(is.vector(attr(fs, "fsb")[[1]]))
  expect_true(is.matrix(attr(fs, "scoring_matrix")[[1]]))
  fp <- attr(fs, "fs_pattern")[[1]]
  expect_equal(fp$label, rep("x1+x2+x3", nrow(fs)))
  expect_equal(
    fp$pat,
    matrix(TRUE, nrow = 3L, ncol = 1L, dimnames = list(hs_ind, "x1+x2+x3"))
  )
  # values unchanged relative to the lavPredict-based canonical reference
  pars <- lavInspect(fit, "est")
  acov <- attr(lavPredict(fit, acov = TRUE), "acov")
  ref <- R2spa:::compute_lav_fs_matrices(
    acov[[1]], psi = pars$psi, alpha = pars$alpha
  )
  expect_equal(fst[[1]], ref$fsT, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(attr(fs, "fsL")[[1]], ref$fsL, tolerance = 1e-8,
               ignore_attr = TRUE)
})
