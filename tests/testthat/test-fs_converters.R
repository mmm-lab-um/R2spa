# Tests for fs_to_group_list() converter

# --- Hand-built fixtures (no lavaan dependency) ---

make_unified_fixture <- function() {
  # Two-group unified data frame
  df <- data.frame(
    fs_visual = c(0.1, -0.2, 0.3, 0.4, -0.1, 0.2),
    fs_visual_se = c(0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
    visual_by_fs_visual = c(1, 1, 1, 1, 1, 1),
    ev_fs_visual = c(0.01, 0.01, 0.01, 0.02, 0.02, 0.02),
    group = c("VHS", "VHS", "VHS", "VHS", "GW", "GW")
  )
  attr(df, "fsT") <- list(
    VHS = matrix(0.01, 1, 1, dimnames = list("fs_visual", "fs_visual")),
    GW = matrix(0.02, 1, 1, dimnames = list("fs_visual", "fs_visual"))
  )
  attr(df, "fsL") <- list(
    VHS = matrix(1, 1, 1, dimnames = list("fs_visual", "visual")),
    GW = matrix(1, 1, 1, dimnames = list("fs_visual", "visual"))
  )
  attr(df, "fsb") <- list(VHS = 0, GW = 0)
  attr(df, "scoring_matrix") <- list(
    VHS = matrix(1, 1, 1),
    GW = matrix(1, 1, 1)
  )
  attr(df, "group_col") <- "group"
  df
}

make_list_fixture <- function() {
  # Two-group list of data frames
  out <- list(
    VHS = data.frame(
      fs_visual = c(0.1, -0.2, 0.3, 0.4),
      fs_visual_se = c(0.15, 0.15, 0.15, 0.15),
      visual_by_fs_visual = c(1, 1, 1, 1),
      ev_fs_visual = c(0.01, 0.01, 0.01, 0.01)
    ),
    GW = data.frame(
      fs_visual = c(-0.1, 0.2),
      fs_visual_se = c(0.15, 0.15),
      visual_by_fs_visual = c(1, 1),
      ev_fs_visual = c(0.02, 0.02)
    )
  )
  attr(out[[1]], "fsT") <- matrix(0.01, 1, 1,
    dimnames = list("fs_visual", "fs_visual")
  )
  attr(out[[2]], "fsT") <- matrix(0.02, 1, 1,
    dimnames = list("fs_visual", "fs_visual")
  )
  attr(out[[1]], "fsL") <- matrix(1, 1, 1,
    dimnames = list("fs_visual", "visual")
  )
  attr(out[[2]], "fsL") <- matrix(1, 1, 1,
    dimnames = list("fs_visual", "visual")
  )
  attr(out[[1]], "fsb") <- 0
  attr(out[[2]], "fsb") <- 0
  attr(out[[1]], "scoring_matrix") <- matrix(1, 1, 1)
  attr(out[[2]], "scoring_matrix") <- matrix(1, 1, 1)
  # Outer-level attributes
  attr(out, "fsT") <- list(
    VHS = attr(out[[1]], "fsT"),
    GW = attr(out[[2]], "fsT")
  )
  attr(out, "fsL") <- list(
    VHS = attr(out[[1]], "fsL"),
    GW = attr(out[[2]], "fsL")
  )
  attr(out, "fsb") <- list(VHS = 0, GW = 0)
  attr(out, "scoring_matrix") <- list(
    VHS = attr(out[[1]], "scoring_matrix"),
    GW = attr(out[[2]], "scoring_matrix")
  )
  out
}

# --- Multi-group: unified -> list ---

test_that("multi-group unified converts to list", {
  unified <- make_unified_fixture()
  result <- fs_to_group_list(unified)

  expect_type(result, "list")
  expect_equal(length(result), 2L)
  expect_equal(names(result), c("VHS", "GW"))
  expect_equal(nrow(result$VHS), 4L)
  expect_equal(nrow(result$GW), 2L)
  expect_false("group" %in% names(result$VHS))
  expect_false("group" %in% names(result$GW))
})

test_that("multi-group unified -> list preserves per-group attributes", {
  unified <- make_unified_fixture()
  result <- fs_to_group_list(unified)

  fst_vhs <- attr(result$VHS, "fsT")
  fst_gw <- attr(result$GW, "fsT")
  expect_equal(fst_vhs[1, 1], 0.01, tolerance = 1e-10)
  expect_equal(fst_gw[1, 1], 0.02, tolerance = 1e-10)

  for (ak in c("fsL", "fsb", "scoring_matrix")) {
    expect_true(!is.null(attr(result$VHS, ak)))
    expect_true(!is.null(attr(result$GW, ak)))
  }
})

test_that("multi-group unified -> list carries outer list attributes", {
  unified <- make_unified_fixture()
  result <- fs_to_group_list(unified)

  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix")) {
    outer <- attr(result, ak)
    expect_type(outer, "list")
    expect_equal(names(outer), c("VHS", "GW"))
  }
})

# --- Multi-group: list -> unified ---

test_that("multi-group list converts to unified", {
  lst <- make_list_fixture()
  result <- fs_to_group_list(lst)

  expect_s3_class(result, "data.frame")
  expect_true("group" %in% names(result))
  expect_equal(nrow(result), 6L)
  expect_equal(unique(result$group), c("VHS", "GW"))
})

test_that("multi-group list -> unified preserves list-valued attributes", {
  lst <- make_list_fixture()
  result <- fs_to_group_list(lst)

  fst <- attr(result, "fsT")
  expect_type(fst, "list")
  expect_equal(names(fst), c("VHS", "GW"))
  expect_equal(fst$VHS[1, 1], 0.01, tolerance = 1e-10)
  expect_equal(fst$GW[1, 1], 0.02, tolerance = 1e-10)
})

# --- Round-trip ---

test_that("unified -> list -> unified round-trip", {
  unified <- make_unified_fixture()
  lst <- fs_to_group_list(unified)
  back <- fs_to_group_list(lst)

  expect_equal(nrow(back), nrow(unified))
  expect_equal(sum(back$group == "VHS"), 4L)
  expect_equal(sum(back$group == "GW"), 2L)

  # Compare data values (ignore row ordering)
  vhs_orig <- unified[unified$group == "VHS", names(unified) != "group"]
  vhs_back <- back[back$group == "VHS", names(back) != "group"]
  expect_equal(as.data.frame(vhs_back), as.data.frame(vhs_orig),
               ignore_attr = TRUE)

  # Compare attributes
  fst_orig <- attr(unified, "fsT")
  fst_back <- attr(back, "fsT")
  expect_equal(fst_back$VHS, fst_orig$VHS, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(fst_back$GW, fst_orig$GW, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("list -> unified -> list round-trip", {
  lst <- make_list_fixture()
  unified <- fs_to_group_list(lst)
  back <- fs_to_group_list(unified)

  expect_equal(nrow(back$VHS), nrow(lst$VHS))
  expect_equal(nrow(back$GW), nrow(lst$GW))
  expect_equal(as.data.frame(back$VHS), as.data.frame(lst$VHS),
               ignore_attr = TRUE)
  expect_equal(as.data.frame(back$GW), as.data.frame(lst$GW),
               ignore_attr = TRUE)
})

# --- Single-group edge case ---

test_that("single-group unified converts to plain data frame", {
  df <- data.frame(
    fs_visual = c(0.1, -0.2, 0.3),
    fs_visual_se = c(0.15, 0.15, 0.15),
    visual_by_fs_visual = c(1, 1, 1),
    ev_fs_visual = c(0.01, 0.01, 0.01),
    group = c("VHS", "VHS", "VHS")
  )
  attr(df, "fsT") <- list(
    VHS = matrix(0.01, 1, 1, dimnames = list("fs_visual", "fs_visual"))
  )
  attr(df, "fsL") <- list(
    VHS = matrix(1, 1, 1, dimnames = list("fs_visual", "visual"))
  )
  attr(df, "fsb") <- list(VHS = 0)
  attr(df, "scoring_matrix") <- list(VHS = matrix(1, 1, 1))

  result <- fs_to_group_list(df)
  expect_s3_class(result, "data.frame")
  expect_false("group" %in% names(result))
  expect_equal(nrow(result), 3L)
  # Attributes should be bare matrices, not lists
  expect_true(is.matrix(attr(result, "fsT")))
  expect_true(is.matrix(attr(result, "fsL")))
})

# --- Error handling ---

test_that("treats data frame without group column as single-group", {
  no_grp_df <- data.frame(fs_visual = 1:3)
  attr(no_grp_df, "fsT") <- list(diag(0.1))
  attr(no_grp_df, "fsL") <- list(diag(1))
  res <- fs_to_group_list(no_grp_df)
  expect_true(is.data.frame(res))
  expect_false("group" %in% names(res))
})

test_that("errors on unnamed list", {
  bad_lst <- list(data.frame(fs_visual = 1:3))
  expect_error(fs_to_group_list(bad_lst),
               "must be named")
})

test_that("errors on unsupported input type", {
  expect_error(fs_to_group_list("not_valid"),
               "must be data frame or list")
})

# --- Integration with get_fs() ---

test_that("fs_to_group_list() works on real get_fs() multigroup output", {
  hs_model <- "visual =~ x1 + x2 + x3"
  fit <- lavaan::cfa(hs_model,
                     data = HolzingerSwineford1939,
                     group = "school")
  fs_unified <- get_fs(fit)

  # Convert to list
  fs_list <- fs_to_group_list(fs_unified)
  expect_type(fs_list, "list")
  expect_equal(length(fs_list), 2L)

  # Compare with get_fs(..., format = "list")
  fs_list_direct <- get_fs(fit, format = "list")
  expect_equal(nrow(fs_list$VHS), nrow(fs_list_direct$VHS))
  expect_equal(nrow(fs_list$GW), nrow(fs_list_direct$GW))

  # Round-trip: list->unified path always uses "group" as column name
  fs_back <- fs_to_group_list(fs_list)
  expect_s3_class(fs_back, "data.frame")
  expect_true("group" %in% names(fs_back))
})

test_that("integration: get_fs() multigroup output uses original variable name", {
  hs_model <- "visual =~ x1 + x2 + x3"
  fit <- lavaan::cfa(hs_model,
                     data = HolzingerSwineford1939,
                     group = "school")
  fs_unified <- get_fs(fit)

  expect_true("school" %in% names(fs_unified))
  expect_equal(sort(unique(fs_unified$school)), c("Grant-White", "Pasteur"))
  expect_equal(attr(fs_unified, "group_col"), "school")

  # Convert to list — should split on "school" not "group"
  fs_list <- fs_to_group_list(fs_unified)
  expect_type(fs_list, "list")
  expect_equal(names(fs_list), c("Pasteur", "Grant-White"))
})

# --- Missing data: nested per-pattern attributes + fs_pattern ---
#
# VHS has two observed-indicator patterns, GW is complete (k = 1).

fs1 <- matrix(0.01, 1, 1, dimnames = list("fs_visual", "fs_visual"))
fs2 <- matrix(0.03, 1, 1, dimnames = list("fs_visual", "fs_visual"))
fs3 <- matrix(0.02, 1, 1, dimnames = list("fs_visual", "fs_visual"))
ld <- matrix(1, 1, 1, dimnames = list("fs_visual", "visual"))
sm <- matrix(1, 1, 1)
missing_pat_vhs <- cbind(
  `x1+x2` = c(x1 = TRUE, x2 = TRUE, x3 = FALSE),
  `x3` = c(x1 = FALSE, x2 = FALSE, x3 = TRUE)
)
missing_pat_gw <- matrix(
  TRUE, nrow = 3L, ncol = 1L,
  dimnames = list(c("x1", "x2", "x3"), "x1+x2+x3")
)

make_unified_missing_fixture <- function() {
  df <- data.frame(
    fs_visual = c(0.1, -0.2, 0.3, 0.4, -0.1, 0.2),
    fs_visual_se = c(0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
    visual_by_fs_visual = c(1, 1, 1, 1, 1, 1),
    ev_fs_visual = c(0.01, 0.01, 0.01, 0.02, 0.02, 0.02),
    group = c("VHS", "VHS", "VHS", "VHS", "GW", "GW")
  )
  attr(df, "fsT") <- list(VHS = list(`x1+x2` = fs1, `x3` = fs2), GW = fs3)
  attr(df, "fsL") <- list(VHS = list(`x1+x2` = ld, `x3` = ld), GW = ld)
  attr(df, "fsb") <- list(VHS = list(`x1+x2` = 0, `x3` = 0), GW = 0)
  attr(df, "scoring_matrix") <- list(
    VHS = list(`x1+x2` = sm, `x3` = sm),
    GW = sm
  )
  attr(df, "fs_pattern") <- list(
    VHS = list(
      label = c("x1+x2", "x1+x2", "x3", "x3"),
      pat = missing_pat_vhs
    ),
    GW = list(
      label = c("x1+x2+x3", "x1+x2+x3"),
      pat = missing_pat_gw
    )
  )
  attr(df, "group_col") <- "group"
  df
}

make_list_missing_fixture <- function() {
  out <- list(
    VHS = data.frame(
      fs_visual = c(0.1, -0.2, 0.3, 0.4),
      fs_visual_se = c(0.15, 0.15, 0.15, 0.15),
      visual_by_fs_visual = c(1, 1, 1, 1),
      ev_fs_visual = c(0.01, 0.01, 0.01, 0.01)
    ),
    GW = data.frame(
      fs_visual = c(-0.1, 0.2),
      fs_visual_se = c(0.15, 0.15),
      visual_by_fs_visual = c(1, 1),
      ev_fs_visual = c(0.02, 0.02)
    )
  )
  attr(out$VHS, "fsT") <- list(`x1+x2` = fs1, `x3` = fs2)
  attr(out$GW, "fsT") <- fs3
  attr(out$VHS, "fsL") <- list(`x1+x2` = ld, `x3` = ld)
  attr(out$GW, "fsL") <- ld
  attr(out$VHS, "fsb") <- list(`x1+x2` = 0, `x3` = 0)
  attr(out$GW, "fsb") <- 0
  attr(out$VHS, "scoring_matrix") <- list(`x1+x2` = sm, `x3` = sm)
  attr(out$GW, "scoring_matrix") <- sm
  attr(out$VHS, "fs_pattern") <- list(
    label = c("x1+x2", "x1+x2", "x3", "x3"),
    pat = missing_pat_vhs
  )
  attr(out$GW, "fs_pattern") <- list(
    label = c("x1+x2+x3", "x1+x2+x3"),
    pat = missing_pat_gw
  )
  out
}

test_that("missing-data unified -> list preserves nested and fs_pattern", {
  unified <- make_unified_missing_fixture()
  result <- fs_to_group_list(unified)

  expect_type(result, "list")
  expect_equal(names(result), c("VHS", "GW"))
  # VHS stays nested per pattern, GW stays a plain matrix
  expect_equal(attr(result$VHS, "fsT"), attr(unified, "fsT")$VHS)
  expect_true(is.matrix(attr(result$GW, "fsT")))
  expect_equal(attr(result$VHS, "fsb"), attr(unified, "fsb")$VHS)
  expect_equal(attr(result$VHS, "fs_pattern"), attr(unified, "fs_pattern")$VHS)
  expect_equal(attr(result$GW, "fs_pattern"), attr(unified, "fs_pattern")$GW)
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix", "fs_pattern")) {
    expect_equal(names(attr(result, ak)), c("VHS", "GW"))
  }
})

test_that("missing-data list -> unified preserves nested and fs_pattern", {
  lst <- make_list_missing_fixture()
  result <- fs_to_group_list(lst)

  expect_s3_class(result, "data.frame")
  expect_true("group" %in% names(result))
  expect_equal(nrow(result), 6L)
  expect_equal(attr(result, "fsT")$VHS, attr(lst$VHS, "fsT"))
  expect_equal(attr(result, "fsT")$GW, fs3)
  expect_equal(attr(result, "fsb")$VHS, attr(lst$VHS, "fsb"))
  expect_equal(attr(result, "fsb")$GW, 0)
  expect_equal(attr(result, "fs_pattern")$VHS, attr(lst$VHS, "fs_pattern"))
  expect_equal(attr(result, "fs_pattern")$GW, attr(lst$GW, "fs_pattern"))
})

test_that("missing-data unified <-> list round-trips both directions", {
  unified <- make_unified_missing_fixture()
  lst <- make_list_missing_fixture()

  back_u <- fs_to_group_list(fs_to_group_list(unified))
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix", "fs_pattern")) {
    expect_equal(attr(back_u, ak), attr(unified, ak))
  }
  expect_equal(as.data.frame(back_u), as.data.frame(unified),
               ignore_attr = TRUE)

  back_l <- fs_to_group_list(fs_to_group_list(lst))
  expect_type(back_l, "list")
  for (g in c("VHS", "GW")) {
    for (ak in c("fsT", "fsL", "fsb", "scoring_matrix", "fs_pattern")) {
      expect_equal(attr(back_l[[g]], ak), attr(lst[[g]], ak))
    }
    expect_equal(as.data.frame(back_l[[g]]), as.data.frame(lst[[g]]),
                 ignore_attr = TRUE)
  }
  expect_equal(names(back_l), c("VHS", "GW"))
})

test_that("single-group missing-data unified unwraps to nested values", {
  df <- data.frame(
    fs_visual = c(0.1, 0.2),
    fs_visual_se = c(0.15, 0.15),
    visual_by_fs_visual = c(1, 1),
    ev_fs_visual = c(0.01, 0.02)
  )
  attr(df, "fsT") <- list(list(`x1+x2` = fs1, `x3` = fs2))
  attr(df, "fsL") <- list(list(`x1+x2` = ld, `x3` = ld))
  attr(df, "fsb") <- list(list(`x1+x2` = 0, `x3` = 0))
  attr(df, "scoring_matrix") <- list(list(`x1+x2` = sm, `x3` = sm))
  attr(df, "fs_pattern") <- list(list(
    label = c("x1+x2", "x3"),
    pat = missing_pat_vhs
  ))

  result <- fs_to_group_list(df)
  expect_s3_class(result, "data.frame")
  expect_false("group" %in% names(result))
  pats <- attr(result, "fsT")
  expect_type(pats, "list")
  expect_equal(names(pats), c("x1+x2", "x3"))
  expect_equal(attr(result, "fs_pattern")$label, c("x1+x2", "x3"))
  expect_equal(attr(result, "fs_pattern")$pat, missing_pat_vhs)
})
