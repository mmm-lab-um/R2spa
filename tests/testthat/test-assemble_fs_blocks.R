# Hand-built block tests for assemble_fs_blocks()
# No lavaan dependency — only base R and the internal helper.

make_test_fs <- function(n, q, lv_names, scores = NULL) {
  if (is.null(scores)) {
    scores <- matrix(runif(n * q, -1, 1), nrow = n, ncol = q)
  }
  colnames(scores) <- lv_names
  fsL <- diag(q)
  colnames(fsL) <- lv_names
  rownames(fsL) <- lv_names
  attr(scores, "fsL") <- fsL
  attr(scores, "fsb") <- setNames(rep(0, q), lv_names)
  attr(scores, "scoring_matrix") <- diag(q)
  scores
}

make_test_block <- function(case_idx, n, q, lv_names, diff_attrs = FALSE,
                            offset = 0, pat_label = NULL, pat = NULL) {
  fs <- make_test_fs(length(case_idx), q, lv_names)
  mult <- if (diff_attrs) (1 + offset) else 1
  fsT <- diag(mult * 0.1, nrow = q)
  rownames(fsT) <- colnames(fsT) <- paste0("fs_", lv_names)
  list(
    case_idx = case_idx,
    fs = fs,
    fsT = fsT,
    pat_label = pat_label,
    pat = pat
  )
}

make_blocks_simple <- function() {
  setNames(list(list(
    make_test_block(1:5, n = 5, q = 1, lv_names = "visual")
  )), "")
}

make_blocks_two_patterns <- function() {
  lv <- "visual"
  setNames(list(list(
    make_test_block(
      c(1L, 2L, 3L, 4L, 5L), 5, 1, lv, diff_attrs = FALSE,
      pat_label = "x1+x2",
      pat = c(x1 = TRUE, x2 = TRUE, x3 = FALSE)
    ),
    make_test_block(
      c(6L, 7L, 8L), 3, 1, lv, diff_attrs = FALSE,
      pat_label = "x3",
      pat = c(x1 = FALSE, x2 = FALSE, x3 = TRUE)
    )
  )), "")
}

make_blocks_diff_attrs <- function() {
  lv <- "visual"
  setNames(list(list(
    make_test_block(c(1L, 2L, 3L), 3, 1, lv, diff_attrs = FALSE, offset = 0),
    make_test_block(c(4L, 5L), 2, 1, lv, diff_attrs = TRUE, offset = 1)
  )), "")
}

make_blocks_two_group <- function() {
  lv <- "visual"
  list(
    "VHS" = list(make_test_block(1:4, 4, 1, lv)),
    "GW" = list(make_test_block(1:3, 3, 1, lv))
  )
}

make_blocks_multi_factor <- function() {
  lvs <- c("visual", "speed")
  setNames(list(list(make_test_block(1:3, 3, 2, lvs))), "")
}

test_that("single block assembles (unified)", {
  blks <- make_blocks_simple()
  res <- assemble_fs_blocks(blks, format = "unified")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
  expect_named(res, c("fs_visual", "fs_visual_se", "visual_by_fs_visual",
                      "ev_fs_visual"), ignore.order = FALSE)
  expect_false("group" %in% names(res))
})

test_that("single block assembles (list)", {
  blks <- make_blocks_simple()
  res <- assemble_fs_blocks(blks, format = "list")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
  expect_false("group" %in% names(res))
})

test_that("two blocks same group ordered correctly (unified)", {
  blks <- make_blocks_two_patterns()
  res <- assemble_fs_blocks(blks, format = "unified")
  expect_equal(nrow(res), 8L)
  expect_equal(names(res),
               c("fs_visual", "fs_visual_se", "visual_by_fs_visual",
                 "ev_fs_visual"))
  expect_false("group" %in% names(res))
  expect_false(any(is.na(res$fs_visual)))
})

test_that("two blocks same group ordered correctly (list)", {
  blks <- make_blocks_two_patterns()
  res <- assemble_fs_blocks(blks, format = "list")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 8L)
  expect_false(any(is.na(res$fs_visual)))
})

test_that("differing block attributes become a per-pattern named list", {
  blks <- make_blocks_diff_attrs()
  expect_no_message({
    res <- assemble_fs_blocks(blks, format = "unified")
  })
  fst <- attr(res, "fsT")
  pats <- fst[[1]]
  expect_type(pats, "list")
  # hand-built blocks lack pat_label, so labels fall back to pattern_<i>
  expect_equal(names(pats), c("pattern_1", "pattern_2"))
  # each block's fsT is preserved under its own pattern's entry
  expect_equal(pats[["pattern_1"]][1, 1], 0.1, tolerance = 1e-10)
  expect_equal(pats[["pattern_2"]][1, 1], 0.2, tolerance = 1e-10)
  expect_equal(nrow(res), 5L)
})

test_that("single block keeps plain-matrix attributes and one-label fs_pattern", {
  blks <- make_blocks_simple()
  res <- assemble_fs_blocks(blks, format = "unified")
  expect_true(is.matrix(attr(res, "fsT")[[1]]))
  expect_true(is.matrix(attr(res, "fsL")[[1]]))
  fp <- attr(res, "fs_pattern")[[1]]
  expect_equal(fp$label, rep("pattern_1", 5))
  expect_null(fp$pat)
})

test_that("fs_pattern records per-row pattern membership and pat matrix", {
  blks <- make_blocks_two_patterns()
  res <- assemble_fs_blocks(blks, format = "unified")
  fp <- attr(res, "fs_pattern")[[1]]
  expect_equal(fp$label, c(rep("x1+x2", 5L), rep("x3", 3L)))
  expect_equal(fp$pat, cbind(
    `x1+x2` = c(x1 = TRUE, x2 = TRUE, x3 = FALSE),
    `x3` = c(x1 = FALSE, x2 = FALSE, x3 = TRUE)
  ))
  expect_equal(rownames(fp$pat), c("x1", "x2", "x3"))
  expect_equal(colnames(fp$pat), c("x1+x2", "x3"))
})

test_that("two groups unified shape", {
  blks <- make_blocks_two_group()
  res <- assemble_fs_blocks(blks, format = "unified")
  expect_equal(nrow(res), 7L)
  expect_true("group" %in% names(res))
  expect_equal(sort(unique(res$group)), c("GW", "VHS"))
  expect_equal(sum(res$group == "VHS"), 4L)
  expect_equal(sum(res$group == "GW"), 3L)
  expect_equal(attr(res, "group_col"), "group")
})

test_that("two groups unified uses group_col when supplied", {
  blks <- make_blocks_two_group()
  res <- assemble_fs_blocks(blks, format = "unified", group_col = "school")
  expect_equal(nrow(res), 7L)
  expect_true("school" %in% names(res))
  expect_false("group" %in% names(res))
  expect_equal(sort(unique(res$school)), c("GW", "VHS"))
  expect_equal(sum(res$school == "VHS"), 4L)
  expect_equal(sum(res$school == "GW"), 3L)
  expect_equal(attr(res, "group_col"), "school")
})

test_that("two groups list shape", {
  blks <- make_blocks_two_group()
  res <- assemble_fs_blocks(blks, format = "list", group_col = "school")
  expect_type(res, "list")
  expect_equal(length(res), 2L)
  expect_equal(names(res), c("VHS", "GW"))
  expect_true("school" %in% names(res[[1]]))
  expect_true("school" %in% names(res[[2]]))
  expect_equal(nrow(res[[1]]), 4L)
  expect_equal(nrow(res[[2]]), 3L)
})

test_that("list format carries group-level attributes on outer list", {
  blks <- make_blocks_two_group()
  res <- assemble_fs_blocks(blks, format = "list", group_col = "school")
  for (ak in c("fsT", "fsL", "fsb", "scoring_matrix", "fs_pattern")) {
    outer_attr <- attr(res, ak)
    expect_type(outer_attr, "list")
    expect_equal(names(outer_attr), c("VHS", "GW"))
  }
})

test_that("multi-factor block assembles correctly", {
  blks <- make_blocks_multi_factor()
  res <- assemble_fs_blocks(blks, format = "unified")
  expect_equal(nrow(res), 3L)
  fs_cols <- grep("^fs_", names(res), value = TRUE)[-grep("_se$", names(res))]
  se_cols <- grep("_se$", names(res), value = TRUE)
  ev_cols <- grep("^ev_", names(res), value = TRUE)
  expect_equal(length(fs_cols), 2L)
  expect_equal(length(se_cols), 2L)
  expect_equal(length(ev_cols), 2L)
})

test_that("block attributes carried through (single factor)", {
  blks <- make_blocks_simple()
  res <- assemble_fs_blocks(blks, format = "unified")
  fst <- attr(res, "fsT")
  fsl <- attr(res, "fsL")
  fsb <- attr(res, "fsb")
  fsm <- attr(res, "scoring_matrix")
  expect_type(fst, "list")
  expect_true(is.matrix(fst[[1]]))
  expect_true(is.matrix(fsl[[1]]))
  expect_true(is.vector(fsb[[1]]))
  expect_true(is.matrix(fsm[[1]]))
})

test_that("unified attributes keyed by group label", {
  blks <- make_blocks_two_group()
  res <- assemble_fs_blocks(blks, format = "unified")
  fst <- attr(res, "fsT")
  expect_equal(names(fst), c("VHS", "GW"))
})

test_that("single-factor column naming follows convention", {
  blks <- make_blocks_simple()
  res <- assemble_fs_blocks(blks, format = "list")
  names_vec <- names(res)
  expect_true(any(grepl("^fs_", names_vec)))
  expect_true(any(grepl("_se$", names_vec)))
  expect_true(any(grepl("_by_", names_vec)))
  expect_true(any(grepl("^ev_", names_vec)))
})
