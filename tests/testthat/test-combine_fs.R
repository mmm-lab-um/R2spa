# combine_fs(): combine several get_fs() results into one block-diagonal
# multi-factor result. Each input is one block; the combined fsL/fsT/psi are
# block-diagonal across the inputs (every cross-block entry is exactly zero).
# The result is downstream-transparent (feeds tspa() via the multi-factor
# derivation) and carries a per-row (per_obs) attribute layout.
#
# Conventions: every block-diagonal result is checked against an independent
# control built from the per-input resolve_fs_per_row() blocks and
# block_diag(), so the cross-block-zero invariant is asserted, not assumed.
# mirt is Suggests-only: its fixtures are guarded and its tests skip.
library(lavaan)

hs <- HolzingerSwineford1939
n_hs <- nrow(hs)

# ---- two single-factor lavaan results (same rows) -------------------------
fit_vis <- cfa("vis =~ x1 + x2 + x3", data = hs)
fit_qua <- cfa("qua =~ x4 + x5 + x6", data = hs)
f_vis <- get_fs(fit_vis)
f_qua <- get_fs(fit_qua)
rv_vis <- R2spa:::resolve_fs_per_row(f_vis)
rv_qua <- R2spa:::resolve_fs_per_row(f_qua)
comb_vq <- combine_fs(list(f_vis, f_qua))

# a differently-named 200-row lavaan result (for the row-count-mismatch guard)
f_qua_short <- get_fs(cfa("qua =~ x4 + x5 + x6", data = hs[1:200, ]))

# ---- mirt fixtures (Suggests-only) ----------------------------------------
have_mirt <- requireNamespace("mirt", quietly = TRUE)
if (have_mirt) {
  set.seed(7)
  mirt_items <- data.frame(
    lapply(1:6, function(k) rbinom(n_hs, 1, sample(0.3:0.7))),
    check.names = FALSE
  )
  colnames(mirt_items) <- paste0("y", 1:6)
  # one 2-factor fit (a single k x k block with internal structure)
  fmirt2f <- get_fs(suppressWarnings(mirt::mirt(mirt_items, 2)))
  # two separate 1-factor fits (two 1 x 1 blocks); both default to "F1", so
  # the second is renamed to disambiguate the latent names
  fmirt1fa <- get_fs(suppressWarnings(mirt::mirt(mirt_items[, 1:3, drop = FALSE], 1)))
  fmirt1fb <- get_fs(suppressWarnings(mirt::mirt(mirt_items[, 4:6, drop = FALSE], 1)))
  # rename a single-factor mirt result's factor (only the bits combine_fs reads:
  # the per-row fsL/fsT/fsb dimnames and the fs_<f> score column)
  rename_mirt_factor <- function(f, new) {
    L <- attr(f, "fsL"); T <- attr(f, "fsT"); B <- attr(f, "fsb")
    for (i in seq_along(L)) {
      dimnames(L[[i]]) <- list(paste0("fs_", new), new)
      dimnames(T[[i]]) <- list(paste0("fs_", new), paste0("fs_", new))
      if (!is.null(B[[i]])) names(B[[i]]) <- paste0("fs_", new)
    }
    attr(f, "fsL") <- L; attr(f, "fsT") <- T; attr(f, "fsb") <- B
    names(f)[names(f) == "fs_F1"] <- paste0("fs_", new)
    f
  }
  fmirt1fb <- rename_mirt_factor(fmirt1fb, "G2")
  # a 50-observation 1-factor fit (for the id-join union guard)
  fmirt1f_short <- get_fs(
    suppressWarnings(mirt::mirt(mirt_items[1:50, 1:3, drop = FALSE], 1))
  )
}

# ===========================================================================
test_that("T1: two lavaan single-factor results combine to a block-diagonal result", {
  expect_true(is.data.frame(comb_vq))
  expect_equal(nrow(comb_vq), n_hs)
  expect_true(isTRUE(attr(comb_vq, "per_obs")))
  expect_equal(
    colnames(comb_vq),
    c("fs_vis", "fs_qua", "fs_vis_se", "fs_qua_se",
      "vis_by_fs_vis", "vis_by_fs_qua", "qua_by_fs_vis", "qua_by_fs_qua",
      "ev_fs_vis", "ecov_fs_qua_fs_vis", "ev_fs_qua")
  )

  # per-row fsL / fsT are block-diagonal: cross entries exactly zero
  Lr <- attr(comb_vq, "fsL")[[1]]
  Tr <- attr(comb_vq, "fsT")[[1]]
  expect_equal(Lr[1, 2], 0); expect_equal(Lr[2, 1], 0)
  expect_equal(Tr[1, 2], 0); expect_equal(Tr[2, 1], 0)

  # within-block entries are the inputs' own
  b1 <- rv_vis$blocks[[rv_vis$pattern_idx[1]]]
  b2 <- rv_qua$blocks[[rv_qua$pattern_idx[1]]]
  expect_equal(Lr[1, 1], b1$fsL[1, 1])
  expect_equal(Tr[1, 1], b1$fsT[1, 1])
  expect_equal(Tr[2, 2], b2$fsT[1, 1])

  # independent block_diag control for the whole per-row fsT
  expect_equal(Tr, R2spa:::block_diag(b1$fsT, b2$fsT), ignore_attr = TRUE)

  # cross-block derived columns are exactly zero
  expect_equal(comb_vq[["vis_by_fs_qua"]], rep(0, n_hs))
  expect_equal(comb_vq[["qua_by_fs_vis"]], rep(0, n_hs))
  expect_equal(comb_vq[["ecov_fs_qua_fs_vis"]], rep(0, n_hs))

  # scores match the inputs
  expect_equal(unname(comb_vq$fs_vis), unname(as.numeric(f_vis$fs_vis)))
  expect_equal(unname(comb_vq$fs_qua), unname(as.numeric(f_qua$fs_qua)))

  # psi is block-diagonal, matching block_diag of the inputs' psi
  expect_equal(
    attr(comb_vq, "psi"),
    R2spa:::block_diag(R2spa:::combine_fs_psi(f_vis),
                       R2spa:::combine_fs_psi(f_qua)),
    ignore_attr = TRUE
  )
})

# ===========================================================================
test_that("cbind-style separate arguments equal the list form", {
  expect_identical(combine_fs(f_vis, f_qua), comb_vq)
  # a list first-arg mixed with a variadic arg is rejected
  expect_error(combine_fs(list(f_vis), f_qua), "not both")
})

# ===========================================================================
test_that("T2: a 2-factor mirt result + a lavaan result give a 3 x 3 block diagonal", {
  skip_if_not_installed("mirt")
  comb <- combine_fs(list(fmirt2f, f_vis))
  n <- n_hs
  expect_equal(nrow(comb), n)
  expect_equal(ncol(comb), 3 + 3 + 3 * 3 + 3 * 4 / 2)  # 3 scores, se, ld, ev

  Lr <- attr(comb, "fsL")[[1]]
  Tr <- attr(comb, "fsT")[[1]]
  # dimnames in global latent order (mirt F1, F2 first, then lavaan vis)
  expect_equal(rownames(Lr), c("fs_F1", "fs_F2", "fs_vis"))
  expect_equal(colnames(Lr), c("F1", "F2", "vis"))
  # the mirt 2 x 2 block is the input's own (internal structure preserved)
  rm_ <- R2spa:::resolve_fs_per_row(fmirt2f)
  bm <- rm_$blocks[[rm_$pattern_idx[1]]]
  expect_equal(Lr[1:2, 1:2], bm$fsL[1:2, 1:2], tolerance = 1e-12)
  expect_equal(Tr[1:2, 1:2], bm$fsT[1:2, 1:2], tolerance = 1e-12)
  # every mirt<->lavaan cross entry is exactly zero
  expect_equal(Lr[c(1, 2), 3], c(0, 0), ignore_attr = TRUE)
  expect_equal(Lr[3, c(1, 2)], c(0, 0), ignore_attr = TRUE)
  expect_equal(Tr[c(1, 2), 3], c(0, 0), ignore_attr = TRUE)
  expect_equal(Tr[3, c(1, 2)], c(0, 0), ignore_attr = TRUE)
  # psi block-diagonal: mirt block on the diagonal, lavaan block, zero cross
  psi <- attr(comb, "psi")
  expect_equal(psi[3, c(1, 2)], c(0, 0), ignore_attr = TRUE)
  expect_equal(psi[c(1, 2), 3], c(0, 0), ignore_attr = TRUE)
})

# ===========================================================================
test_that("T3: two separate single-factor mirt fits give two 1 x 1 blocks", {
  skip_if_not_installed("mirt")
  # disambiguated (F1, G2): two 1 x 1 blocks, zero cross
  comb <- combine_fs(list(fmirt1fa, fmirt1fb))
  expect_equal(rownames(attr(comb, "fsL")[[1]]), c("fs_F1", "fs_G2"))
  expect_equal(attr(comb, "fsL")[[1]][1, 2], 0)
  expect_equal(attr(comb, "fsL")[[1]][2, 1], 0)
  expect_equal(attr(comb, "psi")[1, 2], 0)

  # same-named fits (both "F1") are rejected: duplicate latent names
  expect_error(combine_fs(list(fmirt1fa, fmirt1fa)), "duplicate latent name")
})

# ===========================================================================
test_that("T4: id-join unions the rows; absent rows get an all-NA block", {
  skip_if_not_installed("mirt")
  # mirt on 50 obs (ids 1:50), lavaan on all 301 (ids 1:301)
  fm <- fmirt1f_short
  fm$uid <- seq_len(nrow(fm))
  fl <- f_vis
  fl$uid <- seq_len(nrow(fl))
  comb <- combine_fs(list(fm, fl), id = "uid")
  expect_equal(nrow(comb), n_hs)  # union of 1:50 and 1:301
  # row 100 is absent from the mirt input
  Lr <- attr(comb, "fsL")[[100]]
  expect_true(all(is.na(Lr[1, 1])))   # mirt (F1) block all-NA
  expect_false(all(is.na(Lr[2, 2])))  # lavaan (vis) block present
  expect_true(is.na(comb$fs_F1[100]))
  expect_false(is.na(comb$fs_vis[100]))
  # row 10 is present in both
  expect_false(is.na(comb$fs_F1[10]))
  expect_false(is.na(comb$fs_vis[10]))
  expect_false(all(is.na(attr(comb, "fsL")[[10]][1, 1])))
})

# ===========================================================================
test_that("T5: id = NULL with mismatched row counts errors", {
  expect_error(
    combine_fs(list(f_vis, f_qua_short)),
    "same number of rows"
  )
})

# ===========================================================================
test_that("T6: duplicate latent names across inputs error", {
  f_vis2 <- get_fs(fit_vis)  # same model -> same latent name "vis"
  expect_error(combine_fs(list(f_vis, f_vis2)), "duplicate latent name")
})

# ===========================================================================
test_that("T7: scoring_matrix is always emitted, item-level where the input has one", {
  # lavaan-only: both inputs carry item-level scoring matrices -> p_total = 6
  sm <- attr(comb_vq, "scoring_matrix")
  expect_true(is.list(sm))
  expect_equal(length(sm), n_hs)
  expect_equal(dim(sm[[1]]), c(2L, 6L))  # 2 scores x 6 items
  # block-structured: each latent's row is populated only on its own items
  expect_false(any(is.na(sm[[1]][1, 1:3])))  # vis row, vis items
  expect_true(all(is.na(sm[[1]][1, 4:6])))   # vis row, qua items
  expect_true(all(is.na(sm[[1]][2, 1:3])))   # qua row, vis items
  expect_false(any(is.na(sm[[1]][2, 4:6])))  # qua row, qua items

  # mirt (no item model) leaves its rows NA; lavaan rows are populated
  skip_if_not_installed("mirt")
  comb <- combine_fs(list(fmirt1fa, f_vis))  # 301 rows each
  smm <- attr(comb, "scoring_matrix")
  expect_equal(dim(smm[[1]]), c(2L, 3L))  # F1 (NA) + vis (3 items)
  expect_true(all(is.na(smm[[1]][1, ])))   # mirt F1 row all-NA
  expect_false(any(is.na(smm[[1]][2, ])))  # lavaan vis row populated
})

# ===========================================================================
test_that("T8: the combined result feeds tspa() via the multi-factor derivation", {
  fit <- tspa(model = "qua ~ vis", data = comb_vq)
  co <- coef(fit)
  expect_true("qua~vis" %in% names(co))
  # the vis latent variance equals the block-diagonal psi diagonal (vis block)
  expect_equal(co[["vis~~vis"]], R2spa:::combine_fs_psi(f_vis)[1, 1],
               tolerance = 1e-4)
  # the derived measurement inputs come from the attributes (not explicit)
  expect_true(!is.null(attr(fit, "tspa_args")$fsT))
})

# ===========================================================================
test_that("T9: a multi-group (list-of-groups) result is rejected, not misread as inputs", {
  # mimic a multi-group get_fs() result: a (named) list of per-group data
  # frames carrying a group_col attribute on the list
  mg <- list(g1 = f_vis, g2 = f_qua)
  attr(mg, "group_col") <- "grp"
  expect_error(combine_fs(mg), "multi-group")
  # a plain list of single-group results (no group_col) still combines
  expect_identical(combine_fs(list(f_vis, f_qua)), comb_vq)
})

# ===========================================================================
test_that("T10: duplicate-latent message lists names cleanly (no deparse())", {
  f2a <- get_fs(cfa("a =~ x1 + x2 + x3\nb =~ x4 + x5 + x6", data = hs))
  f2b <- get_fs(cfa("a =~ x4 + x5 + x6\nb =~ x1 + x2 + x3", data = hs))
  err <- tryCatch(combine_fs(list(f2a, f2b)), error = identity)
  msg <- conditionMessage(err)
  expect_match(msg, "across inputs: a, b")
  expect_false(grepl("c(", msg, fixed = TRUE))
})

# ===========================================================================
test_that("T11: factor id columns align by label regardless of level order", {
  a <- f_vis; b <- f_qua
  lab <- as.character(seq_len(n_hs))
  a$fid <- factor(lab)                    # forward level order
  b$fid <- factor(lab, levels = rev(lab)) # reversed level order
  comb <- combine_fs(list(a, b), id = "fid")
  expect_equal(nrow(comb), n_hs)
  # aligned by label: row i carries observation i's scores from both inputs
  expect_equal(unname(comb$fs_vis), unname(as.numeric(a$fs_vis)))
  expect_equal(unname(comb$fs_qua), unname(as.numeric(b$fs_qua)))
})

# ===========================================================================
test_that("T12: NA or duplicated id values are rejected", {
  a <- f_vis; b <- f_qua
  a$uid <- seq_len(n_hs); b$uid <- seq_len(n_hs)
  a_na <- a; a_na$uid[5L] <- NA
  expect_error(combine_fs(list(a_na, b), id = "uid"),
               "NA values in the id column")
  a_dup <- a; a_dup$uid[6L] <- a_dup$uid[7L]
  expect_error(combine_fs(list(a_dup, b), id = "uid"),
               "duplicated values in the id column")
})
