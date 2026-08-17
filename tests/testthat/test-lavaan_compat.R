library(lavaan)

## Golden partable canary (PLAN 04 drift detector) --------------------------
##
## If a future lavaan release changes the partable layout, exactly these
## tests fail with a column diff; the fix is confined to R/lavaan_compat.R.
## The pinned table was generated live against lavaan 0.7-2 (2026-07-16).

canon_mod <- "ind60 =~ x1 + x2 + x3
              dem60 =~ y1 + y2 + y3 + y4
              dem60 ~ ind60"
canon_fit <- sem(canon_mod, data = PoliticalDemocracy, std.lv = TRUE)

mg_mod <- "visual =~ x1 + x2 + x3
           speed =~ x7 + x8 + x9
           visual ~ speed"
mg_fit <- sem(mg_mod, data = HolzingerSwineford1939, group = "school",
              std.lv = TRUE)

test_that("Golden partable canary: canonical layout is pinned", {
  tsp_layout_reset()
  got <- tsp_partable_read(canon_fit)
  expected <- data.frame(
    lhs = c("ind60", "ind60", "ind60", "dem60", "dem60", "dem60", "dem60",
            "dem60", "x1", "x2", "x3", "y1", "y2", "y3", "y4", "ind60",
            "dem60"),
    op = c("=~", "=~", "=~", "=~", "=~", "=~", "=~", "~", "~~", "~~", "~~",
           "~~", "~~", "~~", "~~", "~~", "~~"),
    rhs = c("x1", "x2", "x3", "y1", "y2", "y3", "y4", "ind60", "x1", "x2",
            "x3", "y1", "y2", "y3", "y4", "ind60", "dem60"),
    value = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1, 1),
    free = c(rep(1L, 15), 0L, 0L),
    group = rep(1L, 17),
    block = rep(1L, 17),
    label = rep("", 17),
    user = c(rep(1L, 8), rep(0L, 9)),
    ustart = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1, 1),
    exo = rep(0L, 17),
    stringsAsFactors = FALSE
  )
  expect_identical(got, expected)
})

test_that("tsp_partable_read matches the lavInspect list view", {
  li <- as.data.frame(lavInspect(canon_fit, what = "list"))
  got <- tsp_partable_read(canon_fit)
  cols <- c("lhs", "op", "rhs", "group", "block", "label", "user", "ustart",
            "exo")
  expect_identical(got[, cols], li[, cols])
  # 0.7-2: the raw view's `free` column is the 1-based position in the free
  # estimate vector (0 = fixed); the canonical df carries the 0/1 flag.
  expect_identical(got$free, as.integer(li$free != 0))
  fixed <- got$free == 0
  expect_identical(got$value[fixed], li$start[fixed])
  expect_true(all(is.na(got$value[!fixed])))
})

test_that("tsp_partable_read works on a multigroup fit", {
  got <- tsp_partable_read(mg_fit)
  li <- as.data.frame(lavInspect(mg_fit, what = "list"))
  cols <- c("lhs", "op", "rhs", "group", "block", "label", "user", "ustart",
            "exo")
  expect_equal(nrow(got), nrow(li))
  expect_identical(got[, cols], li[, cols])
  expect_identical(got$free, as.integer(li$free != 0))
  fixed <- got$free == 0
  expect_identical(got$value[fixed], li$start[fixed])
  expect_true(all(is.na(got$value[!fixed])))
})

## Wrapper A/B ---------------------------------------------------------------
## Each wrapper must return exactly what the lavaan calls used to return at
## the migrated call sites (Phase 1 = no behavior change).

test_that("tsp_model_matrices matches lavTech(what = 'est')", {
  expect_identical(tsp_model_matrices(canon_fit),
                   lavTech(canon_fit, what = "est"))
  expect_identical(tsp_model_matrices(mg_fit),
                   lavTech(mg_fit, what = "est"))
})

test_that("tsp_free_matrices matches lavTech(what = 'free')", {
  expect_identical(tsp_free_matrices(canon_fit),
                   lavTech(canon_fit, what = "free"))
  expect_identical(tsp_free_matrices(mg_fit),
                   lavTech(mg_fit, what = "free"))
})

test_that("tsp_partable_mats matches lavTech(what = 'partable', list.by.group = TRUE)", {
  expect_identical(tsp_partable_mats(canon_fit),
                   lavTech(canon_fit, what = "partable", list.by.group = TRUE))
  expect_identical(tsp_partable_mats(mg_fit),
                   lavTech(mg_fit, what = "partable", list.by.group = TRUE))
})

test_that("tsp_nobs / tsp_ngroups / tsp_norig match slots and lavInspect", {
  for (fit in list(canon_fit, mg_fit)) {
    expect_identical(tsp_nobs(fit), unlist(fit@Data@nobs))
    expect_identical(tsp_nobs(fit), lavInspect(fit, what = "nobs"))
    expect_identical(tsp_ngroups(fit), fit@Data@ngroups)
    expect_identical(tsp_ngroups(fit), lavInspect(fit, what = "ngroups"))
    expect_identical(tsp_norig(fit), unlist(fit@Data@norig))
    expect_identical(tsp_norig(fit), lavInspect(fit, what = "norig"))
  }
})

test_that("grand_standardized_solution output is unchanged by the wrapper", {
  # Single-group: must match lavaan::standardizedSolution() as before.
  got <- suppressMessages(grandStandardizedSolution(canon_fit))
  lav <- subset(standardizedSolution(canon_fit), op == "~")
  expect_equal(got$est.std, lav$est.std)
  expect_equal(got$se, lav$se, tolerance = 1e-7)
  # The returned frame keeps its historical columns (incl. `exo`).
  expect_true(all(c("lhs", "op", "rhs", "exo", "group", "block", "label",
                    "est.std", "se") %in% names(got)))
  # Multigroup path runs through the wrappers and stays finite; the MG
  # hand-calculation A/B lives in test-grandStandardizedSolution.R.
  mg <- grandStandardizedSolution(mg_fit)
  expect_true(all(is.finite(mg$est.std)))
  expect_true(all(is.finite(mg$se)))
})

## Layout probing ------------------------------------------------------------

test_that("tsp_resolve_layout errors loudly on an unknown layout", {
  ver <- as.character(packageVersion("lavaan"))
  pt <- tsp_partable_raw(canon_fit)
  expect_error(
    tsp_resolve_layout(pt[!names(pt) %in% c("lhs", "op", "rhs")]),
    "partable layout not supported.*R2spa is tested up to lavaan"
  )
  # missing both fixed-indicator and fixed-value columns
  expect_error(
    tsp_resolve_layout(pt[!names(pt) %in% c("free", "fix", "start")]),
    "partable layout not supported"
  )
  err <- tryCatch(tsp_resolve_layout(pt[!names(pt) %in% c("free", "fix",
                                                          "start")]),
                  error = conditionMessage)
  expect_match(err, ver, fixed = TRUE)
  expect_match(err, "0.7-2", fixed = TRUE)
})

test_that("tsp_layout is memoized per lavaan version", {
  tsp_layout_reset()
  l1 <- tsp_layout(canon_fit)
  expect_type(l1, "list")
  # once resolved, the probe is cached: even a non-fit argument must not
  # re-probe (and must not error)
  expect_identical(l1, tsp_layout("not a lavaan fit"))
  tsp_layout_reset()
  expect_error(tsp_layout("not a lavaan fit"))
})
