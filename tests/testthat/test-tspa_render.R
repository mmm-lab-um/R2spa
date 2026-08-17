library(lavaan)

## Shared canonical setups ----------------------------------------------------

se_sf1 <- data.frame(ind60 = 0.1213615, dem60 = 0.6756472)
se_sf2 <- data.frame(visual = c(0.3391326, 0.3118280),
                     speed = c(0.2786875, 0.2740507))
se_sf3 <- data.frame(ind60 = 0.1213615, dem60 = 0.6756472,
                     dem65 = 0.5724405)

fs_s1 <- get_fs(PoliticalDemocracy, "ind60 =~ x1 + x2 + x3",
                format = "list")
fs_s2 <- get_fs(PoliticalDemocracy, "dem60 =~ y1 + y2 + y3 + y4",
                format = "list")
fs_sg <- cbind(fs_s1, fs_s2)

mod2g <- "visual =~ x1 + x2 + x3\nspeed =~ x7 + x8 + x9"
mod4g <- "visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6\nspeed =~ x7 + x8 + x9"
fs2g <- get_fs(HolzingerSwineford1939, model = mod2g, std.lv = TRUE,
               group = "school", format = "list")
fs4g <- get_fs(HolzingerSwineford1939, model = mod4g, std.lv = TRUE,
               group = "school", format = "list")
fs2g_long <- do.call(rbind, fs2g)
fs4g_long <- do.call(rbind, fs4g)

fsu_2f <- get_fs(PoliticalDemocracy,
                 "ind60 =~ x1 + x2 + x3\ndem60 =~ y1 + y2 + y3 + y4",
                 std.lv = TRUE)

## Pinned rendered format -------------------------------------------------------
## The reference builders below are verbatim copies of the pre-cutover
## string-append builders, kept in the tests as the format specification
## (the Phase 2 A/B gate is now frozen in place). lavaan-estimate drift
## affects both sides equally, so the comparisons are stable across lavaan
## versions.

ref_sf <- function(model, se) {
  if (nrow(se) != 0) {
    ev <- se^2
    var <- names(se)
    len <- length(se)
    group <- nrow(se)
    fs <- paste0("fs_", var)
    latent_var <- lapply(seq_len(len), function(x) {
      paste0(
        var[x], "=~ c(",
        paste0(rep(1, group), collapse = ", "), ") * ", fs[x], "\n"
      )
    })
    error_constraint <- lapply(seq_len(len), function(x) {
      paste0(
        fs[x], "~~ c(", paste(ev[[x]], collapse = ", "), ") * ", fs[x], "\n"
      )
    })
    paste0(
      c("# latent variables (indicated by factor scores)",
        paste(latent_var, collapse = ""),
        "# constrain the errors",
        paste(error_constraint, collapse = ""),
        "# structural model",
        model),
      collapse = "\n"
    )
  }
}

ref_mf <- function(model, data, fsT, fsL, fsb) {
  if (is.list(fsT) && length(fsT) > 1) {
    if (!is.list(fsL) || length(fsL) != length(fsT)) {
      stop("'fsL' must be a list of the same length as 'fsT' for a ",
           "multigroup model.")
    }
    ngroup <- length(fsT)
    fsL1 <- fsL[[1]]
    fsT_in <- !upper.tri(fsT[[1]])
  } else if (is.list(fsL) && length(fsL) > 1) {
    if (!is.list(fsT) || length(fsT) != length(fsL)) {
      stop("'fsT' must be a list of the same length as 'fsL' for a ",
           "multigroup model.")
    }
    ngroup <- length(fsL)
    fsL1 <- fsL[[1]]
    fsT_in <- !upper.tri(fsT[[1]])
  } else {
    ngroup <- 1
    fsL1 <- if (is.list(fsL)) fsL[[1]] else fsL
    fsT_in <- !upper.tri(if (is.list(fsT)) fsT[[1]] else fsT)
  }
  var <- colnames(fsL1)
  nvar <- length(var)
  fs <- rownames(fsL1)

  # latent variables
  loadings_mat <- matrix(unlist(fsL), ncol = ngroup)
  loadings <- apply(loadings_mat, 1, function(x) {
    paste0("c(", paste0(x, collapse = ", "), ") * ")
  }) |>
    paste0(fs)
  loadings_list <- split(loadings, factor(rep(var, each = nvar),
                                           levels = var))
  loadings_c <- lapply(loadings_list, function(x) {
    paste0(x, collapse = " + ")
  })
  latent_var_str <- paste("# latent variables (indicated by factor scores)\n",
                           var, "=~", loadings_c)
  # error variances
  ev_rhs <- fs[col(fsT_in)[fsT_in]]
  ev_lhs <- fs[row(fsT_in)[fsT_in]]
  errors_mat <- matrix(unlist(fsT), ncol = ngroup)[as.vector(fsT_in), ,
                                                    drop = FALSE]
  errors <- apply(errors_mat, 1, function(x) {
    paste0("c(", paste0(x, collapse = ", "), ")")
  })
  error_constraint_str <- paste0("# constrain the errors\n",
                                  ev_lhs, " ~~ ", errors, " * ", ev_rhs)
  if (!is.null(fsb)) {
    # intercepts
    intercepts_mat <- matrix(unlist(fsb), ncol = ngroup)
    intercepts <- split(intercepts_mat, rep(seq_len(nrow(intercepts_mat)), ngroup))
    intercept_constraint <- paste0("# constrain the intercepts\n",
                                    fs, " ~ ", intercepts, " * 1")
  } else {
    intercept_constraint <- ""
  }

  paste0(c(
    latent_var_str,
    error_constraint_str,
    intercept_constraint,
    "# structural model",
    model
  ),
  collapse = "\n")
}

test_that("SF renderer reproduces the pinned format character-for-character", {
  expect_identical(
    tspa_sf("dem60 ~ ind60", data.frame(x = 1), se_sf1),
    ref_sf("dem60 ~ ind60", se_sf1)
  )
  expect_identical(
    tspa_sf("visual ~ speed", data.frame(x = 1), se_sf2),
    ref_sf("visual ~ speed", se_sf2)
  )
  expect_identical(
    tspa_sf("dem60 ~ ind60\ndem65 ~ ind60 + dem60", data.frame(x = 1), se_sf3),
    ref_sf("dem60 ~ ind60\ndem65 ~ ind60 + dem60", se_sf3)
  )
  # verbatim user model with leading newline, comments, trailing blank lines
  m5 <- "\n# comment\ndem60 ~ ind60\n\n# trailing\n\n"
  expect_identical(tspa_sf(m5, data.frame(x = 1), se_sf3), ref_sf(m5, se_sf3))
})

test_that("MF renderer reproduces the pinned format character-for-character", {
  m <- "visual ~ speed"
  expect_identical(
    tspa_mf(m, fs2g, attr(fs2g, "fsT"), attr(fs2g, "fsL"), NULL),
    ref_mf(m, fs2g, attr(fs2g, "fsT"), attr(fs2g, "fsL"), NULL)
  )
  m4 <- "visual ~ speed\ntextual ~ visual + speed"
  expect_identical(
    tspa_mf(m4, fs4g, attr(fs4g, "fsT"), attr(fs4g, "fsL"), NULL),
    ref_mf(m4, fs4g, attr(fs4g, "fsT"), attr(fs4g, "fsL"), NULL)
  )
  m2 <- "dem60 ~ ind60"
  expect_identical(
    tspa_mf(m2, NULL, attr(fsu_2f, "fsT"), attr(fsu_2f, "fsL"), NULL),
    ref_mf(m2, NULL, attr(fsu_2f, "fsT"), attr(fsu_2f, "fsL"), NULL)
  )
})

## Schema construction ---------------------------------------------------------

test_that("SF schema: per-group rows, se^2 values, generated labels", {
  sch <- tspa_schema_sf("visual ~ speed", se_sf2)
  expect_equal(names(sch),
               c("lhs", "op", "rhs", "value", "free", "group", "label",
                 "kind"))
  # one verbatim user row
  expect_equal(sum(sch$kind == "user"), 1)
  expect_identical(sch$rhs[sch$kind == "user"], "visual ~ speed")
  # struct rows: one fixed loading per latent per group
  struct <- sch[sch$kind == "struct", , drop = FALSE]
  expect_equal(nrow(struct), 2L * 2L)
  expect_equal(struct$lhs, c("visual", "visual", "speed", "speed"))
  expect_equal(struct$value, rep(1, 4))
  expect_equal(struct$free, rep(0L, 4))
  expect_true(all(grepl("^__r2spa_ld[0-9]+__$", struct$label)))
  # error rows: se^2 per group
  errors <- sch[sch$kind == "error_var", , drop = FALSE]
  expect_equal(nrow(errors), 2L * 2L)
  g1 <- errors$value[errors$group == 1]
  g2 <- errors$value[errors$group == 2]
  expect_equal(g1, c(0.3391326^2, 0.2786875^2), tolerance = 1e-15)
  expect_equal(g2, c(0.3118280^2, 0.2740507^2), tolerance = 1e-15)
})

test_that("MF schema: per-group value routing is explicit", {
  sch <- tspa_schema_mf("visual ~ speed", attr(fs2g, "fsT"),
                        attr(fs2g, "fsL"), NULL)
  T_list <- attr(fs2g, "fsT")
  L_list <- attr(fs2g, "fsL")
  # struct rows carry the per-group loadings
  struct <- sch[sch$kind == "struct", , drop = FALSE]
  for (g in 1:2) {
    got <- struct$value[struct$group == g &
                          struct$lhs == "visual" &
                          struct$rhs == "fs_visual"]
    expect_equal(got, L_list[[g]][1, 1], tolerance = 1e-15)
    got2 <- struct$value[struct$group == g &
                           struct$lhs == "speed" &
                           struct$rhs == "fs_visual"]
    expect_equal(got2, L_list[[g]][1, 2], tolerance = 1e-15)
  }
  # error rows follow the lower triangle in column-major order (one row per
  # entry per group), with the per-group values of that entry
  errors <- sch[sch$kind %in% c("error_var", "error_cov"), , drop = FALSE]
  e1 <- errors[errors$group == 1, , drop = FALSE]
  e2 <- errors[errors$group == 2, , drop = FALSE]
  expect_equal(e1$lhs, c("fs_visual", "fs_speed", "fs_speed"))
  expect_equal(e1$rhs, c("fs_visual", "fs_visual", "fs_speed"))
  expect_equal(e1$kind, c("error_var", "error_cov", "error_var"))
  expect_equal(e1$value, c(T_list[[1]][1, 1], T_list[[1]][2, 1],
                           T_list[[1]][2, 2]), tolerance = 1e-15)
  expect_equal(e2$value, c(T_list[[2]][1, 1], T_list[[2]][2, 1],
                           T_list[[2]][2, 2]), tolerance = 1e-15)
  # every label lives in the generated namespace
  expect_true(all(grepl("^__r2spa_(ld|ev|int)[0-9]+__$",
                        sch$label[!is.na(sch$label)])))
  expect_false(any(grepl("^__r2spa_", sch$rhs)))
})

## tspa() end-to-end: rendered model pinned, fit usable, row order frozen ------

test_that("tspa() single-group fit carries the pinned model and fits cleanly", {
  fit <- tspa("dem60 ~ ind60", data = fs_sg, se_fs = c(ind60 = 0.1213615,
                                                      dem60 = 0.6756472))
  expect_identical(attr(fit, "tspaModel"), ref_sf("dem60 ~ ind60", se_sf1))
  expect_setequal(names(coef(fit)),
                  c("dem60~ind60", "ind60~~ind60", "dem60~~dem60"))
  v <- vcov(fit)
  expect_equal(dim(v), c(3, 3))
  expect_equal(v, t(v), tolerance = 1e-12)
  expect_true(all(is.finite(v)))
  # 7 parameter rows: 4 fixed + 3 free
  expect_equal(nrow(parameterestimates(fit)), 7L)
})

test_that("tspa() multigroup fit carries the pinned model and fits cleanly", {
  fit <- tspa("visual ~ speed", data = fs4g_long,
              se_fs = data.frame(visual = c(0.3391326, 0.3118280),
                                 speed = c(0.2786875, 0.2740507)),
              group = "school")
  # se_fs without fsT/fsL goes through the single-factor (SF) builder
  expect_identical(attr(fit, "tspaModel"),
                   ref_sf("visual ~ speed", se_sf2))
  v <- vcov(fit)
  expect_equal(v, t(v), tolerance = 1e-12)
  expect_true(all(is.finite(v)))
})

test_that("Multigroup parameter row order is the frozen pre-cutover order", {
  fit <- tspa("visual ~ speed", data = fs2g_long, fsT = attr(fs2g, "fsT"),
              fsL = attr(fs2g, "fsL"), group = "school")
  expect_identical(attr(fit, "tspaModel"),
                   ref_mf("visual ~ speed", fs2g, attr(fs2g, "fsT"),
                          attr(fs2g, "fsL"), NULL))
  pv <- lavaan::lavaanify(attr(fit, "tspaModel"), ngroup = 2)
  one_group <- data.frame(
    lhs   = c("visual", "visual", "speed", "speed",
              "fs_visual", "fs_visual", "fs_speed",
              "visual", "visual", "speed"),
    op    = c("=~", "=~", "=~", "=~",
              "~~", "~~", "~~",
              "~", "~~", "~~"),
    rhs   = c("fs_visual", "fs_speed", "fs_visual", "fs_speed",
              "fs_visual", "fs_speed", "fs_speed",
              "speed", "visual", "speed"),
    block = rep(1L, 10),
    group = 1,
    stringsAsFactors = FALSE
  )
  expected <- rbind(one_group, one_group)
  expected$group <- rep(1:2, each = nrow(one_group))
  expected$block <- expected$group
  expect_identical(as.data.frame(pv)[c("lhs", "op", "rhs", "block", "group")],
                   expected)
})

## tspa() contract --------------------------------------------------------------

test_that("tspa() attributes are intact and complete", {
  fit <- tspa("visual ~ speed", data = fs4g_long, fsT = attr(fs4g, "fsT"),
              fsL = attr(fs4g, "fsL"), group = "school")
  expect_true(inherits(fit, "lavaan"))
  expect_identical(attr(fit, "tspaModel"),
                   tspa_mf("visual ~ speed", fs4g, attr(fs4g, "fsT"),
                           attr(fs4g, "fsL"), NULL))
  expect_identical(attr(fit, "fsT"), attr(fs4g, "fsT"))
  expect_identical(attr(fit, "fsL"), attr(fs4g, "fsL"))
  expect_type(attr(fit, "tspa_call"), "language")
  # tspa() attaches no attributes beyond a plain lavaan fit of the same
  # model
  plain <- sem(attr(fit, "tspaModel"), data = fs4g_long, group = "school")
  expect_equal(
    sort(setdiff(names(attributes(fit)), names(attributes(plain)))),
    c("fsL", "fsT", "tspaModel", "tspa_call")
  )
  # single-group fit carries no fsT/fsL attributes
  fit_sg <- tspa("dem60 ~ ind60", data = fs_sg,
                 se_fs = c(ind60 = 0.1213615, dem60 = 0.6756472))
  expect_null(attr(fit_sg, "fsT"))
  expect_null(attr(fit_sg, "fsL"))
  expect_type(attr(fit_sg, "tspa_call"), "language")
  plain_sg <- sem(attr(fit_sg, "tspaModel"), data = fs_sg)
  expect_equal(
    sort(setdiff(names(attributes(fit_sg)), names(attributes(plain_sg)))),
    c("tspaModel", "tspa_call")
  )
})

## Interaction (product-score) auto-alias ---------------------------------------

test_that("Ambiguous product-score candidates are a clear error", {
  amb <- data.frame(fs_a = rnorm(50), fs_b = rnorm(50), fs_c = rnorm(50),
                    fs_ab = rnorm(50), fs_bc = rnorm(50), fs_y = rnorm(50))
  amb[["fs_a:fs_bc"]] <- amb[["fs_a"]] * amb[["fs_bc"]]
  amb[["fs_ab:fs_c"]] <- amb[["fs_ab"]] * amb[["fs_c"]]
  se_amb <- c(abc = 0.1, a = 0.2, ab = 0.5, bc = 0.6, c = 0.3, y = 0.1)
  expect_error(
    tspa("y ~ abc", data = amb, se_fs = se_amb),
    "Cannot determine which product-score column"
  )
})
