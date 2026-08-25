# compute_fs_prod() -- product factor-score indicators (double-mean-centered)
#
# Successor of the quarantined get_fs_int(): for each requested pair of
# distinct latents it appends fs_a:fs_b (the DMC product of the score
# columns), fs_a:fs_b_se (per-row; pattern-resolved under FIML), and
# fs_a:fs_b_ld (per-row implied loading gamma). The pure-matrix helpers
# fs_prod_se2() / fs_prod_gamma() carry the two formulas; the joint-normal
# derivation is in the roxygen @details of compute_fs_prod().
library(lavaan)

# ---------------------------------------------------------------------------
# Group 1: helper unit tests
# ---------------------------------------------------------------------------

test_that("fs_prod_se2() matches the joint-normal moment formula (independent MC)", {
  # Nonzero c = T[1, 2] and nonzero tau_ij: the terms the old
  # separate-single-factor formula dropped.
  L <- matrix(c(0.9, 0.15, 0.1, 0.85), 2)
  Tm <- matrix(c(0.05, 0.01, 0.01, 0.08), 2)
  psi <- matrix(c(1.2, 0.3, 0.3, 1.0), 2)
  n <- 2e5
  set.seed(4321)
  xi <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = psi)
  e <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Tm)
  # u = L_i (xi - alpha) e_j + L_j (xi - alpha) e_i + e_i e_j, alpha = 0
  # (xi %*% L[, k] is the n-vector of L[k, ] xi_row for every row)
  a1 <- as.numeric(xi %*% L[, 1L, drop = FALSE])
  a2 <- as.numeric(xi %*% L[, 2L, drop = FALSE])
  u <- a1 * e[, 2] + a2 * e[, 1] + e[, 1] * e[, 2]
  expect_equal(var(u), fs_prod_se2(L, Tm, psi, 1, 2), tolerance = 0.02)
  # The product is symmetric in the pair
  expect_equal(fs_prod_se2(L, Tm, psi, 2, 1), fs_prod_se2(L, Tm, psi, 1, 2))
})

test_that("fs_prod_se2() reduces to the legacy separate-models formula for diagonal L/T", {
  L <- diag(c(0.8, 0.9))
  Tm <- diag(c(0.05, 0.08))
  psi <- matrix(c(1.2, 0.3, 0.3, 1.0), 2)
  expect_equal(
    fs_prod_se2(L, Tm, psi, 1, 2),
    0.8^2 * 1.2 * 0.08 + 0.9^2 * 1.0 * 0.05 + 0.05 * 0.08,
    tolerance = 1e-12
  )
})

test_that("fs_prod_gamma() is L[i, i] L[j, j] + L[i, j] L[j, i] and symmetric in the pair", {
  L <- matrix(c(0.9, 0.15, 0.1, 0.85), 2) # non-symmetric
  expect_equal(fs_prod_gamma(L, 1, 2), L[1, 1] * L[2, 2] + L[1, 2] * L[2, 1])
  expect_equal(fs_prod_gamma(L, 1, 2), fs_prod_gamma(L, 2, 1))
})

test_that("parse_product_spec() accepts string, list, and 2-column matrix/data-frame forms", {
  lv <- c("x", "m", "z")
  pairs_xm_xz <- list(c("x", "m"), c("x", "z"))
  expect_equal(parse_product_spec("x:m + x:z", lv), pairs_xm_xz)
  expect_equal(parse_product_spec(list(c("x", "m"), c("x", "z")), lv),
               pairs_xm_xz)
  # column-major: rows (x, m) and (x, z)
  expect_equal(parse_product_spec(matrix(c("x", "x", "m", "z"), ncol = 2), lv),
               pairs_xm_xz)
  expect_equal(
    parse_product_spec(data.frame(a = c("x", "x"), b = c("m", "z")), lv),
    pairs_xm_xz
  )
})

test_that("parse_product_spec() rejects malformed, same-factor, unknown, and duplicate specs", {
  lv <- c("x", "m", "z")
  expect_error(parse_product_spec("x:x", lv), "same-factor")
  expect_error(parse_product_spec("a:b", lv), "Unknown latent name")
  expect_error(parse_product_spec("a:b", lv), "a, b")
  expect_error(parse_product_spec("x", lv), "not of the form 'a:b'")
  expect_error(parse_product_spec("a:b:c", lv), "not of the form 'a:b'")
  # reversed order is the same product: (x, m) + (m, x) is a duplicate
  expect_error(parse_product_spec("x:m + m:x", lv), "Duplicated")
  expect_error(parse_product_spec(matrix(c("x", "m", "z"), ncol = 3), lv),
               "2 columns")
  expect_error(parse_product_spec(NULL, lv), "character string")
})

# ---------------------------------------------------------------------------
# Group 2: population MC validation (no estimation error involved)
# ---------------------------------------------------------------------------

test_that("population MC: Var(product error) and gamma hold for both scoring methods", {
  n <- 2e5
  set.seed(2116)
  psi_pop <- matrix(c(1.0, 0.45, 0.45, 1.3), 2) # non-unit + correlated
  rownames(psi_pop) <- colnames(psi_pop) <- c("f1", "f2")
  # 6 items: x1-x3 on f1 (0.9/0.8/0.7), m1-m3 on f2 (0.85/0.75/0.65),
  # plus a 0.1 cross-loading of m1 (row 4) on f1.
  lambda_pop <- matrix(c(0.9, 0.0, 0.8, 0.0, 0.7, 0.0,
                         0.1, 0.85, 0.0, 0.75, 0.0, 0.65), nrow = 6)
  rownames(lambda_pop) <- c("x1", "x2", "x3", "m1", "m2", "m3")
  colnames(lambda_pop) <- c("f1", "f2")
  # One off-diagonal error covariance: x1 ~~ m1 = 0.2
  theta_pop <- diag(6)
  rownames(theta_pop) <- colnames(theta_pop) <- rownames(lambda_pop)
  theta_pop[1, 4] <- theta_pop[4, 1] <- 0.2
  xi <- MASS::mvrnorm(n, mu = rep(0, 2), Sigma = psi_pop)
  delta <- MASS::mvrnorm(n, mu = rep(0, 6), Sigma = theta_pop)
  y <- xi %*% t(lambda_pop) + delta
  rownames(y) <- paste0("obs", seq_len(n))

  for (method in c("regression", "Bartlett")) {
    fs_pop <- compute_fscore(as.data.frame(y), lambda_pop, theta_pop,
                             psi = psi_pop, method = method,
                             fs_matrices = TRUE)
    L_pop <- attr(fs_pop, "fsL")
    T_pop <- attr(fs_pop, "fsT")
    # Product measurement error: raw product minus E[fs1 fs2 | xi]
    # (latent means are zero, so the score intercepts b are zero).
    # xi %*% L_pop[, k] is the n-vector of L_pop[k, ] xi_row.
    u <- as.numeric(fs_pop[, 1] * fs_pop[, 2]) -
      as.numeric(xi %*% L_pop[, 1L, drop = FALSE]) *
      as.numeric(xi %*% L_pop[, 2L, drop = FALSE])
    expect_equal(var(u), fs_prod_se2(L_pop, T_pop, psi_pop, 1, 2),
                 tolerance = 0.02,
                 info = paste("Var(u) vs formula,", method))
    # gamma: the coefficient of xi1 xi2 in E[fs1 fs2 | xi]
    dd <- data.frame(
      p = fs_pop[, 1] * fs_pop[, 2],
      ix = xi[, 1] * xi[, 2],
      i1 = xi[, 1]^2,
      i2 = xi[, 2]^2
    )
    expect_equal(coef(lm(p ~ ix + i1 + i2, dd))[["ix"]],
                 fs_prod_gamma(L_pop, 1, 2), tolerance = 0.02,
                 info = paste("gamma via lm,", method))
  }
})

# ---------------------------------------------------------------------------
# Group 3: integration with get_fs() (exact, no MC)
# ---------------------------------------------------------------------------

pd_data <- PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")]
pd_model <- "ind60 =~ x1 + x2 + x3
             dem60 =~ y1 + y2 + y3 + y4"
fit_pd <- cfa(pd_model, data = pd_data) # no std.lv: latent variances != 1
fit_pd_std <- cfa(pd_model, data = pd_data, std.lv = TRUE)

# Hand computation of the two per-pair quantities from a get_fs() result's
# own attributes (deliberately NOT calling the fs_prod_* helpers): the
# formula se2 = tau1 s22 + tau2 s11 + s11 s22 + c^2 + 2 tau12 c and
# gamma = L[1, 1] L[2, 2] + L[1, 2] L[2, 1] for the pair (ind60, dem60).
prod_handcheck <- function(fs, se_col, ld_col) {
  L <- attr(fs, "fsL")[[1]]
  Tm <- attr(fs, "fsT")[[1]]
  P <- attr(fs, "psi")[[1]]
  tau1 <- as.numeric(L[1, ] %*% P %*% L[1, ])
  tau2 <- as.numeric(L[2, ] %*% P %*% L[2, ])
  tau12 <- as.numeric(L[1, ] %*% P %*% L[2, ])
  se2_hand <- tau1 * Tm[2, 2] + tau2 * Tm[1, 1] + Tm[1, 1] * Tm[2, 2] +
    Tm[1, 2]^2 + 2 * tau12 * Tm[1, 2]
  expect_equal(unique(se_col), sqrt(se2_hand))
  expect_equal(se_col, rep(sqrt(se2_hand), nrow(fs)))
  expect_equal(ld_col,
               rep(L[1, 1] * L[2, 2] + L[1, 2] * L[2, 1], nrow(fs)))
}

test_that("get_fs(product = ) appends the double-prefixed product columns with hand-checked values", {
  for (method in c("regression", "Bartlett")) {
    fs <- get_fs(fit_pd, product = "ind60:dem60", method = method)
    expect_true(all(
      c("fs_ind60:fs_dem60", "fs_ind60:fs_dem60_se", "fs_ind60:fs_dem60_ld") %in%
        names(fs)
    ))
    prod_handcheck(fs, fs[["fs_ind60:fs_dem60_se"]],
                   fs[["fs_ind60:fs_dem60_ld"]])
    # DMC identity: component-center on the sample score means, then center
    # the product on its own sample mean.
    x <- fs[["fs_ind60"]]
    yv <- fs[["fs_dem60"]]
    p <- (x - mean(x)) * (yv - mean(yv))
    p <- p - mean(p)
    expect_equal(fs[["fs_ind60:fs_dem60"]], p)
  }
})

test_that("get_fs() product entry points agree; the lavaan entry point equals compute_fs_prod()", {
  for (method in c("regression", "Bartlett")) {
    fs_lav <- get_fs(fit_pd, method = method, product = "ind60:dem60")
    expect_identical(fs_lav,
                     compute_fs_prod(get_fs(fit_pd, method = method),
                                     "ind60:dem60"))
    # The data-frame entry point fits the same cfa() and must agree exactly.
    fs_df <- get_fs(pd_data, model = pd_model, method = method,
                    product = "ind60:dem60")
    expect_identical(fs_df, fs_lav)
  }
})

test_that("get_fs(format = 'list', product = ) gives the same product columns", {
  fs_u <- get_fs(fit_pd, method = "regression", product = "ind60:dem60")
  fs_l <- get_fs(fit_pd, method = "regression", format = "list",
                 product = "ind60:dem60")
  new_cols <- c("fs_ind60:fs_dem60", "fs_ind60:fs_dem60_se",
                "fs_ind60:fs_dem60_ld")
  expect_equal(fs_l[, new_cols], fs_u[, new_cols])
})

test_that("compute_fs_prod() preserves every input attribute and leaves the existing columns unchanged", {
  fs_in <- get_fs(fit_pd, method = "regression")
  fs_out <- compute_fs_prod(fs_in, "ind60:dem60")
  for (an in setdiff(names(attributes(fs_in)), "names")) {
    expect_identical(attr(fs_out, an), attr(fs_in, an),
                     info = paste("attribute", an))
  }
  expect_equal(fs_out[, names(fs_in)], fs_in, ignore_attr = TRUE)
  expect_equal(names(fs_out),
               c(names(fs_in), "fs_ind60:fs_dem60",
                 "fs_ind60:fs_dem60_se", "fs_ind60:fs_dem60_ld"))
})

test_that("product accepts list and 2-column matrix forms; reversed pairs are duplicates", {
  fs_str <- get_fs(fit_pd, method = "regression", product = "ind60:dem60")
  fs_list <- get_fs(fit_pd, method = "regression",
                    product = list(c("ind60", "dem60")))
  expect_identical(fs_list, fs_str)
  fs_mat <- get_fs(fit_pd, method = "regression",
                   product = matrix(c("ind60", "dem60"), ncol = 2))
  expect_identical(fs_mat, fs_str)
  # (a, b) and (b, a) are the same product
  expect_error(get_fs(fit_pd, product = "ind60:dem60 + dem60:ind60"),
               "Duplicated")
})

test_that("std.lv = TRUE fit: product columns present and hand-check passes", {
  fs <- get_fs(fit_pd_std, product = "ind60:dem60")
  expect_true(all(
    c("fs_ind60:fs_dem60", "fs_ind60:fs_dem60_se",
      "fs_ind60:fs_dem60_ld") %in% names(fs)
  ))
  prod_handcheck(fs, fs[["fs_ind60:fs_dem60_se"]],
                 fs[["fs_ind60:fs_dem60_ld"]])
})

# ---------------------------------------------------------------------------
# Group 4: FIML (per-pattern) behavior
# ---------------------------------------------------------------------------

test_that("FIML: product SE/ld are pattern-resolved and NA on fully-missing rows", {
  set.seed(1334)
  dd <- pd_data
  dd$x1[!rbinom(nrow(dd), 1, 0.4)] <- NA
  # One fully-missing row (no scorable pattern), placed in the middle:
  # assemble_fs_blocks() sizes the output by max(case_idx), so a
  # fully-missing row at the END of the data is silently dropped by
  # get_fs() (reported as a package bug); a middle row is kept with
  # NA score/label.
  na_row <- as.data.frame(lapply(pd_data, function(col) NA))
  dd <- rbind(dd[1:30, , drop = FALSE], na_row, dd[31:nrow(dd), , drop = FALSE])
  fit_m <- suppressWarnings(cfa(pd_model, data = dd, missing = "fiml"))
  fs <- get_fs(fit_m, product = "ind60:dem60")

  # Per-pattern blocks: >= 2 patterns
  T_pats <- attr(fs, "fsT")[[1]]
  expect_type(T_pats, "list")
  expect_gte(length(T_pats), 2L)
  expect_true(all(vapply(T_pats, is.matrix, logical(1))))

  se_col <- fs[["fs_ind60:fs_dem60_se"]]
  ld_col <- fs[["fs_ind60:fs_dem60_ld"]]
  p_col <- fs[["fs_ind60:fs_dem60"]]
  expect_gte(length(unique(se_col[!is.na(se_col)])), 2L)

  labels <- attr(fs, "fs_pattern")[[1]]$label
  L_pats <- attr(fs, "fsL")[[1]]
  P <- attr(fs, "psi")[[1]]
  for (pl in names(T_pats)) {
    rows_p <- which(labels == pl)
    L_pl <- L_pats[[pl]]
    T_pl <- T_pats[[pl]]
    i <- match("ind60", colnames(L_pl))
    j <- match("dem60", colnames(L_pl))
    expect_equal(se_col[rows_p],
                 rep(sqrt(fs_prod_se2(L_pl, T_pl, P, i, j)),
                     length(rows_p)),
                 info = paste("SE for pattern", pl))
    expect_equal(ld_col[rows_p],
                 rep(fs_prod_gamma(L_pl, i, j), length(rows_p)),
                 info = paste("ld for pattern", pl))
  }

  # The fully-missing row: NA in all three product columns
  na_rows <- which(is.na(labels))
  expect_gte(length(na_rows), 1L)
  expect_true(all(is.na(p_col[na_rows])))
  expect_true(all(is.na(se_col[na_rows])))
  expect_true(all(is.na(ld_col[na_rows])))
})

# ---------------------------------------------------------------------------
# Group 5: end-to-end tspa() with auto-alias (adapted from the quarantined
# int_setup(), seed 2116, 4-factor correlated model, Bartlett, std.lv = TRUE)
# ---------------------------------------------------------------------------

prod_setup <- function(n = 500) {
  set.seed(2116)
  cov_xmz_ey <- matrix(c(1, 0.1, 0.15, 0,
                         0.1, 1, 0.12, 0,
                         0.15, 0.12, 1, 0,
                         0, 0, 0, 0.481351), nrow = 4)
  eta <- as.data.frame(MASS::mvrnorm(n, mu = rep(0, 4),
                                     Sigma = cov_xmz_ey))
  names(eta) <- c("x", "m", "z", "ey")
  lk <- list(x = c(0.9, 0.8, 0.7), m = c(0.85, 0.75, 0.65),
             z = c(0.8, 0.7, 0.6), y = c(0.75, 0.7, 0.65))
  obs <- setNames(lapply(c("x", "m", "z"), function(v0) {
    eta[[v0]] %*% t(lk[[v0]]) + rnorm(n * 3)
  }), c("x", "m", "z"))
  etay <- 0.3 * eta$x + 0.4 * eta$m + 0.2 * eta$z +
    0.1 * eta$x * eta$m + 0.15 * eta$x * eta$z +
    0.12 * eta$m * eta$z + eta$ey
  obs$y <- etay %*% t(lk$y) + rnorm(n * 3)
  df <- as.data.frame(do.call(cbind, obs[c("x", "m", "z", "y")]))
  names(df) <- c(paste0("x", 1:3), paste0("m", 1:3), paste0("z", 1:3),
                 paste0("y", 1:3))
  df
}

test_that("tspa() end-to-end on get_fs(product = ) output (auto-alias, manual-rename equivalence)", {
  df <- prod_setup()
  fs <- get_fs(df,
               model = "x =~ x1 + x2 + x3
                        m =~ m1 + m2 + m3
                        z =~ z1 + z2 + z3
                        y =~ y1 + y2 + y3",
               std.lv = TRUE, method = "Bartlett",
               product = "x:m + x:z + m:z")
  for (pr in c("fs_x:fs_m", "fs_x:fs_z", "fs_m:fs_z")) {
    expect_true(all(paste0(pr, c("", "_se", "_ld")) %in% names(fs)))
  }
  se <- c(y = fs[1, "fs_y_se"], x = fs[1, "fs_x_se"],
          m = fs[1, "fs_m_se"], z = fs[1, "fs_z_se"],
          xm = fs[1, "fs_x:fs_m_se"], xz = fs[1, "fs_x:fs_z_se"],
          mz = fs[1, "fs_m:fs_z_se"])
  expect_true(all(is.finite(se)))

  m <- "y ~ x + m + z + xm + xz + mz"
  fit1 <- tspa(m, data = fs, se_fs = se)
  expect_s4_class(fit1, "lavaan")
  expect_true(all(is.finite(unname(coef(fit1)))))
  expect_true(all(is.finite(vcov(fit1))))
  # The auto-aliased model names (not the `:`-separated data columns) are
  # used in the rendered model (the SF renderer writes `fs_v ~~ se^2 * fs_v`).
  m0 <- attr(fit1, "tspaModel")
  expect_match(m0, "fs_xm ~~ ", fixed = TRUE)
  expect_match(m0, "fs_xz ~~ ", fixed = TRUE)
  expect_match(m0, "fs_mz ~~ ", fixed = TRUE)
  expect_no_match(m0, "fs_x:fs_m", fixed = TRUE)

  # The manual pre-rename workaround gives a bit-identical fit
  fs2 <- fs
  fs2[["fs_xm"]] <- fs2[["fs_x:fs_m"]]
  fs2[["fs_xz"]] <- fs2[["fs_x:fs_z"]]
  fs2[["fs_mz"]] <- fs2[["fs_m:fs_z"]]
  fit2 <- tspa(m, data = fs2, se_fs = se)
  expect_identical(coef(fit1), coef(fit2))
  expect_identical(vcov(fit1), vcov(fit2))
})

# ---------------------------------------------------------------------------
# Group 6: error cases
# ---------------------------------------------------------------------------

test_that("compute_fs_prod() rejects multi-group, merMod, and provenance-less inputs", {
  # Multi-group via the data-frame entry point
  hs_mg <- HolzingerSwineford1939[c("x1", "x2", "x3", "x4", "x5", "x6",
                                    "school")]
  mg_model <- "visual =~ x1 + x2 + x3
               textual =~ x4 + x5 + x6"
  expect_error(
    get_fs(hs_mg, model = mg_model, group = "school",
           product = "visual:textual"),
    "multi-group"
  )
  fit_mg <- cfa(mg_model, data = hs_mg, group = "school")
  # Multi-group unified result
  expect_error(compute_fs_prod(get_fs(fit_mg), "visual:textual"),
               "single-group")
  # Multi-group list result
  expect_error(compute_fs_prod(get_fs(fit_mg, format = "list"),
                               "visual:textual"),
               "single-group")
  # merMod result (3-D fsT/fsL)
  fs_mer <- get_fs(lme4::lmer(Reaction ~ Days + (Days | Subject),
                              lme4::sleepstudy))
  expect_error(compute_fs_prod(fs_mer, "u0:u1"), "merMod")
  # A plain data frame is not a get_fs() result
  expect_error(compute_fs_prod(data.frame(x = 1:5), "a:b"),
               "not a usable get_fs")
  # Plain-matrix fsT/fsL attributes but no fs_pattern provenance
  fake <- data.frame(fs_a = rnorm(5), fs_b = rnorm(5))
  attr(fake, "fsT") <- diag(2)
  attr(fake, "fsL") <- diag(2)
  expect_error(compute_fs_prod(fake, "a:b"), "fs_pattern")
})

test_that("get_fs(product = ) validation errors: same-factor, unknown latent, local = TRUE, NULL", {
  expect_error(get_fs(fit_pd, product = "ind60:ind60"), "same-factor")
  expect_error(get_fs(fit_pd, product = "nope:dem60"), "nope")
  expect_error(get_fs(pd_data, model = pd_model, local = TRUE,
                      product = "ind60:dem60"),
               "local = TRUE")
  expect_error(compute_fs_prod(get_fs(fit_pd), NULL), "character string")
})
