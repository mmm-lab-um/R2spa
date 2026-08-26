# tspa(product = TRUE) -- opt-in product-score auto-compute
#
# With product = TRUE, tspa() computes the double-mean-centered product
# indicators for every model latent that names the product of two of the
# model's factor scores -- by concatenation (`xm` for `x` and `m`) or in
# lavaan's interaction syntax (`x:m`, rendered under the concatenated name)
# -- and wires them into the stage-2 measurement model: the single-factor
# (score-scale) path
# joins the (pooled) product SE into se_fs (loading 1, like every other
# single-factor latent); the multi-factor path adds a fixed loading
# gamma and fixed error variance se_P^2 (from the (pooled) fsL/fsT and
# the psi attribute). v1: single-group only; corrected_se = TRUE
# rejected; attribute-less data without the product columns rejected.
library(lavaan)

# ---------------------------------------------------------------------------
# Setup: 4-factor correlated model (the Phase-1 end-to-end scenario)
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

prod_model <- "x =~ x1 + x2 + x3
               m =~ m1 + m2 + m3
               z =~ z1 + z2 + z3
               y =~ y1 + y2 + y3"

# ---------------------------------------------------------------------------
# Group 1: tspa_product_latents() detection
# ---------------------------------------------------------------------------

test_that("tspa_product_latents() detects concatenated product latents", {
  prods <- tspa_product_latents(
    "y ~ x + m + z + xm + xz + mz",
    c("y", "fs_x", "x1"), c("y", "x", "m", "z")
  )
  expect_equal(prods$v, c("xm", "xz", "mz"))
  expect_equal(sort(prods$a), c("m", "x", "x"))
  expect_equal(sort(prods$b), c("m", "z", "z"))
})

test_that("tspa_product_latents() is orientation-agnostic and ignores non-product tokens", {
  # reversed concatenation: `mx` is the product of (m, x)
  prods <- tspa_product_latents("y ~ x + m + mx",
                                c("y", "fs_x"), c("y", "x", "m"))
  expect_equal(nrow(prods), 1L)
  expect_equal(prods$v, "mx")
  expect_equal(sort(unname(prods$a)), "m")
  expect_equal(sort(unname(prods$b)), "x")
  # no product latents -> NULL
  expect_null(tspa_product_latents("y ~ x + m", c("y"), c("y", "x", "m")))
  # an unknown (non-concatenation) model variable is ignored, not an error
  expect_null(tspa_product_latents("y ~ x + m + w", c("y"), c("y", "x", "m")))
  # numeric values, loading labels, and the c() wrapper are not latents
  expect_null(tspa_product_latents(
    "y ~ 1 + c(b1, b2) * x + 0.5 * m", c("y"), c("y", "x", "m")))
})

test_that("tspa_product_latents() strips comments and rejects ambiguous candidates", {
  prods <- tspa_product_latents(
    paste0("# the product xm\n", "y ~ x + m + xm"),
    c("y"), c("y", "x", "m")
  )
  expect_equal(prods$v, "xm")
  # `abc` = (a, bc) and (ab, c): ambiguous
  expect_error(
    tspa_product_latents("y ~ abc", c("y"), c("y", "a", "bc", "ab", "c")),
    "Cannot determine which factor-score pair"
  )
})

test_that("tspa_product_latents() detects a:b interaction-syntax tokens", {
  prods <- tspa_product_latents(
    "y ~ x + m + z + x:m + x:z + m:z",
    c("y", "fs_x", "x1"), c("y", "x", "m", "z")
  )
  expect_equal(prods$tok, c("x:m", "x:z", "m:z"))
  expect_equal(prods$v, c("xm", "xz", "mz"))
  expect_equal(sort(prods$a), c("m", "x", "x"))
  expect_equal(sort(prods$b), c("m", "z", "z"))
})

test_that("tspa_product_latents() leaves non-product a:b tokens alone and rejects conflicts", {
  # `x:g` with g not a factor score: a genuine lavaan interaction,
  # passed through (not claimed)
  expect_null(tspa_product_latents("y ~ x + m + x:g",
                                   c("y"), c("y", "x", "m")))
  # `x:x` is a squared term in lavaan, not a same-factor product
  expect_null(tspa_product_latents("y ~ x + m + x:x",
                                   c("y"), c("y", "x", "m")))
  # a statement label (`b1:`) is not a product token
  expect_null(tspa_product_latents("b1: y ~ x + m",
                                   c("y"), c("y", "x", "m")))
  # the same pair named twice (`x:m` and `xm`) is an error
  expect_error(
    tspa_product_latents("y ~ x + m + x:m + xm", c("y"), c("y", "x", "m")),
    "names the same factor-score pair more than once"
  )
  # a render name colliding with another model variable is an error
  # (`xm` is a regular latent here: its score column is present)
  expect_error(
    tspa_product_latents("y ~ x + m + xm + x:m", c("y", "fs_xm"),
                         c("y", "x", "m", "xm")),
    "collides with another variable"
  )
})

test_that("tspa_rewrite_product_toks() rewrites whole model variables only", {
  prods <- tspa_product_latents(
    "y ~ x + m + z + x:m + x:z", c("y"), c("y", "x", "m", "z"))
  expect_identical(
    tspa_rewrite_product_toks("y ~ x + m + z + x:m + x:z", prods),
    "y ~ x + m + z + xm + xz"
  )
  # a longer token containing the product token is untouched
  prods2 <- tspa_product_latents(
    "y ~ xx:m + x:m", c("y"), c("y", "x", "xx", "m"))
  expect_identical(
    tspa_rewrite_product_toks("y ~ xx:m + x:m", prods2),
    "y ~ xxm + xm"
  )
  # concatenated tokens (form b) need no rewrite
  prods3 <- tspa_product_latents("y ~ x + m + xm", c("y"), c("y", "x", "m"))
  expect_identical(
    tspa_rewrite_product_toks("y ~ x + m + xm", prods3),
    "y ~ x + m + xm"
  )
})

# ---------------------------------------------------------------------------
# Group 2: single-factor path -- auto-compute equals the manual workflow
# ---------------------------------------------------------------------------

test_that("sf: product = TRUE auto-compute equals the manual se_fs workflow", {
  df <- prod_setup()
  for (method in c("regression", "Bartlett")) {
    fs_plain <- get_fs(df, model = prod_model, std.lv = TRUE,
                       method = method)
    fs_prod <- get_fs(df, model = prod_model, std.lv = TRUE,
                      method = method,
                      product = "x:m + x:z + m:z")
    se_reg <- c(y = fs_plain[1, "fs_y_se"], x = fs_plain[1, "fs_x_se"],
                m = fs_plain[1, "fs_m_se"], z = fs_plain[1, "fs_z_se"])
    se_man <- c(se_reg, xm = fs_prod[1, "fs_x:fs_m_se"],
                xz = fs_prod[1, "fs_x:fs_z_se"],
                mz = fs_prod[1, "fs_m:fs_z_se"])
    m <- "y ~ x + m + z + xm + xz + mz"
    fit_auto <- suppressWarnings(
      tspa(m, data = fs_plain, se_fs = se_reg, product = TRUE)
    )
    fit_man <- suppressWarnings(tspa(m, data = fs_prod, se_fs = se_man))
    expect_identical(attr(fit_auto, "tspaModel"),
                     attr(fit_man, "tspaModel"),
                     info = paste("model string,", method))
    expect_equal(coef(fit_auto), coef(fit_man),
                 info = paste("coefs,", method))
    expect_equal(vcov(fit_auto), vcov(fit_man),
                 info = paste("vcov,", method))
  }
})

test_that("sf: pre-existing product columns (either orientation) are used as-is", {
  df <- prod_setup()
  # the user computed the pair in the opposite orientation
  fs_rev <- get_fs(df, model = prod_model, std.lv = TRUE,
                   method = "Bartlett", product = "m:x")
  se_reg <- c(y = fs_rev[1, "fs_y_se"], x = fs_rev[1, "fs_x_se"],
              m = fs_rev[1, "fs_m_se"], z = fs_rev[1, "fs_z_se"])
  se_man <- c(se_reg, xm = fs_rev[1, "fs_m:fs_x_se"])
  m <- "y ~ x + m + z + xm"
  fit_auto <- tspa(m, data = fs_rev, se_fs = se_reg, product = TRUE)
  fit_man <- tspa(m, data = fs_rev, se_fs = se_man)
  expect_equal(coef(fit_auto), coef(fit_man))
  expect_equal(vcov(fit_auto), vcov(fit_man))
  # no second (duplicated) product column was created: only the
  # original reversed-orientation triple is present
  dnames <- names(attr(fit_auto, "tspa_args")$data)
  expect_equal(sum(grepl(":", dnames)), 3L)
  expect_false("fs_x:fs_m" %in% dnames)
})

test_that("sf: an explicit se_fs entry for the product latent wins (no override)", {
  df <- prod_setup()
  fs_plain <- get_fs(df, model = prod_model, std.lv = TRUE,
                     method = "Bartlett")
  se_reg <- c(y = fs_plain[1, "fs_y_se"], x = fs_plain[1, "fs_x_se"],
              m = fs_plain[1, "fs_m_se"], z = fs_plain[1, "fs_z_se"])
  m <- "y ~ x + m + z + xm"
  fit_a <- tspa(m, data = fs_plain, se_fs = se_reg, product = TRUE)
  se_a <- attr(fit_a, "tspa_args")$se_fs[["xm"]]
  # an explicit (deliberately different) value is kept verbatim
  fit_b <- tspa(m, data = fs_plain,
                se_fs = c(se_reg, xm = 0.11), product = TRUE)
  expect_equal(attr(fit_b, "tspa_args")$se_fs[["xm"]], 0.11)
  expect_false(isTRUE(all.equal(se_a, 0.11)))
})

# ---------------------------------------------------------------------------
# Group 3: multi-factor path -- product rows with gamma and se_P^2
# ---------------------------------------------------------------------------

test_that("mf: product latent gets the implied loading gamma and se_P^2", {
  df <- prod_setup()
  for (method in c("regression", "Bartlett")) {
    fs <- get_fs(df, model = prod_model, std.lv = TRUE, method = method)
    fit <- suppressWarnings(
      tspa("y ~ x + m + z + xm", data = fs, product = TRUE)
    )
    L <- attr(fs, "fsL")[[1]]
    Tm <- attr(fs, "fsT")[[1]]
    P <- attr(fs, "psi")[[1]]
    i <- match("x", colnames(L))
    j <- match("m", colnames(L))
    g <- fs_prod_gamma(L, i, j)
    se2 <- fs_prod_se2(L, Tm, P, i, j)
    m0 <- attr(fit, "tspaModel")
    expect_match(m0, paste0("xm =~ ", as.character(g), " * fs_xm"),
                 fixed = TRUE,
                 info = paste("loading,", method))
    expect_match(m0,
                 paste0("fs_xm ~~ ", as.character(se2), " * fs_xm"),
                 fixed = TRUE,
                 info = paste("error variance,", method))
    expect_true(all(is.finite(unname(coef(fit)))),
                info = paste("coefs finite,", method))
    expect_true(all(is.finite(vcov(fit))),
                info = paste("vcov finite,", method))
  }
  # regression scoring: gamma differs from 1 (the Bartlett value)
  fs_reg <- get_fs(df, model = prod_model, std.lv = TRUE,
                   method = "regression")
  L <- attr(fs_reg, "fsL")[[1]]
  g_reg <- fs_prod_gamma(L, 1L, 2L)
  expect_false(isTRUE(all.equal(g_reg, 1)))
})

test_that("mf: product = FALSE (default) is a no-op (the product latent is not created)", {
  df <- prod_setup()
  fs <- get_fs(df, model = prod_model, std.lv = TRUE, method = "Bartlett")
  # without product the product latent is not in the data: lavaan rejects
  expect_error(
    suppressWarnings(tspa("y ~ x + m + z + xm", data = fs))
  )
  # a model without the product latent fits, and no product column is
  # created in the working data
  fit_ok <- suppressWarnings(
    tspa("y ~ x + m + z", data = fs, product = FALSE)
  )
  expect_false(any(grepl("fs_xm", names(attr(fit_ok, "tspa_args")$data))))
  expect_no_match(attr(fit_ok, "tspaModel"), "xm", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# Group 4: FIML (per-pattern) behavior
# ---------------------------------------------------------------------------

test_that("sf FIML: product SE is pooled from the per-row column (reduce = mean)", {
  set.seed(1334)
  pd <- PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")]
  pd_model <- "ind60 =~ x1 + x2 + x3
               dem60 =~ y1 + y2 + y3 + y4"
  dd <- pd
  dd$x1[!rbinom(nrow(dd), 1, 0.4)] <- NA
  fit_m <- suppressWarnings(cfa(pd_model, data = dd, missing = "fiml"))
  fs_fiml <- get_fs(fit_m, method = "regression")
  se_reg <- c(ind60 = fs_fiml[1, "fs_ind60_se"],
              dem60 = fs_fiml[1, "fs_dem60_se"])
  fit <- suppressWarnings(
    tspa("dem60 ~ ind60 + ind60dem60", data = fs_fiml, se_fs = se_reg,
         product = TRUE)
  )
  # per-row product SE column was created (one value per pattern)
  dat <- attr(fit, "tspa_args")$data
  se_col <- paste0("fs_", sort(c("ind60", "dem60")), collapse = ":fs_")
  # canonical (sorted) pair order: fs_dem60:fs_ind60_se
  expect_true("fs_dem60:fs_ind60_se" %in% names(dat))
  expect_gte(length(unique(dat[["fs_dem60:fs_ind60_se"]], na.rm = TRUE)),
             2L)
  # the model used the per-group mean of the per-row column
  se_used <- attr(fit, "tspa_args")$se_fs[["ind60dem60"]]
  expect_equal(as.numeric(se_used),
               mean(dat[["fs_dem60:fs_ind60_se"]], na.rm = TRUE))
  expect_true(all(is.finite(unname(coef(fit)))))
})

test_that("mf FIML: product rows use the pooled matrices; fit is finite", {
  set.seed(1334)
  pd <- PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")]
  pd_model <- "ind60 =~ x1 + x2 + x3
               dem60 =~ y1 + y2 + y3 + y4"
  dd <- pd
  dd$x1[!rbinom(nrow(dd), 1, 0.4)] <- NA
  fit_m <- suppressWarnings(cfa(pd_model, data = dd, missing = "fiml"))
  fs_fiml <- get_fs(fit_m, method = "regression")
  fit <- suppressWarnings(
    tspa("dem60 ~ ind60 + ind60dem60", data = fs_fiml, product = TRUE)
  )
  # gamma/se_P^2 are finite and the loading is the pooled-matrix value
  m0 <- attr(fit, "tspaModel")
  ld_str <- sub("ind60dem60 =~ ", "",
                regmatches(m0, regexpr("ind60dem60 =~ [^ *]+", m0)))
  expect_true(is.finite(as.numeric(ld_str)))
  expect_true(all(is.finite(unname(coef(fit)))))
})

# ---------------------------------------------------------------------------
# Group 5: replay via tspa_args
# ---------------------------------------------------------------------------

test_that("product fits replay identically via tspa_args (idempotent)", {
  df <- prod_setup()
  fs_plain <- get_fs(df, model = prod_model, std.lv = TRUE,
                     method = "Bartlett")
  se_reg <- c(y = fs_plain[1, "fs_y_se"], x = fs_plain[1, "fs_x_se"],
              m = fs_plain[1, "fs_m_se"], z = fs_plain[1, "fs_z_se"])
  m <- "y ~ x + m + z + xm + xz + mz"
  fit_sf <- tspa(m, data = fs_plain, se_fs = se_reg, product = TRUE)
  args_sf <- attr(fit_sf, "tspa_args")
  expect_true("product" %in% names(args_sf))
  expect_true(isTRUE(args_sf$product))
  fit_sf_re <- do.call(tspa, args_sf)
  expect_identical(coef(fit_sf_re), coef(fit_sf))
  expect_identical(vcov(fit_sf_re), vcov(fit_sf))

  fit_mf <- suppressWarnings(
    tspa("y ~ x + m + z + xm", data = fs_plain, product = TRUE)
  )
  args_mf <- attr(fit_mf, "tspa_args")
  fit_mf_re <- do.call(tspa, args_mf)
  expect_identical(coef(fit_mf_re), coef(fit_mf))
  expect_identical(vcov(fit_mf_re), vcov(fit_mf))
})

# ---------------------------------------------------------------------------
# Group 6: rejections
# ---------------------------------------------------------------------------

test_that("product = TRUE rejects multigroup, corrected_se, and attribute-less data", {
  # multigroup (mf path)
  hs <- HolzingerSwineford1939[c("x1", "x2", "x3", "x4", "x5", "x6",
                                 "school")]
  mg_model <- "visual =~ x1 + x2 + x3
               textual =~ x4 + x5 + x6"
  fit_mg <- cfa(mg_model, data = hs, group = "school")
  fs_mg <- get_fs(fit_mg, std.lv = TRUE)
  expect_error(
    tspa("visual ~ textual + visualtextual", data = fs_mg, product = TRUE,
         group = "school"),
    "multigroup"
  )

  # corrected_se + product
  pd <- PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")]
  pd_model <- "ind60 =~ x1 + x2 + x3
               dem60 =~ y1 + y2 + y3 + y4"
  fit_pd <- cfa(pd_model, data = pd)
  fs <- get_fs(fit_pd)
  expect_error(
    tspa("dem60 ~ ind60 + ind60dem60", data = fs, product = TRUE,
         corrected_se = TRUE,
         vfsLT = attr(get_fs(fit_pd, vfsLT = TRUE), "vfsLT"),
         fsT = attr(fs, "fsT"), fsL = attr(fs, "fsL")),
    "not supported with 'product = TRUE'"
  )

  # attribute-less (cbind'd) data without the product columns
  cb <- cbind(fs[, c("fs_ind60", "fs_ind60_se", "fs_dem60", "fs_dem60_se")])
  expect_error(
    tspa("dem60 ~ ind60 + ind60dem60", data = cb, product = TRUE),
    "lacks the stage-1 attributes"
  )
  # ... but a cbind'd frame that already carries the product columns works
  fs_prod <- get_fs(fit_pd, product = "ind60:dem60")
  cb2 <- cbind(
    fs[, c("fs_ind60", "fs_ind60_se", "fs_dem60", "fs_dem60_se")],
    fs_prod[, c("fs_ind60:fs_dem60", "fs_ind60:fs_dem60_se")]
  )
  fit_cb <- suppressWarnings(
    tspa("dem60 ~ ind60 + ind60dem60", data = cb2, product = TRUE)
  )
  expect_true(all(is.finite(unname(coef(fit_cb)))))

  # local-mode (complete data) results carry block-diagonal fsL/fsT/psi,
  # so auto-compute works and reduces to the separate-models formula
  # (gamma = lambda_x * lambda_m, c = tau_ab = 0)
  df <- prod_setup()
  fs_local <- get_fs(df, model = prod_model, std.lv = TRUE, local = TRUE)
  fit_local <- suppressWarnings(
    tspa("y ~ x + m + z + xm", data = fs_local, product = TRUE)
  )
  L <- attr(fs_local, "fsL")
  if (is.list(L)) L <- L[[1]]
  Tm <- attr(fs_local, "fsT")
  if (is.list(Tm)) Tm <- Tm[[1]]
  P <- attr(fs_local, "psi")
  if (is.list(P) && length(P) == 1L) P <- P[[1]]
  g_local <- L[1L, 1L] * L[2L, 2L]
  expect_match(attr(fit_local, "tspaModel"),
               paste0("xm =~ ", as.character(g_local), " * fs_xm"),
               fixed = TRUE)
  expect_true(all(is.finite(unname(coef(fit_local)))))

  # product must be a single logical
  expect_error(
    tspa("dem60 ~ ind60", data = fs, product = "ind60:dem60"),
    "must be a single TRUE/FALSE value"
  )
})

# ---------------------------------------------------------------------------
# Group 7: lavaan interaction syntax (`a:b` product latents)
# ---------------------------------------------------------------------------

test_that("sf: a:b syntax fits identically to the concatenated syntax", {
  df <- prod_setup()
  for (method in c("regression", "Bartlett")) {
    fs_plain <- get_fs(df, model = prod_model, std.lv = TRUE,
                       method = method)
    se_reg <- c(y = fs_plain[1, "fs_y_se"], x = fs_plain[1, "fs_x_se"],
                m = fs_plain[1, "fs_m_se"], z = fs_plain[1, "fs_z_se"])
    fit_colon <- suppressWarnings(
      tspa("y ~ x + m + z + x:m + x:z + m:z", data = fs_plain,
           se_fs = se_reg, product = TRUE)
    )
    fit_concat <- suppressWarnings(
      tspa("y ~ x + m + z + xm + xz + mz", data = fs_plain,
           se_fs = se_reg, product = TRUE)
    )
    expect_identical(attr(fit_colon, "tspaModel"),
                     attr(fit_concat, "tspaModel"),
                     info = paste("model string,", method))
    expect_equal(coef(fit_colon), coef(fit_concat),
                 info = paste("coefs,", method))
    expect_equal(vcov(fit_colon), vcov(fit_concat),
                 info = paste("vcov,", method))
  }
})

test_that("sf: an explicit product se_fs may be keyed by the a:b token", {
  df <- prod_setup()
  fs_prod <- get_fs(df, model = prod_model, std.lv = TRUE,
                    method = "Bartlett", product = "x:m")
  se_reg <- c(y = fs_prod[1, "fs_y_se"], x = fs_prod[1, "fs_x_se"],
              m = fs_prod[1, "fs_m_se"], z = fs_prod[1, "fs_z_se"])
  # keyed by the model token (as.data.frame() stores it as `x.m`)
  fit_tok <- suppressWarnings(tspa(
    "y ~ x + m + z + x:m", data = fs_prod,
    se_fs = c(se_reg, "x:m" = 0.11), product = TRUE))
  # keyed by the render name
  fit_nm <- suppressWarnings(tspa(
    "y ~ x + m + z + xm", data = fs_prod,
    se_fs = c(se_reg, xm = 0.11), product = TRUE))
  expect_identical(attr(fit_tok, "tspaModel"), attr(fit_nm, "tspaModel"))
  expect_equal(coef(fit_tok), coef(fit_nm))
  # the renamed entry is what the schema (and the replay) sees
  expect_equal(attr(fit_tok, "tspa_args")$se_fs[["xm"]], 0.11)
  expect_false("x.m" %in% colnames(attr(fit_tok, "tspa_args")$se_fs))
})

test_that("mf: a:b syntax fits identically to the concatenated syntax", {
  df <- prod_setup()
  for (method in c("regression", "Bartlett")) {
    fs <- get_fs(df, model = prod_model, std.lv = TRUE, method = method)
    fit_colon <- suppressWarnings(
      tspa("y ~ x + m + z + x:m + x:z + m:z", data = fs, product = TRUE))
    fit_concat <- suppressWarnings(
      tspa("y ~ x + m + z + xm + xz + mz", data = fs, product = TRUE))
    expect_identical(attr(fit_colon, "tspaModel"),
                     attr(fit_concat, "tspaModel"),
                     info = paste("model string,", method))
    expect_equal(coef(fit_colon), coef(fit_concat),
                 info = paste("coefs,", method))
    # the generated model carries the render name, not the interaction token
    m0 <- attr(fit_colon, "tspaModel")
    expect_no_match(m0, "x:m", fixed = TRUE)
    expect_match(m0, "xm =~ ", fixed = TRUE)
  }
})

test_that("a:b product fits replay identically via tspa_args (idempotent)", {
  df <- prod_setup()
  fs_plain <- get_fs(df, model = prod_model, std.lv = TRUE,
                     method = "Bartlett")
  se_reg <- c(y = fs_plain[1, "fs_y_se"], x = fs_plain[1, "fs_x_se"],
              m = fs_plain[1, "fs_m_se"], z = fs_plain[1, "fs_z_se"])
  fit_sf <- tspa("y ~ x + m + z + x:m + x:z + m:z", data = fs_plain,
                 se_fs = se_reg, product = TRUE)
  fit_sf_re <- do.call(tspa, attr(fit_sf, "tspa_args"))
  expect_identical(coef(fit_sf_re), coef(fit_sf))
  expect_identical(vcov(fit_sf_re), vcov(fit_sf))
  fit_mf <- suppressWarnings(
    tspa("y ~ x + m + z + x:m", data = fs_plain, product = TRUE))
  fit_mf_re <- do.call(tspa, attr(fit_mf, "tspa_args"))
  expect_identical(coef(fit_mf_re), coef(fit_mf))
  expect_identical(vcov(fit_mf_re), vcov(fit_mf))
})

# ---------------------------------------------------------------------------
# Group 8: pool_se_col() / pool_se_fs() refactor regression
# ---------------------------------------------------------------------------

test_that("pool_se_fs() matches the hand computation (refactor regression)", {
  set.seed(7)
  n <- 60
  g <- factor(rep(c("a", "b"), each = n / 2))
  dat <- data.frame(
    group = g,
    fs_v_se = c(runif(n / 2, 0.1, 0.3), runif(n / 2, 0.2, 0.5)),
    fs_w_se = c(runif(n / 2, 0.1, 0.4), runif(n / 2, 0.3, 0.6))
  )
  out <- pool_se_fs(dat, c("v", "w"), "mean", "group")
  expect_equal(out[["v"]], c(mean(dat$fs_v_se[g == "a"]),
                             mean(dat$fs_v_se[g == "b"])))
  expect_equal(out[["w"]], c(mean(dat$fs_w_se[g == "a"]),
                             mean(dat$fs_w_se[g == "b"])))
  expect_equal(rownames(out), c("a", "b"))
  # no group column: a named vector
  out_sg <- pool_se_fs(dat, c("v", "w"), "median", NULL)
  expect_equal(out_sg, c(v = median(dat$fs_v_se), w = median(dat$fs_w_se)))
  # pool_se_col() agrees with pool_se_fs() column-by-column
  expect_equal(pool_se_col(dat, "fs_v_se", "mean", "group"), out[["v"]])
  expect_equal(pool_se_col(dat, "fs_v_se", "median"),
               median(dat$fs_v_se))
})
