# Product factor-score indicators (latent interactions)

``` r

library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(R2spa)
```

## Product factor-score indicators

[`compute_fs_prod()`](https://mmm-lab-um.github.io/R2spa/reference/compute_fs_prod.md)
creates double-mean-centered (DMC) product indicators for pairs of
first-order factor scores, together with their standard errors and
model-implied loadings. It is the successor of the (now removed) student
function `get_fs_int()`: the column conventions are kept (for latent
names `a` and `b` the columns are `fs_a:fs_b`, `fs_a:fs_b_se`, and
`fs_a:fs_b_ld`), while the standard error now uses the general
joint-model formula computed from the stage-1 `fsL`/`fsT`/`psi`
attributes of a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result, instead of the separate-single-factor special case the old
function used (see
[`?compute_fs_prod`](https://mmm-lab-um.github.io/R2spa/reference/compute_fs_prod.md)
for the derivation).

The product columns can be obtained in two equivalent ways:

- `get_fs(..., product = "a:b + c:d")` — factor scores and product
  columns in one call
- `compute_fs_prod(get_fs(...), product = "a:b + c:d")` — the same
  columns, computed from an existing
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result

and fed into stage 2 two ways:

- manually, with the product standard errors listed in `se_fs`
- `tspa(..., product = TRUE)` —
  [`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
  detects the product latents (a model latent named by concatenating two
  factor-score names, `xm`, or in lavaan’s interaction syntax, `x:m`),
  computes the missing product columns on the fly from the data’s
  stage-1 attributes, and wires the product standard errors in
  automatically; it also fixes the product indicators’ measurement-error
  covariances, computed from the stage-1 attributes (products sharing a
  factor score have correlated measurement errors)

Differences from `get_fs_int()`:

- the product pairs are always specified explicitly (there is no
  all-pairs default)
- the input must be a single-group
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result — the `fsL`/`fsT`/`psi` attributes carry the stage-1
  information. The old separate-measurement-models workflow (with a
  user-supplied `lat_var`) is not needed: a joint model, with or without
  `std.lv`, covers both cases, because the latent (co)variances are read
  from the model itself
- same-factor products (`"x:x"`) are not supported (v1)

## Illustrative example

We will simulate a dataset with four first-order latent variables: `x`,
`m`, `z`, `y`, where `y` is driven by the main effects and the three
pairwise latent interactions.

``` r

set.seed(1211)
## Sample size:
num_obs <- 5000

# Structural Parameters:
## gamma_x = 0.3, gamma_m = 0.4, gamma_z = 0.2
## gamma_xm = 0.1, gamma_xz = 0.15, gamma_mz = 0.12

# Correlation between latent variables:
## cor_xm = 0.1, cor_xz = 0.15, cor_zm = 0.12

# Data simulation:
## x, z, m, ey
## The residual variance is chosen so that the population variance of the
## `y` latent is exactly 1: under the correlations above,
## Var(0.3x + 0.4m + 0.2z + 0.1xm + 0.15xz + 0.12mz) = 0.41109961, so the
## residual variance is 1 - 0.41109961 = 0.58890039.
cov_xmz_ey <- matrix(c(1, 0.1, 0.15, 0,
                       0.1, 1, 0.12, 0,
                       0.15, 0.12, 1, 0,
                       0, 0, 0, 0.58890039), nrow = 4)
eta <- as.data.frame(
  MASS::mvrnorm(num_obs,
    mu = rep(0, 4), Sigma = cov_xmz_ey,
    empirical = FALSE
  )
)
names(eta) <- c("x", "m", "z", "ey")

# xm, xz, mz
eta <- eta |>
  transform(
    xm = x * m,
    xz = x * z,
    mz = m * z
  )

# y
etay <- 0.3 * eta$x + 0.4 * eta$m + 0.2 * eta$z +
  0.1 * eta$xm + 0.15 * eta$xz + 0.12 * eta$mz + eta$ey

# Observed Indicators
lambda_x <- c(0.9, 0.8, 0.7)
lambda_m <- c(0.85, 0.75, 0.65)
lambda_z <- c(0.8, 0.7, 0.6)
lambda_y <- c(0.75, 0.7, 0.65)

x_obs <- eta$x %*% t(lambda_x) + rnorm(num_obs * length(lambda_x))
m_obs <- eta$m %*% t(lambda_m) + rnorm(num_obs * length(lambda_m))
z_obs <- eta$z %*% t(lambda_z) + rnorm(num_obs * length(lambda_z))
y_obs <- etay %*% t(lambda_y) + rnorm(num_obs * length(lambda_y))

# Dataset: raw score
df <- cbind(x_obs, m_obs, z_obs, y_obs)
df <- as.data.frame(df)
names(df) <- c(
  paste0("x", 1:3), paste0("m", 1:3),
  paste0("z", 1:3), paste0("y", 1:3)
)
```

### Product indicators from a joint model

With the latent variances fixed to 1 (`std.lv = TRUE`), the factor
scores and the three product columns come from one call.
Double-mean-centering is used.

``` r

fs_dat <- get_fs(df, model = "x =~ x1 + x2 + x3
                              m =~ m1 + m2 + m3
                              z =~ z1 + z2 + z3
                              y =~ y1 + y2 + y3",
                 std.lv = TRUE,
                 method = "Bartlett",
                 product = "x:m + x:z + m:z")
head(fs_dat[c("fs_x", "fs_m", "fs_x:fs_m",
              "fs_x:fs_m_se", "fs_x:fs_m_ld")])
#>         fs_x       fs_m   fs_x:fs_m fs_x:fs_m_se fs_x:fs_m_ld
#> 1  1.5511355  2.4451070  3.66912919     1.201306            1
#> 2 -2.0135397 -0.5405614  0.96487878     1.201306            1
#> 3 -1.4696740  1.3977934 -2.17786370     1.201306            1
#> 4 -0.2264411 -0.9138223  0.08336388     1.201306            1
#> 5  1.9367880 -0.6619337 -1.40558826     1.201306            1
#> 6  0.9512282  1.7866185  1.57591888     1.201306            1
```

`fs_x:fs_m_se` is the per-row standard error (constant here, complete
data) and `fs_x:fs_m_ld` the implied loading (1 for Bartlett scores from
a joint model). The same columns can also be computed after the fact,
here for a subset of the pairs:

``` r

fs_dat2 <- compute_fs_prod(
  get_fs(df, model = "x =~ x1 + x2 + x3
                      m =~ m1 + m2 + m3
                      z =~ z1 + z2 + z3
                      y =~ y1 + y2 + y3",
         std.lv = TRUE,
         method = "Bartlett"),
  product = "x:m + m:z"
)
head(fs_dat2[c("fs_x:fs_m", "fs_m:fs_z",
               "fs_x:fs_m_se", "fs_m:fs_z_ld")])
#>     fs_x:fs_m   fs_m:fs_z fs_x:fs_m_se fs_m:fs_z_ld
#> 1  3.66912919  0.01977866     1.201306            1
#> 2  0.96487878  0.38755226     1.201306            1
#> 3 -2.17786370  0.90395574     1.201306            1
#> 4  0.08336388 -3.00806208     1.201306            1
#> 5 -1.40558826 -1.03335845     1.201306            1
#> 6  1.57591888 -2.77363791     1.201306            1
```

### 2S-PA with product factor-score indicators (manual)

``` r

tspa_fit <- tspa("
  y ~ x + m + z + xm + xz + mz
 ", data = fs_dat,
  se_fs = c(
    y  = fs_dat[1, "fs_y_se"],
    x  = fs_dat[1, "fs_x_se"],
    m  = fs_dat[1, "fs_m_se"],
    z  = fs_dat[1, "fs_z_se"],
    xm = fs_dat[1, "fs_x:fs_m_se"],
    xz = fs_dat[1, "fs_x:fs_z_se"],
    mz = fs_dat[1, "fs_m:fs_z_se"]
  ))
coef(tspa_fit)[c("y~x", "y~m", "y~z", "y~xm", "y~xz", "y~mz")]
#>       y~x       y~m       y~z      y~xm      y~xz      y~mz 
#> 0.2924104 0.4154073 0.2287085 0.1272880 0.1677773 0.1144332
```

The product score columns (`fs_x:fs_m`, etc.) are aliased to the plain
model names (`fs_xm`, etc.) automatically by
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md).

We report the unstandardized structural estimates: they are on the scale
of the latent variables and of their DMC products, and agree with the
simulated values (0.3, 0.4, 0.2 for the main effects; 0.1, 0.15, 0.12
for the interactions) for two reasons. First, the residual variance
above was chosen so that the population variance of the `y` latent is
exactly 1, so the `std.lv = TRUE` score scale is the simulated scale.
Second, the product indicators’ measurement errors are correlated —
products sharing a factor score (e.g. `x:m` and `x:z`) both contain that
score’s error — and
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) fixes
those error covariances in the stage-2 model (computed from the stage-1
`fsL`/`fsT`/`psi` attributes), so the interaction estimates are not
attenuated. The interaction terms are deliberately not standardized:
[`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
rescales each coefficient by the standard deviation of its predictor,
but for the product indicator that standard deviation is the SD of the
DMC product, not SD(x)·SD(m), so the conventional standardized
interaction effect b·SD(x)·SD(m)/SD(y) is not what
[`standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
reports for those rows.

### 2S-PA with `product = TRUE` (auto-compute)

With `product = TRUE`,
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
recognizes `x:m`, `x:z`, and `m:z` as the products of the score pairs
`x` + `m`, `x` + `z`, and `m` + `z` (rendered under the concatenated
names `xm`, `xz`, and `mz`, since in the generated model `x:m` would be
parsed as an interaction of the latent variables), computes the missing
product columns from the data’s stage-1 attributes, and joins the
product standard errors into `se_fs` — so only the non-product standard
errors need to be supplied:

``` r

fs_np <- get_fs(df, model = "x =~ x1 + x2 + x3
                             m =~ m1 + m2 + m3
                             z =~ z1 + z2 + z3
                             y =~ y1 + y2 + y3",
                std.lv = TRUE,
                method = "Bartlett")
tspa_auto <- tspa("
  y ~ x + m + z + x:m + x:z + m:z
 ", data = fs_np,
  se_fs = c(
    y = fs_np[1, "fs_y_se"],
    x = fs_np[1, "fs_x_se"],
    m = fs_np[1, "fs_m_se"],
    z = fs_np[1, "fs_z_se"]
  ),
  product = TRUE)
all.equal(coef(tspa_fit), coef(tspa_auto))
#> [1] TRUE
```

`TRUE` — the auto-computed fit (interaction-syntax model) is identical
to the manual one (concatenated names). With `se_fs` omitted entirely,
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) takes
the multi-factor route and the product latents get their implied
loadings and error variances from the stage-1 attributes, as in the next
section.

### Latent variances not fixed to 1

The examples so far fixed the latent variances to 1. Without `std.lv`,
the latent variances are estimated — and no separate latent-variance
input is needed for the product indicators: the standard-error formula
reads the latent (co)variances from the joint model’s `psi` attribute.
With `se_fs` omitted, `tspa(product = TRUE)` is fully automatic:

``` r

fs_np2 <- get_fs(df, model = "x =~ x1 + x2 + x3
                              m =~ m1 + m2 + m3
                              z =~ z1 + z2 + z3
                              y =~ y1 + y2 + y3",
                 method = "Bartlett")
tspa_fit2 <- tspa("
  y ~ x + m + z + x:m + x:z + m:z
 ", data = fs_np2, product = TRUE)
names.oplv <- c("y~x", "y~m", "y~z", "y~xm", "y~xz", "y~mz")
cbind(
  "std.lv = TRUE"  = coef(tspa_fit)[names.oplv],
  "std.lv = FALSE" = coef(tspa_fit2)[names.oplv]
)
#>      std.lv = TRUE std.lv = FALSE
#> y~x      0.2924104      0.2525789
#> y~m      0.4154073      0.3772234
#> y~z      0.2287085      0.2153787
#> y~xm     0.1272880      0.1292186
#> y~xz     0.1677773      0.1766308
#> y~mz     0.1144332      0.1266503
```

The unstandardized coefficients are on the estimated latent scale (the
latent variances are estimated, not fixed to 1), and the main effects
agree with the `std.lv = TRUE` example on the standardized scale — a
rescaling of the same model. (As before, the interaction terms are read
unstandardized: standardizing them by the SD of the product indicator is
not the conventional b·SD(x)·SD(m)/SD(y) form.)
