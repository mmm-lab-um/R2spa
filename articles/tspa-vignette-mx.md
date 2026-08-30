# 2S-PA with OpenMx and IRT (mirt)

``` r

library(lavaan)
library(R2spa)
library(OpenMx)
library(mirt)
```

Two-stage path analysis (2S-PA) can be fit with `OpenMx` as well as the
default `lavaan` backend.
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
takes a `lavaan`-style structural model string plus the (per-row)
factor-score loading / error information from
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) and
returns a fitted `OpenMx` RAM model.

## With a Linear Factor Model

The example is from <https://lavaan.ugent.be/tutorial/sem.html>.

``` r

# CFA
my_cfa <- "
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
"
(fs_dat <- get_fs(PoliticalDemocracy, model = my_cfa, std.lv = TRUE,
                  format = "list")) |> head()
#>     fs_ind60   fs_dem60 fs_ind60_se fs_dem60_se ind60_by_fs_ind60
#> 1 -0.8101568 -1.3119114   0.1859987   0.3012901         0.9553858
#> 2  0.1888466 -1.3644831   0.1859987   0.3012901         0.9553858
#> 3  1.0960931  1.3107705   0.1859987   0.3012901         0.9553858
#> 4  1.8702043  1.4849083   0.1859987   0.3012901         0.9553858
#> 5  1.2446060  0.9193277   0.1859987   0.3012901         0.9553858
#> 6  0.3348621  0.4886331   0.1859987   0.3012901         0.9553858
#>   ind60_by_fs_dem60 dem60_by_fs_ind60 dem60_by_fs_dem60 ev_fs_ind60
#> 1        0.05816994        0.01834111         0.8688889  0.03459552
#> 2        0.05816994        0.01834111         0.8688889  0.03459552
#> 3        0.05816994        0.01834111         0.8688889  0.03459552
#> 4        0.05816994        0.01834111         0.8688889  0.03459552
#> 5        0.05816994        0.01834111         0.8688889  0.03459552
#> 6        0.05816994        0.01834111         0.8688889  0.03459552
#>   ecov_fs_dem60_fs_ind60 ev_fs_dem60
#> 1            0.004017388  0.09077571
#> 2            0.004017388  0.09077571
#> 3            0.004017388  0.09077571
#> 4            0.004017388  0.09077571
#> 5            0.004017388  0.09077571
#> 6            0.004017388  0.09077571
```

``` r

tspa_fit <- tspa(model = "dem60 ~ ind60", data = fs_dat,
                 fsT = attr(fs_dat, "fsT"),
                 fsL = attr(fs_dat, "fsL"))
parameterEstimates(tspa_fit)
#>         lhs op      rhs   est    se     z pvalue ci.lower ci.upper
#> 1     ind60 =~ fs_ind60 0.955 0.000    NA     NA    0.955    0.955
#> 2     ind60 =~ fs_dem60 0.058 0.000    NA     NA    0.058    0.058
#> 3     dem60 =~ fs_ind60 0.018 0.000    NA     NA    0.018    0.018
#> 4     dem60 =~ fs_dem60 0.869 0.000    NA     NA    0.869    0.869
#> 5  fs_ind60 ~~ fs_ind60 0.035 0.000    NA     NA    0.035    0.035
#> 6  fs_ind60 ~~ fs_dem60 0.004 0.000    NA     NA    0.004    0.004
#> 7  fs_dem60 ~~ fs_dem60 0.091 0.000    NA     NA    0.091    0.091
#> 8     dem60  ~    ind60 0.460 0.113 4.089      0    0.240    0.681
#> 9     ind60 ~~    ind60 1.000 0.169 5.900      0    0.668    1.332
#> 10    dem60 ~~    dem60 0.788 0.150 5.267      0    0.495    1.081
```

### Using OpenMx

The measurement information is read straight from the `fsL` / `fsT`
attributes of the
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
output. The score columns (`fs_ind60`, `fs_dem60`) and the latent names
(`ind60`, `dem60`) must already be present in `data` / the model string.

``` r

tspa_mx <- tspa_mx_model(model = "dem60 ~ ind60", data = fs_dat,
                         fsL = attr(fs_dat, "fsL"),
                         fsT = attr(fs_dat, "fsT"))
#> Running m1 with 7 parameters
# Run OpenMx
tspa_mx_fit <- mxRun(tspa_mx)
#> Running m1 with 7 parameters
#> Warning: In model 'm1' Optimizer returned a non-zero status code 5. The Hessian
#> at the solution does not appear to be convex. See ?mxCheckIdentification for
#> possible diagnosis (Mx status RED).
# Summarize the results
summary(tspa_mx_fit)
#> Summary of m1 
#>  
#> The Hessian at the solution does not appear to be convex. See ?mxCheckIdentification for possible diagnosis (Mx status RED). 
#>  
#> free parameters:
#>        name matrix   row      col    Estimate Std.Error A
#> 1 m1.A[4,3]      A dem60    ind60  0.46046539 0.1126082  
#> 2 m1.S[3,3]      S ind60    ind60  1.00000048 0.1694787  
#> 3 m1.S[4,4]      S dem60    dem60  0.78797115 0.1495929  
#> 4 m1.M[1,1]      M     1 fs_ind60 -0.07735089        NA  
#> 5 m1.M[1,2]      M     1 fs_dem60 -0.02589978        NA  
#> 6 m1.M[1,3]      M     1    ind60  0.08049424        NA  
#> 7 m1.M[1,4]      M     1    dem60 -0.01264573        NA  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              7                    143              393.7496
#>    Saturated:              5                    145                    NA
#> Independence:              4                    146                    NA
#> Number of observations/statistics: 75/150
#> 
#> 
#> ** Information matrix is not positive definite (not at a candidate optimum).
#>   Be suspicious of these results. At minimum, do not trust the standard errors.
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       107.7496               407.7496                 409.4212
#> BIC:      -223.6512               423.9720                 401.9099
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 12:16:34 
#> Wall clock time: 0.03143764 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

#### With Mean Structure

``` r

my_cfa <- "
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    x1 + y1 ~ 0
    ind60 + dem60 ~ 1
"
fs_dat <- get_fs(PoliticalDemocracy, model = my_cfa,
                 meanstructure = TRUE,
                 format = "list")
```

lavaan

``` r

tspa_fit <- tspa(model = "dem60 ~ ind60\ndem60 + ind60 ~ 1", data = fs_dat,
                 fsT = attr(fs_dat, "fsT"),
                 fsL = attr(fs_dat, "fsL"),
                 fsb = attr(fs_dat, "fsb"))
parameterEstimates(tspa_fit, standardized = TRUE)
#>         lhs op      rhs    est    se      z pvalue ci.lower ci.upper std.lv
#> 1     ind60 =~ fs_ind60  0.955 0.000     NA     NA    0.955    0.955  0.640
#> 2     ind60 =~ fs_dem60  0.182 0.000     NA     NA    0.182    0.182  0.122
#> 3     dem60 =~ fs_ind60  0.006 0.000     NA     NA    0.006    0.006  0.012
#> 4     dem60 =~ fs_dem60  0.869 0.000     NA     NA    0.869    0.869  1.819
#> 5  fs_ind60 ~~ fs_ind60  0.016 0.000     NA     NA    0.016    0.016  0.016
#> 6  fs_ind60 ~~ fs_dem60  0.006 0.000     NA     NA    0.006    0.006  0.006
#> 7  fs_dem60 ~~ fs_dem60  0.398 0.000     NA     NA    0.398    0.398  0.398
#> 8  fs_ind60 ~1           0.193 0.000     NA     NA    0.193    0.193  0.193
#> 9  fs_dem60 ~1          -0.203 0.000     NA     NA   -0.203   -0.203 -0.203
#> 10    dem60  ~    ind60  1.439 0.352  4.089  0.000    0.749    2.129  0.460
#> 11    dem60 ~1          -1.810 1.794 -1.009  0.313   -5.327    1.706 -0.865
#> 12    ind60 ~1           5.054 0.079 64.155  0.000    4.900    5.209  7.547
#> 13    ind60 ~~    ind60  0.449 0.076  5.900  0.000    0.300    0.598  1.000
#> 14    dem60 ~~    dem60  3.453 0.656  5.267  0.000    2.168    4.738  0.788
#>    std.all
#> 1    0.973
#> 2    0.061
#> 3    0.019
#> 4    0.918
#> 5    0.036
#> 6    0.072
#> 7    0.101
#> 8    0.294
#> 9   -0.102
#> 10   0.460
#> 11  -0.865
#> 12   7.547
#> 13   1.000
#> 14   0.788
```

OpenMx

``` r

tspa_mx <- tspa_mx_model(model = "dem60 ~ ind60\ndem60 + ind60 ~ 1",
                         data = fs_dat)
#> Running m1 with 5 parameters
# We recently supported dropping fsL, fsT, and fsb when they can
# be derived from the get_fs() result's attributes. So we skip
# tspa_mx_model(...,     fsL = attr(fs_dat, "fsL"),
#                        fsT = attr(fs_dat, "fsT"),
#                        fsb = attr(fs_dat, "fsb"))
# Run OpenMx
tspa_mx_fit <- mxRun(tspa_mx)
#> Running m1 with 5 parameters
# Summarize the results
summary(tspa_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix   row   col   Estimate  Std.Error A
#> 1 m1.A[4,3]      A dem60 ind60  1.4393156 0.35196271  
#> 2 m1.S[3,3]      S ind60 ind60  0.4485415 0.07601938  
#> 3 m1.S[4,4]      S dem60 dem60  3.4532713 0.65558361  
#> 4 m1.M[1,3]      M     1 ind60  5.2585434 0.07878429  
#> 5 m1.M[1,4]      M     1 dem60 -2.3798659 1.86521218  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              5                    145              444.4392
#>    Saturated:              5                    145                    NA
#> Independence:              4                    146                    NA
#> Number of observations/statistics: 75/150
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       154.4392               454.4392                 455.3088
#> BIC:      -181.5966               466.0266                 450.2680
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 12:16:35 
#> Wall clock time: 0.0167954 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

The two routes fit the same model, so the path, variance, and covariance
estimates (and their SEs) agree. They differ only in the unidentifiable
split of the mean structure between the corrected latents and their
(observed) factor-score indicators:
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) fixes
the exogenous latent mean at zero and lets the factor-score mean
estimate the data value, whereas
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
fixes the score residual means at zero and estimates the latent means
(the latent mean carries the data value). Read the means accordingly —
compare the two routes on the covariance quantities, not on how the mean
is split.

### Compare to joint model

``` r

jreg <- sem(
    "
    # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    # latent regression
    dem60 ~ ind60
    ",
    data = PoliticalDemocracy, std.lv = TRUE)
coef(jreg)["dem60~ind60"]
#> dem60~ind60 
#>   0.5187309
```

## Combined with IRT

Example from Lai & Hsiao (2022, Psychological Methods)

### Not accounting for error

``` r

# Simulate data with mirt
set.seed(1235)
num_obs <- 1000
# Simulate theta
eta <- MASS::mvrnorm(num_obs, mu = c(0, 0), Sigma = diag(c(1, 1 - 0.5^2)),
                     empirical = TRUE)
th1 <- eta[, 1]
th2 <- -1 + 0.5 * th1 + eta[, 2]
# items and response data
a1 <- matrix(1, 10)
d1 <- matrix(rnorm(10))
a2 <- matrix(runif(10, min = 0.5, max = 1.5))
d2 <- matrix(rnorm(10))
dat1 <- simdata(a = a1, d = d1, N = num_obs, itemtype = "2PL", Theta = th1)
dat2 <- simdata(a = a2, d = d2, N = num_obs, itemtype = "2PL", Theta = th2)
# Factor scores
mod1 <- mirt(dat1, model = 1, itemtype = "Rasch", verbose = FALSE)
mod2 <- mirt(dat2, model = 1, itemtype = "2PL", verbose = FALSE)
fs1 <- fscores(mod1, full.scores.SE = TRUE)
fs2 <- fscores(mod2, full.scores.SE = TRUE)
lm(fs2[, 1] ~ fs1[, 1])  # attenuated coefficient
#> 
#> Call:
#> lm(formula = fs2[, 1] ~ fs1[, 1])
#> 
#> Coefficients:
#> (Intercept)     fs1[, 1]  
#>  -3.730e-06    2.812e-01
```

### Factor-score intercepts via get_fs()

The `mirt` fit above is also scored directly by
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md),
which returns the EAP scores together with their per-row standard
errors, implied loadings and error variances (and, via
[`fs_indiv()`](https://mmm-lab-um.github.io/R2spa/reference/fs_indiv.md),
the definition variables) — so the hand-rolled `fs_eap` frame in the
next section is optional. A factor prior mean can be given with
`prior_mean` (mirt’s default is a zero-mean factor). When the mean is
not zero, the per-row factor-score intercepts stored in the `fsb`
attribute are non-zero: each `fsb_i` is the prior mean scaled by that
row’s shrinkage factor (the posterior variance `Vpost_i`), so the
intercepts — and the score means — track the prior mean, and vary from
row to row.

``` r

fs  <- get_fs(mod1)
fse <- get_fs(mod1, prior_mean = c(F1 = 2))
head(fs, 3)                                   # score / SE / loading / error variance
#>        fs_F1  fs_F1_se F1_by_fs_F1  ev_fs_F1
#> 1  0.4717843 0.4708481   0.6414382 0.2216979
#> 2  0.8269502 0.4766043   0.6197787 0.2271516
#> 3 -0.8717728 0.4707905   0.6416367 0.2216437
round(c(fsb.r1.default = unname(attr(fs, "fsb")[[1]]),
        fsb.r2.default = unname(attr(fs, "fsb")[[2]]),
        fsb.r1.prior2  = unname(attr(fse, "fsb")[[1]]),
        fsb.r2.prior2  = unname(attr(fse, "fsb")[[2]])), 4)
#> fsb.r1.default fsb.r2.default  fsb.r1.prior2  fsb.r2.prior2 
#>         0.0000         0.0000         0.8356         0.9396
round(c(score.mean.default = mean(fs$fs_F1),
        score.mean.prior2  = mean(fse$fs_F1)), 4)
#> score.mean.default  score.mean.prior2 
#>             0.0000             0.7596
```

Under mirt’s default zero-mean prior every `fsb` is zero (the scores are
mean zero); with `prior_mean = c(F1 = 2)` they are non-zero and differ
by row. These `fsb` values are exactly what
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
reads as its score-intercept row (“With Mean Structure” above).

[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
results can now also be passed to
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
directly: with all of `se_fs`/`fsL`/`fsT`/`fsb` omitted, the measurement
inputs are derived from the result’s attributes — the per-row loadings
and error variances become definition variables referencing the result’s
own `*_by_*`/`ev_*`/`ecov_*` columns, and the `int_fs_*` score-intercept
columns are appended from the `fsb` attribute (the
[`fs_indiv()`](https://mmm-lab-um.github.io/R2spa/reference/fs_indiv.md)
round-trip, now automatic). The minimal model that exercises the
score-intercept row, on the nonzero-prior and the zero-mean results:

``` r

# All measurement inputs (per-row loadings, error variances, and the
# int_fs_* intercept columns) are derived from the get_fs() result:
tspa_mx_fit <- tspa_mx_model("F1 ~ 1", data = fse)
#> Running m1 with 2 parameters
summary(tspa_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix row col    Estimate  Std.Error A
#> 1 m1.S[2,2]      S  F1  F1  0.80251349 0.06993011  
#> 2 m1.M[1,2]      M   1  F1 -0.08697211 0.03855771  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              2                    998              2276.339
#>    Saturated:              2                    998                    NA
#> Independence:              2                    998                    NA
#> Number of observations/statistics: 1000/1000
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:        280.339               2280.339                 2280.351
#> BIC:      -4617.601               2290.155                 2283.802
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 12:16:36 
#> Wall clock time: 0.02426338 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
# The zero-mean result derives the same way (its fsb attribute is all-zero,
# so the appended int_fs_F1 column is all-zero):
tspa_mx_fit <- tspa_mx_model("F1 ~ 1", data = fs)
#> Running m1 with 2 parameters
summary(tspa_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix row col      Estimate  Std.Error A
#> 1 m1.S[2,2]      S  F1  F1  9.637746e-01 0.07128714  
#> 2 m1.M[1,2]      M   1  F1 -8.414985e-06 0.03899255  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              2                    998              2395.021
#>    Saturated:              2                    998                    NA
#> Independence:              2                    998                    NA
#> Number of observations/statistics: 1000/1000
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       399.0212               2399.021                 2399.033
#> BIC:     -4498.9185               2408.837                 2402.485
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 12:16:36 
#> Wall clock time: 0.02253413 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

Under `prior_mean = c(F1 = 2)` the appended `int_fs_F1` column is
nonzero and row-varying; for the zero-mean result it is all-zero. For a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result, then, the
[`fs_indiv()`](https://mmm-lab-um.github.io/R2spa/reference/fs_indiv.md) +
hand-rolled character matrices of the next section — and the explicit
`attr(..., "fsL")`-style arguments used elsewhere in this vignette — are
optional sugar, not a requirement; they remain the route when the
per-row quantities come from outside
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md).

### 2S-PA with EAP (shrinkage) scores

Because EAP (shrinkage) scores are used, we use the following properties
for these scores ($`\tilde \eta_i`$ for person $`i`$):

- Standard error (SE, from software output) =
  $`\sqrt{V(\eta) (1 - \rho_{\tilde \eta, i})}`$, where
  $`\rho_{\tilde \eta, i}`$ is the reliability of $`\tilde \eta_i`$
- Loading of $`\tilde \eta_i`$ on $`\eta`$ is equal to
  $`\rho_{\tilde \eta, i}`$ = $`1 - \text{SE}^2 / V(\eta)`$
- Standard error of measurement for $`\tilde \eta_i`$ is equal to
  $`\rho_{\tilde \eta, i} \text{SE}`$
- Based on the above, the total variance of $`\tilde \eta_i`$ is also
  $`\rho_{\tilde \eta, i}`$.

``` r

# Combine into data set
fs_eap <- cbind(fs1, fs2) |>
    as.data.frame() |>
    setNames(c("fs1", "se_fs1", "fs2", "se_fs2")) |>
    # Compute reliability and error variances
    within(expr = {
        rel_fs1 <- 1 - se_fs1^2
        rel_fs2 <- 1 - se_fs2^2
        ev_fs1 <- se_fs1^2 * (1 - se_fs1^2)
        ev_fs2 <- se_fs2^2 * (1 - se_fs2^2)
    })
```

The per-person reliability and error variance are supplied as
*definition variables*: each `fsL` / `fsT` cell holds the name of the
`data` column that carries that person’s value. The structural latent
names (`f1`, `f2`) are kept distinct from the score columns (`fs1`,
`fs2`).

``` r

cross_load <- matrix(c("rel_fs1", NA, NA, "rel_fs2"), nrow = 2) |>
    `dimnames<-`(list(c("fs1", "fs2"), c("f1", "f2")))
err_cov <- matrix(c("ev_fs1", NA, NA, "ev_fs2"), nrow = 2) |>
    `dimnames<-`(rep(list(c("fs1", "fs2")), 2))
tspa_mx <- tspa_mx_model(model = "f2 ~ f1\nf2 + f1 ~ 1",
                         data = fs_eap, fsL = cross_load, fsT = err_cov)
#> Running m1 with 5 parameters
# Run OpenMx
tspa_mx_fit <- mxRun(tspa_mx)
#> Running m1 with 5 parameters
# Summarize the results
summary(tspa_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix row col      Estimate  Std.Error A
#> 1 m1.A[4,3]      A  f2  f1  0.5004173870 0.05514387  
#> 2 m1.S[3,3]      S  f1  f1  0.9050845033 0.06795171  
#> 3 m1.S[4,4]      S  f2  f2  0.7782328697 0.07524372  
#> 4 m1.M[1,3]      M   1  f1 -0.0012615895 0.03807609  
#> 5 m1.M[1,4]      M   1  f2 -0.0004635611 0.04094092  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              5                   1995              4644.151
#>    Saturated:              5                   1995                    NA
#> Independence:              4                   1996                    NA
#> Number of observations/statistics: 1000/2000
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       654.1514               4654.151                 4654.212
#> BIC:     -9136.8204               4678.690                 4662.810
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 12:16:37 
#> Wall clock time: 0.1747611 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

The same model fits in one call from a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result: with the measurement arguments omitted,
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
derives the definition-variable matrices (and the `int_fs_*` intercept
columns) from the result’s own `*_by_*`/`ev_*`/`ecov_*` columns and
`fsb` attribute, as in the previous section. The hand-rolled route stays
when the per-person quantities come from outside
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) —
as here, where the per-person reliabilities and error variances are the
hand-computed EAP formulas above, not a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result’s implied-loadings / error-covariance columns.

### Joint model reference

The gold-standard comparison is a joint (one-step) MIMIC/SEM fit. With
2PL items this uses the polychoric-WLS approximation here; a full FIML
fit is computationally prohibitive at this scale.

``` r

dat <- cbind(dat1, dat2)
colnames(dat) <- paste0("i", 1:20)
wls_fit <- sem("
f1 =~ i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 + i10
f2 =~ i11 + i12 + i13 + i14 + i15 + i16 + i17 + i18 + i19 + i20
f2 ~ f1
", data = dat, ordered = TRUE, std.lv = TRUE)
coef(wls_fit)["f2~f1"]
#>     f2~f1 
#> 0.5490816
```

### Multidimensional Measurement Model

When a multidimensional model is used, for EAP scores we use the
following properties generalized from the unidimensional case.
Specifically, let $`V(\boldsymbol \eta)`$ = $`\boldsymbol{\psi}`$.

- Loading of $`\tilde{\boldsymbol \eta}_i`$ on $`\boldsymbol \eta`$ is
  $`\boldsymbol{\Lambda}_{\tilde \eta, i}`$, which is also the
  reliability matrix (i.e., the squared covariance between
  $`\tilde{\boldsymbol \eta}_i`$ and $`\boldsymbol \eta`$).
- Error matrix ($`\mathbf{E}`$, from software output) =
  $`(\mathbf{I} - \boldsymbol{\Lambda}_{\tilde \eta, i})\boldsymbol{\psi}`$.
  Therefore,
  $`\boldsymbol{\Lambda}_{\tilde \eta, i} = \mathbf{I} - \mathbf{E} \boldsymbol{\psi}^{-1}`$.
- Error variance of $`\tilde \eta_i`$ is
  $`\boldsymbol{\Lambda}_{\tilde \eta, i} \mathbf{E}`$.

``` r

dat <- cbind(dat1, dat2)
colnames(dat) <- paste0("Item_", 1:20)
# Multidimensional IRT
mod <- "
F1 = 1-10
F2 = 11-20
COV = F1*F2
CONSTRAIN = (1-10, a1)
"
mfit <- mirt(dat, model = mod, itemtype = "2PL", verbose = FALSE)
```

``` r

# Factor scores, plus the per-row loadings and error covariances and
# score intercepts get_fs() attaches (the same full-factor-covariance
# correction the EAP section above derives from the hand-computed
# reliabilities and error variances).
fs_dat <- get_fs(mfit)
```

``` r

# Structural model on the two latent dimensions (the mirt factor names, kept
# distinct from the score columns). All measurement inputs are derived from
# the get_fs() result's attributes: the per-row loadings and error
# covariances become definition variables referencing its own
# *_by_*/ev_*/ecov_* columns, and the int_fs_* score-intercept columns are
# appended from the fsb attribute (the fs_indiv() round-trip, now automatic):
tspa_mx2 <- tspa_mx_model(model = "F2 ~ F1\nF2 + F1 ~ 1", data = fs_dat)
#> Running m1 with 5 parameters
# Run OpenMx
tspa_mx2_fit <- mxRun(tspa_mx2)
#> Running m1 with 5 parameters
# Summarize the results
summary(tspa_mx2_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix row col      Estimate  Std.Error A
#> 1 m1.A[4,3]      A  F2  F1  4.753223e-01 0.05204783  
#> 2 m1.S[3,3]      S  F1  F1  1.000064e+00 0.07394681  
#> 3 m1.S[4,4]      S  F2  F2  7.740987e-01 0.07359209  
#> 4 m1.M[1,3]      M   1  F1 -9.402292e-07 0.03972034  
#> 5 m1.M[1,4]      M   1  F2  5.419113e-05 0.04081000  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              5                   1995              4327.057
#>    Saturated:              5                   1995                    NA
#> Independence:              4                   1996                    NA
#> Number of observations/statistics: 1000/2000
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:        337.057               4337.057                 4337.117
#> BIC:      -9453.915               4361.596                 4345.716
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 12:16:43 
#> Wall clock time: 0.2071025 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```
