# Corrected Standard Errors

This vignette shows how to obtain Taylor-series corrected standard
errors in 2S-PA, accounting for the uncertainty in the weights for
obtaining the factor scores. It is described in [this
paper](https://doi.org/10.1007/s00181-020-01942-z).

The second-stage estimation treats as known the loading matrix and the
error covariance matrix of the factor scores as indicators of the true
latent variables. However, those are estimates based on the first-stage
measurement model. When the sample size is large and/or the reliability
of the factor scores is high, the uncertainty in the estimated loading
and error covariance matrix is generally negligible. Otherwise, a
first-order correction on the standard errors of the second-stage
estimation is possible, as illustrated below.

## First-order correction for SE

``` math
\hat V_{\gamma, c} = \hat V_{\gamma} + \boldsymbol{J}_\boldsymbol{\gamma}(\hat{\boldsymbol{\theta}}) \hat V_{\theta} \boldsymbol{J}_\boldsymbol{\gamma}(\hat{\boldsymbol{\theta}})^\top,
```

where $`\boldsymbol{J}_\boldsymbol{\gamma}`$ is the Jacobian matrix of
$`\hat{\boldsymbol{\gamma}}`$ with respect to $`\boldsymbol{\theta}`$,
or

``` math
\hat V_{\gamma, c} = \hat V_{\gamma} + (\boldsymbol{H}_\gamma)^{-1} \left(\frac{\partial^2 \ell}{\partial \theta \partial \gamma^\top}\right) \hat V_{\theta} \left(\frac{\partial^2 \ell}{\partial \theta \partial \gamma^\top}\right)^\top (\boldsymbol{H}_\gamma)^{-1},
```

where $`V_{\gamma}`$ is the naive covariance matrix of the structural
parameter estimates $`\hat{\boldsymbol{\gamma}}`$ assuming the
measurement model parameters $`\boldsymbol{\theta}`$ (the score loadings
and error variances) are known, $`\boldsymbol{H}_\gamma`$ is the Hessian
matrix of the log-likelihood $`\ell`$ with respect to
$`\hat{\boldsymbol{\gamma}}`$, and $`V_{\theta}`$ can be obtained in the
first-stage measurement model analysis.

``` r

library(lavaan)
library(R2spa)
library(numDeriv)
library(boot)
```

This example is based on the `PoliticalDemocracy` dataset, from the
lavaan tutorial (<https://lavaan.ugent.be/tutorial/sem.html>), also used
as an example in Bollen (1989).

## Separate Measurement Models

``` r

model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4

  # regressions
    dem60 ~ ind60
'
```

``` r

cfa_ind60 <- cfa("ind60 =~ x1 + x2 + x3", data = PoliticalDemocracy)
# Regression factor scores
fs1 <- get_fs(cfa_ind60, vfsLT = TRUE, format = "list")
# Delta method variance of (loading, error variance)
vldev1 <- attr(fs1, which = "vfsLT")
cfa_dem60 <- cfa("dem60 =~ y1 + y2 + y3 + y4",
                 data = PoliticalDemocracy)
# Regression factor scores
fs2 <- get_fs(cfa_dem60, vfsLT = TRUE, format = "list")
# Delta method variance of (loading, error variance)
vldev2 <- attr(fs2, which = "vfsLT")
fs_dat <- data.frame(
  fs_ind60 = fs1$fs_ind60,
  fs_dem60 = fs2$fs_dem60
)
# Combine sampling variance of loading and error variance
# Note: loadings first, then error variance
vldev <- block_diag(vldev1, vldev2)[c(1, 3, 2, 4), c(1, 3, 2, 4)]
```

``` r

# 2S-PA
# Assemble loadings
ld <- block_diag(attr(fs1, which = "fsL"),
                 attr(fs2, which = "fsL"))
ev <- block_diag(attr(fs1, which = "fsT"),
                 attr(fs2, which = "fsT"))
tspa_fit <- tspa(model = "dem60 ~ ind60",
                 data = fs_dat,
                 fsL = ld,
                 fsT = ev)
# Unadjusted covariance
vcov(tspa_fit)
#>              d60~60 i60~~6 d60~~6
#> dem60~ind60   0.131              
#> ind60~~ind60 -0.001  0.006       
#> dem60~~dem60 -0.006  0.000  0.477
# Adjusted covariance matrix
(vc_cor <- vcov_corrected(tspa_fit, vfsLT = vldev, which_free = c(1, 4, 5, 7)))
#>              d60~60 i60~~6 d60~~6
#> dem60~ind60   0.133              
#> ind60~~ind60 -0.001  0.006       
#> dem60~~dem60  0.001  0.000  0.518
# Corrected standard errors
sqrt(diag(vc_cor))
#>  dem60~ind60 ind60~~ind60 dem60~~dem60 
#>   0.36415174   0.07585532   0.71996885
```

The `engine` argument of
[`vcov_corrected()`](https://mmm-lab-um.github.io/R2spa/reference/vcov_corrected.md)
(and of `tspa(corrected_se = TRUE)`) selects how the Jacobian is
evaluated. The default `"analytic"` evaluates it refit-free and
deterministically as an influence-function closed form (PLAN 16,
sections 2.4 and 4.3), $`J = -H^{-1} C`$ with $`H`$ and $`C`$ obtained
by central-differencing the analytic log-likelihood score; it covers
single- and multi-group models, saturated and restricted (df \> 0)
structural models, and mean-structure models, and is a pure function of
the base fit + `vfsLT` (bit-reproducible, no refits). `"fd"` (central
finite differences, one stage-2 refit on each side of each free element)
is retained as the A/B reference and the transparent fallback for a
shape the analytic form cannot handle. The two agree to the
finite-difference noise floor whenever `"analytic"` applies.

Standardized solution

The corrected standard error of the `dem60 ~ ind60` coefficient, derived
by hand as the standard error of the standardized coefficient plus a
first-order delta-method term:

``` r

keep_dem_ind <- function(std) {
  std[std$lhs == "dem60" & std$rhs == "ind60" & std$op == "~", , drop = FALSE]
}
tspa_est_std <- function(ld_ev) {
  ld1 <- ld
  ev1 <- ev
  ld1[c(1, 4)] <- ld_ev[1:2]
  ev1[c(1, 4)] <- ld_ev[3:4]
  tfit <- tspa(model = "dem60 ~ ind60",
               data = fs_dat,
               fsL = ld1,
               fsT = ev1)
  keep_dem_ind(standardizedSolution(tfit))$est.std
}
tspa_est_std(c(ld[c(1, 4)], ev[c(1, 4)]))
#> [1] 0.4531693
# Jacobian of the standardized coefficient w.r.t. (loadings, error variances)
jac_std <- numDeriv::jacobian(tspa_est_std,
                              x = c(ld[c(1, 4)], ev[c(1, 4)]))
ss0 <- keep_dem_ind(standardizedSolution(tspa_fit))
sqrt(ss0$se^2 + jac_std %*% vldev %*% t(jac_std))
#>           [,1]
#> [1,] 0.1013887
```

The same value comes from the one-call option. Because the two factors
were fitted in separate models, only the two loadings and the two error
variances are free (positions 1, 4, 5, 7), so we pass the matching 4 x 4
submatrix of `vldev` together with `which_free`:

``` r

tspa_fit_corr <- tspa(model = "dem60 ~ ind60",
                      data = fs_dat,
                      fsL = ld,
                      fsT = ev,
                      vfsLT = vldev,
                      which_free = c(1, 4, 5, 7),
                      corrected_se = TRUE)
sc <- keep_dem_ind(standardizedSolution(tspa_fit_corr))
sc$est.std              # unchanged point estimate
#> [1] 0.4531693
sc$se                   # corrected standard SE (matches the hand-derived value)
#> [1] 0.1013887
attr(tspa_fit_corr, "tspa_corrected")   # TRUE
#> [1] TRUE
```

Compare to joint structural and measurement model (unstandardized
solution):

``` r

sem_fit <- sem(model, data = PoliticalDemocracy)
# Larger standard error
(vc_j <- 
  vcov(sem_fit)[c("dem60~ind60", "ind60~~ind60", "dem60~~dem60"),
                c("dem60~ind60", "ind60~~ind60", "dem60~~dem60")])
#>               dem60~ind60  ind60~~ind60  dem60~~dem60
#> dem60~ind60   0.143355743 -0.0033833569  0.0648673299
#> ind60~~ind60 -0.003383357  0.0075268147 -0.0001098306
#> dem60~~dem60  0.064867330 -0.0001098306  0.7655699561
sqrt(diag(vc_j))
#>  dem60~ind60 ind60~~ind60 dem60~~dem60 
#>   0.37862348   0.08675722   0.87496855
```

Bootstrap Standard Errors

``` r

run_tspa <- function(df, inds) {
  # local = TRUE: score each latent from its own (separate) measurement model
  fs_dat <- get_fs(object = df[inds, ],
                   model = "ind60 =~ x1 + x2 + x3
                            dem60 =~ y1 + y2 + y3 + y4",
                   local = TRUE,
                   se = "none", test = "none")
  # tspa() reads fsL/fsT/fsb and the per-score SEs from the get_fs() result
  tspa_fit <- tspa(model = "dem60 ~ ind60",
                   data = fs_dat,
                   test = "none")
  coef(tspa_fit)
}
boo <- boot(PoliticalDemocracy, statistic = run_tspa, R = 1999)
```

``` r

# Use MAD to downweigh outlying replications
boo$t |>
    apply(MARGIN = 2, FUN = mad) |>
    setNames(c("dem60~ind60", "ind60~~ind60", "dem60~~dem60"))
#>  dem60~ind60 ind60~~ind60 dem60~~dem60 
#>   0.32521386   0.07271009   0.87676911
```

The standard errors seem to diverge slightly among methods.
SE(`dem60~ind60`) with bootstrap is the lowest. SE(`ind60~~ind60`) with
the joint model is particularly higher than the other two.
SE(`dem60~~dem60`) is the lowest with the corrected 2S-PA. It should be
pointed out that the joint model and 2S-PA are different estimators.

``` r

run_sem <- function(df, inds) {
  sem_fit <- sem(model, data = df[inds, ], se = "none", test = "none")
  coef(sem_fit)
}
boo <- boot(PoliticalDemocracy, statistic = run_sem, R = 999)
```

## Joint Measurement Model

``` r

cfa_joint <- cfa("ind60 =~ x1 + x2 + x3
                  dem60 =~ y1 + y2 + y3 + y4",
                 data = PoliticalDemocracy)
# Factor score
fs_joint <- get_fs_lavaan(cfa_joint, vfsLT = TRUE)
# Delta method variance
vldev_joint <- attr(fs_joint, which = "vfsLT")
```

``` r

tspa_fit2 <- tspa(model = "dem60 ~ ind60",
                  data = fs_joint,
                  fsT = attr(fs_joint, "fsT"),
                  fsL = attr(fs_joint, "fsL"))
# Unadjusted covariance
vcov(tspa_fit2)
#>              d60~60 i60~~6 d60~~6
#> dem60~ind60   0.124              
#> ind60~~ind60 -0.001  0.006       
#> dem60~~dem60 -0.006  0.000  0.430
# Adjusted covariance
(vc2_cor <- vcov_corrected(tspa_fit2, vfsLT = vldev_joint))
#>              d60~60 i60~~6 d60~~6
#> dem60~ind60   0.130              
#> ind60~~ind60 -0.001  0.006       
#> dem60~~dem60 -0.009  0.000  0.493
# Corrected standard errors
sqrt(diag(vc2_cor))
#>  dem60~ind60 ind60~~ind60 dem60~~dem60 
#>   0.36049259   0.07654825   0.70229676
```

Bootstrap Standard Errors

``` r

run_tspa2 <- function(df, inds) {
  fs_joint <- get_fs(df[inds, ],
                     model = "ind60 =~ x1 + x2 + x3
                              dem60 =~ y1 + y2 + y3 + y4",
                     se = "none", test = "none")
  # tspa() reads fsL/fsT/fsb and the per-score SEs from the get_fs() result
  tspa_fit2 <- tspa(model = "dem60 ~ ind60",
                    data = fs_joint,
                    test = "none")
  coef(tspa_fit2)
}
boo2 <- boot(PoliticalDemocracy, statistic = run_tspa2, R = 1999)
# run_tspa2b <- function(df, inds) {
#   fsb_joint <- get_fs(df[inds, ],
#                       model = "ind60 =~ x1 + x2 + x3
#                                dem60 =~ y1 + y2 + y3 + y4",
#                       se = "none", test = "none", method = "Bartlett")
#   tspa_fit2b <- tspa(model = "dem60 ~ ind60",
#                      data = fsb_joint,
#                      fsT = attr(fsb_joint, "fsT"),
#                      fsL = diag(2),
#                      test = "none")
#   coef(tspa_fit2b)
# }
# boo2b <- boot(PoliticalDemocracy, statistic = run_tspa2b, R = 4999)
# run_sem2 <- function(df, inds) {
#   sem_fit <- sem(model, data = df[inds, ])
#   coef(sem_fit)[c(6, 14, 15)]
# }
# boo2j <- boot(PoliticalDemocracy, statistic = run_sem2, R = 4999)
```

``` r

# Use MAD to downweigh outlying replications
boo2$t |>
    apply(MARGIN = 2, FUN = mad) |>
    setNames(c("dem60~ind60", "ind60~~ind60", "dem60~~dem60"))
#>  dem60~ind60 ind60~~ind60 dem60~~dem60 
#>   0.34364943   0.07851415   0.88168143
```

### With Bartlett’s Method

``` r

# Factor score
fsb_joint <- get_fs(PoliticalDemocracy,
                    model = "ind60 =~ x1 + x2 + x3
                             dem60 =~ y1 + y2 + y3 + y4",
                    method = "Bartlett",
                    vfsLT = TRUE)
vldevb_joint <- attr(fsb_joint, which = "vfsLT")
```

``` r

tspa_fit2b <- tspa(model = "dem60 ~ ind60",
                   data = fsb_joint,
                   fsT = attr(fsb_joint, "fsT"),
                   fsL = diag(2) |>
                       `dimnames<-`(list(c("fs_ind60", "fs_dem60"),
                                         c("ind60", "dem60"))))
# Unadjusted covariance matrix
vcov(tspa_fit2b)
#>              d60~60 i60~~6 d60~~6
#> dem60~ind60   0.124              
#> ind60~~ind60 -0.001  0.006       
#> dem60~~dem60 -0.006  0.000  0.430
# Adjusted covariance matrix
(vc2b_cor <- vcov_corrected(
    tspa_fit2b,
    # With Bartlett scoring the loadings are fixed, so propagate the error
    # variances only (positions 5 and 7)
    vfsLT = vldevb_joint[c(5, 7), c(5, 7)],
    which_free = c(5, 7)))
#>              d60~60 i60~~6 d60~~6
#> dem60~ind60   0.124              
#> ind60~~ind60 -0.001  0.006       
#> dem60~~dem60 -0.006  0.000  0.448
# Corrected standard errors
sqrt(diag(vc2b_cor))
#>  dem60~ind60 ind60~~ind60 dem60~~dem60 
#>   0.35271754   0.07634634   0.66959805
```

## Multiple Groups

``` r

# Multigroup, three-factor example
mod <- "
  # latent variables
    visual =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed =~ x7 + x8 + x9
"
# Factor scores based on partial invariance
fit <- cfa(mod,
           data = HolzingerSwineford1939,
           std.lv = TRUE,
           group = "school",
           group.equal = c("loadings", "intercepts"),
           group.partial = c("visual=~x2", "x7~1"))
fs_dat <- get_fs_lavaan(fit, vfsLT = TRUE)
vldev <- attr(fs_dat, which = "vfsLT")
```

``` r

tspa_fit <- tspa(model = "visual ~~ textual + speed
                          textual ~~ speed",
                 data = fs_dat,
                 group = "school",
                 fsL = attr(fs_dat, which = "fsL"),
                 fsT = attr(fs_dat, which = "fsT"))
# Unadjusted covariance of the structural parameters (the 12 factor
# variances/covariances, 6 per group); the 6 factor-score intercepts (fs_*~1)
# are omitted.
keep <- grepl("~~", rownames(vcov(tspa_fit)))
vcov(tspa_fit)[keep, keep]
#>                     visual~~textual visual~~speed textual~~speed visual~~visual
#> visual~~textual         0.012032756   0.004143016    0.003562702    0.009059744
#> visual~~speed           0.004143016   0.015289623    0.005749677    0.006286186
#> textual~~speed          0.003562702   0.005749677    0.012308431    0.002169626
#> visual~~visual          0.009059744   0.006286186    0.002169626    0.026249333
#> textual~~textual        0.007226793   0.002111028    0.004878946    0.003126897
#> speed~~speed            0.001464756   0.006962597    0.006774549    0.001505415
#> visual~~textual.g2      0.000000000   0.000000000    0.000000000    0.000000000
#> visual~~speed.g2        0.000000000   0.000000000    0.000000000    0.000000000
#> textual~~speed.g2       0.000000000   0.000000000    0.000000000    0.000000000
#> visual~~visual.g2       0.000000000   0.000000000    0.000000000    0.000000000
#> textual~~textual.g2     0.000000000   0.000000000    0.000000000    0.000000000
#> speed~~speed.g2         0.000000000   0.000000000    0.000000000    0.000000000
#>                     textual~~textual speed~~speed visual~~textual.g2
#> visual~~textual          0.007226793  0.001464756        0.000000000
#> visual~~speed            0.002111028  0.006962597        0.000000000
#> textual~~speed           0.004878946  0.006774549        0.000000000
#> visual~~visual           0.003126897  0.001505415        0.000000000
#> textual~~textual         0.016702352  0.001425195        0.000000000
#> speed~~speed             0.001425195  0.032202255        0.000000000
#> visual~~textual.g2       0.000000000  0.000000000        0.011424159
#> visual~~speed.g2         0.000000000  0.000000000        0.005834597
#> textual~~speed.g2        0.000000000  0.000000000        0.006283615
#> visual~~visual.g2        0.000000000  0.000000000        0.008537899
#> textual~~textual.g2      0.000000000  0.000000000        0.007689097
#> speed~~speed.g2          0.000000000  0.000000000        0.003681750
#>                     visual~~speed.g2 textual~~speed.g2 visual~~visual.g2
#> visual~~textual          0.000000000       0.000000000       0.000000000
#> visual~~speed            0.000000000       0.000000000       0.000000000
#> textual~~speed           0.000000000       0.000000000       0.000000000
#> visual~~visual           0.000000000       0.000000000       0.000000000
#> textual~~textual         0.000000000       0.000000000       0.000000000
#> speed~~speed             0.000000000       0.000000000       0.000000000
#> visual~~textual.g2       0.005834597       0.006283615       0.008537899
#> visual~~speed.g2         0.020016361       0.008699815       0.010689107
#> textual~~speed.g2        0.008699815       0.016930572       0.004219633
#> visual~~visual.g2        0.010689107       0.004219633       0.021628068
#> textual~~textual.g2      0.002940790       0.006708957       0.003370422
#> speed~~speed.g2          0.017174234       0.011969242       0.005282811
#>                     textual~~textual.g2 speed~~speed.g2
#> visual~~textual             0.000000000     0.000000000
#> visual~~speed               0.000000000     0.000000000
#> textual~~speed              0.000000000     0.000000000
#> visual~~visual              0.000000000     0.000000000
#> textual~~textual            0.000000000     0.000000000
#> speed~~speed                0.000000000     0.000000000
#> visual~~textual.g2          0.007689097     0.003681750
#> visual~~speed.g2            0.002940790     0.017174234
#> textual~~speed.g2           0.006708957     0.011969242
#> visual~~visual.g2           0.003370422     0.005282811
#> textual~~textual.g2         0.017541486     0.002565923
#> speed~~speed.g2             0.002565923     0.055832832
# Adjusted covariance
vcov_corrected(tspa_fit, vfsLT = vldev)[keep, keep]
#>                     visual~~textual visual~~speed textual~~speed visual~~visual
#> visual~~textual        1.653582e-02  4.136616e-03   2.336324e-03   0.0046883897
#> visual~~speed          4.136616e-03  3.426164e-02   6.309753e-03   0.0008815580
#> textual~~speed         2.336324e-03  6.309753e-03   1.950154e-02   0.0026658224
#> visual~~visual         4.688390e-03  8.815580e-04   2.665822e-03   0.0503725941
#> textual~~textual       6.491077e-03  1.913223e-03   4.460499e-03   0.0041860447
#> speed~~speed           1.728255e-03  1.053080e-03   3.422944e-03   0.0036590340
#> visual~~textual.g2     1.170824e-04  5.276775e-05  -2.019370e-05  -0.0019412372
#> visual~~speed.g2       5.360119e-05  3.321316e-05   3.940979e-06  -0.0004494630
#> textual~~speed.g2     -7.386367e-06 -1.044665e-05   2.627310e-05   0.0005816077
#> visual~~visual.g2     -1.897690e-05 -1.372446e-04  -4.954036e-04  -0.0095441033
#> textual~~textual.g2   -1.709957e-04 -1.074418e-04  -1.663762e-04  -0.0007817335
#> speed~~speed.g2       -2.699262e-04 -3.577680e-04  -1.063824e-04   0.0028169442
#>                     textual~~textual  speed~~speed visual~~textual.g2
#> visual~~textual         6.491077e-03  0.0017282555       1.170824e-04
#> visual~~speed           1.913223e-03  0.0010530796       5.276775e-05
#> textual~~speed          4.460499e-03  0.0034229443      -2.019370e-05
#> visual~~visual          4.186045e-03  0.0036590340      -1.941237e-03
#> textual~~textual        1.738704e-02  0.0018292298       1.174471e-04
#> speed~~speed            1.829230e-03  0.0545698797       6.050726e-04
#> visual~~textual.g2      1.174471e-04  0.0006050726       1.486744e-02
#> visual~~speed.g2        4.329687e-05 -0.0003043030       5.921234e-03
#> textual~~speed.g2       2.061742e-05 -0.0008122445       5.326339e-03
#> visual~~visual.g2      -3.171240e-05  0.0013761982       1.179427e-02
#> textual~~textual.g2    -2.256062e-04  0.0003306567       8.528362e-03
#> speed~~speed.g2         4.328671e-05 -0.0136523556       2.263978e-03
#>                     visual~~speed.g2 textual~~speed.g2 visual~~visual.g2
#> visual~~textual         5.360119e-05     -7.386367e-06     -0.0000189769
#> visual~~speed           3.321316e-05     -1.044665e-05     -0.0001372446
#> textual~~speed          3.940979e-06      2.627310e-05     -0.0004954036
#> visual~~visual         -4.494630e-04      5.816077e-04     -0.0095441033
#> textual~~textual        4.329687e-05      2.061742e-05     -0.0000317124
#> speed~~speed           -3.043030e-04     -8.122445e-04      0.0013761982
#> visual~~textual.g2      5.921234e-03      5.326339e-03      0.0117942702
#> visual~~speed.g2        2.731199e-02      1.020675e-02      0.0080868149
#> textual~~speed.g2       1.020675e-02      1.979243e-02      0.0036942385
#> visual~~visual.g2       8.086815e-03      3.694239e-03      0.0517510845
#> textual~~textual.g2     3.019288e-03      6.413910e-03      0.0057608769
#> speed~~speed.g2         1.568174e-02      1.374092e-02      0.0004145391
#>                     textual~~textual.g2 speed~~speed.g2
#> visual~~textual           -0.0001709957   -2.699262e-04
#> visual~~speed             -0.0001074418   -3.577680e-04
#> textual~~speed            -0.0001663762   -1.063824e-04
#> visual~~visual            -0.0007817335    2.816944e-03
#> textual~~textual          -0.0002256062    4.328671e-05
#> speed~~speed               0.0003306567   -1.365236e-02
#> visual~~textual.g2         0.0085283624    2.263978e-03
#> visual~~speed.g2           0.0030192875    1.568174e-02
#> textual~~speed.g2          0.0064139098    1.374092e-02
#> visual~~visual.g2          0.0057608769    4.145391e-04
#> textual~~textual.g2        0.0187518537    1.971515e-03
#> speed~~speed.g2            0.0019715145    9.585557e-02
```

The correction can be requested in place on
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md), so
[`vcov()`](https://rdrr.io/r/stats/vcov.html) (and hence the SEs derived
from it) and
[`standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
on the returned fit report the corrected standard errors (the point
estimates are unchanged):

``` r

tspa_fit_corr <- tspa(model = "visual ~~ textual + speed
                               textual ~~ speed",
                      data = fs_dat,
                      group = "school",
                      fsL = attr(fs_dat, which = "fsL"),
                      fsT = attr(fs_dat, which = "fsT"),
                      vfsLT = vldev,
                      corrected_se = TRUE)
attr(tspa_fit_corr, "tspa_corrected")   # TRUE
#> [1] TRUE
round(cbind(unadjusted = sqrt(diag(vcov(tspa_fit))),
            corrected    = sqrt(diag(vcov(tspa_fit_corr)))), 4)
#>                     unadjusted corrected
#> visual~~textual         0.1097    0.1286
#> visual~~speed           0.1237    0.1851
#> textual~~speed          0.1109    0.1396
#> visual~~visual          0.1620    0.2244
#> textual~~textual        0.1292    0.1319
#> speed~~speed            0.1794    0.2336
#> fs_visual~1             0.0682    0.0682
#> fs_textual~1            0.0751    0.0751
#> fs_speed~1              0.0646    0.0646
#> visual~~textual.g2      0.1069    0.1219
#> visual~~speed.g2        0.1415    0.1653
#> textual~~speed.g2       0.1301    0.1407
#> visual~~visual.g2       0.1471    0.2275
#> textual~~textual.g2     0.1324    0.1369
#> speed~~speed.g2         0.2363    0.3096
#> fs_visual~1.g2          0.0648    0.0648
#> fs_textual~1.g2         0.0776    0.0776
#> fs_speed~1.g2           0.0920    0.0920
```

### Corrected grand-standardized solution (multigroup)

[`grandStandardizedSolution()`](https://mmm-lab-um.github.io/R2spa/reference/grand_standardized_solution.md)
pools across groups (grand mean / grand SD) rather than standardizing
each group separately as
[`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
does. Feeding it the corrected fit’s covariance yields the corrected
grand-standardized SEs; we use a structural path model so there is a
regression to standardize:

``` r

path_mg <- "textual ~ visual + speed
            visual  ~ speed"
tspa_path  <- tspa(path_mg, data = fs_dat, group = "school",
                   fsL = attr(fs_dat, "fsL"), fsT = attr(fs_dat, "fsT"))
tspa_path_corr <- tspa(path_mg, data = fs_dat, group = "school",
                       fsL = attr(fs_dat, "fsL"), fsT = attr(fs_dat, "fsT"),
                       vfsLT = vldev, corrected_se = TRUE)
gr0 <- grandStandardizedSolution(tspa_path)       # uncorrected
gr1 <- grandStandardizedSolution(tspa_path_corr)  # corrected (SE only)
keep <- gr1$op == "~"
round(cbind(grand_std = gr1$est.std[keep],
            se_plain  = gr0$se[keep],
            se_corr   = gr1$se[keep]), 4)
#>      grand_std se_plain se_corr
#> [1,]    0.4143   0.0970  1.8119
#> [2,]    0.2100   0.1232  0.3635
#> [3,]    0.4020   0.1355  1.2168
#> [4,]    0.5244   0.1341 21.2261
#> [5,]    0.0683   0.1116 16.7520
#> [6,]    0.4658   0.0878  2.1946
# SE-only: identical point estimates, corrected SEs >= uncorrected
all(gr1$est.std[keep] == gr0$est.std[keep])
#> [1] TRUE
all(gr1$se[keep] >= gr0$se[keep])
#> [1] TRUE
```

``` r

keep_dem_ind <- function(std) {
  std[std$lhs %in% c("visual", "textual", "speed") &
    std$rhs %in% c("visual", "textual", "speed") & 
      std$lhs != std$rhs, , drop = FALSE]
}
tspa_est_std <- function(ld_ev) {
  fsL1 <- attr(fs_dat, which = "fsL")
  fsT1 <- attr(fs_dat, which = "fsT")
  fsL1[[1]][] <- ld_ev[1:9]
  fsL1[[2]][] <- ld_ev[10:18]
  fsT1[[1]][] <- ld_ev[c(19:21, 20, 22:23, 21, 23, 24)]
  fsT1[[2]][] <- ld_ev[c(25:27, 26, 28:29, 27, 29, 30)]
  tfit <- tspa(model = "visual ~~ textual + speed
                        textual ~~ speed",
               data = fs_dat,
               group = "school",
               fsL = fsL1,
               fsT = fsT1) 
  keep_dem_ind(standardizedSolution(tfit))$est.std
}
ld_ev0 <- c(attr(fs_dat, which = "fsL")[[1]],
            attr(fs_dat, which = "fsL")[[2]],
            attr(fs_dat, which = "fsT")[[1]][c(1:3, 5:6, 9)],
            attr(fs_dat, which = "fsT")[[2]][c(1:3, 5:6, 9)])
tspa_est_std(ld_ev0)
#> [1] 0.4939729 0.3427430 0.3334265 0.5425498 0.5412629 0.3472992
jac <- numDeriv::jacobian(tspa_est_std, x = ld_ev0)
v0 <- lavInspect(tspa_fit, "vcov.std.lv")[c(1:3, 10:12), c(1:3, 10:12)]
sqrt(diag(jac %*% vldev %*% t(jac)) / diag(v0))
#>    visual~~textual      visual~~speed     textual~~speed visual~~textual.g2 
#>          1.1139120          1.3925355          0.9611655          1.4150576 
#>   visual~~speed.g2  textual~~speed.g2 
#>          1.6964048          0.6082279
```
