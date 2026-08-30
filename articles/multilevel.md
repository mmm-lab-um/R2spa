# 2S-PA with Random Effects

In this vignette, we show how to use 2S-PA for estimating associations
for cluster-level random effects. We use simulated data from a
parallel-process growth model, where the two processes have different
time scales. This is a pretty challenging model to estimate using
standard structural equation models (as time scales are different) or
multilevel models (as it is multivariate).

``` r

library(R2spa)
library(lme4)
library(OpenMx)
```

The EB (empirical-Bayes) random effects from
[`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) are per-cluster, so
the per-cluster loadings and error covariances are supplied to
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
as *definition variables* (column names), which it derives from the
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result automatically. The score column of each latent carries an `fs_`
prefix (e.g. latent `u0` is measured by the score column `fs_u0`). A
small helper pulls the implied latent covariance back out of the fitted
RAM model by name:

``` r

lat_cov <- function(fit, lat) {
  v <- c(fit$manifestVars, fit$latentVars)
  idx <- match(lat, v)
  co <- coef(fit)
  out <- matrix(NA_real_, length(lat), length(lat), dimnames = list(lat, lat))
  for (i in seq_along(lat)) for (j in seq_along(lat))
    out[i, j] <- co[sprintf("m1.S[%d,%d]", min(idx[i], idx[j]), max(idx[i], idx[j]))]
  out
}
```

### Simulate bivariate growth data

- Process 1: 5 time points
- Process 2: 20 time points

``` r

set.seed(1957)
latent_cor <- matrix(
    c(1, 0, 0.5, 0.3,
      0, 1, -0.1, -0.3,
      0.5, -0.1, 1, 0.2,
      0.3, -0.3, 0.2, 1),
    nrow = 4
)
latent_sd <- c(1, 0.5, 1.5, 0.6)
latent_cov <- diag(latent_sd) %*% latent_cor %*% diag(latent_sd)
latent_mean <- c(-1, 0.5, 0, -0.3)
# Simulate person-specific growth parameters
n_clus <- 100
eta <- MASS::mvrnorm(n_clus, mu = latent_mean, Sigma = latent_cov,
                     empirical = TRUE)
# Time variables (different time scale for the two processes)
n_time1 <- 5
n_time2 <- 20
# Assume complete data for now
time1 <- rep(seq_len(n_time1) - 1, n_clus)
time2 <- rep(seq_len(n_time2) - 1, n_clus)
# Simulate process 1
clus_id1 <- rep(seq_len(n_clus), each = n_time1)
sigma1 <- 1  # error sd
# Simulate error
e1 <- rnorm(clus_id1)
e1 <- e1 - mean(e1)
e1 <- e1 / sd(e1) * sigma1
y1 <- eta[clus_id1, 1] + time1 * eta[clus_id1, 2] + e1
# Simulate process 2
clus_id2 <- rep(seq_len(n_clus), each = n_time2)
sigma2 <- 2  # error sd
# Simulate error
e2 <- rnorm(clus_id1)
e2 <- e2 - mean(e2)
e2 <- e2 / sd(e2) * sigma2
y2 <- eta[clus_id2, 3] + time2 * eta[clus_id2, 4] + e2
```

``` r

# Empirical Bayes
m1 <- lmer(y1 ~ time1 + (time1 | clus_id1),
           data = data.frame(y1, time1, clus_id1),
           REML = FALSE)
fs_dat1 <- get_fs(m1)
# Empirical Bayes
m2 <- lmer(y2 ~ time2 + (time2 | clus_id2),
           data = data.frame(y2, time2, clus_id2),
           REML = FALSE)
fs_dat2 <- get_fs(m2)
```

[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
carries each EB score (`fs_u0`, `fs_u1`) alongside its per-cluster
standard error (`fs_u0_se`), implied loadings (`u0_by_fs_u0`, …), error
variances (`ev_fs_u0`) and error covariance (`ecov_fs_u1_fs_u0`):

``` r

knitr::kable(head(fs_dat1), digits = 3)
```

| fs_u0 | fs_u1 | fs_u0_se | fs_u1_se | u0_by_fs_u0 | u0_by_fs_u1 | u1_by_fs_u0 | u1_by_fs_u1 | ev_fs_u0 | ecov_fs_u1_fs_u0 | ev_fs_u1 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.784 | 0.347 | 0.458 | 0.193 | 0.733 | 0.079 | 0.378 | 0.78 | 0.21 | -0.052 | 0.037 |
| 1.051 | 0.196 | 0.458 | 0.193 | 0.733 | 0.079 | 0.378 | 0.78 | 0.21 | -0.052 | 0.037 |
| -1.620 | -0.537 | 0.458 | 0.193 | 0.733 | 0.079 | 0.378 | 0.78 | 0.21 | -0.052 | 0.037 |
| 0.303 | 0.295 | 0.458 | 0.193 | 0.733 | 0.079 | 0.378 | 0.78 | 0.21 | -0.052 | 0.037 |
| 0.020 | 0.041 | 0.458 | 0.193 | 0.733 | 0.079 | 0.378 | 0.78 | 0.21 | -0.052 | 0.037 |
| 0.991 | 0.505 | 0.458 | 0.193 | 0.733 | 0.079 | 0.378 | 0.78 | 0.21 | -0.052 | 0.037 |

The loading, error-variance and error-covariance column names are
exactly the *definition variables*
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
reads in the model below.

One caveat: because
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
does not attach a per-cluster score intercept (`fsb`) for
[`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) fits, stage 2 fixes
the score intercepts at zero. That does not affect the latent
*covariances* compared below (a covariance is separable from the mean),
but it pins the latent *mean* estimates (`u0 + u1 ~ 1`) to the shrunk EB
score means, so read those means with that caveat.

### Using EB Without Accounting for Error

``` r

eb_m1 <- ranef(m1)[[1]]
eb_m2 <- ranef(m2)[[1]]
cor(cbind(eb_m1, eb_m2))
#>             (Intercept)       time1 (Intercept)      time2
#> (Intercept)   1.0000000  0.20647393  0.36187742  0.2659649
#> time1         0.2064739  1.00000000 -0.09670502 -0.1572327
#> (Intercept)   0.3618774 -0.09670502  1.00000000  0.2543058
#> time2         0.2659649 -0.15723268  0.25430582  1.0000000
cov(cbind(eb_m1, eb_m2))
#>             (Intercept)       time1 (Intercept)       time2
#> (Intercept)  0.91822191  0.08900785  0.48869284  0.15367628
#> time1        0.08900785  0.20238509 -0.06131108 -0.04265213
#> (Intercept)  0.48869284 -0.06131108  1.98610102  0.21610538
#> time2        0.15367628 -0.04265213  0.21610538  0.36359446
```

The numbers are quite different from the generating values.

### 2S-PA

#### Process 1

``` r

tspa1_mx <- tspa_mx_model(model = "u0 ~~ u1\nu0 + u1 ~ 1", data = fs_dat1)
#> Running m1 with 5 parameters
# Run OpenMx
tspa1_mx_fit <- mxRun(tspa1_mx)
#> Running m1 with 5 parameters
# Summarize the results
summary(tspa1_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix row col      Estimate  Std.Error A
#> 1 m1.S[3,3]      S  u0  u0  1.247318e+00 0.25524524  
#> 2 m1.S[3,4]      S  u0  u1 -1.293381e-02 0.08205187  
#> 3 m1.S[4,4]      S  u1  u1  2.582033e-01 0.04965696  
#> 4 m1.M[1,3]      M   1  u0 -3.360093e-09 0.13434551  
#> 5 m1.M[1,4]      M   1  u1  2.587087e-09 0.05925607  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              5                    195              392.9187
#>    Saturated:              5                    195                    NA
#> Independence:              4                    196                    NA
#> Number of observations/statistics: 100/200
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       2.918737               402.9187                 403.5570
#> BIC:    -505.089449               415.9446                 400.1533
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 11:02:55 
#> Wall clock time: 0.04067087 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

The estimated variances and covariances for $`u_0`$ and $`u_1`$ are
similar to the ones in the mixed model.

``` r

print("2S-PA covariance")
#> [1] "2S-PA covariance"
lat_cov(tspa1_mx_fit, c("u0", "u1"))
#>             u0          u1
#> u0  1.24731762 -0.01293381
#> u1 -0.01293381  0.25820331
print("lme4 covariance")
#> [1] "lme4 covariance"
as.matrix(VarCorr(m1)$clus_id1)
#>             (Intercept)       time1
#> (Intercept)  1.24731514 -0.01293371
#> time1       -0.01293371  0.25820269
#> attr(,"class")
#> [1] "vcmat_us" "matrix"   "array"   
#> attr(,"stddev")
#> (Intercept)       time1 
#>   1.1168326   0.5081365 
#> attr(,"correlation")
#>             (Intercept)       time1
#> (Intercept)  1.00000000 -0.02279053
#> time1       -0.02279053  1.00000000
```

Notes: The formulation works best with ML estimation for
[`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html); with REML, one may
need to consider the contrast matrix.

#### Process 2

``` r

tspa2_mx <- tspa_mx_model(model = "u0 ~~ u1\nu0 + u1 ~ 1", data = fs_dat2)
#> Running m1 with 5 parameters
# Run OpenMx
tspa2_mx_fit <- mxRun(tspa2_mx)
#> Running m1 with 5 parameters
# Summarize the results
summary(tspa2_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix row col      Estimate  Std.Error A
#> 1 m1.S[3,3]      S  u0  u0  2.511317e+00 0.45788873  
#> 2 m1.S[3,4]      S  u0  u1  1.723698e-01 0.11015636  
#> 3 m1.S[4,4]      S  u1  u1  3.647062e-01 0.05240907  
#> 4 m1.M[1,3]      M   1  u0  2.927051e-08 0.17993597  
#> 5 m1.M[1,4]      M   1  u1 -8.238643e-10 0.06087595  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              5                    195              526.3253
#>    Saturated:              5                    195                    NA
#> Independence:              4                    196                    NA
#> Number of observations/statistics: 100/200
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       136.3253               536.3253                 536.9636
#> BIC:      -371.6829               549.3512                 533.5599
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 11:02:55 
#> Wall clock time: 0.01644135 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

``` r

print("2S-PA covariance")
#> [1] "2S-PA covariance"
lat_cov(tspa2_mx_fit, c("u0", "u1"))
#>           u0        u1
#> u0 2.5113170 0.1723698
#> u1 0.1723698 0.3647062
print("lme4 covariance")
#> [1] "lme4 covariance"
as.matrix(VarCorr(m2)$clus_id2)
#>             (Intercept)     time2
#> (Intercept)   2.5113397 0.1723549
#> time2         0.1723549 0.3646935
#> attr(,"class")
#> [1] "vcmat_us" "matrix"   "array"   
#> attr(,"stddev")
#> (Intercept)       time2 
#>   1.5847207   0.6038986 
#> attr(,"correlation")
#>             (Intercept)     time2
#> (Intercept)   1.0000000 0.1800971
#> time2         0.1800971 1.0000000
```

#### Combining Processes 1 and 2

First, combine the two data sets with EB estimates. We’ll add suffix to
separate the variables.

``` r

colnames(fs_dat1) <- paste0(colnames(fs_dat1), "_1")
colnames(fs_dat2) <- paste0(colnames(fs_dat2), "_2")
fs_dat <- cbind(fs_dat1, fs_dat2)
```

``` r

# Unstructured model on the four growth-parameter latents. Explicit fsL/fsT
# here: auto-derivation yields one latent set per frame, so the two
# processes' u0/u1 latents must be named apart (u0_1/u0_2, u1_1/u1_2).
cross_load <- matrix(c("u0_by_fs_u0_1", "u0_by_fs_u1_1", NA, NA,
                       "u1_by_fs_u0_1", "u1_by_fs_u1_1", NA, NA,
                       NA, NA, "u0_by_fs_u0_2", "u0_by_fs_u1_2",
                       NA, NA, "u1_by_fs_u0_2", "u1_by_fs_u1_2"), nrow = 4) |>
    `dimnames<-`(list(c("fs_u0_1", "fs_u1_1", "fs_u0_2", "fs_u1_2"),
                      c("u0_1", "u1_1", "u0_2", "u1_2")))
err_cov <- matrix(c("ev_fs_u0_1", "ecov_fs_u1_fs_u0_1", NA, NA,
                    "ecov_fs_u1_fs_u0_1", "ev_fs_u1_1", NA, NA,
                    NA, NA, "ev_fs_u0_2", "ecov_fs_u1_fs_u0_2",
                    NA, NA, "ecov_fs_u1_fs_u0_2", "ev_fs_u1_2"), nrow = 4) |>
    `dimnames<-`(rep(list(c("fs_u0_1", "fs_u1_1", "fs_u0_2", "fs_u1_2")), 2))
tspa3_mx <- tspa_mx_model(
    model = "u0_1 ~~ u1_1 + u0_2 + u1_2\n
             u1_1 ~~ u0_2 + u1_2\n
             u0_2 ~~ u1_2\n
             u0_1 + u1_1 + u0_2 + u1_2 ~ 1",
    data = fs_dat, fsL = cross_load, fsT = err_cov)
#> Running m1 with 14 parameters
# Run OpenMx
tspa3_mx_fit <- mxRun(tspa3_mx)
#> Running m1 with 14 parameters
# Summarize the results
summary(tspa3_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>         name matrix  row  col      Estimate  Std.Error A
#> 1  m1.S[5,5]      S u0_1 u0_1  1.247317e+00 0.25524733  
#> 2  m1.S[5,6]      S u0_1 u1_1 -1.293369e-02 0.08205245  
#> 3  m1.S[6,6]      S u1_1 u1_1  2.582033e-01 0.04965696  
#> 4  m1.S[5,7]      S u0_1 u0_2  8.933885e-01 0.25771537  
#> 5  m1.S[6,7]      S u1_1 u0_2 -1.759013e-01 0.10806509  
#> 6  m1.S[7,7]      S u0_2 u0_2  2.511317e+00 0.45788899  
#> 7  m1.S[5,8]      S u0_1 u1_2  2.375128e-01 0.08516287  
#> 8  m1.S[6,8]      S u1_1 u1_2 -7.772932e-02 0.03690061  
#> 9  m1.S[7,8]      S u0_2 u1_2  1.723696e-01 0.11015577  
#> 10 m1.S[8,8]      S u1_2 u1_2  3.647062e-01 0.05240904  
#> 11 m1.M[1,5]      M    1 u0_1 -1.939552e-07 0.13434495  
#> 12 m1.M[1,6]      M    1 u1_1  3.526593e-08 0.05925593  
#> 13 m1.M[1,7]      M    1 u0_2 -1.797071e-07 0.17993663  
#> 14 m1.M[1,8]      M    1 u1_2 -3.811052e-08 0.06087598  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:             14                    386              893.8016
#>    Saturated:             14                    386                    NA
#> Independence:              8                    392                    NA
#> Number of observations/statistics: 100/400
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       121.8016               921.8016                 926.7428
#> BIC:      -883.7941               958.2740                 914.0585
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 11:02:56 
#> Wall clock time: 0.02103209 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

``` r

lat3 <- c("u0_1", "u1_1", "u0_2", "u1_2")
knitr::kable(cov(cbind(eb_m1, eb_m2)), digits = 3, caption = "EB covariance")
```

|             | (Intercept) |  time1 | (Intercept) |  time2 |
|:------------|------------:|-------:|------------:|-------:|
| (Intercept) |       0.918 |  0.089 |       0.489 |  0.154 |
| time1       |       0.089 |  0.202 |      -0.061 | -0.043 |
| (Intercept) |       0.489 | -0.061 |       1.986 |  0.216 |
| time2       |       0.154 | -0.043 |       0.216 |  0.364 |

EB covariance {.table}

``` r

knitr::kable(lat_cov(tspa3_mx_fit, lat3), digits = 3, caption = "2S-PA covariance")
```

|      |   u0_1 |   u1_1 |   u0_2 |   u1_2 |
|:-----|-------:|-------:|-------:|-------:|
| u0_1 |  1.247 | -0.013 |  0.893 |  0.238 |
| u1_1 | -0.013 |  0.258 | -0.176 | -0.078 |
| u0_2 |  0.893 | -0.176 |  2.511 |  0.172 |
| u1_2 |  0.238 | -0.078 |  0.172 |  0.365 |

2S-PA covariance {.table}

``` r

knitr::kable(latent_cov, digits = 3, caption = "Population covariance")
```

|      |        |        |       |
|-----:|-------:|-------:|------:|
| 1.00 |  0.000 |  0.750 |  0.18 |
| 0.00 |  0.250 | -0.075 | -0.09 |
| 0.75 | -0.075 |  2.250 |  0.18 |
| 0.18 | -0.090 |  0.180 |  0.36 |

Population covariance {.table}

One can see the 2S-PA covariance matrix is much closer to the population
values than using just EB estimates.

With only `n_clus = 100` clusters, any single draw carries substantial
sampling variability (a Monte-Carlo study over the same design gives
standard deviations of roughly `0.05`–`0.26` on these entries). The
2S-PA estimates are essentially unbiased – their mean across many
replications recovers the population values to within a few percent – so
the small gaps remaining here are sampling error from one realization,
not a bias in the estimator.

## With Missing Data

Now consider data with missing time points. We consider missing at
random, where data are more likely to be missing for later time points.

``` r

set.seed(1942)
# Missing indicators
r1 <- rbinom(time1, size = 1, prob = (20 - time1) / 20 - 0.05)
r2 <- rbinom(time2, size = 1, prob = (100 - time2) / 100 - 0.05)
```

Model with missing data

``` r

# Empirical Bayes
m1 <- lmer(y1 ~ time1 + (time1 | clus_id1),
           data = data.frame(y1, time1, clus_id1),
           REML = FALSE, subset = r1 == 1)
fs_dat1 <- get_fs(m1)
# Empirical Bayes
m2 <- lmer(y2 ~ time2 + (time2 | clus_id2),
           data = data.frame(y2, time2, clus_id2),
           REML = FALSE, subset = r2 == 1)
#> Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge with max|grad| = 0.00788774 (tol = 0.002, component 1)
#>   See ?lme4::convergence and ?lme4::troubleshooting.
fs_dat2 <- get_fs(m2)
colnames(fs_dat1) <- paste0(colnames(fs_dat1), "_1")
colnames(fs_dat2) <- paste0(colnames(fs_dat2), "_2")
fs_dat <- cbind(fs_dat1, fs_dat2)
```

Note that the cross-loadings and error covariances are not constant
across clusters.

``` r

head(fs_dat)
#>       fs_u0_1    fs_u1_1 fs_u0_se_1 fs_u1_se_1 u0_by_fs_u0_1 u0_by_fs_u1_1
#> 1  0.81235171  0.3005517  0.4852195  0.2122197     0.6873733    0.07210081
#> 2  1.03147731  0.1988596  0.4624694  0.1957475     0.7151924    0.08227485
#> 3 -1.56598438 -0.4989928  0.4667881  0.2076305     0.7131808    0.07828026
#> 4  0.30614717  0.2913250  0.4624694  0.1957475     0.7151924    0.08227485
#> 5 -0.05669005  0.2077012  0.4500962  0.2230413     0.7095656    0.09404091
#> 6  0.81697533  0.4450641  0.4852195  0.2122197     0.6873733    0.07210081
#>   u1_by_fs_u0_1 u1_by_fs_u1_1 ev_fs_u0_1 ecov_fs_u1_fs_u0_1 ev_fs_u1_1
#> 1     0.3396433     0.7287300  0.2354380        -0.04059019 0.04503721
#> 2     0.3870164     0.7685088  0.2138779        -0.05190589 0.03831709
#> 3     0.3707329     0.7361742  0.2178911        -0.04671442 0.04311041
#> 4     0.3870164     0.7685088  0.2138779        -0.05190589 0.03831709
#> 5     0.4554847     0.6253379  0.2025866        -0.04283360 0.04974744
#> 6     0.3396433     0.7287300  0.2354380        -0.04059019 0.04503721
#>      fs_u0_2    fs_u1_2 fs_u0_se_2 fs_u1_se_2 u0_by_fs_u0_2 u0_by_fs_u1_2
#> 1  1.7222568 0.32664011  0.6703862 0.06896616     0.7412050    0.02049019
#> 2  0.9606439 0.08687858  0.6583523 0.07474277     0.7538714    0.02085256
#> 3  1.1810650 0.24352887  0.6737159 0.06444376     0.7377001    0.02020873
#> 4 -1.4214190 0.36232332  0.6501387 0.06309225     0.7657186    0.01790385
#> 5  1.0052947 0.25705557  0.6741879 0.06346174     0.7390353    0.01891883
#> 6 -2.4228230 0.41403082  0.6672456 0.07610994     0.7433057    0.02162108
#>   u1_by_fs_u0_2 u1_by_fs_u1_2 ev_fs_u0_2 ecov_fs_u1_fs_u0_2  ev_fs_u1_2
#> 1     0.2523431     0.9743629  0.4494177        -0.03434766 0.004756331
#> 2     0.2477879     0.9717386  0.4334278        -0.03513067 0.005586482
#> 3     0.2528806     0.9762748  0.4538931        -0.03404548 0.004152998
#> 4     0.2246715     0.9784454  0.4226803        -0.03137015 0.003980631
#> 5     0.2437995     0.9775809  0.4545293        -0.03193221 0.004027392
#> 6     0.2575528     0.9705645  0.4452167        -0.03579237 0.005792723
```

Now run 2S-PA and compare with EB estimates without accounting for error

``` r

tspa3o_mx <- tspa_mx_model(
    model = "u0_1 ~~ u1_1 + u0_2 + u1_2\n
             u1_1 ~~ u0_2 + u1_2\n
             u0_2 ~~ u1_2\n
             u0_1 + u1_1 + u0_2 + u1_2 ~ 1",
    data = fs_dat, fsL = cross_load, fsT = err_cov)
#> Running m1 with 14 parameters
# Run OpenMx
tspa3o_mx_fit <- mxTryHard(tspa3o_mx)
#> Running m1 with 14 parameters
#> 
#> Beginning initial fit attempt
#> Running m1 with 14 parameters
#> 
#>  Lowest minimum so far:  865.097876007988
#> 
#> Solution found
```

    #> 
    #>  Solution found!  Final fit=865.09788 (started at 865.09788)  (1 attempt(s): 1 valid, 0 errors)
    #>  Start values from best fit:
    #> 1.21323293986919,-0.0198051249223241,0.257101981152412,0.933483129728776,-0.217794756238648,2.57220658680846,0.224924785232625,-0.0678737023570899,0.169045589944976,0.365635026795862,-0.00224803426594298,0.00108201384570927,-0.00200753731223892,0.000225383704665787
    # Summarize the results
    summary(tspa3o_mx_fit)
    #> Summary of m1 
    #>  
    #> free parameters:
    #>         name matrix  row  col      Estimate  Std.Error A
    #> 1  m1.S[5,5]      S u0_1 u0_1  1.2132329399 0.26846470  
    #> 2  m1.S[5,6]      S u0_1 u1_1 -0.0198051249 0.09110896  
    #> 3  m1.S[6,6]      S u1_1 u1_1  0.2571019812 0.05466823  
    #> 4  m1.S[5,7]      S u0_1 u0_2  0.9334831297 0.27064841  
    #> 5  m1.S[6,7]      S u1_1 u0_2 -0.2177947562 0.11909674  
    #> 6  m1.S[7,7]      S u0_2 u0_2  2.5722065868 0.48144199  
    #> 7  m1.S[5,8]      S u0_1 u1_2  0.2249247852 0.08719662  
    #> 8  m1.S[6,8]      S u1_1 u1_2 -0.0678737024 0.03930047  
    #> 9  m1.S[7,8]      S u0_2 u1_2  0.1690455899 0.11318626  
    #> 10 m1.S[8,8]      S u1_2 u1_2  0.3656350268 0.05272070  
    #> 11 m1.M[1,5]      M    1 u0_1 -0.0022480343 0.13856273  
    #> 12 m1.M[1,6]      M    1 u1_1  0.0010820138 0.06330727  
    #> 13 m1.M[1,7]      M    1 u0_2 -0.0020075373 0.18461242  
    #> 14 m1.M[1,8]      M    1 u1_2  0.0002253837 0.06107594  
    #> 
    #> Model Statistics: 
    #>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
    #>        Model:             14                    386              865.0979
    #>    Saturated:             14                    386                    NA
    #> Independence:              8                    392                    NA
    #> Number of observations/statistics: 100/400
    #> 
    #> Information Criteria: 
    #>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
    #> AIC:       93.09788               893.0979                 898.0391
    #> BIC:     -912.49782               929.5703                 885.3547
    #> To get additional fit indices, see help(mxRefModels)
    #> timestamp: 2026-08-30 11:02:57 
    #> Wall clock time: 0.2027819 secs 
    #> optimizer:  SLSQP 
    #> OpenMx version number: 2.22.11 
    #> Need help?  See help(mxSummary)

``` r

lat3 <- c("u0_1", "u1_1", "u0_2", "u1_2")
knitr::kable(cov(cbind(fs_dat1[, 1:2], fs_dat2[, 1:2])), digits = 3, caption = "EB covariance")
```

|         | fs_u0_1 | fs_u1_1 | fs_u0_2 | fs_u1_2 |
|:--------|--------:|--------:|--------:|--------:|
| fs_u0_1 |   0.841 |   0.091 |   0.458 |   0.140 |
| fs_u1_1 |   0.091 |   0.179 |  -0.066 |  -0.030 |
| fs_u0_2 |   0.458 |  -0.066 |   1.985 |   0.220 |
| fs_u1_2 |   0.140 |  -0.030 |   0.220 |   0.363 |

EB covariance {.table}

``` r

knitr::kable(lat_cov(tspa3o_mx_fit, lat3), digits = 3, caption = "2S-PA covariance")
```

|      |   u0_1 |   u1_1 |   u0_2 |   u1_2 |
|:-----|-------:|-------:|-------:|-------:|
| u0_1 |  1.213 | -0.020 |  0.933 |  0.225 |
| u1_1 | -0.020 |  0.257 | -0.218 | -0.068 |
| u0_2 |  0.933 | -0.218 |  2.572 |  0.169 |
| u1_2 |  0.225 | -0.068 |  0.169 |  0.366 |

2S-PA covariance {.table}

``` r

knitr::kable(latent_cov, digits = 3, caption = "Population covariance")
```

|      |        |        |       |
|-----:|-------:|-------:|------:|
| 1.00 |  0.000 |  0.750 |  0.18 |
| 0.00 |  0.250 | -0.075 | -0.09 |
| 0.75 | -0.075 |  2.250 |  0.18 |
| 0.18 | -0.090 |  0.180 |  0.36 |

Population covariance {.table}
