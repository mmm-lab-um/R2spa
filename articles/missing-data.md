# 2S-PA with Missing Data

With missing data, two things change in two-stage path analysis. Stage 1
must be fitted with a missing-data method (`missing = "fiml"` in
`lavaan`), and each case is then scored on its *own* observed
indicators, so the score standard errors, implied loadings, and error
variances are **per-row** quantities rather than one pooled matrix. And
stage 2 has to use those per-row values:
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) (the
`lavaan` backend) reduces the per-pattern values to a single
representative matrix (`reduce =`), whereas
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
(the `OpenMx` backend) carries the per-row values through definition
variables, and its raw-data FIML evaluation can also handle rows whose
scores are dropped as `NA`. The stage-2 fits in this vignette therefore
use
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md).

For illustration, we create two missing data patterns on
`PoliticalDemocracy`: (1) MCAR missingness in a single item (`x1`, about
30% of cases); (2) block missingness, where the last 10 cases are
missing all four `dem65` indicators — a factor with no observed
indicator at all.

## Single-Factor Model

``` r

set.seed(2225)
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(R2spa)
library(OpenMx)
#> To take full advantage of multiple cores, use:
#>   mxOption(key='Number of Threads', value=parallel::detectCores()) #now
#>   Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #before library(OpenMx)
data(PoliticalDemocracy)
pd2 <- PoliticalDemocracy
# Add MCAR missing data to x1 (about 30% of the 75 cases)
pd2[!rbinom(75, size = 1, prob = 0.7), "x1"] <- NA
```

Factor scores with
[`lavaan::lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html):

``` r

fit <- cfa("ind60 =~ x1 + x2 + x3", data = pd2, missing = "fiml")
fs_lavaan <- lavPredict(fit, method = "Bartlett", se = TRUE, fsm = TRUE)
# From R2spa
fs <- R2spa::get_fs(fit, method = "Bartlett")
# One SE per missing-data pattern
round(table(pattern = attr(fs, "fs_pattern")[[1]]$label,
            se = fs$fs_ind60_se), 4)
#>           se
#> pattern    0.149802246676676 0.179196081489146
#>   x1+x2+x3                48                 0
#>   x2+x3                    0                27
# Compare factor scores
plot(fs_lavaan, fs$fs_ind60)
abline(a = 0, b = 1)
```

![](missing-data_files/figure-html/unnamed-chunk-2-1.png)

Because each case is scored conditionally on its own observed
indicators,
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
returns per-row standard errors and error variances: the 48 complete
cases (pattern `"x1+x2+x3"`) get SE 0.150, and the 27 cases observed
only on `x2` and `x3` get SE 0.179. The
`fsT`/`fsL`/`fsb`/`scoring_matrix` attributes are correspondingly
*lists* with one entry per pattern, and the `fs_pattern` attribute maps
each row to its pattern. The plot confirms that the `R2spa` scores
reproduce
[`lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html)
exactly.

### Multiple group

``` r

mg_fit <- cfa("ind60 =~ x1 + x2 + x3",
  data = cbind(pd2, group = rep(1:2, c(40, 35))),
  missing = "fiml",
  group = "group"
)
mg_fs <- R2spa::get_fs(mg_fit)
head(mg_fs, 10)
#>        fs_ind60 fs_ind60_se ind60_by_fs_ind60 ev_fs_ind60 group
#> 1  -0.629900908  0.06327411         0.9907469 0.004003613     1
#> 2  -0.003598338  0.06481599         0.9902860 0.004201113     1
#> 3   0.531624709  0.06327411         0.9907469 0.004003613     1
#> 4   1.093155269  0.06481599         0.9902860 0.004201113     1
#> 5   0.751323885  0.06481599         0.9902860 0.004201113     1
#> 6   0.032310438  0.06481599         0.9902860 0.004201113     1
#> 7  -0.001195823  0.06327411         0.9907469 0.004003613     1
#> 8  -0.081264245  0.06481599         0.9902860 0.004201113     1
#> 9   0.087672535  0.06327411         0.9907469 0.004003613     1
#> 10  0.158894033  0.06327411         0.9907469 0.004003613     1
# Both groups contain the same two patterns
lapply(attr(mg_fs, "fs_pattern"), function(p) table(p$label))
#> $`1`
#> 
#> x1+x2+x3    x2+x3 
#>       25       15 
#> 
#> $`2`
#> 
#> x1+x2+x3    x2+x3 
#>       23       12
```

The multigroup result is a single data frame with a `group` column, and
each group’s `fsT`/`fsL`/`fsb`/`scoring_matrix` attribute is itself a
per-pattern list.

## Multiple-Factor Model

The measurement model constrains the two factors’ loadings to be equal
(labels `a`–`d`) and adds covariances between the corresponding
indicators:

``` r

# Make last 10 cases completely missing on y5-y8 (all dem65 indicators)
pd2[66:75, c("y5", "y6", "y7", "y8")] <- NA
fit2 <- cfa("
  dem60 =~ a * y1 + b * y2 + c * y3 + d * y4
  dem65 =~ a * y5 + b * y6 + c * y7 + d * y8
  y1 ~~ y5
  y2 ~~ y6
  y3 ~~ y7
  y4 ~~ y8
", data = pd2, missing = "fiml")
```

This block-missing pattern is the hard case: for the 10 affected rows,
`dem65` has no observed indicator at all.

### Regression scores

With regression (EBM) scores, the `dem65` score is still *defined* for
those rows — regression scoring uses the cross-factor covariances, so
`dem65` is imputed from the `dem60` indicators. But the per-row
measurement information is degenerate: the implied loading of the score
on `dem65` is 0, i.e. the score carries no information about `dem65`, so
measurement-error correction is impossible. We therefore drop those
`dem65` scores; OpenMx’s raw-data FIML then fits those rows from the
structural model alone.

``` r

# Obtain tidy-ed factor scores data
fs_dat <- R2spa::augment_lav_predict(fit2)
# Per-row definition-variable matrices: each cell names the `fs_dat` column
# carrying that case's loadings / error variances (the exact, per-row
# correction that `tspa()`'s `reduce = ` pooling only approximates).
cross_load <- matrix(c(
  "dem60_by_fs_dem60", "dem60_by_fs_dem65",
  "dem65_by_fs_dem60", "dem65_by_fs_dem65"
), nrow = 2) |>
  `dimnames<-`(list(c("fs_dem60", "fs_dem65"), c("dem60", "dem65")))
err_cov <- matrix(c(
  "ev_fs_dem60", "ecov_fs_dem60_fs_dem65",
  "ecov_fs_dem60_fs_dem65", "ev_fs_dem65"
), nrow = 2) |>
  `dimnames<-`(rep(list(c("fs_dem60", "fs_dem65")), 2))
# Drop the non-correctable dem65 scores of the 10 all-missing rows;
# OpenMx's raw-data FIML then fits those rows as usual.
fs_dat[66:75, "fs_dem65"] <- NA
tspa_mx <- tspa_mx_model(
  "dem65 ~ dem60; dem65 + dem60 ~ 1",
  data = fs_dat,
  fsL = cross_load, fsT = err_cov,
  fsb = c(fs_dem60 = "int_fs_dem60", fs_dem65 = "int_fs_dem65")
)
#> Running m1 with 5 parameters
# Run OpenMx
tspa_mx_fit <- mxRun(tspa_mx)
#> Running m1 with 5 parameters
# Summarize the results
summary(tspa_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix   row   col      Estimate Std.Error A
#> 1 m1.A[4,3]      A dem65 dem60  9.125862e-01 0.0715752  
#> 2 m1.S[3,3]      S dem60 dem60  4.559668e+00 0.8414167  
#> 3 m1.S[4,4]      S dem65 dem65  5.369811e-01 0.2310376  
#> 4 m1.M[1,3]      M     1 dem60 -4.568894e-08 0.2621047  
#> 5 m1.M[1,4]      M     1 dem65 -5.709867e-08 0.1421203  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              5                    135              404.9387
#>    Saturated:              5                    135                    NA
#> Independence:              4                    136                    NA
#> Number of observations/statistics: 75/140
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       134.9387               414.9387                 415.8083
#> BIC:      -177.9222               426.5261                 410.7675
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 11:02:50 
#> Wall clock time: 0.0296402 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

The per-row measurement values reach OpenMx as *definition variables*:
`fsL` and `fsT` are character matrices whose cells name the `fs_dat`
column holding that case’s value, so the stage-2 correction is per-row
rather than pooled. `fsb` does the same for the per-row measurement
intercepts (the `int_*` columns; a FIML fit estimates a mean structure,
so the intercepts vary with the missing-data pattern).

### Bartlett scores

Bartlett scoring does not use the cross-factor covariances, so for the
10 rows with no observed `dem65` indicator the score is `NA` by
convention — in contrast to the EBM case above, where the score was
defined but non-correctable and had to be dropped.

``` r

# Obtain tidy-ed factor scores data (Bartlett weights)
fsb_dat <- augment_lav_predict(fit2, method = "Bartlett")
# fsb_dat[66:75, "fs_dem65"] is already NA (no observed indicator);
# stated explicitly to mirror the EBM case
fsb_dat[66:75, "fs_dem65"] <- NA
tspab_mx <- tspa_mx_model(
  "dem65 ~ dem60; dem65 + dem60 ~ 1",
  data = fsb_dat,
  fsL = cross_load, fsT = err_cov,
  fsb = c(fs_dem60 = "int_fs_dem60", fs_dem65 = "int_fs_dem65")
)
#> Running m1 with 5 parameters
# Run OpenMx
tspab_mx_fit <- mxRun(tspab_mx)
#> Running m1 with 5 parameters
# Summarize the results
summary(tspab_mx_fit)
#> Summary of m1 
#>  
#> free parameters:
#>        name matrix   row   col      Estimate  Std.Error A
#> 1 m1.A[4,3]      A dem65 dem60  9.125862e-01 0.07157514  
#> 2 m1.S[3,3]      S dem60 dem60  4.559668e+00 0.84138958  
#> 3 m1.S[4,4]      S dem65 dem65  5.369811e-01 0.23103931  
#> 4 m1.M[1,3]      M     1 dem60 -4.578959e-08 0.26210338  
#> 5 m1.M[1,4]      M     1 dem65 -5.290055e-08 0.14212079  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              5                    135              536.4213
#>    Saturated:              5                    135                    NA
#> Independence:              4                    136                    NA
#> Number of observations/statistics: 75/140
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:      266.42127               546.4213                 547.2908
#> BIC:      -46.43962               558.0087                 542.2500
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-08-30 11:02:50 
#> Wall clock time: 0.01694989 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

### Compared to Joint Model

The gold-standard reference is the joint model, which fits the
measurement and structural parts simultaneously under FIML:

``` r

jfit <- sem("
  dem60 =~ a * y1 + b * y2 + c * y3 + d * y4
  dem65 =~ a * y5 + b * y6 + c * y7 + d * y8
  y1 ~~ y5
  y2 ~~ y6
  y3 ~~ y7
  y4 ~~ y8
  dem65 ~ dem60
", data = pd2, missing = "fiml")
summary(jfit)
#> lavaan 0.7-2 ended normally after 53 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        29
#>   Number of equality constraints                     3
#> 
#>   Number of observations                            75
#>   Number of missing patterns                         2
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                27.499
#>   Degrees of freedom                                18
#>   P-value (Chi-square)                           0.070
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Observed
#>   Observed information based on                Hessian
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem60 =~                                            
#>     y1         (a)    1.000                           
#>     y2         (b)    1.337    0.153    8.738    0.000
#>     y3         (c)    1.182    0.132    8.942    0.000
#>     y4         (d)    1.334    0.136    9.777    0.000
#>   dem65 =~                                            
#>     y5         (a)    1.000                           
#>     y6         (b)    1.337    0.153    8.738    0.000
#>     y7         (c)    1.182    0.132    8.942    0.000
#>     y8         (d)    1.334    0.136    9.777    0.000
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem65 ~                                             
#>     dem60             0.913    0.074   12.357    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>  .y1 ~~                                               
#>    .y5                0.627    0.395    1.587    0.113
#>  .y2 ~~                                               
#>    .y6                1.238    0.765    1.618    0.106
#>  .y3 ~~                                               
#>    .y7                1.534    0.695    2.207    0.027
#>  .y4 ~~                                               
#>    .y8                0.200    0.551    0.363    0.716
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .y1                5.465    0.298   18.348    0.000
#>    .y2                4.256    0.442    9.629    0.000
#>    .y3                6.563    0.397   16.534    0.000
#>    .y4                4.453    0.379   11.733    0.000
#>    .y5                5.181    0.313   16.561    0.000
#>    .y6                2.782    0.403    6.899    0.000
#>    .y7                6.249    0.367   17.043    0.000
#>    .y8                3.996    0.383   10.420    0.000
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .y1                2.093    0.490    4.276    0.000
#>    .y2                6.502    1.252    5.195    0.000
#>    .y3                5.443    1.038    5.242    0.000
#>    .y4                2.689    0.716    3.756    0.000
#>    .y5                2.528    0.551    4.588    0.000
#>    .y6                3.677    0.827    4.446    0.000
#>    .y7                3.386    0.760    4.453    0.000
#>    .y8                2.628    0.702    3.745    0.000
#>     dem60             4.560    1.038    4.394    0.000
#>    .dem65             0.537    0.258    2.085    0.037
```

``` r

# The three estimates of the dem65 ~ dem60 path:
c(EBM.mx = coef(tspa_mx_fit)["m1.A[4,3]"],
  Bartlett.mx = coef(tspab_mx_fit)["m1.A[4,3]"],
  joint = coef(jfit)["dem65~dem60"])
#>      EBM.mx.m1.A[4,3] Bartlett.mx.m1.A[4,3]     joint.dem65~dem60 
#>             0.9125862             0.9125862             0.9125862
```

Both 2S-PA fits (regression and Bartlett scores) recover the
`dem65 ~ dem60` path at 0.913, identical to the joint-model estimate to
four decimal places.

### Multiple group

The same pipeline applies to multigroup fits:
[`augment_lav_predict()`](https://mmm-lab-um.github.io/R2spa/reference/augment_lav_predict.md)
returns a list with one tidy-ed data frame per group, each carrying that
group’s per-row values. All 10 block-missing rows fall in group 2, so
only group 2 contains the pattern with no observed `dem65` indicator.
(`lavaan` warns that a single label used in a multigroup fit imposes
cross-group equality — that is intentional here: `a`–`d` are constrained
equal in both groups.)

``` r

mg_fit2 <- cfa("
  dem60 =~ a * y1 + b * y2 + c * y3 + d * y4
  dem65 =~ a * y5 + b * y6 + c * y7 + d * y8
  y1 ~~ y5
  y2 ~~ y6
  y3 ~~ y7
  y4 ~~ y8
",
  data = cbind(pd2, group = rep(1:2, c(40, 35))),
  missing = "fiml", group = "group"
)
#> Warning: lavaan->lav_model_pt():  
#>    using a single label per parameter in a multiple group setting implies 
#>    imposing equality constraints across all the groups; If this is not 
#>    intended, either remove the label(s), or use a vector of labels (one for 
#>    each group); See the Multiple groups section in the man page of 
#>    model.syntax.
mg_out <- get_fs(mg_fit2)
head(mg_out)
#>     fs_dem60   fs_dem65 fs_dem60_se fs_dem65_se dem60_by_fs_dem60
#> 1 -2.8000738 -2.3481803   0.4602081    0.484165         0.7056551
#> 2 -3.0162173 -1.6560141   0.4602081    0.484165         0.7056551
#> 3  2.6733722  3.0172082   0.4602081    0.484165         0.7056551
#> 4  2.4571381  2.4920427   0.4602081    0.484165         0.7056551
#> 5  1.3378821  1.3982512   0.4602081    0.484165         0.7056551
#> 6  0.1326308 -0.4178734   0.4602081    0.484165         0.7056551
#>   dem60_by_fs_dem65 dem65_by_fs_dem60 dem65_by_fs_dem65 ev_fs_dem60
#> 1         0.2749477         0.2331707          0.681784   0.2117915
#> 2         0.2749477         0.2331707          0.681784   0.2117915
#> 3         0.2749477         0.2331707          0.681784   0.2117915
#> 4         0.2749477         0.2331707          0.681784   0.2117915
#> 5         0.2749477         0.2331707          0.681784   0.2117915
#> 6         0.2749477         0.2331707          0.681784   0.2117915
#>   ecov_fs_dem65_fs_dem60 ev_fs_dem65 group
#> 1              0.1496662   0.2344158     1
#> 2              0.1496662   0.2344158     1
#> 3              0.1496662   0.2344158     1
#> 4              0.1496662   0.2344158     1
#> 5              0.1496662   0.2344158     1
#> 6              0.1496662   0.2344158     1
```
