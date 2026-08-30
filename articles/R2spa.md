# Two-Stage Path Analysis (2S-PA) Model Examples

The example is from <https://lavaan.ugent.be/tutorial/sem.html>.

Load packages

``` r

library(lavaan)
library(R2spa)
```

## Why two-stage path analysis

When the constructs in a path model are measured with error, plugging
composite or factor scores into a regression as if they were error-free
biases the path coefficients (Cole & Preacher, 2014). A single joint
structural equation model, estimating the measurement and structural
models together, is the gold standard, but it can be hard to specify, is
prone to convergence and admissibility problems in small samples, and
forces a whole-model refit whenever the measurement *or* the structural
part changes (Rosseel, 2012).

Two-stage path analysis (2S-PA; Lai & Hsiao, 2022) is a practical
alternative that splits the problem:

- **Stage 1 — score each construct** with an appropriate psychometric
  model (a `lavaan` CFA, an `lme4` mixed model, a `mirt` IRT model, or
  even an EFA) and record each score’s measurement error / reliability.
- **Stage 2 — path-analyze the scores**, correcting for measurement
  error by fixing the loading(s) and error variance (and covariances) to
  values implied in Stage 1. Observation-specific values can be used via
  *definition variables* (currently available in *OpenMx*).

Because the stages are separate, you can tune and diagnose each on its
own, mix estimators across constructs (one construct from a CFA, another
from an IRT model), and avoid the small-sample convergence problems of
one large joint model. Simulations find 2S-PA tracks the joint model in
large samples but converges more often, with less standard-error bias
and better Type-I error / CI coverage in small samples and with
categorical items (Lai & Hsiao, 2022; Lai et al., 2023); the same
two-stage idea also covers empirical-Bayes random slopes from mixed
models (Lai & Liu, 2026).

`R2spa` implements both stages:
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) is
stage 1 (scores plus their loadings, error covariances and intercepts)
and [`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) /
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
are stage 2.

## How the measurement error is supplied

[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) needs
each score’s measurement error, in one of two forms:

- **`se_fs`** — one standard error per score (a scalar, or a list / data
  frame of per-group scalars for multigroup). This treats the error
  *variance* as constant across observations and assumes each score is
  an indicator of a single latent (no cross-loadings, no correlated
  score errors). It is the simple path: suitable for a unidimensional
  construct scored from continuous, roughly equally reliable items.
- **`fsT`** (with **`fsL`** and **`fsb`**) — the general form
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  returns as attributes (`attr(fs, "fsT")`, `attr(fs, "fsL")`,
  `attr(fs, "fsb")`). `fsT` is the score error *variance-covariance* (so
  it also carries correlated errors and cross-loadings) and `fsL` the
  implied loading of each latent on the score. Use this when scores are
  composite (several indicators per latent, or a score loading on
  several latents), and note it is the only form that carries
  per-observation measurement properties, which the `mirt` (IRT) and
  `merMod` (mixed model) paths need because loadings and error variances
  vary by person / cluster (see the `multiple-factors`,
  `tspa-vignette-mx` and `multilevel` vignettes).

Each construct below is a single-factor score from continuous items, so
`se_fs` is the natural form: the single-factor example passes it
explicitly, while the multi-factor and multigroup examples omit it,
since [`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
reads those standard errors from the
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result automatically. The `multiple-factors` vignette shows the
`fsT`/`fsL` form. The `OpenMx` backend
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
likewise omits the measurement inputs on a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result, deriving them from the result’s attributes (per-row quantities
as definition variables — the exact, non-pooled route; see the
`tspa-vignette-mx` vignette).

## Reading the `get_fs()` output

[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
returns one row per observation. For a construct called `<name>` it
contributes these columns, with the matching matrix-form quantities
attached as attributes on the returned data frame:

| Column(s) / attribute | Meaning |
|----|----|
| `fs_<name>` | the factor score |
| `fs_<name>_se` | that score’s standard error — the measurement-error term [`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) corrects for |
| `ev_fs_<name>` | the score’s error variance |
| `ecov_fs_<a>_fs_<b>` | error covariance between scores `a` and `b` (only from a shared multi-factor model) |
| `<indicator>_by_fs_<name>` | model-implied loading — the mean of the indicator on the score |
| `attr(fs, "fsT")` | the score error variance-covariance (carries correlated errors and cross-loadings) |
| `attr(fs, "fsL")` | the implied loading of each latent on the score |
| `attr(fs, "fsb")` | the scoring-equation intercept, one per score (named by the score, e.g. `fs_ind60`; `0` for a complete-data single-factor CFA, varying per row when the measurement quantities do) |

The attributes are what
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
consumes, and `tspa(fsb = ...)` uses `fsb` to fix the stage-2 intercept
constraints. `fs_indiv(include_intercept = TRUE)` materializes the `fsb`
values as `int_fs_<f>` columns, which is how the
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
examples wire up the score-intercept definition variables
(`fsb = c(fs_ind60 = "int_fs_ind60", ...)`).

## Mean (sum) scores

[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) can
also score a construct as its **mean score** — the plain mean of the
items that load on it, with no latent distribution — a simpler
alternative to the default regression (EB) scores. Request it with
`method = "mean"` (the item-to-factor split is auto-derived from the
estimated loadings, or supplied via `sum_items`):

``` r

fs_mean <- get_fs(PoliticalDemocracy,
                  model = "ind60 =~ x1 + x2 + x3
                           dem60 =~ y1 + y2 + y3 + y4",
                  method = "mean")
head(fs_mean)   # fs_ind60 == rowMeans(x1, x2, x3); constant fs_<name>_se per factor
```

    ##   fs_ind60 fs_dem60 fs_ind60_se fs_dem60_se ind60_by_fs_ind60 ind60_by_fs_dem60
    ## 1 3.545951 1.458333   0.2725027    1.014152          1.666078                 0
    ## 2 4.671723 1.145833   0.2725027    1.014152          1.666078                 0
    ## 3 5.813729 8.874997   0.2725027    1.014152          1.666078                 0
    ## 4 6.707119 9.224997   0.2725027    1.014152          1.666078                 0
    ## 5 5.752078 7.499999   0.2725027    1.014152          1.666078                 0
    ## 6 4.853819 6.041666   0.2725027    1.014152          1.666078                 0
    ##   dem60_by_fs_ind60 dem60_by_fs_dem60 ev_fs_ind60 ecov_fs_dem60_fs_ind60
    ## 1                 0          1.235413  0.07425774                      0
    ## 2                 0          1.235413  0.07425774                      0
    ## 3                 0          1.235413  0.07425774                      0
    ## 4                 0          1.235413  0.07425774                      0
    ## 5                 0          1.235413  0.07425774                      0
    ## 6                 0          1.235413  0.07425774                      0
    ##   ev_fs_dem60
    ## 1    1.028505
    ## 2    1.028505
    ## 3    1.028505
    ## 4    1.028505
    ## 5    1.028505
    ## 6    1.028505

A mean score is a composite of items, so its measurement error is
carried by the `fsT`/`fsL`/`fsb` form (the cross terms are exactly zero
for a plain CFA).
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) reads
those from the result’s attributes:

``` r

tspa_mean_fit <- tspa(model = "dem60 ~ ind60",
                      data = fs_mean,
                      fsT = attr(fs_mean, "fsT"),
                      fsL = attr(fs_mean, "fsL"),
                      fsb = attr(fs_mean, "fsb"))
summary(tspa_mean_fit, standardized = TRUE)
```

    ## lavaan 0.7-2 ended normally after 22 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                         3
    ## 
    ##   Number of observations                            75
    ## 
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                                 0.000
    ##   Degrees of freedom                                 2
    ##   P-value (Chi-square)                           1.000
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Expected
    ##   Information saturated (h1) model          Structured
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
    ##   ind60 =~                                                              
    ##     fs_ind60          1.666                               1.116    0.971
    ##     fs_dem60          0.000                               0.000    0.000
    ##   dem60 =~                                                              
    ##     fs_ind60          0.000                               0.000    0.000
    ##     fs_dem60          1.235                               2.583    0.931
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
    ##   dem60 ~                                                               
    ##     ind60             1.319    0.369    3.577    0.000    0.422    0.422
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
    ##  .fs_ind60 ~~                                                           
    ##    .fs_dem60          0.000                               0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
    ##    .fs_ind60          4.468                               4.468    3.890
    ##    .fs_dem60          5.184                               5.184    1.868
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
    ##    .fs_ind60          0.074                               0.074    0.056
    ##    .fs_dem60          1.029                               1.029    0.134
    ##     ind60             0.449    0.078    5.779    0.000    1.000    1.000
    ##    .dem60             3.591    0.704    5.101    0.000    0.822    0.822

`method = "mean"` is for complete-data fits: it errors when stage 1
retains missing data (e.g. `missing = "fiml"`) and is not supported
together with `corrected_fsT`, `vfsLT`, `reliability`, `prior_mean`, or
`prior_cov`.

## Single group, single factor

``` r

model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4

  # regressions
    dem60 ~ ind60
'
```

To call
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md), a data
frame of factor scores is needed for all latent variables. To get this
data frame, apply
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) to
all latent variables and specify model parameters as their respective
definitions. Combine factor scores for all latent variables using
[`cbind()`](https://rdrr.io/r/base/cbind.html) so that it can be used in
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) for
model building.

``` r

fs_dat_ind60 <- get_fs(PoliticalDemocracy, 
                       model = "ind60 =~ x1 + x2 + x3")
fs_dat_dem60 <- get_fs(PoliticalDemocracy, 
                       model = "dem60 =~ y1 + y2 + y3 + y4")
fs_dat <- cbind(fs_dat_ind60, fs_dat_dem60)
head(fs_dat)
```

    ##     fs_ind60 fs_ind60_se ind60_by_fs_ind60 ev_fs_ind60   fs_dem60 fs_dem60_se
    ## 1 -0.5261683   0.1213615         0.9657673  0.01472862 -2.7487224   0.6756472
    ## 2  0.1436527   0.1213615         0.9657673  0.01472862 -3.0360803   0.6756472
    ## 3  0.7143559   0.1213615         0.9657673  0.01472862  2.6718589   0.6756472
    ## 4  1.2399257   0.1213615         0.9657673  0.01472862  2.9936997   0.6756472
    ## 5  0.8319080   0.1213615         0.9657673  0.01472862  1.9242932   0.6756472
    ## 6  0.2123845   0.1213615         0.9657673  0.01472862  0.9922798   0.6756472
    ##   dem60_by_fs_dem60 ev_fs_dem60
    ## 1         0.8868049   0.4564991
    ## 2         0.8868049   0.4564991
    ## 3         0.8868049   0.4564991
    ## 4         0.8868049   0.4564991
    ## 5         0.8868049   0.4564991
    ## 6         0.8868049   0.4564991

To build a Two-Stage Path Analysis model, simply call the
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
function with model = regressions, data = the combined factor score data
frame using
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md),
and specify standard error as either a list or a data frame. Values for
standard error can be found at column named `fs_[variable name]_se`. For
example, the standard error of latent variable ind60 can be found at
column `fs_ind60_se` in the `fs_dat` data frame.

``` r

tspa_fit <- tspa(model = "dem60 ~ ind60", 
                 data = fs_dat, 
                 se_fs = list(ind60 = 0.1213615, dem60 = 0.6756472))
```

We recently added support for
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) to
directly grep the `fs_[variable name]_se` columns from the
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
output, so the `se_fs` argument can be omitted:

``` r

tspa_fit <- tspa(model = "dem60 ~ ind60", 
                 data = fs_dat)
summary(tspa_fit, standardized = TRUE)
```

    ## lavaan 0.7-2 ended normally after 17 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                         3
    ## 
    ##   Number of observations                            75
    ## 
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                                 0.000
    ##   Degrees of freedom                                 0
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Expected
    ##   Information saturated (h1) model          Structured
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
    ##   ind60 =~                                                              
    ##     fs_ind60          1.000                               0.645    0.983
    ##   dem60 =~                                                              
    ##     fs_dem60          1.000                               1.891    0.942
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
    ##   dem60 ~                                                               
    ##     ind60             1.329    0.332    4.000    0.000    0.453    0.453
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
    ##    .fs_ind60          0.015                               0.015    0.034
    ##    .fs_dem60          0.456                               0.456    0.113
    ##     ind60             0.416    0.070    5.914    0.000    1.000    1.000
    ##    .dem60             2.842    0.543    5.235    0.000    0.795    0.795

To view the Two-Stage Path Analysis model, use
`attributes([model name])$tspaModel`. Function
[`cat()`](https://rdrr.io/r/base/cat.html) can help tidy up the model
output. In the output, the values of error constraints are computed by
squaring standard errors from the previous section.

``` r

cat(attributes(tspa_fit)$tspaModel)
```

    ## # latent variables (indicated by factor scores)
    ## ind60 =~ 1 * fs_ind60
    ## dem60 =~ 1 * fs_dem60
    ## 
    ## # constrain the errors
    ## fs_ind60 ~~ 0.0147286194470875 * fs_ind60
    ## fs_dem60 ~~ 0.456499146148539 * fs_dem60
    ## 
    ## # structural model
    ## dem60 ~ ind60

## Single group, multiple factors

``` r

model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # # residual correlations
  #   y1 ~~ y5
  #   y2 ~~ y4 + y6
  #   y3 ~~ y7
  #   y4 ~~ y8
  #   y6 ~~ y8
'
```

To call
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md), a data
frame of factor scores is needed for all latent variables. To get this
data frame, apply
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) to
all latent variables and specify model parameters as their respective
definitions. Combine factor scores for all latent variables using
[`cbind()`](https://rdrr.io/r/base/cbind.html) to call
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) for
model building.

``` r

fs_dat_ind60 <- get_fs(PoliticalDemocracy, 
                       model = "ind60 =~ x1 + x2 + x3")
fs_dat_dem60 <- get_fs(PoliticalDemocracy, 
                       model = "dem60 =~ y1 + y2 + y3 + y4")
fs_dat_dem65 <- get_fs(PoliticalDemocracy, 
                       model = "dem65 =~ y5 + y6 + y7 + y8")
fs_dat <- cbind(fs_dat_ind60, fs_dat_dem60, fs_dat_dem65)
```

We recently made
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) to
support a multi-factor model input, but computing the factor scores
separately for each latent variable via the argument `local = TRUE`.

``` r

# This is equivalent to the previous chunk:
fs_dat <- get_fs(PoliticalDemocracy, 
                 model = "ind60 =~ x1 + x2 + x3
                          dem60 =~ y1 + y2 + y3 + y4
                          dem65 =~ y5 + y6 + y7 + y8",
                 local = TRUE)
head(fs_dat)  # check that the error covariances and cross-loadings are zero.
```

    ##     fs_ind60   fs_dem60  fs_dem65 fs_ind60_se fs_dem60_se fs_dem65_se
    ## 1 -0.5261683 -2.7487224 -1.371719   0.1213615   0.6756472   0.5724405
    ## 2  0.1436527 -3.0360803 -0.950851   0.1213615   0.6756472   0.5724405
    ## 3  0.7143559  2.6718589  2.738012   0.1213615   0.6756472   0.5724405
    ## 4  1.2399257  2.9936997  1.785091   0.1213615   0.6756472   0.5724405
    ## 5  0.8319080  1.9242932  1.544704   0.1213615   0.6756472   0.5724405
    ## 6  0.2123845  0.9922798 -1.050841   0.1213615   0.6756472   0.5724405
    ##   ind60_by_fs_ind60 ind60_by_fs_dem60 ind60_by_fs_dem65 dem60_by_fs_ind60
    ## 1         0.9657673                 0                 0                 0
    ## 2         0.9657673                 0                 0                 0
    ## 3         0.9657673                 0                 0                 0
    ## 4         0.9657673                 0                 0                 0
    ## 5         0.9657673                 0                 0                 0
    ## 6         0.9657673                 0                 0                 0
    ##   dem60_by_fs_dem60 dem60_by_fs_dem65 dem65_by_fs_ind60 dem65_by_fs_dem60
    ## 1         0.8868049                 0                 0                 0
    ## 2         0.8868049                 0                 0                 0
    ## 3         0.8868049                 0                 0                 0
    ## 4         0.8868049                 0                 0                 0
    ## 5         0.8868049                 0                 0                 0
    ## 6         0.8868049                 0                 0                 0
    ##   dem65_by_fs_dem65 ev_fs_ind60 ecov_fs_dem60_fs_ind60 ev_fs_dem60
    ## 1         0.8998252  0.01472862                      0   0.4564991
    ## 2         0.8998252  0.01472862                      0   0.4564991
    ## 3         0.8998252  0.01472862                      0   0.4564991
    ## 4         0.8998252  0.01472862                      0   0.4564991
    ## 5         0.8998252  0.01472862                      0   0.4564991
    ## 6         0.8998252  0.01472862                      0   0.4564991
    ##   ecov_fs_dem65_fs_ind60 ecov_fs_dem65_fs_dem60 ev_fs_dem65
    ## 1                      0                      0   0.3276882
    ## 2                      0                      0   0.3276882
    ## 3                      0                      0   0.3276882
    ## 4                      0                      0   0.3276882
    ## 5                      0                      0   0.3276882
    ## 6                      0                      0   0.3276882

To build a Two-Stage Path Analysis model, simply call the
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
function with model = regressions (for all predictors) and data = the
factor score data frame created by combining results from
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md).
Each score’s standard error is available in the column named
`fs_[variable name]_se`.[^1] For example, the standard error of latent
variable ind60 is in column `fs_ind60_se` of the `fs_dat` data frame.
Because the combined frame carries each score’s `fs_[variable name]_se`
column from
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md),
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) reads
those standard errors from the columns automatically, so the `se_fs`
argument can be omitted (passing it explicitly as a list or a data frame
still works). In the output model, the values of error constraints are
computed by squaring the standard errors.

``` r

tspa_3var_fit <- tspa(model = "dem60 ~ ind60
                               dem65 ~ ind60 + dem60", 
                      data = fs_dat)
cat(attributes(tspa_3var_fit)$tspaModel)
```

    ## # latent variables (indicated by factor scores)
    ## ind60 =~ 0.965767270434244 * fs_ind60 + 0 * fs_dem60 + 0 * fs_dem65
    ## # latent variables (indicated by factor scores)
    ## dem60 =~ 0 * fs_ind60 + 0.886804906625876 * fs_dem60 + 0 * fs_dem65
    ## # latent variables (indicated by factor scores)
    ## dem65 =~ 0 * fs_ind60 + 0 * fs_dem60 + 0.899825222245714 * fs_dem65
    ## # constrain the errors
    ## fs_ind60 ~~ 0.0147286194470875 * fs_ind60
    ## # constrain the errors
    ## fs_dem60 ~~ 0 * fs_ind60
    ## # constrain the errors
    ## fs_dem65 ~~ 0 * fs_ind60
    ## # constrain the errors
    ## fs_dem60 ~~ 0.456499146148539 * fs_dem60
    ## # constrain the errors
    ## fs_dem65 ~~ 0 * fs_dem60
    ## # constrain the errors
    ## fs_dem65 ~~ 0.327688166837143 * fs_dem65
    ## # constrain the intercepts
    ## fs_ind60 ~ 0 * 1
    ## # constrain the intercepts
    ## fs_dem60 ~ 0 * 1
    ## # constrain the intercepts
    ## fs_dem65 ~ 0 * 1
    ## # structural model
    ## dem60 ~ ind60
    ##                                dem65 ~ ind60 + dem60

## Multigroup, multiple factors

``` r

model <- ' 
  # latent variable definitions
    visual =~ x1 + x2 + x3
    speed =~ x7 + x8 + x9

  # regressions
    visual ~ speed
'
```

To call
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md), a data
frame of factor scores is needed for all latent variables for each
group. To get this data frame, you can use the above strategy to compute
factor scores separately for each latent variable with the
`local = TRUE` argument. Alternatively, you can use the default
(`local = FALSE`), which computes factor scores for all latent variables
simultaneously, but resulting in correlated measurement errors and
cross-loadings with the regression method for factor scores (see more
discussion in the
[`Scoring Matrices: lavaan CFA and lme4`](https://mmm-lab-um.github.io/R2spa/articles/scoring-matrices.md)
vignette).

``` r

hs_mod <- '
visual =~ x1 + x2 + x3
speed =~ x7 + x8 + x9
'

# get factor scores (joint model)
fs_hs <- get_fs(HolzingerSwineford1939, model = hs_mod, group = "school")
head(fs_hs)  # with cross-loadings and correlated measurement errors
```

    ##      fs_visual      fs_speed fs_visual_se fs_speed_se visual_by_fs_visual
    ## 1 -0.818816609  0.0004278866    0.3483615   0.2590721           0.6513207
    ## 2 -0.021268888  0.5125828715    0.3483615   0.2590721           0.6513207
    ## 3 -0.494499099 -0.7827848448    0.3483615   0.2590721           0.6513207
    ## 4  0.385267108 -0.2971447381    0.3483615   0.2590721           0.6513207
    ## 5 -0.637443617  0.1213326700    0.3483615   0.2590721           0.6513207
    ## 6  0.008040415  0.6167608726    0.3483615   0.2590721           0.6513207
    ##   visual_by_fs_speed speed_by_fs_visual speed_by_fs_speed ev_fs_visual
    ## 1         0.06263767          0.1017434         0.6168097    0.1213557
    ## 2         0.06263767          0.1017434         0.6168097    0.1213557
    ## 3         0.06263767          0.1017434         0.6168097    0.1213557
    ## 4         0.06263767          0.1017434         0.6168097    0.1213557
    ## 5         0.06263767          0.1017434         0.6168097    0.1213557
    ## 6         0.06263767          0.1017434         0.6168097    0.1213557
    ##   ecov_fs_speed_fs_visual ev_fs_speed  school
    ## 1              0.02238692  0.06711836 Pasteur
    ## 2              0.02238692  0.06711836 Pasteur
    ## 3              0.02238692  0.06711836 Pasteur
    ## 4              0.02238692  0.06711836 Pasteur
    ## 5              0.02238692  0.06711836 Pasteur
    ## 6              0.02238692  0.06711836 Pasteur

To build a Two-Stage Path Analysis model, simply call the
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
function with model = the regression relation, data = the combined
factor score data frame, and `group` = the grouping column (here
`"school"`). Each score’s standard error is available in the column
named `fs_[variable name]_se`. For example, the standard error of the
multigroup variable visual is in column `fs_visual_se` of the `fs_hs`
data frame. Because `fs_hs` comes from
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md),
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) reads
the per-group standard errors from those columns automatically, so the
`se_fs` argument can be omitted (passing it explicitly as a data frame
still works). The [`unique()`](https://rdrr.io/r/base/unique.html)
function is handy for inspecting the per-group standard errors: for
example, `unique(fs_hs$fs_visual_se)` shows the standard error of visual
in each group.

Function
[`standardizedsolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
enables user to view a table of the standardized estimate, standard
error, z score, p-value, and confidence interval bounds for each
multigroup regression relation.

``` r

# tspa model
tspa_fit <- tspa(model = "visual ~ speed",
                 data = fs_hs,
                 group = "school"
                 # group.equal = "regressions"
                 )
stdsol <- standardizedsolution(tspa_fit)
subset(stdsol, subset = op == "~")
```

    ##       lhs op   rhs group est.std    se     z pvalue ci.lower ci.upper
    ## 10 visual  ~ speed     1   0.335 0.115 2.920  0.003    0.110    0.560
    ## 24 visual  ~ speed     2   0.500 0.094 5.332  0.000    0.316    0.684

To view the Two-Stage Path Analysis model for multigroup, use
`attributes([model name])$tspaModel`. Function
[`cat()`](https://rdrr.io/r/base/cat.html) can help tidy up the model
output. In the output, the values of error constraints are computed by
squaring standard errors from the previous section.

``` r

cat(attributes(tspa_fit)$tspaModel)
```

    ## # latent variables (indicated by factor scores)
    ## visual =~ c(0.651320651151504, 0.630164103758516) * fs_visual + c(0.062637668106974, 0.0735865303140717) * fs_speed
    ## # latent variables (indicated by factor scores)
    ## speed =~ c(0.101743373278373, 0.163519554160675) * fs_visual + c(0.616809744932919, 0.757903096858229) * fs_speed
    ## # constrain the errors
    ## fs_visual ~~ c(0.121355749245672, 0.110660002585844) * fs_visual
    ## # constrain the errors
    ## fs_speed ~~ c(0.0223869195387284, 0.0276266319224011) * fs_visual
    ## # constrain the errors
    ## fs_speed ~~ c(0.0671183590966231, 0.0713803547181482) * fs_speed
    ## # constrain the intercepts
    ## fs_visual ~ c(0, 0) * 1
    ## # constrain the intercepts
    ## fs_speed ~ c(0, 0) * 1
    ## # structural model
    ## visual ~ speed

## Related vignettes

- [`Multiple factors`](https://mmm-lab-um.github.io/R2spa/articles/multiple-factors.md)
  — a score that is a composite of several indicators: the `fsT`/`fsL`
  path and
  [`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
  cross-loadings.
- [`2S-PA with Random Effects`](https://mmm-lab-um.github.io/R2spa/articles/multilevel.md)
  — stage 1 from an `lme4` mixed model, feeding empirical-Bayes random
  effects into 2S-PA.
- [`Scoring Matrices: lavaan CFA and lme4`](https://mmm-lab-um.github.io/R2spa/articles/scoring-matrices.md)
  — the per-observation scoring matrices
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  attaches, and the empirical-Bayes identity behind them.
- [`Correction to Measurement Error`](https://mmm-lab-um.github.io/R2spa/articles/correction-error.md)
  — the measurement-error correction factor and its effect on standard
  errors.
- [`EFA Scores`](https://mmm-lab-um.github.io/R2spa/articles/efa-score.md)
  — stage 1 from an exploratory factor analysis (EFA) rather than a CFA.
- [`Linear Growth Modeling with 2S-PA`](https://mmm-lab-um.github.io/R2spa/articles/tspa-growth-vignette.md)
  — 2S-PA for longitudinal / growth-curve measurement models.
- [`2S-PA with OpenMx and IRT (mirt)`](https://mmm-lab-um.github.io/R2spa/articles/tspa-vignette-mx.md)
  — stage 2 via
  [`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
  (OpenMx), including the `mirt` (IRT) path.

## References

- Cole, D. A., & Preacher, K. J. (2014). Manifest variable path
  analysis: Potentially serious and misleading consequences due to
  uncorrected measurement error. *Psychological Methods, 19*(2),
  300–315. <https://doi.org/10.1037/a0033805>
- Lai, M. H. C., & Hsiao, Y.-Y. (2022). Two-stage path analysis with
  definition variables: An alternative framework to account for
  measurement error. *Psychological Methods, 27*(4), 568–588.
  <https://doi.org/10.1037/met0000410>
- Lai, M. H. C., Tse, W. W.-Y., Zhang, G., Li, Y., & Hsiao, Y.-Y.
  (2023). Correcting for unreliability and partial invariance: A
  two-stage path analysis approach. *Structural Equation Modeling: A
  Multidisciplinary Journal, 30*(2), 258–271.
  <https://doi.org/10.1080/10705511.2022.2125397>
- Lai, M. H. C., & Liu, S. (2026). A two-stage approach to account for
  measurement error when using empirical Bayes estimates of random
  slopes. *Psychological Methods.* Advance online publication.
  <https://doi.org/10.1037/met0000838>
- Rosseel, Y. (2012). *lavaan*: An R package for structural equation
  modeling. *Journal of Statistical Software, 48*(2), 1–36.
  <https://doi.org/10.18637/jss.v048.i02>

[^1]: The text inside `[]` should be replaced with the actual variable
    name.
