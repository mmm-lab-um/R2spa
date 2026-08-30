# Multi-Factor Measurement Model

``` r

library(lavaan)
library(R2spa)
```

The example is from <https://lavaan.ugent.be/tutorial/sem.html>.

## Factor score

For a CFA model with multiple latent factors, even when each indicator
only loads on one factor, the resulting factor scores generally are
weighted composites of **ALL** indicators. Consider the regression
score, which has the form

``` math
\tilde{\boldsymbol \eta} = \mathbf{A}(\mathbf{y} - \hat{\boldsymbol \mu}) + \boldsymbol{\alpha}
```

where
$`\mathbf{A} = \boldsymbol{\Psi}\boldsymbol{\Lambda}^\top \hat{\boldsymbol{\Sigma}}^{-1}`$
is a $`q`$$`\times`$$`p`$ matrix,
$`\hat{\boldsymbol \mu} = \boldsymbol{\nu} + \boldsymbol{\Lambda} \boldsymbol{\alpha}`$
and
$`\boldsymbol{\Sigma} = \boldsymbol{\Lambda} \boldsymbol{\Psi}\boldsymbol{\Lambda}^\top + \boldsymbol{\Theta}`$
are the model-implied means and covariances of the indicators
$`\mathbf{y}`$, and $`\boldsymbol{\alpha}`$ and $`\boldsymbol{\Psi}`$
are the latent means and latent covariances.

Therefore, assuming that the model is correctly specified such that
$`\mathbf{y} = \boldsymbol{\nu} + \boldsymbol{\Lambda} \boldsymbol{\eta} + \boldsymbol{\varepsilon}`$,

``` math
  \begin{aligned}
  \tilde{\boldsymbol \eta} & = \mathbf{A}(\boldsymbol{\Lambda} \boldsymbol{\eta} + \boldsymbol{\varepsilon} - \boldsymbol{\Lambda} \boldsymbol{\alpha}) + \boldsymbol{\alpha} \\
  & = (\mathbf{I} - \mathbf{A}\boldsymbol{\Lambda}) \boldsymbol{\alpha} +
  \mathbf{A}\boldsymbol{\Lambda} \boldsymbol{\eta} + \mathbf{A}\boldsymbol{\varepsilon}.
  \end{aligned}
```

If we consider $`\tilde{\boldsymbol \eta}`$ as indicators of
$`\boldsymbol \eta`$, we can see that
$`\boldsymbol{\nu}_\tilde{\boldsymbol{\eta}} = (\mathbf{I} - \mathbf{A}\boldsymbol{\Lambda}) \boldsymbol{\alpha}`$
is the intercept,
$`\boldsymbol{\Lambda}_\tilde{\boldsymbol{\eta}} = \mathbf{A}\boldsymbol{\Lambda}`$
is the $`q`$$`\times`$$`q`$ loading matrix, and
$`\boldsymbol{\Theta}_\tilde{\boldsymbol{\eta}} = \mathbf{A}\boldsymbol{\varepsilon}`$
is the error covariance matrix.

We can see that $`\boldsymbol{\Lambda}_\tilde{\boldsymbol{\eta}}`$ is
generally not diagonal, as the following shows:

``` r

# CFA
my_cfa <- "
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
"
cfa_fit <- cfa(model = my_cfa,
               data  = PoliticalDemocracy,
               std.lv = TRUE)
# A matrix
pars <- lavInspect(cfa_fit, what = "est")
lambda_mat <- pars$lambda
psi_mat <- pars$psi
sigma_mat <- cfa_fit@implied$cov[[1]]
ginvsigma <- MASS::ginv(sigma_mat)
alambda <- psi_mat %*% crossprod(lambda_mat, ginvsigma %*% lambda_mat)
alambda
#>            ind60      dem60
#> ind60 0.95538579 0.01834111
#> dem60 0.05816994 0.86888887
```

We can also use
[`R2spa::get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md):

``` r

(fs_dat <- get_fs(PoliticalDemocracy, model = my_cfa, std.lv = TRUE)) |> head()
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

Therefore, cross-loadings must be accounted for in stage 2. Because
`fs_dat` comes from
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md),
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) reads
the implied loading matrix `fsL` (which here carries the cross-loadings)
and the score error covariance `fsT` from the result’s attributes, so
the `fsT`/`fsL` arguments no longer need to be passed (they can still be
passed explicitly if preferred).

``` r

tspa_fit <- tspa(model = "dem60 ~ ind60", data = fs_dat)
cat(attr(tspa_fit, "tspaModel"))
#> # latent variables (indicated by factor scores)
#> ind60 =~ 0.955385785456618 * fs_ind60 + 0.0581699407736301 * fs_dem60
#> # latent variables (indicated by factor scores)
#> dem60 =~ 0.018341113191188 * fs_ind60 + 0.868888868390615 * fs_dem60
#> # constrain the errors
#> fs_ind60 ~~ 0.0345955179271982 * fs_ind60
#> # constrain the errors
#> fs_dem60 ~~ 0.00401738814566846 * fs_ind60
#> # constrain the errors
#> fs_dem60 ~~ 0.0907757082269276 * fs_dem60
#> # constrain the intercepts
#> fs_ind60 ~ 0 * 1
#> # constrain the intercepts
#> fs_dem60 ~ 0 * 1
#> # structural model
#> dem60 ~ ind60
parameterestimates(tspa_fit)
#>         lhs op      rhs   est    se     z pvalue ci.lower ci.upper
#> 1     ind60 =~ fs_ind60 0.955 0.000    NA     NA    0.955    0.955
#> 2     ind60 =~ fs_dem60 0.058 0.000    NA     NA    0.058    0.058
#> 3     dem60 =~ fs_ind60 0.018 0.000    NA     NA    0.018    0.018
#> 4     dem60 =~ fs_dem60 0.869 0.000    NA     NA    0.869    0.869
#> 5  fs_ind60 ~~ fs_ind60 0.035 0.000    NA     NA    0.035    0.035
#> 6  fs_ind60 ~~ fs_dem60 0.004 0.000    NA     NA    0.004    0.004
#> 7  fs_dem60 ~~ fs_dem60 0.091 0.000    NA     NA    0.091    0.091
#> 8  fs_ind60 ~1          0.000 0.000    NA     NA    0.000    0.000
#> 9  fs_dem60 ~1          0.000 0.000    NA     NA    0.000    0.000
#> 10    dem60  ~    ind60 0.460 0.113 4.089      0    0.240    0.681
#> 11    ind60 ~~    ind60 1.000 0.169 5.900      0    0.668    1.332
#> 12    dem60 ~~    dem60 0.788 0.150 5.267      0    0.495    1.081
#> 13    ind60 ~1          0.000 0.000    NA     NA    0.000    0.000
#> 14    dem60 ~1          0.000 0.000    NA     NA    0.000    0.000
```

The non-diagonal `fsL` (cross-loadings) is what makes the score-based
model recover the joint-model structure, because each score is a
weighted composite of all the indicators.

## Three-factor model example

``` r

# CFA
cfa_3fac <-  "
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
"
cfa_3fac_fit <- cfa(model = cfa_3fac,
                    data  = PoliticalDemocracy,
                    std.lv = TRUE)
# A matrix
pars <- lavInspect(cfa_3fac_fit, what = "est")
lambda_mat <- pars$lambda
psi_mat <- pars$psi
sigma_mat <- cfa_3fac_fit@implied$cov[[1]]
ginvsigma <- MASS::ginv(sigma_mat)
alambda <- psi_mat %*% crossprod(lambda_mat, ginvsigma %*% lambda_mat)
alambda
#>             ind60        dem60      dem65
#> ind60  0.95064774 -0.005967124 0.02951139
#> dem60 -0.02069724  0.533020047 0.41603787
#> dem65  0.09868528  0.401095944 0.49133377
```

We can also use
[`R2spa::get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md):

``` r

(fs_dat_3fac <- get_fs(PoliticalDemocracy, model = cfa_3fac, std.lv = TRUE)) |>
  head()
#>     fs_ind60   fs_dem60    fs_dem65 fs_ind60_se fs_dem60_se fs_dem65_se
#> 1 -0.7990475 -1.1571745 -1.13720655   0.1844193   0.2422541   0.2270976
#> 2  0.2152486 -1.0238236 -0.80871922   0.1844193   0.2422541   0.2270976
#> 3  1.1028297  1.3890842  1.46389950   0.1844193   0.2422541   0.2270976
#> 4  1.8585004  1.3163888  1.43045560   0.1844193   0.2422541   0.2270976
#> 5  1.2432985  0.9522026  1.03348589   0.1844193   0.2422541   0.2270976
#> 6  0.3067994  0.1109206  0.05023217   0.1844193   0.2422541   0.2270976
#>   ind60_by_fs_ind60 ind60_by_fs_dem60 ind60_by_fs_dem65 dem60_by_fs_ind60
#> 1         0.9506477       -0.02069724        0.09868528      -0.005967124
#> 2         0.9506477       -0.02069724        0.09868528      -0.005967124
#> 3         0.9506477       -0.02069724        0.09868528      -0.005967124
#> 4         0.9506477       -0.02069724        0.09868528      -0.005967124
#> 5         0.9506477       -0.02069724        0.09868528      -0.005967124
#> 6         0.9506477       -0.02069724        0.09868528      -0.005967124
#>   dem60_by_fs_dem60 dem60_by_fs_dem65 dem65_by_fs_ind60 dem65_by_fs_dem60
#> 1           0.53302         0.4010959        0.02951139         0.4160379
#> 2           0.53302         0.4010959        0.02951139         0.4160379
#> 3           0.53302         0.4010959        0.02951139         0.4160379
#> 4           0.53302         0.4010959        0.02951139         0.4160379
#> 5           0.53302         0.4010959        0.02951139         0.4160379
#> 6           0.53302         0.4010959        0.02951139         0.4160379
#>   dem65_by_fs_dem65 ev_fs_ind60 ecov_fs_dem60_fs_ind60 ev_fs_dem60
#> 1         0.4913338   0.0340105           0.0003881641  0.05868703
#> 2         0.4913338   0.0340105           0.0003881641  0.05868703
#> 3         0.4913338   0.0340105           0.0003881641  0.05868703
#> 4         0.4913338   0.0340105           0.0003881641  0.05868703
#> 5         0.4913338   0.0340105           0.0003881641  0.05868703
#> 6         0.4913338   0.0340105           0.0003881641  0.05868703
#>   ecov_fs_dem65_fs_ind60 ecov_fs_dem65_fs_dem60 ev_fs_dem65
#> 1            0.005026024             0.05337525   0.0515733
#> 2            0.005026024             0.05337525   0.0515733
#> 3            0.005026024             0.05337525   0.0515733
#> 4            0.005026024             0.05337525   0.0515733
#> 5            0.005026024             0.05337525   0.0515733
#> 6            0.005026024             0.05337525   0.0515733
```

As before,
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) reads
the cross-loadings (`fsL`) and the score error covariance (`fsT`) from
the [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result.

``` r

tspa_fit_3fac <- tspa(model = "dem60 ~ ind60
                               dem65 ~ ind60 + dem60",
              data = fs_dat_3fac)
cat(attr(tspa_fit_3fac, "tspaModel"))
#> # latent variables (indicated by factor scores)
#> ind60 =~ 0.950647742844845 * fs_ind60 + -0.0206972362902852 * fs_dem60 + 0.0986852834195762 * fs_dem65
#> # latent variables (indicated by factor scores)
#> dem60 =~ -0.00596712444175296 * fs_ind60 + 0.533020047181513 * fs_dem60 + 0.401095943690547 * fs_dem65
#> # latent variables (indicated by factor scores)
#> dem65 =~ 0.0295113941487123 * fs_ind60 + 0.416037872255943 * fs_dem60 + 0.491333767947423 * fs_dem65
#> # constrain the errors
#> fs_ind60 ~~ 0.0340104951546672 * fs_ind60
#> # constrain the errors
#> fs_dem60 ~~ 0.000388164070398798 * fs_ind60
#> # constrain the errors
#> fs_dem65 ~~ 0.00502602420598861 * fs_ind60
#> # constrain the errors
#> fs_dem60 ~~ 0.0586870313996517 * fs_dem60
#> # constrain the errors
#> fs_dem65 ~~ 0.0533752457536214 * fs_dem60
#> # constrain the errors
#> fs_dem65 ~~ 0.051573299369388 * fs_dem65
#> # constrain the intercepts
#> fs_ind60 ~ 0 * 1
#> # constrain the intercepts
#> fs_dem60 ~ 0 * 1
#> # constrain the intercepts
#> fs_dem65 ~ 0 * 1
#> # structural model
#> dem60 ~ ind60
#>                                dem65 ~ ind60 + dem60
parameterestimates(tspa_fit_3fac)
#>         lhs op      rhs    est    se      z pvalue ci.lower ci.upper
#> 1     ind60 =~ fs_ind60  0.951 0.000     NA     NA    0.951    0.951
#> 2     ind60 =~ fs_dem60 -0.021 0.000     NA     NA   -0.021   -0.021
#> 3     ind60 =~ fs_dem65  0.099 0.000     NA     NA    0.099    0.099
#> 4     dem60 =~ fs_ind60 -0.006 0.000     NA     NA   -0.006   -0.006
#> 5     dem60 =~ fs_dem60  0.533 0.000     NA     NA    0.533    0.533
#> 6     dem60 =~ fs_dem65  0.401 0.000     NA     NA    0.401    0.401
#> 7     dem65 =~ fs_ind60  0.030 0.000     NA     NA    0.030    0.030
#> 8     dem65 =~ fs_dem60  0.416 0.000     NA     NA    0.416    0.416
#> 9     dem65 =~ fs_dem65  0.491 0.000     NA     NA    0.491    0.491
#> 10 fs_ind60 ~~ fs_ind60  0.034 0.000     NA     NA    0.034    0.034
#> 11 fs_ind60 ~~ fs_dem60  0.000 0.000     NA     NA    0.000    0.000
#> 12 fs_ind60 ~~ fs_dem65  0.005 0.000     NA     NA    0.005    0.005
#> 13 fs_dem60 ~~ fs_dem60  0.059 0.000     NA     NA    0.059    0.059
#> 14 fs_dem60 ~~ fs_dem65  0.053 0.000     NA     NA    0.053    0.053
#> 15 fs_dem65 ~~ fs_dem65  0.052 0.000     NA     NA    0.052    0.052
#> 16 fs_ind60 ~1           0.000 0.000     NA     NA    0.000    0.000
#> 17 fs_dem60 ~1           0.000 0.000     NA     NA    0.000    0.000
#> 18 fs_dem65 ~1           0.000 0.000     NA     NA    0.000    0.000
#> 19    dem60  ~    ind60  0.448 0.114  3.937  0.000    0.225    0.671
#> 20    dem65  ~    ind60  0.146 0.069  2.112  0.035    0.010    0.281
#> 21    dem65  ~    dem60  0.913 0.073 12.435  0.000    0.769    1.057
#> 22    ind60 ~~    ind60  1.000 0.169  5.902  0.000    0.668    1.332
#> 23    dem60 ~~    dem60  0.799 0.153  5.224  0.000    0.499    1.099
#> 24    dem65 ~~    dem65  0.026 0.043  0.620  0.535   -0.057    0.110
#> 25    ind60 ~1           0.000 0.000     NA     NA    0.000    0.000
#> 26    dem60 ~1           0.000 0.000     NA     NA    0.000    0.000
#> 27    dem65 ~1           0.000 0.000     NA     NA    0.000    0.000
```

Compare to SEM:

``` r

sem_3fac <- sem("
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # structural model
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
  ",
  data = PoliticalDemocracy
)
standardizedsolution(sem_3fac)
#>      lhs op   rhs est.std    se      z pvalue ci.lower ci.upper
#> 1  ind60 =~    x1   0.920 0.023 39.823  0.000    0.874    0.965
#> 2  ind60 =~    x2   0.973 0.017 58.858  0.000    0.941    1.006
#> 3  ind60 =~    x3   0.872 0.031 28.034  0.000    0.811    0.933
#> 4  dem60 =~    y1   0.845 0.039 21.698  0.000    0.769    0.921
#> 5  dem60 =~    y2   0.760 0.054 14.142  0.000    0.655    0.866
#> 6  dem60 =~    y3   0.705 0.063 11.225  0.000    0.582    0.828
#> 7  dem60 =~    y4   0.860 0.036 23.650  0.000    0.789    0.931
#> 8  dem65 =~    y5   0.803 0.046 17.602  0.000    0.714    0.893
#> 9  dem65 =~    y6   0.783 0.049 15.918  0.000    0.687    0.879
#> 10 dem65 =~    y7   0.819 0.043 19.122  0.000    0.735    0.903
#> 11 dem65 =~    y8   0.847 0.038 22.389  0.000    0.773    0.921
#> 12 dem60  ~ ind60   0.448 0.102  4.393  0.000    0.248    0.648
#> 13 dem65  ~ ind60   0.146 0.070  2.071  0.038    0.008    0.283
#> 14 dem65  ~ dem60   0.913 0.048 19.120  0.000    0.819    1.006
#> 15    x1 ~~    x1   0.154 0.042  3.636  0.000    0.071    0.238
#> 16    x2 ~~    x2   0.053 0.032  1.634  0.102   -0.010    0.116
#> 17    x3 ~~    x3   0.240 0.054  4.417  0.000    0.133    0.346
#> 18    y1 ~~    y1   0.286 0.066  4.348  0.000    0.157    0.415
#> 19    y2 ~~    y2   0.422 0.082  5.166  0.000    0.262    0.582
#> 20    y3 ~~    y3   0.503 0.089  5.676  0.000    0.329    0.676
#> 21    y4 ~~    y4   0.261 0.063  4.173  0.000    0.138    0.383
#> 22    y5 ~~    y5   0.355 0.073  4.842  0.000    0.211    0.499
#> 23    y6 ~~    y6   0.387 0.077  5.024  0.000    0.236    0.538
#> 24    y7 ~~    y7   0.329 0.070  4.696  0.000    0.192    0.467
#> 25    y8 ~~    y8   0.283 0.064  4.416  0.000    0.157    0.408
#> 26 ind60 ~~ ind60   1.000 0.000     NA     NA    1.000    1.000
#> 27 dem60 ~~ dem60   0.799 0.091  8.737  0.000    0.620    0.978
#> 28 dem65 ~~ dem65   0.026 0.046  0.579  0.562   -0.063    0.116
```

The [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html) fit above
is the **gold-standard** joint model: it estimates the measurement model
and the structural paths together, treating the indicators as observed,
so no separate measurement-error correction is needed.
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) reaches
the same structural estimates a different way — by correcting the
stage-1 *scores* for their measurement error — so the two agree when the
sample is large. The simulation studies find 2S-PA tracks the joint
model closely while staying more reliable, with less biased standard
errors and better coverage, in small samples and with unreliable or
categorical items (Lai & Hsiao, 2022; Lai et al., 2023).
`PoliticalDemocracy` is a real data set, so there is no “true value” to
compare against here; the point of putting the two side by side is that
they target the same structural model.

## Local vs joint scores

The [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
calls above all score from a *joint* multi-factor measurement model: one
fit with all the latents estimated together. With freely correlated
factors, each score is a weighted composite of **all** the indicators
(the off-diagonals of the `alambda` matrices above), so the stage-2
model has to account for the cross-loadings explicitly. The canonical
2S-PA setup (Lai & Hsiao, 2022) instead scores each latent from its
**own** local measurement model. `get_fs(..., local = TRUE)` does this
in one call and merges the per-construct results into the usual
multi-factor layout:

``` r

(fs_local <- get_fs(PoliticalDemocracy,
                    model = cfa_3fac,
                    std.lv = TRUE,
                    local = TRUE)) |> head()
#>     fs_ind60   fs_dem60   fs_dem65 fs_ind60_se fs_dem60_se fs_dem65_se
#> 1 -0.7883166 -1.2889577 -0.7194375   0.1818264   0.3168311   0.3002329
#> 2  0.2152236 -1.4237084 -0.4987011   0.1818264   0.3168311   0.3002329
#> 3  1.0702633  1.2529141  1.4360292   0.1818264   0.3168311   0.3002329
#> 4  1.8576831  1.4038349  0.9362424   0.1818264   0.3168311   0.3002329
#> 5  1.2463824  0.9023585  0.8101644   0.1818264   0.3168311   0.3002329
#> 6  0.3181990  0.4653096 -0.5511437   0.1818264   0.3168311   0.3002329
#>   ind60_by_fs_ind60 ind60_by_fs_dem60 ind60_by_fs_dem65 dem60_by_fs_ind60
#> 1         0.9657673                 0                 0                 0
#> 2         0.9657673                 0                 0                 0
#> 3         0.9657673                 0                 0                 0
#> 4         0.9657673                 0                 0                 0
#> 5         0.9657673                 0                 0                 0
#> 6         0.9657673                 0                 0                 0
#>   dem60_by_fs_dem60 dem60_by_fs_dem65 dem65_by_fs_ind60 dem65_by_fs_dem60
#> 1         0.8868049                 0                 0                 0
#> 2         0.8868049                 0                 0                 0
#> 3         0.8868049                 0                 0                 0
#> 4         0.8868049                 0                 0                 0
#> 5         0.8868049                 0                 0                 0
#> 6         0.8868049                 0                 0                 0
#>   dem65_by_fs_dem65 ev_fs_ind60 ecov_fs_dem60_fs_ind60 ev_fs_dem60
#> 1         0.8998252  0.03306085                      0    0.100382
#> 2         0.8998252  0.03306085                      0    0.100382
#> 3         0.8998252  0.03306085                      0    0.100382
#> 4         0.8998252  0.03306085                      0    0.100382
#> 5         0.8998252  0.03306085                      0    0.100382
#> 6         0.8998252  0.03306085                      0    0.100382
#>   ecov_fs_dem65_fs_ind60 ecov_fs_dem65_fs_dem60 ev_fs_dem65
#> 1                      0                      0   0.0901398
#> 2                      0                      0   0.0901398
#> 3                      0                      0   0.0901398
#> 4                      0                      0   0.0901398
#> 5                      0                      0   0.0901398
#> 6                      0                      0   0.0901398
```

The merged result has the same layout as the joint output, with the
cross terms *exactly* zero: the `fsT` and `fsL` attributes are
block-diagonal, the `psi` attribute is diagonal, and the off-diagonal
`_by_` loading columns and all `ecov_*` columns are zero. Each score is
a *pure* per-construct score — the cross-factor structure is not
estimated by design:

``` r

attr(fs_local, "fsL")[[1]]
#>              ind60     dem60     dem65
#> fs_ind60 0.9657673 0.0000000 0.0000000
#> fs_dem60 0.0000000 0.8868049 0.0000000
#> fs_dem65 0.0000000 0.0000000 0.8998252
sapply(c("fs_ind60", "fs_dem60", "fs_dem65"),
       function(col) max(abs(fs_local[[col]] - fs_dat_3fac[[col]])))
#>   fs_ind60   fs_dem60   fs_dem65 
#> 0.06198296 0.53354046 0.76242550
```

So local scores are not the same as joint scores: with freely correlated
factors they differ (the maximum score differences above), because the
joint fit scores every latent from all the indicators — even the joint
fit’s own per-factor estimates shift under free correlation. With the
factors constrained uncorrelated (e.g. `ind60 ~~ 0 * dem60` in the joint
model), the likelihood factorizes and the two agree to optimizer
tolerance.

Because the merged result carries the usual `fsT`/`fsL` attributes (with
the zero cross terms), it feeds
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
directly with no explicit `fsT`/`fsL` arguments — the stage-2 model
picks up the zero cross-loadings and score error covariances from the
attributes:

``` r

tspa_fit_local <- tspa(model = "dem60 ~ ind60
                        dem65 ~ ind60 + dem60",
                       data = fs_local)
parameterestimates(tspa_fit_local)
#>         lhs op      rhs   est    se      z pvalue ci.lower ci.upper
#> 1     ind60 =~ fs_ind60 0.966 0.000     NA     NA    0.966    0.966
#> 2     ind60 =~ fs_dem60 0.000 0.000     NA     NA    0.000    0.000
#> 3     ind60 =~ fs_dem65 0.000 0.000     NA     NA    0.000    0.000
#> 4     dem60 =~ fs_ind60 0.000 0.000     NA     NA    0.000    0.000
#> 5     dem60 =~ fs_dem60 0.887 0.000     NA     NA    0.887    0.887
#> 6     dem60 =~ fs_dem65 0.000 0.000     NA     NA    0.000    0.000
#> 7     dem65 =~ fs_ind60 0.000 0.000     NA     NA    0.000    0.000
#> 8     dem65 =~ fs_dem60 0.000 0.000     NA     NA    0.000    0.000
#> 9     dem65 =~ fs_dem65 0.900 0.000     NA     NA    0.900    0.900
#> 10 fs_ind60 ~~ fs_ind60 0.033 0.000     NA     NA    0.033    0.033
#> 11 fs_ind60 ~~ fs_dem60 0.000 0.000     NA     NA    0.000    0.000
#> 12 fs_ind60 ~~ fs_dem65 0.000 0.000     NA     NA    0.000    0.000
#> 13 fs_dem60 ~~ fs_dem60 0.100 0.000     NA     NA    0.100    0.100
#> 14 fs_dem60 ~~ fs_dem65 0.000 0.000     NA     NA    0.000    0.000
#> 15 fs_dem65 ~~ fs_dem65 0.090 0.000     NA     NA    0.090    0.090
#> 16 fs_ind60 ~1          0.000 0.000     NA     NA    0.000    0.000
#> 17 fs_dem60 ~1          0.000 0.000     NA     NA    0.000    0.000
#> 18 fs_dem65 ~1          0.000 0.000     NA     NA    0.000    0.000
#> 19    dem60  ~    ind60 0.453 0.113  4.000  0.000    0.231    0.675
#> 20    dem65  ~    ind60 0.129 0.072  1.787  0.074   -0.012    0.271
#> 21    dem65  ~    dem60 0.898 0.077 11.705  0.000    0.748    1.049
#> 22    ind60 ~~    ind60 1.000 0.169  5.914  0.000    0.669    1.331
#> 23    dem60 ~~    dem60 0.795 0.152  5.235  0.000    0.497    1.092
#> 24    dem65 ~~    dem65 0.071 0.047  1.528  0.126   -0.020    0.163
#> 25    ind60 ~1          0.000 0.000     NA     NA    0.000    0.000
#> 26    dem60 ~1          0.000 0.000     NA     NA    0.000    0.000
#> 27    dem65 ~1          0.000 0.000     NA     NA    0.000    0.000
```

`model` may also be a character vector (or a named list) of complete
single-factor model strings, each element fit verbatim — the escape
hatch for anything the strict string grammar rejects, e.g. a
within-factor residual covariance (`y1 ~~ y4`). `vfsLT` (and hence
`tspa(corrected_se = TRUE)`), `prior_cov`, and `reliability` are not
supported in `local` mode (v1); see
[`?get_fs`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
(Details) for the full contract.
