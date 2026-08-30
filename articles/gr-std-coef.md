# Grand Standardized Coefficients

This vignette shows how to obtain grand standardized estimates for the
latent variables for multi-group `lavaan` objects, using the
[`R2spa::grandStandardizedSolution()`](https://mmm-lab-um.github.io/R2spa/reference/grand_standardized_solution.md)
function.

``` r

library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(R2spa)
```

The example is from <https://lavaan.ugent.be/tutorial/sem.html>.

## Single Group

For single groups, standardized solution can be obtained by first
obtaining the latent variable covariance matrix:

``` math
  \begin{aligned}
    \mathbf{\eta} & = \mathbf{\alpha} + \mathbf{\Gamma} \mathbf{X} + \mathbf{B} \mathbf{\eta} + \mathbf{\zeta} \\
    (\mathbf{I} - \mathbf{B}) \mathbf{\eta} & = \mathbf{\alpha} + \mathbf{\Gamma} \mathbf{X} + \mathbf{\zeta} \\
    (\mathbf{I} - \mathbf{B}) \mathrm{Var}(\mathbf{\eta}) (\mathbf{I} - \mathbf{B})^\top & = \mathrm{Var}(\mathbf{\Gamma} \mathbf{X}) + \mathbf{\Psi} \\
    \mathrm{Var}(\mathbf{\eta}) & = (\mathbf{I} - \mathbf{B})^{-1} [\mathrm{Var}(\mathbf{\Gamma} \mathbf{X}) + \mathbf{\Psi}] {(\mathbf{I} - \mathbf{B})^\top}^{-1}
  \end{aligned}
```

``` r

myModel <- '
   # latent variables
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
   # regressions
     dem60 ~ ind60
'
fit <- sem(model = myModel,
           data  = PoliticalDemocracy)
# Latent variable covariances
lavInspect(fit, "cov.lv")
#>       ind60 dem60
#> ind60 0.449      
#> dem60 0.646 4.382
# Using R2spa:::veta() (internal helper; the ::: is deliberate since it is not exported)
fit_est <- lavInspect(fit, "est")
R2spa:::veta(fit_est$beta, psi = fit_est$psi)
#>           ind60     dem60
#> ind60 0.4485415 0.6455929
#> dem60 0.6455929 4.3824834
```

The standardized estimates in the $`\mathbf{B}`$ matrix are obtained as

``` math
\mathbf{B}_s = \mathbf{S}_\eta^{-1} \mathbf{B} \mathbf{S}_\eta^{-1}
```

where $`\mathbf{S}_\eta`$ is a diagonal matrix containing the square
root of the diagonal elements of $`\mathrm{Var}(\mathbf{\eta})`$

``` r

S_eta <- diag(
  sqrt(diag(
    R2spa:::veta(fit_est$beta, psi = fit_est$psi)
  ))
)
solve(S_eta) %*% fit_est$beta %*% S_eta
#>           [,1] [,2]
#> [1,] 0.0000000    0
#> [2,] 0.4604657    0
# Compare to lavaan
standardizedSolution(fit)[8, ]
#>     lhs op   rhs est.std  se     z pvalue ci.lower ci.upper
#> 8 dem60  ~ ind60    0.46 0.1 4.593      0    0.264    0.657
```

## Multiple Groups

``` r

reg <- '
  # latent variable definitions
    visual =~ x1 + x2 + x3
    speed =~ x7 + x8 + x9

  # regressions
    visual ~ c(b1, b1) * speed
'
reg_fit <- sem(reg, data = HolzingerSwineford1939,
               group = "school",
               group.equal = c("loadings", "intercepts"))
```

### Separate standardization by group

``` r

standardizedSolution(reg_fit, type = "std.lv") |>
  subset(subset = label == "b1")
#>       lhs op   rhs group label est.std    se     z pvalue ci.lower ci.upper
#> 7  visual  ~ speed     1    b1   0.369 0.073 5.032      0    0.226    0.513
#> 30 visual  ~ speed     2    b1   0.496 0.086 5.768      0    0.327    0.664
# Compare to doing it by hand
reg_fit_est <- lavInspect(reg_fit, what = "est")
# Group 1:
S_eta1 <- diag(
  sqrt(diag(
    R2spa:::veta(reg_fit_est[[1]]$beta, psi = reg_fit_est[[1]]$psi)
  ))
)
solve(S_eta1) %*% reg_fit_est[[1]]$beta %*% S_eta1
#>      [,1]      [,2]
#> [1,]    0 0.3694845
#> [2,]    0 0.0000000
# Group 2:
S_eta2 <- diag(
  sqrt(diag(
    R2spa:::veta(reg_fit_est[[2]]$beta, psi = reg_fit_est[[2]]$psi)
  ))
)
solve(S_eta2) %*% reg_fit_est[[2]]$beta %*% S_eta2
#>      [,1]      [,2]
#> [1,]    0 0.4958876
#> [2,]    0 0.0000000
```

### Grand standardization

For each group we have the group-specific covariance matrix
$`\mathrm{Var}(\mathbf{\eta}_g)`$. The group-specific latent means are

``` math
\mathrm{E}(\mathbf{\eta}_g) = (\mathbf{I} - \mathbf{B})^{-1}[\mathbf{\alpha} + \Gamma \mathrm{E}(\mathbf{X})]
```

The grand mean is
$`\mathrm{E}(\mathbf{\eta}) = \sum_{g = 1}^G n_g \mathrm{E}(\mathbf{\eta}_g) / N`$

The grand covariance matrix is:

``` math
\mathrm{Var}(\mathbf{\eta}) = \frac{1}{N} \sum_{g = 1}^G n_g \left\{\mathrm{Var}(\mathbf{\eta}_g) + [\mathrm{E}(\mathbf{\eta}_g) - \mathrm{E}(\mathbf{\eta})][\mathrm{E}(\mathbf{\eta}_g) - \mathrm{E}(\mathbf{\eta})]^\top \right\}
```

So we need to involve the mean as well

``` r

# Latent means
lavInspect(reg_fit, what = "mean.lv")  # lavaan
#> $Pasteur
#> visual  speed 
#>      0      0 
#> 
#> $`Grant-White`
#> visual  speed 
#> -0.204 -0.167
R2spa:::eeta(reg_fit_est[[1]]$beta, alpha = reg_fit_est[[1]]$alpha)
#>        intercept
#> visual         0
#> speed          0
R2spa:::eeta(reg_fit_est[[2]]$beta, alpha = reg_fit_est[[2]]$alpha)
#>         intercept
#> visual -0.2036667
#> speed  -0.1666424
# Grand covariance
ns <- lavInspect(reg_fit, what = "nobs")
R2spa:::veta_grand(ns, beta_list = lapply(reg_fit_est, `[[`, "beta"),
                   psi_list = lapply(reg_fit_est, `[[`, "psi"),
                   alpha_list = lapply(reg_fit_est, `[[`, "alpha"))
#>           visual     speed
#> visual 0.5497595 0.2085522
#> speed  0.2085522 0.4059395
```

The function
[`grandStandardizedSolution()`](https://mmm-lab-um.github.io/R2spa/reference/grand_standardized_solution.md)
automates the computation of grand standardized coefficients, with the
asymptotic standard error obtained using the delta method:

``` r

grandStandardizedSolution(reg_fit)
#>       lhs op   rhs exo group block label est.std    se     z pvalue ci.lower
#> 7  visual  ~ speed   0     1     1    b1   0.431 0.073 5.867      0    0.287
#> 30 visual  ~ speed   0     2     2    b1   0.431 0.073 5.867      0    0.287
#>    ci.upper
#> 7     0.575
#> 30    0.575
```
