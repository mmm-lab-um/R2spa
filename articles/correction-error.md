# Correction to Measurement Error

When computing the reliability or the standard errors of factor scores,
we typically condition on the estimates of the model parameters, aka
assuming that they are known. However, when sample size is small,
ignoring those can be problematic as the model parameters are imprecise.
This vignette provides details on how to obtain more realistic estimates
of standard errors of factor scores.

In short: by default
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
computes score error variances *conditional on the estimated* parameters
— i.e. it treats the scoring weights as known. That is the within-sample
error variance, but in small samples, or weakly-identified models (where
the weights are estimated with substantial error), it understates the
score’s error variance. `get_fs(..., corrected_fsT = TRUE)` adds the
delta-method correction term for the sampling error of the scoring
weights, inflating `fsT` accordingly. The effect is negligible in large,
well-identified samples (the weight-sampling variance tends to zero) and
grows as identification gets poorer. Scope/constraints: `corrected_fsT`
applies to `lavaan` regression / Bartlett scores; it is ignored for
`merMod` objects and is not supported when the data contain missing
values.

## Unidimensional Case

Consider the conditional variance of the factor scores,
$`Var(\tilde \eta \mid \eta)`$, with
$`\tilde \eta = \boldsymbol{a}^\top \boldsymbol{y}`$. From the factor
model, we have

``` math
Var(\tilde \eta \mid \eta) = Var(\boldsymbol{a}^\top \boldsymbol{\lambda} \eta + \boldsymbol{a}^\top \boldsymbol{\varepsilon} \mid \eta) = \eta^2 Var(\boldsymbol{a}^\top \boldsymbol{\lambda}) + Var(\boldsymbol{a}^\top \boldsymbol{\varepsilon}).
```

The second term is the error component. Here is a small simulation:

``` r

#' Load package and set seed
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
library(R2spa)
set.seed(191254)

#' Fixed conditions
num_obs <- 100
lambda <- c(.3, .5, .7, .9)
theta <- c(.5, .4, .5, .7)

#' Simulate true score, and standardized
eta <- rnorm(num_obs)
eta <- (eta - mean(eta))
eta <- eta / sqrt(mean(eta^2))

#' Function for simulating y
simy <- function(eta, lambda) {
    t(tcrossprod(lambda, eta) +
          rnorm(length(lambda) * length(eta), sd = sqrt(theta)))
}
```

``` r

#' Simulation
nsim <- 2500
tfsy_sim <- fsy_sim <- matrix(NA, nrow = num_obs, ncol = nsim)
# Also save the scoring matrix
a_sim <- matrix(NA, nrow = length(lambda), ncol = nsim)
for (i in seq_len(nsim)) {
    y <- simy(eta = eta, lambda = lambda)
    tfsy_sim[, i] <- R2spa::compute_fscore(
      y, lambda = lambda, theta = diag(theta), psi = matrix(1)
    )
    fsy <- R2spa::get_fs(y, std.lv = TRUE)
    fsy_sim[, i] <- fsy[, 1]
    a_sim[, i] <- attr(fsy, which = "scoring_matrix")[[1]]
}
```

``` r

#' Average conditional variance
# a known
apply(tfsy_sim, 1, var) |> mean()
#> [1] 0.1871828
# compare to theoretical value
true_a <- R2spa:::compute_a_reg(lambda, psi = matrix(1), theta = diag(theta))
true_a %*% diag(theta) %*% t(true_a)
#>           [,1]
#> [1,] 0.1893211
# a estimated
apply(fsy_sim, 1, var) |> mean() 
#> [1] 0.2078126
```

The standard errors are larger when using sample estimates of weights.

When $`\boldsymbol{a}`$ is known, the error variance is simply
$`\boldsymbol{a}^\top \boldsymbol{\Theta} \boldsymbol{a}`$. However,
when $`\boldsymbol{a}`$ is unknown and estimated as
$`\tilde{\boldsymbol{a}}`$ from the data,

``` math
Var(\tilde{\boldsymbol{a}}^\top \boldsymbol{\varepsilon}) = E[Var(\tilde{\boldsymbol{a}}^\top \boldsymbol{\varepsilon} \mid \tilde{\boldsymbol{a}})] + Var[E(\tilde{\boldsymbol{a}}^\top \boldsymbol{\varepsilon} \mid \tilde{\boldsymbol{a}})].
```

Under the assumption of $`E(\boldsymbol{\varepsilon})`$ =
$`\boldsymbol{0}`$, the second term is 0. The first term is an
expectation of a quadratic form,

``` math
E(\tilde{\boldsymbol{a}}^\top \boldsymbol{\Theta} \tilde{\boldsymbol{a}}) = E(\tilde{\boldsymbol{a}}^\top) \boldsymbol{\Theta} E(\tilde{\boldsymbol{a}}) + \mathrm{tr}(\boldsymbol{\Theta} \boldsymbol{V}_{\tilde{\boldsymbol{a}}})
```

where $`\boldsymbol{V}_{\tilde{\boldsymbol{a}}}`$ is the covariance
matrix of $`\tilde{\boldsymbol{a}}`$. We can see this using the
simulated data

``` r

correct_fac <- sum(diag(diag(theta) %*% cov(t(a_sim))))
true_a %*% diag(theta) %*% t(true_a) + correct_fac
#>           [,1]
#> [1,] 0.2100938
```

### Correction Factor

We can approximate $`\boldsymbol{V}_{\tilde{\boldsymbol{a}}}`$ using the
[delta method](https://en.wikipedia.org/wiki/Delta_method), and the
other terms in the above using sample estimates

``` math
\tilde{\boldsymbol{a}}^\top \hat{\boldsymbol{\Theta}} \tilde{\boldsymbol{a}} + \mathrm{tr}(\hat{\boldsymbol{\Theta}} \hat{\boldsymbol{V}}_{\tilde{\boldsymbol{a}}}).
```

The delta method estimate,
$`\hat{\boldsymbol{V}}_{\tilde{\boldsymbol{a}}}`$, is

``` math
\boldsymbol{J}_{\tilde{\boldsymbol{a}}}(\boldsymbol{\theta})^\top \hat{\boldsymbol{V}}_{\tilde{\boldsymbol{\theta}}} \boldsymbol{J}_{\tilde{\boldsymbol{a}}}(\boldsymbol{\theta}),
```

where $`\boldsymbol{J}_{\tilde{\boldsymbol{a}}}(\boldsymbol{\theta})`$
is the Jacobian matrix, and
$`\hat{\boldsymbol{V}}_{\tilde{\boldsymbol{\theta}}}`$ is the sampling
variance of the sample estimator of $`\boldsymbol{\theta}`$. The
Jacobian can be approximated using finite difference with
[`lavaan::lav_func_jacobian_complex()`](https://rdrr.io/pkg/lavaan/man/lav_func.html).

``` r

y <- simy(eta = eta, lambda = lambda)
cfa_fit <- cfa("f =~ y1 + y2 + y3 + y4",
               data = setNames(data.frame(y), paste0("y", 1:4)),
               std.lv = TRUE)
# Jacobian
R2spa:::compute_a(coef(cfa_fit), lavobj = cfa_fit)
#> [[1]]
#>         [,1]      [,2]      [,3]      [,4]
#> f 0.09890619 0.3246767 0.4500339 0.3134454
jac_a <- lavaan::lav_func_jacobian_complex(
    \(x) R2spa:::compute_a(x, lavobj = cfa_fit)[[1]],
    coef(cfa_fit)
)
sum(diag(lavInspect(cfa_fit, what = "est")$theta %*%
             jac_a %*% vcov(cfa_fit) %*% t(jac_a)))
#> [1] 0.02152978
fsy <- R2spa::get_fs(y, std.lv = TRUE)
attr(fsy, which = "fsT")
#> [[1]]
#>           fs_f1
#> fs_f1 0.2014365
# Using `get_fs(..., corrected_fsT = TRUE)
fsy2 <- R2spa::get_fs(y, std.lv = TRUE, corrected_fsT = TRUE)
attr(fsy2, which = "fsT")
#> [[1]]
#>           fs_f1
#> fs_f1 0.2229663
```

## Multidimensional Case

``` r

#' Set seed
set.seed(251329)

#' Fixed conditions
num_obs <- 100
lambda <- Matrix::bdiag(list(c(.3, .5, .7, .9), c(.7, .6, .7))) |>
    as.matrix()
theta <- c(.5, .4, .5, .7, .7, .8, .5)
psi <- matrix(c(1, -.4, -.4, 1), nrow = 2)

#' Simulate true score (exact means and covariances)
eta <- MASS::mvrnorm(num_obs, mu = rep(0, 2), Sigma = psi, empirical = TRUE)

#' Function for simulating y
simy <- function(eta, lambda) {
    tcrossprod(eta, lambda) +
        MASS::mvrnorm(num_obs, mu = rep(0, length(theta)), Sigma = diag(theta))
}
```

``` r

#' Simulation
nsim <- 2500
tfsy_sim <- fsy_sim <- array(NA, dim = c(num_obs, ncol(lambda), nsim))
# Also save the scoring matrix
a_sim <- array(NA, dim = c(ncol(lambda), nrow(lambda), nsim))
for (i in seq_len(nsim)) {
    y <- simy(eta = eta, lambda = lambda)
    tfsy_sim[, , i] <- R2spa::compute_fscore(
      y, lambda = lambda, theta = diag(theta), psi = psi
    )
    fsy <- R2spa::get_fs(
        data.frame(y) |> setNames(paste0("y", 1:7)),
        model = "f1 =~ y1 + y2 + y3 + y4\nf2 =~ y5 + y6 + y7",
        std.lv = TRUE)
    fsy_sim[, , i] <- as.matrix(fsy[, 1:2])
    a_sim[, , i] <- attr(fsy, which = "scoring_matrix")[[1]]
}
```

``` r

#' Average conditional variance
# a known
apply(tfsy_sim, MARGIN = 1, FUN = \(x) cov(t(x))) |> 
    rowMeans() |> matrix(ncol = 2, nrow = 2)
#>            [,1]       [,2]
#> [1,]  0.1789934 -0.0483855
#> [2,] -0.0483855  0.2022312
# compare to theoretical value
true_a <- R2spa:::compute_a_reg(lambda, psi = psi, theta = diag(theta))
true_a %*% diag(theta) %*% t(true_a)
#>             [,1]        [,2]
#> [1,]  0.18076106 -0.04855749
#> [2,] -0.04855749  0.20339709
# a estimated
apply(fsy_sim, MARGIN = 1, FUN = \(x) cov(t(x))) |> 
    rowMeans() |> matrix(ncol = 2, nrow = 2)
#>             [,1]        [,2]
#> [1,]  0.20175239 -0.05254471
#> [2,] -0.05254471  0.23427145
```

``` r

correct_fac11 <- sum(diag(diag(theta) %*% cov(t(a_sim[1, , ]))))
correct_fac22 <- sum(diag(diag(theta) %*% cov(t(a_sim[2, , ]))))
correct_fac21 <- sum(diag(diag(theta) %*%
                              cov(t(a_sim[2, , ]), t(a_sim[1, , ]))))
correct_fac <- matrix(c(correct_fac11, correct_fac21,
                        correct_fac21, correct_fac22), nrow = 2)
true_a %*% diag(theta) %*% t(true_a) + correct_fac
#>             [,1]        [,2]
#> [1,]  0.20114155 -0.05375924
#> [2,] -0.05375924  0.23429418
```

The $`i`$, $`j`$ element of the correction matrix can be approximated by
$`\mathrm{tr}(\hat{\boldsymbol{\Theta}} \hat{\boldsymbol{V}}_{{\tilde{\boldsymbol{a}}_i}{\tilde{\boldsymbol{a}}_j}})`$,
with
$`\hat{\boldsymbol{V}}_{{\tilde{\boldsymbol{a}}_i}{\tilde{\boldsymbol{a}}_j}})`$
obtained by

``` math
\boldsymbol{J}_{\tilde{\boldsymbol{a}}_i}(\boldsymbol{\theta})^\top \hat{\boldsymbol{V}}_{\tilde{\boldsymbol{\theta}}} \boldsymbol{J}_{\tilde{\boldsymbol{a}}_j}(\boldsymbol{\theta}),
```

``` r

y <- simy(eta = eta, lambda = lambda)
cfa_fit <- cfa("f1 =~ y1 + y2 + y3 + y4\nf2 =~ y5 + y6 + y7",
               data = setNames(data.frame(y), paste0("y", 1:7)),
               std.lv = TRUE)
# Jacobian
R2spa:::compute_a(coef(cfa_fit), lavobj = cfa_fit)
#> [[1]]
#>           [,1]        [,2]        [,3]        [,4]        [,5]        [,6]
#> f1  0.06546889  0.22696068  0.38977223  0.33101521 -0.05012587 -0.02226577
#> f2 -0.01227313 -0.04254719 -0.07306865 -0.06205377  0.37397828  0.16612010
#>           [,7]
#> f1 -0.05314801
#> f2  0.39652577
jac_a1 <- lavaan::lav_func_jacobian_complex(
    \(x) R2spa:::compute_a(x, lavobj = cfa_fit)[[1]][1, ],
    coef(cfa_fit)
)
jac_a2 <- lavaan::lav_func_jacobian_complex(
    \(x) R2spa:::compute_a(x, lavobj = cfa_fit)[[1]][2, ],
    coef(cfa_fit)
)
# Correction[1, 1]
sum(diag(lavInspect(cfa_fit, what = "est")$theta %*%
             jac_a1 %*% vcov(cfa_fit) %*% t(jac_a1)))
#> [1] 0.01551285
# Correction[2, 2]
sum(diag(lavInspect(cfa_fit, what = "est")$theta %*%
             jac_a2 %*% vcov(cfa_fit) %*% t(jac_a2)))
#> [1] 0.02207038
# Correction[2, 1]
sum(diag(lavInspect(cfa_fit, what = "est")$theta %*%
             jac_a2 %*% vcov(cfa_fit) %*% t(jac_a1)))
#> [1] -0.004873356
fsy <- R2spa::get_fs(
    data.frame(y) |> setNames(paste0("y", 1:7)),
    model = "f1 =~ y1 + y2 + y3 + y4\nf2 =~ y5 + y6 + y7",
    std.lv = TRUE)
attr(fsy, which = "fsT")
#> [[1]]
#>             fs_f1       fs_f2
#> fs_f1  0.16691097 -0.05653148
#> fs_f2 -0.05653148  0.19891917
# Using `get_fs(..., corrected_fsT = TRUE)
fsy2 <- R2spa::get_fs(
    data.frame(y) |> setNames(paste0("y", 1:7)),
    model = "f1 =~ y1 + y2 + y3 + y4\nf2 =~ y5 + y6 + y7",
    std.lv = TRUE,
    corrected_fsT = TRUE)
attr(fsy2, which = "fsT")
#> [[1]]
#>             fs_f1       fs_f2
#> fs_f1  0.18242382 -0.06140484
#> fs_f2 -0.06140484  0.22098955
```
