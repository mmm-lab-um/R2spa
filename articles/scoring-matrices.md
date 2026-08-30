# Scoring Matrices: lavaan CFA and lme4

## Overview

The **scoring matrix** is the linear operator that converts the modelled
observed data into the factor / empirical-Bayes (EB) scores reported by
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md).
Together with the score error-covariance (`fsT`), implied loadings
(`fsL`) and implied intercepts (`fsb`), it fully specifies how each
reported score is constructed from the model.
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
exposes the scoring matrix as the `scoring_matrix` attribute alongside
these three:

- for **lavaan** CFA models it is a named list of one (score x item)
  matrix per group;
- for **lme4** models it is a named list of one (random-effect x
  observation) matrix per cluster.

In this vignette we reconstruct the scores **by hand** from the stored
scoring matrix and the modelled data, demonstrating the identity in both
settings. (The vignette [2S-PA with Random
Effects](https://mmm-lab-um.github.io/R2spa/articles/multilevel.md)
shows the downstream 2S-PA use of the `ev_`/`ecov_` columns that
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
derives from `fsT`.)

``` r

library(R2spa)
library(lavaan)
library(lme4)
```

## lavaan CFA

We fit a single-factor CFA on
[`lavaan::PoliticalDemocracy`](https://rdrr.io/pkg/lavaan/man/PoliticalDemocracy.html)
with a mean structure, so that the observed-item intercepts `nu` and the
factor intercept `alpha` are estimated parameters:

``` r

cfa_fit <- cfa("f =~ y1 + y2 + y3 + y4",
               data = lavaan::PoliticalDemocracy,
               meanstructure = TRUE)
fs <- get_fs(cfa_fit)
head(fs)
#>         fs_f   fs_f_se f_by_fs_f   ev_fs_f
#> 1 -2.7487224 0.6756472 0.8868049 0.4564991
#> 2 -3.0360803 0.6756472 0.8868049 0.4564991
#> 3  2.6718589 0.6756472 0.8868049 0.4564991
#> 4  2.9936997 0.6756472 0.8868049 0.4564991
#> 5  1.9242932 0.6756472 0.8868049 0.4564991
#> 6  0.9922798 0.6756472 0.8868049 0.4564991
```

The `scoring_matrix` attribute is a named list keyed by group label (an
empty string for single-group fits); each entry is the matrix `S` with
one row per factor and one column per item:

``` r

sm <- attr(fs, "scoring_matrix")
S <- sm[[1]]
S
#>        [,1]     [,2]      [,3]      [,4]
#> f 0.2298863 0.112681 0.1071893 0.2788005
# Same results as in lavaan
attr(lavPredict(cfa_fit, fsm = TRUE), "fsm")
#> [[1]]
#>          y1       y2        y3        y4
#> f 0.2298863 0.112681 0.1071893 0.2788005
```

### Reconstructing the scores by hand

[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
computes the regression scores as

`f = S %*% (y' - (lambda %*% alpha + nu))' + alpha`,

where `lambda` is the item-loading matrix, `nu` the item intercepts, and
`alpha` the factor intercept. From `lavInspect(cfa_fit, "est")` we can
therefore recompute the scores directly from the raw data:

``` r

est <- lavInspect(cfa_fit, "est")
Lambda <- as.matrix(est$lambda)
Theta  <- as.matrix(est$theta)
nu     <- as.numeric(est$nu)
alpha  <- as.numeric(est$alpha)
y      <- lavInspect(cfa_fit, "data")
fs_hand <- S %*% (t(y) - (as.vector(Lambda %*% alpha) + nu)) + as.vector(alpha)
max(abs(fs_hand - fs$fs_f))
#> [1] 0
```

The maximum deviation from the reported `fs_f` column is at the level of
machine precision. The same `S` also generates the companion attributes,
so all of them are redundant views of the same linear operator:

- `fsL = S %*% Lambda`: implied loadings of the items on the factor
  scores (i.e. the factor score’s regression on the items, in loading
  form);
- `fsb = alpha - fsL %*% alpha`: implied intercepts of the factor
  scores;
- `fsT = S %*% Theta %*% t(S)`: error variance-covariance of the
  regression scores.

``` r

max(abs(attr(fs, "fsL")[[1]] - S %*% Lambda))
#> [1] 0
max(abs(attr(fs, "fsb")[[1]] - (alpha - attr(fs, "fsL")[[1]] %*% alpha)))
#> [1] 0
max(abs(attr(fs, "fsT")[[1]] - S %*% Theta %*% t(S)))
#> [1] 0
```

## lme4

We fit the classic `sleepstudy` growth model with a random intercept and
slope for `Subject`, using ML estimation (the 2S-PA convention; see the
Notes there). Each cluster (subject) obtains a single 2-dimensional EB
estimate of the random effect:

``` r

lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
             REML = FALSE)
fs1 <- get_fs(lme1)
head(fs1)
#>        fs_u0      fs_u1 fs_u0_se fs_u1_se u0_by_fs_u0 u0_by_fs_u1 u1_by_fs_u0
#> 1   2.815789  9.0755068 9.434845 1.858658   0.7366482  0.03978511   0.7195288
#> 2 -40.047855 -8.6441517 9.434845 1.858658   0.7366482  0.03978511   0.7195288
#> 3 -38.432497 -5.5134706 9.434845 1.858658   0.7366482  0.03978511   0.7195288
#> 4  22.831765 -4.6586649 9.434845 1.858658   0.7366482  0.03978511   0.7195288
#> 5  21.549514 -2.9444450 9.434845 1.858658   0.7366482  0.03978511   0.7195288
#> 6   8.815409 -0.2351823 9.434845 1.858658   0.7366482  0.03978511   0.7195288
#>   u1_by_fs_u1 ev_fs_u0 ecov_fs_u1_fs_u0 ev_fs_u1
#> 1   0.8287253  89.0163        -11.46683 3.454609
#> 2   0.8287253  89.0163        -11.46683 3.454609
#> 3   0.8287253  89.0163        -11.46683 3.454609
#> 4   0.8287253  89.0163        -11.46683 3.454609
#> 5   0.8287253  89.0163        -11.46683 3.454609
#> 6   0.8287253  89.0163        -11.46683 3.454609
```

Its `scoring_matrix` attribute is a named list with one entry per
cluster; each entry `S_j` has one row per random effect and one column
per observation in the cluster:

``` r

sm1 <- attr(fs1, "scoring_matrix")
length(sm1)
#> [1] 18
names(sm1)
#>  [1] "308" "309" "310" "330" "331" "332" "333" "334" "335" "337" "349" "350"
#> [13] "351" "352" "369" "370" "371" "372"
knitr::kable(sm1[["308"]], digits = 3,
             caption = "Scoring matrix for subject 308")
```

|       |      1 |      2 |      3 |      4 |     5 |     6 |     7 |      8 |      9 |     10 |
|:------|-------:|-------:|-------:|-------:|------:|------:|------:|-------:|-------:|-------:|
| fs_u0 |  0.215 |  0.184 |  0.152 |  0.121 | 0.089 | 0.058 | 0.026 | -0.005 | -0.036 | -0.068 |
| fs_u1 | -0.031 | -0.024 | -0.016 | -0.008 | 0.000 | 0.008 | 0.016 |  0.024 |  0.032 |  0.039 |

Scoring matrix for subject 308 {.table style="width:100%;"}

### Reconstructing the scores by hand

The theoretical EB score of cluster `j` is

`u_hat_j = S_j %*% (y_j - X_j %*% beta)`,

with

`S_j = (G^{-1} + Z_j' Z_j / sigma^2)^{-1} Z_j' / sigma^2`,

where `y_j` is the cluster’s response vector, `X_j` its rows of the
fixed-effects design, `Z_j` its rows of the random-effects design,
`beta = fixef(.)` the fixed-effect estimates, and `G`, `sigma^2` the
random-effects and residual variances. `R2spa` reconstructs the scaled
matrix `D = G / sigma^2` from the LME4 theta via `get_D()`; because
`Kz = Z_j' Z_j` is unscaled, `sigma^2` cancels in `S_j`. As with lavaan,
`S_j` also determines the companion attributes:

- `fsL_j = S_j %*% Z_j`: loading of the EB estimate on the random
  effect, i.e. the operator of `E[u_hat_j | u_j]`;
- `fsT_j = sigma^2 %*% S_j %*% t(S_j)`: conditional (measurement-error)
  covariance of the EB scores, `Var(u_hat_j | u_j)`.

Reconstructing the scores from the raw data per cluster:

``` r

mf   <- model.frame(lme1)
y    <- model.response(mf)
X    <- as.matrix(lme1@pp$X)
beta <- fixef(lme1)
idx_list <- split(seq_len(nrow(mf)), lme1@flist[[1]])
fs_hand1 <- t(vapply(names(sm1), function(lv) {
  idx <- idx_list[[lv]]
  sm1[[lv]] %*% (y[idx] - as.numeric(X[idx, ] %*% beta))
}, numeric(2)))
max(abs(fs_hand1 - as.matrix(fs1[, 1:2])))
#> [1] 1.350031e-13
max(abs(fs_hand1 - as.matrix(ranef(lme1)[[1]])))
#> [1] 1.350031e-13
```

Both the
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
output and
[`lme4::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) are
recovered to machine precision.

Under `method = "ML"` the same per-cluster identity holds with
prior-free OLS weights `S_j = (Z_j' Z_j)^+ Z_j'` (no random-effects
prior enters `S_j`), so the reconstructed scores are the per-cluster OLS
(ML) estimates of the random effects rather than the EB estimates shown
here.

### Notes

- Only the **first random-effect term** is scored (`[[1]]` convention,
  matching `ranef(lme1)[[1]]`).
- Cluster sizes may differ (unbalanced designs); the per-cluster list
  accommodates differing matrix widths.
- The fixed and random designs need not coincide (`Z != X`); the scoring
  identity is based on the random design `Z`.
- The identity holds under both ML and REML estimation; the vignette
  uses ML by 2S-PA convention.

## Comparison

``` r

knitr::kable(
  rbind(
    c("scoring matrix `S`",
      "regression (default):
      S = Psi Lambda' (Lambda Psi Lambda' + Theta)^{-1};
      \"Bartlett\":
      S = (Lambda' Theta^{-1} Lambda)^{-1} Lambda' Theta^{-1}",
      "S_j = (G^{-1} + Z_j' Z_j / sigma^2)^{-1} Z_j' / sigma^2 (EB
      weights)"),
    c("shape & orientation",
      "score x item (q x p), shared within a group",
      "num_re x n_j per cluster j (p x n_j)"),
    c("container",
      "named list of one matrix per group",
      "named list of one matrix per cluster
      (list, not array: n_j can differ)"),
    c("score reconstruction",
      "f = S %*% (y' - (Lambda %*% alpha + nu))' + alpha
      (all cases in one product)",
      "u_hat_j = S_j %*% (y_j - X_j %*% beta)
      (one product per cluster)")
  ),
  col.names = c("property", "lavaan CFA", "lme4"),
  escape = FALSE
)
```

| property           | lavaan CFA            | lme4 |
|:-------------------|:----------------------|:-----|
| scoring matrix `S` | regression (default): |      |

      S = Psi Lambda' (Lambda Psi Lambda' + Theta)^{-1};
      "Bartlett":
      S = (Lambda' Theta^{-1} Lambda)^{-1} Lambda' Theta^{-1} |S_j = (G^{-1} + Z_j' Z_j / sigma^2)^{-1} Z_j' / sigma^2 (EB
      weights)   |

\|shape & orientation \|score x item (q x p), shared within a group
\|num_re x n_j per cluster j (p x n_j) \| \|container \|named list of
one matrix per group \|named list of one matrix per cluster (list, not
array: n_j can differ) \| \|score reconstruction \|f = S %*% (y’ -
(Lambda %*% alpha + nu))’ + alpha (all cases in one product) \|u_hat_j =
S_j %*% (y_j - X_j %*% beta) (one product per cluster) \|
