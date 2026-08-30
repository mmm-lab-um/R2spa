# Compute factor scores

Compute factor scores

## Usage

``` r
compute_fscore(
  y,
  lambda,
  theta,
  psi = NULL,
  nu = NULL,
  alpha = NULL,
  method = c("regression", "Bartlett", "mean"),
  center_y = TRUE,
  acov = FALSE,
  fs_matrices = FALSE,
  sum_items = NULL
)
```

## Arguments

- y:

  An N x p matrix where each row is a response vector. If there is only
  one observation, it should be a matrix of one row.

- lambda:

  A p x q matrix of factor loadings.

- theta:

  A p x p matrix of unique variance-covariances.

- psi:

  A q x q matrix of latent factor variance-covariances.

- nu:

  A vector of length p of measurement intercepts.

- alpha:

  A vector of length q of latent means.

- method:

  A character string indicating the method for computing factor scores:
  `"regression"` (default), `"Bartlett"`, or `"mean"` (sum scores, i.e.
  the plain uncentered item means; see
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  for the full description).

- center_y:

  Logical indicating whether `y` should be mean-centered. Default to
  `TRUE`. Ignored for `method = "mean"`, whose scores are raw item
  means.

- acov:

  Logical indicating whether the asymptotic covariance matrix of factor
  scores should be returned as an attribute.

- fs_matrices:

  Logical indicating whether covariances of the error portion of factor
  scores (`fsT`), factor score loading matrix (\\L\\; `fsL`) and
  intercept vector (\\b\\; `fsb`) should be returned. The loading and
  intercept matrices are the implied loadings and intercepts by the
  model when using the factor scores as indicators of the latent
  variables. If `TRUE`, these matrices will be added as attributes.

- sum_items:

  For `method = "mean"` only: a named list mapping factor names to the
  items included in each factor's sum score. `NULL` (default)
  auto-derives the assignment from the loadings. See
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  for the full description.

## Value

An N x p matrix of factor scores.

## Examples

``` r
library(lavaan)
fit <- cfa(" ind60 =~ x1 + x2 + x3
             dem60 =~ y1 + y2 + y3 + y4 ",
           data = PoliticalDemocracy)
fs_lavaan <- lavPredict(fit, method = "Bartlett")
# Using R2spa::compute_fscore()
est <- lavInspect(fit, what = "est")
fs_hand <- compute_fscore(lavInspect(fit, what = "data"),
                          lambda = est$lambda,
                          theta = est$theta,
                          psi = est$psi,
                          method = "Bartlett")
fs_hand - fs_lavaan  # same scores
#>       ind60 dem60
#>  [1,]     0     0
#>  [2,]     0     0
#>  [3,]     0     0
#>  [4,]     0     0
#>  [5,]     0     0
#>  [6,]     0     0
#>  [7,]     0     0
#>  [8,]     0     0
#>  [9,]     0     0
#> [10,]     0     0
#> [11,]     0     0
#> [12,]     0     0
#> [13,]     0     0
#> [14,]     0     0
#> [15,]     0     0
#> [16,]     0     0
#> [17,]     0     0
#> [18,]     0     0
#> [19,]     0     0
#> [20,]     0     0
#> [21,]     0     0
#> [22,]     0     0
#> [23,]     0     0
#> [24,]     0     0
#> [25,]     0     0
#> [26,]     0     0
#> [27,]     0     0
#> [28,]     0     0
#> [29,]     0     0
#> [30,]     0     0
#> [31,]     0     0
#> [32,]     0     0
#> [33,]     0     0
#> [34,]     0     0
#> [35,]     0     0
#> [36,]     0     0
#> [37,]     0     0
#> [38,]     0     0
#> [39,]     0     0
#> [40,]     0     0
#> [41,]     0     0
#> [42,]     0     0
#> [43,]     0     0
#> [44,]     0     0
#> [45,]     0     0
#> [46,]     0     0
#> [47,]     0     0
#> [48,]     0     0
#> [49,]     0     0
#> [50,]     0     0
#> [51,]     0     0
#> [52,]     0     0
#> [53,]     0     0
#> [54,]     0     0
#> [55,]     0     0
#> [56,]     0     0
#> [57,]     0     0
#> [58,]     0     0
#> [59,]     0     0
#> [60,]     0     0
#> [61,]     0     0
#> [62,]     0     0
#> [63,]     0     0
#> [64,]     0     0
#> [65,]     0     0
#> [66,]     0     0
#> [67,]     0     0
#> [68,]     0     0
#> [69,]     0     0
#> [70,]     0     0
#> [71,]     0     0
#> [72,]     0     0
#> [73,]     0     0
#> [74,]     0     0
#> [75,]     0     0
```
