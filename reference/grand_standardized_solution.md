# Grand Standardized Solution

Grand standardized solution of a two-stage path analysis model.

## Usage

``` r
grand_standardized_solution(
  object,
  model_list = NULL,
  se = TRUE,
  acov_par = NULL,
  free_list = NULL,
  level = 0.95
)

grandStandardizedSolution(
  object,
  model_list = NULL,
  se = TRUE,
  acov_par = NULL,
  free_list = NULL,
  level = 0.95
)
```

## Arguments

- object:

  An object of class lavaan.

- model_list:

  A list of string variable describing the structural path model, in
  `lavaan` syntax.

- se:

  A Boolean variable. If TRUE, standard errors for the grand
  standardized parameters will be computed.

- acov_par:

  An asymptotic variance-covariance matrix for a fitted model object;
  defaults to `vcov(object)`. Supplying the covariance of an
  SE-corrected
  [`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) fit
  (produced with `tspa(corrected_se = TRUE)`) carries the first-order
  correction through to the grand-standardized standard errors reported
  by this function.

- free_list:

  A list of model matrices that indicate the position of the free
  parameters in the parameter vector.

- level:

  The confidence level required.

## Value

A data frame (class `lavaan.data.frame`) with one row per structural
(`~`) parameter in partable order: `lhs`, `op`, `rhs`, `exo`, `group`,
`block`, `label`, the grand standardized estimate `est.std`, and, when
`se = TRUE`, `se`, `z`, `pvalue`, `ci.lower`, and `ci.upper`.

## Details

For a **multigroup** fit the grand standardized solution is computed
from the grand (pooled) means and variance-covariances and is **not**
equivalent to
[`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html),
which matches it only in the single-group case (the function emits that
message when it applies). When the input fit is an SE-corrected
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) fit
(produced with `tspa(corrected_se = TRUE)`, attribute
`tspa_corrected = TRUE`) — or when `acov_par = vcov(corrected_fit)` is
supplied — the reported standard errors are the first-order corrected
grand-standardized standard errors, while the point estimates
(`est.std`) are unchanged.

Structural slopes may be free or user-fixed. A user-fixed slope is
reported alongside the free ones: its `est.std` is the user value
rescaled by the (grand) SD ratio, and its `se` is the first-order delta
approximation, both matching
[`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
(which also reports a delta SE for fixed slopes). A structural
regression whose outcome or predictor is not part of the beta matrix
(e.g. observed-on-observed) is still rejected with an error. The
delta-method SEs are propagated over the free `beta`/`psi` (and `alpha`,
single group) elements as marked by the lavaan free-position matrices; a
free parameter whose block matrix carries no mark contributes no
first-order term (in lavaan 0.7-x the free-position matrices do not
mark, e.g., some observed-variable intercept blocks, whose omission only
matters when their sampling uncertainty affects the group means pooled
into the grand covariance).

## See also

- `vignette("Grand Standardized Coefficients", package = "R2spa")` for
  the grand-standardization workflow.

## Examples

``` r
library(lavaan)

## A single-group, two-factor example
mod1 <- '
   # latent variables
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
   # regressions
     dem60 ~ ind60
'
fit1 <- sem(model = mod1,
          data  = PoliticalDemocracy)
grand_standardized_solution(fit1)
#> The grand standardized solution is equivalent to lavaan::standardizedSolution() for a model with a single group.
#>     lhs op   rhs exo group block label est.std  se     z pvalue ci.lower
#> 8 dem60  ~ ind60   0     1     1          0.46 0.1 4.593      0    0.264
#>   ci.upper
#> 8    0.657

## A single-group, three-factor example
mod2 <- '
    # latent variables
      ind60 =~ x1 + x2 + x3
      dem60 =~ y1 + y2 + y3 + y4
      dem65 =~ y5 + y6 + y7 + y8
    # regressions
      dem60 ~ ind60
      dem65 ~ ind60 + dem60
'
fit2 <- sem(model = mod2,
            data  = PoliticalDemocracy)
grand_standardized_solution(fit2)
#> The grand standardized solution is equivalent to lavaan::standardizedSolution() for a model with a single group.
#>      lhs op   rhs exo group block label est.std    se      z pvalue ci.lower
#> 12 dem60  ~ ind60   0     1     1         0.448 0.102  4.393  0.000    0.248
#> 13 dem65  ~ ind60   0     1     1         0.146 0.070  2.071  0.038    0.008
#> 14 dem65  ~ dem60   0     1     1         0.913 0.048 19.120  0.000    0.819
#>    ci.upper
#> 12    0.648
#> 13    0.283
#> 14    1.006

## A multigroup, two-factor example
mod3 <- '
  # latent variable definitions
    visual =~ x1 + x2 + x3
    speed =~ x7 + x8 + x9
  # regressions
    visual ~ c(b1, b1) * speed
'
fit3 <- sem(mod3, data = HolzingerSwineford1939,
            group = "school",
            group.equal = c("loadings", "intercepts"))
grand_standardized_solution(fit3)
#>       lhs op   rhs exo group block label est.std    se     z pvalue ci.lower
#> 7  visual  ~ speed   0     1     1    b1   0.431 0.073 5.867      0    0.287
#> 30 visual  ~ speed   0     2     2    b1   0.431 0.073 5.867      0    0.287
#>    ci.upper
#> 7     0.575
#> 30    0.575

## A multigroup, three-factor example
mod4 <- '
  # latent variable definitions
    visual =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed =~ x7 + x8 + x9

  # regressions
    visual ~ c(b1, b1) * textual + c(b2, b2) * speed
'
fit4 <- sem(mod4, data = HolzingerSwineford1939,
            group = "school",
            group.equal = c("loadings", "intercepts"))
grand_standardized_solution(fit4)
#>       lhs op     rhs exo group block label est.std    se     z pvalue ci.lower
#> 10 visual  ~ textual   0     1     1    b1   0.419 0.073 5.704      0    0.275
#> 11 visual  ~   speed   0     1     1    b2   0.324 0.078 4.145      0    0.171
#> 46 visual  ~ textual   0     2     2    b1   0.419 0.073 5.704      0    0.275
#> 47 visual  ~   speed   0     2     2    b2   0.324 0.078 4.145      0    0.171
#>    ci.upper
#> 10    0.563
#> 11    0.477
#> 46    0.563
#> 47    0.477
```
