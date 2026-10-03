
# R2spa

<!-- badges: start -->

[![R-CMD-check](https://github.com/mmm-lab-um/R2spa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mmm-lab-um/R2spa/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`R2spa` is a free and open-source R package that performs two-stage path
analysis (2S-PA). With 2S-PA, researchers can perform path analysis by
first obtaining factor scores and then adjusting for measurement errors
using estimates of observation-specific reliability or standard error of
those factor scores. As a viable alternative to SEM, 2S-PA has been
shown to give equally-good estimates as SEM in relatively simple models
and large sample sizes, as well as to give more accurate parameter
estimates, has better control of Type I error rates, and has
substantially less convergence problems in more complex models or small
sample sizes.

## Installation

This package is still in developmental stage and can be installed on
GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("mmm-lab-um/R2spa")
```

## Documentation

The package website, <https://mmm-lab-um.github.io/R2spa/>, hosts the
full reference and a series of worked articles. Start with the
[Two-Stage Path Analysis (2S-PA) Model
Examples](https://mmm-lab-um.github.io/R2spa/articles/R2spa.html)
getting-started article.

## Example

The canonical workflow scores each latent from one multi-factor
measurement model (`get_fs(..., local = TRUE)`), then fits the
structural path with `tspa()`, which reads each score’s standard error
from the `fs_<name>_se` columns automatically (no hard-coded `se_fs`
needed):

``` r
library(lavaan)
library(R2spa)

# One multi-factor measurement model
meas_model <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
'
```

``` r
# Stage 1: score each latent from its own measurement model (local = TRUE)
fs_dat <- get_fs(PoliticalDemocracy, model = meas_model, local = TRUE)
# get_fs() gives a data frame of factor scores and their standard errors
head(fs_dat)
#>     fs_ind60   fs_dem60  fs_dem65 fs_ind60_se fs_dem60_se fs_dem65_se
#> 1 -0.5261683 -2.7487224 -1.371719   0.1213615   0.6756472   0.5724405
#> 2  0.1436527 -3.0360803 -0.950851   0.1213615   0.6756472   0.5724405
#> 3  0.7143559  2.6718589  2.738012   0.1213615   0.6756472   0.5724405
#> 4  1.2399257  2.9936997  1.785091   0.1213615   0.6756472   0.5724405
#> 5  0.8319080  1.9242932  1.544704   0.1213615   0.6756472   0.5724405
#> 6  0.2123845  0.9922798 -1.050841   0.1213615   0.6756472   0.5724405
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
#> 1         0.8998252  0.01472862                      0   0.4564991
#> 2         0.8998252  0.01472862                      0   0.4564991
#> 3         0.8998252  0.01472862                      0   0.4564991
#> 4         0.8998252  0.01472862                      0   0.4564991
#> 5         0.8998252  0.01472862                      0   0.4564991
#> 6         0.8998252  0.01472862                      0   0.4564991
#>   ecov_fs_dem65_fs_ind60 ecov_fs_dem65_fs_dem60 ev_fs_dem65
#> 1                      0                      0   0.3276882
#> 2                      0                      0   0.3276882
#> 3                      0                      0   0.3276882
#> 4                      0                      0   0.3276882
#> 5                      0                      0   0.3276882
#> 6                      0                      0   0.3276882
```

``` r
# Stage 2: fit the structural path; se_fs is read from the fs_*_se columns
tspa_fit <- tspa(
  model = "dem60 ~ ind60
          dem65 ~ ind60 + dem60",
  data = fs_dat
)
```

Because the latent constructs have no intrinsic scale, the
**standardized** coefficients are usually the parameters of substantive
interest (see the *Two-Stage Path Analysis* vignette,
`vignette("R2spa")`, for the reasoning and references):

``` r
standardizedSolution(tspa_fit) |>
  subset(op %in% c("~", "~~") & !grepl("^fs_", lhs))
#>      lhs op   rhs est.std    se      z pvalue ci.lower ci.upper
#> 19 dem60  ~ ind60   0.453 0.101  4.480  0.000    0.255    0.651
#> 20 dem65  ~ ind60   0.129 0.073  1.771  0.076   -0.014    0.272
#> 21 dem65  ~ dem60   0.898 0.049 18.314  0.000    0.802    0.994
#> 22 ind60 ~~ ind60   1.000 0.000     NA     NA    1.000    1.000
#> 23 dem60 ~~ dem60   0.795 0.092  8.668  0.000    0.615    0.974
#> 24 dem65 ~~ dem65   0.071 0.047  1.514  0.130   -0.021    0.164
```

> **Optional dependencies.** The exact (non-pooled) stage-2 route,
> `tspa_mx_model()`, requires `OpenMx`, and the item-response-theory
> (`mirt`) scoring path requires `mirt` (both `Suggests`); see the
> *2S-PA with OpenMx and IRT (mirt)* vignette,
> `vignette("tspa-vignette-mx")`.

This package is based upon work supported by the National Science
Foundation under Grant No. 2141790.

<!-- `devtools::build_readme()` -->
