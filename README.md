
# R2spa

<!-- badges: start -->

[![R-CMD-check](https://github.com/mmm-lab-um/R2spa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mmm-lab-um/R2spa/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`R2spa` implements **two-stage path analysis (2S-PA)** for estimating
relations among constructs measured with error. It separates measurement
model estimation from structural model estimation:

1.  **Stage 1:** Fit appropriate measurement models and obtain factor
    scores, together with estimates of their measurement error.
2.  **Stage 2:** Estimate a path model using the factor scores as
    indicators of the latent variables, with constraints that account
    for measurement error.

This differs from using factor scores in an ordinary regression, which
treats them as error-free. It also differs from joint structural
equation modeling (SEM), which estimates the measurement and structural
parameters simultaneously. The two-stage approach lets researchers
evaluate measurement models separately and use the resulting scores in
different structural models.

Simulation studies have found that 2S-PA can perform similarly to joint
SEM, though results depend on the measurement model and analysis
conditions.

## Installation

The package is under development. Install it from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("mmm-lab-um/R2spa")
```

## A lavaan example

The `PoliticalDemocracy` data from lavaan contain indicators of
industrialization in 1960 (`ind60`) and political democracy in 1960 and
1965 (`dem60` and `dem65`). We first obtain a score for each construct,
then estimate the relations among the three constructs.

``` r
library(lavaan)
library(R2spa)

# Measurement models for the three constructs
meas_model <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
'
```

``` r
# Stage 1: fit a separate measurement model for each construct and obtain scores
fs_dat <- get_fs(PoliticalDemocracy, model = meas_model, local = TRUE)
```

``` r
# Stage 2: estimate the structural paths, accounting for score measurement error
tspa_fit <- tspa(
  model = "dem60 ~ ind60
           dem65 ~ ind60 + dem60",
  data = fs_dat
)
```

`local = TRUE` fits a separate single-factor model for each construct.
`tspa()` uses the score measurement information returned by `get_fs()`
and returns a lavaan fit. See the worked vignette for score details and
interpretation.

The standardized path coefficients express associations in
latent-variable standard deviation units. They are useful when comparing
results across different choices of latent-variable scale:

``` r
standardizedSolution(tspa_fit) |>
  subset(op == "~")
#>      lhs op   rhs est.std    se      z pvalue ci.lower ci.upper
#> 19 dem60  ~ ind60   0.453 0.101  4.480  0.000    0.255    0.651
#> 20 dem65  ~ ind60   0.129 0.073  1.771  0.076   -0.014    0.272
#> 21 dem65  ~ dem60   0.898 0.049 18.314  0.000    0.802    0.994
```

## Learn more

Start with the [worked 2S-PA
examples](https://mmm-lab-um.github.io/R2spa/articles/R2spa.html) for a
step-by-step introduction. The [package
website](https://mmm-lab-um.github.io/R2spa/) also includes the function
reference and articles on:

- [Scores from a multi-factor measurement
  model](https://mmm-lab-um.github.io/R2spa/articles/multiple-factors.html),
  including correlated score errors and score cross-loadings.
- [Categorical items and observation-specific measurement
  error](https://mmm-lab-um.github.io/R2spa/articles/tspa-vignette-mx.html),
  using `mirt` for item response models and `tspa_mx_model()` for stage
  2.
- [Missing
  data](https://mmm-lab-um.github.io/R2spa/articles/missing-data.html)
  and [random effects from mixed
  models](https://mmm-lab-um.github.io/R2spa/articles/multilevel.html).
- [Standard errors that account for stage-1 estimation
  uncertainty](https://mmm-lab-um.github.io/R2spa/articles/corrected-se.html).

Observation-specific measurement error is a central feature of the
original 2S-PA framework. `tspa_mx_model()` uses OpenMx definition
variables to retain these individual-specific values. `tspa()` uses
lavaan and pools varying measurement quantities within each group. The
OpenMx and IRT examples require the optional packages `OpenMx` and
`mirt`, respectively.

## Methodological background

- Lai, M. H. C., & Hsiao, Y.-Y. (2022). Two-stage path analysis with
  definition variables: An alternative framework to account for
  measurement error. *Psychological Methods, 27*(4), 568–588.
  <https://doi.org/10.1037/met0000410>
- Lai, M. H. C., Tse, W. W.-Y., Zhang, G., Li, Y., & Hsiao, Y.-Y.
  (2023). Correcting for unreliability and partial invariance: A
  two-stage path analysis approach. *Structural Equation Modeling: A
  Multidisciplinary Journal, 30*(2), 258–271.
  <https://doi.org/10.1080/10705511.2022.2125397>
- Vermunt, J. K. (2025). Stepwise estimation of latent variable models:
  An overview of approaches. *Statistical Modelling, 25*(6), 530–551.
  <https://doi.org/10.1177/1471082X251355693>

## Acknowledgment

This package is based upon work supported by the National Science
Foundation under Grant No. 2141790.

<!-- `devtools::build_readme()` -->
