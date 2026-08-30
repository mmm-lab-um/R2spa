# Linear Growth Modeling with Two-Stage Path Analysis

``` r

library(lavaan)
library(R2spa)
```

We use the example from [this
tutorial](https://thechangelab.stanford.edu/tutorials/growth-modeling/measurement-invariance-2nd-order-growth-model/)
from Chapter 14 of the book *Growth Modeling* (Grimm, Ram & Estabrook,
2017) to demonstrate how to perform linear growth modeling with
two-stage path analysis (2S-PA)

``` r

# Load data
eclsk <- read.table(
  "https://raw.githubusercontent.com/LRI-2/Data/main/GrowthModeling/ECLS_Science.dat",
  na.strings = "."
)
names(eclsk) <- c("id", "s_g3", "r_g3", "m_g3", "s_g5", "r_g5", "m_g5", "s_g8",
                  "r_g8", "m_g8", "st_g3", "rt_g3", "mt_g3", "st_g5", "rt_g5",
                  "mt_g5", "st_g8", "rt_g8", "mt_g8")
```

## Joint structural equation model: Latent growth with strict invariance model

In the
[tutorial](https://quantdev.ssri.psu.edu/tutorials/growth-modeling-chapter-14-modeling-change-latent-variables-measured-continuous),
the authors first performed longitudinal invariance testing and found
support for strict invariance. They moved on to fit a latent growth
model based on the strict invariance model, as shown below.

``` r

jsem_growth_mod <- "
# factor loadings (with constraints)
eta1 =~ 15.1749088 * s_g3 + l2 * r_g3 +  l3 * m_g3
eta2 =~ 15.1749088 * s_g5 + l2 * r_g5 +  l3 * m_g5
eta3 =~ 15.1749088 * s_g8 + l2 * r_g8 +  l3 * m_g8

# factor variances
eta1 ~~ psi * eta1
eta2 ~~ psi * eta2
eta3 ~~ psi * eta3

# unique variances/covariances 
s_g3 ~~ u1 * s_g3 + s_g5 + s_g8
s_g5 ~~ u1 * s_g5 + s_g8
s_g8 ~~ u1 * s_g8
r_g3 ~~ u2 * r_g3 + r_g5 + r_g8
r_g5 ~~ u2 * r_g5 + r_g8
r_g8 ~~ u2 * r_g8
m_g3 ~~ u3 * m_g3 + m_g5 + m_g8
m_g5 ~~ u3 * m_g5 + m_g8
m_g8 ~~ u3 * m_g8

# observed variable intercepts
s_g3 ~ i1 * 1
s_g5 ~ i1 * 1
s_g8 ~ i1 * 1
r_g3 ~ i2 * 1
r_g5 ~ i2 * 1
r_g8 ~ i2 * 1
m_g3 ~ i3 * 1
m_g5 ~ i3 * 1
m_g8 ~ i3 * 1

# latent basis model
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + start(.5) * eta2 + 1 * eta3

i ~~ start(.8) * i
s ~~ start(.5) * s
i ~~ start(0) * s

i ~ 0 * 1
s ~ 1
"
jsem_growth_fit <- lavaan(jsem_growth_mod,
                          data = eclsk,
                          meanstructure = TRUE,
                          estimator = "ML",
                          fixed.x = FALSE)
```

## 2S-PA: Latent growth with strict invariance model

With 2S-PA, we first specify the strict invariance model, directly
adopted from the
[tutorial](https://quantdev.ssri.psu.edu/tutorials/growth-modeling-chapter-14-modeling-change-latent-variables-measured-continuous).
Based on the strict invariance model, we then obtain the factor scores.

``` r

strict_mod <- "
# factor loadings
eta1 =~ 15.1749088 * s_g3 + l2 * r_g3 + l3 * m_g3
eta2 =~ 15.1749088 * s_g5 + l2 * r_g5 + l3 * m_g5
eta3 =~ 15.1749088 * s_g8 + l2 * r_g8 + l3 * m_g8

# factor variances/covariances
eta1 ~~ 1 * eta1 + eta2 + eta3
eta2 ~~ eta2 + eta3
eta3 ~~ eta3

# unique variances/covariances
s_g3 ~~ u1 * s_g3 + s_g5 + s_g8
s_g5 ~~ u1 * s_g5 + s_g8
s_g8 ~~ u1 * s_g8
r_g3 ~~ u2 * r_g3 + r_g5 + r_g8
r_g5 ~~ u2 * r_g5 + r_g8
r_g8 ~~ u2 * r_g8
m_g3 ~~ u3 * m_g3 + m_g5 + m_g8
m_g5 ~~ u3 * m_g5 + m_g8
m_g8 ~~ u3 * m_g8

# latent variable intercepts
eta1 ~ 0 * 1
eta2 ~ 1
eta3 ~ 1

# observed variable intercepts
s_g3 ~ i1 * 1
s_g5 ~ i1 * 1
s_g8 ~ i1 * 1
r_g3 ~ i2 * 1
r_g5 ~ i2 * 1
r_g8 ~ i2 * 1
m_g3 ~ i3 * 1
m_g5 ~ i3 * 1
m_g8 ~ i3 * 1
"
```

``` r

# Get factor scores
fs_dat <- get_fs(eclsk, model = strict_mod)
```

In the second stage, we use the factor scores obtained from the strict
invariance model to model the growth trajectory across time points. The
growth model is the same as the “latent basis model” in the joint
structural equation model.

``` r

# Growth model
tspa_growth_mod <- "
i =~ 1 * eta1 + 1 * eta2 + 1 * eta3
s =~ 0 * eta1 + start(.5) * eta2 + 1 * eta3

# factor variances
eta1 ~~ psi * eta1
eta2 ~~ psi * eta2
eta3 ~~ psi * eta3

i ~~ start(.8) * i
s ~~ start(.5) * s
i ~~ start(0) * s

i ~ 1
s ~ 1
"
```

``` r

# Fit the growth model
tspa_growth_fit <- tspa(tspa_growth_mod, fs_dat,
                        fsT = attr(fs_dat, "fsT"),
                        fsL = attr(fs_dat, "fsL"),
                        fsb = attr(fs_dat, "fsb"),
                        estimator = "ML")
```

For objects extracted from
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md),
one can omit the `fsT`, `fsL`, and `fsb` arguments, as the
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md)
function will automatically extract the necessary information from the
`fs_dat` object.

``` r

tspa_growth_fit <- tspa(tspa_growth_mod, fs_dat, estimator = "ML")
```

## Comparison

The joint structural equation model and 2S-PA yielded comparable
estimates of the mean and variances of the intercept and slope factors.

``` r

parameterestimates(tspa_growth_fit) |>
  subset(subset = lhs %in% c("i", "s") & op %in% c("~1", "~~"))
#>    lhs op rhs label    est    se       z pvalue ci.lower ci.upper
#> 28   i ~~   i        0.944 0.051  18.592  0.000    0.845    1.044
#> 29   s ~~   s        0.118 0.016   7.298  0.000    0.086    0.149
#> 30   i ~~   s       -0.076 0.020  -3.757  0.000   -0.116   -0.036
#> 31   i ~1            0.001 0.035   0.016  0.987   -0.068    0.069
#> 32   s ~1            1.993 0.019 107.413  0.000    1.957    2.030
parameterestimates(jsem_growth_fit) |>
  subset(subset = lhs %in% c("i", "s") & op  %in% c("~1", "~~"))
#>    lhs op rhs label    est    se      z pvalue ci.lower ci.upper
#> 46   i ~~   i        0.941 0.052 18.138      0    0.839    1.042
#> 47   s ~~   s        0.117 0.017  6.826      0    0.083    0.150
#> 48   i ~~   s       -0.074 0.020 -3.646      0   -0.115   -0.034
#> 49   i ~1            0.000 0.000     NA     NA    0.000    0.000
#> 50   s ~1            1.991 0.024 84.586      0    1.945    2.038
```
