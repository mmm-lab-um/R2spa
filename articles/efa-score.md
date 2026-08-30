# EFA Scores

This vignette demonstrates how to use the
[`compute_fscore()`](https://mmm-lab-um.github.io/R2spa/reference/compute_fscore.md)
function to calculate factor scores based on exploratory factor
analysis, and compare the results to those calculated by the
[`psych::factor.scores()`](https://rdrr.io/pkg/psych/man/factor.scores.html)
function.

``` r

library(psych)
data(bfi, package = "psych")
library(lavaan)
#> This is lavaan 0.7-2
#> lavaan is FREE software! Please report any bugs.
#> 
#> Attaching package: 'lavaan'
#> The following object is masked from 'package:psych':
#> 
#>     cor2cov
```

## Example: Big Five Inventory

``` r

# Complete-case (listwise) correlation — the same rows used to compute the
# scores below, so stage 1 (EFA) and stage 2 (2S-PA) share one data basis.
cc <- complete.cases(bfi[, 1:25])
corr_bfi <- cor(bfi[cc, 1:25])
# EFA (Target rotation)
target_mat_bfi <- matrix(0, nrow = 25, ncol = 5)
target_mat_bfi[1:5, 1] <- NA
target_mat_bfi[6:10, 2] <- NA
target_mat_bfi[11:15, 3] <- NA
target_mat_bfi[16:20, 4] <- NA
target_mat_bfi[21:25, 5] <- NA
fa_target_bfi <- psych::fa(
    corr_bfi, n.obs = sum(cc),
    nfactors = 5,
    rotate = "targetQ", Target = target_mat_bfi,
    scores = "Bartlett",
    n.rotations = 1)  # <-- avoids the psych 2.6.5 faRotations() bug
#> Loading required namespace: GPArotation
# Factor correlations
fa_target_bfi$Phi |>
    (`[`)(paste0("MR", 1:5), paste0("MR", 1:5)) |>
    knitr::kable(digits = 2, caption = "EFA Factor Correlation")
```

|     |   MR1 |   MR2 |   MR3 |   MR4 |  MR5 |
|:----|------:|------:|------:|------:|-----:|
| MR1 |  1.00 |  0.22 |  0.35 | -0.08 | 0.14 |
| MR2 |  0.22 |  1.00 |  0.28 | -0.20 | 0.19 |
| MR3 |  0.35 |  0.28 |  1.00 | -0.21 | 0.16 |
| MR4 | -0.08 | -0.20 | -0.21 |  1.00 | 0.01 |
| MR5 |  0.14 |  0.19 |  0.16 |  0.01 | 1.00 |

EFA Factor Correlation {.table}

``` r

# Correlation with sum scores
bfi |>
    transform(A = (7 - A1) + A2 + A3 + A4 + A5,
           C = C1 + C2 + C3 + (7 - C4) + (7 - C5),
           E = (7 - E1) + (7 - E2) + E3 + E4 + E5,
           N = N1 + N2 + N3 + N4 + N5,
           O = O1 + (7 - O2) + O3 + O4 + (7 - O5)) |>
    subset(select = A:O) |>
    cor(use = "complete") |>
    knitr::kable(digits = 2, caption = "Sum scores")
```

|     |     A |     C |     E |     N |     O |
|:----|------:|------:|------:|------:|------:|
| A   |  1.00 |  0.26 |  0.47 | -0.19 |  0.14 |
| C   |  0.26 |  1.00 |  0.27 | -0.23 |  0.19 |
| E   |  0.47 |  0.27 |  1.00 | -0.23 |  0.22 |
| N   | -0.19 | -0.23 | -0.23 |  1.00 | -0.08 |
| O   |  0.14 |  0.19 |  0.22 | -0.08 |  1.00 |

Sum scores {.table}

Hand calculate the Bartlett scores using weights

``` r

# Bartlett score for first person
bscores <-
    psych::factor.scores(bfi[, 1:25], f = fa_target_bfi,
                         method = "Bartlett")
fa_target_bfi$weights
#>             MR4          MR3          MR2          MR1          MR5
#> A1  0.014402183  0.092380968  0.035746321 -0.193565413 -0.045036549
#> A2  0.041124721 -0.058105557  0.015488614  0.382773898  0.012937569
#> A3  0.041432241 -0.006587634 -0.014415176  0.443581386  0.000701728
#> A4  0.008071086 -0.004947141  0.071824635  0.193766446 -0.110175036
#> A5 -0.009604888  0.065022160 -0.028184518  0.282085524  0.003273756
#> C1  0.027090659 -0.017220587  0.250848487 -0.015579885  0.092474527
#> C2  0.064163619 -0.054953490  0.366653319  0.038376054  0.025965380
#> C3  0.019205794 -0.038889507  0.255423628  0.037393639 -0.053014339
#> C4  0.040554906  0.022465420 -0.363740725  0.038737702 -0.013671015
#> C5  0.053243933 -0.048611615 -0.294985973  0.042119161  0.106815302
#> E1  0.003629865 -0.232100888  0.063304564  0.030924582  0.002797333
#> E2  0.068569876 -0.388082053  0.012996445  0.085579476  0.053356965
#> E3  0.024907048  0.190395745 -0.027277362  0.059977519  0.168130445
#> E4 -0.012652895  0.317206605 -0.013257326  0.096760720 -0.162301617
#> E5  0.032448621  0.194592910  0.117003697 -0.045875151  0.089532244
#> N1  0.395908950  0.187207604  0.032047073 -0.151290626 -0.099970010
#> N2  0.308558388  0.096554724  0.033254576 -0.099289054  0.008145985
#> N3  0.260704614 -0.021915505 -0.005781887  0.066188637  0.025631532
#> N4  0.183264898 -0.196804612 -0.068885085  0.128908778  0.135782260
#> N5  0.137867460 -0.080178115  0.010361065  0.125979385 -0.081435235
#> O1  0.005371816  0.042302579  0.011644648 -0.022922677  0.308623067
#> O2  0.040541202  0.017885382 -0.021159771  0.073057487 -0.280501573
#> O3  0.017091505  0.087387231 -0.025874628 -0.003339764  0.474002153
#> O4  0.048235950 -0.120850280 -0.018694335  0.099355799  0.237795617
#> O5  0.020079980  0.034390843  0.001426520  0.025223473 -0.336536454
# Calculation by hand
y1 <- scale(bfi[, 1:25])[1, ]  # z-score
crossprod(fa_target_bfi$weights, as.matrix(y1))
#>            [,1]
#> MR4 -0.37637333
#> MR3  0.03664292
#> MR2 -1.64227837
#> MR1 -1.07426561
#> MR5 -2.25597237
# Compare to results from psych::fa()
bscores$scores[1, ]
#>         MR4         MR3         MR2         MR1         MR5 
#> -0.37637333  0.03664292 -1.64227837 -1.07426561 -2.25597237
```

``` r

# Covariance of Bartlett scores
cov(bscores$scores, use = "complete") |>
    (`[`)(paste0("MR", 1:5), paste0("MR", 1:5)) |>
    knitr::kable(digits = 2, caption = "With Bartlett scores")
```

|     |   MR1 |   MR2 |   MR3 |   MR4 |  MR5 |
|:----|------:|------:|------:|------:|-----:|
| MR1 |  1.34 |  0.20 |  0.29 | -0.06 | 0.13 |
| MR2 |  0.20 |  1.30 |  0.27 | -0.20 | 0.18 |
| MR3 |  0.29 |  0.27 |  1.27 | -0.21 | 0.13 |
| MR4 | -0.06 | -0.20 | -0.21 |  1.17 | 0.01 |
| MR5 |  0.13 |  0.18 |  0.13 |  0.01 | 1.42 |

With Bartlett scores {.table}

### Using `compute_fscore()` and perform a two-stage analysis

Use
[`R2spa::compute_fscore()`](https://mmm-lab-um.github.io/R2spa/reference/compute_fscore.md)

``` r

# Obtain error covariances
yc <- scale(bfi[, 1:25])
yc <- yc[complete.cases(yc), ]
lam <- fa_target_bfi$loadings
colnames(lam) <- c("N", "E", "C", "A", "O")
phi <- fa_target_bfi$Phi
th <- diag(fa_target_bfi$uniquenesses)
# # scoring weights
# a <- solve(crossprod(lam, solve(th, lam)), t(solve(th, lam)))
# ecov_fs <- a %*% th %*% t(a)
# dimnames(ecov_fs) <- rep(list(c("N", "E", "C", "A", "O")), 2)
# Two-stage analysis
library(R2spa)
bfi_fs <- compute_fscore(yc, lambda = lam, theta = th,
                         method = "Bartlett", center_y = FALSE,
                         fs_matrices = TRUE)
head(bfi_fs)
#>                 N           E            C           A          O
#> 61617 -0.37637333  0.03664292 -1.642278370 -1.07426561 -2.2559724
#> 61618  0.06602659  0.60969483 -0.764512534 -0.19811184 -0.3023409
#> 61620  0.66383177  0.37175492 -0.004350127 -0.95821234  0.2931325
#> 61621 -0.14486826 -0.04522127 -1.366864455  0.02201861 -1.4712146
#> 61622 -0.37429555  0.53847538 -0.074453722 -1.02844224 -1.0240800
#> 61623  0.17203114  1.46173826  1.798527481  0.10919228  0.5072230
# Scoring matrix
attr(bfi_fs, which = "scoring_matrix")
#>          [,1]        [,2]         [,3]         [,4]         [,5]        [,6]
#> N  0.01440218  0.04112472  0.041432241  0.008071086 -0.009604888  0.02709066
#> E  0.09238097 -0.05810556 -0.006587634 -0.004947141  0.065022160 -0.01722059
#> C  0.03574632  0.01548861 -0.014415176  0.071824635 -0.028184518  0.25084849
#> A -0.19356541  0.38277390  0.443581386  0.193766446  0.282085524 -0.01557988
#> O -0.04503655  0.01293757  0.000701728 -0.110175036  0.003273756  0.09247453
#>          [,7]        [,8]        [,9]       [,10]        [,11]       [,12]
#> N  0.06416362  0.01920579  0.04055491  0.05324393  0.003629865  0.06856988
#> E -0.05495349 -0.03888951  0.02246542 -0.04861161 -0.232100888 -0.38808205
#> C  0.36665332  0.25542363 -0.36374073 -0.29498597  0.063304564  0.01299645
#> A  0.03837605  0.03739364  0.03873770  0.04211916  0.030924582  0.08557948
#> O  0.02596538 -0.05301434 -0.01367102  0.10681530  0.002797333  0.05335696
#>         [,13]       [,14]       [,15]       [,16]        [,17]        [,18]
#> N  0.02490705 -0.01265289  0.03244862  0.39590895  0.308558388  0.260704614
#> E  0.19039575  0.31720661  0.19459291  0.18720760  0.096554724 -0.021915505
#> C -0.02727736 -0.01325733  0.11700370  0.03204707  0.033254576 -0.005781887
#> A  0.05997752  0.09676072 -0.04587515 -0.15129063 -0.099289054  0.066188637
#> O  0.16813044 -0.16230162  0.08953224 -0.09997001  0.008145985  0.025631532
#>         [,19]       [,20]        [,21]       [,22]        [,23]       [,24]
#> N  0.18326490  0.13786746  0.005371816  0.04054120  0.017091505  0.04823595
#> E -0.19680461 -0.08017812  0.042302579  0.01788538  0.087387231 -0.12085028
#> C -0.06888509  0.01036107  0.011644648 -0.02115977 -0.025874628 -0.01869433
#> A  0.12890878  0.12597939 -0.022922677  0.07305749 -0.003339764  0.09935580
#> O  0.13578226 -0.08143523  0.308623067 -0.28050157  0.474002153  0.23779562
#>         [,25]
#> N  0.02007998
#> E  0.03439084
#> C  0.00142652
#> A  0.02522347
#> O -0.33653645
# Error covariance
attr(bfi_fs, which = "fsT")
#>              fs_N         fs_E         fs_C         fs_A         fs_O
#> fs_N  0.161187532 -0.007284208  0.008122536  0.025896466  0.007956845
#> fs_E -0.007284208  0.263478478 -0.008953332 -0.066595623 -0.029653071
#> fs_C  0.008122536 -0.008953332  0.297667780 -0.014661536 -0.015195505
#> fs_A  0.025896466 -0.066595623 -0.014661536  0.326546394 -0.007596324
#> fs_O  0.007956845 -0.029653071 -0.015195505 -0.007596324  0.438466397
```

Recover factor covariances with 2S-PA

``` r

ts_fit <- tspa("",
               data = data.frame(bscores$scores) |>
                   setNames(c("fs_N", "fs_E", "fs_C", "fs_A", "fs_O")),
               fsT = attr(bfi_fs, which = "fsT"),
               fsL = diag(5) |>
                   `dimnames<-`(list(c("fs_N", "fs_E", "fs_C", "fs_A", "fs_O"),
                                     c("N", "E", "C", "A", "O"))))
lavInspect(ts_fit, what = "cor.lv") |>
    (`[`)(c("A", "C", "E", "N", "O"), c("A", "C", "E", "N", "O")) |>
    knitr::kable(digits = 2, caption = "With Bartlett scores and 2S-PA")
```

|     |     A |     C |     E |     N |    O |
|:----|------:|------:|------:|------:|-----:|
| A   |  1.00 |  0.22 |  0.35 | -0.09 | 0.14 |
| C   |  0.22 |  1.00 |  0.28 | -0.20 | 0.19 |
| E   |  0.35 |  0.28 |  1.00 | -0.21 | 0.16 |
| N   | -0.09 | -0.20 | -0.21 |  1.00 | 0.00 |
| O   |  0.14 |  0.19 |  0.16 |  0.00 | 1.00 |

With Bartlett scores and 2S-PA {.table}
