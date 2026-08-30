# First-order correction of sampling covariance for 2S-PA estimates

First-order correction of sampling covariance for 2S-PA estimates

## Usage

``` r
vcov_corrected(
  tspa_fit,
  vfsLT,
  which_free = NULL,
  engine = c("analytic", "fd"),
  ...
)
```

## Arguments

- tspa_fit:

  A fit from
  [`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) with
  `fsT` and `fsL` supplied (multi-factor measurement model, single- or
  multi-group), so that it carries the `fsT`, `fsL`, and `tspa_args`
  attributes. A fit corrected in place via `tspa(corrected_se = TRUE)`
  (attribute `tspa_corrected = TRUE`) is rejected.

- vfsLT:

  The sampling covariance matrix of the free `fsL`/`fsT` elements, taken
  (or sub-matrixed) from the `vfsLT` attribute of a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result fitted with `vfsLT = TRUE`. Its row and column order is the
  same `fsL`-then-`fsT` order described by `which_free`. When
  `which_free = NULL` (the default) this is the full square matrix
  returned by
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md);
  when `which_free` is a subset of length `k`, this **must** be the
  matching `k x k` principal submatrix, i.e.
  `vfsLT_full[which_free, which_free]`.

- which_free:

  An optional numeric vector of positions selecting which `fsL`/`fsT`
  elements to treat as free (and therefore to propagate through the
  Jacobian). Positions run over the `fsL` matrix (column-major order)
  followed by the lower-triangular part of `fsT` including the diagonal
  (column-major order). For a two-factor model `fsL` and `fsT` are both
  2 x 2: `fsL` occupies positions 1:4 and `fsT` positions 5:7, so the
  two error variances are 5 and 7 and the error covariance between the
  two factor scores (the `[2, 1]` element of `fsT`) is position 6. For a
  multigroup fit the positions run per group, in group order: group 1's
  full `fsL` (column-major), then group 2's full `fsL`, ..., then group
  1's lower-triangular `fsT` (column-major), then group 2's, i.e. all
  loadings across groups first, then all error-variance elements
  (matching the order of the `vfsLT` attribute from
  `get_fs(vfsLT = TRUE)`). `NULL` (the default) treats every per-group
  `fsL`/`fsT` element as free. A non-`NULL` `which_free` of length `k`
  requires `vfsLT` to be the matching `k x k` principal submatrix (see
  `vfsLT`).

- engine:

  The engine used to evaluate the Jacobian `J = d(thetahat)/d(eta)`.
  `"analytic"` (the default) uses a refit-free, deterministic
  influence-function closed form (PLAN 16, sections 2.4 and 4.3):
  `J = -H^{-1} C`, with `H` (the log-likelihood Hessian over the free
  params) and `C` (the cross-derivative w.r.t. the fixed `fsL`/`fsT`
  entries) obtained by central-differencing the analytic log-likelihood
  score. It covers single- and multi-group models, saturated and
  restricted (df \> 0) structural models, and mean-structure models, and
  is a pure function of the base fit + `vfsLT` (bit-reproducible, no
  refits). `"fd"` uses central finite differences (one stage-2 refit on
  each side of each free element) and is retained as the A/B reference;
  the analytic path transparently falls back to it only for a shape it
  cannot handle (an unrecognised free parameter or unequal per-group
  free-param counts). The two agree to the finite-difference noise floor
  whenever `"analytic"` applies.

- ...:

  Currently not used.

## Value

A corrected covariance matrix in the same dimension as `vcov(tspa_fit)`
(symmetric).

## Details

`vcov_corrected()` applies the first-order (delta-method) two-stage
approximation: the stage-2 covariance `vcov(tspa_fit)` is augmented by
`J %*% vfsLT %*% t(J)`, where `J` is the Jacobian of the stage-2
estimates with respect to the selected `fsL`/`fsT` free elements. With
the default `engine = "analytic"`, `J` is evaluated refit-free via a
closed form (see the `engine` argument below); with `engine = "fd"` it
is instead estimated by central differences, one full stage-2 refit on
each side of each free element while reusing the base fit's
coefficients, so the cost is `2 x (number of free elements)` stage-2
refits.

The correction is **partial by design**. It propagates only the
*sampling* uncertainty of the stage-1 estimates of `fsL`
(loadings/cross-loadings) and `fsT` (score error variance-covariance,
including its off-diagonal). It does **not** account for the sampling
uncertainty of the factor-score *values* (the `fs_<name>` columns fitted
in stage 2), their standard errors (`se_fs`), or the score intercepts
(`fsb`); these are held fixed. The two-stage approximation also treats
the stages as independent, so the cross-covariance between the stage-1
estimates and the stage-2 data is not modelled.

The same correction can also be requested in place via
`tspa(..., corrected_se = TRUE, vfsLT = <matrix>)`;
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) then
sets the `tspa_corrected = TRUE` attribute on the returned fit, and
[`standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
on that fit reports corrected standard SEs. Such an already-corrected
fit is rejected by `vcov_corrected()` (the correction is never applied
twice).

The `engine` argument selects how the Jacobian is evaluated. The default
`"analytic"` is refit-free and deterministic: an influence-function
closed form (PLAN 16, sections 2.4 and 4.3), `J = -H^{-1} C`, with `H`
(the log-likelihood Hessian over the free params) and `C` (the
cross-derivative w.r.t. the fixed `fsL`/`fsT` entries) obtained by
central-differencing the analytic log-likelihood score. It covers
single- and multi-group models, saturated and restricted (df \> 0)
structural models, and mean-structure models, and is a pure function of
the base fit + `vfsLT` (bit-reproducible, no refits). `"fd"` (central
finite differences, one stage-2 refit per side of each free element) is
retained as the A/B reference; the analytic path falls back to it only
for a shape it cannot handle. The two agree to the finite-difference
noise floor whenever `"analytic"` applies.

## See also

- `vignette("Corrected Standard Errors", package = "R2spa")` for the
  corrected-SE workflow.

- `vignette("Correction to Measurement Error", package = "R2spa")` for
  the underlying error correction.

## Examples

``` r
library(lavaan)

# Two-factor model, Bartlett scoring. The fsL/fsT free elements run over
# positions 1:7 (fsL -> 1:4, fsT -> 5:7); the two error variances are 5
# and 7, the score-error covariance is 6.
fs <- get_fs(PoliticalDemocracy,
             model = "ind60 =~ x1 + x2 + x3
                      dem60 =~ y1 + y2 + y3 + y4",
             method = "Bartlett", vfsLT = TRUE)
fit <- tspa("dem60 ~ ind60", data = fs,
            fsT = attr(fs, "fsT"), fsL = attr(fs, "fsL"))
vfsLT <- attr(fs, "vfsLT")

# Propagate only the two error variances (positions 5, 7): with a
# non-NULL which_free, pass the matching principal submatrix of vfsLT.
vc <- vcov_corrected(fit, vfsLT = vfsLT[c(5, 7), c(5, 7)],
                     which_free = c(5, 7))
sqrt(diag(vc))   # corrected standard errors
#>  dem60~ind60 ind60~~ind60 dem60~~dem60 
#>   0.35271754   0.07634634   0.66959805 
```
