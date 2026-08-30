# Two-Stage Path Analysis

Fit a two-stage path analysis (2S-PA) model.

## Usage

``` r
tspa(
  model,
  data,
  reliability = NULL,
  se = "standard",
  se_fs = NULL,
  fsT = NULL,
  fsL = NULL,
  fsb = NULL,
  reduce = c("mean", "median"),
  vfsLT = NULL,
  corrected_se = FALSE,
  which_free = NULL,
  product = FALSE,
  engine = "analytic",
  ...
)
```

## Arguments

- model:

  A string variable describing the structural path model, in `lavaan`
  syntax.

- data:

  A data frame containing factor scores. When `data` is a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result and the corresponding arguments are omitted, the measurement
  inputs are derived from it: the `fsT`/`fsL`/`fsb` attributes for a
  multi-factor fit, or, for a
  [`cbind()`](https://rdrr.io/r/base/cbind.html)'d frame whose
  attributes [`cbind()`](https://rdrr.io/r/base/cbind.html) drops, the
  `fs_<v>`/`fs_<v>_se` score columns for a single-factor fit (see
  `Details`).

- reliability:

  A numeric vector representing the reliability indexes of each latent
  factor. Currently `tspa()` does not support the reliability argument.
  Please use `se`.

- se:

  Deprecated to avoid conflict with the argument of the same name in
  [`lavaan::lavaan()`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

- se_fs:

  A numeric vector representing the standard errors of each factor score
  variable for single-group 2S-PA. A list or data frame storing the
  standard errors of each group in each latent factor for multigroup
  2S-PA. An explicit `se_fs` always wins over derivation from a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result and, even when empty
  ([`list()`](https://rdrr.io/r/base/list.html)), suppresses the
  multi-factor (attribute) derivation; the single-factor `se_fs` is
  derived from the data's `fs_<v>_se` columns only when this argument is
  omitted.

- fsT:

  An error variance-covariance matrix of the factor scores, which can be
  obtained from the output of
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  using [`attr()`](https://rdrr.io/r/base/attr.html) with the argument
  `which = "fsT"`. When a group was fitted with missing data
  (`missing = "fiml"`), the attribute carries per-pattern values (a
  per-group list of one matrix per observed-indicator pattern); for a
  `merMod` fit it is a 3-D per-cluster array; for a `mirt` fit it is a
  per-observation list (one matrix per row, marked `mirt_per_obs`).
  Values of these per-unit shapes are reduced to a single representative
  per-group matrix by `reduce`; the pooled (not the nested/per-cluster)
  matrix is attached to the returned fit as the `fsT` attribute. When
  omitted, `fsT` is derived from the `fsT` attribute of a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result passed as `data` (see `Details`); an explicit `fsT` always
  wins.

- fsL:

  A matrix of loadings and cross-loadings from the latent variables to
  the factor scores `fs`, which can be obtained from the output of
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  using [`attr()`](https://rdrr.io/r/base/attr.html) with the argument
  `which = "fsL"`. For details see the Multi-Factor Measurement Model
  vignette:
  `vignette("Multi-Factor Measurement Model", package = "R2spa")`. As
  with `fsT`, per-pattern (FIML missing data), per-cluster (merMod), and
  per-observation (mirt) values are supported and reduced per group by
  `reduce`; the pooled (not the nested) matrix is attached to the
  returned fit as the `fsL` attribute. When omitted, `fsL` is derived
  from the `fsL` attribute of a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result passed as `data` (see `Details`); an explicit `fsL` always
  wins.

- fsb:

  A vector of intercepts for the factor scores `fs`, which can be
  obtained from the output of
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  using [`attr()`](https://rdrr.io/r/base/attr.html) with the argument
  `which = "fsb"`. As with `fsT`, per-pattern (FIML missing data),
  per-cluster (merMod), and per-observation (mirt) values are supported
  and reduced per group by `reduce`; the pooled (not the
  nested/per-cluster) vector is used for the stage-2 intercept
  constraints. When omitted, `fsb` is derived from the `fsb` attribute
  of a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result passed as `data` (which may be absent — `merMod` results carry
  none); an explicit `fsb` always wins.

- reduce:

  Controls how per-unit `fsL`/`fsT`/`fsb` from a group fitted with
  missing data (`missing = "fiml"`; per-pattern lists of matrices),
  per-cluster values from a `merMod` fit (3-D arrays), or
  per-observation values from a `mirt` fit, are collapsed to a single
  representative value per group for stage 2. A no-op when the per-unit
  quantities are constant within the group (e.g. complete single-group
  data). A group (or the whole data) with no scorable rows at all (every
  row missing) is an error: there is nothing to pool. One of `"mean"`
  (the default) or `"median"`. With `"mean"` the pooled `fsT` is a
  convex combination of the per-unit (positive semi-definite) matrices
  and so remains positive semi-definite; with `"median"` the reduction
  is element-wise and need not be, in which case a warning is emitted
  when the pooled `fsT` is not positive semi-definite. The pooled (not
  the nested/per-cluster) `fsT` and `fsL` are what get attached to the
  returned fit.

- vfsLT:

  The sampling covariance matrix of the free `fsL`/`fsT` elements, taken
  from the `vfsLT` attribute of a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result fitted with `vfsLT = TRUE`. Required when
  `corrected_se = TRUE`; ignored otherwise.

- corrected_se:

  A logical; when `TRUE`, the stage-2 covariance of the returned fit is
  replaced by the first-order (delta-method) correction of
  [`vcov_corrected()`](https://mmm-lab-um.github.io/R2spa/reference/vcov_corrected.md)
  and the `tspa_corrected` attribute is set to `TRUE`. Requires a
  multi-factor fit (both `fsT` and `fsL` supplied) and `vfsLT`.
  Supported for single-group and multigroup fits (the multigroup
  `fsL`/`fsT` are the per-group list attributes from
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)).
  Default `FALSE` (the returned fit is unchanged).

- which_free:

  An optional numeric vector of positions selecting which `fsL`/`fsT`
  free elements to propagate through the corrected covariance (see
  [`vcov_corrected()`](https://mmm-lab-um.github.io/R2spa/reference/vcov_corrected.md));
  used only when `corrected_se = TRUE`.

- product:

  A logical; when `TRUE`, `tspa()` automatically computes the
  double-mean-centered product indicators (via
  [`compute_fs_prod()`](https://mmm-lab-um.github.io/R2spa/reference/compute_fs_prod.md))
  for every latent in `model` that names the product of two of the
  model's factor scores — by concatenation (the latent `xm` is the
  product of the scores of `x` and `m`) or in lavaan's interaction
  syntax (`x:m`, rendered under the concatenated name) — and
  incorporates them into the stage-2 measurement model, so
  `tspa("y ~ x + m + x:m", data = get_fs(...), product = TRUE)` needs no
  pre-computed product columns. In the single-factor (score-scale) path
  the product latent loads 1 on its indicator with error variance the
  (per-group pooled, `reduce`) product SE, like every other
  single-factor latent; in the multi-factor path it loads with the
  implied loading `gamma` and error variance `se_P^2`, both evaluated at
  the (pooled) `fsL`/`fsT` with the `psi` attribute. When two product
  latents share a factor score (e.g. `x:m` and `x:z` share `x`), their
  indicators' measurement errors are correlated, and `tspa()` fixes
  those error covariances in the stage-2 model, estimated from the
  stage-1 `fsL`/`fsT`/`psi` (the Isserlis expansion of the joint-normal
  score-error moments — the same joint-normality assumptions the
  [`compute_fs_prod()`](https://mmm-lab-um.github.io/R2spa/reference/compute_fs_prod.md)
  product SE rests on); with a single product latent there is nothing to
  fix. Single-group models only (v1); not supported with
  `corrected_se = TRUE`. Default `FALSE`, which leaves the manual
  workflow
  ([`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  with `product` set, or
  [`compute_fs_prod()`](https://mmm-lab-um.github.io/R2spa/reference/compute_fs_prod.md),
  up front, the product SE in `se_fs`) unchanged.

- engine:

  The Jacobian engine passed to
  [`vcov_corrected()`](https://mmm-lab-um.github.io/R2spa/reference/vcov_corrected.md)
  when `corrected_se = TRUE`: `"analytic"` (the default, a refit-free,
  deterministic influence-function closed form covering single- and
  multi-group, saturated and restricted, and mean-structure models) or
  `"fd"` (finite differences, retained as the A/B reference). Ignored
  unless `corrected_se = TRUE`.

- ...:

  Additional arguments passed to
  [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html). See
  [`lavOptions`](https://rdrr.io/pkg/lavaan/man/lavOptions.html) for a
  complete list.

## Value

An object of class `lavaan` carrying the following `R2spa`-specific
attributes: `tspaModel`, the stage-2 model syntax actually fitted;
`tspa_call`, the matched `tspa()` call; and `tspa_args`, the argument
list used to build the stage-2 model (captured at fit time as evaluated
values), which lets the fit be re-evaluated without its original
environment, e.g. by the
[`vcov_corrected()`](https://mmm-lab-um.github.io/R2spa/reference/vcov_corrected.md)
SE correction. On a fit whose measurement inputs were derived from
`data` (see `Details`), `tspa_args` carries the *resolved* values
(post-pooling plain matrices, the derived `se_fs`), so a replay via
`do.call(tspa, attr(fit, "tspa_args"))` re-passes them explicitly, skips
the derivation, and cannot double-pool. When `fsT`/`fsL` are supplied
(multi-factor measurement model), the (possibly reduced) matrices are
also attached as the `fsT`/`fsL` attributes, and `pooled_fs` records the
`reduce` method used when per-unit values were collapsed. When
`corrected_se = TRUE`, the returned fit additionally carries
`tspa_corrected = TRUE` and its covariance is the first-order corrected
matrix, so [`vcov()`](https://rdrr.io/r/stats/vcov.html), `se()`, and
[`standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
on it report the corrected standard errors, for multigroup fits as well
as single-group ones.

## Details

When the factor-score attributes are heterogeneous across the units of a
group, `tspa()` first re-expresses them as long-form,
individual-specific values and then reduces them to a single
representative set per group; that pooled set is what feeds the stage-2
model and is attached to the returned fit. The heterogeneous cases are
per-pattern values from a group fitted with missing data
(`missing = "fiml"`), where `fsL`/`fsT`/`fsb` are per-group lists of one
matrix/vector per observed-indicator pattern, and per-cluster values
from a `merMod` fit, where `fsL`/`fsT` are 3-D arrays (one slice per
cluster). Pooling, rather than fitting each pattern as its own tiny
stage-2 sub-group, keeps small (possibly near-empty) patterns from
making the measurement model under-identified or numerically fragile.
The default `reduce` is `"mean"`, a convex combination of positive
semi-definite matrices, so the pooled `fsT` stays positive
semi-definite; the opt-in `"median"` trades that guarantee for
robustness and emits a warning when the pooled `fsT` is not positive
semi-definite. For homogeneous inputs the reduction is a no-op, so
complete-data behavior is unchanged.

### Derivation from a [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) result

`tspa()` silently derives its measurement inputs from a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result passed as `data`, so the canonical call is
`tspa(model, data = get_fs(...))`. Two derivations exist, in this order,
and explicit arguments always win over them:

- *Multi-factor* — fires when `fsT`, `fsL`, and `se_fs` are all omitted
  and `data` carries both `fsT` and `fsL` attributes; those attributes
  (and the `fsb` attribute when present) become the measurement inputs.
  Derivation is provenance-gated: `data` must actually resolve as a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result (it carries the `fs_pattern` attribute, or a `merMod` 3-D /
  `mirt` per-observation shape), so a hand-rolled data frame with plain
  matrix `fsT`/`fsL` attributes but no such provenance is not derived
  and is rejected with an informative error. When both forms are
  available and nothing is passed, this (attribute) form wins.

- *Single-factor* — fires when `fsT` is still `NULL` after that and
  `se_fs` was omitted; `se_fs` is then built from the data's own
  `fs_<v>` score columns that carry a matching numeric `fs_<v>_se`
  column (simple names only; product-score columns `fs_<v1>:<v2>` are
  not derived). This is the path taken by a
  [`cbind()`](https://rdrr.io/r/base/cbind.html)'d
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result, because [`cbind()`](https://rdrr.io/r/base/cbind.html) drops
  the attributes and leaves only the score/SE columns.

A supplied `se_fs` (even an empty
[`list()`](https://rdrr.io/r/base/list.html)) suppresses the
multi-factor derivation, keeping the single-factor path.

In the single-factor derivation with a group column, the derived `se_fs`
carries one row per group, in the order the group column's values first
appear in the data (lavaan's `group =` order for both character and
factor columns); a within-group-constant SE column keeps its constant
and a varying column (FIML missing data) is reduced by `reduce`,
idempotent with the FIML pooling that may run on the result. The group
column is the data's `group_col` attribute (if that column exists), else
the `group =` argument (if it names a column of the data), else a
literal `group` column. A
[`cbind()`](https://rdrr.io/r/base/cbind.html)'d multi-group frame with
no such group signal (no `group_col` attribute —
[`cbind()`](https://rdrr.io/r/base/cbind.html) drops it — and no
`group =` argument) is silently fitted as single-group with row-mean
SEs, so pass `group =` to control the stage-2 grouping.

The multi-factor derivation also picks up the `fsb` (intercept)
attribute, so a derived fit includes the stage-2 intercept constraints.
On `std.lv = TRUE` data the derived intercepts are zero (the default
`regression`/`Bartlett` scores), so the derived fit's free estimates and
standard errors equal the corresponding values of the explicit no-`fsb`
form (the intercept rows are simply fixed); the differences are that its
`tspaModel` string carries an extra intercept block and that, for
multigroup fits (where lavaan enables a mean structure), the no-`fsb`
form additionally estimates the factor-score intercepts freely. On data
with nonzero `fsb` (e.g. a `mirt` fit) even the free estimates differ
from the no-`fsb` form.

### Product-score auto-compute (`product = TRUE`)

With `product = TRUE`, a model latent that names the product of two of
the model's factor scores is treated as a latent interaction measured by
the double-mean-centered product indicator of the two scores. The pair
may be named by concatenation (`xm` for `x` and `m`) or in lavaan's
interaction syntax (`x:m`); the interaction-syntax form is rewritten to
the concatenated render name, because in the generated model `x:m` would
be parsed by lavaan as an interaction of the (latent) variables. An
`a:b` token whose parts are not both factor scores (e.g. `x:g` with `g`
an observed covariate) is not claimed and is passed through to lavaan as
an ordinary interaction. The product columns are computed on the fly
when absent
([`compute_fs_prod()`](https://mmm-lab-um.github.io/R2spa/reference/compute_fs_prod.md)
from the data's stage-1 attributes; pre-existing `fs_a:fs_b` columns —
in either orientation — are used as-is, so a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result with `product` set works too) and the product is wired into the
stage-2 measurement model: in the single-factor path the product SE
joins `se_fs` (per-group pooled by `reduce`, the same convention as the
score SEs; an explicit product SE may be keyed by either the render name
`xm` or the token `x:m`), and in the multi-factor path the product
latent gets a fixed loading `gamma` and fixed error variance `se_P^2`
from the (pooled) `fsL`/`fsT`/`psi`. Rejected with an informative error
(v1): multigroup models, `corrected_se = TRUE`, and data without stage-1
attributes that lacks the product columns (e.g. a
[`cbind()`](https://rdrr.io/r/base/cbind.html)'d
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result — compute the product columns up front or pass the
un-[`cbind()`](https://rdrr.io/r/base/cbind.html)'d result). A model
variable matching two different factor-score pairs, the same pair named
twice (`x:m` and `xm`), or a render name colliding with another model
variable, is an error.

### OpenMx route and mean structure

The OpenMx counterpart is
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md).
Both routes fit the same model (the covariance/structural quantities
agree to optimizer tolerance); they differ only in the unidentifiable
mean split between the corrected latents and their (observed)
factor-score indicators: `tspa()` fixes the exogenous latent mean at
zero and estimates the factor-score mean, whereas
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
does the opposite. Compare the two on the covariance quantities, not on
how the mean is split (see
[`tspa_mx_model()`](https://mmm-lab-um.github.io/R2spa/reference/tspa_mx_model.md)
under "Mean structure").

### Re-fitting with `lavaan::update()`

The returned fit carries a self-contained `@call`, so `lavaan::update()`
works on it as it does on a hand-written
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html) fit. In
particular `update(fit, meanstructure = TRUE)` re-fits the stage-2 model
with a mean structure, equivalent to calling
`tspa(..., meanstructure = TRUE)` from the start.

## See also

- `vignette("Two-Stage Path Analysis (2S-PA) Model Examples", package = "R2spa")`
  for an end-to-end walkthrough.

- `vignette("Multi-Factor Measurement Model", package = "R2spa")` for
  multi-factor measurement models.

- `vignette("Linear Growth Modeling with Two-Stage Path Analysis", package = "R2spa")`
  for growth models.

- `vignette("2S-PA with Missing Data", package = "R2spa")` for missing
  data (`missing = "fiml"`).

- `vignette("Product factor-score indicators (latent interactions)", package = "R2spa")`
  for `product = TRUE`.

- `vignette("Corrected Standard Errors", package = "R2spa")` for
  `corrected_se = TRUE`.

## Examples

``` r
library(lavaan)

# single-group, two-factor example, factor scores obtained separately
# get factor scores
fs_dat_ind60 <- get_fs(object = PoliticalDemocracy,
                       model = "ind60 =~ x1 + x2 + x3")
fs_dat_dem60 <- get_fs(object = PoliticalDemocracy,
                       model = "dem60 =~ y1 + y2 + y3 + y4")
fs_dat <- cbind(fs_dat_ind60, fs_dat_dem60)
# tspa model
tspa(model = "dem60 ~ ind60", data = fs_dat,
     se_fs = c(ind60 = fs_dat_ind60[1, "fs_ind60_se"],
               dem60 = fs_dat_dem60[1, "fs_dem60_se"]))
#> lavaan 0.7-2 ended normally after 17 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         3
#> 
#>   Number of observations                            75
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0

# the same fit with se_fs derived from the fs_<v>_se columns of the
# cbind()ed get_fs() results (no explicit se_fs)
tspa(model = "dem60 ~ ind60", data = fs_dat)
#> lavaan 0.7-2 ended normally after 17 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         3
#> 
#>   Number of observations                            75
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0

# single-group, three-factor example
mod2 <- "
  # latent variables
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
"
fs_dat2 <- get_fs(PoliticalDemocracy, model = mod2, std.lv = TRUE)
tspa(model = "dem60 ~ ind60
              dem65 ~ ind60 + dem60",
     data = fs_dat2,
     fsT = attr(fs_dat2, "fsT"),
     fsL = attr(fs_dat2, "fsL"))
#> lavaan 0.7-2 ended normally after 24 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         6
#> 
#>   Number of observations                            75
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0

# the same fit with the measurement inputs derived from the get_fs()
# result (no explicit fsT/fsL)
tspa(model = "dem60 ~ ind60
              dem65 ~ ind60 + dem60",
     data = fs_dat2)
#> lavaan 0.7-2 ended normally after 24 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         6
#> 
#>   Number of observations                            75
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 3
#>   P-value (Chi-square)                           1.000

# multigroup, two-factor example
mod3 <- "
  # latent variables
    visual =~ x1 + x2 + x3
    speed =~ x7 + x8 + x9
"
fs_dat3 <- get_fs(HolzingerSwineford1939, model = mod3, std.lv = TRUE,
                  group = "school")
tspa(model = "visual ~ speed",
     data = fs_dat3,
     fsT = attr(fs_dat3, "fsT"),
     fsL = attr(fs_dat3, "fsL"),
     group = "school")
#> lavaan 0.7-2 ended normally after 28 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        10
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0
#>   Test statistic for each group:
#>     Pasteur                                      0.000
#>     Grant-White                                  0.000

# the same fit with the measurement inputs derived from the get_fs()
# result (no explicit fsT/fsL)
tspa(model = "visual ~ speed",
     data = fs_dat3,
     group = "school")
#> lavaan 0.7-2 ended normally after 28 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         6
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 4
#>   P-value (Chi-square)                           1.000
#>   Test statistic for each group:
#>     Pasteur                                      0.000
#>     Grant-White                                  0.000

# multigroup, three-factor example
mod4 <- "
  # latent variables
    visual =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed =~ x7 + x8 + x9
"
fs_dat4 <- get_fs(HolzingerSwineford1939, model = mod4, std.lv = TRUE,
                  group = "school")
tspa(model = "visual ~ speed
              textual ~ visual + speed",
     data = fs_dat4,
     fsT = attr(fs_dat4, "fsT"),
     fsL = attr(fs_dat4, "fsL"),
     group = "school")
#> lavaan 0.7-2 ended normally after 39 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        18
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0
#>   Test statistic for each group:
#>     Pasteur                                      0.000
#>     Grant-White                                  0.000

# get factor scores
fs_dat_visual <- get_fs(object = HolzingerSwineford1939,
                        model = "visual =~ x1 + x2 + x3",
                        group = "school",
                        format = "list")
fs_dat_speed <- get_fs(object = HolzingerSwineford1939,
                       model = "speed =~ x7 + x8 + x9",
                       group = "school",
                       format = "list")
fs_hs <- cbind(do.call(rbind, fs_dat_visual),
               do.call(rbind, fs_dat_speed))

# tspa model
tspa(model = "visual ~ speed",
     data = fs_hs,
     se_fs = data.frame(visual = c(0.3391326, 0.311828),
                        speed = c(0.2786875, 0.2740507)),
     group = "school",
     group.equal = "regressions")
#> lavaan 0.7-2 ended normally after 19 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        10
#>   Number of equality constraints                     1
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.014
#>   Degrees of freedom                                 1
#>   P-value (Chi-square)                           0.907
#>   Test statistic for each group:
#>     Pasteur                                      0.010
#>     Grant-White                                  0.003

# manually adding equality constraints on the regression coefficients
tspa(model = "visual ~ c(b1, b1) * speed",
     data = fs_hs,
     se_fs = list(visual = c(0.3391326, 0.311828),
                  speed = c(0.2786875, 0.2740507)),
     group = "school")
#> lavaan 0.7-2 ended normally after 19 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        10
#>   Number of equality constraints                     1
#> 
#>   Number of observations per group:                   
#>     Pasteur                                        156
#>     Grant-White                                    145
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.014
#>   Degrees of freedom                                 1
#>   P-value (Chi-square)                           0.907
#>   Test statistic for each group:
#>     Pasteur                                      0.010
#>     Grant-White                                  0.003

# Missing data (FIML): per-pattern fsL/fsT are pooled within the group
data("HolzingerSwineford1939", package = "lavaan")
hs <- HolzingerSwineford1939
set.seed(1334)
hs$x2[!rbinom(nrow(hs), 1, 0.4)] <- NA
hs$x8[!rbinom(nrow(hs), 1, 0.4)] <- NA
mod_fin <- "
  visual =~ x1 + x2 + x3
  speed  =~ x7 + x8 + x9
"
fit_fin <- suppressWarnings(cfa(mod_fin, data = hs, missing = "fiml"))
fs_fin <- get_fs(fit_fin)
tspa("visual ~ speed", data = fs_fin,
     fsT = attr(fs_fin, "fsT"), fsL = attr(fs_fin, "fsL"),
     reduce = "mean")
#> lavaan 0.7-2 ended normally after 16 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         3
#> 
#>   Number of observations                           301
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0
# opt-in element-wise reduction (may lose positive semi-definiteness)
suppressWarnings(tspa("visual ~ speed", data = fs_fin,
     fsT = attr(fs_fin, "fsT"), fsL = attr(fs_fin, "fsL"),
     reduce = "median"))
#> lavaan 0.7-2 ended normally after 17 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         3
#> 
#>   Number of observations                           301
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0

# merMod: per-cluster fsL/fsT are pooled (one value per cluster)
library(lme4)
lmod <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fs_mer <- get_fs(lmod)
tspa("u1 ~ u0", data = fs_mer,
     fsT = attr(fs_mer, "fsT"), fsL = attr(fs_mer, "fsL"))
#> lavaan 0.7-2 ended normally after 51 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         3
#> 
#>   Number of observations                            18
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0

# Product-score auto-compute (opt-in): the model's `x:m` latent is the
# product of the `x` and `m` scores (rendered under the concatenated
# name `xm`), computed on the fly from the data's stage-1 attributes
# (no pre-computed product columns needed)
set.seed(2116)
covx <- matrix(c(1, 0.4, 0.4, 1), 2)
eta <- as.data.frame(MASS::mvrnorm(500, rep(0, 2), covx))
names(eta) <- c("x", "m")
lk <- list(x = c(0.9, 0.8, 0.7), m = c(0.85, 0.75, 0.65),
           y = c(0.75, 0.7, 0.65))
etay <- 0.5 * eta$x + 0.4 * eta$m + 0.3 * eta$x * eta$m
obs <- setNames(lapply(c("x", "m"), function(v0) {
  eta[[v0]] %*% t(lk[[v0]]) + rnorm(1500)
}), c("x", "m"))
obs$y <- etay %*% t(lk$y) + rnorm(1500)
df <- as.data.frame(do.call(cbind, obs))
names(df) <- c(paste0("x", 1:3), paste0("m", 1:3), paste0("y", 1:3))
fs_prod <- get_fs(df, model = "x =~ x1 + x2 + x3
                    m =~ m1 + m2 + m3
                    y =~ y1 + y2 + y3", std.lv = TRUE)
tspa("y ~ x + m + x:m", data = fs_prod, product = TRUE)
#> Warning: lavaan->lav_object_post_check():  
#>    some estimated lv variances are negative
#> lavaan 0.7-2 ended normally after 38 iterations
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                        11
#> 
#>   Number of observations                           500
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 3
#>   P-value (Chi-square)                           1.000
```
