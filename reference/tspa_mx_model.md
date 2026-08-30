# Two-Stage Path Analysis (OpenMx)

Fit a two-stage path analysis (2S-PA) model in OpenMx. This is the
OpenMx counterpart of
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md): the
structural path model is expressed on the (true) latent factors, each
factor score is a single indicator with known measurement error, and —
unlike `lavaan::tspa()` — the measurement quantities (loadings / error
variances / intercepts) may be fixed per-group constants *or*
per-observation definition-variable columns. The per-row
definition-variable form is the exact (non-pooled) correction that
`lavaan::tspa(reduce = )` only approximates.

## Usage

``` r
tspa_mx_model(
  model,
  data,
  se_fs = NULL,
  fsL = NULL,
  fsT = NULL,
  fsb = NULL,
  ...
)
```

## Arguments

- model:

  A character string describing the structural path model in `lavaan`
  syntax, using the **latent** (factor) names. Phase 1 restricts every
  variable in `model` to a corrected latent (one that has a factor
  score). Latent variances are added automatically, so do not declare
  them here.

- data:

  A data frame carrying the factor-score columns (`fs_<latent>`) and,
  for definition-variable entries, the per-observation columns they
  reference. A
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result works directly: with `se_fs`/`fsL`/`fsT`/`fsb` omitted, the
  measurement inputs are derived from its attributes (see Details), and
  the `int_fs_*` intercept columns are appended automatically —
  [`fs_indiv()`](https://mmm-lab-um.github.io/R2spa/reference/fs_indiv.md)
  is no longer needed to obtain them.
  [`fs_indiv()`](https://mmm-lab-um.github.io/R2spa/reference/fs_indiv.md)
  on a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result produces the equivalent fully explicit table.
  Definition-variable columns must be numeric and free of `NA`.

- se_fs:

  A named numeric vector of standard errors (one per latent) for the
  single-score-per-latent case; implies fixed unit loadings and error
  variances `se_fs^2`. An explicit `se_fs` always wins over derivation;
  when omitted (along with `fsL`, `fsT`, and `fsb`), the measurement
  inputs are derived from a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result passed as `data` (see Details).

- fsL:

  A `q x p` loading matrix including cross-loadings: rows = score names
  (`fs_<latent>`), columns = latent names. The matrix must be uniformly
  numeric (every cell a fixed loading) or uniformly character (every
  cell a definition-variable column name); mixing fixed values and
  column names in one matrix is not supported. Or omitted, in which case
  the value is derived from the `fsL` attribute of a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result passed as `data` (see Details); an explicit `fsL` always wins.

- fsT:

  A `q x q` error variance-covariance matrix over the score names; the
  lower triangle (incl. diagonal) is used, and every score must have an
  error variance (a complete diagonal). The matrix must be uniformly
  numeric (every cell fixed) or uniformly character (every cell a
  definition-variable column name); mixing fixed values and column names
  in one matrix is not supported. Or omitted, in which case the value is
  derived from the `fsT` attribute of a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result passed as `data` (see Details); an explicit `fsT` always wins.

- fsb:

  A vector of score intercepts (length `q`, named by score, either
  order), uniformly numeric (every entry fixed) or uniformly character
  (every entry a definition-variable column name); mixing fixed values
  and column names is not supported. `NULL` (default) fixes all score
  intercepts at zero. Or omitted, in which case the value is derived
  from the `fsb` attribute of a
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result passed as `data` (see Details); the derivation omits it — fixed
  zero intercepts — when the result carries no `fsb` attribute.

- ...:

  Additional arguments passed on to
  [`OpenMx::mxRun()`](https://rdrr.io/pkg/OpenMx/man/mxRun.html) (e.g.
  `intervals = TRUE`).

## Value

A fitted `OpenMx` `MxModel`.
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`summary()`](https://rdrr.io/r/base/summary.html) work as usual.

## Details

The internal design is a single-level RAM model (no sub-model, no
`umx`): the lavaan structural string is parsed with
[`lavaan::lavaanify()`](https://rdrr.io/pkg/lavaan/man/model.syntax.html),
the corrected latents are given auto latent variances, the score
indicators and their errors/intercepts are attached per
[tspa()'s](https://www.rdocumentation.org/packages/R2spa) schema, and
the whole thing is fit with `mxFitFunctionML()` (raw-data FIML).

### Auto-derivation from a [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) result

`tspa_mx_model()` derives its measurement inputs from a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result passed as `data`, so the canonical call is
`tspa_mx_model(model, data = get_fs(...))`. Derivation fires only when
all of `se_fs`, `fsL`, `fsT`, and `fsb` are omitted; explicit arguments
always win. Derivation is provenance-gated: `data` must resolve as a
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result, so a hand-rolled frame with plain matrix `fsT`/`fsL` attributes
but no such provenance is rejected with an informative error.

Attributes are dispatched by shape: constant quantities (plain matrices
— complete-data, `local = TRUE`, and `format = "list"` results; a plain
`fsb` vector) become fixed numeric matrices, while per-row quantities
(`mirt_per_obs`/`per_obs`-marked results), per-cluster quantities
(`merMod` results: 3-D `fsL`/`fsT` arrays, one row per cluster), and
per-pattern quantities (single-group FIML, keyed by `fs_pattern$label`)
become definition-variable matrices referencing the result's own
`*_by_*`, `ev_*`, and `ecov_*` columns. A `merMod` result carries no
`fsb` attribute, so its score intercepts stay fixed at zero. A
non-`NULL` `fsb` attribute appends `int_fs_*` intercept columns to a
working copy of `data` (the `fs_indiv(include_intercept = TRUE)`
equivalent); `NULL` keeps the default fixed-zero intercepts. Unlike
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) (whose
`reduce = ` argument pools per-unit quantities), there is no pooling
here — the OpenMx route is exact-or-fail. A multigroup result (the
`group_col` attribute) is refused with the Phase-1 message; a `mirt`
multigroup result (a `group` column, no `group_col` attribute) derives
as a single pooled per-row-corrected fit (no per-group structural
parameters; the `group` column is inert). A `get_fs(product = )` result
derives identically: its extra `fs_a:fs_b` (and `_se`/`_ld`) columns are
inert to derivation, but the `:` they carry is illegal in
[`OpenMx::mxData()`](https://rdrr.io/pkg/OpenMx/man/mxData.html) column
names, so the frame is un-fittable until the product columns are dropped
(a pre-existing limitation that affects the explicit-argument route
too).

### Mean structure

Both routes fit the same model, so the covariance/structural quantities
(paths, latent variances/covariances, and their SEs) agree to optimizer
tolerance. The two differ only in the unidentifiable split of the mean
structure between the corrected latents and their (observed)
factor-score indicators:
[`tspa()`](https://mmm-lab-um.github.io/R2spa/reference/tspa.md) fixes
the exogenous latent mean at zero and lets the factor-score mean
estimate the data value, whereas this OpenMx route fixes the score
residual means at zero and estimates the latent means (the latent mean
carries the data value). Read the means accordingly — compare the two
routes on the covariance quantities, not on how the mean is split
between latent and indicator. (In the no-mean-structure case the
difference is flat: the split is arbitrary and only the total is
identifiable.)

## See also

- `vignette("2S-PA with OpenMx and IRT (mirt)", package = "R2spa")` for
  the OpenMx route and IRT stage 1.

- `vignette("Two-Stage Path Analysis (2S-PA) Model Examples", package = "R2spa")`
  for the lavaan route.

## Examples

``` r
if (FALSE) { # \dontrun{
## Measurement inputs derived from a get_fs() result:
fs <- get_fs(PoliticalDemocracy, "dem60 =~ y1 + y2 + y3 + y4
                                  ind60 =~ x1 + x2 + x3")
# measurement inputs derived from the get_fs() result
tspa_mx_model("dem60 ~ ind60; dem60 + ind60 ~ 1", data = fs)

## Equivalent, fully explicit (per-row definition variables via
## fs_indiv()):
dat <- fs_indiv(fs, include_intercept = TRUE)
tspa_mx_model("dem60 ~ ind60; dem60 + ind60 ~ 1",
  data = dat,
  fsL = matrix(c("ind60_by_fs_ind60", "ind60_by_fs_dem60",
                 "dem60_by_fs_ind60", "dem60_by_fs_dem60"),
               nrow = 2, dimnames = list(c("fs_ind60", "fs_dem60"),
                                         c("ind60", "dem60"))),
  fsT = matrix(c("ev_fs_ind60", "ecov_fs_ind60_fs_dem60", NA,
                 "ev_fs_dem60"),
               nrow = 2, dimnames = list(c("fs_ind60", "fs_dem60"),
                                         c("fs_ind60", "fs_dem60"))),
  fsb = c(fs_ind60 = "int_fs_ind60", fs_dem60 = "int_fs_dem60"))
} # }
```
