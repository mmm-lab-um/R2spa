# Get Factor Scores and the Corresponding Standard Error of Measurement

[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md) is
an S3 generic that extracts factor scores from fitted models. Methods
are available for `data.frame` (fits a CFA internally), `lavaan`
objects, `lmerMod` objects, and fitted `mirt` models (single-group
`SingleGroupClass` and multi-group `MultipleGroupClass`; `mirt` is a
`Suggests` dependency). Multi-group mirt results carry a trailing
`group` column and a per-group (`list`) `psi` attribute.

## Usage

``` r
get_fs_lavaan(
  lavobj,
  method = c("regression", "Bartlett", "ML", "EB", "mean"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  prior_mean = NULL,
  prior_cov = NULL,
  sum_items = NULL,
  ...
)
```

## Arguments

- lavobj:

  A lavaan model object when using `get_fs_lavaan()`.

- method:

  Character. Method for computing factor scores. For `lavaan` and data
  frame objects: `"regression"` (default, consistent with
  [`lavPredict`](https://rdrr.io/pkg/lavaan/man/lavPredict.html)),
  `"Bartlett"`, or `"mean"` (a third, distinct method: sum scores, each
  score being the plain uncentered mean of the items assigned to its
  factor, using no latent distribution), with `"ML"` an alias for
  `"Bartlett"` and `"EB"` an alias for `"regression"`. For `merMod`
  objects: `"EB"` (empirical Bayes, default; identical to the first
  random-effect term's
  [`ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html)()
  estimates) or `"ML"` (a prior-free, per-cluster OLS estimate of the
  random effects, using no random-effects prior, analogous to Bartlett
  scores for `lavaan` objects). The `"ML"`/`"EB"` aliases apply to the
  lavaan path only; for `merMod` objects the two strings are distinct
  methods. Bartlett scores have more desirable properties than
  regression scores and may be preferred for 2S-PA. `method = "mean"`
  takes the item-to-factor assignment from `sum_items` (auto-derived
  from the estimated loadings when `NULL`); it errors when the model was
  fitted with missing data retained (e.g. FIML/digamma), and is not
  supported together with `corrected_fsT`, `vfsLT`, `reliability`,
  `prior_mean`, or `prior_cov`.

- corrected_fsT:

  Logical. Whether to correct for the sampling error in the factor score
  weights when computing the error variance estimates of factor scores.
  Currently ignored for `merMod` objects.

- vfsLT:

  Logical. Whether to return the covariance matrix of `fsT` and `fsL`,
  returned as attribute `vfsLT`; used for second-order SE correction of
  2S-PA results. Currently ignored for `merMod` objects.

- reliability:

  Logical. Whether to return the reliability of factor scores. Available
  only for single-factor lavaan models; for multi-factor models a
  warning is issued and no `reliability` attribute is returned.

- prior_mean:

  An optional numeric vector of length `q` (the number of latent
  variables) giving fixed external prior means for the latent variables.
  `NULL` (default) uses the lavaan-estimated (group-specific) latent
  means. Non-NULL values are treated as fixed external priors shared
  across all lavaan groups. For `mirt` `SingleGroupClass` objects it
  instead sets the factor prior mean used for the EAP scores; the
  factor-score intercepts (`fsb`) then vary per observation as
  `Vpost_i %*% solve(psi) %*% prior_mean`, i.e. the latent mean scaled
  by the per-observation shrinkage factor (zero when
  `prior_mean = NULL`), where `psi` is the mirt model's estimated factor
  covariance. For `mirt` `MultipleGroupClass` objects a non-NULL
  `prior_mean` (length `q`) is applied as the factor prior mean to every
  group (mirt's per-group EAP is otherwise centred on a zero-mean
  standard-normal prior); each observation's regression form uses the
  factor covariance of its own group. Only supported for lavaan objects
  with regression (EB) scoring (and for mirt); `reliability = TRUE` is
  not supported together with user-supplied `prior_mean`/`prior_cov`,
  and `prior_cov` is not supported for mirt. Conceptually similar to the
  `mean` argument of
  [`mirt::fscores()`](https://philchalmers.github.io/mirt/reference/fscores.html).

- prior_cov:

  An optional numeric `q x q` covariance matrix (a scalar or 1 x 1
  matrix is accepted when `q = 1`) giving fixed external prior
  covariance for the latent variables. `NULL` (default) uses the
  lavaan-estimated (group-specific) latent covariance. Non-NULL values
  must be finite, symmetric and positive definite; when `q > 1` the
  matrix must be named (row and column names matching the latent
  variable names), so its entries map unambiguously onto the model's
  latent variables. Values are treated as fixed external priors shared
  across all lavaan groups. Only supported for lavaan objects with
  regression (EB) scoring; `reliability = TRUE` is not supported
  together with user-supplied `prior_mean`/`prior_cov`. With
  `corrected_fsT = TRUE` or `vfsLT = TRUE` the supplied covariance is
  treated as fixed, i.e. no sampling uncertainty from the prior itself
  is propagated. Conceptually similar to the `cov` argument of
  [`mirt::fscores()`](https://philchalmers.github.io/mirt/reference/fscores.html).

- sum_items:

  A named list mapping each factor name to the item names that make up
  its sum score, e.g.
  `list(ind60 = c("x1", "x2", "x3"), dem60 = c("y1", "y2", "y3", "y4"))`.
  `NULL` (default) auto-derives the assignment from the estimated
  loadings, which requires each indicator to load on exactly one factor
  and every factor to have at least one item. A supplied list must cover
  all model factors, and each item may belong to only one sum. Only used
  for `lavaan` and data frame objects with `method = "mean"`.

- ...:

  additional arguments passed to
  [`cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html) (when `object` is a
  data frame). See
  [`lavOptions`](https://rdrr.io/pkg/lavaan/man/lavOptions.html) for a
  complete list.

## Value

A data frame containing the factor scores (with prefix `"fs_"`), the
standard errors (with suffix `"_se"`), the implied loadings of indicator
`_by_` factor scores, and the error variance-covariance of the factor
scores (with prefix `"ev_"` or `"ecov_"`). For multi-group lavaan models
in `"unified"` format, a `group` column is included. The following
attributes are attached: \* `fsT`: error covariance of factor scores
(matrix or named list by group) \* `fsL`: loading matrix of factor
scores (matrix or named list by group) \* `fsb`: intercepts of factor
scores (vector or named list by group); with `method = "mean"` the
intercept is the mean of the factor's item intercepts, the measurement
intercept of the score regressed on the uncentered latent (same
`E[fs] - fsL %*% alpha` convention as the other methods; equals the
score's column mean for models without a mean structure) \*
`scoring_matrix`: weights for computing factor scores from the observed
data, as a named list. For lavaan models: one score x item matrix per
group; with `method = "mean"` the weights are the item-mean weights, so
`S %*% y` reproduces the raw scores exactly (no centering offset). For
`merMod` models: one `num_re` x `n_j` matrix per cluster, where
`S_j %*% (y_j - X_j %*% beta)` with `y_j`/`X_j` the cluster's rows of
the model response and the fixed-effects design reproduces the cluster's
EB scores for method `"EB"` and the per-cluster OLS (ML) scores for
method `"ML"`. \* `psi`: effective (prior-adjusted) covariance matrix of
the latent variables (`q x q`), group-level (not per-pattern), and a
point estimate only (no sampling SEs of the latents are attached).
Mirrors the `fsT` shape: a named list keyed by group label for
`"unified"` output; a direct attribute on each group data frame (plus a
list-valued attribute on the outer list) for `"list"` output; for
`merMod` objects a single `q x q` matrix. With `prior_cov` supplied it
equals the prior (shared across groups), otherwise the per-group lavaan
estimate. For `merMod` objects the matrix is the first random-effects
term's `VarCorr`, with dimnames renamed to match the `fsL` column names
(`u0`/`u1`/..., or the legacy `u0_eb`/`u1_eb` names). \* `alpha`:
effective (prior-adjusted) means of the latent variables (a named vector
of length `q`), with the same group nesting and point-estimate semantics
as `psi`. With `prior_mean` supplied it equals the prior, otherwise the
per-group lavaan estimate; a named zero vector (`0` per latent) when the
model has no (estimated) mean structure. For `merMod` objects a named
zero vector (random effects are mean zero). \* `fs_pattern`: for lavaan
models, a named list by group of `list(label, pat)` entries. `label` is
a character vector with one entry per case in the group giving that
case's observed-indicator pattern name (`NA` for cases whose indicators
are all missing); `pat` is a logical matrix with rows = indicators and
one column per pattern, the columns being named by pattern name.

        For a lavaan group without missing data, its `fsT`/`fsL`/`fsb`/
        `scoring_matrix` elements are the plain matrix/vector for the whole
        group. When a group's cases split into multiple observed-indicator
        patterns (missing data), each such element is instead a named list
        with one entry per pattern; the pattern name is the observed
        indicator names joined with `"+"` in indicator order (e.g.
        `"x1+x3"`).

        With `local = TRUE` and missing data (e.g. `missing = "fiml"`),
        the `fsT`/`fsL`/`fsb`/`scoring_matrix` attributes are instead
        per-row lists (one entry per data row) and the result carries a
        `per_obs = TRUE` attribute (the same convention as mirt's
        `mirt_per_obs`).

        Note: for a single-group lavaan fit in `"unified"` format, the
        per-group attribute wrappers (`fsT`, `fsL`, `fsb`,
        `scoring_matrix`, `psi`, and `alpha`) are each a one-element list
        named with the empty string `""`; `x[[""]]` does not match in R
        list subsetting, so read these attributes positionally (e.g.
        `attr(fs, "fsT")[[1]]`, `attr(fs, "psi")[[1]]`) rather than by
        name.

## Details

`get_fs_lavaan()` is superseded by
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md). It
is retained for backward compatibility and delegates to
`get_fs(object, format = "list")` internally. New code should call
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
directly.

## See also

- `vignette("Two-Stage Path Analysis (2S-PA) Model Examples", package = "R2spa")`
  for end-to-end stage-1/stage-2 examples.

- `vignette("Scoring Matrices: lavaan CFA and lme4", package = "R2spa")`
  for the scoring-matrix internals.

- `vignette("EFA Scores", package = "R2spa")` for EFA-based factor
  scores.

- `vignette("2S-PA with Missing Data", package = "R2spa")` for
  `missing = "fiml"`.

## Examples

``` r
library(lavaan)
get_fs(PoliticalDemocracy[c("x1", "x2", "x3")])
#>          fs_f1  fs_f1_se f1_by_fs_f1   ev_fs_f1
#> 1  -0.52616832 0.1213615   0.9657673 0.01472862
#> 2   0.14365274 0.1213615   0.9657673 0.01472862
#> 3   0.71435592 0.1213615   0.9657673 0.01472862
#> 4   1.23992565 0.1213615   0.9657673 0.01472862
#> 5   0.83190803 0.1213615   0.9657673 0.01472862
#> 6   0.21238453 0.1213615   0.9657673 0.01472862
#> 7   0.11880855 0.1213615   0.9657673 0.01472862
#> 8   0.11322703 0.1213615   0.9657673 0.01472862
#> 9   0.25617279 0.1213615   0.9657673 0.01472862
#> 10  0.37112496 0.1213615   0.9657673 0.01472862
#> 11  0.67281395 0.1213615   0.9657673 0.01472862
#> 12  0.56885577 0.1213615   0.9657673 0.01472862
#> 13  1.31369791 0.1213615   0.9657673 0.01472862
#> 14  0.22042629 0.1213615   0.9657673 0.01472862
#> 15  0.57849228 0.1213615   0.9657673 0.01472862
#> 16  0.37805983 0.1213615   0.9657673 0.01472862
#> 17  0.05734046 0.1213615   0.9657673 0.01472862
#> 18 -0.01609202 0.1213615   0.9657673 0.01472862
#> 19  0.88923616 0.1213615   0.9657673 0.01472862
#> 20  1.11445897 0.1213615   0.9657673 0.01472862
#> 21  0.94657339 0.1213615   0.9657673 0.01472862
#> 22  0.90122770 0.1213615   0.9657673 0.01472862
#> 23  0.58409450 0.1213615   0.9657673 0.01472862
#> 24  0.64089192 0.1213615   0.9657673 0.01472862
#> 25  0.91021968 0.1213615   0.9657673 0.01472862
#> 26 -0.89660969 0.1213615   0.9657673 0.01472862
#> 27 -0.13195991 0.1213615   0.9657673 0.01472862
#> 28 -0.52968769 0.1213615   0.9657673 0.01472862
#> 29 -0.81799629 0.1213615   0.9657673 0.01472862
#> 30 -1.27199371 0.1213615   0.9657673 0.01472862
#> 31 -0.32096024 0.1213615   0.9657673 0.01472862
#> 32 -1.16780103 0.1213615   0.9657673 0.01472862
#> 33 -0.12295473 0.1213615   0.9657673 0.01472862
#> 34 -0.04285945 0.1213615   0.9657673 0.01472862
#> 35 -0.34323505 0.1213615   0.9657673 0.01472862
#> 36 -0.60541633 0.1213615   0.9657673 0.01472862
#> 37  0.17688718 0.1213615   0.9657673 0.01472862
#> 38 -0.55066055 0.1213615   0.9657673 0.01472862
#> 39 -1.05988219 0.1213615   0.9657673 0.01472862
#> 40 -0.04138802 0.1213615   0.9657673 0.01472862
#> 41 -0.12611837 0.1213615   0.9657673 0.01472862
#> 42 -0.60322892 0.1213615   0.9657673 0.01472862
#> 43 -0.11057176 0.1213615   0.9657673 0.01472862
#> 44 -1.06423085 0.1213615   0.9657673 0.01472862
#> 45 -1.08354999 0.1213615   0.9657673 0.01472862
#> 46 -0.84009484 0.1213615   0.9657673 0.01472862
#> 47 -1.14678213 0.1213615   0.9657673 0.01472862
#> 48 -0.57578976 0.1213615   0.9657673 0.01472862
#> 49  0.07186692 0.1213615   0.9657673 0.01472862
#> 50  0.14682421 0.1213615   0.9657673 0.01472862
#> 51  0.35871830 0.1213615   0.9657673 0.01472862
#> 52 -0.43403195 0.1213615   0.9657673 0.01472862
#> 53  0.44603111 0.1213615   0.9657673 0.01472862
#> 54  0.26352000 0.1213615   0.9657673 0.01472862
#> 55  0.55051165 0.1213615   0.9657673 0.01472862
#> 56  0.23453122 0.1213615   0.9657673 0.01472862
#> 57  0.27968138 0.1213615   0.9657673 0.01472862
#> 58  0.70960640 0.1213615   0.9657673 0.01472862
#> 59  0.25227978 0.1213615   0.9657673 0.01472862
#> 60  1.18849297 0.1213615   0.9657673 0.01472862
#> 61  0.21104946 0.1213615   0.9657673 0.01472862
#> 62 -1.16516281 0.1213615   0.9657673 0.01472862
#> 63 -0.85560065 0.1213615   0.9657673 0.01472862
#> 64  0.13398476 0.1213615   0.9657673 0.01472862
#> 65 -0.07912189 0.1213615   0.9657673 0.01472862
#> 66 -0.27146711 0.1213615   0.9657673 0.01472862
#> 67 -0.04417217 0.1213615   0.9657673 0.01472862
#> 68 -1.33425662 0.1213615   0.9657673 0.01472862
#> 69 -0.38720750 0.1213615   0.9657673 0.01472862
#> 70 -0.55355511 0.1213615   0.9657673 0.01472862
#> 71 -0.72242623 0.1213615   0.9657673 0.01472862
#> 72  0.30607449 0.1213615   0.9657673 0.01472862
#> 73  0.77707950 0.1213615   0.9657673 0.01472862
#> 74  0.06847481 0.1213615   0.9657673 0.01472862
#> 75 -0.11052927 0.1213615   0.9657673 0.01472862

# Multiple factors
get_fs(PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")],
       model = " ind60 =~ x1 + x2 + x3
                 dem60 =~ y1 + y2 + y3 + y4 ")
#>       fs_ind60    fs_dem60 fs_ind60_se fs_dem60_se ind60_by_fs_ind60
#> 1  -0.54258816 -2.74640573   0.1245694   0.6307323         0.9553858
#> 2   0.12647664 -2.85646114   0.1245694   0.6307323         0.9553858
#> 3   0.73408891  2.74401728   0.1245694   0.6307323         0.9553858
#> 4   1.25253604  3.10856431   0.1245694   0.6307323         0.9553858
#> 5   0.83355267  1.92455641   0.1245694   0.6307323         0.9553858
#> 6   0.22426801  1.02292332   0.1245694   0.6307323         0.9553858
#> 7   0.12517739  1.00406461   0.1245694   0.6307323         0.9553858
#> 8   0.11783867 -0.37216403   0.1245694   0.6307323         0.9553858
#> 9   0.25175134 -1.24897911   0.1245694   0.6307323         0.9553858
#> 10  0.39938631  2.85267059   0.1245694   0.6307323         0.9553858
#> 11  0.67497777  1.41959595   0.1245694   0.6307323         0.9553858
#> 12  0.56462020  1.08769844   0.1245694   0.6307323         0.9553858
#> 13  1.31236592  1.54090232   0.1245694   0.6307323         0.9553858
#> 14  0.23246021  1.77370863   0.1245694   0.6307323         0.9553858
#> 15  0.58638481  2.45676871   0.1245694   0.6307323         0.9553858
#> 16  0.38404785  2.35887573   0.1245694   0.6307323         0.9553858
#> 17  0.05076465  0.04034088   0.1245694   0.6307323         0.9553858
#> 18 -0.01747337 -1.86718064   0.1245694   0.6307323         0.9553858
#> 19  0.90920762  3.61477756   0.1245694   0.6307323         0.9553858
#> 20  1.12553557  0.88355273   0.1245694   0.6307323         0.9553858
#> 21  0.97202590  3.62673300   0.1245694   0.6307323         0.9553858
#> 22  0.87820036 -3.02428925   0.1245694   0.6307323         0.9553858
#> 23  0.57540754 -1.51695438   0.1245694   0.6307323         0.9553858
#> 24  0.66221224  2.76341635   0.1245694   0.6307323         0.9553858
#> 25  0.92358281  2.00507336   0.1245694   0.6307323         0.9553858
#> 26 -0.89353051 -0.92008050   0.1245694   0.6307323         0.9553858
#> 27 -0.13984744 -1.19025576   0.1245694   0.6307323         0.9553858
#> 28 -0.53828496 -1.01247764   0.1245694   0.6307323         0.9553858
#> 29 -0.80834865  0.10456709   0.1245694   0.6307323         0.9553858
#> 30 -1.25324343 -0.71847055   0.1245694   0.6307323         0.9553858
#> 31 -0.33373641 -1.61401581   0.1245694   0.6307323         0.9553858
#> 32 -1.17441075 -3.27250363   0.1245694   0.6307323         0.9553858
#> 33 -0.12409974 -1.17530231   0.1245694   0.6307323         0.9553858
#> 34 -0.04239173 -0.53796274   0.1245694   0.6307323         0.9553858
#> 35 -0.34010528  0.74552889   0.1245694   0.6307323         0.9553858
#> 36 -0.58953870  1.61018662   0.1245694   0.6307323         0.9553858
#> 37  0.17453657 -0.28144814   0.1245694   0.6307323         0.9553858
#> 38 -0.54457243  0.37694690   0.1245694   0.6307323         0.9553858
#> 39 -1.05196602 -0.62919501   0.1245694   0.6307323         0.9553858
#> 40 -0.05504697 -0.03346842   0.1245694   0.6307323         0.9553858
#> 41 -0.12364358 -0.38394102   0.1245694   0.6307323         0.9553858
#> 42 -0.59058710  1.35347275   0.1245694   0.6307323         0.9553858
#> 43 -0.11968796  0.89227782   0.1245694   0.6307323         0.9553858
#> 44 -1.07176064 -2.08481096   0.1245694   0.6307323         0.9553858
#> 45 -1.09139097 -2.07944291   0.1245694   0.6307323         0.9553858
#> 46 -0.83287255  1.59590721   0.1245694   0.6307323         0.9553858
#> 47 -1.14519896 -1.53352201   0.1245694   0.6307323         0.9553858
#> 48 -0.56115378  2.08138051   0.1245694   0.6307323         0.9553858
#> 49  0.06493340 -1.04044137   0.1245694   0.6307323         0.9553858
#> 50  0.15671638  1.72618633   0.1245694   0.6307323         0.9553858
#> 51  0.34626130 -1.24967043   0.1245694   0.6307323         0.9553858
#> 52 -0.45158373 -2.31742576   0.1245694   0.6307323         0.9553858
#> 53  0.43233465 -1.07533341   0.1245694   0.6307323         0.9553858
#> 54  0.25779725 -0.02904676   0.1245694   0.6307323         0.9553858
#> 55  0.51730650 -2.78207923   0.1245694   0.6307323         0.9553858
#> 56  0.20104991 -2.49001474   0.1245694   0.6307323         0.9553858
#> 57  0.25318620 -2.52145147   0.1245694   0.6307323         0.9553858
#> 58  0.72354623  1.86717109   0.1245694   0.6307323         0.9553858
#> 59  0.24619740 -0.93321102   0.1245694   0.6307323         0.9553858
#> 60  1.21681210  3.19853937   0.1245694   0.6307323         0.9553858
#> 61  0.18167599 -3.15685030   0.1245694   0.6307323         0.9553858
#> 62 -1.16605067 -3.41334680   0.1245694   0.6307323         0.9553858
#> 63 -0.86491026 -3.11864398   0.1245694   0.6307323         0.9553858
#> 64  0.10990059 -0.47238885   0.1245694   0.6307323         0.9553858
#> 65 -0.07376176  2.95292007   0.1245694   0.6307323         0.9553858
#> 66 -0.28782931 -1.96509718   0.1245694   0.6307323         0.9553858
#> 67 -0.02508160  2.96218478   0.1245694   0.6307323         0.9553858
#> 68 -1.31843215 -1.59567027   0.1245694   0.6307323         0.9553858
#> 69 -0.40462357 -1.79146161   0.1245694   0.6307323         0.9553858
#> 70 -0.55568363 -1.01578892   0.1245694   0.6307323         0.9553858
#> 71 -0.71308015  0.08818212   0.1245694   0.6307323         0.9553858
#> 72  0.31014319  1.70765911   0.1245694   0.6307323         0.9553858
#> 73  0.79092897  1.86102556   0.1245694   0.6307323         0.9553858
#> 74  0.08770237  3.12885767   0.1245694   0.6307323         0.9553858
#> 75 -0.14138149 -2.41398025   0.1245694   0.6307323         0.9553858
#>    ind60_by_fs_dem60 dem60_by_fs_ind60 dem60_by_fs_dem60 ev_fs_ind60
#> 1           0.181827       0.005867694         0.8688887  0.01551752
#> 2           0.181827       0.005867694         0.8688887  0.01551752
#> 3           0.181827       0.005867694         0.8688887  0.01551752
#> 4           0.181827       0.005867694         0.8688887  0.01551752
#> 5           0.181827       0.005867694         0.8688887  0.01551752
#> 6           0.181827       0.005867694         0.8688887  0.01551752
#> 7           0.181827       0.005867694         0.8688887  0.01551752
#> 8           0.181827       0.005867694         0.8688887  0.01551752
#> 9           0.181827       0.005867694         0.8688887  0.01551752
#> 10          0.181827       0.005867694         0.8688887  0.01551752
#> 11          0.181827       0.005867694         0.8688887  0.01551752
#> 12          0.181827       0.005867694         0.8688887  0.01551752
#> 13          0.181827       0.005867694         0.8688887  0.01551752
#> 14          0.181827       0.005867694         0.8688887  0.01551752
#> 15          0.181827       0.005867694         0.8688887  0.01551752
#> 16          0.181827       0.005867694         0.8688887  0.01551752
#> 17          0.181827       0.005867694         0.8688887  0.01551752
#> 18          0.181827       0.005867694         0.8688887  0.01551752
#> 19          0.181827       0.005867694         0.8688887  0.01551752
#> 20          0.181827       0.005867694         0.8688887  0.01551752
#> 21          0.181827       0.005867694         0.8688887  0.01551752
#> 22          0.181827       0.005867694         0.8688887  0.01551752
#> 23          0.181827       0.005867694         0.8688887  0.01551752
#> 24          0.181827       0.005867694         0.8688887  0.01551752
#> 25          0.181827       0.005867694         0.8688887  0.01551752
#> 26          0.181827       0.005867694         0.8688887  0.01551752
#> 27          0.181827       0.005867694         0.8688887  0.01551752
#> 28          0.181827       0.005867694         0.8688887  0.01551752
#> 29          0.181827       0.005867694         0.8688887  0.01551752
#> 30          0.181827       0.005867694         0.8688887  0.01551752
#> 31          0.181827       0.005867694         0.8688887  0.01551752
#> 32          0.181827       0.005867694         0.8688887  0.01551752
#> 33          0.181827       0.005867694         0.8688887  0.01551752
#> 34          0.181827       0.005867694         0.8688887  0.01551752
#> 35          0.181827       0.005867694         0.8688887  0.01551752
#> 36          0.181827       0.005867694         0.8688887  0.01551752
#> 37          0.181827       0.005867694         0.8688887  0.01551752
#> 38          0.181827       0.005867694         0.8688887  0.01551752
#> 39          0.181827       0.005867694         0.8688887  0.01551752
#> 40          0.181827       0.005867694         0.8688887  0.01551752
#> 41          0.181827       0.005867694         0.8688887  0.01551752
#> 42          0.181827       0.005867694         0.8688887  0.01551752
#> 43          0.181827       0.005867694         0.8688887  0.01551752
#> 44          0.181827       0.005867694         0.8688887  0.01551752
#> 45          0.181827       0.005867694         0.8688887  0.01551752
#> 46          0.181827       0.005867694         0.8688887  0.01551752
#> 47          0.181827       0.005867694         0.8688887  0.01551752
#> 48          0.181827       0.005867694         0.8688887  0.01551752
#> 49          0.181827       0.005867694         0.8688887  0.01551752
#> 50          0.181827       0.005867694         0.8688887  0.01551752
#> 51          0.181827       0.005867694         0.8688887  0.01551752
#> 52          0.181827       0.005867694         0.8688887  0.01551752
#> 53          0.181827       0.005867694         0.8688887  0.01551752
#> 54          0.181827       0.005867694         0.8688887  0.01551752
#> 55          0.181827       0.005867694         0.8688887  0.01551752
#> 56          0.181827       0.005867694         0.8688887  0.01551752
#> 57          0.181827       0.005867694         0.8688887  0.01551752
#> 58          0.181827       0.005867694         0.8688887  0.01551752
#> 59          0.181827       0.005867694         0.8688887  0.01551752
#> 60          0.181827       0.005867694         0.8688887  0.01551752
#> 61          0.181827       0.005867694         0.8688887  0.01551752
#> 62          0.181827       0.005867694         0.8688887  0.01551752
#> 63          0.181827       0.005867694         0.8688887  0.01551752
#> 64          0.181827       0.005867694         0.8688887  0.01551752
#> 65          0.181827       0.005867694         0.8688887  0.01551752
#> 66          0.181827       0.005867694         0.8688887  0.01551752
#> 67          0.181827       0.005867694         0.8688887  0.01551752
#> 68          0.181827       0.005867694         0.8688887  0.01551752
#> 69          0.181827       0.005867694         0.8688887  0.01551752
#> 70          0.181827       0.005867694         0.8688887  0.01551752
#> 71          0.181827       0.005867694         0.8688887  0.01551752
#> 72          0.181827       0.005867694         0.8688887  0.01551752
#> 73          0.181827       0.005867694         0.8688887  0.01551752
#> 74          0.181827       0.005867694         0.8688887  0.01551752
#> 75          0.181827       0.005867694         0.8688887  0.01551752
#>    ecov_fs_dem60_fs_ind60 ev_fs_dem60
#> 1             0.005632564   0.3978232
#> 2             0.005632564   0.3978232
#> 3             0.005632564   0.3978232
#> 4             0.005632564   0.3978232
#> 5             0.005632564   0.3978232
#> 6             0.005632564   0.3978232
#> 7             0.005632564   0.3978232
#> 8             0.005632564   0.3978232
#> 9             0.005632564   0.3978232
#> 10            0.005632564   0.3978232
#> 11            0.005632564   0.3978232
#> 12            0.005632564   0.3978232
#> 13            0.005632564   0.3978232
#> 14            0.005632564   0.3978232
#> 15            0.005632564   0.3978232
#> 16            0.005632564   0.3978232
#> 17            0.005632564   0.3978232
#> 18            0.005632564   0.3978232
#> 19            0.005632564   0.3978232
#> 20            0.005632564   0.3978232
#> 21            0.005632564   0.3978232
#> 22            0.005632564   0.3978232
#> 23            0.005632564   0.3978232
#> 24            0.005632564   0.3978232
#> 25            0.005632564   0.3978232
#> 26            0.005632564   0.3978232
#> 27            0.005632564   0.3978232
#> 28            0.005632564   0.3978232
#> 29            0.005632564   0.3978232
#> 30            0.005632564   0.3978232
#> 31            0.005632564   0.3978232
#> 32            0.005632564   0.3978232
#> 33            0.005632564   0.3978232
#> 34            0.005632564   0.3978232
#> 35            0.005632564   0.3978232
#> 36            0.005632564   0.3978232
#> 37            0.005632564   0.3978232
#> 38            0.005632564   0.3978232
#> 39            0.005632564   0.3978232
#> 40            0.005632564   0.3978232
#> 41            0.005632564   0.3978232
#> 42            0.005632564   0.3978232
#> 43            0.005632564   0.3978232
#> 44            0.005632564   0.3978232
#> 45            0.005632564   0.3978232
#> 46            0.005632564   0.3978232
#> 47            0.005632564   0.3978232
#> 48            0.005632564   0.3978232
#> 49            0.005632564   0.3978232
#> 50            0.005632564   0.3978232
#> 51            0.005632564   0.3978232
#> 52            0.005632564   0.3978232
#> 53            0.005632564   0.3978232
#> 54            0.005632564   0.3978232
#> 55            0.005632564   0.3978232
#> 56            0.005632564   0.3978232
#> 57            0.005632564   0.3978232
#> 58            0.005632564   0.3978232
#> 59            0.005632564   0.3978232
#> 60            0.005632564   0.3978232
#> 61            0.005632564   0.3978232
#> 62            0.005632564   0.3978232
#> 63            0.005632564   0.3978232
#> 64            0.005632564   0.3978232
#> 65            0.005632564   0.3978232
#> 66            0.005632564   0.3978232
#> 67            0.005632564   0.3978232
#> 68            0.005632564   0.3978232
#> 69            0.005632564   0.3978232
#> 70            0.005632564   0.3978232
#> 71            0.005632564   0.3978232
#> 72            0.005632564   0.3978232
#> 73            0.005632564   0.3978232
#> 74            0.005632564   0.3978232
#> 75            0.005632564   0.3978232

# Local per-construct scoring: each latent is scored from its own
# single-factor model (the canonical 2S-PA stage 1) and the results are
# merged into the usual multi-factor layout
get_fs(PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4",
                            "y5", "y6", "y7", "y8")],
       model = " ind60 =~ x1 + x2 + x3
                 dem60 =~ y1 + y2 + y3 + y4
                 dem65 =~ y5 + y6 + y7 + y8 ",
       local = TRUE)
#>       fs_ind60    fs_dem60    fs_dem65 fs_ind60_se fs_dem60_se fs_dem65_se
#> 1  -0.52616832 -2.74872236 -1.37171890   0.1213615   0.6756472   0.5724405
#> 2   0.14365274 -3.03608028 -0.95085100   0.1213615   0.6756472   0.5724405
#> 3   0.71435592  2.67185886  2.73801194   0.1213615   0.6756472   0.5724405
#> 4   1.23992565  2.99369974  1.78509094   0.1213615   0.6756472   0.5724405
#> 5   0.83190803  1.92429320  1.54470402   0.1213615   0.6756472   0.5724405
#> 6   0.21238453  0.99227979 -1.05084086   0.1213615   0.6756472   0.5724405
#> 7   0.11880855  0.99227979 -0.53096192   0.1213615   0.6756472   0.5724405
#> 8   0.11322703 -0.21586296  1.39783099   0.1213615   0.6756472   0.5724405
#> 9   0.25617279 -1.44378438 -0.39762599   0.1213615   0.6756472   0.5724405
#> 10  0.37112496  2.92254844  3.03518364   0.1213615   0.6756472   0.5724405
#> 11  0.67281395  1.34957735  1.69157706   0.1213615   0.6756472   0.5724405
#> 12  0.56885577  0.99227979  1.24608800   0.1213615   0.6756472   0.5724405
#> 13  1.31369791  1.34957735  2.00369398   0.1213615   0.6756472   0.5724405
#> 14  0.22042629  1.84912937 -0.70338102   0.1213615   0.6756472   0.5724405
#> 15  0.57849228  2.31552345  2.88831061   0.1213615   0.6756472   0.5724405
#> 16  0.37805983  2.40746412  0.80499951   0.1213615   0.6756472   0.5724405
#> 17  0.05734046 -0.15715192 -1.19771390   0.1213615   0.6756472   0.5724405
#> 18 -0.01609202 -2.10674549 -0.98511812   0.1213615   0.6756472   0.5724405
#> 19  0.88923616  3.60483453  3.47431108   0.1213615   0.6756472   0.5724405
#> 20  1.11445897  0.63497829 -0.05469352   0.1213615   0.6756472   0.5724405
#> 21  0.94657339  3.60483453  3.62118412   0.1213615   0.6756472   0.5724405
#> 22  0.90122770 -3.39337795 -2.62734059   0.1213615   0.6756472   0.5724405
#> 23  0.58409450 -1.81938756 -1.27737266   0.1213615   0.6756472   0.5724405
#> 24  0.64089192  2.65451561  2.28230961   0.1213615   0.6756472   0.5724405
#> 25  0.91021968  1.97337239  2.21580445   0.1213615   0.6756472   0.5724405
#> 26 -0.89660969 -0.77027212 -1.38899894   0.1213615   0.6756472   0.5724405
#> 27 -0.13195991 -1.27818082 -0.39126438   0.1213615   0.6756472   0.5724405
#> 28 -0.52968769 -0.99453217 -1.91311530   0.1213615   0.6756472   0.5724405
#> 29 -0.81799629  0.42472744 -0.03742971   0.1213615   0.6756472   0.5724405
#> 30 -1.27199371 -0.45058999 -2.42453443   0.1213615   0.6756472   0.5724405
#> 31 -0.32096024 -1.72743303 -0.98362660   0.1213615   0.6756472   0.5724405
#> 32 -1.16780103 -3.22096319 -2.92108665   0.1213615   0.6756472   0.5724405
#> 33 -0.12295473 -1.24467171 -0.98362660   0.1213615   0.6756472   0.5724405
#> 34 -0.04285945 -0.53007648 -2.27361739   0.1213615   0.6756472   0.5724405
#> 35 -0.34323505  0.75106691  0.37463233   0.1213615   0.6756472   0.5724405
#> 36 -0.60541633  1.82492000 -0.19951669   0.1213615   0.6756472   0.5724405
#> 37  0.17688718 -0.28513805 -1.00185417   0.1213615   0.6756472   0.5724405
#> 38 -0.55066055  0.38898019 -1.59180427   0.1213615   0.6756472   0.5724405
#> 39 -1.05988219 -0.56048618 -0.75214950   0.1213615   0.6756472   0.5724405
#> 40 -0.04138802  0.03253769 -1.06763415   0.1213615   0.6756472   0.5724405
#> 41 -0.12611837 -0.22523213 -1.34262601   0.1213615   0.6756472   0.5724405
#> 42 -0.60322892  1.39285768  0.84762975   0.1213615   0.6756472   0.5724405
#> 43 -0.11057176  0.92901951  0.27384517   0.1213615   0.6756472   0.5724405
#> 44 -1.06423085 -1.84123563  0.31144974   0.1213615   0.6756472   0.5724405
#> 45 -1.08354999 -1.82951509 -1.61654871   0.1213615   0.6756472   0.5724405
#> 46 -0.84009484  1.87568627  1.66585009   0.1213615   0.6756472   0.5724405
#> 47 -1.14678213 -1.46208990 -2.62734059   0.1213615   0.6756472   0.5724405
#> 48 -0.57578976  2.39046791  2.06987729   0.1213615   0.6756472   0.5724405
#> 49  0.07186692 -0.86170657 -0.99922629   0.1213615   0.6756472   0.5724405
#> 50  0.14682421  1.81245195  0.80546890   0.1213615   0.6756472   0.5724405
#> 51  0.35871830 -1.46208990 -0.53813741   0.1213615   0.6756472   0.5724405
#> 52 -0.43403195 -2.26050229 -1.18649882   0.1213615   0.6756472   0.5724405
#> 53  0.44603111 -1.17574904  0.50699008   0.1213615   0.6756472   0.5724405
#> 54  0.26352000  0.02760093  1.69204606   0.1213615   0.6756472   0.5724405
#> 55  0.55051165 -3.03608028 -2.62734059   0.1213615   0.6756472   0.5724405
#> 56  0.23453122 -2.68411252 -1.27737266   0.1213615   0.6756472   0.5724405
#> 57  0.27968138 -2.67878261 -1.92858516   0.1213615   0.6756472   0.5724405
#> 58  0.70960640  1.84912937  1.54470402   0.1213615   0.6756472   0.5724405
#> 59  0.25227978 -1.15540938 -1.03962578   0.1213615   0.6756472   0.5724405
#> 60  1.18849297  3.03011868  3.47431108   0.1213615   0.6756472   0.5724405
#> 61  0.21104946 -3.39337795 -2.77421362   0.1213615   0.6756472   0.5724405
#> 62 -1.16516281 -3.39337795 -2.92108665   0.1213615   0.6756472   0.5724405
#> 63 -0.85560065 -3.10602002 -1.23014270   0.1213615   0.6756472   0.5724405
#> 64  0.13398476 -0.44957067 -0.24588734   0.1213615   0.6756472   0.5724405
#> 65 -0.07912189  3.03011868  3.32743805   0.1213615   0.6756472   0.5724405
#> 66 -0.27146711 -1.87744937  0.34647922   0.1213615   0.6756472   0.5724405
#> 67 -0.04417217  3.03011868  2.88681895   0.1213615   0.6756472   0.5724405
#> 68 -1.33425662 -1.23074578 -0.94504866   0.1213615   0.6756472   0.5724405
#> 69 -0.38720750 -1.65235182 -2.33359452   0.1213615   0.6756472   0.5724405
#> 70 -0.55355511 -0.99453217 -0.39275590   0.1213615   0.6756472   0.5724405
#> 71 -0.72242623  0.33138994 -0.30824367   0.1213615   0.6756472   0.5724405
#> 72  0.30607449  1.79828639  1.98383160   0.1213615   0.6756472   0.5724405
#> 73  0.77707950  1.81469038  1.72470940   0.1213615   0.6756472   0.5724405
#> 74  0.06847481  3.22923146  3.18205667   0.1213615   0.6756472   0.5724405
#> 75 -0.11052927 -2.44305891 -2.33508605   0.1213615   0.6756472   0.5724405
#>    ind60_by_fs_ind60 ind60_by_fs_dem60 ind60_by_fs_dem65 dem60_by_fs_ind60
#> 1          0.9657673                 0                 0                 0
#> 2          0.9657673                 0                 0                 0
#> 3          0.9657673                 0                 0                 0
#> 4          0.9657673                 0                 0                 0
#> 5          0.9657673                 0                 0                 0
#> 6          0.9657673                 0                 0                 0
#> 7          0.9657673                 0                 0                 0
#> 8          0.9657673                 0                 0                 0
#> 9          0.9657673                 0                 0                 0
#> 10         0.9657673                 0                 0                 0
#> 11         0.9657673                 0                 0                 0
#> 12         0.9657673                 0                 0                 0
#> 13         0.9657673                 0                 0                 0
#> 14         0.9657673                 0                 0                 0
#> 15         0.9657673                 0                 0                 0
#> 16         0.9657673                 0                 0                 0
#> 17         0.9657673                 0                 0                 0
#> 18         0.9657673                 0                 0                 0
#> 19         0.9657673                 0                 0                 0
#> 20         0.9657673                 0                 0                 0
#> 21         0.9657673                 0                 0                 0
#> 22         0.9657673                 0                 0                 0
#> 23         0.9657673                 0                 0                 0
#> 24         0.9657673                 0                 0                 0
#> 25         0.9657673                 0                 0                 0
#> 26         0.9657673                 0                 0                 0
#> 27         0.9657673                 0                 0                 0
#> 28         0.9657673                 0                 0                 0
#> 29         0.9657673                 0                 0                 0
#> 30         0.9657673                 0                 0                 0
#> 31         0.9657673                 0                 0                 0
#> 32         0.9657673                 0                 0                 0
#> 33         0.9657673                 0                 0                 0
#> 34         0.9657673                 0                 0                 0
#> 35         0.9657673                 0                 0                 0
#> 36         0.9657673                 0                 0                 0
#> 37         0.9657673                 0                 0                 0
#> 38         0.9657673                 0                 0                 0
#> 39         0.9657673                 0                 0                 0
#> 40         0.9657673                 0                 0                 0
#> 41         0.9657673                 0                 0                 0
#> 42         0.9657673                 0                 0                 0
#> 43         0.9657673                 0                 0                 0
#> 44         0.9657673                 0                 0                 0
#> 45         0.9657673                 0                 0                 0
#> 46         0.9657673                 0                 0                 0
#> 47         0.9657673                 0                 0                 0
#> 48         0.9657673                 0                 0                 0
#> 49         0.9657673                 0                 0                 0
#> 50         0.9657673                 0                 0                 0
#> 51         0.9657673                 0                 0                 0
#> 52         0.9657673                 0                 0                 0
#> 53         0.9657673                 0                 0                 0
#> 54         0.9657673                 0                 0                 0
#> 55         0.9657673                 0                 0                 0
#> 56         0.9657673                 0                 0                 0
#> 57         0.9657673                 0                 0                 0
#> 58         0.9657673                 0                 0                 0
#> 59         0.9657673                 0                 0                 0
#> 60         0.9657673                 0                 0                 0
#> 61         0.9657673                 0                 0                 0
#> 62         0.9657673                 0                 0                 0
#> 63         0.9657673                 0                 0                 0
#> 64         0.9657673                 0                 0                 0
#> 65         0.9657673                 0                 0                 0
#> 66         0.9657673                 0                 0                 0
#> 67         0.9657673                 0                 0                 0
#> 68         0.9657673                 0                 0                 0
#> 69         0.9657673                 0                 0                 0
#> 70         0.9657673                 0                 0                 0
#> 71         0.9657673                 0                 0                 0
#> 72         0.9657673                 0                 0                 0
#> 73         0.9657673                 0                 0                 0
#> 74         0.9657673                 0                 0                 0
#> 75         0.9657673                 0                 0                 0
#>    dem60_by_fs_dem60 dem60_by_fs_dem65 dem65_by_fs_ind60 dem65_by_fs_dem60
#> 1          0.8868049                 0                 0                 0
#> 2          0.8868049                 0                 0                 0
#> 3          0.8868049                 0                 0                 0
#> 4          0.8868049                 0                 0                 0
#> 5          0.8868049                 0                 0                 0
#> 6          0.8868049                 0                 0                 0
#> 7          0.8868049                 0                 0                 0
#> 8          0.8868049                 0                 0                 0
#> 9          0.8868049                 0                 0                 0
#> 10         0.8868049                 0                 0                 0
#> 11         0.8868049                 0                 0                 0
#> 12         0.8868049                 0                 0                 0
#> 13         0.8868049                 0                 0                 0
#> 14         0.8868049                 0                 0                 0
#> 15         0.8868049                 0                 0                 0
#> 16         0.8868049                 0                 0                 0
#> 17         0.8868049                 0                 0                 0
#> 18         0.8868049                 0                 0                 0
#> 19         0.8868049                 0                 0                 0
#> 20         0.8868049                 0                 0                 0
#> 21         0.8868049                 0                 0                 0
#> 22         0.8868049                 0                 0                 0
#> 23         0.8868049                 0                 0                 0
#> 24         0.8868049                 0                 0                 0
#> 25         0.8868049                 0                 0                 0
#> 26         0.8868049                 0                 0                 0
#> 27         0.8868049                 0                 0                 0
#> 28         0.8868049                 0                 0                 0
#> 29         0.8868049                 0                 0                 0
#> 30         0.8868049                 0                 0                 0
#> 31         0.8868049                 0                 0                 0
#> 32         0.8868049                 0                 0                 0
#> 33         0.8868049                 0                 0                 0
#> 34         0.8868049                 0                 0                 0
#> 35         0.8868049                 0                 0                 0
#> 36         0.8868049                 0                 0                 0
#> 37         0.8868049                 0                 0                 0
#> 38         0.8868049                 0                 0                 0
#> 39         0.8868049                 0                 0                 0
#> 40         0.8868049                 0                 0                 0
#> 41         0.8868049                 0                 0                 0
#> 42         0.8868049                 0                 0                 0
#> 43         0.8868049                 0                 0                 0
#> 44         0.8868049                 0                 0                 0
#> 45         0.8868049                 0                 0                 0
#> 46         0.8868049                 0                 0                 0
#> 47         0.8868049                 0                 0                 0
#> 48         0.8868049                 0                 0                 0
#> 49         0.8868049                 0                 0                 0
#> 50         0.8868049                 0                 0                 0
#> 51         0.8868049                 0                 0                 0
#> 52         0.8868049                 0                 0                 0
#> 53         0.8868049                 0                 0                 0
#> 54         0.8868049                 0                 0                 0
#> 55         0.8868049                 0                 0                 0
#> 56         0.8868049                 0                 0                 0
#> 57         0.8868049                 0                 0                 0
#> 58         0.8868049                 0                 0                 0
#> 59         0.8868049                 0                 0                 0
#> 60         0.8868049                 0                 0                 0
#> 61         0.8868049                 0                 0                 0
#> 62         0.8868049                 0                 0                 0
#> 63         0.8868049                 0                 0                 0
#> 64         0.8868049                 0                 0                 0
#> 65         0.8868049                 0                 0                 0
#> 66         0.8868049                 0                 0                 0
#> 67         0.8868049                 0                 0                 0
#> 68         0.8868049                 0                 0                 0
#> 69         0.8868049                 0                 0                 0
#> 70         0.8868049                 0                 0                 0
#> 71         0.8868049                 0                 0                 0
#> 72         0.8868049                 0                 0                 0
#> 73         0.8868049                 0                 0                 0
#> 74         0.8868049                 0                 0                 0
#> 75         0.8868049                 0                 0                 0
#>    dem65_by_fs_dem65 ev_fs_ind60 ecov_fs_dem60_fs_ind60 ev_fs_dem60
#> 1          0.8998252  0.01472862                      0   0.4564991
#> 2          0.8998252  0.01472862                      0   0.4564991
#> 3          0.8998252  0.01472862                      0   0.4564991
#> 4          0.8998252  0.01472862                      0   0.4564991
#> 5          0.8998252  0.01472862                      0   0.4564991
#> 6          0.8998252  0.01472862                      0   0.4564991
#> 7          0.8998252  0.01472862                      0   0.4564991
#> 8          0.8998252  0.01472862                      0   0.4564991
#> 9          0.8998252  0.01472862                      0   0.4564991
#> 10         0.8998252  0.01472862                      0   0.4564991
#> 11         0.8998252  0.01472862                      0   0.4564991
#> 12         0.8998252  0.01472862                      0   0.4564991
#> 13         0.8998252  0.01472862                      0   0.4564991
#> 14         0.8998252  0.01472862                      0   0.4564991
#> 15         0.8998252  0.01472862                      0   0.4564991
#> 16         0.8998252  0.01472862                      0   0.4564991
#> 17         0.8998252  0.01472862                      0   0.4564991
#> 18         0.8998252  0.01472862                      0   0.4564991
#> 19         0.8998252  0.01472862                      0   0.4564991
#> 20         0.8998252  0.01472862                      0   0.4564991
#> 21         0.8998252  0.01472862                      0   0.4564991
#> 22         0.8998252  0.01472862                      0   0.4564991
#> 23         0.8998252  0.01472862                      0   0.4564991
#> 24         0.8998252  0.01472862                      0   0.4564991
#> 25         0.8998252  0.01472862                      0   0.4564991
#> 26         0.8998252  0.01472862                      0   0.4564991
#> 27         0.8998252  0.01472862                      0   0.4564991
#> 28         0.8998252  0.01472862                      0   0.4564991
#> 29         0.8998252  0.01472862                      0   0.4564991
#> 30         0.8998252  0.01472862                      0   0.4564991
#> 31         0.8998252  0.01472862                      0   0.4564991
#> 32         0.8998252  0.01472862                      0   0.4564991
#> 33         0.8998252  0.01472862                      0   0.4564991
#> 34         0.8998252  0.01472862                      0   0.4564991
#> 35         0.8998252  0.01472862                      0   0.4564991
#> 36         0.8998252  0.01472862                      0   0.4564991
#> 37         0.8998252  0.01472862                      0   0.4564991
#> 38         0.8998252  0.01472862                      0   0.4564991
#> 39         0.8998252  0.01472862                      0   0.4564991
#> 40         0.8998252  0.01472862                      0   0.4564991
#> 41         0.8998252  0.01472862                      0   0.4564991
#> 42         0.8998252  0.01472862                      0   0.4564991
#> 43         0.8998252  0.01472862                      0   0.4564991
#> 44         0.8998252  0.01472862                      0   0.4564991
#> 45         0.8998252  0.01472862                      0   0.4564991
#> 46         0.8998252  0.01472862                      0   0.4564991
#> 47         0.8998252  0.01472862                      0   0.4564991
#> 48         0.8998252  0.01472862                      0   0.4564991
#> 49         0.8998252  0.01472862                      0   0.4564991
#> 50         0.8998252  0.01472862                      0   0.4564991
#> 51         0.8998252  0.01472862                      0   0.4564991
#> 52         0.8998252  0.01472862                      0   0.4564991
#> 53         0.8998252  0.01472862                      0   0.4564991
#> 54         0.8998252  0.01472862                      0   0.4564991
#> 55         0.8998252  0.01472862                      0   0.4564991
#> 56         0.8998252  0.01472862                      0   0.4564991
#> 57         0.8998252  0.01472862                      0   0.4564991
#> 58         0.8998252  0.01472862                      0   0.4564991
#> 59         0.8998252  0.01472862                      0   0.4564991
#> 60         0.8998252  0.01472862                      0   0.4564991
#> 61         0.8998252  0.01472862                      0   0.4564991
#> 62         0.8998252  0.01472862                      0   0.4564991
#> 63         0.8998252  0.01472862                      0   0.4564991
#> 64         0.8998252  0.01472862                      0   0.4564991
#> 65         0.8998252  0.01472862                      0   0.4564991
#> 66         0.8998252  0.01472862                      0   0.4564991
#> 67         0.8998252  0.01472862                      0   0.4564991
#> 68         0.8998252  0.01472862                      0   0.4564991
#> 69         0.8998252  0.01472862                      0   0.4564991
#> 70         0.8998252  0.01472862                      0   0.4564991
#> 71         0.8998252  0.01472862                      0   0.4564991
#> 72         0.8998252  0.01472862                      0   0.4564991
#> 73         0.8998252  0.01472862                      0   0.4564991
#> 74         0.8998252  0.01472862                      0   0.4564991
#> 75         0.8998252  0.01472862                      0   0.4564991
#>    ecov_fs_dem65_fs_ind60 ecov_fs_dem65_fs_dem60 ev_fs_dem65
#> 1                       0                      0   0.3276882
#> 2                       0                      0   0.3276882
#> 3                       0                      0   0.3276882
#> 4                       0                      0   0.3276882
#> 5                       0                      0   0.3276882
#> 6                       0                      0   0.3276882
#> 7                       0                      0   0.3276882
#> 8                       0                      0   0.3276882
#> 9                       0                      0   0.3276882
#> 10                      0                      0   0.3276882
#> 11                      0                      0   0.3276882
#> 12                      0                      0   0.3276882
#> 13                      0                      0   0.3276882
#> 14                      0                      0   0.3276882
#> 15                      0                      0   0.3276882
#> 16                      0                      0   0.3276882
#> 17                      0                      0   0.3276882
#> 18                      0                      0   0.3276882
#> 19                      0                      0   0.3276882
#> 20                      0                      0   0.3276882
#> 21                      0                      0   0.3276882
#> 22                      0                      0   0.3276882
#> 23                      0                      0   0.3276882
#> 24                      0                      0   0.3276882
#> 25                      0                      0   0.3276882
#> 26                      0                      0   0.3276882
#> 27                      0                      0   0.3276882
#> 28                      0                      0   0.3276882
#> 29                      0                      0   0.3276882
#> 30                      0                      0   0.3276882
#> 31                      0                      0   0.3276882
#> 32                      0                      0   0.3276882
#> 33                      0                      0   0.3276882
#> 34                      0                      0   0.3276882
#> 35                      0                      0   0.3276882
#> 36                      0                      0   0.3276882
#> 37                      0                      0   0.3276882
#> 38                      0                      0   0.3276882
#> 39                      0                      0   0.3276882
#> 40                      0                      0   0.3276882
#> 41                      0                      0   0.3276882
#> 42                      0                      0   0.3276882
#> 43                      0                      0   0.3276882
#> 44                      0                      0   0.3276882
#> 45                      0                      0   0.3276882
#> 46                      0                      0   0.3276882
#> 47                      0                      0   0.3276882
#> 48                      0                      0   0.3276882
#> 49                      0                      0   0.3276882
#> 50                      0                      0   0.3276882
#> 51                      0                      0   0.3276882
#> 52                      0                      0   0.3276882
#> 53                      0                      0   0.3276882
#> 54                      0                      0   0.3276882
#> 55                      0                      0   0.3276882
#> 56                      0                      0   0.3276882
#> 57                      0                      0   0.3276882
#> 58                      0                      0   0.3276882
#> 59                      0                      0   0.3276882
#> 60                      0                      0   0.3276882
#> 61                      0                      0   0.3276882
#> 62                      0                      0   0.3276882
#> 63                      0                      0   0.3276882
#> 64                      0                      0   0.3276882
#> 65                      0                      0   0.3276882
#> 66                      0                      0   0.3276882
#> 67                      0                      0   0.3276882
#> 68                      0                      0   0.3276882
#> 69                      0                      0   0.3276882
#> 70                      0                      0   0.3276882
#> 71                      0                      0   0.3276882
#> 72                      0                      0   0.3276882
#> 73                      0                      0   0.3276882
#> 74                      0                      0   0.3276882
#> 75                      0                      0   0.3276882

# Vector form: one complete single-factor model string per latent, fit
# verbatim (e.g. a within-factor residual covariance the strict string
# grammar rejects)
get_fs(PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")],
       model = c("ind60 =~ x1 + x2 + x3",
                 "dem60 =~ y1 + y2 + y3 + y4
                  y1 ~~ y4"),
       local = TRUE)
#> Warning: lavaan->lav_object_post_check():  
#>    the covariance matrix of the residuals of the observed variables (theta) 
#>    is not positive definite ; use lavInspect(fit, "theta") to investigate.
#>       fs_ind60     fs_dem60 fs_ind60_se fs_dem60_se ind60_by_fs_ind60
#> 1  -0.52616832 -3.482363561   0.1213615          NA         0.9657673
#> 2   0.14365274 -4.226253608   0.1213615          NA         0.9657673
#> 3   0.71435592  3.043357833   0.1213615          NA         0.9657673
#> 4   1.23992565  3.876514686   0.1213615          NA         0.9657673
#> 5   0.83190803  3.681706567   0.1213615          NA         0.9657673
#> 6   0.21238453  2.617461300   0.1213615          NA         0.9657673
#> 7   0.11880855  2.617461300   0.1213615          NA         0.9657673
#> 8   0.11322703 -0.733925249   0.1213615          NA         0.9657673
#> 9   0.25617279 -1.906371531   0.1213615          NA         0.9657673
#> 10  0.37112496  4.606997293   0.1213615          NA         0.9657673
#> 11  0.67281395  2.193926473   0.1213615          NA         0.9657673
#> 12  0.56885577  2.617461300   0.1213615          NA         0.9657673
#> 13  1.31369791  2.193926473   0.1213615          NA         0.9657673
#> 14  0.22042629  1.667491994   0.1213615          NA         0.9657673
#> 15  0.57849228  4.221173427   0.1213615          NA         0.9657673
#> 16  0.37805983  2.052992030   0.1213615          NA         0.9657673
#> 17  0.05734046 -0.358098888   0.1213615          NA         0.9657673
#> 18 -0.01609202 -2.254446012   0.1213615          NA         0.9657673
#> 19  0.88923616  4.861883741   0.1213615          NA         0.9657673
#> 20  1.11445897  3.041000291   0.1213615          NA         0.9657673
#> 21  0.94657339  4.861883741   0.1213615          NA         0.9657673
#> 22  0.90122770 -3.802718655   0.1213615          NA         0.9657673
#> 23  0.58409450 -1.510555964   0.1213615          NA         0.9657673
#> 24  0.64089192  3.769919095   0.1213615          NA         0.9657673
#> 25  0.91021968  2.420942322   0.1213615          NA         0.9657673
#> 26 -0.89660969  0.369255631   0.1213615          NA         0.9657673
#> 27 -0.13195991 -1.458001288   0.1213615          NA         0.9657673
#> 28 -0.52968769 -2.091861670   0.1213615          NA         0.9657673
#> 29 -0.81799629  1.642330703   0.1213615          NA         0.9657673
#> 30 -1.27199371  0.684402047   0.1213615          NA         0.9657673
#> 31 -0.32096024 -1.272511149   0.1213615          NA         0.9657673
#> 32 -1.16780103 -3.356384626   0.1213615          NA         0.9657673
#> 33 -0.12295473 -0.022775870   0.1213615          NA         0.9657673
#> 34 -0.04285945 -0.869845651   0.1213615          NA         0.9657673
#> 35 -0.34323505  0.052334920   0.1213615          NA         0.9657673
#> 36 -0.60541633  0.688161013   0.1213615          NA         0.9657673
#> 37  0.17688718 -0.912540974   0.1213615          NA         0.9657673
#> 38 -0.55066055  1.163563921   0.1213615          NA         0.9657673
#> 39 -1.05988219 -0.871426961   0.1213615          NA         0.9657673
#> 40 -0.04138802 -2.169358240   0.1213615          NA         0.9657673
#> 41 -0.12611837 -2.619661501   0.1213615          NA         0.9657673
#> 42 -0.60322892  1.636501394   0.1213615          NA         0.9657673
#> 43 -0.11057176  0.007466061   0.1213615          NA         0.9657673
#> 44 -1.06423085 -1.795830233   0.1213615          NA         0.9657673
#> 45 -1.08354999 -1.724444561   0.1213615          NA         0.9657673
#> 46 -0.84009484  1.669122251   0.1213615          NA         0.9657673
#> 47 -1.14678213 -1.934090918   0.1213615          NA         0.9657673
#> 48 -0.57578976  2.273570206   0.1213615          NA         0.9657673
#> 49  0.07186692 -1.294362204   0.1213615          NA         0.9657673
#> 50  0.14682421  2.004363896   0.1213615          NA         0.9657673
#> 51  0.35871830 -1.934090918   0.1213615          NA         0.9657673
#> 52 -0.43403195 -2.933191503   0.1213615          NA         0.9657673
#> 53  0.44603111 -1.069288009   0.1213615          NA         0.9657673
#> 54  0.26352000 -0.745155896   0.1213615          NA         0.9657673
#> 55  0.55051165 -4.226253608   0.1213615          NA         0.9657673
#> 56  0.23453122 -3.963936887   0.1213615          NA         0.9657673
#> 57  0.27968138 -4.649788562   0.1213615          NA         0.9657673
#> 58  0.70960640  1.667491994   0.1213615          NA         0.9657673
#> 59  0.25227978 -1.283394346   0.1213615          NA         0.9657673
#> 60  1.18849297  3.374103646   0.1213615          NA         0.9657673
#> 61  0.21104946 -3.802718655   0.1213615          NA         0.9657673
#> 62 -1.16516281 -3.802718655   0.1213615          NA         0.9657673
#> 63 -0.85560065 -3.058828608   0.1213615          NA         0.9657673
#> 64  0.13398476 -0.182343617   0.1213615          NA         0.9657673
#> 65 -0.07912189  3.374103646   0.1213615          NA         0.9657673
#> 66 -0.27146711 -2.346117556   0.1213615          NA         0.9657673
#> 67 -0.04417217  3.374103646   0.1213615          NA         0.9657673
#> 68 -1.33425662 -3.166360868   0.1213615          NA         0.9657673
#> 69 -0.38720750 -2.612944604   0.1213615          NA         0.9657673
#> 70 -0.55355511 -2.091861670   0.1213615          NA         0.9657673
#> 71 -0.72242623 -1.395712591   0.1213615          NA         0.9657673
#> 72  0.30607449  1.189863436   0.1213615          NA         0.9657673
#> 73  0.77707950  1.868752533   0.1213615          NA         0.9657673
#> 74  0.06847481  5.257699189   0.1213615          NA         0.9657673
#> 75 -0.11052927 -2.710754128   0.1213615          NA         0.9657673
#>    ind60_by_fs_dem60 dem60_by_fs_ind60 dem60_by_fs_dem60 ev_fs_ind60
#> 1                  0                 0          1.144843  0.01472862
#> 2                  0                 0          1.144843  0.01472862
#> 3                  0                 0          1.144843  0.01472862
#> 4                  0                 0          1.144843  0.01472862
#> 5                  0                 0          1.144843  0.01472862
#> 6                  0                 0          1.144843  0.01472862
#> 7                  0                 0          1.144843  0.01472862
#> 8                  0                 0          1.144843  0.01472862
#> 9                  0                 0          1.144843  0.01472862
#> 10                 0                 0          1.144843  0.01472862
#> 11                 0                 0          1.144843  0.01472862
#> 12                 0                 0          1.144843  0.01472862
#> 13                 0                 0          1.144843  0.01472862
#> 14                 0                 0          1.144843  0.01472862
#> 15                 0                 0          1.144843  0.01472862
#> 16                 0                 0          1.144843  0.01472862
#> 17                 0                 0          1.144843  0.01472862
#> 18                 0                 0          1.144843  0.01472862
#> 19                 0                 0          1.144843  0.01472862
#> 20                 0                 0          1.144843  0.01472862
#> 21                 0                 0          1.144843  0.01472862
#> 22                 0                 0          1.144843  0.01472862
#> 23                 0                 0          1.144843  0.01472862
#> 24                 0                 0          1.144843  0.01472862
#> 25                 0                 0          1.144843  0.01472862
#> 26                 0                 0          1.144843  0.01472862
#> 27                 0                 0          1.144843  0.01472862
#> 28                 0                 0          1.144843  0.01472862
#> 29                 0                 0          1.144843  0.01472862
#> 30                 0                 0          1.144843  0.01472862
#> 31                 0                 0          1.144843  0.01472862
#> 32                 0                 0          1.144843  0.01472862
#> 33                 0                 0          1.144843  0.01472862
#> 34                 0                 0          1.144843  0.01472862
#> 35                 0                 0          1.144843  0.01472862
#> 36                 0                 0          1.144843  0.01472862
#> 37                 0                 0          1.144843  0.01472862
#> 38                 0                 0          1.144843  0.01472862
#> 39                 0                 0          1.144843  0.01472862
#> 40                 0                 0          1.144843  0.01472862
#> 41                 0                 0          1.144843  0.01472862
#> 42                 0                 0          1.144843  0.01472862
#> 43                 0                 0          1.144843  0.01472862
#> 44                 0                 0          1.144843  0.01472862
#> 45                 0                 0          1.144843  0.01472862
#> 46                 0                 0          1.144843  0.01472862
#> 47                 0                 0          1.144843  0.01472862
#> 48                 0                 0          1.144843  0.01472862
#> 49                 0                 0          1.144843  0.01472862
#> 50                 0                 0          1.144843  0.01472862
#> 51                 0                 0          1.144843  0.01472862
#> 52                 0                 0          1.144843  0.01472862
#> 53                 0                 0          1.144843  0.01472862
#> 54                 0                 0          1.144843  0.01472862
#> 55                 0                 0          1.144843  0.01472862
#> 56                 0                 0          1.144843  0.01472862
#> 57                 0                 0          1.144843  0.01472862
#> 58                 0                 0          1.144843  0.01472862
#> 59                 0                 0          1.144843  0.01472862
#> 60                 0                 0          1.144843  0.01472862
#> 61                 0                 0          1.144843  0.01472862
#> 62                 0                 0          1.144843  0.01472862
#> 63                 0                 0          1.144843  0.01472862
#> 64                 0                 0          1.144843  0.01472862
#> 65                 0                 0          1.144843  0.01472862
#> 66                 0                 0          1.144843  0.01472862
#> 67                 0                 0          1.144843  0.01472862
#> 68                 0                 0          1.144843  0.01472862
#> 69                 0                 0          1.144843  0.01472862
#> 70                 0                 0          1.144843  0.01472862
#> 71                 0                 0          1.144843  0.01472862
#> 72                 0                 0          1.144843  0.01472862
#> 73                 0                 0          1.144843  0.01472862
#> 74                 0                 0          1.144843  0.01472862
#> 75                 0                 0          1.144843  0.01472862
#>    ecov_fs_dem60_fs_ind60 ev_fs_dem60
#> 1                       0   -1.016191
#> 2                       0   -1.016191
#> 3                       0   -1.016191
#> 4                       0   -1.016191
#> 5                       0   -1.016191
#> 6                       0   -1.016191
#> 7                       0   -1.016191
#> 8                       0   -1.016191
#> 9                       0   -1.016191
#> 10                      0   -1.016191
#> 11                      0   -1.016191
#> 12                      0   -1.016191
#> 13                      0   -1.016191
#> 14                      0   -1.016191
#> 15                      0   -1.016191
#> 16                      0   -1.016191
#> 17                      0   -1.016191
#> 18                      0   -1.016191
#> 19                      0   -1.016191
#> 20                      0   -1.016191
#> 21                      0   -1.016191
#> 22                      0   -1.016191
#> 23                      0   -1.016191
#> 24                      0   -1.016191
#> 25                      0   -1.016191
#> 26                      0   -1.016191
#> 27                      0   -1.016191
#> 28                      0   -1.016191
#> 29                      0   -1.016191
#> 30                      0   -1.016191
#> 31                      0   -1.016191
#> 32                      0   -1.016191
#> 33                      0   -1.016191
#> 34                      0   -1.016191
#> 35                      0   -1.016191
#> 36                      0   -1.016191
#> 37                      0   -1.016191
#> 38                      0   -1.016191
#> 39                      0   -1.016191
#> 40                      0   -1.016191
#> 41                      0   -1.016191
#> 42                      0   -1.016191
#> 43                      0   -1.016191
#> 44                      0   -1.016191
#> 45                      0   -1.016191
#> 46                      0   -1.016191
#> 47                      0   -1.016191
#> 48                      0   -1.016191
#> 49                      0   -1.016191
#> 50                      0   -1.016191
#> 51                      0   -1.016191
#> 52                      0   -1.016191
#> 53                      0   -1.016191
#> 54                      0   -1.016191
#> 55                      0   -1.016191
#> 56                      0   -1.016191
#> 57                      0   -1.016191
#> 58                      0   -1.016191
#> 59                      0   -1.016191
#> 60                      0   -1.016191
#> 61                      0   -1.016191
#> 62                      0   -1.016191
#> 63                      0   -1.016191
#> 64                      0   -1.016191
#> 65                      0   -1.016191
#> 66                      0   -1.016191
#> 67                      0   -1.016191
#> 68                      0   -1.016191
#> 69                      0   -1.016191
#> 70                      0   -1.016191
#> 71                      0   -1.016191
#> 72                      0   -1.016191
#> 73                      0   -1.016191
#> 74                      0   -1.016191
#> 75                      0   -1.016191

# Multiple-group
hs_model <- ' visual  =~ x1 + x2 + x3 '
fit <- cfa(hs_model,
           data = HolzingerSwineford1939,
           group = "school")
get_fs(HolzingerSwineford1939, hs_model, group = "school")
#>        fs_visual fs_visual_se visual_by_fs_visual ev_fs_visual      school
#> 1   -0.821165191    0.3391326           0.6734826   0.11501089     Pasteur
#> 2   -0.124009418    0.3391326           0.6734826   0.11501089     Pasteur
#> 3   -0.370072089    0.3391326           0.6734826   0.11501089     Pasteur
#> 4    0.440928618    0.3391326           0.6734826   0.11501089     Pasteur
#> 5   -0.691389016    0.3391326           0.6734826   0.11501089     Pasteur
#> 6   -0.110032619    0.3391326           0.6734826   0.11501089     Pasteur
#> 7   -0.904127845    0.3391326           0.6734826   0.11501089     Pasteur
#> 8   -0.031747573    0.3391326           0.6734826   0.11501089     Pasteur
#> 9   -0.439478981    0.3391326           0.6734826   0.11501089     Pasteur
#> 10  -0.938939050    0.3391326           0.6734826   0.11501089     Pasteur
#> 11  -0.436821880    0.3391326           0.6734826   0.11501089     Pasteur
#> 12   0.305033497    0.3391326           0.6734826   0.11501089     Pasteur
#> 13   0.522076263    0.3391326           0.6734826   0.11501089     Pasteur
#> 14  -0.090367931    0.3391326           0.6734826   0.11501089     Pasteur
#> 15   0.526276771    0.3391326           0.6734826   0.11501089     Pasteur
#> 16  -0.226580678    0.3391326           0.6734826   0.11501089     Pasteur
#> 17  -0.582016192    0.3391326           0.6734826   0.11501089     Pasteur
#> 18   0.017040431    0.3391326           0.6734826   0.11501089     Pasteur
#> 19   0.563052459    0.3391326           0.6734826   0.11501089     Pasteur
#> 20   0.746621910    0.3391326           0.6734826   0.11501089     Pasteur
#> 21   0.234672405    0.3391326           0.6734826   0.11501089     Pasteur
#> 22   1.157487518    0.3391326           0.6734826   0.11501089     Pasteur
#> 23  -0.162272449    0.3391326           0.6734826   0.11501089     Pasteur
#> 24  -0.556027059    0.3391326           0.6734826   0.11501089     Pasteur
#> 25  -0.321443540    0.3391326           0.6734826   0.11501089     Pasteur
#> 26   0.153141050    0.3391326           0.6734826   0.11501089     Pasteur
#> 27   0.696234416    0.3391326           0.6734826   0.11501089     Pasteur
#> 28  -0.020961039    0.3391326           0.6734826   0.11501089     Pasteur
#> 29   0.532601236    0.3391326           0.6734826   0.11501089     Pasteur
#> 30  -0.727687585    0.3391326           0.6734826   0.11501089     Pasteur
#> 31  -0.676719580    0.3391326           0.6734826   0.11501089     Pasteur
#> 32  -1.120216393    0.3391326           0.6734826   0.11501089     Pasteur
#> 33  -0.313631732    0.3391326           0.6734826   0.11501089     Pasteur
#> 34  -0.187091845    0.3391326           0.6734826   0.11501089     Pasteur
#> 35  -0.887709484    0.3391326           0.6734826   0.11501089     Pasteur
#> 36  -0.760795908    0.3391326           0.6734826   0.11501089     Pasteur
#> 37   0.556943532    0.3391326           0.6734826   0.11501089     Pasteur
#> 38  -0.458666570    0.3391326           0.6734826   0.11501089     Pasteur
#> 39   0.514741536    0.3391326           0.6734826   0.11501089     Pasteur
#> 40   0.373009089    0.3391326           0.6734826   0.11501089     Pasteur
#> 41  -0.528550562    0.3391326           0.6734826   0.11501089     Pasteur
#> 42  -0.865864795    0.3391326           0.6734826   0.11501089     Pasteur
#> 43  -1.182344640    0.3391326           0.6734826   0.11501089     Pasteur
#> 44  -0.435334517    0.3391326           0.6734826   0.11501089     Pasteur
#> 45   0.306520860    0.3391326           0.6734826   0.11501089     Pasteur
#> 46   0.821604565    0.3391326           0.6734826   0.11501089     Pasteur
#> 47   1.213927875    0.3391326           0.6734826   0.11501089     Pasteur
#> 48  -0.851887996    0.3391326           0.6734826   0.11501089     Pasteur
#> 49  -0.085053749    0.3391326           0.6734826   0.11501089     Pasteur
#> 50  -0.508885873    0.3391326           0.6734826   0.11501089     Pasteur
#> 51   0.502467638    0.3391326           0.6734826   0.11501089     Pasteur
#> 52   0.284732253    0.3391326           0.6734826   0.11501089     Pasteur
#> 53   0.202677755    0.3391326           0.6734826   0.11501089     Pasteur
#> 54  -0.335953502    0.3391326           0.6734826   0.11501089     Pasteur
#> 55   0.556410369    0.3391326           0.6734826   0.11501089     Pasteur
#> 56  -0.058746970    0.3391326           0.6734826   0.11501089     Pasteur
#> 57  -0.066932487    0.3391326           0.6734826   0.11501089     Pasteur
#> 58   0.554230368    0.3391326           0.6734826   0.11501089     Pasteur
#> 59  -0.321761185    0.3391326           0.6734826   0.11501089     Pasteur
#> 60  -0.421834819    0.3391326           0.6734826   0.11501089     Pasteur
#> 61   0.345476529    0.3391326           0.6734826   0.11501089     Pasteur
#> 62   0.194809883    0.3391326           0.6734826   0.11501089     Pasteur
#> 63  -0.207870208    0.3391326           0.6734826   0.11501089     Pasteur
#> 64  -0.441658981    0.3391326           0.6734826   0.11501089     Pasteur
#> 65   0.102070958    0.3391326           0.6734826   0.11501089     Pasteur
#> 66   0.311198487    0.3391326           0.6734826   0.11501089     Pasteur
#> 67   0.676364229    0.3391326           0.6734826   0.11501089     Pasteur
#> 68   0.297858262    0.3391326           0.6734826   0.11501089     Pasteur
#> 69  -1.055487128    0.3391326           0.6734826   0.11501089     Pasteur
#> 70  -0.737997019    0.3391326           0.6734826   0.11501089     Pasteur
#> 71  -1.576099236    0.3391326           0.6734826   0.11501089     Pasteur
#> 72   0.534360181    0.3391326           0.6734826   0.11501089     Pasteur
#> 73  -0.105888156    0.3391326           0.6734826   0.11501089     Pasteur
#> 74   0.266237302    0.3391326           0.6734826   0.11501089     Pasteur
#> 75  -0.352427927    0.3391326           0.6734826   0.11501089     Pasteur
#> 76  -0.334783784    0.3391326           0.6734826   0.11501089     Pasteur
#> 77   0.133588508    0.3391326           0.6734826   0.11501089     Pasteur
#> 78  -1.035662965    0.3391326           0.6734826   0.11501089     Pasteur
#> 79   0.762507108    0.3391326           0.6734826   0.11501089     Pasteur
#> 80  -0.260699265    0.3391326           0.6734826   0.11501089     Pasteur
#> 81  -0.329095893    0.3391326           0.6734826   0.11501089     Pasteur
#> 82   0.752413211    0.3391326           0.6734826   0.11501089     Pasteur
#> 83   0.149268188    0.3391326           0.6734826   0.11501089     Pasteur
#> 84  -0.208880471    0.3391326           0.6734826   0.11501089     Pasteur
#> 85  -1.078285998    0.3391326           0.6734826   0.11501089     Pasteur
#> 86   0.306043760    0.3391326           0.6734826   0.11501089     Pasteur
#> 87   0.349677056    0.3391326           0.6734826   0.11501089     Pasteur
#> 88   0.165686549    0.3391326           0.6734826   0.11501089     Pasteur
#> 89   0.077307606    0.3391326           0.6734826   0.11501089     Pasteur
#> 90  -0.077401396    0.3391326           0.6734826   0.11501089     Pasteur
#> 91  -0.081863485    0.3391326           0.6734826   0.11501089     Pasteur
#> 92   0.106748566    0.3391326           0.6734826   0.11501089     Pasteur
#> 93  -0.211593616    0.3391326           0.6734826   0.11501089     Pasteur
#> 94  -0.926665153    0.3391326           0.6734826   0.11501089     Pasteur
#> 95  -0.739484382    0.3391326           0.6734826   0.11501089     Pasteur
#> 96   0.570387167    0.3391326           0.6734826   0.11501089     Pasteur
#> 97  -0.913642554    0.3391326           0.6734826   0.11501089     Pasteur
#> 98   0.547484887    0.3391326           0.6734826   0.11501089     Pasteur
#> 99  -0.602850599    0.3391326           0.6734826   0.11501089     Pasteur
#> 100  0.225794270    0.3391326           0.6734826   0.11501089     Pasteur
#> 101  0.620447015    0.3391326           0.6734826   0.11501089     Pasteur
#> 102  0.158885005    0.3391326           0.6734826   0.11501089     Pasteur
#> 103 -0.127938344    0.3391326           0.6734826   0.11501089     Pasteur
#> 104 -0.420347455    0.3391326           0.6734826   0.11501089     Pasteur
#> 105  1.327978307    0.3391326           0.6734826   0.11501089     Pasteur
#> 106  0.181843348    0.3391326           0.6734826   0.11501089     Pasteur
#> 107 -0.148932224    0.3391326           0.6734826   0.11501089     Pasteur
#> 108  0.612373626    0.3391326           0.6734826   0.11501089     Pasteur
#> 109 -0.066558798    0.3391326           0.6734826   0.11501089     Pasteur
#> 110 -0.420880619    0.3391326           0.6734826   0.11501089     Pasteur
#> 111  1.127036295    0.3391326           0.6734826   0.11501089     Pasteur
#> 112  0.237591068    0.3391326           0.6734826   0.11501089     Pasteur
#> 113  0.853758689    0.3391326           0.6734826   0.11501089     Pasteur
#> 114 -0.143618023    0.3391326           0.6734826   0.11501089     Pasteur
#> 115  0.475206679    0.3391326           0.6734826   0.11501089     Pasteur
#> 116 -0.670554590    0.3391326           0.6734826   0.11501089     Pasteur
#> 117  0.022672257    0.3391326           0.6734826   0.11501089     Pasteur
#> 118  0.302002707    0.3391326           0.6734826   0.11501089     Pasteur
#> 119  0.151392125    0.3391326           0.6734826   0.11501089     Pasteur
#> 120 -0.475300449    0.3391326           0.6734826   0.11501089     Pasteur
#> 121 -0.346740056    0.3391326           0.6734826   0.11501089     Pasteur
#> 122 -0.078888759    0.3391326           0.6734826   0.11501089     Pasteur
#> 123  1.197237913    0.3391326           0.6734826   0.11501089     Pasteur
#> 124  0.539243306    0.3391326           0.6734826   0.11501089     Pasteur
#> 125  0.867258388    0.3391326           0.6734826   0.11501089     Pasteur
#> 126  0.592287901    0.3391326           0.6734826   0.11501089     Pasteur
#> 127 -0.500540901    0.3391326           0.6734826   0.11501089     Pasteur
#> 128 -0.361193954    0.3391326           0.6734826   0.11501089     Pasteur
#> 129  0.626883588    0.3391326           0.6734826   0.11501089     Pasteur
#> 130 -0.437514518    0.3391326           0.6734826   0.11501089     Pasteur
#> 131  0.695972854    0.3391326           0.6734826   0.11501089     Pasteur
#> 132  0.424715775    0.3391326           0.6734826   0.11501089     Pasteur
#> 133 -0.203725744    0.3391326           0.6734826   0.11501089     Pasteur
#> 134 -0.441499507    0.3391326           0.6734826   0.11501089     Pasteur
#> 135  0.735619838    0.3391326           0.6734826   0.11501089     Pasteur
#> 136  0.783874697    0.3391326           0.6734826   0.11501089     Pasteur
#> 137  0.565709540    0.3391326           0.6734826   0.11501089     Pasteur
#> 138  0.258425494    0.3391326           0.6734826   0.11501089     Pasteur
#> 139  0.861093397    0.3391326           0.6734826   0.11501089     Pasteur
#> 140 -0.059757233    0.3391326           0.6734826   0.11501089     Pasteur
#> 141 -0.920340689    0.3391326           0.6734826   0.11501089     Pasteur
#> 142  0.845629236    0.3391326           0.6734826   0.11501089     Pasteur
#> 143  1.227427574    0.3391326           0.6734826   0.11501089     Pasteur
#> 144  1.054223601    0.3391326           0.6734826   0.11501089     Pasteur
#> 145 -1.246596805    0.3391326           0.6734826   0.11501089     Pasteur
#> 146 -0.473120468    0.3391326           0.6734826   0.11501089     Pasteur
#> 147 -0.560171503    0.3391326           0.6734826   0.11501089     Pasteur
#> 148 -0.365394462    0.3391326           0.6734826   0.11501089     Pasteur
#> 149  0.084744422    0.3391326           0.6734826   0.11501089     Pasteur
#> 150  0.910676146    0.3391326           0.6734826   0.11501089     Pasteur
#> 151  1.094189533    0.3391326           0.6734826   0.11501089     Pasteur
#> 152 -0.013149231    0.3391326           0.6734826   0.11501089     Pasteur
#> 153 -0.166472976    0.3391326           0.6734826   0.11501089     Pasteur
#> 154  0.008695459    0.3391326           0.6734826   0.11501089     Pasteur
#> 155 -0.094989494    0.3391326           0.6734826   0.11501089     Pasteur
#> 156 -0.457123143    0.3391326           0.6734826   0.11501089     Pasteur
#> 157 -0.915287109    0.3118280           0.6990509   0.09723667 Grant-White
#> 158  0.035963597    0.3118280           0.6990509   0.09723667 Grant-White
#> 159  0.355636604    0.3118280           0.6990509   0.09723667 Grant-White
#> 160 -0.387353871    0.3118280           0.6990509   0.09723667 Grant-White
#> 161 -0.622393942    0.3118280           0.6990509   0.09723667 Grant-White
#> 162  0.195944561    0.3118280           0.6990509   0.09723667 Grant-White
#> 163  1.353023831    0.3118280           0.6990509   0.09723667 Grant-White
#> 164 -0.341506254    0.3118280           0.6990509   0.09723667 Grant-White
#> 165 -0.199493575    0.3118280           0.6990509   0.09723667 Grant-White
#> 166 -0.689869149    0.3118280           0.6990509   0.09723667 Grant-White
#> 167 -0.463929554    0.3118280           0.6990509   0.09723667 Grant-White
#> 168 -0.423001505    0.3118280           0.6990509   0.09723667 Grant-White
#> 169  0.279743296    0.3118280           0.6990509   0.09723667 Grant-White
#> 170 -0.916908219    0.3118280           0.6990509   0.09723667 Grant-White
#> 171  0.589344501    0.3118280           0.6990509   0.09723667 Grant-White
#> 172  0.191474701    0.3118280           0.6990509   0.09723667 Grant-White
#> 173  0.935275715    0.3118280           0.6990509   0.09723667 Grant-White
#> 174  0.393715904    0.3118280           0.6990509   0.09723667 Grant-White
#> 175  0.086569994    0.3118280           0.6990509   0.09723667 Grant-White
#> 176  0.555606898    0.3118280           0.6990509   0.09723667 Grant-White
#> 177 -0.558217193    0.3118280           0.6990509   0.09723667 Grant-White
#> 178  0.766715894    0.3118280           0.6990509   0.09723667 Grant-White
#> 179  0.115548801    0.3118280           0.6990509   0.09723667 Grant-White
#> 180  0.901249191    0.3118280           0.6990509   0.09723667 Grant-White
#> 181  0.174316971    0.3118280           0.6990509   0.09723667 Grant-White
#> 182 -0.078980322    0.3118280           0.6990509   0.09723667 Grant-White
#> 183 -0.581882977    0.3118280           0.6990509   0.09723667 Grant-White
#> 184 -0.661179262    0.3118280           0.6990509   0.09723667 Grant-White
#> 185 -0.245341176    0.3118280           0.6990509   0.09723667 Grant-White
#> 186 -0.195801662    0.3118280           0.6990509   0.09723667 Grant-White
#> 187 -0.281221524    0.3118280           0.6990509   0.09723667 Grant-White
#> 188 -0.293909378    0.3118280           0.6990509   0.09723667 Grant-White
#> 189 -0.604192958    0.3118280           0.6990509   0.09723667 Grant-White
#> 190 -0.738437335    0.3118280           0.6990509   0.09723667 Grant-White
#> 191 -0.304109345    0.3118280           0.6990509   0.09723667 Grant-White
#> 192  0.104931733    0.3118280           0.6990509   0.09723667 Grant-White
#> 193 -0.025781487    0.3118280           0.6990509   0.09723667 Grant-White
#> 194 -0.897318824    0.3118280           0.6990509   0.09723667 Grant-White
#> 195 -0.892560027    0.3118280           0.6990509   0.09723667 Grant-White
#> 196 -0.078402465    0.3118280           0.6990509   0.09723667 Grant-White
#> 197 -0.379063934    0.3118280           0.6990509   0.09723667 Grant-White
#> 198 -0.324926380    0.3118280           0.6990509   0.09723667 Grant-White
#> 199 -0.684299797    0.3118280           0.6990509   0.09723667 Grant-White
#> 200 -0.304109345    0.3118280           0.6990509   0.09723667 Grant-White
#> 201  0.622793169    0.3118280           0.6990509   0.09723667 Grant-White
#> 202 -0.152835419    0.3118280           0.6990509   0.09723667 Grant-White
#> 203 -0.421902013    0.3118280           0.6990509   0.09723667 Grant-White
#> 204 -0.060883872    0.3118280           0.6990509   0.09723667 Grant-White
#> 205 -0.303298790    0.3118280           0.6990509   0.09723667 Grant-White
#> 206  0.425021826    0.3118280           0.6990509   0.09723667 Grant-White
#> 207  0.131478875    0.3118280           0.6990509   0.09723667 Grant-White
#> 208 -0.914998172    0.3118280           0.6990509   0.09723667 Grant-White
#> 209  0.324226132    0.3118280           0.6990509   0.09723667 Grant-White
#> 210 -0.086170767    0.3118280           0.6990509   0.09723667 Grant-White
#> 211  0.428424818    0.3118280           0.6990509   0.09723667 Grant-White
#> 212  0.188465179    0.3118280           0.6990509   0.09723667 Grant-White
#> 213 -0.306958111    0.3118280           0.6990509   0.09723667 Grant-White
#> 214  0.581736955    0.3118280           0.6990509   0.09723667 Grant-White
#> 215 -0.743485051    0.3118280           0.6990509   0.09723667 Grant-White
#> 216  0.263580524    0.3118280           0.6990509   0.09723667 Grant-White
#> 217  0.178786831    0.3118280           0.6990509   0.09723667 Grant-White
#> 218 -0.064832113    0.3118280           0.6990509   0.09723667 Grant-White
#> 219  0.499976414    0.3118280           0.6990509   0.09723667 Grant-White
#> 220 -0.092839593    0.3118280           0.6990509   0.09723667 Grant-White
#> 221 -0.263702931    0.3118280           0.6990509   0.09723667 Grant-White
#> 222 -0.983966325    0.3118280           0.6990509   0.09723667 Grant-White
#> 223  1.434912536    0.3118280           0.6990509   0.09723667 Grant-White
#> 224 -1.037582228    0.3118280           0.6990509   0.09723667 Grant-White
#> 225 -0.047569849    0.3118280           0.6990509   0.09723667 Grant-White
#> 226  1.084767775    0.3118280           0.6990509   0.09723667 Grant-White
#> 227  0.092011181    0.3118280           0.6990509   0.09723667 Grant-White
#> 228 -0.562687052    0.3118280           0.6990509   0.09723667 Grant-White
#> 229 -0.304919900    0.3118280           0.6990509   0.09723667 Grant-White
#> 230  1.038920175    0.3118280           0.6990509   0.09723667 Grant-White
#> 231 -0.789854287    0.3118280           0.6990509   0.09723667 Grant-White
#> 232 -0.602282928    0.3118280           0.6990509   0.09723667 Grant-White
#> 233 -0.894630830    0.3118280           0.6990509   0.09723667 Grant-White
#> 234  0.532614527    0.3118280           0.6990509   0.09723667 Grant-White
#> 235 -0.548955945    0.3118280           0.6990509   0.09723667 Grant-White
#> 236 -0.221514635    0.3118280           0.6990509   0.09723667 Grant-White
#> 237 -0.095849115    0.3118280           0.6990509   0.09723667 Grant-White
#> 238 -0.122235502    0.3118280           0.6990509   0.09723667 Grant-White
#> 239  0.892276864    0.3118280           0.6990509   0.09723667 Grant-White
#> 240  0.328439663    0.3118280           0.6990509   0.09723667 Grant-White
#> 241  1.217391042    0.3118280           0.6990509   0.09723667 Grant-White
#> 242  0.574513902    0.3118280           0.6990509   0.09723667 Grant-White
#> 243  0.160168762    0.3118280           0.6990509   0.09723667 Grant-White
#> 244  0.654909662    0.3118280           0.6990509   0.09723667 Grant-White
#> 245 -0.509777155    0.3118280           0.6990509   0.09723667 Grant-White
#> 246  1.201493560    0.3118280           0.6990509   0.09723667 Grant-White
#> 247  0.584874625    0.3118280           0.6990509   0.09723667 Grant-White
#> 248  0.075142371    0.3118280           0.6990509   0.09723667 Grant-White
#> 249  0.550976266    0.3118280           0.6990509   0.09723667 Grant-White
#> 250 -0.886308302    0.3118280           0.6990509   0.09723667 Grant-White
#> 251  0.552075757    0.3118280           0.6990509   0.09723667 Grant-White
#> 252  1.415972940    0.3118280           0.6990509   0.09723667 Grant-White
#> 253  0.298650301    0.3118280           0.6990509   0.09723667 Grant-White
#> 254 -0.143028906    0.3118280           0.6990509   0.09723667 Grant-White
#> 255  0.245195154    0.3118280           0.6990509   0.09723667 Grant-White
#> 256  0.247072593    0.3118280           0.6990509   0.09723667 Grant-White
#> 257  0.817322291    0.3118280           0.6990509   0.09723667 Grant-White
#> 258  0.651771976    0.3118280           0.6990509   0.09723667 Grant-White
#> 259  1.338875623    0.3118280           0.6990509   0.09723667 Grant-White
#> 260 -1.160005528    0.3118280           0.6990509   0.09723667 Grant-White
#> 261  0.163306449    0.3118280           0.6990509   0.09723667 Grant-White
#> 262 -0.387353871    0.3118280           0.6990509   0.09723667 Grant-White
#> 263 -0.517128372    0.3118280           0.6990509   0.09723667 Grant-White
#> 264  0.065103160    0.3118280           0.6990509   0.09723667 Grant-White
#> 265 -0.115438510    0.3118280           0.6990509   0.09723667 Grant-White
#> 266  0.094049376    0.3118280           0.6990509   0.09723667 Grant-White
#> 267  0.396725409    0.3118280           0.6990509   0.09723667 Grant-White
#> 268  0.672356312    0.3118280           0.6990509   0.09723667 Grant-White
#> 269  1.165974090    0.3118280           0.6990509   0.09723667 Grant-White
#> 270 -0.483518949    0.3118280           0.6990509   0.09723667 Grant-White
#> 271  0.035024877    0.3118280           0.6990509   0.09723667 Grant-White
#> 272  0.741974248    0.3118280           0.6990509   0.09723667 Grant-White
#> 273 -0.170386603    0.3118280           0.6990509   0.09723667 Grant-White
#> 274 -0.205873481    0.3118280           0.6990509   0.09723667 Grant-White
#> 275  0.714777307    0.3118280           0.6990509   0.09723667 Grant-White
#> 276 -0.620772831    0.3118280           0.6990509   0.09723667 Grant-White
#> 277 -0.313626938    0.3118280           0.6990509   0.09723667 Grant-White
#> 278 -0.157466035    0.3118280           0.6990509   0.09723667 Grant-White
#> 279  0.118269386    0.3118280           0.6990509   0.09723667 Grant-White
#> 280  0.101111673    0.3118280           0.6990509   0.09723667 Grant-White
#> 281 -0.625403463    0.3118280           0.6990509   0.09723667 Grant-White
#> 282  0.486638761    0.3118280           0.6990509   0.09723667 Grant-White
#> 283 -0.178676540    0.3118280           0.6990509   0.09723667 Grant-White
#> 284  0.274013189    0.3118280           0.6990509   0.09723667 Grant-White
#> 285 -0.316347523    0.3118280           0.6990509   0.09723667 Grant-White
#> 286 -0.026752814    0.3118280           0.6990509   0.09723667 Grant-White
#> 287  0.245323318    0.3118280           0.6990509   0.09723667 Grant-White
#> 288 -0.356336853    0.3118280           0.6990509   0.09723667 Grant-White
#> 289 -0.581594057    0.3118280           0.6990509   0.09723667 Grant-White
#> 290  0.263002667    0.3118280           0.6990509   0.09723667 Grant-White
#> 291 -0.864680712    0.3118280           0.6990509   0.09723667 Grant-White
#> 292 -0.377964443    0.3118280           0.6990509   0.09723667 Grant-White
#> 293 -0.112717909    0.3118280           0.6990509   0.09723667 Grant-White
#> 294  0.114449326    0.3118280           0.6990509   0.09723667 Grant-White
#> 295  0.001287274    0.3118280           0.6990509   0.09723667 Grant-White
#> 296  0.597634438    0.3118280           0.6990509   0.09723667 Grant-White
#> 297 -0.252531637    0.3118280           0.6990509   0.09723667 Grant-White
#> 298 -0.472901881    0.3118280           0.6990509   0.09723667 Grant-White
#> 299 -0.187255397    0.3118280           0.6990509   0.09723667 Grant-White
#> 300 -0.542415283    0.3118280           0.6990509   0.09723667 Grant-White
#> 301  0.358774274    0.3118280           0.6990509   0.09723667 Grant-White
# Or without the model
get_fs(HolzingerSwineford1939[c("school", "x4", "x5", "x6")],
       group = "school")
#>             fs_f1  fs_f1_se f1_by_fs_f1   ev_fs_f1      school
#> 1    0.3074500370 0.2999315   0.8833584 0.08995892     Pasteur
#> 2   -0.7746062892 0.2999315   0.8833584 0.08995892     Pasteur
#> 3   -1.5843019574 0.2999315   0.8833584 0.08995892     Pasteur
#> 4    0.2739579120 0.2999315   0.8833584 0.08995892     Pasteur
#> 5    0.1440153923 0.2999315   0.8833584 0.08995892     Pasteur
#> 6   -1.0440895948 0.2999315   0.8833584 0.08995892     Pasteur
#> 7    1.0507357396 0.2999315   0.8833584 0.08995892     Pasteur
#> 8    0.1041882698 0.2999315   0.8833584 0.08995892     Pasteur
#> 9    0.7750146375 0.2999315   0.8833584 0.08995892     Pasteur
#> 10   0.4822117444 0.2999315   0.8833584 0.08995892     Pasteur
#> 11  -0.4511886490 0.2999315   0.8833584 0.08995892     Pasteur
#> 12   0.3522691973 0.2999315   0.8833584 0.08995892     Pasteur
#> 13   0.0657041070 0.2999315   0.8833584 0.08995892     Pasteur
#> 14   0.3259750264 0.2999315   0.8833584 0.08995892     Pasteur
#> 15   1.3008341323 0.2999315   0.8833584 0.08995892     Pasteur
#> 16  -0.2804588573 0.2999315   0.8833584 0.08995892     Pasteur
#> 17  -0.3604017581 0.2999315   0.8833584 0.08995892     Pasteur
#> 18  -0.7502293722 0.2999315   0.8833584 0.08995892     Pasteur
#> 19   1.2600468631 0.2999315   0.8833584 0.08995892     Pasteur
#> 20   0.2874908636 0.2999315   0.8833584 0.08995892     Pasteur
#> 21  -1.1394826729 0.2999315   0.8833584 0.08995892     Pasteur
#> 22  -0.3791151940 0.2999315   0.8833584 0.08995892     Pasteur
#> 23   1.0444979094 0.2999315   0.8833584 0.08995892     Pasteur
#> 24  -0.6248930150 0.2999315   0.8833584 0.08995892     Pasteur
#> 25  -0.3689426673 0.2999315   0.8833584 0.08995892     Pasteur
#> 26  -0.5663524842 0.2999315   0.8833584 0.08995892     Pasteur
#> 27  -0.9568515617 0.2999315   0.8833584 0.08995892     Pasteur
#> 28  -0.8137619455 0.2999315   0.8833584 0.08995892     Pasteur
#> 29  -0.5028199109 0.2999315   0.8833584 0.08995892     Pasteur
#> 30  -0.3054100713 0.2999315   0.8833584 0.08995892     Pasteur
#> 31  -0.3728773637 0.2999315   0.8833584 0.08995892     Pasteur
#> 32  -0.8529175744 0.2999315   0.8833584 0.08995892     Pasteur
#> 33   0.4101382620 0.2999315   0.8833584 0.08995892     Pasteur
#> 34   0.1848026386 0.2999315   0.8833584 0.08995892     Pasteur
#> 35  -0.9680814159 0.2999315   0.8833584 0.08995892     Pasteur
#> 36  -0.1436070439 0.2999315   0.8833584 0.08995892     Pasteur
#> 37  -0.0126071782 0.2999315   0.8833584 0.08995892     Pasteur
#> 38  -0.3531066368 0.2999315   0.8833584 0.08995892     Pasteur
#> 39   1.0665717701 0.2999315   0.8833584 0.08995892     Pasteur
#> 40   0.7993915544 0.2999315   0.8833584 0.08995892     Pasteur
#> 41   0.1110975387 0.2999315   0.8833584 0.08995892     Pasteur
#> 42  -0.2735495637 0.2999315   0.8833584 0.08995892     Pasteur
#> 43  -0.7167372327 0.2999315   0.8833584 0.08995892     Pasteur
#> 44  -0.6594424814 0.2999315   0.8833584 0.08995892     Pasteur
#> 45   0.9507364688 0.2999315   0.8833584 0.08995892     Pasteur
#> 46  -0.0478281107 0.2999315   0.8833584 0.08995892     Pasteur
#> 47  -1.7573348450 0.2999315   0.8833584 0.08995892     Pasteur
#> 48  -0.4620326417 0.2999315   0.8833584 0.08995892     Pasteur
#> 49  -1.0305566597 0.2999315   0.8833584 0.08995892     Pasteur
#> 50   0.7618675383 0.2999315   0.8833584 0.08995892     Pasteur
#> 51  -1.3414986833 0.2999315   0.8833584 0.08995892     Pasteur
#> 52  -0.0761397743 0.2999315   0.8833584 0.08995892     Pasteur
#> 53  -1.8231705136 0.2999315   0.8833584 0.08995892     Pasteur
#> 54   0.0094666323 0.2999315   0.8833584 0.08995892     Pasteur
#> 55   0.0436302463 0.2999315   0.8833584 0.08995892     Pasteur
#> 56  -0.0001315726 0.2999315   0.8833584 0.08995892     Pasteur
#> 57  -0.6571393750 0.2999315   0.8833584 0.08995892     Pasteur
#> 58   0.6121542642 0.2999315   0.8833584 0.08995892     Pasteur
#> 59  -1.0957208458 0.2999315   0.8833584 0.08995892     Pasteur
#> 60  -0.9289257761 0.2999315   0.8833584 0.08995892     Pasteur
#> 61   0.9553426588 0.2999315   0.8833584 0.08995892     Pasteur
#> 62   0.3647448029 0.2999315   0.8833584 0.08995892     Pasteur
#> 63  -1.1003270330 0.2999315   0.8833584 0.08995892     Pasteur
#> 64   0.4259742697 0.2999315   0.8833584 0.08995892     Pasteur
#> 65   0.1689666035 0.2999315   0.8833584 0.08995892     Pasteur
#> 66   0.6586050419 0.2999315   0.8833584 0.08995892     Pasteur
#> 67  -0.2952375720 0.2999315   0.8833584 0.08995892     Pasteur
#> 68  -0.4745082474 0.2999315   0.8833584 0.08995892     Pasteur
#> 69  -0.5613604418 0.2999315   0.8833584 0.08995892     Pasteur
#> 70   0.5632119383 0.2999315   0.8833584 0.08995892     Pasteur
#> 71  -0.7769093955 0.2999315   0.8833584 0.08995892     Pasteur
#> 72  -1.1434174113 0.2999315   0.8833584 0.08995892     Pasteur
#> 73  -0.7808440920 0.2999315   0.8833584 0.08995892     Pasteur
#> 74   0.9270309952 0.2999315   0.8833584 0.08995892     Pasteur
#> 75  -0.5426470306 0.2999315   0.8833584 0.08995892     Pasteur
#> 76  -1.1065648385 0.2999315   0.8833584 0.08995892     Pasteur
#> 77   0.6233841094 0.2999315   0.8833584 0.08995892     Pasteur
#> 78  -1.7010974136 0.2999315   0.8833584 0.08995892     Pasteur
#> 79  -0.0013773330 0.2999315   0.8833584 0.08995892     Pasteur
#> 80  -1.2772946275 0.2999315   0.8833584 0.08995892     Pasteur
#> 81  -1.0344913561 0.2999315   0.8833584 0.08995892     Pasteur
#> 82   0.4345151789 0.2999315   0.8833584 0.08995892     Pasteur
#> 83  -1.5280645151 0.2999315   0.8833584 0.08995892     Pasteur
#> 84  -0.6101143231 0.2999315   0.8833584 0.08995892     Pasteur
#> 85  -1.5122284832 0.2999315   0.8833584 0.08995892     Pasteur
#> 86   1.2011204798 0.2999315   0.8833584 0.08995892     Pasteur
#> 87  -0.3258523027 0.2999315   0.8833584 0.08995892     Pasteur
#> 88   0.6687775411 0.2999315   0.8833584 0.08995892     Pasteur
#> 89  -1.6382363147 0.2999315   0.8833584 0.08995892     Pasteur
#> 90   0.3062042720 0.2999315   0.8833584 0.08995892     Pasteur
#> 91  -0.4745082474 0.2999315   0.8833584 0.08995892     Pasteur
#> 92  -0.6798846937 0.2999315   0.8833584 0.08995892     Pasteur
#> 93  -1.1865077366 0.2999315   0.8833584 0.08995892     Pasteur
#> 94  -0.5011882981 0.2999315   0.8833584 0.08995892     Pasteur
#> 95   0.7362448336 0.2999315   0.8833584 0.08995892     Pasteur
#> 96   0.4934415622 0.2999315   0.8833584 0.08995892     Pasteur
#> 97   0.9661866515 0.2999315   0.8833584 0.08995892     Pasteur
#> 98   1.7212764661 0.2999315   0.8833584 0.08995892     Pasteur
#> 99  -0.8199997483 0.2999315   0.8833584 0.08995892     Pasteur
#> 100  0.7369163271 0.2999315   0.8833584 0.08995892     Pasteur
#> 101 -0.5403439270 0.2999315   0.8833584 0.08995892     Pasteur
#> 102  0.0525570079 0.2999315   0.8833584 0.08995892     Pasteur
#> 103  0.4973762860 0.2999315   0.8833584 0.08995892     Pasteur
#> 104  0.5434412113 0.2999315   0.8833584 0.08995892     Pasteur
#> 105  1.4015048920 0.2999315   0.8833584 0.08995892     Pasteur
#> 106  0.5338429790 0.2999315   0.8833584 0.08995892     Pasteur
#> 107  1.5005470281 0.2999315   0.8833584 0.08995892     Pasteur
#> 108 -0.4353526184 0.2999315   0.8833584 0.08995892     Pasteur
#> 109  1.7269399974 0.2999315   0.8833584 0.08995892     Pasteur
#> 110 -0.1863115671 0.2999315   0.8833584 0.08995892     Pasteur
#> 111  0.7431541299 0.2999315   0.8833584 0.08995892     Pasteur
#> 112  0.3345159128 0.2999315   0.8833584 0.08995892     Pasteur
#> 113  0.3111963144 0.2999315   0.8833584 0.08995892     Pasteur
#> 114  0.6750153713 0.2999315   0.8833584 0.08995892     Pasteur
#> 115 -1.3822859470 0.2999315   0.8833584 0.08995892     Pasteur
#> 116  0.4299090164 0.2999315   0.8833584 0.08995892     Pasteur
#> 117  0.4368182853 0.2999315   0.8833584 0.08995892     Pasteur
#> 118 -0.6334339242 0.2999315   0.8833584 0.08995892     Pasteur
#> 119 -0.9153928519 0.2999315   0.8833584 0.08995892     Pasteur
#> 120 -0.2662544424 0.2999315   0.8833584 0.08995892     Pasteur
#> 121 -0.1238362896 0.2999315   0.8833584 0.08995892     Pasteur
#> 122 -0.1987871727 0.2999315   0.8833584 0.08995892     Pasteur
#> 123  1.5676284956 0.2999315   0.8833584 0.08995892     Pasteur
#> 124 -0.2906313821 0.2999315   0.8833584 0.08995892     Pasteur
#> 125  0.7125393874 0.2999315   0.8833584 0.08995892     Pasteur
#> 126  0.1324998831 0.2999315   0.8833584 0.08995892     Pasteur
#> 127  1.1488177518 0.2999315   0.8833584 0.08995892     Pasteur
#> 128  0.5559168169 0.2999315   0.8833584 0.08995892     Pasteur
#> 129  0.8572606191 0.2999315   0.8833584 0.08995892     Pasteur
#> 130 -0.9789254060 0.2999315   0.8833584 0.08995892     Pasteur
#> 131  1.4416206448 0.2999315   0.8833584 0.08995892     Pasteur
#> 132  0.4542859333 0.2999315   0.8833584 0.08995892     Pasteur
#> 133 -1.3845890506 0.2999315   0.8833584 0.08995892     Pasteur
#> 134 -0.2883282757 0.2999315   0.8833584 0.08995892     Pasteur
#> 135  0.4430560881 0.2999315   0.8833584 0.08995892     Pasteur
#> 136  1.2089899229 0.2999315   0.8833584 0.08995892     Pasteur
#> 137  1.1942112109 0.2999315   0.8833584 0.08995892     Pasteur
#> 138  0.6013102486 0.2999315   0.8833584 0.08995892     Pasteur
#> 139  0.2371053667 0.2999315   0.8833584 0.08995892     Pasteur
#> 140  1.0053422804 0.2999315   0.8833584 0.08995892     Pasteur
#> 141  0.8095640810 0.2999315   0.8833584 0.08995892     Pasteur
#> 142 -1.4408264861 0.2999315   0.8833584 0.08995892     Pasteur
#> 143  1.3821199900 0.2999315   0.8833584 0.08995892     Pasteur
#> 144  2.7284791267 0.2999315   0.8833584 0.08995892     Pasteur
#> 145  0.0123440494 0.2999315   0.8833584 0.08995892     Pasteur
#> 146  0.8277032180 0.2999315   0.8833584 0.08995892     Pasteur
#> 147  0.7862444827 0.2999315   0.8833584 0.08995892     Pasteur
#> 148 -1.1325734026 0.2999315   0.8833584 0.08995892     Pasteur
#> 149  1.7660956264 0.2999315   0.8833584 0.08995892     Pasteur
#> 150 -0.3712457509 0.2999315   0.8833584 0.08995892     Pasteur
#> 151  1.8944065607 0.2999315   0.8833584 0.08995892     Pasteur
#> 152  0.6098511578 0.2999315   0.8833584 0.08995892     Pasteur
#> 153  0.2654170028 0.2999315   0.8833584 0.08995892     Pasteur
#> 154 -1.1434174113 0.2999315   0.8833584 0.08995892     Pasteur
#> 155 -0.6005160981 0.2999315   0.8833584 0.08995892     Pasteur
#> 156  0.3562038937 0.2999315   0.8833584 0.08995892     Pasteur
#> 157 -0.3952560297 0.3152173   0.8801489 0.09936192 Grant-White
#> 158 -0.6339724772 0.3152173   0.8801489 0.09936192 Grant-White
#> 159  0.2006240287 0.3152173   0.8801489 0.09936192 Grant-White
#> 160 -0.3424279836 0.3152173   0.8801489 0.09936192 Grant-White
#> 161  0.4351919667 0.3152173   0.8801489 0.09936192 Grant-White
#> 162  0.3115124621 0.3152173   0.8801489 0.09936192 Grant-White
#> 163  2.1561291304 0.3152173   0.8801489 0.09936192 Grant-White
#> 164 -0.2901419159 0.3152173   0.8801489 0.09936192 Grant-White
#> 165 -0.0836930350 0.3152173   0.8801489 0.09936192 Grant-White
#> 166 -0.0180739626 0.3152173   0.8801489 0.09936192 Grant-White
#> 167 -0.3482023390 0.3152173   0.8801489 0.09936192 Grant-White
#> 168 -2.1021597427 0.3152173   0.8801489 0.09936192 Grant-White
#> 169 -0.6455211880 0.3152173   0.8801489 0.09936192 Grant-White
#> 170 -1.4615522765 0.3152173   0.8801489 0.09936192 Grant-White
#> 171  1.0262088017 0.3152173   0.8801489 0.09936192 Grant-White
#> 172 -1.0506495591 0.3152173   0.8801489 0.09936192 Grant-White
#> 173  0.4308707217 0.3152173   0.8801489 0.09936192 Grant-White
#> 174  0.9582254779 0.3152173   0.8801489 0.09936192 Grant-White
#> 175 -0.2535530327 0.3152173   0.8801489 0.09936192 Grant-White
#> 176  1.4214142722 0.3152173   0.8801489 0.09936192 Grant-White
#> 177 -0.9540052265 0.3152173   0.8801489 0.09936192 Grant-White
#> 178  1.0029529231 0.3152173   0.8801489 0.09936192 Grant-White
#> 179  1.2184135369 0.3152173   0.8801489 0.09936192 Grant-White
#> 180 -1.2498710090 0.3152173   0.8801489 0.09936192 Grant-White
#> 181 -0.5198466305 0.3152173   0.8801489 0.09936192 Grant-White
#> 182 -0.0471041875 0.3152173   0.8801489 0.09936192 Grant-White
#> 183  0.4393404762 0.3152173   0.8801489 0.09936192 Grant-White
#> 184 -1.0311730096 0.3152173   0.8801489 0.09936192 Grant-White
#> 185 -0.9502259332 0.3152173   0.8801489 0.09936192 Grant-White
#> 186 -0.1141763435 0.3152173   0.8801489 0.09936192 Grant-White
#> 187 -0.4004884068 0.3152173   0.8801489 0.09936192 Grant-White
#> 188  0.1050635903 0.3152173   0.8801489 0.09936192 Grant-White
#> 189 -0.3354112770 0.3152173   0.8801489 0.09936192 Grant-White
#> 190 -1.5556596125 0.3152173   0.8801489 0.09936192 Grant-White
#> 191 -0.8442006782 0.3152173   0.8801489 0.09936192 Grant-White
#> 192  0.0780283916 0.3152173   0.8801489 0.09936192 Grant-White
#> 193  0.2011659712 0.3152173   0.8801489 0.09936192 Grant-White
#> 194 -2.5263954621 0.3152173   0.8801489 0.09936192 Grant-White
#> 195 -0.6914909578 0.3152173   0.8801489 0.09936192 Grant-White
#> 196  2.0234378663 0.3152173   0.8801489 0.09936192 Grant-White
#> 197  0.9733807823 0.3152173   0.8801489 0.09936192 Grant-White
#> 198  0.4204058784 0.3152173   0.8801489 0.09936192 Grant-White
#> 199 -0.9260589225 0.3152173   0.8801489 0.09936192 Grant-White
#> 200 -0.5669003212 0.3152173   0.8801489 0.09936192 Grant-White
#> 201  0.1202188947 0.3152173   0.8801489 0.09936192 Grant-White
#> 202  0.6210804306 0.3152173   0.8801489 0.09936192 Grant-White
#> 203 -1.3421940527 0.3152173   0.8801489 0.09936192 Grant-White
#> 204  0.1625820976 0.3152173   0.8801489 0.09936192 Grant-White
#> 205 -0.0323180992 0.3152173   0.8801489 0.09936192 Grant-White
#> 206  0.2444403061 0.3152173   0.8801489 0.09936192 Grant-White
#> 207 -0.7443190039 0.3152173   0.8801489 0.09936192 Grant-White
#> 208 -0.4379884220 0.3152173   0.8801489 0.09936192 Grant-White
#> 209 -1.8529784508 0.3152173   0.8801489 0.09936192 Grant-White
#> 210 -0.8256352699 0.3152173   0.8801489 0.09936192 Grant-White
#> 211  1.2003900978 0.3152173   0.8801489 0.09936192 Grant-White
#> 212 -0.3328743349 0.3152173   0.8801489 0.09936192 Grant-White
#> 213  0.0199679685 0.3152173   0.8801489 0.09936192 Grant-White
#> 214  1.6758280024 0.3152173   0.8801489 0.09936192 Grant-White
#> 215 -0.7190680830 0.3152173   0.8801489 0.09936192 Grant-White
#> 216 -0.2050463207 0.3152173   0.8801489 0.09936192 Grant-White
#> 217  1.9621400656 0.3152173   0.8801489 0.09936192 Grant-White
#> 218 -0.9345287129 0.3152173   0.8801489 0.09936192 Grant-White
#> 219 -0.3534347427 0.3152173   0.8801489 0.09936192 Grant-White
#> 220 -1.9580925486 0.3152173   0.8801489 0.09936192 Grant-White
#> 221 -1.3602175025 0.3152173   0.8801489 0.09936192 Grant-White
#> 222  0.0859562303 0.3152173   0.8801489 0.09936192 Grant-White
#> 223 -0.2340765190 0.3152173   0.8801489 0.09936192 Grant-White
#> 224  0.6780569596 0.3152173   0.8801489 0.09936192 Grant-White
#> 225 -0.4295186317 0.3152173   0.8801489 0.09936192 Grant-White
#> 226 -0.6920329003 0.3152173   0.8801489 0.09936192 Grant-White
#> 227 -0.7158307215 0.3152173   0.8801489 0.09936192 Grant-White
#> 228 -0.1960345878 0.3152173   0.8801489 0.09936192 Grant-White
#> 229 -0.3676788793 0.3152173   0.8801489 0.09936192 Grant-White
#> 230  1.7742566022 0.3152173   0.8801489 0.09936192 Grant-White
#> 231 -0.6792418741 0.3152173   0.8801489 0.09936192 Grant-White
#> 232  0.0760333654 0.3152173   0.8801489 0.09936192 Grant-White
#> 233  1.4989512714 0.3152173   0.8801489 0.09936192 Grant-White
#> 234 -0.7881352546 0.3152173   0.8801489 0.09936192 Grant-White
#> 235 -1.3564381627 0.3152173   0.8801489 0.09936192 Grant-White
#> 236 -0.3424279836 0.3152173   0.8801489 0.09936192 Grant-White
#> 237  1.1080670102 0.3152173   0.8801489 0.09936192 Grant-White
#> 238  1.0476803416 0.3152173   0.8801489 0.09936192 Grant-White
#> 239  1.0224294726 0.3152173   0.8801489 0.09936192 Grant-White
#> 240  0.3823639473 0.3152173   0.8801489 0.09936192 Grant-White
#> 241  0.5644730554 0.3152173   0.8801489 0.09936192 Grant-White
#> 242  0.7880342609 0.3152173   0.8801489 0.09936192 Grant-White
#> 243  0.4308707217 0.3152173   0.8801489 0.09936192 Grant-White
#> 244  0.4038355230 0.3152173   0.8801489 0.09936192 Grant-White
#> 245 -0.4437627683 0.3152173   0.8801489 0.09936192 Grant-White
#> 246  0.0812657691 0.3152173   0.8801489 0.09936192 Grant-White
#> 247  0.1483379252 0.3152173   0.8801489 0.09936192 Grant-White
#> 248  0.5392221863 0.3152173   0.8801489 0.09936192 Grant-White
#> 249  0.5359848088 0.3152173   0.8801489 0.09936192 Grant-White
#> 250 -0.1702417405 0.3152173   0.8801489 0.09936192 Grant-White
#> 251  0.4361030987 0.3152173   0.8801489 0.09936192 Grant-White
#> 252  2.2284336635 0.3152173   0.8801489 0.09936192 Grant-White
#> 253  1.3648069594 0.3152173   0.8801489 0.09936192 Grant-White
#> 254  0.6513529802 0.3152173   0.8801489 0.09936192 Grant-White
#> 255  1.6943933841 0.3152173   0.8801489 0.09936192 Grant-White
#> 256 -0.1574506784 0.3152173   0.8801489 0.09936192 Grant-White
#> 257 -0.2768089380 0.3152173   0.8801489 0.09936192 Grant-White
#> 258  1.8047399107 0.3152173   0.8801489 0.09936192 Grant-White
#> 259 -0.0328600509 0.3152173   0.8801489 0.09936192 Grant-White
#> 260  0.2154100812 0.3152173   0.8801489 0.09936192 Grant-White
#> 261  0.6358664831 0.3152173   0.8801489 0.09936192 Grant-White
#> 262 -0.4412257995 0.3152173   0.8801489 0.09936192 Grant-White
#> 263  0.2438983636 0.3152173   0.8801489 0.09936192 Grant-White
#> 264  0.9739226982 0.3152173   0.8801489 0.09936192 Grant-White
#> 265  1.0061903098 0.3152173   0.8801489 0.09936192 Grant-White
#> 266  0.9596785882 0.3152173   0.8801489 0.09936192 Grant-White
#> 267  2.1052961106 0.3152173   0.8801489 0.09936192 Grant-White
#> 268  0.9501249128 0.3152173   0.8801489 0.09936192 Grant-White
#> 269  1.1403346218 0.3152173   0.8801489 0.09936192 Grant-White
#> 270 -0.4152744951 0.3152173   0.8801489 0.09936192 Grant-White
#> 271  0.5026333389 0.3152173   0.8801489 0.09936192 Grant-White
#> 272  0.0051818802 0.3152173   0.8801489 0.09936192 Grant-White
#> 273 -0.3096184562 0.3152173   0.8801489 0.09936192 Grant-White
#> 274 -0.3857023185 0.3152173   0.8801489 0.09936192 Grant-White
#> 275 -0.8874750488 0.3152173   0.8801489 0.09936192 Grant-White
#> 276 -1.1134004701 0.3152173   0.8801489 0.09936192 Grant-White
#> 277 -0.2283021636 0.3152173   0.8801489 0.09936192 Grant-White
#> 278  0.1678144746 0.3152173   0.8801489 0.09936192 Grant-White
#> 279 -1.1662284537 0.3152173   0.8801489 0.09936192 Grant-White
#> 280 -0.4527744745 0.3152173   0.8801489 0.09936192 Grant-White
#> 281 -0.6952703137 0.3152173   0.8801489 0.09936192 Grant-White
#> 282  1.1655855175 0.3152173   0.8801489 0.09936192 Grant-White
#> 283 -0.4908164323 0.3152173   0.8801489 0.09936192 Grant-White
#> 284  0.4541265645 0.3152173   0.8801489 0.09936192 Grant-White
#> 285 -0.7591050564 0.3152173   0.8801489 0.09936192 Grant-White
#> 286 -0.4623281857 0.3152173   0.8801489 0.09936192 Grant-White
#> 287  1.3363186770 0.3152173   0.8801489 0.09936192 Grant-White
#> 288 -0.7823609350 0.3152173   0.8801489 0.09936192 Grant-White
#> 289  0.1140753232 0.3152173   0.8801489 0.09936192 Grant-White
#> 290 -0.2611117177 0.3152173   0.8801489 0.09936192 Grant-White
#> 291 -0.5849237710 0.3152173   0.8801489 0.09936192 Grant-White
#> 292 -1.4087242662 0.3152173   0.8801489 0.09936192 Grant-White
#> 293 -0.2430882519 0.3152173   0.8801489 0.09936192 Grant-White
#> 294 -0.1760160958 0.3152173   0.8801489 0.09936192 Grant-White
#> 295 -0.7448609198 0.3152173   0.8801489 0.09936192 Grant-White
#> 296 -0.1341948089 0.3152173   0.8801489 0.09936192 Grant-White
#> 297 -0.7480982973 0.3152173   0.8801489 0.09936192 Grant-White
#> 298 -0.9345287129 0.3152173   0.8801489 0.09936192 Grant-White
#> 299  0.8873740285 0.3152173   0.8801489 0.09936192 Grant-White
#> 300 -0.0566578363 0.3152173   0.8801489 0.09936192 Grant-White
#> 301  0.5830384728 0.3152173   0.8801489 0.09936192 Grant-White

# Fixed external latent prior (shared across groups) for regression scores;
# conceptually similar to mirt::fscores(mean, cov)
fit <- cfa("visual =~ x1 + x2 + x3",
           data = HolzingerSwineford1939,
           group = "school", group.equal = c("loadings", "intercepts"))
get_fs(fit, prior_mean = c(visual = -0.12), prior_cov = 0.33)
#>         fs_visual fs_visual_se visual_by_fs_visual ev_fs_visual      school
#> 1   -0.8661444157    0.2460880           0.7578508   0.06055928     Pasteur
#> 2   -0.1702694303    0.2460880           0.7578508   0.06055928     Pasteur
#> 3   -0.3296920230    0.2460880           0.7578508   0.06055928     Pasteur
#> 4    0.2799417046    0.2460880           0.7578508   0.06055928     Pasteur
#> 5   -0.7047116341    0.2460880           0.7578508   0.06055928     Pasteur
#> 6   -0.1351731231    0.2460880           0.7578508   0.06055928     Pasteur
#> 7   -0.7538261380    0.2460880           0.7578508   0.06055928     Pasteur
#> 8   -0.1903321424    0.2460880           0.7578508   0.06055928     Pasteur
#> 9   -0.4470125320    0.2460880           0.7578508   0.06055928     Pasteur
#> 10  -0.8320432414    0.2460880           0.7578508   0.06055928     Pasteur
#> 11  -0.3236677575    0.2460880           0.7578508   0.06055928     Pasteur
#> 12   0.1876927826    0.2460880           0.7578508   0.06055928     Pasteur
#> 13   0.5737488481    0.2460880           0.7578508   0.06055928     Pasteur
#> 14  -0.2474880790    0.2460880           0.7578508   0.06055928     Pasteur
#> 15   0.4584285626    0.2460880           0.7578508   0.06055928     Pasteur
#> 16  -0.1522102556    0.2460880           0.7578508   0.06055928     Pasteur
#> 17  -0.5081722451    0.2460880           0.7578508   0.06055928     Pasteur
#> 18  -0.0770187318    0.2460880           0.7578508   0.06055928     Pasteur
#> 19   0.5627157155    0.2460880           0.7578508   0.06055928     Pasteur
#> 20   0.4062682214    0.2460880           0.7578508   0.06055928     Pasteur
#> 21  -0.0459497806    0.2460880           0.7578508   0.06055928     Pasteur
#> 22   0.8444544127    0.2460880           0.7578508   0.06055928     Pasteur
#> 23  -0.2424655900    0.2460880           0.7578508   0.06055928     Pasteur
#> 24  -0.4640496724    0.2460880           0.7578508   0.06055928     Pasteur
#> 25  -0.4620663771    0.2460880           0.7578508   0.06055928     Pasteur
#> 26   0.2408515785    0.2460880           0.7578508   0.06055928     Pasteur
#> 27   0.7381939298    0.2460880           0.7578508   0.06055928     Pasteur
#> 28  -0.1301675700    0.2460880           0.7578508   0.06055928     Pasteur
#> 29   0.6148660991    0.2460880           0.7578508   0.06055928     Pasteur
#> 30  -0.7508377521    0.2460880           0.7578508   0.06055928     Pasteur
#> 31  -0.5723405940    0.2460880           0.7578508   0.06055928     Pasteur
#> 32  -1.0235601413    0.2460880           0.7578508   0.06055928     Pasteur
#> 33  -0.3377198338    0.2460880           0.7578508   0.06055928     Pasteur
#> 34  -0.1311524028    0.2460880           0.7578508   0.06055928     Pasteur
#> 35  -0.6344987619    0.2460880           0.7578508   0.06055928     Pasteur
#> 36  -0.8220321431    0.2460880           0.7578508   0.06055928     Pasteur
#> 37   0.4453918847    0.2460880           0.7578508   0.06055928     Pasteur
#> 38  -0.2765365489    0.2460880           0.7578508   0.06055928     Pasteur
#> 39   0.5075633242    0.2460880           0.7578508   0.06055928     Pasteur
#> 40   0.1305302102    0.2460880           0.7578508   0.06055928     Pasteur
#> 41  -0.4520180851    0.2460880           0.7578508   0.06055928     Pasteur
#> 42  -0.6816299704    0.2460880           0.7578508   0.06055928     Pasteur
#> 43  -0.8681210674    0.2460880           0.7578508   0.06055928     Pasteur
#> 44  -0.3557587429    0.2460880           0.7578508   0.06055928     Pasteur
#> 45   0.1556017972    0.2460880           0.7578508   0.06055928     Pasteur
#> 46   0.7893256009    0.2460880           0.7578508   0.06055928     Pasteur
#> 47   0.8364266019    0.2460880           0.7578508   0.06055928     Pasteur
#> 48  -0.6465336632    0.2460880           0.7578508   0.06055928     Pasteur
#> 49  -0.0007985378    0.2460880           0.7578508   0.06055928     Pasteur
#> 50  -0.5643330410    0.2460880           0.7578508   0.06055928     Pasteur
#> 51   0.4794897372    0.2460880           0.7578508   0.06055928     Pasteur
#> 52  -0.0038410533    0.2460880           0.7578508   0.06055928     Pasteur
#> 53   0.2448656552    0.2460880           0.7578508   0.06055928     Pasteur
#> 54  -0.3487496446    0.2460880           0.7578508   0.06055928     Pasteur
#> 55   0.5938049244    0.2460880           0.7578508   0.06055928     Pasteur
#> 56  -0.1442027026    0.2460880           0.7578508   0.06055928     Pasteur
#> 57   0.1255515586    0.2460880           0.7578508   0.06055928     Pasteur
#> 58   0.5286211770    0.2460880           0.7578508   0.06055928     Pasteur
#> 59  -0.2745396395    0.2460880           0.7578508   0.06055928     Pasteur
#> 60  -0.3788234628    0.2460880           0.7578508   0.06055928     Pasteur
#> 61   0.3250726897    0.2460880           0.7578508   0.06055928     Pasteur
#> 62   0.3270931788    0.2460880           0.7578508   0.06055928     Pasteur
#> 63  -0.3808472736    0.2460880           0.7578508   0.06055928     Pasteur
#> 64  -0.5121962794    0.2460880           0.7578508   0.06055928     Pasteur
#> 65   0.2889948716    0.2460880           0.7578508   0.06055928     Pasteur
#> 66   0.0984425466    0.2460880           0.7578508   0.06055928     Pasteur
#> 67   0.6248874975    0.2460880           0.7578508   0.06055928     Pasteur
#> 68   0.3671950391    0.2460880           0.7578508   0.06055928     Pasteur
#> 69  -0.8490803818    0.2460880           0.7578508   0.06055928     Pasteur
#> 70  -0.7528412973    0.2460880           0.7578508   0.06055928     Pasteur
#> 71  -1.0897051443    0.2460880           0.7578508   0.06055928     Pasteur
#> 72   0.4153147447    0.2460880           0.7578508   0.06055928     Pasteur
#> 73  -0.0439193340    0.2460880           0.7578508   0.06055928     Pasteur
#> 74   0.2639096626    0.2460880           0.7578508   0.06055928     Pasteur
#> 75  -0.2615029538    0.2460880           0.7578508   0.06055928     Pasteur
#> 76  -0.1933138925    0.2460880           0.7578508   0.06055928     Pasteur
#> 77  -0.0599815914    0.2460880           0.7578508   0.06055928     Pasteur
#> 78  -0.7157075652    0.2460880           0.7578508   0.06055928     Pasteur
#> 79   0.6740086372    0.2460880           0.7578508   0.06055928     Pasteur
#> 80  -0.1331526341    0.2460880           0.7578508   0.06055928     Pasteur
#> 81  -0.3407251556    0.2460880           0.7578508   0.06055928     Pasteur
#> 82   0.7111187976    0.2460880           0.7578508   0.06055928     Pasteur
#> 83  -0.0178625717    0.2460880           0.7578508   0.06055928     Pasteur
#> 84  -0.2905952611    0.2460880           0.7578508   0.06055928     Pasteur
#> 85  -0.9182712197    0.2460880           0.7578508   0.06055928     Pasteur
#> 86   0.0974407701    0.2460880           0.7578508   0.06055928     Pasteur
#> 87   0.2097524120    0.2460880           0.7578508   0.06055928     Pasteur
#> 88   0.1014648044    0.2460880           0.7578508   0.06055928     Pasteur
#> 89   0.1937339841    0.2460880           0.7578508   0.06055928     Pasteur
#> 90  -0.1221397592    0.2460880           0.7578508   0.06055928     Pasteur
#> 91  -0.0258668029    0.2460880           0.7578508   0.06055928     Pasteur
#> 92   0.2318356131    0.2460880           0.7578508   0.06055928     Pasteur
#> 93  -0.2073659610    0.2460880           0.7578508   0.06055928     Pasteur
#> 94  -0.8039696544    0.2460880           0.7578508   0.06055928     Pasteur
#> 95  -0.7207503119    0.2460880           0.7578508   0.06055928     Pasteur
#> 96   0.6289012316    0.2460880           0.7578508   0.06055928     Pasteur
#> 97  -0.8851954015    0.2460880           0.7578508   0.06055928     Pasteur
#> 98   0.1074485544    0.2460880           0.7578508   0.06055928     Pasteur
#> 99  -0.5512930413    0.2460880           0.7578508   0.06055928     Pasteur
#> 100  0.1265297555    0.2460880           0.7578508   0.06055928     Pasteur
#> 101  0.6710099590    0.2460880           0.7578508   0.06055928     Pasteur
#> 102 -0.1131337515    0.2460880           0.7578508   0.06055928     Pasteur
#> 103 -0.2224095137    0.2460880           0.7578508   0.06055928     Pasteur
#> 104 -0.4109144482    0.2460880           0.7578508   0.06055928     Pasteur
#> 105  0.9758067324    0.2460880           0.7578508   0.06055928     Pasteur
#> 106  0.2017448590    0.2460880           0.7578508   0.06055928     Pasteur
#> 107 -0.5112180824    0.2460880           0.7578508   0.06055928     Pasteur
#> 108  0.5276160865    0.2460880           0.7578508   0.06055928     Pasteur
#> 109 -0.2685492537    0.2460880           0.7578508   0.06055928     Pasteur
#> 110 -0.2625014085    0.2460880           0.7578508   0.06055928     Pasteur
#> 111  0.8966047962    0.2460880           0.7578508   0.06055928     Pasteur
#> 112  0.0964423154    0.2460880           0.7578508   0.06055928     Pasteur
#> 113  0.7441979376    0.2460880           0.7578508   0.06055928     Pasteur
#> 114 -0.2645285334    0.2460880           0.7578508   0.06055928     Pasteur
#> 115  0.5065718556    0.2460880           0.7578508   0.06055928     Pasteur
#> 116 -0.6615908301    0.2460880           0.7578508   0.06055928     Pasteur
#> 117 -0.0178559281    0.2460880           0.7578508   0.06055928     Pasteur
#> 118  0.4584488203    0.2460880           0.7578508   0.06055928     Pasteur
#> 119  0.2538952425    0.2460880           0.7578508   0.06055928     Pasteur
#> 120 -0.4349776229    0.2460880           0.7578508   0.06055928     Pasteur
#> 121 -0.4089142248    0.2460880           0.7578508   0.06055928     Pasteur
#> 122 -0.0900487738    0.2460880           0.7578508   0.06055928     Pasteur
#> 123  0.8845595870    0.2460880           0.7578508   0.06055928     Pasteur
#> 124  0.5837768824    0.2460880           0.7578508   0.06055928     Pasteur
#> 125  0.7211332177    0.2460880           0.7578508   0.06055928     Pasteur
#> 126  0.3751959485    0.2460880           0.7578508   0.06055928     Pasteur
#> 127 -0.5883995374    0.2460880           0.7578508   0.06055928     Pasteur
#> 128 -0.5021715592    0.2460880           0.7578508   0.06055928     Pasteur
#> 129  0.4142993540    0.2460880           0.7578508   0.06055928     Pasteur
#> 130 -0.4209424903    0.2460880           0.7578508   0.06055928     Pasteur
#> 131  0.7191466084    0.2460880           0.7578508   0.06055928     Pasteur
#> 132  0.3862357247    0.2460880           0.7578508   0.06055928     Pasteur
#> 133 -0.2895934845    0.2460880           0.7578508   0.06055928     Pasteur
#> 134 -0.2665085069    0.2460880           0.7578508   0.06055928     Pasteur
#> 135  0.3069899433    0.2460880           0.7578508   0.06055928     Pasteur
#> 136  0.5687164015    0.2460880           0.7578508   0.06055928     Pasteur
#> 137  0.6860604822    0.2460880           0.7578508   0.06055928     Pasteur
#> 138  0.1395631194    0.2460880           0.7578508   0.06055928     Pasteur
#> 139  0.8103834537    0.2460880           0.7578508   0.06055928     Pasteur
#> 140 -0.0539506900    0.2460880           0.7578508   0.06055928     Pasteur
#> 141 -0.6475321179    0.2460880           0.7578508   0.06055928     Pasteur
#> 142  0.8073781319    0.2460880           0.7578508   0.06055928     Pasteur
#> 143  0.8133618820    0.2460880           0.7578508   0.06055928     Pasteur
#> 144  0.7652388467    0.2460880           0.7578508   0.06055928     Pasteur
#> 145 -0.9844397998    0.2460880           0.7578508   0.06055928     Pasteur
#> 146 -0.3697938833    0.2460880           0.7578508   0.06055928     Pasteur
#> 147 -0.5553034537    0.2460880           0.7578508   0.06055928     Pasteur
#> 148 -0.3868512736    0.2460880           0.7578508   0.06055928     Pasteur
#> 149  0.0332790570    0.2460880           0.7578508   0.06055928     Pasteur
#> 150  0.7943311540    0.2460880           0.7578508   0.06055928     Pasteur
#> 151  0.8444577267    0.2460880           0.7578508   0.06055928     Pasteur
#> 152 -0.0058210268    0.2460880           0.7578508   0.06055928     Pasteur
#> 153 -0.1271453123    0.2460880           0.7578508   0.06055928     Pasteur
#> 154 -0.0529522353    0.2460880           0.7578508   0.06055928     Pasteur
#> 155 -0.3969028952    0.2460880           0.7578508   0.06055928     Pasteur
#> 156 -0.5152016012    0.2460880           0.7578508   0.06055928     Pasteur
#> 157 -1.0105865497    0.2051125           0.8500167   0.04207115 Grant-White
#> 158 -0.1639869772    0.2051125           0.8500167   0.04207115 Grant-White
#> 159  0.1464895899    0.2051125           0.8500167   0.04207115 Grant-White
#> 160 -0.6395921201    0.2051125           0.8500167   0.04207115 Grant-White
#> 161 -0.6703227143    0.2051125           0.8500167   0.04207115 Grant-White
#> 162  0.0102507587    0.2051125           0.8500167   0.04207115 Grant-White
#> 163  1.0423360918    0.2051125           0.8500167   0.04207115 Grant-White
#> 164 -0.6283440905    0.2051125           0.8500167   0.04207115 Grant-White
#> 165 -0.4987987089    0.2051125           0.8500167   0.04207115 Grant-White
#> 166 -0.8637203668    0.2051125           0.8500167   0.04207115 Grant-White
#> 167 -0.7965090547    0.2051125           0.8500167   0.04207115 Grant-White
#> 168 -0.6909779191    0.2051125           0.8500167   0.04207115 Grant-White
#> 169 -0.3181655436    0.2051125           0.8500167   0.04207115 Grant-White
#> 170 -0.9272746177    0.2051125           0.8500167   0.04207115 Grant-White
#> 171  0.2985310257    0.2051125           0.8500167   0.04207115 Grant-White
#> 172 -0.0855505439    0.2051125           0.8500167   0.04207115 Grant-White
#> 173  0.6588751951    0.2051125           0.8500167   0.04207115 Grant-White
#> 174  0.0729074909    0.2051125           0.8500167   0.04207115 Grant-White
#> 175 -0.0189387402    0.2051125           0.8500167   0.04207115 Grant-White
#> 176  0.2018321995    0.2051125           0.8500167   0.04207115 Grant-White
#> 177 -0.4659538877    0.2051125           0.8500167   0.04207115 Grant-White
#> 178  0.4414633016    0.2051125           0.8500167   0.04207115 Grant-White
#> 179 -0.0560401262    0.2051125           0.8500167   0.04207115 Grant-White
#> 180  0.5241774642    0.2051125           0.8500167   0.04207115 Grant-White
#> 181 -0.1718988642    0.2051125           0.8500167   0.04207115 Grant-White
#> 182 -0.2853207224    0.2051125           0.8500167   0.04207115 Grant-White
#> 183 -0.8688727113    0.2051125           0.8500167   0.04207115 Grant-White
#> 184 -0.9388206576    0.2051125           0.8500167   0.04207115 Grant-White
#> 185 -0.5100467336    0.2051125           0.8500167   0.04207115 Grant-White
#> 186 -0.8555069865    0.2051125           0.8500167   0.04207115 Grant-White
#> 187 -0.4437787572    0.2051125           0.8500167   0.04207115 Grant-White
#> 188 -0.4343257699    0.2051125           0.8500167   0.04207115 Grant-White
#> 189 -0.7432841403    0.2051125           0.8500167   0.04207115 Grant-White
#> 190 -0.7879993982    0.2051125           0.8500167   0.04207115 Grant-White
#> 191 -0.3941879956    0.2051125           0.8500167   0.04207115 Grant-White
#> 192 -0.3199834895    0.2051125           0.8500167   0.04207115 Grant-White
#> 193 -0.4933237017    0.2051125           0.8500167   0.04207115 Grant-White
#> 194 -0.9658941954    0.2051125           0.8500167   0.04207115 Grant-White
#> 195 -0.8320939831    0.2051125           0.8500167   0.04207115 Grant-White
#> 196 -0.2093229080    0.2051125           0.8500167   0.04207115 Grant-White
#> 197 -0.6344168671    0.2051125           0.8500167   0.04207115 Grant-White
#> 198 -0.6179935846    0.2051125           0.8500167   0.04207115 Grant-White
#> 199 -0.7715761206    0.2051125           0.8500167   0.04207115 Grant-White
#> 200 -0.3941879956    0.2051125           0.8500167   0.04207115 Grant-White
#> 201  0.3572309424    0.2051125           0.8500167   0.04207115 Grant-White
#> 202 -0.5292066503    0.2051125           0.8500167   0.04207115 Grant-White
#> 203 -0.6946349755    0.2051125           0.8500167   0.04207115 Grant-White
#> 204  0.0254538599    0.2051125           0.8500167   0.04207115 Grant-White
#> 205 -0.4358439616    0.2051125           0.8500167   0.04207115 Grant-White
#> 206  0.2945742102    0.2051125           0.8500167   0.04207115 Grant-White
#> 207 -0.2321169775    0.2051125           0.8500167   0.04207115 Grant-White
#> 208 -0.9725876401    0.2051125           0.8500167   0.04207115 Grant-White
#> 209  0.3085588739    0.2051125           0.8500167   0.04207115 Grant-White
#> 210 -0.2941530317    0.2051125           0.8500167   0.04207115 Grant-White
#> 211 -0.1001329721    0.2051125           0.8500167   0.04207115 Grant-White
#> 212 -0.0365804603    0.2051125           0.8500167   0.04207115 Grant-White
#> 213 -0.5733012352    0.2051125           0.8500167   0.04207115 Grant-White
#> 214 -0.0143824212    0.2051125           0.8500167   0.04207115 Grant-White
#> 215 -0.9597985152    0.2051125           0.8500167   0.04207115 Grant-White
#> 216 -0.0244349120    0.2051125           0.8500167   0.04207115 Grant-White
#> 217 -0.0760975616    0.2051125           0.8500167   0.04207115 Grant-White
#> 218 -0.1500023184    0.2051125           0.8500167   0.04207115 Grant-White
#> 219  0.5348030818    0.2051125           0.8500167   0.04207115 Grant-White
#> 220 -0.3826402167    0.2051125           0.8500167   0.04207115 Grant-White
#> 221 -0.2090019893    0.2051125           0.8500167   0.04207115 Grant-White
#> 222 -0.8165911327    0.2051125           0.8500167   0.04207115 Grant-White
#> 223  0.7592328168    0.2051125           0.8500167   0.04207115 Grant-White
#> 224 -0.9126692810    0.2051125           0.8500167   0.04207115 Grant-White
#> 225 -0.4473900064    0.2051125           0.8500167   0.04207115 Grant-White
#> 226  0.8352517957    0.2051125           0.8500167   0.04207115 Grant-White
#> 227 -0.1928767219    0.2051125           0.8500167   0.04207115 Grant-White
#> 228 -0.5617551904    0.2051125           0.8500167   0.04207115 Grant-White
#> 229 -0.3525320296    0.2051125           0.8500167   0.04207115 Grant-White
#> 230  0.8240037710    0.2051125           0.8500167   0.04207115 Grant-White
#> 231 -0.8913916692    0.2051125           0.8500167   0.04207115 Grant-White
#> 232 -0.7885971676    0.2051125           0.8500167   0.04207115 Grant-White
#> 233 -0.5586976375    0.2051125           0.8500167   0.04207115 Grant-White
#> 234  0.6351589644    0.2051125           0.8500167   0.04207115 Grant-White
#> 235 -0.7305179240    0.2051125           0.8500167   0.04207115 Grant-White
#> 236 -0.3352112332    0.2051125           0.8500167   0.04207115 Grant-White
#> 237 -0.3336701331    0.2051125           0.8500167   0.04207115 Grant-White
#> 238 -0.6496199683    0.2051125           0.8500167   0.04207115 Grant-White
#> 239  0.8267404051    0.2051125           0.8500167   0.04207115 Grant-White
#> 240 -0.1278042793    0.2051125           0.8500167   0.04207115 Grant-White
#> 241  0.9632789856    0.2051125           0.8500167   0.04207115 Grant-White
#> 242  0.4709508157    0.2051125           0.8500167   0.04207115 Grant-White
#> 243 -0.3072172682    0.2051125           0.8500167   0.04207115 Grant-White
#> 244  0.5372417006    0.2051125           0.8500167   0.04207115 Grant-White
#> 245 -0.8077570794    0.2051125           0.8500167   0.04207115 Grant-White
#> 246  0.6451902857    0.2051125           0.8500167   0.04207115 Grant-White
#> 247  0.2027297182    0.2051125           0.8500167   0.04207115 Grant-White
#> 248 -0.2412261375    0.2051125           0.8500167   0.04207115 Grant-White
#> 249  0.3341142151    0.2051125           0.8500167   0.04207115 Grant-White
#> 250 -1.0476879357    0.2051125           0.8500167   0.04207115 Grant-White
#> 251  0.3304571588    0.2051125           0.8500167   0.04207115 Grant-White
#> 252  0.9842797517    0.2051125           0.8500167   0.04207115 Grant-White
#> 253 -0.0490469273    0.2051125           0.8500167   0.04207115 Grant-White
#> 254 -0.2236073210    0.2051125           0.8500167   0.04207115 Grant-White
#> 255 -0.3732083989    0.2051125           0.8500167   0.04207115 Grant-White
#> 256  0.0756441249    0.2051125           0.8500167   0.04207115 Grant-White
#> 257  0.5865115385    0.2051125           0.8500167   0.04207115 Grant-White
#> 258  0.3201295564    0.2051125           0.8500167   0.04207115 Grant-White
#> 259  0.9070176878    0.2051125           0.8500167   0.04207115 Grant-White
#> 260 -1.0808342452    0.2051125           0.8500167   0.04207115 Grant-White
#> 261 -0.0901051239    0.2051125           0.8500167   0.04207115 Grant-White
#> 262 -0.6395921201    0.2051125           0.8500167   0.04207115 Grant-White
#> 263 -0.5885060704    0.2051125           0.8500167   0.04207115 Grant-White
#> 264 -0.4291716864    0.2051125           0.8500167   0.04207115 Grant-White
#> 265 -0.2950505554    0.2051125           0.8500167   0.04207115 Grant-White
#> 266  0.0278924788    0.2051125           0.8500167   0.04207115 Grant-White
#> 267  0.0239374023    0.2051125           0.8500167   0.04207115 Grant-White
#> 268  0.6615889257    0.2051125           0.8500167   0.04207115 Grant-White
#> 269  0.8598867146    0.2051125           0.8500167   0.04207115 Grant-White
#> 270 -0.7578894770    0.2051125           0.8500167   0.04207115 Grant-White
#> 271 -0.3884132391    0.2051125           0.8500167   0.04207115 Grant-White
#> 272  0.6920197706    0.2051125           0.8500167   0.04207115 Grant-White
#> 273 -0.2698178670    0.2051125           0.8500167   0.04207115 Grant-White
#> 274 -0.5492869892    0.2051125           0.8500167   0.04207115 Grant-White
#> 275  0.4177259013    0.2051125           0.8500167   0.04207115 Grant-White
#> 276 -0.7536346463    0.2051125           0.8500167   0.04207115 Grant-White
#> 277 -0.6617884202    0.2051125           0.8500167   0.04207115 Grant-White
#> 278 -0.3969246297    0.2051125           0.8500167   0.04207115 Grant-White
#> 279 -0.1430091195    0.2051125           0.8500167   0.04207115 Grant-White
#> 280 -0.2293574349    0.2051125           0.8500167   0.04207115 Grant-White
#> 281 -0.6213526306    0.2051125           0.8500167   0.04207115 Grant-White
#> 282  0.3578287118    0.2051125           0.8500167   0.04207115 Grant-White
#> 283 -0.2749931200    0.2051125           0.8500167   0.04207115 Grant-White
#> 284 -0.1822264666    0.2051125           0.8500167   0.04207115 Grant-White
#> 285 -0.5748194269    0.2051125           0.8500167   0.04207115 Grant-White
#> 286 -0.2235844174    0.2051125           0.8500167   0.04207115 Grant-White
#> 287 -0.1071261710    0.2051125           0.8500167   0.04207115 Grant-White
#> 288 -0.4559243005    0.2051125           0.8500167   0.04207115 Grant-White
#> 289 -0.8308738066    0.2051125           0.8500167   0.04207115 Grant-White
#> 290 -0.1004327263    0.2051125           0.8500167   0.04207115 Grant-White
#> 291 -0.8655383128    0.2051125           0.8500167   0.04207115 Grant-White
#> 292 -0.6380739235    0.2051125           0.8500167   0.04207115 Grant-White
#> 293 -0.3820195437    0.2051125           0.8500167   0.04207115 Grant-White
#> 294 -0.0523830649    0.2051125           0.8500167   0.04207115 Grant-White
#> 295 -0.4851120654    0.2051125           0.8500167   0.04207115 Grant-White
#> 296  0.3037062787    0.2051125           0.8500167   0.04207115 Grant-White
#> 297 -0.5188790479    0.2051125           0.8500167   0.04207115 Grant-White
#> 298 -0.4939461137    0.2051125           0.8500167   0.04207115 Grant-White
#> 299 -0.3181672777    0.2051125           0.8500167   0.04207115 Grant-White
#> 300 -0.9081129669    0.2051125           0.8500167   0.04207115 Grant-White
#> 301  0.3636017293    0.2051125           0.8500167   0.04207115 Grant-White

# Product-score indicator for the ind60 x dem60 interaction (single-group
# lavaan models only, v1); see compute_fs_prod() for the derivation
get_fs(PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")],
       model = " ind60 =~ x1 + x2 + x3
                 dem60 =~ y1 + y2 + y3 + y4 ",
       product = "ind60:dem60")
#>       fs_ind60    fs_dem60 fs_ind60_se fs_dem60_se ind60_by_fs_ind60
#> 1  -0.54258816 -2.74640573   0.1245694   0.6307323         0.9553858
#> 2   0.12647664 -2.85646114   0.1245694   0.6307323         0.9553858
#> 3   0.73408891  2.74401728   0.1245694   0.6307323         0.9553858
#> 4   1.25253604  3.10856431   0.1245694   0.6307323         0.9553858
#> 5   0.83355267  1.92455641   0.1245694   0.6307323         0.9553858
#> 6   0.22426801  1.02292332   0.1245694   0.6307323         0.9553858
#> 7   0.12517739  1.00406461   0.1245694   0.6307323         0.9553858
#> 8   0.11783867 -0.37216403   0.1245694   0.6307323         0.9553858
#> 9   0.25175134 -1.24897911   0.1245694   0.6307323         0.9553858
#> 10  0.39938631  2.85267059   0.1245694   0.6307323         0.9553858
#> 11  0.67497777  1.41959595   0.1245694   0.6307323         0.9553858
#> 12  0.56462020  1.08769844   0.1245694   0.6307323         0.9553858
#> 13  1.31236592  1.54090232   0.1245694   0.6307323         0.9553858
#> 14  0.23246021  1.77370863   0.1245694   0.6307323         0.9553858
#> 15  0.58638481  2.45676871   0.1245694   0.6307323         0.9553858
#> 16  0.38404785  2.35887573   0.1245694   0.6307323         0.9553858
#> 17  0.05076465  0.04034088   0.1245694   0.6307323         0.9553858
#> 18 -0.01747337 -1.86718064   0.1245694   0.6307323         0.9553858
#> 19  0.90920762  3.61477756   0.1245694   0.6307323         0.9553858
#> 20  1.12553557  0.88355273   0.1245694   0.6307323         0.9553858
#> 21  0.97202590  3.62673300   0.1245694   0.6307323         0.9553858
#> 22  0.87820036 -3.02428925   0.1245694   0.6307323         0.9553858
#> 23  0.57540754 -1.51695438   0.1245694   0.6307323         0.9553858
#> 24  0.66221224  2.76341635   0.1245694   0.6307323         0.9553858
#> 25  0.92358281  2.00507336   0.1245694   0.6307323         0.9553858
#> 26 -0.89353051 -0.92008050   0.1245694   0.6307323         0.9553858
#> 27 -0.13984744 -1.19025576   0.1245694   0.6307323         0.9553858
#> 28 -0.53828496 -1.01247764   0.1245694   0.6307323         0.9553858
#> 29 -0.80834865  0.10456709   0.1245694   0.6307323         0.9553858
#> 30 -1.25324343 -0.71847055   0.1245694   0.6307323         0.9553858
#> 31 -0.33373641 -1.61401581   0.1245694   0.6307323         0.9553858
#> 32 -1.17441075 -3.27250363   0.1245694   0.6307323         0.9553858
#> 33 -0.12409974 -1.17530231   0.1245694   0.6307323         0.9553858
#> 34 -0.04239173 -0.53796274   0.1245694   0.6307323         0.9553858
#> 35 -0.34010528  0.74552889   0.1245694   0.6307323         0.9553858
#> 36 -0.58953870  1.61018662   0.1245694   0.6307323         0.9553858
#> 37  0.17453657 -0.28144814   0.1245694   0.6307323         0.9553858
#> 38 -0.54457243  0.37694690   0.1245694   0.6307323         0.9553858
#> 39 -1.05196602 -0.62919501   0.1245694   0.6307323         0.9553858
#> 40 -0.05504697 -0.03346842   0.1245694   0.6307323         0.9553858
#> 41 -0.12364358 -0.38394102   0.1245694   0.6307323         0.9553858
#> 42 -0.59058710  1.35347275   0.1245694   0.6307323         0.9553858
#> 43 -0.11968796  0.89227782   0.1245694   0.6307323         0.9553858
#> 44 -1.07176064 -2.08481096   0.1245694   0.6307323         0.9553858
#> 45 -1.09139097 -2.07944291   0.1245694   0.6307323         0.9553858
#> 46 -0.83287255  1.59590721   0.1245694   0.6307323         0.9553858
#> 47 -1.14519896 -1.53352201   0.1245694   0.6307323         0.9553858
#> 48 -0.56115378  2.08138051   0.1245694   0.6307323         0.9553858
#> 49  0.06493340 -1.04044137   0.1245694   0.6307323         0.9553858
#> 50  0.15671638  1.72618633   0.1245694   0.6307323         0.9553858
#> 51  0.34626130 -1.24967043   0.1245694   0.6307323         0.9553858
#> 52 -0.45158373 -2.31742576   0.1245694   0.6307323         0.9553858
#> 53  0.43233465 -1.07533341   0.1245694   0.6307323         0.9553858
#> 54  0.25779725 -0.02904676   0.1245694   0.6307323         0.9553858
#> 55  0.51730650 -2.78207923   0.1245694   0.6307323         0.9553858
#> 56  0.20104991 -2.49001474   0.1245694   0.6307323         0.9553858
#> 57  0.25318620 -2.52145147   0.1245694   0.6307323         0.9553858
#> 58  0.72354623  1.86717109   0.1245694   0.6307323         0.9553858
#> 59  0.24619740 -0.93321102   0.1245694   0.6307323         0.9553858
#> 60  1.21681210  3.19853937   0.1245694   0.6307323         0.9553858
#> 61  0.18167599 -3.15685030   0.1245694   0.6307323         0.9553858
#> 62 -1.16605067 -3.41334680   0.1245694   0.6307323         0.9553858
#> 63 -0.86491026 -3.11864398   0.1245694   0.6307323         0.9553858
#> 64  0.10990059 -0.47238885   0.1245694   0.6307323         0.9553858
#> 65 -0.07376176  2.95292007   0.1245694   0.6307323         0.9553858
#> 66 -0.28782931 -1.96509718   0.1245694   0.6307323         0.9553858
#> 67 -0.02508160  2.96218478   0.1245694   0.6307323         0.9553858
#> 68 -1.31843215 -1.59567027   0.1245694   0.6307323         0.9553858
#> 69 -0.40462357 -1.79146161   0.1245694   0.6307323         0.9553858
#> 70 -0.55568363 -1.01578892   0.1245694   0.6307323         0.9553858
#> 71 -0.71308015  0.08818212   0.1245694   0.6307323         0.9553858
#> 72  0.31014319  1.70765911   0.1245694   0.6307323         0.9553858
#> 73  0.79092897  1.86102556   0.1245694   0.6307323         0.9553858
#> 74  0.08770237  3.12885767   0.1245694   0.6307323         0.9553858
#> 75 -0.14138149 -2.41398025   0.1245694   0.6307323         0.9553858
#>    ind60_by_fs_dem60 dem60_by_fs_ind60 dem60_by_fs_dem60 ev_fs_ind60
#> 1           0.181827       0.005867694         0.8688887  0.01551752
#> 2           0.181827       0.005867694         0.8688887  0.01551752
#> 3           0.181827       0.005867694         0.8688887  0.01551752
#> 4           0.181827       0.005867694         0.8688887  0.01551752
#> 5           0.181827       0.005867694         0.8688887  0.01551752
#> 6           0.181827       0.005867694         0.8688887  0.01551752
#> 7           0.181827       0.005867694         0.8688887  0.01551752
#> 8           0.181827       0.005867694         0.8688887  0.01551752
#> 9           0.181827       0.005867694         0.8688887  0.01551752
#> 10          0.181827       0.005867694         0.8688887  0.01551752
#> 11          0.181827       0.005867694         0.8688887  0.01551752
#> 12          0.181827       0.005867694         0.8688887  0.01551752
#> 13          0.181827       0.005867694         0.8688887  0.01551752
#> 14          0.181827       0.005867694         0.8688887  0.01551752
#> 15          0.181827       0.005867694         0.8688887  0.01551752
#> 16          0.181827       0.005867694         0.8688887  0.01551752
#> 17          0.181827       0.005867694         0.8688887  0.01551752
#> 18          0.181827       0.005867694         0.8688887  0.01551752
#> 19          0.181827       0.005867694         0.8688887  0.01551752
#> 20          0.181827       0.005867694         0.8688887  0.01551752
#> 21          0.181827       0.005867694         0.8688887  0.01551752
#> 22          0.181827       0.005867694         0.8688887  0.01551752
#> 23          0.181827       0.005867694         0.8688887  0.01551752
#> 24          0.181827       0.005867694         0.8688887  0.01551752
#> 25          0.181827       0.005867694         0.8688887  0.01551752
#> 26          0.181827       0.005867694         0.8688887  0.01551752
#> 27          0.181827       0.005867694         0.8688887  0.01551752
#> 28          0.181827       0.005867694         0.8688887  0.01551752
#> 29          0.181827       0.005867694         0.8688887  0.01551752
#> 30          0.181827       0.005867694         0.8688887  0.01551752
#> 31          0.181827       0.005867694         0.8688887  0.01551752
#> 32          0.181827       0.005867694         0.8688887  0.01551752
#> 33          0.181827       0.005867694         0.8688887  0.01551752
#> 34          0.181827       0.005867694         0.8688887  0.01551752
#> 35          0.181827       0.005867694         0.8688887  0.01551752
#> 36          0.181827       0.005867694         0.8688887  0.01551752
#> 37          0.181827       0.005867694         0.8688887  0.01551752
#> 38          0.181827       0.005867694         0.8688887  0.01551752
#> 39          0.181827       0.005867694         0.8688887  0.01551752
#> 40          0.181827       0.005867694         0.8688887  0.01551752
#> 41          0.181827       0.005867694         0.8688887  0.01551752
#> 42          0.181827       0.005867694         0.8688887  0.01551752
#> 43          0.181827       0.005867694         0.8688887  0.01551752
#> 44          0.181827       0.005867694         0.8688887  0.01551752
#> 45          0.181827       0.005867694         0.8688887  0.01551752
#> 46          0.181827       0.005867694         0.8688887  0.01551752
#> 47          0.181827       0.005867694         0.8688887  0.01551752
#> 48          0.181827       0.005867694         0.8688887  0.01551752
#> 49          0.181827       0.005867694         0.8688887  0.01551752
#> 50          0.181827       0.005867694         0.8688887  0.01551752
#> 51          0.181827       0.005867694         0.8688887  0.01551752
#> 52          0.181827       0.005867694         0.8688887  0.01551752
#> 53          0.181827       0.005867694         0.8688887  0.01551752
#> 54          0.181827       0.005867694         0.8688887  0.01551752
#> 55          0.181827       0.005867694         0.8688887  0.01551752
#> 56          0.181827       0.005867694         0.8688887  0.01551752
#> 57          0.181827       0.005867694         0.8688887  0.01551752
#> 58          0.181827       0.005867694         0.8688887  0.01551752
#> 59          0.181827       0.005867694         0.8688887  0.01551752
#> 60          0.181827       0.005867694         0.8688887  0.01551752
#> 61          0.181827       0.005867694         0.8688887  0.01551752
#> 62          0.181827       0.005867694         0.8688887  0.01551752
#> 63          0.181827       0.005867694         0.8688887  0.01551752
#> 64          0.181827       0.005867694         0.8688887  0.01551752
#> 65          0.181827       0.005867694         0.8688887  0.01551752
#> 66          0.181827       0.005867694         0.8688887  0.01551752
#> 67          0.181827       0.005867694         0.8688887  0.01551752
#> 68          0.181827       0.005867694         0.8688887  0.01551752
#> 69          0.181827       0.005867694         0.8688887  0.01551752
#> 70          0.181827       0.005867694         0.8688887  0.01551752
#> 71          0.181827       0.005867694         0.8688887  0.01551752
#> 72          0.181827       0.005867694         0.8688887  0.01551752
#> 73          0.181827       0.005867694         0.8688887  0.01551752
#> 74          0.181827       0.005867694         0.8688887  0.01551752
#> 75          0.181827       0.005867694         0.8688887  0.01551752
#>    ecov_fs_dem60_fs_ind60 ev_fs_dem60 fs_ind60:fs_dem60 fs_ind60:fs_dem60_se
#> 1             0.005632564   0.3978232        0.84766209            0.4836628
#> 2             0.005632564   0.3978232       -1.00378073            0.4836628
#> 3             0.005632564   0.3978232        1.37184752            0.4836628
#> 4             0.005632564   0.3978232        3.25108370            0.4836628
#> 5             0.005632564   0.3978232        0.96171401            0.4836628
#> 6             0.005632564   0.3978232       -0.41309615            0.4836628
#> 7             0.005632564   0.3978232       -0.51681895            0.4836628
#> 8             0.005632564   0.3978232       -0.68636044            0.4836628
#> 9             0.005632564   0.3978232       -0.95693729            0.4836628
#> 10            0.005632564   0.3978232        0.49681244            0.4836628
#> 11            0.005632564   0.3978232        0.31569057            0.4836628
#> 12            0.005632564   0.3978232       -0.02836862            0.4836628
#> 13            0.005632564   0.3978232        1.37972256            0.4836628
#> 14            0.005632564   0.3978232       -0.23018846            0.4836628
#> 15            0.005632564   0.3978232        0.79810673            0.4836628
#> 16            0.005632564   0.3978232        0.26341602            0.4836628
#> 17            0.005632564   0.3978232       -0.64045724            0.4836628
#> 18            0.005632564   0.3978232       -0.60987920            0.4836628
#> 19            0.005632564   0.3978232        2.64407816            0.4836628
#> 20            0.005632564   0.3978232        0.35196489            0.4836628
#> 21            0.005632564   0.3978232        2.88277329            0.4836628
#> 22            0.005632564   0.3978232       -3.29843704            0.4836628
#> 23            0.005632564   0.3978232       -1.51537211            0.4836628
#> 24            0.005632564   0.3978232        1.18746300            0.4836628
#> 25            0.005632564   0.3978232        1.20934616            0.4836628
#> 26            0.005632564   0.3978232        0.17961487            0.4836628
#> 27            0.005632564   0.3978232       -0.47605091            0.4836628
#> 28            0.005632564   0.3978232       -0.09750365            0.4836628
#> 29            0.005632564   0.3978232       -0.72703179            0.4836628
#> 30            0.005632564   0.3978232        0.25791336            0.4836628
#> 31            0.005632564   0.3978232       -0.10384929            0.4836628
#> 32            0.005632564   0.3978232        3.20075831            0.4836628
#> 33            0.005632564   0.3978232       -0.49665041            0.4836628
#> 34            0.005632564   0.3978232       -0.61969996            0.4836628
#> 35            0.005632564   0.3978232       -0.89606344            0.4836628
#> 36            0.005632564   0.3978232       -1.59177245            0.4836628
#> 37            0.005632564   0.3978232       -0.69162812            0.4836628
#> 38            0.005632564   0.3978232       -0.84778002            0.4836628
#> 39            0.005632564   0.3978232        0.01938665            0.4836628
#> 40            0.005632564   0.3978232       -0.64066279            0.4836628
#> 41            0.005632564   0.3978232       -0.59503329            0.4836628
#> 42            0.005632564   0.3978232       -1.44184867            0.4836628
#> 43            0.005632564   0.3978232       -0.74930004            0.4836628
#> 44            0.005632564   0.3978232        1.59191320            0.4836628
#> 45            0.005632564   0.3978232        1.62698008            0.4836628
#> 46            0.005632564   0.3978232       -1.97169243            0.4836628
#> 47            0.005632564   0.3978232        1.11368267            0.4836628
#> 48            0.005632564   0.3978232       -1.81047968            0.4836628
#> 49            0.005632564   0.3978232       -0.71006453            0.4836628
#> 50            0.005632564   0.3978232       -0.37198346            0.4836628
#> 51            0.005632564   0.3978232       -1.07521763            0.4836628
#> 52            0.005632564   0.3978232        0.40400663            0.4836628
#> 53            0.005632564   0.3978232       -1.10740902            0.4836628
#> 54            0.005632564   0.3978232       -0.64999330            0.4836628
#> 55            0.005632564   0.3978232       -2.08169280            0.4836628
#> 56            0.005632564   0.3978232       -1.14312236            0.4836628
#> 57            0.005632564   0.3978232       -1.28090185            0.4836628
#> 58            0.005632564   0.3978232        0.70847948            0.4836628
#> 59            0.005632564   0.3978232       -0.87225925            0.4836628
#> 60            0.005632564   0.3978232        3.24951627            0.4836628
#> 61            0.005632564   0.3978232       -1.21602903            0.4836628
#> 62            0.005632564   0.3978232        3.33763021            0.4836628
#> 63            0.005632564   0.3978232        2.05484206            0.4836628
#> 64            0.005632564   0.3978232       -0.69442094            0.4836628
#> 65            0.005632564   0.3978232       -0.86031770            0.4836628
#> 66            0.005632564   0.3978232       -0.07689257            0.4836628
#> 67            0.005632564   0.3978232       -0.71680146            0.4836628
#> 68            0.005632564   0.3978232        1.46127786            0.4836628
#> 69            0.005632564   0.3978232        0.08236246            0.4836628
#> 70            0.005632564   0.3978232       -0.07804786            0.4836628
#> 71            0.005632564   0.3978232       -0.70538605            0.4836628
#> 72            0.005632564   0.3978232       -0.11288628            0.4836628
#> 73            0.005632564   0.3978232        0.82943391            0.4836628
#> 74            0.005632564   0.3978232       -0.36809690            0.4836628
#> 75            0.005632564   0.3978232       -0.30121299            0.4836628
#>    fs_ind60:fs_dem60_ld
#> 1             0.8311908
#> 2             0.8311908
#> 3             0.8311908
#> 4             0.8311908
#> 5             0.8311908
#> 6             0.8311908
#> 7             0.8311908
#> 8             0.8311908
#> 9             0.8311908
#> 10            0.8311908
#> 11            0.8311908
#> 12            0.8311908
#> 13            0.8311908
#> 14            0.8311908
#> 15            0.8311908
#> 16            0.8311908
#> 17            0.8311908
#> 18            0.8311908
#> 19            0.8311908
#> 20            0.8311908
#> 21            0.8311908
#> 22            0.8311908
#> 23            0.8311908
#> 24            0.8311908
#> 25            0.8311908
#> 26            0.8311908
#> 27            0.8311908
#> 28            0.8311908
#> 29            0.8311908
#> 30            0.8311908
#> 31            0.8311908
#> 32            0.8311908
#> 33            0.8311908
#> 34            0.8311908
#> 35            0.8311908
#> 36            0.8311908
#> 37            0.8311908
#> 38            0.8311908
#> 39            0.8311908
#> 40            0.8311908
#> 41            0.8311908
#> 42            0.8311908
#> 43            0.8311908
#> 44            0.8311908
#> 45            0.8311908
#> 46            0.8311908
#> 47            0.8311908
#> 48            0.8311908
#> 49            0.8311908
#> 50            0.8311908
#> 51            0.8311908
#> 52            0.8311908
#> 53            0.8311908
#> 54            0.8311908
#> 55            0.8311908
#> 56            0.8311908
#> 57            0.8311908
#> 58            0.8311908
#> 59            0.8311908
#> 60            0.8311908
#> 61            0.8311908
#> 62            0.8311908
#> 63            0.8311908
#> 64            0.8311908
#> 65            0.8311908
#> 66            0.8311908
#> 67            0.8311908
#> 68            0.8311908
#> 69            0.8311908
#> 70            0.8311908
#> 71            0.8311908
#> 72            0.8311908
#> 73            0.8311908
#> 74            0.8311908
#> 75            0.8311908
```
