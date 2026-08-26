# R2spa 0.0.4

## New Features
- Add function `get_fs_int()` for estimating interaction effects in 2S-PA (#82).
- Add computation of reliability function in `get_fscore()` (#81)
- Add functions of obtaining tidy-ed factor scores data for `get_fscore()` (#79)
- `get_fs()` now supports multi-group `mirt` models (`MultipleGroupClass`):
  per-observation factor scores are extracted from the whole fit, carry a
  trailing `group` column (the model's group levels, `NA` for
  completely-missing rows) and a per-group (\code{list}) `psi` attribute, with
  each observation using its own group's factor covariance.
- Re-integrated `grand_standardized_solution()` / legacy alias
  `grandStandardizedSolution()` (from `.quarantine/`): grand-standardized
  path coefficients with delta-method SEs for (multi)group `lavaan` fit objects.
- `tspa()` now accepts **multi-factor mirt** factor scores — single- and
  multi-group `get_fs()` output (per-observation `fsL`/`fsT` lists) — by
  reducing them to one representative set per group via `reduce` (default
  `"mean"`), the same per-unit pooling path as FIML/merMod; completely-missing
  mirt rows are dropped from the reduction.
- `tspa()` now accepts a `get_fs()` result as `data` directly: when
  `fsT`/`fsL` (multi-factor) or `se_fs` (single-factor) are omitted they are
  derived from the result's `fsT`/`fsL`/`fsb` attributes or its
  `fs_<v>`/`fs_<v>_se` score columns; explicit arguments always take
  precedence.
- `tspa_mx_model()` now accepts a `get_fs()` result as `data` directly:
  when `se_fs`/`fsL`/`fsT`/`fsb` are all omitted, the measurement inputs
  are derived from the result's attributes — constant quantities (complete
  data, `local = TRUE`, `format = "list"`) become fixed numerics, while
  per-row (mirt, `local = TRUE` FIML) and per-pattern (single-group FIML)
  quantities become definition variables referencing the result's own
  `*_by_*`/`ev_*`/`ecov_*` columns, with the `int_fs_*` score-intercept
  columns appended automatically from the `fsb` attribute. Explicit
  arguments always take precedence, a non-`get_fs()` frame fails fast with
  an actionable message instead of the old misleading error, and there is
  no pooling (the OpenMx route is exact-or-fail, unlike
  `tspa(reduce = )`); multigroup remains refused in Phase 1.

- `get_fs()` gains `local = TRUE`: each latent is scored from its own local
  measurement model (per-construct stage 1, the canonical 2S-PA setup)
  instead of the single joint multi-factor model; the merged result carries
  the usual multi-factor attributes with exactly-zero cross-terms
  (block-diagonal `fsT`/`fsL`/`psi`, zero `ecov_*` columns) and feeds
   `tspa()` directly. `model` may be a single string (strict per-latent
   `=~` grammar) or a character vector / named list of complete single-factor
   model strings (the escape hatch). `vfsLT` (hence
   `tspa(corrected_se = TRUE)`), `prior_cov`, and `reliability` are not
   supported in `local` mode (v1).
- New exported `compute_fs_prod()` + `get_fs(product = )`: double-mean-centered
   product factor-score indicators for pairs of distinct latents in a
   single-group `lavaan` model (v1). For each requested pair (`product` as an
   `"a:b + c:d"` string, a list of length-2 name pairs, or a 2-column
   matrix/data frame) the result gains `fs_a:fs_b` (the product indicator),
   `fs_a:fs_b_se` (per-row standard error) and `fs_a:fs_b_ld` (implied
   loading). The SE uses the general joint-model formula
   `se_P^2 = tau_a s_b^2 + tau_b s_a^2 + s_a^2 s_b^2 + c^2 + 2 tau_ab c`
   (derivation in `?compute_fs_prod`), so correlated factors, cross-loadings,
   and error covariances are handled — the separate-single-factor special
   case it reduces to. Under FIML the SE/loading resolve per
   observed-indicator pattern; the output feeds `tspa()` directly via the
   existing product-score auto-alias.
- `get_fs()` gains a `product` argument forwarding to `compute_fs_prod()`
   (single-group `lavaan` models only; rejected with `local = TRUE` and for
   multi-group models, v1).
- `tspa()` gains an opt-in `product` argument (default `FALSE`): when
    `TRUE`, the double-mean-centered product indicators for every model
    latent that names the product of two of the model's factor scores — by
    concatenation (`xm` for `x` and `m`) or in lavaan's interaction syntax
    (`x:m`, rendered under the concatenated name) — are computed on the fly
    and incorporated into the stage-2 measurement model — in the
    single-factor (score-scale) path the product SE joins `se_fs`
    (per-group pooled, loading 1, as before), and in the multi-factor path
    the product latent gets a fixed implied loading `gamma` and fixed error
    variance `se_P^2` from the (pooled) `fsL`/`fsT` and the `psi` attribute.
    An `a:b` token whose parts are not both factor scores (e.g. `x:g` with
    `g` an observed covariate) is not claimed and passes through to lavaan
    as an ordinary interaction; an explicit product SE may be keyed by
    either the render name or the `a:b` token. Single-group models only
    (v1); not supported with `corrected_se = TRUE`; naming the same pair
    twice (`x:m` and `xm`) or a render name colliding with another model
    variable is an error; pre-existing product columns (either orientation)
    are used as-is, so the manual `get_fs(product = )` + `se_fs` workflow is
    unchanged.

## Breaking changes
- The quarantined student function `get_fs_int()` (and its test file) is
  removed; its functionality is superseded by `get_fs(product = )` /
  `compute_fs_prod()`, which keep its column-naming convention
  (`fs_a:fs_b`, `fs_a:fs_b_se`, `fs_a:fs_b_ld`) but use the corrected
  joint-model SE formula.

## Improvements
- The stage-2 model string attached to `tspa()` fits (the `tspaModel`
  attribute) is now rendered with normalized operator spacing (`lhs =~ rhs`,
  `lhs ~~ rhs`) and without `c()` for single-value fixed statements (e.g.
  `ind60 =~ 1 * fs_ind60` instead of `ind60=~ c(1) * fs_ind60`); multigroup
  statements keep `c(v1, v2)`; fitted models, parameter tables, and covariance
  are unchanged.
- `tspa(corrected_se = TRUE)` now supports multigroup fits, and
  `grandStandardizedSolution()` threads `vcov(object)` (or an explicit
  `acov_par = vcov(corrected_fit)`) so a corrected fit yields corrected
  grand-standardized SEs with unchanged point estimates; the
   re-integrated `corrected-se` vignette documents the in-place single-
   and multigroup correction.
- `grandStandardizedSolution()` now reports user-fixed structural slopes
  (`~ k*var`): its `est.std` is the user value rescaled by the (grand) SD
  ratio and its `se` is the first-order delta approximation, matching
  `lavaan::standardizedSolution()` (which also reports a delta SE for fixed
  slopes). Previously a single fixed structural path rejected the whole
  solution. Free slopes are unchanged (still anchored by free position); a
  structural regression outside the beta matrix (e.g. observed-on-observed)
  is still rejected.
- Add examples for `get_fs_int()` (#82).
- Update naming for `get_fscore()` (#79)
    * Rename `vc` to `ev` (error variance-covariance) for better consistency

 ## Bug Fixes
- `tspa_mx_model()` no longer mis-specifies off-diagonal factor-score
  covariances when `lavaanify()` presents a `~~` row with the (lhs, rhs)
  pair reversed relative to the score order: the definition-variable lookup
  now falls back to the transposed triangle, so a lower-triangle-only `fsT`
  (the documented and the auto-derived convention) is honored. Previously
  the `c(1)` definition-variable sentinel leaked into the model as a fixed
  unit covariance between the scores, and such fits (e.g. multi-factor mirt
  per-row models, single-group FIML per-pattern models) aborted with
  "implied covariance not positive definite".
- Fix a bug in the `se_fs` argument in `tspa()` (#90).
- `grandStandardizedSolution()` now assigns standardized estimates and SEs
  to partable rows by their global free position instead of assuming the
  structural paths follow the beta matrix column-major order; the two agree
  only by accident (e.g. single-predictor models), so models with multiple
  endogenous variables per group could report `est.std`/`se`/CIs on the
  wrong rows.

## Documentation
- Updated vignettes for:
    * tspa-growth-vignette (#50)
    * get_fs_int-vignette (#82)
    * reliability (#81)
     * missing-data (#79)
     * 2S-PA with OpenMx and IRT (mirt): the multidimensional example now
       uses the auto-derived measurement inputs (no hand-rolled
       `cross_load`/`err_cov` matrices)

## Other
- General code clean-up (#82).

# R2spa 0.0.3

- Add function `tspa_plot()` for bivariate and residual plots (#23)

- `get_fs()` gains argument `corrected_fsT` for computing corrected error estimates (#50)

- New function `vcov_corrected()` for computing corrected SEs (#39)

- New function `get_fs_lavaan()` for computing factor scores and relevant matrices directly from a `lavaan` output (#61)

- Initial support for 2S-PA with *OpenMx* with `tspa_mx()`

- Update naming of relevant matrices when computing factor scores:
    * `fsT`: error covariance of factor scores
    * `fsL`: loading matrix of factor scores
    * `fsb`: intercepts of factor scores
    * `scoring_matrix`: weights for computing factor scores from items

- New vignettes for:
    * Corrected error variance of factor scores (#50)
    * Corrected standard errors incorporating uncertainty in measurement parameters of factor scores (#39)
    * Using 2S-PA with EFA scores
    * Using 2S-PA with OpenMx and definition variables (PR #57)
    * Latent interaction with categorical indicators (#27)
    * Growth modeling

- Better error messages for `tspa()` (#53)

- Support mean structure and growth model (#36, #19)

- Clean up code with `lintr` (#33)

# R2spa 0.0.2

- Use `pkgdown` to create website, with GitHub action (#22)

- `get_fs()` now returns a list with multi-group models (#29).

- New function `grandStandardizedSolution()` computes standardized solution based on grand mean and grand SD (#13).

- `tspa()` gains argument `vc` and `cross_loadings`, which is useful for factor scores obtained from multi-factor models (#7). See `vignette("Multi-Factor Measurement Model")`. 

# R2spa 0.0.1

- **Work-In-Progress!**

- 0.0.1 version
