# R2spa 0.0.4

## New Features
- `get_fs()` gains a `brmsfit` method: Bayesian Gaussian mixed models (brms)
  are now a supported stage-1 input alongside lavaan, lme4, and mirt. The
  per-level posterior-mean random effects of the (single) random-effects term
  are the factor scores, with `fsL`/`fsT`/`scoring_matrix` from the same EB/ML
  formulas the `merMod` method uses, generalised to a `p x p` random-effects
  covariance `D` reconstructed from the posterior of the term's `sd_`/`cor_`
  hyperparameters. `method = "EB"` (default) returns the posterior-mean random
  effects; `method = "ML"` a prior-free per-cluster OLS estimate. v1 requires
  the Gaussian family and a single random grouping factor (random slopes
  supported; `p > 2` terms via the generic Cholesky). `brms`/`posterior`/
  `reformulas` are `Suggests`-only and guarded by `require_brms()`; the fixed
  and random designs are built directly from the model formula + data
  (bit-identical to `lme4::getME`), so no throwaway fit is required.
- `vcov_corrected()` and `tspa(corrected_se = TRUE)` gain an `engine`
  argument. `"fd"` (the default) is the original central finite-difference
  Jacobian; `"analytic"` evaluates the same Jacobian refit-free and
  deterministically via a closed form for the saturated single-group path
  model (PLAN 16; free structural regressions and latent variances only),
  transparently falling back to `"fd"` otherwise — multigroup, a model that
  is not exactly saturated (df > 0), or a saturated model with other free
  parameters (e.g. free latent covariances or means). The two engines agree
  to the finite-difference noise floor whenever `"analytic"` applies, and
  the analytic result is bit-reproducible (no refits, no optimizer jitter).
- `tspa_mx_model()` now auto-derives its measurement inputs from
  `get_fs.merMod()` results as well: with `se_fs`/`fsL`/`fsT`/`fsb` all
  omitted, the per-cluster 3-D `fsL`/`fsT` array attributes become
  definition-variable matrices referencing the result's own per-cluster
  `*_by_*`/`ev_*`/`ecov_*` columns (one row per cluster, exact
  non-pooled correction); `merMod` results carry no `fsb`, so the score
  intercepts stay fixed at zero.
- Product factor-score indicators now carry fixed measurement-error
  covariances in the stage-2 model: products sharing a factor score (e.g.
  `xm` and `xz`) have correlated measurement errors (the shared score's
  error enters both), and `tspa()` computes the error covariance for every
  pair of product latents from the stage-1 `fsL`/`fsT`/`psi` attributes
  (new pure-matrix helper `fs_prod_ecov()`: a diagonal score-error matrix —
  disjoint-indicator CFAs — reduces the shared-factor pair to `tau_jl
  s_i^2`, the non-shared pair's latent covariance times the shared score's
  error variance; local-mode separate-models results reduce to zero) and
  emits it as a fixed `fs_v1 ~~ fs_v2` statement. Works on both the
  single-factor (`se_fs`) and the multi-factor (`fsT`/`fsL`) paths, for both
  `tspa(product = TRUE)` and the manual workflow (product columns
  pre-computed with `get_fs(product = )` / `compute_fs_prod()` and listed in
  `se_fs`), removing the first-order (n-independent) attenuation of the
  product coefficients that the previously free covariance absorbed. Data
  lacking the stage-1 attributes (e.g. a `cbind()`'d `get_fs()` result) with
  two or more product latents is rejected with an informative error (remedy:
  pass the un-`cbind()`'d result). Single-group models (v1).
- Re-integrated the product factor-score vignette from `.quarantine/` as
  `product-factor-scores` (`vignettes/product-factor-scores.Rmd`): the
  rewritten `get_fs_int-vignette.Rmd` on `compute_fs_prod()` /
  `get_fs(product = )` / `tspa(product = )` (including the fixed
  product-indicator error covariances above; the example DGP is scaled so
  the `y` latent has population variance exactly 1, so the estimates are
  read directly against the simulated coefficients).
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
- `augment_lav_predict()`'s factor-score standard-error columns are renamed
  from the legacy `se_fs_*` to the canonical `fs_*_se` used by `get_fs()`,
  `get_fs_lavaan()`, and `fs_indiv()` (#85), so all per-row APIs now agree
  on column names. The documented `ev_*`/`ecov_*` ordering difference
  (upper-triangular in `augment_lav_predict()`, lower-triangular elsewhere)
  is intentional and unchanged.

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
- Documented the mean-structure convention difference between `tspa()`
  (lavaan, which fixes the exogenous latent mean at zero and estimates the
  factor-score mean) and `tspa_mx_model()` (OpenMx, which fixes the score
  residual means at zero and estimates the latent means): the two routes fit
  the same model, so only the unidentifiable mean split differs — compare on
  the covariance quantities. Noted in both functions' docs, the OpenMx
  vignette, and a regression test.
- Updated vignettes for:
    * tspa-growth-vignette (#50)
    * missing-data (#79) — re-integrated and modernized: the per-row
      scoring of missing-data fits is explained, and the stage-2 examples
      use `tspa_mx_model()` (per-row definition variables) with
      `tspa(reduce = )` pooling contrasted
    * 2S-PA with OpenMx and IRT (mirt): the multidimensional example now
      uses the auto-derived measurement inputs (no hand-rolled
      `cross_load`/`err_cov` matrices)
    * efa-score: complete-case (listwise) correlation instead of FIML
      `lavCor`, so stage 1 (EFA) and stage 2 share one data basis;
      `n.rotations = 1` pins `psych::fa()` against the psych 2.6.5
      `faRotations()` bug

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
