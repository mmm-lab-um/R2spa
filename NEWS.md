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
