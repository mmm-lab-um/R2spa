# R2spa - AGENTS.md

## Overview
`R2spa` — an R package implementing **two-stage path analysis (2S-PA)**, a latent-variable /
structural-equation-modeling (SEM) technique. Stage 1 extracts factor scores and their
observation-specific SEs/reliability from a `lavaan` CFA, `lme4` mixed model, or `mirt` IRT
fit. Stage 2 feeds those scores and error-variance estimates into a `lavaan::sem()` path model
that corrects for measurement error (or an OpenMx model via `tspa_mx_model()`). Core
machinery: factor-score computation, schema-driven stage-2 model assembly. Both stage-2
entrypoints auto-derive the measurement inputs (`fsT`/`fsL`/`fsb`/`se_fs`) from a `get_fs()`
result when omitted (PLANs 13/15). The first-order delta-method SE correction (`vcov_corrected()`, also exposed
in-place via `tspa(corrected_se = TRUE)`), multigroup grand standardization, and the
`tspa_mx_model()` OpenMx path have all been **re-integrated** into `R/` (2026-08). The
latent-interaction function `get_fs_int()` was **removed** and replaced by the
joint-model `compute_fs_prod()` / `get_fs(product = )` (2026-08, branch
`rejoin/fs-prod`); its column conventions (`fs_a:fs_b`, `fs_a:fs_b_se`,
`fs_a:fs_b_ld`) are kept, and `tspa()` fixes the product indicators'
measurement-error covariances in the stage-2 model.

## Repository Facts
- ~8,000 lines of R across 12 files in `R/`; 26 test files in `tests/testthat/`.
- `.quarantine/` — quarantined consumers of `get_fs()`/`tspa()` (`tests/` — only `_snaps/`
  remains; `vignettes/`; the `R/` subdirectory was deleted when its last files were
  re-integrated or removed), excluded from the package build via `^\.quarantine$` in `.Rbuildignore`.
  Not part of the build, tests, or docs — never modify package code to match it, and don't
  "fix" it unless working on re-integration. Re-integration: `git mv` files back, then
   `document()` → `test()` → `check()`. Its test files are self-contained with provenance
   headers (plan: `archive/PLAN_QUARANTINE.md`). Re-integration log: `tspa_corrected_se.R` (`vcov_corrected()`),
    `grandStandardizedSolution.R`, and `tspa_mx.R`/`tspa_mx_model()` are re-integrated into
    `R/` (2026-08). The `corrected-se.Rmd` vignette (and its
    `boo_separate.RDS`/`boo_joint.RDS` fixtures, now shared with the corrected-SE tests) was
    re-integrated to `vignettes/` when `tspa(corrected_se)` gained multigroup support and the
    corrected grand-standardized SE path (2026-08). `get_fs_int.R` (latent interaction)
    and its test file were **deleted** (2026-08, branch `rejoin/fs-prod`) — superseded by
     `R/compute_fs_prod.R` / `get_fs(product = )` with the same column conventions; the
     quarantined `get_fs_int-vignette.Rmd` vignette was rewritten on top of
     `compute_fs_prod()` / `get_fs(product = )` / `tspa(product = TRUE)` (2026-08-26) and
     re-integrated as `vignettes/product-factor-scores.Rmd` (2026-08-27, with the
     product-indicator error-covariance fix); `categorical-interaction.Rmd` remains
     stale. The `tspa()` `tspa_args` attribute
   (self-contained evaluated argument list) and the `tsp_set_vcov()` lavaan-compat boundary
   (the single `@vcov[["vcov"]]` write behind `corrected_se`) are staying-code features
   consumed by the now-package-internal `vcov_corrected()`.
- `legacy/` — `tspa_plot.R` (diagnostic plotting, removed from package). `archive/` — 
  `tspa-plot-vignette.Rmd` + completed plan files (`PLAN_01` … `PLAN_15`,
  `PLAN_QUARANTINE`). Both directories are ignored for development.
- **Actively developed** — intensive 2026-08 re-integration + plan work (PLAN 06–15; see
  `STATUS.md` for the full issue log). Version 0.0.4 is still "developmental". Suite
  ~4,050 expectations passing, 0 fail (2026-08-27); last full `R CMD check` 0/0/0 was
  2026-08-26 — re-run after code changes.
- Target dev environment: Linux (WSL/Ubuntu-like), R 4.6.1.
- No `TODO`/`FIXME`/`HACK` markers in the codebase.
- No `library()`/`require()` in function bodies — only in roxygen `@examples` blocks.

## Build / Test Lifecycle — Exact Order
This package uses `devtools` + `roxygen2` + `testthat` (edition 3). Never skip or reorder.

1. `devtools::load_all()` — after any `R/*.R` edit.
2. `devtools::document()` — **mandatory** after any roxygen tag change (`@export`, `@importFrom`,
   `@param`, new/removed function). Regenerates `NAMESPACE` and `man/*.Rd`.
3. `devtools::test()` — after every functional change, always after step 2 if exports changed.
4. `devtools::check()` — before marking a change complete. Expect slowness (vignettes pull
   all `Suggests`). Expected finding: 0 errors / 0 warnings / 0 NOTEs (`R/tspa_mx.R` was
   re-integrated 2026-08, so `OpenMx` is imported in `NAMESPACE`; the former
   `'OpenMx' in Imports but not imported from` NOTE is obsolete).

**Hard rules, no exceptions:**
- **Never hand-edit `NAMESPACE` or `man/*.Rd`** — they are machine-generated. Edit roxygen in
  `R/` and rerun `devtools::document()`.
- **Never add `library()`/`require()` in a function body.** Namespace every external call
  (`lavaan::sem()`) or import via `@importFrom pkg fun` in the function's roxygen block, then
  run `devtools::document()`.
- New function from an already-imported package requires a matching `@importFrom` tag before
  running `devtools::document()`.

## Project Layout
- `R/` — all package logic, one module per file. No `src/`. See prioritized reference below.
- `tests/testthat/` — edition 3 (`Config/testthat/edition: 3` in `DESCRIPTION`). Tests validate
  numerically against `lavaan` reference fits and hand-calculated matrices using
  `expect_equal(..., tolerance = ...)`, `ignore_attr = TRUE`. Attribute-carrying data frames are
  intentional — don't "fix" them without verifying.
- `vignettes/` — many `.Rmd` with cached `.RDS`/`.csv` fixtures. Don't regenerate casually; they
  back specific vignette narratives.
- `DESCRIPTION` — `Imports:` `lavaan`, `lme4`, `MASS`, `OpenMx`. `Suggests:` `boot`,
   `DiagrammeR`, `ggplot2`, `knitr`, `lintr`, `magrittr`, `Matrix`, `mirt`, `numDeriv`,
   `psych`, `rmarkdown`, `testthat (>= 3.0.0)`, `tidyr`, `umx`. `OpenMx` is consumed by
   `R/tspa_mx.R` (re-integrated 2026-08) and imported in `NAMESPACE`.
- `.github/workflows/` — `R-CMD-check.yaml` and `pkgdown.yaml`. Catch failures locally first.

## Dependency Rules
- **No heavyweight frameworks** (`tidyverse`, `dplyr`, `purrr`) when base R or existing imports
   (`MASS`; `Matrix` (Suggests) for merMod ops) suffice. `Imports` is deliberately minimal (4:
   `lavaan`, `lme4`, `MASS`, `OpenMx`) — keep it. `OpenMx` is consumed by `R/tspa_mx.R`
   (re-integrated 2026-08); dropping it needs that OpenMx path settled first.
- No `%>%` in `R/` — `magrittr` is `Suggests`-only. Use base-R `|>` if piping is needed.
- Any new `Suggests`/`Imports` requires a `DESCRIPTION` edit and justification. Prefer existing deps.

## Naming / Style
- Exported: `snake_case` (current inventory: `get_fs`, `tspa`, `get_fs_lavaan`, `get_fs_lmer`,
  `compute_fscore`, `compute_fs_prod`, `fs_indiv`, `augment_lav_predict`, `fs_to_group_list`,
  `block_diag`, `tspa_mx_model`, `vcov_corrected`,
  `grand_standardized_solution`/`grandStandardizedSolution`
  — the legacy CamelCase alias pair is kept in sync). S3 methods on `get_fs()`:
  `data.frame`, `default`, `lavaan`, `merMod`, and `mirt`'s `SingleGroupClass`/
  `MultipleGroupClass` (`mirt` is `Suggests`-only, guarded by `require_mirt()`).
- Column conventions from `get_fs()` — downstream functions parse by name:
  - `fs_<name>` score | `<name>_se` SE | `ev_<name>` error variance
  - `ecov_<name1>_<name2>` error covariance | `<indicator>_by_<name>` implied loading
  - `fs_a:fs_b` DMC product indicator | `fs_a:fs_b_se` per-row SE | `fs_a:fs_b_ld` implied
    loading (from `get_fs(product = )`/`compute_fs_prod()`; `tspa_sf_alias()` maps
    `fs_a:fs_b` to the `fs_ab` model name)
  - Attributes: `fsT` (error cov), `fsL` (loadings), `fsb` (intercepts), `scoring_matrix`
    (lavaan: per-group score×item matrices; `merMod`: named list of per-cluster
    `num_re × n_j` matrices), `psi`/`alpha` (effective latent covariance/mean —
    `prior_cov`/`prior_mean` if supplied, else the model estimate; mirt MG: per-group
    list), `fs_pattern` (per-case observed-indicator pattern; missing data),
    `group_col` (merMod cluster / mirt group column name), `mirt_per_obs`/`per_obs`
    markers (per-row attribute lists), `pooled_fs` (`tspa(reduce = )` marker)
- Roxygen: markdown (`Roxygen: list(markdown = TRUE)` in `DESCRIPTION`).
- Internal helpers: `snake_case`, co-located with their caller — don't move to shared `utils.R`.

## Prioritized `R/` File Reference
1. **`get_fscore.R`** (~750 lines) — `get_fs()` S3 generic, legacy wrappers (`get_fs_lavaan()`,
   `get_fs_lmer()`), `fs_to_group_list()`, helpers (`augment_fs()`, `assemble_fs_blocks()`;
   `check_blocks_identical()` was deleted with PLAN 06 — per-pattern blocks are kept, not
   dropped).
2. **`get_fs_methods.R`** (~2,000 lines) — S3 methods (`get_fs.data.frame()`,
   `get_fs.default()`, `get_fs.lavaan()`, `get_fs.merMod()`, `mirt`'s
   `get_fs.SingleGroupClass()`/`get_fs.MultipleGroupClass()`), block builders
   (`get_fs_blocks.lavaan()`, `get_fs_blocks.merMod()`); `local = TRUE` per-construct scoring
   internals (`get_fs_local()`, `merge_local_fs()`, …; PLAN 14); `mirt` helpers
   (`mirt_full_cov()`, `mirt_group_pars()`, `require_mirt()`).
3. **`get_fscore_math.R`** (~760 lines) — `compute_fscore()`, `augment_lav_predict()`,
   `compute_a*`, `compute_fspars()`, `correct_evfs()`, `compute_evfs()`, `compute_ldfs()`,
   `compute_fsrel()`. Pure math, no S3. Touchpoint for SE bugs, missing data, multigroup.
4. **`tspa.R`** (~1,630 lines) — `tspa()` entrypoint; auto-derivation of the measurement
   inputs (PLAN 13: with `fsT`/`fsL` omitted, the multi-factor inputs are derived from a
   `get_fs()` result's attributes via `derive_sf_se_fs()`/`fs_group_order()`, and a
   single-factor `se_fs` from the `fs_<v>`/`fs_<v>_se` columns — explicit arguments always
   win), per-unit pooling (PLAN 09: `reduce = "mean"/"median"` collapses per-pattern /
   per-cluster attributes via `is_per_unit_fs()`/`pool_per_unit()`/`pool_se_fs()`); owned
   partable stage-2 schema
   (`tspa_schema_sf()`/`tspa_schema_mf()` → `tspa_render()`), product-score auto-alias
    (`tspa_sf_alias()`), opt-in `product = TRUE` auto-compute (detects model
    latents naming the product of two factor scores — concatenated (`xm`) or
    lavaan interaction syntax (`x:m`, rewritten to the render name
    `xm` via `tspa_rewrite_product_toks()`; a non-score `a:b` passes
    through to lavaan) — via `tspa_product_latents()`, computes missing
    DMC product columns via `compute_fs_prod()` with
    `tspa_ensure_product_cols()`, joins the pooled product SE into `se_fs`
    on the sf path (an explicit product SE may be keyed by the `a:b`
    token, stored check.names()-ed as `x.m`), emits fixed
    `gamma`/`se_P^2` product rows on the mf path via `prods` in
    `tspa_schema_mf()`, and fixes the product-indicator
    measurement-error covariances — pairs of product latents sharing a
    factor score have correlated errors (the shared score's error enters
    both) — on both the sf and mf paths, `product = TRUE` and manual flows
    alike, via `tspa_prod_ecov()` → `fs_prod_ecov()` (the manual flow,
    product columns pre-computed and listed in `se_fs`, is detected via
    `tspa_sf_alias()`'s `prod_map`)); `tspa_sf()`/`tspa_mf()` emit the
    model string fed to `lavaan::sem()`.
5. **`lavaan_compat.R`** (~390 lines) — `tsp_*` wrappers, the only file that reads lavaan
   internals (layout/partable probing, tested-up-to version canary). Consumed by its own
   canary tests (`test-lavaan_compat.R`) plus in-package: `tsp_set_vcov()` (the in-place
   `tspa(corrected_se = )` covariance overwrite) and `tsp_beta_names()` (fixed-slope
   reporting in `grandStandardizedSolution()`, PLAN 12).
6. **`helper.R`** / **`globals.R`** — `block_diag()`, NSE NOTE suppression. Low risk.
7. **`tspa_corrected_se.R`** — `vcov_corrected()`, the first-order (delta-method) corrected-SE
   path (re-integrated 2026-08-23). Central-difference stage-2 Jacobian `J` over the free
   `fsL`/`fsT` elements at `h0 = 1e-5` (refits via `do.call(tspa, tspa_args)`), returning
   `vcov(fit) + J %*% vfsLT %*% t(J)`; in-place via `tspa(corrected_se = TRUE, vfsLT = ...)`,
   which overwrites the fit covariance through the `tsp_set_vcov()` lavaan-compat boundary
    (`fit@vcov[["vcov"]]` only, `est.std` unchanged → `standardizedSolution()` reports
    corrected std SEs). Single- **and** multi-group: `which_free` positions run per group
    (group 1's `fsL`, group 2's, … then all groups' lower-tri `fsT`) in the same order as the
    `vfsLT` attribute; a double-correction guard rejects an already-`tspa_corrected` fit.
    `grandStandardizedSolution()` threads the corrected covariance, so a corrected fit reports
    corrected grand-standardized SEs (point estimates unchanged). Helpers:
    `tsp_tri2full_colmajor()`, `check_refit_convergence()`.
8. **`compute_fs_prod.R`** — `compute_fs_prod()` (+ `get_fs(product = )`), the
   double-mean-centered product factor-score indicators for single-group lavaan
   models (v1). Attribute-driven: per-pattern `fsL`/`fsT` via
   `resolve_fs_per_row()`, shared `psi` attribute; joint-normal SE formula
   `se_P^2 = tau_a s_b^2 + tau_b s_a^2 + s_a^2 s_b^2 + c^2 + 2 tau_ab c`
   (derivation in roxygen `@details`; pure-matrix helpers `fs_prod_se2()`,
   `fs_prod_ecov()` — the measurement-error covariance between two product
   indicators, `tau_ik c_jl + tau_il c_jk + tau_jk c_il + tau_jl c_ik +
   c_ik c_jl + c_il c_jk` with `tau_uv = L_u psi L_v'`, `c_uv = T[u,v]`
   (a diagonal `T` reduces the shared-factor pair to `tau_jl s_i^2`) — and
   `fs_prod_gamma()`, spec parser `parse_product_spec()` co-located).
   Replaces the removed quarantined `get_fs_int()`.
9. **`fs_indiv.R`** (~670 lines) — exported `fs_indiv()`: re-derives the individual-specific
   columns per row (`_se`, `<lvs>_by_<lv>_*`, `ev_*`/`ecov_*`, per-pattern intercepts) from the
   row's `fsL`/`fsT`/`fsb` via the shared value-only engine `fs_row_cols()` (also used by
   `augment_lav_predict()`, so SEs are pattern-consistent); per-row (`mirt`) dispatch via
   `resolve_per_obs()`.
10. **`grandStandardizedSolution.R`** (~320 lines) — `grand_standardized_solution()` (+
   legacy CamelCase alias kept in sync), multigroup grand standardization; reports
   user-fixed structural slopes alongside free ones (PLAN 12: `out_idx` free-position anchor
   + β-dimnames fallback for fixed cells). Threads the corrected covariance, so a
   `tspa(corrected_se = TRUE)` fit reports corrected grand-standardized SEs.
11. **`tspa_mx.R`** (~610 lines) — `tspa_mx_model()`, the OpenMx stage-2 route (exact, no
   pooling); PLAN 15 `tspa_mx_derive_measurement()` auto-derives the measurement inputs from
   a `get_fs()` result (per-row/per-pattern quantities become definition-variable matrices
   over the result's own `_by_`/`ev_`/`ecov_` columns; `int_fs_*` intercept columns from the
   `fsb` attribute); `tspa_mx_defvar_col()` handles `lavaanify()`'s reversed `~~` orientation.

`.quarantine/R/` no longer exists (both of its files were deleted in 2026-08):
`get_fs_int.R` (latent interaction) was removed and replaced by
`R/compute_fs_prod.R` (2026-08, branch `rejoin/fs-prod`);
`tspa_corrected_se.R` (`vcov_corrected()`), `grandStandardizedSolution.R`
(multigroup standardization), and `tspa_mx.R`/`tspa_mx_model()` (OpenMx) were re-integrated
to `R/` earlier in 2026-08. `.quarantine/vignettes/` holds the still-stale
`categorical-interaction.Rmd` and `reliability.Rmd` (+ `sim_results_reliability.RDS`);
the rewritten product vignette was re-integrated as
`vignettes/product-factor-scores.Rmd` (2026-08-27).

## General Instruction
Trust and follow the rules above exactly. Never hand-edit `NAMESPACE` or `man/*.Rd`, never call
`library()`/`require()` in function bodies, always `devtools::document()` before testing when
roxygen changed, prefer existing `Imports`/base R over new dependencies.
