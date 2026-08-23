# R2spa - AGENTS.md

## Overview
`R2spa` — an R package implementing **two-stage path analysis (2S-PA)**, a latent-variable /
structural-equation-modeling (SEM) technique. Stage 1 extracts factor scores and their
observation-specific SEs/reliability from a `lavaan` CFA (or `lme4` mixed model). Stage 2 feeds
those scores and error-variance estimates into a `lavaan::sem()` path model that corrects for
measurement error. Core machinery: factor-score computation, schema-driven stage-2 model
assembly. The first-order delta-method SE correction (`vcov_corrected()`, also exposed
in-place via `tspa(corrected_se = TRUE)`), multigroup grand standardization, and the
`tspa_mx_model()` OpenMx path have all been **re-integrated** into `R/` (2026-08). Only
latent interaction (`get_fs_int()`) remains **quarantined** in `.quarantine/` while the
remaining contracts settle (see `archive/PLAN_QUARANTINE.md`).

## Repository Facts
- ~2,500 lines of R across 7 files in `R/`; 9 test files in `tests/testthat/`.
- `.quarantine/` — quarantined consumers of `get_fs()`/`tspa()` (`R/`, `tests/`,
  `vignettes/`), excluded from the package build via `^\.quarantine$` in `.Rbuildignore`.
  Not part of the build, tests, or docs — never modify package code to match it, and don't
  "fix" it unless working on re-integration. Re-integration: `git mv` files back, then
   `document()` → `test()` → `check()`. Its test files are self-contained with provenance
   headers (plan: `archive/PLAN_QUARANTINE.md`). Re-integration log: `tspa_corrected_se.R` (`vcov_corrected()`),
    `grandStandardizedSolution.R`, and `tspa_mx.R`/`tspa_mx_model()` are re-integrated into
    `R/` (2026-08); only `get_fs_int.R` remains. The `corrected-se.Rmd` vignette (and its
    `boo_separate.RDS`/`boo_joint.RDS` fixtures, now shared with the corrected-SE tests) was
    re-integrated to `vignettes/` when `tspa(corrected_se)` gained multigroup support and the
    corrected grand-standardized SE path (2026-08). The `tspa()` `tspa_args` attribute
   (self-contained evaluated argument list) and the `tsp_set_vcov()` lavaan-compat boundary
   (the single `@vcov[["vcov"]]` write behind `corrected_se`) are staying-code features
   consumed by the now-package-internal `vcov_corrected()`.
- `legacy/` — `tspa_plot.R` (diagnostic plotting, removed from package). `archive/` — 
  `tspa-plot-vignette.Rmd` + completed plan files (`PLAN_01` … `PLAN_04`, `PLAN_QUARANTINE`).
  Both directories are ignored for development.
- **Unmaintained for ~2 years** — last commit Nov 2024 (version 0.0.4, still "developmental").
  Do not assume CI is currently green; verify before trusting existing behavior as a spec.
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
  `compute_fscore`, `augment_lav_predict`, `fs_to_group_list`, `block_diag`, `tspa_mx_model`,
  `vcov_corrected`, `grand_standardized_solution`/`grandStandardizedSolution` — the legacy
  CamelCase alias pair is kept in sync; quarantined: `get_fs_int` only).
- Column conventions from `get_fs()`/`get_fs_int()` (the latter is quarantined — the
  conventions remain the spec for re-integration) — downstream functions parse by name:
  - `fs_<name>` score | `<name>_se` SE | `ev_<name>` error variance
  - `ecov_<name1>_<name2>` error covariance | `<indicator>_by_<name>` implied loading
  - Attributes: `fsT` (error cov), `fsL` (loadings), `fsb` (intercepts), `scoring_matrix`
    (lavaan: per-group score×item matrices; `merMod`: named list of per-cluster
    `num_re × n_j` matrices)
- Roxygen: markdown (`Roxygen: list(markdown = TRUE)` in `DESCRIPTION`).
- Internal helpers: `snake_case`, co-located with their caller — don't move to shared `utils.R`.

## Prioritized `R/` File Reference
1. **`get_fscore.R`** (~510 lines) — `get_fs()` S3 generic, legacy wrappers (`get_fs_lavaan()`,
   `get_fs_lmer()`), `fs_to_group_list()`, helpers (`augment_fs()`, `check_blocks_identical()`,
   `assemble_fs_blocks()`).
2. **`get_fs_methods.R`** (~610 lines) — S3 methods (`get_fs.data.frame()`, `get_fs.default()`,
   `get_fs.lavaan()`, `get_fs.merMod()`), block builders (`get_fs_blocks.lavaan()`,
   `get_fs_blocks.merMod()`). Future `get_fs.mirt()` here.
3. **`get_fscore_math.R`** (~540 lines) — `compute_fscore()`, `augment_lav_predict()`,
   `compute_a*`, `compute_fspars()`, `correct_evfs()`, `compute_evfs()`, `compute_ldfs()`,
   `compute_fsrel()`. Pure math, no S3. Touchpoint for SE bugs, missing data, multigroup.
4. **`tspa.R`** (~530 lines) — `tspa()` entrypoint; owned partable stage-2 schema
   (`tspa_schema_sf()`/`tspa_schema_mf()` → `tspa_render()`), product-score auto-alias
   (`tspa_sf_alias()`); `tspa_sf()`/`tspa_mf()` emit the model string fed to `lavaan::sem()`.
5. **`lavaan_compat.R`** (~275 lines) — `tsp_*` wrappers, the only file that reads lavaan
   internals (layout/partable probing, tested-up-to version canary). Currently consumed only
   by its own canary tests (`test-lavaan_compat.R`); its package consumers are quarantined.
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

Quarantined (in `.quarantine/R/`, do not build against): `get_fs_int.R` (latent interaction)
**only**. `tspa_corrected_se.R` (`vcov_corrected()`), `grandStandardizedSolution.R`
(multigroup standardization), and `tspa_mx.R`/`tspa_mx_model()` (OpenMx) were re-integrated
to `R/` in 2026-08.

## General Instruction
Trust and follow the rules above exactly. Never hand-edit `NAMESPACE` or `man/*.Rd`, never call
`library()`/`require()` in function bodies, always `devtools::document()` before testing when
roxygen changed, prefer existing `Imports`/base R over new dependencies.
