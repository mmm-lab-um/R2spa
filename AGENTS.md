# R2spa - AGENTS.md

## Overview
`R2spa` — an R package implementing **two-stage path analysis (2S-PA)**, a latent-variable /
structural-equation-modeling (SEM) technique. Stage 1 extracts factor scores and their
observation-specific SEs/reliability from a `lavaan` CFA (or `lme4` mixed model). Stage 2 feeds
those scores and error-variance estimates into a `lavaan::sem()` path model that corrects for
measurement error. Core machinery: factor-score computation, schema-driven stage-2 model
assembly. The delta-method SE correction (`vcov_corrected`), grand standardization, latent
interaction, and the OpenMx path are currently **quarantined** in `.quarantine/` while the
`get_fs()`/`tspa()` output contracts settle (see `archive/PLAN_QUARANTINE.md`).

## Repository Facts
- ~2,500 lines of R across 7 files in `R/`; 9 test files in `tests/testthat/`.
- `.quarantine/` — quarantined consumers of `get_fs()`/`tspa()` (`R/`, `tests/`,
  `vignettes/`), excluded from the package build via `^\.quarantine$` in `.Rbuildignore`.
  Not part of the build, tests, or docs — never modify package code to match it, and don't
  "fix" it unless working on re-integration. Re-integration: `git mv` files back, then
  `document()` → `test()` → `check()`. Its test files are self-contained with provenance
  headers (plan: `archive/PLAN_QUARANTINE.md`).
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
   all `Suggests`). Expected finding: 1 NOTE — `'OpenMx' in DESCRIPTION Imports but not
   imported from` (OpenMx stays in Imports until `R/tspa_mx.R` is re-integrated; the check
   is otherwise expected to be 0 errors / 0 warnings).

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
  `psych`, `rmarkdown`, `testthat (>= 3.0.0)`, `tidyr`, `umx`. Note: `OpenMx` is in
  `Imports` but imported from nowhere in staying code (its only consumer,
  `.quarantine/R/tspa_mx.R`, is quarantined) → the one expected check NOTE.
- `.github/workflows/` — `R-CMD-check.yaml` and `pkgdown.yaml`. Catch failures locally first.

## Dependency Rules
- **No heavyweight frameworks** (`tidyverse`, `dplyr`, `purrr`) when base R or existing imports
  (`MASS`; `Matrix` (Suggests) for merMod ops) suffice. `Imports` is deliberately minimal (4:
  `lavaan`, `lme4`, `MASS`, `OpenMx`) — keep it. `OpenMx` stays until `tspa_mx.R` is
  re-integrated; dropping it needs the quarantined OpenMx path settled first.
- No `%>%` in `R/` — `magrittr` is `Suggests`-only. Use base-R `|>` if piping is needed.
- Any new `Suggests`/`Imports` requires a `DESCRIPTION` edit and justification. Prefer existing deps.

## Naming / Style
- Exported: `snake_case` (current inventory: `get_fs`, `tspa`, `get_fs_lavaan`, `get_fs_lmer`,
  `compute_fscore`, `augment_lav_predict`, `fs_to_group_list`, `block_diag`; quarantined:
  `get_fs_int`, `tspa_mx_model`, `vcov_corrected`,
  `grand_standardized_solution`/`grandStandardizedSolution` — the legacy CamelCase alias pair
  must be kept in sync on re-integration).
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

Quarantined (in `.quarantine/R/`, do not build against): `get_fs_int.R` (latent
interaction), `tspa_corrected_se.R` (`vcov_corrected()`, sensitive to `fsT`/`fsL`
shapes), `grandStandardizedSolution.R` (multigroup standardization, Jacobian SEs),
`tspa_mx.R` (`tspa_mx_model()`, OpenMx path).

## General Instruction
Trust and follow the rules above exactly. Never hand-edit `NAMESPACE` or `man/*.Rd`, never call
`library()`/`require()` in function bodies, always `devtools::document()` before testing when
roxygen changed, prefer existing `Imports`/base R over new dependencies.
