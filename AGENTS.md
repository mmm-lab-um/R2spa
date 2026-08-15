# R2spa - AGENTS.md

## Overview
`R2spa` — an R package implementing **two-stage path analysis (2S-PA)**, a latent-variable /
structural-equation-modeling (SEM) technique. Stage 1 extracts factor scores and their
observation-specific SEs/reliability from a `lavaan` CFA (or `lme4` mixed model). Stage 2 feeds
those scores and error-variance estimates into a `lavaan::sem()` path model that corrects for
measurement error. Core machinery: factor-score computation, delta-method SE correction,
lavaan/OpenMx model assembly.

## Repository Facts
- ~2,200 lines of R across 10 files in `R/`; 8 test files in `tests/testthat/`.
- `legacy/` — `tspa_plot.R` (diagnostic plotting, removed from package). `archive/` — 
  `tspa-plot-vignette.Rmd`. Both directories are ignored for development.
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
4. `devtools::check()` — before marking a change complete. Expect slowness (OpenMx compiles
  vendored code; vignettes pull all `Suggests`).

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
- `DESCRIPTION` — `Imports:` `lavaan`, `lme4`, `MASS`, `Matrix`, `OpenMx`. `Suggests:` `boot`,
  `DiagrammeR`, `ggplot2`, `knitr`, `magrittr`, `mirt`, `numDeriv`, `psych`, `rmarkdown`,
  `testthat (>= 3.0.0)`, `tidyr`, `umx`.
- `.github/workflows/` — `R-CMD-check.yaml` and `pkgdown.yaml`. Catch failures locally first.

## Dependency Rules
- **No heavyweight frameworks** (`tidyverse`, `dplyr`, `purrr`) when base R or existing imports
  (`Matrix`, `MASS`) suffice. `Imports` is deliberately 5 packages — keep it.
- No `%>%` in `R/` — `magrittr` is `Suggests`-only. Use base-R `|>` if piping is needed.
- Any new `Suggests`/`Imports` requires a `DESCRIPTION` edit and justification. Prefer existing deps.

## Naming / Style
- Exported: `snake_case` (`get_fs`, `tspa`, `get_fs_int`, `vcov_corrected`), except legacy
  CamelCase aliases kept for back-compat (`grandStandardizedSolution` ↔
  `grand_standardized_solution` — keep in sync).
- Column conventions from `get_fs()`/`get_fs_int()` — downstream functions parse by name:
  - `fs_<name>` score | `<name>_se` SE | `ev_<name>` error variance
  - `ecov_<name1>_<name2>` error covariance | `<indicator>_by_<name>` implied loading
  - Attributes: `fsT` (error cov), `fsL` (loadings), `fsb` (intercepts), `scoring_matrix`
- Roxygen: markdown (`Roxygen: list(markdown = TRUE)` in `DESCRIPTION`).
- Internal helpers: `snake_case`, co-located with their caller — don't move to shared `utils.R`.

## Prioritized `R/` File Reference
1. **`get_fscore.R`** (~420 lines) — `get_fs()` S3 generic, legacy wrappers (`get_fs_lavaan()`,
   `get_fs_lmer()`), `fs_to_group_list()`, helpers (`augment_fs()`, `check_blocks_identical()`,
   `assemble_fs_blocks()`).
2. **`get_fs_methods.R`** (~320 lines) — S3 methods (`get_fs.data.frame()`, `get_fs.lavaan()`,
   `get_fs.merMod()`), block builders (`get_fs_blocks.*`, `get_D()`). Future `get_fs.mirt()` here.
3. **`get_fscore_math.R`** (~470 lines) — `compute_fscore()`, `augment_lav_predict()`,
   `compute_a*`, `compute_fspars()`, `correct_evfs()`, `compute_evfs()`, `compute_ldfs()`,
   `compute_fsrel()`. Pure math, no S3. Touchpoint for SE bugs, missing data, multigroup.
4. **`tspa.R`** (~300 lines) — `tspa()` entrypoint; dispatches to `tspa_sf()`/`tspa_mf()`.
5. **`tspa_corrected_se.R`** (~80 lines) — `vcov_corrected()`. Sensitive to `fsT`/`fsL` shapes.
6. **`grandStandardizedSolution.R`** (~260 lines) — multigroup standardization, Jacobian SEs.
7. **`get_fs_int.R`** (~110 lines) — `get_fs_int()`, latent interaction (double-mean centering).
8. **`tspa_mx.R`** (~230 lines) — `tspa_mx_model()`, OpenMx path. No test coverage; manual-test
   via `tspa-vignette-mx.Rmd` if changed.
9. **`helper.R`** / **`globals.R`** — `block_diag()`, NSE NOTE suppression. Low risk.

## General Instruction
Trust and follow the rules above exactly. Never hand-edit `NAMESPACE` or `man/*.Rd`, never call
`library()`/`require()` in function bodies, always `devtools::document()` before testing when
roxygen changed, prefer existing `Imports`/base R over new dependencies.
