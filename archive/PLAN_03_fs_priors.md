# Plan: user-supplied latent priors (`prior_mean` / `prior_cov`) for `get_fs()`

## Objective

- Finalize and later implement a lavaan-only `get_fs()` option for user-supplied latent
  prior means and covariance for regression/EB factor scores, using the `mean` / `cov`
  arguments of `mirt::fscores()` as the conceptual API reference (R2spa's new arguments
  are named `prior_mean` / `prior_cov`).
- Support the intended multiple-group use case by applying a shared prior across groups so
  factor scores are comparable, while keeping `tspa()` and related downstream functions
  compatible.

## Important Details

- User explicitly suggested referencing the `mean` and `cov` arguments in
  `mirt::fscores()`; the R2spa API mirrors that design but the new arguments are named
  `prior_mean` and `prior_cov` (explicit `prior_` prefix, per user decision).
- User-confirmed scope decisions:
  - Initial implementation is lavaan-only.
  - Support a shared prior only; do not implement group-specific list priors yet.
  - Support `corrected_fsT = TRUE` and `vfsLT = TRUE` with custom priors by treating the
    supplied covariance as fixed.
  - `reliability = TRUE` should error when either `prior_mean` or `prior_cov` is non-NULL.
- API shape:
  - Add `prior_mean = NULL` and `prior_cov = NULL` to:
    - `get_fs.data.frame()`
    - `get_fs.lavaan()`
    - `get_fs_lavaan()` wrapper
  - Do not add formal `prior_mean` / `prior_cov` arguments to `get_fs.merMod()` or
    `get_fs_lmer()`; instead reject non-NULL `prior_mean` or `prior_cov` passed through
    `...`.
  - In `get_fs.data.frame()`, place `prior_mean = NULL, prior_cov = NULL` before `...` so
    they are not forwarded to `lavaan::cfa()`.
- Prior semantics:
  - `prior_mean = NULL` and `prior_cov = NULL` preserve current behavior exactly.
  - `prior_mean` and `prior_cov` may be supplied independently.
  - Non-NULL values are treated as fixed external priors.
  - The supplied prior is shared across all lavaan groups.
  - `prior_mean` is a numeric vector of length q, where q is the number of latent
    variables.
  - `prior_cov` is a numeric q x q matrix; for q = 1, also accept a scalar or 1 x 1
    matrix.
  - If `prior_mean` names are supplied, they must match latent-variable names and are
    reordered to model order.
  - If `prior_cov` row/column names are supplied, they must match latent-variable names
    and are reordered to model order.
  - If names are absent, assign model latent-variable names internally.
  - Validate `prior_cov` as finite, square, symmetric, and positive definite.
- Method restrictions:
  - Allow `prior_mean` / `prior_cov` only for `method = "regression"` or alias `"EB"`.
  - Error if `prior_mean` or `prior_cov` is non-NULL and `method = "Bartlett"` or alias
    `"ML"`.
  - Error if `reliability = TRUE` and either `prior_mean` or `prior_cov` is non-NULL.
  - Allow `corrected_fsT = TRUE` with custom priors.
  - Allow `vfsLT = TRUE` with custom priors.
  - When a custom `prior_cov` is supplied to corrected-SE/vcov machinery, treat it as
    fixed; do not propagate sampling uncertainty from the external prior itself.
- Mathematical behavior:
  - `compute_fscore()` already accepts separate `psi` and `alpha` inputs, so the core
    scoring math already supports user-supplied priors.
  - For regression/EB scoring, external priors are applied as:
    - `eta ~ N(prior_mean, prior_cov)`
    - fitted lavaan `lambda`, `theta`, and `nu` remain group- or pattern-specific.
  - Derived outputs must be internally consistent:
    - `fs`
    - `fsL`
    - `fsT`
    - `fsb`
    - `scoring_matrix`
    - `se_*`
    - `ev_*`
    - `ecov_*`
  - The implied score representation is:
    - `score = fsb + fsL eta + eps`
    - `fsb = prior_mean - fsL %*% prior_mean`
    - `Var(eps) = A theta A'`
    - where `A` is the regression scoring matrix computed using the supplied
      `prior_cov`.
- Math-layer plan:
  - Add optional `psi_override = NULL` to internal functions in
    `/home/marklai/R2spa/R/get_fscore_math.R`:
    - `compute_fspars()`
    - `compute_a()`
    - `compute_evfs()`
    - `compute_ldfs()`
    - `compute_grad_ld_evfs()`
    - `vcov_ld_evfs()`
    - `correct_evfs()`
  - In `compute_fspars()`, after replacing estimated coefficients, if `psi_override` is
    non-NULL, set the group's `mat$psi <- psi_override` before scoring matrix computation.
  - When `psi_override = NULL`, existing behavior remains unchanged.
  - `compute_fsrel()` does not need changes for this release because `reliability = TRUE`
    errors when custom priors are present.
- Lavaan method wiring:
  - `get_fs.lavaan()` should normalize method, validate priors, then pass values to:
    - `correct_evfs()` when `corrected_fsT = TRUE`
    - `get_fs_blocks.lavaan()`
    - `vcov_ld_evfs()` when `vfsLT = TRUE`
  - `get_fs_blocks.lavaan()` should accept `prior_mean = NULL, prior_cov = NULL`.
  - Its internal `prepare_fs()` helper should use:
    - `psi_use <- if (is.null(prior_cov)) est$psi else prior_cov`
    - `alpha_use <- if (is.null(prior_mean)) est$alpha else prior_mean`
  - This applies to:
    - complete single-group data
    - complete multiple-group data
    - each missing-data observed pattern
- Backward compatibility:
  - `get_fs_lavaan()` should gain `prior_mean = NULL, prior_cov = NULL` after existing
    formal arguments and pass them through to `get_fs()`.
  - `get_fs_lmer()` should not gain formal arguments; merMod rejection should happen
    through `get_fs.merMod()` inspecting `...`.
  - Matrix input currently dispatches to `get_fs.data.frame()`, so `prior_mean` /
    `prior_cov` supplied via the generic should work through `...`.
- `mirt::fscores()` investigation facts:
  - mirt version observed: 1.46.1
  - `mirt::fscores()` signature includes:
    - `object`
    - `method = "EAP"`
    - `mean = NULL`
    - `cov = NULL`
    - `covdata = NULL`
    - many other arguments
  - `mirt::fscores()` passes the user arguments as `gmean = mean, gcov = cov` to
    `fscores.internal()`.
  - Documentation confirms `mean` is a custom latent-variable mean vector and `cov` is a
    custom latent covariance matrix, with NULL using values from the fitted object.
- Documentation plan:
  - Update shared `get_fs()` Roxygen block in `/home/marklai/R2spa/R/get_fscore.R`.
  - Add `@param prior_mean` and `@param prior_cov`.
  - Document:
    - NULL uses lavaan's estimated group-specific latent means/covariances.
    - Non-NULL values are fixed external priors shared across groups.
    - Only regression/EB scoring is supported.
    - `corrected_fsT` and `vfsLT` are supported but treat supplied covariance as fixed.
    - `reliability = TRUE` is unsupported with user-supplied priors.
    - The design is conceptually similar to `mirt::fscores(mean, cov)`.
  - Add a small lavaan example using `lavaan::cfa()` and
    `get_fs(fit, prior_mean = ..., prior_cov = ...)`.
  - Run `devtools::document()` after Roxygen changes.
  - Do not hand-edit `man/*.Rd` or `NAMESPACE`; no new exports are planned.
- Testing plan:
  - Add a new test file:
    - `/home/marklai/R2spa/tests/testthat/test-get_fs_priors.R`
  - Test areas:
    - Single-group custom `prior_mean` + `prior_cov` matches manual `compute_fscore()`
      with overridden `alpha` and `psi`.
    - Single-group custom `prior_mean` only matches manual `compute_fscore()` using
      model-estimated `psi` and supplied `alpha`.
    - Single-group custom `prior_cov` only matches manual `compute_fscore()` using
      supplied `psi` and model-estimated `alpha`.
    - Returned attributes match manual computation:
      - `fs`
      - `fsT`
      - `fsL`
      - `fsb`
      - `scoring_matrix`
    - Named inputs in non-model order are correctly reordered.
    - Invalid input errors:
      - wrong-length `prior_mean`
      - unknown names
      - non-square `prior_cov`
      - asymmetric `prior_cov`
      - non-positive-definite `prior_cov`
      - `method = "Bartlett"` with non-NULL `prior_mean` / `prior_cov` errors.
      - `reliability = TRUE` with non-NULL `prior_mean` / `prior_cov` errors.
    - `corrected_fsT = TRUE` with custom `prior_cov`:
      - prior-based `fsT` equals uncorrected prior `fsT` plus
        `correct_evfs(..., psi_override = prior_cov)`.
    - `corrected_fsT = TRUE` with custom `prior_mean` only:
      - correction still uses model-estimated `psi`.
    - `vfsLT = TRUE` with custom `prior_cov`:
      - `attr(fs, "vfsLT")` exists with expected dimensions.
      - `vcov_corrected()` runs on a simple `tspa()` fit using prior-adjusted
        `fsT` / `fsL`.
    - Data-frame entry point:
      - `get_fs(data, model, group = ..., prior_mean = ..., prior_cov = ...)` matches
        explicit `lavaan::cfa()` then `get_fs(fit, ...)`.
      - `prior_mean` / `prior_cov` are not accidentally forwarded to `cfa()`.
    - Multiple groups:
      - shared prior applied to all groups.
      - with group-equal measurement parameters, identical raw indicator values across
        groups produce the same regression score under the same prior.
      - `format = "list"` and `format = "unified"` preserve per-group attributes.
      - `fs_to_group_list()` still works.
    - Missing data:
      - `get_fs_blocks.lavaan()` with custom priors matches manual `compute_fscore()`
        calls per observed pattern.
    - Non-lavaan methods:
      - `get_fs()` on merMod with non-NULL `prior_mean` or `prior_cov` errors.
      - `get_fs_lmer()` with non-NULL `prior_mean` or `prior_cov` errors.
  - Verification lifecycle:
    - After implementation:
      - `devtools::load_all()`
      - `devtools::document()`
      - `devtools::test()`
      - `devtools::check()`
    - Expected:
      - existing tests remain green
      - new prior tests pass
      - `man/get_fs.Rd` regenerated
      - no manual `NAMESPACE` changes needed
  - Out of scope for first release:
    - group-specific `prior_mean` / `prior_cov` lists keyed by lavaan group labels
    - reliability computation under user-supplied priors
    - propagating sampling uncertainty from an externally estimated prior covariance
    - extending `prior_mean` / `prior_cov` to merMod
    - extending `augment_lav_predict()` / OpenMx-facing helpers
    - changing `compute_lav_fs_matrices()` for this feature

## Work State

### Completed

- Investigated `mirt::fscores()` as the API reference:
  - confirmed mirt version 1.46.1
  - confirmed `mean = NULL` and `cov = NULL` arguments
  - confirmed `mirt::fscores()` forwards them as `gmean = mean, gcov = cov`
  - confirmed documentation that `mean` and `cov` override fitted-object latent
    means/covariances when non-NULL.
- Renamed the planned R2spa arguments to `prior_mean` / `prior_cov` (per user decision;
  keeps `mirt::fscores(mean, cov)` as the conceptual reference while disambiguating from
  lavaan group-specific estimates).
- Read relevant R2spa sources and tests:
  - `/home/marklai/R2spa/R/get_fscore.R`
  - `/home/marklai/R2spa/R/get_fs_methods.R`
  - `/home/marklai/R2spa/R/get_fscore_math.R`
  - `/home/marklai/R2spa/R/tspa.R`
  - `/home/marklai/R2spa/R/tspa_corrected_se.R`
  - `/home/marklai/R2spa/R/get_fs_int.R`
  - `/home/marklai/R2spa/tests/testthat/test-get_fscore.R`
  - `/home/marklai/R2spa/tests/testthat/test-compute_fscore.R`
  - `/home/marklai/R2spa/tests/testthat/test-lavPredict_equivalence.R`
- Confirmed the core mathematical feasibility:
  - `compute_fscore()` already accepts separate `psi` and `alpha`.
  - User-supplied priors can be injected in `get_fs_blocks.lavaan()` / `prepare_fs()`
    before calling `compute_fscore()`.
- Confirmed current data-frame and lavaan method behavior:
  - `get_fs.data.frame()` passes `...` to `lavaan::cfa()`.
  - `get_fs.lavaan()` currently ignores `...`, so explicit `prior_mean` / `prior_cov`
    arguments can be added without leaking into `cfa()`.
- Checked repository state:
  - `git status --short` showed a clean worktree.
  - recent commits include `00bf670`, `9c60ff8`, `5f72883`, `db36781`, and `258b673`.
- Asked scope/tradeoff questions and received user answers:
  - shared prior only
  - support `corrected_fsT` / `vfsLT`
  - `reliability = TRUE` should error initially with custom priors
  - new argument names are `prior_mean` / `prior_cov`
- Produced a complete implementation plan covering API, math wiring, validation,
  restrictions, docs, tests, verification lifecycle, and out-of-scope items.

### Completed — implementation (2026-08-16, archived as-is)

- API: `prior_mean = NULL, prior_cov = NULL` added to `get_fs.data.frame()` (before
  `...`, so not forwarded to `lavaan::cfa()`) and `get_fs.lavaan()` in
  `R/get_fs_methods.R`, and to the `get_fs_lavaan()` wrapper in `R/get_fscore.R`.
- Validation: new `validate_fs_priors()` helper in `R/get_fs_methods.R`: finite
  values, q-length (named inputs reordered to model order; unnamed used in model
  order), `prior_cov` finite/square/symmetric/positive definite (scalar or 1 x 1
  accepted when q = 1).
- Method restrictions: priors error for `method = "Bartlett"`/`"ML"` and for
  `reliability = TRUE`; `corrected_fsT` / `vfsLT` accepted (prior covariance treated
  as fixed).
- Math wiring: `psi_override = NULL` threaded through `compute_fspars()` (sets the
  group's `mat$psi` after estimated-coefficient replacement), `compute_a()`,
  `compute_evfs()`, `compute_ldfs()`, `compute_grad_ld_evfs()`, `vcov_ld_evfs()`,
  `correct_evfs()` in `R/get_fscore_math.R`.
- Scoring: `get_fs_blocks.lavaan()` accepts the priors; its `prepare_fs()` helper
  uses `psi_use <- if (is.null(prior_cov)) est$psi else prior_cov` and
  `alpha_use <- if (is.null(prior_mean)) est$alpha else prior_mean` for complete
  data, each missing-data pattern, and every group (shared prior). All derived
  outputs (`fs`, `fsT`, `fsL`, `fsb`, `scoring_matrix`, `se_*`, `ev_*`, `ecov_*`)
  are recomputed from the prior-based scoring matrix.
- merMod: `get_fs.merMod()` rejects non-NULL `prior_mean`/`prior_cov` passed through
  `...` (covers the `get_fs_lmer()` path too).
- Docs: `@param prior_mean` / `@param prior_cov` + multi-group example in the
  `get_fs()` Roxygen block (`R/get_fscore.R`); `man/get_fs.Rd` and
  `man/get_fs_lavaan.Rd` regenerated via `devtools::document()` (idempotent re-run).
- Tests: `tests/testthat/test-get_fs_priors.R` — 22 `test_that` blocks covering
  manual-`compute_fscore()` equivalence (both priors / one prior each way), attribute
  consistency, name reordering, all planned validation errors, `corrected_fsT` and
  `vfsLT` with fixed priors, `vcov_corrected()` on a prior-adjusted `tspa()` fit,
  data.frame/matrix/`get_fs_lavaan()` entry points, q = 1 scalar forms, multi-group
  (per-group equivalence, unified vs list, identical-group-data invariance), shared
  prior with missing data, and merMod rejection.
- Verification (per AGENTS.md lifecycle, 2026-08-16):
  - `devtools::load_all()` — OK.
  - `devtools::document()` — OK, idempotent (second run changes nothing).
  - `devtools::test()` — **658 pass, 0 fail, 0 warn, 0 skip** (596 at the PLAN 03
    [mermod] close; +62 expectations from the new file, verified by A/B run).
  - `devtools::check()` (as-cran default) — **0 errors, 2 WARNINGs, 3 NOTEs**:
    identical to the pre-existing baseline (S3-consistency for `get_fs.lavaan()` /
    `get_fs.merMod()`; undocumented `fsm` / `...` Rd args; `.lintr` hidden file,
    unused `Matrix` Import, top-level files). New roxygen example runs
    (`checking examples ... OK`); all 14 vignettes rebuild; tests OK under check.

### Active

- (none — plan complete; archived 2026-08-16)

### Blocked

- (none)

## Relevant Files

- `/home/marklai/R2spa/R/get_fscore.R`: user-facing `get_fs()` documentation,
  `get_fs_lavaan()` wrapper, and shared Roxygen block where `@param prior_mean` /
  `@param prior_cov` should be added.
- `/home/marklai/R2spa/R/get_fs_methods.R`: contains `get_fs.data.frame()`,
  `get_fs.lavaan()`, `get_fs.merMod()`, and `get_fs_blocks.lavaan()`; main API injection
  point and where prior validation/wiring belongs.
- `/home/marklai/R2spa/R/get_fscore_math.R`: contains `compute_fscore()`, scoring-matrix
  code, `correct_evfs()`, `vcov_ld_evfs()`, and related internal helpers; needs
  `psi_override` support for corrected-SE/vcov compatibility.
- `/home/marklai/R2spa/R/tspa.R`: downstream consumer of `fs`, `fsT`, `fsL`, and `fsb`;
  used to verify that prior-adjusted scores remain compatible with two-stage path
  analysis.
- `/home/marklai/R2spa/R/tspa_corrected_se.R`: `vcov_corrected()` consumes `fsL`, `fsT`,
  and `vfsLT`; relevant for testing `vfsLT = TRUE` with a fixed custom covariance.
- `/home/marklai/R2spa/R/get_fs_int.R`: downstream interaction-indicator function that
  consumes factor scores and SEs; relevant for compatibility verification.
- `/home/marklai/R2spa/tests/testthat/test-get_fscore.R`: existing `get_fs()` tests;
  useful reference for expected output structure and backward compatibility.
- `/home/marklai/R2spa/tests/testthat/test-compute_fscore.R`: existing
  `compute_fscore()` tests; useful reference for manual comparisons with overridden
  `alpha` / `psi`.
- `/home/marklai/R2spa/tests/testthat/test-lavPredict_equivalence.R`: existing
  equivalence tests around `lavPredict()` and fs matrices; useful baseline for ensuring
  no unintended changes.
- `/home/marklai/R2spa/tests/testthat/test-get_fs_priors.R`: planned new test file for
  custom `prior_mean` / `prior_cov` behavior.
- `/home/marklai/R2spa/man/get_fs.Rd`: generated documentation target for Roxygen
  changes; should not be edited manually.
- `/home/marklai/R2spa/NAMESPACE`: likely unchanged because no new exported functions are
  planned.
