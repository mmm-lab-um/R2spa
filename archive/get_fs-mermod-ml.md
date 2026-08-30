# PLAN: `method = "ML"` for `get_fs.merMod()`

Status: **DONE** (2026-08-18; implemented per this plan, committed as `aa5be3e` on
`refactor/core`; suite 1211 pass at that point). Follow-ups (2026-08-18): first-term
Zden fold perf change + multi-term `clen` restore in the EB branch landed in `b730b53`
(suite 1415); the pre-existing lme4-2.x `get_D()` convention work item
(`get_d-mermod-lme4-2x.md`) is DONE in `966ff9b` (bit-identical rework +
`test-lme4_compat.R` canary). Remaining follow-ups: STATUS #5 `correct_evfs()`
perf, F1–F4 deferrals. Approved decisions at bottom.
Sequencing: this is Plan A of the pair (see also `get_fs-lavaan-mean.md` = Plan B);
execute A first, sequentially, single agent.

## Context (for a fresh session)

- Repo `/home/marklai/R2spa`. Test baseline verified 2026-08-17: **707 pass / 0 fail /
  0 warn / 0 skip** (`devtools::test()`).
- Lifecycle (AGENTS.md, never reorder): `devtools::load_all()` → `devtools::document()`
  (roxygen changed) → `devtools::test()` → `devtools::check()`.
  Never hand-edit `NAMESPACE`/`man/*.Rd`; no `library()`/`require()` in function bodies;
  namespace external calls.
- merMod code path: `get_fs.merMod()` (`R/get_fs_methods.R:486-562`) →
  `get_fs_blocks.merMod()` (`R/get_fs_methods.R:373-460`). EB method only today.
  Legacy wrapper `get_fs_lmer()` (`R/get_fscore.R:498-513`).

## User request

Support `method = "ML"` for merMod: the prior-free (ML) estimate of u0, u1, ... per
cluster — OLS of the cluster's residuals on its random-effects design block —
analogous to Bartlett scores for lavaan (estimator that does not use the prior D).

## Verified math (R 4.6.1, session 2026-08-17)

Reference fit: `lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)`.
For cluster j with row index set `{j}`: `r_j = y_j − X_j β`, `Z_j` = the cluster's
rows of `lme4::getME(object, "Z")` sliced to the first RE term's `num_re` columns,
`Kz_j = Z_j' Z_j` (p × p, p = `num_re`):

- **Score** `û_ml,j = (Z'Z)⁺ Z' r_j` — the ML estimate of u_j with no prior on u
  (≡ OLS of `r_j` on `Z_j`). Verified identical (to ~1e-15) to hand `lm(r_j ~ Z_j)`
  and to `solve(Kz) %*% crossprod(zj, r)`; e.g. subject 308 → (−7.2125, 11.2965)
  vs ranef EB (2.2591, 9.1994) — correctly different.
- **fsL_j = I_p**: from `r_j = Z_j u_j + e_j`, score = `u_j + Kz⁺ Z_j' e_j`.
- **fsT_j = σ² Kz_j⁺**: Var(score error) = σ² Kz⁺ Z'Z Kz⁺ = σ² Kz⁺ (Penrose
  condition A⁺ A A⁺ = A⁺). Verified diag: single-RE case gives se = σ/√n_j
  (25.5918/√10 = 8.0928 for the 10-obs cluster).
- **fsb = NULL** (as in EB; E[u] = 0 is the model-implied center).
- Prior covariance D is **not used** — that is the defining property of "ML" here
  (`get_D()` stays EB-only).

Rank deficiency (RESOLVED, user-confirmed): use `MASS::ginv(Kz_j)` always — exact at
full rank, minimum-norm solution when singular; matches `compute_a_reg()`'s
unconditional `MASS::ginv`; no separate error path. `fsT_j = s^2 * ginv(Kz_j)` holds
via the Penrose identity in both cases.

## Design

- `get_fs.merMod(object, method = c("EB", "ML"), ...)`; pass `method` into
  `get_fs_blocks.merMod(object, method = c("EB", "ML"), legacy_names = FALSE)`.
  `match.arg` at the top of `get_fs.merMod()`. The per-cluster loop branches:
  EB keeps the existing code untouched; ML adds the OLS branch.
- ML residuals: `y <- model.response(model.frame(object))`,
  `X <- as.matrix(lme4::getME(object, "X"))`, `beta <- lme4::fixef(object)`
  (same access pattern the existing tests use at
  `tests/testthat/test-get_fscore.R:431-436`). Per cluster `r_j = y[idx] − X[idx, ] %*% beta`.
- Block fields for ML: `fs = û_ml,j` (1 × p row, cols `re_names`),
  `fsL = diag(p)` (dimnames fs_names / re_names), `fsT = s^2 * ginv(Kz_j)`
  (dimnames fs_names), `fsb = NULL`, `scoring_matrix = ginv(Kz_j) %*% t(zj)`
  (p × n_j; reproduces the score via `S_j (y_j − X_j β)`, mirroring the EB
  contract documented in `get_fs()`'s `@return` and `vignettes/scoring-matrices.Rmd`).
- **Column names are unchanged from EB** (`fs_u0`, `fs_u0_se`, `u0_by_fs_u0`,
  `ev_`/`ecov_`; `legacy_names = TRUE` composes exactly as for EB, since the rename
  is off the `_eb` suffix in `re_names`). One row per cluster, in canonical level
  order — all existing invariants (`fsL_arr`/`fsT_arr` attributes, per-cluster
  `scoring_matrix` list) carry over without change.
- Gaussian lmer only — same implicit constraint as EB today (`stats::sigma()`).

### Defect found while planning (fix included)

`get_fs_lmer()` captures `method` as a formal (`R/get_fscore.R:500`) but **never
forwards it** — `get_fs_lmer(lme, method = "ML")` would silently return EB scores.
Fix: forward `method` into the `get_fs()` call; extend the doc from
"Currently only `EB`" to the real choices `c("EB", "ML")`. (`corrected_fsT`/`vfsLT`/
`fsm` remain documented no-ops for merMod — no change.)

## Files

- `R/get_fs_methods.R` — `get_fs.merMod()` (signature + dispatch),
  `get_fs_blocks.merMod()` (ML branch; any helper co-located).
- `R/get_fscore.R` — `get_fs_lmer()` method forwarding + `@param method`;
  `get_fs()` generic `@param method` paragraph gains the merMod sentence
  ("for `merMod` objects the choices are EB (empirical Bayes, default) and ML
  (prior-free per-cluster OLS)").
- `tests/testthat/test-get_fscore.R` — new tests in the merMod section.
- `man/` regenerated by `devtools::document()` (NAMESPACE unchanged — no new exports).

## Tests (new)

1. ML scores equal hand-per-cluster `(Z'Z)⁺ Z' r` (sleepstudy `lme1`;
   also cross-checked against `lm()` on each cluster's residuals; tol 1e-10).
 2. `fsL` is identity per cluster (check dimnames `fs_u*` / `u*`);
   `fsT` equals `σ² ginv(Kz_j)` per cluster (tol 1e-10; full rank, so ginv ≡ solve).
3. Scoring identity: `sm[[j]] %*% (y_j − X_j β)` reproduces the ML score (tol 1e-10).
4. `_se` columns = `sqrt(diag(fsT))` (single-RE model: `σ/√n_j`).
5. Column-name layout identical to EB defaults and to `legacy_names = TRUE`
   (ML composes with legacy names).
6. Single-RE model `(1 | Subject)` works (p = 1 path).
7. Unbalanced clusters work (drop rows as in the existing test at :484).
8. `get_fs(lme1, method = "regression")` → clean error listing the EB/ML choices.
9. `get_fs_lmer(lme1, method = "ML")` ≡ `get_fs(lme1, method = "ML")` (regression
   test for the forwarding fix); default (no method) unchanged ≡ EB/ranef.
10. Singular-`Kz` cluster (random slope on a within-cluster constant predictor)
    returns finite minimum-norm scores via ginv, no crash (approved fallback).

## Verification

- Lifecycle per AGENTS.md. Expect: baseline 707 + new tests, 0 fail.
- `devtools::check()` — expect the one pre-existing OpenMx NOTE only.

## Risks / notes

- No impact on lavaan path, `tspa()`, or quarantined code.
- The `scoring-matrices` vignette documents the EB identity; if it is built
  (check `vignettes/`), add one ML sentence or leave as-is (it is EB-specific —
  verify before deciding; do not regenerate cached fixtures).

## Subagent roles (dispatched by the lead agent)

The lead agent implements and owns the lifecycle; the specialists are dispatched
via the task tool (subagent types `r-architect`, `r-doc`, `r-tester`). Their
outputs are advisory; the lead reviews, applies, and runs the lifecycle. Subagents
do not hand-edit `man/`, do not run `document()`/`test()`/`check()` themselves, and
change nothing in `.quarantine/`.

1. **`r-architect`** — at plan start: brief with the verified math, the per-cluster
   block fields, the `get_fs_lmer()` forwarding fix, and the file:line map; ask it
   to (a) confirm the ML branch fits `get_fs_blocks.merMod()` without disturbing
   the EB path and the "one row per cluster" assembly exception
   (`R/get_fs_methods.R:509-516`), (b) confirm the attribute shapes
   (`fsL_arr`/`fsT_arr`, per-cluster `scoring_matrix`) flow through
   `augment_fs()` and `tspa()` unchanged, (c) sanity-check the ginv-always choice
   and the residual access pattern (`model.frame`/`getME("X")`/`fixef`).
   Re-dispatch before `check()` for a diff review (naming, cross-file consistency).
2. **`r-tester`** — after the code is drafted: implement/extend the tests above per
   `tests/testthat/` conventions (testthat edition 3; independent references via
   `lme4::getME()`, `model.frame()`, `lm()` on cluster residuals;
   `expect_equal(..., tolerance = 1e-10)`). Invite it to add edge cases: factor
   level ordering, `Z ≠ X`, unbalanced clusters, p = 1, `legacy_names` composition.
3. **`r-doc`** — after code + tests, before `devtools::document()`: produce the
   exact roxygen wording for (a) `get_fs()`'s `@param method` (new merMod sentence:
   choices EB (empirical Bayes, default) and ML (prior-free per-cluster OLS));
   (b) the corrected `get_fs_lmer()` `@param method` (currently documented
   "only `EB`" while silently no-op; after the fix it forwards EB/ML);
   (c) any `@return`/`@details` deltas (identity `fsL`, per-cluster `σ²(Z'Z)⁻¹`
   `fsT`, ML `scoring_matrix` identity). The lead applies the text verbatim to
   `R/*.R`, then runs `document()`.

Order: `r-architect` → (lead code in parallel with `r-tester` drafting tests) →
`r-doc` → lead runs `load_all() → document() → test() → check()`.

## Approved decisions (2026-08-17)

- Method string: `"ML"`; merMod choices become `c("EB", "ML")`.
- Rank-deficient `Kz_j`: `MASS::ginv` always (minimum-norm fallback; matches
  `compute_a_reg()`'s ginv usage).
- Column names/attributes unchanged from EB; ML composes with `legacy_names`.
- `get_fs_lmer()` forwarding defect is fixed as part of this plan.
- Execution: sequential with Plan B (this plan first), single lead agent, with the
  subagent dispatch above.
