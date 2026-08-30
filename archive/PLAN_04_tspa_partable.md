# Plan: Robust redesign of `tspa()` stage-2 model construction (partable- and `sam()`-informed, syntax-rendering)

Replace the fragile syntax-string *appending* in `tspa()`'s stage-2 assemblers
(`tspa_sf` :209-248 / `tspa_mf` :250-323 in `R/tspa.R`) with **structured, table-driven
model construction**, informed by lavaan 0.7-2's parameter-table (partable)
machinery and `lavaan::sam()` ("Structural After Measurement").

**Explicit design goal (user-directed): contain the maintenance cost of future
partable format changes.** The format has already drifted once (0.6.x → 0.7-2),
the current release is 3 weeks old, and the package declares no lavaan version
bound. The design therefore keeps partable *off* the estimation path
(decision 2 below) and confines every lavaan-internal touch point to one
compat module plus canary tests, so that a future format change surfaces as
one localized test failure with a table diff — not as a
user-discovered estimation bug.

Draft/review copy of this plan: `.opencode/plans/tspa-partable-sam.md`.

## Progress

- [x] **Phase 1 — compat module + canaries (no behavior change)** — done
      2026-08-16: `R/lavaan_compat.R` (contract table, `tsp_layout()` +
      `tsp_resolve_layout()`, `tsp_partable_read()`, `tsp_model_matrices()`,
      `tsp_free_matrices()`, `tsp_partable_mats()`, `tsp_nobs()`,
      `tsp_ngroups()`, `tsp_norig()`, `tsp_layout_reset()`); 5
      grandStandardizedSolution sites + `tspa_corrected_se.R` @Data@ngroups
      site migrated; `tests/testthat/test-lavaan_compat.R` (golden canary +
      wrapper A/B + layout error/memoization — 40 expectations). Suite:
      **698 pass, 0 fail** (658 baseline + 40); `document()` idempotent,
      NAMESPACE/man unchanged; scoped `R CMD check`: only the pre-existing
      2 WARNING / 2 NOTE items. Drift finding: 0.7-2 partable `free` column
      is the 1-based *position* in the free estimate vector (0 = fixed),
      not a 0/1 flag; `partable()` and `lavInspect(what="list")` share
      identical values/modes (primary = the latter, the legacy view).
- [x] **Phase 2 — schema + renderer + cutover** — done
      2026-08-16 (cont.): `R/tspa.R` — schema constructors
      (`tspa_schema_sf()`, `tspa_schema_mf()`; 8 columns
      `lhs, op, rhs, value, free, group, label, kind`, one row per
      (statement term, group), user model carried as a single verbatim
      `raw` row), statement assembler `tspa_statements()`, and
      `tspa_render(sch, style)` reproducing the legacy strings
      **character-for-character** (per-path spacing quirks, bare-value
      SG intercepts, blank line before `# structural model` when no
      intercepts). `tspa_sf`/`tspa_mf` bodies cut over onto schema +
      render; verbatim legacy copies kept as `tspa_sf_legacy`/
      `tspa_mf_legacy` behind `tspa_env$render` (default `"schema"`) for
      A/B until Phase 3. Product-score auto-aliasing `tspa_sf_alias()`
      (`fs_a:fs_b` → `fs_v` in data + `se_fs` names; >1 candidate →
      clear error; `:` is invalid in lavaan variable names). A/B:
      **10/10 case strings byte-identical** (SF SG/MG/3fac/comments,
      MF SG list/growth/plain-matrix, MF MG 2fac/3fac); interaction
      alias end-to-end verified (auto vs manual rename → identical
      coef/vcov/tspaModel). New `tests/testthat/test-tspa_render.R`
      (schema rows, render byte-identity incl. trailing-newline model,
      alias incl. ambiguity error, contract/attrs, 67 expectations).
      Suite: **765 pass, 0 fail**; `document()` idempotent (NAMESPACE
      unchanged); scoped `R CMD check`: same 2 WARNING / 2 NOTE as
      baseline, tests OK under check.
- [x] **Phase 3 — hardening & docs** — done
      2026-08-17 (cont.): CI lavaan axis in
      `.github/workflows/R-CMD-check.yaml` (`lavaan-baseline` job:
      0.7-2 pinned from the CRAN archive + full check;
      `lavaan-dev-tripwire` job: dev lavaan from
      `Lavaan/lavaan@master` + only the scoped subset
      `test-tspa*.R` / `test-lavaan_compat.R` /
      `test-get_fs_int.R`); controlled re-knit of `vignettes/
      R2spa.Rmd` + `vignettes/multiple-factors.Rmd` — verification-only,
      all 5 printed `tspaModel` blocks byte-identical (checked in the knit
      artifacts); `vignettes/get_fs_int-vignette.Rmd` manual rename
      workaround removed — auto-alias reproduces the old fits
      bit-identically (model/coef/vcov/standardizedSolution, both
      data paths). Migration fallback removed (`tspa_env`,
      `tspa_sf_legacy`, `tspa_mf_legacy`); the A/B format guarantee is
      frozen in `test-tspa_render.R` as in-test pinned reference builders
      (`ref_sf()`/`ref_mf()`) + a `lavaanify()` parameter row-order
      golden. Suite: **763 pass, 0 fail**; `document()` idempotent.
      Full `devtools::check()` in the final (fallback-free) state:
      **0 errors, 2 WARNING / 3 NOTE** = the pre-existing baseline
      items, no new findings; all 16 vignettes rebuild; tests OK.
      STATUS.md updated; plan archived as
      `archive/PLAN_04_tspa_partable.md`.
      Deviation note: the draft's "full check before fallback removal"
      ordering was collapsed into a single full check of the final
      fallback-free state — the fallback was A/B safety net only,
      already satisfied in Phase 2 (10/10 byte-identity), so checking
      the dead path added no information.

### Implementation notes (as they diverge from the draft)

- **Baseline count**: the "596 pass" baseline in this plan predates the PLAN 03
  (fs priors) commit (`32fc817`). The true pre-Phase-1 suite is **658 pass,
  0 fail** (2026-08-16, R 4.6.1, lavaan 0.7-2). Phase gates use "identical
  numbers to the baseline at phase start", i.e. 658.
- **Renderer output**: the renderer reproduces the legacy string
  **character-for-character** (legacy statement order: measurement →
  errors (var/cov interleaved in lower-tri column-major order) → intercepts →
  user structural block). Rationale: the Phase 2 A/B gate requires
  "estimates, vcov, AND parameter row order identical" vs the legacy
  string-append model; lavaan's partable/parameter ordering follows statement
  order (auto-added free parameters last), so any reorder (e.g. "user block
  first") would move `~` rows in `lavaanify()`/`parameterestimates()` and
  break the row-order gate. The §3 "user block verbatim → error-var →
  error-cov → intercept" ordering is kept as the **schema row organization**
  (user rows first in the table, `kind`-marked), not the rendered order.
  Consequence: the Phase 3 vignette re-knit is verification-only.
- **Canonical partable df** (`tsp_partable_read()`): the plan's 10 columns
  (`lhs, op, rhs, value, free, group, block, label, user, ustart`) plus a
  trailing `exo` column — `grand_standardized_solution()`'s returned frame
  must keep its current `exo` column (user-visible; existing tests consume it).
- **`lavTech(what = "est"/"free")` shape** verified on 0.7-2: a **flat** list
  of 6 matrices (`lambda, theta, psi, beta, nu, alpha`) repeated per group
  (not nested). `lavInspect(what = "est")` is nested per group with
  `group.block` element names — the compat fallback flattens to the lavTech
  shape by reassigning per-group block names.

## Scope

1. New internal compat module `R/lavaan_compat.R` — the only file allowed to
   read lavaan internals (partable/list views, `@Data@` slots, `lavTech` views).
2. R2spa-owned stage-2 model schema (frozen by R2spa, not by lavaan) +
   `tspa_render()` renderer to a lavaan model **syntax string** (decisions 2).
3. Cutover of `tspa_sf` / `tspa_mf` internals onto schema + renderer, with the
   legacy string builder kept as an internal fallback until Phase 3 check.
4. Canary/anchor tests (golden partable, golden estimates, renderer A/B,
   `tspa()`-contract tests).
5. CI lavaan version matrix (0.7-2 baseline + dev snapshot, decision 3).
6. Vignette re-knit of the two files that print `tspaModel`
   (`vignettes/R2spa.Rmd`, `vignettes/multiple-factors.Rmd`) and the
   `get_fs_int` vignette (interaction-rename workaround removed by an
   automatic alias).

Out of scope: `DESCRIPTION` changes (decision 1); any `sam()` integration or
partable-write path (decisions 2, 4); STATUS follow-up F1 (tspa native
unified/list score shapes) — deferred to its own plan (decision 5).

## Decisions (2026-08-16, user)

1. **Version contract: NO change to `DESCRIPTION`** — `lavaan` stays unbounded.
   Drift defense rests on: the compat module's dependency-contract table, its
   "layout not supported" error (must name the tested-up-to version), the
   canary tests, and the CI dev slot. The code already de-facto requires
   ≥ 0.7-2 (uses `lavTech`); that floor stays undeclared.
2. **Partable in the estimation path: SYNTAX RENDERING ONLY.** Partable is
   consumed read-only by the compat module and used in canary/contract tests.
   No partable renderer, no `tsp_partable_write()`.
3. **CI matrix: 0.7-2 (pinned baseline) + latest dev** (lavaan website source
   / r-universe). Non-baseline slot runs only the scoped test subset
   (`test-tspa*.R` + compat tests). No 0.6.x slot — current code does not work
   on 0.6.x anyway (`lavTech` is a 0.7-x export), so there is no compatibility
   to preserve.
4. **`sam()` spike: SKIPPED.** `sam()` stays reference material (Evidence
   below); no phase depends on it. If partable-native construction is ever
   reconsidered, `sam(model = <parameter table>, se = "twostep")` and the
   `lav_partable_*` toolkit are the entry points to re-verify.
5. **F1 (STATUS follow-up, tspa native unified/list shapes): DEFERRED** to its
   own plan. The schema below (`group` column, list-valued `fsT` ⇒ per-group
   rows) makes it a small follow-up.

## Evidence (live probes, lavaan 0.7-2, 2026-08-16)

- **Environment**: R 4.6.1; **lavaan 0.7-2, packaged 2026-07-16** (~3 weeks
  old). Installed build has **no `NEWS.md`** and ships
  `understanding_lavaan_internals.R` → major restructure release; the dev line
  is active, so "changes partable format again" is a realistic near-term risk.
- **Partable drift is proven, not hypothetical**: 0.6.x partable had a `fix`
  column; 0.7-2 dropped it. In 0.7-2, `partable(fit)` and
  `lavInspect(fit, what = "list")` (class `lavaan.data.frame`) share columns:
  `id, lhs, op, rhs, user, block, group, free, ustart, exo, label, plabel, start, est, se`
  — with `user`, `block`, `plabel`, `est`, `se` new vs 0.6.x.
- `partable(object)` takes a single argument (no `which`/column filter).
- New 0.7-x accessor family:
  `lavTech(object, what = "free", add_labels, add_class, list.by.group,
  drop.list_single_group, ...)`. R2spa already uses it
  (`R/grandStandardizedSolution.R:84,96,101` + `tests/testthat/
  test-grandStandardizedSolution.R`) with `what = "est"/"free"/"partable"`.
- `lavInspect()` semantics changed in 0.7-2: `what = "est"/"free"` return
  `lavaan.list` of **matrices** (lambda/theta/psi/beta), no `value` argument
  (positional name-based access errors); `what = "list"` returns the
  `lavaan.data.frame` partable view; `drop.list.single.group` controls list
  flattening for single-group fits.
- **`coef(fit)` remains a named vector and `vcov()` is name-indexed** in 0.7-2
  → name-based point/SE access (`coef(fit)["x3~~x4"]`) is the stable anchor
  for tests (no positional row subsetting).
- **`lavaan::sam()`** (new):
  `sam(model, data, aux, cmd = "sem", se = "twostep", mm_list,
  mm_args = list(bounds = "wide.zerovar"), struc_args =
  list(estimator = "ML"), sam_method = "local", local_options =
  list(M.method = "ML", lambda.correction = TRUE, alpha.correction = 0L,
  twolevel.method = "h1"), global_options, bootstrap =
  list(R = 1000L, type = "ordinary", show.progress = FALSE),
  output = "lavaan")`. **`model` accepts a parameter table (e.g. output of
  `lavParTable()`) in addition to lavaan syntax**; `se = "twostep"` (default)
  accounts for stage-1 measurement uncertainty.
- **`lav_partable_*` toolkit** (12 exports, 0.7-2; arguments not captured —
  not needed per decision 4): `lav_partable_df, _merge, _attributes,
  _constraints_ceq, _constraints_ciq, _constraints_def, _labels,
  _unrestricted, _from_lm, _independence, _add, _complete, _npar, _ndat`.
- `startvalues()` / `lavPrep` are **not exported** in 0.7-2 (pre-fit parTable
  access, if ever needed, via `lavModelInfo(...)$parTable` / `lavParTable`).
- **R2spa's current lavaan-internal surface** (grep-verified; line numbers
  pre-PLAN 03-landing for `get_fs_methods.R`/`get_fscore_math.R`):
  - `R/grandStandardizedSolution.R` — `lavTech(what="est")` (:84),
    `unlist(object@Data@nobs)` (:87),
    `lavInspect(what="list")` + `subset(..., op == "~")` →
    `[, c("lhs","op","rhs","exo","group","block","label")]` (:98-100) —
    **R2spa already consumes a partable-equivalent table for reading**,
    `lavTech(what="free")` (:96), `lavTech(what="partable",
    list.by.group = TRUE)` → per-group `$beta` position vectors (:101-114).
   - `R/tspa_corrected_se.R` — `tspa_fit@Data@ngroups` slot (:23, per commit
     `5f72883`); `update_tspa()` `eval(call)` re-fit contract :69-80;
     `coef(update_tspa(...))` Jacobian path :41-65.
  - `R/get_fs_methods.R` — `lavInspect(what="est"/"data")` in
    `get_fs_blocks.lavaan()`; `unlist(object@Data@norig)` slot;
    `lavInspect(what="est")` in `get_fs.lavaan()` (reliability/priors paths).
  - `R/get_fscore_math.R` — `lavInspect(what="ngroups"/"est"
    (drop.list.single.group = FALSE)/"meanstructure"/"free"/"implied"/
    "norig")` across `augment_lav_predict()`, `compute_fspars()`,
    `correct_evfs()`, and the reliability helpers.
- **No version gating exists anywhere** (no `packageVersion()` in `R/`, no
  bound in `DESCRIPTION`) — existing robustness posture is direct slot access
  (commit `5f72883`, PLAN 02 Step 5) with no declared contract.
- **Baseline green on 0.7-2**: full suite **596 pass** (STATUS.md,
  2026-08-16); `test(filter = "tspa")` 39 pass.
- **Repo state at plan time**: PLAN 03 (user priors for `get_fs()`) is
  mid-implementation, **uncommitted**: `R/get_fs_methods.R`, `R/get_fscore.R`,
  `R/get_fscore_math.R`, regenerated `man/get_fs.Rd` + `man/get_fs_lavaan.Rd`,
  new `tests/testthat/test-get_fs_priors.R`, plus `_PLAN_03.md`. This plan
  touches almost none of those files — **PLAN 03 must land first** so the
  diffs stay disentangled (STATUS.md already flags that its commit scope
  "must be decided explicitly").
- `get_fscore_math.R`/`tspa.R` sensitive numerics and the `vcov_corrected()`
  re-fit path are the two places where a change must be provably inert
  (golden estimates gate, below).

## Maintenance cost: partable format drift (the user question)

Three candidate architectures, costed against "lavaan changes partable format
again" (already happened once, 0.6.x → 0.7-2, within one release window; the
current release is 3 weeks old and dev-active):

| Option | Drift coupling on the ESTIMATION path | Correctness (name collisions, MG per-group values, user labels) | One-time effort | Ongoing cost per future lavaan release |
|---|---|---|---|---|
| **A. Keep string appending** (status quo) | Low — lavaan model *syntax* is the most stable lavaan API (decades) | Weak — the original motivation of this plan: `fs_` name collisions, ambiguous `# constrain the errors` blocks, MG per-group value routing | 0 | Low, but breakage = ad-hoc bugfixes discovered by users |
| **B. Raw partable as the construction path** (build the model as a partable df for `sem()`/`sam()`) | **High** — estimation coupled to a schema that changed once already, on a 3-week-old release whose `lav_partable_*` toolkit is early-lifecycle | Strong — structured rows; 0.7-2 `user` column distinguishes user-written rows from program-added rows | Medium | **Highest**: every partable column/semantics change is an R2spa release; failure mode can be *silent mis-estimation* between drift and detection (malformed writes are validated lazily by the parser) |
| **C. Own schema → render to stable syntax; partable confined to a compat module + tests** — CHOSEN (decision 2) | Low-medium — the only format coupling lives in ONE compat module (plus tests); the estimation path depends on lavaan *syntax* (stable) + R2spa's own frozen schema | Strong — injected rows carried in an R2spa-owned table with a `kind` column; user rows verbatim; generated label namespace prevents collisions | Medium-high (~200–300 LOC compat + renderer + tests) | **Low and bounded**: next drift first surfaces as a *single canary test failure* with a column diff; fix is localized to the compat module; estimation behavior provably unchanged by A/B golden-estimate tests |

Cost model for "lavaan changes partable again" under the chosen design (C):

1. **Detection**: the canary test
   `tsp_partable_read(canonical_fit) == expected literal table` fails loudly
   with a column diff in CI (the dev slot runs a current lavaan) — before any
   user sees it.
2. **Fix scope**: column-name / fixed-value mapping inside the compat module
   only (e.g. fixed-value detection `free == 0` / `start` non-NA now vs
   legacy `fix`; `user` column present vs inferred). No changes to schema,
   renderer, `tspa()` contract, or user-visible output.
3. **Validation**: golden-estimate tests (point estimates + vcov on
   canonical SG/MG models, name-addressed via `coef()`) confirm estimates did
   not move.
4. **Estimated effort**: one focused PR, typically < 1 day, plus CI re-run.

Under option B the same event can additionally corrupt the *construction*
path (wrong columns written, constraints silently dropped): detection depends
on estimate-level tests and the fix can span `tspa.R` — strictly more
expensive, with a wider silent-failure window. Option A avoids drift entirely
but keeps the correctness bugs that motivate this plan.

Costs the chosen design **always** pays (priced once): ~2× CI check time on
the dev slot (OpenMx compile dominant; mitigated by the scoped test subset), a
"lavaan internals we depend on" table kept current at the top of the compat
module, and the vignette re-knit in Phase 3.

## Design

### 1. Compat module — `R/lavaan_compat.R` (new; internal; unexported)

The ONLY file that reads lavaan internals. Header doc carries the
**"Lavaan dependency contract" table** — every internal consumed (partable/
list columns with meaning, `@Data@` slots, `lavTech`/`lavInspect` view
shapes), each annotated "verified against lavaan 0.7-2 (2026-07-16)". This
table is the audit artifact for the next drift (decision 1 makes it the
declared contract in place of a `DESCRIPTION` bound).

Functions (all internal, `snake_case`, co-located here per the single-boundary
rule):

- `tsp_layout()` — lazy, memoized probe of one canonical fit's
  `partable()`/`lavInspect(what = "list")`: resolves the fixed-indicator
  column (`free` now; `fix` historically), the fixed-value column (`start`
  now; `fix` historically), the user-row flag (`user` column when present,
  else inferred from `block`), and group/label columns. Unknown shape →
  **explicit error naming the installed version and the tested-up-to
  version** ("lavaan X.Y.Z partable layout not supported; R2spa tested up to
  0.7-2") — never silent guessing.
- `tsp_partable_read(fit)` → canonical df
  `[lhs, op, rhs, value, free, group, block, label, user, ustart]`.
  **Read-only by design** (no `tsp_partable_write()` — decision 2); used by
  canary/contract tests and by `grand_standardized_solution()`'s row
  identification.
- `tsp_model_matrices(fit)` → per-group lambda/theta/psi/beta/alpha lists
  (wraps `lavTech(what = "est")` with `lavInspect` fallback).
- `tsp_free_matrices(fit)` (wraps `lavTech(what = "free")`).
- `tsp_partable_mats(fit)` (wraps `lavTech(what = "partable",
  list.by.group = TRUE)`) — the per-group `$beta` position vectors
  `grand_standardized_solution()` currently extracts at :101-114.
- `tsp_nobs(fit)`, `tsp_ngroups(fit)`, `tsp_norig(fit)` — direct slot access
  with `lavInspect` fallback; consolidates the ad-hoc slot sites from commit
  `5f72883`.

**Capability probing, not `packageVersion()` string gating** (survives
patch/dev releases; the installed build has no NEWS and is a fast-moving dev
line).

Phase-1 migration (only sites that touch drifting surfaces):
`R/grandStandardizedSolution.R` (:84/:87/:96/:98/:101 → the five wrappers),
`R/tspa_corrected_se.R:23` → `tsp_ngroups()`. The `get_fs_methods.R` /
`get_fscore_math.R` sites (`what = "est"/"data"/"free"/"meanstructure"/
"implied"`) are migrated **only where the corresponding wrapper exists**
(`tsp_model_matrices`/`tsp_free_matrices`); the rest (`"data"`,
`"meanstructure"`, `"implied"`) stay on plain `lavInspect()` and are recorded
in the contract table as known coupling — deliberately deferred to keep
Phase 1 reviewable and the numerics-sensitive math layer untouched.

### 2. R2spa-owned stage-2 model schema (frozen by R2spa, not by lavaan)

Data.frame, one row per model statement:

```
lhs | op | rhs | value | free | group | label | kind
```

- `kind ∈ {user, struct, error_var, error_cov, intercept}` — injected rows
  are marked; the user `model` string's statements are carried **verbatim**
  (`kind = "user"`/`"struct"`) — no re-parsing of user syntax, and user
  rows can never be mistaken for (or collide with) injected ones.
- `group = NULL` ⇒ applies to all groups; multigroup construction emits
   per-group rows from the per-group `fsT`/`fsL`/`fsb` lists (the current
  `tspa_mf()` multigroup value routing, `R/tspa.R:148-185` + :255-311, made
  explicit and unit-testable).
- Construction inputs: `model` (user string), `fsT`, `fsL`, `fsb`, `se_fs`
  (single-factor path). This replaces the inline string building in
  `tspa_sf`/`tspa_mf` (currently `paste0("# constrain the errors\n",
  ev_lhs, " ~~ ", errors, " * ", ev_rhs)` at `R/tspa.R:301-302` and
  `paste0("# constrain the intercepts\n", fs, " ~ ", intercepts, " * 1")` at
  :307-308; the single-factor path additionally builds the
  `var =~ c(...) * fs_` / `fs_ ~~ c(...) * fs_` blocks at :221-244 with
  per-group `c(...)` value vectors that the renderer must reproduce
  verbatim).
- **Interaction** (`get_fs_int`): product scores keep `:`-separated *data*
  column names (e.g. `fs_x:fs_m`); the schema aliases them to generated model
  names — removes the manual rename workaround shown in
  `vignettes/get_fs_int-vignette.Rmd` (assigning the `fs_x:fs_m` column to a
  plain name such as `fs_xm`); old-spelling user models keep working because
  the alias is automatic, not mandatory.

### 3. Renderer — `tspa_render(schema_df)` (R/tspa.R, internal; the only
renderer, decision 2)

- Deterministic statement order: user block verbatim → error-variance rows →
  error-covariance rows → intercept rows; multigroup statements use lavaan's
  `group(` qualification exactly as the current `tspa_mf()` does (verified in
  the Phase 2 A/B step, not re-invented here).
- Fixed values per row via `*` syntax — the same patterns the current code
  already emits: `fs_a ~~ <v>*fs_a` (error variance), `fs_a ~ <v>*1`
  (intercept), `fs_a ~~ <v>*fs_b` (error covariance). The A/B equivalence
  test validates them end-to-end (no separate format spike — decision 4).
- **Generated label namespace** (e.g. `__r2spa_ev1__`) for every injected
  constraint, so user labels/variable names cannot collide with R2spa's
  constraints — closes the core correctness gap of option A without any
  partable dependency.
- The rendered string is attached as `attr(fit, "tspaModel")` exactly as
  today (two vignettes print it — Phase 3 re-knits them).
- No partable renderer, no `sem(model = <partable df>)` path (decision 2).

### 4. Canary / anchor tests

New files: `tests/testthat/test-lavaan_compat.R`, `tests/testthat/
test-tspa_render.R`.

- **Golden partable canary** (the drift detector):
  `tsp_partable_read(canonical_fit)` equals an expected literal table pinned
  in the test. Format change ⇒ this one test fails with a column diff ⇒
  compat-module-only fix (cost model step 2).
- **Golden estimates**: canonical SG/MG models (1 factor; 2 factors + error
  covariance) — point estimates, SEs, `vcov()` vs hand-calculated references
  (existing `test-tspa.R` pattern), rows addressed **by name** through
  `coef(fit)["x3~~x4"]` — no positional row subsetting (the 0.7-2
  `lavInspect` change makes positional access the fragile kind).
- **Renderer A/B**: `tspa_render(schema)` re-fit ≡ legacy string-append
  model — estimates, vcov, and parameter row order identical. Run during
  migration (Phase 2), kept as a permanent regression.
- **Contract tests**: `tspa()` attributes intact (`tspaModel` :200, `fsT`
  /`fsL` :202-203, `tspa_call` :205 — `tspa()` attaches no other attrs);
  `fsL`/`fsT` still re-callable (the `vcov_corrected()`
  `eval(call)` path, `R/tspa_corrected_se.R:69-80`); `vcov_corrected()`
  end-to-end on an MG fit; `grand_standardized_solution()` A/B vs
  pre-change output (its existing MG test already hand-calculates via
  `lavTech`).
- Existing suites untouched: `test-tspa.R`, `test-grandStandardizedSolution.
  R`, `test-get_fs_int.R`, priors tests must remain green through all phases.

### 5. CI version matrix (`.github/workflows/R-CMD-check.yaml`) — decision 3

- Lavaan axis: **0.7-2 (pinned baseline) + latest dev** (lavaan website
  source / r-universe). No 0.6.x slot (decision 3 rationale).
- Non-baseline slot runs only `tests/testthat/test-tspa*.R` +
  `test-lavaan_compat.R` + `test-get_fs_int.R` (scoped subset keeps the
  OpenMx-heavy check affordable at ~2× baseline time).
- The dev slot is the "lavaan changes partable again" tripwire: it fails on
  the canary test with the column diff (cost model step 1).

## `tspa()` change specifics (invariants)

- **Public API, dispatch, and attributes UNCHANGED.** In particular
  `fsL`/`fsT` remain named `tspa()` arguments and
  `attr(fit, "tspa_call")` (`R/tspa.R:205`) remains `match.call()` —
  `vcov_corrected()` / `update_tspa()`'s `eval(call)` re-fit contract
  (`R/tspa_corrected_se.R:69-80`) is bit-identical (the hardest constraint in
  the plan).
- `tspa_sf` / `tspa_mf` **bodies** replaced by schema construction +
  `tspa_render()`; function names and signatures stay (the test suite calls
  `tspa_mf()` directly — `test-tspa.R:387-422`, including the `lavaanify()`
  row/order tests).
- Cutover safety: legacy string builder kept as an internal
  `render = "string"` fallback inside `tspa()` until the Phase 3 full check
  passes; removed only at completion.
- Multigroup detection / validation (`R/tspa.R:148-185`, `tspa()` attribute
  validation from commit `5f72883`) unchanged.
- `tspa_mx.R` (OpenMx path) untouched.

## Vignettes

- **`vignettes/R2spa.Rmd`** (prints `tspaModel` at :62, :107, :161) and
  **`vignettes/multiple-factors.Rmd`** (:78, :146): the only files whose
  printed text changes when the renderer's statement order/labeling changes.
  **Controlled re-knit of these two only** in Phase 3 (AGENTS.md: cached
  `.RDS`/`.csv` fixtures back specific narratives — do not regenerate all
  vignettes). If `tspa_render()` output matches the legacy string
  character-for-character, re-knit is verification-only (cached fixtures
  should reproduce unchanged hashes where applicable).
- **`vignettes/get_fs_int-vignette.Rmd`**: replace the manual rename of the
  `fs_x:fs_m` product-score column with the automatic schema alias; keep
  estimates identical.
- All other vignettes (14/14 currently build, STATUS.md): untouched.

## Verification (AGENTS.md lifecycle, per phase — exact order)

At every phase: `devtools::load_all()` → `devtools::document()` (only if
roxygen touched; Phases 1–2 add no exports, so NAMESPACE is unchanged and
`document()` must come out idempotent) → `devtools::test()` → scoped check.
Phased gates:

1. **Phase 1 (compat module + canaries, no behavior change)**
   - Gate: **full suite 596 pass, 0 fail** (identical numbers to STATUS.md
     baseline — a no-behavior-change phase must not move them);
     `document()` idempotent; no new exports; the five grandStandardized
     internal sites behave identically (existing MG A/B test covers them).
   - Scoped check: `R CMD build --no-build-vignettes .` +
     `R CMD check --ignore-vignettes --no-manual <tarball>` (the scoping
     trick documented in PLAN 02/STATUS.md).
2. **Phase 2 (schema + renderer + cutover)**
   - Gate: A/B identity (estimates/vcov/row order, SG + MG + 2-factor +
     interaction); all `test-tspa*.R`, `test-grandStandardizedSolution.R`,
     `test-get_fs_int.R` green; `vcov_corrected()` MG run finite/symmetric.
3. **Phase 3 (hardening & docs)**
   - CI matrix live; vignette re-knit (2 + 1 files, above); STATUS.md note;
     plan archived to `archive/` on completion.
   - Full `devtools::check()` (builds all vignettes): expect the pre-existing
     **2 WARNING / 3 NOTE baseline identical** (S3 consistency ×2 WARNINGs;
     `.lintr`, unused `Matrix` import, top-level files NOTEs) with **no new
      findings**; `DESCRIPTION` unchanged per decision 1.

## Risks

- **0.7-2 APIs (`lavTech`, `lavaan.data.frame`) are 3 weeks old** and may
  themselves shift before stabilization; the user chose NOT to pin a
  `DESCRIPTION` bound (decision 1) → hedged by: the compat module as the only
  coupling point, capability probing (not version strings), the compat error
  naming the tested-up-to version, the CI dev slot, and syntax-first
  rendering (the estimation path never parses lavaan's table format).
- **`vcov_corrected()` re-fit path** (slow, per-Jacobian `eval(call)`) must
  stay bit-identical → Phase 2 A/B gates it explicitly; any `tspa()`
  argument-name or `tspa_call` change is a release-breaker.
- **Vignettes print `tspaModel`** → renderer changes alter printed text →
  controlled re-knit of only the 2 affected files (plus the interaction
  vignette's workaround removal).
- **PLAN 03 uncommitted work** in shared repo state (3 R files, 2 Rd files,
  1 new test file) → **land it first**; this plan's files barely overlap
  (`tspa.R`, `tspa_corrected_se.R`, `grandStandardizedSolution.R`, new files),
  but `tests/` and `man/` share the repo and an interleaved commit history
  would muddy the A/B baseline.
- **Math-layer coupling stays**: `get_fscore_math.R`'s
  `lavInspect(what = "data"/"meanstructure"/"implied")` sites are documented
  in the contract table but not wrapped in Phase 1 — if a future lavaan
  changes *those* views, detection falls to the existing equivalence tests
  (`test-lavPredict_equivalence.R`, `test-compute_fscore.R`) rather than the
  canary. Acceptable: those are the most stable `lavInspect` views and the
  math layer is the numerics-sensitive touchpoint AGENTS.md warns to keep
  frozen.

## Files touched

- `R/lavaan_compat.R` — **new** internal compat module (contract table,
  `tsp_layout()`, `tsp_partable_read()`, `tsp_model_matrices()`,
  `tsp_free_matrices()`, `tsp_partable_mats()`, `tsp_nobs()`,
  `tsp_ngroups()`, `tsp_norig()`).
- `R/tspa.R` — schema constructor + `tspa_render()` (internal, co-located);
  `tspa_sf`/`tspa_mf` bodies rewritten onto them; temporary
  `render = "string"` fallback.
- `R/grandStandardizedSolution.R` — five internal-touch sites migrated onto
  the compat module.
- `R/tspa_corrected_se.R` — `@Data@ngroups` slot site → `tsp_ngroups()`
  (contract-only change).
- `tests/testthat/test-lavaan_compat.R` — **new**: golden partable canary,
  wrapper A/B tests.
- `tests/testthat/test-tspa_render.R` — **new**: schema construction,
  renderer A/B, `tspa()` contract tests (attributes, re-callability),
  interaction-alias test.
- `.github/workflows/R-CMD-check.yaml` — lavaan axis (0.7-2 + dev, scoped
  non-baseline slot).
- `vignettes/R2spa.Rmd`, `vignettes/multiple-factors.Rmd`,
  `vignettes/get_fs_int-vignette.Rmd` — controlled Phase-3 re-knit /
  workaround removal (no narrative changes).
- `STATUS.md` — note added at Phase 3 (new open item for the deferred F1
  cross-reference, if not already present).
- **Unchanged by design**: `DESCRIPTION` (decision 1), `NAMESPACE` (no new
  exports), `R/tspa_mx.R`, `R/get_fscore*.R` / `R/get_fs_methods.R`
  (PLAN 03's territory), all other vignettes and their fixtures.
