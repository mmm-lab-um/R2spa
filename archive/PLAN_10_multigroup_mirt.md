# PLAN 10 — Multi-group mirt `get_fs()` support (`MultipleGroupClass`)

**Status: IMPLEMENTED** (2026-08-22, `31b3809` — all phases P0–P4 complete).
See `STATUS.md` Closed **#18**. Follow-up **F7** (resolved). Built on PLAN 08
(single-group mirt `get_fs()`) and the ψ fix (`fc16f81`). Final suite **3358
pass** / 0 fail; `R CMD check` 0/0/0. `mirt` stays a `Suggests`-only dependency
(namespaced `mirt::` calls, `require_mirt()` guard).

> **Implementation outcome / deviations discovered during P0.** (1) `mirt::multipleGroup()`
> requires `model` as a **number** of factors (not the string `"1"`) and `group` as a
> factor/character **vector** of length `nrow` — the string form leaves `model$x` a 0-row
> template → `subscript out of bounds`. (2) The standard identifiable form is
> `invariance = "slopes"` (metric): it fixes **every** group's factor to `(0, I)`, so the
> per-group EAP is on a common standard-normal scale and the group mean is carried by the
> item intercepts — hence `alpha = 0` by default (identical to the single-group path), and a
> per-group `psi` is the estimated factor covariance (coincides with `diag(q)` under metric
> invariance; differs only under `free_mean(s)`/`free_var`). (3) mirt **drops completely-
> missing rows** from every extraction (`N`/`group`/`rowID` are scorable-length); the full-row
> `group` column is rebuilt as `groupNames[group][keep]`, `NA` for the missing rows. The
> `format = "list"` (per-group data frames) is not implemented — `mirt` output is always a
> single per-observation data frame + `group` column, like the single-group path.

## 1. Scope
Implement `get_fs.MultipleGroupClass()` — today a `stop()` stub at
`R/get_fs_methods.R:1078` ("Multi-group mirt models are not supported … fit or
extract a single group first") — so a fitted multi-group IRT model
(`mirt::multipleGroup` → S4 `MultipleGroupClass`) yields per-observation
`fs_<factor>` scores with **per-row** loadings / SEs / error covariances, plus the
group structure (a `group` column in `"unified"`; per-group data frames in
`"list"`) and **per-group** `psi`/`alpha`. It reuses the SG `mirt_per_obs`
per-row engine and `fs_indiv()`'s per-row branch, so `fs_indiv(get_fs(<mg mirt
fit>))` works with truly row-specific, group-correct values.

## 2. What already exists (verified)
- `get_fs.SingleGroupClass()` (`fc16f81`): per-row `fsL = I − Vpost_i·Ψ⁻¹`,
  `fsT = fsL·Vpost_i`, `fsb`; **flat** per-row list attrs (`fsL`/`fsT`/`fsb`,
  length `n`) + the `mirt_per_obs` marker; `psi = mirt_full_cov(fit)` (the full
  estimated factor covariance); no `group` column. Single-group only.
- `get_fs.MultipleGroupClass()`: `stop()` stub. S3 dispatch on the S4 class
  already works for `SingleGroupClass` (exported + `S3method`), so the
  `MultipleGroupClass` method dispatches the same way.
- `mirt_full_cov(fit)` (`R/get_fs_methods.R:896`): reads a **SingleGroupClass**'s
  `coef(fit)$GroupPars` `COV_ij` → full `q×q`. Reused per group below.
- `fs_indiv()` per-row branch `resolve_per_obs()` (`R/fs_indiv.R:249`): reads
  **only** the per-row `fsL`/`fsT`/`fsb` list attrs (length `n`); does **not**
  read `psi`/`alpha`; mints one block per row and calls
  `make_resolved(…, group_col = NULL, group_vals = NULL)`. So the mirt per-row
  contract is independent of the latent-moment attributes and currently ignores
  any group column (see §9 OQ).

## 3. mirt multi-group data plumbing (VERIFIED empirically, mirt 1.46.1)
- **Fit:** `mirt::multipleGroup(items, model, group = <group-vector>, invariance =
  ..., method = "EM")`. `model` = `1` (1-factor), a character, or a
  `mirt.model()` object — **a raw matrix is NOT accepted by `multipleGroup`**
  (unlike `mirt::mirt`), so ≥2 factors need `mirt.model()`. `invariance` =
  `""` (configural / separate per-group), `c("slopes")` (metric),
  `c("slopes","intercepts","free_var","free_means")` (scalar), etc. The first
  factor level is the reference group.
- **Scores/SEs:** `fscores(object, full.scores = TRUE, full.scores.SE = TRUE)` →
  an `n×(2q)` **matrix** (cols `<factor>`, `<factor>.SE`), **no group column**,
  rows in the original (global) data order.
- **Per-obs posterior cov:** `fscores(object, full.scores = TRUE, return.acov =
  TRUE)` → list, one `q×q` matrix per observation, global order. **Works for
  `MultipleGroupClass`** (verified: length `n`, correct dims).
- **Group membership:** `mirt::extract.mirt(object, "group")` → the grouping
  factor in global order (length `n`). This is the only group source (the
  `fscores` output carries none).
- **Per-group `ψ`:** `mirt_full_cov(mirt::extract.group(object, g))` —
  `extract.group` returns a real `SingleGroupClass`, so the existing helper runs
  unchanged per group (verified: per-group `MEAN`/`COV` present on the
  extracted object).
- **NOT extractable directly on the multi-group object:**
  `extract.mirt(object, "theta")` and `("lambda")` error (`Could not extract
  element`) — per-group (Λ, θ) live in the per-group sub-models, not one table.
  `extract.group()` is the clean per-group route.
- **`coef(object)` (multi-group):** a list whose data-frame form is **wide** —
  one column per `(<g>.<item>.<par>)` plus `(<g>.GroupPars.MEAN_i)` and
  `(<g>.GroupPars.COV_ij)`. (Direct per-group `ψ` source; see D2.)
- `extract.mirt(object, c("nfact","factorNames","itemnames"))` behave as in SG.

## 4. Key modeling decisions
- **D1 — score from the whole fit, not per extracted group.** `fscores(mfit)`
  scores every observation with its group's parameters and honors cross-group
  invariance constraints (shared Λ/δ); re-scoring `extract.group(g)` separately
  matches only for `invariance=""` (configural). → Use `fscores(mfit)` for
  scores/SEs/`acov`; use `extract.group(g)` **only** for per-group `ψ`.
- **D2 — per-group `ψ`.** For an observation in group `g`, use
  `Ψ_g = mirt_full_cov(extract.group(mfit, g))`. The SG scalar
  `mirt_full_cov(fit)` becomes a per-group map built once (`n_groups`
  `extract.group` calls), then indexed by each row's group. Perf fallback: parse
  the wide `coef()` `COV_ij` columns once (skip per-group `extract.group`).
- **D3 — latent mean / `alpha`.** Mirror SG: `alpha = 0` (default) or a user
  `prior_mean` (length `q`, **shared across groups** — the `prior_*` "shared
  across groups" convention). Do **not** default `alpha` to the per-group mirt
  means (`<g>.GroupPars.MEAN_i`); that is a distinct modeling choice (§8).
- **D4 — reuse the SG per-row engine verbatim.** The per-row
  `compute_lav_fs_matrices(Vpost_i, Ψ_{g_i}, alpha, "regression")` loop is the
  SG loop with `Ψ` switched per row by group; column naming, `fs_row_cols`, `se`
  / loadings / ev layout, and the all-NA missing-row convention are all unchanged.
- **D5 — output grouping.** `"unified"` (default): one data frame, one row per
  observation in global order, a trailing `group` column; **flat** per-row
  `fsL`/`fsT`/`fsb` list attrs (length `n`, as SG); `fs_pattern = list(label =
  seq_len(n), pat = NULL)`; `mirt_per_obs = TRUE`; `psi` = a **named list by
  group** (each `q×q`, mirroring the lavaan-MG `psi` shape) and `alpha` = the
  shared `q`-vector (OQ3). `"list"`: the unified df split by `group` into
  per-group dfs (group column dropped per df), each carrying its per-row attrs +
  that group's `psi`/`alpha`, plus list-valued outer attrs (mirror
  `assemble_fs_blocks()`'s list branch / `fs_to_group_list()`).
- **D6 — `scoring_matrix`.** Not attached (matches SG mirt), §8.

## 5. `get_fs.MultipleGroupClass()` contract
`get_fs.MultipleGroupClass(object, prior_mean = NULL,
format = c("unified","list"), ...)` (parallel to SG; `prior_cov` /
`corrected_fsT` / `vfsLT` / `reliability` / `method` unsupported → rejected like
the SG path).
1. `require_mirt()`; assert `inherits(object, "MultipleGroupClass")`.
2. `q`/`fn` from `extract.mirt` (same guards as SG); `fs_names = paste0("fs_", fn)`.
3. `grp = extract.mirt(object, "group")` (length `n`); `group_labels =
   levels(as.factor(grp))`.
4. `alpha` = `0` or `validate_fs_priors(prior_mean, NULL, fn)$mean` (D3).
5. `full = fscores(object, full.scores=TRUE, full.scores.SE=TRUE)`;
   `acov = fscores(object, full.scores=TRUE, return.acov=TRUE)`.
6. Reconcile completely-missing rows **globally** exactly as SG does (OQ4):
   scorable rows = complement of `extract.mirt(object, "completely_missing")`
   in `1:n`; `acov[[k]]` = `k`-th scorable row.
7. `psi_g <- setNames(lapply(group_labels, function(g) mirt_full_cov(extract.group(object, g))), group_labels)` (D2).
8. Per-row loop over scorable obs `i` (group `g_i`):
   `m_i = compute_lav_fs_matrices(acov[[k]], psi_g[[g_i]], alpha, "regression")`
   → per-row `fsL/fsT/fsb` (SG body, `Ψ` per group); all-NA row for the rest.
9. Assemble the df (fs/se/ld/ev via `fs_row_cols` + SG column-naming), add the
   `group` column, attach flat per-row `fsT/fsL/fsb` + `fs_pattern` +
   `mirt_per_obs=TRUE` + per-group `psi`/`alpha` (D5). `"list"`: split +
   per-group attrs (§5 / D5).

## 6. Phases
- **P0 — Spike hardening** (partly done in-planning). ✅ class / `fscores` /
  `acov` / `extract.group` / wide `coef` / `group`-vector verified (1-D). TODO:
  (a) build a working **identified 2-factor** multi-group fixture (the 2-D
  `multipleGroup` construction hit a mirt-internal shape error in the spike —
  right `mirt.model()` / data-gen / anchoring; §9 OQ1); (b) confirm the **global**
  completely-missing → scorable-row mapping for `acov` (OQ4); (c) A/B that
  `fscores(mfit)` honors invariance vs per-group SG (D1); (d) confirm `grp`
  aligns row-for-row with `fscores`/`acov` in all cases.
- **P1 — per-group ψ helper.** `mirt_group_psi(mfit)` (named list via
  `extract.group`+`mirt_full_cov`, internal) or a group-arg on `mirt_full_cov`;
  unit-test against the wide `coef()` `COV_ij` columns.
- **P2 — implement `get_fs.MultipleGroupClass`** per §5 ("unified" + "list"),
  replacing the stub.
- **P3 — consumer + docs.** Confirm `fs_indiv()` round-trips the multi-group
  mirt output (per-row branch unchanged; optional: carry the `group` column
  through `resolve_per_obs`, §9 OQ). Update `get_fs()` roxygen (drop the
  "Multi-group mirt models are not supported" line; methods list; add a
  multi-group mirt `@examples`), `NEWS`, and a vignette section (extend
  `tspa-vignette-mx.Rmd`'s mirt part).
- **P4 — tests + lifecycle** (§7): `document()` → `test()` → `check()`.

## 7. Verification / test plan (mirror `test-get_fs_mirt.R`, PLAN 08 rigor)
New `tests/testthat/test-get_fs_mirt_mg.R` (or extend the mirt file):
- **Per-group ψ:** `attr(fs, "psi")[[g]] == mirt_full_cov(extract.group(mfit, g))`
  ∀ g; names == groups; `alpha` shared.
- **Score/SE identity:** `fs_<f>`, `_<f>_se` == `fscores(mfit, full.scores=TRUE)`
  cols `[<f>]`, `["SE_<f>"]`.
- **Per-row algebra, per group:** unidim `<f>_by_fs_<f> == 1 − SE²`,
  `ev == (1−SE²)·SE²`; 2-D off-diagonal `fsL`/`fsT` vs the hand regression form
  **with that group's `Ψ`**.
- **Cross-validation vs SG (the key net):** for each group g,
  `fs_mg[group == g, ]` structure columns **==** `get_fs(extract.group(mfit, g))`
  structure columns (score/se/ld/ev + per-row matrices) — the ψ-fix validation
  pattern.
- **Grouping/structure:** `"unified"` trailing `group` col (length `n`); per-row
  list attrs length `n`; `mirt_per_obs`; `fs_pattern`; `"list"` → per-group dfs
  with matching rows/attrs + per-group `psi`/`alpha`; `fs_to_group_list`
  round-trip.
- **Invariance:** a metric fit (shared slopes) scores per group, with per-group
  `Ψ` from each group's estimated covariance (a `free_var` case) — distinct from
  a configural fit.
- **`fs_indiv(fit)`:** 1:1 per-row values on the multi-group output; NA-row
  convention; the `group` column behavior per D5/OQ.
- **Guards:** S3 dispatch of the S4 `MultipleGroupClass`; `require_mirt()`
  guard; unsupported args rejected.
- **Independence:** the tree without the method still passes 3298 (bisect), as
  in PLAN 08.
- `devtools::document()` (roxygen2 8.1.0 pinned). No new `NAMESPACE` export
  needed if `mirt_group_psi` stays internal; the S3 method registers via
  `@export` + the already-present `S3method(get_fs, MultipleGroupClass)`.
  `check()` 0 errors / 0 warnings.

## 8. Out of scope / follow-ups
- `scoring_matrix` for mirt (SG omits it) — D6.
- `alpha` = per-group mirt means (letting `fsb` capture between-group location).
- `prior_cov` for mirt (SG supports `prior_mean` only).
- `corrected_fsT` / `vfsLT` / `reliability` for mirt (unsupported in SG too).
- 2S-PA consumption: `tspa()` on multi-group mirt scores (a downstream F1/pooled
  item; `tspa()` currently rejects `mirt_per_obs`).
- mirt ≥10-factor `COV_ij` naming (mirt_full_cov's `nchar != 2` guard) — inherit.

## 9. Risks / open questions
- **OQ1 (fitting, blocks 2-D tests only):** a clean identified **2-factor**
  multi-group test model. The spike's 2-D `multipleGroup` failed with a
  mirt-internal shape error; resolution is a fixture detail (proper
  `mirt.model()` / data-gen / anchoring), not an architecture change.
- **OQ2 (group column name):** mirt takes the group as an input vector, so there
  is no user-named column to echo. Use the fixed name `"group"` (uniform); check
  it doesn't collide with an item also named `group`.
- **OQ3 (`alpha` shape):** attach `alpha` as a per-group list (entries equal the
  shared `q`-vector) to mirror `psi`'s lavaan-MG shape, or a plain scalar —
  decide during P3 (keep it simple; entries identical).
- **OQ4 (missing-data alignment across groups):** SG reconciles `acov`
  (scorable) vs `full` (all) via `completely_missing`. For multi-group the
  mapping is global; verify mirt's multi-group `completely_missing` + `acov`
  lengths align in global row order. If mirt's multi-group missing handling
  differs, fall back to per-group scoring for the missing cases.
- **Perf:** `n_groups × extract.group` is a one-time O(n_groups) cost (small);
  the per-row loop and a single `fscores(mfit)` call scale as SG.
- **Consumer drops the group column:** `resolve_per_obs` sets
  `group_col = NULL`, so `fs_indiv()` currently would drop a `group` column.
  Non-blocking (the per-row values are group-agnostic once `Ψ` is right), but a
  P3 nicety: pass it through for parity with lavaan-MG `fs_indiv`.
