# PLAN 11 — mirt per-obs multi-factor `tspa()` (close PLAN 09 §8 gap)

Closes the deferred item in `archive/PLAN_09_per_unit_pooled_fs.md` §8:
"*mirt per-obs `tspa()` — the per-obs `fsL`/`fsT` list is not pooled
(`resolve_fs_per_row` supports it, but `is_per_unit_fs` deliberately does not
trigger on it); remains unsupported.*" Scope = **both** single-group and
multi-group **multi-factor** mirt. Baseline: `refactor/core` @ `481716e`
(branch `mirt-mf-tspa`). `load_all()`-only (no shared-library install).

## 0. Task framing (from the ask)
`get_fs()` already supports mirt (`.SingleGroupClass`, 1-D/`q`-factor;
`.MultipleGroupClass`, per-obs trailing `group` column, PLAN 08/10). The
`get_fs()` output is per-observation: `fsL`/`fsT`/`fsb` are **flat lists of n
matrices/vectors** (`R/get_fs_methods.R:1263-1269`) + `mirt_per_obs = TRUE`.
`fs_indiv()` already handles that shape. But `tspa()` does **not** pool it and
mis-reads the n-row list as n "groups". Make multi-factor mirt poolable so
`stage-1 mirt fit → get_fs() → tspa()` works, with and without completely-
missing rows.

### 0a. Conflict check (2026-08-23)
Baseline drifted from `2177783` to `481716e` ("multigroup corrected_se"). That
commit edits `R/tspa.R` in two regions only — the `corrected_se` roxygen
(top) and the corrected-se block (removed the `tsp_ngroups > 1` guard) —
neither overlaps this plan's edit regions (`is_per_unit_fs`, the call site,
`pool_per_unit`). No other plan file changed. **No conflict; re-anchored to
content.**

## 1. Current state / feasibility (analysis)
**Two precise gaps, both in `R/tspa.R` (the resolver is already done):**
1. **Trigger** — `is_per_unit_fs(fsT, fsL)` is `FALSE` for the mirt flat
   per-row list ⇒ `pool_per_unit()` is never invoked.
2. **Group** — `resolve_per_obs()` (`R/fs_indiv.R:306`) hard-codes
   `group_col = NULL, group_vals = NULL`, discarding the `group` column that
   `get_fs.MultipleGroupClass()` emits (`R/get_fs_methods.R:1256-1258`), so
   MG mirt would pool all rows into **one** group.

**Reusable (no new math):** `resolve_fs_per_row()` → `resolve_per_obs()`
already mints one block per row; `fs_row_cols()` yields per-row se/ld/ev
(mirt-compatible); `pool_per_unit()`'s within-group reduction loop is
shape-agnostic and just needs per-row `group_vals`.

**Already supported (unchanged, no work):** single-factor mirt (SG & MG) via
the `se_fs` path — `pool_se_fs()` reduces the per-row `fs_<v>_se` columns by
the `group` column with `na.rm`. `corrected_se`/`vfsLT` is N/A for mirt (no
`vfsLT` source; `corrected_se = TRUE` errors as designed).

## 2. Design — marker-gated trigger + pool-side group synthesis
Do **NOT** touch `resolve_per_obs()`. It is shared with `fs_indiv()`, whose
MG-mirt behavior is pinned to *drop* the group column
(`tests/testthat/test-get_fs_mirt_multigroup.R:94`
`expect_identical(setdiff(names(fs1), "group"), names(fs_indiv(fs1)))`).
Carrying `group` in the resolver would flip that invariant. Keep all changes
inside `R/tspa.R`.

**Trigger: marker-gated, not shape-gated.** A bare flat list of matrices is
ambiguous; the `mirt_per_obs` attribute is authoritative and only set by the
two mirt `get_fs()` methods. Single-factor mirt has `fsT = NULL`, so
`is.list(fsT) && length(fsT) > 1` is `FALSE` ⇒ it stays on the `se_fs` path.

## 3. `R/tspa.R` changes (owned by **r-architect**)
- `is_per_unit_fs(fsT, fsL)` → `is_per_unit_fs(fsT, fsL, mirt_per_obs =
  FALSE)`; add `|| (mirt_per_obs && is.list(fsT) && length(fsT) > 1L)`.
- Call site: `is_per_unit_fs(fsT, fsL, mirt_per_obs = isTRUE(attr(data,
  "mirt_per_obs")))`.
- `pool_per_unit()`: after `resolved <- resolve_fs_per_row(fs)`, synthesize
  effective `g_col`/`g_vals` (fall back to the mirt `group` column when
  `mirt_per_obs` set + `"group" %in% names(fs)`), then use `g_col`/`g_vals` in
  the multigroup branch; `glabs` from the factor **levels** when `mirt_mg`
  (guarantees a K-entry list == lavaan `group=` levels). SG mirt (no `group`)
  falls through to the single-group reduction unchanged.

## 4. Documentation (owned by **r-doc** — roxygen only, `R/tspa.R`)
- `@param fsT` / `@param fsL` / `@param reduce`: extend the "per-unit
  (per-pattern FIML / per-cluster merMod)" wording to include **per-observation
  mirt** (SG & MG multi-factor), reduced per group by `reduce`.
- Internal comments on `is_per_unit_fs` (marker-gated path) / `pool_per_unit`
  (group synthesis) and the call-site comment.
- `NEWS.md`: one bullet under the current top version.
- No new `man/`; `man/tspa.Rd` regenerates. `NAMESPACE` unchanged. No
  `DESCRIPTION`/new dep (`mirt` already `Suggests`).

## 5. Tests (owned by **r-tester** — `tests/testthat/test-tspa_pooled.R` only)
Gate the new block behind `skip_if_not_installed("mirt")`; build a 2-factor
SG `mirt::mirt` fit and the existing 2-factor `mirt::multipleGroup` fixture (2
groups, one completely-missing row).
1. **SG multi-factor:** `tspa("F2 ~ F1", data = fs, fsT = attr(fs,"fsT"),
   fsL = attr(fs,"fsL"))` fits; pooled `fsT`/`fsL` == scorable per-row matrix
   mean (hand A/B); model string == a manual `tspa_mf` fed the hand-pooled
   matrices.
2. **MG multi-factor:** `tspa("F2 ~ F1", data = fs, group = "group", fsT =
   ..., fsL = ...)`; one plain `q×q` per group; per-group per-row means match
   hand computation; **group order == `levels(fs$group)`**; `attr(fit,"fsT")`
   is a K-named list.
3. **NA exclusion:** the completely-missing row (`group = NA`) is absent from
   every group's reduction (hand A/B over scorable rows only).
4. **Update existing guard `test-tspa_pooled.R:437-444`:** keep
   `expect_false(is_per_unit_fs(fsT_flat, fsL_flat))` (default marker ⇒ shape
   alone insufficient); add
   `expect_true(is_per_unit_fs(fsT_flat, fsL_flat, mirt_per_obs = TRUE))`.
5. **`reduce = "median"`** exercised once on mirt (PSD guard reachable).

## 6. Verification lifecycle (AGENTS.md order)
1. `devtools::load_all()`. 2. `devtools::document()` (roxygen only; expect
`man/tspa.Rd` regen, `NAMESPACE` byte-identical). 3. `devtools::test()`
(new mirt block + all existing). 4. `devtools::check()` — expect **0/0/0**
(mirt tests conditional).

## 7. Subagent role assignments
- **r-architect** — §3; dependency/contract audit vs `481716e`.
- **r-doc** — §4. **r-tester** — §5.

## 8. Out of scope
- Single-factor mirt (already works via `se_fs`; optional assertion only).
- `corrected_se`/`vfsLT` for mirt (no `vfsLT` source).
- Cluster-size / per-observation weighting (unweighted `reduce = mean`
  default, same policy as merMod, PLAN 09 §8).
- Free-structural-per-pattern mirt sub-stage-2 (not applicable).

## 9. Risks / regression safety
- **Complete-data lavaan/merMod:** `g_vals`/`g_col` are `resolved$*` on every
  non-mirt path (identity); `test-tspa.R` / `test-tspa_render.R` pass
  untouched.
- **`fs_indiv()` MG-mirt:** untouched (`resolve_per_obs` unchanged) ⇒
  `test-get_fs_mirt_multigroup.R:94` stays green.
- **Single-factor mirt:** `fsT = NULL` ⇒ `is_per_unit_fs` `FALSE` ⇒ `se_fs`
  path unchanged.
- **K-group mismatch:** `glabs` from mirt factor **levels** ⇒ K-entry pooled
  list == lavaan `group=` levels (avoids the `length(fsT) != ngroup` schema
  error a naive `unique()` would risk).
- **PSD:** `mean` of PSD per-row `fsT` is PSD; `median` guarded.
- **roxygen2 8.1.0** pinned; never hand-edit `NAMESPACE`/`man`.
