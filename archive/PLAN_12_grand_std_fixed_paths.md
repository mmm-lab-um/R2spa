# PLAN 12 — Report fixed structural paths in `grand_standardized_solution()`

**Date:** 2026-08-24
**Owner roles:** `r-architect` (code), `r-tester` (tests), `r-doc` (roxygen/NEWS)
**Status:** implemented (2026-08-24) — decisions D1 (real delta SE) + D2 (slopes-only scope) delivered; R1/R2 verifications closed (§9).
**Blocked by / relates to:** none. Builds on the uncommitted grand-std row-mapping bug fix + `tsp_partable_positions()`/`tsp_check_free_positions()` compat helpers.
**Lava canary (verified 0.7-2):** `op=="~"` is slope-only (means are `op=="~1"`, excluded upstream by `i_struct`). `lavInspect(fit,"est")$beta` **does** carry row/col variable dimnames (SG flat block list; MG one list per group in group order) — the fixed-slope anchor. `lavTech(fit,"est")$beta` strips dimnames, so the dimname source must stay `lavInspect`.

---

## 1. Problem

`grand_standardized_solution()` (`R/grandStandardizedSolution.R`) refuses models in which any
structural **outcome** is fully fixed:

```
R/grandStandardizedSolution.R:126-129
  out_positions <- p$free[p$op == "~" & p$user == 1L & p$rhs %in% out_names]
  ...
  if (any(out_positions == 0))
    stop("grand_standardized_solution(): every structural parameter must be free; ",
         "fixed structural paths are rejected.")
```

`out_positions` are the partable **free positions** of the outcome's `~` rows; a user-fixed slope
(e.g. `dem60 ~ 1.5*ind60`) has `free=0`, which trips the guard and aborts the **whole** solution.

This is inconsistent with `lavaan::standardizedSolution()`, which reports both `est.std` **and**
a nonzero delta-method SE for fixed slopes.

Verified on `lavaan 0.7-2`, `PoliticalDemocracy`, CFA `ind60/dem60 + ind60 ~ dem60`:

| fit | `estimated` | `std.err` (lavaan) | `est.std` (lavaan) |
|---|---|---|---|
| free slope (`dem60 ~ ind60`) | 1.448 | 0.329 | 0.443 |
| fixed slope (`dem60 ~ 1.5*ind60`) | 1.5 (fixed) | **0.057** | **0.473** |
| (2S-PA corrected, free, from `vcov_corrected`) | 1.448 | 0.386 | 0.444 |

- The plain-SEM grand-std value is `lavaan est × sqrt(Var(ind60))/Var(dem60)`:
  `1.5 × 9.76/31.87 = 0.473` ✓, `1.448 × … = 0.443` ✓ — the ratio rescaling is correct.
- The fixed slope's `std.err 0.057` is a **delta** SE propagated through the ratio
  `c × sd_x/sd_y` with `V(x,y)` (involving `dem60~~dem60` **and** `ind60~~dem60`).

## 2. Root cause (why free-only)

The function standardizes the **entire** structural coefficient matrix `β` (free **and** fixed),
not just free positions:

- `std_beta_est()` `:217-293` — every `β` cell `(i,j)` gets a `tmp_std_beta` value using
  `sd_y = sqrt(ψ_{y,y})`, `sd_x = sqrt(ψ_{x,x})`, and (cross) `cov_xy`.
- `grand_std_beta_est()` `:262-292` — same, over the grand SD (pooled `N_g`).
- The SE Jacobians `J` (grand :239-248) and `A` (plain :235-236) are built **over all free β
  (`b_idx`) and all free ψ (`psi_idx`) for every cell**, including fixed ones.

So the plain-SEM `estimated`/`est.std` **and** the first-order delta SE for a fixed slope are
**already computed** in `tmp_std_beta` / `tmp_acov_std_beta`. Only the output **selection** is
free-only:

```
R/grandStandardizedSolution.R:144-148   est.std anchor
  for (i in seq_along(out_rows[[b]]))
    out_idx[[b]] <- which(beta_free[[b]] == out_positions[[b]][i])     # fixed → 0 → NA

R/grandStandardizedSolution.R:167-169   se anchor
  tmp_se[[b]]  <- sqrt(tmp_acov_std_beta)
  beta_se[[b]] <- tmp_se[[b]][out_idx[[b]],  ]
```

A fixed cell has no free-β position, so `which(...) == integer(0)` → `out_idx` slot is `NA` →
nothing is reported. The guard at `:126-129` exists to fail fast on this gap.

## 3. Scope & decisions (confirmed)

- **D1 — SE for fixed slopes = the real delta (matching `standardizedSolution()`'s `std.err`).**
  Report the first-order delta SE; do **not** emit a bare ratio rescale and do **not** `NA` it.
- **D2 — Scope = structural slopes (op `"~"`) only.** No partial fix of other `user=1` `free=0`
  entries (intercepts/means). A defensive guard keeps the "must be free" rejection for **non-slope**
  rows. The `lhs ∈ out_names` filter already confines output to structural rows.
  - Consequence: a **fixed intercept or mean** (e.g. `ind60 ~ 1*1`) still reaches the guard (op `~`,
    user 1, free 0) → still rejected. Correct under D2.
  - Consequence: a **fixed regression coefficient** (`dem60 ~ 1.5*ind60`) now reports
    `estimated=1.5`, `est.std=0.473`, `se=0.057` (delta, matches `standardizedSolution()`).

## 4. Approach (r-architect)

The only code paths that need to change:

1. **Select by cell identity, not free position.** Replace `out_idx` (free-position → within-matrix
   row) with a **cell `(row, col)` resolved by name** from the partable row:
   `out_lnm = out_rows[[b]]$lhs` (target row in `β`), `out_rnm = out_rows[[b]]$rhs` (input col).
   Every `β` cell is present in `tmp_std_beta`/`tmp_acov_std_beta` (free or fixed), so no new
   computation is needed — only a correct row/col into the full matrices.
   - Keep the existing free-position path as the primary anchor for free cells (unchanged); add the
     name-based cell index as a fallback that is used **iff** the free-position anchor is `NA`.
   - This preserves bit-identical output for all-free models (regression).

2. **Names of the `β`/`ψ` matrices.** `lavInspect(fit, "est")$beta[[b]]` arrays carry **no
   dimnames** (verified on plain-SEM SG + MG; consistent with the plan-04 audit of `@Model`). The
   cell must therefore be located via dimname sources that actually exist:
   - **Primary (verify in step 1):** `lavInspect(fit, "est")$beta` dimnames if the install
     populates them, else
   - **Fallback:** `fit@Model`-level parameter names, or the partable row index into
     `lavInspect(fit, "est.mat")` (the full `est` matrix is 1-indexed over all free+fixed params and
     carries `parameter` dimnames). `est.mat` is the robust route: `which(rownames(est.mat) ==
     paste(lhs, rhs, sep = op))` gives the global row; combine with the `beta_free` column layout
     to recover `(row, col)`.
   - **Decision point (verify before coding):** confirm which dimname source is stable across
     plain-SEM SG, plain-SEM MG, and 2S-PA-corrected SG/MG. If `est.mat` dimnames are stable (they
     are — they're the `parameter` strings), build a small helper
     `tsp_beta_names(fit, b)` returning `list(lhs = cnames, rhs = rnames)` for group `b`, co-located
     in `R/lavaan_compat.R` (lavaan-internals boundary) with a canary test. This keeps the
     dimname-probing out of `grandStandardizedSolution.R` and consistent with the existing
     `tsp_partable_positions`/`tsp_check_free_positions` helpers there.

3. **Rewrite the guard.** `:126-129` becomes: keep rejecting on any **non-slope** fixed output row;
   for **slope** rows, no longer require `free>0`. i.e.
   - drop the blanket `any(out_positions == 0)` stop;
   - the only hard failure becomes: a fixed row that is not an outcome slope (defensive, per D2).
   - Document in roxygen (D2): fixed structural slopes **are** reported; all other fixed outcome
     entries (intercepts/means) are not.

4. **`which_free` / `J` / `A` are unchanged.** They already cover every `β` cell. No touch to
   `std_beta_est`/`grand_std_beta_est`.

5. **2S-PA path (corrected in-place / standalone).** The same `out_idx`/`J`/`A` code serves the
   delta SE over the corrected `vcov` (passed in as `v`). A fixed slope's fixed-β cell has no free
   β, but the Jacobian already differentiates it w.r.t. free β and free ψ, so the corrected SE for a
   fixed slope is valid and is what `J %*% v %*% t(J)` produces. **No separate branch** needed; the
   name-based selection fix applies identically. (Verify against the plain `standardizedSolution()`
   std.err with a corrected `vcov` for at least one SG case; see test §5.)

## 5. Test plan (r-tester)

New `tests/testthat/test-*` (extend the existing `grandStandardizedSolution` / `lavaan_compat`
test files; do **not** create a new file unless conventions demand it):

1. **Fixed slope, plain SEM, SG:** `ind60 ~ dem60` model with `dem60 ~ 1.5*ind60`. Assert
   `grand_standardized_solution(fit)` returns a row with `estimated == 1.5` (or the user value),
   `est.std ≈ standardizedSolution(fit)` est.std (`0.473`, tolerate 1e-6), `se ≈ std.err`
   (`0.057`, tolerate ~5e-3 for first-order delta). Guard: `se > 0` and finite.
2. **Fixed + free mix, SG:** same model **without** the fixed slope (free `dem60 ~ ind60`) plus any
   other fixed structural. Assert the free row is unchanged (regression vs pre-change values) and
   the fixed row matches `standardizedSolution()` as in (1). Guards: `identical` free-row values to
   the current output (snapshot) + fixed-row values to `standardizedSolution()` reference.
3. **Fixed slope, plain SEM, MG:** group model with one fixed structural slope in one group (via
   `group.equal` or explicit per-group). Assert the fixed row's `est.std`/`se` match
   hand-calculated grand ratio on that group's grand SD, and that the non-fixed groups are
   unaffected. (Use `group.equal="regressions"` / `group.free` — **not** the colon-statement group
   prefixes, which are non-functional in this lavaan build; see cross-group work.)
4. **Fixed slope, 2S-PA corrected SG:** run `tspa(..., corrected_se=TRUE)`, assert the fixed slope's
   `se` equals a hand-built `sqrt(J v J')` at that cell (reuse `vcov_corrected` internals via a
   small test harness), and `est.std` equals `1.5 × sd_ind/sd_dem` on the corrected variances. This
   is the delta-vs-delta consistency check (matches plain `standardizedSolution()` std.err on the
   corrected fit).
5. **Regression (all-free unchanged):** the free-slope case from §1 without the fix produces
   byte-identical `est.std`/`se` to the values captured before the change (guards against the
   name-based anchor perturbing the common path).
6. **`tsp_beta_names` canary:** asserts the dimname source used by the new helper is stable (row/col
   order matches the partable/label order) on the SG free-slope and SG fixed-slope fits.
7. **Defensive guard (D2):** a model with a **fixed intercept** on an outcome (`dem60 ~ 1*1` fixed,
   e.g. `0*1` or a fixed mean) still errors, with the message updated to name the offending
   non-slope row.

Tolerance guidance: `est.std` vs `standardizedSolution()` `tolerance=1e-6` (same ratio, exact);
`se` vs `std.err` `tolerance ~5e-3` (first-order delta, documented as an approximation).

## 6. Docs (r-doc)

- **Roxygen `R/grandStandardizedSolution.R`** `@details` (currently `:16-17`): replace
  "The implementation currently assumes every structural parameter is free; lavaan models that
  fix or constrain one or more structural parameters are rejected."
  with: "Structural **slopes** may be freely estimated or fixed by the user. A fixed slope is
  reported alongside the free ones: its `est.std` is the user value rescaled by the (grand) SD
  ratio, and its `se` is the first-order delta approximation (matching
  `lavaan::standardizedSolution()`, which also reports a delta SE for fixed slopes). Non-slope
  outcome entries (intercepts/means) must remain free; a fixed one is rejected."
  The `@details` already states the ratio formula; keep it.
- **`man/grand_standardized_solution.Rd`:** regenerate via `devtools::document()` (never hand-edit).
- **`NEWS.md`:** add a bullet under the current dev heading:
  "`grand_standardized_solution()` now reports user-fixed structural slopes (est.std + delta se,
  matching `lavan::standardizedSolution()`), instead of rejecting the whole solution when any
  structural parameter is fixed."

## 7. Acceptance

- `devtools::load_all()` → `devtools::document()` → `devtools::test()` → `devtools::check()` all
  green (0 errors / 0 warnings; expected NOTEs unchanged).
- `sem("ind60 ~ a1*x1 + ...; dem60 ~ b1*x5; dem60 ~ 1.5*ind60", PoliticalDemocracy)` now yields, from
  `grand_standardized_solution()`, an `est.std ≈ 0.473`, `se ≈ 0.057` for the fixed path — no error.
- All-free models reproduce previous `est.std`/`se` exactly.

## 8. Out of scope

- Partial fix of **intercepts/means** on outcomes (D2).
- `lavaan::sam()` integration (separate SAM line; method choice is `sam_method`, not `M.method`).
- Changing `which_free`/`J`/`A` (already correct; no touch).
- Multi-group free-position mapping (already handled by the uncommitted `tsp_partable_positions`
  fix; not part of this plan).

## 9. Verification log (closed 2026-08-24)

**R1 — β dimname source: RESOLVED.** `lavInspect(fit, what="est")$beta` is a `lavaan.matrix`
with row/col **variable dimnames**: single group → flat block list (`est$beta`), multiple groups
→ one list per group in group order (`est[[g]]$beta`). `lavTech(fit, what="est")$beta` has the
same dims but **no** dimnames. So the dimname anchor is `lavInspect(est)$beta`, exposed by
`tsp_beta_names()` (handles the SG-flat / MG-nested shape via `tsp_ngroups`).

**R2 — `op=="~"` slope-only: RESOLVED.** Structural means/intercepts carry `op=="~1"`, **not**
`"~"`, so `i_struct <- which(partable$op == "~")` excludes them upstream (verified on free-LV-mean
and free-observed-mean models). Every `~` row is a beta block (observed-on-observed included —
beta rows/columns expand to include observed endogenous/exogenous variables).

**Column-major correction.** The col-major position of cell `(lhs = row r, rhs = col c)` in an
`nrow × ncol` matrix is `(c - 1) * nrow + r` (stride = **nrow**, multiplier = `c-1`). Verified
equal to the free-position anchor (`which(beta_free == pos)`) on a free slope (both = 2 for
`dem60 ~ ind60`).

**D2 guard scope.** Since every `~` row resolves to a beta cell, the `out_idx == 0` guard is a
drift/NA safeguard, not reachable from a real `~` row — tested by mocking `tsp_beta_names`
(`testthat::local_mocked_bindings`, testthat 3.3.2) to drop the predictor.

**Live acceptance (lavaan 0.7-2, R 4.6.1).** SG fixed `dem60 ~ 1.5*ind60`: `est.std 0.473`,
`se 0.0569` = `standardizedSolution()` exactly. Mixed free+fixed (3-factor): all rows match. MG
fixed `0.5*speed`: identical grand-std `est.std`/`se` in both groups. 2S-PA corrected SG fixed:
`se 0.0423` = `standardizedSolution(corrected)`. All-free models bit-identical (regression).
Targeted + full testthat suite: 0 failures (only the pre-existing `tspa_mx` negative-lv-variance
warning).
