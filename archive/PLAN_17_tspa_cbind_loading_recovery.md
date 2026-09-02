# PLAN 17 — `tspa()` recovers the factor-score loading from a cbind'd `get_fs()` result

**Date:** 2026-08-31
**Owner roles:** `r-architect` (code), `r-tester` (tests), `r-doc` (roxygen/NEWS/vignette pointer)
**Status:** **implemented** (2026-09-01) — branch `plan17-cbind-loading`, worktree
`/home/marklai/R2spa-wt-plan17` (forked from `origin/devel` @ `9b67694`). Full suite
4346 pass / 0 fail; `R CMD check` (`--as-cran --no-manual`) 0/0/0.
**Scope decision (2026-09-01, user-confirmed):** **narrow / derive-only** — the loading is
recovered only when `se_fs` is *auto-derived* (the exact reported footgun), NOT for an
explicit `se_fs`. This preserves the documented "explicit `se_fs` ⇒ unit loading" contract
(the OpenMx `tspa_mx_model()` path documents the same) and leaves the OpenMx path and all
explicit-`se_fs` goldens untouched. It refines §4 **D1** ("when the single-factor path runs")
to "when the single-factor `se_fs` is auto-derived." Known edge (documented in `@details`):
`do.call(tspa, attr(fit, "tspa_args"))` re-passes the derived `se_fs` as explicit, so a
replayed derived single-factor fit uses the unit loading (the original fit is unaffected;
`vcov_corrected()`, multi-factor only, is unaffected).
**Blocked by / relates to:** PLAN 13 (the auto-derivation; it already *documents* that
`cbind()` drops attributes and the single-factor path then relies on the `fs_<v>`/`fs_<v>_se`
columns). This plan closes the statistical gap that leaves: for a **shrinkage score** the
single-factor fallback silently uses a unit loading, which biases the structural coefficient.

---

## 1. Problem

`get_fs()` attaches the measurement inputs as **attributes** (`fsL`/`fsT`/`fsb`) *and* as
plain columns (`fs_<v>`, `fs_<v>_se`, `<v>_by_fs_<v>`, `ev_fs_<v>`). `tspa()` auto-derives the
inputs from a `get_fs()` result (PLAN 13): the **multi-factor** path reads the `fsL`/`fsT`
attributes (loading = `fsL`), the **single-factor** path derives `se_fs` from the
`fs_<v>`/`fs_<v>_se` columns and hard-codes a **unit loading** (1).

`cbind()` on a data frame **drops attributes** but keeps the columns. So the canonical flow
"add a distal, then fit a structural model"

```r
fs <- get_fs(lmer_fit)               # carries fsL (= w), fsT attributes
fit <- tspa("d ~ u0", data = cbind(fs, d = d))   # cbind drops fsL/fsT
```

loses the `fsL`/`fsT` attributes; the multi-factor derivation is skipped and the single-factor
path runs with **loading 1**. For a **shrinkage score** (lme4 `method="EB"`, brmsfit, mirt —
`fsL = w < 1`, the posterior-variance-derived loading `1 − Vpost/ψ`) that unit loading is wrong
and **biases the structural coefficient upward** (and the latent variance to `w²·ψ`).

## 2. Verified facts (lavaan 0.7-2, R 4.6.1, devel tree)

- `cbind(fs, d)` keeps `fs_u0`, `fs_u0_se`, `u0_by_fs_u0`, `ev_fs_u0`; `fsL`/`fsT` become
  `NULL` (confirmed).
- The implied-loading column equals the attribute: `u0_by_fs_u0 == attr(fs,"fsL")[,]` exactly
  (sleepstudy: both `0.859310`). Column names from `create_fsL_names()`
  (R/get_fscore_math.R:70, `<lv>_by_<fs>`) and `ev_`/`ecov_` from `create_fsT_names()` (:62-68).
- Single random-intercept DGP (`y~1+(1|sj)`, `d = u + v`, true β=1), `method="EB"`:
  | route | `d ~ u0` | note |
  |---|---|---|
  | single-stage SEM (gold) | 1.1047 | recovers the true model |
  | `tspa("u0 ~~ u0", data = fs)` (pure result) | ψ = 0.7694 = lme4 ψ | multi-factor path, `fsL = w`, **correct** |
  | `tspa("d ~ u0", data = cbind(fs, d))` | **1.4175** | unit-loading fallback, **+28% biased** |
  | `tspa("d ~ u0", …, fsL, fsT)` (explicit) | 1.1253 | `fsL = w`, **correct** |
- Over 150 reps the same pattern holds: unit-loading cbind ≈ +21–28% vs the SEM; `fsL = w`
  ≈ +0.3%. The bias is first-order (∝ `1/w`), not a small-sample artifact.

## 3. Root cause

R/tspa.R:500-536. The multi-factor derivation (D4) fires only when `attr(data,"fsT")` and
`attr(data,"fsL")` are non-`NULL` **and** resolve as a `get_fs()` result (provenance gate,
:506-527). When `cbind()` has dropped them, control falls to the single-factor derivation
(:528-536, `derive_sf_se_fs()`, R/tspa.R:1532), which builds `se_fs` from the `fs_<v>_se`
columns and emits a **fixed unit loading** in the schema: `tspa_schema_sf()`
(R/tspa.R:1604) writes `tspa_row(var[k], "=~", fs[k], 1, g, "struct", lab)` (:1613). The
`<v>_by_fs_<v>` column that carries the correct loading is present in `data` but never read.

## 4. Scope & decisions (confirmed with the user 2026-08-31)

- **D1 — Auto-recover the loading from the columns.** When the single-factor path runs and the
  data carries the `<v>_by_fs_<v>` implied-loading column(s), read the per-group loading from
  them (pooled exactly like the SE, via the same `reduce`) and use it in the schema instead of
  `1`. The result is then **identical** to passing `fsL`/`fsT` explicitly — no user action needed
  to cbind a distal. (Chosen over "warn-only" and over "warn now / recover later".)
- **D2 — Single-factor path only.** The multi-factor path already uses `fsL` from the
  attributes; nothing to change there. The fix lives entirely in the `se_fs`/single-factor
  branch.
- **D3 — Backward compatible by construction.** No `<v>_by_fs_<v>` column (a hand-rolled
  `fs_<v>`/`fs_<v>_se` frame) → loading stays `1` (today's behavior, unchanged). Column present
  → use it. When the implied loading is `1` (e.g. an unbiased CFA regression score) this is a
  no-op; when it is `w < 1` (shrinkage scores) it is a **correction**. No currently-*correct*
  call changes.
- **D4 — The error-variance row is unchanged.** The `ev_fs_<v>` column equals `fs_<v>_se²`
  (both are `fsT`), so `tspa_row(fs[k], "~~", fs[k], se[g,k]^2, …)` stays as-is; only the
  loading row (`:1613`) changes from `1` to the recovered value.
- **D5 — Branch.** Implement on **`devel`** (core `R/tspa.R`; the fixture is the lme4
  `method="EB"` score, no `brms` needed). Cherry-pick to `get_fs-brms`, where the same fix
  covers the brmsfit scores.
- **D6 — Out of scope (v1).** Recovering the *off-diagonal* error covariance (`ecov_*`) and the
  score intercepts (`int_fs_*`) for a **multi-factor** cbind'd result — that path needs the full
  `fsL`/`fsT` matrix, not the single-factor `se_fs`; it is already served by passing `fsL`/`fsT`
  explicitly. v1 is the single-factor (one latent + one distal) case, which is the reported
  footgun.

## 5. Approach (`r-architect`)

All in `R/tspa.R`. No new exports, no `DESCRIPTION` change, no `NAMESPACE` change. Helpers stay
co-located.

1. **`derive_sf_se_fs()` (R/tspa.R:1532) → also recover the loading.** For each latent `v` it
   already finds, detect the `<v>_by_fs_<v>` column (via `create_fsL_names()` naming, so it stays
   in lockstep with `get_fs()`). Return `list(se = <pooled se data frame (as today)>, ld =
   <pooled loading data frame>)`, where `ld` is pooled per group by the **same** `reduce` used
   for the SE (constant in the balanced case; mean/median pool in the per-pattern/per-cluster
   case). No `<v>_by_fs_<v>` column → `ld` is a data frame of `1`s (D3).
   - Reuse the group-column resolution already in the function (:1549-1560); factor the SE and
     loading pooling into one shared helper so both use the identical group order (avoids a
     silent row-order mismatch).
2. **Capture at the call site (R/tspa.R:528-536).** When `derive_sf_se_fs()` returns a non-`NULL`
   `se`, set `se_fs <- derived$se` **and** `sf_ld <- derived$ld` (a new local; `NULL` when the
   multi-factor path was used).
3. **Thread to the schema.** `tspa_sf()` (R/tspa.R:1870; called at :819) gains an `ld = NULL`
   argument (default → unit loading, preserving every existing caller) and forwards it to
   `tspa_schema_sf()`.
4. **`tspa_schema_sf()` (R/tspa.R:1604) accepts `ld`.** Default `ld` = a unit matrix sized to
   `se` (so the existing `tspa_sf(model, data, se, prod_ecov)` signature and all its callers are
   unchanged). The loading row becomes
   `tspa_row(var[k], "=~", fs[k], ld[g, k], g, "struct", lab)` — replacing the literal `1` at
   :1613. (The `error_var`/`error_cov`/`intercept` rows are untouched, D4.)
5. **`tspa_args` replay (refined for the derive-only scope).** The resolved `ld` is NOT
    recorded in `tspa_args`; `tspa_args` records the *derived* `se_fs`, so a replay via
    `do.call(tspa, attr(fit,"tspa_args"))` re-passes that `se_fs` as explicit and therefore
    falls back to the unit loading (it does not re-derive the loading). Documented edge: the
    self-contained-replay guarantee does not extend to a derived single-factor fit.
    `vcov_corrected()` is unaffected (it requires the multi-factor `fsT`/`fsL` path).

## 6. Tests (`r-tester`, edition 3, `tests/testthat/`)

New `test-tspa_cbind_loading.R`:

- **T1 — The fix (A/B).** lme4 `method="EB"` random-intercept fit; the structural coefficient
  from `tspa("d ~ u0", data = cbind(fs, d=d))` (derived, no explicit args) equals the one from
  the **explicit** `tspa("d ~ u0", …, fsL=attr(fs,"fsL"), fsT=attr(fs,"fsT"))` to
  `tolerance = 1e-8` (bit-exact, both reach the same MLE).
- **T2 — Regression guard (the footgun).** Same fit; the derived structural coefficient is
  within, say, 0.15 SE of the single-stage SEM reference (`sem("u =~ 1*y1 + … + 1*ym; d ~ u")`),
  i.e. the +21–28% unit-loading bias is gone. Pin the corrected value.
- **T3 — Latent-variance consistency.** `tspa("u0 ~~ u0", data = cbind(fs, d=d))` now recovers
  the lme4 random-effect variance to `tolerance = 1e-6` (was `w²·ψ` before).
- **T4 — No-op for unit implied loading.** A lavaan CFA single-factor result whose
  `<v>_by_fs_<v>` is `1`: derived fit bit-equals the pre-fix (unit-loading) fit — confirms D3
  (no regression for scores that already have a unit loading).
- **T5 — Hand-rolled frame preserved.** A hand-rolled `data.frame(fs_v = …, fs_v_se = …)` with
  **no** `<v>_by_fs_<v>` column → derived fit equals the explicit `se_fs` unit-loading fit
  (today's behavior, D3).
- **T6 — Precedence / no double application.** Explicit `se_fs` still suppresses derivation
  (PLAN 13 D3); explicit `fsL`/`fsT` still win (D4/multi-factor); a multi-factor *attribute*
  result (no cbind) is untouched by this change (A/B vs current).
- **T7 — Per-row / multigroup.** An unbalanced or multigroup `get_fs()` result cbind'd with a
  distal: the pooled loading uses the same `reduce` order as the pooled SE (group-order canary,
  cf. PLAN 13's character/factor level-order checks).

## 7. Docs (`r-doc`)

- Roxygen `@details` (R/tspa.R derivation subsection): one paragraph — a cbind'd `get_fs()`
  result (attributes dropped) now recovers the per-latent loading from the `<v>_by_fs_<v>`
  column, so a single-factor structural fit matches the explicit `fsL`/`fsT` form; without that
  column the loading stays `1`.
- `NEWS.md` bullet (bug fix): single-factor 2S-PA on a cbind'd shrinkage-score `get_fs()` result
  no longer silently uses a unit loading (which biased the structural coefficient).
- Vignette pointer: one line in `vignettes/R2spa.Rmd` (the canonical cbind flow) noting the
  loading is now recovered automatically.

## 8. Verification / exit criteria

- `devtools::document()` → `devtools::test()` (0 fail; the 7 new blocks green; no existing
  single-factor golden shifts except the intentional cbind'd-shrinkage-score correction, which is
  A/B-gated against the explicit-`fsL` form in T1/T2).
- `devtools::check()` (as-cran, `--no-manual`) 0/0/0.
- Cross-branch: cherry-pick to `get_fs-brms`; the identical T1/T2 pass there with a **brmsfit**
  `method="EB"` result (closing the loop on the original brmsfit observation).
