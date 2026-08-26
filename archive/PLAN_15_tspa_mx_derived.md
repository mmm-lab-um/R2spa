# PLAN 15 — `tspa_mx_model()` auto-derives measurement inputs from a `get_fs()` result

**Date:** 2026-08-25
**Owner roles:** `r-architect` (code), `r-tester` (tests), `r-doc` (roxygen/NEWS/vignettes)
**Status:** implemented (2026-08-25) — decisions D1–D8 delivered as specified; see §9.
**Blocked by / relates to:** none. Builds on PLAN 13 (the `tspa()` derivation whose
gate semantics this mirrors) and the per-row `fs_indiv()`/`resolve_fs_per_row()`
machinery (PLAN 07/11).

---

## 1. Problem

`tspa()` accepts `get_fs()` results directly since PLAN 13, but its OpenMx
counterpart `tspa_mx_model()` does not. Today the caller must re-supply the
measurement inputs that `get_fs()` already attached:

- constant-quantity results: `fsL = attr(fs, "fsL"), fsT = attr(fs, "fsT"),
  fsb = attr(fs, "fsb")` (the length-1-list unwrap at R/tspa_mx.R:123-132 was
  written for exactly this shape, but the caller still has to pass it);
- per-row results (mirt, FIML): hand-built **character** definition-variable
  matrices plus a data frame from `fs_indiv(fs, include_intercept = TRUE)` —
  the roxygen example (R/tspa_mx.R:49-68) is the documented route, and it is
  entirely manual.

Verified runtime facts (lavaan 0.7-2, R 4.6.1):

- `tspa_mx_model("dem60 ~ ind60", data = fs)` on a `get_fs()` result errors with
  the **misleading** `"'fsL' rows must be named by the factor-score names."`
  (`tspa_mx_spec()` falls into the else-branch with `fsL = NULL`).
- Explicit attributes fit fine: `tspa_mx_model(..., fsL = attr(fs, "fsL"),
  fsT = attr(fs, "fsT"), fsb = attr(fs, "fsb"))` runs (smoke-verified on the
  PoliticalDemocracy 2-factor CFA).
- Per-row list attributes can't even be passed explicitly —
  `tspa_mx_unwrap()` rejects length-`n` lists with the wrong message
  ("not supported yet (Phase 1 is single-group)" — it is per-row, not
  multigroup).

## 2. Why this is cheap (current state)

- A `get_fs()` result already contains **every** definition-variable column the
  character matrices would reference: `fs_<v>`, `fs_<v>_se`,
  `<latent_j>_by_fs_<v>` (all `q^2`), `ev_fs_<v>`, `ecov_fs_<a>_fs_<b>`. Only
  `int_fs_<v>` needs the `fs_indiv(include_intercept = TRUE)` round-trip — and
  its values are just the per-row `fsb` attribute, so they can be appended
  directly (same names, same values, same NA convention).
- The spec layer already supports mixed fixed/defvar cells
  (`tspa_mx_cells()`, R/tspa_mx.R:134-145) and a character `fsb` vector
  (:192-197); the data-contract guards (:96-115) — score columns present,
  defvar columns present, defvar columns NA-free — already produce the exact
  errors the per-row path needs (e.g. unscored rows).
- `resolve_fs_per_row()` (R/fs_indiv.R:217) + `first_pattern_value()`
  (R/fs_indiv.R:411) are the existing provenance/name-resolution internals
  `tspa()`'s PLAN 13 gate reuses — same here, no new shape logic.
- Attribute shapes (verified; see also PLAN 13 §2 table):

  | Shape | `fsT`/`fsL` form after `tspa_mx_unwrap()` |
  |---|---|
  | lavaan SG complete (incl. `local = TRUE` compact, R/get_fs_methods.R:768-814) | plain matrix |
  | `format = "list"` SG | plain matrix |
  | lavaan SG FIML (joint) | named list, one matrix per observed pattern; rows map via `fs_pattern$label` |
  | `local = TRUE` FIML / mirt SG | per-row list (length `n`) + `per_obs` / `mirt_per_obs` marker |
  | mirt MG / lavaan MG | per-row list + `group` column + `group_col` attr (MG) / per-group list |
  | hand-rolled frame | plain matrices, **no** `fs_pattern` → provenance gate fails |

- No pooling semantics exist in `tspa_mx_model()` and none are wanted: the
  OpenMx route's value is the exact (non-pooled) correction. Constant
  quantities → fixed numbers; non-constant → definition variables. That is the
  entire dispatch.
- **Product-score results** (`get_fs(..., product = )`, branch
  `rejoin/fs-prod`): v1 of `compute_fs_prod()` rejects per-obs inputs, so a
  product result is a plain-matrix (complete data) or per-pattern (FIML)
  result **plus** the extra `fs_a:fs_b`/`_se`/`_ld` columns — it flows through
  the D2/D3 branches unchanged (derivation only reads the `fsT`/`fsL`/`fsb`
  attributes and the standard `_by_`/`ev_`/`ecov_` column names; the product
  columns are ignored). The MG signal is the `group_col` attribute, not list
  length (the same disambiguation `compute_fs_prod.R` uses) — consistent with
  D7. One A/B test pins that a product result derives correctly.

## 3. Scope & decisions (confirmed with the user 2026-08-25)

- **D1 — Explicit wins, silent derivation.** If any of `se_fs`/`fsL`/`fsT`/
  `fsb` is given, current behavior is byte-identical. Derivation fires only
  when all four are `NULL` (PLAN 13 semantics; every such call errors today,
  so no currently-working call changes behavior).
- **D2 — Constant → fixed numeric.** Unwrapped `fsL`/`fsT` plain matrices
  (and `fsb` a plain vector) are passed as fixed numeric matrices. `fsb`
  attribute `NULL` (no mean structure) → default fixed-zero behavior.
- **D3 — Non-constant → definition variables.** Per-row list attributes
  (`per_obs`/`mirt_per_obs` marker) or per-pattern named lists (SG FIML,
  keyed by `fs_pattern$label`) → character matrices referencing the frame's
  own columns. **Intercept:** `fsb` attr `NULL` → omit (fixed 0); per-row /
  per-pattern `fsb` → append `int_fs_<v>` columns to a *working copy* of
  `data` (names/values exactly `fs_indiv(include_intercept = TRUE)`) and
  point `fsb` at them. No pooling, no `reduce`.
- **D4 — Provenance gate, PLAN 13 style.** Derivation only when both
  `attr(data, "fsT")` and `attr(data, "fsL")` are present **and**
  `resolve_fs_per_row(data)` succeeds (reused via `tryCatch` for its
  informative errors). A hand-rolled frame with plain matrix attributes but no
  provenance is **not** derived; it falls to the D5 fail-fast carrying the
  gate's message.
- **D5 — Fail-fast replaces the misleading error.** All-NULL on a
  non-`get_fs()` frame becomes an actionable "no measurement inputs — supply
  `se_fs`, or `fsL` and `fsT`, or pass a `get_fs()` result as `data`" error
  (PLAN 13 wording adapted), instead of `"'fsL' rows must be named..."`.
- **D6 — Existing guards unchanged.** Score-column presence, defvar-column
  presence, and the NA-free defvar guard (R/tspa_mx.R:96-115) run on the
  (possibly `int_fs_`-augmented) frame and give the existing clear errors —
  including for unscored rows in per-row results.
- **D7 — Multigroup: clear refusal.** Length-`>1` list attributes with a
  `group_col` attribute (unified MG; per-row MG carries the same marker) →
  the existing "Multigroup … not supported yet (Phase 1 is single-group)"
  message, emitted by derivation rather than the wrong unwrap error.
  List-format MG (a named list, not a data frame) keeps the existing
  `is.data.frame(data)` error (:76).
- **D8 — No replay concerns.** `tspa_mx_model()` returns a self-contained
  `MxModel` (RAM built from the derived spec); nothing equivalent to
  `tspa_args` is needed.

## 4. Approach (r-architect) — `R/tspa_mx.R` only

1. **Capture** `se_fs_given <- !is.null(se_fs)` and `fsb_given <- !is.null(fsb)`
   before any coercion (mirrors R/tspa.R:349-352).
2. **New internal `tspa_mx_derive_measurement(data)`**, co-located after
   `tspa_mx_unwrap()`, returning `list(fsL, fsT, fsb, data)` (`data` possibly
   augmented with `int_fs_*` columns) or `stop()ping`:
   - attr presence: both `fsT` and `fsL` non-NULL, else signal "no
     attributes" (caller → D5);
   - provenance: `tryCatch({ resolve_fs_per_row(data); TRUE }, error = ...)`
     (D4) — a failure carries the gate message to D5;
   - unwrap both via `tspa_mx_unwrap()`; `fsb` unwrapped when a list;
   - dispatch:
     - both matrices → D2: return as fixed numeric (`fsb` vector or NULL);
     - list with `attr(data, "group_col")` present (or a trailing group
       column + per-row marker, i.e. MG per-row) → D7 stop;
     - per-row (marker present) → D3 per-row: names from
       `first_pattern_value(attr_L)` dimnames (row = score, col = latent);
       `L_char[s, v] <- paste0(v, "_by_", s)`; `T_char` lower-tri (incl. diag)
       `ev_`/`ecov_` names, upper-tri NA (the spec uses the lower triangle);
       `fsb` NULL → omit; per-row `fsb` list → append one `int_fs_<s>` column
       per row to a **copy** of `data` (all-NA elements stay NA → D6 guard)
       and return character `fsb`.
     - per-pattern named list (no marker; SG FIML) → D3 per-pattern: same
       `L_char`/`T_char` (columns are pattern-constant, hence complete);
       `fsb` per-pattern named list → append `int_fs_*` by mapping rows to
       `attr(data, "fs_pattern")$label` values.
   - column-existence pre-check is **not** needed separately: the existing
     dv_cols guard (:104-108) reports missing columns with the standard
     message.
3. **Wire into `tspa_mx_model()`** between the existing arg-guards (:79-90)
   and `tspa_mx_spec()` (:92): when `!se_fs_given && is.null(fsL) &&
   is.null(fsT) && !fsb_given` → call the helper; on gate failure store the
   message; then the D5 fail-fast fires whenever no measurement inputs ended
   up (replacing the fall-through to `tspa_mx_spec()`'s misleading error).
   The existing guards (:79-90) and spec/`lav_to_mx_ram` are untouched.

## 5. Test plan (r-tester) — `tests/testthat/test-tspa_mx_derived.R`

A/B pattern (derived ≡ explicit, `coef()` + covariance at tight tolerance),
mirroring `test-tspa_derived.R`:

1. **SG complete, 2-factor joint CFA** (PoliticalDemocracy): derived ≡
   explicit `attr()` matrices (D2).
2. **`local = TRUE` compact** (3-factor, R2spa.Rmd model): derived ≡
   explicit attributes; and derived ≡ manual character-matrix route (block-
   diagonal, zero cross-terms) — pins D2/D3 agreement on the same data.
3. **SG FIML joint** (`local = FALSE`, NAs injected): per-pattern defvar
   route ≡ manual character-matrix + `fs_indiv()`-data route (D3 per-pattern).
   **Narrowed (implementation finding V3c):** full per-row/per-pattern
   q≥2 defvar models abort in OpenMx 2.22.11 ("implied covariance not
   positive definite") for both the derived **and** the documented manual
   route (off-diagonal defvar cells trigger it; seeds 7/42/99/123/1334,
   SLSQP and CSOLNP) — a pre-existing OpenMx RAM limitation, not a
   derivation defect. Pin instead: derived model **string** ≡ manual model
   string (identical), appended `int_fs_*` ≡ `fs_indiv()` values, and a
   q=1 / diagonal-only per-row fit as the end-to-end numerical A/B.
4. **Local FIML** (`per_obs`) and **mirt SG** (`mirt_per_obs`): per-row
   defvar ≡ manual route (mirt control = the roxygen-example matrices).
   **mirt MG documented behavior:** a `MultipleGroupClass` result carries no
   `group_col` attribute (only a `group` column), so D7 does not fire — it
   derives as a single **pooled per-row-corrected** fit (the exact per-row
   measurement quantities; no per-group structural parameters; the `group`
   column is inert). Pin with one A/B (≡ explicit character matrices).
5. **Intercepts:** mean-structure CFA (nonzero constant `fsb`) → numeric path
   ≡ explicit `fsb`; mirt `prior_mean = c(F1 = 2)` → `int_fs_*` appended,
   derived ≡ explicit `fsb = c(fs_F1 = "int_fs_F1")` on
   `fs_indiv(fs, include_intercept = TRUE)` data (D3 intercept).
6. **Explicit-wins (D1):** `se_fs` given → se_fs path (no derivation); one of
   `fsL`/`fsT` → existing xor error; fully explicit → unchanged fit.
7. **Fail-fast (D4/D5/D7):** hand-rolled frame (plain matrix attrs, no
   `fs_pattern`) → D5 message incl. the gate text; plain non-`get_fs` frame
   all-NULL → D5; unified MG result → D7 message.
8. **NA handling (D6):** per-row result with a completely-missing row →
   existing "Definition-variable column(s) contain NA" error.
9. **Product-score edge:** `get_fs(..., product = "ind60:dem60")` result
   (complete data) → derived ≡ explicit-attribute control (product columns
   inert to derivation).
10. **Regression:** `test-tspa_mx.R` untouched and green; current suite
    baseline (re-measured at start — product-score work landed after 3720)
    preserved.

## 6. Docs (r-doc)

1. **Roxygen `tspa_mx_model`:** new `@details` subsection "Auto-derivation
   from `get_fs()` results" (D1–D7 in user terms; emphasize *no pooling* —
   per-row quantities go through definition variables, constant ones are
   fixed); `@param data` — a `get_fs()` result works directly;
   `fs_indiv()` no longer needed for the intercepts; `@param se_fs`/`fsL`/
   `fsT`/`fsb` — "or omitted to derive from a `get_fs()` result";
   `@examples` — the existing `\dontrun` example reduced to the direct call
   plus one explicit-attribute contrast. Then `devtools::document()`.
2. **`NEWS.md`:** 0.0.4 New Features entry.
3. **`vignettes/tspa-vignette-mx.Rmd`:** "Factor-score intercepts via
   `get_fs()`" section extended to the direct
   `tspa_mx_model("f2 ~ f1\nf2 + f1 ~ 1", data = fs)` call (zero- and
   nonzero-`fsb` via `prior_mean`); EAP section: note the hand-rolled
   `fs_eap` + character matrices can be replaced by passing the `get_fs()`
   result directly (keep the hand-rolled narrative as the pedagogical base).
   Re-knit.
4. **Optional:** one-line pointer in `vignettes/R2spa.Rmd` (mx line of
   "Related vignettes" or "How the measurement error is supplied") if it
   reads naturally.

## 7. Acceptance

- `tspa_mx_model(model, fs)` runs with zero measurement args on: lavaan SG
  complete (joint and `local = TRUE`), SG FIML (joint and local), mirt SG
  (default and `prior_mean`); results ≡ explicit-args A/B controls.
- All D1–D8 decisions pinned by `test-tspa_mx_derived.R`.
- Suite green (baseline 3720/0, 1 pre-existing warn); `R CMD check` 0/0/0.
- Roxygen/NEWS/vignette updated; `document()` clean; vignettes re-knit.

## 8. Out of scope

- Multigroup (the Phase-1 guard stands; D7 only improves the message path).
- `se_fs`-column derivation (the OpenMx `se_fs` path is constant-by-design;
  per-row data uses the defvar route).
- cbind'd frames without attributes (pass the `get_fs()` result itself).
- Pooling / `reduce` of per-unit quantities in the OpenMx route — by design
  (exact or fail).

## 9. Verification log (closed 2026-08-25)

**Hook.** R/tspa_mx.R (purely additive, +198 lines): `se_fs_given`/`fsb_given`
captured at the top of `tspa_mx_model()` (D1); derivation + D5 fail-fast wired
between the arg-guards and `tspa_mx_spec()`; new internal
`tspa_mx_derive_measurement()` co-located after `tspa_mx_unwrap()`, returning
`list(fsL, fsT, fsb, data, prov_err)` and dispatching: matrices → D2;
`group_col` attr → D7 stop; `per_obs`/`mirt_per_obs` → D3 per-row (character
matrices + `int_fs_*` appended to a working copy); `fs_pattern$label` → D3
per-pattern (row-wise label mapping).

**Deviations (all minor, all in the helper's internals, none in the contract):**

1. **Length-1 unwrap is local to the helper.** Unified results wrap *every*
   attribute (incl. `fs_pattern`) in a `list("" = …)`; the helper unwraps
   length-1 lists itself and dispatches longer lists only on the `group_col`
   marker, so `tspa_mx_unwrap()`'s multigroup message is preserved for
   genuine explicit MG inputs.
2. **`first_pattern_value()` is called on the already length-1-unwrapped
   `L`/`T`** (a unified SG attribute is double-wrapped: `""` then the
   per-row/per-pattern list).
3. **mirt MG derives as pooled per-row (not refused).** A
   `MultipleGroupClass` result carries no `group_col` attribute (only a
   `group` column), so D7 does not fire; it fits as a single pooled
   per-row-corrected model (exact per-row measurement quantities, no
   per-group structural parameters, `group` column inert). Documented in the
   roxygen `@details` and pinned by a dedicated A/B test.

**Findings (pre-existing, not derivation defects):**

- **V3c — OpenMx 2.22.11 q≥2 off-diagonal-defvar abort.** Full per-row/
  per-pattern q≥2 defvar models abort with "implied covariance not positive
  definite" for BOTH the derived and the documented manual character-matrix
  routes (seeds 7/42/99/123/1334, SLSQP and CSOLNP; HS 2/3-factor and
  synthetic; diagonal-only controls fit fine). Test item 3 is pinned at
  string identity + `int_fs_*` identity + a q=1 end-to-end A/B accordingly;
  the vignette demonstrates the mirt derived route on q=1 results for the
  same reason.
- **Product frames are un-fittable by OpenMx as-is**: `mxData()` rejects the
  `:` in `fs_a:fs_b` column names ("is illegal because it contains the ':'
  character") — the explicit-argument route fails identically (byte-identical
  message). Pinned: both routes error the same; with the product columns
  dropped, derived ≡ explicit bit-exact and the model string is unchanged
  with/without them (inertness). One-line roxygen note added.
- **Fixed numeric `fsb` is a no-op in the RAM model** (pre-existing, by
  design — R/tspa_mx.R comment: a fixed/absent mean leaves the observed
  score mean at its data value; only *character* `fsb` emits a score-mean
  line). Item 5a pins the D2 numeric dispatch + A/B equivalence, not a mean
  constraint.

**A/B results (r-architect verification, all bit-exact, max|Δcoef| = 0):**
V1 SG complete 2-factor (D2) · V2 `local = TRUE` compact 3-factor (D2) ·
V3 SG FIML joint (structure: model string identical; `int_fs_*` ≡
`fs_indiv()`; q=1 numerical A/B) · V4 mirt SG per-row (default and
`prior_mean = c(F1 = 2)` auto-append) · V5 error paths (D5 + gate text, D5
base, D7 MG message, D6 NA guard) · V6 explicit-wins (se_fs path, xor
error) — plus r-tester's 70 expectations (13 blocks): items 1–9 of §5 all
pass, A/Bs pinned with `expect_identical` on coefficients and 1e-10
tolerances on path/variance values.

**Gate.** `devtools::test()`: **FAIL 0 | WARN 1 | PASS 3940** (baseline
before the plan: 3870 / 0 fail; the 1 warn is the pre-existing
`test-tspa_mx.R` "some estimated lv variances are negative").
`R CMD check`: **0 errors / 0 warnings / 0 notes** (3m 19s, vignettes
re-knit). `document()` clean; NAMESPACE unchanged; only
`man/tspa_mx_model.Rd` regenerated.
