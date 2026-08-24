# PLAN 13 — `tspa()` auto-derives measurement inputs from a `get_fs()` result

**Date:** 2026-08-24
**Owner roles:** `r-architect` (code), `r-tester` (tests), `r-doc` (roxygen/NEWS/vignettes)
**Status:** draft
**Blocked by / relates to:** none. Builds on PLAN 09 (per-unit pooling: `is_per_unit_fs()` /
`pool_per_unit()`) and the self-contained `tspa_args` replay record (7f704dc).

---

## 1. Problem

The canonical 2S-PA usage is `tspa(model, data = get_fs(...))`, but today the caller must
re-supply the measurement inputs that `get_fs()` already attached to the data frame:
`fsT = attr(fs, "fsT"), fsL = attr(fs, "fsL")` (multi-factor) or `se_fs = list(...)`
(single-factor). `tspa()` ignores those attributes on `data`.

Verified runtime facts (lavaan 0.7-2, R 4.6.1):

- `tspa(model, fs_dat)` on cbind'd single-factor data: `tspa_sf()` returns invisible `NULL` on
  an empty `se_fs` (R/tspa.R:912-916) → lavaan errors `"model is NULL or not a valid type for
  it!"`.
- `tspa(model, fs2)` on a unified multi-factor `get_fs()` result: `fsT`/`fsL` are `NULL`,
  `se_fs` empty → same confusing lavaan error.
- Same for FIML and merMod results.

The boilerplate appears in every multi-factor example: R/tspa.R:138-143, :151-157, :166-173,
:214-227; vignettes/multiple-factors.Rmd.

## 2. Why this is cheap (current state)

- `get_fs()` already attaches everything: `fsT`/`fsL`/`fsb` (R/get_fscore_math.R:329-348,
  R/get_fscore.R:437-440, :444-592; R/get_fs_methods.R:216-238, :382-383, :450-456, :789-806,
  :1066-1072, :1263-1269) plus `fs_pattern`, `group_col`, the `mirt_per_obs` marker,
  `scoring_matrix`, and per-row `fs_<v>_se` columns.
- **The per-unit path already treats the data's attributes as the source of truth.**
  `pool_per_unit()` (R/tspa.R:471-620) calls `resolve_fs_per_row(data)` (R/fs_indiv.R:217-240)
  and reads `attr(data, "fsT"/"fsL"/"fsb"/"fs_pattern"/"group_col"/"mirt_per_obs")`. The explicit
  `fsT`/`fsL` arguments serve only as detection gates (`is_per_unit_fs()`, R/tspa.R:281-282,
  :448-460) and a `have_int` flag (R/tspa.R:283). For FIML/merMod/mirt shapes, derivation
  changes nothing semantically.
- `cbind()` drops the attributes (verified) — the canonical single-factor vignette flow
  (vignettes/R2spa.Rmd:119-134) therefore has **no** attributes and must rely on the
  `fs_<v>_se` columns. `se_fs` derivation from columns is essential, not optional.
- Per-group se reduction already exists: `pool_se_fs()` (R/tspa.R:635-654), today invoked only
  when an explicit `se_fs` was supplied (gated on `nrow(se_fs) > 0`, R/tspa.R:349-364).
- `tspa_args` (R/tspa.R:398-410) records the **in-scope resolved values** — replay via
  `do.call(tspa, attr(fit, "tspa_args"))` (consumed by `vcov_corrected()`,
  R/tspa_corrected_se.R:206-239) re-passes them explicitly, so derivation is skipped on replay
  and no double-pooling occurs (a pooled plain matrix fails `is_per_unit_fs()`). No change
  needed there; document it.

Attribute shapes (verified):

| Shape | `fsT`/`fsL` form | Notes |
|---|---|---|
| SG unified | length-1 list named `""` | assembled R/get_fscore.R:451-453, :584-589; already accepted by the schema (R/tspa.R:724-748) |
| MG unified | named list, one matrix per group, stage-2 group order | labels from `@Data@group.label` (R/get_fs_methods.R:216-238); plus `group` column + `group_col` attr (R/get_fscore.R:581; read at R/tspa.R:345) |
| `format = "list"` | same lists as outer-list attributes | `attr(fs, "fsT")` on the named list works (R/get_fscore.R:557-567) |
| merMod | 3-D arrays (per cluster), no `fsb` | R/get_fs_methods.R:652, :789-790 |
| mirt | per-row list of matrices + `mirt_per_obs` marker | R/get_fs_methods.R:1066-1072, :1254-1258, :1263-1269 |

## 3. Scope & decisions (confirmed with the user 2026-08-24)

- **D1 — Silent auto-derivation; no new flag.** Derivation fires only for argument values that
  are `NULL`, and every such call errors today — no currently-working call changes behavior.
  A flag would protect against a regression that does not exist and would permanently add an
  argument + second code path against the package's deliberately minimal surface.
- **D2 — Explicit arguments always win.** Any non-`NULL` `fsT`/`fsL`/`fsb`/`se_fs` suppresses
  the corresponding derivation entirely. (A mismatch warning between args and attributes is a
  follow-up at most, not v1.)
- **D3 — `se_fs` given ⇒ no multi-factor derivation.** A user who has a multi-factor result but
  explicitly passes `se_fs` today gets the single-factor (diagonal, unit-loading) model and
  keeps getting exactly that.
- **D4 — Both forms available, nothing passed ⇒ multi-factor wins** (the richer, statistically
  correct form for composite scores; the data's own attributes define it).
- **D5 — `vfsLT` stays explicit-only** (set only by `get_fs(..., vfsLT = TRUE)`,
  R/get_fs_methods.R:450-456; carries `which_free` sub-matrixing semantics). `group=` stays
  **explicit-only** for multigroup (no auto-forwarding from the group column in v1).
- **D6 — The SF latent set comes from the data columns, not the model string.** Verified
  (lavaan 0.7-2): `lavaan::lavaanify("visual ~ c(b1, b1) * speed")` **errors** although
  `sem()` accepts the same syntax — and that exact syntax is documented in this package's own
  examples (R/tspa.R:196-200). Model-string parsing is not viable; derive from the package's
  own column-naming convention (AGENTS.md: downstream functions parse by name).
- **D7 — Product-score columns (`fs_a:fs_b`) are out of scope** for derivation (`get_fs_int` is
  quarantined); column derivation is restricted to simple `fs_<v>` / `fs_<v>_se` pairs.

## 4. Approach (r-architect)

All in `R/tspa.R`. No new exports, no `DESCRIPTION` change, no schema/renderer/pooling changes
(frozen/byte-gated in `tests/testthat/test-tspa_render.R`). Helpers stay co-located (AGENTS.md).

1. **Capture `se_fs_given <- !is.null(se_fs)` BEFORE the coercion** at R/tspa.R:253-255
   (after coercion, a missing `se_fs` is a 0×0 data frame, not `NULL`). Placement: after the
   `reliability`/`se` guards (:245-251) and the `xor` both-or-neither check (:256-258).
2. **Derivation A (multi-factor)** — fires iff `fsT` is `NULL` **and** `fsL` is `NULL` **and**
   `!se_fs_given`:
   - Require both `attr(data, "fsT")` and `attr(data, "fsL")` non-`NULL` (data may be a unified
     data frame *or* a named list — outer-list attributes work for both, verified).
   - **Provenance gate:** `resolve_fs_per_row(data)` (R/fs_indiv.R:217) must succeed. Reuse its
     informative errors (missing `fs_pattern` :314-321; list input :401-431; merMod 3-D /
     row-count :466-485; mirt :249-308) instead of duplicating shape rules. A hand-rolled df
     with plain `fsT`/`fsL` matrix attributes but no `fs_pattern` (and no 3-D/mirt shape) is
     therefore **not** derived — it falls through to the improved error (§4.4).
   - `fsT <- attr(data, "fsT"); fsL <- attr(data, "fsL"); fsb <- attr(data, "fsb")` (may be
     `NULL` — merMod has none).
   - Everything downstream is untouched: per-unit detection + pooling (:281-309 — already reads
     the data), name check (:316-330), `group=` check (:332-334), schema, fit attributes
     (:391-397), `tspa_args` (:398-410).
3. **Derivation B (single-factor)** — fires iff after A, `fsT` is still `NULL` **and**
   `!se_fs_given`:
   - Latent set: `v` for each column `fs_<v>` that also has a matching numeric `fs_<v>_se`
     column (simple names only; `fs_a:fs_b` excluded per D7). Order = first appearance of the
     score columns (equals cbind order in the vignette flows — the order users write today,
     vignettes/R2spa.Rmd:130-132, :179-183).
   - Group handling: reuse the `group_col` logic at R/tspa.R:345-348. No group column →
     single-row `se_fs` (a named numeric vector is fine; the coercion at :253-255 normalizes
     it). Group column → one row per group: constant within the group ⇒ the constant; varying
     ⇒ `pool_se_fs()` reduction by `reduce` (idempotent with the existing FIML pooling at
     :349-364, which then re-derives the same values).
   - Group order: factor **levels** for factors, first appearance for character vectors
     (mirrors :580-595). Align `pool_se_fs()`'s current `unique()` (R/tspa.R:643) to that
     convention as part of this change (its only caller, :361, currently fires only for
     character group columns in practice — latent-bug-proofing, not a behavior fix).
   - A new co-located helper next to `pool_se_fs()` builds the derived `se_fs` data frame
     (~40 lines, base R only).
4. **Fail-fast error** — if after both derivations `fsT` is `NULL` and `se_fs` is empty:
   `stop()` with a clear message naming the three accepted forms — explicit `se_fs`
   (single-factor), explicit `fsT`/`fsL` (multi-factor), or a `get_fs()` result as `data`.
   Replaces lavaan's "model is NULL or not a valid type for it!".
5. **Validation that the design adds on top of existing checks:**
   - Provenance gate → hand-rolled attributes rejected with informative errors.
   - Existing name check (:316-330) runs on derived values (always consistent — both sides come
     from the same frame) and on explicit values (unchanged).
   - Latent-name mismatch between the model string and the data fails in lavaan with its
     standard "missing observed variables" error — same contract as explicit args today; no
     new pre-check (model parsing not viable, D6).

## 5. Test plan (r-tester)

New `tests/testthat/test-tspa_derived.R` (self-contained fixtures; repo testthat-3 conventions;
fixture patterns noted per case):

1. **SG MF unified** (length-1 `""`-list attrs): derived fit ≡ explicit-args fit — identical
   `tspaModel` string, `coef()`/`vcov()` equal, fit attributes equal.
2. **SG MF list format** (plain df, direct attributes).
3. **MG MF** unified and list format, `group = "school"` (pattern: test-tspa.R:306-350).
4. **SG SF from `cbind()`** (no attributes): derived `se_fs` values equal the per-column
   constants; fit ≡ explicit-`se_fs` fit (pattern: test-tspa.R:44-53).
5. **MG SF from cbind'd MG frames**: per-group derived `se_fs` ≡ explicit per-group
   `data.frame` (pattern: test-tspa.R:246-255).
6. **FIML SG + MG MF**: derived ≡ explicit (reuse the hand-computed weighted-mean reference
   from test-tspa_pooled.R:78-87); replay via `do.call(tspa, attr(fit, "tspa_args"))`
   reproduces the estimates/vcov (the replay passes resolved plain values, so `is_per_unit_fs()`
   is false on replay — assert estimates/vcov identity, not the `pooled_fs` marker).
7. **merMod (SG)** and **mirt (SG + MG, `skip_if_not_installed("mirt")`)**: derived ≡ explicit
   (pattern: test-tspa_pooled.R:510-554).
8. **Precedence:** attributes present + *differing* explicit args ⇒ fit equals the
   explicit-args-only fit (strip the attributes from a copy for the control).
9. **`se_fs`-given suppression (D3):** multi-factor frame + explicit `se_fs` ⇒ SF path,
   identical to today's explicit behavior.
10. **Provenance gate:** hand-rolled df with `fsT`/`fsL` attrs but no `fs_pattern` ⇒ informative
    error naming explicit `fsT`/`fsL` (not "model is NULL"); `fsT` without `fsL` ⇒ informative
    error; nothing available ⇒ the new fail-fast error.
11. **`corrected_se = TRUE`** with derived `fsT`/`fsL` + explicit `vfsLT`: equals the
    explicit-args corrected fit (SG + MG); the `tspa_args` of a derived fit are fully resolved
    (non-`NULL` `fsT`/`fsL`; `se_fs` where applicable) and self-contained.
12. **`model = ""`** with derived attrs (the efa-score pattern, vignettes/efa-score.Rmd:116-122)
    still works — confirms no model-string parsing anywhere.

**Existing tests that must stay green unchanged:** all of `test-tspa.R`, `test-tspa_render.R`
(byte-identical renderer gate), `test-tspa_pooled.R`, `test-tspa_corrected_se.R`,
`test-tspa_mx.R`, and the `get_fs`/`fs_indiv` suites — derivation fires only on `NULL`-arg calls
that error today, and no schema/renderer/pooling internals change.

## 6. Docs (r-doc)

All roxygen in `R/tspa.R` is owned by r-doc; r-doc runs `devtools::document()` after editing
(never hand-edit `man/` or `NAMESPACE`).

- **Roxygen `R/tspa.R`:**
  - `@param data` — note that when `data` is a `get_fs()` result and the corresponding
    arguments are omitted, the measurement inputs are derived from its attributes
    (multi-factor) or its `fs_<v>`/`fs_<v>_se` columns (single-factor).
  - `@param se_fs` / `@param fsT` / `@param fsL` / `@param fsb` — precedence note (explicit
    wins; the derivation conditions per D1–D4).
  - `@details` — a derivation subsection: the two derivations, the provenance gate, and the
    `cbind()` caveat (a cbind'd frame has no attributes, so its `se_fs` is derived from the
    `fs_<v>_se` columns).
  - `@return` — note that on a derived fit `tspa_args` carries the **resolved** values (replay
    via `do.call(tspa, ...)` / `vcov_corrected()` re-passes them explicitly; no double-pooling).
  - `@examples` — add the minimal forms (`tspa(model, fs)`; MG: `tspa(model, fs, group = ...)`)
    alongside the existing explicit-args examples (keep those — explicit usage remains fully
    supported).
- **Vignettes** (after tests are green; re-knit at check time):
  - `vignettes/R2spa.Rmd` + `vignettes/multiple-factors.Rmd`: drop the
    `fsT = attr(...), fsL = attr(...)` boilerplate from the canonical multi-factor examples.
  - `multilevel`, `tspa-vignette-mx`, `corrected-se` vignettes: **no change** (explicit args
    remain fully supported; their narratives show the explicit form deliberately).
- **`NEWS.md`** — bullet under `# R2spa 0.0.4` → `## New Features`:
  `tspa()` now accepts a `get_fs()` result as `data` directly: when `fsT`/`fsL` (multi-factor)
  or `se_fs` (single-factor) are omitted they are derived from the result's attributes /
  `fs_<v>_se` columns; explicit arguments always take precedence.

## 7. Acceptance

- `devtools::load_all()` → `devtools::document()` → `devtools::test()` → `devtools::check()`
  all green (0 errors / 0 warnings / 0 NOTEs).
- `tspa("dem60 ~ ind60", fs_dat)` (cbind'd single-factor) and
  `tspa("dem60 ~ ind60\ndem65 ~ ind60 + dem60", fs2)` (unified 3-factor; MG with
  `group = "school"`) run with **no** explicit measurement arguments and equal the
  explicit-args fits (estimates, `vcov()`, `tspaModel` string).
- Every existing explicit-arg call is byte-identical (renderer gate `test-tspa_render.R`;
  pinned fits in `test-tspa.R`).

## 8. Out of scope

- Auto-derivation of `vfsLT` (explicit-only by D5).
- Auto-forwarding `group=` from the data's group column (D5; revisit as a follow-up).
- Product-score (`fs_a:fs_b`, `get_fs_int`) derivation (D7; the quarantined flow keeps the
  explicit `se_fs` + existing alias mechanism, R/tspa.R:930-968).
- `tspa_mx_model()` (separate entrypoint with its own `se_fs`/`fsL`/`fsT` precedence,
  R/tspa_mx.R:156-208; single-group Phase 1).
- A mismatch warning between explicit args and data attributes (follow-up at most).
- `se_fs` derivation on the multi-factor path (the `fsT` diagonal already carries the error
  variances there; `se_fs` is SF-only by construction).

## 9. Verification log

(to be filled during implementation)
