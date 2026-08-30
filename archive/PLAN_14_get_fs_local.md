# PLAN 14 — Per-construct ("local") stage-1 scoring: `get_fs(..., local = TRUE)`

**Date:** 2026-08-24
**Owner roles:** `r-architect` (code), `r-tester` (tests), `r-doc` (roxygen/NEWS/vignettes)
**Status:** implemented (2026-08-24) — decisions D1–D6 delivered as specified; see §9.
**Blocked by / relates to:** none (benefits from PLAN 13: a local-mode result feeds
`tspa(model, fs)` directly with no explicit arguments once the `per_obs` marker OR-checks
land).

---

## 1. Problem

Two ways to get multiple factor-score sets today:

- **(a) Joint:** `get_fs(data, model = <multi-factor CFA>)` — one joint fit; all scores
  computed from it. With the default freely-correlated factors, each score is a weighted
  composite of **all** indicators (cross-factor items leak into every score — verified:
  joint `fsL` off-diagonals 0.0059/0.1818, joint `fsT` off-diagonal 0.0056 on the 3-factor
  `PoliticalDemocracy` example), and the per-factor estimates themselves shift under free
  correlation (e.g. `psi_ind60` 0.4485 joint vs 0.4455 local).
- **(b) Canonical 2S-PA (Lai & Hsiao 2022):** call `get_fs()` once per single-factor model
  and `cbind()` the results — the pattern the vignettes teach. But `cbind()` **drops the
  `fsT`/`fsL` attributes**, which locks cbind-ed users into the single-factor `se_fs`
  stage-2 path (R/tspa.R:336-384): they can never reach the multi-factor stage-2 machinery
  (cross-loading structure, `corrected_se`, `tspa_mx_model()`'s `fsL` argument).

The requested feature: one call that scores each latent from its **own local model**
(analogous to lavaan local models / `lavaan::sam()`-style per-model estimation), restricted
to the measurement stage. This makes the canonical pattern (b) a single call with full
multi-factor output.

## 2. Why this is feasible (current state)

All new code composes the existing per-fit machinery one fit at a time; no new dependencies.

| Stage | Location | Role |
|---|---|---|
| Entry point | R/get_fscore.R:99-132 (`get_fs.data.frame`) | Only the data.frame method owns the model string; `local` belongs here |
| Fit → blocks | R/get_fs_methods.R:263-498 (`get_fs.lavaan`); :134-239 (`get_fs_blocks.lavaan`); :148-207 (`prepare_fs`, per-missing-pattern blocks) | Each local fit is an ordinary single-factor `lavaan` object flowing through this unchanged code |
| Math | R/get_fscore_math.R:268-351 (`compute_fscore`; `fsL` :331, `fsT = a Θ aᵀ` :346, `fsb` :339-343); :467-472 (regression `a`); :474-494 (Bartlett) | Pure per-fit math; a 1-factor local fit gives a `1×p_k` `a` row automatically |
| Assembly / naming | R/get_fscore.R:400-442 (`augment_fs` — the column-name spec); :444-592 (`assemble_fs_blocks`); :287-398 (`fs_to_group_list`) | The merged object must emit exactly this layout so `fs_to_group_list()` round-trips unchanged |
| SE paths | R/get_fs_methods.R:316-330 (missing-data guard for `corrected_fsT`/`reliability`/`vfsLT`); :361-371 (`correct_evfs`); :450-456 (`vfsLT`) | Per-fit; work on each local fit independently **except `vfsLT`** (§3 D4) |
| Downstream (name/attr-parsing spec) | R/tspa.R:336-384 (sf vs mf branch) + PLAN 13 derivation; R/fs_indiv.R:217-240 (`resolve_fs_per_row`); R/tspa_mx.R | A merged object matching the joint layout needs **zero** changes here except the one-line marker generalizations (§4.5) |

Verified facts (lavaan 0.7-2, R 4.6.1):

- FIML: row order is preserved across separate fits on the same data (301/301 rows, same
  group labels and per-group counts in MG); each local fit has its **own** pattern set
  (over different item sets) — pattern sets are not comparable across latents, so a
  per-pattern *merged* attribute is impossible; **per-row merging is the only
  FIML-consistent representation**. Rows with all of a factor's items missing get NA
  scores (the existing all-NA block convention).
- `corrected_fsT` and `reliability` are computable per single-factor local fit
  (local `reliability` = 0.9607411, matching the existing SG pin in
  test-get_fscore.R:1087-1092).
- Joint vs local scores: free correlation ⇒ differ (max |diff| 0.369 regression / 0.242
  Bartlett); factors constrained uncorrelated (`ind60 ~~ 0 * dem60`) ⇒ agree to ~1.4e-7
  (the likelihood factorizes, so per-factor MLEs coincide).
- `cfa()` accepts structural `~` paths, so today's joint mode tolerates them; a fitted
  `lavaan` object does **not** retain its syntax string → `local` can only work from the
  data-frame entry point.
- The merged free elements for `vfsLT` would need **cross-latent** sampling covariances
  (estimates from separate fits on the same observations are correlated) — a block-diagonal
  assembly of per-local `vfsLT`s would be **wrong**, not just incomplete.

## 3. Scope & decisions (confirmed with the user 2026-08-24: full scope)

- **D1 — Hybrid input: strict-grammar string split + parse-free vector/list escape hatch.**
  `local = TRUE` on the data-frame entry point, `model`:
  - a **single string** → strict-grammar split into per-latent models (§4.2);
  - a **character vector of length ≥ 2, or a named list of strings** → each element is a
    *complete single-factor model string*, fit verbatim as `cfa(element, data, group, ...)`
    (no parsing at all); each element must define exactly one latent; latent order = vector
    order; latent names must be unique across elements. This form accepts any lavaan syntax
    valid for a 1-factor CFA (within-factor residual covariances, fixed values, thresholds
    where the backend supports them) — the documented escape hatch;
  - `model = NULL` → the auto `f1 =~ <all items>` model (R/get_fscore.R:113-119) is
    trivially local; run the normal single-fit path (no-op).
  - `local = TRUE` on a `lavaan`/`merMod`/`mirt` object → **error** (no syntax string to
    split; a fitted object is a joint model by construction).
- **D2 — Merged output layout is byte-identical to the joint layout** (same columns, same
  order, same attribute shapes) with **exactly-zero** cross-terms (`_by_` off-diagonals,
  `ecov_*`, `fsT`/`fsL`/`psi` off-diagonals) encoding "no shared measurement model".
  Downstream (`tspa()`, `fs_indiv()`, `fs_to_group_list()`, `tspa_mx_model()`) is untouched
  except the one-line marker generalizations (§4.5).
- **D3 — `per_obs = TRUE` marker** on FIML local results (per-row attribute lists, same
  convention as mirt's `mirt_per_obs`); the three consumption sites OR-check
  `per_obs || mirt_per_obs`.
- **D4 — `vfsLT = TRUE` is rejected** with a documented message (no cross-latent sampling
  covariances across separate fits ⇒ no `tspa(corrected_se = TRUE)`, no corrected
  grand-standardized SEs from a local stage 1). `prior_cov` rejected in v1 (a q×q prior
  cannot be reduced to per-local 1×1 priors without silently dropping off-diagonals;
  follow-up: per-latent priors). `reliability = TRUE` rejected in v1 (per-latent attribute
  shape is new contract surface for an attribute `tspa()` already deprecates; follow-up).
- **D5 — Supported argument interactions:** `group` (forwarded to every local `cfa()`),
  `std.lv` (forwarded; joint `std.lv = TRUE` ≡ per-local `std.lv = TRUE`), `missing` and
  all other `...` (forwarded), `method` (regression/Bartlett/mean all work per local fit),
  `sum_items` (auto-derivation is trivially satisfied per local model; a user-supplied
  `sum_items` is sliced per latent), `corrected_fsT = TRUE` (per-fit `correct_evfs()`
  needs no cross-latent quantity), `prior_mean` (validated once against all latent names,
  then sliced per local fit).
- **D6 — Cross-factor structure is *not estimated by design*.** The merged `psi`/`fsT`/`fsL`
  off-diagonals are zeros carrying the independence assumption, documented as such.

## 4. Approach (r-architect)

All new code co-located in `R/get_fscore.R` (AGENTS.md helper rule), except the three
one-line marker checks. No `DESCRIPTION` change (parser is base R; `block_diag()` already
exists at R/helper.R:4). No new exports.

1. **`local = FALSE` formal** on `get_fs.data.frame()` (R/get_fscore.R:99-112; placed before
   `...` so it is never forwarded), documented on the generic (R/get_fscore.R:11-118).
   Guards: error on fitted-object input (`local = TRUE` with a `lavaan`/`merMod`/`mirt`
   object); `model = NULL` ⇒ normal single-fit path (trivially local).
2. **`split_local_models(model)`** (new internal) — strict, error-first:
   - strip `#` comments; `;` → statement separator; a line ending in `+` continues;
   - every statement must match `lhs =~ i1 + i2 + ...` where `lhs` is exactly **one** bare
     identifier and the RHS is a `+`-separated list of bare identifiers;
   - anything else is an error **naming the offending line** and pointing to the
     alternatives (joint mode — `local = FALSE` — or the vector form): multi-latent LHS
     (`a b =~ c`), **any** `~~` (latent-latent covariance/correlation, cross-factor
     residual covariance, and within-factor residual covariance — v1 limitation, the
     vector form works around it), `~` structural, `\|~`, thresholds `$`, labels, `c()`,
     fixed values, `0*`;
   - cross-checks: every item appears in exactly one statement; every latent appears once;
     latent order = statement order.
3. **Orchestration** in `get_fs.data.frame()` with `local = TRUE`: for each latent (model
   order), `cfa(model_k, data, group = group, ...)` → the **existing method unchanged**
   `get_fs(fit_k, method = ..., corrected_fsT = ..., format = "list", prior_mean =
   <per-latent slice>, sum_items = <per-latent slice>, ...)`, reusing per fit: the
   missing-data guard (R/get_fs_methods.R:316-330), `prepare_fs` pattern handling
   (:148-207), `augment_fs` naming (R/get_fscore.R:400-442), `assemble_fs_blocks`
   (:444-592). Then `merge_local_fs()` (§4.4).
4. **`merge_local_fs()`** (new internal) — per-row merge, the single representation that is
   consistent under FIML:
   - resolve each per-latent result to **per-row blocks** with `resolve_fs_per_row()`
     (R/fs_indiv.R:217-240): one `(fsL_k, fsT_k, fsb_k)` per row (constant within a group
     when complete; pattern-specific under FIML; all-NA block for unscorable rows);
   - per row, **block-diagonalize** across latents (`block_diag()` from R/helper.R:4, with
     an all-NA 1×1 block where that latent couldn't score the row — the existing NA-row
     convention);
   - when every row of a group carries an identical block (complete data), materialize the
     compact per-group matrix instead of the list;
   - **columns**: order identical to the joint output (R/get_fscore.R:407-432):
     `fs_v1 … fs_vq | fs_v1_se … fs_vq_se | <v_j>_by_fs_v_i (q², latent-outer/score-inner)
     | ev_fs_vi / ecov_fs_vi_fs_vj (q(q+1)/2, i-outer j≤i)`; off-diagonal loading columns
     and all `ecov_*` columns exactly 0;
   - **attributes**: `fsT`/`fsL`/`psi` block-diagonal per group (SG: length-1 list named
     `""` — the existing convention, R/get_fscore.R:178-184; MG: named list keyed by group
     label), off-diagonals 0 (`psi` off-diagonal 0 = the documented independence
     assumption); `fsb` concatenated length-`q` (per group); `scoring_matrix` `q×p` per
     group with row k = local fit k's `a_k` (1×p_k) **zero-padded to all items** (item
     order = first appearance across local models) — `S %*% y` reproduces the merged
     scores exactly; `alpha` concatenated local latent means; `fs_pattern`: complete data →
     `list(label = <all-items label>, pat = p×1 all-TRUE)` per group; FIML →
     `list(label = seq_len(n), pat = NULL)` per group (the mirt per-row convention,
     R/get_fs_methods.R:1069) **+ `per_obs = TRUE` marker**;
   - `format = "unified"`/`"list"` both work (MG adds the trailing `group` column and
     `group_col` attribute, as `assemble_fs_blocks` already does, R/get_fscore.R:571-589).
5. **One-line marker generalizations** (`per_obs || mirt_per_obs`):
   - R/fs_indiv.R:222 (`resolve_fs_per_row` dispatch);
   - R/tspa.R:281-282 (`is_per_unit_fs()` call site in `tspa()` — after PLAN 13 this is
     what lets a FIML local result pool through the derived path);
   - R/tspa.R:480-485 (`pool_per_unit()`'s mirt `mirt_mg` recovery from the data's group
     column).

## 5. Test plan (r-tester)

New `tests/testthat/test-get_fs_local.R` (~15-25 `test_that` units; house conventions from
`test-get_fscore.R` / `test-tspa_pooled.R`):

1. **Parser** (via `R2spa:::split_local_models`): valid splits (comments, `;`, trailing-`+`
   continuation, 3+ latents); each error class from §4.2 (multi-latent LHS, each `~~`
   flavor, `~`, `|~`, `$`, labels, `c()`, fixed values, duplicate item, duplicate latent);
   vector/list form: verbatim fit, one-latent validation, duplicate latent names across
   elements, order preservation.
2. **SG equivalence to the canonical pattern**: local output's score/`_se` columns
   `identical()` to `cbind(get_fs(m_k1), get_fs(m_k2))` (the vignette workflow); merged
   attributes equal block-diagonals of the per-local attributes; `scoring_matrix %*% y`
   reproduces the scores (1e-10, mirroring the merMod identity tests at
   test-get_fscore.R:428-447).
3. **Layout pin**: column set/order equals the joint `get_fs()` layout for the same
   latents; off-diagonal `_by_` and `ecov_` columns exactly 0; `fsT`/`fsL`/`psi`
   block-diagonal with correct dimnames; SG `""`-wrap convention; `fs_to_group_list()`
   round-trip (both directions).
4. **Joint A/B**: with `~~ 0*`-constrained factors, local ≡ joint scores/`fsT`/`fsL` to
   1e-5; with free correlation, assert **not** equal (locks the semantics, guards against
   accidental joint reuse).
5. **MG**: `group` column, per-group block-diagonal attributes equal per-local per-group
   values, row counts; `group.equal` forwarding.
6. **FIML**: per-row list shapes (length = nrow); per-row blocks equal each local fit's
   per-row blocks (cross-checked via `fs_indiv()` on each local output —
   R/fs_indiv.R:77-166); all-NA-row convention; `per_obs` marker; `fs_indiv(local_output)`
   row-equal to the cbind of per-local `fs_indiv()` outputs; downstream pooling
   `tspa(reduce = "mean")` works — and (PLAN 13) `tspa(model, local_fs)` with **no**
   explicit `fsT`/`fsL` works on the FIML local result.
7. **Downstream (complete data)**: `tspa()` mf path on local output fits; `tspaModel`
   carries the zero cross-loading/covariance rows; `se_fs`-only path on the same data;
   `tspa_mx_model()` accepts the merged `fsL` (R/tspa_mx.R:70).
8. **Guards**: `vfsLT`/`prior_cov`/`reliability` errors with the documented messages;
   `local = TRUE` on a fitted `lavaan` object errors; `model = NULL` no-op (≡ plain
   single-factor `get_fs`).
9. **`corrected_fsT = TRUE`**: merged blocks equal per-local corrected values (A/B against
   two separate `corrected_fsT` calls).
10. **`std.lv`**: local + `std.lv` ≡ per-local `std.lv` (diagonal `psi` = 1 per group).

## 6. Docs (r-doc)

- **Roxygen** (r-doc owns all roxygen; runs `devtools::document()` after editing):
  - `@param local` on the generic (R/get_fscore.R:11-118);
  - `@details` subsection "Local per-construct scoring": the canonical 2S-PA stage 1 is
    per-construct models; local mode produces *pure* per-construct scores with exactly-zero
    cross-terms; **local ≠ joint when factors are correlated** (verified numbers: max
    score diff 0.369 regression / 0.242 Bartlett; uncorrelated ⇒ agreement ~1e-7 — the
    bridge case); the `vector`/list escape hatch; the v1 rejections (`vfsLT`, `prior_cov`,
    `reliability`) with one-line reasons;
  - `@return` extension: `per_obs` marker on FIML results (per-row attribute lists);
  - `@examples`: a 3-factor `local = TRUE` call (PoliticalDemocracy) + one vector-form call
    (small/fast; examples run under R CMD check).
- **Vignette** (`vignettes/multiple-factors.Rmd`): new "Local vs joint scores" subsection
  with the verified local-vs-joint numbers, the zero-cross-term semantics, and the
  `tspa(model, local_fs)` one-call flow (post-PLAN 13).
- **`NEWS.md`**: bullet under `# R2spa 0.0.4` → `## New Features`:
  `get_fs()` gains `local = TRUE`: each latent is scored from its own local measurement
  model (per-construct stage 1, the canonical 2S-PA setup) instead of the single joint
  multi-factor model; the merged result carries the usual multi-factor attributes with
  exactly-zero cross-terms (block-diagonal `fsT`/`fsL`, zero `ecov_*` columns) and feeds
  `tspa()` directly. `vfsLT`/`prior_cov`/`reliability` are not supported in `local` mode
  (v1).

## 7. Acceptance

- `devtools::load_all()` → `devtools::document()` → `devtools::test()` →
  `devtools::check()` all green (0 errors / 0 warnings / 0 NOTEs).
- `get_fs(PoliticalDemocracy, model = <3-factor string>, local = TRUE)` works; its result
  equals the per-construct cbind pattern for scores/SEs and feeds
  `tspa("dem60 ~ ind60\ndem65 ~ ind60 + dem60", fs_local)` with **no** explicit
  `fsT`/`fsL` (PLAN 13 derivation), producing a stage-2 model with zero cross-loadings.
- All existing `get_fs()`/`tspa()`/`fs_indiv()` behavior unchanged (full suite green).

## 8. Out of scope

- `local = TRUE` on fitted objects (`lavaan`/`merMod`/`mirt`) — the syntax string is the
  spec; error.
- String form: only single-latent `=~` definitions (bare identifiers, `+` lists, `#`
  comments, `;`, trailing-`+` continuation). **No** `~~` of any kind (v1), no `~`, `|~`,
  thresholds, labels, `c()`, fixed values. The vector/list form is the documented escape
  hatch (any 1-factor lavaan syntax, validated to define exactly one latent).
- **No `vfsLT`** → no `tspa(corrected_se = TRUE)`, no corrected grand-standardized SEs from
  a local stage 1 (D4).
- **No `prior_cov`** (per-latent priors are a follow-up); `prior_mean` is supported.
- **No `reliability` attribute** (computable per local fit; attribute-contract design
  deferred).
- Cross-factor structure is not estimated by design (D6).
- No changes to `augment_lav_predict()`, `get_fs_lavaan()`/`get_fs_lmer()` wrappers, or any
  quarantined code. Quarantine conflict check: none — `.quarantine/R/get_fs_int.R`
  operates on single-factor `get_fs()` results via the same column conventions the merged
  object preserves per-factor; its re-integration is unaffected.

## 9. Verification log (closed 2026-08-24)

**Hook.** `local = FALSE` formal on `get_fs.data.frame()` (R/get_fs_methods.R, before
`...`); fitted-object guard in the generic (R/get_fscore.R); new internals in
R/get_fs_methods.R: `get_fs_local()`, `split_local_models()`, `parse_local_statement()`,
`local_model_syntax_error()`, `merge_local_fs()`; marker OR-checks in R/fs_indiv.R
(`resolve_fs_per_row` dispatch + `resolve_per_obs()` group-col recovery) and R/tspa.R
(`is_per_unit_fs` call site, `pool_per_unit` mirt_mg).

**R1 — FIML `scoring_matrix` placement (found in implementation).** Each local latent's
per-pattern `a`-matrix spans only the pattern's *observed* items; placed at observed-item
positions, `NA` at that latent's missing items, `0` at other latents' items. Verified on
empty-pattern and partial-missing edge cases.

**R2 — pre-existing `resolve_group_blocks` bug (found in implementation, fixed at root).**
The single-pattern (plain-matrix) branch in R/fs_indiv.R routed fully-NA/empty-pattern rows
to the value block instead of an all-NA block (the multi-pattern branch already handled
it). Required for local FIML; benefits `fs_indiv()` for all inputs. Full suite: zero
regressions.

**R3 — parser bugs (found by r-tester, fixed by r-architect).** (a) `local_model_syntax_error()`
had arity `(line, txt, detail)` while five call sites passed 2–3 detail strings → raw R
`unused argument` errors; fixed with `...` + concatenation (the seven formerly-clean classes
stay byte-identical). (b) `parse_local_statement("a =~")`: `strsplit` drops the trailing
empty field → RHS `NA` → raw `missing value where TRUE/FALSE needed`; fixed by treating
missing/`NA`/zero-length RHS as the empty-RHS case. All 16 error classes now emit the full
intended message (stable prefix `"unsupported model syntax on line"` + detail + pointer
sentence), pinned in the tests (14 prefix pins + 7 class-specific detail fragments).

**A/B numbers.** local ≡ joint(uncorrelated, `ind60 ~~ 0 * dem60`) at 1.16e-5 (tests use
tolerance 1e-4); free correlation: max score diff observed ~1.52 in the test (plan §1's
0.369/0.242 are the PoliticalDemocracy 3-factor magnitudes). Scoring-matrix identity exact
for regression/Bartlett: `score = S %*% (y − colMeans(y))` (per-column centering — the naive
row-wise `y − colMeans(y)` is wrong under R's column-major recycling); `method = "mean"` →
`S %*% y` raw. FIML per-row matrices do not naively satisfy the identity (same as the joint
FIML path) — not pinned.

**Smoke/edge matrix (18/18).** SG equivalence to the cbind pattern (identical scores/SEs);
block-diagonal attrs; layout pin; MG (per-group block-diagonals, `group.equal`); FIML
(per-row lists, `per_obs`, all-NA-row blocks, listwise default → clear row-count error,
`fs_indiv()` dispatch, derived pooled `tspa()` ≡ cbind-of-locals control); vector form
(verbatim fit incl. within-factor `y1 ~~ y4`); all four guards; `corrected_fsT` A/B;
`std.lv` (psi diagonal 1); `model = NULL` no-op; downstream `tspa()` derived (zero
cross-loadings in `tspaModel`), `se_fs`-only, `tspa_mx_model()` accepts the merged `fsL`.

**Tests.** `test-get_fs_local.R`: 18 blocks / 113 expectations, ~4 s (full suite ~74 s).
Full suite **3720 pass / 0 fail** (3607 baseline + 113; 1 pre-existing warn in
`test-tspa_mx.R`). `R CMD check` **0 / 0 / 0** (vignettes re-knit; the new
"Local vs joint scores" section in `multiple-factors.Rmd` verified live: local fit +
derived `tspa()` fit, zero cross-terms printed).
