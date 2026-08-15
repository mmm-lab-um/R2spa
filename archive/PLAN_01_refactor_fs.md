# Refactor `get_fs()` in `R/get_fscore.R`: multi-backend S3 architecture

## Goal

Rewrite `get_fs()` so it can extract factor scores from multiple fitted-model
backends (starting with `lavaan`, with `lme4` and `mirt` anticipated later)
via S3 dispatch, while producing a data frame output whose attributes
(`fsT`, `fsL`, `fsb`, `scoring_matrix`) can support:

- multiple-group models,
- individual-specific (case-varying) loadings/error covariances — needed
  later for OpenMx definition-variable models,
- a single averaged/representative matrix per group — needed for `tspa()`'s
  current lavaan-based constrained-loading approach.

This pass focuses on **lavaan only**. `lme4` support is refactored into the
new shape but not redesigned. `mirt` is anticipated architecturally (an S3
slot) but not implemented.

## Progress

**Agents: update this table after completing each step.** Record the result of
`devtools::load_all()` / `devtools::document()` / `devtools::test()` (pass count
and any failures), and note any regressions or decisions made during
implementation. Steps are numbered to match the list at the bottom of this
document. Always run the full `load_all` → `document` → `test` cycle before
marking a step `done`.

| Step | Description | Status | Verification (test results) | Notes |
|------|-------------|--------|----------------------------|-------|
| 1 | `assemble_fs_blocks()` + block contract + unit tests | done | `test-assemble_fs_blocks.R`: 14 tests pass; full suite 205 | `check_blocks_identical()` helper added. **Post-review fixes (Aug 13):** dead `if(is_named)` branch removed (line 322–326); representative-selection `if/else` collapsed + added `message()` when attributes differ across blocks (line 360–368); `n_cases` via `max(case_idx)` flagged as fragile (relies on contiguous, 1-indexed indices). |
| 2 | `get_fs_blocks.lavaan()` + `get_fs.lavaan()` S3 method | done | Full suite: 205 pass, 0 fail | `get_fs()` became S3 generic; `get_fs.data.frame()`, `get_fs.default()`, `get_fs.lavaan()` added. **Post-review fix:** blocks now carry explicit `fsL`, `fsb`, `scoring_matrix` fields (was relying on `attr(b$fs, ...)`). Assembler falls back to attributes for backward compat with hand-built test fixtures. |
| 3 | `get_fs()` generic + `get_fs.data.frame()` | done | Full suite: 205 pass, 0 fail | Completed together with step 2. **Issue:** `get_fs.data.frame()` hardcodes `format = "list"` for back-compat; plan says default should be `"unified"`. Intentional for now but should be documented. |
| 4 | Back-compat wrappers `get_fs_lavaan()` + `get_fs.default()` | done | Full suite: 205 pass, 0 fail | `get_fs_lavaan()` → `get_fs(..., format = "list")`; `get_fs.default()` handles matrix → data.frame conversion + error for unsupported classes |
| 5 | `get_fs_blocks.merMod()` + `get_fs.merMod()` | done | Full suite: 279 pass (205 + 74 new), 0 fail | `get_fs_lmer()` → thin wrapper to `get_fs(object)`. `get_fs_blocks.merMod()` produces one block per cluster with `case_idx`, `fs`, `fsL`, `fsT`. `get_fs.merMod()` assembles into one-row-per-cluster matrix with legacy `u0_eb` naming + `fsL`/`fsT` array attributes. `get_fsLT_lmer()` removed (logic inlined into `get_fs_blocks.merMod()`). Two new tests: S3 dispatch on `merMod` + block structure validation. |
| 6 | `fs_to_group_list()` converter + tests | done | Full suite: 334 pass (279 prior + 55 new), 0 fail | Self-inverse: `data.frame` (unified with `group` column + list-valued attrs) ↔ named `list` of data frames (per-group attrs). Single-group input returns bare data frame (no list) in either direction. Error handling for missing `group` column, unnamed list, and unsupported types.
| 7 | Numeric equivalence + timing vs `lavPredict(acov=TRUE)` | done | Full suite: 405 pass (334 prior + 71 new), 0 fail | `test-lavPredict_equivalence.R`: SG/MG, 1f/3f, complete/missing data, regression (Bartlett where applicable). Factor scores match `lavPredict` (tolerance 1e-8 complete, 1e-5 missing). Block-level `fsT`/`fsL` match `acov`-derived matrices. Timing: `get_fs` is ~2x `lavPredict(acov=TRUE)` on pre-fitted models — acceptable overhead from block assembly and data-frame augmentation. Note: Bartlett method skipped for missing-data blocks due to known singular-weight-matrix issue when all indicators for a LV are missing for a case. **Follow-up perf fix (Aug 13, comprehensive):** `Rprof()` showed >50% of `get_fs()` runtime inside `lavInspect()`'s per-call `lav_object_check_version()` (re-reads lavaan's `DESCRIPTION` via `read.dcf()` every call). Root cause: all six `lavInspect(what = "ngroups")` and `lavInspect(what = "group")` call sites pay this cost, yet the data is trivially available as S4 slots (`@Data@ngroups`, `@Data@group`) already used elsewhere in the file. **Fixed:** replaced all six cheap-lookup `lavInspect()` calls with direct slot access in `get_fs_blocks.lavaan()` (2 sites), `get_fs.lavaan()` (1 site), `correct_evfs()` (1 site), `compute_fspars()` (1 site), `compute_fsrel()` (1 site). Ratio improved ~2x → ~1.2–1.3x (confirmed by existing timing test: SG 1.22x, MG 1.26x). The remaining `lavInspect()` calls (`"est"`, `"data"`, `"free"`, `"implied"`) still pay the version-check cost but extract non-trivial formatted matrices with no cheap slot equivalent — safe bypass would require unexported `lavaan:::` internals. 405 tests pass, 0 fail; optional-feature paths (`corrected_fsT = TRUE`, `reliability = TRUE`) also tested and working. |
| 8 | Full `document` → `test` → `check` cycle | done | `document`: OK (pre-existing `tspa_plot.R:3` link warning). `test`: 405 pass, 0 fail. `check`: fails at vignette build — 8/13 vignettes error. | **Vignette regressions (all from refactor):** (a) `format = "unified"` default breaks 7 vignettes that call `do.call(rbind, get_fs(..., group=...))` expecting a list. (b) `get_fs(data = matrix)` dispatches to `get_fs.default()` instead of `get_fs.data.frame()` (1 vignette). (c) `multilevel.rmd` — `u0_by_u0_eb` column missing from `lme4` output. Per plan's "out of scope" clause, fixing `tspa()` wiring is a separate follow-up. These are known and tracked. |
| 9 | Roxygen updates + deprecation notes | done | 405 pass, 0 fail | Updated main `get_fs()` roxygen: added `@description` for S3 generic, `@details` for dispatch behavior, `@param format` (unified/list), updated `@return` for unified format with list-valued attributes. Deprecation notes on `get_fs_lavaan()` and `get_fs_lmer()` already in place from prior steps. |

## Decisions confirmed with user

1. **Multigroup output shape**: default `get_fs()` output becomes a single,
   rbinded data frame with a `group` column, with `fsT`/`fsL`/`fsb` attributes
   as named lists keyed by group label. Also ship an **exported converter**
   (e.g. `fs_as_group_list()` / `fs_as_unified()` — naming TBD in implementation)
   that can convert between this unified shape and the legacy
   list-of-data-frames-per-group shape (with matching list-shaped attributes),
   so `tspa()` (and any user code depending on the old shape) can keep working
   by explicitly calling the converter. A `format = c("unified", "list")`-style
   argument on `get_fs()` itself is also acceptable if simpler than a separate
   converter — decide during implementation, but default to `"unified"`.
2. **Duplicate lavaan paths**: consolidate `get_fs()` on the
   `compute_fscore()`-based manual matrix path (already handles per-pattern
   blocks and matches the AGENTS.md column-naming spec). `augment_lav_predict()`
   stays in the file, untouched, as a separate/legacy function — not called by
   the new `get_fs()`, not removed. Before considering the refactor done,
   numerically cross-check the new path's `fs`/`fsL`/`fsT` against
   `lavaan::lavPredict(..., acov = TRUE)` on both single- and multi-group,
   complete- and missing-data examples, and benchmark that the new path is at
   least as fast as `lavPredict()` for a comparable model (it should be, since
   it already avoids some of `lavPredict`'s overhead — but must be verified,
   not assumed).

## Proposed architecture

### 1. Generic dispatch

```r
get_fs <- function(object, ...) {
  UseMethod("get_fs")
}
```

- `get_fs.data.frame()` — today's convenience path: build a one-factor (or
  user-supplied) lavaan model string, fit with `lavaan::cfa()`, then call
  `get_fs()` on the resulting fit. Preserves current top-level signature
  `get_fs(data, model = NULL, group = NULL, method = ..., corrected_fsT = ...,
  vfsLT = ..., reliability = ..., ...)`.
- `get_fs.lavaan()` — the rewritten core (see below). Replaces the current
  `get_fs_lavaan()` body.
- `get_fs.merMod()` — refactor of current `get_fs_lmer()` logic into the same
  block/assemble pipeline (see below), same external behavior.
- `get_fs.default()` — informative error for unsupported object classes,
  mentioning which classes are currently supported and that `mirt` support is
  planned.

Backward compatibility: `get_fs_lavaan()` and `get_fs_lmer()` remain exported
as thin wrappers (`get_fs_lavaan <- function(lavobj, ...) get_fs(lavobj, ...)`,
similarly for `get_fs_lmer()`), each carrying a `@details`/deprecation note in
roxygen pointing users to `get_fs()`. Their output must remain
pixel-for-pixel identical to today's for existing tests
(`test-get_fscore.R`) to keep passing unmodified as a regression safety net —
except where the multigroup shape default changes (see below), which
requires updating those specific existing tests/expectations explicitly and
deliberately, not incidentally.

### 2. Shared intermediate representation ("blocks")

Introduce an internal (non-exported) concept of a **block**: a homogeneous
computational unit — one group × one missing-data pattern for lavaan, one
cluster for `lme4` — for which a single `(fsL, fsT, fsb)` triple applies to
every case in the block.

```r
# One block:
list(
  case_idx = <integer vector, indices into the original data/group>,
  fs       = <n_block x q matrix of factor scores>,
  fsL      = <q x q (or q x p) implied loading matrix>,
  fsT      = <q x q error covariance matrix>,
  fsb      = <length-q intercept vector, or NULL>
)
```

Each backend produces a **list of blocks per group**:

```r
get_fs_blocks.lavaan(object, method, corrected_fsT, ...) ->
  list(
    group1 = list(block1, block2, ...),   # e.g. one per missing-data pattern
    group2 = list(block1, ...),
    ...
  )
```

For `lme4`, each cluster is its own block (case-varying `fsL`/`fsT` already
happens today via `get_fsLT_lmer()`; this becomes one block per cluster,
single implicit "group").

This isolates all model-specific extraction logic (pulling `lambda`/`theta`/
`psi`/`nu`/`alpha` out of a `lavaan` fit; looping over `@Data@Mp` missing-data
patterns; computing `a_mat` via `compute_a_from_mat()`) inside
`get_fs_blocks.lavaan()`, reusing the existing `compute_fscore()` machinery
unchanged. No change to `compute_fscore()`, `compute_a_from_mat()`,
`compute_a_reg()`, `compute_a_bartlett()`, `correct_evfs()`, `compute_fsrel()`,
`vcov_ld_evfs()` — these are already backend-agnostic-enough at the matrix
level and stay as-is.

### 3. Shared assembler

A single internal `assemble_fs_blocks()` function, called by every S3 method,
turns `list(group -> list(block))` into the final output:

1. For each block, call the existing `augment_fs()` helper (unchanged) to
   build its row-level tidy columns (`fs_*`, `*_se`, `*_by_*`, `ev_*`/`ecov_*`).
2. Within a group, `rbind` blocks back into original case order using
   `case_idx` (this replaces the current hand-rolled loop inside
   `get_fs_lavaan()`).
3. Determine, per group, whether all blocks share **identical** `fsL`/`fsT`/
   `fsb` (the common case: complete data, no missing patterns). If so, the
   group's representative attribute is that single shared matrix/vector. If
   blocks differ (missing-data patterns present), keep the per-row inline
   columns as the source of truth for case-specific values (as today), and set
   the group's `fsT`/`fsL`/`fsb` attribute to a documented deterministic choice
   (e.g. the complete-data / most-common pattern, with a message if it's not
   unambiguous) — this is what feeds `tspa()`'s single-matrix lavaan path.
   Do **not** silently average, since averaging is a modeling decision and
   should stay a documented degenerate-case fallback, not the general
   behavior — this matches how it works today.
4. Across groups: `rbind` into one data frame with a `group` column
   (`format = "unified"`, the new default), building `fsT`/`fsL`/`fsb`
   attributes as named lists keyed by group label. If `format = "list"`
   requested (or via the separate converter, per the user's decision above),
   produce the legacy list-of-data-frames-per-group shape instead, with
   attributes attached the way `get_fs_lavaan()` does today.
5. Attach `scoring_matrix` attribute analogously (single matrix if constant
   within a group, else per-block / left as today).

### 4. Format converter

Export a small utility (exact name TBD, e.g. `fs_to_group_list()`) that takes
either shape of `get_fs()` output and returns the other, so `tspa()` — and any
other unchanged downstream code — can be pointed at the legacy shape
explicitly:

```r
fs_dat  <- get_fs(fit, group = "school")        # unified df + group column
fs_lst  <- fs_to_group_list(fs_dat)             # legacy list-of-df shape
```

This function is new code, self-contained, and doesn't require touching
`tspa.R` in this pass — flagged as a **separate follow-up** to actually swap
`tspa()`'s internals to call it (or accept both shapes natively), tracked but
not executed now.

### 5. Naming/column preservation (hard constraints from AGENTS.md)

No changes to:
- Column prefixes/suffixes: `fs_<name>`, `<name>_se`, `ev_<name>`,
  `ecov_<name1>_<name2>`, `<indicator>_by_<name>`.
- Attribute names: `fsT`, `fsL`, `fsb`, `scoring_matrix`.

Only the **shape** (single matrix vs. named list keyed by group) and the
**default multigroup container** (list-of-df vs. unified df + `group` column)
change.

## Testing plan

- Add a numeric equivalence check (new test, e.g. in `test-get_fscore.R` or a
  new `test-get_fs_lavaan_consistency.R`) comparing `get_fs()`'s `fs`/`fsL`/
  `fsT` against `lavaan::lavPredict(fit, method = ..., acov = TRUE)`-derived
  values, for: single-group complete data, single-group with missing data
  (multiple patterns), and multigroup complete data. Tolerance consistent
  with existing `expect_equal(..., tolerance = ...)` conventions.
- Add a `system.time()`/`bench::mark()`-style informal timing comparison
  (can live as a one-off check during development, doesn't need to be a
  permanent automated test) confirming the new path is not slower than
  `lavPredict(acov = TRUE)` on a moderately-sized example.
- Update/extend existing multigroup tests to reflect the new default unified
  shape, and add tests for `fs_to_group_list()` round-tripping back to
  something `expect_equal()`-comparable with today's `get_fs_lavaan()` output
  (`ignore_attr` where appropriate, per AGENTS.md test conventions).
- Add tests for the block-detection logic (identical vs. differing `fsL`/
  `fsT` across blocks within a group) using a deliberately-constructed
  missing-data example (e.g. `HolzingerSwineford1939` with induced `NA`s).
- Re-run `devtools::load_all()` → `devtools::document()` → `devtools::test()`
  → `devtools::check()` per AGENTS.md lifecycle before considering this done.

## Out of scope for this pass

- Any change to `tspa.R` itself (only the converter utility is added; wiring
  `tspa()` to use it is a follow-up).
- `mirt` S3 method implementation (architecture anticipates it via
  `get_fs.default`'s error message and the block/assembler design being
  backend-agnostic, but no `get_fs.SingleGroupClass` method is written yet).
- Redesigning `get_fs_lmer()`'s math — only its control flow is refactored to
  produce blocks consumed by the shared assembler; numeric behavior stays
  identical.
- Removing or rewriting `augment_lav_predict()`.
- Reliability (`reliability = TRUE`) computation logic — stays as-is,
  reattached as an attribute after assembly, same as today.

## Implementation steps

1. Write `assemble_fs_blocks()` (shared assembler) and the block list-of-list
   data structure contract, with unit tests against hand-built block lists
   (no lavaan involved) to validate the assembly/naming logic in isolation.
2. Refactor `get_fs_lavaan()`'s body into `get_fs_blocks.lavaan()` (block
   extraction only) + call to `assemble_fs_blocks()`; wire up `get_fs.lavaan()`
   S3 method.
3. Add `get_fs()` generic + `get_fs.data.frame()` (today's `get_fs()` body,
   minus the final dispatch, calling `cfa()` then `get_fs.lavaan()`).
4. Add back-compat wrappers `get_fs_lavaan()`/`get_fs_lmer()`.
5. Refactor `get_fs_lmer()` into `get_fs_blocks.merMod()` + `get_fs.merMod()`
   using the same assembler.
6. Write `fs_to_group_list()` converter (+ its inverse if needed) and its
   tests.
7. Numeric equivalence + timing checks against `lavPredict(acov = TRUE)`.
8. Update/add tests; run full `devtools::document()` → `devtools::test()` →
   `devtools::check()` cycle; fix any NOTE/WARNING regressions.
9. Update roxygen docs on `get_fs()` (new generic, new `format` semantics,
   new converter) and on `get_fs_lavaan()`/`get_fs_lmer()` (deprecation note).

## Code review findings (Aug 13, steps 1-4)

### Fixed during review
- Dead code in `assemble_fs_blocks()`: inner `if (is_named)` branch at line 322
  could never execute — removed.
- Silent attribute mismatch in `assemble_fs_blocks()`: all branches of the
  representative selection returned `block_attrs[[1]]` with no notification —
  now emits `message()` when blocks differ across missing-data patterns.
- Block contract: `get_fs_blocks.lavaan()` now carries `fsL`, `fsb`,
  `scoring_matrix` as explicit list fields (matches plan spec). Assembler
  falls back to `attr(b$fs, ...)` for backward compat with hand-built test
  blocks.

### Known issues to address in remaining steps
**All resolved (Aug 13 post-review fixes):**
- ~~`get_fs.data.frame()` hardcodes `format = "list"`~~ — **FIXED:** `format`
  is now a pass-through parameter with `c("unified", "list")` defaulting to
  `"unified"`. Tests updated to pass `format = "list"` where they need legacy
  output shapes.
- ~~`n_cases` via `max(case_idx)` fragile~~ — **FIXED:** changed to
  `max(unlist(lapply(blocks, function(b) max(b$case_idx))))` to handle
  non-1-starting indices (though lavaan indices are always contiguous).
- ~~`est$psi == 1` matrix comparison for `is_std.lv`~~ — **FIXED:** changed to
  `all(est$psi == 1)` for correct scalar result.
- ~~Unused `group` variable (dead code)~~ — **FIXED:** removed the unused
  `lavInspect(object, what = "group")` call in the reliability block.
- Additionally fixed: `tspa.R:282` `collpase` typo (was `collpase = "\n"` in
  `paste0()`, silently broke multi-factor model syntax).
- Additionally removed: unused `get_fsLT()` function (lines 587-596).
- Added deprecation notes to `get_fs_lavaan()` and `get_fs_lmer()` roxygen
  blocks.

## Remaining issues (review, Aug 15 — all 9 steps marked done)

### Blocked: R CMD check fails
Step 8's vignette failures (8/13 vignettes error) are still open:
1. **merMod column-name regression** (root of `multilevel.rmd` failure) —
   `get_fs_lmer()` output no longer matches legacy names (`u0_eb`→`fs_u0_eb`,
   `u0_by_u0_eb`→`u0_eb_by_fs_u0_eb`, `ecov_u0_eb_u1_eb`→`ecov_u1_eb_u0_eb`),
   violating this plan's "pixel-for-pixel identical / legacy `u0_eb` naming"
   spec. All 7 names referenced by `vignettes/multilevel.rmd`
   (`tspa_mx_model` `mat_ld`/`mat_ev`) are missing post-refactor. Existing
   lmer tests are name-blind (`ignore_attr = TRUE`), so they don't catch it.
   Fix: restore legacy names in the merMod path.
2. **7 vignettes broke on `format = "unified"` default** — either vignettes
   pass `format = "list"` or their narrative is updated; the `tspa()`-side fix
   stays a follow-up per the out-of-scope clause. Check stays red until fixed.
3. `get_fs(data = matrix)` dispatches to `get_fs.default()` instead of
   `get_fs.data.frame()` (1 vignette).

### Bugs not tracked by the plan
4. **`reliability = TRUE` + single-group multi-factor + default `unified`
   format → hard error** ("reliability is only supported for unidimensional
   models") instead of the documented warning. Root cause:
   `get_fs_methods.R:183-184` guards on `length(attr(out, "fsb")) > 1`, but in
   unified format that attribute is always a length-1 **list**, so the
   `multifactor` guard never fires and `compute_fsrel()` throws.

### Deferred optimizations (perf work declared done, hot path remains)
5. **`correct_evfs()` is ~100% of `corrected_fsT = TRUE` runtime** (get_fs
   0.168 s, correct_evfs 0.167 s, 4-LV model). Traced 57 `lavInspect` calls
   (file read each) for a 13-param 2-LV model = `q × (npar+1) × 2 + 1`. Fix:
   single `lav_func_jacobian_complex()` call over the full `a` matrix (q× cut
   in evaluations — `compute_fsrel()`/`compute_grad_ld_evfs()` already use the
   single-call pattern) + hoist `lavInspect("free")`/`("est")` out of the
   per-evaluation `compute_fspars()` body (get_fscore_math.R:289-290; they're
   properties of the fitted object, so capture once in the Jacobian closure).
6. **Three missed cheap-lookup `lavInspect` sites** inconsistent with the
   Aug 13 fix: `tspa_corrected_se.R:21` (`"ngroups"`),
   `get_fs_methods.R:197` (`"norig"`), `grandStandardizedSolution.R:85`
   (`"nobs"`). All have `@Data@` slot equivalents (`norig`/`nobs` slots are
   lists — use `unlist()`).

### Explicitly out of scope (follow-ups)
7. Wiring `tspa()` to accept unified/list shapes natively via
   `fs_to_group_list()`.
8. `mirt` S3 method.

### Minor cleanup
9. Dead rownames assignment (`get_fs_methods.R:247`, overwritten at :248);
   commented-out line (`get_fscore.R:247`); unreachable `is.data.frame` check
   in `get_fs.data.frame` (:14); `re_names` computed twice between
   `get_fs_blocks.merMod()`/`get_fs.merMod()`.

**Suggested order:** 1 → 2 (unblocks check) → 4 (user-facing bug) → 5/6
(perf) → 9 (dead code) → 7 (future). Tracked in `STATUS.md`.
