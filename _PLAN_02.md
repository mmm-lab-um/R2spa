# Plan: `method` aliases, merMod score naming, and STATUS.md items 2/3/6/7

## Progress

- [x] **Step 1 — `method` aliases ("ML" / "EB")** — done 2026-08-15.
  - `normalize_fs_method()` helper added at `R/get_fs_methods.R:5`; called as
    the first line of `get_fs.lavaan()`.
  - `method` formal updated to `c("regression", "Bartlett", "ML", "EB")` in
    `get_fs.data.frame()`, `get_fs.lavaan()`, and `get_fs_lavaan()`.
  - Aliases documented in the shared `@param method` (one roxygen edit
    covers `get_fs()`, `get_fs.data.frame()`, `get_fs.lavaan()` via
    `@rdname get_fs`; `get_fs_lavaan()` picks it up via `@inherit get_fs`).
    `man/get_fs.Rd` / `man/get_fs_lavaan.Rd` regenerated via
    `devtools::document()` — not hand-edited.
  - 4 new tests (5 assertions) in `tests/testthat/test-get_fscore.R`:
    `EB` ≡ regression, `ML` ≡ Bartlett (data.frame and fitted-lavaan paths,
    incl. `get_fs_lavaan()`), invalid-method error listing all 4 choices.
  - `R/get_fscore_math.R` and the lower-level `compute_fscore()` /
    `augment_lav_predict()` machinery untouched, per plan.
  - Verified: `load_all()` OK, `document()` OK (2 Rd files), `test()`
    425 pass / 0 fail / 0 warn. `R CMD check --ignore-vignettes`: 2 WARNINGs
    / 3 NOTEs, A/B-verified (via `git stash`) to be pre-existing and
    identical on the pre-change tree.
   - **Committed in `db36781` (2026-08-15)** — by explicit decision that
     commit also includes the PLAN 01 style reformat (interleaved in the
     same files). Uncommitted PLAN 01 leftovers remain: `.lintr`
     (untracked) and `Suggests: lintr` in `DESCRIPTION`.
- [x] **Step 2 — merMod score naming (`fs_u0` vs. `u0_eb`) +
      `legacy_names`** — done 2026-08-15.
  - `get_fs_blocks.merMod()` / `get_fs.merMod()` take `legacy_names`
    (default `FALSE`); base names get the `"_eb"` suffix **only** when
    `legacy_names = TRUE`.
  - Default (`fs_u0`-style) output: scores `fs_u0`/`fs_u1`, SEs
    `fs_u0_se`/`fs_u1_se`, loadings `u0_by_fs_u0`, evs `ev_fs_u0`,
    **`ecov_fs_u1_fs_u0`** — i.e. the natural `augment_fs()` row/col
    iteration order. **Per user decision (2026-08-15)** the default is left
    in that natural (u-later-first) order rather than coerced to
    `ecov_fs_u0_fs_u1`, since no consumer pins it — avoids a reordering pass
    on the new-default path.
  - `legacy_names = TRUE`: post-processing rename in `get_fs.merMod()`
    (helpers `rename_legacy_fs_cols()` + `reorder_ecov_col()`) reproduces
    `u0_eb` / `u0_by_u0_eb` / `ev_u0_eb` / **`ecov_u0_eb_u1_eb` (ascending)**
    — exactly what `vignettes/multilevel.rmd` hardcodes, so the legacy path
    *does* reorder. This closes STATUS item 1 (merMod column-name
    regression).
  - `get_fs_lmer()` defaults `legacy_names = TRUE` (back-compat for
    `multilevel.rmd`).
  - Fixed the overwritten `fsT_j` rownames bug (single correct `fs_-`prefixed
    assignment) and the duplicated `re_names`/`fs_names` recompute (now
    derived from `colnames/rownames(blocks[[1]]$fsL)`) — STATUS item 7's
    first + fourth sub-bullets, both covered as a side effect.
   - 3 new tests in `tests/testthat/test-get_fscore.R` pin the default,
     legacy, and `get_fs_lmer()` column-name sets.
   - **Revised 2026-08-15 (review decision — compat relaxed):** the legacy
     path is **name-compatible, not byte-compatible**, with the true
     pre-refactor `get_fs_lmer()` output (d1faa7e, which had 9 columns and
     no `_se`/attributes/subject row names). Decision: keep the extra
     `_se` columns, `fsL`/`fsT` array attributes, and NULL row names
     (all current consumers work — `vignettes/multilevel.rmd` /
     `tspa_mx_model` verified end-to-end) and **document the delta** in the
     roxygen of `get_fs_lmer()` and `get_fs.merMod()` instead of dropping
     columns for byte-compat.
   - **Lints deferred (user decision 2026-08-15):** 2 new style lints in
     `rename_legacy_fs_cols()` (`R/get_fs_methods.R:392` indentation, `:397`
     commented-code) and all pre-existing lints are intentionally not
     addressed this pass.
   - Verified: `load_all()` OK, `document()` OK (2 Rd, idempotent), `test()`
     **428 pass / 0 fail / 0 warn**. Full `devtools::check()` (incl. vignette
     build) reserved for the later STATUS-closing pass.
- [x] Step 3 — vignette breakage on `format = "unified"` (STATUS item 2) —
      done 2026-08-16.
    - Root causes (verified 2026-08-16 on the post-Steps 1–2 tree; re-run
      justified because `R/` changed):
      - `corrected-se.Rmd` lines 270–288 (`group length is 0 but data length
        > 0`): `tspa_mf()` assumed list `fsT` ⇒ list `fsL`; the vignette
        legitimately passes list-valued unified `fsT` with a plain identity
        `fsL` (Bartlett). `fsL[[1]]` on a matrix extracted a column, so
        `colnames()` was empty and the split failed cryptically.
      - `tspa-vignette-mx.Rmd` lines 61–79 (`incorrect number of
        dimensions`): direct matrix indexing/arithmetic on the now
        list-valued `fsT`/`fsL` attributes.
      - `multilevel.rmd` already passed (Step 2 `legacy_names` default).
    - `R/tspa.R`: `tspa()` now validates `fsT`/`fsL` group-count consistency
      before the name check (plain matrix = 1 group, so a length-1 list may
      be mixed with a plain matrix for a single-group model), gives a clear
      `must be a list of the same length` error for true mismatches, and
      derives `multigroup` from the validated `fsT` shape; `tspa_mf()`
      accepts any matrix/length-1-list combination for single-group and
      keeps the multigroup mismatch guards.
    - `vignettes/tspa-vignette-mx.Rmd`: `format = "list"` on the two
      single-group `get_fs()` calls (direct attribute arithmetic needs the
      plain-matrix shape; `tspa_mx_model()` does not use the attributes).
    - 6 new assertions in `tests/testthat/test-tspa.R`: single-group
      list/matrix mixed shapes numerically equal the plain-matrix fit;
      multigroup shape mismatch is a clear error in both directions.
    - Verified 2026-08-16: `load_all()` OK; `document()` idempotent;
      `test()` **434 pass / 0 fail / 0 warn** (428 baseline + 6 new).
      Clean full vignette rebuild: **13/13 build, 0 failures** —
      `corrected-se.Rmd` and `tspa-vignette-mx.Rmd` pass including the
      previously untested chunks after their first failure. Scoped
      `R CMD check --ignore-vignettes --no-manual`: 2 WARNINGs, 3 NOTEs —
      identical to the pre-existing baseline documented in STATUS.md
      (S3 consistency, undocumented `fsm`/`...`, `.lintr`, unused `Matrix`
      import, missing `stats::model.frame`).
    - STATUS.md item 2 closes in the final STATUS pass (per the
      Verification section) after Step 6.
- [x] Step 4 — `get_fs(data = matrix)` dispatch (STATUS item 3) — done
      2026-08-16.
    - No code change needed: `get_fs.data.frame()` already re-dispatches
      matrices (`R/get_fs_methods.R` — `is.matrix(data)` →
      `as.data.frame(data)` → `get_fs(data, ...)`), which was already in
      place in the working tree (STATUS.md simply not updated).
    - Live check 2026-08-16: `get_fs(as.matrix(PoliticalDemocracy), ...)`
      is `identical()` to `get_fs(PoliticalDemocracy, ...)` incl.
      `fsT`/`fsL`/`fsb` attributes; unsupported classes still hit the
      clear `get_fs.data.frame()` stop message.
    - Added regression test in `tests/testthat/test-get_fscore.R`:
      `get_fs(as.matrix(...), single_model, format = "list")` is
      `identical()` to the precomputed data-frame result.
    - Verified: `test()` **435 pass / 0 fail / 0 warn** (434 after
      Step 3 + 1 new).
    - STATUS.md item 3 closes in the final STATUS pass as
      already-resolved.
- [x] Step 5 — cheap-lookup `lavInspect()` sites (STATUS item 6) — done
      2026-08-16.
    - Replaced exactly the 3 sites listed in STATUS item 6:
      - `R/tspa_corrected_se.R:21` — `lavInspect(tspa_fit, "ngroups")` →
        `tspa_fit@Data@ngroups`.
      - `R/get_fs_methods.R:227` (line number post-Step 2 shifts; plan
        said :197) — `lavInspect(object, "norig")` →
        `unlist(object@Data@norig)`.
      - `R/grandStandardizedSolution.R:85` — `lavInspect(object, "nobs")`
        → `unlist(object@Data@nobs)`.
    - Per-site behavior check (plan §5 requirement that surrounding logic
      behaves identically):
      - Slot values proven equal to `lavInspect()` for SG and MG fits:
        `@Data@ngroups` identical; `unlist(@Data@nobs)` /
        `unlist(@Data@norig)` `all.equal()` — both `nobs`- and `norig`
        -slot list form unlists to the same named vector `lavInspect()`
        returns.
      - `grand_standardized_solution()` MG path: covered by
        `test-grandStandardizedSolution.R`, which A/Bs the function's
        `est.std`/`se` against a hand calculation that still uses
        `lavInspect(fit3, "nobs")` — passes.
      - `norig` reliability line (no unit test — only SG reliability
        tested): live A/B 2026-08-16 on an MG single-factor CFA —
        pre-edit `lavInspect("norig")` path vs. post-edit slot path give
        identical per-group + overall reliability, and the actual
        `get_fs(..., reliability = TRUE)` attribute matches both.
      - `vcov_corrected()` (no unit test): runs end-to-end on an MG
        `tspa()` fit — 18×18 finite symmetric matrix.
    - Out of scope (noted, not fixed): `R/get_fscore_math.R:128`
      (`lavInspect(lavobj, what = "ngroups")` in `augment_lav_predict()`)
      is a fourth cheap lookup missing from the STATUS item 6 list —
      candidate for a future step, deliberately not expanded here.
    - Verified: `test()` **435 pass / 0 fail / 0 warn** (no new tests —
      no roxygen tag changes, so no `document()` regeneration needed).
    - STATUS.md item 6 closes in the final STATUS pass.
- [x] Step 6 — dead-code cleanup (STATUS item 7) — done 2026-08-16.
    - Removed the two remaining STATUS item 7 pieces:
      - `R/get_fscore.R` (`augment_fs()`) — deleted commented-out
        `# fs_se[is.nan(fs_se)] <- 0` line.
      - `R/get_fs_methods.R` (`get_fs.data.frame()`) — deleted the
        unreachable `if (!is.data.frame(data)) data <- as.data.frame(data)`
        guard (S3 dispatch on `data.frame` class already guarantees it).
      - The other two sub-bullets (overwritten `fsT` rownames, duplicate
        `re_names` computation) were already fixed in Step 2.
    - Verified: load_all + `document()` idempotent + `test()`
      **435 pass / 0 fail / 0 warn**.
    - Final full verification (plan §Verification item 4): `R CMD build .`
      (knits all vignettes) + `R CMD check --no-manual` → all **13/13
      vignettes rebuild** (`checking package vignettes ... OK`); check
      status **2 WARNINGs, 3 NOTEs** — the same 5 pre-existing baseline
      items, no new findings from any PLAN 02 step.
    - Post-check cleanup: removed `R2spa.Rcheck/`, rebuilt tarball
      restored (tracked file), `tests/testthat/_problems/` artifact removed.
    - STATUS.md final pass done: items 1–3, 6, 7 closed (2026-08-16);
      items 4, 5 left open; verification state refreshed.

## Scope

1. `get_fs()` (data.frame + lavaan path): accept `method = "ML"` (alias for
   `"Bartlett"`) and `method = "EB"` (alias for `"regression"`), in addition
   to the existing `"regression"`/`"Bartlett"` spellings.
2. `get_fs.lmer()`/`get_fs.merMod()`: new default `fs_u0`-style random-effect
   score naming (dropping the legacy `_eb` suffix, now that `fs_` signals
   "factor score" and the method is separately documented), with a
   `legacy_names` argument for old-style `u0_eb`-style output. **Revised
   2026-08-15 (review):** the legacy path only promises name-compatible
   output (same legacy column names, same order) — not byte-compatible
   (keeps the extra `_se` columns, `fsL`/`fsT` attributes, NULL row names;
   delta documented on `get_fs_lmer()`/`get_fs.merMod()`). This also
   resolves STATUS.md item **1** (merMod column-name regression), so item 1
   will be closed as a side effect even though it wasn't explicitly listed.
3. STATUS.md item **2** — vignette breakage under `format = "unified"`
   default.
4. STATUS.md item **3** — `get_fs(data = matrix)` dispatch.
5. STATUS.md item **6** — remaining expensive `lavInspect()` cheap-lookup
   sites.
6. STATUS.md item **7** — dead-code cleanup.
7. Explicitly **out of scope**: STATUS.md items 4 (reliability guard) and 5
   (`correct_evfs()` perf).

## 1. `method` aliases ("ML" / "EB")

Add a small helper (in `R/get_fs_methods.R`, near the other lavaan-path
helpers) that normalizes user-facing method strings before they reach the
existing `"regression"`/`"Bartlett"` machinery in `get_fscore_math.R`
(left untouched — that file is the sensitive numerics touchpoint per
AGENTS.md):

```r
normalize_fs_method <- function(method) {
  method <- match.arg(method, c("regression", "Bartlett", "ML", "EB"))
  switch(method, ML = "Bartlett", EB = "regression", method)
}
```

- Call `method <- normalize_fs_method(method)` as the first line inside
  `get_fs.lavaan()` (the convergence point for both the `data.frame` and
  `lavaan` dispatch paths — `get_fs.data.frame()` just forwards `method`
  through `get_fs(fit, method = method, ...)`, which re-dispatches to
  `get_fs.lavaan()`). Everything downstream (`correct_evfs()`,
  `get_fs_blocks.lavaan()`, `vcov_ld_evfs()`, `compute_fsrel()`) keeps calling
  `match.arg(method, c("regression", "Bartlett"))` unchanged and will receive
  an already-normalized value.
- Update the formal default / documented choices to
  `method = c("regression", "Bartlett", "ML", "EB")` in `get_fs()` (roxygen
  in `R/get_fscore.R`), `get_fs.data.frame()`, `get_fs.lavaan()`
  (`R/get_fs_methods.R`), and the legacy `get_fs_lavaan()` wrapper
  (`R/get_fscore.R`) — cosmetic/documentation only, since the actual
  validation happens once in `get_fs.lavaan()`.
- Update the `@param method` roxygen text at each of those sites to mention
  the aliases, e.g. *"`'ML'` is an alias for `'Bartlett'`; `'EB'` is an alias
  for `'regression'`."*
- No change to `compute_fscore()`, `augment_lav_predict()`,
  `compute_lav_fs_matrices()`, etc. — those remain `"regression"`/`"Bartlett"`
  only, consistent with them being lower-level exported utilities, not the
  `get_fs()` user surface this task targets.

## 2. merMod score naming (`fs_u0` vs. `u0_eb`)

Current flow: `get_fs_blocks.merMod()` builds `re_names <- paste0("u", k, "_eb")`
and uses them as the pre-`augment_fs()` column names; `augment_fs()` then
prepends `"fs_"` unconditionally, yielding `fs_u0_eb` today (not the
pre-refactor `u0_eb`, and not the desired new `fs_u0`).

New design:
- Add `legacy_names = FALSE` to `get_fs_blocks.merMod()` and
  `get_fs.merMod()`.
- Base names become `paste0("u", seq_len(num_re) - 1)` (i.e. `"u0"`, `"u1"`,
  ...). When `legacy_names = TRUE`, append `"_eb"` as before.
- Default (`legacy_names = FALSE`) result, via the existing shared
  `augment_fs()` prefixing, matches the lavaan naming convention already
  documented in AGENTS.md:
  - score: `fs_u0`
  - loading: `u0_by_fs_u0` (pattern `<indicator>_by_<name>`)
  - error variance/covariance: `ev_fs_u0`, **`ecov_fs_u1_fs_u0`** — left in
    `augment_fs()`'s natural row/col iteration order (u-later-first).
    **Revised 2026-08-15 (user decision):** originally planned to coerce to
    `ecov_fs_u0_fs_u1` (u0-first); since nothing pins the default-mode ecov
    order, keeping the natural order is simpler (no reordering pass). Only the
    LEGACY path reorders (below), because `vignettes/multilevel.rmd` hardcodes
    the ascending `ecov_u0_eb_u1_eb`.
- `legacy_names = TRUE`: reproduces the documented pre-refactor names
  (`u0_eb`, `u0_by_u0_eb`, `ecov_u0_eb_u1_eb`) via an explicit post-processing
  rename step in `get_fs.merMod()` after `augment_fs()` runs (since
  `augment_fs()`'s `"fs_"` prefixing is shared with the lavaan path and won't
  be special-cased there). This step will:
  - strip the leading `"fs_"` that `augment_fs()` added wherever it was
    applied to an `"u<k>_eb"`-based name (scores, SEs, the `<name>` side of
    `_by_`, and the `ev_`/`ecov_` suffixes), and
  - reorder each `ecov_A_B` pair by ascending numeric `u<k>` index (extracted
    via regex) so two-random-effect models render as `ecov_u0_eb_u1_eb`
    rather than the index-order-dependent `ecov_u1_eb_u0_eb` that
    `augment_fs()`'s row/col iteration currently produces. This is a
    self-contained rename table, not a change to `augment_fs()` itself.
   - **Caveat — moot as shipped:** the initial worry was that a *lexicographic*
     sort would misorder `u2` vs `u10`; the implemented `reorder_ecov_col()`
     instead compares the **extracted integers** (`as.integer(m[2])` vs
     `as.integer(m[3])`), so it is correct for any number of random effects,
     not just `u0..u9`.
- `get_fs_lmer()` (the deprecated back-compat wrapper in `R/get_fscore.R`)
  will call through with `legacy_names = TRUE` by default, so old code
  calling `get_fs_lmer()` keeps getting the exact old column names — this is
  what actually closes STATUS.md item 1. `get_fs()`/`get_fs.merMod()` called
  directly default to the new `fs_u0`-style names.
- Fix the dead/overwritten-attribute bug along the way: remove the
  overwritten `rownames(fsT_j) <- colnames(fsT_j) <- paste0("fs_", re_names)`
  line at `get_fs_methods.R:247` and keep only the single correct assignment
  using the (new) `fs_names <- paste0("fs_", re_names)` — this is STATUS item
  7's first sub-bullet.
- Fix the duplicate `re_names` computation (STATUS item 7's fourth
  sub-bullet): in `get_fs.merMod()`, derive `re_names`/`fs_names` for the
  `fsL_arr`/`fsT_arr` attribute construction from
  `colnames(blocks[[1]]$fsL)` / `rownames(blocks[[1]]$fsL)` instead of
  recomputing from `object@cnms`.
- Update roxygen for `get_fs.merMod()` and `get_fs_lmer()` to document the
  new default naming and the `legacy_names` argument. **Revised 2026-08-15
  (review decision):** also document the legacy-path delta vs. the true
  pre-refactor `get_fs_lmer()` output (d1faa7e): legacy mode additionally
  returns `_se` score-error columns, per-cluster `fsL`/`fsT` array
  attributes, and NULL row names (true legacy: no `_se` columns, no
  attributes, ranef subject-ID row names). Byte-compat explicitly relaxed
  to name-compat; no columns dropped.

No existing test pins the exact merMod column names (confirmed by search —
tests use `ignore_attr = TRUE` and slice `[, 1:2]`), so this is safe to
change without breaking `devtools::test()`. Will add a couple of small new
tests: one confirming default `fs_u0`-style names, one confirming
`legacy_names = TRUE` reproduces `u0_eb`-style names and that
`get_fs_lmer()` still returns those by default.

## 3. STATUS.md item 2 — vignette breakage on `format = "unified"`

> **Verified 2026-08-15 (pre-Step 1 baseline — do not re-verify):**
> `R CMD build` on the pre-Step 1 working tree fails on exactly **3 of 13**
> vignettes — `corrected-se.Rmd`, `multilevel.rmd`,
> `tspa-vignette-mx.Rmd` — not the 8/13 failure count claimed in STATUS.md
> item 2. Treat these three as the ground-truth failure set for this task;
> only re-run the build if vignettes or `R/` code change. Practical note:
> `devtools::check(args = "--ignore-vignettes")` still knits vignettes in the
> *build* phase (and fails there); to scope a check, use
> `R CMD build --no-build-vignettes .` then
> `R CMD check --ignore-vignettes --no-manual <tarball>`.

A static read-through of every vignette calling `get_fs(..., group = ...)`
(`R2spa.Rmd`, `corrected-se.Rmd`, `reliability.Rmd`) did **not** turn up code
that actually breaks under the unified default — `R2spa.Rmd` uses `cbind()`
and `$`-column access (fine with unified output), `corrected-se.Rmd`
immediately overwrites the `get_fs()` result with `get_fs_lavaan()`, and
`reliability.Rmd` only reads the `reliability` attribute. This conflicts with
STATUS.md's claim of "7 vignettes break" / `check()` failing at 8/13
vignette builds — that failure count plausibly comes substantially from item
1 (merMod naming breaking `multilevel.rmd`'s `tspa_mx_model` usage), which
task 2 above fixes.

Since this can only be resolved authoritatively by actually building the
vignettes, the implementation step here is:
1. Run `devtools::document()` then attempt to build/knit each vignette
   (e.g. via `tools::buildVignettes(dir = ".")` or `devtools::check()`,
   accepting it will be slow) to get the current, real list of failures and
   error messages — not trusting the STATUS.md count as ground truth per
   AGENTS.md's "unmaintained, verify before trusting" guidance.
2. For any vignette that fails **specifically because of list-vs-unified
   handling** of a multi-group `get_fs()` result, fix it by either adding
   `format = "list"` to the `get_fs()` call (minimal-diff fix, preserves the
   vignette's narrative) or, if the narrative is explicitly about the new
   unified shape, updating the downstream code to use `$group`
   filtering / `fs_to_group_list()` instead.
3. Re-run the build to confirm each fixed vignette now succeeds, and note
   any remaining failures unrelated to items 1/2 (e.g. pre-existing/unrelated
   breakage) rather than silently patching around them.

## 4. STATUS.md item 3 — `get_fs(data = matrix)` dispatch

Reading `get_fs.default()` (`R/get_fs_methods.R`) shows it already has:
```r
get_fs.default <- function(data, ...) {
  if (is.matrix(data)) {
    data <- as.data.frame(data)
    return(get_fs(data, ...))
  }
  stop(...)
}
```
which re-dispatches a converted matrix into `get_fs.data.frame()` — this
looks like it already fixes the issue described in STATUS.md (possibly
already fixed in a prior uncommitted change and STATUS.md just wasn't
updated). Implementation step: add a small regression test
(`get_fs(as.matrix(...))` produces the same result as
`get_fs(as.data.frame(...))`) to lock in the behavior; if it passes as
expected, close item 3 in STATUS.md as already-resolved rather than making
further code changes. If it turns out **not** to work as expected once
actually run, fix `get_fs.default()`/add a `get_fs.matrix()` method as
needed.

## 5. STATUS.md item 6 — cheap-lookup `lavInspect()` sites

Replace with direct slot access (`unlist()` where the slot is a per-group
list), matching the pattern already used elsewhere in the same files:

| File | Current | Replacement |
|---|---|---|
| `R/tspa_corrected_se.R:21` | `lavInspect(tspa_fit, what = "ngroups")` | `tspa_fit@Data@ngroups` |
| `R/get_fs_methods.R:197` | `lavInspect(object, what = "norig")` | `unlist(object@Data@norig)` |
| `R/grandStandardizedSolution.R:85` | `lavInspect(object, what = "nobs")` | `unlist(object@Data@nobs)` |

Will double check each call site's surrounding logic (e.g.
`grandStandardizedSolution.R`'s `if (length(ns) == 1) ns <- NULL` guard)
still behaves identically with the slot-based value before/after, using
existing tests (`test-grandStandardizedSolution.R`) plus the corrected-se and
tspa tests as the safety net.

## 6. STATUS.md item 7 — remaining dead-code cleanup

- `R/get_fscore.R:247` — delete the commented-out
  `# fs_se[is.nan(fs_se)] <- 0` line inside `augment_fs()`.
- `R/get_fs_methods.R:14` (`get_fs.data.frame()`) — remove the unreachable
  `if (!is.data.frame(data)) data <- as.data.frame(data)` guard (S3 dispatch
  on `data.frame` already guarantees this).
- The other two sub-bullets (rownames double-assignment, duplicate
  `re_names` computation) are handled as part of task 2 above.

## Verification (per AGENTS.md build/test lifecycle — exact order)

1. `devtools::load_all()`
2. `devtools::document()` (roxygen changes: `method` param docs,
   `legacy_names` param docs)
3. `devtools::test()` — expect the existing 420 tests to still pass, plus
   new tests for: `method = "ML"/"EB"` aliasing, merMod default vs.
   `legacy_names` naming, and the matrix-dispatch regression test.
4. `devtools::check()` — including the vignette build investigation from
   task 3; report final pass/fail state.
5. Update `STATUS.md`:
   - Move items 1 (via task 2 side effect), 2, 3, 6, 7 to **Closed** with
     today's date, and a short note per item on what changed.
   - Update the **Verification state** section with fresh
     `test()`/`document()`/`check()` results.
   - Leave items 4 and 5 in **Open** (explicitly out of scope for this pass).

## Files touched

- `R/get_fscore.R` — `method` doc/default updates, `get_fs_lavaan()`
  doc/default, `get_fs_lmer()` (`legacy_names = TRUE` passthrough),
  `augment_fs()` dead-comment removal.
- `R/get_fs_methods.R` — `normalize_fs_method()` helper, `get_fs.data.frame()`
  dead-code removal + method doc, `get_fs.lavaan()` normalization call +
  method doc, `get_fs_blocks.merMod()` / `get_fs.merMod()` renaming rework +
  `legacy_names`, `get_fs_methods.R:197` lavInspect fix.
- `R/tspa_corrected_se.R` — `lavInspect("ngroups")` fix.
- `R/grandStandardizedSolution.R` — `lavInspect("nobs")` fix.
- `tests/testthat/test-get_fscore.R` (or a new small test file) — new
  regression tests noted above.
- Vignette file(s) as determined by the task-3 build investigation
  (unknown count until vignettes are actually built).
- `STATUS.md` — close out resolved items, refresh verification state.
