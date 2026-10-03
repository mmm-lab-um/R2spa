# PLAN 18 — Usability first batch: P1 correctness + doc/vignette fixes + onboarding

Status: **implemented** (2026-10-03)
Source: usability review of the package (API ease-of-use + vignette coverage). This is the
first implementation batch. The full S3 score-result class refactor is **deferred** (see end).

## Motivation

The review found:

1. `tspa()` silently replaces explicit measurement inputs with auto-derived values,
   violating the documented "explicit arguments always win" contract.
2. `get_fs()` drops trailing fully-missing cases, and row reordering/subsetting detaches
   the row-specific attributes (which live in data-frame attributes, not columns).
3. Several vignettes and help pages contradict the current implementation (mean-score
   `vfsLT`, `std.lv` "bias", formulas, `vignette()` lookups, loading direction, EAP SE).
4. The onboarding path (README, intro vignette, pkgdown) teaches a manual workflow the
   current API no longer needs, and has no curated navigation.

Scope decisions (2026-10-03):

- Prioritized first batch: P1 API correctness + all documentation/vignette corrections +
  onboarding improvements.
- Row safety: **detect + document** (runtime guard + prominent docs). The S3
  score-result class is a follow-up plan.
- Explicit-input precedence: **fix the code** to honor explicit arguments; update the
  tests that assert the old replacement behavior.
- Multilevel vignette: **fix the `e2` bug and re-run** the simulation, updating the
  numerical narrative.

## Ground rules (from AGENTS.md — no exceptions)

1. After any `R/*.R` edit: `devtools::load_all()`.
2. After any roxygen tag change: `devtools::document()` (regenerates NAMESPACE + man/).
3. After every functional change: `devtools::test()`.
4. Before marking the plan complete: `devtools::check()` — expected 0 errors / 0 warnings
   / 0 NOTEs (`--no-manual` on this LaTeX-less machine).
5. **Never hand-edit `NAMESPACE` or `man/*.Rd`.** Edit roxygen in `R/`, rerun `document()`.
6. **Never add `library()`/`require()` in a function body.** Namespace calls
   (`lavaan::sem()`) or `@importFrom` + `document()`.
7. Don't touch `.quarantine/`, `legacy/`, or `data-raw/`.
8. Don't regenerate vignette RDS fixtures **except** the one marked in Phase 5 item 7
   (`sim_scoring_se.RDS`).
9. No new `Imports`/`Suggests`; no tidyverse/dplyr/purrr; no `%>%` in `R/`.
10. End each phase with the full test suite green before starting the next.
11. Add a `NEWS.md` entry at the end (Phase 5 item 18).

Verification shorthand used below: "probed" means the behavior was confirmed in an
in-memory R session against the 2026-10-03 `devel` checkout; re-probe where noted, since
line numbers may drift.

---

## Phase 1 — `tspa()` explicit-input precedence (P1)

**Problem.** `man/tspa.Rd` promises explicit `se_fs`/`fsT`/`fsL`/`fsb` win over
derivation, but three paths violate it:

1. **`fsb`** — `R/tspa.R:513–536`: the derivation block runs when `fsT`, `fsL`, `se_fs`
   are all omitted and unconditionally does `fsb <- attr(data, "fsb")`, overwriting a
   supplied `fsb`. Probed: two-factor fit with `fsb = c(fs_a = 2, fs_b = 3)` produced
   resolved `tspa_args$fsb` of zeros from the data attribute.
2. **Per-unit `fsT`/`fsL`** — `R/tspa.R:595–603`: when the supplied inputs have a
   poolable per-unit shape, `pool_per_unit(data, reduce, ...)` re-reads the data's
   attributes and discards the supplied matrices. Probed: doubling every supplied
   per-pattern `fsT` left the fitted/pooled ratio at 1.
3. **`se_fs`** — `R/tspa.R:697–710`: when the per-row `fs_<v>_se` columns vary within a
   group, an explicit `se_fs` is replaced by `pool_se_fs()` values. This is asserted by
   `tests/testthat/test-tspa_pooled.R:206–235` — i.e. the bug is codified in a test.

**Changes**

- Track supplied-ness before derivation: `fsb_given <- !missing(fsb)`; use `missing()`
  for `fsT`/`fsL` too where the NULL-default convention allows. Keep the existing
  NULL-vs-`list()` convention for `se_fs` (an explicit empty `list()` suppresses the
  multi-factor derivation — see `R/tspa.R:494–495` and the roxygen at `R/tspa.R:66–67`).
- Derive `fsb` at `R/tspa.R:536` only when `!fsb_given`.
- Thread explicit values into pooling: give `pool_per_unit()` (`R/tspa.R:1032`) optional
  `fsT`/`fsL`/`fsb` override arguments. After `resolve_fs_per_row()`, overwrite each
  resolved block's `fsT`/`fsL`/`fsb` from the matching override element, then let the
  existing per-row accumulation pool the overridden values:
  - lavaan per-pattern: override is a per-group list of one matrix/vector per
    pattern, in the same pattern-label order as the resolved blocks;
  - `merMod`: override is a 3-D array, slice `j` per cluster in resolved block order;
  - mirt per-obs: override is a flat list, one matrix/vector per row.
- Restrict the within-group `se_fs` pooling (`R/tspa.R:697–710`) to the **derived** case
  (`se_fs` omitted). An explicit `se_fs` is used as-is (this is the existing
  complete-data behavior; FIML simply no longer overrides it).
- Roxygen: after the fix the "explicit arguments always win" claim
  (`R/tspa.R:29–39`, `man/tspa.Rd`) is true; re-read the single-factor loading note
  (`R/tspa.R:49–64`) to confirm it still matches (derived `se_fs` recovers
  shrinkage loadings; explicit `se_fs` keeps unit loadings).

**Tests**

- Rewrite `tests/testthat/test-tspa_pooled.R:206–235`: explicit
  `se_fs = c(visual = 0.35)` must now produce the model string rendered from 0.35
  (the test currently asserts the pooled replacement — invert the expectation and keep
  the `tspa_render(tspa_schema_sf(...))` byte-comparison technique with the explicit
  value).
- Add a sibling test: `se_fs` **omitted** on the same FIML multigroup input → pooled
  per-group means (the derived path keeps its current behavior).
- New test (explicit `fsb`): two-factor fit, `tspa(..., fsT = attr(fs, "fsT"),
  fsL = attr(fs, "fsL"), fsb = c(fs_a = 2, fs_b = 3))` → `attr(fit, "tspa_args")$fsb`
  equals `c(fs_a = 2, fs_b = 3)` and the rendered `tspaModel` carries the supplied
  intercepts.
- New test (explicit per-unit `fsT`/`fsL`): on a per-pattern (FIML) or per-cluster
  (`merMod`) result, pass `fsT`/`fsL` equal to the attributes **doubled** → the pooled
  `attr(fit, "fsT")` is doubled relative to the no-argument fit.
- Full suite green (other FIML pooling tests in `test-tspa_pooled.R` must keep passing —
  they exercise the derived path).

---

## Phase 2 — Row preservation in `get_fs()` (P1)

**Problem.** Output length per group is inferred from the largest scorable case index:
`n_cases <- max(unlist(lapply(blocks, function(b) max(b$case_idx))))` at
`R/get_fscore.R:599` (duplicated at `R/get_fscore_math.R:183` in
`augment_lav_predict()`). A fully-missing case after the last scorable case disappears.
Probed: 302 input rows with row 302 fully missing → 301 output rows. The existing
test documents the workaround: `tests/testthat/test-compute_fs_prod.R:243–253`
("a fully-missing row at the END of the data is silently dropped by `get_fs()`
(reported as a package bug)").

**Changes**

- `assemble_fs_blocks()` (`R/get_fscore.R:579`): add argument `group_n = NULL` — a named
  integer vector of per-group case counts in the same order as `blocks_by_group`.
  Use `n_cases <- group_n[group_labels[g]]` when provided, else the current
  `max(case_idx)` fallback. Unscrenable rows stay the existing NA scaffold: NA score/SE
  columns and NA `fs_pattern` label (the documented NA-row convention).
- `get_fs.lavaan` (`R/get_fs_methods.R:1260–1280`): pass
  `group_n = vapply(lavInspect(object, "data"), nrow, integer(1))` — the per-group row
  counts of the model frame are the ground truth.
- `augment_lav_predict()` (`R/get_fscore_math.R:183`): `n_cases <- nrow(fs_g)` —
  `lavPredict()` keeps one row per case with lavaan's NA convention for unscoreable
  factors, and `fs_g` is already in hand.
- Verify and leave untouched: `merge_local_fs()` compact path
  (`R/get_fs_methods.R:933`, `case_idx = seq_len(length(rows_g))` already covers all
  rows); `get_fs.merMod` (`R/get_fs_methods.R:1652`, hand-rolled one-row-per-cluster
  assembly, clusters are never dropped); mirt (fully-missing rows are "minted" as
  all-NA blocks by `get_fs.SingleGroupClass()`).

**Tests**

- New test (place with the existing FIML tests; check `test-get_fs.R` /
  `test-compute_fs_prod.R` conventions): leading, middle, **and trailing** fully-missing
  rows → `nrow(output) == nrow(input)`; unscoreable rows have NA scores, NA SEs, NA
  pattern labels, and NA product columns (for a `product =` result).
- Multigroup variant: trailing fully-missing row in one group → per-group row counts
  match the input per group.
- `augment_lav_predict()` equivalent: trailing fully-missing row → NA row present.
- Update `test-compute_fs_prod.R:243–253`: the trailing fully-missing row is now legal;
  replace the workaround (keep a middle row too) and assert row preservation.

---

## Phase 3 — Row-reorder detection + docs; mirt multigroup identity (P1)

### 3A. Row safety: detect + document (no S3 class)

**Problem.** `get_fs()` output is an ordinary data frame; row-specific quantities
(`fs_pattern` labels, per-obs attribute lists) are positional. Subset/reorder the rows
and the attributes no longer describe the right observations. `fs_indiv()`
(`R/fs_indiv.R:359–393` unified path; `:282–355` per-obs path) and `compute_fs_prod()`
trust the stored positional labels and silently use another observation's measurement
quantities. Probed: reversing a 301-row FIML result before `fs_indiv()` gave mismatched
SEs on 2 rows (with a random shuffle: 134).

**Changes**

- Add an internal alignment check co-located in `R/fs_indiv.R` (snake_case, e.g.
  `check_fs_row_alignment(fs, resolved)`): for each resolved block, compare the carried
  `fs_<v>_se` columns for the block's rows against the per-block SE values the existing
  `fs_row_cols()` loop already computes; beyond tolerance (relative ~1e-8), `stop()`
  with an actionable message, e.g.: "the row-specific attributes of this score result
  (fs_pattern / per-row fsL/fsT) do not match the data columns — the rows appear to
  have been reordered or subset after scoring. Keep the original row order, or re-score
  the subset with get_fs() before calling fs_indiv()/compute_fs_prod()."
- Wire the check into `fs_indiv()` (both the unified lavaan path and the `per_obs`
  mirt path) and into `compute_fs_prod()` (`R/compute_fs_prod.R` entry, after
  `resolve_fs_per_row()`). Skip the check when the carried `_se` columns are absent
  (hand-rolled input) to avoid false positives.
- New test: shuffled (and reversed) FIML result → `fs_indiv()` errors with the message;
  an unmodified result passes unchanged. Same for `compute_fs_prod()`.

**Docs**

- `get_fs()` roxygen: new "Row safety" details note — the output is an ordinary data
  frame; row-specific quantities are positional and do **not** follow subsetting or
  reordering; add a stable ID column **before** scoring if you need joins, and re-score
  a subset with `get_fs()` rather than subsetting a scored result.
- `fs_indiv()` roxygen: pointer to the same note.
- `vignettes/R2spa.Rmd`, "Reading the `get_fs()` output" section: short callout with the
  same advice.

### 3B. mirt multigroup identity

**Problem.** A multigroup mirt result carries a literal `group` column but no
`group_col` attribute (`R/get_fs_methods.R:2151–2166`), so the same stage-1 result
behaves differently downstream: `fs_indiv()` drops the group column (tested),
`combine_fs()` accepts it despite its single-group contract, and `tspa_mx_model()`
treats it as one structural group (documented, tested).

**Changes (minimal — guard + docs; no behavior change for currently-valid inputs)**

- `combine_fs()` (`R/combine_fs.R:142–148`): extend the existing `group_col` guard to
  also reject inputs with `"group" %in% names(x)` — same fail-fast message.
- New test: `combine_fs()` on a mirt MG input errors.
- `fs_indiv()` roxygen: document that mirt MG results carry no `group_col` attribute, so
  the `group` column is dropped from the output (matches
  `test-get_fs_mirt_multigroup.R:128–137` — keep that test).
- `tspa_mx_model()` roxygen: keep the existing single-structural-group note
  (`R/tspa_mx.R:43–46`).

---

## Phase 4 — Small validation + reference-doc fixes (roxygen; then `document()`)

1. **mirt `...` rejected** — `R/get_fs_methods.R:1830` (`get_fs.SingleGroupClass`) and
   `:2000` (`get_fs.MultipleGroupClass`): non-empty `...` is currently accepted and
   ignored (probed: `prior_cov = matrix(2)`, `product = "F1:F2"`, `method = "ML"` all
   silently no-op). Error naming the unsupported options passed. Follow the `merMod`
   prior-dots guard pattern (`R/get_fs_methods.R:1638–1644`). Test: each of the three
   probe calls errors.
2. **SE validation** —
   - `tspa()`: extend the `se_fs` check (`R/tspa.R:500–503`) to require finite **and
     nonnegative** values (negative SEs are currently squared into positive error
     variances; probed: `a = -0.2, b = -0.3` rendered fixed variances 0.04/0.09).
   - `tspa_mx_model()`: numeric `se_fs` path (`R/tspa_mx.R:572–578`) and SE-column path
     (`R/tspa_mx.R:533–550`) check `is.na` only — use `is.finite`.
   - Tests: negative / NA / Inf `se_fs` error on both routes.
3. **`grand_standardized_solution()` `model_list` docs** — the help says "a list of
   string variable describing the structural path model, in lavaan syntax", but the
   implementation takes a list of estimated matrices (e.g. `lavInspect(fit, "est")`;
   per-group list for MG — see `R/grandStandardizedSolution.R:119,145–147` and the test
   at `test-grandStandardizedSolution.R:89–99`). Fix the `@param` roxygen to describe
   the matrix layout, and add early validation with an actionable error when a supplied
   list doesn't look like estimated matrices.
4. **`tspa()` `reliability`/`se` help** (`R/tspa.R:153–157`): `reliability` says
   "Please use `se`" → should say use `se_fs`; clarify `se` (character values are
   forwarded to `lavaan::sem()`; numeric score-SE usage is deprecated in favor of
   `se_fs`).
5. **`compute_fscore()` return docs** (`R/get_fscore_math.R` roxygen): claims an
   N × p score matrix; with p indicators and q factors the result is N × q. Fix.
6. **`get_fs.default` error** (`R/get_fs_methods.R:1130–1135`): still says mirt support
   is "planned" although mirt methods exist — update the message (verify exact text).
7. **`block_diag()`** (`R/helper.R`): document the square-matrix and all-or-none
   dimname restrictions in roxygen (a single-matrix call currently fails because it is
   coerced to a list of scalars).
8. **`vignette()` lookups** — title-based lookups fail (`vignette()` resolves by
   topic/basename). Fix to topics: `R/get_fscore.R:292–296` (`"R2spa"`, `"scoring-matrices"`,
   `"efa-score"`, `"missing-data"`, `"sim-scoring-se"`), `R/tspa.R:297–303`, and
   `R/tspa_corrected_se.R:101–103` (`"corrected-se"`, `"correction-error"`). Grep
   `vignette("` across `R/` for stragglers. (The README's `vignette("R2spa")` is the
   correct pattern.)
9. **`tspa()` "Re-fitting with `lavaan::update()`"** (`R/tspa.R:135–142`): add that
   `update()` re-runs the stored lavaan model, so the R2spa correction
   (`corrected_se = TRUE`) and R2spa attributes (`tspa_args`, `tspa_corrected`, …) are
   not re-applied by `update()`.
10. Run `devtools::document()`, then `devtools::test()`.

---

## Phase 5 — Vignette, README, pkgdown content fixes

Verify each item against the code before editing; fix the text, keep executed numbers
unless the code changed.

1. `vignettes/R2spa.Rmd:164–166` — says `method = "mean"` is "not supported together
   with `corrected_fsT`, `vfsLT`, …", but `vfsLT` **is** supported for mean scores
   (`R/get_fs_methods.R:1160–1182`; tested at `test-get_fs_mean.R:370–404`, including
   the hand-verified loading-variance entry). Fix the sentence.
2. `vignettes/R2spa.Rmd:168` — heading "Single group, single factor" actually contains
   two separately scored constructs; rename (e.g. "Single group, two separate
   single-factor measurement models").
3. `vignettes/R2spa.Rmd:313` — `standardizedsolution()` → `standardizedSolution()`.
4. Loading-direction wording — `vignettes/R2spa.Rmd:124` and
   `vignettes/scoring-matrices.Rmd:96–98` (and the `fs_indiv()`/`get_fs()` roxygen):
   both the `_by_` columns and `fsL = S Λ` (rows = score names, columns = latent
   names; see `R/get_fs_methods.R:1538–1543`) are the **score-on-latent** loading —
   the coefficient of the latent in the score equation `score = fsb + fsL·η + error`
   (verified: the `fsL` value matches the regression of the score on the latent, not
   the reverse). The `_by_` column name `<latent>_by_fs_<score>` reads lavaan-style
   (as if latent-on-score) but carries the score-on-latent value. Fix the wording to
   say score-on-latent.
5. `vignettes/R2spa.Rmd:332–340` related vignettes — add: missing-data, corrected-se,
   product-factor-scores, sim-scoring-se.
6. `vignettes/R2spa.Rmd` — add the Phase 3A row-safety callout in "Reading the
   `get_fs()` output".
7. `vignettes/sim-scoring-se.Rmd`:
   - `:23`, `:32` — `vignette()` topic names (see Phase 4 item 8).
   - `:91–96` — rewrite the "Why not `std.lv = TRUE`?" note. Fixing the endogenous
     latent's residual variance is a **valid latent-scale identification choice**, not a
     misspecification; the unstandardized slope on that scale need not equal the
     unit-variance-scale slope (0.5/√(1 − 0.5²) ≈ 0.577 — probed: marker and
     residual-unit identification give the same **standardized** path). Do not teach
     that `std.lv` biases the path.
   - `:109–117` and `:204–208` — mean/sum scores: the item **weights** are fixed, but
     the implied score loadings `fsL = M Λ` remain free functions of the **estimated**
     item loadings (`test-get_fs_mean.R:391–404` verifies the nonzero loading-variance
     entry). So the corrected SE for mean scores propagates loadings **and** error
     variances. Bartlett: identity `fsL` is fixed, only error variances are free.
   - **Re-run the seeded sim loop** (chunk `sim-loop`, `set.seed(20260830)` — "a few
     minutes" per the vignette; this is the only fixture regeneration this plan
     allows, refreshing `vignettes/sim_scoring_se.RDS`) and confirm the cached tables
     still support the (corrected) narrative; update numbers/prose if they drift.
8. `vignettes/gr-std-coef.Rmd:56–70` — displayed equation
   `B_s = S_η⁻¹ B S_η⁻¹` must be `B_s = S_η⁻¹ B S_η` (the executed code at `:70`,
   `solve(S_eta) %*% beta %*% S_eta`, is correct).
9. `vignettes/multiple-factors.Rmd:42` — "Aε is the error covariance matrix" → it is
   the error **vector**; its covariance is A Θ Aᵀ.
10. `vignettes/tspa-vignette-mx.Rmd`:
    - `:219–234` — EAP measurement error: the code uses variance ρ·SE² (SE = √ρ·SE);
      the prose gives ρ·SE. Fix the prose to match the code and state the
      unit-prior-variance assumption behind the "total variance equals reliability"
      statement.
    - `:181–184` — "default-prior scores are mean zero" should distinguish the
      prior/model mean (zero) from the realized sample mean.
11. `vignettes/corrected-se.Rmd` —
    - Add an early "what is and is not corrected" box: the correction is partial — it
      propagates only the sampling uncertainty of the stage-1 `fsL`/`fsT` estimates;
      score values, `se_fs`, and `fsb` are held fixed, and the stages are treated as
      independent (source: `R/tspa_corrected_se.R:14–22`).
    - Label the separate-model example's block-diagonal `vfsLT` as an approximation
      (both models use the same subjects; the cross-model sampling covariance is
      ignored).
    - Remove the internal "PLAN 16" pointer from user-facing text (`:97`); describe
      the engines without the plan reference.
12. `vignettes/missing-data.Rmd` —
    - Add a worked **pooled `tspa()` route** section (before or alongside the OpenMx
      route): FIML stage 1 → `tspa()` with `reduce = "mean"` pooling the per-pattern
      `fsT`/`fsL` → inspect the pooled quantities → compare the pooled lavaan estimate
      with the exact OpenMx route and the joint model. For the block-missing example,
      the 10 non-correctable `dem65` scores must be dropped first (same convention as
      the OpenMx route, `:140–142`) — say so explicitly.
    - State where `vfsLT`/`corrected_fsT`/`reliability` are unavailable under missing
      data (`R/get_fs_methods.R:1208–1216`).
    - Fix the `VignetteIndexEntry` label (`:2–7`: "missing-data" → "2S-PA with Missing
      Data").
13. `vignettes/gr-std-coef.Rmd:2–5` — `VignetteIndexEntry` "gr-std-coef" → "Grand
    Standardized Coefficients".
14. `vignettes/multilevel.rmd`:
    - `:85` — **bug**: `e2 <- rnorm(clus_id1)` generates 500 errors recycled over the
      2,000 process-2 observations. → `e2 <- rnorm(clus_id2)`.
    - Re-run the seeded chunks (deterministic under the same seeds) and refresh the
      `:227–234` narrative with the new covariance tables (the unbiasedness claim
      should hold; report the replication SD honestly if the numbers move).
    - `:21–25` — OpenMx is `Suggests`-only and loaded unconditionally here; guard it
      like `vignettes/missing-data.Rmd:17–19` (`requireNamespace` + conditional
      `library` + reader-facing skip notice).
    - Note that `tspa_mx_model()` already returns a **fitted** model
      (`R/tspa_mx.R:266–267`); the subsequent `mxRun`/`mxTryHard` calls are refits or
      intentional retries — make that explicit in the prose (or drop the redundant
      refits, keeping `mxTryHard` where it is a deliberate convergence retry).
15. Accessibility/navigation —
    - `toc: true` for the long articles: corrected-se, tspa-vignette-mx, missing-data.
    - One-line prerequisites / optional-dependency note at the top of each article that
      uses `Suggests` (OpenMx, mirt, boot, numDeriv), including what is skipped when
      the dependency is absent.
16. `README.Rmd` —
    - Lead with the current canonical workflow instead of separate fits + `cbind()` +
      hard-coded SEs (`:48–68`): `get_fs(..., local = TRUE)` on one multi-factor model
      → `tspa(model, data = fs)` (auto-derived `se_fs`; see `:202–212` of
      `vignettes/R2spa.Rmd` for the evaluated version) → `standardizedSolution()`
      subset of the structural rows.
    - Add a link to the pkgdown site and the getting-started article, plus a short
      optional-dependencies note (OpenMx/mirt for the exact route and IRT).
    - `devtools::build_readme()` after editing.
17. `_pkgdown.yml` — currently only `url` + template. Add curated groups: articles
    organized by learner task (Getting started / Scoring strategies / Missing data &
    groups / Inference (SEs & standardization) / Extensions (products, IRT, random
    effects, growth, EFA) / Theory & simulation evidence) and a grouped reference
    section.
18. `NEWS.md` — one entry summarizing: explicit measurement inputs are now always
    honored in `tspa()`; fully-missing rows are preserved in `get_fs()` output;
    row-reorder guard in `fs_indiv()`/`compute_fs_prod()`; unsupported `mirt` options
    and negative SEs are now rejected.

---

## Phase 6 — Final verification

1. `devtools::document()` (after all roxygen changes) → confirm NAMESPACE/man diffs are
   only the intended ones.
2. `devtools::test()` — full suite green (~4,500 expectations; one existing
   negative-latent-variance warning is expected).
3. `devtools::check()` — expect 0 errors / 0 warnings / 0 NOTEs (`--no-manual`).
4. Knit every edited vignette to confirm clean evaluation — especially
   `multilevel.rmd`, `missing-data.Rmd`, `sim-scoring-se.Rmd`.
5. Re-run the review probes; all must now pass:
   - two-factor fit with `fsb = c(fs_a = 2, fs_b = 3)` → `tspa_args$fsb` = c(2, 3);
   - trailing fully-missing row → `nrow(output) == nrow(input)`;
   - reordered FIML result → `fs_indiv()` errors with the row-safety message.

---

## Deferred (next plan — do not start here)

- S3 score-result class with metadata-aware `[`/`subset` and stable case IDs — the full
  row-safety fix (Phase 3A here is the interim detect + document version).
- Replay fidelity: persist the recovered `sf_ld` in `tspa_args` so
  `do.call(tspa, attr(fit, "tspa_args"))` reproduces a derived single-factor fit
  (currently falls back to unit loadings; documented at `R/tspa.R:59–64`).
- R2spa-aware refit helper that preserves `corrected_se` through `update()`.
- Missing tutorials: multigroup `mirt`, external lavaan `prior_cov`/`prior_mean`,
  ID-aligned mixed-method `combine_fs()` (shuffled IDs, missing subjects), and a
  diagnostics/reporting guide.
- pkgdown stable-vs-development publishing policy (the published site shows 0.0.4 while
  the repo is 0.0.5; the workflow deploys only `main`/`master` pushes).

## Out of scope / do not touch

- `.quarantine/`, `legacy/`, `archive/` (except this plan file).
- Vignette RDS fixtures except `sim_scoring_se.RDS` (Phase 5 item 7).
- New `Imports`/`Suggests`; no tidyverse/dplyr/purrr; no `%>%` in `R/`.
