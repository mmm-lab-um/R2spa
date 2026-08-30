# Plan: Quarantine get_fs()-/tspa()-dependent code (`.quarantine/`)

Create `.quarantine/{R,tests,vignettes}` and move all files that **consume** `get_fs()` /
`tspa()` (and are therefore at risk while those two are being revised/refactored — active
work: PLAN 04, uncommitted) out of the package. The core files being refactored stay.

**Decisions (user-confirmed 2026-08-17):**
- `R/lavaan_compat.R` + `tests/testthat/test-lavaan_compat.R` **stay** in the package
  (PLAN-04 lavaan-drift canary layer; independent of get_fs()/tspa()).
- Core-workflow vignettes (`R2spa.Rmd`, `multiple-factors.Rmd`, `correction-error.Rmd`
  + its RDS fixtures, `tspa-growth-vignette.Rmd`, `efa-score.Rmd`, `scoring-matrices.Rmd`)
  **stay** in `vignettes/` as the working spec for the refactor.
- No `DESCRIPTION` change: `OpenMx` remains in `Imports` until the OpenMx path is re-integrated.

**Baseline state** (verify before starting): `git status` should show the uncommitted
PLAN-04 work (modified `R/tspa.R`, `R/get_fscore.R`, `R/get_fs_methods.R`,
`R/tspa_corrected_se.R`, `R/grandStandardizedSolution.R`, `tests/testthat/test-tspa.R`,
several vignettes; untracked `R/lavaan_compat.R`, `tests/testthat/test-lavaan_compat.R`,
`tests/testthat/test-tspa_render.R`). All line numbers below are from that working tree
(HEAD `32fc817` + uncommitted). **If numbers have drifted, trust the verbatim anchor
strings, not the line numbers.**

**Hard rules (AGENTS.md):** never hand-edit `NAMESPACE`/`man/*.Rd`; never `library()`/
`require()` in function bodies; after any roxygen change run `devtools::document()`;
test with `devtools::test()` after functional changes; `devtools::check()` before done
(slow — OpenMx compiles, vignettes pull `Suggests`).

---

## Progress

- [x] **Phase 1 — move whole files** (R/, test files, vignettes + fixtures) — done 2026-08-17 (17 staged renames + 8 untracked html)
- [x] **Phase 2 — extract embedded test blocks** (5 in-place edits, 2 new quarantined test files, 2 appends) — done 2026-08-17
- [x] **Phase 3 — hygiene** (`.Rbuildignore`, orphaned `@importFrom`, dangling roxygen links) — done 2026-08-17 (plan-file already covered by `^_PLAN_.*\.md$`; only `^\.quarantine$` added)
- [x] **Phase 4 — regenerate + verify** (`load_all` → `document` → `test` → `check`) — done 2026-08-17: tests 707 pass / 0 fail; check = 0 error, 0 warning, 1 NOTE ("OpenMx in Imports not imported from" — expected consequence of the user decision to keep OpenMx in Imports until re-integration)
- [x] **Re-integration (partial) — grand standardization** — done 2026-08-22, branch `rejoin/grand-std-sol` (parent `126a63a`): `git mv` of `R/grandStandardizedSolution.R` + `tests/testthat/test-grandStandardizedSolution.R` back (pure rename; zero code edits); `document()` → `test()` (3377 pass / 0 fail) → `check()` (0 errors / 0 warnings / 0 notes, as-cran). `vignettes/gr-std-coef.Rmd` deliberately not restored in this commit (reintroduced separately in a follow-up commit — see Re-integration notes).

---

## Verified dependency facts (why this set moves)

- No remaining (staying) `R/` file **calls** any quarantined function. `get_fs`/`tspa`
  core files reference quarantined functions only in roxygen examples/comments.
- `R/tspa.R:tspa_sf_alias()` (product-score column auto-alias) has **no hard dependency**
  on `get_fs_int` (pure column-name handling for `fs_a:fs_b` columns) → **stays**.
- Confirmed consumers of quarantined functions in-package:
  - `get_fs_int`: `tests/testthat/test-tspa_render.R:374` (in `int_setup()`).
  - `tspa_mx_model`: `tests/testthat/test-tspa.R:179-235`,
    `tests/testthat/test-get_fscore.R:684`.
  - `vcov_corrected`: `tests/testthat/test-tspa_render.R:322-344`,
    `tests/testthat/test-get_fs_priors.R:175-205`.
  - `grandStandardizedSolution`/`grand_standardized_solution`:
    `tests/testthat/test-grandStandardizedSolution.R` (whole file),
    `tests/testthat/test-lavaan_compat.R:107-121`.
  - `tsp_*` wrappers of `R/lavaan_compat.R`: only `R/tspa_corrected_se.R:23` and
    `R/grandStandardizedSolution.R:84-99` (both move) — so the compat layer becomes
    unused-but-kept internal code (user decision: keep).
- `R/tspa_mx.R` is the **only** OpenMx consumer in `R/` (all 6 `importFrom(OpenMx, ...)`
  entries come from its roxygen).
- `pnorm`/`qnorm`/`tail`/`combn` are used **only** by quarantined files
  (`grandStandardizedSolution.R`, `get_fs_int.R`) → their `importFrom` entries vanish
  on regeneration; no remaining bare call sites.
- **Orphaned-import gotcha:** `importFrom(lavaan, vcov)` is currently declared *only* in
  `R/grandStandardizedSolution.R:19`, but `R/get_fscore_math.R:427,518` make bare
  `vcov(fit)` calls. Phase 3 must re-declare it or those calls break.
- `R/globals.R` (`utils::globalVariables("op")`) exists for
  `grandStandardizedSolution.R:96` (`subset(..., op == "~")`). After the move the
  suppression is vestigial but harmless → **keep**, no action.
- No `tests/testthat/_snaps/` snapshot files exist; nothing to move there.
- `vignettes/.gitignore` ignores `*.html` and `*.R` → the `.html` build artifacts are
  **untracked**: plain `mv`, not `git mv`. `.Rmd`/`.rmd`/`.RDS` are tracked: `git mv`.
- Untracked `R2spa.Rcheck/` dir (leftover check run) — do not touch.

---

## Phase 1 — Move whole files

### 1.1 Setup

```sh
mkdir -p .quarantine/R .quarantine/tests .quarantine/vignettes
```

### 1.2 R/ → .quarantine/R/ (4 files; all exports they carry leave NAMESPACE in Phase 4)

```sh
git mv R/get_fs_int.R               .quarantine/R/
git mv R/tspa_mx.R                  .quarantine/R/
git mv R/tspa_corrected_se.R        .quarantine/R/
git mv R/grandStandardizedSolution.R .quarantine/R/
```

Export inventory leaving the package: `get_fs_int`, `tspa_mx_model`, `vcov_corrected`,
`grand_standardized_solution`, `grandStandardizedSolution`; imports: `OpenMx` (6 fns),
`stats` (`pnorm`, `qnorm`), `utils` (`combn`, `tail`), `lavaan` (`vcov` — re-declared in
Phase 3; `lav_func_jacobian_complex` — drops safely, only namespaced `lavaan::…` calls
remain in staying code), plus the S3-less internal helpers they define
(`check_inputs`, `make_mx_ld/vc/int`, `update_tspa`, `veta`, `eeta`, `veta_grand`,
`.fill_matrix_list`, `.combine_est`, `std_beta_est`, `grand_std_beta_est` — all of which
are used only by their own file, verify no other references before moving if in doubt:
`grep -rn "veta\|eeta\|make_mx\|update_tspa" R/`).

### 1.3 tests/testthat/ → .quarantine/tests/ (whole files)

```sh
git mv tests/testthat/test-get_fs_int.R             .quarantine/tests/
git mv tests/testthat/test-grandStandardizedSolution.R .quarantine/tests/
```

### 1.4 vignettes/ → .quarantine/vignettes/ (8 vignettes + their fixtures)

```sh
git mv vignettes/get_fs_int-vignette.Rmd .quarantine/vignettes/
git mv vignettes/categorical-interaction.Rmd .quarantine/vignettes/
git mv vignettes/reliability.Rmd            .quarantine/vignettes/
git mv vignettes/corrected-se.Rmd           .quarantine/vignettes/
git mv vignettes/tspa-vignette-mx.Rmd       .quarantine/vignettes/
git mv vignettes/missing-data.Rmd           .quarantine/vignettes/
git mv vignettes/multilevel.rmd             .quarantine/vignettes/
git mv vignettes/gr-std-coef.Rmd            .quarantine/vignettes/
# fixtures (each referenced by exactly the vignette above it)
git mv vignettes/sim_results_reliability.RDS    .quarantine/vignettes/
git mv vignettes/boo_joint.RDS                  .quarantine/vignettes/
git mv vignettes/boo_separate.RDS               .quarantine/vignettes/
# untracked .html builds (gitignored) — plain mv
mv vignettes/get_fs_int-vignette.html  .quarantine/vignettes/
mv vignettes/categorical-interaction.html .quarantine/vignettes/
mv vignettes/reliability.html          .quarantine/vignettes/
mv vignettes/corrected-se.html         .quarantine/vignettes/
mv vignettes/tspa-vignette-mx.html     .quarantine/vignettes/
mv vignettes/missing-data.html         .quarantine/vignettes/
mv vignettes/multilevel.html           .quarantine/vignettes/
mv vignettes/gr-std-coef.html          .quarantine/vignettes/
```

(`*.html` names may differ slightly — if a `mv` reports "No such file", just skip it;
they are regenerable and gitignored.)

Fixture-to-vignette mapping (verified via `readRDS` grep):
`corrected-se.Rmd` ← `boo_separate.RDS`, `boo_joint.RDS`;
`reliability.Rmd` ← `sim_results_reliability.RDS`.
`sim_correction-error.RDS` / `sim_correction-error-multi.RDS` / `sim_output.csv`
belong to **staying** vignettes — do not move.

### 1.5 Phase 1 check

```sh
git status --short   # 4 R/ + 2 tests/ + 8 vignettes + 3 RDS staged as renames (R)
ls .quarantine/R .quarantine/tests .quarantine/vignettes
```

---

## Phase 2 — Extract embedded test blocks

Quarantined test files go to `.quarantine/tests/`. Every extracted block gets a
provenance comment header, e.g.:

```r
# Quarantined with R/tspa_mx.R (see QUARANTINE_PLAN.md).
# Extracted from tests/testthat/test-tspa.R, "# Compare to Mx" block
# (lines 159-239 as of 2026-08-17) plus its setup (lines 116-157).
```

Extracted blocks reference objects defined earlier in their source files — copy the
minimal setup listed below so each quarantined file is self-contained when the code is
re-integrated.

### 2.1 New `.quarantine/tests/test-tspa_mx.R`

Append, in this order:

1. **From `test-tspa.R`:**
   - Setup copy: lines 116-157 — anchors: `cfa_3var1 <- '` … through the `tspa_3var
     <- tspa(` call ending in `se_fs = c(ind60 = 0.1213615, dem60 = 0.6756472,
     dem65 = 0.5724405) )`. (Includes `fs_dat_3var` and `tspa_3var`, both of which the
     remaining part of `test-tspa.R` also uses — **keep** them there; a copy goes here.)
   - Block copy: lines 159-239 — anchor start `# Compare to Mx`, anchor end: the lone
     `}` closing `if (requireNamespace("umx", quietly = TRUE)) {` immediately before
     the blank line + `########## Testing section #############`.
2. **From `test-get_fscore.R`** (provenance header noting the extra setup):
   - Setup copy: `hs_model_2 <- ' visual =~ x1 + x2 + x3` definition (line 240), then
     lines 637-642 (`hs <- HolzingerSwineford1939`, `set.seed(1334)`, the 3 NA-injection
     lines), then lines 650-654 (`cfa_fiml <- cfa(` … `a2 <- augment_lav_predict(
     cfa_fiml, method = "Bartlett")`).
   - Block copy: lines 675-700 — anchor start `if (requireNamespace("umx", quietly =
     TRUE)) {` containing `lcov_umx <- umxLav2RAM(`, anchor end: file end (line 700,
     the closing `}` of the if-block; this is the last line of the file).

Then **edit in place**:
- `test-tspa.R`: delete lines 159-239 (inclusive).
- `test-get_fscore.R`: delete lines 675-700 (inclusive; the `cfa_fiml` def at 650-653
  and `a2` at 654 **stay** — used by the remaining test at ~line 662).

### 2.2 New `.quarantine/tests/test-vcov_corrected.R`

1. **From `test-tspa_render.R`:**
   - Setup copy: line 17 — `mod2g <- "visual =~ x1 + x2 + x3\nspeed =~ x7 + x8 + x9"`.
   - Block copy: lines 322-344 — the single
     `test_that("tspa_call stays re-callable and vcov_corrected() runs on an MG fit", {`
     block ending in `})` just before `## Interaction (product-score) auto-alias`.
2. **From `test-get_fs_priors.R`:**
   - Setup copy: lines 1-18 — `library(lavaan)`, `library(lme4)`, the
     `########## Single-group example ##########` line, `prior_model <- '`,
     `prior_fit <- cfa(prior_model, data = PoliticalDemocracy)`, `pm <- c(ind60 =
     0.3, dem60 = -0.4)`, `pc <- matrix(c(1.2, 0.25, 0.25, 0.8), 2, 2, …)`.
   - Block copy: lines 175-205 — the
     `test_that("vcov_corrected() works with prior-adjusted factor scores", {` block
     ending in `})` just before `test_that("data.frame entry point matches lavaan
     method entry point"`.

Then **edit in place**:
- `test-tspa_render.R`: delete lines 322-344 (inclusive).
- `test-get_fs_priors.R`: delete lines 175-205 (inclusive).

### 2.3 Append to `.quarantine/tests/test-get_fs_int.R`

From `test-tspa_render.R`, the product-score section:
- Lines 346-405: the section header (`## Interaction (product-score) auto-alias
  ------…`), `int_setup <- function(n = 500) {` (calls `get_fs_int` at its line 374),
  and `test_that("Product-score columns are auto-aliased (no manual rename needed", {`.
- Lines 419-427: `test_that("tspa_sf_alias is a no-op when the score column already
  exists", {` … `})` (end of file).

### 2.4 Keep (do NOT move): `test-tspa_render.R` lines 407-417

`test_that("Ambiguous product-score candidates are a clear error", {` is
self-contained (manual data, no `get_fs_int`) and tests **core** `tspa.R` aliasing —
it stays. After deletions in 2.2/2.3, the remaining file should read: …line 320
`})`, blank, the section header (346-347), the ambiguous-candidates test (407-417).
(Adjust the section header if you prefer, e.g. keep it as-is.)

### 2.5 Append to `.quarantine/tests/test-grandStandardizedSolution.R`

From `test-lavaan_compat.R`, the wrapper A/B block for grandSS:
- Setup copy: lines 9-18 — `canon_mod <- …`, `canon_fit <- sem(…)`, `mg_mod <- …`,
  `mg_fit <- sem(…)`.
- Block copy: lines 107-121 — `test_that("grand_standardized_solution output is
  unchanged by the wrapper", {` … `})`.

Then **edit in place**:
- `test-lavaan_compat.R`: delete lines 107-121 (inclusive). Keep the file's
  `canon_fit`/`mg_fit` setup (lines 9-18) — the remaining wrapper canaries use them.

### 2.6 Phase 2 check

Every remaining test file must parse with no dangling references:

```sh
Rscript -e 'for (f in list.files("tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)) parse(f)'
# quick eyeball:
grep -n "tspa_mx\|vcov_corrected\|get_fs_int\|grandStandardized\|mxMatrix\|umxLav2RAM" tests/testthat/*.R
# expected: only the kept "Ambiguous product-score candidates" test area (no matches at all)
```

---

## Phase 3 — Hygiene

1. **`.Rbuildignore`** — add (one line, keep existing entries intact):

   ```
   ^\.quarantine$
   ^QUARANTINE_PLAN\.md$
   ```

2. **Re-declare the orphaned `vcov` import** (required — `R/get_fscore_math.R:427,518`
   call `vcov()` bare; the only `@importFrom lavaan vcov` was in the removed
   `grandStandardizedSolution.R`): in `R/get_fscore.R`, extend the existing roxygen
   import block of `get_fs()` (lines 86-88: `@importFrom lavaan cfa sem`,
   `@importFrom lavaan lavInspect lavTech coef`, `@importFrom stats setNames`) with:

   ```
   #' @importFrom lavaan vcov
   ```

3. **Dangling roxygen links** (removed Rd topics would be check warnings):
   - `R/get_fscore.R` `@param vfsLT` (lines ~35-38): text references
     `[vcov_corrected()]` and `[tspa()]` — reword to drop the `vcov_corrected()`
     cross-ref (e.g. "returned as attribute `vfsLT`; used for second-order SE
     correction of 2S-PA results" — no link).
   - `R/get_fscore_math.R:101` (roxygen comment around the `fs_matrices`/scoring
     output docs): mentions `[tspa_mx_model()]` — reword to plain `tspa_mx_model()`
     without link brackets.
   - Grep for stragglers in **staying** files only:
     `grep -rn "vcov_corrected\|tspa_mx_model\|grandStandardizedSolution\|get_fs_int" R/`
     — remaining hits should be comments only (`R/tspa.R:188,482` mention
     `get_fs_int` in code comments — fine, no roxygen link; leave or reword at
     discretion).

---

## Phase 4 — Regenerate + verify (AGENTS.md ordering, no skips)

1. `Rscript -e 'devtools::load_all()'`

2. `Rscript -e 'devtools::document()'` — then **inspect the diff** (never hand-edit):

   - `NAMESPACE`: exports reduced to
     `augment_lav_predict`, `block_diag`, `compute_fscore`, `fs_to_group_list`,
     `get_fs`, `get_fs_lavaan`, `get_fs_lmer`, `tspa` (+ the 4 `get_fs` S3 methods);
     `importFrom(OpenMx, …)` block gone; `stats` block reduced to `setNames`;
     `utils` block gone entirely; `lavaan` block loses nothing **unless** the Phase 3
     step 2 was missed (`vcov` must remain).
   - `man/`: `get_fs_int.Rd`, `tspa_mx_model.Rd`, `vcov_corrected.Rd`,
     `grand_standardized_solution.Rd` removed; `get_fs.Rd` updated (new
     `@importFrom`/rewritten `vfsLT` text).

3. `Rscript -e 'devtools::test()'` — expect all remaining files green
   (`test-assemble_fs_blocks`, `test-compute_fscore`, `test-fs_converters`,
   `test-get_fscore`, `test-get_fs_priors`, `test-lavaan_compat`,
   `test-lavPredict_equivalence`, `test-tspa`, `test-tspa_render`). Quarantined files
   are outside `tests/testthat/` and do not run.

4. `Rscript -e 'devtools::check()'` — slow (OpenMx compiles; 6 remaining vignettes
   build). Acceptable findings: pre-existing WARNING/NOTE items on the uncommitted
   baseline (see PLAN 04 notes: 2 WARNING / 2 NOTE). Anything **new** (broken vignette,
   link warnings, no-visible-global) must be fixed.

5. `git status` / `git diff --stat` sanity: working tree should show only the
   quarantines (staged renames) + the small Phase-2/3 edits; then stop and let the
   user commit (do not commit unless asked).

---

## Expected final state

```
.quarantine/
  R/           get_fs_int.R  tspa_mx.R  tspa_corrected_se.R  grandStandardizedSolution.R
  tests/       test-tspa_mx.R  test-vcov_corrected.R  test-get_fs_int.R
               test-grandStandardizedSolution.R
  vignettes/   get_fs_int-vignette.Rmd  categorical-interaction.Rmd  reliability.Rmd
               corrected-se.Rmd  tspa-vignette-mx.Rmd  missing-data.Rmd
               multilevel.rmd  gr-std-coef.Rmd
               sim_results_reliability.RDS  boo_joint.RDS  boo_separate.RDS
               *.html (untracked builds)
R/             get_fscore.R  get_fs_methods.R  get_fscore_math.R  tspa.R
               lavaan_compat.R  helper.R  globals.R
vignettes/     R2spa.Rmd  multiple-factors.Rmd  correction-error.Rmd  tspa-growth-vignette.Rmd
               efa-score.Rmd  scoring-matrices.Rmd  (RDS fixtures: sim_correction-error*.RDS,
               sim_output.csv, *.html)
```

## Re-integration notes (for later)

- Restore by `git mv`-ing files back into `R/`, `tests/testthat/`, `vignettes/` once
  `get_fs()`/`tspa()` contracts are settled, then: `document()` → `test()` → `check()`.
- The quarantined test files are self-contained (setup copied) with provenance
  comments; re-run them as-is against the new contracts.
- If `get_fs()`/`tspa()` output contracts change, the quarantined code (especially
  `get_fs_int()` column conventions and `vcov_corrected()` attribute usage) must be
  re-verified before unquarantining.
- `tspa.R:tspa_sf_alias()` continues to work for hand-built `fs_a:fs_b` columns; a
  minimal core test for aliasing (ambiguous-candidates case) is retained in
  `test-tspa_render.R`.
- **2026-08-22 — `grandStandardizedSolution.R` re-joined** (R file + test file;
  branch `rejoin/grand-std-sol`, merged into `refactor/core`): the quarantined code
  was contract-stable, so a pure `git mv` restored it with zero edits and the
  self-contained test file (incl. the Phase 2.5 wrapper A/B block) passed
  unmodified. Still quarantined: `R/get_fs_int.R`, `R/tspa_corrected_se.R` (+ their
  test files, and the stale `.quarantine/tests/test-tspa_mx.R` leftover from the
  `bcd42a3` rewrite-style re-integration) and 5 vignettes + 3 RDS fixtures
  (`vignettes/gr-std-coef.Rmd` was reintroduced separately in a follow-up commit).
