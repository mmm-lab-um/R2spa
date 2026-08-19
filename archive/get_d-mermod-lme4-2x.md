# Plan: `get_D()` — lme4-2.x RE covariance convention (merMod EB)

**Status: DONE** (2026-08-18, committed `966ff9b` on `refactor/core`).
**Key correction after step 1 (baseline):** the "missing `sigma^2` factor"
premise of issue 2 was a **misread** — the b730b53 clen restore had already
made the call-site parse exact under lme4 2.x, and the EB formulas carry the
explicit `sigma^2` scale themselves (the "1378.179 vs 1.435" probe was
exactly `sigma^2 = 960.4319` apart). Baseline measured per-term identity
`s^2 * tcrossprod(L_i) == VarCorr(term)` at max diff **0** on all fixtures,
and the final change is **numerically bit-identical** (RDS-verified, all
8 EB/ML outputs on 4 fixtures). What actually landed: (1) `get_D(object)`
self-contained per-term parse from `@cnms` split — the convention is
documented at the definition (R/get_fs_methods.R:~601) and the function can't
be mis-called; (2) new `tests/testthat/test-lme4_compat.R` canary (the
planned step 7) — per-term identity at 1e-15, `get_D == VarCorr[[1]]/s^2`
at 1e-15, warning-free multi-term parse; (3) term-1-only attribute reference
tests for the 2+1/2+2 fixtures at 1e-12 (in `test-get_fscore.R`); (4) roxygen
accuracy fixes (EB = first-term ranef; `corrected_fsT`/`format` declared-but-
unused notes for merMod). Suite 2339 pass / 0 fail; check 0/0/1 (OpenMx).
**Opened: 2026-08-18** — surfaced while landing the first-term Z fold perf change
in `get_fs_blocks.merMod()`; findings empirically verified (17-digit identities) by
r-architect on PRE-change code, i.e. all issues below are pre-existing.

Related landed work (same session, committed as `b730b53` on `refactor/core`):
- First-term Z fold (`Zden`) replacing the `as.matrix(getME("Z"))` densification —
  ~13–16× on large fixtures, output bitwise identical.
- Multi-term `clen` restore in the EB branch (`attr(theta, "clen") <-
  lengths(object@cnms, use.names = FALSE)` before `get_D`) — removes the 2+1
  `vec2mlist` replacement-length warning and the 2+2 hard crash
  (non-conformable 3×3 `D` for `num_re = 2`); first-term block bit-identical.
- New tests: multi-term fold regression (2+1) and 2+2 EB crash regression in
  `tests/testthat/test-get_fscore.R`.

## Problem

`get_D(theta)` in `R/get_fs_methods.R` (~:571) decodes lme4's random-effects
Cholesky parameters under the **lme4-1.x** convention:

```r
get_D <- function(theta) {
  L_D <- lme4::vec2mlist(theta, symm = FALSE)[[1]]
  tcrossprod(L_D)
}
```

Two independent issues on lme4 >= 2.x (2.0.6 installed):

1. **Multi-term theta parsing** — PARTIALLY ADDRESSED by the landed clen restore.
   lme4 2.x no longer attaches `"clen"` to `@theta` for multi-term fits; without
   it `vec2mlist()` falls back to a single-block parse of the mixed theta
   (2+1 → fractional block size → warning, first-term block still bit-correct;
   2+2 → 3+3 = 6 = 3·4/2 triangular-number coincidence → one mixed 3×3 block
   → crash). The clen restore in `get_fs_blocks.merMod` fixes both for this
   caller, but `get_D` itself still assumes a well-formed single/first block
   and stays fragile for any future caller. (Single-term theta lengths are
   triangular numbers, so the fallback was always exact there.)

2. **Scale/fill convention — THE WORK ITEM.** lme4 2.x stores the per-term
   Cholesky of `VarCorr / sigma^2`. Exact (17-digit) identity verified for
   multi-term, single-term p=2, and p=1 fixtures:

   ```
   VarCorr(term) == sigma^2 * tcrossprod(fill(upper-tri of term's theta params))
   ```

   Current `get_D` applies `tcrossprod` to the `vec2mlist(symm = FALSE)` parse
   (1.x convention) **without the sigma^2 factor**, so the EB
   `fsT` / `fsL` / `scoring_matrix` attributes returned for merMod fits are
   numerically inconsistent with lme4 2.x's `VarCorr`.

   Caveat to nail down first: the current single-term EB suite is
   self-consistent with the 1.x convention (sleepstudy `(Days|Subject)`
   asserts `fsT == sigma^2 * solve(Kz)` at 1e-12 and passes), so step 1 below
   must document exactly which quantities disagree with `VarCorr`-derived
   expectations before any value is touched. (r-architect's probe: e.g. a
   1378.179 vs 1.435 gap on a sleepstudy single-RE variance relative to
   `VarCorr/sigma^2`.)

## Scope of impact

- `get_fs(mobj, method = "EB")` for `merMod`: `fsT`, `fsL`, `scoring_matrix`
  attribute values (per-cluster arrays). Scores (`fs_u*`) are UNAFFECTED (they
  come from `getME("b")`); `method = "ML"` UNAFFECTED (no `D`); lavaan path
  UNAFFECTED.
- Stage 2 (`R/tspa.R`) consumes `fsT`/`fsL`/`ev_*`/`ecov_*` — latent SE
  correction for EB merMod fits changes numerically.
- Vignettes: check `vignettes/scoring-matrices.Rmd` and any merMod EB numeric
  narratives for stale values.

## Proposed approach

1. **Baseline characterization (no code change).** For fixtures (sleepstudy
   p=2; sleepstudy p=1; 2+1 crossed; 2+2 crossed), compute:
   current `D` vs the 2.x-convention per-term matrix
   `sigma^2 * tcrossprod(L_term)`, pinning down exactly what `vec2mlist`
   returns vs what `get_D` needs under 2.x (read lme4 2.0.6
   `VarCorr.merMod`/`varCorr` source to get the fill/upper-tri direction
   exactly right). Record per-quantity diffs (fsT/fsL/scoring_matrix) in this
   file.
2. **Rework `get_D`.** Parse per-term blocks with the clen idiom (decide:
   restore-at-call-site as today, or change `get_D` to take `object` and do
   the restore internally so the function can't be mis-called). Keep output
   bitwise for single-term iff step 1 shows equivalence; otherwise accept the
   value change and update references.
3. **Update tests** (`tests/testthat/test-get_fscore.R`): re-point merMod EB
   `fsT`/`fsL`/`scoring_matrix` references at `VarCorr`-derived values
   (1e-12); add multi-term EB reference tests (2+1, 2+2) using per-term
   `VarCorr`. No tolerance bumps beyond 1e-12.
4. **Stage-2 sweep.** Find and re-run every consumer of merMod EB `fsT`
   (grep tests + vignette sources for `method = "EB"` merMod usage;
   `tspa()` smoke on a merMod EB fit before/after).
5. **Docs.** If `get_fs`'s `@return` fsT wording or the scoring-matrices
   vignette describe EB SE semantics that change, update roxygen/vignette
   (then `document()`; vignette re-render only for the affected fixture).
6. **Lifecycle.** `load_all` → `document` → `test` → `check` (expect 0/0/1,
   the OpenMx NOTE).
7. **Version canary (optional, nice).** `R/lavaan_compat.R` sets the pattern of
   a tested-to version canary; consider an lme4 version canary asserting the
   2.x identity on the installed lme4 so a future lme4 convention change fails
   loudly instead of silently.

## Exit criteria

- `get_D` output equals `sigma^2 * tcrossprod(L_term)` for each term at
  machine precision on all fixtures (17-digit check recorded).
- Full suite green with updated references; 0 failures, 0 suite warnings.
- `devtools::check()`: 0 errors / 0 warnings / 1 NOTE (OpenMx).
- This file updated to DONE with the measured before/after table.

## Risks

- (Resolved) EB merMod stage-2 latent SEs shift for users — the baseline
  showed NO value shift (bit-identical), so the commit message records
  "no user-visible value change"; vignettes checked and confirmed
  bit-unchanged (no re-render).
- (Resolved) lme4 1.x vs 2.x dual convention: source doc comment records
  that both eras store the Cholesky of the scaled covariance (`get_D(θ)·σ²
  ≡ VarCorr` parity per `archive/PLAN_03_mermod_scoring_matrix.md`); the
  `test-lme4_compat.R` canary (step 7) is the drift guard.
- (Resolved 2026-08-18: the Z-fold + clen changes are committed as `b730b53`;
  this work item now diffs cleanly against it.)
- Noted, not fixed (pre-existing, logged here for the backlog): multi-term
  merMod EB attributes are the term-1-only conditional quantities while
  lme4's full random-effect posterior is joint (posterior-cov gap up to
  ~4.5e-2 on the 2+2 fixture); `get_fs.merMod()` declares `corrected_fsT` /
  `vfsLT` / `format` but never routes them (now documented as ignored /
  not used).

## References

- `get_D`: `R/get_fs_methods.R` (function `get_D`, now `get_D(object)`, ~:598)
- convention doc comment: `R/get_fs_methods.R` ~:601-639
- lme4 canary: `tests/testthat/test-lme4_compat.R`
- New multi-term/term-1 attribute tests: `tests/testthat/test-get_fscore.R`
  (~:939–:1080; term-1 reference helper `check_eb_term1_attrs()`)
- Z-fold plan/spec: session plan `get_fs-mermod-ml.md` (sibling file)
- r-architect findings report, 2026-08-18 (this session)

## Completion record (2026-08-18, `966ff9b`)

Baseline table (step 1) — per-term identity `s^2 * tcrossprod(L_i) ==
VarCorr(x)[[term]]` measured at **max diff 0** on all fixtures (sleepstudy
p=2 REML + ML, `(1 | Subject)` p=1, 2+1, 2+2), confirming the b730b53 clen
restore already made the first-term value exact. Representative pre-change
values: p=1 `VarCorr = 1378.179` vs `get_D = 1.43492`, ratio = `sigma^2 =
960.4319` exactly (the "discrepancy" that motivated this item was the
expected scale factor, not a parse artifact).

| Step | Outcome |
|---|---|
| 1. baseline | identity at 0; single-term drift verdict: **none** |
| 2. rework | `get_D(object)` — `@cnms` block split, lower-tri column-major fill, `tcrossprod`; first-term SCALED contract documented in-source |
| 3. tests | canary `test-lme4_compat.R` (3 tests, 9 expectations @1e-15 + warning-free parse); term-1 attribute reference tests on 2+1/2+2 at 1e-12; stale comment updated |
| 4. stage-2 sweep | no R/ or test consumer routed merMod EB `fsT`/`fsL` into `tspa()`; vignette `scoring-matrices.Rmd` values bit-unchanged (no re-render) |
| 5. docs | roxygen accuracy fixes (EB first-term ranef; `corrected_fsT`/`format` declared-but-unused for merMod) |
| 6. lifecycle | suite 2339 pass / 0 fail / 0 warn / 0 skip; check 0/0/1 (OpenMx) |
| 7. canary | landed (test-only, ~95 lines) |

Exit-criteria note: criterion 1 as written (`get_D == s^2 * tcrossprod(L)`)
is the lme4-side `VarCorr` identity; `get_D` returns the SCALED first term
(`VarCorr[[1]] / s^2`) by contract — the canary pins both sides at 1e-15.
