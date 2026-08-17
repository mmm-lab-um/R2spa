# R2spa — Open Issues & Follow-ups

Tracks remaining work from completed plans (see `archive/`). Update status as
items are resolved; move finished items to the **Closed** section with the
date and commit/PR reference.

**Last updated:** 2026-08-17 (PLAN 04 (tspa partable) — lavaan compat module + R2spa-owned stage-2 model schema/renderer, product-score auto-alias, CI lavaan axis; baseline check findings cleared — S3 `get_fs()` arg rename, Rd `fsm`/`...` docs, `Matrix` → Suggests, top-level `.Rbuildignore` exclusions; full check 0 errors / 0 warnings / 0 notes)

## Open

| # | Priority | Issue | Where | Source | Status |
|---|----------|-------|-------|--------|--------|
| 4 | P2 | **`reliability = TRUE` + single-group multi-factor + default `unified` format → hard error** ("reliability is only supported for unidimensional models") instead of the documented warning. Guard at `get_fs_methods.R:183-184` checks `length(attr(out, "fsb")) > 1`, but in unified format `fsb` is always a length-1 list, so `multifactor` never fires. | `R/get_fs_methods.R` (`get_fs.lavaan()` reliability block) | PLAN 01, remaining issue 4 | open |
| 5 | P2 | **`correct_evfs()` is ~100% of `corrected_fsT = TRUE` runtime** (0.167 s of 0.168 s, 4-LV model). `q × (npar+1) × 2 + 1` `lavInspect` file-reads per call (57 traced on a 13-param model). Fix: one `lav_func_jacobian_complex()` over the full `a` matrix (q× fewer evaluations) + hoist `lavInspect("free")`/`("est")` out of the per-evaluation `compute_fspars()` body. | `R/get_fscore_math.R` (`correct_evfs()` :366-401, `compute_fspars()` :289-290) | PLAN 01, remaining issue 5 | open |

## Follow-ups (future work, explicitly out of scope for their plans)

| # | Item | Source | Status |
|---|------|--------|--------|
| F1 | Wire `tspa()` to accept unified/list factor-score shapes natively (via `fs_to_group_list()` or dual support) — currently relies on attributes; the `format` default change shifts multigroup output to list-valued attributes that `tspa()` must keep parsing. | PLAN 01, remaining issue 7 | deferred |
| F2 | `mirt` S3 method (`get_fs.SingleGroupClass`) — architecture (blocks + assembler) is backend-agnostic and ready. | PLAN 01, remaining issue 8 | deferred |

## Closed

| # | Issue | Closed | Reference |
|---|-------|--------|-----------|
| 1 | **merMod column-name regression** — restored pre-refactor `u0_eb`-style *column names* in the merMod path via `legacy_names` switch; `get_fs_lmer()` defaults `legacy_names = TRUE` so `vignettes/multilevel.rmd` (`tspa_mx_model`) works unchanged; default `get_fs()`/`get_fs.merMod()` now use `fs_u0`-style names. Legacy output is name-compatible (not byte-identical) with the pre-refactor result — extra `_se` columns, `fsL`/`fsT` attributes, NULL row names (delta documented on `get_fs_lmer()`/`get_fs.merMod()`). Fixed the related overwritten-`fsT`-rownames + duplicate-`re_names` bugs (item 7 sub-bullets). | 2026-08-15 | PLAN 02, Step 2 |
| 2 | **Vignette breakage on `format = "unified"`** (verified failure set 3/13: `corrected-se.Rmd`, `multilevel.rmd`, `tspa-vignette-mx.Rmd`) — `multilevel.rmd` fixed by item 1; `R/tspa.R` now validates `fsT`/`fsL` group-count consistency (plain matrix = 1 group, so single-group length-1 list attributes may be mixed with plain matrices, e.g. Bartlett identity `fsL`) with a clear mismatch error, and `tspa_mf()` accepts all single-group shape combinations; `tspa-vignette-mx.Rmd` uses `format = "list"` for its direct attribute arithmetic. 6 new regression tests. **13/13 vignettes build** (verified 2026-08-16). | 2026-08-16 | PLAN 02, Step 3 |
| 3 | **`get_fs(data = matrix)` dispatch** — already resolved in the working tree: `get_fs.data.frame()` re-dispatches matrices (`is.matrix(data)` → `as.data.frame(data)` → `get_fs()`). Live check: `get_fs(as.matrix(...))` is `identical()` to `get_fs(data.frame)` incl. attributes. Locked in with a regression test in `test-get_fscore.R`. | 2026-08-16 | PLAN 02, Step 4 |
| 6 | **Missed cheap-lookup `lavInspect` sites** — all 3 replaced with direct slot access: `tspa_corrected_se.R` (`@Data@ngroups`), `get_fs_methods.R` (`unlist(@Data@norig)`), `grandStandardizedSolution.R` (`unlist(@Data@nobs)`). Slot values proven equal to `lavInspect()` output (SG + MG); MG standardization path A/B'd by existing tests; MG reliability `norig` path A/B'd live. One further cheap lookup noted out of scope: `get_fscore_math.R:128` (`"ngroups"` in `augment_lav_predict()`). | 2026-08-16 | PLAN 02, Step 5 |
| 7 | **Dead-code cleanup** — remaining pieces removed: commented-out `# fs_se[is.nan(fs_se)] <- 0` line in `augment_fs()` (`get_fscore.R`) and the unreachable `is.data.frame()` guard in `get_fs.data.frame()` (`get_fs_methods.R`); the merMod `fsT_j` rownames + duplicate `re_names` sub-bullets were fixed in Step 2. | 2026-08-16 | PLAN 02, Step 6 |
| 8 | **merMod `scoring_matrix` + Z-design fix + vignette** — `get_fs.merMod()` now emits a per-cluster `scoring_matrix` (named list, one `num_re × n_j` matrix per cluster; score = `S_j %*% (y_j − X_j β)`), and `get_fs_blocks.merMod()` builds `Kz` from the random-effects design `Z` (`lme4::getME(object, "Z")`) instead of the fixed-effects `pp$X` — fixing the `Z ≠ X` crash (e.g. `Reaction ~ Days + (1 \| Subject)` → "non-conformable arguments"). Numerically inert where the old code worked (`Z == X` ⇒ `crossprod(zj) == crossprod(xj)` exactly). Roxygen (`get_fs`/`get_fs.merMod`/`get_fs_lmer`) + `AGENTS.md` attribute listing document the new attribute for both backends; new vignette `vignettes/scoring-matrices.Rmd` (lavaan CFA + lme4, hand-reconstruction of scores, comparison table). 5 new `test_that` blocks (score identity vs `ranef()`, structure, Z≠X regression, unbalanced clusters, legacy `get_fs_lmer()`). | 2026-08-16 | PLAN 03 |
| 9 | **merMod cluster-name robustness + efficiency pass** — `get_fs_blocks.merMod()` named blocks by *first appearance* in the data (`unique(as.character(object@flist[[1]]))`) while every lme4 structure it consumes (`split()`/factor levels, the `Z` columns, the `b` random-effect vector, `ranef()` row order) follows the *canonical level order* — blocks were mislabeled whenever appearance order ≠ level order (shuffled rows, reversed factor levels, non-monotonic numeric cluster ids; values were always correct, only the labels wrong). Now names/rows come from `levels(as.factor(object@flist[[1]]))`, Z is sliced by level index `(j−1)·num_re+1:num_re`, and EB scores use `lme4::getME(object, "b")` reshaped level-major (bit-identical to `ranef()`, ~7× faster) — also dropping the unused `model.frame(object)` call. Verified: lme4 2.0.6 exposes no `pp$Z` (refclass `merPredD` has only `Zt`), so exported `getME()` is both the only and the stable choice; direct-slot access is not faster (sub-µs either way); per-cluster sparse `Zt[rws, idx]` slicing benchmarked and rejected (~55 µs/call vs 0.2 ms full densification at K=300). 3 new regression tests (shuffled rows, reversed factor levels, non-monotonic numeric ids); each confirmed to fail against the appearance-order naming. | 2026-08-16 | PLAN 03, follow-up |
| 10 | **User-supplied latent priors `prior_mean`/`prior_cov` for `get_fs()`** (lavaan-only) — new arguments on `get_fs.data.frame()`, `get_fs.lavaan()`, and the `get_fs_lavaan()` wrapper (placed before `...` so they are never forwarded to `lavaan::cfa()`); non-NULL priors are treated as **fixed external priors shared across all lavaan groups**. Semantics: NULL preserves current behavior exactly; priors may be supplied independently; `prior_mean` is a length-q vector, `prior_cov` a q×q matrix (scalar/1×1 for q = 1); named inputs validated against latent names and reordered to model order; `prior_cov` validated finite/square/symmetric/positive-definite. Restrictions: regression/EB only (Bartlett/ML error), `reliability = TRUE` errors, merMod + `get_fs_lmer()` reject via `...`. Math: `psi_override` threaded through `compute_fspars()`/`compute_a()`/`compute_evfs()`/`compute_ldfs()`/`compute_grad_ld_evfs()`/`vcov_ld_evfs()`/`correct_evfs()`, so `corrected_fsT = TRUE` and `vfsLT = TRUE` treat the supplied covariance as fixed (no prior sampling uncertainty propagated); `prepare_fs()` in `get_fs_blocks.lavaan()` uses the overridden `psi`/`alpha` for complete data, every missing-data pattern, and every group. All derived outputs (`fs`, `fsT`, `fsL`, `fsb`, `scoring_matrix`, `se_*`, `ev_*`, `ecov_*`) are recomputed from the prior-based scoring matrix — verified equivalent to manual `compute_fscore()` calls. Roxygen docs + multi-group example; 22 new `test_that` blocks in `test-get_fs_priors.R` (equivalence, attributes, reordering, all validation errors, corrected-SE/vcov, `vcov_corrected()` on a prior-adjusted `tspa()` fit, data.frame/matrix/legacy entry points, q = 1 forms, multi-group incl. identical-group invariance, missing data, merMod rejection). `fs_to_group_list()` round-trip verified on prior-adjusted unified output. Out of scope (unchanged): group-specific priors, reliability under priors, merMod priors, `augment_lav_predict()`/OpenMx-facing helpers. | 2026-08-16 | PLAN 03 (fs priors) |

## Verification state (2026-08-16 — PLAN 02, PLAN 03 [mermod], PLAN 03 [fs priors])

- `devtools::test()`: **596 pass, 0 fail, 0 warn** (435 at the PLAN 02
  close; +151 expectations from PLAN 03's 5 new `test_that` blocks in
  `test-get_fscore.R`; +10 from the 3 robustness regression tests in the
  PLAN 03 follow-up). Also re-run inside the full check (`checking tests
  ... OK`).
- `devtools::document()`: OK, idempotent (pre-existing `tspa_plot.R` link
  warning — file now in `legacy/`).
- PLAN 03 full `devtools::check()` (2026-08-16, builds vignettes):
  `checking package vignettes ... OK` — **all 14 vignettes rebuild**;
  `checking tests ... OK`; status **0 errors, 2 WARNINGs, 4 NOTEs**.
  Caveat: `devtools::check()` defaults to `cran = TRUE` ⇒ `R CMD check
  --as-cran`, which adds an as-cran-only "top-level files" NOTE (a plain
  `R CMD check` skips that step entirely). A/B-verified in a clean
  worktree at `9c60ff8`: that NOTE **already fired on the pre-existing**
  top-level files (`PERF_FIX_SUMMARY.md`, `STATUS.md`, `archive/`,
  `dependency_analysis.md`, `legacy/`) before this change; PLAN 03 only
  adds `_PLAN_03_mermod_scoring_matrix.md` to that already-triggered
  list. The other findings are all pre-existing baseline items:
  S3-consistency for `get_fs.lavaan()`/`get_fs.merMod()`, undocumented
  `fsm`/`...` Rd args (2 WARNINGs); `.lintr` hidden file, unused
  `Matrix` Import, `stats::model.frame` not imported (3 NOTEs).
   Under the plain-`R CMD check` measurement used for the recorded
   baseline, the current tree reproduces exactly the same 2W/3N.
- PLAN 03 follow-up full `devtools::check()` (2026-08-16, after the
  robustness/efficiency fixes; as-cran default): **0 errors, 2 WARNINGs,
  3 NOTEs** — tests OK, all 14 vignettes rebuild. Strictly better than the
  PLAN 03 baseline (2W/4N): removing the unused `model.frame(object)` call
  in `get_fs_blocks.merMod()` resolved the pre-existing
  "`stats::model.frame` not imported" NOTE. The remaining findings are all
  pre-existing baseline items (S3 consistency for `get_fs.lavaan()`/
  `get_fs.merMod()`, undocumented `fsm`/`...` Rd args; `.lintr`, unused
  `Matrix` Import, top-level files incl. the PLAN 03 plan file).
- PLAN 03 (fs priors) verification (2026-08-16): `devtools::test()`
  **658 pass, 0 fail, 0 warn, 0 skip** (596 at the PLAN 03 [mermod] close;
  +62 expectations from the 22 new `test_that` blocks in
  `test-get_fs_priors.R`, A/B-verified by running the suite with the new file
  removed). `devtools::document()` OK and idempotent (regenerates
  `man/get_fs.Rd`, `man/get_fs_lavaan.Rd`; re-run changes nothing). Full
  `devtools::check()` (as-cran default; `checking examples ... OK` incl. the
  new multi-group prior example; all 14 vignettes rebuild; tests OK under
  check): **0 errors, 2 WARNINGs, 3 NOTEs** — byte-for-byte the pre-existing
  baseline items (S3-consistency for `get_fs.lavaan()`/`get_fs.merMod()`;
  undocumented `fsm`/`...` Rd args; `.lintr` hidden file, unused `Matrix`
   Import, top-level files). No new check findings introduced.
- PLAN 04 (tspa partable / schema renderer) verification (2026-08-17):
  `devtools::test()` **763 pass, 0 fail, 0 warn, 0 skip** (658 at the PLAN 04
  start; +40 expectations from the new `R/lavaan_compat.R` compat module and
  its golden canary `tests/testthat/test-lavaan_compat.R`; +65 from
  `tests/testthat/test-tspa_render.R` — pinned-format renderer tests, schema
  construction, `tspa()` contract/attributes, product-score auto-alias incl.
  ambiguity errors). Phase-2 A/B gate: schema renderer output
  character-for-character identical to the legacy string-append builders in
  10/10 canonical cases (SF SG/MG/3-predictor/verbatim user models incl.
  comments and trailing newlines; MF SG 2- and 3-factor, list-attr and plain
  matrix `fsT`/`fsL`, growth with intercepts); product-score auto-alias
  bit-identical to the old manual-rename workaround (model, coef, vcov,
  standardizedSolution; both vignette data paths). `DESCRIPTION` and
  `NAMESPACE` unchanged by design (no new exports; lavaan bound
  intentionally undeclared — drift defense is the compat module's
  dependency-contract table, the "layout not supported" error naming the
  tested-up-to version, the canary tests, and the new CI lavaan axis in
  `.github/workflows/R-CMD-check.yaml`: pinned 0.7-2 full-check job +
  lavaan-dev tripwire running only the `test-tspa*.R` /
  `test-lavaan_compat.R` / `test-get_fs_int.R` subset). Vignettes:
  `R2spa.Rmd` + `multiple-factors.Rmd` re-knit verification-only (all 5
  printed `tspaModel` blocks unchanged, confirmed in the knit artifacts);
  `get_fs_int-vignette.Rmd` now relies on the automatic product-score alias
  (manual rename workaround removed; estimates identical). The cutover
  fallback (`tspa_env$render = "string"` + verbatim legacy builders) was
  removed at plan completion; the pinned format it guaranteed is frozen in
  `test-tspa_render.R` as in-test reference builders `ref_sf()`/`ref_mf()`.
  Full `devtools::check()` in the final state (as-cran default; all 16
  vignettes rebuild; `checking tests ... OK`; vignette output re-build OK):
  **0 errors, 2 WARNINGs, 3 NOTEs** — byte-for-byte the pre-existing
   baseline items (S3-consistency for `get_fs.lavaan()`/`get_fs.merMod()`;
   undocumented `fsm`/`...` Rd args; `.lintr` hidden file, unused `Matrix`
   Import, top-level files). No new check findings introduced.
- Check-cleanup verification (2026-08-17, after PLAN 04) — **all
  pre-existing baseline findings cleared**:
  1. **S3 generic/method consistency WARNING** (the `get_fs()` one):
     first formal renamed `data` → `object` in the `get_fs()` generic
     **and all four methods** (`tools:::checkS3methods` compares the
     first positional arg symmetrically across every method, so the
     generic-only rename would have just moved the warning onto
     `.data.frame`/`.default`). Named call sites updated for identical
     behavior: `R/tspa.R` `@examples` (4), `tests/testthat/test-tspa.R`
     (2), `vignettes/R2spa.Rmd` (7), `vignettes/corrected-se.Rmd` (2),
     `vignettes/categorical-interaction.Rmd` (1), `README.Rmd`/`README.md`
     (2 each). One canonical `@param object` kept on the generic; 3
     stale `@param object` lines removed from method blocks.
  2. **Undocumented Rd args WARNING**: `@param fsm` ("Currently not
     used." — signature-only arg) added to `get_fs.merMod`; `@param ...`
     ("Additional arguments, passed on to `get_fs()`") added to
     `get_fs_lmer`.
  3. **Unused `Matrix` Import NOTE**: `Matrix` moved `Imports` →
     `Suggests` in `DESCRIPTION` (still exercised by
     `vignettes/correction-error.Rmd`, so no "suggested but not used"
     flip).
  4. **`.lintr` + top-level-files NOTEs**: `.Rbuildignore` now excludes
     the dev-only top-level entries (`.lintr`, `archive/`, `legacy/`,
     `PERF_FIX_SUMMARY.md`, `STATUS.md`, `dependency_analysis.md`,
     `^_PLAN_.*\.md$` plan files).
  `devtools::document()` regenerated only `man/get_fs.Rd`,
  `man/get_fs_lavaan.Rd`, `man/get_fs_lmer.Rd`, `man/tspa.Rd`;
  `NAMESPACE` byte-identical. `devtools::test()`: **763 pass, 0 fail,
  0 warn, 0 skip** (count unchanged). Final full `devtools::check()`
  (as-cran default; all 16 vignettes rebuild; `checking tests ... OK`;
  vignette re-build OK; run 2026-08-17): **0 errors, 0 warnings,
  0 notes** — first fully clean check in the recorded history. The only
  intermediate finding was the untracked top-level `_PLAN_QUARANTINE.md`
  (future quarantine plan doc, now excluded by the
  `^_PLAN_.*\.md$` build-exclusion).
- Vignette build history: exactly **3 of 13 failed** on the pre-PLAN 02
  Step-1 tree (`corrected-se.Rmd`, `multilevel.rmd`, `tspa-vignette-mx.Rmd`
  — the "7/8 of 13" reports were wrong); **13/13 build on 2026-08-16**
  after PLAN 02 Steps 2–3; **14/14 on 2026-08-16** after PLAN 03's
  `scoring-matrices.Rmd`. Caveat: `devtools::check(args =
  "--ignore-vignettes")` still knits vignettes in the *build* phase; scope
  with `R CMD build --no-build-vignettes .` + `R CMD check
  --ignore-vignettes --no-manual <tarball>`.

## Notes

- PLAN 01/02 code is committed (258b673..5f72883; plan files archived in
  9c60ff8). PLAN 03 (merMod `scoring_matrix`) code is committed (00bf670);
  its plan file is archived as `archive/PLAN_03_mermod_scoring_matrix.md`.
- PLAN 03 (fs priors) is committed (`32fc817`); its plan file is archived as
  `archive/PLAN_03_fs_priors.md`.
- Working-tree changes for **PLAN 04 (tspa partable) + check cleanup**
  are **uncommitted** as of 2026-08-17: new `R/lavaan_compat.R` (the only
  file that reads lavaan internals), `tests/testthat/test-lavaan_compat.R`,
  `tests/testthat/test-tspa_render.R`, `archive/PLAN_04_tspa_partable.md`
  (archived plan), and the untracked future-work plan
  `_PLAN_QUARANTINE.md` (quarantine the `get_fs()`/`tspa()`-consuming code
  into `.quarantine/` while those two are revised; phases not started;
  excluded from the package build via `^_PLAN_.*\.md$`). Modified:
  `R/tspa.R` (stage-2 schema + `tspa_render()` + product-score auto-alias),
  `R/get_fscore.R` + `R/get_fs_methods.R` (S3 first-arg `data` → `object`
  rename + `@param fsm`/`...` docs), `R/grandStandardizedSolution.R` +
  `R/tspa_corrected_se.R` (migrated onto the compat module),
  `DESCRIPTION` (`Matrix` moved `Imports` → `Suggests`),
  `.Rbuildignore` (top-level dev-file exclusions),
  `.github/workflows/R-CMD-check.yaml` (lavaan axis),
  `tests/testthat/test-tspa.R` + 4 vignettes + `README.Rmd`/`README.md`
  (named-arg call sites / auto-alias), and the regenerated
  `man/get_fs.Rd`, `man/get_fs_lavaan.Rd`, `man/get_fs_lmer.Rd`,
  `man/tspa.Rd`. `NAMESPACE` unchanged. Commit scope must be decided
  explicitly (all listed files are PLAN 04 + check-cleanup work).
- Suggested order (all plans complete): 4 (user-facing bug) →
  5 (perf) → F1 (future).
