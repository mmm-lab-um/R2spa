# R2spa — Open Issues & Follow-ups

Tracks remaining work from completed plans (see `archive/`). Update status as
items are resolved; move finished items to the **Closed** section with the
date and commit/PR reference.

**Last updated:** 2026-08-22 — **multi-group mirt `get_fs()` (31b3809 / PLAN 10 / F7)**: `get_fs.MultipleGroupClass()` implemented — whole-fit per-observation scores + trailing `group` column + a per-group (`list`) `psi` attr, reusing the single-group per-row `mirt_per_obs` engine; fully-missing rows reconciled via `completely_missing`; new `test-get_fs_mirt_multigroup.R`. Suite **3358 pass** / 0 fail, `R CMD check` 0/0/0. **mirt `get_fs()` ψ fix (fc16f81)**: `get_fs.SingleGroupClass()` now uses the model's **full** estimated factor covariance (`mirt_full_cov()`, between-factor covariances included) as `psi`, not the unit-variance quadrature prior — fixing a ≥2-factor latent-covariance bug that only surfaced off the diagonal (1-D invariant). **P0/P1/P2 vignette review** (d5fac0d / 4d44a56 / 8867f6f): motivation + Crossref-verified references; `get_fs()` column glossary + 2S-PA framing (P1); vignette cross-link hub via relative `.html` links + `\VignetteIndexEntry`/title normalization (P2 — `?vignette=` does not resolve under rmarkdown 2.31, so the hub uses plain sibling links); correction-error scope note. New plan `archive/PLAN_10_multigroup_mirt.md` (follow-up **F7**: multi-group mirt `get_fs()`). Suite **3298 pass** / 0 fail. Prior: 2026-08-20 — **PLAN 09**: `tspa()` per-unit `fsL`/`fsT` pooling — new `reduce = c("mean", "median")` collapses FIML per-pattern and merMod per-cluster `fsT`/`fsL`/`fsb` to one representative value per group (replaces PLAN 06 §4's "not supported" `stop()`; enables missing-data + merMod 2S-PA). Prior: 2026-08-19 two features landed on `refactor/core`: **PLAN 07** `b9ca3c9` — new exported `fs_indiv()` per-row API + `psi`/`alpha` latent-moment attributes on `get_fs()` + `augment_lav_predict()` refactored onto the shared engine — resolves F4, closed #14; then the **`mirt` `SingleGroupClass` `get_fs()` method** (follow-up F2) + per-row `mirt` branch in `fs_indiv()` — closed #15; plus the `correct_evfs()` performance fix closing the last tracked Open item #5. Suite **3176 pass** / 0 fail, check 0/0/1 OpenMx NOTE. Prior context: WI-1/WI-2 `10b5b0d`/`966ff9b` (2339), PLAN 06 `7c51d23`, `method="ML"`/`"mean"` + merMod speedup `aa5be3e`/`b730b53`)

## Open

_No tracked open issues as of 2026-08-19. The sole open item, `#5`
(`correct_evfs()` performance, P2), was closed that day — see Closed. Remaining
work is future-scope (Follow-ups F1/F5/F6) or the quarantined stage-2
re-integration (`vcov_corrected()`, `tspa_mx()`, `get_fs_int()`, grand
standardization — see `archive/PLAN_QUARANTINE.md`).*

## Follow-ups (future work, explicitly out of scope for their plans)

| # | Item | Source | Status |
|---|------|--------|--------|
| F1 | Wire `tspa()` to accept unified/list factor-score shapes natively (via `fs_to_group_list()` or dual support) — currently relies on attributes; the `format` default change shifts multigroup output to list-valued attributes that `tspa()` must keep parsing. | PLAN 01, remaining issue 7 | deferred |
| F2 | `mirt` S3 method (`get_fs.SingleGroupClass`) — **resolved 2026-08-19**: per-row regression-form `fsL = I − Vpost·Ψ⁻¹` (univariate = `1 − SE²`); `mirt` stays a `Suggests` dependency; `get_fs.MultipleGroupClass` errors. See #15. | PLAN 01, remaining issue 8 | resolved 2026-08-19 (#15) |
| F3 | **Resolved 2026-08-20 via the pooling variant (#16)** — `tspa(reduce = "mean"/"median")` converts per-pattern (FIML) / per-cluster (merMod) `fsL`/`fsT`/`fsb` to long-form individual-specific values, then reduces them to one representative value per group feeding the usual schema (small patterns are pooled, not fit as tiny sub-groups). The **sub-group** variant below remains deferred. | *sub-group variant (deferred):* `tspa()` per-pattern stage 2 via lavaan sub-groups `(group × pattern)` from a synthesized pattern column off `fs_pattern` (extend `tspa_schema_mf()` group indices to (group, pattern) pairs; user `group.equal` still applies). Hazard: lavaan's level ordering for `group = c(a, b)` must be derived + pinned canary-style (cf. `R/lavaan_compat.R`); tiny sub-groups (n ≲ 5) may not converge. Design: PLAN 06 §7a. | PLAN 06, §7a | pooling resolved 2026-08-20 (#16); sub-group deferred |
| F4 | Per-row column API (user-requested): user-level function appending individual-specific `fsL`/`fsT`/`fsb` values as columns to the factor-score data frame, keyed by `fs_pattern`'s per-row labels (e.g. `ld_fs_visual_x1`, `ev_fs_visual`) so each case carries the matrices matching its own pattern. **Resolved 2026-08-19** by the exported `fs_indiv()` (see #14). | PLAN 06, §7b | resolved 2026-08-19 (#14) |
| F5 | mirt-EAP per-row SE variant: the per-row `mirt` SE is currently mirt's `fscores(..., full.scores.SE = TRUE)`; also expose the alternative `sqrt(diag(Vpost_i))` (regression-SE convention) as an option/flag so both conventions are available. | plan `archive/PLAN_08_mirt_fs.md`, follow-up | deferred |
| F6 | mirt `psi` = the estimated factor covariance. The **single-group** path is handled: `get_fs.SingleGroupClass()` uses `mirt_full_cov()` (reads `COV_ij` / factor variances) as `psi` — **resolved 2026-08-22 (fc16f81)** (supersedes the `psi = diag(q)` wording in #15 / PLAN 08 for the SG path). The **multi-group** path uses the **group-specific** covariance (per-group `mirt_group_pars()`) — **resolved 2026-08-22 (31b3809, F7)**; both SG and MG paths are done. | plan `archive/PLAN_08_mirt_fs.md` | SG resolved (fc16f81); MG resolved (31b3809) |
| F7 | ~~multi-group mirt `get_fs()` (`MultipleGroupClass`)~~ — **resolved 2026-08-22 (31b3809)**: `get_fs.MultipleGroupClass()` implemented (see Closed #18) — per-group `psi` via `mirt_full_cov(mirt::extract.group(object, k))` (new `mirt_group_pars()`), whole-fit scores/SEs/`return.acov`, trailing `group` column + per-group (`list`) `psi` attr, reusing the SG per-row `mirt_per_obs` list attrs. | plan `archive/PLAN_10_multigroup_mirt.md` | **done 2026-08-22 (31b3809)** |

## Closed

| # | Issue | Closed | Reference |
|---|-------|--------|-----------|
| 18 | **Multi-group mirt `get_fs()` (`MultipleGroupClass`) (PLAN 10)** — `get_fs.MultipleGroupClass()` replaces the prior "not supported" stub. Per-observation scores/SEs + posterior covariance come from `mirt::fscores()` on the whole fit; each observation's regression-form `(fsL, fsT, fsb)` is built by the shared `compute_lav_fs_matrices()` using its own group's factor covariance (new co-located `mirt_group_pars()`, reusing `mirt_full_cov()` on `mirt::extract.group(object, k)`). Output = the single-group per-row column set + a trailing `group` column (the model's `groupNames`; `NA` for completely-missing rows, matching the all-NA row convention). `attr("psi")` is a named list (one `q x q` per group) — the one structural difference from the single-group result. Completely-missing rows are reconciled against the full row set via `mirt::extract.mirt(object, "completely_missing")` (mirt drops them from every extraction; `N`/`group`/`rowID` are scorable-length). `fs_indiv()` dispatches on the unchanged `mirt_per_obs` marker (per-row `fsL`/`fsT`/`fsb`), dropping the `group` column. No `NAMESPACE`/`Imports` change. Docs: `get_fs()` description + `prior_mean`/`format` params; NEWS entry. New `test-get_fs_mirt_multigroup.R` (9 units); the single-group dispatch test no longer asserts the removed guard. Suite **3358 pass** / 0 fail; `R CMD check` 0 errors / 0 warnings / 0 notes. | 2026-08-22 | `31b3809` |
| 17 | **mirt `get_fs()` ψ fix — use the estimated factor covariance** — `get_fs.SingleGroupClass()` set `psi = diag(q)` (mirt's **unit-variance quadrature** prior, not an estimate); with ≥2 correlated factors the shrinkage loadings / error covariances (`fsL = I − Vpost·Ψ⁻¹`, `fsT = fsL·Vpost`) silently dropped the between-factor covariances (1-D invariant: mirt fixes `COV_11 = 1`). New co-located `mirt_full_cov()` (R/get_fs_methods.R) reads the model's estimated `COV_ij` from `coef(fit)$GroupPars` and supplies it as `psi`. Verified: `get_fs()` now matches the previously-correct hand-rolled full-covariance score matrices exactly (1-D and 2-D); the lavaan path already used `est$psi` (unchanged); new correlated-2-factor mirt regression test in `test-get_fs_mirt.R`. Roxygen/comments corrected. Suite **3298 pass** / 0 fail; check 0 errors / 0 warnings. | 2026-08-22 | `fc16f81` |
| 16 | **`tspa()` per-unit `fsL`/`fsT` pooling (PLAN 09)** — new `tspa(reduce = c("mean", "median"), ...)` collapses per-unit factor-score attributes to one representative value per group before the usual stage-2 schema. Detection `is_per_unit_fs()` (3-D arrays for merMod; a per-group list-of-matrices for multiple missing-data patterns), pooling `pool_per_unit()` (reuses `resolve_fs_per_row()` + `fs_row_cols()` from `R/fs_indiv.R` for long-form per-row values; reduces each column by `mean`/`median` with `na.rm`; reassembles one `fsT`/`fsL`/`fsb` per group or single), and single-factor `pool_se_fs()`. Replaces PLAN 06 §4's per-pattern `stop()` — FIML and merMod factor scores now fit (merMod `tspa()`, previously broken at `upper.tri`, now works). `mean` keeps the pooled `fsT` positive semi-definite (convex); `median` is opt-in with a PSD guard warning. Pooled values + `pooled_fs` marker attached to the fit; complete-data paths are a no-op (byte-identical, identity-tested). No new `Imports`; no `NAMESPACE` change (internal helpers) + `@param reduce`/`fsT`/`fsL`/`fsb`/`@details`/`@examples`. Out of scope: mirt per-obs (not pooled), the sub-group variant, cluster-size weighting. 12 new `test_that` (66 expectations) in `test-tspa_pooled.R` + 3 stale PLAN 06 units in `test-tspa.R` rewritten. Suite **3247 pass** / 0 fail, 0 skip; check 0/0/1 OpenMx NOTE. | 2026-08-20 | this commit + `archive/PLAN_09_per_unit_pooled_fs.md` |
| 5 | **`correct_evfs()` dominated `corrected_fsT = TRUE` runtime** — the per-row `lav_func_jacobian_complex()` loop did `q × (npar+1) × 2 + 1` `lavInspect()` file-reads per call (57 traced on a 13-param 2-factor fit). Two numerically-**exact** optimizations in `R/get_fscore_math.R`: (a) one `lav_func_jacobian_complex()` over the **full** `p × c` `a` matrix instead of one call per row — row `i`'s Jacobian is the column-major slice `J[i + p*(0:(c-1)), ]`, entry-identical to the per-row call (complex-step perturbations are parameter-local and linear in the output entries), giving `q×` fewer `a`-matrix evaluations; and (b) `compute_fspars()`/`compute_a()` now accept pre-fetched `frees`/`mats`, hoisting the `lavInspect("free")`/`("est")` file-reads out of the per-evaluation loop (`psi_override` semantics unchanged). Verified **bit-identical** (`identical()`, `all.equal` tolerance 0) across 1/2/3-factor SG, prior-adjusted, and 2-group MG fits under both `regression` and `Bartlett`. `corrected_fsT = TRUE` runtime (regression, median) — 1-factor 0.011→0.005 s, 2-factor 0.033→0.007 s, 3-factor **0.076→0.010 s (~8×)**, 2-group **0.198→0.029 s (~7×)**; the win scales with factor count. No `NAMESPACE`/roxygen change (internal, undocumented helpers). | 2026-08-19 | `#5` (this commit) |
| 15 | **`mirt` Item-Response `get_fs()` support (`SingleGroupClass`)** — new `get_fs.SingleGroupClass()` (S4 mirt fit, univariate & multi-factor) + a `get_fs.MultipleGroupClass()` guard that errors with a "fit/extract the single group first" message. Scores are the EAP posterior means (`fs_<factor>`, named by the mirt factor names via `mirt::extract.mirt(m, "factorNames")`); per-observation `fsL`/`fsT` are **per-row lists** of `q × q` regression-form matrices from `compute_lav_fs_matrices(Vpost_i, diag(q), 0, "regression")` where `Vpost_i` is mirt's `fscores(..., return.acov = TRUE)` posterior covariance (`diag(Vpost_i) = SE²`). Unidimensional: `F1_by_fs_F1 = 1 − SE²`, `ev = (1 − SE²)·SE²` (machine precision); multi-factor off-diagonals give the shrinkage loadings / error covariances. `fs_pattern = list(label = seq_len(n), pat = NULL)`, marker `mirt_per_obs = TRUE`, `psi = diag(q)`, `alpha = 0` (mirt's default unit/zero prior; `prior_*` unsupported here). `fs_indiv()` gained a **first-dispatch** per-row branch (`resolve_per_obs()`) that mints one block per row — the `fs_indiv()` body is otherwise unchanged, so a mirt output round-trips through it like lavaan/merMod. Missing data: no-scorable-indicator rows → all-NA block (NA score/SE/ev/intercept), matching lavaan. `mirt` stays in `DESCRIPTION: Suggests` (namespaced `mirt::` calls, `require_mirt()` guard). New `tests/testthat/test-get_fs_mirt.R` (72 expectations). Resolves follow-up F2. | 2026-08-19 | this commit + `archive/PLAN_08_mirt_fs.md` |
| 14 | **`fs_indiv()` per-row API + effective latent `psi`/`alpha` + `augment_lav_predict()` reconcile (PLAN 07)** — new exported `fs_indiv(fs)` re-derives, per row, the individual-specific values (`_se`, `<lvs>_by_<lv>_*`, `ev_*`/`ecov_*`, and the per-pattern `fsb` intercepts) from the pattern's `fsL`/`fsT`/`fsb` via a shared **value-only** engine `fs_row_cols()` that `augment_lav_predict()` now also uses (SEs are therefore always pattern-consistent); returns a data frame with exactly the factor-structure columns (score columns kept, the other `get_fs()` extras dropped) + `group_col`/`id_vals`/`block_label` attributes, under a one-block-per-group contract (lavaan: complete-data `fsT`; all-missing patterns → all-NA block). Adds effective latent-moment attributes to **every** get_fs() backend: `psi` = `prior_cov` if supplied else the lavaan/RE-term covariance, `alpha` = `prior_mean` if supplied else the lavaan mean (merMod: named zero vector) — point estimates only, carried through `fs_to_group_list()`. Refactors `augment_lav_predict()` onto `augment_fs2()` (drops the ~54-line duplicate `fsL`/`fsT`/`fsb` derivation). Resolves follow-up F4; documents F5 (mirt-EAP SE variant, out of scope). | 2026-08-19 | `b9ca3c9` + `archive/PLAN_07_fs_indiv_and_latent_moments.md` |
| 13 | **`get_D()` lme4-2.x RE-covariance convention (merMod EB)** — `get_D(object)` is now self-contained: splits `@theta` by `@cnms` block lengths (the same idiom lme4's `mkVarCorr` uses; no dependence on the `"clen"` attribute, absent on `@theta` since lme4 2.x even for single-term fits) and returns the SCALED first RE-term covariance `VarCorr(x)[[1]]/sigma(x)^2`; the lme4 2.0.6 convention is documented in-source (R/get_fs_methods.R ~:601). Baseline (step 1) showed the `b730b53` clen restore had already made the value exact under 2.x — the plan's "missing `sigma^2` factor" premise was a misread (the motivating probe gap of 1378.179 vs 1.435 was exactly `sigma^2`; the EB formulas carry the explicit scale) — so the change is numerically **bit-identical** (RDS-verified, all 8 EB/ML outputs on 4 fixtures; no user-visible value change, no stage-2 SE shift, no vignette re-knit). New `tests/testthat/test-lme4_compat.R` canary (per-term `s^2 * tcrossprod(L) == VarCorr` at 1e-15, `get_D == VarCorr[[1]]/s^2`, warning-free multi-term parse) guards future lme4 theta drift — needed because scores are D-free via `getME("b")`, so ranef-identity tests alone would stay green through a convention change. Term-1-only EB `fsT`/`fsL`/`scoring_matrix` reference tests for the 2+1/2+2 fixtures at 1e-12 (joint-vs-term-1 posterior gap recorded in-test, no multi-term score-identity pin by design); roxygen accuracy fixes (EB scores documented as first-term `ranef()`; `corrected_fsT`/`format` declared-but-ignored notes for merMod). | 2026-08-18 | `966ff9b` + `archive/get_d-mermod-lme4-2x.md` |
| 4 | **`reliability = TRUE` + single-group multi-factor + default `unified` format → hard error** — the `get_fs.lavaan()` reliability guard now derives dimensionality from the model (`nrow(est$psi)` for SG, `nrow(est[[1]]$psi)` for MG — format-agnostic; the PLAN 06 missing-data stop for reliability still fires upstream, so per-pattern fits cannot misfire) instead of the `fsb` attribute-shape test that never fired for SG in the unified format. SG multi-factor now warns ("Computation of reliability for a multi-factor model is not currently supported.") and omits the attribute in BOTH formats; all other cases verified unchanged (SG unidimensional pin .9607411, MG unidimensional per-group + overall values/names, non-standardized warning, priors/mean rejections). 4 new tests in `test-get_fscore.R` (SG 2-factor unified + list — the unified case is the pre-fix-failing regression net — MG 2-factor both formats, MG 1-factor guard); `@param reliability` wording updated (roxygen only). | 2026-08-18 | `10b5b0d` (PLAN 01, remaining issue 4) |
| 12 | **Per-pattern factor-score attributes for lavaan missing-data fits** — `assemble_fs_blocks()` now keeps one `fsT`/`fsL`/`fsb`/`scoring_matrix` value per observed-indicator pattern: a k=1 group keeps a plain matrix/vector (complete-data values/shapes unchanged, regression-tested), and a group with k≥2 patterns gets a named list keyed by pattern label (observed indicators joined with `"+"` in indicator order). New per-group `fs_pattern` attribute = `list(label, pat)` (per-case pattern label — `NA` for cases with all indicators missing — plus a named logical p×k indicator-by-pattern matrix) so a future API (follow-up F4) can index per-row columns. `prepare_fs()` carries `pat_label`/`pat` through the blocks; the "blocks have differing fsT/fsL/fsb attributes" message and `check_blocks_identical()` deleted (nothing is dropped anymore). `fs_to_group_list()` round-trips the nested shapes and `fs_pattern` in both directions. SE paths (`corrected_fsT`/`reliability`/`vfsLT`) now `stop()` explicitly on multi-pattern data (previously a cryptic dimension error deep in `compute_fspars()`/`correct_evfs()`); `tspa()` rejects nested per-pattern `fsT`/`fsL`/`fsb` attributes with an explicit error (both the single-group k>1 "misread as k groups" trap and multigroup covered; per-pattern stage 2 is follow-up F3). Tests: per-pattern contract + `fs_pattern` content in `test-assemble_fs_blocks.R`; new `test-get_fs_missing.R` (SG/MG, pattern matrices matched against `lavPredict(acov=TRUE)` **by pattern label**, `fs_pattern` vs raw NA positions, `format="list"`, complete-data regression guard); SE-path errors in `test-get_fscore.R`; guard tests in `test-tspa.R` (incl. real `get_fs()` missing-data output); `fs_pattern` round-trips in `test-fs_converters.R`. Implementation note: a pattern's `fsT`/`fsL` equal lavaan's raw `acov` entry only via the package's canonical mapping `compute_lav_fs_matrices()` (`fsT = (I − AΨ⁻¹)A`, pinned in `test-lavPredict_equivalence.R`) — the tests assert that reference, not raw acov equality. Quarantined consumers (`vcov_corrected`, `tspa_mx_model`, `grandStandardizedSolution`) read `fsT`/`fsL` as single matrices and need adapting at re-integration. Also: `^\.git$` added to `.Rbuildignore` (worktree `.git` *file* leaked into built packages, adding a hidden-file NOTE). | 2026-08-18 | `archive/PLAN_06_per_pattern_fs_attrs.md` (implementation `bb64d2e`, merged into `refactor/core` in `7c51d23`) |
| 11 | **Quarantine of `get_fs()`/`tspa()`-consuming code** — while those two contracts are being revised, every in-package consumer moved to `.quarantine/{R,tests,vignettes}/` (excluded from the build via `^\.quarantine$` in `.Rbuildignore`): `R/get_fs_int.R`, `R/tspa_mx.R`, `R/tspa_corrected_se.R`, `R/grandStandardizedSolution.R`; full test files `test-get_fs_int.R`, `test-grandStandardizedSolution.R`; 8 vignettes + 3 RDS fixtures. Embedded blocks extracted into 2 new self-contained quarantined test files (`test-tspa_mx.R`: the Mx comparison block + umx/OpenMx missing-data block with copied setups; `test-vcov_corrected.R`: the MG `vcov_corrected()` + prior-adjusted `vcov_corrected()` tests with copied setups) and appended to the 2 quarantined files (product-score auto-alias section → `test-get_fs_int.R`; grandSS wrapper A/B → `test-grandStandardizedSolution.R`). Kept in-package by decision: `R/lavaan_compat.R` (now consumed only by its own canary tests — its only package consumers were the two quarantined files) and the 6 core vignettes. `NAMESPACE` shrinks to 8 exports + 4 `get_fs` S3 methods (no `OpenMx`/`utils`; `stats` = `setNames`; `lavaan::vcov` re-declared on `get_fs()` for the bare `vcov()` calls in `get_fscore_math.R`). Roxygen links to quarantined topics reworded. `tspa()` product-score auto-alias (`tspa_sf_alias`) retained (no `get_fs_int` dependency) with the ambiguous-candidates core test kept in `test-tspa_render.R`. Quarantined tests are self-contained (setups copied) with provenance headers for re-integration (`git mv` back → `document()` → `test()` → `check()`). `OpenMx` kept in `DESCRIPTION: Imports` until re-integration. | 2026-08-17 | `archive/PLAN_QUARANTINE.md` |
| 1 | **merMod column-name regression** — restored pre-refactor `u0_eb`-style *column names* in the merMod path via `legacy_names` switch; `get_fs_lmer()` defaults `legacy_names = TRUE` so `vignettes/multilevel.rmd` (`tspa_mx_model`) works unchanged; default `get_fs()`/`get_fs.merMod()` now use `fs_u0`-style names. Legacy output is name-compatible (not byte-identical) with the pre-refactor result — extra `_se` columns, `fsL`/`fsT` attributes, NULL row names (delta documented on `get_fs_lmer()`/`get_fs.merMod()`). Fixed the related overwritten-`fsT`-rownames + duplicate-`re_names` bugs (item 7 sub-bullets). | 2026-08-15 | PLAN 02, Step 2 |
| 2 | **Vignette breakage on `format = "unified"`** (verified failure set 3/13: `corrected-se.Rmd`, `multilevel.rmd`, `tspa-vignette-mx.Rmd`) — `multilevel.rmd` fixed by item 1; `R/tspa.R` now validates `fsT`/`fsL` group-count consistency (plain matrix = 1 group, so single-group length-1 list attributes may be mixed with plain matrices, e.g. Bartlett identity `fsL`) with a clear mismatch error, and `tspa_mf()` accepts all single-group shape combinations; `tspa-vignette-mx.Rmd` uses `format = "list"` for its direct attribute arithmetic. 6 new regression tests. **13/13 vignettes build** (verified 2026-08-16). | 2026-08-16 | PLAN 02, Step 3 |
| 3 | **`get_fs(data = matrix)` dispatch** — already resolved in the working tree: `get_fs.data.frame()` re-dispatches matrices (`is.matrix(data)` → `as.data.frame(data)` → `get_fs()`). Live check: `get_fs(as.matrix(...))` is `identical()` to `get_fs(data.frame)` incl. attributes. Locked in with a regression test in `test-get_fscore.R`. | 2026-08-16 | PLAN 02, Step 4 |
| 6 | **Missed cheap-lookup `lavInspect` sites** — all 3 replaced with direct slot access: `tspa_corrected_se.R` (`@Data@ngroups`), `get_fs_methods.R` (`unlist(@Data@norig)`), `grandStandardizedSolution.R` (`unlist(@Data@nobs)`). Slot values proven equal to `lavInspect()` output (SG + MG); MG standardization path A/B'd by existing tests; MG reliability `norig` path A/B'd live. One further cheap lookup noted out of scope: `get_fscore_math.R:128` (`"ngroups"` in `augment_lav_predict()`). | 2026-08-16 | PLAN 02, Step 5 |
| 7 | **Dead-code cleanup** — remaining pieces removed: commented-out `# fs_se[is.nan(fs_se)] <- 0` line in `augment_fs()` (`get_fscore.R`) and the unreachable `is.data.frame()` guard in `get_fs.data.frame()` (`get_fs_methods.R`); the merMod `fsT_j` rownames + duplicate `re_names` sub-bullets were fixed in Step 2. | 2026-08-16 | PLAN 02, Step 6 |
| 8 | **merMod `scoring_matrix` + Z-design fix + vignette** — `get_fs.merMod()` now emits a per-cluster `scoring_matrix` (named list, one `num_re × n_j` matrix per cluster; score = `S_j %*% (y_j − X_j β)`), and `get_fs_blocks.merMod()` builds `Kz` from the random-effects design `Z` (`lme4::getME(object, "Z")`) instead of the fixed-effects `pp$X` — fixing the `Z ≠ X` crash (e.g. `Reaction ~ Days + (1 \| Subject)` → "non-conformable arguments"). Numerically inert where the old code worked (`Z == X` ⇒ `crossprod(zj) == crossprod(xj)` exactly). Roxygen (`get_fs`/`get_fs.merMod`/`get_fs_lmer`) + `AGENTS.md` attribute listing document the new attribute for both backends; new vignette `vignettes/scoring-matrices.Rmd` (lavaan CFA + lme4, hand-reconstruction of scores, comparison table). 5 new `test_that` blocks (score identity vs `ranef()`, structure, Z≠X regression, unbalanced clusters, legacy `get_fs_lmer()`). | 2026-08-16 | PLAN 03 |
| 9 | **merMod cluster-name robustness + efficiency pass** — `get_fs_blocks.merMod()` named blocks by *first appearance* in the data (`unique(as.character(object@flist[[1]]))`) while every lme4 structure it consumes (`split()`/factor levels, the `Z` columns, the `b` random-effect vector, `ranef()` row order) follows the *canonical level order* — blocks were mislabeled whenever appearance order ≠ level order (shuffled rows, reversed factor levels, non-monotonic numeric cluster ids; values were always correct, only the labels wrong). Now names/rows come from `levels(as.factor(object@flist[[1]]))`, Z is sliced by level index `(j−1)·num_re+1:num_re`, and EB scores use `lme4::getME(object, "b")` reshaped level-major (bit-identical to `ranef()`, ~7× faster) — also dropping the unused `model.frame(object)` call. Verified: lme4 2.0.6 exposes no `pp$Z` (refclass `merPredD` has only `Zt`), so exported `getME()` is both the only and the stable choice; direct-slot access is not faster (sub-µs either way); per-cluster sparse `Zt[rws, idx]` slicing benchmarked and rejected (~55 µs/call vs 0.2 ms full densification at K=300). 3 new regression tests (shuffled rows, reversed factor levels, non-monotonic numeric ids); each confirmed to fail against the appearance-order naming. | 2026-08-16 | PLAN 03, follow-up |
| 10 | **User-supplied latent priors `prior_mean`/`prior_cov` for `get_fs()`** (lavaan-only) — new arguments on `get_fs.data.frame()`, `get_fs.lavaan()`, and the `get_fs_lavaan()` wrapper (placed before `...` so they are never forwarded to `lavaan::cfa()`); non-NULL priors are treated as **fixed external priors shared across all lavaan groups**. Semantics: NULL preserves current behavior exactly; priors may be supplied independently; `prior_mean` is a length-q vector, `prior_cov` a q×q matrix (scalar/1×1 for q = 1); named inputs validated against latent names and reordered to model order; `prior_cov` validated finite/square/symmetric/positive-definite. Restrictions: regression/EB only (Bartlett/ML error), `reliability = TRUE` errors, merMod + `get_fs_lmer()` reject via `...`. Math: `psi_override` threaded through `compute_fspars()`/`compute_a()`/`compute_evfs()`/`compute_ldfs()`/`compute_grad_ld_evfs()`/`vcov_ld_evfs()`/`correct_evfs()`, so `corrected_fsT = TRUE` and `vfsLT = TRUE` treat the supplied covariance as fixed (no prior sampling uncertainty propagated); `prepare_fs()` in `get_fs_blocks.lavaan()` uses the overridden `psi`/`alpha` for complete data, every missing-data pattern, and every group. All derived outputs (`fs`, `fsT`, `fsL`, `fsb`, `scoring_matrix`, `se_*`, `ev_*`, `ecov_*`) are recomputed from the prior-based scoring matrix — verified equivalent to manual `compute_fscore()` calls. Roxygen docs + multi-group example; 22 new `test_that` blocks in `test-get_fs_priors.R` (equivalence, attributes, reordering, all validation errors, corrected-SE/vcov, `vcov_corrected()` on a prior-adjusted `tspa()` fit, data.frame/matrix/legacy entry points, q = 1 forms, multi-group incl. identical-group invariance, missing data, merMod rejection). `fs_to_group_list()` round-trip verified on prior-adjusted unified output. Out of scope (unchanged): group-specific priors, reliability under priors, merMod priors, `augment_lav_predict()`/OpenMx-facing helpers. | 2026-08-16 | PLAN 03 (fs priors) |

## Verification state (2026-08-16/17 — PLAN 02, PLAN 03 [mermod], PLAN 03 [fs priors], PLAN 04 + check cleanup, QUARANTINE)

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
- QUARANTINE verification (2026-08-17): `devtools::test()` **707 pass,
  0 fail, 0 warn, 0 skip** (763 pre-quarantine; the delta is exactly the
  56 expectations in the quarantined blocks: Mx comparison + umx/OpenMx
  missing-data in `test-tspa.R`/`test-get_fscore.R`, `vcov_corrected()` in
  `test-tspa_render.R`/`test-get_fs_priors.R`, product-score auto-alias in
  `test-tspa_render.R`, grandSS wrapper A/B in `test-lavaan_compat.R`).
  `devtools::document()` idempotent — re-run changes nothing; `NAMESPACE`
  and `man/` diffs are purely the removed exports (4 Rd files deleted,
  `get_fs.Rd`/`get_fs_lavaan.Rd`/`augment_lav_predict.Rd` updated for the
  reworded `vfsLT` param / de-linked `tspa_mx_model()`). Full
  `devtools::check()`: **0 errors, 0 warnings, 1 NOTE** — the sole NOTE is
  "'OpenMx' in DESCRIPTION Imports but not imported from anywhere", the
   expected direct consequence of keeping `OpenMx` in `Imports` until the
   OpenMx path is re-integrated. No new findings beyond that expected NOTE.
 - PLAN 06 (per-pattern missing-data attributes) verification (2026-08-18):
   implemented in worktree `plan06` (no shared-library install;
   `load_all()` only), then merged into `refactor/core`. Worktree
   `devtools::test()`: **886 pass, 0 fail, 0 warn, 0 skip**. Merged-tree
   `devtools::test()` (with the in-flight `method="ML"`/`"mean"` + merMod
   speedup commits `aa5be3e`/`b730b53` present): **1594 pass, 0 fail,
   0 warn, 0 skip**. The only merge conflict was a docs block in the
   `get_fs()` roxygen `@return` (both sides' attribute docs in the same
   block); resolved as the union and both Rd files regenerated via
   `devtools::document()` (`NAMESPACE` untouched). Full
   `devtools::check()` in the merged tree: **0 errors, 0 warnings,
   1 NOTE** — the sole NOTE is the expected OpenMx baseline item.
- WI-1/WI-2 verification (2026-08-18): `10b5b0d` (reliability guard, this
  table's #4) + `966ff9b` (lme4-2.x `get_D()`, Closed #13). Suites:
  1594 → **1604** after WI-1 (+4 `test_that` blocks, +10 expectations) →
  **2339 pass / 0 fail / 0 warn / 0 skip** after WI-2 (new
  `tests/testthat/test-lme4_compat.R` canary +9; +726 per-cluster term-1
  attribute expectations on the 2+1/2+2 fixtures). Both
  `devtools::document()` idempotent (`NAMESPACE` unchanged; only
  `man/get_fs.Rd` + `man/get_fs_lavaan.Rd` regenerated, from roxygen
  accuracy wording). merMod `fsT`/`fsL`/`scoring_matrix` values: **no
  change** (bit-identical RDS-verified; no stage-2 SE shift,
  `vignettes/scoring-matrices.Rmd` needs no re-knit). Full
   `devtools::check()` at each step: **0 errors, 0 warnings, 1 NOTE** (the
   expected OpenMx item).
- PLAN 07 + mirt `get_fs()` verification (2026-08-19): `devtools::test()`
  **3104 pass / 0 fail / 0 warn / 0 skip** at the PLAN 07 (#14) close (2339 →
  +765 from `test-fs_indiv.R` + `test-get_fs_latent_moments.R`), and **3176
  pass / 0 fail / 0 warn / 0 skip** after the mirt `#15` method (+72 from
  `test-get_fs_mirt.R`). Independently green at each commit (bisect-checked:
  the PLAN 07 tree alone runs its full 3104 with the mirt test removed).
  `augment_lav_predict()` A/B: 1–4-factor lavaan fits identical to the old
  inline implementation after the `augment_fs2()` reconcile. Unidimensional
  mirt: `F1_by_fs_F1 == 1 − SE²` and `ev == (1 − SE²)·SE²` to ≈1e-16,
  off-diagonals match the regression-form algebra by hand; the mirt per-row
  branch is dormant unless the `mirt_per_obs` marker is set (never on
  lavaan/merMod). `devtools::document()` run with the pinned **roxygen2 8.1.0**
  (the sandbox's global 7.3.1 would churn `RoxygenNote` + `importFrom`
  splitting) — `NAMESPACE`/`man/` carry only the intended `fs_indiv` export /
  mirt-method changes. Real bug fixed en route: `resolve_group_blocks()` no
  longer `stop()`s on all-missing-pattern rows (returns an all-NA block,
  preserving `nrow`). Full `devtools::check()` (`--no-manual`; the PDF-manual
  build errors are the environmental missing `pdflatex`, not a package issue):
  **0 errors, 0 warnings, 1 NOTE** — the sole NOTE is the expected OpenMx
  `Imports` baseline item.
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
- **PLAN 04 (tspa partable) + check cleanup is committed** —
  `7c173ab` (schema renderer + compat module + render tests + CI lavaan
  axis) and `0d782b4` (S3 `get_fs()` arg rename, Rd `fsm`/`...` docs,
  `Matrix` → `Suggests`, top-level `.Rbuildignore` exclusions); plan
  archived as `archive/PLAN_04_tspa_partable.md`.
- **QUARANTINE is committed** in `6dddf5f` (renames into
  `.quarantine/{R,tests,vignettes}/`, embedded-block deletions, the 2 new
  self-contained quarantined test files, hygiene changes: `.Rbuildignore`,
  `NAMESPACE`, `man/`, `vfsLT` reword + `@importFrom lavaan vcov`); docs
  refreshed in `4e4a805`. Plan archived as `archive/PLAN_QUARANTINE.md`.
  `OpenMx` stays in `DESCRIPTION: Imports` until the OpenMx path is
  re-integrated (the sole expected check NOTE).
- **PLAN 07 (`fs_indiv()` + latent `psi`/`alpha`) is committed** in `b9ca3c9`;
  plan archived as `archive/PLAN_07_fs_indiv_and_latent_moments.md` (resolves
  F4, closed #14).
- **The `mirt` `SingleGroupClass` `get_fs()` method is committed** (follow-up
  F2, closed #15); recorded in `archive/PLAN_08_mirt_fs.md`. `mirt` remains a
  `Suggests`-only dependency (no `DESCRIPTION` `Imports` change); per-obs
  `fs_indiv()` support rides on the `mirt_per_obs` marker. Follow-ups F5
  (mirt-EAP SE variant) + F6 (2-factor `psi` unit-by-convention) documented.
- Suggested order (all plans complete): F1 / F3 / F5 / F6 (future; #5 perf, F2,
  and F4 resolved 2026-08-19).
- 2026-08-18 work-item plans archived under `archive/` (previously tracked only
  via gitignored `.opencode/plans/` paths): `get_d-mermod-lme4-2x.md` (lme4
  `get_D()`, now #13 / `966ff9b`), `get_fs-lavaan-mean.md` (`method = "mean"`,
  lavaan), `get_fs-mermod-ml.md` (`method = "ML"`, merMod; the two landed
  together in `aa5be3e`, with the merMod speedup in `b730b53`).
- **Vignette doc hygiene (2026-08-20)** — doc-only cleanup of shipped
  vignettes; no `R/`/roxygen/test/check impact (no `document()` needed).
  (a) `correction-error.Rmd` unwraps the `eval=FALSE` sim blocks'
  `scoring_matrix` attribute with `[[1]]` — under the default
  `format = "unified"` a single-group `scoring_matrix` attribute is a
  length-1 list, so the bare value no longer assigns into the numeric
  `a_sim` array for a reader running the block (verified: `[[1]]` yields the
  plain 1×3 / 2×6 matrix). (b) `multiple-factors.Rmd` normalizes
  `standardizedSolution` → canonical `standardizedsolution` (×3; lavaan 0.7.2
  exports both casings). (c) `R2spa.Rmd` drops the unused `a/b/c` loading
  labels from the example-1 model string so it matches the CFA the code
  actually fits. The 3 edited vignettes re-knit clean (R 4.6.1, lavaan 0.7.2,
  R2spa 0.0.4). Deferred (tracked for the multilevel re-integration, see
  `archive/PLAN_QUARANTINE.md`): the dead `(multilevel.html)` cross-ref in
  `vignettes/scoring-matrices.Rmd:35-37`.
