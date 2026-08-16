# R2spa — Open Issues & Follow-ups

Tracks remaining work from completed plans (see `archive/`). Update status as
items are resolved; move finished items to the **Closed** section with the
date and commit/PR reference.

**Last updated:** 2026-08-16 (PLAN 02 Steps 3–6 + final verification)

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

## Verification state (2026-08-16 — full PLAN 02)

- `devtools::test()`: **435 pass, 0 fail, 0 warn** (425 after Step 1; 428
  after Step 2; 434 after Step 3's tspa shape tests; 435 after Step 4's
  matrix-dispatch test).
- `devtools::document()`: OK, idempotent (pre-existing `tspa_plot.R` link
  warning — file now in `legacy/`).
- Full `R CMD build .` + `R CMD check --no-manual` (2026-08-16):
  `checking package vignettes ... OK` — **all 13 vignettes rebuild**;
  status **2 WARNINGs, 3 NOTEs** — identical to the pre-existing baseline
  (S3-consistency for `get_fs.lavaan()`/`get_fs.merMod()`, undocumented
  `fsm`/`...` Rd args, `.lintr` hidden file, unused `Matrix` Import,
  `stats::model.frame` not imported). An A/B `git stash` check around
  Step 1 showed the same 2W/3N on the pre-PLAN 02 tree; a 2026-08-16
  scoped re-check after Steps 3–5 confirmed the same 5 items verbatim.
- Vignette build history: exactly **3 of 13 failed** on the pre-Step 1
  tree (`corrected-se.Rmd`, `multilevel.rmd`, `tspa-vignette-mx.Rmd` —
  the "7/8 of 13" reports were wrong); **13/13 build on 2026-08-16**
  after Steps 2–3. Caveat: `devtools::check(args = "--ignore-vignettes")`
  still knits vignettes in the *build* phase; scope with
  `R CMD build --no-build-vignettes .` +
  `R CMD check --ignore-vignettes --no-manual <tarball>`.

## Notes

- Working-tree refactor from PLAN 01 **plus** all of PLAN 02 (Steps 1–6)
  are **uncommitted** as of 2026-08-16; the edits sit in shared files, so
  commit scope must be decided explicitly.
- Suggested order (PLAN 02 complete): 4 (user-facing bug) → 5 (perf) →
  F1 (future).
