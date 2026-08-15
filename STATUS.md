# R2spa — Open Issues & Follow-ups

Tracks remaining work from completed plans (see `archive/`). Update status as
items are resolved; move finished items to the **Closed** section with the
date and commit/PR reference.

**Last updated:** 2026-08-15 (from review of PLAN 01)

## Open

| # | Priority | Issue | Where | Source | Status |
|---|----------|-------|-------|--------|--------|
| 1 | P1 | **merMod column-name regression** — `get_fs_lmer()` output renamed vs legacy: `u0_eb`→`fs_u0_eb`, `u0_by_u0_eb`→`u0_eb_by_fs_u0_eb`, `ecov_u0_eb_u1_eb`→`ecov_u1_eb_u0_eb`; breaks `vignettes/multilevel.rmd` (`tspa_mx_model` `mat_ld`/`mat_ev`) and downstream user code. Violates PLAN 01's back-compat spec; tests are name-blind (`ignore_attr = TRUE`) so they don't catch it. Restore legacy names in the merMod path. | `R/get_fs_methods.R` (`get_fs_blocks.merMod()` / `augment_fs` interaction) | PLAN 01, remaining issue 1 | open |
| 2 | P1 | **7 vignettes break on `format = "unified"` default** — vignettes that expected a list from `get_fs(..., group = ...)` (e.g. `do.call(rbind, ...)`) now get a unified data frame. Either pass `format = "list"` or update the vignette narrative. `R CMD check` stays red until fixed. | `vignettes/` | PLAN 01, remaining issue 2 | open |
| 3 | P2 | **`get_fs(data = matrix)` dispatch** — goes to `get_fs.default()` instead of `get_fs.data.frame()` (1 vignette depends on the data.frame path). | `R/get_fs_methods.R` | PLAN 01, remaining issue 3 | open |
| 4 | P2 | **`reliability = TRUE` + single-group multi-factor + default `unified` format → hard error** ("reliability is only supported for unidimensional models") instead of the documented warning. Guard at `get_fs_methods.R:183-184` checks `length(attr(out, "fsb")) > 1`, but in unified format `fsb` is always a length-1 list, so `multifactor` never fires. | `R/get_fs_methods.R` (`get_fs.lavaan()` reliability block) | PLAN 01, remaining issue 4 | open |
| 5 | P2 | **`correct_evfs()` is ~100% of `corrected_fsT = TRUE` runtime** (0.167 s of 0.168 s, 4-LV model). `q × (npar+1) × 2 + 1` `lavInspect` file-reads per call (57 traced on a 13-param model). Fix: one `lav_func_jacobian_complex()` over the full `a` matrix (q× fewer evaluations) + hoist `lavInspect("free")`/`("est")` out of the per-evaluation `compute_fspars()` body. | `R/get_fscore_math.R` (`correct_evfs()` :366-401, `compute_fspars()` :289-290) | PLAN 01, remaining issue 5 | open |
| 6 | P3 | **Missed cheap-lookup `lavInspect` sites** (each pays the `read.dcf` version-check cost): `tspa_corrected_se.R:21` (`"ngroups"`), `get_fs_methods.R:197` (`"norig"`), `grandStandardizedSolution.R:85` (`"nobs"`). All have `@Data@` slot equivalents; `norig`/`nobs` slots are lists — use `unlist()`. | 3 files | PLAN 01, remaining issue 6 | open |
| 7 | P3 | **Dead-code cleanup** — rownames assignment at `get_fs_methods.R:247` overwritten at :248; commented-out line `get_fscore.R:247`; unreachable `is.data.frame` check in `get_fs.data.frame` (:14); `re_names` computed twice between `get_fs_blocks.merMod()`/`get_fs.merMod()`. | `R/get_fs_methods.R`, `R/get_fscore.R` | PLAN 01, remaining issue 9 | open |

## Follow-ups (future work, explicitly out of scope for their plans)

| # | Item | Source | Status |
|---|------|--------|--------|
| F1 | Wire `tspa()` to accept unified/list factor-score shapes natively (via `fs_to_group_list()` or dual support) — currently relies on attributes; the `format` default change shifts multigroup output to list-valued attributes that `tspa()` must keep parsing. | PLAN 01, remaining issue 7 | deferred |
| F2 | `mirt` S3 method (`get_fs.SingleGroupClass`) — architecture (blocks + assembler) is backend-agnostic and ready. | PLAN 01, remaining issue 8 | deferred |

## Closed

| # | Issue | Closed | Reference |
|---|-------|--------|-----------|
| — | *(none yet)* | | |

## Verification state (2026-08-15)

- `devtools::test()`: 420 pass, 0 fail.
- `devtools::document()`: OK (pre-existing `tspa_plot.R` link warning — file
  now in `legacy/`).
- `devtools::check()`: **FAILS at vignette build** (8/13 vignettes) — open
  items 1 and 2.

## Notes

- Working-tree refactor from PLAN 01 is **uncommitted** as of 2026-08-15.
- Suggested order: 1 → 2 (unblocks `check`) → 4 (user-facing bug) → 5/6
  (perf) → 7 (dead code) → F1 (future).
