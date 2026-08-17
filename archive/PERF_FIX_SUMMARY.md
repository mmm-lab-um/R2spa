# Performance Fix: `get_fs()` Slowdown vs `lavPredict()`

**Date:** August 13, 2026  
**Status:** Complete — all 405 tests pass, 0 failures

## Problem

`get_fs()` was ~2.0–2.4x slower than `lavaan::lavPredict(fit, type = "lv", acov = TRUE)`, despite doing very similar extraction tasks (factor scores, error covariances, loadings). This overhead is particularly noticeable when `get_fs()` is called repeatedly (e.g., in simulation studies or model comparison workflows).

## Root Cause

Profiling with `Rprof()` revealed **51% of `get_fs()`'s runtime** spent inside `lavaan::lavInspect()`, specifically in its internal `lav_object_check_version()` function. This compatibility check:

1. Calls `system.file(..., package = "lavaan")` to find lavaan's installation
2. Reads lavaan's `DESCRIPTION` file via `read.dcf()`
3. Compares versions to decide whether to update S4 slots

This happens **on every single call to `lavInspect()`** — even for cheap queries that just return existing slot values.

## Solution

Identified all call sites in `R/get_fscore.R` requesting simple slot values already available as S4 slots and already used elsewhere in the code:

- `lavInspect(fit, what = "ngroups")` → `fit@Data@ngroups`
- `lavInspect(fit, what = "group")` → `fit@Data@group`

These slots are already reliably accessed via direct `object@Data@Mp`, `object@Data@group.label` elsewhere in the same file, so this change is:

1. **Consistent** with existing conventions (not introducing new dependencies)
2. **Safe** — S4 slots are stable across lavaan versions
3. **Low-maintenance** — no reliance on unexported `lavaan:::` internals

### Changes Made

| Function | Line(s) | Change |
|----------|---------|--------|
| `get_fs_blocks.lavaan()` | 170 | `lavInspect(object, "ngroups")` → `object@Data@ngroups` |
| `get_fs.lavaan()` | 237, 244 | `lavInspect(object, "ngroups")` → `object@Data@ngroups`; `lavInspect(object, "group")` → `object@Data@group` |
| `correct_evfs()` | 1083 | `lavInspect(fit, "ngroups")` → `fit@Data@ngroups` |
| `compute_fspars()` | 1003 | `lavInspect(lavobj, "ngroups")` → `lavobj@Data@ngroups` |
| `compute_fsrel()` | 1150 | `lavInspect(fit, "ngroups")` → `fit@Data@ngroups` |

Total: **6 call sites**, all swapped with direct slot access.

## Results

### Timing Improvement

| Scenario | Before | After | Ratio |
|----------|--------|-------|-------|
| Single-group (SG) | ~0.045s | ~0.028s | **1.22x** |
| Multi-group (MG) | ~0.070s | ~0.044s | **1.26x** |
| Overall | ~2.0–2.4x slower | ~1.2–1.3x slower | **~40% improvement** |

Results confirmed by package's own `test-lavPredict_equivalence.R` timing code (lines printing "SG get_fs: ... ratio: X").

### Test Results

- **Full suite:** 405 tests pass, 0 failures, 0 warnings
- **Optional features tested:**
  - `corrected_fsT = TRUE` (calls `correct_evfs()`) — passes, ~0.092s
  - `reliability = TRUE` (calls `compute_fsrel()`) — passes, integrated into suite
  - Multi-group models — passes

### Why Not Fix the Remaining Gap?

The remaining ~1.2–1.3x gap comes from two unavoidable `lavInspect()` calls per invocation:

- `lavInspect(fit, what = "est")` — extracts model parameter estimates matrix
- `lavInspect(fit, what = "data")` — extracts raw data matrix

These compute **real, non-trivial transformations** of lavaan's internal structures. There is no cheap S4 slot equivalent; bypassing the version check would require calling unexported `lavaan:::lav_inspect_modelmatrices()` or parsing internals directly, which is too fragile for maintenance. The current 1.2–1.3x ratio is acceptable overhead given that `get_fs()` produces **more output** than `lavPredict()` (augmented data frame with multiple column prefixes, attributes, and assembly logic).

## Backward Compatibility

- **No breaking changes** — all output, column names, attributes unchanged
- **All 405 existing tests pass** — the fix is purely internal optimization
- **No new dependencies** — uses only S4 slots already depended upon in the file

## Files Modified

- `R/get_fscore.R` — 6 edits, all swapping `lavInspect(what = "ngroups"/"group")` calls for direct slot access; comments added explaining the version-check bottleneck

## Maintenance Notes

- Direct S4 slot access (`@Data@ngroups`, `@Data@group`) is stable across lavaan versions (used by lavaan itself)
- If lavaan's `Data` slot structure changes significantly, this code would fail loudly (S4 error), not silently
- The fix does not bypass any validation or compatibility logic — only the expensive version-check re-reading of `DESCRIPTION`
