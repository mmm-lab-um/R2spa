# Plan: `get_fs()` preserves `fsT`/`fsL`/`fsb` per missing-data pattern

## Context

With missing data, lavaan partitions each group's cases into one **block per
distinct observed-indicator pattern** (`@Data@Mp`), scored on only that
pattern's indicators. `assemble_fs_blocks()` (R/get_fscore.R:348-470) currently
keeps only block 1's matrices as the group-level `fsT`/`fsL`/`fsb`/
`scoring_matrix` attributes and prints
"blocks have differing fsT/fsL/fsb attributes ...". The per-row columns are
already pattern-correct (filled from each case's own block, R/get_fscore.R:368-380);
only the group-level attributes are lossy.

Goal: preserve per-pattern matrices in the attributes, store the pattern/case
membership so a future API can append individual-specific loadings/error
covariances as columns, and make downstream unsupported combinations fail
explicitly. **`tspa()` does not gain missing-data support in this change** —
it gets an explicit error, and the per-pattern stage-2 design is recorded as a
future plan (Section 7).

Confirmed decisions (user):
- tspa(): do NOT support missing data now; explicit `stop()` + future-plan note.
- Attribute shape follows `lavPredict(acov = TRUE)` nesting (verified on
  lavaan 0.7.2): per-group list; group value is a plain matrix when the group
  has one pattern, a per-pattern list when it has k≥2. We name the patterns
  (lavaan leaves them unnamed).
- Pattern membership **must be stored** (per-case pattern label + pattern
  definitions) for the future column-append API.
- `corrected_fsT`/`reliability`/`vfsLT` + missing patterns → explicit
  informative `stop()` (today: cryptic "incorrect number of dimensions").

## 1. Attribute contract

For lavaan objects, all four attributes (`fsT`, `fsL`, `fsb`, `scoring_matrix`)
plus a new `fs_pattern`:

- Group with **k = 1** pattern (no missing data): identical to today — group
  element is a plain matrix/vector. All complete-data outputs unchanged.
- Group with **k ≥ 2** patterns: group element is a **named list of k
  matrices**, one per pattern. Pattern name = observed indicators joined with
  `"+"` in indicator order (e.g. `"x1+x2+x3"`, `"x4"`). Unique per pattern.
- New **`fs_pattern`**: named list per group; each element
  `list(label = <character vector, one entry per case in the group, its
  pattern name>, pat = <named logical p × k matrix: rows = indicators,
  cols = pattern names>)`. For complete groups: single label (all indicators),
  single all-TRUE column. This mirrors lavaan's `@Data@Mp` (`pat`, `case.idx`)
  in public-attribute form and is what a future per-row column API indexes on.

Applies uniformly to `format = "unified"` and `"list"`, single- and
multi-group. The message and `check_blocks_identical()` (R/get_fscore.R:338-346)
are deleted — nothing is dropped anymore.

## 2. `R/get_fs_methods.R`

- `prepare_fs()` (lines 142-195): on each pattern block, carry
  `pat_label = paste0(colnames(y)[pat_m], collapse = "+")` and
  `pat = setNames(pat_m, colnames(y))`.
- `get_fs.lavaan()` (~line 294): before the math, if
  `any(vapply(object@Data@Mp, function(m) !is.null(m), logical(1)))` and
  (`corrected_fsT` || `reliability` || `vfsLT`) → `stop()` with a clear
  message that these SE paths are not supported with multiple missing-data
  patterns (replace the cryptic failure deep in `compute_fspars()`/
  `correct_evfs()`).

## 3. `R/get_fscore.R`

- `assemble_fs_blocks()` (lines 348-470), per group:
  - if `length(blocks) > 1`: for each `ak ∈ attr_keys`,
    `setNames(lapply(blocks, `[[`, ak), vapply(blocks, `[[`, "pat_label"))`
    (fall back to `paste0("pattern_", i)` when a hand-built block lacks
    `pat_label` so `test-assemble_fs_blocks.R` fixtures work);
    build per-row `label` vector from `case_idx` + `pat_label`; group's
    `fs_pattern` = `list(label = ..., pat = do.call(cbind, pat list))`.
  - else: exactly today's behavior (first/only block's matrices) plus
    `fs_pattern` with the single all-observed label.
  - delete the message block (lines 395-412) and `check_blocks_identical()`.
  - add `"fs_pattern"` to `attr_keys` (line 358).
- `fs_to_group_list()` (lines 186-294): add `"fs_pattern"` to `attr_keys`
  (line 187); generic pass-through already handles nested values — verified by
  reading both directions (element attrs at lines 231-238, reconstruction at
  lines 274-289).
- Roxygen for `get_fs` `@return` (lines 75-86) and `fs_to_group_list`
  `@return`: document the nested shape and `fs_pattern`.

`augment_fs()` and per-row columns: unchanged.

## 4. `R/tspa.R` — guard only (no missing-data support)

- After the existing `fsT`/`fsL` shape checks (lines 141-160): detect a
  per-group value that is itself a list of matrices (i.e. a group with k≥2
  patterns — a length-1 list element that is a list, or any group element
  `is.list()` and not a bare vector for `fsb`) → `stop()` with an explicit
  message: stage-2 `tspa()` does not yet support groups with multiple
  missing-data patterns; use complete data or see the planned support
  (Section 7). This is required because left unchecked, a single-group
  k=3-pattern fit would make `attr(fs, "fsT")` a length-3 list, silently
  misread as 3 groups (lines 148-165), and a multigroup fit would fail
  cryptically at `T_list[[g]][i, j]` (R/tspa.R:333).
- `tspa` roxygen `@param fsT/fsL/fsb`: note the current limitation.

Everything else in `tspa()` (schema, render, sf/mf paths) unchanged.

## 5. Tests

- `test-assemble_fs_blocks.R`:
  - replace "representative fsT from first block when attributes differ"
    (lines 103-109) with: per-pattern named-list contract — `attr(res, "fsT")`
    group element is a named list whose matrices equal the blocks' `fsT`;
    labels correct; **no message** (`expect_no_message`);
  - new: `fs_pattern` attribute content (per-row labels, `pat` matrix) for
    the two-block fixture; k=1 case keeps plain-matrix shape.
- New `test-get_fs_missing.R` (NA-injected `HolzingerSwineford1939`,
  `set.seed(1334)`, NAs on cols 7:9 — the existing test recipe):
  - single-group: `attr(fs, "fsT")` group element = named list; each pattern's
    matrix equals `lavPredict(fit, acov = TRUE)` acov pattern matrix, **matched
    by pattern label** (lavaan's pattern order is case-count-descending and
    internal — never assume index alignment);
  - multigroup (`group = "school"`): both groups nested; labels correct;
  - `fs_pattern` consistency: per-row label agrees with which raw indicators
    are NA for that row;
  - complete-data regression guard: k=1 groups keep the plain-matrix shape
    and all current attribute values.
- `test-get_fs_lavaan_SE_paths.R` or extend `test-get_fscore.R`:
  `corrected_fsT`/`reliability`/`vfsLT` on multi-pattern data →
  `expect_error(regexp = ...)` with the new informative message.
- `test-tspa.R`: passing nested (multi-pattern) `fsT`/`fsL` to `tspa()` →
  `expect_error(regexp = ...)`; single-group k>1 attribute also rejected
  (the "misread as 3 groups" trap); all existing complete-data tspa tests
  untouched.
- `test-fs_converters.R`: unified↔list round-trip preserving `fs_pattern`
  and nested per-pattern lists (both directions, single- and multi-group).
- `test-lavPredict_equivalence.R` / `test-get_fs_priors.R`: no changes needed
  (they call `get_fs_blocks.lavaan()` directly; block-level assertions still
  hold).

## 6. Verification lifecycle (AGENTS.md order)

1. `devtools::load_all()` after R edits.
2. `devtools::document()` (roxygen changed) — NAMESPACE/man regenerated, never
   hand-edited.
3. `devtools::test()`.
4. `devtools::check()` — expect 0 errors / 0 warnings, the one known
   OpenMx NOTE.

## 7. Future plan (recorded for a later change, NOT part of this one)

**a) `tspa()` per-pattern stage 2.** When a group has k≥2 patterns, fit stage 2
with lavaan sub-groups `(group × pattern)`: synthesize an internal pattern
column per row from `fs_pattern` (e.g. `__r2spa_pat__`), call
`sem(..., group = c(<user group var>, "__r2spa_pat__"))`, and fix each
sub-group's loadings/error (co)variances/intercepts to its own pattern's
matrices (extend `tspa_schema_mf()` group indices to (group, pattern) pairs).
Structural paths estimated per sub-group; user `group.equal` still applies.
Known hazard to verify first: lavaan's exact level ordering for
`group = c(a, b)` (lavaan 0.7.2) — the schema value-order mapping must be
derived empirically and pinned with a canary-style test
(cf. R/lavaan_compat.R). Limitation to document: tiny pattern sub-groups
(n ≲ 5) may fail to converge; no automatic merging.

**b) Per-row column API (user-requested).** A user-level function that appends
individual-specific `fsL`/`fsT`/`fsb` values as columns to the factor-score
data frame, keyed by `fs_pattern`'s per-row labels (e.g. `ld_fs_visual_x1`,
`ev_fs_visual`), so each case carries the matrices matching its own pattern.

## Risks / notes

- Complete-data outputs are byte-compatible (k=1 path unchanged) — the main
  regression surface is the nested attribute consumption in
  `test-tspa.R`/`test-tspa_render.R`, all of which use complete data and are
  expected to pass untouched.
- `fs_pattern` row order must match the group data-frame row order
  (`case_idx` are group-local row numbers in `assemble_fs_blocks`).
- merMod path and quarantined code untouched; the quarantined consumers
  (`vcov_corrected`, `tspa_mx_model`, `grandStandardizedSolution`) read
  `attr(..., "fsT"/"fsL")` as single matrices and will need adaptation at
  re-integration (note for `.quarantine` re-integration plan).
- Vignettes: none assert the old message or missing-data attribute shapes
  (grepped) — no fixture regeneration needed.
