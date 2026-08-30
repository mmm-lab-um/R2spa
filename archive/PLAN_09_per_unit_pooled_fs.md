# PLAN 09 — `tspa()` per-unit factor-score pooling (FIML per-pattern + merMod per-cluster)

**Status: PLAN** (to be implemented by assigned subagents)

Supersedes the `stop()` gate added in **PLAN 06 §4** ("`tspa()` does not yet
support groups with multiple missing-data patterns") and realizes the pooling
variant of follow-up **F3** (per-pattern stage 2, `STATUS.md`) — the robust
alternative to F3's sub-group design. Builds on **PLAN 07**'s `fs_indiv()` /
`resolve_fs_per_row()` / `fs_row_cols()` (`R/fs_indiv.R`), which already do the
"expand per-pattern/per-cluster/merMod attributes into long-form
individual-specific values" half of this work.

## 0. Task framing (from the ask)

1. **Feasibility** of `tspa()` supporting *direct* per-pattern `fsL`/`fsT` (FIML
   + merMod); the flagged hazard is small sample size for some patterns.
2. **Regardless of (1): implement** support by first converting to long-form
   individual-specific values, then reducing to a single representative value
   (the mean) that feeds the existing stage-2 machinery.
3. **Explore** whether the median is an alternative to the mean.

Confirmed user decisions (asked 2026-08-20):
- **API:** new argument `reduce = c("mean", "median")` on `tspa()`, default
  `"mean"`; a no-op unless per-unit heterogeneity is present.
- **merMod:** unweighted per-cluster (one weight per row of the `get_fs().merMod`
  output; do NOT pull cluster sizes `n_j`).
- **median:** shipped now, opt-in, with a PSD guard (warn when the pooled `fsT`
  is not positive semi-definite).
- **Scope:** multi-factor `fsT`/`fsL`/`fsb` for FIML per-pattern **and** merMod
  per-cluster (core) **+** single-factor `se_fs` for FIML. mirt per-obs is
  **out of scope** (unchanged behavior — do not alter its path).

## 1. Feasibility verdict (Task 1 — analysis only, no code)

Direct per-pattern staging is feasible **only for FIML**, by treating the
observed-indicator pattern as the stage-2 lavaan `group` (`group =
c(<user group>, "__r2spa_pat__")`), passing a pattern-length list for
`fsT`/`fsL`/`fsb` (the schema already renders per-group `c(...)` values), and
pooling structural parameters across patterns.

Key nuance on the flagged hazard:
- Small patterns hurt **only if structural / latent-variance / latent-mean
  parameters are estimated separately per pattern** (under-identified in tiny
  strata → non-convergence).
- If measurement parameters stay **fixed per pattern** (they are — taken
  verbatim from `fsT`/`fsL`) and structural parameters are **equality-constrained
  (pooled)**, small patterns are harmless — they contribute fixed-measurement
  observations to the pooled MLE.

Why we still pick pooling (Task 2) over the direct path:
1. The owned schema (`tspa_schema_mf` / `tspa_render`) treats the user model as a
   verbatim "raw" row and injects only per-group measurement values. Pooling
   structural paths across pattern-groups would require label-equality surgery
   on arbitrary user syntax (or `group.equal`, which constrains by parameter
   *type* — imprecise / surprising semantics).
2. `2^k` pattern strata, some near-empty → per-group setup overhead and numeric
   fragility.
3. No generalization: mirt is per-observation (cannot be a `group`); merMod is
   per-cluster ("one group per cluster" is the small-group pathology itself).

The **pooling** approach reuses every stage-2 mechanism, sidesteps the
small-sample issue by pooling information, and is the path implemented here.

## 2. Core design (Task 2)

The "convert to long-form individual-specific values" step **already exists**:
- `resolve_fs_per_row(fs)` (`R/fs_indiv.R:217`) dispatches on input shape
  (lavaan unified / lavaan list / merMod 3-D `fsT`/`fsL` / mirt per-obs marker)
  and returns `list(n, scores, pattern_idx, blocks, group_col, group_vals,
  id_vals, legacy)`, where `blocks` is one `(fsL, fsT, fsb)` value per
  pattern (lavaan) or per cluster (merMod) and `pattern_idx` maps each row to a
  block. All-missing (NA) rows are routed to a dedicated all-NA block.
- `fs_row_cols(fs, fsL, fsT, fsb)` (`R/fs_indiv.R:182`) returns an `n × K`
  name-free numeric matrix per block: `[1:q]` = per-row SE,
  `[(q+1):(q+q^2)]` = `c(as.matrix(fsL))` (column-major, per latent: loadings of
  that latent on each score), `[(q+q^2+1):(q+q^2+q(q+1)/2)]` = `fsT` lower
  triangle (row-major, i-outer / j<=i-inner, the r2spa column order), and
  (if `fsb`) the last `q` = `fsb`.

**Reduce = per-observation-equal** (the "long-form" reading): repeat each
unit's matrices once per member row, then average. For a group with patterns
`p1..pm` and member counts `c1..cm` this equals the case-count-weighted pattern
mean (for the mean); for the median it is the element-wise median of that same
expanded multiset. `na.rm = TRUE` in the reduction drops the all-NA rows the
engine already emits.

**PSD:** the **mean** of PS covariance matrices (a convex combination) is
guaranteed PS → the pooled `fsT` stays a valid covariance and the stage-2
measurement block (fixed `fs_i ~~ fs_j` errors) stays proper. The **median** is
element-wise and has no such guarantee → a non-PSD pooled `fsT` can make the
stage-2 model improper / non-convergent. Hence: `mean` default; `median` opt-in
with a PSD warning.

## 3. `R/tspa.R` — code (owned by **r-architect**)

New signature:
```r
tspa <- function(model, data, reliability = NULL, se = "standard",
                 se_fs = NULL, fsT = NULL, fsL = NULL, fsb = NULL,
                 reduce = c("mean", "median"), ...) {
```
`reduce <- match.arg(reduce)` near the top.

### 3.1 Internal helpers (co-located in `R/tspa.R`, `snake_case`, no `library()`)

- **`is_per_unit_fs(fsT, fsL)`** → `TRUE` when per-unit heterogeneity is present
  and poolable:
  - **merMod:** `is.array(x) && length(dim(x)) == 3L` for either `fsT`/`fsL`; OR
  - **lavaan per-pattern:** `is.list(x)` and ∃ element `e` with
    `is.list(e) && length(e) > 0L && all(vapply(e, is.matrix, logical(1)))`.
  - Do **NOT** trigger on mirt's flat per-obs list (a length-n list of *bare*
    matrices, not nested per group) — keep that path out of scope.

- **`pool_per_unit(fs, reduce, have_int)`** → `list(fsT, fsL, fsb)`:
  - `resolved <- resolve_fs_per_row(fs)` (errors clearly if `fs` isn't a
    resolvable `get_fs()` result — the "data must be the get_fs() result"
    constraint).
  - `q <- ncol(resolved$blocks[[1]]$fsT)`; assemble per-row `se_mat (×q)`,
    `ld_mat (×q^2)`, `ev_mat (×q(q+1)/2)`, `int_mat (×q)` exactly as
    `fs_indiv()`'s block loop (`R/fs_indiv.R:94-113`); use `fs_row_cols(resolved$scores[rows, ,
    drop = FALSE], blk$fsL, blk$fsT, if (have_int) blk$fsb else NULL)`.
  - `reduce_fn <- if (reduce == "median") stats::median else mean`; for every
    block column `reduced[,k] <- reduce_fn(mat[,k], na.rm = TRUE)`.
  - Reassemble (name-free, no dependence on `legacy` column names):
    - `fsT`: take the row-major lower-triangle (`ev`) reduced values into a
      symmetric `q × q` (position `(i,i)` and `(i,j), j<=i`, mirror to `(j,i)`)
      with rownames/colnames = `rownames(blocks[[1]]$fsT)`.
    - `fsL`: `matrix(reduced_ld, nrow = q, ncol = q)` (column-major) with
      rownames = `rownames(blocks[[1]]$fsL)` (scores), colnames =
      `colnames(blocks[[1]]$fsL)` (latents).
    - `fsb`: reduced `int_mat` column (named vector), or `NULL` when `!have_int`.
  - **Per-group vs single:** if `resolved$group_vals` is non-null (lavaan MG),
    reduce *within each group* (`which(group_vals == g)`) and return **named
    lists** (one single matrix/vector per group, named by group label).
    Otherwise return **a single** matrix/vector (SG FIML, merMod, mirt). The
    returned shape is exactly what the existing schema already accepts
    (multigroup list or single matrix).
  - **PSD guard** (single- and per-group `fsT`, cheap):
    `emin <- min(eigen(T, symmetric = TRUE, only.values = TRUE)$values)`; if
    `!is.finite(emin) || emin < -.Machine$double.eps^0.5` → `warning()`
    ("pooled fsT is not positive semi-definite; consider reduce = \"mean\"").
    (Always satisfied for `reduce = "mean"`; can fire for `"median"`.)

- **`pool_se_fs(data, se_names, reduce, group_col)`** → per-group `se_fs`
  (single-factor FIML): `reduce_fn(data[[paste0("fs_", v, "_se")]])` per latent
  `v`, `na.rm = TRUE`; per-group when `group_col` non-null (a data frame: one
  row per group × latent cols) else a numeric vector. Reads the materialized se
  columns directly (robust to the `cbind`-combined `@examples` frames).

### 3.2 Wiring

**Multi-factor path** (`is.null(fsT)` is FALSE, ~`R/tspa.R:227`):
1. Keep the both-or-none check (`xor(is.null(fsT), is.null(fsL))`, `:147-165`)
   — and note `3`-D arrays / nested lists still satisfy "both provided".
2. **Replace** the per-pattern hard stop (`:166-199`) with: if
   `is_per_unit_fs(fsT, fsL)`, then
   `pooled <- pool_per_unit(data, reduce, have_int = !is.null(fsb))`; set
   `fsT <- pooled$fsT; fsL <- pooled$fsL; if (!is.null(fsb)) fsb <- pooled$fsb`.
   Do this **before** `multigroup <- ...` (`:201`) and the name-match (`:207-221`)
   so the existing detection/`tspa_mf` run on clean per-group (or single)
   matrices. Collapsing the 3-D merMod array to a 2-D `fsT` also fixes the
   `upper.tri()` call in `tspa_schema_mf` (`:364`) — **this enables `tspa()` for
   merMod factor scores** (currently broken).
3. If not per-unit, keep the both-or-none error and (optionally) a residual
   backstop. Otherwise unchanged.
4. Attach the **pooled** `fsT`/`fsL` + a marker (e.g. `attr(fit, "pooled_fs") <-
   <reduce or TRUE>`) on the returned fit (extend `:247-250`) for transparency.

**Single-factor path** (`is.null(fsT)` TRUE, `:233`): after
`se_fs <- as.data.frame(as.list(se_fs))` (`:144-146`), detect within-group se
variation as the FIML signal (a complete group's `fs_<v>_se` column is
constant within the group ⇒ no false positive):
`vapply over se_fs columns: ∃ group with `sum(!is.na(unique(col))) > 1``.
If `group_col <- attr(data, "group_col")` (or `"group"`) is in `names(data)` and
variation is present → `se_fs <- pool_se_fs(data, colnames(se_fs), reduce,
group_col)`. Else leave `se_fs` as-is (unchanged homogeneous behavior). Requires
`data` to carry the `fs_*_se` columns (i.e. the `get_fs()` result); if not,
leave `se_fs` as-is (no change) rather than erroring.

**Unchanged:** the `data` frame passed to `sem()` (extra `ev_*`/`_by_*`/`_se`
columns are already ignored by lavaan), the user still passes `fsT`/`fsL` to
select the mf path, MG still requires `group =` in `...` (`:223-225`), `se`/
`reliability` handling, `tspa_schema_*`/`tspa_render`.

**No new imports** (base `mean`/`median`/`eigen`; `stats::median`, `stats::setNames`
already in scope). No `library()`/`require()` in bodies. No new `NAMESPACE`
exports (helpers stay internal; `fs_indiv`/`resolve_fs_per_row`/`fs_row_cols`
are package-internal and callable directly).

## 4. Documentation (owned by **r-doc** — roxygen only, `R/tspa.R`)

- **`@param reduce`** — full prose (r-architect adds a minimal `@param reduce`
  line so R CMD check passes; r-doc refines it): the reduction of per-unit
  `fsL`/`fsT`/`fsb` (FIML per-pattern) or per-cluster (merMod) to a single
  representative per-group value; no-op for homogeneous inputs; `mean` preserves
  positive semi-definiteness, `median` may not (a warning is emitted when it
  does).
- **`@param fsT`/`fsL`/`fsb`** — remove "Must not contain per-pattern lists
  (groups fitted with missing data): `tspa()` does not yet support that." and
  state that per-pattern (missing-data) and merMod per-cluster values are now
  accepted and pooled per group via `reduce` (the pooled — not the nested —
  values are attached to the returned fit).
- **`@details`** — short paragraph: per-pattern/per-cluster values → long-form
  individual-specific → pooled per group; why (small patterns), and the
  mean-vs-median PSD tradeoff + that the pooled values are attached as fit
  attributes.
- **`@examples`** — add examples that actually run under `devtools::run_examples` /
  `R CMD check`:
  - a **single- or multi-group FIML** two-factor example (NA-inject two
    indicators with a fixed seed, fit with `missing = "fiml"`, `get_fs(...)`,
    `tspa(..., fsT = attr(...,"fsT"), fsL = attr(...,"fsL"), reduce = "mean")`);
    and a `reduce = "median"` invocation to surface the argument.
  - a **merMod** example: `library(lme4); lmod <- lmer(Reaction ~ Days +
    (Days | Subject), sleepstudy); fs_mer <- get_fs(lmod); tspa("u1 ~ u0",
    data = fs_mer, fsT = attr(fs_mer,"fsT"), fsL = attr(fs_mer,"fsL"))`.
  - keep existing examples intact.

## 5. Tests (owned by **r-tester** — `tests/testthat/test-tspa_pooled.R` only)

Reuse fixtures from `tests/testthat/test-fs_indiv.R` (`mk_fiml_hs`,
`fit_fiml`/`fs_fiml`, `mk_fiml_all`/`fs_fiml_all`, merMod `lmod`/`fs_mer`,
helpers `rowcols`, `sg_attr`, `symm_from_lower`). FIML **multi-factor** needs a
two-factor CFA with NA-injection (extend the single-factor fixture; keep a fixed
seed). Cases:
- **FIML mf SG + MG:** `tspa()` runs without error; the attached
  `attr(fit,"fsT")` equals the case-count-weighted (per-row) mean of the
  per-pattern `fsT` (rebuild with `symm_from_lower`), is PSD, and matches a
  hand computation; same for `fsL`. For MG, one pooled `fsT` per group.
- **merMod:** `tspa()` runs (previously broke at `upper.tri`); pooled
  `fsL`/`fsT` `== colMeans`/mean of the 3-D per-cluster slices (unweighted);
  single group.
- **Single-factor FIML `se_fs` path:** `tspa()` on `se_fs`-style input with
  FIML data pools `fs_<v>_se` (== `mean(..., na.rm = TRUE)`); a complete-data
  single-factor input is byte-identical to the pre-change path.
- **Identity guards:** single-group complete data and multigroup complete data
  → **identical** stage-2 model string/fit to the legacy un-pooled path (the
  pool is a no-op when values are constant within group).
- **`reduce = "median"`:** differs from mean only when values vary; the PSD
  warning path fires on a **contrived** non-PSD per-unit set (e.g. hand-built
  nested `fsT` feeding `tspa(reduce = "median")`), mean does not warn.
- **All-missing (NA-pattern) rows** are excluded via `na.rm`, not dropped (pooled
  value equals the hand mean over non-NA rows only).
- **mirt per-obs unchanged:** passing a per-obs mirt `fsL`/`fsT` list is *not*
  routed into pooling by `is_per_unit_fs` (out of scope) — assert it does not
  hit the pooler (keeps existing behavior / error path).
- **`data` must be the get_fs result:** a `data` frame without resolvable
  attributes + a nested `fsT` → informative error (not a silent mis-pool).

Keep edition-3 `expect_equal(..., tolerance = ..., ignore_attr = TRUE)` style.
Do NOT edit anything under `R/` or `man/`.

## 6. Verification lifecycle (AGENTS.md order; owner noted)

1. **r-architect:** `devtools::load_all()`, `devtools::document()` (roxygen
   changed), `devtools::test()` (confirm no regressions in
   `test-tspa.R`/`test-tspa_render.R` — all complete-data, expected to pass
   untouched).
2. **r-doc:** `devtools::document()`, then run the new `@examples`
   (`devtools::run_examples()` / knit the example snippets) to confirm they run.
3. **r-tester:** `devtools::test()` (targeted `test-tspa_pooled` first, then the
   full suite; report pass/fail counts).
4. **Orchestrator:** `devtools::check()` → expect **0 errors / 0 warnings / 1
   NOTE** (the known OpenMx `Imports` baseline). Then update `STATUS.md`:
   close **F3** as delivered via the `reduce`-pooling variant, and log that the
   *sub-group* (per-pattern-free-structural) variant remains future/out of
   scope with its `group = c(a, b)` lavaan canary caveat (former F3 text).

## 7. Subagent role assignments

| Subagent | Owns | Files it may edit | Verifies |
|---|---|---|---|
| **r-architect** | Core strategy, cross-file consistency, deps, build lifecycle | `R/tspa.R` (code + signature + minimal `@param reduce`) | `load_all`+`document`+`test` (no regressions) |
| **r-doc** | roxygen2 documentation maintenance | roxygen in `R/tspa.R` (`@param reduce`, `@param fsT/fsL/fsb`, `@details`, `@examples`) + regenerated `man/` | `document()` + examples run |
| **r-tester** | testthat implementation | `tests/testthat/test-tspa_pooled.R` (new) | `devtools::test()` |

**Execution order (serial, to respect the document→test→check lifecycle and the
shared `R/tspa.R`):** r-architect → r-doc → r-tester → orchestrator final check.
(r-doc and r-tester both finish after r-architect; r-tester may author the test
file in parallel with r-doc since they touch different files, but the full
`devtools::test()` + final `check()` run serially once to avoid `man/`/NAMESPACE
write contention.)

## 8. Out of scope (follow-ups)

- **Sub-group per-pattern stage 2** (former **F3**): free-structural-per-pattern
  lavaan sub-groups — deferred (fragile for tiny patterns; needs
  `group = c(a,b)` lavaan canary). The delivered pooling is the robust default.
- **mirt per-obs `tspa()`** — the per-obs `fsL`/`fsT` list is *not* pooled
  (`resolve_fs_per_row` supports it, but `is_per_unit_fs` deliberately does not
  trigger on it); remains unsupported, as today.
- **cluster-size weighting for merMod** (chosen: unweighted per-cluster).
- Vignette regeneration not required (no narrative asserts the removed error);
  a vignette *addendum* showing the FIML/merMod 2S-PA flow is a natural follow-up.

## 9. Risks / notes

- **Complete-data regression surface:** complete SG/MG fits must be byte-identical
  (the pool is a no-op when values are constant within a group); the identity
  tests guard this. `test-tspa.R` / `test-tspa_render.R` use complete data and
  must pass untouched.
- **`data` contract:** pooling requires `data` to be the (unmodified) `get_fs()`
  result carrying `fs_pattern` (lavaan) or the 3-D `fsT`/`fsL` (merMod). A
  stripped/subset/`cbind` frame that isn't resolvable → informative error (mf)
  or fall-through to unchanged `se_fs` (sf).
- **PSD of the median:** documented + guarded by the warning; `mean` is the safe
  default. `mean` of PS matrices is PS (convex), so the mf measurement block
  stays proper for the default path.
- **roxygen2 is 8.1.0** in this environment (no `RoxygenNote`/`importFrom`
  churn). Never hand-edit `NAMESPACE`/`man`.
- **quarantine untouched:** `vcov_corrected`/`tspa_mx`/`grandStandardizedSolution`
  remain out of the build; they read `fsT`/`fsL` as single matrices and will need
  the same "pool first" treatment at re-integration (note for that plan).
