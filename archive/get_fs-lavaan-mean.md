# PLAN: `method = "mean"` (sum-score) for `get_fs.lavaan()` / `get_fs.data.frame()`

Status: **DONE 2026-08-18** — implemented, tested (final suite 1211 pass / 0
fail), documented, `R CMD check` clean for code/docs/tests (0 errors / 0 real
warnings / 1 pre-existing OpenMx NOTE; full vignette build currently blocked by
a transient GitHub rate-limit on `tspa-growth-vignette.Rmd`'s ECLS download —
unrelated to the code). Approved decisions at bottom; C1–C8 corrections
applied.
Sequencing: this is Plan B of the pair (execute after `get_fs-mermod-ml.md` = Plan A),
sequentially, single agent.

## Context (for a fresh session)

- Repo `/home/marklai/R2spa`. Test baseline verified 2026-08-17: **707 pass /
  0 fail / 0 warn / 0 skip**.
- Lifecycle (AGENTS.md, never reorder): `load_all()` → `document()` (roxygen changed
  here) → `test()` → `check()`. Never hand-edit `NAMESPACE`/`man/*.Rd`.
- lavaan code path: `get_fs.lavaan()` (`R/get_fs_methods.R:251-371`) →
  `get_fs_blocks.lavaan()` (`R/get_fs_methods.R:129-227`) → `compute_fscore()`
  (`R/get_fscore_math.R:236-286`) → `compute_a_from_mat()`
  (`R/get_fscore_math.R:361-381`). `normalize_fs_method()` at
  `R/get_fs_methods.R:5-8` maps `ML→Bartlett`, `EB→regression` (lavaan path only).

## User request

A scoring method for lavaan CFAs (and the data.frame entry point) that returns
sum scores by averaging each factor's items: `fsL` = mean of loadings,
`fsb` = mean of intercepts, `fsT` computed from the CFA theta matrix. Must support
multiple dimensions via a user-specifiable item→sum assignment. No missing data:
error for now.

Method name (RESOLVED, user-confirmed): `"mean"` (distinct from regression/Bartlett
/ ML / EB; `ML` keeps meaning the Bartlett alias on the lavaan path).

Score constant (RESOLVED, user-confirmed): **(a) raw sum score**:
`fs_k = (1/|I_k|) Σ_{i∈I_k} y_i` (NO centering — the one `get_fs()` output whose
raw column is not model-mean-centered; document prominently),
`fsb = M ν + M λ α` (mean of model-implied item means; with α = 0: exactly the mean
of the item intercepts, as requested). The (score, fsb) pair is exactly
self-consistent: E[fs] = fsL α + fsb. Implementation: the "mean" branch of
`compute_fscore()` uses `a_mat = M` but skips centering (no `+ alpha` on the score,
`center_y` ignored) and substitutes the fsb formula `M ν + fsL α` in the shared tail
(fsL = M λ, fsT = M θ M', scoring_matrix = M unchanged).

Item→sum assignment (RESOLVED, user-confirmed): `sum_items`, a **named list,
factor → item names**, e.g.
`sum_items = list(ind60 = c("x1","x2","x3"), dem60 = c("y1","y2","y3","y4"))`.
`sum_items = NULL` auto-derives from the estimated loadings (an indicator is
assigned to the unique factor with a nonzero loading; >1 → error naming the item(s);
0 → the item is unused). Strict coverage: the assignment must cover every factor of
the model (a factor with no items → error); unknown factor/item names → error.
Subset sums are out of scope (they would change the stage-2 `fsL` column contract).

## Architect corrections (r-architect review 2026-08-17, applied to design)

- **C1 (orientation)**: `M` is **q × p (factor × item)**, row k holding `1/|I_k|`
  on factor k's items — the shared tail (`a_mat %*% y1c`, `a_mat %*% lambda`,
  `a_mat %*% theta %*% t(a_mat)`) requires q×p. `compute_a_mean()` builds
  `M <- matrix(0, ncol(lambda), nrow(lambda), dimnames = list(colnames(lambda),
  rownames(lambda)))`; `M[k, match(items, rownames(lambda))] <- 1/length(items)`.
- **C2**: `compute_fscore()` must forward `sum_items` into its
  `compute_a_from_mat(...)` call (else user assignment silently ignored).
- **C3 (corrected identity)**: the self-consistent pair is
  `fsb = E[fs] = M(ν + λ α) = fsL·α + M ν` (mean of item intercepts when α = 0;
  for no-mean models `ν` falls back to sample item means so `fsb` equals the
  score column mean exactly). Do NOT state "E[fs] = fsL α + fsb" (double counts
  M λα).
- **C4 (stage-2 caveat)**: with mean scores, stage-2 latent means are pinned to
  `fsb / L` (not stage-1 α); path coefficients are location-invariant so this
  does not affect 2S-PA paths; regression/Bartlett scores WITH fsb supplied
  instead pin a population-scale wrong mean — mean scores behave better in
  stage 2. Document this.
- **C5**: add a comment at the `do.call(compute_a_from_mat, ...)` in
  `compute_fspars()` (R/get_fscore_math.R:~327): the routing relies on the
  element NAMES of `mat[c("lambda","psi","theta")]`; do not drop them.
- **C6**: post-Plan-A test baseline is **1014** (not 707); `get_fs_lavaan()`
  signature is at `R/get_fscore.R:140`.
- **C7 (found in pre-check diff review)**: explicit `sum_items` was indexed
  POSITIONALLY (`sum_items[[k]]`) against model factor position k, so a user
  list named in a different order than the model factors silently swapped the
  per-factor scores and produced phantom `fsL` cross-loadings (no error; passed
  the green suite because all tests passed model order). FIX: at the end of the
  explicit `else` branch, normalize to model order before the positional loops:
  `sum_items <- sum_items[match(lv_names, names(sum_items))]`. Regression test
  (reversed-order list must equal auto-derived output at tol 0) added to
  test-get_fs_mean.R.
- **C8 (user directive, implemented 2026-08-18)**: for `method = "mean"`,
  `fsb` is the MEASUREMENT INTERCEPT of the score regressed on the uncentered
  latent = **the mean of the factor's item intercepts, `M ν`** — NOT
  `M(ν + λα) = E[fs]`. This unifies the mean method with the package's
  `fsb = E[fs] − fsL·α` convention (regression/Bartlett already use it) and lets
  stage-2 (`tspa()`) recover the endogenous latent mean instead of silently
  returning 0. `compute_fscore()` fsb line changed to `as.numeric(a_mat %*% nu)`
  for "mean" (α = 0 → identical to before, so no behavior change for models
  without a mean structure). @return wording + two meanstructure assertions
  updated to the `M ν` formula; new deterministic test with α = 3 pins
  `fsb ≠ E[fs]` (a plain CFA's latent mean is un-identified, so fitted fixtures
  all have α = 0 and could not pin this).

## Verified math (R 4.6.1 / lavaan session 2026-08-17)

Scoring matrix `M` (q × p, factor × item): row k holds `1/|I_k|` for the items
of factor k. With `L = M λ`, `T = M θ M'`:

- `fsL = M λ` — diagonal entry per factor = mean of its items' loadings
  (e.g. PoliticalDemocracy ind60: 1.672353 = mean of 1, 2.19, 1.82). Verified.
- `fsT = M θ M'` — for one factor: `(1/|I_k|²)` × sum of all elements of that
  item set's θ sub-block (i.e. variance of the item-mean residual incl. residual
  covariances); off-diagonal across factors picks up cross-factor residual
  covariances (zero in plain CFAs). Verified: single-factor θ-block sum/9 =
  0.073436.
- The existing `compute_fscore()` tail already computes
  `fsL = a λ`, `fsT = a θ a'`, `fsb = α − fsL α`,
  `scoring_matrix = a` for any `a` — so a "mean" `a`-builder plugs in with a
  branch in `compute_a_from_mat()` (`R/get_fscore_math.R:361-381`) and no other
  math code. For "mean", the branch additionally skips centering (raw score, no
  `+ α`) and the fsb formula is `M ν + fsL α` instead of `α − fsL α` (approved
  raw-score convention above).
- No-mean-structure CFA (all package examples): `est$alpha`/`est$nu` are NULL →
  `compute_fscore()` substitutes `α = 0`, `ν = colMeans(y)`.

## Design

- Add `"mean"` to the `match.arg` choice vectors of: `get_fs.lavaan()` (:253),
  `get_fs.data.frame()` (:100), `get_fs_lavaan()` (`R/get_fscore.R:140`),
  `normalize_fs_method()` (identity — no aliasing; note `ML` must keep mapping to
  Bartlett), `get_fs_blocks.lavaan()` (:137), and `compute_fscore()` /
  `compute_a_from_mat()` (`R/get_fscore_math.R`). `compute_fscore()` forwards
  `sum_items` into its `compute_a_from_mat(...)` call (C2).
- `compute_a_from_mat(method = "mean", lambda, theta, sum_items)` →
  `compute_a_mean(lambda, sum_items)`: build `M` (**q × p**) via
  `M[k, match(items, rownames(lambda))] <- 1/length(items)` (C1).
  `sum_items = NULL` → auto-derive from `lambda`: an indicator is assigned to the
  unique factor with a nonzero loading; 0 nonzero loadings → skip the item;
  >1 → `stop()` naming the item(s) and requesting explicit `sum_items`; a factor
  left with no items → `stop()`. Explicit input: must be a named list; unknown
  factor names → `stop()`; factors missing from the list → `stop()` (strict full
  coverage); unknown item names → `stop()`; an item in two sums → `stop()`; a
  factor with zero items → `stop()`. All stops use `call. = FALSE`.
- New parameter `sum_items` added to: `get_fs.lavaan()`, `get_fs.data.frame()`,
  `get_fs_lavaan()` (pass-through), `get_fs_blocks.lavaan()`, `prepare_fs()`,
  `compute_fscore()` (new arg, default NULL). One roxygen `@param` on the generic;
  `@rdname`/`@inheritParams` wiring so `man/get_fs.Rd` carries it exactly once.
  (Plan A already edited `get_fs()`'s `@param method`; this plan appends the
  "mean" sentence to the same paragraph — keep it readable as one block.)
- **Guard rails** (hard `stop()` with one clear message each, inserted in
  `get_fs.lavaan()` immediately after the `inherits` check — before
  `format <- match.arg()`, the prior validation, `correct_evfs`, the blocks
  builder, `vfsLT`, and the `reliability` machinery, which only accept
  regression/Bartlett). Missing-data gate:
  `!all(vapply(object@Data@Mp, is.null, logical(1)))` (only models fitted with
  missing data retained have a non-NULL pattern; listwise-deleted fits pass —
   correct semantics). Second gate: `corrected_fsT`/`vfsLT`/`reliability`/
   `prior_mean`/`prior_cov` non-default → single message listing the offending
   arguments ("not supported together with: ...").
- `augment_lav_predict()` is left at its two methods (it wraps `lavPredict()`,
  which has no mean scoring) — note this in its `@param method` if cheap.
- Multigroup: `M` is re-derived PER GROUP from each group's `est$lambda`
  (architect C-b: a loading fixed to zero in one group only would make a
  once-derived M silently wrong; per-group derivation also yields correct
  per-group auto-derivation errors). Explicit `sum_items` is model-level and
  used as-is; per-group `fsL = M λ_g` etc. fall out of the existing loop. The
  `prepare_fs` inner closure captures `sum_items` — no signature change.

## Files

- `R/get_fscore_math.R` — `compute_fscore()` (+ `@param method`, new arg),
  `compute_a_from_mat()` mean branch + co-located `compute_a_mean()` helper.
- `R/get_fs_methods.R` — `normalize_fs_method()`, `get_fs.data.frame()`,
  `get_fs.lavaan()` (guard rails + new arg), `get_fs_blocks.lavaan()` (match.arg,
  missing-data error, thread the assignment).
- `R/get_fscore.R` — `get_fs()` `@param method` paragraph (lavaan sentence +
  note the method is lavaan-only) + new `@param` for the assignment;
  `get_fs_lavaan()` signature.
- `tests/testthat/test-compute_fscore.R` — math-level tests (scores, M, fsL, fsT,
  fsb against hand matrices from `lavInspect`).
- `tests/testthat/test-get_fs_mean.R` — **new** API test file (keeps
  `test-get_fscore.R` from growing; named per repo convention).
- `man/` regenerated by `devtools::document()` (NAMESPACE unchanged — no new exports).

## Tests (new)

Math level (`test-compute_fscore.R`):
1. Single factor: raw scores (= item means), `fsL` = mean of loadings,
   `fsT` = θ-block sum/9, `fsb = M ν + M λ α` — all hand-computed from
   `lavInspect(..., "est")`.
2. Two factors (PoliticalDemocracy ind60/dem60): per-factor raw scores,
   block-diagonal fsT (zero off-diagonal), full `M λ`/`M θ M'` equality.
3. With mean structure (`ind60 ~ 1; x1~1; ...`): fsb = mean of item intercepts
   (α = 0 case) / full `M λ α + M ν` (α ≠ 0 case); scores stay raw (uncentered).

API level (`test-get_fs_mean.R`):
4. `get_fs(fit, method = "mean")` output column layout matches the other methods'
   shape (fs_*, _se, _by_, ev_/ecov_, attributes fsT/fsL/fsb/scoring_matrix).
5. Multigroup CFA (HolzingerSwineford1939, group = "school"): unified + list
   formats, per-group attributes.
6. `get_fs(data.frame, model = ..., method = "mean")` works end-to-end.
7. `get_fs_lavaan(fit, method = "mean")` wrapper.
8. Missing data (inject NAs) → informative error.
9. `corrected_fsT = TRUE` / `vfsLT = TRUE` / `reliability = TRUE` /
   `prior_mean` / `prior_cov` + mean → clean errors.
10. Cross-loading model without explicit assignment → error naming the ambiguous
    item; with explicit assignment → works.
11. Assignment validation: unknown item name / unknown factor name → error;
    strict full-coverage rule (a factor omitted from the assignment → error).
12. Smoke: `tspa(model = "dem60 ~ ind60", data = fs, fsT = attr(fs, "fsT"),
    fsL = attr(fs, "fsL"))` runs (name contract intact).

## Verification

- Lifecycle per AGENTS.md. DONE 2026-08-18: final suite **1211 pass / 0 fail /
  0 warn / 0 skip** (707 + 307 Plan A + 190 Plan B + 7 fsb-convention tests);
  `R CMD check` = 0 errors / 0 warnings / 1 pre-existing OpenMx NOTE for all
  non-vignette steps (full vignette build temporarily blocked by a GitHub
  rate-limit on an ECLS download in a pre-existing vignette). NAMESPACE diff
  empty; man/ regenerated by document() only.
- `devtools::check()` — expect the one pre-existing OpenMx NOTE only; roxygen
  warnings for the new `@param` must be 0 (new arg documented on every signature
  that carries it).

## Risks / notes

- Per the approved raw-score convention, the "mean" score is the only `get_fs()`
  output whose raw column is not model-mean-centered — document that prominently
  in `@param method` and the `@return`/details so it is not mistaken for
  regression-style scores.
- Stage-2 contract is unchanged: `tspa()` reads the same attributes; no tspa
  edits needed. The raw score / `fsb = M(ν + λ α)` pair is exactly
  self-consistent: `fsb = E[fs] = fsL·α + M ν` (C3 — the ONLY scoring method
  whose stored `fsb` equals the score's actual mean). Caveat (C4): with mean
  scores the stage-2 latent mean is pinned to `fsb / L`, not stage-1 `α`;
  paths are location-invariant (no 2S-PA impact).
- Quarantined code (`.quarantine/`) is never touched; `augment_lav_predict()`
  users see no change.
- Roxygen trap: `compute_fscore()`'s `@param method` is stale today ("Currently,
  only 'regression' is supported.") — this plan must fix it, not extend it
  blindly.

## Subagent roles (dispatched by the lead agent)

The lead agent implements and owns the lifecycle (Plan B runs only after
Plan A is fully checked). Specialists are dispatched via the task tool (types
`r-architect`, `r-doc`, `r-tester`); outputs are advisory — the lead reviews,
applies, and runs the lifecycle. Subagents do not hand-edit `man/`, do not run
`document()`/`test()`/`check()` themselves, and change nothing in `.quarantine/`.

1. **`r-architect`** — at plan start: brief with the verified math, the approved
   raw-score conventions (no centering; fsb = M ν + M λ α), the `sum_items`
   contract, and the file:line map; ask it to (a) recommend the cleanest
   decomposition of the "mean" branch — where in `compute_fscore()` /
   `compute_a_from_mat()` the centering skip and the fsb swap belong, keeping the
   shared tail (fsL/fsT/scoring_matrix) untouched; (b) review the
   `sum_items` threading through 6 signatures and the auto-derivation/validation
   rules in `compute_a_mean()`; (c) confirm the guard-rail gate placement in
   `get_fs.lavaan()` (before `correct_evfs`/`vcov_ld_evfs`/`compute_fsrel`, which
   only accept regression/Bartlett); (d) confirm `normalize_fs_method()` keeps
   `ML→Bartlett` and maps `"mean"` as identity; (e) re-verify the stage-2
   self-consistency claim against `tspa_schema_mf()` (`R/tspa.R:282-349`).
   Re-dispatch before `check()` for a diff review.
2. **`r-tester`** — after the code is drafted: implement the test lists above —
   math level in `test-compute_fscore.R`, API level in the NEW file
   `tests/testthat/test-get_fs_mean.R` — per repository conventions (edition 3;
   hand references from `lavInspect`; `expect_equal(..., tolerance = ...)`;
   `ignore_attr = TRUE` where attributes are intentional). Reference values to
   anchor on: PoliticalDemocracy single-factor `fsT` = 0.073436 (θ-block sum/9),
   mean loading 1.672353. Invite it to extend the validation/error-case matrix
   (auto-derivation ambiguity, strict coverage, missing data, wrapped
   `get_fs_lavaan()`, multigroup shapes).
3. **`r-doc`** — after code + tests, before `devtools::document()`: produce the
   exact roxygen wording for (a) `get_fs()`'s `@param method` — one coherent
   paragraph covering all five methods and per-class availability, disambiguating
   that `"ML"` = Bartlett alias for lavaan objects but prior-free OLS for
   `merMod`, and that `"mean"` is lavaan-only; (b) the new `@param sum_items`
   (named-list shape with example, NULL auto-derivation, strict full coverage,
   error conditions); (c) the `compute_fscore()` `@param method` fix (add
   "mean"; drop the stale "only 'regression'" sentence) + `@param center_y`
   note ("ignored for method 'mean'"); (d) `@return`/`@details` deltas (raw
   uncentered scores; per-factor fsb semantics; mean-only error cases);
   (e) note in `augment_lav_predict()`'s `@param method` that it intentionally
   stays at regression/Bartlett. Wiring: `sum_items` documented once, inherited
   into every `@rdname get_fs` method page. The lead applies verbatim to `R/*.R`,
   then runs `document()`.

Order: `r-architect` → (lead code in parallel with `r-tester` drafting tests) →
`r-doc` → lead runs `load_all() → document() → test() → check()`.

## Approved decisions (2026-08-17)

- Method string: `"mean"` (lavaan + data.frame + `get_fs_lavaan()`).
- Score: raw item mean, NO centering; `fsb = M(ν + λ α) = E[fs]` (mean of
  intercepts when α = 0); `fsL = M λ`; `fsT = M θ M'`; `scoring_matrix = M`
  (q × p); `M` re-derived per group from each group's loadings.
- Assignment: `sum_items` named list (factor → items); `NULL` = auto-derive from
  loadings (unique nonzero loading per indicator); strict full-factor coverage;
  validation errors for unknown names / ambiguity / uncovered factors /
  item-in-two-sums.
- No missing data: hard error. `corrected_fsT`/`vfsLT`/`reliability`/
  `prior_mean`/`prior_cov` + "mean": hard error.
- Test placement: math in `test-compute_fscore.R`, API in new
  `test-get_fs_mean.R`.
- Execution: sequential after Plan A, single lead agent, with the subagent
  dispatch above.
