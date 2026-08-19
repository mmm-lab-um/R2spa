# Plan: `fs_indiv()` (individual-specific fsL/fsT data frame) + latent `psi`/`alpha` attributes on `get_fs()`

## Context

Two requested capabilities, plus a redundancy the user asked to surface:

1. **`fs_indiv()`** — a function that converts a `get_fs()` result's `fsL`/`fsT` (and
   `fsb`) attributes into a **per-row data frame** (one row per observation for lavaan,
   one row per cluster for `merMod`) carrying **individual-specific** loadings, SEs, and
   error variances/covariances — analogous in *shape* to `mirt::fscores(full.scores.SE = TRUE)`.
2. **`get_fs()` latent-moment attributes** — also return the **effective (prior-adjusted)**
   latent means (`alpha`) and covariance (`psi`) as attributes.
3. **Reconcile redundancy** with the existing exported `augment_lav_predict()` (staying code,
   `R/get_fscore_math.R:116`), which already returns a per-observation SE/ev/ld data frame.

> Note: `augment_lav_predict()` is **staying** code (`R/get_fscore_math.R`), exported
> (`NAMESPACE:7`) and covered by staying tests; only its *consumer* (`tspa_mx.R`) + vignettes
> are quarantined. (The user referred to it as "in quarantine"; it is not.)

### Key technical findings (verified numerically on lavaan 0.x)

- For a **linear CFA**, an observation's `fsL`/`fsT`/`fsb` depend only on its
  **observed-indicator pattern**, not on response values. So "individual-specific"
  = **pattern-resolved** (rows differ only under missing data / across groups) — matches
  mirt's linear-model behavior. `get_fs()` **already** emits pattern-resolved per-row
  `_se`/`ev_`/`ecov_`/`<latent>_by_fs_<score>` columns (verified: per-pattern SE
  0.37829 / 0.40327 / 0.39417, each `== sqrt(diag fsT_pattern)`, 301 rows, `all.equal` vs
  `compute_lav_fs_matrices` per pattern).
- **SE value**: `get_fs() $fsT = AΘAᵀ` **equals** `augment_lav_predict()`'s
  `(I − VΨ⁻¹)V` exactly (both 0.1551/0.0255 on the HS1939 2-factor regression fit), where
  `V = lavPredict(acov=TRUE) == Vpost = (Ψ⁻¹ + ΛᵀΘ⁻¹Λ)⁻¹` (posterior variance). We keep
  `sqrt_or_na(diag(AΘAᵀ))` (R2spa's measurement-error SE = `get_fs()`'s own `_se` column).
  This is **not** mirt's EAP SE `sqrt(diag(Vpost))` (larger); mirt-EAP is documented as
  out-of-scope / default-off.
- **`augment_lav_predict()` overlaps Request 1 ~90%** (identical per-row values). Differences:
  it takes a **lavaan fit** (re-runs `lavPredict`), uses `se_*`-prefix + `int_*` columns +
  char-matrix `ld/ev/int` attributes (OpenMx path), and supports **only** regression/Bartlett
  (no `mean`/priors/`corrected_fsT`).

### Confirmed decisions (user)

- **API**: `fs_indiv(fs)` is a **standalone** exported function (keeps `get_fs()` contract
  unchanged; reusable on cached `fs`); `get_fs()` itself is not given a new argument.
- **Latent moments**: return the **effective / prior-adjusted** `psi`/`alpha`
  (a zero `alpha` vector when the model has no mean structure), **group-level**, shape mirroring
  `fsT`. **Point estimates only** (no sampling SEs of the latents).
- **Redundancy**: **Option 1** — `fs_indiv()` is the canonical per-row engine;
  **refactor `augment_lav_predict()` onto it** (thin wrapper preserving its public contract).
  `fs_indiv()` **must support `get_fs.merMod()`** (one row per cluster).
- **Name**: `fs_indiv()` (renameable; alternatives `fs_to_indiv()`, `fs_individual()`).

## 1. Feature 1 — `R/fs_indiv.R`: `fs_indiv(fs, include_intercept = FALSE, ...)`

**Output** — a single long `data.frame`, `nrow(== input get_fs() rows)`, reusing `get_fs()`
naming exactly:

`fs_<f>` · `fs_<f>_se` · `<latent_j>_by_fs_<f>` (q²) · `ev_fs_<f>` (q) ·
`ecov_fs_<a>_fs_<b>` (q choose 2, lower-tri, matching `augment_fs`) ·
[optional `int_fs_<f>` (q) when `include_intercept = TRUE`] · trailing `group` column
(unified MG, named via `attr(fs,"group_col")` or `"group"`) · subject/cluster-id column (merMod).

**Shared per-row engine (the core DRY).** Factor the per-block column generation out of
`augment_fs()` (`R/get_fscore.R:350`) and `augment_fs2()` (`R/get_fscore_math.R:10`) into one
value-only helper, e.g. `fs_row_cols(fs, fsL, fsT, fsb)` → `matrix(se, loadings, ev, ecov [, int])`.
All three of `augment_fs`, `augment_fs2`, and `fs_indiv` call it; each keeps its **own**
naming/layout (r2spa `*_se` vs legacy `se_*`/`int_*`/upper-tri), so values stay equal and
naming stays per-consumer.

**Row resolution** (per input shape):
- **lavaan `unified`/`list`** (`R/get_fscore.R:171,240,392` shapes): for each group `g`,
  read `label = attr(fs,"fs_pattern")[[g]]$label` and the group's `fsT`/`fsL`/`fsb` attributes.
  Attribute a **plain matrix/vector** (single/complete pattern) → every row uses it; a
  **named list by pattern label** → row `i` uses `attr[[label_i]]` (lookup by the
  `"+"`-joined pattern label). Emit pattern-specific values per row, place in original row order.
- **merMod** (`get_fs.merMod`, `R/get_fs_methods.R:660-741`): one row per cluster; per-cluster
  matrices come from the 2×2×n_clus arrays `attr(fs,"fsT")[,,j]` / `attr(fs,"fsL")[,,j]`
  (dim-3 = cluster id = `names(blocks)` = subject levels). Preserve `legacy_names` translation
  (reuse `rename_legacy_fs_cols`, `R/get_fs_methods.R:749`).
- Internal helper `resolve_fs_per_row(fs)` (co-located) drives the SG/MG/unified/list/merMod
  dispatch. Values are cross-checked against `get_fs()`'s own per-row columns in tests.

## 2. Feature 2 — latent `psi`/`alpha` attributes on `get_fs()`

- **Effective values** from the same `psi_use`/`alpha_use` already used in scoring
  (`R/get_fs_methods.R:148-150`): `psi_effective = prior_cov else est$psi`;
  `alpha_effective = prior_mean else est$alpha else rep(0, q)` (a named zero vector when the
  model has no mean structure — matches `compute_fscore`'s `alpha` default).
- **Group-level** (not per-pattern); `psi` q×q, `alpha` length q.
- **Shape mirrors `fsT`:** `unified` → named list by group label; `list` → per-group-df
  attribute + outer-list attribute. Attach **post-`assemble_fs_blocks`** in
  `get_fs.lavaan()` (next to the `vfsLT`/`reliability` attachments,
  `R/get_fs_methods.R:391-438`), branching on `format`. **`merMod`** (`get_fs.merMod`):
  `psi = VarCorr(object)[[1]]` **with dimnames renamed to `re_names`** (`u0`/`u1`/…) so it
  aligns with the `fsL` column names; `alpha = rep(0, q)` named by `re_names` (RE are mean-0).
- **Round-trip:** add `"psi"`, `"alpha"` to `attr_keys` in `fs_to_group_list()`
  (`R/get_fscore.R:241, 402`) so the attributes survive the unified↔list conversion.
- **Docs:** extend the roxygen `@return` on `get_fs()` with the two new attributes.
- **Non-breaking:** `tspa()` does not read these attrs (R/tspa.R:129); purely additive.

## 3. Reconcile `augment_lav_predict()` (Option 1)

- **Keep the public contract exactly**: lavaan-fit input, `drop_list_single`, legacy column
  layout (`fs_*`, `se_*`, `<ind>_by_fs_*`, `ev_*`/`ecov_*` **upper-tri** order, `int_*`), and
  the char-matrix `ld`/`ev`/`int` attributes (kept for the quarantined OpenMx path + vignettes).
- **Refactor the value source**: source per-pattern `fsL`/`fsT`/`fsb` from the **canonical**
  `compute_fscore`/`get_fs` blocks (the same engine `get_fs()` uses) instead of re-deriving via
  `lavPredict(acov=TRUE)` + `compute_lav_fs_matrices()`; emit through the **shared engine** (§1)
  with legacy naming + char-matrix attrs. Result: **one source of truth** for fsL/fsT/SE;
  `augment_lav_predict(fit)` ≡ `fs_indiv(get_fs(fit, method))` in **values**, differing only in
  naming / `int_*` / char-matrix attrs.
- **Guards:** staying tests `tests/testthat/test-get_fscore.R:1169-1221` (complete + missing,
  SG + MG, regression + Bartlett) must stay green. `compute_lav_fs_matrices()` remains a small
  tested helper (used by `test-get_fs_missing.R:72,138,182` and
  `test-lavPredict_equivalence.R:42,139,276`) — no longer the per-row *producer*.

## 4. Files touched

| File | Change |
| --- | --- |
| **new** `R/fs_indiv.R` | `fs_indiv()` (exported) + `resolve_fs_per_row()` + shared engine `fs_row_cols()` |
| `R/get_fscore.R` | `augment_fs()` → calls shared engine; `fs_to_group_list()` `attr_keys` += `psi`/`alpha` |
| `R/get_fscore_math.R` | `augment_lav_predict()` refactor onto shared engine + canonical fsL/fsT source; `augment_fs2()` → shared engine |
| `R/get_fs_methods.R` | `get_fs.lavaan()` + `get_fs.merMod()` attach `psi`/`alpha`; `get_fs()` roxygen `@return` |
| **new** `tests/testthat/test-fs_indiv.R` | (r-tester) |
| **new** `tests/testthat/test-get_fs_latent_moments.R` | (r-tester) |

`NAMESPACE` + `man/*` regenerated by `devtools::document()` (new `@export` for `fs_indiv`,
new `@param`/`@return`). **No `DESCRIPTION` change** — no new deps (base R + existing
`lavaan`/`lme4`/`MASS` only). `mirt` stays in `Suggests` and is used **only** as a conceptual
doc reference — never in `R/` code or tests (keeps check fast / no version coupling).

## 5. Execution model (subagents)

Charters from `~/.config/opencode/agents/`. Permissions: **@R-Architect** may edit
`R/` (logic + roxygen) and `DESCRIPTION`; **denies** `man/*`, `NAMESPACE`. **@R-Tester** may
edit **only** `tests/testthat/*`. **@R-Doc** may edit **only** `R/*` (roxygen). None
hand-edit `NAMESPACE`/`man/*`.

### @R-Architect (code + strategy + final gate)
- Implement `R/fs_indiv.R`: `fs_indiv()`, `resolve_fs_per_row()`, shared `fs_row_cols()`.
  (May carry a correct minimal roxygen header incl. `@export` so the function loads and
  `document()` works; roxygen *content* polished by @R-Doc next.)
- Implement the shared-engine refactor of `augment_fs()` / `augment_fs2()`.
- Implement `psi`/`alpha` attachment in `get_fs.lavaan()` / `get_fs.merMod()`; add
  `psi`/`alpha` to `fs_to_group_list()` `attr_keys`.
- Refactor `augment_lav_predict()` onto the shared engine + canonical fsL/fsT source
  (preserving its public contract).
- Guard the `Imports`/`Suggests` boundary; confirm no new dep.
- **After** delegating @R-Doc and @R-Tester, **verify both ran** (`document()` then `test()`),
  then run `devtools::check()` as the **final gate**: 0 errors / 0 warnings + the 1 existing
  OpenMx NOTE; treat any new WARNING/NOTE vs baseline as **blocking**.

### @R-Tester (tests/ only)
- **`tests/testthat/test-fs_indiv.R`** (mirrors new `R/fs_indiv.R`):
  - SG complete: n rows; per-row `fs_*_se == sqrt_or_na(diag(fsT))`; `ev`/`ecov`/`<ld>` equal
    the single `fsT`/`fsL`; matches `get_fs()`'s own per-row columns.
  - MG: per-group resolution; `group` column present and aligned.
  - **FIML missing data: different patterns → different per-row SE**; each row's SE
    `== sqrt(diag fsT_pattern)` (the individualized F4 behavior).
  - `method="mean"` + `prior_cov`/`prior_mean` + `corrected_fsT`: values flow through the attrs.
  - **merMod**: one row per cluster; row-j SE `== sqrt(diag(fsT_arr[,,j]))`; `legacy_names`
    variant; subject-id column.
  - `unified` vs `list` input → identical row values.
  - **Equivalence** `all.equal()` of `fs_indiv()` values vs `augment_lav_predict()` values
    (ignoring naming / `int_*` / char-attr deltas).
- **`tests/testthat/test-get_fs_latent_moments.R`**:
  - SG no-mean-structure: `alpha` attr `== rep(0,q)` named; `psi == est$psi`.
  - With priors: `psi == prior_cov`, `alpha == prior_mean`; `prior_cov`-only → `psi` = prior,
    `alpha` = `est$alpha`/zeros.
  - MG: per-group `psi`/`alpha == est[[g]]`; priors → all groups == prior; shape mirrors `fsT`.
  - merMod: `psi == VarCorr(object)[[1]]` (renamed to re_names); `alpha` zeros.
  - `fs_to_group_list()` round-trips `psi`/`alpha`.
- Numeric idiom per suite: `expect_equal(..., tolerance = <small>)`, `ignore_attr = TRUE`
  when comparing attribute-bearing objects, `expect_error()` for any validation. Prefer existing
  fixtures (`HolzingerSwineford1939`, `PoliticalDemocracy`, `sleepstudy`). Run
  `devtools::load_all()` **then** `devtools::test()` after writing/editing (and re-run after
  any architect/doc change lands). If a test reveals a real bug → report to @R-Architect, do not
  silently edit the test.

### @R-Doc (roxygen in R/ only)
- Roxygen for `fs_indiv()`: title, `@description`, `@param` (incl. `include_intercept`),
  `@return` (object + the per-row columns), `@examples` (existing datasets), `@export`.
- Update `get_fs()` `@return` to document the new `psi`/`alpha` attributes (incl. group nesting
  and the "group-level, effective/prior-adjusted, zero `alpha` when no mean structure" semantics).
- Update `augment_lav_predict()` roxygen to reflect the refactor (canonical source) while
  preserving its documented contract (input, `se_`/`int_` naming, char-matrix `ld/ev/int`).
- Preserve the repo's tag set/order; markdown is on (`Roxygen: list(markdown = TRUE)`).
- Run `devtools::document()` and **confirm only the intended `man/*.Rd` + `NAMESPACE` changed**;
  run `devtools::test()` to confirm nothing broke; hand the final `check()` to @R-Architect.

### Sequencing (respecting dependencies + edit scopes)
1. **@R-Architect** — all `R/` logic + roxygen skeleton (incl. `@export fs_indiv`).
2. **@R-Doc** — polish/own roxygen blocks → `devtools::document()`.
3. **@R-Tester** — write the two test files → `load_all()` + `test()`.
4. **@R-Architect** — verify doc+test ran; final `devtools::check()` gate.
Steps 2 and 3 touch disjoint scopes (`R/` roxygen vs `tests/`) and **may run in parallel**, but
default to sequential (2 → 3) to avoid interleaving `document()`/`test()` on divergent state.

### Lifecycle (unchanged, no reordering)
`devtools::load_all()` → **`devtools::document()`** (mandatory: new `@export`, `@param`,
`@return`) → `devtools::test()` → `devtools::check()` (expected: 0 err / 0 warn + the 1
existing OpenMx NOTE).

## 6. Risks / notes
- **Numeric equivalence**: the `compute_lav_fs_matrices →` canonical-source swap is safe because
  per-pattern `fsT` from `compute_fscore` is already validated `== compute_lav_fs_matrices(acov,…)`.
  Guard with the new equivalence test + existing staying tests.
- **Naming parity**: shared engine is value-only; each caller applies its own names
  (`get_fs`/`fs_indiv` → `*_se`; `augment_lav_predict` → `se_*`/`int_*`/upper-tri).
- **merMod latent naming**: `psi` dimnames renamed to `re_names` to align with the score/`fsL`
  columns (documented choice, low risk).
- **No new deps**; `mirt` never in code/tests (Suggests, doc-only reference).

## 7. Open minor decisions (confirm during implementation, non-blocking)
- Public name `fs_indiv()` (vs `fs_to_indiv()` / `fs_individual()`).
- `fs_indiv` **keeps** a `group` column (MG) and a subject-id column (merMod) in the long table.
- Latent-moment tests in a **new** `test-get_fs_latent_moments.R` (vs extending
  `test-get_fs_priors.R`).
