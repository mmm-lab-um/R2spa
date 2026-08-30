# PLAN 08 — mirt Item-Response `get_fs()` support (`SingleGroupClass`)

**Status: COMPLETE** (2026-08-19) — implemented in this commit; `mirt` remains a
`Suggests`-only dependency (namespaced `mirt::` calls, `require_mirt()` guard).
Closes follow-up **F2** (mirt S3 method, tracked from PLAN 01 / `STATUS.md`).
Builds on PLAN 07's `fs_indiv()` (the consumer of the per-row `mirt` attrs).

## 1. Scope
Give `get_fs()` a `mirt` backend so a fitted **single-group** IRT model
(`mirt::mirt` → `SingleGroupClass`) yields `fs_<factor>` scores with
**per-observation** loadings / error variances / SEs, and — because `fs_indiv()`
drives on the row-level `fsL`/`fsT` — so `fs_indiv(get_fs(<mirt fit>))` works end
to end with truly row-specific values (a different contract than the
lavaan/merMod one-block-per-group contract). Multi-group mirt
(`MultipleGroupClass`) is **out of scope**: the method errors with a "fit/extract
the single group first" message.

## 2. Key modeling decision (per-row regression loading, NOT the identity)
The implied loading for observation `i` is the **regression/shrinkage** form
`fsL_i = I − Vpost_i·Ψ⁻¹`, **not** the identity. This reuses the package's
canonical engine `R/get_fscore_math.R::compute_lav_fs_matrices(acov = Vpost_i,
psi, alpha, "regression")` row by row, yielding (per row):

- `fsL_i = I − Vpost_i·Ψ⁻¹`  →  unidimensional (unit latent variance): **`1 − SE²`**
- `fsT_i = fsL_i · Vpost_i`     →  unidimensional: **`(1 − SE²)·SE²`** (error var)
- `fsb_i = alpha − fsL_i·alpha` →  zero here (mirt `alpha = 0`)

This matches the **merMod EB** path (also shrinkage: `fsL_j = DKz −
DKz·inv_W·DKz`, ≠ `I`) and deliberately differs from merMod `ML`/`Bartlett`
(both `fsL = I`, Bartlett `fsT = SE²` by construction). Confirmed against the
hand formula to machine precision.

## 3. mirt data plumbing (verified empirically on `HolzingerSwineford1939`)
- Fit: `mirt::mirt(dat, n)` — **do NOT pass `itemtype`** (default is the
  single-itemtype recycle error: *"itemtype must have length 1 or be the same
  length as items"*).
- EAP means + SEs: `fscores(object, full.scores = TRUE, full.scores.SE = TRUE)`;
  score columns are the **factor names** from `mirt::extract.mirt(object,
  "factorNames")` (e.g. `F1`), SE columns `<f>.SE`.
- Per-observation posterior covariance: `fscores(object, full.scores = TRUE,
  return.acov = TRUE)` → a length-n **list of `q × q` matrices `Vpost_i`**;
  `diag(Vpost_i)` = the per-factor SE² for row `i`; off-diagonals = latent error
  covariances.
- Latent prior: `psi = diag(q)` (unit variance), `alpha = rep(0, q)` — mirt's
  **default** factor prior. `prior_cov`/`prior_mean` are not threaded here (no
  user-facing priors for the mirt path yet). ⚠ `object@Model$Theta` is the
  **quadrature NODE grid**, not the factor covariance — do not use for `psi`.

## 4. `get_fs()` method contract
`get_fs.SingleGroupClass(object, format = c("unified", "list"), ...)`:
- `fs_<f>` columns = EAP means (named by the mirt factor names).
- Per-observation attrs: `fsT` / `fsL` = **lists** (length `n`, each a `q × q`
  with rownames/colnames = factor names); `fsb` = named zero vector length `q`.
- `fs_pattern = list(label = seq_len(n), pat = NULL)`; marker
  `attr(out, "mirt_per_obs") = TRUE` (set only here — the `fs_indiv()`
  dispatch key); `scoring_matrix` not attached.
- `psi` = `diag(q)` (rownames/colnames = factor names), `alpha` = named zero
  vector.
- `MultipleGroupClass`: `stop()` with
  "Multi-group mirt models are not supported by get_fs() yet; fit/extract the
  single-group model (object@Model) on the resulting SingleGroupClass object."

## 5. `fs_indiv()` per-row branch
`R/fs_indiv.R::resolve_fs_per_row()` checks the `mirt_per_obs` marker **first**
and routes to a new `resolve_per_obs()`, which mints **one block per row**
(each the row's own `fsL`/`fsT`, plus an all-NA `fsb` for all-missing rows). The
shared `fs_indiv()` block loop then emits per-row values. **`fs_indiv()` itself
is byte-for-byte unchanged** apart from this dispatch branch (the branch is
dormant unless the marker is set, so the lavaan/merMod paths are untouched).
Missing data: a row with no scorable indicator gets an all-NA block from
`resolve_group_blocks()` → NA score/SE/ev/intercept (same convention as lavaan).

## 6. Verification
- `tests/testthat/test-get_fs_mirt.R` — 72 expectations: unidimensional
  `F1_by_fs_F1 == 1 − SE²` (1e-16); `SE == sqrt(1−SE²)` identity; `q = 2`
  off-diagonal `fsT`/`fsL` vs hand formula + mirt SE identity; structural
  contracts (per-row lists, `fs_pattern`, marker, `psi`/`alpha`); `fs_indiv()`
  1:1 match (score/se/ev, 2-factor loadings, `nrow`, column set, NA rows);
  multi-factor (2 factors); `MultipleGroupClass` error; S3 dispatch of an S4
  object; `require_mirt()` guard.
- Independence: the tree **without** the mirt method still passes its full
  3104 (bisect-checked at the PLAN 07 commit); adding mirt → **3176 / 0 / 0 / 0**.
- `devtools::document()` run with pinned **roxygen2 8.1.0** (the sandbox's
  global 7.3.1 churns `RoxygenNote` + `importFrom` splitting). `NAMESPACE`:
  +`S3method(get_fs, SingleGroupClass)` / `+S3method(get_fs, MultipleGroupClass)`;
  `man/get_fs.Rd` gains the mirt aliases/usage/@description line.
- `devtools::check()` (`--no-manual`; PDF-manual build errors are the
  environmental no-`pdflatex`, not a package issue): **0 errors / 0 warnings /
  1 NOTE** — the expected OpenMx `Imports` baseline item.

## 7. Out of scope (follow-ups logged in `STATUS.md`)
- **F5** — expose the alternative mirt per-row SE `sqrt(diag(Vpost_i))`
  (regression-SE convention) alongside the `fscores(..., full.scores.SE = TRUE)`
  value.
- **F6** — multi-factor mirt: `psi` is the unit matrix under mirt's default unit
  prior (not a parameter); thread a non-unit covariance through if mirt ever
  exposes one (mirrors lavaan `prior_cov`).
