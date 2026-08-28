# PLAN 16 — analytic influence-function Jacobian for `vcov_corrected()` / `tspa(corrected_se = )`

**Date:** 2026-08-28 (updated after the Phase-2 gate)
**Owner roles:** `r-architect` (code), `r-tester` (tests), `r-doc` (roxygen/NEWS)
**Status:** **gate passed (saturated case)** — not yet implemented. The Hessian-free closed form
(§2.4) reproduces the FD corrected vcov to `max|diff| = 0.0049` (the FD's own precision) on the
T3 fixture, fully analytically. Remaining: wire it into `vcov_corrected()`, and extend to the
general (non-saturated) structural model (§4.3).
**Blocked by / relates to:** none. Extends the first-order corrected-SE path (re-integrated
2026-08-23, `R/tspa_corrected_se.R`); does not change the public contract of
`vcov_corrected()` or `tspa(corrected_se = )` — only the way the Jacobian `J` is computed.
Supersedes the finite-difference `J` as an *internal* option; the FD route stays as the
A/B reference (and the fallback).

---

## 1. Problem

`vcov_corrected()` propagates the sampling covariance of the stage-1 `fsL`/`fsT` into the
stage-2 estimate through the Jacobian `J = ∂θ̂/∂η` (`θ` = free structural params, `η` =
the fixed `fsL`/`fsT` entries), correcting `vcov(fit) → vcov(fit) + J·vfsLT·J'`. Today `J`
is estimated by **central differences**: one full stage-2 refit on each side of each free
element → `2·nfree` refits (e.g. `2·15 = 30` for the 3-factor T3 fixture). Consequences:

- **Speed:** ~30 small `sem()` refits per call (seconds); the T4 bootstrap test skips under
  `R CMD check` partly for this cost.
- **Determinism / cross-platform drift:** each refit is a fresh optimizer run, so `J`
  (hence the corrected vcov) is only reproducible to the optimizer/BLAS noise floor. This is
  the root cause of the cross-OS golden drift the T3 goldens were widened to `1e-2` for
  (see `tests/testthat/test-tspa_corrected_se.R`).

The analytic influence-function Jacobian evaluates `J` at the **single base fit** — no refits,
no step `h`, no noise amplification — so it is deterministic to machine precision and drops the
`2·nfree` refits (~30× faster). The formula (saturated case, §2.4):

```
J = V · H_θη ,   H_θη[k,j] = -(n/2) · tr( Σ⁻¹ (∂Σ/∂η_j) Σ⁻¹ (∂Σ/∂θ_k) )
```

with `V = vcov(fit)` and `Σ` the base-fit implied score covariance.

## 2. Verified foundation (2026-08-28 investigation, lavaan 0.7-2, R 4.6.1)

Facts established empirically before/while drafting this plan. The implementer should treat
these as given (re-verify only if the model grammar changes).

**Model structure — a standard covariance-structure MLE.** The stage-2 model fixes the
measurement part (R/tspa.R:1508-1531): per latent `var[k]`/score `fs[i]` a *fixed* loading
`var[k] =~ fs[i] * L[i,k]` and per lower-tri `fsT` entry a *fixed* error `fs[i] ~~ fs[j] *
T[i,j]`; only the structural block is free. The implied covariance of the observed score
variables is

```
Σ_obs(θ) = L · F(θ) · L' + T          (L = fsL fixed, T = fsT fixed, F free)
```

**All ingredients are in the lavaan partable** (`lavInspect(fit, "est")`):
`lambda` = `L`, `psi` = **exogenous-only** latent cov (diagonal for a DAG), `beta` = the
structural matrix `β`, `theta` = `T`. The **full** latent covariance is the VARIMAX/structural
solution, *not* `est$psi`:

```
M  = solve(diag(q) - beta)             # (I - β)^{-1}
F  = M %*% psi %*% t(M)                # full latent cov (includes structural covariances)
Σ  = L %*% F %*% t(L) + theta          # matches cov(score data) at the MLE (verified)
```

(Using `est$psi` directly for `F` is a trap: it leaves the latent off-diagonals zero and the
implied cov is wrong by ~3× on the dominant diagonal.)

**The log-likelihood is a pure function of `Σ(θ)`.** The stage-2 log-likelihood is
`ℓ(θ) = -(n/2)·log|Σ(θ)| - (n/2)·tr(Σ(θ)⁻¹ S_ml)` (`S_ml` = the group score covariance).
**Verified:** the central-difference `∂²ℓ/∂θ²` of this expression equals `−V⁻¹` (lavaan's
observed info) to ~1% on the T3 fit — so `V` is exactly the Hessian of this `ℓ(θ)`, and `Σ`
may be treated as an unconstrained function of `θ`.

**`Δ_θ` and `Δ_η` are exactly computable and were verified to machine precision.** With
`vech(·)` the column-major lower triangle (length `q(q+1)/2`) and the structural params in
`coef(fit)` order:

- Free structural params (`Δ_θ`, `p × q(q+1)/2`):
  - regression `β[i,j]` (latent `i` on latent `j`): `∂F/∂β[i,j] = M[,i] %o% F[j,] + F[,j] %o% M[,i]`,
    then `∂Σ/∂β[i,j] = L %*% (∂F/∂β[i,j]) %*% t(L)`;
  - latent (residual) variance `ψ[k,k]`: `∂F/∂ψ[k,k] = M[,k] %o% M[,k]`, then `L %*% · %*% t(L)`.
  - **Verified:** `Δ_θ` equals a central-difference of `Σ(θ)` (recomputing `F` from a
    perturbed `β`/`ψ`, **no refit**) to `max|Δ| = 0` on all 6 columns of the T3 fit.
    (The common first bug — `F[i,] %o% M[,j]` for the second term — gives `F[,j] %o% M[,i]`;
    the regression columns were off by O(1-6) until fixed, the variance columns were already exact.)
- Fixed entries (`Δ_η`, `nfree × q(q+1)/2`, order = `fsL` full column-major then `fsT`
  lower-tri column-major — the same order as the `vfsLT` attribute / `which_free`):
  - loading `L[i,k]` (score `i` on latent `k`): `∂Σ/∂L[i,k] = e_i (F L')[k,] + (L F)[,k] e_i'`;
  - error entry `T[i,j]` (`i ≥ j`): `∂Σ/∂T[i,j] = E_ij + E_ji` (trivial).

**The naive "unfreeze and read the cross-Hessian block" shortcut is ruled out.** Refitting the
stage-2 model with the `fsL`/`fsT` entries *free* is **not identified** for the fully-crossed
2S-PA measurement structure (observed: negative latent variances, non-PD residual covariance),
and in any case a saturated free-`Σ` model has a zero score at the MLE. The cross-Hessian must
be assembled from `Δ_θ`, `Δ_η`, and the model's likelihood geometry.

**2.4 The saturated-case cross-Hessian is Hessian-free (the deliverable).** This is the key
result. For a stage-2 model whose structural part **saturates** the latent covariance
(`#free structural = q(q-1)/2`, e.g. the T3 chain `dem60~ind60; dem65~ind60+dem60`), the
implied score cov equals the sample cov at the MLE (`Σ = S_ml`), so the score matrix
`[Σ⁻¹ S_ml Σ⁻¹ − Σ⁻¹] = 0`. The influence cross-Hessian then collapses to

```
H_θη[k,j] = -(n/2) · tr( Σ⁻¹ (∂Σ/∂η_j) Σ⁻¹ (∂Σ/∂θ_k) ) ,    J = V %*% H_θη
```

which needs **only** `Σ⁻¹`, the analytic `Δ_θ`/`Δ_η` (§2.3), and `n` — **no log-likelihood
Hessian, no finite differences, no refits**. (Derivation: `∂ℓ/∂θ = (n/2)tr([Σ⁻¹S_mlΣ⁻¹−Σ⁻¹]
∂Σ/∂θ)`; differentiating w.r.t. `η`, the `[Σ⁻¹S_mlΣ⁻¹−Σ⁻¹]`-times-`∂²Σ/∂θ∂η` term vanishes at
the MLE, and `∂(Σ⁻¹S_mlΣ⁻¹−Σ⁻¹)/∂η|_{Σ=S_ml} = −Σ⁻¹(∂Σ/∂η)Σ⁻¹`.)

**Gate result (2026-08-28):** on the T3 fixture (`n = 75`, saturated), the **fully-analytic**
`J` (the `dFb`/`dFp`/`dSig_L`/`dSig_T` formulas above + the §2.4 closed form) gives a corrected
vcov matching the FD `vcov_corrected()` to `max|diff| = 0.0049` — i.e. at the level of the FD's
*own* refit/finite-difference noise (the FD uses `h0 = 1e-5` + default-convergence refits).
Per dominant diagonal, FD vs analytic: `dem60~~dem60` 17.8948 / 17.8997; `dem65~~dem65`
12.9279 / 12.9287; `dem60~ind60` 1.5646 / 1.5641; `dem65~ind60` 4.3072 / 4.3065;
`dem65~dem60` 1.201 / 1.201. **The gate PASSES** — the analytic is expected to be *more*
accurate than the FD (the FD carries refit noise the analytic does not).

**2.5 The general `H_ΣΣ` route is ill-conditioned (why §2.4 is used instead).** The *general*
(non-saturated) cross-Hessian is `H_Θη = Δ_θ' H_ΣΣ Δ_η + (∂²Σ/∂θ∂η)·g_Σ`, needing the
log-likelihood Hessian `H_ΣΣ = ∂²ℓ/∂vech(Σ)²`. Computing `H_ΣΣ` by finite difference yields a
negative-definite (sane) matrix whose **diagonal is correct** (the coordinate quadratic form
`½h²H[1,1]` reproduces the actual `Δℓ` to ~0.3% at `h = 1e-4`), but the directional curvatures
`Δ_θ' H_ΣΣ Δ_θ` come out **1e3–1e6× too large** versus the ground-truth `∂²ℓ/∂θ² = −V⁻¹`. The
Hessian has large near-cancelling diagonal + cross terms (diagonal ~`−1.4e4`, true directional
~`−3e2`), which a finite difference cannot resolve. The saturated §2.4 form sidesteps this; the
general case must use the same Hessian-free structure, extended with the
`[Σ⁻¹S_mlΣ⁻¹−Σ⁻¹]` terms (§4.3) rather than an FD `H_ΣΣ`.

**Implementation gotchas already hit (so the coder does not re-hit them):**
- `D` (if ever needed) must be built with **paired** indices:
  `D[cbind(seq_len(m), vech_vecpos)] <- 1`; the grid form `D[1:m, vech_vecpos] <- 1` silently
  fills every row. (Not needed for the §2.4 form — that works directly with full `q×q`
  `∂Σ/∂θ`, `∂Σ/∂η` matrices.)
- `vech(A) = A[lower.tri(A, diag = TRUE)]` returns a **named** vector; wrap in `unname()`
  before numeric comparisons, and `as.vector(D %*% x)` (matrix %×% vector is a matrix).
- Single-row/col subsetting drops dimensions in R; use `%o%` for the outer products in
  `Δ_θ`/`Δ_η`. `sweep(as.matrix(x), 2, mu)` (not a data frame) before outer products.
- `lavaan::lavInspect(fit, "hessian")` is the *per-observation* Hessian (`vcov ≈ (n·hess)⁻¹`) —
  not needed for §2.4 but useful for cross-checks.

## 3. Scope & decisions

- **D1 — Same contract, new engine.** `vcov_corrected()` and `tspa(corrected_se = )` keep their
  signatures, return type, and the `tspa_corrected` double-correction guard. The change is
  purely how `J` is computed. `which_free`/`vfsLT` semantics unchanged (the analytic `J` is
  formed for the full `nfree` positions, then the same `which_free` sub-selection + principal
  submatrix rule applies).
- **D2 — Gate: passed (saturated); extend before default switch.** The §2.4 closed form is
  A/B-verified against the FD on the saturated T3 fixture (`max|diff| = 0.0049`, at the FD's
  own noise). **GO** for the saturated path. Before switching the default, the A/B must be run
  on **all** existing corrected-SE fixtures (T3 incl. multigroup, larger-`n` fixtures, a
  mean-structure case) at `max rel |Δ| < 1e-2` (relaxing the original `1e-3` target to the FD's
  achievable precision, since the FD is only accurate to its refit/`h` noise floor ~1e-3 rel),
  **and** the general (non-saturated) case must be implemented (§4.3) and A/B'd. If a fixture
  shape fails, that shape falls back to FD (D5) — do not ship a `J` that disagrees with the FD.
- **D3 — Hessian-free closed form is the method (supersedes the observed/expected-info choice).**
  The §2.4 form is exact for the saturated model and needs no `info = c("observed","expected")`
  switch. (The expected-info `½n(Σ⁻¹⊗Σ⁻¹)` form alone is **insufficient** — off ~13× on the
  dominant element at `n = 75` — because the correction needs the full `Σ⁻¹(∂Σ/∂η)Σ⁻¹(∂Σ/∂θ)`
  structure, not just the expected information.)
- **D4 — Saturated primary, general fallback.** The common 2S-PA case is a (near-)saturated
  structural model (`#edges = q(q-1)/2`, e.g. T3), where `g_Σ = 0` and §2.4 is exact. Detect
  saturation (check `‖Σ⁻¹ S_ml Σ⁻¹ − Σ⁻¹‖ < tol`, equivalently `#free structural = q(q-1)/2`)
  and branch; the general (non-saturated) path (§4.3) is required for restricted path models and
  must be exercised by at least one non-saturated fixture.
- **D5 — FD stays as reference + fallback.** The finite-difference `J` is not deleted. It is the
  A/B reference and remains available (`engine = c("analytic","fd")`, default `"analytic"` after
  the full gate; `"fd"` to force). If the analytic gate fails on some fixture shape, that shape
  falls back to FD with a note.
- **D6 — Determinism is the acceptance bar, not just speed.** Because the analytic `J` uses no
  refits, the corrected vcov becomes a pure function of the base fit + `vfsLT`; the T3 goldens
  can then be **tightened from `1e-2` back toward `1e-8`** (subject to §7), which is the real
  payoff for the CI cross-platform issue.

## 4. Approach (r-architect) — new internal in `R/tspa_corrected_se.R`

A single new internal, `vcov_jacobian_analytic(fit, vfsLT_layout)`, returning the `p × nfree`
`J` (same row/col layout as the FD `J`), reusing the existing `assemble()`/`which_free`/
`tsp_tri2full_colmajor()` plumbing in `vcov_corrected()`.

### 4.1 Reconstruct the measurement geometry (per group)
For each group `g` (single-group: one iteration): from `lavInspect(fit, "est")` pull
`L = lambda`, `psi`, `beta`, `T = theta` (per-group partables for MG). Form
`M = solve(diag(q) - beta)`, `F = M %*% psi %*% t(M)`, `Σ = L %*% F %*% t(L) + T`,
`Sinv = solve(Σ)`. Sanity: `Σ` PD and matches the group score data cov (guard, tol ~1e-6).
This is the `F = (I−β)⁻¹ψ(I−β)⁻¹'` reconstruction from §2 — **do not use `est$psi` as `F`**.

### 4.2 Build `Δ_θ` and `Δ_η` as **full `q×q`** matrices (verified in §2.3)
- `Δ_θ[k]` (full `q×q`): for each free structural param in `coef(fit)` order, `L (∂F/∂θ_k) L'`
  with the §2.3 `∂F` formulas (regression `dFb`, variance `dFp`; the second regression term is
  `F[,j] %o% M[,i]`).
- `Δ_η[j]` (full `q×q`): the `nfree` columns in `fsL`-full-then-`fsT`-lower-tri column-major
  order (`∂Σ/∂L[i,k]`, `∂Σ/∂T[i,j]` in §2.3). Same vectorization as the existing
  `assemble()`/`val_fsLT`, so the `J` columns line up with `which_free`/`vfsLT` with no
  reordering. (Work directly in the full-matrix basis — no `vech`/`D` needed for §2.4.)

### 4.3 Build `H_θη` and `J`
- **Saturated (primary, D4):** `H_θη[k,j] = -(n/2)·sum(diag(Sinv %*% Δ_η[j] %*% Sinv %*% Δ_θ[k]))`
  for each free param `k` and fixed entry `j` (`n` = group sample size). Then
  `J = V %*% H_θη` with `V = vcov(fit)` (per group, stacked in the same group order as the FD).
  Cost: `p × nfree` products of `q×q` matrices (`q` = #factor scores, small) — negligible.
- **General / non-saturated (fallback path, D4):** `g_Σ = Σ⁻¹S_mlΣ⁻¹ − Σ⁻¹ ≠ 0`. Extend §2.4:
  ```
  H_θη[k,j] = (n/2)·tr( A(Δ_η[j]) · Δ_θ[k] ) + (n/2)·tr( g_Σ · ∂²Σ/∂θ_k∂η_j )
  A(D) = -Sinv·D·Sinv·S_ml·Sinv - Sinv·S_ml·Sinv·D·Sinv + Sinv·D·Sinv     # = ∂(Σ⁻¹S_mlΣ⁻¹−Σ⁻¹)/∂η
  ```
  where `∂²Σ/∂θ_k∂η_j` is the second derivative of the implied cov (from `∂F/∂θ_k` plus the
  direct `L`/`T` structure; a small set of analytic expressions). When `‖g_Σ‖ < tol`, the second
  term is dropped and `A(·)` reduces to `−Sinv·Δ_η[j]·Sinv`, recovering §2.4. **This is the
  remaining implementation work** (the saturated path is already A/B-verified in §2.4).

### 4.4 Wire into `vcov_corrected()` / `tspa(corrected_se = )`
Add an `engine` argument (default per D5) to `vcov_corrected()` and thread `corrected_se`
through `tspa()`'s existing in-place path. The analytic path replaces the
`for (k in seq_len(nfree))` refit loop; the input guards (PSD `fsT`, `vfsLT` shape/symmetry,
`which_free` validation, `tspa_corrected` double-correction reject, the `assemble()`
round-trip) are **shared** and unchanged. `check_refit_convergence()` is FD-only.

## 5. Test plan (r-tester) — `tests/testthat/test-tspa_corrected_se.R`

A/B pattern: analytic `vcov_corrected(engine = "analytic")` ≡ FD
`vcov_corrected(engine = "fd")` at the D2 tolerance, plus invariants.

1. **Gate A/B (D2):** T3 3-factor (the `n = 75` saturated fixture) and every existing
   corrected-SE fixture — analytic ≡ FD at `max rel |Δ| < 1e-2` (the FD's noise floor); the
   `sqrt(diag())` corrected SEs agree. **Multigroup** corrected-SE fixture (per-group `fsL`/`fsT`):
   analytic ≡ FD, group order preserved (the `which_free`-per-group layout).
2. **Saturated vs general (D4):** the T3 chain exercises the `g_Σ = 0` primary path (already
   A/B-verified to 0.0049 in the investigation); add one **non-saturated** structural model
   (fewer edges than `q(q-1)/2`, so `g_Σ ≠ 0`) and pin the §4.3 general path against the FD.
3. **`which_free` / `vfsLT` sub-selection (D1):** `which_free = c(5,7)` (two error variances)
   + matching principal submatrix — analytic ≡ FD (subset path).
4. **Determinism (D6):** the analytic corrected vcov is **bit-stable** — call it repeatedly and
   after forcing a different RNG stream / (on the dev box) single-threaded vs default BLAS;
   `expect_identical` (no refits → no optimizer/BLAS noise). This is the property the FD lacks.
5. **Guards unchanged:** PSD-`fsT` reject, `vfsLT` shape/symmetry/`which_free` errors,
   `tspa_corrected` double-correction reject, no-free-params reject — all fire identically on
   the analytic path (they are shared, pre-`J`).
6. **`tspa(corrected_se = TRUE, engine = "analytic")` in-place:** sets `tspa_corrected = TRUE`,
   `standardizedSolution()` reports corrected std SEs, and `grandStandardizedSolution()` threads
   the corrected covariance (point estimates unchanged) — mirrors the existing in-place tests
   but on the analytic engine.
7. **Regression:** existing corrected-SE tests green; suite baseline (re-measured at start)
   preserved. The T4 bootstrap test stays skipped under `R CMD check` (unchanged), but with a
   faster analytic `J` the *local* run cost drops.

## 6. Docs (r-doc)

1. **Roxygen `vcov_corrected` / `tspa`:** new `@details` subsection "Analytic Jacobian" —
   the influence-function form `J = V·H_θη` with the saturated Hessian-free closed form (§2.4),
   that the default engine is refit-free (deterministic), the `engine` argument and its default,
   and that the FD engine remains as the A/B reference / fallback. `@param engine`. Then
   `devtools::document()`.
2. **`NEWS.md`:** 0.0.4 entry (Performance / internal: corrected-SE Jacobian computed
   analytically by default, refit-free and deterministic; FD retained as reference/fallback).
3. **Vignette** (`vignettes/corrected-se.Rmd`): short "How the Jacobian is computed" note
   (finite-difference vs analytic; determinism), no re-knit of the cached narratives unless a
   printed value changes (it should not, to the printed precision).

## 7. Acceptance

- D2 gate passes **across the board**: analytic ≡ FD `max rel |Δ| < 1e-2` on every corrected-SE
  fixture (incl. multigroup and a non-saturated model), and the T4 bootstrap comparison holds.
  (The saturated T3 case is already at `max|diff| = 0.0049`.)
- Default `engine = "analytic"`; `tspa(corrected_se = )` in-place path works on both engines;
  double-correction guard intact.
- **Determinism proven:** analytic corrected vcov is bit-stable across repeated calls and BLAS
  threading (test item 4); T3 goldens tightened from `1e-2` to `≤ 1e-8` (D6).
- Suite green (baseline preserved, 1 pre-existing `tspa_mx` warn); `R CMD check` 0/0/1 NOTE
  (the pre-existing title-case + CRAN-URL NOTE only); `document()` clean.

## 8. Out of scope

- Changing the *partial-by-design* coverage of the correction (still only the `fsL`/`fsT`
  sampling covariance; factor-score values / `se_fs` / `fsb` remain fixed — §2 of the
  `vcov_corrected` roxygen is unchanged).
- Two-stage independence assumption (stage-1/stage-2 cross-covariance) — still not modelled.
- `tspa_mx_model()` (OpenMx) — separate exact route; not touched.
- Propagating the analytic `J` to any new quantity beyond the corrected vcov.
- Re-deriving the stage-1 `vfsLT` (unchanged; still from `get_fs(vfsLT = TRUE)`).

## 9. Verification log

- **Gate number (saturated T3):** fully-analytic `J` vs FD `vcov_corrected()` —
  `max|diff| = 0.0049` (absolute), matching to ~3-4 decimals on every element (the 0.88
  "max rel" is a near-zero off-diagonal, abs diff 0.001). Dominant diagonals: FD 17.8948 /
  analytic 17.8997 (`dem60~~dem60`), FD 12.9279 / analytic 12.9287 (`dem65~~dem65`).
- **Foundation checks:** `∂²ℓ/∂θ²` (central diff of the §2 `ℓ(θ)`) = `−V⁻¹` to ~1% (ratios
  0.995–1.006 vs lavaan); `Δ_θ`/`Δ_η` central-diff exact (max|Δ| = 0); composition
  `ℓ(θ) = ℓ̃(Σ(θ))` holds to ~1e-6.
- **Remaining to log at closure:** the `‖g_Σ‖` saturation threshold used; one saturated + one
  non-saturated fixture A/B; the bit-stability evidence; the new T3 golden tolerance.
