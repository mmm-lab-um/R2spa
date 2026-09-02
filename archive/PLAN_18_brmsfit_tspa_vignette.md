# PLAN 18 — Vignette: 2S-PA with `brms` fit factor scores

**Date:** 2026-08-31
**Branch:** `get_fs-brms` (the `get_fs.brmsfit()` method lives here; `devel` is brms-free).
**Owner roles:** `r-doc` (vignette), `r-architect` (technical review), `r-tester` (fixtures).
**Status:** implemented (2026-09-02; vignette at `vignettes/tspa-brms.Rmd`, fixture
`vignettes/sim_brms_ri_vignette.rds`; PLAN 17 cherry-picked into this branch first as a
prerequisite; `R CMD check` 0/0/0)
**Relates to:** `vignettes/tspa-brms-random-sigma.Rmd` (the location-scale / random-`sigma`
case, already built). This is the **general** brmsfit 2S-PA companion: ordinary random
intercepts/slopes scored as factor scores and corrected in `tspa()`. Cross-ref PLAN 17 (the
`cbind()` loading-recovery fix) for the structural-model flow.

---

## 1. Purpose

Document the full **brms → 2S-PA** workflow that is *not* the location-scale case: a Gaussian
mixed model with a **single random-effects term** (random intercept, and by extension a random
slope) is scored by `get_fs.brmsfit()` into a factor score + its SE + implied loading, and the
score is corrected for measurement error in `tspa()`. The vignette makes three things explicit
that the random-sigma vignette does not:

1. The **scoring convention**: the brmsfit score is the **posterior mean** of the random effect
   (an EAP-analog / Bayesian EB score), and the implied loading is
   `fsL = 1 − Vpost/ψ = w` (the shrinkage weight, derived from the **posterior variance**), with
   `fsT = fsL·Vpost`. This is the *same* convention as `mirt` and the lme4 `method="EB"` path —
   `fsL` is **not** `1`.
2. The two scoring **methods** available to `get_fs.brmsfit()`: `method="EB"` (posterior mean,
   the default) and `method="ML"` (the prior-free group-mean / Bartlett analog).
3. The **`cbind()` footgun** when a distal is added for a structural model (cross-ref PLAN 17),
   and the robust way to fit one.

## 2. Scope

- **In:** one brms Gaussian model, `y ~ … + (1 | g)` (and a brief `(1 + x | g)` random-slope
  note); `get_fs(fit, method = "EB"/"ML")`; the stage-1 output (score, `fsL`, `fsT`, SE, the
  `u0_by_fs_u0` / `ev_fs_u0` columns); `tspa()` **latent-variance** model (`u0 ~~ u0`) and a
  **structural** model (`d ~ u0`); a **validation against the lme4 analog** (the frequentist
  gold standard); a pointer to the exact per-subject `tspa_mx_model()` route.
- **Out (v1, defer to the existing random-sigma vignette or a follow-up):** location-scale /
  random-`sigma` models (already in `tspa-brms-random-sigma.Rmd`); **multiple** random-effects
  grouping factors (unsupported by `get_fs.brmsfit()`); multigroup brms; `corrected_se` /
  grand-standardization on brms results.
- **Companion, not a rewrite:** cross-link from and to `tspa-brms-random-sigma.Rmd` and the
  `R2spa.Rmd` hub.

## 3. Vignette outline (`r-doc`)

File: `vignettes/tspa-brms.Rmd` (suggested; short and general, distinct from the random-sigma
title). YAML index entry "2S-PA with brms fit factor scores". Standard
`rmarkdown::html_vignette`, `knitr::opts_chunk$set(collapse = TRUE, comment = "#>")`.

1. **Intro.** 2S-PA framing (stage 1 `get_fs()`, stage 2 `tspa()`); the motivation: a latent
   trait measured by repeated / clustered observations, scored from a brms fit, corrected for
   score error before relating it to a distal outcome (Cole & Preacher, 2014 — same refs as the
   existing vignettes).
2. **The data.** A seeded random-intercept DGP (`y_ij = b0 + u_j + e_ij`, `d_j = β u_j + v_j`),
   or `lme4::sleepstudy` (real, no DGP assumptions). State `n` subjects × `m` reps and the
   brms model `y ~ 1 + (1 | g)` (and a one-line `(1 + x | g)` variant note).
3. **Stage 1 — `get_fs.brmsfit()`.** Fit `brm(..., chains, iter, seed)`; `fs <- get_fs(fit)`.
   Show the output columns and attributes; **work the identity** `fsL = 1 − Vpost/ψ = w` and
   `fsT = fsL·Vpost` on the numbers (so the reader sees the loading is the shrinkage weight, not
   1). Sub-section **EB vs ML**: `get_fs(fit, method="EB")` (posterior mean) vs
   `method="ML"` (prior-free group mean), with a sentence on when to prefer each.
4. **Stage 2a — latent-variance model.** `tspa("u0 ~~ u0", data = fs)` — recovers the random
   effect's variance (compare to `lme4::VarCorr` / the brms posterior of `sd_g`). This uses the
   **derived** form directly (no distal, so no `cbind`, so no footgun — it works today).
5. **Stage 2b — structural model.** Relate the trait to a distal `d`. Show the **robust**
   explicit form `tspa("d ~ u0", data = cbind(fs, d=d), fsL = attr(fs,"fsL"),
   fsT = attr(fs,"fsT"))` (always correct), then a callout box: *adding a column drops the
   attributes, so the derived `tspa("d ~ u0", data = cbind(fs, d=d))` form relies on PLAN 17 to
   recover the loading; until then, pass `fsL`/`fsT` explicitly.* (This keeps the vignette
   correct whether or not PLAN 17 has landed.)
6. **Validation against lme4.** Fit the same model in `lme4` (`lmer`, `REML=FALSE`), score with
   `get_fs(mer_fit, method="EB")`, and show the brmsfit and lme4 2S-PA estimates agree (both
   recover the latent variance and the structural coefficient), contrasting with the naive
   (uncorrected) regression. This is the empirical payoff and ties brms to the familiar lme4.
7. **The exact route (pointer).** One short paragraph: `tspa_mx_model()` honours the
   person-specific `fsT` (no pooling); cross-ref `vignettes/tspa-vignette-mx.Rmd`.
8. **Scope and notes.** `brms` is a `Suggests` dep (guarded `require_brms()`); single RE term
   only; EB-only vs ML; posterior-based scoring (no design matrix); links out.
9. **References.** Cole & Preacher (2014) + the brms / lme4 references (match the existing
   vignettes' reference style).

## 4. Fixtures (`r-tester`)

- A seeded DGP + a brms random-intercept fit is cheap enough to run **live** for a small
  (`n ≤ ~60`, `chains = 2, iter = ~2500`) vignette; if build time is a concern, follow the
  random-sigma vignette's pattern: run the fit once, `saveRDS` the `get_fs()` result (and the
  distal) to `vignettes/sim_brms_ri_vignette.rds`, and `readRDS` it in the body (the code that
  builds the fixture sits in a collapsed/hidden chunk).
- The lme4 analog is trivial to run live (no MCMC).
- Pin the seed; keep the fixture small so `R CMD check` vignette build stays fast.
- If `brms` is unavailable the vignette chunk should be guarded (the package already
  `require_brms()`-guards `get_fs.brmsfit()`); note the vignette build assumes `brms` installed
  (it is a `Suggests` and pulled for the vignette build).

## 5. Lifecycle (exact order, per AGENTS.md)

- **No R changes** in this plan → `devtools::document()` is a no-op (run it only if PLAN 17
  lands in the same branch and changes roxygen).
- `devtools::load_all()` → build the vignette chunk by chunk (verify the brms fit / `get_fs` /
  `tspa` all run and the numbers in §3.6 agree) → `devtools::test()` (fixtures load; no existing
  test touched) → `devtools::check()` (as-cran, `--no-manual`): the new vignette builds with
  0 errors/warnings.
- Update the `R2spa.Rmd` hub with a relative `.html` link to the new vignette (the P2 hub
  convention: plain sibling links, since `?vignette=` does not resolve under rmarkdown 2.31).

## 6. Exit criteria

- `vignettes/tspa-brms.Rmd` builds under `devtools::check()` with 0 errors / 0 warnings.
- The EB-vs-ML, `fsL = 1 − Vpost/ψ`, lme4-validation, and structural-model sections all render
  with the numbers the text claims (spot-check the identity on the displayed values).
- Cross-links present: to `tspa-brms-random-sigma.Rmd`, `tspa-vignette-mx.Rmd`, and from the
  `R2spa.Rmd` hub.
- Suite unchanged (0 fail); `R CMD check` 0/0/0 on the `get_fs-brms` tree.

## 7. Implementation notes — deviations from the outline

Two premises in the outline did not hold for the single-factor (one latent + one observed
distal) structural model and were corrected in the vignette:

1. **§3.5 explicit-vs-derived is inverted.** The outline recommended the *explicit*
   `fsL`/`fsT` form as "always correct" and the *derived* `cbind()` form as the footgun. In
   practice the explicit **multi-factor** `fsL`/`fsT` form *fails* on a brms (or lme4) result
   plus an observed distal — it routes through `fs_indiv()`, which requires the `fs_pattern`
   attribute that neither `get_fs.brmsfit()` nor `get_fs.merMod()` carries (only the lavaan
   path sets it). The **derived single-factor** form is the one that works, and with PLAN 17 it
   recovers the implied loading. So the vignette uses the derived form as the correct route and
   frames the footgun as *loading recovery*, not the `cbind()` itself. This is why PLAN 17 had
   to be cherry-picked into this branch before the vignette could be correct.
2. **"Naive is biased" does not apply.** For a random-intercept model whose distal is generated
   from the *true* random effect, the naive OLS of the distal on the EB score is **unbiased in
   expectation** (it equals `beta` exactly, because the EB attenuation and its noise cancel in
   `Cov/Var`). The real, demonstrable contrast is the **loading**: a unit-loading correction
   (the explicit `se_fs` form, or the pre-PLAN-17 derived form) over-corrects and biases the
   coefficient (~1.28 here), while the PLAN-17-derived form recovers the shrinkage-weight
   loading and is correct (~1.09), matching the lme4-EB analog (~1.08). `n` was raised 60 → 200
   so the estimates sit closer to the true values.
