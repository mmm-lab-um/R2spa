# PLAN 03 — `scoring_matrix` for `get_fs.merMod()` + documentation vignette

Status: **DONE** (2026-08-16; landed as `00bf670` + the PLAN 03 follow-up —
see STATUS.md Closed #8 and #9). Approved decisions at bottom.

## Context (for a fresh session)

- Repo `/home/marklai/R2spa`, branch `refactor/core`, working tree clean at `9c60ff8`
  (PLAN 02 complete, archived as `archive/PLAN_02_backward_compatible.md`).
- Tests: **435 pass / 0 fail**. Check baseline: **2 WARNINGs + 3 NOTEs** (pre-existing).
  Vignettes **13/13 build**. `VignetteBuilder: knitr` — a new `.Rmd` in `vignettes/` is
  auto-discovered; no DESCRIPTION or registry edit needed.
- `R2spa_0.0.4.tar.gz` is **tracked in git** — never delete; after any `R CMD build`,
  restore with `git checkout -- R2spa_0.0.4.tar.gz`.
- Lifecycle (AGENTS.md, never reorder): `devtools::load_all()` → after roxygen changes
  `devtools::document()` → `devtools::test()` → `devtools::check()`.
  Never hand-edit `NAMESPACE` or `man/*.Rd`. No `library()`/`require()` in function bodies.
- Check log pitfall: read `R2spa.Rcheck/00check.log` **before** removing `R2spa.Rcheck/`.

## Verified math (all confirmed at machine precision this session)

Re-verify before editing if desired (Step 0); key reference values from the
`lme1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)` fit:

- `G <- as.matrix(VarCorr(lme1)[[1]])` = `[[612.100, 9.604],[9.604, 35.072]]`;
  `sigma(lme1)` = 25.5918 (`sigma^2` = 654.94).
- **`get_D()` is correct — do not "fix" it.**
  `R2spa:::get_D(lme1@theta) * sigma(lme1)^2 ≡ G` (max diff **0**).
  `get_D()` (`R/get_fs_methods.R:292-295`) reconstructs the **scaled** RE covariance
  `G/σ²` from lme4's theta via `lme4::vec2mlist(theta, symm = FALSE)[[1]]` +
  `tcrossprod()`. The `symm = FALSE` argument is essential (with `symm = TRUE` the
  mirrored triangle breaks the `tcrossprod` structure). σ² cancels against the
  unscaled `Kz = Z'Z`, which is why the formulas below work with scaled `D`.
- Current `fsL ≡ S %*% Z` (loading of the EB estimate `û` on the latent `u`,
  i.e. `E[û | u] = (S·Z) u`): max diff **1.3e-14**. Keep the formula.
- Current `fsT ≡ sigma^2 * S %*% t(S) = Var(û | u)` (score measurement-error
  covariance, matching lavaan `fsT` semantics): max diff **2.0e-13**. Keep the formula.
- **The one real defect:** `get_fs_blocks.merMod()` builds `Kz` from the
  **fixed-effects** design `object@pp$X` (line 255) instead of the random design `Z`.
  Numerically harmless only when `Z == X` (sleepstudy, hence uncaught). For any
  `Z ≠ X` model it is wrong and **crashes**:
  `get_fs(lmer(Reaction ~ Days + (1 | Subject), sleepstudy))` →
  "non-conformable arguments" (verified).
- Scoring matrix identity: with `S_j = solve(solve(G) + crossprod(z_j)/sigma^2) %*%
  t(z_j)/sigma^2` (p×n_j), `S_j %*% (y_j - X_j %*% beta) ≡ ranef[.][[1]][j, ]` —
  max diff **2.9e-13** (Z==X), **3.4e-13** (Z≠X, extra fixed covariate),
  **5.9e-14** (p=1 `(1|Subject)`). `y_j`/`X_j` are the cluster's rows of
  `model.response(model.frame(f))` / `as.matrix(f@pp$X)`; `beta <- fixef(f)`.
- Authoring shortcut: `S_j ≡ inv_W %*% D %*% t(z_j)` using the **existing** loop
  variables (`inv_W = solve(DKz + diag(nrow(Kz)))`, `D = get_D(theta)`,
  `DKz = D %*% crossprod(z_j)`): max diff to the direct formula **3.1e-16**.
- Extracting `Z`: `lme4::getME(object, "Z")` (exported) returns the **full**
  n × (p·G) design matrix in this lme4 version (`object@pp$Z` does not exist —
  `'Z' is not a valid field ... "merPredD"`). Full Z for `lme1` is 180×36
  (18 subjects × 2 RE params); for the p=1 model it is 180×18. Per-cluster slice:
  rows = `case_idx[[j]]` (from `split(seq_len(nrow(model.frame(f))), f@flist[[1]])`),
  columns = `(level_j - 1)*p + 1:p` where
  `level_j <- as.integer(f@flist[[1]])` (factor level index of the cluster of row).
  Verified: rows 1:10 (subject 308) → cols 1:2; rows 11:20 (subject 309) → cols 3:4;
  values equal `model.matrix(~ 1 + Days)` per cluster (only `colnames` differ —
  full Z uses repeated subject labels; drop them and set our own).
- lavaan side (for the vignette): `compute_fscore()`
  (`R/get_fscore_math.R:236-286`): `meany = lambda %*% alpha + nu`;
  `fs = t(a_mat %*% (t(y) - meany) + alpha)`; attribute `scoring_matrix = a_mat`
  (q×p, score×item); `fsL = a_mat %*% lambda`; `fsb = alpha - fsL %*% alpha`;
  `fsT = a_mat %*% theta %*% t(a_mat)`. So scores reconstruct by hand as
  `S %*% (t(y) - (nu + Lambda %*% alpha)) + alpha` from `lavaan::lavInspect(fit, "est")`.
- Consumers of the attributes: `augment_fs()` builds `_se`/`ev_`/`ecov_` columns from
  `fsT` (R/get_fscore.R:281-294); `get_fs.merMod()` packs `fsL_j`/`fsT_j` into 3-way
  array attributes (R/get_fs_methods.R:353-375). `vignettes/multilevel.rmd` feeds
  `ev_u0_eb`-style columns into `tspa_mx_model()` — **no cache chunks, no hard-coded
  narrative numbers** (prose is relative), so value-affecting changes re-render
  cleanly (none expected here — Step 1 is a no-op on Z==X models).
- No test pins merMod `fsT`/`fsL`/SE *values* — only score-identity vs `ranef()`
  (test-get_fscore.R:363-368), column naming (:392-417), and block structure
  (:379-391).

## Decisions (user-confirmed)

1. **Container**: merMod output attribute `scoring_matrix` = **named list per cluster**
   (names = cluster levels, same as `names(blocks)`). Not a 3-way array (unbalanced
   clusters); matches the lavaan convention (per-group named list).
2. **Scope**: add the attribute **and** fix the Z-design bug (`Kz` from `Z`, not
   `pp$X`). `get_D()`, `fsL_j`, `fsT_j` formulas stay unchanged.
3. **Plan file first**: this file is written; implementation is a separate step
   (user chose "write plan file, stop") — but if the executing session is this
   session, proceed through Steps 0-7.

## Steps

- [x] **0. (optional) Re-verify** the identities above with a short `Rscript`
      before editing (key checks: `get_D·σ² ≡ G`; `S(y_j−X_jβ) ≡ ranef`; current
      crash on `(1|Subject)`).
- [x] **1. Code — `R/get_fs_methods.R`**
  - `get_fs_blocks.merMod()` (lines 242-291):
    - After `case_idx <- split(seq_len(nrow(mf)), object@flist[[1]])`:
      ```
      Zmat <- as.matrix(lme4::getME(object, "Z"))
      gids <- as.integer(object@flist[[1]])
      ```
      Guard: `stopifnot(ncol(Zmat) >= num_re * n_clus)` (multi-RE-term models carry a
      wider Z; we slice the **first bar's** block — inherits the existing `[[1]]`
      convention used by `ranef(object)[[1]]`, `cnms[[1]]`, `flist[[1]]`).
    - Per cluster `j` with `idx <- case_idx[[j]]`:
      `zj <- Zmat[idx, (gids[idx[1]] - 1) * num_re + seq_len(num_re), drop = FALSE]`
    - Replace `xj <- as.matrix(Xlist[[j]])` / `Kz <- crossprod(xj)` with
      `Kz <- crossprod(zj)`. Delete `Xlist` (line 255) — it had no other use
      (scores come from `ranef()`; `_se`/`ev_` columns come from `fsT`).
      `mf` is still needed for `case_idx`.
    - `D`, `inv_W`, `fsL_j`, `fsT_j` lines unchanged.
    - New block field:
      `scoring_matrix_j <- inv_W %*% D %*% t(zj)` (p×n_j);
      `rownames(scoring_matrix_j) <- fs_names`;
      `colnames(scoring_matrix_j) <- as.character(seq_len(nrow(zj)))`;
      set `scoring_matrix = scoring_matrix_j` in `blocks[[j]]` (replaces `NULL`, :286).
  - `get_fs.merMod()` (lines 316-378): after `out` is assembled (near the existing
    fsL/fsT array attribute code, :353-375):
      ```
      attr(out, "scoring_matrix") <- setNames(
        lapply(blocks, function(b) b$scoring_matrix), names(blocks))
      ```
- [x] **2. Roxygen** (edit roxygen only, then `devtools::document()`):
  - `get_fs` generic docs, attribute bullets (`R/get_fscore.R` ~:43-59):
    `scoring_matrix` is a named list — lavaan: one matrix per group (score×item);
    merMod: one p×n_j matrix per cluster, `score = S_j %*% (y_j - X_j beta)`.
  - `get_fs.merMod` `@rdname get_fs` block (`R/get_fs_methods.R` :297-315) and the
    legacy `get_fs_lmer()` docs (`R/get_fscore.R` ~:430-470, note its "carries ..."
    list at ~:452): mention the new `scoring_matrix` attribute.
  - `AGENTS.md` attribute listing (~line 66): one-line note — merMod
    `scoring_matrix` is a per-cluster named list (lavaan: per-group matrices).
- [x] **3. Tests** — append `test_that` blocks to `tests/testthat/test-get_fscore.R`
      (existing style: fit model locally per block, `lme1 <- lmer(Reaction ~ Days +
      (Days | Subject), sleepstudy)`, `get_fs(lme1)`):
  1. **Score identity** (core): `sm <- attr(fs, "scoring_matrix")`; for all 18
     clusters `j` with `idx <- which(flist level == j)`, `y_j` from
     `model.response(model.frame(lme1))`, `X_j` from `as.matrix(lme1@pp$X)`:
     `expect_equal(t(sm[[lvl]] %*% (y_j - X_j %*% fixef(lme1))),
     ranef(lme1)[[1]][j, ], tolerance = 1e-10)`.
  2. **Structure**: `expect_named(sm, as.character(levels(sleepstudy$Subject)))`;
     `length(sm) == 18`; each entry 2×10 with rownames `fs_u0`/`fs_u1`.
  3. **Z ≠ X regression**: `lmer(Reaction ~ Days + (1 | Subject), sleepstudy)` —
     `get_fs()` no longer errors; `sm` entries 1×10; identity holds (tol 1e-10).
     (Currently errors "non-conformable arguments" — this test pins the fix.)
  4. **Unbalanced**: `sleepstudy` subset with 2 rows removed (different subjects);
     identity holds for every cluster; entry dims differ (e.g. ×9 vs ×10).
  5. **Legacy**: `get_fs_lmer(lme1)` — `attr(out, "scoring_matrix")` present and
     equal to the non-legacy value (column *names* differ, attribute doesn't).
  - All 435 existing tests must stay green: Step 1 changes nothing numerically
    wherever the old code worked (when `Z == X`, `crossprod(zj) == crossprod(xj)`
    exactly), so sleepstudy outputs and all vignette values are untouched.
- [x] **4. New vignette — `vignettes/scoring-matrices.Rmd`**
  (YAML per existing convention, e.g. `efa-score.Rmd`: `title: "Scoring Matrices: lavaan CFA and lme4"`,
  `author: "Mark Lai"`, `date: "`r Sys.Date()`"`, `output: rmarkdown::html_vignette`,
  `vignette: >` with `%\VignetteIndexEntry{scoring-matrices}`,
  `%\VignetteEngine{knitr::rmarkdown}`, `%\VignetteEncoding{UTF-8}`;
  setup chunk `knitr::opts_chunk$set(collapse = TRUE, comment = "#>")`):
  1. **Overview**: the scoring matrix is the linear operator converting modelled
     observed data into the reported factor/EB scores; `get_fs()` exposes it as the
     `scoring_matrix` attribute alongside `fsL`/`fsT`/`fsb`.
  2. **lavaan CFA**: fit `cfa("f =~ y1+y2+y3+y4", lavaan::PoliticalDemocracy)`;
     show `get_fs()` head + `attr(fs, "scoring_matrix")`; **reconstruct scores by
     hand** `S %*% (t(y) - (nu + Lambda %*% alpha)) + alpha` from
     `lavaan::lavInspect(fit, "est")` and demonstrate agreement with the `fs_f`
     column (~1e-8); explain `fsL = S·Λ` (implied loadings of factors on scores),
     `fsb = α − fsL·α` (implied intercepts), `fsT = S·Θ·S'` (score error variance).
  3. **lme4**: fit `lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
     REML = FALSE)` (2S-PA convention from `multilevel.rmd`); show `get_fs()` head +
     `attr(fs, "scoring_matrix")` (named list, one 2×10 matrix per subject);
     **reconstruct scores** `S_j %*% (y_j - X_j %*% beta)` per cluster and demonstrate
     agreement with `fs_u0`/`fs_u1` (~1e-10) and with `lme4::ranef()`; state formulas:
     `S_j = (G⁻¹ + Z_j'Z_j/σ²)⁻¹ Z_j'/σ²`; `D ≡ G/σ²` reconstructed from theta
     (`get_D`), why σ² cancels against `Kz = Z_j'Z_j`; `fsL = S_j Z_j` (loading of `û`
     on `u`, = `E[û|u]` operator); `fsT = σ² S_j S_j' = Var(û | u)` (score
     measurement-error covariance). Notes: single RE term (first `[[1]]` convention);
     unbalanced clusters handled (per-cluster list); `Z ≠ X` supported.
  4. **Comparison table** (`knitr::kable`): rows = `scoring_matrix` / shape &
     orientation / container (per-group list vs per-cluster list) / score-
     reconstruction formula; columns = lavaan CFA / lme4.
  - No new dependencies (lavaan, lme4, knitr, kable already used in vignettes);
    no cache chunks (fits take seconds at build time).
- [x] **5. Document & test**: `devtools::document()` → `devtools::test()`
      (expect 435 existing + ~5 new, 0 fail).
- [x] **6. Check & status**: `devtools::check()` (full — OpenMx vendored compile makes
      it slow; vignettes now 14/14). Read `R2spa.Rcheck/00check.log` before removing
      `R2spa.Rcheck/`. Update `STATUS.md`: add item (merMod `scoring_matrix` +
      vignette), mark completed with date; bump verification-section vignette count
      13 → 14.
- [x] **7. Report** diff summary to user. Commit **only if explicitly requested**.

## Risks / notes

- `as.matrix(Zmat)` is dense n × (p·G): fine for typical fits (sleepstudy 180×36);
  for very large clustered data this is the main memory consideration. No action.
- Multi-RE-term models: first-term-only convention inherited (documented in roxygen +
  vignette); not a new limitation.
- REML vs ML: the score identity holds either way; vignette uses `REML = FALSE`.
- Naming: `colnames` of the per-cluster `S_j` are cluster-position labels
  `"1".."n_j"` (order = model-frame row order, which lme4 sorts by cluster).
- Do not touch: `get_D()`, `fsL_j`/`fsT_j` formula lines, lavaan path,
  `assemble_fs_blocks()`, `fs_to_group_list()`.

## Approved decisions (2026-08-16)

- Container: **named list per cluster** (user selected "Named list (Recommended)").
- Fix scope: **full fix** — Z-design bug fix included; user corrected the prior
  misreading of `get_D()` (it is correct — scaled RE covariance `G/σ²`) and of
  `fsL` (it is the loading of `û` on `u`: `fsL = S %*% Z`).
- Workflow: plan file written first; implementation deferred unless the user says to
  proceed in-session.
