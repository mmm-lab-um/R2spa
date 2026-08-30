# Salvaged drafts from deleted stale branches (2026-08-27)

These files were the only unique content on remote branches that were deleted on
2026-08-27 as part of the pre-release branch cleanup. Source branches (all on
`origin`, gone now) and the commits they pointed at:

| Directory | Branch (tip) | Context |
|---|---|---|
| `branch-marklhc-issue83/` | `marklhc/issue83` (`9d7b84a`) | Open issue #83 (factor score reliability after alignment). `alignment-reliability.rmd` = working notes; `compare-identification.*` = identification-comparison vignette + cached fixtures (also carried by the deleted `identification` branch, of which it is a superset). |
| `branch-marklhc-issue87/` | `marklhc/issue87` (`bc336e1`) | Open issue #87 (probit regression example with Mx). `vignettes/probit.rmd` = 93-line draft. |
| `branch-openmx/` | `openmx` (`d677541`) | 2023 categorical-interaction exploration that preceded the current `compute_fs_prod()` work. Salvaged: the three `Simulation_*.R` scripts, the `vignettes/test_cat_int.Rmd` draft, and the small input/summary data files. **Not** salvaged: 20 machine-generated `SimDesign-results_*/results-row-*.rds` output dumps (reproducible by running the scripts), and its two 2-line edits to then-existing vignettes (`categorical-interaction.Rmd` typo fix — already present in `.quarantine/vignettes/categorical-interaction.Rmd`; `tspa-vignette-mx.Rmd` trailing newline — file since rewritten). |
| `branch-vignette/` | `vignette` (`5908960`) | Early vignette drafts. Salvaged: `invariance-factor-score.Rmd`, `three_var.Rmd`. **Not** salvaged: its 109-line `efa-score.Rmd` draft, superseded by the current `vignettes/efa-score.Rmd` (131 lines). |

Everything here is build-ignored via `^archive$` in `.Rbuildignore`. These are
drafts, not package code — don't edit package code to match them.
