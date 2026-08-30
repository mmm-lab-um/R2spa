# Get Factor Scores and the Corresponding Scoring Matrices for Mixed-Effect Models

`get_fs_lmer()` is superseded by
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md). It
is retained for backward compatibility and delegates to `get_fs(object)`
internally. New code should call
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
directly.

## Usage

``` r
get_fs_lmer(
  object,
  method = c("EB"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  fsm = FALSE,
  legacy_names = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted model object of class
  [lme4::lmerMod](https://rdrr.io/pkg/lme4/man/merMod-class.html).

- method:

  `"EB"` (empirical Bayes, default) or `"ML"` (prior-free per-cluster
  OLS), forwarded to
  [`get_fs.merMod()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md).

- corrected_fsT:

  Currently not used.

- vfsLT:

  Currently not used.

- fsm:

  Currently not used.

- legacy_names:

  Logical. Passed to
  [`get_fs.merMod()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md).
  Defaults to `TRUE` so `get_fs_lmer()` keeps returning the pre-refactor
  `u0_eb`-style *column names* (in the legacy column order). Note the
  legacy output is name-compatible, not byte-identical, with the
  pre-refactor result: it additionally carries score-error columns
  (`u0_eb_se`, ...), per-cluster `fsL`/`fsT` array attributes, a
  per-cluster `scoring_matrix` list attribute (see
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)),
  and has NULL row names (the pre-refactor output had no `_se` columns,
  no attributes, and used the ranef subject IDs as row names).

- ...:

  Additional arguments, passed on to
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md).
