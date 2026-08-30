# Compute product factor-score indicators (double-mean-centered)

`compute_fs_prod()` takes a single-group
[`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
result and appends, for each requested pair of distinct latent
variables, the double-mean-centered (DMC) product indicator of the two
factor scores together with its per-row standard error and its implied
loading. It is the re-integrated successor of the quarantined
`get_fs_int()` student function: the column-naming conventions are kept
(the product column joins the two full score column names, so for latent
names `a` and `b` the columns are `fs_a:fs_b`, `fs_a:fs_b_se`, and
`fs_a:fs_b_ld` — e.g. the product of `fs_ind60` and `fs_dem60` is
`fs_ind60:fs_dem60`), while the standard error is computed from the
stage-1 `fsL`/`fsT`/`psi` attributes — the general joint-model formula
below, instead of the separate-single-factor special case the old
function used.

## Usage

``` r
compute_fs_prod(fs, product)
```

## Arguments

- fs:

  A single-group
  [`get_fs()`](https://mmm-lab-um.github.io/R2spa/reference/get_fs.md)
  result: the unified data frame (`format = "unified"`, the default) or
  the single-group data frame of `format = "list"`. Rejected with an
  informative error (v1: single-group lavaan models only): multi-group
  results (a named list of per-group data frames, or a unified data
  frame with a group column), `merMod` results (3-D `fsT` attribute),
  and per-observation results (the `per_obs`/`mirt_per_obs` markers —
  local-mode FIML merges and `mirt` models).

- product:

  A character string of the form `"a:b + c:d"` (pairs of distinct latent
  names, `+`-separated) or a list of length-2 character vectors of
  latent names (a 2-column matrix or data frame of names is coerced to a
  list). Same-factor products (`"x:x"`) and duplicated pairs are
  rejected (v1).

## Value

The input with, per requested pair in order, three appended columns:
`fs_a:fs_b` the DMC product indicator, `fs_a:fs_b_se` the per-row
standard error (constant for complete data; one value per
observed-indicator pattern under FIML; `NA` for unscorable rows), and
`fs_a:fs_b_ld` the per-row implied loading `gamma` (pattern-specific
under FIML because `fsL` varies by pattern). All existing columns and
attributes are untouched.

## Details

Let \\fs_k = L_k \xi + b_k + e_k\\ for \\k = 1, \ldots, q\\, where \\L\\
is the block's `fsL` (\\q \times q\\, rows = scores, cols = latents),
\\T\\ is the block's `fsT` (score error covariance), \\\Psi\\ is the
latent (co)variance (the `psi` attribute), and \\e\\ is the score error.
The derivation relies on: \\E\[e\] = 0\\ (absorbed in \\b\\); \\Cov(e,
\xi) = 0\\ elementwise for ALL latents (true for both regression and
Bartlett scoring, where \\e\\ is linear in the indicator errors); and
joint normality of \\(\xi, e)\\.

For a pair of distinct latents \\(a, b)\\, the double-mean centered
(DMC) product indicator is

\$\$P = (fs_a - mi_a)(fs_b - mi_b) - mi_P\$\$

where \\mi_a\\ and \\mi_b\\ are the SAMPLE means of the score columns
and \\mi_P\\ the sample mean of the (component-centered) product. Its
conditional expectation is

\$\$E\[P \mid \xi\] = (L_a (\xi - \alpha))(L_b (\xi - \alpha))\$\$

(the centering constants cancel in the error), so the measurement error
is

\$\$u = P - E\[P \mid \xi\] = L_a (\xi - \alpha) e_b + L_b (\xi -
\alpha) e_a + e_a e_b\$\$

Under joint normality (Isserlis), for \\\tau_k = L_k \Psi L_k'\\ (a
scalar: row \\k\\ of \\L\\ times \\\Psi\\ times its transpose),
\\\tau\_{ab} = L_a \Psi L_b'\\, \\s_k^2 = T\[k, k\]\\, and \\c = T\[a,
b\]\\:

\$\$se_P^2 = \tau_a s_b^2 + \tau_b s_a^2 + s_a^2 s_b^2 + c^2 + 2
\tau\_{ab} c\$\$

(the \\c^2\\ coefficient is 1, not 2: \\Var(e_a e_b) = E\[e_a^2
e_b^2\] - c^2 = (s_a^2 s_b^2 + 2 c^2) - c^2 = s_a^2 s_b^2 + c^2\\).

The implied loading is

\$\$\gamma = L\[a, a\] L\[b, b\] + L\[a, b\] L\[b, a\]\$\$

— the coefficient of \\\xi_a \xi_b\\ in \\E\[P \mid \xi\]\\ (Bartlett
joint model: \\\gamma = 1\\; separate single-factor models: \\\gamma =
\lambda_a \lambda_b\\).

Sanity check: with diagonal \\L\\ and \\T\\ (separate single-factor
measurement models, \\\lambda_k = L\[k, k\]\\) the formula reduces to
the classic \\\lambda_a^2 \Psi\[a, a\] s_b^2 + \lambda_b^2 \Psi\[b, b\]
s_a^2 + s_a^2 s_b^2\\ — the formula the old (quarantined) `get_fs_int()`
had. The terms \\c^2\\ and \\2 \tau\_{ab} c\\ are exactly what the old
code dropped; they are nonzero in joint models (correlated factors,
cross-loadings, error covariances).

Per-row values are pattern-resolved: under FIML each observed-indicator
pattern carries its own `fsL`/`fsT` block, so `fs_a:fs_b_se` and
`fs_a:fs_b_ld` take one value per pattern (the implied loading varies by
pattern because `fsL` does), while the latent (co)variance `psi` is
group-level and shared. For complete data the standard error is constant
across rows. Rows with no scorable pattern (fully-missing rows) get `NA`
in all three columns.

## References

The joint-normal moment expansion (Isserlis's theorem) used to derive
`se_P^2` and the measurement-error covariances is:

Isserlis, L. (1918). On a formula for the product-moment coefficient of
any order of a normal frequency distribution in any number of variables.
*Biometrika*, 12(1-2), 134-139.
[doi:10.1093/biomet/12.1-2.134](https://doi.org/10.1093/biomet/12.1-2.134)
.

## See also

- `vignette("Product factor-score indicators (latent interactions)", package = "R2spa")`
  for the full latent-interaction workflow.

## Examples

``` r
library(lavaan)
fit <- cfa(
  "ind60 =~ x1 + x2 + x3
   dem60 =~ y1 + y2 + y3 + y4",
  data = PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")],
  std.lv = TRUE
)
fs <- get_fs(fit, product = "ind60:dem60")
head(fs[c("fs_ind60", "fs_dem60", "fs_ind60:fs_dem60",
          "fs_ind60:fs_dem60_se", "fs_ind60:fs_dem60_ld")])
#>     fs_ind60   fs_dem60 fs_ind60:fs_dem60 fs_ind60:fs_dem60_se
#> 1 -0.8101568 -1.3119114         0.6045906            0.3449699
#> 2  0.1888466 -1.3644831        -0.7159414            0.3449699
#> 3  1.0960931  1.3107705         0.9784631            0.3449699
#> 4  1.8702043  1.4849083         2.3188185            0.3449699
#> 5  1.2446060  0.9193277         0.6859375            0.3449699
#> 6  0.3348621  0.4886331        -0.2946387            0.3449699
#>   fs_ind60:fs_dem60_ld
#> 1             0.831191
#> 2             0.831191
#> 3             0.831191
#> 4             0.831191
#> 5             0.831191
#> 6             0.831191

# The same columns, computed directly from a get_fs() result
fs2 <- compute_fs_prod(get_fs(fit), product = "ind60:dem60")
head(fs2[c("fs_ind60:fs_dem60", "fs_ind60:fs_dem60_se")])
#>   fs_ind60:fs_dem60 fs_ind60:fs_dem60_se
#> 1         0.6045906            0.3449699
#> 2        -0.7159414            0.3449699
#> 3         0.9784631            0.3449699
#> 4         2.3188185            0.3449699
#> 5         0.6859375            0.3449699
#> 6        -0.2946387            0.3449699
```
