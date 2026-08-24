# R2spa's single boundary to lavaan internals.
#
# DESIGN RULE (PLAN 04): this file is the ONLY place in the package that
# reads lavaan internals (partable/list views, `@Data@` slots, `lavTech`
# views). Every other file consumes the `tsp_*` wrappers below, so a future
# lavaan partable/format drift surfaces here (and in the canary tests) as one
# localized failure with a table diff — never as a user-discovered
# estimation bug. Capability probing (column inspection) is used in place of
# `packageVersion()` string gating, so patch/dev releases are handled
# without edits; the error path names the installed and tested-up-to
# versions.
#
# ---------------------------------------------------------------------------
# Lavaan dependency contract (audit artifact for the next format drift)
#
# The package declares no lavaan version bound in DESCRIPTION (PLAN 04,
# decision 1); this table is the declared contract. Every entry was
# verified live against lavaan 0.7-2 (packaged 2026-07-16) on R 4.6.1.
#
# Consumed here (R/lavaan_compat.R):
#   View / slot                                Semantics (0.7-2)
#   ----------------------------------------   --------------------------------
#   partable(fit) / lavInspect(what = "list")  parameter table; columns:
#     id, lhs, op, rhs, user, block, group,
#     free, ustart, exo, label, plabel, start,
#     est, se
#     - lhs/op/rhs  statement identifiers (`=~`/`~~`/`~`/`|`/`~~|` ...)
#     - user        1 = row written in the model syntax; 0 = added by
#                   lavaan (0.7-2; absent in 0.6.x)
#     - block       parameter block index (1 = lambda, 2 = theta, ...)
#     - group       1-based group index
#     - free        0.7-2: 1-based position of the parameter in the free
#                   estimate vector; 0 = fixed (0.6.x used a `fix` column
#                   instead — drift within the 0.7 line: 0.6.x `free` was a
#                   plain 0/1 flag)
#     - ustart      user-specified start value (start()); NA otherwise
#     - exo         1 = involves an exogenous variable
#     - label/plabel  parameter labels (auto `.pN.` when unnamed)
#     - start       fixed value for fixed rows; for free rows 0.7-2 fills
#                   the final estimate (data-dependent — not goldenable)
#     - est/se      point estimate / SE (excluded from the canonical df)
#   raw `free` column of that same view
#     per-row 1-based position of the parameter in the free estimate vector
#     (0 = fixed); consumed via tsp_partable_positions() by
#     R/grandStandardizedSolution.R (partable-row <-> matrix-position
#     mapping), canary-covered in test-lavaan_compat.R
#   fit@Data@nobs / @Data@ngroups / @Data@norig
#                                          per-group nobs list / group count /
#                                          per-group norig list (`lavaanData`
#                                          S4 slot); values proven equal to
#                                          lavInspect(what = <name>) (STATUS
#                                          item 6)
#   fit@optim$converged / lavInspect("converged")
#                                          optimizer convergence flag; scalar
#                                          logical (MG fits optimize the
#                                          combined free-parameter vector in
#                                          a single call; verified 0.7-2,
#                                          SG + MG); read via tsp_converged()
#   lavTech(fit, what = "est" | "free")        flat list of 6 matrices
#                                             (lambda, theta, psi, beta, nu,
#                                             alpha) repeated per group; no
#                                             dimnames
#   lavTech(fit, what = "partable",
#           list.by.group = TRUE)              named list of lists, one per
#                                             group; per-group free-position
#                                             matrices (`$beta` etc.; 1s mark
#                                             free parameters). Single-group:
#                                             unnamed length-1 list.
#   lavInspect(fit, what = "est")              0.7-2 fallback for lavTech:
#                                             lavaan.list nested per group,
#                                             elements named `<group>.<block>`
#   coef(fit) / vcov(fit)                      named vector / name-indexed
#                                             matrix — the stable name-based
#                                             anchor for tests (no positional
#                                             row subsetting; 0.7-2 changed
#                                             lavInspect list-access
#                                             semantics)
#   fit@vcov                                     S4 slot; list of
#                                                 {se = method string,
#                                                 information = string,
#                                                 vcov = p x p covariance
#                                                 of the free parameters}.
#                                                 The only field written
#                                                 by the package:
#                                                 fit@vcov[["vcov"]], via
#                                                 tsp_set_vcov() for the
#                                                 corrected-SE path (R/
#                                                 tspa_corrected_se.R);
#                                                 canary-covered. The
#                                                 method-string keys
#                                                 (fit@vcov$se,
#                                                 @Options$se) stay
#                                                 intact, so se() /
#                                                 standardizedSolution()
#                                                 dispatch is unchanged
#
# Known but NOT wrapped (deliberately deferred, PLAN 04 §1 — most stable
# views; detection falls to the existing equivalence tests rather than the
# canary):
#   - lavInspect(what = "est"/"data")      R/get_fs_methods.R
#   - unlist(fit@Data@norig)               R/get_fs_methods.R
#   - lavInspect(what = "data"/"meanstructure"/"implied"/"free"/"nobs"/
#     "ngroups"/"orig")                    R/get_fscore_math.R
#   - fit@implied$cov / fit@implied$mean   vignettes (user-facing examples)
#   - lavaan::lav_func_jacobian_complex    R/get_fscore_math.R (correct_evfs,
#     compute_grad_ld_evfs, compute_fsrel — purely algebraic closures, no
#     optimizer boundary => complex steps valid),
#     R/grandStandardizedSolution.R (re-integrated, same shape).
#     NOT used by the corrected-SE path: R/tspa_corrected_se.R
#     computes its stage-2 Jacobian by explicit central differences (stepped
#     refits through the optimizer silently degrade; complex literals die in
#     the model-string parser).
#
# ---------------------------------------------------------------------------
#
# Layout state is memoized per lavaan version (first probe wins).

tsp_compat_env <- new.env(parent = emptyenv())

# Lavaan release up to which the partable layout was verified live.
tsp_lavaan_tested_up_to <- "0.7-2"

tsp_unsupported_layout <- function() {
  stop(
    "lavaan ", as.character(utils::packageVersion("lavaan")),
    " partable layout not supported; R2spa is tested up to lavaan ",
    tsp_lavaan_tested_up_to, call. = FALSE
  )
}

tsp_has_lavTech <- function() {
  exists("lavTech", envir = asNamespace("lavaan"), inherits = FALSE)
}

# Raw partable view as a plain data.frame (class stripped for stable
# golden comparison). Primary source is `lavInspect(what = "list")` — the
# view the migrated call sites always consumed (keeps column storage modes
# bit-identical to the pre-Phase-1 output); `partable()` is the fallback.
tsp_partable_raw <- function(fit) {
  raw <- tryCatch(
    lavInspect(fit, what = "list"),
    error = function(e) lavaan::partable(fit)
  )
  as.data.frame(raw)
}

# Resolve column meanings from an actual partable view (capability probing,
# never version strings). Memoized by `tsp_layout()`.
tsp_resolve_layout <- function(pt) {
  must <- c("lhs", "op", "rhs", "group", "block", "label")
  if (!all(must %in% names(pt))) tsp_unsupported_layout()
  if ("free" %in% names(pt)) {
    fixed_flag <- "free" # 0.7-2: free = 1 means free
  } else if ("fix" %in% names(pt)) {
    fixed_flag <- "fix" # 0.6.x: the fix column carried the fixed value;
    # a legacy fixed value of exactly 0 is indistinguishable from a free
    # row (best-effort path; the active 0.7-x path is canary-protected)
  } else {
    tsp_unsupported_layout()
  }
  if ("start" %in% names(pt)) {
    value_col <- "start"
  } else if ("fix" %in% names(pt)) {
    value_col <- "fix" # 0.6.x: the fix column carried fixed values
  } else {
    tsp_unsupported_layout()
  }
  list(
    fixed_flag = fixed_flag,
    value_col = value_col,
    user_col = if ("user" %in% names(pt)) "user" else NULL,
    ustart_col = if ("ustart" %in% names(pt)) "ustart" else NULL,
    exo_col = if ("exo" %in% names(pt)) "exo" else NULL
  )
}

# Lazy, memoized layout probe: resolves the fixed-indicator column
# (`free` now, `fix` historically), the fixed-value column (`start` now,
# `fix` historically), the user-row flag, and the group/label/value columns
# from one partable view. Unknown shape -> explicit error naming the
# installed and tested-up-to versions (never silent guessing).
tsp_layout <- function(fit) {
  ver <- as.character(utils::packageVersion("lavaan"))
  if (identical(tsp_compat_env$version, ver) &&
      !is.null(tsp_compat_env$layout)) {
    return(tsp_compat_env$layout)
  }
  layout <- tsp_resolve_layout(tsp_partable_raw(fit))
  tsp_compat_env$layout <- layout
  tsp_compat_env$version <- ver
  layout
}

tsp_layout_reset <- function() {
  tsp_compat_env$layout <- NULL
  tsp_compat_env$version <- NULL
}

# Canonical partable read (read-only by design — PLAN 04, decision 2): one
# row per parameter, in partable order (statement order; lavaan-added
# parameters last). Columns:
#   lhs, op, rhs        statement identifiers
#   value               fixed value where the row is fixed, NA where free
#   free                1 = free, 0 = fixed
#   group               1-based group index
#   block               parameter block index
#   label               parameter label
#   user                1 = user-written row (NA if the layout has no
#                       `user` column)
#   ustart              user-specified start value (NA if absent)
#   exo                 exogenous flag (kept for grand_standardized_solution()
#                       row identification, whose output carries `exo`)
# The lavaan.data.frame class is intentionally dropped.
tsp_partable_read <- function(fit) {
  lo <- tsp_layout(fit)
  pt <- tsp_partable_raw(fit)
  if (lo$fixed_flag == "free") {
    fixed <- pt$free == 0
  } else {
    fixed <- pt$fix != 0
  }
  value <- pt[[lo$value_col]]
  value[!fixed] <- NA
  # `free` is 1 for free parameters in both layouts (0.7-2 `free` = 1 means
  # free; 0.6.x `fix` = 1 means fixed, so free = !(fix != 0)).
  out <- data.frame(
    lhs = pt$lhs,
    op = pt$op,
    rhs = pt$rhs,
    value = value,
    free = as.integer(!fixed),
    stringsAsFactors = FALSE
  )
  out$group <- pt$group
  out$block <- pt$block
  out$label <- pt$label
  out$user <- if (!is.null(lo$user_col)) pt[[lo$user_col]] else NA
  out$ustart <- if (!is.null(lo$ustart_col)) pt[[lo$ustart_col]] else NA
  out$exo <- if (!is.null(lo$exo_col)) pt[[lo$exo_col]] else NA
  out
}

# Global position of each partable row in the free-estimate vector, 1-based
# (0 for fixed rows). On the 0.7-x line the raw partable `free` column
# carries exactly these positions; this is the row <-> matrix-position
# anchor grand_standardized_solution() uses to assign standardized
# estimates to partable rows: free-position matrices alone do not order the
# partable rows, because a matrix block is traversed column-major while
# partable rows follow model-statement order (they coincide only by
# accident, e.g. single-predictor models). A layout where `free` is a bare
# 0/1 flag (0.6.x) fails the permutation check below and errors loudly
# rather than silently mislabeling rows downstream.
tsp_partable_positions <- function(fit) {
  pt <- tsp_partable_raw(fit)
  if (!"free" %in% names(pt)) tsp_unsupported_layout()
  p <- suppressWarnings(as.numeric(pt$free))
  tsp_check_free_positions(p)
  p
}

# Guard for the 0.7-x `free`-column semantics (see tsp_partable_positions):
# the per-row positions must be a 0-extended 1..npars permutation. A bare
# 0/1 flag layout (0.6.x) fails this and errors loudly rather than
# silently mislabeling standardized-solution rows downstream.
tsp_check_free_positions <- function(p) {
  nfree <- sum(p > 0)
  if (nfree > 0 && !all(sort(p[p > 0]) == seq_len(nfree))) {
    stop(
      "lavaan ", as.character(utils::packageVersion("lavaan")),
      " partable `free` column does not carry 1..npars positions; R2spa is ",
      "tested up to lavaan ", tsp_lavaan_tested_up_to, call. = FALSE
    )
  }
  invisible(TRUE)
}

# Per-group block matrices (`lambda`, `theta`, `psi`, `beta`, `nu`,
# `alpha`) repeated in one flat list — the `lavTech(what = "est")` shape —
# with a `lavInspect(what = "est")` fallback (0.7-2: nested per-group
# lavaan.list with `<group>.<block>` element names, flattened here).
tsp_model_matrices <- function(fit) {
  if (tsp_has_lavTech()) {
    return(lavaan::lavTech(fit, what = "est"))
  }
  tsp_flatten_grouped_est(lavInspect(fit, what = "est"))
}

tsp_free_matrices <- function(fit) {
  if (tsp_has_lavTech()) {
    return(lavaan::lavTech(fit, what = "free"))
  }
  tsp_flatten_grouped_est(lavInspect(fit, what = "free"))
}

# Flatten a nested per-group lavaan.list (elements named `<group>.<block>`)
# to the flat `lavTech` shape, rebuilding plain block names in the same
# within-group order.
tsp_flatten_grouped_est <- function(li) {
  out <- unlist(li, recursive = FALSE)
  names(out) <- unlist(lapply(li, names), use.names = FALSE)
  out
}

# Per-group free-position matrices from
# `lavTech(what = "partable", list.by.group = TRUE)` (named list of lists;
# `$beta` etc. hold 1s at free positions). Not reconstructible from the
# partable df alone, so a missing `lavTech` is a hard error here (the
# installed build must be a 0.7-x with `lavTech`).
tsp_partable_mats <- function(fit) {
  if (tsp_has_lavTech()) {
    return(lavaan::lavTech(fit, what = "partable", list.by.group = TRUE))
  }
  stop(
    "lavaan ", as.character(utils::packageVersion("lavaan")),
    " does not provide lavTech(); R2spa requires the 0.7-x lavaan line"
    , call. = FALSE
  )
}

# Per-group observation counts (vector form; equal to
# `lavInspect(what = "nobs")`).
tsp_nobs <- function(fit) {
  tryCatch(
    unlist(fit@Data@nobs),
    error = function(e) lavInspect(fit, what = "nobs")
  )
}

# Group count (equal to `lavInspect(what = "ngroups")`).
tsp_ngroups <- function(fit) {
  tryCatch(
    fit@Data@ngroups,
    error = function(e) lavInspect(fit, what = "ngroups")
  )
}

# Per-group original observation counts (vector form; equal to
# `lavInspect(what = "norig")`).
tsp_norig <- function(fit) {
  tryCatch(
    unlist(fit@Data@norig),
    error = function(e) lavInspect(fit, what = "norig")
  )
}

# Optimizer convergence flag (scalar logical). MG fits optimize the combined
# free-parameter vector in a single call, so one flag covers all groups
# (verified 0.7-2, SG + MG).
tsp_converged <- function(fit) {
  tryCatch(
    fit@optim$converged,
    error = function(e) lavInspect(fit, what = "converged")
  )
}

# Overwrite the fitted p x p covariance of the free parameters: the lavaan
# @vcov slot is a list of {se = method string, information = string,
# vcov = matrix} and only the matrix element is touched, so vcov(), se(),
# and standardizedSolution() report the corrected values while the SE-method
# dispatch keys (@vcov$se, @Options$se) stay intact. Consumed by the
# corrected-SE path (R/tspa_corrected_se.R); canary-covered slot write.
tsp_set_vcov <- function(fit, m) {
  fit@vcov[["vcov"]] <- m
  fit
}
