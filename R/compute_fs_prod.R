# Product factor-score indicators (double-mean-centered latent interactions).
#
# compute_fs_prod() appends, to a single-group get_fs() result, the
# double-mean-centered (DMC) product indicator of each requested pair of
# distinct latents, together with its per-row standard error and implied
# loading. The derivation is in the roxygen @details of
# compute_fs_prod(); the pure-matrix helpers fs_prod_se2(),
# fs_prod_ecov(), and fs_prod_gamma() carry the SE, error-covariance, and
# implied-loading formulas for unit tests (fs_prod_ecov() is also consumed
# by tspa()). Per-row values are
# resolved through resolve_fs_per_row() (R/fs_indiv.R), so complete data
# (one block) and FIML (one block per observed-indicator pattern) share the
# same code path.

#' Compute product factor-score indicators (double-mean-centered)
#'
#' @description
#' `compute_fs_prod()` takes a single-group [get_fs()] result and appends,
#' for each requested pair of distinct latent variables, the
#' double-mean-centered (DMC) product indicator of the two factor scores
#' together with its per-row standard error and its implied loading. It is
#' the re-integrated successor of the quarantined `get_fs_int()` student
#' function: the column-naming conventions are kept (the product column
#' joins the two full score column names, so for latent names `a` and `b`
#' the columns are `fs_a:fs_b`, `fs_a:fs_b_se`, and `fs_a:fs_b_ld` — e.g.
#' the product of `fs_ind60` and `fs_dem60` is `fs_ind60:fs_dem60`), while
#' the standard error is computed from the stage-1 `fsL`/`fsT`/`psi`
#' attributes — the general joint-model formula below, instead of the
#' separate-single-factor special case the old function used.
#'
#' @details
#' Let \eqn{fs_k = L_k \xi + b_k + e_k} for \eqn{k = 1, \ldots, q}, where
#' \eqn{L} is the block's `fsL` (\eqn{q \times q}, rows = scores, cols =
#' latents), \eqn{T} is the block's `fsT` (score error covariance),
#' \eqn{\Psi} is the latent (co)variance (the `psi` attribute), and
#' \eqn{e} is the score error. The derivation relies on: \eqn{E[e] = 0}
#' (absorbed in \eqn{b}); \eqn{Cov(e, \xi) = 0} elementwise for ALL latents
#' (true for both regression and Bartlett scoring, where \eqn{e} is linear
#' in the indicator errors); and joint normality of \eqn{(\xi, e)}.
#'
#' For a pair of distinct latents \eqn{(a, b)}, the double-mean centered
#' (DMC) product indicator is
#'
#' \deqn{P = (fs_a - mi_a)(fs_b - mi_b) - mi_P}
#'
#' where \eqn{mi_a} and \eqn{mi_b} are the SAMPLE means of the score
#' columns and \eqn{mi_P} the sample mean of the (component-centered)
#' product. Its conditional expectation is
#'
#' \deqn{E[P \mid \xi] = (L_a (\xi - \alpha))(L_b (\xi - \alpha))}
#'
#' (the centering constants cancel in the error), so the measurement error
#' is
#'
#' \deqn{u = P - E[P \mid \xi] = L_a (\xi - \alpha) e_b + L_b (\xi - \alpha) e_a + e_a e_b}
#'
#' Under joint normality (Isserlis), for \eqn{\tau_k = L_k \Psi L_k'}
#' (a scalar: row \eqn{k} of \eqn{L} times \eqn{\Psi} times its
#' transpose), \eqn{\tau_{ab} = L_a \Psi L_b'}, \eqn{s_k^2 = T[k, k]}, and
#' \eqn{c = T[a, b]}:
#'
#' \deqn{se_P^2 = \tau_a s_b^2 + \tau_b s_a^2 + s_a^2 s_b^2 + c^2 + 2 \tau_{ab} c}
#'
#' (the \eqn{c^2} coefficient is 1, not 2: \eqn{Var(e_a e_b) = E[e_a^2 e_b^2] -
#' c^2 = (s_a^2 s_b^2 + 2 c^2) - c^2 = s_a^2 s_b^2 + c^2}).
#'
#' The implied loading is
#'
#' \deqn{\gamma = L[a, a] L[b, b] + L[a, b] L[b, a]}
#'
#' — the coefficient of \eqn{\xi_a \xi_b} in \eqn{E[P \mid \xi]} (Bartlett
#' joint model: \eqn{\gamma = 1}; separate single-factor models:
#' \eqn{\gamma = \lambda_a \lambda_b}).
#'
#' Sanity check: with diagonal \eqn{L} and \eqn{T} (separate single-factor
#' measurement models, \eqn{\lambda_k = L[k, k]}) the formula reduces to
#' the classic \eqn{\lambda_a^2 \Psi[a, a] s_b^2 + \lambda_b^2 \Psi[b, b] s_a^2 +
#' s_a^2 s_b^2} — the formula the old (quarantined) `get_fs_int()` had.
#' The terms \eqn{c^2} and \eqn{2 \tau_{ab} c} are exactly what the old
#' code dropped; they are nonzero in joint models (correlated factors,
#' cross-loadings, error covariances).
#'
#' Per-row values are pattern-resolved: under FIML each
#' observed-indicator pattern carries its own `fsL`/`fsT` block, so
#' `fs_a:fs_b_se` and `fs_a:fs_b_ld` take one value per pattern (the
#' implied loading varies by pattern because `fsL` does), while the latent
#' (co)variance `psi` is group-level and shared. For complete data the
#' standard error is constant across rows. Rows with no scorable pattern
#' (fully-missing rows) get `NA` in all three columns.
#'
#' @param fs A single-group [get_fs()] result: the unified data frame
#'        (`format = "unified"`, the default) or the single-group data
#'        frame of `format = "list"`. Rejected with an informative error
#'        (v1: single-group lavaan models only): multi-group results (a
#'        named list of per-group data frames, or a unified data frame
#'        with a group column), `merMod` results (3-D `fsT` attribute),
#'        and per-observation results (the `per_obs`/`mirt_per_obs`
#'        markers — local-mode FIML merges and `mirt` models).
#' @param product A character string of the form `"a:b + c:d"` (pairs of
#'        distinct latent names, `+`-separated) or a list of length-2
#'        character vectors of latent names (a 2-column matrix or data
#'        frame of names is coerced to a list). Same-factor products
#'        (`"x:x"`) and duplicated pairs are rejected (v1).
#' @return The input with, per requested pair in order, three appended
#'        columns: `fs_a:fs_b` the DMC product indicator, `fs_a:fs_b_se`
#'        the per-row standard error (constant for complete data; one
#'        value per observed-indicator pattern under FIML; `NA` for
#'        unscorable rows), and `fs_a:fs_b_ld` the per-row implied loading
#'        `gamma` (pattern-specific under FIML because `fsL` varies by
#'        pattern). All existing columns and attributes are untouched.
#' @export
#'
#' @references
#' The joint-normal moment expansion (Isserlis's theorem) used to derive
#' `se_P^2` and the measurement-error covariances is:
#'
#' Isserlis, L. (1918). On a formula for the product-moment coefficient of any
#' order of a normal frequency distribution in any number of variables.
#' \emph{Biometrika}, 12(1-2), 134-139. \doi{10.1093/biomet/12.1-2.134}.
#'
#' @examples
#' library(lavaan)
#' fit <- cfa(
#'   "ind60 =~ x1 + x2 + x3
#'    dem60 =~ y1 + y2 + y3 + y4",
#'   data = PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")],
#'   std.lv = TRUE
#' )
#' fs <- get_fs(fit, product = "ind60:dem60")
#' head(fs[c("fs_ind60", "fs_dem60", "fs_ind60:fs_dem60",
#'           "fs_ind60:fs_dem60_se", "fs_ind60:fs_dem60_ld")])
#'
#' # The same columns, computed directly from a get_fs() result
#' fs2 <- compute_fs_prod(get_fs(fit), product = "ind60:dem60")
#' head(fs2[c("fs_ind60:fs_dem60", "fs_ind60:fs_dem60_se")])
compute_fs_prod <- function(fs, product) {
  # v1 scope: single-group lavaan results only.
  if (is.list(fs) && !is.data.frame(fs)) {
    stop(
      "A multi-group list of per-group data frames is not supported ",
      "(v1: single-group only).",
      call. = FALSE
    )
  }
  T_attr <- attr(fs, "fsT")
  if (is.array(T_attr) && length(dim(T_attr)) == 3L) {
    stop("merMod results are not supported (v1: lavaan models only).",
         call. = FALSE)
  }
  if (isTRUE(attr(fs, "per_obs")) || isTRUE(attr(fs, "mirt_per_obs"))) {
    stop(
      "Per-observation results (local-mode FIML merges, mirt) are not ",
      "supported (v1: lavaan models only).",
      call. = FALSE
    )
  }
  # Unified multi-group output is a data frame with a group-column
  # attribute. A single-group list-format result with missing data also
  # carries a list `fsT` (length > 1), but that list is keyed by pattern
  # label and has no group-column attribute, so the multi-group signal is
  # the group_col attribute, not the list length alone.
  grp_col <- attr(fs, "group_col")
  if (is.list(T_attr) && length(T_attr) > 1L &&
      !is.null(grp_col) && grp_col %in% names(fs)) {
    stop("Multi-group results are not supported (v1: single-group only).",
         call. = FALSE)
  }
  resolved <- tryCatch(
    resolve_fs_per_row(fs),
    error = function(e) {
      stop(
        "Input is not a usable get_fs() result: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  # The latent (co)variance: a plain matrix (list-format single group) or
  # the length-1 list named "" (unified single group).
  psi <- fs_psi_matrix(attr(fs, "psi"))
  L1 <- resolved$blocks[[1L]]$fsL
  lv_names <- colnames(L1)
  if (is.null(lv_names) || anyNA(lv_names)) {
    stop(
      "The input's 'fsL' has no latent-variable column names; is it a ",
      "current get_fs() result?",
      call. = FALSE
    )
  }
  score_cols <- paste0("fs_", lv_names)
  missing_cols <- setdiff(score_cols, names(fs))
  if (length(missing_cols) > 0L) {
    stop(
      "The input is missing score column(s) ",
      paste(missing_cols, collapse = ", "),
      " (expected 'fs_<latent>' columns for the latents ",
      paste(lv_names, collapse = ", "), ").",
      call. = FALSE
    )
  }
  pairs <- parse_product_spec(product, lv_names)
  pidx <- resolved$pattern_idx
  out <- fs
  for (pair in pairs) {
    a <- pair[1L]
    b <- pair[2L]
    i <- match(a, lv_names)
    j <- match(b, lv_names)
    # Per-block values (one entry per observed-indicator pattern; the
    # all-NA block appended for fully-missing rows yields NA naturally).
    se2_by_block <- vapply(
      resolved$blocks,
      function(blk) fs_prod_se2(blk$fsL, blk$fsT, psi, i, j),
      numeric(1)
    )
    ld_by_block <- vapply(
      resolved$blocks,
      function(blk) fs_prod_gamma(blk$fsL, i, j),
      numeric(1)
    )
    se_col <- sqrt_or_na(se2_by_block[pidx])
    ld_col <- ld_by_block[pidx]
    # DMC product: component-center on the sample score means, then
    # center the product on its own sample mean (NA rows stay NA).
    x <- as.numeric(fs[[paste0("fs_", a)]])
    y <- as.numeric(fs[[paste0("fs_", b)]])
    p <- (x - mean(x, na.rm = TRUE)) * (y - mean(y, na.rm = TRUE))
    p <- p - mean(p, na.rm = TRUE)
    nm <- paste0("fs_", a, ":fs_", b)
    out[[nm]] <- p
    out[[paste0(nm, "_se")]] <- se_col
    out[[paste0(nm, "_ld")]] <- ld_col
  }
  out
}

# Parse the `product` argument into an ordered list of length-2 latent-name
# pairs. Accepts a character string of the form "a:b + c:d" (terms split
# on '+', each part trimmed, pairs split on ':') or a list of length-2
# character vectors (a 2-column matrix/data frame of names is coerced to
# a list). Errors on malformed terms (naming the offending term), unknown
# latent names (listing them), same-factor pairs, and duplicated pairs.
parse_product_spec <- function(product, lv_names) {
  if (is.null(product)) {
    stop(
      "'product' must be a character string of the form 'a:b + c:d' ",
      "(latent names) or a list of length-2 latent-name pairs.",
      call. = FALSE
    )
  }
  if (is.matrix(product) || is.data.frame(product)) {
    # Checked FIRST: a character matrix also satisfies is.character()
    # (typeof ignores dim attributes), and an array may satisfy
    # is.list(), so the matrix/data-frame form must win over both.
    if (ncol(product) != 2L) {
      stop(
        "A matrix or data frame 'product' must have 2 columns (one ",
        "latent name per column).",
        call. = FALSE
      )
    }
    pairs <- lapply(seq_len(nrow(product)), function(r) {
      vals <- vapply(
        seq_len(2L),
        function(k) as.character(product[r, k]),
        character(1L)
      )
      if (anyNA(vals)) {
        stop("A matrix or data frame 'product' must not contain NA.",
             call. = FALSE)
      }
      trimws(vals)
    })
  } else if (is.character(product)) {
    # A character vector of length > 1 is a sequence of 'a:b' terms
    # (the '+' separator is optional between elements).
    terms <- trimws(unlist(strsplit(product, "+", fixed = TRUE)))
    terms <- terms[terms != ""]
    if (length(terms) == 0L) {
      stop("'product' must name at least one 'a:b' pair.", call. = FALSE)
    }
    pairs <- lapply(terms, function(term) {
      parts <- trimws(strsplit(term, ":", fixed = TRUE)[[1L]])
      if (length(parts) != 2L) {
        stop(
          "Product term '", term, "' is not of the form 'a:b' (exactly ",
          "two latent names separated by ':').",
          call. = FALSE
        )
      }
      parts
    })
  } else if (is.list(product)) {
    pairs <- lapply(product, function(p) {
      if (!is.character(p) || length(p) != 2L) {
        stop(
          "Each element of 'product' must be a length-2 character vector ",
          "of latent names.",
          call. = FALSE
        )
      }
      trimws(p)
    })
  } else {
    stop(
      "'product' must be a character string of the form 'a:b + c:d' ",
      "(latent names) or a list of length-2 latent-name pairs.",
      call. = FALSE
    )
  }
  if (length(pairs) == 0L) {
    stop("'product' must name at least one 'a:b' pair.", call. = FALSE)
  }
  all_nm <- unlist(pairs, use.names = FALSE)
  unknown <- setdiff(all_nm, lv_names)
  if (length(unknown) > 0L) {
    stop(
      "Unknown latent name(s) in 'product': ",
      paste(unknown, collapse = ", "),
      ". The model's latent variables are: ",
      paste(lv_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  for (p in pairs) {
    if (p[1L] == p[2L]) {
      stop(
        "'", p[1L], ":", p[2L], "': same-factor products are not ",
        "supported (v1).",
        call. = FALSE
      )
    }
  }
  pair_lbl <- vapply(pairs, function(p) paste(p, collapse = ":"),
                     character(1L))
  # (a, b) and (b, a) are the same product: dedupe on the sorted pair.
  keys <- vapply(pairs, function(p) paste(sort(p), collapse = ":"),
                 character(1L))
  if (any(duplicated(keys))) {
    dup_lbl <- unique(pair_lbl[duplicated(keys) |
                                 duplicated(keys, fromLast = TRUE)])
    stop(
      "Duplicated pair(s) in 'product': ",
      paste(dup_lbl, collapse = ", "), ".",
      call. = FALSE
    )
  }
  pairs
}

# se_P^2 for the DMC product of scores i and j (row indices of the block's
# fsL, column indices of the shared psi):
#   tau_a s_b^2 + tau_b s_a^2 + s_a^2 s_b^2 + c^2 + 2 tau_ab c
# with tau_k = L_k psi L_k', tau_ab = L_a psi L_b', s_k^2 = T[k, k],
# c = T[a, b]. Plain matrices in, scalar out; an all-NA block yields NA
# (NA propagates through %*%).
fs_prod_se2 <- function(L, T, psi, i, j) {
  Li <- as.numeric(L[i, , drop = TRUE])
  Lj <- as.numeric(L[j, , drop = TRUE])
  tau_i <- as.numeric(Li %*% psi %*% Li)
  tau_j <- as.numeric(Lj %*% psi %*% Lj)
  tau_ij <- as.numeric(Li %*% psi %*% Lj)
  s_i2 <- T[i, i]
  s_j2 <- T[j, j]
  c_ij <- T[i, j]
  tau_i * s_j2 + tau_j * s_i2 + s_i2 * s_j2 + c_ij^2 + 2 * tau_ij * c_ij
}

# Cov(u_ij, u_kl) for the measurement errors of the DMC product indicators
# of the score pairs (i, j) and (k, l) (row indices of the block's fsL,
# column indices of the shared psi), where
#   u_kl = X_k e_l + X_l e_k + e_k e_l - c_kl,  X_k = L_k (xi - alpha),
# tau_uv = L_u psi L_v' (a scalar: row u of L times psi times row v'),
# c_uv = T[u, v] (so c_uu = T[u, u]), by the Isserlis expansion of the
# joint-normal (xi, e) moments (E[e] = 0, Cov(e, xi) = 0):
#   tau_ik c_jl + tau_il c_jk + tau_jk c_il + tau_jl c_ik
#     + c_ik c_jl + c_il c_jk
# Symmetric in the two pairs (i, j) <-> (k, l); the diagonal case
# (i, j) = (k, l) reduces exactly to fs_prod_se2(). Nonzero whenever the two
# products share a factor score (or, more generally, when a score-error
# moment links the two pairs) — the error covariance tspa() must fix
# between such product indicators in the stage-2 model. Plain matrices in,
# scalar out; an all-NA block yields NA (NA propagates through %*%).
fs_prod_ecov <- function(L, T, psi, i, j, k, l) {
  Li <- as.numeric(L[i, , drop = TRUE])
  Lj <- as.numeric(L[j, , drop = TRUE])
  Lk <- as.numeric(L[k, , drop = TRUE])
  Ll <- as.numeric(L[l, , drop = TRUE])
  tau_ik <- as.numeric(Li %*% psi %*% Lk)
  tau_il <- as.numeric(Li %*% psi %*% Ll)
  tau_jk <- as.numeric(Lj %*% psi %*% Lk)
  tau_jl <- as.numeric(Lj %*% psi %*% Ll)
  c_ik <- T[i, k]
  c_il <- T[i, l]
  c_jk <- T[j, k]
  c_jl <- T[j, l]
  tau_ik * c_jl + tau_il * c_jk + tau_jk * c_il + tau_jl * c_ik +
    c_ik * c_jl + c_il * c_jk
}

# Implied loading: gamma = L[i, i] L[j, j] + L[i, j] L[j, i], the
# coefficient of xi_i xi_j in E[P | xi] = (L_i (xi - alpha))(L_j (xi - alpha)).
fs_prod_gamma <- function(L, i, j) {
  L[i, i] * L[j, j] + L[i, j] * L[j, i]
}

# Unwrap a get_fs() result's 'psi' attribute to the single-group latent
# (co)variance matrix: a plain matrix (list-format single group) or the
# length-1 list named "" (unified single group). Errors informatively on
# anything else (multi-group per-group lists, missing attribute, ...).
fs_psi_matrix <- function(psi_attr) {
  if (is.matrix(psi_attr)) {
    return(psi_attr)
  }
  if (is.list(psi_attr) && length(psi_attr) == 1L &&
      identical(names(psi_attr), "") && is.matrix(psi_attr[[1L]])) {
    return(psi_attr[[1L]])
  }
  stop(
    "The input carries no usable single-group 'psi' attribute (latent ",
    "covariance); is it a current single-group get_fs() result?",
    call. = FALSE
  )
}
