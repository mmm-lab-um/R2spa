# Individual-specific (per-row) factor-score definition quantities.
#
# fs_indiv() expands the per-block fsL/fsT/fsb attributes of a get_fs()
# result into one long data frame (one row per input row), reusing
# get_fs() column naming. Row resolution (resolve_fs_per_row()) dispatches
# on the input shape: lavaan "unified" (single data frame with list-valued
# attributes), lavaan "list" (named list of per-group data frames with
# direct attributes), and merMod (one row per cluster with 3-D
# per-cluster arrays).
#
# All per-row consumers -- augment_fs() (R/get_fscore.R), augment_fs2()
# (R/get_fscore_math.R), fs_indiv(), and the mirt per-obs paths
# (R/get_fs_methods.R) -- source their per-row se/loadings/error values
# from the shared value-only helper fs_row_cols() defined below, so there
# is one source of truth for those quantities. The r2spa column naming is
# shared through fs_row_colnames(), the naming twin of fs_row_cols();
# augment_fs2() keeps its own legacy se_*/upper-tri naming on top.

#' Individual-specific (per-row) factor-score definition quantities
#'
#' @description
#' `fs_indiv()` converts a [get_fs()] result into a single long data frame
#' (one row per input row) carrying, for every observation, the
#' individual-specific standard errors, implied loadings, and error
#' variances and covariances of the factor scores, reusing `get_fs()`'s
#' column naming. In shape this resembles `mirt::fscores(full.scores.SE =
#' TRUE)` (a per-observation table pairing each score with its standard
#' errors), but the quantities are read off the `fsL`/`fsT`/`fsb` attributes
#' of a `get_fs()` result rather than estimated from a new model.
#'
#' For lavaan models the per-row values are pattern-resolved: an
#' observation's `fsL`/`fsT`/`fsb` depend only on its observed-indicator
#' pattern, not on its response values, so rows within a group differ only
#' when the data contain missing values (rows across groups may always
#' differ). For `merMod` models there is one row per cluster, each carrying
#' that cluster's own `fsL`/`fsT`.
#'
#' @param fs A factor-score object as returned by [get_fs()] (or
#'        [get_fs_lavaan()]/[get_fs_lmer()]): a unified data frame, a named
#'        list of per-group data frames (`format = "list"`), or a `merMod`
#'        result (one row per cluster).
#' @param include_intercept Logical. When `TRUE`, also emit the factor
#'        score intercepts (the `fsb` attribute) as `int_fs_<f>` columns
#'        (`q` of them). Ignored for `merMod` results, which have no `fsb`
#'        (random effects are mean zero).
#' @param ... Currently unused.
#' @return A data frame with `nrow()` equal to the number of rows of the
#'        input [get_fs()] result (one row per observation for lavaan
#'        models, one row per cluster for `merMod` models). Columns, in
#'        order:
#'        * the input's factor-score columns, kept unmodified (`fs_<f>`; for
#'          `merMod` models `fs_u<k>`, or the legacy `u<k>_eb` name when
#'          [get_fs.merMod()] was run with `legacy_names = TRUE`);
#'        * the per-observation standard errors (`fs_<f>_se`; the square
#'          root of the per-row `fsT` diagonal, `NA` where that entry is
#'          negative or non-finite);
#'        * the implied loadings of the latents on the factor scores
#'          (`<latent_j>_by_fs_<f>`, `q^2` of them);
#'        * the error variances and error covariances of the factor scores
#'          (`ev_fs_<f>`, `q` of them, and `ecov_fs_<a>_fs_<b>`, `q choose
#'          2` of them), in the same lower-triangular order [get_fs()] uses;
#'        * optionally (when `include_intercept = TRUE`) the score
#'          intercepts (`int_fs_<f>`, `q` of them);
#'        * a trailing `group` column for multi-group lavaan results (the
#'          input's own group-column name in `"unified"` format, or `group`
#'          in `"list"` format); and
#'        * a trailing `id` column holding the cluster/subject id for `merMod`
#'          results.
#' @export
#' @examples
#' library(lavaan)
#' fit <- cfa("visual =~ x1 + x2 + x3", data = HolzingerSwineford1939)
#' # For a lavaan result the per-row se / loading / ev columns are the same
#' # ones get_fs() already carries; include_intercept = TRUE adds the score
#' # intercepts (int_fs_*) from the fsb attribute (non-NULL under mean
#' # scoring or with prior_mean/prior_cov).
#' fs_indiv(get_fs(fit, method = "mean"), include_intercept = TRUE)
#'
#' # merMod: one row per cluster, with a trailing id column holding the
#' # cluster (Subject) level
#' library(lme4)
#' lmod <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' fs_indiv(get_fs(lmod))
fs_indiv <- function(fs, include_intercept = FALSE, ...) {
  resolved <- resolve_fs_per_row(fs)
  ref_T <- resolved$blocks[[1L]]$fsT
  ref_L <- resolved$blocks[[1L]]$fsL
  q <- ncol(ref_T)
  has_int <- include_intercept && !is.null(resolved$blocks[[1L]]$fsb)

  n <- resolved$n
  k_ld <- q * q
  k_ev <- q * (q + 1L) / 2
  k_int <- if (has_int) q else 0L

  se_mat <- matrix(NA_real_, nrow = n, ncol = q)
  ld_mat <- matrix(NA_real_, nrow = n, ncol = k_ld)
  ev_mat <- matrix(NA_real_, nrow = n, ncol = k_ev)
  int_mat <- if (has_int) matrix(NA_real_, nrow = n, ncol = k_int) else NULL

  # Row indices per block, precomputed once: `which(pattern_idx == b)` inside
  # the loop is O(n x n_blocks) -- quadratic for per-cluster (merMod) and
  # per-row (mirt per-obs) inputs, where every row is its own block. split()
  # is O(n). The explicit factor levels keep one (possibly empty) entry per
  # block, so block b always maps to rows_by_block[[b]].
  rows_by_block <- split(
    seq_len(n),
    factor(resolved$pattern_idx, levels = seq_along(resolved$blocks))
  )

  for (b in seq_along(resolved$blocks)) {
    blk <- resolved$blocks[[b]]
    rows_b <- rows_by_block[[b]]
    vals <- fs_row_cols(
      resolved$scores[rows_b, , drop = FALSE],
      blk$fsL,
      blk$fsT,
      if (has_int) blk$fsb else NULL
    )
    # NOTE: 'drop' is not a valid argument in an assignment sub-expression
    # (extraction only); row-slice assignment below uses the default.
    se_mat[rows_b, ] <- vals[, seq_len(q), drop = FALSE]
    ld_mat[rows_b, ] <-
      vals[, (q + 1L):(q + k_ld), drop = FALSE]
    ev_mat[rows_b, ] <-
      vals[, (q + k_ld + 1L):(q + k_ld + k_ev), drop = FALSE]
    if (has_int) {
      int_mat[rows_b, ] <-
        vals[, (q + k_ld + k_ev + 1L):ncol(vals), drop = FALSE]
    }
  }

  # Column naming: the r2spa convention, shared with augment_fs() and the
  # mirt per-obs paths via fs_row_colnames() (the naming twin of
  # fs_row_cols()).
  fs_names <- rownames(ref_L)   # "fs_<factor>" (merMod: "fs_u<k>" / legacy)
  nm <- fs_row_colnames(ref_L, ref_T)
  se_nm <- nm$se
  ld_nm <- nm$ld
  ev_nm <- nm$ev
  int_nm <- if (has_int) paste0("int_", fs_names) else character(0)

  cols <- list(resolved$scores, se_mat, ld_mat, ev_mat)
  if (has_int) {
    cols[[length(cols) + 1L]] <- int_mat
  }
  out <- do.call(cbind, cols)
  colnames(out) <- c(
    colnames(resolved$scores),
    se_nm,
    ld_nm,
    ev_nm,
    int_nm
  )
  # Legacy u<k>_eb-style input (get_fs.merMod(legacy_names = TRUE)):
  # translate the new column names with the same (idempotent) rules the
  # original output used.
  if (resolved$legacy) {
    colnames(out) <- rename_legacy_fs_cols(colnames(out))
  }
  if (!is.null(resolved$group_vals)) {
    out[[resolved$group_col]] <- resolved$group_vals
  }
  if (!is.null(resolved$id_vals)) {
    out[["id"]] <- resolved$id_vals
  }
  out
}

# Value-only shared engine for per-row se / loadings / error terms.
#
# Returns an n x K numeric matrix with K column blocks, all computed from
# the (pattern-level) fsL/fsT/fsb and repeated for the `nrow(fs)` rows:
#   1. se    q columns:            sqrt_or_na(diag(fsT))
#   2. lds   q^2 columns:          c(fsL) (column-major: per latent,
#                                    loadings of that latent on each score)
#   3. evs   q*(q+1)/2 columns:    fsT[i, j] in i-outer / j<=i-inner order
#                                    (row-major lower triangle, the order
#                                    get_fs() columns are named in)
#   4. int   q columns:            fsb (only when fsb is non-NULL)
# The matrix carries no names and no attributes: every consumer applies
# its own column naming/layout (r2spa naming in augment_fs()/fs_indiv(),
# se_*/int_*/upper-tri legacy naming in augment_fs2()/augment_lav_publish()).
fs_row_cols <- function(fs, fsL, fsT, fsb = NULL) {
  n <- nrow(fs)
  q <- nrow(fsT)
  se <- unname(sqrt_or_na(diag(fsT)))
  lds <- unname(c(as.matrix(fsL)))
  # i-outer / j<=i-inner (row-major lower triangle) -- the order the
  # r2spa ev/ecov columns are named in (and the explicit fs_ev[i, j]
  # order of get_fs()); fsT[lower.tri(fsT, diag = TRUE)] would give
  # COLUMN-major triangle order, which differs for q >= 3.
  evs <- rep(NA_real_, q * (q + 1L) / 2)
  count <- 1L
  for (i in seq_len(q)) {
    for (j in seq_len(i)) {
      evs[count] <- fsT[i, j]
      count <- count + 1L
    }
  }
  row_vals <- c(se, lds, evs)
  if (!is.null(fsb)) {
    row_vals <- c(row_vals, unname(as.numeric(fsb)))
  }
  # nrow/ncol must both be given: with only nrow, matrix() requires the data
  # length to be a multiple of n, which holds only for single-row patterns.
  matrix(row_vals, nrow = n, ncol = length(row_vals), byrow = TRUE)
}

# The r2spa per-row column NAMES, the naming twin of fs_row_cols()'s value
# layout. Returns list(se, ld, ev):
#   se:  <score>_se, one per score (row order of fsT)
#   ld:  <latent>_by_<score>, q^2 of them, column-major per latent
#        (matching c(as.matrix(fsL)))
#   ev:  ev_<score> / ecov_<score_i>_<score_j>, q*(q+1)/2 of them,
#        i-outer / j<=i-inner on fsT (row-major lower triangle)
# Shared by every r2spa-named consumer -- augment_fs() (R/get_fscore.R),
# fs_indiv() above, and the mirt per-obs paths (R/get_fs_methods.R) -- so
# the naming scheme has exactly one implementation.
fs_row_colnames <- function(fsL, fsT) {
  lv_names <- colnames(fsL)
  # Canonical score names derived from the latent names (the same rule
  # score_column_names() uses to locate the score columns); equal to the
  # fsL rownames for every real get_fs() result.
  ld_fs_names <- paste0("fs_", lv_names)
  fs_row <- rownames(fsT)
  fs_col <- colnames(fsT)
  q <- ncol(fsT)
  se_nm <- paste0(fs_row, "_se")
  ld_nm <- character(q * q)
  count <- 1L
  for (j in seq_len(q)) {
    ld_nm[count:(count + q - 1L)] <- paste(lv_names[j], ld_fs_names,
                                           sep = "_by_")
    count <- count + q
  }
  ev_nm <- character(q * (q + 1L) / 2L)
  count <- 1L
  for (i in seq_len(q)) {
    for (j in seq_len(i)) {
      ev_nm[count] <- if (i == j) {
        paste0("ev_", fs_row[i])
      } else {
        paste0("ecov_", fs_row[i], "_", fs_col[j])
      }
      count <- count + 1L
    }
  }
  list(se = se_nm, ld = ld_nm, ev = ev_nm)
}

# Resolve a get_fs() result into (a) the score data frame, (b) one
# (fsL, fsT, fsb) value block per pattern (lavaan) or per cluster
# (merMod), and (c) an integer vector mapping each row to its block.
# Dispatches on the input shape:
#   - named list of data frames        -> lavaan "list" format (per-group
#     attributes sit directly on each group data frame)
#   - data frame with 3-D fsT/fsL      -> merMod (one row per cluster)
#   - data frame otherwise             -> lavaan "unified" (or
#     "list"-format single group: direct matrix/vector attributes)
resolve_fs_per_row <- function(fs) {
  # Per-observation input: the `mirt_per_obs` (mirt) or `per_obs` (local
  # FIML merge, PLAN 14) marker is authoritative and is checked FIRST.
  # Each row carries its own fsL/fsT, so one block is minted per row. The
  # markers are never set on a lavaan/merMod result, so the branches below
  # are byte-identical to before.
  if (is.data.frame(fs) &&
      (isTRUE(attr(fs, "mirt_per_obs")) || isTRUE(attr(fs, "per_obs")))) {
    return(resolve_per_obs(fs))
  }
  if (is.list(fs) && !is.data.frame(fs)) {
    return(resolve_lavaan_list(fs))
  }
  if (!is.data.frame(fs)) {
    stop(
      "'fs' must be a data frame or a named list of data frames as ",
      "returned by get_fs().",
      call. = FALSE
    )
  }
  T_attr <- attr(fs, "fsT")
  if (is.array(T_attr) && length(dim(T_attr)) == 3L) {
    return(resolve_mer_mod(fs))
  }
  resolve_lavaan_unified(fs)
}

# mirt per-observation result: one block per row (each row's own fsL/fsT),
# pattern identity = the row itself. The shared fs_indiv() block loop then
# emits per-row se / loadings / error terms from these per-row blocks. The
# completely-missing rows were minted by get_fs.SingleGroupClass() as all-NA
# blocks, so fs_row_cols() yields NA se/loadings/error terms for them (the
# NA-row convention), with scores already NA in the carried-through score
# columns.
resolve_per_obs <- function(fs) {
  L_attr <- attr(fs, "fsL")
  T_attr <- attr(fs, "fsT")
  b_attr <- attr(fs, "fsb")
  if (is.null(L_attr) || is.null(T_attr) || !is.list(L_attr) || !is.list(T_attr)) {
    stop(
      "'mirt_per_obs' input is missing its per-row 'fsL'/'fsT' list ",
      "attributes; is it a get_fs() (mirt) result?",
      call. = FALSE
    )
  }
  n <- nrow(fs)
  if (length(L_attr) != n || length(T_attr) != n) {
    stop(
      "Per-row 'fsL'/'fsT' lists (length ", length(L_attr), " / ",
      length(T_attr), ") do not match the number of rows (", n, ").",
      call. = FALSE
    )
  }
  ref_L <- L_attr[[1L]]
  # Per-row intercepts: get_fs() (mirt) attaches fsb as a per-row list (each
  # element a q-vector, or all-NA for a completely-missing row). Fall back to a
  # shared constant vector for the (legacy) single-value case.
  b_is_list <- is.list(b_attr)
  if (b_is_list && length(b_attr) != n) {
    stop(
      "Per-row 'fsb' list (length ", length(b_attr), ") does not match the ",
      "number of rows (", n, ").",
      call. = FALSE
    )
  }
  # Completely-missing rows were minted with an all-NA fsT; give them an
  # all-NA intercept too (the lavaan NA-row convention: a row with no scorable
  # pattern has NA score/SE/ev AND NA intercept, not the shared constant).
  b_na <- if (is.null(b_attr)) {
    NULL
  } else {
    rep(NA_real_, if (b_is_list) length(b_attr[[1L]]) else length(b_attr))
  }
  blocks <- lapply(seq_len(n), function(i) {
    Ti <- T_attr[[i]]
    if (all(is.na(Ti))) {
      fsb_i <- b_na
    } else if (b_is_list) {
      fsb_i <- b_attr[[i]]
    } else {
      fsb_i <- b_attr
    }
    list(
      fsL = L_attr[[i]],
      fsT = Ti,
      fsb = fsb_i
    )
  })
  # A per-row result may carry a `group_col` attribute naming its group
  # column (local FIML multi-group output, PLAN 14); mirt per-obs results
  # carry none (their group column, if any, is literally named "group" and
  # is recovered downstream in pool_per_unit()), so the behavior there is
  # unchanged.
  g_col <- attr(fs, "group_col")
  if (!is.null(g_col) && !(g_col %in% names(fs))) {
    g_col <- NULL
  }
  g_vals <- if (is.null(g_col)) {
    NULL
  } else {
    as.character(fs[[g_col]])
  }
  make_resolved(
    fs_scores_df(fs, ref_L),
    n, seq_len(n), blocks,
    group_col = g_col, group_vals = g_vals, id_vals = NULL
  )
}

# lavaan "unified" output (and the "list"-format single-group equivalent,
# whose attributes sit directly on the data frame).
resolve_lavaan_unified <- function(fs) {
  T_attr <- attr(fs, "fsT")
  fp <- attr(fs, "fs_pattern")
  if (is.null(fp)) {
    stop(
      "Input data frame has no 'fs_pattern' attribute; ",
      "is it a get_fs() result?",
      call. = FALSE
    )
  }
  grp_col <- attr(fs, "group_col")
  has_group <- !is.null(grp_col) && grp_col %in% names(fs)
  n <- nrow(fs)

  # Unified single group: the attribute values are wrapped in a length-1
  # list named "" (a per-pattern list is never a length-1 list named "" --
  # a single pattern stays a plain matrix, and pattern labels are
  # non-empty joins of indicator names).
  unified_sg <- is.list(T_attr) && length(T_attr) == 1L &&
    !is.null(names(T_attr)) && names(T_attr) == "" &&
    is.list(fp) && length(fp) == 1L &&
    !is.null(names(fp)) && names(fp) == ""

  if (unified_sg) {
    r <- resolve_group_blocks(
      fp[[1L]]$label,
      T_attr[[1L]],
      attr(fs, "fsL")[[1L]],
      attr(fs, "fsb")[[1L]]
    )
    # Name the score columns from a single-pattern fsL (the per-group
    # attribute is a per-pattern LIST when the fit has missing data).
    scores_df <- fs_scores_df(fs, r$blocks[[1L]]$fsL)
    make_resolved(scores_df, n, r$pattern_idx, r$blocks,
                  group_col = NULL, group_vals = NULL)
  } else if (has_group) {
    # Unified multi-group: attributes are named lists keyed by group
    # label; the group column carries the same labels in group order.
    group_labels <- unique(fs[[grp_col]])
    blocks <- list()
    pattern_idx <- integer(n)
    group_vals <- character(n)
    offset <- 0L
    L_ref <- NULL
    for (g in group_labels) {
      rows_g <- which(fs[[grp_col]] == g)
      r <- resolve_group_blocks(
        fp[[g]]$label,
        T_attr[[g]],
        attr(fs, "fsL")[[g]],
        attr(fs, "fsb")[[g]]
      )
      blocks <- c(blocks, r$blocks)
      pattern_idx[rows_g] <- r$pattern_idx + offset
      offset <- offset + length(r$blocks)
      group_vals[rows_g] <- g
      if (is.null(L_ref)) {
        L_ref <- r$blocks[[1L]]$fsL
      }
    }
    make_resolved(fs_scores_df(fs, L_ref),
                  n, pattern_idx, blocks,
                  group_col = grp_col, group_vals = group_vals)
  } else {
    # "list"-format single group: a plain data frame whose fsT/fsL/fsb
    # attributes are the pattern values directly (or a named list of them
    # when the fit has missing data).
    r <- resolve_group_blocks(
      fp$label,
      T_attr,
      attr(fs, "fsL"),
      attr(fs, "fsb")
    )
    make_resolved(fs_scores_df(fs, first_pattern_value(attr(fs, "fsL"))),
                  n, r$pattern_idx, r$blocks,
                  group_col = NULL, group_vals = NULL)
  }
}

# A group-level fsL/fsT/fsb attribute is a plain matrix/vector for a
# single observed pattern and a named list of them (one per pattern) for
# multi-pattern fits. Column naming only needs the model structure, which
# is identical across patterns, so the first pattern's value suffices.
first_pattern_value <- function(x) {
  if (is.list(x) && !is.matrix(x) && length(x) > 0) x[[1L]] else x
}

# lavaan "list" format: a named list of per-group data frames with
# per-group attributes sitting directly on each element.
resolve_lavaan_list <- function(fs) {
  group_labels <- names(fs)
  if (is.null(group_labels) || !all(nzchar(group_labels))) {
    stop(
      "List input must be named by group label (get_fs(format = 'list')).",
      call. = FALSE
    )
  }
  n <- sum(vapply(fs, nrow, integer(1)))
  first_L1 <- first_pattern_value(attr(fs[[group_labels[1L]]], "fsL"))
  nm_l1 <- score_column_names(first_L1, colnames(fs[[group_labels[1L]]]))
  blocks <- list()
  pattern_idx <- integer(n)
  group_vals <- character(n)
  offset_row <- 0L
  offset_blk <- 0L
  for (g in group_labels) {
    df_g <- fs[[g]]
    if (!is.data.frame(df_g)) {
      stop(
        "List element '", g, "' is not a data frame.",
        call. = FALSE
      )
    }
    fp_g <- attr(df_g, "fs_pattern")
    if (is.null(fp_g)) {
      stop(
        "Group '", g, "' has no 'fs_pattern' attribute; ",
        "is the input a get_fs() result?",
        call. = FALSE
      )
    }
    if (!all(nm_l1$score_nm %in% colnames(df_g))) {
      stop(
        "Group '", g, "' is missing the score column(s) ",
        paste(nm_l1$score_nm[nm_l1$score_nm %notin% colnames(df_g)],
              collapse = ", "), ".",
        call. = FALSE
      )
    }
    r <- resolve_group_blocks(
      fp_g$label,
      attr(df_g, "fsT"),
      attr(df_g, "fsL"),
      attr(df_g, "fsb")
    )
    rows_g <- offset_row + seq_along(fp_g$label)
    blocks <- c(blocks, r$blocks)
    pattern_idx[rows_g] <- r$pattern_idx + offset_blk
    offset_blk <- offset_blk + length(r$blocks)
    offset_row <- offset_row + nrow(df_g)
    group_vals[rows_g] <- g
  }
  scores_df <- do.call(
    rbind,
    lapply(fs, function(df_g) df_g[, nm_l1$score_nm, drop = FALSE])
  )
  make_resolved_with_legacy(scores_df, n, pattern_idx, blocks,
                            nm_l1$legacy,
                            group_col = "group", group_vals = group_vals)
}

# merMod result: one row per cluster; per-cluster fsL/fsT come from the
# 3-D array attributes (dim 3 = cluster, in canonical factor-level order,
# matching the row order).
resolve_mer_mod <- function(fs) {
  T_attr <- attr(fs, "fsT")
  L_attr <- attr(fs, "fsL")
  if (!is.array(L_attr) || length(dim(L_attr)) != 3L) {
    stop(
      "3-D 'fsT' attribute without a matching 3-D 'fsL'; ",
      "is the input a get_fs() (merMod) result?",
      call. = FALSE
    )
  }
  b_attr <- attr(fs, "fsb")
  n_clus <- dim(T_attr)[3L]
  cluster_names <- dimnames(T_attr)[[3L]]
  if (nrow(fs) != n_clus || length(cluster_names) != n_clus) {
    stop(
      "merMod input rows (", nrow(fs), ") do not match the per-cluster ",
      "attributes (", n_clus, " clusters).",
      call. = FALSE
    )
  }
  blocks <- vector("list", n_clus)
  for (j in seq_len(n_clus)) {
    # matrix() enforces 2-D geometry on the subscripted slice: the default
    # drop = TRUE turns a q = 1 slice of a 1 x 1 x N attribute into a bare
    # scalar (ncol() NULL), which breaks ncol()/diag() downstream. (With
    # drop = FALSE a length-1 trailing dim would survive as a 3-D q x q x 1
    # "array", which breaks the same ops.)
    blocks[[j]] <- list(
      fsL = matrix(L_attr[,, j], nrow = dim(L_attr)[1L],
                   dimnames = dimnames(L_attr)[1:2]),
      fsT = matrix(T_attr[,, j], nrow = dim(T_attr)[1L],
                   dimnames = dimnames(T_attr)[1:2]),
      fsb = if (is.null(b_attr)) {
        NULL
      } else if (is.list(b_attr) && length(b_attr) == n_clus) {
        b_attr[[j]]
      } else {
        b_attr
      }
    )
  }
  nm_l1 <- score_column_names(blocks[[1L]]$fsL, colnames(fs))
  make_resolved_with_legacy(fs[, nm_l1$score_nm, drop = FALSE],
                            n_clus, seq_len(n_clus), blocks, nm_l1$legacy,
                            group_col = NULL, group_vals = NULL,
                            id_vals = cluster_names)
}

# One group: given the per-case pattern labels and the group's fsT/fsL/fsb
# attribute values (each either a plain matrix/vector for the single
# pattern, or a named list of them with one entry per pattern label),
# build the value blocks and the row -> block index.
resolve_group_blocks <- function(labels_g, T_g, L_g, b_g) {
  blocks <- list()
  pattern_idx <- integer(length(labels_g))
  if (is.list(T_g)) {
    pat_labels <- names(T_g)
    m <- match(labels_g, pat_labels)
    # Rows whose observed-indicator pattern is entirely missing are kept by
    # get_fs() with NA score/SE/ev, but they have no entry in the per-pattern
    # fsL/fsT/fsb attributes (lavaan cannot score a fully-NA pattern). Preserve
    # them in the long table with all-NA values by routing those rows to a
    # dedicated all-NA block appended after the real (named) blocks.
    na_rows <- is.na(m)
    for (p in seq_along(pat_labels)) {
      rows_p <- which(m == p)
      blocks[[length(blocks) + 1L]] <- list(
        fsL = L_g[[p]],
        fsT = T_g[[p]],
        fsb = b_g[[p]]
      )
      pattern_idx[rows_p] <- length(blocks)
    }
    if (any(na_rows)) {
      ref_L <- L_g[[1L]]
      ref_T <- T_g[[1L]]
      blocks[[length(blocks) + 1L]] <- list(
        fsL = matrix(NA_real_, nrow(ref_L), ncol(ref_L),
                     dimnames = dimnames(ref_L)),
        fsT = matrix(NA_real_, nrow(ref_T), ncol(ref_T),
                     dimnames = dimnames(ref_T)),
        fsb = if (is.null(b_g) || is.null(b_g[[1L]])) {
          NULL
        } else {
          rep(NA_real_, length(b_g[[1L]]))
        }
      )
      pattern_idx[na_rows] <- length(blocks)
    }
  } else {
    blocks[[1L]] <- list(fsL = L_g, fsT = T_g, fsb = b_g)
    pattern_idx[] <- 1L
    # A single observed pattern can still coexist with fully-NA (empty) rows:
    # lavaan drops an all-missing case from the pattern set but get_fs() keeps
    # it with an NA score and a NA pattern label. Route those rows to a
    # dedicated all-NA block (the same convention as the multi-pattern list
    # branch above) instead of the single value block.
    na_rows <- is.na(labels_g)
    if (any(na_rows)) {
      blocks[[length(blocks) + 1L]] <- list(
        fsL = matrix(NA_real_, nrow(L_g), ncol(L_g),
                     dimnames = dimnames(L_g)),
        fsT = matrix(NA_real_, nrow(T_g), ncol(T_g),
                     dimnames = dimnames(T_g)),
        fsb = if (is.null(b_g) || length(b_g) == 0L) {
          NULL
        } else {
          rep(NA_real_, length(b_g))
        }
      )
      pattern_idx[na_rows] <- length(blocks)
    }
  }
  list(blocks = blocks, pattern_idx = pattern_idx)
}

# Score column names from a block's fsL: canonical "fs_<name>" first
# (lavaan and non-legacy merMod output); bare "u<k>[_eb]" names second
# (legacy_names = TRUE merMod output).
score_column_names <- function(L1, available) {
  canon_nm <- paste0("fs_", colnames(L1))
  legacy_nm <- colnames(L1)
  if (all(canon_nm %in% available)) {
    return(list(score_nm = canon_nm, legacy = FALSE))
  }
  if (all(legacy_nm %in% available)) {
    return(list(score_nm = legacy_nm, legacy = TRUE))
  }
  stop(
    "Could not locate the factor score columns in the input (expected: ",
    paste(canon_nm, collapse = ", "), " or ",
    paste(legacy_nm, collapse = ", "), ").",
    call. = FALSE
  )
}
`%notin%` <- Negate(`%in%`)

# Score data frame for a data-frame input.
fs_scores_df <- function(fs, L1) {
  fs[, score_column_names(L1, colnames(fs))$score_nm, drop = FALSE]
}

make_resolved <- function(scores_df, n, pattern_idx, blocks, group_col,
                          group_vals, id_vals = NULL) {
  make_resolved_with_legacy(
    scores_df, n, pattern_idx, blocks,
    legacy = FALSE,
    group_col = group_col,
    group_vals = group_vals,
    id_vals = id_vals
  )
}

# Assemble the resolved structure.
make_resolved_with_legacy <- function(scores_df, n, pattern_idx, blocks,
                                      legacy, group_col, group_vals,
                                      id_vals = NULL) {
  list(
    n = n,
    scores = scores_df,
    pattern_idx = pattern_idx,
    blocks = blocks,
    group_col = group_col,
    group_vals = group_vals,
    id_vals = id_vals,
    legacy = legacy
  )
}
