#' Get Factor Scores and the Corresponding Standard Error of Measurement
#'
#' @description
#' `get_fs()` is an S3 generic that extracts factor scores from fitted models.
#' Methods are available for `data.frame` (fits a CFA internally), `lavaan`
#' objects, `lmerMod` objects, and fitted `mirt` models (single-group
#' `SingleGroupClass` and multi-group `MultipleGroupClass`; `mirt` is a
#' `Suggests` dependency). Multi-group mirt results carry a trailing `group`
#' column and a per-group (`list`) `psi` attribute.
#'
#' @details
#' When `object` is a data frame and `model` is supplied as a lavaan syntax string,
#' the function internally calls [lavaan::cfa()] and then dispatches to the
#' `lavaan` method. When `object` is a fitted model object, the appropriate S3
#' method is called directly.
#'
#' `get_fs()` replaced `get_fs_lavaan()` and `get_fs_lmer()`, which are now
#' thin wrappers retained for backward compatibility.
#'
#' @param object A data frame, a fitted [lavaan] model object, or a fitted
#'        [lme4::lmer] model object (`merMod`).
#' @param model An optional string specifying the measurement model
#'              in \code{lavaan} syntax. Only used when `object` is a data frame.
#'              See \code{\link[lavaan]{model.syntax}} for more information.
#' @param group Character. Name of the grouping variable for multiple group
#'              analysis, which is passed to \code{\link[lavaan]{cfa}}.
#'              Only used when `object` is a data frame.
#' @param method Character. Method for computing factor scores. For
#'               `lavaan` and data frame objects: `"regression"` (default,
#'               consistent with \code{\link[lavaan]{lavPredict}}),
#'               `"Bartlett"`, or `"mean"` (a third, distinct method: sum
#'               scores, each score being the plain uncentered mean of the
#'               items assigned to its factor, using no latent
#'               distribution), with `"ML"` an alias for `"Bartlett"` and
#'               `"EB"` an alias for `"regression"`. For `merMod` objects:
#'               `"EB"` (empirical Bayes, default; identical to the
#'               first random-effect term's \code{\link[lme4]{ranef}}()
#'               estimates) or `"ML"` (a prior-free,
#'               per-cluster OLS estimate of the random effects, using no
#'               random-effects prior, analogous to Bartlett scores for
#'               `lavaan` objects). The `"ML"`/`"EB"` aliases apply to the
#'               lavaan path only; for `merMod` objects the two strings are
#'               distinct methods. Bartlett scores have more desirable
#'               properties than regression scores and may be preferred for
#'               2S-PA. `method = "mean"` takes the item-to-factor
#'               assignment from `sum_items` (auto-derived from the
#'               estimated loadings when `NULL`); it errors when the model
#'               was fitted with missing data retained (e.g. FIML/digamma),
#'               and is not supported together with `corrected_fsT`,
#'               `vfsLT`, `reliability`, `prior_mean`, or `prior_cov`.
#' @param corrected_fsT Logical. Whether to correct for the sampling
#'                      error in the factor score weights when computing
#'                      the error variance estimates of factor scores.
#'                      Currently ignored for `merMod` objects.
#' @param vfsLT Logical. Whether to return the covariance matrix of `fsT`
#'              and `fsL`, returned as attribute `vfsLT`; used for
#'              second-order SE correction of 2S-PA results. Currently
#'              ignored for `merMod` objects.
#' @param reliability Logical. Whether to return the reliability of factor
#'                    scores. Available only for single-factor lavaan models;
#'                    for multi-factor models a warning is issued and no
#'                    `reliability` attribute is returned.
#' @param prior_mean An optional numeric vector of length `q` (the number of
#'        latent variables) giving fixed external prior means for the latent
#'        variables. `NULL` (default) uses the lavaan-estimated (group-specific)
#'        latent means. Non-NULL values are treated as fixed external priors
#'        shared across all lavaan groups. For `mirt` `SingleGroupClass`
#'        objects it instead sets the factor prior mean used for the EAP
#'        scores; the factor-score intercepts (`fsb`) then vary per observation
 #'        as `Vpost_i %*% solve(psi) %*% prior_mean`, i.e. the latent mean
 #'        scaled by the per-observation shrinkage factor (zero when
 #'        `prior_mean = NULL`), where `psi` is the mirt model's estimated
 #'        factor covariance. For `mirt` `MultipleGroupClass` objects a
 #'        non-NULL `prior_mean` (length `q`) is applied as the factor prior
 #'        mean to every group (mirt's per-group EAP is otherwise centred on a
 #'        zero-mean standard-normal prior); each observation's regression form
 #'        uses the factor covariance of its own group.
#'        Only supported for lavaan objects with regression (EB) scoring (and
#'        for mirt); `reliability = TRUE` is not supported together with
#'        user-supplied `prior_mean`/`prior_cov`, and `prior_cov` is not
#'        supported for mirt. Conceptually similar to the `mean` argument of
#'        `mirt::fscores()`.
#' @param prior_cov An optional numeric `q x q` covariance matrix (a scalar or
#'        1 x 1 matrix is accepted when `q = 1`) giving fixed external prior
#'        covariance for the latent variables. `NULL` (default) uses the
#'        lavaan-estimated (group-specific) latent covariance. Non-NULL values
#'        must be finite, symmetric and positive definite, and are treated as
#'        fixed external priors shared across all lavaan groups. Only supported
#'        for lavaan objects with regression (EB) scoring;
#'        `reliability = TRUE` is not supported together with user-supplied
#'        `prior_mean`/`prior_cov`. With `corrected_fsT = TRUE` or
#'        `vfsLT = TRUE` the supplied covariance is treated as fixed, i.e. no
#'        sampling uncertainty from the prior itself is propagated.
#'        Conceptually similar to the `cov` argument of `mirt::fscores()`.
#' @param sum_items A named list mapping each factor name to the item names
#'        that make up its sum score, e.g.
#'        `list(ind60 = c("x1", "x2", "x3"), dem60 = c("y1", "y2", "y3",
#'        "y4"))`. `NULL` (default) auto-derives the assignment from the
#'        estimated loadings, which requires each indicator to load on
#'        exactly one factor and every factor to have at least one item. A
#'        supplied list must cover all model factors, and each item may
#'        belong to only one sum. Only used for `lavaan` and data frame
#'        objects with `method = "mean"`.
#' @param format Output format when `object` is a lavaan or merMod object.
#'        `"unified"` (default) returns a single data frame; for multiple
#'        groups it carries a `group` column and attributes `fsT`, `fsL`,
#'        `fsb`, and `scoring_matrix` are named lists keyed by group
#'        label. `"list"` returns the legacy shape: a named list of data
 #'        frames (one per group) with per-group matrix attributes. Use
 #'        [fs_to_group_list()] to convert between the two. For `mirt`
 #'        `SingleGroupClass` and `MultipleGroupClass` objects `format` is
 #'        accepted but the output is always a single per-observation data
 #'        frame; the multi-group result additionally carries a trailing
 #'        `group` column (the model's group levels, `NA` for
 #'        completely-missing rows) and a per-group (`list`) `psi` attribute.
#' @param ... additional arguments passed to \code{\link[lavaan]{cfa}}
#'            (when `object` is a data frame). See \code{\link[lavaan]{lavOptions}}
#'            for a complete list.
#' @return A data frame containing the factor scores (with prefix `"fs_"`),
#'         the standard errors (with suffix `"_se"`), the implied loadings
#'         of indicator `_by_` factor scores, and the error variance-covariance
#'         of the factor scores (with prefix `"ev_"` or `"ecov_"`).
#'         For multi-group lavaan models in `"unified"` format, a `group` column
#'         is included. The following attributes are attached:
#'         * `fsT`: error covariance of factor scores (matrix or named list by group)
#'         * `fsL`: loading matrix of factor scores (matrix or named list by group)
#'         * `fsb`: intercepts of factor scores (vector or named list by
#'           group); with `method = "mean"` the intercept is the mean of the
#'           factor's item intercepts, the measurement intercept of the score
#'           regressed on the uncentered latent (same `E[fs] - fsL %*% alpha`
#'           convention as the other methods; equals the score's column mean
#'           for models without a mean structure)
#'         * `scoring_matrix`: weights for computing factor scores from the
#'           observed data, as a named list. For lavaan models: one
#'           score x item matrix per group; with `method = "mean"` the
#'           weights are the item-mean weights, so `S %*% y` reproduces the
#'           raw scores exactly (no centering offset). For `merMod` models:
#'           one `num_re` x `n_j` matrix per cluster, where
#'           `S_j %*% (y_j - X_j %*% beta)` with `y_j`/`X_j` the cluster's
#'           rows of the model response and the fixed-effects design
#'           reproduces the cluster's EB scores for method `"EB"` and the
#'           per-cluster OLS (ML) scores for method `"ML"`.
#'         * `psi`: effective (prior-adjusted) covariance matrix of the
#'           latent variables (`q x q`), group-level (not per-pattern), and
#'           a point estimate only (no sampling SEs of the latents are
#'           attached). Mirrors the `fsT` shape: a named list keyed by group
#'           label for `"unified"` output; a direct attribute on each group
#'           data frame (plus a list-valued attribute on the outer list) for
#'           `"list"` output; for `merMod` objects a single `q x q` matrix.
#'           With `prior_cov` supplied it equals the prior (shared across
#'           groups), otherwise the per-group lavaan estimate. For `merMod`
#'           objects the matrix is the first random-effects term's
#'           `VarCorr`, with dimnames renamed to match the `fsL` column
#'           names (`u0`/`u1`/..., or the legacy `u0_eb`/`u1_eb` names).
#'         * `alpha`: effective (prior-adjusted) means of the latent
#'           variables (a named vector of length `q`), with the same group
#'           nesting and point-estimate semantics as `psi`. With
#'           `prior_mean` supplied it equals the prior, otherwise the
#'           per-group lavaan estimate; a named zero vector (`0` per latent)
#'           when the model has no (estimated) mean structure. For `merMod`
#'           objects a named zero vector (random effects are mean zero).
#'         * `fs_pattern`: for lavaan models, a named list by group of
#'           `list(label, pat)` entries. `label` is a character vector with
#'           one entry per case in the group giving that case's
#'           observed-indicator pattern name (`NA` for cases whose indicators
#'           are all missing); `pat` is a logical matrix with rows =
#'           indicators and one column per pattern, the columns being named
#'           by pattern name.
#'
#'         For a lavaan group without missing data, its `fsT`/`fsL`/`fsb`/
#'         `scoring_matrix` elements are the plain matrix/vector for the whole
#'         group. When a group's cases split into multiple observed-indicator
#'         patterns (missing data), each such element is instead a named list
#'         with one entry per pattern; the pattern name is the observed
#'         indicator names joined with `"+"` in indicator order (e.g.
#'         `"x1+x3"`).
#'
#'         Note: for a single-group lavaan fit in `"unified"` format, the
#'         per-group attribute wrappers (`fsT`, `fsL`, `fsb`,
#'         `scoring_matrix`, `psi`, and `alpha`) are each a one-element list
#'         named with the empty string `""`; `x[[""]]` does not match in R
#'         list subsetting, so read these attributes positionally (e.g.
#'         `attr(fs, "fsT")[[1]]`, `attr(fs, "psi")[[1]]`) rather than by
#'         name.
#' @importFrom lavaan cfa sem
#' @importFrom lavaan lavInspect lavTech coef
#' @importFrom lavaan vcov
#' @importFrom stats setNames
#'
#' @export
#'
#' @examples
#' library(lavaan)
#' get_fs(PoliticalDemocracy[c("x1", "x2", "x3")])
#'
#' # Multiple factors
#' get_fs(PoliticalDemocracy[c("x1", "x2", "x3", "y1", "y2", "y3", "y4")],
#'        model = " ind60 =~ x1 + x2 + x3
#'                  dem60 =~ y1 + y2 + y3 + y4 ")
#'
#' # Multiple-group
#' hs_model <- ' visual  =~ x1 + x2 + x3 '
#' fit <- cfa(hs_model,
#'            data = HolzingerSwineford1939,
#'            group = "school")
#' get_fs(HolzingerSwineford1939, hs_model, group = "school")
#' # Or without the model
#' get_fs(HolzingerSwineford1939[c("school", "x4", "x5", "x6")],
#'        group = "school")
#'
#' # Fixed external latent prior (shared across groups) for regression scores;
#' # conceptually similar to mirt::fscores(mean, cov)
#' fit <- cfa("visual =~ x1 + x2 + x3",
#'            data = HolzingerSwineford1939,
#'            group = "school", group.equal = c("loadings", "intercepts"))
#' get_fs(fit, prior_mean = c(visual = -0.12), prior_cov = 0.33)

get_fs <- function(object, ...) {
  UseMethod("get_fs")
}

#' @inherit get_fs
#' @param lavobj A lavaan model object when using [get_fs_lavaan()].
#'
#' @details
#' `get_fs_lavaan()` is superseded by [get_fs()]. It is retained for backward
#' compatibility and delegates to `get_fs(object, format = "list")` internally.
#' New code should call [get_fs()] directly.
#'
#' @export
get_fs_lavaan <- function(
  lavobj,
  method = c("regression", "Bartlett", "ML", "EB", "mean"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  prior_mean = NULL,
  prior_cov = NULL,
  sum_items = NULL,
  ...
) {
  get_fs(
    lavobj,
    method = method,
    corrected_fsT = corrected_fsT,
    vfsLT = vfsLT,
    reliability = reliability,
    format = "list",
    prior_mean = prior_mean,
    prior_cov = prior_cov,
    sum_items = sum_items,
    ...
  )
}

#' Convert Unified Factor Scores to Group List (or Vice Versa)
#'
#' @description
#' `fs_to_group_list()` converts a unified factor-score data frame
#' (single data frame with a `group` column and list-valued attributes)
#' into the legacy list-of-data-frames shape with per-group attributes.
#' It acts as its own inverse: if given a group list, it converts back to
#' the unified shape.
#'
#' @param fs Object returned by [get_fs()]. Either a unified data frame
#'           (with a `group` column and list-valued attributes) or a named
#'           list of data frames (one per group) with per-group attributes.
#'
#' @return If `fs` is a unified data frame, returns a named list of data
#'         frames, one per group, with `fsT`, `fsL`, `fsb`, `scoring_matrix`,
#'         and `fs_pattern` attached as per-group attributes on each element
#'         and as list-valued attributes on the outer list. If `fs` is a
#'         group list, returns a single data frame with a `group` column and
#'         list-valued attributes. A single-group input returns a single
#'         data frame (not a list) in either direction.
#'
#' @export
#'
#' @examples
#' library(lavaan)
#' hs_model <- "visual =~ x1 + x2 + x3"
#' fit <- cfa(hs_model, data = HolzingerSwineford1939, group = "school")
#' fs_unified <- get_fs(fit)                     # unified df
#' fs_list <- fs_to_group_list(fs_unified)       # list-of-df shape
#' fs_back <- fs_to_group_list(fs_list)          # back to unified
#' all.equal(fs_unified, fs_back, check.attributes = FALSE)
fs_to_group_list <- function(fs) {
  # psi/alpha are the group-level latent moments (see get_fs() @return);
  # they are carried through the unified <-> list conversion like fsT.
  attr_keys <- c("fsT", "fsL", "fsb", "scoring_matrix", "fs_pattern",
                 "psi", "alpha")

  if (is.data.frame(fs)) {
    grp_col <- attr(fs, "group_col")
    if (is.null(grp_col)) {
      grp_col <- "group"
    }

    if (!grp_col %in% names(fs)) {
      # Single-group unified result without group column — unwrap attributes
      out <- fs
      for (ak in attr_keys) {
        outer <- attr(fs, ak)
        if (is.list(outer) && length(outer) == 1L) {
          attr(out, ak) <- outer[[1L]]
        } else {
          attr(out, ak) <- outer
        }
      }
      return(out)
    }

    group_labels <- unique(fs[[grp_col]])

    if (length(group_labels) == 1) {
      out <- fs[, !names(fs) %in% grp_col, drop = FALSE]
      for (ak in attr_keys) {
        outer <- attr(fs, ak)
        if (is.list(outer) && length(outer) == 1L) {
          attr(out, ak) <- outer[[1L]]
        } else {
          attr(out, ak) <- outer
        }
      }
      return(out)
    }

    grp_dfs <- split(fs, fs[[grp_col]])
    grp_dfs <- grp_dfs[group_labels]
    for (g in group_labels) {
      grp_dfs[[g]] <- grp_dfs[[g]][,
        !names(grp_dfs[[g]]) %in% grp_col,
        drop = FALSE
      ]
      for (ak in attr_keys) {
        outer <- attr(fs, ak)
        if (is.list(outer) && !is.null(names(outer))) {
          attr(grp_dfs[[g]], ak) <- outer[[g]]
        } else {
          attr(grp_dfs[[g]], ak) <- outer
        }
      }
    }

    # Re-attach outer list-level attributes
    for (ak in attr_keys) {
      outer <- attr(fs, ak)
      if (is.list(outer)) {
        attr(grp_dfs, ak) <- outer
      }
    }
    grp_dfs
  } else if (is.list(fs)) {
    group_labels <- names(fs)
    if (is.null(group_labels) || !all(nzchar(group_labels))) {
      stop("Group list must be named (by group label).", call. = FALSE)
    }

    if (length(group_labels) == 1) {
      out <- fs[[1L]]
      if ("group" %in% names(out)) {
        out <- out[, !names(out) %in% "group", drop = FALSE]
      }
      return(out)
    }

    group_dfs <- lapply(group_labels, function(g) {
      df <- fs[[g]]
      if ("group" %in% names(df)) {
        df <- df[, !names(df) %in% "group", drop = FALSE]
      }
      df[["group"]] <- g
      df
    })
    unified <- do.call(rbind, group_dfs)
    attr(unified, "group_col") <- "group"

    for (ak in attr_keys) {
      per_grp <- lapply(group_labels, function(g) {
        a <- attr(fs[[g]], ak)
        if (!is.null(a)) {
          return(a)
        }
        outer <- attr(fs, ak)
        if (is.list(outer) && !is.null(names(outer))) {
          return(outer[[g]])
        }
        NULL
      })
      if (!all(vapply(per_grp, is.null, logical(1)))) {
        attr(unified, ak) <- setNames(per_grp, group_labels)
      }
    }
    unified
  } else {
    stop("'fs' must be data frame or list.", call. = FALSE)
  }
}

augment_fs <- function(fs, fs_ev) {
  fsL <- attr(fs, "fsL")
  # Column values (se / loadings / lower-tri error terms) come from the
  # shared value-only engine fs_row_cols() (R/fs_indiv.R); this function
  # supplies the r2spa column naming.
  vals <- fs_row_cols(fs, fsL, fs_ev)
  num_lvs <- ncol(fs_ev)
  colnames(fs) <- paste0("fs_", colnames(fs))
  fs_se_names <- paste0(rownames(fs_ev), "_se")
  fs_names <- paste0("fs_", colnames(fsL))
  # Bare name vector (i-outer over the latents), matching the value order
  # of the ld block in fs_row_cols() -- c(as.matrix(fsL)), column-major.
  fs_ld_names <- unlist(lapply(seq_len(ncol(fsL)), function(i) {
    paste(colnames(fsL)[i], fs_names, sep = "_by_")
  }), use.names = FALSE)
  fs_ev_names <- character(num_lvs * (num_lvs + 1) / 2)
  count <- 1
  for (i in seq_len(num_lvs)) {
    for (j in seq_len(i)) {
      if (i == j) {
        fs_ev_names[count] <- paste0("ev_", rownames(fs_ev)[i])
      } else {
        fs_ev_names[count] <- paste0(
          "ecov_",
          rownames(fs_ev)[i],
          "_",
          colnames(fs_ev)[j]
        )
      }
      count <- count + 1
    }
  }
  colnames(vals) <- c(fs_se_names, fs_ld_names, fs_ev_names)
  fs_dat <- cbind(
    as.data.frame(fs),
    vals
  )
  attr(fs_dat, "fsT") <- fs_ev
  attr(fs_dat, "fsL") <- fsL
  attr(fs_dat, "fsb") <- attr(fs, "fsb")
  attr(fs_dat, "scoring_matrix") <- attr(fs, "scoring_matrix")
  fs_dat
}

assemble_fs_blocks <- function(
  blocks_by_group,
  format = c("unified", "list"),
  group_col = NULL
) {
  format <- match.arg(format)
  group_labels <- names(blocks_by_group)
  if (is.null(group_labels) || !all(nzchar(group_labels))) {
    group_labels <- rep("", length(blocks_by_group))
  }
  attr_keys <- c("fsT", "fsL", "fsb", "scoring_matrix", "fs_pattern")

  group_dfs <- vector("list", length(group_labels))
  names(group_dfs) <- group_labels

  for (g in seq_along(group_labels)) {
    blocks <- blocks_by_group[[g]]
    n_cases <- max(unlist(lapply(blocks, function(b) max(b$case_idx))))

    aug_list <- lapply(blocks, function(b) {
      augment_fs(b$fs, b$fsT)
    })

    template_cols <- names(aug_list[[1]])
    grp_df <- as.data.frame(
      matrix(NA, nrow = n_cases, ncol = length(template_cols))
    )
    colnames(grp_df) <- template_cols

    for (bl in seq_along(blocks)) {
      grp_df[blocks[[bl]]$case_idx, ] <- aug_list[[bl]]
    }

    block_attrs <- lapply(blocks, function(b) {
      list(
        fsT = b$fsT,
        fsL = if (!is.null(b$fsL)) b$fsL else attr(b$fs, "fsL"),
        fsb = if (!is.null(b$fsb)) b$fsb else attr(b$fs, "fsb"),
        scoring_matrix = if (!is.null(b$scoring_matrix)) {
          b$scoring_matrix
        } else {
          attr(b$fs, "scoring_matrix")
        }
      )
    })

    # Pattern bookkeeping (mirrors lavaan's @Data@Mp in public form): each
    # block carries the observed-indicator pattern of its cases (pat_label =
    # observed indicator names joined with "+", pat = one named logical
    # column per indicator). Blocks without pat_label/pat (hand-built
    # fixtures) get a positional "pattern_<i>" label and a NULL pat matrix.
    pat_labels <- vapply(seq_along(blocks), function(m) {
      if (!is.null(blocks[[m]]$pat_label)) {
        blocks[[m]]$pat_label
      } else {
        paste0("pattern_", m)
      }
    }, character(1))

    label_vec <- rep(NA_character_, n_cases)
    for (m in seq_along(blocks)) {
      label_vec[blocks[[m]]$case_idx] <- pat_labels[m]
    }

    if (length(blocks) > 1) {
      # One attribute value per observed-indicator pattern, keyed by the
      # pattern label.
      for (ak in attr_keys) {
        attr(grp_df, ak) <- setNames(lapply(block_attrs, `[[`, ak), pat_labels)
      }
      if (all(vapply(blocks, function(b) !is.null(b$pat), logical(1)))) {
        attr(grp_df, "fs_pattern") <- list(
          label = label_vec,
          pat = do.call(cbind, setNames(lapply(blocks, `[[`, "pat"), pat_labels))
        )
      } else {
        attr(grp_df, "fs_pattern") <- list(label = label_vec, pat = NULL)
      }
    } else {
      for (ak in attr_keys) {
        attr(grp_df, ak) <- block_attrs[[1]][[ak]]
      }
      pat1 <- blocks[[1]]$pat
      if (!is.null(pat1)) {
        attr(grp_df, "fs_pattern") <- list(
          label = label_vec,
          pat = matrix(
            pat1,
            ncol = 1L,
            dimnames = list(names(pat1), pat_labels)
          )
        )
      } else {
        attr(grp_df, "fs_pattern") <- list(label = label_vec, pat = NULL)
      }
    }

    group_dfs[[g]] <- grp_df
  }

  if (format == "list") {
    if (length(group_labels) == 1 && group_labels == "") {
      return(group_dfs[[1]])
    }

    if (!is.null(group_col)) {
      for (g in seq_along(group_labels)) {
        if (!group_col %in% names(group_dfs[[g]])) {
          group_dfs[[g]][[group_col]] <- group_labels[g]
        }
      }
    }

    outer_attr_names <- setdiff(
      names(attributes(group_dfs[[1]])),
      c("names", "class", "row.names", "col.names")
    )
    for (ak in outer_attr_names) {
      attr_lst <- setNames(
        lapply(group_dfs, function(df) attr(df, ak)),
        group_labels
      )
      attr(group_dfs, ak) <- attr_lst
    }
    return(group_dfs)
  }

  has_group <- any(nzchar(group_labels))
  col_name <- if (!is.null(group_col)) group_col else "group"
  if (has_group) {
    for (g in seq_along(group_labels)) {
      group_dfs[[g]][[col_name]] <- group_labels[g]
    }
  }
  unified_df <- do.call(rbind, group_dfs)
  rownames(unified_df) <- NULL
  if (has_group) {
    attr(unified_df, "group_col") <- col_name
  }

  for (ak in attr_keys) {
    attr(unified_df, ak) <- setNames(
      lapply(group_dfs, function(df) attr(df, ak)),
      group_labels
    )
  }

  unified_df
}

#' Get Factor Scores and the Corresponding Scoring Matrices for
#' Mixed-Effect Models
#'
#' @description
#' `get_fs_lmer()` is superseded by [get_fs()]. It is retained for backward
#' compatibility and delegates to `get_fs(object)` internally. New code should
#' call [get_fs()] directly.
#'
#' @param object A fitted model object of class [lme4::lmerMod-class].
#' @param method `"EB"` (empirical Bayes, default) or `"ML"` (prior-free
#'        per-cluster OLS), forwarded to [get_fs.merMod()].
#' @param corrected_fsT Currently not used.
#' @param vfsLT Currently not used.
#' @param fsm Currently not used.
#' @param legacy_names Logical. Passed to [get_fs.merMod()]. Defaults to
#'        `TRUE` so `get_fs_lmer()` keeps returning the pre-refactor
#'        `u0_eb`-style *column names* (in the legacy column order).
#'        Note the legacy output is name-compatible, not byte-identical,
#'        with the pre-refactor result: it additionally carries
#'        score-error columns (`u0_eb_se`, ...), per-cluster `fsL`/`fsT`
#'        array attributes, a per-cluster `scoring_matrix` list attribute
#'        (see [get_fs()]), and has NULL row names (the pre-refactor
#'        output had no `_se` columns, no attributes, and used the ranef
#'        subject IDs as row names).
#' @param ... Additional arguments, passed on to [get_fs()].
#'
#' @export
get_fs_lmer <- function(
  object,
  method = c("EB"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  fsm = FALSE,
  legacy_names = TRUE,
  ...
) {
  get_fs(
    object,
    method = method,
    format = "list",
    legacy_names = legacy_names,
    ...
  )
}
