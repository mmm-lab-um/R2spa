#' Get Factor Scores and the Corresponding Standard Error of Measurement
#'
#' @description
#' `get_fs()` is an S3 generic that extracts factor scores from fitted models.
#' Methods are available for `data.frame` (fits a CFA internally), `lavaan`
#' objects, and `lmerMod` objects. Support for `mirt` models is planned.
#'
#' @details
#' When `data` is a data frame and `model` is supplied as a lavaan syntax string,
#' the function internally calls [lavaan::cfa()] and then dispatches to the
#' `lavaan` method. When `data` is a fitted model object, the appropriate S3
#' method is called directly.
#'
#' `get_fs()` replaced `get_fs_lavaan()` and `get_fs_lmer()`, which are now
#' thin wrappers retained for backward compatibility.
#'
#' @param data A data frame, a fitted [lavaan] model object, or a fitted
#'        [lme4::lmer] model object (`merMod`).
#' @param model An optional string specifying the measurement model
#'              in \code{lavaan} syntax. Only used when `data` is a data frame.
#'              See \code{\link[lavaan]{model.syntax}} for more information.
#' @param group Character. Name of the grouping variable for multiple group
#'              analysis, which is passed to \code{\link[lavaan]{cfa}}.
#'              Only used when `data` is a data frame.
#' @param method Character. Method for computing factor scores (options are
#'               "regression" or "Bartlett"; "ML" is an alias for "Bartlett"
#'               and "EB" is an alias for "regression"). Currently, the
#'               default is "regression" to be consistent with
#'               \code{\link[lavaan]{lavPredict}}, but the Bartlett scores
#'               have more desirable properties and may be preferred for
#'               2S-PA.
#' @param corrected_fsT Logical. Whether to correct for the sampling
#'                      error in the factor score weights when computing
#'                      the error variance estimates of factor scores.
#' @param vfsLT Logical. Whether to return the covariance matrix of `fsT`
#'              and `fsL`, which can be used as input for [vcov_corrected()]
#'              to obtain corrected covariances and standard errors for
#'              [tspa()] results. This is currently ignored.
#' @param reliability Logical. Whether to return the reliability of factor
#'                    scores. Available only for single-factor lavaan models.
#' @param format Output format when `data` is a lavaan or merMod object.
#'        `"unified"` (default) returns a single data frame with a `group` column;
#'        for multiple groups, attributes `fsT`, `fsL`, `fsb`, and `scoring_matrix`
#'        are named lists keyed by group label. `"list"` returns the legacy
#'        shape: a named list of data frames (one per group) with per-group
#'        matrix attributes. Use [fs_to_group_list()] to convert between the two.
#' @param ... additional arguments passed to \code{\link[lavaan]{cfa}}
#'            (when `data` is a data frame). See \code{\link[lavaan]{lavOptions}}
#'            for a complete list.
#' @return A data frame containing the factor scores (with prefix `"fs_"`),
#'         the standard errors (with suffix `"_se"`), the implied loadings
#'         of indicator `_by_` factor scores, and the error variance-covariance
#'         of the factor scores (with prefix `"ev_"` or `"ecov_"`).
#'         For multi-group lavaan models in `"unified"` format, a `group` column
#'         is included. The following attributes are attached:
#'         * `fsT`: error covariance of factor scores (matrix or named list by group)
#'         * `fsL`: loading matrix of factor scores (matrix or named list by group)
#'         * `fsb`: intercepts of factor scores (vector or named list by group)
#'         * `scoring_matrix`: weights for computing factor scores from items
#' @importFrom lavaan cfa sem
#' @importFrom lavaan lavInspect lavTech coef
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

get_fs <- function(data, ...) {
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
  method = c("regression", "Bartlett", "ML", "EB"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  ...
) {
  get_fs(
    lavobj,
    method = method,
    corrected_fsT = corrected_fsT,
    vfsLT = vfsLT,
    reliability = reliability,
    format = "list",
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
#'         frames, one per group, with `fsT`, `fsL`, `fsb`, and
#'         `scoring_matrix` attached as per-group attributes on each element
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
  attr_keys <- c("fsT", "fsL", "fsb", "scoring_matrix")

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
  fs_se <- t(as.matrix(sqrt(diag(fs_ev))))
  colnames(fs) <- paste0("fs_", colnames(fs))
  colnames(fs_se) <- paste0(colnames(fs_se), "_se")
  num_lvs <- ncol(fs_ev)
  fs_evs <- rep(NA, num_lvs * (num_lvs + 1) / 2)
  count <- 1
  for (i in seq_len(num_lvs)) {
    for (j in seq_len(i)) {
      fs_evs[count] <- fs_ev[i, j]
      if (i == j) {
        names(fs_evs)[count] <- paste0("ev_", rownames(fs_ev)[i])
      } else {
        names(fs_evs)[count] <- paste0(
          "ecov_",
          rownames(fs_ev)[i],
          "_",
          colnames(fs_ev)[j]
        )
      }
      count <- count + 1
    }
  }
  fsL <- attr(fs, "fsL")
  fs_names <- paste0("fs_", colnames(fsL))
  fs_lds <- lapply(seq_len(ncol(fsL)), function(i) {
    setNames(fsL[, i], paste(colnames(fsL)[i], fs_names, sep = "_by_"))
  })
  fs_lds <- unlist(fs_lds)
  fs_dat <- cbind(
    as.data.frame(fs),
    fs_se,
    t(as.matrix(fs_lds)),
    t(as.matrix(fs_evs))
  )
  attr(fs_dat, "fsT") <- fs_ev
  attr(fs_dat, "fsL") <- fsL
  attr(fs_dat, "fsb") <- attr(fs, "fsb")
  attr(fs_dat, "scoring_matrix") <- attr(fs, "scoring_matrix")
  fs_dat
}

check_blocks_identical <- function(a, b, keys) {
  all(vapply(
    keys,
    function(k) {
      identical(a[[k]], b[[k]])
    },
    logical(1)
  ))
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
  attr_keys <- c("fsT", "fsL", "fsb", "scoring_matrix")

  group_dfs <- vector("list", length(group_labels))
  names(group_dfs) <- group_labels

  for (g in seq_along(group_labels)) {
    grp <- group_labels[g]
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

    if (length(blocks) > 1) {
      all_same <- vapply(
        seq_len(length(blocks))[-1],
        function(i) {
          check_blocks_identical(block_attrs[[i]], block_attrs[[1]], attr_keys)
        },
        logical(1)
      )
      if (!all(all_same)) {
        message(
          "Group '",
          grp,
          "': blocks have differing fsT/fsL/fsb attributes ",
          "(e.g. due to missing-data patterns). Using first block as ",
          "representative for group-level attributes."
        )
      }
    }
    repr <- block_attrs[[1]]

    for (ak in attr_keys) {
      attr(grp_df, ak) <- repr[[ak]]
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
#' @param method Currently only `"EB"` for empirical Bayes.
#' @param corrected_fsT Currently not used.
#' @param vfsLT Currently not used.
#' @param fsm Currently not used.
#' @param legacy_names Logical. Passed to [get_fs.merMod()]. Defaults to
#'        `TRUE` so `get_fs_lmer()` keeps returning the pre-refactor
#'        `u0_eb`-style *column names* (in the legacy column order).
#'        Note the legacy output is name-compatible, not byte-identical,
#'        with the pre-refactor result: it additionally carries
#'        score-error columns (`u0_eb_se`, ...), per-cluster `fsL`/`fsT`
#'        array attributes, and has NULL row names (the pre-refactor
#'        output had no `_se` columns, no attributes, and used the ranef
#'        subject IDs as row names).
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
    format = "list",
    legacy_names = legacy_names,
    ...
  )
}
