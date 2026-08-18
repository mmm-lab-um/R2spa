# S3 methods for get_fs()
#
# Future methods (e.g. get_fs.mirt()) should be added to this file.

normalize_fs_method <- function(method) {
  method <- match.arg(method, c("regression", "Bartlett", "ML", "EB"))
  switch(method, ML = "Bartlett", EB = "regression", method)
}

validate_fs_priors <- function(prior_mean, prior_cov, lv_names) {
  q <- length(lv_names)

  if (!is.null(prior_mean)) {
    prior_nms <- names(prior_mean)
    prior_mean <- as.numeric(prior_mean)
    if (!is.null(prior_nms)) {
      if (!setequal(prior_nms, lv_names)) {
        stop(
          "'prior_mean' names must match the latent variable names (",
          paste(lv_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      prior_mean <- prior_mean[match(lv_names, prior_nms)]
    } else {
      if (length(prior_mean) != q) {
        stop(
          "'prior_mean' must have length ", q,
          " (the number of latent variables).",
          call. = FALSE
        )
      }
    }
    if (any(!is.finite(prior_mean))) {
      stop("'prior_mean' must contain only finite values.", call. = FALSE)
    }
    names(prior_mean) <- lv_names
  }

  if (!is.null(prior_cov)) {
    if (!is.matrix(prior_cov)) {
      prior_cov <- as.matrix(prior_cov)
    }
    if (any(!is.finite(prior_cov))) {
      stop("'prior_cov' must contain only finite values.", call. = FALSE)
    }
    if (nrow(prior_cov) != ncol(prior_cov)) {
      stop("'prior_cov' must be a square matrix.", call. = FALSE)
    }
    if (nrow(prior_cov) != q) {
      stop(
        "'prior_cov' must be a ", q, " x ", q, " matrix (one row and ",
        "column per latent variable).",
        call. = FALSE
      )
    }
    pcv <- prior_cov
    dimnames(pcv) <- NULL
    if (!isTRUE(all.equal(pcv, t(pcv)))) {
      stop("'prior_cov' must be symmetric.", call. = FALSE)
    }
    eig_vals <- eigen(pcv, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig_vals <= 0)) {
      stop("'prior_cov' must be positive definite.", call. = FALSE)
    }
    rn <- rownames(prior_cov)
    cn <- colnames(prior_cov)
    if (q == 1) {
      if (is.null(rn) && !is.null(cn)) {
        rn <- cn
      }
      if (is.null(cn) && !is.null(rn)) {
        cn <- rn
      }
    }
    if (!is.null(rn) || !is.null(cn)) {
      if (!setequal(rn, lv_names) || !setequal(cn, lv_names)) {
        stop(
          "'prior_cov' names must match the latent variable names (",
          paste(lv_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      prior_cov <- prior_cov[match(lv_names, rn),
                             match(lv_names, cn),
                             drop = FALSE]
    }
    rownames(prior_cov) <- colnames(prior_cov) <- lv_names
  }

  list(mean = prior_mean, cov = prior_cov)
}

#' @rdname get_fs
#' @export
get_fs.data.frame <- function(
  object,
  model = NULL,
  group = NULL,
  method = c("regression", "Bartlett", "ML", "EB"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  format = c("unified", "list"),
  prior_mean = NULL,
  prior_cov = NULL,
  ...
) {
  if (is.null(model)) {
    ind_names <- colnames(object)
    if (!is.null(group)) {
      ind_names <- setdiff(ind_names, group)
    }
    model <- paste("f1 =~", paste(ind_names, collapse = " + "))
  }
  fit <- cfa(model, data = object, group = group, ...)
  get_fs(
    fit,
    method = method,
    corrected_fsT = corrected_fsT,
    vfsLT = vfsLT,
    reliability = reliability,
    format = format,
    prior_mean = prior_mean,
    prior_cov = prior_cov
  )
}

get_fs_blocks.lavaan <- function(
  object,
  method,
  add_to_evfs,
  prior_mean = NULL,
  prior_cov = NULL,
  ...
) {
  method <- match.arg(method, c("regression", "Bartlett"))
  est <- lavInspect(object, what = "est")
  y <- lavInspect(object, what = "data")
  miss_pat <- object@Data@Mp

  prepare_fs <- function(y, est, add, mp, method) {
    psi_use <- if (is.null(prior_cov)) est$psi else prior_cov
    alpha_use <- if (is.null(prior_mean)) est$alpha else prior_mean
    if (is.null(mp)) {
      fscore <-
        compute_fscore(
          y,
          lambda = est$lambda,
          theta = est$theta,
          psi = psi_use,
          nu = est$nu,
          alpha = alpha_use,
          method = method,
          fs_matrices = TRUE
        )
      list(
        case_idx = seq_len(nrow(y)),
        fs = fscore,
        fsT = attr(fscore, "fsT") + add,
        fsL = attr(fscore, "fsL"),
        fsb = attr(fscore, "fsb"),
        scoring_matrix = attr(fscore, "scoring_matrix"),
        pat_label = paste0(colnames(y), collapse = "+"),
        pat = setNames(rep_len(TRUE, ncol(y)), colnames(y))
      )
    } else {
      npat <- mp$npatterns
      pats <- mp$pat
      mis_idx <- mp$case.idx
      blocks <- vector("list", npat)
      for (m in seq_len(npat)) {
        idx_m <- mis_idx[[m]]
        pat_m <- pats[m, ]
        fs_m <-
          compute_fscore(
            y[idx_m, pat_m, drop = FALSE],
            lambda = est$lambda[pat_m, , drop = FALSE],
            theta = est$theta[pat_m, pat_m, drop = FALSE],
            psi = psi_use,
            nu = est$nu[pat_m, , drop = FALSE],
            alpha = alpha_use,
            method = method,
            fs_matrices = TRUE
          )
        blocks[[m]] <- list(
          case_idx = idx_m,
          fs = fs_m,
          fsT = attr(fs_m, "fsT") + add,
          fsL = attr(fs_m, "fsL"),
          fsb = attr(fs_m, "fsb"),
          scoring_matrix = attr(fs_m, "scoring_matrix"),
          pat_label = paste0(colnames(y)[pat_m], collapse = "+"),
          pat = setNames(pat_m, colnames(y))
        )
      }
      blocks
    }
  }

  # Direct slot access avoids lavInspect()'s per-call version-compatibility
  # check (lav_object_check_version(), which re-reads lavaan's DESCRIPTION
  # file via read.dcf() every time) -- a major source of overhead when
  # get_fs() is called repeatedly. @Data@ngroups is already relied upon
  # elsewhere in this file (e.g. @Data@Mp, @Data@group.label).
  ngroups <- object@Data@ngroups

  if (ngroups == 1) {
    blocks <- prepare_fs(y, est, add_to_evfs[[1]], miss_pat[[1]], method)
    if (is.null(miss_pat[[1]])) {
      blocks <- list(blocks)
    }
    return(setNames(list(blocks), ""))
  }

  group_labels <- object@Data@group.label
  blocks_by_group <- setNames(vector("list", ngroups), group_labels)

  for (g in seq_len(ngroups)) {
    grp_est <- est[[g]]
    grp_y <- y[[g]]
    grp_add <- add_to_evfs[[g]]
    grp_mp <- miss_pat[[g]]
    blocks_by_group[[g]] <- prepare_fs(grp_y, grp_est, grp_add, grp_mp, method)
    if (is.null(grp_mp)) {
      blocks_by_group[[g]] <- list(blocks_by_group[[g]])
    }
  }

  blocks_by_group
}

#' @rdname get_fs
#' @param format Output format: `"unified"` returns a single data frame with
#'        a `group` column; `"list"` returns a list of data frames per group.
#' @export
get_fs.default <- function(object, ...) {
  if (is.matrix(object)) {
    object <- as.data.frame(object)
    return(get_fs(object, ...))
  }
  stop(
    "get_fs() is not implemented for objects of class '",
    paste(class(object), collapse = "', '"),
    "'. Currently supported: 'data.frame', 'lavaan', ",
    "and 'lmerMod'. Support for 'mirt' models is planned.",
    call. = FALSE
  )
}

#' @rdname get_fs
#' @param format Output format: `"unified"` returns a single data frame with
#'        a `group` column; `"list"` returns a list of data frames per group.
#' @export
get_fs.lavaan <- function(
  object,
  method = c("regression", "Bartlett", "ML", "EB"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  format = c("unified", "list"),
  prior_mean = NULL,
  prior_cov = NULL,
  ...
) {
  method <- normalize_fs_method(method)
  if (!inherits(object, "lavaan")) {
    stop("`object` must be a `lavaan` model object.")
  }
  format <- match.arg(format)

  # The sampling-error SE paths (corrected_fsT/reliability/vfsLT) require
  # the full observed-indicator set, so they cannot be applied to lavaan's
  # per-missing-data-pattern blocks; fail here instead of deep inside
  # compute_fspars()/correct_evfs().
  has_miss_patterns <- any(vapply(
    object@Data@Mp,
    function(m) !is.null(m),
    logical(1)
  ))
  if (has_miss_patterns && (corrected_fsT || reliability || vfsLT)) {
    stop(
      "'corrected_fsT', 'reliability', and 'vfsLT' are not supported when ",
      "the data contain missing values: lavaan scores the cases on their ",
      "observed-indicator patterns, and these SE paths require the full ",
      "indicator set. Fit on complete data or re-run with ",
      "corrected_fsT = FALSE, reliability = FALSE, vfsLT = FALSE.",
      call. = FALSE
    )
  }

  has_priors <- !is.null(prior_mean) || !is.null(prior_cov)
  if (has_priors) {
    if (method == "Bartlett") {
      stop(
        "'prior_mean'/'prior_cov' are only supported for 'regression' ",
        "(EB) scoring.",
        call. = FALSE
      )
    }
    if (reliability) {
      stop(
        "'reliability = TRUE' is not supported with user-supplied ",
        "'prior_mean'/'prior_cov'.",
        call. = FALSE
      )
    }
    est_first <- lavInspect(object, what = "est")
    if (object@Data@ngroups > 1) {
      est_first <- est_first[[1]]
    }
    lv_names <- colnames(est_first$lambda)
    priors <- validate_fs_priors(prior_mean, prior_cov, lv_names)
  } else {
    priors <- list(mean = NULL, cov = NULL)
  }

  if (reliability) {
    corrected_fsT <- TRUE
  }
  if (corrected_fsT) {
    add_to_evfs <- correct_evfs(
      object,
      method = method,
      psi_override = priors$cov
    )
  } else {
    # Direct slot access; see note in get_fs_blocks.lavaan() -- avoids
    # lavInspect()'s expensive per-call version check.
    add_to_evfs <- rep(0, object@Data@ngroups)
  }

  blocks_by_group <- get_fs_blocks.lavaan(
    object,
    method = method,
    add_to_evfs = add_to_evfs,
    prior_mean = priors$mean,
    prior_cov = priors$cov
  )

  group_var <- object@Data@group
  group_col <- if (length(group_var) > 0) group_var else NULL

  out <- assemble_fs_blocks(
    blocks_by_group,
    format = format,
    group_col = group_col
  )

  if (vfsLT) {
    attr(out, "vfsLT") <- vcov_ld_evfs(
      object,
      method = method,
      psi_override = priors$cov
    )
  }
  if (reliability) {
    est <- lavInspect(object, what = "est")
    group_labels <- object@Data@group.label
    ngroups <- object@Data@ngroups
    multifactor <- if (ngroups > 1) {
      length(attr(out, "fsb")[[1]]) > 1
    } else {
      length(attr(out, "fsb")) > 1
    }
    if (multifactor) {
      warning(
        "Computation of reliability for a multi-factor model is not ",
        "currently supported. "
      )
    } else {
      if (ngroups == 1) {
        is_std.lv <- all(est$psi == 1)
        attr(out, "reliability") <-
          compute_fsrel(object, method = method)[[1]]
      } else {
        is_std.lv <- all(unlist(lapply(est, function(x) x$psi)) == 1)
        rels <- compute_fsrel(object, method = method)
        # The norig slot is a per-group list; unlist to the vector form
        # that lavInspect(what = "norig") returns.
        group_n <- unlist(object@Data@norig)
        rels[ngroups + 1] <- sum(unlist(rels) * group_n / sum(group_n))
        attr(out, "reliability") <-
          setNames(rels, c(group_labels, "overall"))
      }
      if (!is_std.lv) {
        warning(
          "Computation of reliability may not be accurate when ",
          "the latent variables are not standardized. "
        )
      }
    }
  }
  out
}

get_fs_blocks.merMod <- function(object, legacy_names = FALSE, ...) {
  num_re <- length(object@cnms[[1]])
  base_names <- paste0("u", seq_len(num_re) - 1)
  re_names <- if (legacy_names) paste0(base_names, "_eb") else base_names
  fs_names <- paste0("fs_", re_names)

  # Grouping factor for the first RE term. lme4 canonicalizes its levels
  # (sorted for atomic vectors, user-specified order for pre-factorized
  # factors) and EVERY structure read below follows that level order --
  # not the row order of the data: split() on a factor, the Z columns,
  # the b/u random-effect vector, and ranef() rows are all indexed by
  # level position. Names must therefore come from the levels, never
  # from first-appearance order in the data.
  f1 <- as.factor(object@flist[[1]])
  n_clus <- nlevels(f1)

  # Row indices per cluster, in level order (split() follows factor
  # levels, so block j corresponds to level j regardless of row order).
  case_idx <- split(seq_len(length(f1)), f1)

  # Random-effects design Z (full n x sum_p matrix from lme4). Its columns
  # are ordered by RE term, then by level index of the first term's
  # grouping factor, so the first term's block for level (i.e. block) j
  # occupies columns (j - 1) * num_re + seq_len(num_re). Multi-term models
  # make Z wider; we slice the first bar only, matching the `[[1]]`
  # convention used for b/cnms/flist below. The random design -- not the
  # fixed design -- determines Kz: with Z != X the fixed-design code
  # produces non-conformable products.
  Zmat <- as.matrix(lme4::getME(object, "Z"))
  D <- get_D(object@theta)
  s <- stats::sigma(object)

  stopifnot(ncol(Zmat) >= num_re * n_clus)

  # EB scores for the first term: getME("b") = crossprod(Lambdat, u),
  # bit-identical to ranef(object)[[1]] values but computed level-major
  # without ranef()'s per-term work (cheaper for multi-term models).
  b <- lme4::getME(object, "b")
  stopifnot(length(b) >= num_re * n_clus)
  u_b <- matrix(
    b[seq_len(num_re * n_clus)],
    nrow = n_clus,
    ncol = num_re,
    byrow = TRUE
  )
  rownames(u_b) <- levels(f1)

  blocks <- vector("list", n_clus)

  for (j in seq_len(n_clus)) {
    idx <- case_idx[[j]]
    zj <- Zmat[idx, (j - 1L) * num_re + seq_len(num_re), drop = FALSE]
    Kz <- crossprod(zj)
    DKz <- D %*% Kz
    inv_W <- solve(DKz + diag(nrow(Kz)))
    fsL_j <- DKz - DKz %*% inv_W %*% DKz
    fsT_j <- s^2 * inv_W %*% DKz %*% D %*% t(inv_W)

    fs_row <- u_b[j, , drop = FALSE]
    colnames(fs_row) <- re_names

    # colnames = indicator/lv names, rownames = fs names (augment_fs convention)
    colnames(fsL_j) <- re_names
    rownames(fsL_j) <- fs_names
    attr(fs_row, "fsL") <- fsL_j

    rownames(fsT_j) <- colnames(fsT_j) <- fs_names

    # Scoring matrix: S_j %*% (y_j - X_j %*% beta) reproduces the EB scores
    # (ranef), where y_j/X_j are the cluster's rows of the model frame and
    # fixed-effects design. See vignettes/scoring-matrices.Rmd.
    scoring_matrix_j <- inv_W %*% D %*% t(zj)
    rownames(scoring_matrix_j) <- fs_names
    colnames(scoring_matrix_j) <- as.character(seq_len(nrow(zj)))

    blocks[[j]] <- list(
      case_idx = idx,
      fs = fs_row,
      fsL = fsL_j,
      fsT = fsT_j,
      fsb = NULL,
      scoring_matrix = scoring_matrix_j
    )
  }

  # Names in canonical level order (matches ranef() row names / @cnms).
  setNames(blocks, levels(f1))
}
get_D <- function(theta) {
  L_D <- lme4::vec2mlist(theta, symm = FALSE)[[1]]
  tcrossprod(L_D)
}

#' @rdname get_fs
#' @param fsm Currently not used.
#' @param format Output format: `"unified"` returns a single data frame with
#'        a `group` column; `"list"` returns a list of data frames per group.
#'        For `merMod` objects there is always a single implicit group, so
#'        `"list"` returns a bare data frame.
#' @param legacy_names Logical. Random-effect score naming convention for
#'        `merMod` objects. `FALSE` (default) uses `fs_u0`-style names
#'        (`fs_u0`/`fs_u1`/..., with loadings `u0_by_fs_u0` and error terms
#'        `ev_fs_u0`, `ecov_fs_u1_fs_u0`). `TRUE` reproduces the pre-refactor
#'        `u0_eb`-style *column names* (`u0_eb`, `u0_by_u0_eb`, `ev_u0_eb`,
#'        `ecov_u0_eb_u1_eb`) in the legacy column order. Note the legacy
#'        output is name-compatible, not byte-identical, with the
#'        pre-refactor `get_fs_lmer()` result: it additionally carries
#'        score-error columns (`u0_eb_se`, ...), per-cluster `fsL`/`fsT`
#'        array attributes, a per-cluster `scoring_matrix` list attribute
#'        (see [get_fs()]), and has NULL row names (the pre-refactor
#'        output had none of these and used the ranef subject IDs as row
#'        names).
#' @export
get_fs.merMod <- function(
  object,
  method = c("EB"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  fsm = FALSE,
  format = c("unified", "list"),
  legacy_names = FALSE,
  ...
) {
  if (!inherits(object, "merMod")) {
    stop("`object` must be an `lmerMod` model object.", call. = FALSE)
  }
  prior_dots <- list(...)[c("prior_mean", "prior_cov")]
  if (!is.null(prior_dots$prior_mean) || !is.null(prior_dots$prior_cov)) {
    stop(
      "'prior_mean'/'prior_cov' are not supported for `lmerMod` objects.",
      call. = FALSE
    )
  }

  blocks <- get_fs_blocks.merMod(object, legacy_names = legacy_names)

  # NOTE: merMod does NOT route through assemble_fs_blocks(). The shared
  # assembler assumes one data row per individual case (nrow = n_cases),
  # filling a template DataFrame by case_idx. merMod's semantics produce
  # one row per cluster (nrow = n_clusters), because each cluster has a
  # single EB estimate shared across its cases. Forcing it through the
  # assembler would produce n_cases rows with repeated values, breaking
  # backward compatibility with ranef()-aligned output. This is an
  # intentional architectural exception.
  aug_list <- lapply(blocks, function(b) {
    augment_fs(b$fs, b$fsT)
  })
  out <- do.call(rbind, aug_list)
  rownames(out) <- NULL

  # Legacy `u<k>_eb`-style names are not produced by augment_fs() (which
  # unconditionally prefixes scores with "fs_" and emits `ecov` columns in
  # row/col iteration order), so translate the column names when requested.
  if (legacy_names) {
    colnames(out) <- rename_legacy_fs_cols(colnames(out))
  }

  # Per-cluster fsL/fsT as array attributes (merMod-specific: each cluster
  # has its own covariance matrix, unlike lavaan groups where attributes are
  # shared within a group). Names come from the first block's fsL.
  n_clus <- length(blocks)
  fsL_1 <- blocks[[1]]$fsL
  re_names <- colnames(fsL_1)
  fs_names <- rownames(fsL_1)
  fsL_arr <- array(
    0,
    dim = c(nrow(fsL_1), ncol(fsL_1), n_clus),
    dimnames = list(fs_names, re_names, names(blocks))
  )
  fsT_arr <- array(
    0,
    dim = c(nrow(fsL_1), ncol(fsL_1), n_clus),
    dimnames = list(fs_names, fs_names, names(blocks))
  )
  for (j in seq_len(n_clus)) {
    fsL_arr[,, j] <- blocks[[j]]$fsL
    fsT_arr[,, j] <- blocks[[j]]$fsT
  }
  attr(out, "fsL") <- fsL_arr
  attr(out, "fsT") <- fsT_arr

  # Per-cluster scoring matrices as a named list (one p x n_j matrix per
  # cluster; list, not array, because cluster sizes can differ).
  attr(out, "scoring_matrix") <- setNames(
    lapply(blocks, function(b) b$scoring_matrix),
    names(blocks)
  )

  out
}

# Translate augment_fs() column names into the pre-refactor `u<k>_eb`-style
# names. augment_fs() unconditionally prefixes scores with "fs_" and emits
# `ecov` columns in matrix row/col iteration order; the legacy convention
# needs neither that prefix (scores/SEs/ev), nor a bare (no "_eb") indicator
# side on loadings, nor the ascending `ecov` order — hence the per-category
# rules below.
rename_legacy_fs_cols <- function(x) {
  vapply(
    x,
    function(nm) {
      if (grepl("_by_", nm, fixed = TRUE)) {
        # <ind>_eb_by_fs_<fs>_eb -> <ind>_by_<fs>_eb
        sub("^(u[0-9]+)_eb_by_fs_(u[0-9]+)_eb$", "\\1_by_\\2_eb", nm)
      } else if (grepl("^ecov_", nm)) {
        # ecov_fs_uA_eb_fs_uB_eb -> ecov_uA_eb_uB_eb, ascending in u<k>
        reorder_ecov_col(
          sub("^ecov_fs_(u[0-9]+)_eb_fs_(u[0-9]+)_eb$",
              "ecov_\\1_eb_\\2_eb",
              nm
          )
        )
      } else if (grepl("^ev_", nm)) {
        # ev_fs_u0_eb -> ev_u0_eb
        sub("^ev_fs_", "ev_", nm)
      } else {
        # scores (fs_u0_eb) and SEs (fs_u0_eb_se): drop the leading "fs_"
        sub("^fs_", "", nm)
      }
    },
    character(1),
    USE.NAMES = FALSE
  )
}

# Reorder `ecov_uA_eb_uB_eb` so the u<k> indices are ascending (A < B).
# Comparing the extracted integers (not the strings) keeps this correct for
# 10+ random effects as well.
reorder_ecov_col <- function(nm) {
  m <- regmatches(nm, regexec("^ecov_u([0-9]+)_eb_u([0-9]+)_eb$", nm))[[1]]
  if (length(m) != 3L) return(nm)
  a <- as.integer(m[2])
  b <- as.integer(m[3])
  if (a <= b) nm else paste0("ecov_u", b, "_eb_u", a, "_eb")
}
