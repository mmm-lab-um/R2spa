# S3 methods for get_fs()
#
# Future methods (e.g. get_fs.mirt()) should be added to this file.

normalize_fs_method <- function(method) {
  method <- match.arg(method, c("regression", "Bartlett", "ML", "EB"))
  switch(method, ML = "Bartlett", EB = "regression", method)
}

#' @rdname get_fs
#' @export
get_fs.data.frame <- function(
  data,
  model = NULL,
  group = NULL,
  method = c("regression", "Bartlett", "ML", "EB"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  format = c("unified", "list"),
  ...
) {
  if (is.null(model)) {
    ind_names <- colnames(data)
    if (!is.null(group)) {
      ind_names <- setdiff(ind_names, group)
    }
    model <- paste("f1 =~", paste(ind_names, collapse = " + "))
  }
  fit <- cfa(model, data = data, group = group, ...)
  get_fs(
    fit,
    method = method,
    corrected_fsT = corrected_fsT,
    vfsLT = vfsLT,
    reliability = reliability,
    format = format
  )
}

get_fs_blocks.lavaan <- function(object, method, add_to_evfs, ...) {
  method <- match.arg(method, c("regression", "Bartlett"))
  est <- lavInspect(object, what = "est")
  y <- lavInspect(object, what = "data")
  miss_pat <- object@Data@Mp

  prepare_fs <- function(y, est, add, mp, method) {
    if (is.null(mp)) {
      fscore <-
        compute_fscore(
          y,
          lambda = est$lambda,
          theta = est$theta,
          psi = est$psi,
          nu = est$nu,
          alpha = est$alpha,
          method = method,
          fs_matrices = TRUE
        )
      list(
        case_idx = seq_len(nrow(y)),
        fs = fscore,
        fsT = attr(fscore, "fsT") + add,
        fsL = attr(fscore, "fsL"),
        fsb = attr(fscore, "fsb"),
        scoring_matrix = attr(fscore, "scoring_matrix")
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
            psi = est$psi,
            nu = est$nu[pat_m, , drop = FALSE],
            alpha = est$alpha,
            method = method,
            fs_matrices = TRUE
          )
        blocks[[m]] <- list(
          case_idx = idx_m,
          fs = fs_m,
          fsT = attr(fs_m, "fsT") + add,
          fsL = attr(fs_m, "fsL"),
          fsb = attr(fs_m, "fsb"),
          scoring_matrix = attr(fs_m, "scoring_matrix")
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
#' @param object A fitted model object. For [get_fs()], the first argument is
#'        named `data` (dispatch argument). Methods accept `object` as their
#'        first parameter name internally.
#' @param format Output format: `"unified"` returns a single data frame with
#'        a `group` column; `"list"` returns a list of data frames per group.
#' @export
get_fs.default <- function(data, ...) {
  if (is.matrix(data)) {
    data <- as.data.frame(data)
    return(get_fs(data, ...))
  }
  stop(
    "get_fs() is not implemented for objects of class '",
    paste(class(data), collapse = "', '"),
    "'. Currently supported: 'data.frame', 'lavaan', ",
    "and 'lmerMod'. Support for 'mirt' models is planned.",
    call. = FALSE
  )
}

#' @rdname get_fs
#' @param object A fitted model object.
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
  ...
) {
  method <- normalize_fs_method(method)
  if (!inherits(object, "lavaan")) {
    stop("`object` must be a `lavaan` model object.")
  }
  format <- match.arg(format)

  if (reliability) {
    corrected_fsT <- TRUE
  }
  if (corrected_fsT) {
    add_to_evfs <- correct_evfs(object, method = method)
  } else {
    # Direct slot access; see note in get_fs_blocks.lavaan() -- avoids
    # lavInspect()'s expensive per-call version check.
    add_to_evfs <- rep(0, object@Data@ngroups)
  }

  blocks_by_group <- get_fs_blocks.lavaan(
    object,
    method = method,
    add_to_evfs = add_to_evfs
  )

  group_var <- object@Data@group
  group_col <- if (length(group_var) > 0) group_var else NULL

  out <- assemble_fs_blocks(
    blocks_by_group,
    format = format,
    group_col = group_col
  )

  if (vfsLT) {
    attr(out, "vfsLT") <- vcov_ld_evfs(object, method = method)
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

  u_eb <- as.data.frame(lme4::ranef(object)[[1]])
  colnames(u_eb) <- re_names

  cluster_levels <- unique(as.character(object@flist[[1]]))
  mf <- model.frame(object)
  case_idx <- split(seq_len(nrow(mf)), object@flist[[1]])

  Xlist <- split(data.frame(object@pp$X), object@flist[[1]])
  D <- get_D(object@theta)
  s <- stats::sigma(object)

  n_clus <- length(Xlist)
  blocks <- vector("list", n_clus)

  for (j in seq_len(n_clus)) {
    xj <- as.matrix(Xlist[[j]])
    Kz <- crossprod(xj)
    DKz <- D %*% Kz
    inv_W <- solve(DKz + diag(nrow(Kz)))
    fsL_j <- DKz - DKz %*% inv_W %*% DKz
    fsT_j <- s^2 * inv_W %*% DKz %*% D %*% t(inv_W)

    fs_row <- as.matrix(u_eb[j, , drop = FALSE])
    colnames(fs_row) <- re_names

    # colnames = indicator/lv names, rownames = fs names (augment_fs convention)
    colnames(fsL_j) <- re_names
    rownames(fsL_j) <- fs_names
    attr(fs_row, "fsL") <- fsL_j

    rownames(fsT_j) <- colnames(fsT_j) <- fs_names

    blocks[[j]] <- list(
      case_idx = case_idx[[j]],
      fs = fs_row,
      fsL = fsL_j,
      fsT = fsT_j,
      fsb = NULL,
      scoring_matrix = NULL
    )
  }

  setNames(blocks, cluster_levels)
}
get_D <- function(theta) {
  L_D <- lme4::vec2mlist(theta, symm = FALSE)[[1]]
  tcrossprod(L_D)
}

#' @rdname get_fs
#' @param object A fitted model object.
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
#'        array attributes, and has NULL row names (the pre-refactor
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
