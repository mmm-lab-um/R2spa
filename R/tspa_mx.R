#' Two-Stage Path Analysis (OpenMx)
#'
#' Fit a two-stage path analysis (2S-PA) model in OpenMx. This is the OpenMx
#' counterpart of [tspa()]: the structural path model is expressed on the
#' (true) latent factors, each factor score is a single indicator with known
#' measurement error, and --- unlike `lavaan::tspa()` --- the measurement
#' quantities (loadings / error variances / intercepts) may be fixed
#' per-group constants \emph{or} per-observation definition-variable columns.
#' The per-row definition-variable form is the exact (non-pooled) correction
#' that `lavaan::tspa(reduce = )` only approximates.
#'
#' The internal design is a single-level RAM model (no sub-model, no `umx`):
#' the lavaan structural string is parsed with [lavaan::lavaanify()], the
#' corrected latents are given auto latent variances, the score indicators and
#' their errors/intercepts are attached per [tspa()'s](
#' https://www.rdocumentation.org/packages/R2spa) schema, and the whole thing
#' is fit with `mxFitFunctionML()` (raw-data FIML).
#'
#' @param model A character string describing the structural path model in
#'   `lavaan` syntax, using the **latent** (factor) names. Phase 1 restricts
#'   every variable in `model` to a corrected latent (one that has a factor
#'   score). Latent variances are added automatically, so do not declare them
#'   here.
#' @param data A data frame carrying the factor-score columns
#'   (`fs_<latent>`) and, for definition-variable entries, the per-observation
#'   columns they reference. [`fs_indiv()`] on a [get_fs()] result produces
#'   exactly this table. Definition-variable columns must be free of `NA`.
#' @param se_fs A named numeric vector of standard errors (one per latent) for
#'   the single-score-per-latent case; implies fixed unit loadings and error
#'   variances `se_fs^2`.
#' @param fsL A `q x p` loading matrix including cross-loadings: rows = score
#'   names (`fs_<latent>`), columns = latent names. Each cell is either a
#'   number (fixed loading) or a character naming a definition-variable column.
#' @param fsT A `q x q` error variance-covariance matrix over the score names;
#'   the lower triangle (incl. diagonal) is used. Each cell is a number
#'   (fixed) or a character naming a definition-variable column.
#' @param fsb A vector of score intercepts (length `q`, named by score, either
#'   order) --- each a number (fixed) or a definition-variable column name.
#'   `NULL` (default) fixes all score intercepts at zero.
#' @param ... Additional arguments passed on to [`OpenMx::mxRun()`]
#'   (e.g. `intervals = TRUE`).
#' @return A fitted `OpenMx` `MxModel`. `coef()`, `vcov()`, and `summary()`
#'   work as usual.
#'
#' @importFrom OpenMx mxModel mxData mxPath mxFitFunctionML mxRun
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Per-row correction straight from a get_fs() result via fs_indiv():
#' fit <- cfa("dem60 =~ y1 + y2 + y3 + y4; ind60 =~ x1 + x2 + x3",
#'            data = PoliticalDemocracy)
#' fs  <- get_fs(PoliticalDemocracy, "dem60 =~ y1 + y2 + y3 + y4
#'                                      ind60 =~ x1 + x2 + x3")
#' dat <- fs_indiv(fs, include_intercept = TRUE)
#' tspa_mx_model("dem60 ~ ind60; dem60 + ind60 ~ 1",
#'   data = dat,
#'   fsL = matrix(c("ind60_by_fs_ind60", "ind60_by_fs_dem60",
#'                  "dem60_by_fs_ind60", "dem60_by_fs_dem60"),
#'                nrow = 2, dimnames = list(c("fs_ind60", "fs_dem60"),
#'                                          c("ind60", "dem60"))),
#'   fsT = matrix(c("ev_fs_ind60", "ecov_fs_ind60_fs_dem60", NA,
#'                  "ev_fs_dem60"),
#'                nrow = 2, dimnames = list(c("fs_ind60", "fs_dem60"),
#'                                          c("fs_ind60", "fs_dem60"))),
#'   fsb = c(fs_ind60 = "int_fs_ind60", fs_dem60 = "int_fs_dem60"))
#' }

tspa_mx_model <- function(model, data, se_fs = NULL, fsL = NULL,
                          fsT = NULL, fsb = NULL, ...) {
  if (!is.character(model)) {
    stop("The structural path model 'model' must be a lavaan syntax string.",
         call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.null(se_fs) && (!is.null(fsL) || !is.null(fsT))) {
    stop("Provide either 'se_fs' or ('fsL' and 'fsT'), not both.",
         call. = FALSE)
  }
  if (xor(is.null(fsL), is.null(fsT))) {
    stop("Provide both 'fsL' and 'fsT', or use 'se_fs'.", call. = FALSE)
  }
  if ((is.list(fsL) && !is.matrix(fsL) && length(fsL) > 1L) ||
      (is.list(fsT) && !is.matrix(fsT) && length(fsT) > 1L)) {
    stop("Multigroup 'fsL'/'fsT' are not supported yet (Phase 1 is single-group).",
         call. = FALSE)
  }

  spec <- tspa_mx_spec(se_fs, fsL, fsT, fsb)

  # Data-contract guards (score + definition-variable columns must be present,
  # and definition-variable columns must be NA-free for OpenMx).
  if (!all(spec$scores %in% names(data))) {
    stop("data is missing factor-score column(s): ",
         paste(setdiff(spec$scores, names(data)), collapse = ", "), ".",
         call. = FALSE)
  }
  dv_cols <- unique(c(unlist(spec$L$coln), unlist(spec$T$coln),
                      unlist(spec$b$coln)))
  dv_cols <- dv_cols[!is.na(dv_cols)]
  if (!all(dv_cols %in% names(data))) {
    stop("data is missing definition-variable column(s): ",
         paste(setdiff(dv_cols, names(data)), collapse = ", "), ".",
         call. = FALSE)
  }
  na_dv <- dv_cols[vapply(dv_cols, function(col) anyNA(data[[col]]), logical(1L))]
  if (length(na_dv)) {
    stop("Definition-variable column(s) contain NA (",
         paste(na_dv, collapse = ", "),
         "): OpenMx definition variables must be complete for every row.",
         call. = FALSE)
  }

  full <- lav_to_mx_ram(tspa_mx_model_string(model, spec), spec, data)
  mxRun(full, ...)
}

# --- Measurement specification ------------------------------------------------

# Unwrap a length-1 list (unified single-group get_fs() attribute shape) so a
# plain matrix is used; any longer list is a (unsupported) multigroup input.
tspa_mx_unwrap <- function(x) {
  if (is.list(x) && !is.matrix(x)) {
    if (length(x) == 1L) return(x[[1L]])
    stop("Length-", length(x), " list 'fsL'/'fsT' is not supported yet (Phase 1 is single-group).",
         call. = FALSE)
  }
  x
}

# Numeric -> values, character -> definition-variable column names, in one
# fixed/defvar cell table. Absence (NA) is preserved in both slots.
tspa_mx_cells <- function(m, arg) {
  if (is.null(dim(m))) m <- matrix(m)
  vals <- array(NA_real_, dim = dim(m), dimnames = dimnames(m))
  coln <- array(NA_character_, dim = dim(m), dimnames = dimnames(m))
  if (is.numeric(m)) vals <- m
  else if (is.character(m)) coln <- m
  else stop("'", arg, "' must be a numeric (fixed) or character (definition-variable) matrix.",
            call. = FALSE)
  list(vals = vals, coln = coln)
}

# Reorder (and align) an error matrix to the score order `S`.
tspa_mx_align_scores <- function(T, S, arg) {
  if (is.null(rownames(T)) || length(S) > 0L && !all(S %in% rownames(T))) {
    stop("'", arg, "' rows must be named by the factor-score names.", call. = FALSE)
  }
  T <- T[rownames = S, , drop = FALSE]
  T <- T[, colnames = S, drop = FALSE]
}

tspa_mx_spec <- function(se_fs, fsL, fsT, fsb) {
  if (!is.null(se_fs)) {
    if (is.null(names(se_fs)) || anyNA(names(se_fs))) {
      stop("'se_fs' must be named by latent name.", call. = FALSE)
    }
    se <- as.numeric(se_fs)
    V <- names(se_fs)
    S <- paste0("fs_", V)
    Lm <- matrix(NA_real_, length(V), length(V), dimnames = list(S, V)); diag(Lm) <- 1
    Tm <- matrix(NA_real_, length(V), length(V), dimnames = list(S, S)); diag(Tm) <- se^2
    L <- tspa_mx_cells(Lm, "fsL")
    T <- tspa_mx_cells(Tm, "fsT")
  } else {
    fsL <- tspa_mx_unwrap(fsL)
    fsT <- tspa_mx_unwrap(fsT)
    if (is.null(rownames(fsL)) || anyNA(rownames(fsL))) {
      stop("'fsL' rows must be named by the factor-score names.", call. = FALSE)
    }
    if (is.null(colnames(fsL)) || anyNA(colnames(fsL))) {
      stop("'fsL' columns must be named by the latent names.", call. = FALSE)
    }
    S <- rownames(fsL)
    V <- colnames(fsL)
    fsT <- tspa_mx_align_scores(fsT, S, "fsT")
    L <- tspa_mx_cells(fsL, "fsL")
    T <- tspa_mx_cells(fsT, "fsT")
  }
  if (!is.null(fsb)) {
    fsb <- tspa_mx_unwrap(fsb)
    if (is.matrix(fsb)) {
      b <- if (nrow(fsb) == 1L) as.matrix(fsb) else t(fsb)
      if (!setequal(colnames(b), S)) {
        stop("'fsb' must be named by all factor-score names.", call. = FALSE)
      }
      b <- b[, S, drop = FALSE]
    } else {
      nm <- names(fsb)
      if (is.null(nm) || anyNA(nm) || !setequal(nm, S)) {
        stop("'fsb' must be a named vector over all factor-score names.",
             call. = FALSE)
      }
      b <- matrix(unname(fsb[S]), nrow = 1L)
    }
    colnames(b) <- S
    b <- tspa_mx_cells(b, "fsb")
  } else {
    b <- list(
      vals = matrix(0, 1L, length(S), dimnames = list(NULL, S)),
      coln = matrix(NA_character_, 1L, length(S), dimnames = list(NULL, S))
    )
  }
  list(scores = S, latents = V, L = L, T = T, b = b)
}

# --- Full lavaan model string (structural + measurement) ---------------------

# cell -> "c(value)" for a fixed value, "c(1)" sentinel for a definition
# variable (overlaid later), NA when the cell is absent.
tspa_mx_cellval <- function(vals, coln, i, j) {
  if (!is.na(vals[i, j])) sprintf("c(%s)", vals[i, j])
  else if (!is.na(coln[i, j])) "c(1)"
  else NA_character_
}

tspa_mx_model_string <- function(model, spec) {
  S <- spec$scores
  V <- spec$latents

  pt1 <- lavaan::lavaanify(model)
  struct_var <- setdiff(unique(c(pt1$lhs, pt1$rhs)), c("", "1", "one"))
  bad <- setdiff(struct_var, V)
  if (length(bad)) {
    stop("Every variable in the structural model must be a corrected latent ",
         "(present in 'fsL'); unexpected: ", paste(bad, collapse = ", "), ".",
         call. = FALSE)
  }
  lines <- character()

  for (k in seq_along(V)) {
    terms <- character()
    for (i in seq_along(S)) {
      cv <- tspa_mx_cellval(spec$L$vals, spec$L$coln, i, k)
      if (!is.na(cv)) terms <- c(terms, paste0(cv, " * ", S[i]))
    }
    if (!length(terms)) {
      stop("Latent '", V[k], "' has no factor-score indicator in 'fsL'.",
           call. = FALSE)
    }
    lines <- c(lines, paste(V[k], "=~", paste(terms, collapse = " + ")))
  }

  for (i in seq_along(S)) {
    for (j in seq_len(i)) {
      cv <- tspa_mx_cellval(spec$T$vals, spec$T$coln, i, j)
      if (!is.na(cv)) lines <- c(lines, paste0(S[i], " ~~ ", cv, " * ", S[j]))
    }
  }
  # Score means are modelled only as per-row definition-variable columns; a
  # fixed/absent mean leaves the observed score mean at its data value.
  for (i in seq_along(S)) {
    if (!is.na(spec$b$coln[1L, i]) && nzchar(spec$b$coln[1L, i])) {
      lines <- c(lines, paste0(S[i], " ~ c(1) * 1"))
    }
  }

  paste(c(model, paste(lines, collapse = "\n")), collapse = "\n")
}

# --- lavaan partable -> OpenMx RAM model --------------------------------------

# One partable row -> (from, to, arrows) for OpenMx.
tspa_mx_op_map <- function(r) {
  if (r$op == "~")        list(f = r$rhs, t = r$lhs, ar = 1L)
  else if (r$op == "~1")  list(f = "one", t = r$lhs, ar = 1L)
  else if (r$op == "=~")  list(f = r$lhs, t = r$rhs, ar = 1L)
  else if (r$op == "~~")  list(f = r$lhs, t = r$rhs, ar = 2L)
  else stop("Unsupported lavaan operator '", r$op, "'.", call. = FALSE)
}

# Sensible starting value when the partable leaves a free parameter without a
# start (a bare start is NA): variance -> 1, mean -> 0, path -> 0.1.
tspa_mx_default_start <- function(op, ustart) {
  if (!is.na(ustart)) return(ustart)
  if (op == "~~") 1 else if (op == "~1") 0 else 0.1
}

# For the measurement rows, the definition-variable column backing a given
# (lhs, op, rhs) cell --- NA if the cell is fixed or absent.
tspa_mx_defvar_col <- function(spec, lhs, op, rhs) {
  si <- match(lhs, spec$scores)
  sj <- match(rhs, spec$scores)
  li <- match(lhs, spec$latents)
  if (op == "=~") {
    return(if (is.na(li) || is.na(sj)) NA_character_ else spec$L$coln[sj, li])
  }
  if (op == "~~") {
    return(if (is.na(si) || is.na(sj)) NA_character_ else spec$T$coln[si, sj])
  }
  if (op == "~1") {
    return(if (is.na(si)) NA_character_ else spec$b$coln[1L, si])
  }
  NA_character_
}

# Partable -> one mxPath per row, overlaying definition-variable labels for the
# measurement cells. Returns the list of mxPath objects.
tspa_mx_paths <- function(pt, spec) {
  lapply(seq_len(nrow(pt)), function(i) {
    r  <- pt[i, , drop = FALSE]
    mm <- tspa_mx_op_map(r)
    dv <- tspa_mx_defvar_col(spec, r$lhs, r$op, r$rhs)
    if (!is.na(dv) && nzchar(dv)) {
      return(mxPath(from = mm$f, to = mm$t, arrows = mm$ar, free = FALSE,
                    values = 0, labels = paste0("data.", dv)))
    }
    # lavaanify() auto-seeds structural-latent variances (and means, when a mean
    # structure is present) as free == 0, ustart == 0; release them.
    if (r$free == 0L && is.numeric(r$ustart) && r$ustart == 0L &&
        r$lhs %in% spec$latents && (mm$ar == 2L && mm$f == mm$t || mm$f == "one")) {
      return(mxPath(from = mm$f, to = mm$t, arrows = mm$ar, free = TRUE,
                    values = if (mm$ar == 2L) 1 else 0))
    }
    free <- r$free > 0L
    args <- list(from = mm$f, to = mm$t, arrows = mm$ar, free = free,
                 values = if (free) tspa_mx_default_start(r$op, r$ustart)
                          else r$ustart)
    lab <- r$label
    if (is.character(lab) && !is.na(lab) && nzchar(lab)) args$labels <- lab
    do.call(mxPath, args)
  })
}

# Convert a fully-specified lavaan model string into a complete OpenMx RAM
# model (data + ML fit function): latents are the `=~` left-hand sides,
# manifests are everything else, and each partable row becomes an mxPath.
lav_to_mx_ram <- function(model_str, spec, data) {
  pt        <- lavaan::lavaanify(model_str)
  latents   <- unique(pt$lhs[pt$op == "=~"])
  all_vars  <- setdiff(unique(c(pt$lhs, pt$rhs)), c(NA_character_, "", "1", "one"))
  manifests <- setdiff(all_vars, latents)
  paths <- tspa_mx_paths(pt, spec)
  # type = "raw" (FIML) requires an explicit mean for every variable; add a
  # free (start 0) mean wherever the structural model omits one.
  need_mean <- setdiff(c(manifests, latents), pt$lhs[pt$op == "~1"])
  if (length(need_mean)) {
    paths <- c(paths, lapply(need_mean, function(v)
      mxPath(from = "one", to = v, free = TRUE, values = 0)))
  }
  mxModel("m1", type = "RAM", manifestVars = manifests, latentVars = latents,
          unlist(paths, recursive = FALSE),
          mxData(observed = data, type = "raw"),
          mxFitFunctionML())
}
