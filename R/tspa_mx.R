#' Two-Stage Path Analysis (OpenMx)
#'
#' Fit a two-stage path analysis (2S-PA) model in OpenMx. This is the OpenMx
#' counterpart of [tspa()]: the structural path model is expressed on the
#' (true) latent factors, each factor score is a single indicator with known
#' measurement error, and --- unlike [tspa()] --- the measurement
#' quantities (loadings / error variances / intercepts) may be fixed
#' constants \emph{or} per-observation definition-variable columns.
#' The per-row definition-variable form is the exact (non-pooled) correction
#' that [tspa()] (with `reduce = `) only approximates.
#'
#' The internal design is a single-level RAM model (no sub-model, no `umx`):
#' the lavaan structural string is parsed with [lavaan::lavaanify()], the
#' corrected latents are given auto latent variances, the score indicators and
#' their errors/intercepts are attached per `tspa()`'s schema, and the whole thing
#' is fit with `mxFitFunctionML()` (raw-data FIML).
#'
#' ## Auto-derivation from a [get_fs()] result
#'
#' `tspa_mx_model()` derives its measurement inputs from a [get_fs()]
#' result passed as `data`, so the canonical call is
#' `tspa_mx_model(model, data = get_fs(...))`. Derivation fires only
#' when all of `se_fs`, `fsL`, `fsT`, and `fsb` are omitted; explicit
#' arguments always win. Derivation is provenance-gated: `data` must
#' resolve as a [get_fs()] result, so a hand-rolled frame with plain
#' matrix `fsT`/`fsL` attributes but no such provenance is rejected
#' with an informative error.
#'
#' Attributes are dispatched by shape: constant quantities (plain
#' matrices --- complete-data, `local = TRUE`, and `format = "list"`
#' results; a plain `fsb` vector) become fixed numeric matrices, while
#' per-row quantities (`mirt_per_obs`/`per_obs`-marked results),
#' per-cluster quantities (`merMod` results: 3-D `fsL`/`fsT` arrays, one
#' row per cluster), and per-pattern quantities (single-group FIML, keyed
#' by `fs_pattern$label`) become definition-variable matrices referencing
#' the result's own `*_by_*`, `ev_*`, and `ecov_*` columns. A `merMod`
#' result carries no `fsb` attribute, so its score intercepts use the default
#' mean handling. A per-unit `fsb` attribute appends `int_fs_*` intercept
#' columns to a working copy of `data` (the `fs_indiv(include_intercept = TRUE)`
#' equivalent); a constant numeric `fsb` is fixed directly. Unlike
#' [tspa()] (whose `reduce = ` argument pools per-unit quantities),
#' there is no pooling here --- the OpenMx route is exact-or-fail. A
#' multigroup result (the `group_col` attribute) is refused with the
#' Phase-1 message; a `mirt` multigroup result (a `group` column, no
#' `group_col` attribute) derives as a single pooled per-row-corrected
#' fit (no per-group structural parameters; the `group` column is
#' inert). A `get_fs(product = )` result derives identically: its extra
#' `fs_a:fs_b` (and `_se`/`_ld`) columns are inert to derivation, but the
#' `:` they carry is illegal in `OpenMx::mxData()` column names, so the
#'   frame is un-fittable until the product columns are dropped (a
#'   pre-existing limitation that affects the explicit-argument route too).
#'
#' ## Mean structure
#'
#' When both routes specify the same mean structure, their covariance and
#' structural quantities (paths, latent variances/covariances, and their SEs)
#' agree to optimizer tolerance. Without a known score intercept, only the
#' total observed mean is identified: [tspa()] and this OpenMx route may
#' report different splits between latent and score means. An explicit numeric
#' `fsb` fixes each score intercept, including zero; a character `fsb` fixes
#' it per observation through a definition variable. With fixed score
#' intercepts and an identified mean structure, compare implied means rather
#' than an endogenous latent's structural intercept alone.
#'
#' @param model A character string describing the structural path model in
#'   `lavaan` syntax, using the **latent** (factor) names. Phase 1 restricts
#'   every variable in `model` to a corrected latent (one that has a factor
#'   score). Latent variances are added automatically, so do not declare them
#'   here.
#' @param data A data frame carrying the factor-score columns
#'   (`fs_<latent>`) and, for definition-variable entries, the per-observation
#'   columns they reference. A [get_fs()] result works directly: with
#'   `se_fs`/`fsL`/`fsT`/`fsb` omitted, the measurement inputs are derived
#'   from its attributes (see Details), and the `int_fs_*` intercept
#'   columns are appended automatically --- `fs_indiv()` is no longer
#'   needed to obtain them. [`fs_indiv()`] on a [get_fs()] result produces
#'   the equivalent fully explicit table. Definition-variable columns must
#'   be numeric and free of `NA`.
#' @param se_fs Standard errors for the single-score-per-latent case, one per
#'   latent, named by latent name. Either (a) a named numeric vector of fixed
#'   SEs (the `k`-th latent's error variance is `se_fs[k]^2`), or (b) a named
#'   character vector whose entries name per-observation SE columns of `data`
#'   (the square root of the measurement error); the model then uses the
#'   column's square, `(data[, col])^2`, as that latent's error variance via an
#'   OpenMx definition variable. In form (b) the fitted model's data gains one
#'   working column per named SE column, `<col>_R2spa_se2` (holding
#'   `data[, col]^2`; it is inert except as a definition-variable source).
#'   Both forms imply fixed unit loadings. An
#'   explicit `se_fs` always wins over derivation; when omitted (along with
#'   `fsL`, `fsT`, and `fsb`), the measurement inputs are derived from a
#'   [get_fs()] result passed as `data` (see Details).
#' @param fsL A `q x p` loading matrix including cross-loadings: rows = score
#'   names (`fs_<latent>`), columns = latent names. The matrix must be
#'   uniformly numeric (every cell a fixed loading) or uniformly character
#'   (every cell a definition-variable column name); mixing fixed values and
#'   column names in one matrix is not supported.
#'   Or omitted, in which case the value is derived from the `fsL` attribute
#'   of a [get_fs()] result passed as `data` (see Details); an explicit `fsL`
#'   always wins.
#' @param fsT A `q x q` error variance-covariance matrix over the score names;
#'   the lower triangle (incl. diagonal) is used, and every score must have
#'   an error variance (a complete diagonal). The matrix must be uniformly
#'   numeric (every cell fixed) or uniformly character (every cell a
#'   definition-variable column name); mixing fixed values and column names in
#'   one matrix is not supported.
#'   Or omitted, in which case the value is derived from the `fsT` attribute
#'   of a [get_fs()] result passed as `data` (see Details); an explicit `fsT`
#'   always wins.
#' @param fsb A vector of score intercepts (length `q`, named by score, either
#'   order), uniformly numeric (every entry fixed) or uniformly character
#'   (every entry a definition-variable column name); mixing fixed values and
#'   column names is not supported.
#'   `NULL` (the default) omits an explicit score-intercept constraint: with
#'   a structural mean model the score intercepts are auto-seeded at zero;
#'   without one, raw-data FIML adds free score means. Supply a numeric zero
#'   to fix an intercept at zero regardless of the structural syntax.
#'   When `fsb` is omitted together with `se_fs`, `fsL`, and `fsT`, its value
#'   is derived from the `fsb` attribute of a [get_fs()] result passed as
#'   `data` (see Details): numeric intercepts stay fixed and per-unit
#'   intercepts use appended `int_fs_*` definition-variable columns. An
#'   explicit `fsb` always wins.
#' @param ... Additional arguments passed on to [`OpenMx::mxRun()`]
#'   (e.g. `intervals = TRUE`).
#' @return A fitted `OpenMx` `MxModel`. `coef()`, `vcov()`, and `summary()`
#'   work as usual.
#'
#' @note `OpenMx` is an optional dependency (a `Suggests`). `tspa_mx_model()`
#'   requires it at run time and stops with a message pointing at
#'   `install.packages("OpenMx")` if it is not installed. The `lavaan`-based
#'   [tspa()] route does not need `OpenMx`.
#'
#' @seealso
#' - `vignette("tspa-vignette-mx", package = "R2spa")` for the OpenMx route and IRT stage 1.
#' - `vignette("R2spa", package = "R2spa")` for the lavaan route.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Measurement inputs derived from a get_fs() result:
#' fs <- get_fs(PoliticalDemocracy, "dem60 =~ y1 + y2 + y3 + y4
#'                                   ind60 =~ x1 + x2 + x3")
#' # measurement inputs derived from the get_fs() result
#' tspa_mx_model("dem60 ~ ind60; dem60 + ind60 ~ 1", data = fs)
#'
#' ## Equivalent, fully explicit (per-row definition variables via
#' ## fs_indiv()):
#' dat <- fs_indiv(fs, include_intercept = TRUE)
#' tspa_mx_model("dem60 ~ ind60; dem60 + ind60 ~ 1",
#'   data = dat,
#'   fsL = matrix(c("ind60_by_fs_ind60", "ind60_by_fs_dem60",
#'                  "dem60_by_fs_ind60", "dem60_by_fs_dem60"),
#'                nrow = 2, dimnames = list(c("fs_ind60", "fs_dem60"),
#'                                          c("ind60", "dem60"))),
#'   fsT = matrix(c("ev_fs_ind60", "ecov_fs_dem60_fs_ind60", NA,
#'                  "ev_fs_dem60"),
#'                nrow = 2, dimnames = list(c("fs_ind60", "fs_dem60"),
#'                                          c("fs_ind60", "fs_dem60"))),
#'   fsb = c(fs_ind60 = "int_fs_ind60", fs_dem60 = "int_fs_dem60"))
#'
#' ## se_fs shorthand for the unit-loading case above: the named entries point
#' ## at per-row SE (sqrt of measurement error) columns of the fs_indiv()
#' ## table; the model uses their square as each score's error variance.
#' tspa_mx_model("dem60 ~ ind60; dem60 + ind60 ~ 1",
#'   data = dat,
#'   se_fs = c(ind60 = "fs_ind60_se", dem60 = "fs_dem60_se"))
#' }

tspa_mx_model <- function(model, data, se_fs = NULL, fsL = NULL,
                          fsT = NULL, fsb = NULL, ...) {
  require_openmx()
  # Record explicit inputs before deriving measurement attributes.
  se_fs_given <- !is.null(se_fs)
  fsb_given <- !is.null(fsb)
  if (!is.character(model)) {
    stop("The structural path model 'model' must be a lavaan syntax string.",
         call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!nrow(data)) {
    stop("'data' must have at least one row.", call. = FALSE)
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

  # Derive only when all measurement inputs are omitted.
  if (!se_fs_given && is.null(fsL) && is.null(fsT) && !fsb_given) {
    derived <- tspa_mx_derive_measurement(data)
    if (is.null(derived$fsL)) {
      stop(
        "No measurement inputs found for the factor scores in 'data'. ",
        "Please supply one of: (1) 'se_fs' (single-factor), (2) 'fsL' and ",
        "'fsT', or (3) a get_fs() result as 'data' (its attributes carry ",
        "the inputs).",
        if (is.null(derived$prov_err)) {
          ""
        } else {
          paste0(
            " The 'fsT'/'fsL' attributes on 'data' were detected but not ",
            "used because the data does not look like a get_fs() result: ",
            derived$prov_err
          )
        },
        call. = FALSE
      )
    }
    fsL <- derived$fsL
    fsT <- derived$fsT
    fsb <- derived$fsb
    data <- derived$data
  }

  # se_fs definition variables: a character entry is a per-row SE (the square
  # root of the measurement error) column; the RAM model needs its square as the
  # error variance, so append a squared working column and point the entry at
  # it (numeric fixed SEs pass through unchanged).
  if (!is.null(se_fs)) {
    resolved_se_fs <- tspa_mx_resolve_se_fs(se_fs, data)
    se_fs <- resolved_se_fs$se_fs
    data <- resolved_se_fs$data
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
  bad_type <- dv_cols[vapply(dv_cols, function(col) !is.numeric(data[[col]]),
                             logical(1L))]
  if (length(bad_type)) {
    stop("Definition-variable column(s) must be numeric (",
         paste(bad_type, collapse = ", "), ").", call. = FALSE)
  }
  na_dv <- dv_cols[vapply(dv_cols, function(col) anyNA(data[[col]]), logical(1L))]
  if (length(na_dv)) {
    stop("Definition-variable column(s) contain NA (",
         paste(na_dv, collapse = ", "),
         "): OpenMx definition variables must be complete for every row.",
         call. = FALSE)
  }

  full <- lav_to_mx_ram(tspa_mx_model_string(model, spec), spec, data)
  OpenMx::mxRun(full, ...)
}

# OpenMx is optional; guard the entry point before using its namespace.
require_openmx <- function() {
  if (!requireNamespace("OpenMx", quietly = TRUE)) {
    stop(
      "'OpenMx' is required for the OpenMx stage-2 route ",
      "(tspa_mx_model()). Install it with install.packages('OpenMx') and ",
      "retry.",
      call. = FALSE
    )
  }
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

# PLAN 15 (D1-D4): derive the measurement inputs (fsL/fsT/fsb) from a
# get_fs() result's attributes for the all-args-NULL case. Dispatch:
#   - both unwrapped attributes plain matrices -> D2 fixed numeric;
#   - group_col attribute present              -> D7 stop (multigroup);
#   - per_obs/mirt_per_obs marker              -> D3 per-row: character
#     definition-variable matrices (+ int_fs_<score> columns appended to a
#     working copy of `data` when the fsb attribute is present);
#   - merMod 3-D fsL/fsT arrays (per-cluster) -> D3 per-row: the same
#     character definition-variable matrices (no fsb attribute -> the
#     default fixed-zero score intercepts);
#   - fs_pattern attribute present (SG FIML)   -> D3 per-pattern: the same
#     character matrices (pattern-constant columns), int_fs_<score> values
#     mapped row-wise through fs_pattern$label.
# Returns list(fsL, fsT, fsb, data, prov_err); fsL is NULL when nothing is
# derivable (prov_err = the D4 gate message, or NULL when the frame simply
# carries no fsT/fsL attributes) -- the caller turns that into the D5
# fail-fast. Stops() for multigroup input (D7) and for shapes it cannot
# handle.
tspa_mx_derive_measurement <- function(data) {
  attr_T <- attr(data, "fsT")
  attr_L <- attr(data, "fsL")
  if (is.null(attr_T) || is.null(attr_L)) {
    # No attributes at all: nothing to derive (D5 base message, caller).
    return(list(fsL = NULL, fsT = NULL, fsb = NULL, data = data,
                prov_err = NULL))
  }

  # D4 provenance gate: derive only from attributes that resolve as a
  # get_fs() result, reusing resolve_fs_per_row()'s informative errors
  # (PLAN 13 style, cf. R/tspa.R). A hand-rolled frame with plain matrix
  # attributes but no get_fs() provenance is NOT derived; its gate message
  # is carried to the caller's D5 fail-fast.
  prov <- tryCatch(
    {
      resolve_fs_per_row(data)
      TRUE
    },
    error = function(e) conditionMessage(e)
  )
  if (!isTRUE(prov)) {
    return(list(fsL = NULL, fsT = NULL, fsb = NULL, data = data,
                prov_err = prov))
  }

  # Unwrap only the unified single-group length-1 wrapper: a longer list
  # (per-pattern, per-row, multigroup) is dispatched below and must keep its
  # list form (tspa_mx_unwrap() would reject it with the multigroup message).
  L <- if (is.list(attr_L) && !is.matrix(attr_L) && length(attr_L) == 1L)
    attr_L[[1L]] else attr_L
  T <- if (is.list(attr_T) && !is.matrix(attr_T) && length(attr_T) == 1L)
    attr_T[[1L]] else attr_T
  # fsb: same length-1-only unwrap (a per-row / per-pattern list stays a
  # list); an all-NULL per-pattern list means "no intercepts" (fixed zero).
  b <- attr(data, "fsb")
  if (is.list(b) && !is.matrix(b)) {
    if (length(b) == 1L) {
      b <- b[[1L]]
    } else if (all(vapply(b, is.null, logical(1L)))) {
      b <- NULL
    }
  }

  if (is.matrix(L) && is.matrix(T)) {
    # D2: constant quantities -> fixed numeric matrices; fsb NULL -> the
    # default fixed-zero intercepts.
    return(list(fsL = L, fsT = T, fsb = b, data = data, prov_err = NULL))
  }

  # D7: the multigroup signal is the group_col attribute, NOT the list
  # length (a single-group FIML result also carries a list fsT, keyed by
  # pattern label -- the same disambiguation compute_fs_prod() uses).
  if (!is.null(attr(data, "group_col"))) {
    stop("Multigroup 'fsL'/'fsT' are not supported yet (Phase 1 is single-group).",
         call. = FALSE)
  }

  # Reference dimnames from the first row/pattern (structural, so valid even
  # for an all-NA first row): rownames = score names, colnames = latent names.
  # Note: from the already length-1-unwrapped L/T -- first_pattern_value()
  # itself only unwraps one list level, and a unified single-group attribute
  # is double-wrapped ("", then the per-row / per-pattern list).
  ref_L <- first_pattern_value(L)
  ref_T <- first_pattern_value(T)
  scores <- rownames(ref_L)
  latents <- colnames(ref_L)
  q <- length(scores)
  p <- length(latents)
  # D3: character definition-variable matrices referencing the frame's own
  # columns (the get_fs() naming, cf. fs_row_colnames()). The spec uses the
  # lower triangle of fsT, so the upper triangle stays NA.
  scores_T <- rownames(ref_T)
  L_char <- matrix(NA_character_, nrow = q, ncol = p,
                   dimnames = list(scores, latents))
  for (j in seq_len(p)) L_char[, j] <- paste0(latents[j], "_by_", scores)
  T_char <- matrix(NA_character_, nrow = q, ncol = q,
                   dimnames = list(scores_T, scores_T))
  for (i in seq_len(q)) {
    T_char[i, i] <- paste0("ev_", scores[i])
    for (j in seq_len(i - 1L)) {
      T_char[i, j] <- paste0("ecov_", scores[i], "_", scores[j])
    }
  }
  b_char <- setNames(paste0("int_", scores), scores)

  is_per_cluster <- is.array(L) && length(dim(L)) == 3L &&
    is.array(T) && length(dim(T)) == 3L
  if (is_per_cluster || isTRUE(attr(data, "per_obs")) ||
      isTRUE(attr(data, "mirt_per_obs"))) {
    # D3 per-row: fsb is a per-row list (one q-vector per row; the legacy
    # single-value case falls back to a shared constant, cf.
    # resolve_per_obs()). All-NA elements stay NA -> the D6 NA guard fires.
    # merMod per-cluster results (3-D fsL/fsT arrays, one row per cluster,
    # no fsb attribute -> fixed-zero score intercepts) take the same path.
    if (!is.null(b) && !is.list(b)) b <- rep(list(b), nrow(data))
    if (!is.null(b)) {
      for (i in seq_len(q)) {
        data[[paste0("int_", scores[i])]] <- vapply(
          b,
          function(v) if (is.null(v)) NA_real_ else as.numeric(v)[i],
          numeric(1L)
        )
      }
    }
    return(list(fsL = L_char, fsT = T_char,
                fsb = if (is.null(b)) NULL else b_char,
                data = data, prov_err = NULL))
  }

  fp <- attr(data, "fs_pattern")
  # Unified single-group results wrap fs_pattern in the same length-1 ""
  # list as the fsT/fsL attributes (cf. resolve_lavaan_unified()).
  if (is.list(fp) && length(fp) == 1L && !is.null(names(fp)) &&
      names(fp) == "") {
    fp <- fp[[1L]]
  }
  if (!is.list(fp) || !is.character(fp$label)) {
    stop(
      "Cannot derive the measurement inputs from 'data': its 'fsL'/'fsT' ",
      "attributes are lists, but the frame is neither a per-row result ",
      "(`per_obs`/`mirt_per_obs` marker) nor a pattern-keyed (SG FIML) ",
      "result with a character 'fs_pattern' label. Pass 'fsL'/'fsT' ",
      "explicitly.",
      call. = FALSE
    )
  }

  # D3 per-pattern (SG FIML joint): the columns are pattern-constant, so the
  # character matrices are complete. int_fs_<score> values map each row to
  # its pattern via fs_pattern$label (NA rows / absent patterns stay NA ->
  # the D6 NA guard fires).
  if (!is.null(b)) {
    labels <- fp$label
    for (i in seq_len(q)) {
      data[[paste0("int_", scores[i])]] <- vapply(
        labels,
        function(lb) {
          if (is.na(lb)) return(NA_real_)
          v <- b[[lb]]
          if (is.null(v)) NA_real_ else as.numeric(v)[i]
        },
        numeric(1L)
      )
    }
  }
  list(fsL = L_char, fsT = T_char,
       fsb = if (is.null(b)) NULL else b_char,
       data = data, prov_err = NULL)
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
  if (is.null(colnames(T)) || length(S) > 0L && !all(S %in% colnames(T))) {
    stop("'", arg, "' columns must be named by the factor-score names.",
         call. = FALSE)
  }
  T <- T[rownames = S, , drop = FALSE]
  T <- T[, colnames = S, drop = FALSE]
}

# se_fs definition variables: a character entry names a per-observation SE
# column of `data` (the square root of the measurement error). The RAM model
# needs the error VARIANCE, and an OpenMx definition variable is assigned the
# referenced column's value directly (no squaring op), so append a squared
# working column (cf. the int_fs_* columns the derivation route appends) and
# point the entry at it. Numeric entries (fixed SEs) pass through unchanged.
tspa_mx_resolve_se_fs <- function(se_fs, data) {
  if (is.null(se_fs) || length(se_fs) == 0L) {
    return(list(se_fs = se_fs, data = data))
  }
  if (is.numeric(se_fs)) {
    return(list(se_fs = se_fs, data = data))
  }
  if (!is.character(se_fs)) {
    stop("'se_fs' must be a named numeric vector (fixed SEs) or a named ",
         "character vector (definition-variable SE columns).", call. = FALSE)
  }
  # Sources are validated against the ORIGINAL data columns (not the working
  # copy, which grows as squared columns are appended), so an internally
  # generated working column can never be mistaken for a user source.
  orig_cols <- names(data)
  data2 <- data
  for (nm in names(se_fs)) {
    col <- unname(se_fs[nm])
    if (is.na(col) || !nzchar(col)) {
      stop("'se_fs' entry '", nm,
           "' must name a non-empty column of 'data'.", call. = FALSE)
    }
    if (!col %in% orig_cols) {
      # a numeric-looking entry that is not a column usually means a mixed
      # fixed/column se_fs (c() coerces the whole vector to character, so the
      # fixed value is now indistinguishable from a column name).
      if (!is.na(suppressWarnings(as.numeric(col)))) {
        stop("'se_fs' entry '", nm, "' has value '", col,
             "', which is not a column of 'data'. A numeric value is treated ",
             "as a fixed SE only when every entry is numeric; a mixed ",
             "fixed/column se_fs is not supported. Use all fixed SEs (numeric) ",
             "or all column names (character), or the fsL/fsT route.",
             call. = FALSE)
      }
      stop("'se_fs' definition-variable column '", col,
           "' is not a column of 'data'.", call. = FALSE)
    }
    if (!is.numeric(data[[col]])) {
      stop("'se_fs' definition-variable column '", col,
           "' must be numeric.", call. = FALSE)
    }
    # Finiteness + non-negativity (not just NA): Inf/NaN would corrupt the
    # squared definition variable, and a negative SE would square into a
    # positive (wrong) error variance.
    if (any(!is.finite(data[[col]]) | data[[col]] < 0)) {
      stop("'se_fs' definition-variable column '", col,
           "' must be finite and non-negative for every row.",
           call. = FALSE)
    }
    # A working column name unique against every existing and already-appended
    # column, so a real data column is never overwritten.
    se2col <- paste0(col, "_R2spa_se2")
    k <- 1L
    while (se2col %in% names(data2)) {
      se2col <- paste0(col, "_R2spa_se2_", k)
      k <- k + 1L
    }
    data2[[se2col]] <- data[[col]]^2
    se_fs[nm] <- se2col
  }
  list(se_fs = se_fs, data = data2)
}

tspa_mx_spec <- function(se_fs, fsL, fsT, fsb) {
  if (!is.null(se_fs)) {
    if (is.null(names(se_fs)) || anyNA(names(se_fs))) {
      stop("'se_fs' must be named by latent name.", call. = FALSE)
    }
    V <- names(se_fs)
    S <- paste0("fs_", V)
    # unit loadings (fixed) -- the single-score-per-latent se_fs convention
    Lm <- matrix(NA_real_, length(V), length(V), dimnames = list(S, V)); diag(Lm) <- 1
    L <- tspa_mx_cells(Lm, "fsL")
    # error variance: a numeric se_fs is a fixed SE (value se^2); a character
    # se_fs (post tspa_mx_resolve_se_fs) points at the per-row squared-SE
    # working columns. Built by hand so fixed and definition-variable cells
    # coexist in one table (tspa_mx_cells() enforces uniform type).
    Tv <- matrix(NA_real_, length(V), length(V), dimnames = list(S, S))
    Tc <- matrix(NA_character_, length(V), length(V), dimnames = list(S, S))
    if (is.numeric(se_fs)) {
      se <- as.numeric(se_fs)
      # Finite and non-negative: NA/Inf/NaN corrupt the fixed variance term,
      # and a negative SE squares into a positive (wrong) error variance.
      if (any(!is.finite(se) | se < 0)) {
        stop("'se_fs' standard errors must be finite and non-negative: ",
             "every latent needs a known factor-score SE.", call. = FALSE)
      }
      diag(Tv) <- se^2
    } else if (is.character(se_fs)) {
      diag(Tc) <- unname(se_fs)
    } else {
      stop("'se_fs' must be numeric (fixed SEs) or character (definition-",
           "variable SE columns).", call. = FALSE)
    }
    T <- list(vals = Tv, coln = Tc)
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
    if (anyNA(diag(fsT))) {
      stop("'fsT' must specify an error variance (diagonal entry) for every ",
           "factor score.", call. = FALSE)
    }
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
    # Absence, not zero: the spec's fixed-value slot must stay NA so
    # tspa_mx_model_string() emits no mean line (an explicit fsb = 0 keeps
    # its user-fixed c(0) row). With a structural mean model, lavaanify()
    # auto-seeds a fixed-zero score mean; otherwise need_mean adds a free one.
    b <- list(
      vals = matrix(NA_real_, 1L, length(S), dimnames = list(NULL, S)),
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
  score_used <- rep(FALSE, length(S))
  for (k in seq_along(V)) {
    terms <- character()
    for (i in seq_along(S)) {
      cv <- tspa_mx_cellval(spec$L$vals, spec$L$coln, i, k)
      if (!is.na(cv)) {
        terms <- c(terms, paste0(cv, " * ", S[i]))
        score_used[i] <- TRUE
      }
    }
    if (!length(terms)) {
      stop("Latent '", V[k], "' has no factor-score indicator in 'fsL'.",
           call. = FALSE)
    }
    lines <- c(lines, paste(V[k], "=~", paste(terms, collapse = " + ")))
  }
  if (any(!score_used)) {
    stop("Factor score(s) ", paste(S[!score_used], collapse = ", "),
         " load on no latent in 'fsL'.", call. = FALSE)
  }

  for (i in seq_along(S)) {
    for (j in seq_len(i)) {
      cv <- tspa_mx_cellval(spec$T$vals, spec$T$coln, i, j)
      if (!is.na(cv)) lines <- c(lines, paste0(S[i], " ~~ ", cv, " * ", S[j]))
    }
  }
  # Score means: an explicit fixed value (a "c(<b>)" row; user-fixed, so
  # the user == 0 auto-seed release in tspa_mx_paths() never re-frees it,
  # even at b = 0) or a per-row definition-variable column (the "c(1)"
  # sentinel, overlaid with the data label in tspa_mx_paths()). An absent
  # fsb (NULL) emits nothing: a structural mean model auto-seeds fixed-zero
  # score means; without one lav_to_mx_ram() adds free raw-data means.
  for (i in seq_along(S)) {
    cv <- tspa_mx_cellval(spec$b$vals, spec$b$coln, 1L, i)
    if (!is.na(cv)) lines <- c(lines, paste0(S[i], " ~ ", cv, " * 1"))
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
    if (is.na(si) || is.na(sj)) return(NA_character_)
    # lavaanify() may present a '~~' row with (lhs, rhs) in either orientation
    # relative to the score order, so a lower-triangle-only fsT (the documented
    # and the derived convention) must be found from both.
    v <- spec$T$coln[si, sj]
    if (is.na(v)) v <- spec$T$coln[sj, si]
    return(v)
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
      return(OpenMx::mxPath(from = mm$f, to = mm$t, arrows = mm$ar,
                            free = FALSE, values = 0,
                            labels = paste0("data.", dv)))
    }
    # lavaanify() auto-seeds structural-latent variances (and means, when a mean
    # structure is present) as user == 0, free == 0, ustart == 0; release them.
    # The user flag keeps a user-declared fixed zero (pathological, but legal
    # lavaan syntax) from being silently re-freed.
    if (r$user == 0L && r$free == 0L && is.numeric(r$ustart) && r$ustart == 0L &&
        r$lhs %in% spec$latents && (mm$ar == 2L && mm$f == mm$t || mm$f == "one")) {
      return(OpenMx::mxPath(from = mm$f, to = mm$t, arrows = mm$ar, free = TRUE,
                            values = if (mm$ar == 2L) 1 else 0))
    }
    free <- r$free > 0L
    args <- list(from = mm$f, to = mm$t, arrows = mm$ar, free = free,
                 values = if (free) tspa_mx_default_start(r$op, r$ustart)
                          else r$ustart)
    lab <- r$label
    if (is.character(lab) && !is.na(lab) && nzchar(lab)) args$labels <- lab
    do.call(OpenMx::mxPath, args)
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
      OpenMx::mxPath(from = "one", to = v, free = TRUE, values = 0)))
  }
  OpenMx::mxModel("m1", type = "RAM", manifestVars = manifests, latentVars = latents,
                  unlist(paths, recursive = FALSE),
                  OpenMx::mxData(observed = data, type = "raw"),
                  OpenMx::mxFitFunctionML())
}
