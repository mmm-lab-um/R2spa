# S3 methods for get_fs()
#
# Future methods (e.g. get_fs.mirt()) should be added to this file.

normalize_fs_method <- function(method) {
  method <- match.arg(
    method,
    c("regression", "Bartlett", "ML", "EB", "mean")
  )
  switch(method, ML = "Bartlett", EB = "regression", method)
}

validate_fs_priors <- function(prior_mean, prior_cov, lv_names) {
  q <- length(lv_names)

  if (!is.null(prior_mean)) {
    # Type check BEFORE as.numeric(): a factor would silently coerce to its
    # integer codes (accepted as a valid prior), a character vector coerces
    # with a warning to NA.
    if (!is.numeric(prior_mean) && !is.logical(prior_mean)) {
      stop("'prior_mean' must be numeric.", call. = FALSE)
    }
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
    # A non-numeric matrix (e.g. character, from a data.frame with factor or
    # character columns) would make the is.finite() check below error with a
    # cryptic base-R message instead of a clean stop.
    if (!is.numeric(prior_cov)) {
      stop("'prior_cov' must be numeric.", call. = FALSE)
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
    # chol() instead of eigen(): it is the standard positive-definiteness
    # test, is faster than a full eigendecomposition, and does not reject a
    # valid covariance whose smallest eigenvalue is within floating-point
    # noise of zero (an eigen() check with `<= 0` can).
    is_pd <- tryCatch({
      chol(pcv)
      TRUE
    }, error = function(e) FALSE)
    if (!is_pd) {
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
    if (q > 1L && (is.null(rn) || is.null(cn))) {
      stop(
        "'prior_cov' must be a named matrix (row and column names matching ",
        "the latent variable names: ", paste(lv_names, collapse = ", "),
        ") when the model has more than one latent variable; an unnamed ",
        "matrix's row and column order is ambiguous.",
        call. = FALSE
      )
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
  local = FALSE,
  method = c("regression", "Bartlett", "ML", "EB", "mean"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  format = c("unified", "list"),
  prior_mean = NULL,
  prior_cov = NULL,
  sum_items = NULL,
  product = NULL,
  ...
) {
  # `local` is a named formal placed before `...` so it is never forwarded
  # to cfa(). It only applies to explicitly supplied multi-statement models:
  # with `model = NULL` the auto single-factor model is trivially local and
  # the normal single-fit path runs unchanged (a no-op).
  model_given <- !is.null(model)
  if (is.null(model)) {
    ind_names <- colnames(object)
    if (!is.null(group)) {
      ind_names <- setdiff(ind_names, group)
    }
    if (length(ind_names) == 0L) {
      stop(
        "get_fs(): no indicator columns available (all columns are the ",
        "group variable); supply 'model' explicitly.",
        call. = FALSE
      )
    }
    model <- paste("f1 =~", paste(ind_names, collapse = " + "))
  }
  if (isTRUE(local) && model_given) {
    # v1 rejections (PLAN 14, D4): these quantities need cross-latent
    # information that separate local fits do not provide.
    if (vfsLT) {
      stop(
        "'vfsLT = TRUE' is not supported with 'local = TRUE': the latents ",
        "are scored from separate local fits, so the cross-latent sampling ",
        "covariances that 'vfsLT' requires do not exist (a block-diagonal ",
        "assembly of the per-local 'vfsLT' values would be wrong, not just ",
        "incomplete). 'tspa(corrected_se = TRUE)' is therefore not available ",
        "from a local stage 1.",
        call. = FALSE
      )
    }
    if (!is.null(prior_cov)) {
      stop(
        "'prior_cov' is not supported with 'local = TRUE' (v1): a q x q ",
        "prior covariance cannot be reduced to the per-latent 1 x 1 priors ",
        "used by the separate local fits without silently dropping its ",
        "off-diagonals.",
        call. = FALSE
      )
    }
    if (reliability) {
      stop(
        "'reliability = TRUE' is not supported with 'local = TRUE' (v1): ",
        "the per-latent 'reliability' attribute would introduce a new ",
        "attribute shape for an attribute 'tspa()' already deprecates.",
        call. = FALSE
      )
    }
    if (!is.null(product)) {
      stop("'product' is not supported with 'local = TRUE' (v1).",
           call. = FALSE)
    }
    return(get_fs_local(
      object,
      model = model,
      group = group,
      method = method,
      corrected_fsT = corrected_fsT,
      format = format,
      prior_mean = prior_mean,
      sum_items = sum_items,
      ...
    ))
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
    prior_cov = prior_cov,
    sum_items = sum_items,
    product = product
  )
}

# ---------------------------------------------------------------------------
# Per-construct ("local") stage-1 scoring (PLAN 14): `get_fs(..., local = TRUE)`.
#
# Each latent is scored from its own local measurement model (the canonical
# two-stage path analysis setup) instead of one joint multi-factor model.
# The merged output reproduces the joint layout (columns and attribute
# shapes) exactly, with exactly-zero cross terms (off-diagonal `_by_`
# columns, `ecov_*` columns, and `fsT`/`fsL`/`psi` off-diagonals) encoding
# "no shared measurement model" (D2/D6).
# ---------------------------------------------------------------------------

# Strict-grammar parser for the local-mode string form: splits a measurement
# model string into per-latent model strings. Only statements of the form
# `<latent> =~ <item1> + <item2> + ...` (bare identifiers) are accepted;
# `#` comments are stripped, `;` separates statements, and a line ending in
# `+` continues onto the next line. Anything else is an error naming the
# offending line and pointing to the alternatives (joint mode, `local = FALSE`,
# or the vector/list form). Returns the per-latent model strings in statement
# order, named by latent.
split_local_models <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop(
      "get_fs(local = TRUE): 'model' must be a single string, a character ",
      "vector, or a list of strings.",
      call. = FALSE
    )
  }
  raw_lines <- strsplit(model, "\n", fixed = TRUE)[[1L]]
  # Per physical line: strip the `#` comment, split on `;`, and record
  # whether the (trimmed, comment-free) line ends in `+` (continuation).
  line_segs <- vector("list", length(raw_lines))
  for (i in seq_along(raw_lines)) {
    lc <- sub("#.*$", "", raw_lines[i])
    segs <- trimws(strsplit(lc, ";", fixed = TRUE)[[1L]])
    line_segs[[i]] <- list(
      segs = segs[segs != ""],
      cont = endsWith(trimws(lc), "+")
    )
  }
  # Assemble statements: a segment starts a statement; a `;` or line end
  # closes it; a trailing `+` at a line end keeps it open for the next line
  # (a duplicated leading `+` on the continuation line is dropped).
  stmts <- list()
  cur <- NULL
  for (i in seq_along(raw_lines)) {
    l <- line_segs[[i]]
    for (s in l$segs) {
      if (is.null(cur)) {
        cur <- list(txt = s, line = i)
      } else if (endsWith(cur$txt, "+") && cur$line < i) {
        cur$txt <- if (startsWith(s, "+")) {
          paste0(cur$txt, " ", substr(s, 2L, nchar(s)))
        } else {
          paste0(cur$txt, " ", s)
        }
      } else {
        # A `;` (or a same-line dangling `+`) closes the open statement.
        stmts[[length(stmts) + 1L]] <- cur
        cur <- list(txt = s, line = i)
      }
    }
    if (!is.null(cur) && !(l$cont && endsWith(cur$txt, "+"))) {
      stmts[[length(stmts) + 1L]] <- cur
      cur <- NULL
    }
  }
  if (!is.null(cur)) {
    if (endsWith(cur$txt, "+")) {
      local_model_syntax_error(
        cur$line,
        cur$txt,
        "a trailing '+' with no items on the following line"
      )
    }
    stmts[[length(stmts) + 1L]] <- cur
  }
  if (length(stmts) == 0L) {
    stop(
      "get_fs(local = TRUE): 'model' contains no statements of the form ",
      "'<latent> =~ <items>'.",
      call. = FALSE
    )
  }
  parsed <- lapply(seq_along(stmts), function(i) {
    parse_local_statement(stmts[[i]]$txt, stmts[[i]]$line)
  })
  latents <- vapply(parsed, `[[`, character(1L), "lhs")
  dup_lv <- latents[duplicated(latents)]
  if (length(dup_lv) > 0L) {
    stop(
      "get_fs(local = TRUE): the latent variable '", dup_lv[1L], "' is ",
      "defined more than once; every latent must be defined by exactly one ",
      "statement.",
      call. = FALSE
    )
  }
  items <- unlist(lapply(parsed, `[[`, "rhs"), use.names = FALSE)
  dup_items <- items[duplicated(items)]
  if (length(dup_items) > 0L) {
    stop(
      "get_fs(local = TRUE): item(s) ",
      paste(dup_items, collapse = ", "),
      " are used in more than one statement; each item may load on exactly ",
      "one latent.",
      call. = FALSE
    )
  }
  setNames(
    vapply(
      parsed,
      function(st) paste0(st$lhs, " =~ ", paste(st$rhs, collapse = " + ")),
      character(1L)
    ),
    latents
  )
}

# Parse one local-mode statement (`lhs =~ i1 + i2 + ...`), erroring (via
# local_model_syntax_error()) with a class-specific detail on anything else.
parse_local_statement <- function(txt, line) {
  if (grepl("~~", txt, fixed = TRUE)) {
    local_model_syntax_error(
      line,
      txt,
      "'~~' (latent-latent covariance/correlation, or residual covariance ",
      "between indicators) is not supported in the local-mode string form; ",
      "the vector form accepts it when it stays within one factor's items"
    )
  }
  if (grepl("|~", txt, fixed = TRUE)) {
    local_model_syntax_error(
      line,
      txt,
      "an ordered-response statement ('|~') is not supported in local mode"
    )
  }
  if (grepl("$", txt, fixed = TRUE)) {
    local_model_syntax_error(
      line,
      txt,
      "thresholds ('$') are not supported in local mode"
    )
  }
  if (!grepl("=~", txt, fixed = TRUE)) {
    if (grepl("~", txt, fixed = TRUE)) {
      local_model_syntax_error(
        line,
        txt,
        "a structural path ('~') is not supported in local mode (measurement ",
        "statements only)"
      )
    }
    local_model_syntax_error(line, txt, "no '=~' operator")
  }
  n_op <- (nchar(txt) - nchar(gsub("=~", "", txt, fixed = TRUE))) / 2L
  if (n_op > 1L) {
    local_model_syntax_error(line, txt, "more than one '=~' operator")
  }
  parts <- strsplit(txt, "=~", fixed = TRUE)[[1L]]
  lhs <- trimws(parts[1L])
  rhs <- trimws(parts[2L])
  if (grepl("\\s", lhs)) {
    local_model_syntax_error(
      line,
      txt,
      "multiple latent names on the left-hand side (only one latent per ",
      "statement)"
    )
  }
  if (!grepl("^[A-Za-z._][A-Za-z0-9._]*$", lhs)) {
    local_model_syntax_error(
      line,
      txt,
      "the left-hand side must be a single bare identifier (labels and fixed ",
      "values are not allowed)"
    )
  }
  # strsplit() drops a trailing empty field, so "a =~" yields a one-element
  # split and parts[2L] is NA (not ""); treat missing/NA/zero-length RHS as
  # the empty-RHS case rather than letting `rhs == ""` hit a missing value.
  if (length(rhs) == 0L || is.na(rhs) || rhs == "") {
    local_model_syntax_error(line, txt, "an empty right-hand side")
  }
  if (grepl("~", rhs, fixed = TRUE)) {
    local_model_syntax_error(
      line,
      txt,
      "'~' is not allowed on the right-hand side"
    )
  }
  if (grepl("*", rhs, fixed = TRUE)) {
    local_model_syntax_error(
      line,
      txt,
      "labels and fixed values ('*') are not allowed on the right-hand side"
    )
  }
  if (grepl("(", rhs, fixed = TRUE) || grepl(")", rhs, fixed = TRUE)) {
    local_model_syntax_error(line, txt, "'c(...)' calls are not allowed")
  }
  # strsplit() drops a trailing empty field, so a dangling '+' would be
  # silently ignored (fitted as if it were not there); reject it explicitly.
  if (endsWith(rhs, "+")) {
    local_model_syntax_error(
      line,
      txt,
      "a dangling '+' operator at the end of the right-hand side"
    )
  }
  toks <- trimws(strsplit(rhs, "+", fixed = TRUE)[[1L]])
  bad <- toks[!grepl("^[A-Za-z._][A-Za-z0-9._]*$", toks)]
  if (length(bad) > 0L) {
    local_model_syntax_error(
      line,
      txt,
      "right-hand side items must be bare identifiers (got: ",
      paste(bad, collapse = ", "), ")"
    )
  }
  list(lhs = lhs, rhs = toks, line = line)
}

# The shared local-mode string-form error: names the offending line and
# points to the alternatives (joint mode or the vector/list form). The call
# sites pass a possibly multi-line detail as separate string arguments (a
# single sentence wrapped across source lines); `...` reassembles them
# verbatim (collapse = "", the fragments carry their own spacing) so the
# message keeps its exact wording.
local_model_syntax_error <- function(line, txt, detail, ...) {
  detail <- paste0(c(detail, ...), collapse = "")
  stop(
    "get_fs(local = TRUE): unsupported model syntax on line ", line,
    " ('", txt, "'): ", detail,
    ". The local-mode string form only accepts statements of the form ",
    "'<latent> =~ <item1> + <item2> + ...'. Fit the joint model instead ",
    "(local = FALSE) or pass one complete single-factor model string per ",
    "latent as a character vector or list.",
    call. = FALSE
  )
}

# Orchestrate local per-construct scoring: resolve the per-latent models,
# fit each with cfa() (all `...` forwarded), score each with the existing
# get_fs.lavaan() path (format = "list" internally), and merge.
get_fs_local <- function(
  object,
  model,
  group = NULL,
  method = c("regression", "Bartlett", "ML", "EB", "mean"),
  corrected_fsT = FALSE,
  format = c("unified", "list"),
  prior_mean = NULL,
  sum_items = NULL,
  ...
) {
  # -- Resolve the per-latent model strings (named in latent order). --
  if (is.list(model)) {
    if (!all(vapply(
      model,
      function(e) is.character(e) && length(e) >= 1L,
      logical(1L)
    ))) {
      stop(
        "get_fs(local = TRUE): the list form of 'model' must contain ",
        "character strings only (one complete single-factor model per ",
        "element).",
        call. = FALSE
      )
    }
    models <- setNames(
      as.list(unlist(lapply(model, paste, collapse = "\n"))),
      names(model)
    )
  } else if (length(model) > 1L) {
    models <- as.list(model)
  } else {
    models <- as.list(split_local_models(model))
  }
  q <- length(models)

  # -- Fit each latent's local model (a vector/list element is fit verbatim;
  #    a parsed statement is a canonical `<latent> =~ <items>` string). --
  fits <- lapply(models, function(mk) cfa(mk, data = object, group = group, ...))

  # -- Validate: each fit defines exactly one latent; names are unique. --
  ests <- lapply(fits, function(fk) lavInspect(fk, what = "est"))
  est1s <- lapply(ests, function(e) if (!is.null(e$psi)) e else e[[1L]])
  n_lv <- vapply(est1s, function(e) nrow(e$psi), integer(1L))
  bad_k <- which(n_lv != 1L)
  if (length(bad_k) > 0L) {
    stop(
      "get_fs(local = TRUE): model element(s) ",
      paste(bad_k, collapse = ", "),
      " define ",
      paste(n_lv[bad_k], collapse = " / "),
      " latent variable(s); every element must define exactly one latent.",
      call. = FALSE
    )
  }
  lv_names <- vapply(
    est1s,
    function(e) colnames(e$lambda)[1L],
    character(1L)
  )
  dup_lv <- lv_names[duplicated(lv_names)]
  if (length(dup_lv) > 0L) {
    stop(
      "get_fs(local = TRUE): latent variable name(s) ",
      paste(dup_lv, collapse = ", "),
      " occur in more than one model element; latent names must be unique.",
      call. = FALSE
    )
  }
  if (!is.null(names(models)) && !all(names(models) == lv_names)) {
    mismatch <- which(names(models) != lv_names)
    stop(
      "get_fs(local = TRUE): model element(s) ",
      paste(mismatch, collapse = ", "),
      " are named ",
      paste(names(models)[mismatch], collapse = " / "),
      " but their fitted latent(s) are ",
      paste(lv_names[mismatch], collapse = " / "), ".",
      call. = FALSE
    )
  }

  # -- prior_mean: validated once against all latent names, then sliced. --
  pm_all <- if (is.null(prior_mean)) {
    NULL
  } else {
    validate_fs_priors(prior_mean, NULL, lv_names)$mean
  }

  # -- sum_items: user-supplied list sliced per latent (auto-derivation is
  #    trivially satisfied per single-factor local model). --
  if (!is.null(sum_items)) {
    if (!is.list(sum_items) || is.null(names(sum_items))) {
      stop(
        "'sum_items' must be a named list mapping factor names to item names.",
        call. = FALSE
      )
    }
    unknown_lv <- setdiff(names(sum_items), lv_names)
    if (length(unknown_lv) > 0L) {
      stop(
        "Unknown factor name(s) in 'sum_items': ",
        paste(deparse(unknown_lv), collapse = ", "),
        ". Local model factors are: ",
        paste(deparse(lv_names), collapse = ", "), ".",
        call. = FALSE
      )
    }
    missing_lv <- setdiff(lv_names, names(sum_items))
    if (length(missing_lv) > 0L) {
      stop(
        "'sum_items' must cover all model factors; no items given for: ",
        paste(deparse(missing_lv), collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # -- Score each local fit through the existing (unchanged) method. --
  fs_list <- lapply(seq_len(q), function(k) {
    get_fs(
      fits[[k]],
      method = method,
      corrected_fsT = corrected_fsT,
      format = "list",
      prior_mean = if (is.null(pm_all)) NULL else pm_all[lv_names[k]],
      sum_items = if (is.null(sum_items)) {
        NULL
      } else {
        setNames(list(sum_items[[lv_names[k]]]), lv_names[k])
      }
    )
  })
  names(fs_list) <- lv_names

  merge_local_fs(
    fs_list,
    latent_names = lv_names,
    format = format,
    group_col = group
  )
}

# Merge per-latent get_fs() results (each format = "list") into the joint
# multi-factor layout. Every row of the data carries one block-diagonal
# (fsL, fsT, fsb) tuple across the latents (an all-NA 1 x 1 block where a
# latent could not score the row); when every row of every group carries an
# identical tuple (complete data) the compact per-group matrix/vector
# attributes are materialized (via assemble_fs_blocks()); otherwise (FIML)
# the flat per-row attribute lists and the `per_obs` marker are attached
# (the mirt per-obs convention). Column order, attribute shapes, and the
# group-column/group_col handling are exactly the joint ones.
merge_local_fs <- function(
  fs_list,
  latent_names,
  format = c("unified", "list"),
  group_col = NULL
) {
  format <- match.arg(format)
  q <- length(fs_list)
  # Positional: callers may pass a named vector (get_fs_local()'s lv_names
  # carries the latent names as element names); dimnames are value-based.
  latent_names <- unname(latent_names)
  fs_names <- paste0("fs_", latent_names)

  # -- Group structure (shared by all per-latent results: same data/group).
  first <- fs_list[[1L]]
  mg <- is.list(first) && !is.data.frame(first)
  group_labels <- if (mg) names(first) else ""
  ngroups <- length(group_labels)
  df_kg <- function(k, g) if (mg) fs_list[[k]][[g]] else fs_list[[k]]

  n_g <- if (mg) vapply(first, nrow, integer(1L)) else nrow(first)
  for (k in seq_len(q)) {
    f <- fs_list[[k]]
    if (mg) {
      if (!identical(names(f), group_labels) ||
          !identical(vapply(f, nrow, integer(1L)), n_g)) {
        stop(
          "get_fs(local = TRUE): the per-latent fits have different group ",
          "labels or per-group row counts; the per-row merge needs every ",
          "local fit to keep the same rows in every group. This happens with ",
          "listwise missing-data handling (the cfa() default), which drops ",
          "different rows in each local fit. Re-run with a FIML method ",
          "(e.g. missing = 'fiml') so all local fits keep every row.",
          call. = FALSE
        )
      }
    } else if (nrow(f) != n_g) {
      stop(
        "get_fs(local = TRUE): the per-latent fits have different row counts ",
        "(", paste(vapply(fs_list, nrow, integer(1L)), collapse = ", "),
        "); the per-row merge needs every local fit to keep the same rows. ",
        "This happens with listwise missing-data handling (the cfa() ",
        "default), which drops different rows in each local fit. Re-run with ",
        "a FIML method (e.g. missing = 'fiml') so all local fits keep every ",
        "row.",
        call. = FALSE
      )
    }
  }
  n <- if (mg) sum(n_g) else n_g

  # -- Resolve each per-latent result to per-row (fsL, fsT, fsb) blocks. --
  res <- lapply(fs_list, resolve_fs_per_row)
  pidx <- lapply(res, `[[`, "pattern_idx")
  blks <- lapply(res, `[[`, "blocks")
  group_vals <- if (mg) res[[1L]]$group_vals else NULL
  if (mg) {
    for (k in seq_len(q)) {
      if (!identical(res[[k]]$group_vals, group_vals)) {
        stop(
          "internal error: per-latent group values differ (local merge).",
          call. = FALSE
        )
      }
    }
  }
  rows_of_group <- if (mg) {
    lapply(group_labels, function(g) which(group_vals == g))
  } else {
    list(seq_len(n))
  }

  # -- Item order: first appearance across the local models (per-latent
  #    indicator names from the per-group pattern attributes). --
  items_by_latent <- vector("list", q)
  for (k in seq_len(q)) {
    fp <- attr(df_kg(k, 1L), "fs_pattern")
    if (is.null(fp) || is.null(fp$pat) || is.null(rownames(fp$pat))) {
      stop(
        "internal error: local merge requires the per-group 'fs_pattern' ",
        "attributes with indicator names.",
        call. = FALSE
      )
    }
    items_by_latent[[k]] <- rownames(fp$pat)
  }
  # unique() keeps first occurrences, i.e. first appearance across the
  # per-latent vectors in latent order (the old setdiff accumulation).
  item_names <- unique(unlist(items_by_latent, use.names = FALSE))
  p_all <- length(item_names)
  item_pos <- lapply(items_by_latent, match, item_names)

  # -- Per-row scores (one column per latent, in latent order). --
  score_mx <- matrix(NA_real_, nrow = n, ncol = q)
  for (k in seq_len(q)) {
    score_mx[, k] <- unname(as.matrix(res[[k]]$scores))
  }
  colnames(score_mx) <- latent_names

  # -- Per-row merged (block-diagonal) blocks. Each local fit defines exactly
  #    one latent (enforced in get_fs_local()), so every per-latent (fsL, fsT,
  #    fsb) block is 1 x 1 and the merged q x q matrices are zero off-diagonal
  #    with per-pattern scalars on the diagonal. A row's merged blocks and
  #    scoring-matrix fragments are determined entirely by its (group,
  #    per-latent pattern) tuple, so one template per distinct tuple is minted
  #    and the reference is shared by every row carrying that tuple (rows keep
  #    distinct list slots; every downstream consumer reads the blocks
  #    read-only). Measured at n = 50k, q = 2, ~50% missing per factor: the
  #    FIML merge drops from ~4.9 s to ~0.05 s, and the compact path (which
  #    reads a single row per group) from ~4.9 s to ~0.15 s.
  g_of_row <- if (mg) match(group_vals, group_labels) else rep(1L, n)
  # Per-latent per-block 1 x 1 scalars (fsL / fsT / fsb) + the unscorable
  # marker (all-NA fsT: the latent could not score the block's rows).
  Lsc <- vector("list", q); Tsc <- vector("list", q)
  Bsc <- vector("list", q); unsc <- vector("list", q)
  for (k in seq_len(q)) {
    Lsc[[k]] <- vapply(blks[[k]], function(b) as.numeric(b$fsL)[1L],
                       numeric(1L))
    Tsc[[k]] <- vapply(blks[[k]], function(b) as.numeric(b$fsT)[1L],
                       numeric(1L))
    Bsc[[k]] <- vapply(blks[[k]], function(b) {
      bb <- b$fsb
      if (is.null(bb)) NA_real_ else as.numeric(bb)[1L]
    }, numeric(1L))
    unsc[[k]] <- vapply(blks[[k]], function(b) all(is.na(b$fsT)), logical(1L))
  }
  diagL_mx <- matrix(NA_real_, n, q)
  diagT_mx <- matrix(NA_real_, n, q)
  B_mx <- matrix(NA_real_, n, q)
  for (k in seq_len(q)) {
    pk_r <- pidx[[k]]
    diagL_mx[, k] <- Lsc[[k]][pk_r]
    diagT_mx[, k] <- Tsc[[k]][pk_r]
    B_mx[, k] <- Bsc[[k]][pk_r]
  }
  # Distinct (group, per-latent pattern) tuples: the key is a string join of
  # the integer fields (non-negative, so the separator is unambiguous).
  # unique()/match() give the dense 1..ntpl encoding in first-appearance
  # order; tpl_row[u] is the first-occurrence row of distinct key u, the
  # template representative.
  key_v <- as.character(g_of_row)
  for (k in seq_len(q)) {
    key_v <- paste0(key_v, "\r", pidx[[k]])
  }
  key_u <- unique(key_v)
  tuple <- match(key_v, key_u)
  ntpl <- length(key_u)
  tpl_row <- match(key_u, key_v)

  # -- Per-row scoring matrices (row k = local fit k's a_k, zero-padded to
  #    all items; all-NA row where that latent could not score the row;
  #    colnames NULL per the joint-path dimname quirk). Precomputed per
  #    (latent, group, pattern column): each row's pattern column plus that
  #    column's observed/missing item positions and the scoring values.
  sm_kg <- vector("list", q * ngroups)
  lab_kg <- vector("list", q * ngroups)
  pat_kg <- vector("list", q * ngroups)
  for (k in seq_len(q)) {
    for (g in seq_len(ngroups)) {
      d <- df_kg(k, g)
      sm_kg[[(k - 1L) * ngroups + g]] <- attr(d, "scoring_matrix")
      fp <- attr(d, "fs_pattern")
      lab_kg[[(k - 1L) * ngroups + g]] <- fp$label
      pat_kg[[(k - 1L) * ngroups + g]] <- fp$pat
    }
  }
  col_of_row <- vector("list", q)
  pos_obs <- vector("list", q)
  pos_miss <- vector("list", q)
  sm_vals <- vector("list", q)
  for (k in seq_len(q)) {
    cor_k <- integer(n)
    po_k <- vector("list", ngroups)
    pm_k <- vector("list", ngroups)
    sv_k <- vector("list", ngroups)
    for (g in seq_len(ngroups)) {
      idx_kg <- (k - 1L) * ngroups + g
      sm <- sm_kg[[idx_kg]]
      if (is.matrix(sm)) {
        # Single-pattern (complete) group: one virtual pattern column.
        po_k[[g]] <- list(item_pos[[k]])
        pm_k[[g]] <- list(integer(0))
        sv_k[[g]] <- list(as.numeric(sm))
        cor_k[rows_of_group[[g]]] <- 1L
      } else {
        pk <- pat_kg[[idx_kg]]
        cols <- colnames(pk)
        lab2col <- seq_along(cols)
        names(lab2col) <- cols
        cor_k[rows_of_group[[g]]] <- lab2col[lab_kg[[idx_kg]]]
        po_g <- vector("list", ncol(pk))
        pm_g <- vector("list", ncol(pk))
        sv_g <- vector("list", ncol(pk))
        for (cc in seq_len(ncol(pk))) {
          obs_c <- rownames(pk)[pk[, cc, drop = FALSE]]
          po_g[[cc]] <- match(obs_c, item_names)
          pm_g[[cc]] <- match(setdiff(items_by_latent[[k]], obs_c), item_names)
          sv_g[[cc]] <- as.numeric(sm[[cols[cc]]])
        }
        po_k[[g]] <- po_g
        pm_k[[g]] <- pm_g
        sv_k[[g]] <- sv_g
      }
    }
    col_of_row[[k]] <- cor_k
    pos_obs[[k]] <- po_k
    pos_miss[[k]] <- pm_k
    sm_vals[[k]] <- sv_k
  }
  # One template per distinct tuple, minted at its first-occurrence row.
  diag_q <- cbind(seq_len(q), seq_len(q))
  dn_L <- list(fs_names, latent_names)
  dn_T <- list(fs_names, fs_names)
  L_tmpl <- vector("list", ntpl)
  T_tmpl <- vector("list", ntpl)
  B_tmpl <- vector("list", ntpl)
  S_tmpl <- vector("list", ntpl)
  for (u in seq_len(ntpl)) {
    r <- tpl_row[u]
    g_r <- g_of_row[r]
    L_u <- matrix(0, q, q, dimnames = dn_L)
    L_u[diag_q] <- diagL_mx[r, ]
    T_u <- matrix(0, q, q, dimnames = dn_T)
    T_u[diag_q] <- diagT_mx[r, ]
    S_u <- matrix(0, q, p_all, dimnames = list(latent_names, NULL))
    for (k in seq_len(q)) {
      if (unsc[[k]][pidx[[k]][r]]) {
        # Unscorable latent (all-NA block): NA across its own items.
        S_u[k, item_pos[[k]]] <- NA_real_
        next
      }
      c_r <- col_of_row[[k]][r]
      S_u[k, pos_obs[[k]][[g_r]][[c_r]]] <- sm_vals[[k]][[g_r]][[c_r]]
      pm_c <- pos_miss[[k]][[g_r]][[c_r]]
      if (length(pm_c) > 0L) {
        S_u[k, pm_c] <- NA_real_
      }
    }
    L_tmpl[[u]] <- L_u
    T_tmpl[[u]] <- T_u
    B_tmpl[[u]] <- setNames(B_mx[r, ], fs_names)
    S_tmpl[[u]] <- S_u
  }
  L_row_list <- L_tmpl[tuple]
  T_row_list <- T_tmpl[tuple]
  B_row_list <- B_tmpl[tuple]
  S_row_list <- S_tmpl[tuple]

  # -- Group-level latent moments (block-diagonal psi, concatenated alpha).
  psi_g <- vector("list", ngroups)
  alpha_g <- vector("list", ngroups)
  names(psi_g) <- names(alpha_g) <- group_labels
  for (g in seq_len(ngroups)) {
    psi_g[[g]] <- block_diag(lapply(seq_len(q), function(k) {
      as.matrix(attr(df_kg(k, g), "psi"))
    }))
    alpha_g[[g]] <- do.call(
      c,
      lapply(seq_len(q), function(k) as.numeric(attr(df_kg(k, g), "alpha")))
    )
  }

  # -- Compact (complete-data) detection: every row of a group carries an
  #    identical per-latent pattern block (across all latents). --
  group_compact <- vapply(seq_len(ngroups), function(g) {
    rows_g <- rows_of_group[[g]]
    all(vapply(seq_len(q), function(k) {
      length(unique(pidx[[k]][rows_g])) == 1L
    }, logical(1L)))
  }, logical(1L))
  per_row_mode <- !all(group_compact)

  if (!per_row_mode) {
    # Compact form: one block per group through the existing assembler
    # (columns, pattern bookkeeping, and group-column handling included).
    blocks_by_group <- setNames(vector("list", ngroups), group_labels)
    for (g in seq_len(ngroups)) {
      rows_g <- rows_of_group[[g]]
      r1 <- rows_g[1L]
      L_g <- L_row_list[[r1]]
      T_g <- T_row_list[[r1]]
      # Per-latent observed-indicator patterns (a single pattern per
      # compact group), merged over all items in first-appearance order.
      pat_merged <- rep(FALSE, p_all)
      for (k in seq_len(q)) {
        fp <- attr(df_kg(k, g), "fs_pattern")
        pat_merged[item_pos[[k]]] <- as.logical(fp$pat[, 1L, drop = TRUE])
      }
      names(pat_merged) <- item_names
      label_g <- paste(item_names[pat_merged], collapse = "+")
      S_g <- matrix(0, nrow = q, ncol = p_all, dimnames = list(latent_names, NULL))
      for (k in seq_len(q)) {
        sm <- attr(df_kg(k, g), "scoring_matrix")
        S_g[k, item_pos[[k]]] <- if (is.matrix(sm)) {
          as.numeric(sm)
        } else {
          as.numeric(sm[[1L]])
        }
      }
      fs_g <- as.data.frame(score_mx[rows_g, , drop = FALSE])
      attr(fs_g, "fsL") <- L_g
      blocks_by_group[[g]] <- list(
        list(
          case_idx = seq_len(length(rows_g)),
          fs = fs_g,
          fsT = T_g,
          fsL = L_g,
          fsb = B_row_list[[r1]],
          scoring_matrix = S_g,
          pat_label = label_g,
          pat = pat_merged
        )
      )
    }
    out <- assemble_fs_blocks(
      blocks_by_group,
      format = format,
      group_col = if (mg) group_col else NULL
    )
  } else {
    # Per-row (FIML) form: flat per-row attribute lists + the `per_obs`
    # marker (the mirt per-obs convention). The output is always a single
    # data frame (a trailing group column + group_col attribute for MG).
    k3 <- q + q * q + q * (q + 1L) / 2L
    # One fs_row_cols() call per distinct tuple: rows sharing a tuple share
    # the (fsL, fsT, fsb) template, hence the same se/loading/error row.
    # fs_row_cols() only reads nrow(fs) (always 1 here): one shared 1-row
    # token.
    dummy_fs <- data.frame(x = 0L)
    vals_u <- matrix(NA_real_, ntpl, k3)
    for (u in seq_len(ntpl)) {
      vals_u[u, ] <- fs_row_cols(
        dummy_fs, L_tmpl[[u]], T_tmpl[[u]], B_tmpl[[u]]
      )[1L, seq_len(k3), drop = FALSE]
    }
    vals <- vals_u[tuple, , drop = FALSE]
    se_nm <- paste0(fs_names, "_se")
    ld_nm <- c(create_fsL_names(latent_names, fs_names))
    ev_nm <- character(q * (q + 1L) / 2L)
    count <- 1L
    for (i in seq_len(q)) {
      for (j in seq_len(i)) {
        ev_nm[count] <- if (i == j) {
          paste0("ev_", fs_names[i])
        } else {
          paste0("ecov_", fs_names[i], "_", fs_names[j])
        }
        count <- count + 1L
      }
    }
    out <- as.data.frame(
      cbind(as.matrix(score_mx), vals),
      check.names = FALSE
    )
    colnames(out) <- c(fs_names, se_nm, ld_nm, ev_nm)
    rownames(out) <- NULL
    if (mg) {
      out[[group_col]] <- group_vals
      attr(out, "group_col") <- group_col
    }
    attr(out, "fsT") <- T_row_list
    attr(out, "fsL") <- L_row_list
    attr(out, "fsb") <- B_row_list
    attr(out, "scoring_matrix") <- S_row_list
    if (mg) {
      attr(out, "fs_pattern") <- setNames(
        lapply(n_g, function(m) list(label = seq_len(m), pat = NULL)),
        group_labels
      )
    } else {
      attr(out, "fs_pattern") <- list(label = seq_len(n), pat = NULL)
    }
    attr(out, "per_obs") <- TRUE
  }

  # Group-level latent moments mirror the fsT shape (see get_fs.lavaan());
  # per-row (FIML) results are always the single-data-frame ("unified")
  # shape, so their psi/alpha are list-valued (SG: length-1 list named "").
  if (format == "unified" || per_row_mode) {
    attr(out, "psi") <- psi_g
    attr(out, "alpha") <- alpha_g
  } else if (ngroups == 1L) {
    attr(out, "psi") <- psi_g[[1L]]
    attr(out, "alpha") <- alpha_g[[1L]]
  } else {
    for (g in seq_len(ngroups)) {
      attr(out[[g]], "psi") <- psi_g[[g]]
      attr(out[[g]], "alpha") <- alpha_g[[g]]
    }
    attr(out, "psi") <- psi_g
    attr(out, "alpha") <- alpha_g
  }
  out
}

get_fs_blocks.lavaan <- function(
  object,
  method,
  add_to_evfs,
  prior_mean = NULL,
  prior_cov = NULL,
  sum_items = NULL,
  ...
) {
  method <- match.arg(method, c("regression", "Bartlett", "mean"))
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
          fs_matrices = TRUE,
          sum_items = sum_items
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
            fs_matrices = TRUE,
            sum_items = sum_items
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
  method = c("regression", "Bartlett", "ML", "EB", "mean"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  reliability = FALSE,
  format = c("unified", "list"),
  prior_mean = NULL,
  prior_cov = NULL,
  sum_items = NULL,
  product = NULL,
  ...
) {
  method <- normalize_fs_method(method)
  if (!inherits(object, "lavaan")) {
    stop("`object` must be a `lavaan` model object.")
  }
  if (method == "mean") {
    # "mean" scores are raw item means: they bypass the corrected-FS-T /
    # vfsLT / reliability / prior machinery, which only supports
    # regression and Bartlett scores.
    # A FIML-style fit carries @Data@Mp even for complete data (a single
    # all-TRUE pattern), so missingness is read from the patterns, not from
    # the slot's presence.
    if (!all(vapply(object@Data@Mp, function(mp) {
      is.null(mp) || all(mp$pat)
    }, logical(1)))) {
      stop(
        "method = 'mean' does not support models fitted with missing ",
        "data; use method = 'regression' or 'Bartlett'.",
        call. = FALSE
      )
    }
    incompatible <- c(
      if (corrected_fsT) "corrected_fsT",
      if (vfsLT) "vfsLT",
      if (reliability) "reliability",
      if (!is.null(prior_mean)) "prior_mean",
      if (!is.null(prior_cov)) "prior_cov"
    )
    if (length(incompatible) > 0) {
      stop(
        "method = 'mean' is not supported together with: ",
        paste(incompatible, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  format <- match.arg(format)

  # The sampling-error SE paths (corrected_fsT/reliability/vfsLT) require
  # the full observed-indicator set, so they cannot be applied to lavaan's
  # per-missing-data-pattern blocks; fail here instead of deep inside
  # compute_fspars()/correct_evfs().
  # A FIML fit stores a (single, all-TRUE) pattern in @Data@Mp even for
  # complete data, so "has missing patterns" must mean "some observed
  # indicator is actually missing", i.e. a non-NULL pattern that is not
  # all-TRUE -- the same test the method = "mean" guard above uses. Using
  # !is.null(m) alone would wrongly lock out these SE paths on complete
  # FIML fits.
  has_miss_patterns <- any(vapply(
    object@Data@Mp,
    function(m) !is.null(m) && !all(m$pat),
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
    prior_cov = priors$cov,
    sum_items = sum_items
  )

  group_var <- object@Data@group
  group_col <- if (length(group_var) > 0) group_var else NULL

  out <- assemble_fs_blocks(
    blocks_by_group,
    format = format,
    group_col = group_col
  )

  # Group-level latent moments (effective / prior-adjusted), attached
  # post-assemble beside the vfsLT/reliability attachments below. The
  # values are the same psi_use/alpha_use the scoring above used
  # (see get_fs_blocks.lavaan(): prior_cov/prior_mean when supplied,
  # otherwise the per-group estimates; a named zero vector when the model
  # has no estimated mean structure). The shape mirrors fsT: unified ->
  # named list by group label; list -> direct attribute on each group data
  # frame plus a list-valued attribute on the outer list (multigroup).
  est_all <- lavInspect(object, what = "est", drop.list.single.group = FALSE)
  ngroups <- object@Data@ngroups
  group_labels <- object@Data@group.label
  if (length(group_labels) != ngroups) {
    # Single-group fits carry @Data@group.label as character(0); use the
    # same "" wrapper convention assemble_fs_blocks() applies to the
    # fsT/fsL attributes (psi/alpha must mirror the fsT shape).
    group_labels <- rep("", ngroups)
  }
  psi_g <- vector("list", ngroups)
  alpha_g <- vector("list", ngroups)
  names(psi_g) <- names(alpha_g) <- group_labels
  for (g in seq_len(ngroups)) {
    est_g <- est_all[[g]]
    psi_g[[g]] <- if (is.null(priors$cov)) est_g$psi else priors$cov
    if (!is.null(priors$mean)) {
      alpha_g[[g]] <- priors$mean
    } else if (is.null(est_g$alpha)) {
      # No (estimated) mean structure: the compute_fscore() zero-alpha
      # convention, named to match the score columns.
      alpha_g[[g]] <- setNames(
        rep(0, ncol(est_g$lambda)),
        colnames(est_g$lambda)
      )
    } else {
      a_est <- est_g$alpha
      alpha_g[[g]] <- setNames(
        as.numeric(a_est),
        if (!is.null(rownames(a_est))) {
          rownames(a_est)
        } else {
          colnames(est_g$lambda)
        }
      )
    }
  }
  if (format == "unified") {
    attr(out, "psi") <- psi_g
    attr(out, "alpha") <- alpha_g
  } else if (ngroups == 1) {
    attr(out, "psi") <- psi_g[[1L]]
    attr(out, "alpha") <- alpha_g[[1L]]
  } else {
    for (g in seq_len(ngroups)) {
      attr(out[[g]], "psi") <- psi_g[[g]]
      attr(out[[g]], "alpha") <- alpha_g[[g]]
    }
    attr(out, "psi") <- psi_g
    attr(out, "alpha") <- alpha_g
  }

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
    # Dimensionality from the model estimates, not from the output
    # attributes: format = "list" carries the SG fsb attribute as a bare
    # vector while "unified" wraps it in a per-group list, so the attr
    # shape cannot be tested uniformly. psi is q x q for the q latent
    # variables (SG: top-level est element; MG: per-group list, and all
    # groups share the latent structure, so group 1 suffices).
    n_latent <- if (ngroups > 1) nrow(est[[1]]$psi) else nrow(est$psi)
    multifactor <- n_latent > 1
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
  if (!is.null(product)) {
    if (ngroups > 1) {
      stop("'product' is not supported for multi-group models (v1); single-group lavaan models only.", call. = FALSE)
    }
    out <- compute_fs_prod(out, product = product)
  }
  out
}

get_fs_blocks.merMod <- function(
  object,
  method = c("EB", "ML"),
  legacy_names = FALSE,
  ...
) {
  method <- match.arg(method)
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

  # Random-effects design Z (sparse n x sum_p matrix from lme4, column-
  # major). First-term fold invariant -- Zden folds the first RE term's
  # columns onto the num_re coefficient columns of each row's own cluster:
  # (a) term-major layout: the first term's columns are the FIRST
  #     num_re * n_clus columns of getME("Z"); level j of its grouping
  #     factor occupies columns (j - 1) * num_re + seq_len(num_re) --
  #     the same layout the per-cluster slicing below relies on;
  # (b) first-term nonzeros have disjoint support per level (every
  #     nonzero lies in a row of its own level's cluster), so row i of
  #     Zden holds the first term's coefficients for i's own cluster;
  # (c) multi-term models make Z wider; only the first term's block is
  #     folded, matching the `[[1]]` convention used for b/cnms/flist
  #     below. The random design -- not the fixed design -- determines
  #     Kz: with Z != X the fixed-design code produces non-conformable
  #     products. Zden is n x num_re (tiny) instead of the
  #     n x (num_re * n_clus) dense matrix as.matrix() used to
  #     allocate (multi-GB on large fits).
  Zsp <- lme4::getME(object, "Z")
  if (!inherits(Zsp, "CsparseMatrix")) {
    stop("`lme4::getME(object, 'Z')` must return a column-compressed ",
         "sparse matrix (got: ", class(Zsp)[1], ").", call. = FALSE)
  }
  stopifnot(ncol(Zsp) >= num_re * n_clus)
  n1 <- Zsp@p[(num_re * n_clus) + 1L]   # nnz in the first p * G columns
  cc0 <- rep(seq_len(num_re * n_clus),
            diff(Zsp@p[seq_len(num_re * n_clus + 1L)]))
  # Fold safety net: the scatter below (m[cbind(i, j)] <- v) silently
  # keeps the LAST value for duplicate (i, j) pairs, so a future change
  # to lme4's Z column layout would mis-score silently instead of
  # erroring. Under the invariants above, a row has at most one nonzero
  # per folded column (its own level's block is the only block with
  # nnz in that row, and each folded column maps from exactly one raw
  # column within the level), so this fires only if the invariant
  # breaks; all-zero rows/levels simply yield no nnz and cannot trip
  # it. Cost: one vectorized O(nnz) pass -- negligible.
  if (any(duplicated(cbind(Zsp@i[seq_len(n1)], (cc0 - 1L) %% num_re)))) {
    stop(
      "internal error: first-term Z fold produced duplicate (row, ",
      "coefficient) scatter pairs; the lme4 Z column layout assumed by ",
      "get_fs_blocks.merMod() (first term's ", num_re * n_clus,
      " columns, level-major) no longer holds.",
      call. = FALSE
    )
  }
  Zden <- matrix(0, nrow = nrow(Zsp), ncol = num_re)
  Zden[cbind(Zsp@i[seq_len(n1)] + 1L, (cc0 - 1L) %% num_re + 1L)] <-
    Zsp@x[seq_len(n1)]
  s <- stats::sigma(object)

  if (method == "EB") {
    # Scaled random-effects covariance of the first RE term; the lme4-2.x
    # theta convention and the first-term-only contract are documented at
    # the definition of get_D().
    D <- get_D(object)
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
  } else {
    # ML (prior-free) scores: per cluster, the MLE of u_j with u treated as
    # fixed is the OLS fit of the cluster's fixed-effects-adjusted residuals
    # on its random-effects design block. Residuals are reconstructed as
    # y_j - X_j %*% beta (no offset()/weights support, same as the EB
    # scoring identity).
    y <- stats::model.response(stats::model.frame(object))
    X <- as.matrix(lme4::getME(object, "X"))
    beta <- lme4::fixef(object)
  }

  blocks <- vector("list", n_clus)

  for (j in seq_len(n_clus)) {
    idx <- case_idx[[j]]
    zj <- Zden[idx, seq_len(num_re), drop = FALSE]
    Kz <- crossprod(zj)

    if (method == "ML") {
      # u_hat_j = (Z'Z)^+ Z' r_j: Bartlett-analog (estimator uses no prior D);
      # fsL = I (score = u_j + (Z'Z)^+ Z' e_j), fsT = sigma^2 (Z'Z)^+ (Penrose:
      # (Z'Z)^+ Z'Z (Z'Z)^+ = (Z'Z)^+). ginv handles rank-deficient Z blocks
      # (e.g. a random slope on a within-cluster constant predictor) with the
      # minimum-norm solution.
      rj <- y[idx] - as.numeric(X[idx, , drop = FALSE] %*% beta)
      # ginv() runs a full SVD per cluster. For the single-coefficient case Kz
      # is a non-negative scalar, so its Moore-Penrose inverse is just 1/Kz
      # (exact, and far cheaper than the SVD); only the multi-coefficient,
      # possibly rank-deficient case keeps ginv's minimum-norm solution. A
      # plain solve() is NOT a safe substitute: it errors only on
      # exactly-singular Kz and silently returns garbage on near-singular Kz.
      Gz <- if (num_re == 1L && Kz[1L, 1L] > 0) {
        matrix(1 / Kz[1L, 1L], 1L, 1L)
      } else {
        MASS::ginv(Kz)
      }
      fs_row <- t(Gz %*% crossprod(zj, rj))
      fsL_j <- diag(num_re)
      fsT_j <- s^2 * Gz
      scoring_matrix_j <- Gz %*% t(zj)
    } else {
      DKz <- D %*% Kz
      inv_W <- solve(DKz + diag(nrow(Kz)))
      fsL_j <- DKz - DKz %*% inv_W %*% DKz
      # tcrossprod(A, B) == A %*% t(B) without materializing the transpose.
      fsT_j <- s^2 * tcrossprod(inv_W %*% DKz %*% D, inv_W)
      fs_row <- u_b[j, , drop = FALSE]
      scoring_matrix_j <- inv_W %*% D %*% t(zj)
    }

    colnames(fs_row) <- re_names

    # colnames = indicator/lv names, rownames = fs names (augment_fs convention)
    colnames(fsL_j) <- re_names
    rownames(fsL_j) <- fs_names
    attr(fs_row, "fsL") <- fsL_j

    rownames(fsT_j) <- colnames(fsT_j) <- fs_names

    # Scoring matrix: S_j %*% (y_j - X_j %*% beta) reproduces the scores
    # (EB / ranef for method "EB", per-cluster OLS for "ML"), where y_j/X_j
    # are the cluster's rows of the model frame and fixed-effects design.
    # See vignettes/scoring-matrices.Rmd.
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
# Reconstruct the SCALED random-effects covariance of the FIRST random-
# effects term of a merMod fit from its @theta parameters.
#
# lme4 >= 2.0 convention (read from the lme4 2.0.6 source; same scale as
# 1.x, which also stored the Cholesky of the scaled covariance):
#   - x@theta packs one block per RE term, in cnms() (formula) order; the
#     block for a term with p coefficients holds p * (p + 1) / 2 entries
#     and is the packed lower-triangular (column-major filled) Cholesky
#     factor L of the SCALED covariance D / sigma^2, where D =
#     VarCorr(x)[[term]] is the unscaled term covariance;
#   - equivalently VarCorr(x)[[term]] == sigma(x)^2 * tcrossprod(L). That
#     is exactly what lme4 2.x implements: mkVarCorr() splits @theta into
#     per-term blocks of length p * (p + 1) / 2 (the same idiom used
#     below), builds one "Covariance.us" object per block (whose
#     getLambda() fills the packed entries into lower.tri(column-major)),
#     and returns sc^2 * tcrossprod(L) for non-GLMM fits.
#   - lme4 2.x no longer attaches the "clen" (block lengths) attribute to
#     @theta -- it is missing even on single-term fits. A fallback parser
#     such as lme4::vec2mlist() therefore re-parses the whole mixed theta
#     as a SINGLE block of (sqrt(8L + 1) - 1) / 2 rows: exact only while
#     L is a triangular number (single term), a replacement-length warning
#     for a 2+1 term split, a mixed 3x3 block for 2+2. The split is thus
#     done explicitly from @cnms here, independent of that attribute.
#
# Returns the first term's p x p scaled covariance (p = the first term's
# number of coefficients), which equals VarCorr(x)[[1]] / sigma(x)^2 for
# LMMs. It must be the first term's, not a block-diagonal combination of
# all terms: get_fs_blocks.merMod() scores only the first term -- its
# clusters (flist[[1]]), Kz block (first-term Z fold) and scores
# (getME(x, "b"), term-major) all follow the flist[[1]] / cnms[[1]]
# convention, and the EB formulas consuming this D carry the explicit
# sigma^2 scale themselves (fsT_j = s^2 * (D Kz + I)^-1 D Kz D ...).
get_D <- function(object) {
  n1 <- length(object@cnms[[1L]])
  blk <- object@theta[seq_len(n1 * (n1 + 1L) / 2L)]
  L_D <- diag(nrow = n1)
  L_D[lower.tri(L_D, diag = TRUE)] <- blk
  tcrossprod(L_D)
}

#' @rdname get_fs
#' @param fsm Currently not used.
#' @param format Currently not used for `merMod` objects: the output is
#'        always a single data frame with one row per cluster (no `group`
#'        column).
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
  method = c("EB", "ML"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  fsm = FALSE,
  format = c("unified", "list"),
  legacy_names = FALSE,
  ...
) {
  method <- match.arg(method)
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

  blocks <- get_fs_blocks.merMod(
    object,
    method = method,
    legacy_names = legacy_names
  )

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
  # Bind at the matrix level rather than via rbind.data.frame: each block is
  # a 1-row, all-numeric data frame, and rbind.data.frame does expensive
  # per-frame name/class/factor matching that degrades as the (large) number
  # of clusters grows. as.matrix() reduces to plain numeric matrices, the
  # matrix rbind is a single C-level allocation, and as.data.frame() runs
  # exactly once. augment_fs() remains the single source of column names.
  out <- as.data.frame(
    do.call(rbind, lapply(aug_list, `as.matrix`)),
    check.names = FALSE
  )
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

  # Group-level latent moments: the (shared) prior covariance of the first
  # random-effects term, with dimnames renamed to the score names
  # (re_names) so they align with the fsL column names; random effects are
  # mean zero, so alpha is the named zero vector.
  psi_re <- as.matrix(lme4::VarCorr(object)[[1L]])
  rownames(psi_re) <- colnames(psi_re) <- re_names
  attr(out, "psi") <- psi_re
  attr(out, "alpha") <- setNames(rep(0, length(re_names)), re_names)

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

# ===========================================================================
# brms (Bayesian mixed models) support
#
# get_fs() method for brms's `brmsfit` object. Mirrors get_fs.merMod(): the
# posterior-mean random effects of the (single) random-effects term are the
# scores, and fsL / fsT / scoring_matrix come from the SAME EB / ML formulas
# get_fs_blocks.merMod() uses, generalised to a p x p random-effects
# covariance D. Point estimates are posterior means (brms's `Estimate`).
#
# brms is Suggests-only (guarded by require_brms()). The fixed design (X),
# random design (Z), and cluster structure are built directly from the model
# formula + data via reformulas (findbars / mkReTrms / nobars) -- bit-identical
# to lme4::getME -- so no throwaway lme4 fit is required. The random-effects
# covariance D is reconstructed from the posterior of the term's sd_ / cor_
# hyperparameters (brms_re_cov()).
#
# v1 scope: Gaussian (normal) family, a SINGLE random-effects term (one
# grouping factor) with any number of coefficients (random slopes; p = 1 or 2
# are the headline cases, p > 2 via the generic Cholesky). Multiple distinct
# grouping factors are deferred. `corrected_fsT` / `vfsLT` / `prior_*` are not
# supported (rejected), matching get_fs.merMod().
# ===========================================================================

require_brms <- function() {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "'brms' is required to extract factor scores from a brms model. ",
      "Install it with install.packages('brms') and retry.",
      call. = FALSE
    )
  }
}

# The plain model formula from a `brmsfit`. `object$formula` is a brmsformula
# (a list), not a bare formula; the main (mu) formula is its `formula` element
# (already a bare formula), so no coercion is needed.
brms_main_formula <- function(object) {
  f <- object$formula
  if (inherits(f, "brmsformula")) f <- f$formula
  f
}

# Reconstruct the (unscaled) p x p random-effects covariance D of a brms
# random-effects term from the posterior means of its sd_ / cor_ hyperparameters.
# brms names the per-term standard deviations `sd_<group>__<coef>` and, for
# p > 1, the `cor_<group>__<co_a>__<co_b>` draws are the lower-triangular
# elements of the Cholesky factor of the term's correlation matrix (brms's
# "generic" correlation parameterisation): L[i, j] (i > j) is the draw named
# `cor_<group>__<co_j>__<co_i>`. D = tcrossprod(diag(sd) %*% L). Coefficient
# names in the parameter strings carry no parentheses ("Intercept", not
# "(Intercept)"), so they are stripped from the lme4-style cnms.
brms_re_cov <- function(draws, group, cnms) {
  p <- length(cnms)
  cnms_p <- gsub("[()]", "", cnms)
  sd_nm <- paste0("sd_", group, "__", cnms_p)
  miss <- setdiff(sd_nm, colnames(draws))
  if (length(miss)) {
    stop("internal error: missing brms sd parameters: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  s <- vapply(sd_nm, function(nm) mean(draws[[nm]]), numeric(1L))
  # Cholesky factor L (lower) of the term's correlation matrix. The brms
  # `cor_<group>__<co_j>__<co_i>` (j < i) draws are the free off-diagonal
  # elements L[i, j] (brms's default "generic" parameterisation); the
  # diagonal is L[i, i] = sqrt(1 - sum_{j < i} L[i, j]^2).
  L <- matrix(0, p, p)
  if (p >= 2L) {
    for (i in 2L:p) {
      for (j in seq_len(i - 1L)) {
        nm <- paste0("cor_", group, "__", cnms_p[j], "__", cnms_p[i])
        if (!nm %in% colnames(draws)) {
          stop("internal error: missing brms cor parameter: ", nm, call. = FALSE)
        }
        L[i, j] <- mean(draws[[nm]])
      }
    }
  }
  for (i in seq_len(p)) {
    L[i, i] <- sqrt(max(0, 1 - sum(L[i, seq_len(i - 1L)]^2)))
  }
  R <- L %*% t(L)
  # D = diag(sd) R diag(sd); the diagonal scaling is built via diag(p) +
  # diag(S) <- s, NOT diag(<vector>), which is ambiguous for a length-1 vector.
  S <- diag(p)
  diag(S) <- s
  D <- S %*% R %*% S
  dimnames(D) <- list(cnms, cnms)
  D
}

# Per-level blocks for the single random-effects term of a brmsfit, reusing the
# EB / ML formulas of get_fs_blocks.merMod() generalised to p x p. One block per
# level (in the term's canonical level order), each carrying the score row, the
# p x p fsL / fsT, and the p x n_j scoring matrix.
get_fs_blocks.brmsfit <- function(object, method = c("EB", "ML"),
                                  legacy_names = FALSE, ...) {
  require_brms()
  method <- match.arg(method)
  f <- brms_main_formula(object)
  data <- object$data
  bars <- reformulas::findbars(f)
  if (length(bars) != 1L) {
    stop(
      "get_fs() for brms models supports a single random-effects term; this ",
      "model has ", length(bars), " (multiple distinct grouping factors are ",
      "not supported yet).",
      call. = FALSE
    )
  }
  rt <- reformulas::mkReTrms(bars, data, calc.lambdat = FALSE)
  group <- names(rt$flist)[[1L]]
  f1 <- rt$flist[[1L]]
  cnms <- rt$cnms[[1L]]
  p <- length(cnms)
  n_clus <- length(levels(f1))
  base_names <- paste0("u", seq_len(p) - 1L)
  re_names <- if (legacy_names) paste0(base_names, "_eb") else base_names
  fs_names <- paste0("fs_", re_names)

  # Random design (level-major, n x (p * n_clus)): row i holds i's own
  # cluster's p coefficients in columns (level(i) - 1) * p + seq_len(p) -- the
  # same invariant get_fs_blocks.merMod() relies on for the lme4 Z fold.
  # rt$Zt is [p * n_clus, n] column-compressed; scatter its nonzeros into
  # Zmat[n, p * n_clus] (= t(Zt)) via direct slot access, which avoids the
  # S3 `t`/`as.matrix` dispatch on the sparse in the package namespace.
  stopifnot(inherits(rt$Zt, "CsparseMatrix"))
  Zt <- rt$Zt
  Zmat <- matrix(0, nrow = ncol(Zt), ncol = nrow(Zt))
  if (length(Zt@x) > 0L) {
    obs_id <- rep(seq_len(ncol(Zt)), diff(Zt@p))
    Zmat[cbind(obs_id, Zt@i + 1L)] <- Zt@x
  }
  stopifnot(dim(Zmat)[1L] == nrow(data), dim(Zmat)[2L] == p * n_clus)
  case_idx <- split(seq_len(nrow(data)), f1)

  # Point estimates (posterior means).
  draws <- posterior::as_draws_df(object)
  sigma <- mean(draws[["sigma"]])
  D <- brms_re_cov(draws, group, cnms)
  Dsc <- D / sigma^2
  # Posterior-mean random effects, one p-vector per level, reindexed to the
  # canonical factor-level order used above. The per-level coefficients
  # (stats::coef, dispatched to brms's method) carry ALL fixed terms, so the
  # random part is the per-level deviation from the fixed effects
  # (brms::fixef), selecting only the term's (cnms) coefficients (the other
  # fixed terms deviate by ~0). Coefficient names are brms-style (no parens).
  co3d <- stats::coef(object)[[group]]
  all_terms <- dimnames(co3d)[[3L]]
  b_full <- as.matrix(co3d[, "Estimate", , drop = FALSE][, 1L, ])
  colnames(b_full) <- all_terms
  b_dev <- b_full - matrix(
    as.numeric(brms::fixef(object)[all_terms, "Estimate", drop = FALSE]),
    nrow = n_clus, ncol = length(all_terms), byrow = TRUE
  )
  b_mat <- b_dev[, match(gsub("[()]", "", cnms), all_terms), drop = FALSE]
  b_mat <- b_mat[match(levels(f1), dimnames(co3d)[[1L]]), , drop = FALSE]

  # The ML path additionally needs the response and fixed-effects design.
  if (method == "ML") {
    fx <- reformulas::nobars(f)
    X <- as.matrix(stats::model.matrix(fx, data))
    y <- stats::model.response(stats::model.frame(fx, data))
    # fixef() row names and model.matrix() column names follow the same
    # fixed-term order; match on the paren-stripped names so "(Intercept)"
    # (model.matrix) aligns with "Intercept" (brms).
    fixef_mat <- brms::fixef(object)
    beta <- as.numeric(fixef_mat[
      match(gsub("[()]", "", colnames(X)), gsub("[()]", "", rownames(fixef_mat))),
      "Estimate", drop = FALSE
    ])
  }

  blocks <- vector("list", n_clus)
  for (j in seq_len(n_clus)) {
    idx <- case_idx[[j]]
    zj <- Zmat[idx, (j - 1L) * p + seq_len(p), drop = FALSE]
    Kz <- crossprod(zj)
    if (method == "ML") {
      # Bartlett analog: prior-free per-cluster OLS of the fixed-effects-
      # adjusted residuals on the random design; (Z'Z)^+ minimum-norm solution.
      rj <- y[idx] - as.numeric(X[idx, , drop = FALSE] %*% beta)
      Gz <- if (p == 1L && Kz[1L, 1L] > 0) {
        matrix(1 / Kz[1L, 1L], 1L, 1L)
      } else {
        MASS::ginv(Kz)
      }
      fs_row <- t(Gz %*% crossprod(zj, rj))
      fsL_j <- diag(p)
      fsT_j <- sigma^2 * Gz
      scoring_j <- Gz %*% t(zj)
    } else {
      DKz <- Dsc %*% Kz
      inv_W <- solve(DKz + diag(p))
      fs_row <- b_mat[j, , drop = FALSE]
      fsL_j <- DKz - DKz %*% inv_W %*% DKz
      fsT_j <- sigma^2 * tcrossprod(inv_W %*% DKz %*% Dsc, inv_W)
      scoring_j <- inv_W %*% Dsc %*% t(zj)
    }
    colnames(fs_row) <- re_names
    colnames(fsL_j) <- re_names
    rownames(fsL_j) <- fs_names
    attr(fs_row, "fsL") <- fsL_j
    rownames(fsT_j) <- colnames(fsT_j) <- fs_names
    rownames(scoring_j) <- fs_names
    colnames(scoring_j) <- as.character(seq_len(nrow(zj)))
    blocks[[j]] <- list(
      case_idx = idx,
      fs = fs_row,
      fsL = fsL_j,
      fsT = fsT_j,
      fsb = NULL,
      scoring_matrix = scoring_j
    )
  }
  attr(blocks, "group") <- group
  attr(blocks, "D") <- D
  setNames(blocks, levels(f1))
}

#' @rdname get_fs
#' @export
get_fs.brmsfit <- function(
  object,
  method = c("EB", "ML"),
  corrected_fsT = FALSE,
  vfsLT = FALSE,
  format = c("unified", "list"),
  legacy_names = FALSE,
  ...
) {
  require_brms()
  method <- match.arg(method)
  if (!inherits(object, "brmsfit")) {
    stop("`object` must be a brms `brmsfit` model object.", call. = FALSE)
  }
  if (!identical(object$family$family, "gaussian")) {
    stop(
      "get_fs() for brms models supports Gaussian (normal) families only ",
      "(this model's family is '", object$family$family, "').",
      call. = FALSE
    )
  }
  prior_dots <- list(...)[c("prior_mean", "prior_cov")]
  if (!is.null(prior_dots$prior_mean) || !is.null(prior_dots$prior_cov)) {
    stop("'prior_mean'/'prior_cov' are not supported for brms objects.",
         call. = FALSE)
  }

  blocks <- get_fs_blocks.brmsfit(object, method = method,
                                  legacy_names = legacy_names)
  group <- attr(blocks, "group")

  # merMod-style assembly (one row per level of the RE term); mirrors
  # get_fs.merMod(), whose only model-specific input here is psi (the term's
  # random-effects covariance rather than lme4::VarCorr()).
  aug_list <- lapply(blocks, function(b) augment_fs(b$fs, b$fsT))
  out <- as.data.frame(
    do.call(rbind, lapply(aug_list, `as.matrix`)),
    check.names = FALSE
  )
  rownames(out) <- NULL
  if (legacy_names) {
    colnames(out) <- rename_legacy_fs_cols(colnames(out))
  }

  n_clus <- length(blocks)
  fsL_1 <- blocks[[1L]]$fsL
  re_names <- colnames(fsL_1)
  fs_names <- rownames(fsL_1)
  fsL_arr <- array(
    0, dim = c(nrow(fsL_1), ncol(fsL_1), n_clus),
    dimnames = list(fs_names, re_names, names(blocks))
  )
  fsT_arr <- array(
    0, dim = c(nrow(fsL_1), ncol(fsL_1), n_clus),
    dimnames = list(fs_names, fs_names, names(blocks))
  )
  for (j in seq_len(n_clus)) {
    fsL_arr[, , j] <- blocks[[j]]$fsL
    fsT_arr[, , j] <- blocks[[j]]$fsT
  }
  attr(out, "fsL") <- fsL_arr
  attr(out, "fsT") <- fsT_arr

  # Group-level latent moments: the term's (posterior-mean) random-effects
  # covariance, dimnamed to the score names; random effects are mean zero.
  psi_re <- attr(blocks, "D")
  rownames(psi_re) <- colnames(psi_re) <- re_names
  attr(out, "psi") <- psi_re
  attr(out, "alpha") <- setNames(rep(0, length(re_names)), re_names)

  attr(out, "scoring_matrix") <- setNames(
    lapply(blocks, function(b) b$scoring_matrix),
    names(blocks)
  )
  out
}

# ===========================================================================
# mirt (Item Response Theory) support
#
# get_fs() methods for mirt's S4 item-response fits. A fitted SingleGroupClass
# is scored by its EAP posterior means; the per-observation EAP posterior
# covariance (Vpost) feeds the shared regression-form matrix engine
# compute_lav_fs_matrices() (R/get_fscore_math.R). The latent covariance psi
# there is the FULL factor covariance the mirt model estimates (mirt_full_cov(),
# between-factor covariances included) -- NOT the unit-variance quadrature prior
# mirt uses to score. The mean is zero by default but can be overridden with
# `prior_mean` (alpha = 0 or prior_mean, shared across the EAP and the
# intercepts), giving the per-observation implied loading / error-covariance /
# intercept:
#   fsL_i = I - Vpost_i %*% solve(psi)         (univariate: 1 - SE^2)
#   fsT_i = fsL_i %*% Vpost_i                 (univariate: (1 - SE^2) * SE^2)
#   fsb_i = (I - fsL_i) %*% alpha = Vpost_i %*% solve(psi) %*% alpha
#         (univariate: SE^2 * alpha;  zero when alpha = 0)
# Because Vpost_i varies per observation, fsL/fsT/fsb are all attached as
# PER-ROW (list) attributes plus the `mirt_per_obs` marker so that fs_indiv()
# (resolve_per_obs(), R/fs_indiv.R) mints one block per row.
# ===========================================================================

# mirt is a Suggests-only dependency; guard every mirt:: call at entry with a
# clear, actionable message (never library()/require() a function body).
require_mirt <- function() {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop(
      "'mirt' is required to extract factor scores from a mirt model. ",
      "Install it with install.packages('mirt') and retry.",
      call. = FALSE
    )
  }
}

# Full factor variance-covariance matrix estimated by a mirt model, q x q,
# dimnames = factor names. Preferred route: coef(simplify = TRUE)$cov, mirt's
# own reconstruction, which fills the lower triangle positionally from the
# GroupPars parameter row (COV_11, COV_21, COV_22, ...) -- no parameter-name
# parsing, so it works for any number of factors. The fallback parses the
# COV_ij names from coef()$GroupPars; that scheme is unambiguous only for
# q < 10 (COV_111 could be factor (11, 1) or (1, 11)), so the fallback keeps
# the guard. simplify omits the `cov` element only for custom-density models,
# where the fallback still applies.
mirt_full_cov <- function(fit) {
  fn <- mirt::extract.mirt(fit, "factorNames")
  q <- length(fn)
  sc <- tryCatch(coef(fit, simplify = TRUE)$cov, error = function(e) NULL)
  if (is.matrix(sc) && all(dim(sc) == q) && all(is.finite(sc))) {
    dimnames(sc) <- list(fn, fn)
    return(sc)
  }
  gp1 <- coef(fit)$GroupPars[1L, , drop = TRUE]
  V <- matrix(0, q, q)
  for (nm in names(gp1)[grepl("^COV_", names(gp1))]) {
    s <- sub("^COV_", "", nm)
    if (nchar(s) != 2L) {
      stop("unsupported mirt factor covariance naming '", nm,
           "': q must be < 10.", call. = FALSE)
    }
    i <- as.integer(substr(s, 1L, 1L))
    j <- as.integer(substr(s, 2L, 2L))
    V[i, j] <- V[j, i] <- as.numeric(gp1[[nm]])
  }
  dimnames(V) <- list(fn, fn)
  V
}

#' @rdname get_fs
#' @param format Currently not used for `mirt` objects: the output is always a
#'        single data frame with one row per observation (no `group` column).
#' @export
get_fs.SingleGroupClass <- function(object, prior_mean = NULL,
                                    format = c("unified", "list"), ...) {
  require_mirt()
  if (!inherits(object, "SingleGroupClass")) {
    stop("`object` must be a mirt `SingleGroupClass` model object.", call. = FALSE)
  }
  format <- match.arg(format)  # accepted but unused: single-group mirt -> one df

  q <- mirt::extract.mirt(object, "nfact")
  fn <- mirt::extract.mirt(object, "factorNames")
  if (length(fn) != q) {
    stop(
      "internal error: nfact (", q, ") does not match the number of factor ",
      "names (", length(fn), ").",
      call. = FALSE
    )
  }
  fs_names <- paste0("fs_", fn)

  # Factor mean + covariance. object@Model$Theta is the quadrature NODE grid,
  # not the (mean, cov) prior, so the mean is not read from the fit -- alpha is
  # `prior_mean` or zero. psi must be the FULL estimated factor covariance
  # (including the between-factor covariances); for a 1-factor model mirt fixes
  # COV_11 = 1, so mirt_full_cov() returns diag(1) and 1-D is unchanged.
  alpha <- if (is.null(prior_mean)) {
    setNames(rep(0, q), fn)
  } else {
    validate_fs_priors(prior_mean, NULL, fn)$mean
  }
  psi <- mirt_full_cov(object)
  # mirt::fscores() prior-mean override, shared by both extraction calls so the
  # scores and the posterior covariances come from the same prior.
  fs_prior <- if (is.null(prior_mean)) list() else list(mean = alpha)

  # EAP posterior means. mirt re-adds completely-missing rows here, so `full`
  # has one row per observation with NA scores for the rows it could not
  # score. full.scores.SE is deliberately NOT requested: it would only add
  # mirt's per-observation posterior-covariance pass on this call (the SE
  # column is not used -- the emitted `_se` columns are the regression-form
  # sqrt(diag(fsT_i)) built from the acov call below).
  full <- do.call(mirt::fscores, c(list(object = object, full.scores = TRUE),
    fs_prior))
  full <- as.data.frame(full)

  # Per-observation EAP posterior covariance (acov). This early-returns BEFORE
  # mirt re-adds completely-missing rows, so the list has one q x q matrix per
  # SCORABLE observation only (named by scorable-row index).
  acov <- do.call(mirt::fscores, c(list(object = object, full.scores = TRUE,
    return.acov = TRUE), fs_prior))
  n <- nrow(full)
  # Plain numeric score matrix (one row per observation); column subsetting on
  # a matrix (not a data frame) keeps the `drop` argument honoured.
  score_mx <- as.matrix(full[, fn])
  if (n < 1 || length(acov) < 1) {
    stop("mirt returned no factor scores; the fit has no scorable observations.")
  }

  # Reconcile row alignment. `full` (score/SE call) includes the
  # completely-missing rows (NA scores); `acov` skips them. extract.mirt(
  # "completely_missing") names those original-data positions, so the scorable
  # rows are their complement -- in which order acov[[k]] is the k-th scorable
  # row (verified: diag(Vpost_k) == SE_row^2 for every scorable row).
  cm <- mirt::extract.mirt(object, "completely_missing")
  if (is.null(cm)) cm <- integer(0)
  keep <- !seq_len(n) %in% cm
  scorsc <- which(keep)
  if (length(acov) != length(scorsc)) {
    stop(
      "internal error: mirt posterior covariances (", length(acov),
      ") do not match the number of scorable rows (", length(scorsc), ").",
      call. = FALSE
    )
  }

  # Per-row regression-form matrices (one source of truth with the lavaan /
  # merMod paths: compute_lav_fs_matrices() with psi = I and the prior mean
  # alpha; alpha = 0 gives all-zero intercepts). fsb_i is a function of the
  # latent mean and the per-row shrinkage factor (Vpost_i %*% solve(psi)).
  fsL_list <- vector("list", n)
  fsT_list <- vector("list", n)
  fsb_list <- vector("list", n)
  naL <- matrix(NA_real_, q, q, dimnames = list(fs_names, fn))
  naT <- matrix(NA_real_, q, q, dimnames = list(fs_names, fs_names))
  naB <- setNames(rep(NA_real_, q), fs_names)
  for (k in seq_len(length(scorsc))) {
    i <- scorsc[k]
    Vpost_i <- as.matrix(acov[[k]])
    m_i <- compute_lav_fs_matrices(Vpost_i, psi, alpha, method = "regression")
    L_i <- m_i$fsL
    T_i <- m_i$fsT
    rownames(L_i) <- fs_names
    colnames(L_i) <- fn
    rownames(T_i) <- colnames(T_i) <- fs_names
    fsL_list[[i]] <- L_i
    fsT_list[[i]] <- T_i
    fsb_list[[i]] <- setNames(as.numeric(m_i$fsb), fs_names)
  }
  # Completely-missing rows: R2spa's NA-row convention (all-NA per-row block,
  # keeping reference dimnames for the column-name resolver).
  for (i in setdiff(seq_len(n), scorsc)) {
    fsL_list[[i]] <- naL
    fsT_list[[i]] <- naT
    fsb_list[[i]] <- naB
  }

  # Per-row se / loadings / error terms -- the shared value-only engine
  # fs_row_cols() (R/fs_indiv.R), applied row by row. fsb is passed as NULL so
  # get_fs() itself emits no intercept columns (fs_indiv() may emit them via
  # include_intercept = TRUE using the attached fsb).
  # fs_row_cols() layout (no intercept here): se (q) + loadings (q^2) +
  # error terms (q*(q+1)/2). The score columns are carried separately in
  # scores_df, so this width excludes them.
  K <- q + q * q + q * (q + 1L) / 2L
  vals <- matrix(NA_real_, nrow = n, ncol = K)
  # fs_row_cols() only reads nrow(fs) (always 1 here): one shared token.
  dummy_fs <- data.frame(x = 0L)
  for (k in seq_len(length(scorsc))) {
    i <- scorsc[k]
    vals[i, ] <- fs_row_cols(
      dummy_fs,
      fsL_list[[i]], fsT_list[[i]], NULL
    )[1L, , drop = FALSE]
  }

  # Assemble the canonical data frame. Column set + order is identical to what
  # fs_indiv(get_fs(m)) emits (shared naming via fs_row_colnames()):
  # fs_<fn> | fs_<fn>_se | <fn_j>_by_fs_<fn> (q^2, column-major per latent)
  # | ev_fs_<fn>/ecov_<a>_<b> (lower-tri row-major, i-outer j<=i). No
  # group/id columns.
  scores_df <- as.data.frame(score_mx)
  colnames(scores_df) <- fs_names
  nm <- fs_row_colnames(fsL_list[[1L]], fsT_list[[1L]])
  out <- as.data.frame(
    cbind(as.matrix(scores_df), vals),
    check.names = FALSE
  )
  colnames(out) <- c(fs_names, nm$se, nm$ld, nm$ev)

  # Per-row attributes + group-level latent moments + the per-obs marker that
  # fs_indiv()'s resolve_per_obs() dispatches on.
  attr(out, "fsT") <- fsT_list
  attr(out, "fsL") <- fsL_list
  attr(out, "fsb") <- fsb_list
  attr(out, "fs_pattern") <- list(label = seq_len(n), pat = NULL)
  attr(out, "psi") <- psi
  attr(out, "alpha") <- alpha
  attr(out, "mirt_per_obs") <- TRUE
  out
}

# Per-group factor (co)variances for a mirt MultipleGroupClass. Returns a list
# with `group_names` (character, one per group, in mirt's code order) and `psi`
# (a named list, one q x q full covariance per group, keyed by label). Each
# group's covariance is the estimated factor covariance of that group's
# single-group model, so this is the natural multi-group generalisation of
# mirt_full_cov(). Under the usual metric invariance every group's factor is
# fixed to (0, I), so all entries coincide; they differ only under free_means /
# free_var or a between-factor covariance model.
mirt_group_pars <- function(object) {
  group_names <- mirt::extract.mirt(object, "groupNames")
  K <- length(group_names)
  psi <- lapply(seq_len(K), function(k) {
    mirt_full_cov(mirt::extract.group(object, k))
  })
  names(psi) <- group_names
  list(group_names = group_names, psi = psi)
}

#' @rdname get_fs
#' @export
get_fs.MultipleGroupClass <- function(object, prior_mean = NULL,
                                      format = c("unified", "list"), ...) {
  require_mirt()
  if (!inherits(object, "MultipleGroupClass")) {
    stop("`object` must be a mirt `MultipleGroupClass` model object.", call. = FALSE)
  }
  format <- match.arg(format)  # accepted but unused: mirt -> one per-obs df + group

  q <- mirt::extract.mirt(object, "nfact")
  fn <- mirt::extract.mirt(object, "factorNames")
  if (length(fn) != q) {
    stop(
      "internal error: nfact (", q, ") does not match the number of factor ",
      "names (", length(fn), ").",
      call. = FALSE
    )
  }
  fs_names <- paste0("fs_", fn)

  # Per-group factor (co)variances (named list, length K) + the group labels in
  # mirt's code order. mirt drops completely-missing rows from every extraction,
  # so all per-observation quantities below are reconciled against the scorable
  # rows exactly as in the single-group path.
  gp <- mirt_group_pars(object)
  psi <- gp$psi
  group_names <- gp$group_names

  # alpha: the prior mean used to build the regression-form intercept. mirt's
  # multi-group EAP is centred on a standard-normal prior per group (the group
  # mean is carried by the item intercepts, not the factor prior), so it is 0
  # by default -- identical to the single-group path -- or the user's
  # prior_mean (a length-q vector applied to every group).
  alpha <- if (is.null(prior_mean)) {
    setNames(rep(0, q), fn)
  } else {
    validate_fs_priors(prior_mean, NULL, fn)$mean
  }
  fs_prior <- if (is.null(prior_mean)) list() else list(mean = alpha)

  # EAP posterior means. mirt re-adds completely-missing rows here, so `full`
  # has one row per observation with NA scores for the rows it could not
  # score. full.scores.SE is deliberately NOT requested: it would only add
  # mirt's per-observation posterior-covariance pass on this call (the SE
  # column is not used -- the emitted `_se` columns are the regression-form
  # sqrt(diag(fsT_i)) built from the acov call below).
  full <- do.call(mirt::fscores, c(list(object = object, full.scores = TRUE),
    fs_prior))
  full <- as.data.frame(full)
  n <- nrow(full)
  # Plain numeric score matrix (one row per observation).
  score_mx <- as.matrix(full[, fn])

  # Per-observation EAP posterior covariance (acov): scorable rows only.
  acov <- do.call(mirt::fscores, c(list(object = object, full.scores = TRUE,
    return.acov = TRUE), fs_prior))

  if (n < 1 || length(acov) < 1) {
    stop("mirt returned no factor scores; the fit has no scorable observations.")
  }

  # Row alignment. `full` (score/SE call) includes the completely-missing rows
  # (NA scores); `acov` skips them. extract.mirt("completely_missing") names
  # those original-data positions, so the scorable rows are their complement, in
  # which order acov[[k]] is the k-th scorable row.
  cm <- mirt::extract.mirt(object, "completely_missing")
  if (is.null(cm)) cm <- integer(0)
  keep <- !seq_len(n) %in% cm
  scorsc <- which(keep)
  if (length(acov) != length(scorsc)) {
    stop(
      "internal error: mirt posterior covariances (", length(acov),
      ") do not match the number of scorable rows (", length(scorsc), ").",
      call. = FALSE
    )
  }

  # Per-row group. extract.mirt("group") returns the group LABEL of each
  # scorable observation (completely-missing rows are dropped, so it has length
  # == scorsc); map it onto the scorable row positions. The per-group psi that
  # each scorable row uses is indexed by matching its label to group_names.
  grp_scor <- as.character(mirt::extract.mirt(object, "group"))
  if (length(grp_scor) != length(scorsc)) {
    stop(
      "internal error: mirt per-observation groups (", length(grp_scor),
      ") do not match the number of scorable rows (", length(scorsc), ").",
      call. = FALSE
    )
  }
  psi_idx <- match(grp_scor, group_names)
  if (any(is.na(psi_idx))) {
    stop("internal error: a scorable group value is not among the model's group names.",
         call. = FALSE)
  }

  # Per-row regression-form matrices -- the shared source of truth with the
  # lavaan / merMod / single-group-mirt paths (compute_lav_fs_matrices()), using
  # each row's own group factor covariance.
  fsL_list <- vector("list", n)
  fsT_list <- vector("list", n)
  fsb_list <- vector("list", n)
  naL <- matrix(NA_real_, q, q, dimnames = list(fs_names, fn))
  naT <- matrix(NA_real_, q, q, dimnames = list(fs_names, fs_names))
  naB <- setNames(rep(NA_real_, q), fs_names)
  for (k in seq_len(length(scorsc))) {
    i <- scorsc[k]
    Vpost_i <- as.matrix(acov[[k]])
    psi_i <- psi[[ psi_idx[k] ]]
    m_i <- compute_lav_fs_matrices(Vpost_i, psi_i, alpha, method = "regression")
    L_i <- m_i$fsL
    T_i <- m_i$fsT
    rownames(L_i) <- fs_names
    colnames(L_i) <- fn
    rownames(T_i) <- colnames(T_i) <- fs_names
    fsL_list[[i]] <- L_i
    fsT_list[[i]] <- T_i
    fsb_list[[i]] <- setNames(as.numeric(m_i$fsb), fs_names)
  }
  # Completely-missing rows: R2spa's NA-row convention (all-NA per-row block).
  for (i in setdiff(seq_len(n), scorsc)) {
    fsL_list[[i]] <- naL
    fsT_list[[i]] <- naT
    fsb_list[[i]] <- naB
  }

  # Per-row se / loadings / error terms -- the shared value-only engine
  # fs_row_cols(), applied row by row (no intercept columns emitted here;
  # fs_indiv() may emit them via include_intercept = TRUE using the fsb attr).
  K <- q + q * q + q * (q + 1L) / 2L
  vals <- matrix(NA_real_, nrow = n, ncol = K)
  # fs_row_cols() only reads nrow(fs) (always 1 here): one shared token.
  dummy_fs <- data.frame(x = 0L)
  for (k in seq_len(length(scorsc))) {
    i <- scorsc[k]
    vals[i, ] <- fs_row_cols(
      dummy_fs,
      fsL_list[[i]], fsT_list[[i]], NULL
    )[1L, , drop = FALSE]
  }

  # Assemble the canonical data frame (identical column set + order to the
  # single-group path, shared naming via fs_row_colnames()), then append the
  # trailing `group` column.
  scores_df <- as.data.frame(score_mx)
  colnames(scores_df) <- fs_names
  nm <- fs_row_colnames(fsL_list[[1L]], fsT_list[[1L]])
  out <- as.data.frame(
    cbind(as.matrix(scores_df), vals),
    check.names = FALSE
  )
  colnames(out) <- c(fs_names, nm$se, nm$ld, nm$ev)

  # The group column: one value per observation; NA for completely-missing rows
  # (which carry no group, mirroring the all-NA row convention).
  g_full <- rep(NA_character_, n)
  g_full[scorsc] <- grp_scor
  out$group <- factor(g_full, levels = group_names)

  # Per-row attributes + group-level latent moments + the per-obs marker that
  # fs_indiv()'s resolve_per_obs() dispatches on. `psi` is a named list (one q x
  # q per group); fs_indiv() does not read it on the per-obs path.
  attr(out, "fsT") <- fsT_list
  attr(out, "fsL") <- fsL_list
  attr(out, "fsb") <- fsb_list
  attr(out, "fs_pattern") <- list(label = seq_len(n), pat = NULL)
  attr(out, "psi") <- psi
  attr(out, "alpha") <- alpha
  attr(out, "mirt_per_obs") <- TRUE
  out
}
