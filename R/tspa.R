#' Two-Stage Path Analysis
#'
#' Fit a two-stage path analysis (2S-PA) model.
#'
#' @param model A string variable describing the structural path model,
#'              in \code{lavaan} syntax.
#' @param data A data frame containing factor scores.
#' @param reliability A numeric vector representing the reliability indexes
#'                    of each latent factor. Currently \code{tspa()} does not
#'                    support the reliability argument. Please use \code{se}.
#' @param se Deprecated to avoid conflict with the argument of the same name
#'           in [lavaan::lavaan()].
#' @param se_fs A numeric vector representing the standard errors of each
#'              factor score variable for single-group 2S-PA. A list or data
#'              frame storing the standard errors of each group in each latent
#'              factor for multigroup 2S-PA.
#' @param fsT An error variance-covariance matrix of the factor scores, which
#'            can be obtained from the output of \code{get_fs()} using
#'            \code{attr()} with the argument \code{which = "fsT"}.
#' @param fsL A matrix of loadings and cross-loadings from the
#'            latent variables to the factor scores \code{fs}, which
#'            can be obtained from the output of \code{get_fs()} using
#'            \code{attr()} with the argument \code{which = "fsL"}.
#'            For details see the multiple-factors vignette:
#'            \code{vignette("multiple-factors", package = "R2spa")}.
#' @param fsb A vector of intercepts for the factor scores \code{fs}, which can
#'            be obtained from the output of \code{get_fs()} using \code{attr()}
#'            with the argument \code{which = "fsb"}.
#' @param ... Additional arguments passed to \code{\link[lavaan]{sem}}. See
#'            \code{\link[lavaan]{lavOptions}} for a complete list.
#' @return An object of class \code{lavaan}, with an attribute \code{tspaModel}
#'         that contains the model syntax.
#'
#' @export
#'
#' @examples
#' library(lavaan)
#'
#' # single-group, two-factor example, factor scores obtained separately
#' # get factor scores
#' fs_dat_ind60 <- get_fs(object = PoliticalDemocracy,
#'                        model = "ind60 =~ x1 + x2 + x3")
#' fs_dat_dem60 <- get_fs(object = PoliticalDemocracy,
#'                        model = "dem60 =~ y1 + y2 + y3 + y4")
#' fs_dat <- cbind(fs_dat_ind60, fs_dat_dem60)
#' # tspa model
#' tspa(model = "dem60 ~ ind60", data = fs_dat,
#'      se_fs = c(ind60 = fs_dat_ind60[1, "fs_ind60_se"],
#'                dem60 = fs_dat_dem60[1, "fs_dem60_se"]))
#'
#' # single-group, three-factor example
#' mod2 <- "
#'   # latent variables
#'     ind60 =~ x1 + x2 + x3
#'     dem60 =~ y1 + y2 + y3 + y4
#'     dem65 =~ y5 + y6 + y7 + y8
#' "
#' fs_dat2 <- get_fs(PoliticalDemocracy, model = mod2, std.lv = TRUE)
#' tspa(model = "dem60 ~ ind60
#'               dem65 ~ ind60 + dem60",
#'      data = fs_dat2,
#'      fsT = attr(fs_dat2, "fsT"),
#'      fsL = attr(fs_dat2, "fsL"))
#'
#' # multigroup, two-factor example
#' mod3 <- "
#'   # latent variables
#'     visual =~ x1 + x2 + x3
#'     speed =~ x7 + x8 + x9
#' "
#' fs_dat3 <- get_fs(HolzingerSwineford1939, model = mod3, std.lv = TRUE,
#'                   group = "school")
#' tspa(model = "visual ~ speed",
#'      data = fs_dat3,
#'      fsT = attr(fs_dat3, "fsT"),
#'      fsL = attr(fs_dat3, "fsL"),
#'      group = "school")
#'
#' # multigroup, three-factor example
#' mod4 <- "
#'   # latent variables
#'     visual =~ x1 + x2 + x3
#'     textual =~ x4 + x5 + x6
#'     speed =~ x7 + x8 + x9
#' "
#' fs_dat4 <- get_fs(HolzingerSwineford1939, model = mod4, std.lv = TRUE,
#'                   group = "school")
#' tspa(model = "visual ~ speed
#'               textual ~ visual + speed",
#'      data = fs_dat4,
#'      fsT = attr(fs_dat4, "fsT"),
#'      fsL = attr(fs_dat4, "fsL"),
#'      group = "school")
#'
#' # get factor scores
#' fs_dat_visual <- get_fs(object = HolzingerSwineford1939,
#'                         model = "visual =~ x1 + x2 + x3",
#'                         group = "school",
#'                         format = "list")
#' fs_dat_speed <- get_fs(object = HolzingerSwineford1939,
#'                        model = "speed =~ x7 + x8 + x9",
#'                        group = "school",
#'                        format = "list")
#' fs_hs <- cbind(do.call(rbind, fs_dat_visual),
#'                do.call(rbind, fs_dat_speed))
#'
#' # tspa model
#' tspa(model = "visual ~ speed",
#'      data = fs_hs,
#'      se_fs = data.frame(visual = c(0.3391326, 0.311828),
#'                         speed = c(0.2786875, 0.2740507)),
#'      group = "school",
#'      group.equal = "regressions")
#'
#' # manually adding equality constraints on the regression coefficients
#' tspa(model = "visual ~ c(b1, b1) * speed",
#'      data = fs_hs,
#'      se_fs = list(visual = c(0.3391326, 0.311828),
#'                   speed = c(0.2786875, 0.2740507)),
#'      group = "school")


tspa <- function(model, data, reliability = NULL, se = "standard",
                 se_fs = NULL, fsT = NULL, fsL = NULL, fsb = NULL, ...) {

  if (!inherits(model, "character")) {
    stop("The structural path model provided is not a string.")
  }

  if (!is.null(reliability)) {
    stop("tspa() currently does not support reliability model")
  }
  if (!is.character(se)) {
    warning("using `se` to set se for factor scores is deprecated. ",
            "use `se_fs` instead.")
  }

  if (!is.data.frame(se_fs)) {
    se_fs <- as.data.frame(as.list(se_fs))
  }
  if (xor(is.null(fsT), is.null(fsL))) {
    stop("Please provide both or none of fsT and fsL.")
  }

  # A plain matrix stands for a single group, so a length-1 list may be mixed
  # with a plain matrix (single-group unified get_fs() output); only differing
  # group counts are an error.
  if (!is.null(fsT)) {
    nT <- if (is.list(fsT)) length(fsT) else 1
    nL <- if (is.list(fsL)) length(fsL) else 1
    if (nT != nL) {
      stop(
        if (nT > nL) {
          "'fsL' must be a list of the same length as 'fsT' for a multigroup model."
        } else {
          "'fsT' must be a list of the same length as 'fsL' for a multigroup model."
        }
      )
    }
  }
  multigroup <- if (!is.null(fsT)) {
    is.list(fsT) && length(fsT) > 1
  } else {
    nrow(se_fs) > 1
  }

  if (!is.null(fsT)) {
    if (multigroup) {
      fs_names <- colnames(fsT[[1]])
      dat_names <- if (is.data.frame(data)) names(data) else names(data[[1]])
    } else {
      fs_names <- if (is.list(fsT)) colnames(fsT[[1]]) else colnames(fsT)
      dat_names <- names(data)
    }
    names_match <- lapply(fs_names, function(x) x %in% dat_names) |> unlist()
    if (any(!names_match)) {
      stop(
        "Names of factor score variables do not match those in the input data."
      )
    }
  }

  if (multigroup && is.null(list(...)[["group"]])) {
    stop("Please specify 'group = ' to fit a multigroup model in lavaan.")
  }

  if (is.null(fsT)) { # single-factor measurement model
    # Product-score columns (get_fs_int: `fs_a:fs_b`) are not valid lavaan
    # variable names; the schema's generated model name for latent `v` is
    # `fs_v`, so a matching product-score column is aliased into a working
    # copy of the data. Manual pre-renames keep working (the alias is a
    # no-op when `fs_v` already exists).
    data <- tspa_sf_alias(data, se_fs)$data
    tspaModel <- tspa_sf(model, data, se_fs)
  } else { # multi-factor measurement model
    tspaModel <- tspa_mf(model, data, fsT, fsL, fsb)
    if (inherits(data, "list")) {
      data <- do.call(rbind, data)
    }
  }

  tspa_fit <- sem(model = tspaModel,
                  data  = data,
                  se = se,
                  ...)
  attr(tspa_fit, "tspaModel") <- tspaModel
  if (!is.null(fsT)) {
    attr(tspa_fit, "fsT") <- fsT
    attr(tspa_fit, "fsL") <- fsL
  }
  attr(tspa_fit, "tspa_call") <- match.call()
  return(tspa_fit)
}

# ---------------------------------------------------------------------------
# Stage-2 model schema (PLAN 04): owned and frozen by R2spa, not by lavaan.
# One row per (statement term, group):
#   lhs   statement left-hand side (NA for verbatim user rows)
#   op    "=~" | "~~" | "~" | "raw" ("raw" rows carry user syntax)
#   rhs   right-hand variable ("1" for intercepts); "raw" rows carry the
#         user's line verbatim
#   value per-group fixed value (NA for "raw" rows)
#   free  0 = fixed (every injected row); NA for "raw" rows
#   group 1-based group index (NA for "raw" rows)
#   label R2spa-generated label for the statement (__r2spa_ldN__ /
#         __r2spa_evN__ / __r2spa_intN__): a generated namespace user
#         labels/variable names cannot collide with; carried as row
#         metadata, never emitted into the syntax
#   kind  "user" | "struct" | "error_var" | "error_cov" | "intercept"
# ---------------------------------------------------------------------------

tspa_row <- function(lhs, op, rhs, value, group, kind, label) {
  data.frame(
    lhs = if (is.null(lhs)) NA_character_ else lhs,
    op = op,
    rhs = rhs,
    value = value,
    free = if (op == "raw") NA_integer_ else 0L,
    group = if (is.null(group)) NA_integer_ else as.integer(group),
    label = label,
    kind = kind,
    stringsAsFactors = FALSE
  )
}

# The user model string, carried verbatim as a single "raw" row (no
# re-parsing and no round-trip through line splitting, so leading/trailing
# newlines survive exactly as written).
tspa_user_rows <- function(model) {
  tspa_row(NA, "raw", paste0(model, collapse = "\n"), NA, NA, "user", NA)
}

# Single-factor (se_fs) schema: per latent, one fixed-loading struct row
# and one error-variance row per group; values follow the se_fs rows
# (groups) in order.
tspa_schema_sf <- function(model, se) {
  var <- colnames(se)
  fs <- paste0("fs_", var)
  ng <- nrow(se)
  rows <- list(tspa_user_rows(model))
  for (k in seq_along(var)) {
    lab <- paste0("__r2spa_ld", k, "__")
    for (g in seq_len(ng)) {
      rows[[length(rows) + 1L]] <- tspa_row(
        var[k], "=~", fs[k], 1, g, "struct", lab
      )
    }
    ev_lab <- paste0("__r2spa_ev", k, "__")
    for (g in seq_len(ng)) {
      rows[[length(rows) + 1L]] <- tspa_row(
        fs[k], "~~", fs[k], se[g, k]^2, g, "error_var", ev_lab
      )
    }
  }
  do.call(rbind, rows)
}

# Multi-factor (fsT/fsL/fsb) schema: per latent, one struct row per score
# term and per group; error rows follow the lower triangle (incl. diagonal)
# of fsT in column-major order — the legacy per-group value routing made
# explicit and unit-testable; per-score intercept rows when fsb is given.
tspa_schema_mf <- function(model, fsT, fsL, fsb) {
  # `fsT`/`fsL` are plain matrices for a single-group model or named lists
  # of them for a multigroup model. Single-group unified get_fs() output
  # carries length-1 list attributes, so either shape is accepted on either
  # side (e.g. list-valued `fsT` with a plain identity `fsL` for Bartlett).
  if (is.list(fsT) && length(fsT) > 1) {
    if (!is.list(fsL) || length(fsL) != length(fsT)) {
      stop("'fsL' must be a list of the same length as 'fsT' for a ",
           "multigroup model.")
    }
    ngroup <- length(fsT)
    L_list <- fsL
    T_list <- fsT
  } else if (is.list(fsL) && length(fsL) > 1) {
    if (!is.list(fsT) || length(fsT) != length(fsL)) {
      stop("'fsT' must be a list of the same length as 'fsL' for a ",
           "multigroup model.")
    }
    ngroup <- length(fsL)
    L_list <- fsL
    T_list <- fsT
  } else {
    ngroup <- 1
    L_list <- if (is.list(fsL)) fsL else list(fsL)
    T_list <- if (is.list(fsT)) fsT else list(fsT)
  }
  fsL1 <- L_list[[1]]
  var <- colnames(fsL1)
  fs <- rownames(fsL1)

  rows <- list(tspa_user_rows(model))
  for (k in seq_along(var)) {
    lab <- paste0("__r2spa_ld", k, "__")
    for (i in seq_along(fs)) {
      for (g in seq_len(ngroup)) {
        rows[[length(rows) + 1L]] <- tspa_row(
          var[k], "=~", fs[i], L_list[[g]][i, k], g, "struct", lab
        )
      }
    }
  }
  ev_count <- 0L
  tri <- which(!upper.tri(T_list[[1]]), arr.ind = TRUE)
  for (k in seq_len(nrow(tri))) {
    i <- tri[k, 1]
    j <- tri[k, 2]
    ev_count <- ev_count + 1L
    lab <- paste0("__r2spa_ev", ev_count, "__")
    kind <- if (i == j) "error_var" else "error_cov"
    for (g in seq_len(ngroup)) {
      rows[[length(rows) + 1L]] <- tspa_row(
        fs[i], "~~", fs[j], T_list[[g]][i, j], g, kind, lab
      )
    }
  }
  if (!is.null(fsb)) {
    B_list <- if (is.list(fsb)) fsb else list(fsb)
    for (i in seq_along(fs)) {
      lab <- paste0("__r2spa_int", i, "__")
      for (g in seq_len(ngroup)) {
        rows[[length(rows) + 1L]] <- tspa_row(
          fs[i], "~", "1", B_list[[g]][i], g, "intercept", lab
        )
      }
    }
  }
  do.call(rbind, rows)
}

# Statements of one schema: rows with the given identity grouped into
# consecutive statements (order of first appearance).
tspa_statements <- function(sch, id) {
  id <- match(id, unique(id))
  starts <- c(1L, which(diff(id) != 0L) + 1L)
  lapply(starts, function(i) sch[id == id[i], , drop = FALSE])
}

# Per-group values of a single-term statement, in group order.
tspa_stmt_values <- function(st) {
  gs <- sort(unique(st$group))
  unlist(lapply(gs, function(g) st$value[st$group == g]))
}

# c(...) value string for a possibly multi-term statement; `terms` are the
# ordered unique rhs values (one row value per group each).
tspa_stmt_cvals <- function(st, terms) {
  trm <- match(st$rhs, terms)
  paste(
    vapply(seq_along(terms), function(k) {
      sub <- st[trm == k, , drop = FALSE]
      paste0("c(", paste(tspa_stmt_values(sub), collapse = ", "), ")")
    }, character(1)),
    collapse = " + "
  )
}

# The single renderer (PLAN 04): schema -> lavaan model syntax string.
# Reproduces the legacy string builders character-for-character, including
# their per-path spacing quirks, so parameter row order, estimates, and
# vcov() are provably unchanged (Phase 2 A/B gate).
tspa_render <- function(sch, style = c("sf", "mf")) {
  style <- match.arg(style)
  user_lines <- sch$rhs[sch$kind == "user"]
  struct <- sch[sch$kind == "struct", , drop = FALSE]
  errors <- sch[sch$kind %in% c("error_var", "error_cov"), , drop = FALSE]
  ints <- sch[sch$kind == "intercept", , drop = FALSE]
  if (style == "sf") {
    latent_var_str <- paste(
      vapply(tspa_statements(struct, struct$lhs), function(st) {
        paste0(st$lhs[1], "=~ ", tspa_stmt_cvals(st, st$rhs[1]),
               " * ", st$rhs[1], "\n")
      }, character(1)),
      collapse = ""
    )
    error_constraint_str <- paste(
      vapply(tspa_statements(errors, paste(errors$lhs, errors$rhs,
                                           sep = "|")),
             function(st) {
               paste0(st$lhs[1], "~~ ", tspa_stmt_cvals(st, st$rhs[1]),
                      " * ", st$rhs[1], "\n")
             }, character(1)),
      collapse = ""
    )
    elems <- c(
      "# latent variables (indicated by factor scores)",
      latent_var_str,
      "# constrain the errors",
      error_constraint_str,
      "# structural model",
      paste(user_lines, collapse = "\n")
    )
  } else {
    latent_var_str <- vapply(
      tspa_statements(struct, struct$lhs),
      function(st) {
        terms <- unique(st$rhs)
        loadings_c <- paste(
          vapply(terms, function(t) {
            tr <- st[st$rhs == t, , drop = FALSE]
            paste0("c(", paste0(tspa_stmt_values(tr), collapse = ", "),
                   ") * ", t)
          }, character(1)),
          collapse = " + "
        )
        paste("# latent variables (indicated by factor scores)\n",
              st$lhs[1], "=~", loadings_c)
      },
      character(1)
    )
    error_constraint_str <- vapply(
      tspa_statements(errors, paste(errors$lhs, errors$rhs, sep = "|")),
      function(st) {
        paste0("# constrain the errors\n", st$lhs[1], " ~~ ",
               tspa_stmt_cvals(st, st$rhs[1]), " * ", st$rhs[1])
      },
      character(1)
    )
    if (nrow(ints) > 0) {
      ng <- length(sort(unique(ints$group)))
      intercept_constraint <- vapply(
        tspa_statements(ints, ints$lhs),
        function(st) {
          vals <- tspa_stmt_values(st)
          intercepts <- if (ng == 1) {
            vals[1]
          } else {
            paste0("c(", paste0(vals, collapse = ", "), ")")
          }
          paste0("# constrain the intercepts\n", st$lhs[1], " ~ ",
                 intercepts, " * 1")
        },
        character(1)
      )
    } else {
      # The legacy builder emitted an empty element here (rendered as a
      # blank line before the structural block) when no intercepts exist.
      intercept_constraint <- ""
    }
    elems <- c(
      latent_var_str,
      error_constraint_str,
      intercept_constraint,
      "# structural model",
      paste(user_lines, collapse = "\n")
    )
  }
  paste0(elems, collapse = "\n")
}

tspa_sf <- function(model, data, se = NULL) {
  if (nrow(se) != 0) {
    return(tspa_render(tspa_schema_sf(model, se), style = "sf"))
  }
}

tspa_mf <- function(model, data, fsT, fsL, fsb) {
  tspa_render(tspa_schema_mf(model, fsT, fsL, fsb), style = "mf")
}

# ---------------------------------------------------------------------------
# Product-score (get_fs_int) auto-alias: the schema's generated model name
# for latent `v` is `fs_v`; a data column `fs_a:fs_b` (a,b latent names in
# se_fs) whose names concatenate to `v` is copied into the working data as
# `fs_v`. Old user models that pre-rename the product-score column keep
# working because the alias is a no-op when `fs_v` already exists.
# ---------------------------------------------------------------------------

tspa_sf_alias <- function(data, se) {
  is_lst <- inherits(data, "list") && !is.data.frame(data)
  dnames <- if (is_lst) names(data[[1]]) else names(data)
  se_names <- colnames(se)
  aliases <- character()
  for (v in se_names) {
    tgt <- paste0("fs_", v)
    if (tgt %in% dnames) next
    cand <- character()
    for (col in dnames) {
      pos <- regexpr(":", col, fixed = TRUE)
      if (pos < 1L) next
      a <- substr(col, 1L, pos - 1L)
      b <- substr(col, pos + 1L, nchar(col))
      # Product-score columns are `fs_a:fs_b` (both parts score names).
      if (!grepl("^fs_", a) || !grepl("^fs_", b)) next
      a <- sub("^fs_", "", a)
      b <- sub("^fs_", "", b)
      if (!(a %in% se_names) || !(b %in% se_names)) next
      if (paste0(a, b) == v || paste0(b, a) == v) cand <- c(cand, col)
    }
    if (length(cand) == 0) next
    if (length(cand) > 1) {
      stop(
        "Cannot determine which product-score column in the input data ",
        "corresponds to the latent variable '", v, "': ",
        paste0("\"", cand, "\"", collapse = " or "),
        ". Rename it to \"", tgt, "\" to disambiguate."
      )
    }
    if (is_lst) {
      for (i in seq_along(data)) data[[i]][[tgt]] <- data[[i]][[cand]]
    } else {
      data[[tgt]] <- data[[cand]]
    }
    aliases <- c(aliases, paste0(tgt, " <- ", cand))
  }
  list(data = data, aliases = aliases)
}
