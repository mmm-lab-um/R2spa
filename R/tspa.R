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
#' fs_dat_ind60 <- get_fs(data = PoliticalDemocracy,
#'                        model = "ind60 =~ x1 + x2 + x3")
#' fs_dat_dem60 <- get_fs(data = PoliticalDemocracy,
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
#' fs_dat_visual <- get_fs(data = HolzingerSwineford1939,
#'                         model = "visual =~ x1 + x2 + x3",
#'                         group = "school",
#'                         format = "list")
#' fs_dat_speed <- get_fs(data = HolzingerSwineford1939,
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

tspa_sf <- function(model, data, se = NULL) {
  if (nrow(se) != 0) {
    ev <- se^2
    var <- names(se)
    len <- length(se)
    group <- nrow(se)
    # col <- colnames(data)
    fs <- paste0("fs_", var)
    tspaModel <- rep(NA, len)
    latent_var <- rep(NA, len)
    error_constraint <- rep(NA, len)

    latent_var <- lapply(seq_len(len), function(x) {
      paste0(
        var[x], "=~ c(",
        paste0(rep(1, group), collapse = ", "), ") * ", fs[x], "\n"
      )
    })

    error_constraint <- lapply(seq_len(len), function(x) {
      paste0(
        fs[x], "~~ c(", paste(ev[[x]], collapse = ", "), ") * ", fs[x], "\n"
      )
    })

    latent_var_str <- paste(latent_var, collapse = "")
    error_constraint_str <- paste(error_constraint, collapse = "")
    tspaModel <- paste0(
      c("# latent variables (indicated by factor scores)",
        latent_var_str,
        "# constrain the errors",
        error_constraint_str,
        "# structural model",
        model),
      collapse = "\n"
    )

    return(tspaModel)
  }
}

tspa_mf <- function(model, data, fsT, fsL, fsb) {
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
    fsL1 <- fsL[[1]]
    fsT_in <- !upper.tri(fsT[[1]])
  } else if (is.list(fsL) && length(fsL) > 1) {
    if (!is.list(fsT) || length(fsT) != length(fsL)) {
      stop("'fsT' must be a list of the same length as 'fsL' for a ",
           "multigroup model.")
    }
    ngroup <- length(fsL)
    fsL1 <- fsL[[1]]
    fsT_in <- !upper.tri(fsT[[1]])
  } else {
    ngroup <- 1
    fsL1 <- if (is.list(fsL)) fsL[[1]] else fsL
    fsT_in <- !upper.tri(if (is.list(fsT)) fsT[[1]] else fsT)
  }
  var <- colnames(fsL1)
  nvar <- length(var)
  fs <- rownames(fsL1)

  # latent variables
  loadings_mat <- matrix(unlist(fsL), ncol = ngroup)
  loadings <- apply(loadings_mat, 1, function(x) {
    paste0("c(", paste0(x, collapse = ", "), ") * ")
  }) |>
    paste0(fs)
  loadings_list <- split(loadings, factor(rep(var, each = nvar),
                                          levels = var))
  loadings_c <- lapply(loadings_list, function(x) {
    paste0(x, collapse = " + ")
  })
  latent_var_str <- paste("# latent variables (indicated by factor scores)\n",
                          var, "=~", loadings_c)
  # error variances
  ev_rhs <- fs[col(fsT_in)[fsT_in]]
  ev_lhs <- fs[row(fsT_in)[fsT_in]]
  errors_mat <- matrix(unlist(fsT), ncol = ngroup)[as.vector(fsT_in), ,
                                                   drop = FALSE]
  errors <- apply(errors_mat, 1, function(x) {
    paste0("c(", paste0(x, collapse = ", "), ")")
  })
  error_constraint_str <- paste0("# constrain the errors\n",
                                 ev_lhs, " ~~ ", errors, " * ", ev_rhs)
  if (!is.null(fsb)) {
    # intercepts
    intercepts_mat <- matrix(unlist(fsb), ncol = ngroup)
    intercepts <- split(intercepts_mat, rep(seq_len(nrow(intercepts_mat)), ngroup))
    intercept_constraint <- paste0("# constrain the intercepts\n",
                                   fs, " ~ ", intercepts, " * 1")
  } else {
    intercept_constraint <- ""
  }

  tspaModel <- paste0(c(
    latent_var_str,
    error_constraint_str,
    intercept_constraint,
    "# structural model",
    model
  ),
  collapse = "\n")

  return(tspaModel)
}
