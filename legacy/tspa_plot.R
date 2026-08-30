#' Diagnostic plots of fitted 2S-PA model
#'
#' @param tspa_fit An object of class [lavaan][lavaan:lavaan-class],
#'                 representing the output generated from the [tspa()]
#'                 function.
#' @param title Character. Set the name of scatter plot. The default value is
#'              "Scatterplot".
#' @param label_x Character. Set the name of the x-axis. The default value is
#'                "fs_" followed by variable names.
#' @param label_y Character. Set the name of the y-axis. The default value is
#'                "fs_" followed by variable names.
#' @param abbreviation Logic input. If `FALSE` is indicated, the group name
#'                     will be shown in full. The default setting is `TRUE`.
#' @param ask Logic input. If `TRUE` is indicated, the user will be asked before
#'            before each plot is generated. The default setting is 'False'.
#' @param fscores_type Character. Set the type of factor score for input.
#'                     The default setting is using factor score from observed
#'                     data (i.e., output from [get_fs()]). If
#'                     `fscore_type = "est"`, then it will use output from
#'                     [lavaan::lavPredict()].
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}. See
#'            \code{\link[graphics]{plot}} for a list.
#'
#' @return A scatterplot between factor scores, and a residual plot.
#'
#' @importFrom grDevices devAskNewPage
#'
#' @export
#'
#' @examples
#' library(lavaan)
#' model <- "
#' # latent variable definitions
#' ind60 =~ x1 + x2 + x3
#' dem60 =~ y1 + a*y2 + b*y3 + c*y4
#' # regressions
#' dem60 ~ ind60
#' "
#' fs_dat_ind60 <- get_fs(
#'   data = PoliticalDemocracy,
#'   model = "ind60 =~ x1 + x2 + x3"
#' )
#' fs_dat_dem60 <- get_fs(
#'   data = PoliticalDemocracy,
#'   model = "dem60 =~ y1 + y2 + y3 + y4"
#' )
#' fs_dat <- cbind(fs_dat_ind60, fs_dat_dem60)
#'
#' tspa_fit <- tspa(
#'   model = "dem60 ~ ind60",
#'   data = fs_dat,
#'   se_fs = list(ind60 = 0.1213615, dem60 = 0.6756472)
#' )
#' tspa_plot(tspa_fit)
tspa_plot <- function(tspa_fit,
                      title = NULL,
                      label_x = NULL,
                      label_y = NULL,
                      abbreviation = TRUE,
                      fscores_type = c("original", "lavaan"),
                      ask = FALSE,
                      ...) {
  fit_pars <- lavaan::parameterestimates(tspa_fit) # parameter estimates
  # HL Note: this may not be the most reliable in extracting exo and
  # endo variables.
  regression_fit <- fit_pars[which(fit_pars$op == "~"), ] # path coefficients
  # latent_dv <- unique(c(unname(t(regression_fit["lhs"])))) # Extract DV
  # latent_iv <- unique(c(unname(t(regression_fit["rhs"])))) # Extract IV
  fscores_type <- match.arg(fscores_type) # Specify the type of factor scores

  # reg_pairs <- names(coef(tspa_fit)[names(coef(tspa_fit)) %in%
  #   apply(expand.grid(latent_dv, latent_iv), 1, paste, collapse = "~")])
  # HL: Perhaps an easier way
  reg_pairs <- paste(regression_fit[["lhs"]], regression_fit[["rhs"]],
                     sep = "~")

  # Type of scatterplot
  if (fscores_type == "original") {
    # factor scores from `get_fs`
    fscores <- lavaan::lavInspect(tspa_fit, what = "data")
  } else {
    # factor scores from estimation using `lavPredict`
    fscores <- lavaan::lavPredict(tspa_fit)
  }

  # Determine multi-group sample
  # HL: To simplify the code, just convert fscores to a list for single groups,
  #     and the for loop can apply equally.
  if (!is.list(fscores)) {
    fscores <- list(fscores)
    g_nums <- list(NULL)
  } else {
    g_nums <- as.list(seq_along(fscores))
  }
  # Ask for abbreviation
  g_names <- names(fscores)
  if (!is.null(g_names) && abbreviation) {
    g_names <- abbreviate(g_names)
  }

  # Prepare for abline
  fscores <- lapply(fscores, as.data.frame)
  # slope <- coef(tspa_fit)[grepl(apply(expand.grid(latent_dv, latent_iv), 1, paste, collapse="~"),
  #                               names(coef(tspa_fit)))]
  # fs_means <- matrix(unlist(lapply(df_latent_scores, apply, 2, mean, na.rm = T)),
  #                    nrow = length(slope),
  #                    byrow = T)
  # intercept <- fs_means[,2] - fs_means[,1]*slope

  # Ask for next plot
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  for (g in seq_along(g_nums)) {
    for (i in seq_along(reg_pairs)) {

      # Scatterplot
      plot_scatter(
        # Note: Use `[[` for list
        fscores_df = fscores[[g]],
        iv = unlist(strsplit(reg_pairs[i], split = "~"))[2],
        dv = unlist(strsplit(reg_pairs[i], split = "~"))[1],
        # ab_slope = slope[g],
        # ab_intercept = intercept[g],
        g_num = g_nums[[g]],
        g_name = g_names[g],
        title = title[i],
        label_x = label_x[i],
        label_y = label_y[i],
        ...
      )

      # Residual Plot
      plot_residual(
        fscores_df = fscores[[g]],
        iv = unlist(strsplit(reg_pairs[i], split = "~"))[2],
        dv = unlist(strsplit(reg_pairs[i], split = "~"))[1],
        tspa_fit = tspa_fit,
        g_num = g_nums[[g]],
        g_name = g_names[g],
        ...
      )
    }
  }
}

# Helper function for the scatter plot
plot_scatter <- function(fscores_df, iv, dv,
                         # ab_slope, ab_intercept,
                         g_num, g_name,
                         title, label_x, label_y,
                         ...) {
  # Group name for multiple-group case
  # HL: The function exists to handle single-group. Make input
  #     consistent to avoid extra conditional statements.
  iv_data <- fscores_df[, paste0("fs_", iv)]
  dv_data <- fscores_df[, paste0("fs_", dv)]

  # ylab
  ylab <- c()
  if (is.null(label_y)) {
    ylab <- paste0("fs_", dv)
  } else {
    ylab <- label_y
  }

  # ylab
  xlab <- c()
  if (is.null(label_x)) {
    xlab <- paste0("fs_", iv)
  } else {
    xlab <- label_x
  }

  # title
  # HL: Simplify to two conditionals
  if (is.null(title)) {
    title <- paste0("Scatterplot")
  }
  if (!is.null(g_name)) {
    title <- paste0(title, " (Group ", g_num, ": ", g_name, ")")
  }

  plot(iv_data,
    dv_data,
    ylab = ylab,
    xlab = xlab,
    main = title,
    pch = 16,
    ...
  )
  # abline(a = ab_intercept, b = ab_slope, lwd = 2)
}

# Helper function for the residual plot
# HL: Incorrect to say `g_num = g_num`, as one cannot assume there
#     `g_num` exists in the parent environment.
plot_residual <- function(fscores_df, iv, dv,
                          tspa_fit,
                          g_num, g_name,
                          ...) {
  # HL: The function exists to handle single-group. Make input
  #     consistent to avoid extra conditional statements.
  iv_data <- fscores_df[, paste0("fs_", iv)]
  dv_data <- fscores_df[, paste0("fs_", dv)]

  predicted_y <- lavaan::lavPredictY(tspa_fit,
    ynames = paste0("fs_", dv), xnames = paste0("fs_", iv),
    assemble = FALSE
  )

  if (!is.null(g_name)) {
    predicted_y <- predicted_y[[g_num]]
    # title
    title <- paste0("Residual Plot", " (Group ", g_num, ": ", g_name, ")")
  } else {
    title <- paste0("Residual Plot")
  }

  resid <- dv_data - predicted_y

  plot(iv_data,
    resid,
    ylab = paste0("Residuals: ", dv),
    xlab = paste0("Fitted values: ", iv),
    main = title,
    pch = 18,
    ...
  )
}
