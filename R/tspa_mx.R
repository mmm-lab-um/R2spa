#' Two-Stage Path Analysis with Definition Variables Using OpenMx
#'
#' Fit a two-stage path analysis (2S-PA) model.
#'
#' @param mx_model A model object of class [OpenMx::MxRAMModel-class], created
#'                 by [OpenMx::mxModel()] or the `umx` package. This
#'                 is the structural model part.
#' @param data A data frame containing factor scores.
#' @param mat_ld A \eqn{p \times p} matrix indicating the loadings of the
#'               factor scores on the latent variables. The *i*th row indicate
#'               loadings of the *i*th factor score variable on the latent
#'               variables. This could be one of the following:
#' * A matrix created by [OpenMx::mxMatrix()] with loading
#'   values.
#' * A named numeric matrix, with rownames and column names
#'   matching the factor score and latent variables.
#' * A named character matrix, with rownames and column names
#'   matching the factor score and latent variables, and the
#'   character values indicating the variable names in `data` for
#'   the corresponding loadings.
#' @param mat_ev Similar to `mat_ld` but for the error variance-covariance
#'               matrix of the factor scores.
#' @param mat_int Similar to `mat_ld` but for the measurement intercept
#'               matrix of the factor scores.
#' @param fs_lv_names A named character vector where each element is
#'                    the name of a factor score variable, and the
#'                    names of the vector are the corresponding name
#'                    of the latent variables.
#' @param ... Additional arguments passed to [OpenMx::mxModel()].
#' @return An object of class [OpenMx::MxModel-class]. Note that the
#'         model has not been run.
#'
#' @importFrom OpenMx mxData mxMatrix mxAlgebraFromString mxExpectationNormal
#'             mxModel mxFitFunctionML
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(mirt)
#' library(umx)
#' library(OpenMx)
#' # Simulate data with mirt
#' set.seed(1324)
#' num_obs <- 100
#' # Simulate theta
#' eta <- MASS::mvrnorm(num_obs, mu = c(0, 0), Sigma = diag(c(1, 1 - 0.5^2)),
#'                      empirical = TRUE)
#' th1 <- eta[, 1]
#' th2 <- -1 + 0.5 * th1 + eta[, 2]
#' # items and response data
#' a1 <- matrix(1, 10)
#' d1 <- matrix(rnorm(10))
#' a2 <- matrix(runif(10, min = 0.5, max = 1.5))
#' d2 <- matrix(rnorm(10))
#' dat1 <- mirt::simdata(a = a1, d = d1,
#'                       N = num_obs, itemtype = "2PL", Theta = th1)
#' dat2 <- mirt::simdata(a = a2, d = d2, N = num_obs,
#'                       itemtype = "2PL", Theta = th2)
#' # Factor scores
#' mod1 <- mirt(dat1, model = 1, itemtype = "Rasch", verbose = FALSE)
#' mod2 <- mirt(dat2, model = 1, itemtype = "2PL", verbose = FALSE)
#' fs1 <- fscores(mod1, full.scores.SE = TRUE)
#' fs2 <- fscores(mod2, full.scores.SE = TRUE)
#' # Combine factor scores and standard errors into data set
#' fs_dat <- as.data.frame(cbind(fs1, fs2))
#' names(fs_dat) <- c("fs1", "se_fs1", "fs2", "se_fs2")
#' # Compute reliability and error variances
#' fs_dat <- within(fs_dat, expr = {
#'   rel_fs1 <- 1 - se_fs1^2
#'   rel_fs2 <- 1 - se_fs2^2
#'   ev_fs1 <- se_fs1^2 * (1 - se_fs1^2)
#'   ev_fs2 <- se_fs2^2 * (1 - se_fs2^2)
#' })
#' # OpenMx model (from umx so that lavaan syntax can be used)
#' fsreg_umx <- umxLav2RAM(
#'   "
#'     fs2 ~ fs1
#'     fs2 + fs1 ~ 1
#'   ",
#'   printTab = FALSE)
#' # Prepare loading and error covariance matrices
#' cross_load <- matrix(c("rel_fs2", NA, NA, "rel_fs1"), nrow = 2) |>
#'   `dimnames<-`(rep(list(c("fs2", "fs1")), 2))
#' err_cov <- matrix(c("ev_fs2", NA, NA, "ev_fs1"), nrow = 2) |>
#'   `dimnames<-`(rep(list(c("fs2", "fs1")), 2))
#' # Create 2S-PA model (with definition variables)
#' tspa_mx <-
#'   tspa_mx_model(fsreg_umx,
#'     data = fs_dat,
#'     mat_ld = cross_load,
#'     mat_ev = err_cov
#'   )
#' # Run OpenMx
#' tspa_mx_fit <- mxRun(tspa_mx)
#' # Summarize the results
#' summary(tspa_mx_fit)
#' }

tspa_mx_model <- function(mx_model, data, mat_ld, mat_ev,
                          mat_int = NULL,
                          fs_lv_names = NULL, ...) {
  str_name <- mx_model$name
  q <- ncol(mx_model$A$values)
  if (!isTRUE(attr(class(mat_ld), which = "package") == "OpenMx")) {
    mv_name <- mx_model$manifestVars
    ld_ind <- match(mv_name, table = colnames(mat_ld))
    if (any(is.na(ld_ind))) {
      stop("`colnames(mat_ld)` must match the names of the ",
           "latent factor variables.")
    }
    if (is.null(fs_lv_names) &&
      !identical(rownames(mat_ld), colnames(mat_ld))) {
      fs_lv_names <- rownames(mat_ld)
      names(fs_lv_names) <- colnames(mat_ld)
    }
    mat_ld <- make_mx_ld(mat_ld[ld_ind, ld_ind])
  }
  if (!is.null(fs_lv_names)) {
    dup_fs <- data[fs_lv_names]
    names(dup_fs) <- names(fs_lv_names)
    data <- cbind(data, dup_fs)
  }
  if (!isTRUE(attr(class(mat_ev), which = "package") == "OpenMx")) {
    mv_name <- mx_model$manifestVars
    if (!is.null(fs_lv_names)) {
      mv_name <- fs_lv_names[match(mv_name, table = names(fs_lv_names))]
    }
    ev_ind <- match(mv_name, table = rownames(mat_ev))
    if (any(is.na(ev_ind))) {
      stop("`rownames(mat_ev)` must match the names of the ",
           "factor score variables.")
    }
    mat_ev <- make_mx_vc(mat_ev[ev_ind, ev_ind])
  }
  if (!is.null(mat_int)) {
    if (!isTRUE(attr(class(mat_int), which = "package") == "OpenMx")) {
      mv_name <- mx_model$manifestVars
      if (!is.null(fs_lv_names)) {
        mv_name <- fs_lv_names[match(mv_name, table = names(fs_lv_names))]
      }
      int_ind <- match(mv_name, table = colnames(mat_int))
      if (any(is.na(int_ind))) {
        stop("`colnames(mat_int)` must match the names of the ",
             "factor score variables.")
      }
      mat_int <- make_mx_int(mat_int[, int_ind, drop = FALSE])
    }
  } else {
    zero_row_vector <- matrix(0, ncol = length(mx_model$manifestVars))
    dimnames(zero_row_vector) <- list(NULL, mx_model$manifestVars)
    mat_int <- make_mx_int(zero_row_vector)
  }
  mxModel(
    "2SPAD",
    mxData(observed = data, type = "raw"),
    mx_model, mat_ld, mat_ev, mat_int,
    mxMatrix(type = "Iden", nrow = q, ncol = q, name = "I"),
    mxAlgebraFromString(paste0("(L %*%", str_name,
                               ".F %*% solve(I - ", str_name, ".A)) %&% ",
                               str_name, ".S + E"), name = "expCov"),
    mxAlgebraFromString(paste0(str_name, ".M %*% t(L %*%",
                               str_name, ".F %*% solve(I - ",
                               str_name, ".A)) + b"), name = "expMean"),
    # mxAlgebraFromString(paste0(str_name, ".M"), name = "expMean"),
    mxExpectationNormal(
      covariance = "expCov", means = "expMean",
      dimnames = mx_model$manifestVars),
    mxFitFunctionML()
  )
}

make_mx_ld <- function(ld_mat) {
  if (is.numeric(ld_mat)) {
    val <- ld_mat
    lab <- NA
  } else if (is.character(ld_mat)) {
    lab <- ld_mat
    lab[!is.na(lab)] <- paste0("data.", lab[!is.na(lab)])
    val <- NA
  } else {
    stop("`ld_mat` must be either a numeric matrix or a character matrix.")
  }
  mxMatrix(
    type = "Full", nrow = nrow(ld_mat), ncol = nrow(ld_mat),
    free = FALSE,
    values = val,
    labels = lab,
    name = "L"
  )
}

make_mx_vc <- function(vc_mat) {
  if (is.numeric(vc_mat)) {
    val <- vc_mat
    lab <- NA
  } else if (is.character(vc_mat)) {
    lab <- vc_mat
    lab[!is.na(lab)] <- paste0("data.", lab[!is.na(lab)])
    val <- NA
  } else {
    stop("`vc_mat` must be either a numeric matrix or a character matrix.")
  }
  mxMatrix(
    type = "Symm", nrow = nrow(vc_mat), ncol = nrow(vc_mat),
    free = FALSE,
    values = val,
    labels = lab,
    name = "E"
  )
}

make_mx_int <- function(int_mat) {
  if (!is.matrix(int_mat) || nrow(int_mat) != 1) {
    stop("`int_mat` must be a 1xN matrix (i.e., a row vector).")
  }
  if (is.numeric(int_mat)) {
    val <- int_mat
    lab <- NA
  } else if (is.character(int_mat)) {
    lab <- int_mat
    lab[!is.na(lab)] <- paste0("data.", lab[!is.na(lab)])
    val <- NA
  } else {
    stop("`int_mat` must be either a numeric matrix or a character matrix.")
  }
  mxMatrix(
    type = "Full", nrow = 1, ncol = length(int_mat),
    free = FALSE,
    values = val,
    labels = lab,
    name = "b"
  )
}