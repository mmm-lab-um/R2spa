# Simple simulation

library(SimDesign)
library(mirt)  # for IRT
library(OpenMx)
library(umx)
library(dplyr)
library(R2spa)
library(ufs)
library(ggplot2)
library(devtools)
set.seed(52422)

# Write into a function

DESIGNFACTOR <- createDesign(
  N = c(250, 1000, 5000),
  beta = list(set1 = c(.2, .3, -.4),
              set2 = c(.15, .30, .22)),
  cov_xm = c(0.2, 0.5, 0.8)
)

GenData <- function (condition, fixed_objects = NULL) {
  num_obs <- condition$N # Sample size
  betas <- unname(unlist(condition$beta))
  cov_xm <- condition$cov_xm # latent correlation: varied

  # Loadings
  lambdax <- 1.7 * c(.6, .7, .8, .5, .5, .3, .9, .7, .6)  # on normal ogive metric
  lambdam <- 1.7 * c(.8, .7, .7, .6, .5)  # on normal ogive metric
  lambday <- c(.5, .7, .9)
  thresx <- matrix(c(1, 1.5, 1.3, 0.5, 2, -0.5, -1, 0.5, 0.8), nrow = 1)  # binary
  thresm <- matrix(c(-2, -2, -0.4, -0.5, 0,
                     -1.5, -0.5, 0, 0, 0.5,
                     -1, 1, 1.5, 0.8, 1), nrow = 3,
                   byrow = TRUE)  # ordinal with 4 categories
  inty <- 1.5  # intercept of y

  # Error variance of Y
  r2y <- crossprod(
    betas,
    matrix(c(1, cov_xm, 0, cov_xm, 1, 0, 0, 0, 1 + cov_xm^2), nrow = 3)
  ) %*% betas
  evar <- as.numeric(1 - r2y)

  #' Simulate latent variables X and M with variances of 1, and the disturbance of
  #' Y
  cov_xm_ey <- matrix(c(1, cov_xm, 0,
                        cov_xm, 1, 0,
                        0, 0, evar), nrow = 3)
  eta <- MASS::mvrnorm(num_obs, mu = rep(0, 3), Sigma = cov_xm_ey,
                       empirical = FALSE)
  # Add product term
  eta <- cbind(eta, eta[, 1] * eta[, 2])

  #' Compute latent y (standardized)
  etay <- inty + eta[, -3] %*% betas + eta[, 3]

  #' Simulate y (continuous)
  y <- t(
    lambday %*% t(etay) +
      rnorm(num_obs * length(lambday), sd = sqrt(1 - lambday^2))
  )

  #' Simulate latent continuous responses for X and M
  xstar <- eta[, 1] %*% t(lambdax) + rnorm(num_obs * length(lambdax))
  mstar <- eta[, 2] %*% t(lambdam) + rnorm(num_obs * length(lambdam))
  #' Obtain categorical items
  x <- vapply(
    seq_along(lambdax),
    FUN = \(i) {
      findInterval(xstar[, i], thresx[, i])
    },
    FUN.VALUE = numeric(num_obs))
  m <- vapply(
    seq_along(lambdam),
    FUN = \(i) {
      findInterval(mstar[, i], thresm[, i])
    },
    FUN.VALUE = numeric(num_obs))

  # Compute Factor Scores
  fsy <- get_fs(data = y, std.lv = TRUE)

  irtx <- mirt(data.frame(x), itemtype = "2PL",
               verbose = FALSE)  # IRT (2-PL) for x
  fsx <- fscores(irtx, full.scores.SE = TRUE)

  irtm <- mirt(data.frame(m), itemtype = "graded",
               verbose = FALSE)  # IRT (GRM) for m
  fsm <- fscores(irtm, full.scores.SE = TRUE)

  # Assemble data
  fs_dat <- data.frame(fsy[c(1, 3, 4)], fsx, fsm) |>
    setNames(c(
      "fsy", "rel_fsy", "ev_fsy",
      "fsx", "se_fsx", "fsm", "se_fsm"
    )) |>
    # Compute reliability; only needed for 2S-PA
    within(expr = {
      rel_fsx <- 1 - se_fsx^2
      rel_fsm <- 1 - se_fsm^2
      ev_fsx <- se_fsx^2 * (1 - se_fsx^2)
      ev_fsm <- se_fsm^2 * (1 - se_fsm^2)
    }) |>
    # Add interaction
    within(expr = {
      fsxm <- fsx * fsm
      ld_fsxm <- rel_fsx * rel_fsm
      ev_fsxm <- rel_fsx^2 * ev_fsm + rel_fsm^2 * ev_fsx + ev_fsm * ev_fsx
    })

  return(fs_dat)
}

# Test it
GenData(condition = DESIGNFACTOR[1,])

extract_res <- function (condition, dat, fixed_objects = NULL) {

  fs_dat = dat

  # Using factor score
  fs_summary <- summary(lm(fsy ~ fsx * fsm, data = fs_dat))

  # Openmx
  fsreg_mx <- mxModel("str",
                      type = "RAM",
                      manifestVars = c("fsy", "fsx", "fsm", "fsxm"),
                      # Path coefficients
                      mxPath(
                        from = c("fsx", "fsm", "fsxm"), to = "fsy",
                        values = 0
                      ),
                      # Variances
                      mxPath(
                        from = c("fsy", "fsx", "fsm", "fsxm"),
                        to = c("fsy", "fsx", "fsm", "fsxm"),
                        arrows = 2,
                        values = c(0.6, 1, 1, 1)
                      ),
                      # Covariances
                      mxPath(
                        from = c("fsx", "fsm", "fsxm"),
                        to = c("fsx", "fsm", "fsxm"),
                        arrows = 2, connect = "unique.pairs",
                        values = 0
                      ),
                      mxPath(
                        from = "one", to = c("fsy", "fsx", "fsm", "fsxm"),
                        values = 0
                      )
  )

  fsreg_umx <- umxLav2RAM(
    "
      fsy ~ fsx + fsm + fsxm
      fsy + fsx + fsm + fsxm ~ 1
      fsx ~~ fsm + fsxm
      fsm ~~ fsxm
    ",
    printTab = FALSE)

  # Loading
  matL <- mxMatrix(
    type = "Diag", nrow = 4, ncol = 4,
    free = FALSE,
    labels = c("data.rel_fsy", "data.rel_fsx", "data.rel_fsm", "data.ld_fsxm"),
    name = "L"
  )
  # Error
  matE <- mxMatrix(
    type = "Diag", nrow = 4, ncol = 4,
    free = FALSE,
    labels = c("data.ev_fsm", "data.ev_fsx", "data.ev_fsm", "data.ev_fsxm"),
    name = "E"
  )

  tspa_mx <-
    R2spa:::tspa_mx_model(fsreg_umx, data = fs_dat,
                          mat_ld = matL, mat_vc = matE)
  tspa_mx_fit <- mxRun(tspa_mx)
  tspa_summary <- mxStandardizeRAMpaths(tspa_mx_fit, SE = TRUE)$m1

  tspa_pe <- round(tspa_summary[tspa_summary$label %in% c("fsx_to_fsy", "fsm_to_fsy", "fsxm_to_fsy"), ]$Raw.Value, 3)
  names(tspa_pe) <- names(round(fs_summary$coefficients[-1,"Estimate"], 3))

  tspa_se <- round(tspa_summary[tspa_summary$label %in% c("fsx_to_fsy", "fsm_to_fsy", "fsxm_to_fsy"), ]$Raw.SE, 3)
  names(tspa_se) <- names(tspa_pe)

  # Extract parameter estimates and standard errors
  paret <- c(fscore = round(fs_summary$coefficients[-1,"Estimate"], 3),
             fscore_se = round(fs_summary$coefficients[-1,"Std. Error"], 3),
             tspa_pe = tspa_pe,
             tspa_se = tspa_se)
  names(paret) <- c("fsx_est", "fsm_est", "fsxm_est",
                    "fsx_se", "fsm_se", "fsxm_se",
                    "tspax_est", "tspam_est", "tspaxm_est",
                    "tspax_se", "tspam_se", "tspaxm_se")
  return(paret)
}

# Simple simulation of raw data frames in a list
results_df <- list()
for (i in seq_len(nrow(DESIGNFACTOR))) {
  temp_result <- c()
  for (j in 1:2) {
    temp_row <- as.data.frame(t(unlist(extract_res(dat = GenData(DESIGNFACTOR[i,])))))
    temp_result <- rbind(temp_result, temp_row)
  }
  results_df[[i]] <- temp_result
  names(results_df)[i] <- paste0("Condition_", i)
}

saveRDS(results_df, "cat_dat_0902.RDS")
