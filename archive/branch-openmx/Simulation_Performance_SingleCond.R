# Simple simulation

library(SimDesign)
library(mirt)  # for IRT
library(OpenMx)
library(umx)
library(R2spa)
library(ufs)
library(ggplot2)
library(tidyverse)
library(dplyr)
devtools::load_all()

# Write into a function

DESIGNFACTOR <- createDesign(
  N = 5000,
  beta = list(set1 = c(.2, .3, -.4),
              set2 = c(.15, .30, .22)),
  cov_xm = 0.8
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
# GenData(condition = DESIGNFACTOR[1,])

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

evaluate_res <- function (condition, results, fixed_objects = NULL) {

  # Population parameter
  pop_par <- rep(unlist(condition$beta), 2)
  names(pop_par) <- c("fsx", "fsm", "fsxm",
                      "tspax", "tspam", "tspaxm")

  # Separate estimates and se
  results_est <- as.data.frame(results[names(results)[grepl("_est", names(results))]])
  results_se <- as.data.frame(results[names(results)[grepl("_se", names(results))]])

  # Helper function: relative SE bias
  rse_bias <- function(est_se, est) {
    est_se <- as.matrix(est_se)
    est <- as.matrix(est)
    est_se <- colMeans(est_se)
    emp_sd <- apply(est, 2L, sd)
    est_se / emp_sd - 1
  }

  # Helper functuion for calculating coverage rate
  coverage_rate <- function(est_se, est, pop) {
    ci_est <- list()
    est_se <- as.matrix(est_se)
    est <- as.matrix(est)
    lo.95 <- est - qnorm(.975)*est_se
    hi.95 <- est + qnorm(.975)*est_se

    for (i in seq_len(length(colnames(est)))) {
      ci_est[[i]] <- cbind(lo.95[,i], hi.95[,i])
      names(ci_est)[i] <- colnames(est)[i]
    }

    temp <- rep(NA, length(pop_par))

    for (i in seq_len(length(ci_est))) {
      temp[i] <- ECR(ci_est[[i]], parameter = pop_par[i])
    }

    return(temp)

  }

  c(
    std_bias = bias(results_est,
                    parameter = pop_par,
                    type = "standardized"),
    coverage = coverage_rate(results_se,
                             results_est,
                             pop = pop_par),
    rmse = RMSE(results_est,
                parameter = pop_par),
    rse_bias = rse_bias(results_se,
                        results_est)
  )
}

sim_cat_singlecond <- runSimulation(design = DESIGNFACTOR,
                              replications = 2000,
                              generate = GenData,
                              analyse = extract_res,
                              summarise = evaluate_res,
                              save = TRUE,
                              save_results = TRUE,
                              filename = "sim_cat_singlecond",
                              parallel = TRUE,
                              ncores = min(4L, parallel::detectCores() - 1))

# Summarize the results

sim_temp <- readRDS("sim_cat_singlecond.rds") %>%
  rename("COVERAGE.fs_x_est" = "coverage1",
         "COVERAGE.fs_m_est" = "coverage2",
         "COVERAGE.fs_xm_est" = "coverage3",
         "COVERAGE.tspa_x_est" = "coverage4",
         "COVERAGE.tspa_m_est" = "coverage5",
         "COVERAGE.tspa_xm_est" = "coverage6",
         "RSEbias.fs_x_se" = "rse_bias.fsx_se",
         "RSEbias.fs_m_se" = "rse_bias.fsm_se",
         "RSEbias.fs_xm_se" = "rse_bias.fsxm_se",
         "RSEbias.tspa_x_se" = "rse_bias.tspax_se",
         "RSEbias.tspa_m_se" = "rse_bias.tspam_se",
         "RSEbias.tspa_xm_se" = "rse_bias.tspaxm_se",
         "STDbias.fs_x_est" = "std_bias.fsx_est",
         "STDbias.fs_m_est" = "std_bias.fsm_est",
         "STDbias.fs_xm_est" = "std_bias.fsxm_est",
         "STDbias.tspa_x_est" = "std_bias.tspax_est",
         "STDbias.tspa_m_est" = "std_bias.tspam_est",
         "STDbias.tspa_xm_est" = "std_bias.tspaxm_est",
         "RMSE.fs_x_est" = "rmse.fsx_est",
         "RMSE.fs_m_est" = "rmse.fsm_est",
         "RMSE.fs_xm_est" = "rmse.fsxm_est",
         "RMSE.tspa_x_est" = "rmse.tspax_est",
         "RMSE.tspa_m_est" = "rmse.tspam_est",
         "RMSE.tspa_xm_est" = "rmse.tspaxm_est") %>%
  mutate(beta = as.character(rep(c(rep(0.2, 1), rep(0.15, 1)), 1)))

sim_results <- sim_temp %>%
  gather("var", "val", STDbias.fs_x_est:RSEbias.tspa_xm_se) %>%
  select(-c(SIM_TIME:WARNINGS)) %>%
  separate(col = var, into = c("stats", "parmet"), sep = "\\.") %>%
  separate(col = parmet, into = c("method", "par", "result"),  sep = "_") %>%
  select(-result) %>%
  spread(stats, val) %>%
  relocate(REPLICATIONS, .after = last_col()) %>%
  mutate(N_lab = as_factor(paste0("italic(N) == ", N)),
         beta_lab = as_factor(paste0("beta == ", beta)),
         cov_xm_lab = as_factor(paste0("Correlation_XM == ", cov_xm)))

write_csv(sim_results, "sim_cat_singlecond.csv")

# Plot the results

sim_plots <- read.csv("sim_cat_singlecond.csv")
# # Bias
# sim_plots %>%
#   ggplot(aes(x = factor(N), y = bias, color = method)) +
#   geom_boxplot() +
#   facet_grid(cor_xm_lab ~ rel_lab, labeller = label_parsed) +
#   labs(x = "Sample Size (N)", y = "Bias")

# Standard Bias
sim_plots %>%
  ggplot(aes(x = factor(N), y = STDbias, color = method)) +
  geom_boxplot() +
  facet_grid(beta_lab ~ cov_xm_lab, labeller = label_parsed) +
  labs(x = "Sample Size (N)", y = "Standardized Bias")

# Relative SE Bias
sim_plots %>%
  ggplot(aes(x = factor(N), y = RSEbias, color = method)) +
  geom_boxplot() +
  facet_grid(beta_lab ~ cov_xm_lab, labeller = label_parsed) +
  labs(x = "Sample Size (N)", y = "Relative SE Bias")

# Coverage rate
sim_plots %>%
  ggplot(aes(x = factor(N), y = COVERAGE, color = method)) +
  geom_boxplot() +
  facet_grid(beta_lab ~ cov_xm_lab, labeller = label_parsed) +
  labs(x = "Sample Size (N)", y = "Coverage Rate (95%)")

# RMSE
sim_plots %>%
  ggplot(aes(x = factor(N), y = RMSE, color = method)) +
  geom_boxplot() +
  facet_grid(beta_lab ~ cov_xm_lab, labeller = label_parsed) +
  labs(x = "Sample Size (N)", y = "RMSE")



