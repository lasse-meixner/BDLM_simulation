library(tidyverse)
library(BLR)
library(bayesm) # for multivariate Inverse Wishart model

# Loading source for data generation
source("generate_data_source.R")

## Load all R scripts in the src/ directory ----
load_src_files <- function(src_path = "../src/") {
  source_files <- list.files(src_path, "*.R$")
  map(paste0(src_path, source_files), source)
}
load_src_files()

# Function to fit Hahn model ----
fit_BLRs <- function(data) {
  # Note: since BLR doesn't allow control of low-level cat statements (typical bs R package SE), I suppress this manually

  invisible(capture.output({ # Suppress output from BLR calls
    # 1. Hahn & Linero
    naive <- BLR(y = data$Y, XR = data$X, XF = as.matrix(data$A))
    fitted_ps <- BLR(y = data$A, XR = data$X)
    hcph <- BLR(y = data$Y, XF = cbind(data$A - fitted_ps$yHat), XR = data$X)
    linero <- BLR(y = data$Y, XF = cbind(data$A, fitted_ps$yHat), XR = data$X)

    # 2. FDML
    fitted_a_step1 <- fitted_ps
    fitted_y_step1 <- BLR(y = data$Y, XR = data$X)
      
    a_res_step2 <- data$A - fitted_a_step1$mu - data$X %*% fitted_a_step1$bR
    y_res_step2 <- data$Y - fitted_y_step1$mu - data$X %*% fitted_y_step1$bR

    fitted_fdml_full <- lm(y_res_step2 ~ a_res_step2)

    n <- dim(data$X)[1] # number of observations
    ix_step1 <- 1:(floor(n/2))
    fitted_a_step1 <- BLR(y = data$A[ix_step1], XR = data$X[ix_step1,])
    fitted_y_step1 <- BLR(y = data$Y[ix_step1], XR = data$X[ix_step1,])

    ix_step2 <- (floor(n/2) + 1):n
    a_res_step2 <- data$A[ix_step2] - fitted_a_step1$mu - data$X[ix_step2,] %*% fitted_a_step1$bR
    y_res_step2 <- data$Y[ix_step2] - fitted_y_step1$mu - data$X[ix_step2,] %*% fitted_y_step1$bR
  })) # end of silenced BLR calls

  fitted_fdml_split <- lm(y_res_step2 ~ a_res_step2)
  # Ensure the return object is not printed
  invisible(list(Naive = naive, HCPH = hcph, Linero = linero, "FDML-Full" = fitted_fdml_full, "FDML-Split" = fitted_fdml_split))
}

fit_mvn_iw_model <- function(data) {

  # reg_data holds data for each regression as separate list
  reg_data <- NULL
  reg_data[[1]] <- list(y = data$Y, X = data$X)
  reg_data[[2]] <- list(y = data$A, X = data$X)

  # Get 1000 draws
  invisible(capture.output({
    draws <- rsurGibbs(
    Data = list(regdata = reg_data),
    Prior = list(
      betabar = rep(0, ncol(data$X)*2), # prior mean (2*P)x1
      A = diag(1/ncol(data$X), ncol(data$X)*2), # prior precision (2*P)x(2*P)
      nu = 4, # IW prior degrees of freedom
      V = diag(1, 2, 2) # IW prior scale matrix 2x2
      ),
    Mcmc = list(R=12000, keep=1, nprint=0))
    }))
  # Return the transformed draws for alpha
  Sigma_draws <- draws$Sigmadraw[-(1:2000), , drop = FALSE] # R x (2*2), remove burn-in
  alpha_draws <- Sigma_draws[, 2] / Sigma_draws[, 4]
}

# Function to extract results for IW model
extract_results_IW <- function(draws, gamma, method_name, additional_results_info) {
  gamma_hat <- mean(draws)
  interval <- unname(quantile(draws, c(0.025, 0.975)))
  squared_error <- (gamma_hat - gamma)^2
  catch <- check_interval(gamma, interval)
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]
  
  table <- data.frame(
    gamma_hat = gamma_hat, 
    squared_error = squared_error,
    LCL = LCL,
    UCL = UCL,
    catch = catch,
    interval_width = interval_width,
    Method = method_name
  )

  # Append additional_results_info (setting parameters) to pass through for downstream analysis
  for (name in names(additional_results_info)) {
    table[[name]] <- additional_results_info[[name]]
  }
  table
}

# Function to extract results for Hahn and Linero models from BLR object ----
extract_results_blr <- function(fit, gamma, method_name, additional_results_info) {
  gamma_hat <- fit$bF[1]
  interval <- get_interval(fit)
  squared_error <- (gamma_hat - gamma)^2
  catch <- check_interval(gamma, interval)
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]
  
  table <- data.frame(
    gamma_hat = gamma_hat, 
    squared_error = squared_error,
    LCL = LCL,
    UCL = UCL,
    catch = catch,
    interval_width = interval_width,
    Method = method_name
  )

  # Append additional_results_info (setting parameters) to pass through for downstream analysis
  for (name in names(additional_results_info)) {
    table[[name]] <- additional_results_info[[name]]
  }
  table
}

# Function to extract results for FDML model from lm object ----
extract_results_lm <- function(fit, gamma, method_name, additional_results_info) {
  gamma_hat <- summary(fit)$coefficients[2, 1]
  interval <- unname(confint(fit)[2,])
  squared_error <- (gamma_hat - gamma)^2
  catch <- check_interval(gamma, interval)
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]
  
  table <- data.frame(
    gamma_hat = gamma_hat, 
    squared_error = squared_error,
    LCL = LCL,
    UCL = UCL,
    catch = catch,
    interval_width = interval_width,
    Method = method_name
  )

  # Append additional_results_info (setting parameters) to pass through for downstream analysis
  for (name in names(additional_results_info)) {
    table[[name]] <- additional_results_info[[name]]
  }
  table
}

# Main simulation function for BRLs for a given setting ----
sim_iter_BLRs <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit_BRL <- fit_BLRs(data)
  
  # extract results for each BRL model
  BRLs_extraction <- lapply(names(fit_BRL), function(model_name) {
    if (model_name %in% c("FDML-Full", "FDML-Split")) {
      extract_results_lm(fit_BRL[[model_name]], data$gamma, model_name, additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
    } else {
      extract_results_blr(fit_BRL[[model_name]], data$gamma, model_name, additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
    }
  })

  # combine results
  do.call(rbind, BRLs_extraction)
}

# Main simulation function for BDML-IW for a given setting ----
sim_iter_BDML_iw <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit_IW <- fit_mvn_iw_model(data)
  
  # extract results for IW model
  IW_extraction <- extract_results_IW(fit_IW, data$gamma, "BDML-IW", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))

  # combine results
  IW_extraction
}