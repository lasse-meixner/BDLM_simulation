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
  hahn <- BLR(y = data$Y, XF = cbind(data$A - fitted_ps$yHat), XR = data$X)
  linero <- BLR(y = data$Y, XF = cbind(data$A, fitted_ps$yHat), XR = data$X)

  # 2. FDML
  fitted_a_step1 <- fitted_ps
  fitted_y_step1 <- BLR(y = data$Y, XR = data$X)
     
  a_res_step2 <- data$A - fitted_a_step1$mu - data$X %*% fitted_a_step1$bR
  y_res_step2 <- data$Y - fitted_y_step1$mu - data$X %*% fitted_y_step1$bR

  fitted_dml_full <- lm(y_res_step2 ~ a_res_step2)

  n <- dim(data$X)[1] # number of observations
  ix_step1 <- 1:(floor(n/2))
  fitted_a_step1 <- BLR(y = data$A[ix_step1], XR = data$X[ix_step1,])
  fitted_y_step1 <- BLR(y = data$Y[ix_step1], XR = data$X[ix_step1,])

  ix_step2 <- (floor(n/2) + 1):n
  a_res_step2 <- data$A[ix_step2] - fitted_a_step1$mu - data$X[ix_step2,] %*% fitted_a_step1$bR
  y_res_step2 <- data$Y[ix_step2] - fitted_y_step1$mu - data$X[ix_step2,] %*% fitted_y_step1$bR
})) # end of silenced BLR calls

  fitted_dml_split <- lm(y_res_step2 ~ a_res_step2)
  # Ensure the return object is not printed
  invisible(list(naive = naive, hahn = hahn, linero = linero, FDML_full = fitted_dml_full, FDML_split = fitted_dml_split))
}

fit_mvn_iw_model <- function(data) {

  get_draw <- function(data){
    # Fit multivariate normal model with inverse Wishart prior for Sigma matrix
    iw_draw <- rmultireg(
      Y = cbind(data$Y, data$A), 
      X = data$X, 
      Bbar = cbind(rep(0, ncol(data$X)), rep(0, ncol(data$X))), # prior means
      A = diag(ncol(data$X))*0.2, # prior precision matrix
      nu = 3, # prior degrees of freedom
      V = matrix(4, 2, 2) # prior scale matrix
    )
    # return the fitted model
    Sigma_draw <- iw_draw$Sigma
    alpha  <- Sigma_draw[1, 2] / Sigma_draw[2, 2]
    # return draw for alpha
    alpha
  }

  # Get 1100 draws
  draws <- replicate(1100, get_draw(data))
  # Return the draws
  draws
  
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

# Main simulation function for a given setting ----
sim_iter_BLRs <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit_BRL <- fit_BLRs(data)
  fit_IW <- fit_mvn_iw_model(data)
  
  # extract results for each BRL model
  BRLs_extraction <- lapply(names(fit_BRL), function(model_name) {
    if (model_name %in% c("FDML_full", "FDML_split")) {
      extract_results_lm(fit_BRL[[model_name]], data$gamma, model_name, additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
    } else {
      extract_results_blr(fit_BRL[[model_name]], data$gamma, model_name, additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
    }
  })
  
  # extract results for IW model
  IW_extraction <- extract_results_IW(fit_IW, data$gamma, "BDML_IW", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))

  # combine results
  do.call(rbind, c(BRLs_extraction, list(IW_extraction)))
}
