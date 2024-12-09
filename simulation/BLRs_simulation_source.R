library(tidyverse)
library(BLR)

## requires set working directory to "programs/Linero2023/BDLM_clean"

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
  # first stage
  fitted_ps <- BLR(y = data$A, XR = data$X)
  # naive 
  naive <- BLR(y = data$Y, XR = data$X, XF = as.matrix(data$A))
  # Hahn second stage
  hahn <- BLR(y = data$Y, XF = cbind(data$A - fitted_ps$yHat), XR = data$X)
  # Linero second stage
  linero <- BLR(y = data$Y, XF = cbind(data$A, fitted_ps$yHat), XR = data$X)

  list(naive = naive, hahn = hahn, linero = linero)
}

# Function to extract results for Hahn and Linero models ----
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

# Main simulation function for a given setting ----
sim_iter_BLRs <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit <- fit_BLRs(data)
  
  # extract results for each model
  extraction <- lapply(names(fit), function(model_name) {
    extract_results_blr(fit[[model_name]], data$gamma, model_name, additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  })

  do.call(rbind, extraction)
}
