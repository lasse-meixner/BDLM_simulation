## Load necessary libraries ----
library(cmdstanr)
library(tidyverse)

## source data generation
source("generate_data_source.R")

## Load all R scripts in the src/ directory ----
load_src_files <- function(src_path = "../src/") {
  source_files <- list.files(src_path, "*.R$")
  map(paste0(src_path, source_files), source)
}
load_src_files()

## Compile the Stan models ----
model_dml_b <- cmdstan_model("../dml_b.stan")
model_dml_b2 <- cmdstan_model("../dml_b2.stan")
model_dml_r2d2 <- cmdstan_model("../dml_r2d2.stan")
model_dml_b2_iw <- cmdstan_model("../dml_b2_iw.stan")

## Function to fit model B ----
fit_model_dml_b <- function(data, seed) {
  N <- nrow(data$X)
  P <- ncol(data$X)
  
  fitted_dml_b <- model_dml_b$sample(
    data = list(K = 2, J = P, N = N, x = data$X, y = cbind(data$Y, data$A)),
    seed = seed,
    chains = 1,
    parallel_chains = 1,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
}

## Function to fit model B2 (hierarchical) ----
fit_model_dml_b2 <- function(data, seed) {
  N <- nrow(data$X)
  P <- ncol(data$X)
  
  model_dml_b2$sample(
    data = list(K = 2, J = P, N = N, x = data$X, y = cbind(data$Y, data$A)),
    seed = seed,
    chains = 1,
    parallel_chains = 1,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
}

## Function to fit model R2D2 ----
fit_model_dml_r2d2 <- function(data, seed) {
  N <- nrow(data$X)
  P <- ncol(data$X)

  model_dml_r2d2$sample(
    data = list(J = P, N = N, x = data$X, y = cbind(data$Y, data$A), b = 0.5),
    seed = seed,
    chains = 1,
    parallel_chains = 1,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
}

## Function to fit model IW-Hierarchical ----
fit_model_dml_b2_iw <- function(data, seed) {
  N <- nrow(data$X)
  P <- ncol(data$X)
  
  model_dml_b2_iw$sample(
    data = list(J = P, N = N, x = data$X, y = cbind(data$Y, data$A)),
    seed = seed,
    chains = 1,
    parallel_chains = 1,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
}

## Function to extract results ----
extract_results_dml <- function(fit, gamma, type, additional_results_info) {
#' Extract results from the DML-B2 model fit
#'
#' @param fit A fitted model object from which to extract results.
#' @param gamma The true value of the parameter gamma, called "alpha" in the paper
#' @param type The type of model used, e.g. "DML_B2"
#' @return A data frame containing the extracted results, including:
#'         - gamma_hat: The estimated value of gamma.
#'         - squared_error: The squared error of the estimate.
#'         - LCL: The lower confidence limit of the interval.
#'         - UCL: The upper confidence limit of the interval.
#'         - catch: Whether the true gamma is within the confidence interval.
#'         - interval_width: The width of the confidence interval.
#'         - Method: The method used ("DML-B2").
#'
#' @details This function extracts the posterior draws of the parameter alpha
#'          from the fitted model, computes summary statistics, and checks
#'          whether the true gamma value is within the computed interval.
  draws <- fit$draws("alpha")
  fit_stat <- data.frame(
    mean = mean(draws),
    q2.5 = quantile(draws, 0.025),
    q97.5 = quantile(draws, 0.975)
  )
  gamma_hat <- fit_stat$mean
  interval <- c(fit_stat$q2.5, fit_stat$q97.5)
  squared_error <- (gamma_hat - gamma)^2
  catch <- check_interval(gamma, interval) # Uses check_interval.R
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
    Method = type # e.g. "DML-B2"
  )

  # Append additional_results_info (setting parameters) to pass through for downstream analysis
  for (name in names(additional_results_info)) {
    table[[name]] <- additional_results_info[[name]]
  }
  table
}

## Main simulation function BDML_b for given setting ----
sim_iter_BDML_b <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit <- fit_model_dml_b(data, seed)
  res <- extract_results_dml(fit, data$gamma, type = "BDML_b", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  res
}

## Main simulation function BDML_b2 for given setting ----
sim_iter_BDML_b2 <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit <- fit_model_dml_b2(data, seed)
  res <- extract_results_dml(fit, data$gamma, type = "BDML_b2", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  res
}

## Main simulation function BDML_r2d2 for given setting ----
sim_iter_BDML_r2d2 <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit_r2d2 <- fit_model_dml_r2d2(data, seed)
  res_r2d2 <- extract_results_dml(fit_r2d2, data$gamma, type = "BDML_r2d2", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  res_r2d2
}

## Main simulation function BDML_b2_iw for given setting ----
sim_iter_BDML_b2_iw <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit_iw <- fit_model_dml_b2_iw(data, seed)
  res_iw <- extract_results_dml(fit_iw, data$gamma, type = "BDML_b2_iw", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  res_iw
}