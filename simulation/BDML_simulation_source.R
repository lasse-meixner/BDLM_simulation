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
model_bdml_lkj <- cmdstan_model("../stan/bdml_lkj.stan")
model_bdml_lkj_hp <- cmdstan_model("../stan/bdml_lkj_hp.stan")
model_bdml_r2d2 <- cmdstan_model("../stan/bdml_r2d2.stan")
model_bdml_iw_hp <- cmdstan_model("../stan/bdml_iw_hp.stan")

## Function to fit model B ----
fit_model_bdml_lkj <- function(data, seed) {
  N <- nrow(data$X)
  P <- ncol(data$X)
  
  fitted_bdml_lkj <- model_bdml_lkj$sample(
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
fit_model_bdml_lkj_hp <- function(data, seed) {
  N <- nrow(data$X)
  P <- ncol(data$X)
  
  model_bdml_lkj_hp$sample(
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
fit_model_bdml_r2d2 <- function(data, seed) {
  N <- nrow(data$X)
  P <- ncol(data$X)

  model_bdml_r2d2$sample(
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
fit_model_bdml_iw_hp <- function(data, seed) {
  N <- nrow(data$X)
  P <- ncol(data$X)
  
  model_bdml_iw_hp$sample(
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
extract_results_bdml <- function(fit, gamma, type, additional_results_info) {
#' Extract results from the DML-B2 model fit
#'
#' @param fit A fitted model object from which to extract results.
#' @param gamma The true value of the parameter gamma, called "alpha" in the paper
#' @param type The type of model used, e.g. "BDML-LKJ-HP"
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

## Main simulation function bdml_lkj for given setting ----
sim_iter_bdml_lkj <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit <- fit_model_bdml_lkj(data, seed)
  res <- extract_results_bdml(fit, data$gamma, type = "BDML-LKJ", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  res
}

## Main simulation function bdml_lkj_hp for given setting ----
sim_iter_bdml_lkj_hp <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit <- fit_model_bdml_lkj_hp(data, seed)
  res <- extract_results_bdml(fit, data$gamma, type = "BDML-LKJ-HP", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  res
}

## Main simulation function bdml_r2d2 for given setting ----
sim_iter_bdml_r2d2 <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit_r2d2 <- fit_model_bdml_r2d2(data, seed)
  res_r2d2 <- extract_results_bdml(fit_r2d2, data$gamma, type = "BDML-R2D2", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  res_r2d2
}

## Main simulation function bdml_iw_hp for given setting ----
sim_iter_bdml_iw_hp <- function(N, P, setting, sigma, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(N, P, setting, sigma)
  fit_iw <- fit_model_bdml_iw_hp(data, seed)
  res_iw <- extract_results_bdml(fit_iw, data$gamma, type = "BDML-IW-HP", additional_results_info = list(setting = setting, sigma = sigma, N = N, P = P))
  res_iw
}