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
  n <- nrow(data$X)
  p <- ncol(data$X)
  
  fitted_bdml_lkj <- model_bdml_lkj$sample(
    data = list(k = 2, p = p, n = n, X = data$X, Y = cbind(data$Y, data$D)),
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
  n <- nrow(data$X)
  p <- ncol(data$X)
  
  model_bdml_lkj_hp$sample(
    data = list(k = 2, p = p, n = n, X = data$X, Y = cbind(data$Y, data$D)),
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
  n <- nrow(data$X)
  p <- ncol(data$X)

  model_bdml_r2d2$sample(
    data = list(p = p, n = n, X = data$X, Y = cbind(data$Y, data$D), b = 0.5),
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
  n <- nrow(data$X)
  p <- ncol(data$X)
  
  model_bdml_iw_hp$sample(
    data = list(p = p, n = n, X = data$X, Y = cbind(data$Y, data$D)),
    seed = seed,
    chains = 1,
    parallel_chains = 1,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
}

## Function to extract results ----
extract_results_bdml <- function(fit, alpha, type, additional_results_info) {
#' Extract results from the DML-B2 model fit
#'
#' @param fit A fitted model object from which to extract results.
#' @param alpha The true value of the parameter alpha, called "alpha" in the paper
#' @param type The type of model used, e.g. "BDML-LKJ-HP"
#' @return A data frame containing the extracted results, including:
#'         - alpha_hat: The estimated value of alpha.
#'         - squared_error: The squared error of the estimate.
#'         - LCL: The lower confidence limit of the interval.
#'         - UCL: The upper confidence limit of the interval.
#'         - catch: Whether the true alpha is within the confidence interval.
#'         - interval_width: The width of the confidence interval.
#'         - Method: The method used ("DML-B2").
#'
#' @details This function extracts the posterior draws of the parameter alpha
#'          from the fitted model, computes summary statistics, and checks
#'          whether the true alpha value is within the computed interval.
  draws <- fit$draws("alpha")
  fit_stat <- data.frame(
    mean = mean(draws),
    q2.5 = quantile(draws, 0.025),
    q97.5 = quantile(draws, 0.975)
  )
  alpha_hat <- fit_stat$mean
  interval <- c(fit_stat$q2.5, fit_stat$q97.5)
  squared_error <- (alpha_hat - alpha)^2
  catch <- check_interval(alpha, interval) # Uses check_interval.R
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]
  
  table <- data.frame(
    alpha_hat = alpha_hat, 
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
sim_iter_bdml_lkj <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit <- fit_model_bdml_lkj(data, seed)
  res <- extract_results_bdml(fit, data$alpha, type = "BDML-LKJ", additional_results_info = list(R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, n = n, p = p))
  res
}

## Main simulation function bdml_lkj_hp for given setting ----
sim_iter_bdml_lkj_hp <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit <- fit_model_bdml_lkj_hp(data, seed)
  res <- extract_results_bdml(fit, data$alpha, type = "BDML-LKJ-HP", additional_results_info = list(R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, n = n, p = p))
  res
}

## Main simulation function bdml_r2d2 for given setting ----
sim_iter_bdml_r2d2 <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit_r2d2 <- fit_model_bdml_r2d2(data, seed)
  res_r2d2 <- extract_results_bdml(fit_r2d2, data$alpha, type = "BDML-R2D2", additional_results_info = list(R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, n = n, p = p))
  res_r2d2
}

## Main simulation function bdml_iw_hp for given setting ----
sim_iter_bdml_iw_hp <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit_iw <- fit_model_bdml_iw_hp(data, seed)
  res_iw <- extract_results_bdml(fit_iw, data$alpha, type = "BDML-IW-HP", additional_results_info = list(R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, n = n, p = p))
  res_iw
}
