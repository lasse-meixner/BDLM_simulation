library(tidyverse)
library(mvtnorm)

# Auxiliary functions for data generation ----
## General function to prepare root data ----

prepare_root_data <- function(N, P) {
  X <- matrix(rnorm(n = N * P), nrow = N)
  theta <- rep(1/sqrt(P), P) # gamma in the paper
  beta <- rnorm(P, 0, 1/sqrt(P))
  A_hat <- as.numeric(X %*% theta)
  A <- A_hat + rnorm(N)
  mu_hat <- as.numeric(X %*% beta)
  list(X = X, theta = theta, beta = beta, A_hat = A_hat, A = A, mu_hat = mu_hat)
}

## Setting specific data tweak functions ----
tweak_data_random <- function(data, N) {
  data$gamma <- rnorm(1) # alpha in paper
  data$w <- rnorm(1) # kappa in paper
  data
}

tweak_data_fixed <- function(data) {
  data$gamma <- 2 # alpha in paper
  data$w <- -0.5 # kappa in paper
  data
}

tweak_data_noisy_fs <- function(data) {
  data <- tweak_data_fixed(data)
  # override theta
  data$theta <- runif(P, 0, 1/sqrt(ncol(data$X)))
  data
}

tweak_data_hahn <- function(data) {
  data$gamma <- 2
  data$w <- -2
  data
}

tweak_data_naive <- function(data) {
  data$gamma <- 2
  data$w <- 0
  data
}

tweak_data_mixed <- function(data, P) {
  data$gamma <- 1
  data$w <- -1
  data$theta_w <- data$theta
  data$theta_w[floor(P/2):P] <- 0 # Set half of the coefficients to zero
  data$A_hat <- as.numeric(data$X %*% data$theta_w)
  data
}

tweak_data_joint <- function(data, N, P) {
  data$delta <- 0.5 # renamed this to delta to avoid maximum confusion
  B <- matrix(rnorm(n = P * 2), nrow = 2)
  Sigma <- matrix(c(1, data$delta, data$delta, 1), byrow = TRUE, nrow = 2)
  errors <- rmvnorm(n = N, mean = rep(0, 2), sigma = Sigma)
  W <- data$X %*% t(B) + errors
  data$Y <- W[, 1]
  data$A <- W[, 2]
  data
}

## Main Function to generate data based on setting ----
generate_data <- function(N, P, setting, sigma) {
#' Generate data based on the specified setting
#'
#' @param N An integer specifying the number of observations.
#' @param P An integer specifying the number of predictors.
#' @param setting A string specifying the data generation setting. 
#'                Options are "random", "fixed", "hahn", "naive", "mixed", "joint".
#' @param sigma A numeric value specifying the standard deviation of the noise in the outcome equation.
#' @return A list containing the generated data, including predictors (X), 
#'         response (Y), and other relevant parameters based on the setting.
#'
#' @details This function generates data according to different settings. 
#'          Each setting tweaks the data in a specific way:
#'          - "random": Random gamma and w values.
#'          - "fixed": Fixed gamma and w values.
#'          - "hahn": Specific gamma and w values for Hahn setting.
#'          - "naive": Specific gamma and w values for Naive setting.
#'          - "mixed": Mixed setting with some coefficients set to zero.
#'          - "joint": Joint setting with correlated errors.
  data <- prepare_root_data(N, P)
  
  if (setting == "random") {
    data <- tweak_data_random(data, N)
  } else if (setting == "fixed") {
    data <- tweak_data_fixed(data)
  } else if (setting == "noisy_fs") {
    data <- tweak_data_noisy_fs(data)
  } else if (setting == "hahn") {
    data <- tweak_data_hahn(data)
  } else if (setting == "naive") {
    data <- tweak_data_naive(data)
  } else if (setting == "mixed") {
    data <- tweak_data_mixed(data, P)
  } else if (setting == "joint") {
    data <- tweak_data_joint(data, N, P)
    return(data)
  }
  if (setting != "joint") {
    data$Y <- data$gamma * data$A + data$w * data$A_hat + data$mu_hat + sigma * rnorm(N)
  }

  # save sigma
  data$sigma <- sigma

  # return data
  data
}