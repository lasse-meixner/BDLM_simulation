# Auxiliary functions for data generation ----
## General function to prepare root data ----

library(tidyverse)
library(mvtnorm)

#------------------------------------------------------------------------------#
# 1) ROOT PREP  
#------------------------------------------------------------------------------#
prepare_root_data <- function(n, p) {
  list(X = matrix(rnorm(n * p), nrow = n), n = n, p = p)
}

#------------------------------------------------------------------------------#
# 2) STRUCTURAL COMPONENTS BUILDER  
#------------------------------------------------------------------------------#
build_structures <- function(data,
                             gamma_fun = function(p) rep(1/sqrt(p), p),
                             beta_rndcomp_fun  = function(p) rnorm(p, 0, 1/sqrt(p))) {
  data$gamma  <- gamma_fun(data$p)
  data$beta_rndcomp   <- beta_rndcomp_fun(data$p)
  data$D_hat  <- as.numeric(data$X %*% data$gamma)
  data$D      <- data$D_hat + rnorm(data$n)
  data$mu_hat <- as.numeric(data$X %*% data$beta_rndcomp)
  data
}

#------------------------------------------------------------------------------#
# 3) OUTCOME BUILDER  
#------------------------------------------------------------------------------#
build_outcome <- function(data, sigma) {
  data$Y     <- data$alpha * data$D +
               data$w * data$D_hat +
               data$mu_hat +
               sigma * rnorm(data$n)
  data$sigma <- sigma
  data
}

#------------------------------------------------------------------------------#
# 4) SETTING-SPECIFIC GENERATORS  
#------------------------------------------------------------------------------#

generate_random_data <- function(n, p, sigma) {
  data <- prepare_root_data(n, p)
  data$alpha <- rnorm(1)
  data$w     <- rnorm(1)
  data <- build_structures(data)
  build_outcome(data, sigma)
}

generate_fixed_data <- function(n, p, sigma) {
  data <- prepare_root_data(n, p)
  data$alpha <- 2
  data$w     <- -0.5
  data <- build_structures(data)
  build_outcome(data, sigma)
}

generate_noisy_fs_data <- function(n, p, sigma) {
  data <- prepare_root_data(n, p)
  data$alpha <- 2
  data$w     <- -0.5
  data <- build_structures(
    data,
    gamma_fun = function(p) runif(p, 0, 1/sqrt(p)),
    beta_rndcomp_fun  = function(p) rnorm(p, 0, 1/sqrt(p))
  )
  build_outcome(data, sigma)
}

generate_hahn_data <- function(n, p, sigma) {
  data <- prepare_root_data(n, p)
  data$alpha <- 2
  data$w     <- -2
  data <- build_structures(data)
  build_outcome(data, sigma)
}

generate_naive_data <- function(n, p, sigma) {
  data <- prepare_root_data(n, p)
  data$alpha <- 2
  data$w     <- 0
  data <- build_structures(data)
  build_outcome(data, sigma)
}

generate_mixed_data <- function(n, p, sigma) {
  data <- prepare_root_data(n, p)
  data$alpha <- 1
  data$w     <- -1
  data <- build_structures(
    data,
    gamma_fun = function(p) {
      θ <- rep(1/sqrt(p), p)
      θ[floor(p/2):p] <- 0
      θ
    }
  )
  # mixed uses A_hat only
  build_outcome(data, sigma)
}

generate_joint_data <- function(n, p, sigma = NULL) {
  data <- prepare_root_data(n, p)
  data$delta <- 0.5
  B <- matrix(rnorm(p * 2), nrow = 2)
  Sigma <- matrix(c(1, data$delta, data$delta, 1), 2, 2)
  errors <- rmvnorm(n, mean = c(0, 0), sigma = Sigma)

  W <- data$X %*% t(B) + errors
  data$Y <- W[,1]
  data$D <- W[,2]
  data$sigma <- NA
  data
}

#------------------------------------------------------------------------------#
# 5) DISPATCHER  
#------------------------------------------------------------------------------#
#' Generate data based on the specified setting
#'
#' @param n       Number of obs
#' @param p       Number of predictors
#' @param setting One of "random", "fixed", "noisy_fs", 
#'                "hahn", "naive", "mixed", "joint"
#' @param sigma   Noise sd (ignored for "joint")
generate_data <- function(n, p, setting, sigma = 1) {
  switch(setting,
    random    = generate_random_data(n, p, sigma),
    fixed     = generate_fixed_data(n, p, sigma),
    noisy_fs  = generate_noisy_fs_data(n, p, sigma),
    hahn      = generate_hahn_data(n, p, sigma),
    naive     = generate_naive_data(n, p, sigma),
    mixed     = generate_mixed_data(n, p, sigma),
    joint     = generate_joint_data(n, p),
    stop("Unknown setting: ", setting)
  )
}
