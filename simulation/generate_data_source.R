# Auxiliary functions for data generation ----
## General function to prepare root data ----

library(tidyverse)
library(mvtnorm)

#------------------------------------------------------------------------------#
# 1) ROOT PREP  
#------------------------------------------------------------------------------#
prepare_root_data <- function(N, P) {
  list(X = matrix(rnorm(N * P), nrow = N), N = N, P = P)
}

#------------------------------------------------------------------------------#
# 2) STRUCTURAL COMPONENTS BUILDER  
#------------------------------------------------------------------------------#
build_structures <- function(dat,
                             theta_fun = function(P) rep(1/sqrt(P), P),
                             beta_fun  = function(P) rnorm(P, 0, 1/sqrt(P))) {
  dat$theta  <- theta_fun(dat$P)
  dat$beta   <- beta_fun(dat$P)
  dat$A_hat  <- as.numeric(dat$X %*% dat$theta)
  dat$A      <- dat$A_hat + rnorm(dat$N)
  dat$mu_hat <- as.numeric(dat$X %*% dat$beta)
  dat
}

#------------------------------------------------------------------------------#
# 3) OUTCOME BUILDER  
#------------------------------------------------------------------------------#
build_outcome <- function(dat, sigma) {
  dat$Y     <- dat$gamma * dat$A +
               dat$w     * dat$A_hat +
               dat$mu_hat +
               sigma * rnorm(dat$N)
  dat$sigma <- sigma
  dat
}

#------------------------------------------------------------------------------#
# 4) SETTING-SPECIFIC GENERATORS  
#------------------------------------------------------------------------------#

generate_random_data <- function(N, P, sigma) {
  dat <- prepare_root_data(N, P)
  dat$gamma <- rnorm(1)
  dat$w     <- rnorm(1)
  dat <- build_structures(dat)
  build_outcome(dat, sigma)
}

generate_fixed_data <- function(N, P, sigma) {
  dat <- prepare_root_data(N, P)
  dat$gamma <- 2
  dat$w     <- -0.5
  dat <- build_structures(dat)
  build_outcome(dat, sigma)
}

generate_noisy_fs_data <- function(N, P, sigma) {
  dat <- prepare_root_data(N, P)
  dat$gamma <- 2
  dat$w     <- -0.5
  dat <- build_structures(
    dat,
    theta_fun = function(P) runif(P, 0, 1/sqrt(P)),
    beta_fun  = function(P) rnorm(P, 0, 1/sqrt(P))
  )
  build_outcome(dat, sigma)
}

generate_hahn_data <- function(N, P, sigma) {
  dat <- prepare_root_data(N, P)
  dat$gamma <- 2
  dat$w     <- -2
  dat <- build_structures(dat)
  build_outcome(dat, sigma)
}

generate_naive_data <- function(N, P, sigma) {
  dat <- prepare_root_data(N, P)
  dat$gamma <- 2
  dat$w     <- 0
  dat <- build_structures(dat)
  build_outcome(dat, sigma)
}

generate_mixed_data <- function(N, P, sigma) {
  dat <- prepare_root_data(N, P)
  dat$gamma <- 1
  dat$w     <- -1
  dat <- build_structures(
    dat,
    theta_fun = function(P) {
      θ <- rep(1/sqrt(P), P)
      θ[floor(P/2):P] <- 0
      θ
    }
  )
  # mixed uses A_hat only
  build_outcome(dat, sigma)
}

generate_joint_data <- function(N, P, sigma = NULL) {
  dat <- prepare_root_data(N, P)
  dat$delta <- 0.5
  B <- matrix(rnorm(P * 2), nrow = 2)
  Sigma <- matrix(c(1, dat$delta, dat$delta, 1), 2, 2)
  errors <- rmvnorm(N, mean = c(0, 0), sigma = Sigma)

  W <- dat$X %*% t(B) + errors
  dat$Y <- W[,1]
  dat$A <- W[,2]
  dat$sigma <- NA
  dat
}

#------------------------------------------------------------------------------#
# 5) DISPATCHER  
#------------------------------------------------------------------------------#
#' Generate data based on the specified setting
#'
#' @param N       Number of obs
#' @param P       Number of predictors
#' @param setting One of "random", "fixed", "noisy_fs", 
#'                "hahn", "naive", "mixed", "joint"
#' @param sigma   Noise sd (ignored for "joint")
generate_data <- function(N, P, setting, sigma = 1) {
  switch(setting,
    random    = generate_random_data(N, P, sigma),
    fixed     = generate_fixed_data(N, P, sigma),
    noisy_fs  = generate_noisy_fs_data(N, P, sigma),
    hahn      = generate_hahn_data(N, P, sigma),
    naive     = generate_naive_data(N, P, sigma),
    mixed     = generate_mixed_data(N, P, sigma),
    joint     = generate_joint_data(N, P),
    stop("Unknown setting: ", setting)
  )
}
