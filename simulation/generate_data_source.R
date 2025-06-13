# Auxiliary functions for data generation ----

library(tidyverse)
library(mvtnorm)

generate_data <- function(n, p, R_Y2, R_D2, rho, alpha) {
  
  # 1a) build the 2×2 covariance for (beta_j,gamma_j)'
  Sigma <- (1/p) * matrix(
    c(R_Y2,
      rho * sqrt(R_Y2 * R_D2),
      rho * sqrt(R_Y2 * R_D2),
      R_D2),
    nrow = 2, byrow = TRUE
  )
  
  # 1b) draw p independent coefficient pairs
  #     (β_j, γ_j)' ~ N(0, Sigma)
  B <- rmvnorm(n = p, mean = c(0, 0), sigma = Sigma)
  beta  <- B[,1]
  gamma <- B[,2]
  
  # 2) residual standard‐deviations
  sigma_eps <- sqrt(1 - R_Y2)
  sigma_V   <- sqrt(1 - R_D2)
  
  # 3) simulate data
  X   <- matrix(rnorm(n * p), nrow = n, ncol = p)        # regressors
  V   <- rnorm(n, mean = 0, sd = sigma_V)                # first‐stage noise
  eps <- rnorm(n, mean = 0, sd = sigma_eps)              # outcome noise
  
  D <- as.numeric(X %*% gamma + V)
  Y <- as.numeric(alpha * D + X %*% beta + eps)
  
  # return a list
  list(
    X     = X,
    D     = D,
    Y     = Y,
    R_Y2  = R_Y2,
    R_D2  = R_D2,
    rho   = rho,
    alpha = alpha
  )
}