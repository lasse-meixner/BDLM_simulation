gen_ridge_data <- function(N, P, tau, omega, gamma) {
  X <- matrix(rnorm(N * P),nrow = N, ncol = P)
  phi <- tau * rnorm(P) / sqrt(P)
  b <- tau * rnorm(P) / sqrt(P)
  beta <- omega * phi + b
  
  A <- as.numeric(X %*% phi) + rnorm(N)
  Y <- as.numeric(X %*% beta) + gamma * A + rnorm(N)
  
  return(list(X, A, Y, phi, b, beta))
}