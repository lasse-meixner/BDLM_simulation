sim_gp <- function(eta, P) {
  K <- exp(-as.matrix(dist(eta))^2)
  m <- MASS::mvrnorm(n = P, mu = rep(0, length(eta)), Sigma = K)
  return(m)
}
