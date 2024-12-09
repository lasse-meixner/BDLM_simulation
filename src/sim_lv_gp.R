sim_lv_gp <- function(N, P, sigma_x) {
  eta <- rnorm(N)
  my_gps <- sim_gp(eta, P)
  
  X <- scale(t(my_gps + sigma_x * rnorm(length(my_gps))))
  
  sim_f <- function(X) {
    K <- exp(-as.matrix(dist(X))^2)
    mu <- MASS::mvrnorm(n = 2, mu = rep(0, nrow(X)), Sigma = K)
    return(t(mu))
  }
  
  fg <- sim_f(X)
  
  A <- fg[,1] + rnorm(nrow(fg))
  Y <- fg[,1] + fg[,2] + A + rnorm(nrow(fg))

  return(list(X = X, Y = Y, A = A, fg = fg))
  
}