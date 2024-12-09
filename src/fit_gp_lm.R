fit_gp_lm <- function(X, Y, Z = NULL, a, b, s_sq, X_test, Z_test) {
  Sigma <- a * exp(-b * as.matrix(dist(X))^2)
  Sigma_test <- a * exp(-b * as.matrix(dist(X_test))^2)
  Sigma_cross <- a * exp(-b * as.matrix(pdist::pdist(X_test, X))^2)
  
  if(!is.null(Z)) {
    num_lin <- ncol(Z)
    for(j in 1:num_lin) {
      Sigma <- Sigma + Z[,j] %*% t(Z[,j]) * 100
      Sigma_test <- Sigma_test + Z_test[,j] %*% t(Z_test[,j]) * 100
      Sigma_cross <- Sigma_cross + Z_test[,j] %*% t(Z[,j]) * 100
    }    
  }
  
  K <- Sigma + s_sq * diag(length(Y))
  
  mu_hat <- Sigma %*% solve(K, Y)
  S_hat  <- Sigma - Sigma %*% solve(K, Sigma)
  
  mu_test <- Sigma_cross %*% solve(K,Y)
  S_test  <- Sigma_test - Sigma_cross %*% solve(K, t(Sigma_cross))
  
  return(list(mu_hat = mu_hat, S_hat = S_hat, 
              mu_test = mu_test, S_test = S_test))
  
}