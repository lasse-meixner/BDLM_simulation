ridge_eb <- function(Y, eigen_decomp, transformed = FALSE) {
  Gamma <- eigen_decomp$vectors
  evals <- eigen_decomp$values
  if(!transformed) {
    Y <- t(Gamma) %*% Y
  }

  neg_loglik <- function(lambda) {
    - 0.5 * sum(log(lambda / (lambda + evals))) +
    0.5 * sum(Y^2 * lambda / (lambda + evals))
  }

  out <- optimize(neg_loglik, c(0.01, 5))
  return(out)

}
