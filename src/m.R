m <- function(lambda, rho) {
  v <- 1 - rho + lambda
  out <- -v + sqrt(v^2 + 4 * lambda * rho)
  out <- out / 2 / rho / lambda
  return(out)
}

mm <- function(lambda,rho) {
  numerator <- -lambda * (1 + rho) - (-1 + rho) * (-1 + rho + sqrt((1+lambda)^2 + 2 * (-1 + lambda) * rho + rho^2))
  denominator <- 2 * lambda^2 * rho * sqrt(lambda^2 + (-1 + rho)^2 + 2 * lambda * (1 + rho))
  return(numerator/denominator)
}