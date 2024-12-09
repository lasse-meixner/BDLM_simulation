theoretical_bias <- function(omega,lambda,eta,r) {
  m <- m(lambda,r)
  numerator <- r * lambda / eta * (1 - lambda * m)
  denominator <- (1 - r) + lambda/eta*r*(1 + (eta - lambda) * m)
  return(omega * numerator / denominator)
}