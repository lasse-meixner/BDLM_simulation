theoretical_mse <- function(omega, lambda, eta, r) {
  
  tdenom <- 1 - r + lambda * r * m(lambda,r) + lambda/eta * (r - lambda * r * m(lambda,r))
  
  out <- 
    (1 + omega^2) * lambda^2 * r * m(lambda,r) / eta + # Written
    lambda^3 * r * mm(lambda,r) / eta + # Written
    r^2 * lambda^2 * omega^2 * (1 - lambda * m(lambda,r))^2 / eta^2 / tdenom + # Written
    lambda^4 * omega^2 /eta/(lambda-eta) * r * mm(lambda,r) + # Written
    -lambda^3 * omega^2 * (r * (m(lambda,r))/(lambda-eta)^2) + #Written
    -lambda^3 * omega^2 * 2 * (r/eta - lambda * r * m(lambda,r)/eta) * (r * (m(lambda,r) + lambda * mm(lambda,r))/eta)/tdenom +
    -lambda^3 * omega^2 * (r/eta - lambda * r  * m(lambda,r)/eta)^2 * (r / eta * (1 + (eta - 2 * lambda) * m(lambda,r) - lambda * (lambda-eta) * mm(lambda,r))) / (tdenom)^2 +
    -lambda^3 * omega^2 * r * (-mm(lambda,r) - m(lambda,r)/(eta-lambda))/(eta-lambda) + 
    r - 2 * lambda * r * m(lambda,r) - lambda^2 * r * mm(lambda,r)
  
  return(out)
}