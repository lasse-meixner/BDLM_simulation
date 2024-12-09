indirect_bias <- function(lambda, eta, r, eval) {
  
  psi_00 <- 1
  psi_10 <- mean(1 / (lambda + eval))
  psi_20 <- mean(1 / (lambda + eval)^2)
  psi_30 <- mean(1 / (lambda + eval)^3)
  
  psi_01 <- mean(eval)
  psi_11 <- psi_00 - lambda * psi_10
  psi_21 <- psi_10 - lambda * psi_20
  psi_31 <- psi_20 - lambda * psi_30
  
  psi_02 <- mean(eval^2)
  psi_12 <- psi_01 - lambda * psi_11
  psi_22 <- psi_11 - lambda * psi_21
  psi_32 <- psi_21 - lambda * psi_31
  
  psi_03 <- mean(eval^3)
  psi_13 <- psi_02 - lambda * psi_12
  psi_23 <- psi_12 - lambda * psi_22
  psi_33 <- psi_22 - lambda * psi_32
  
  p_t <- 1 / (lambda * psi_32 + lambda/eta * psi_33)
  final_numerator <- lambda * ((psi_11 + psi_12/eta) * (-p_t * psi_22 / eta) + 
                                 psi_11 / eta + 
                                 (psi_22 + psi_23 / eta) * psi_22 / eta * p_t)
  
  # final_numerator <- lambda * (psi_11 / eta - psi_22 / eta * (psi_21 + psi_22/eta) / (psi_32 + psi_33/eta))
  final_denominator <- (1-r)/r + lambda * (psi_10 - (psi_21 + psi_22/eta)^2/(psi_32 + psi_33/eta) + psi_11/eta)
  
  final_numerator / final_denominator  
}