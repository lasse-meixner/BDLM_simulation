# tau <- 3
# lambda <- 1
# gamma <- 3
# omega <- 2.5
# P <- 4000
# N <- 3000

# eta = c(0.4, 2),
# gamma = 2,
# omega = 1, r = c(0.25, 0.5, 1, 2),
# lambda = lambda_grid

r <- 2
# tau <- 1
# eta <- r/(tau^2)
eta <- 0.4
tau <- sqrt(r / eta)
lambda <- 0.1
gamma <- 2
omega <- 1
N <- 3000
P <- r * N

ridge_test_direct <- function(tau = 3, lambda = NULL, gamma = 3, omega = 1, N = 300, P = 300) {
    
  tau_sq <- tau^2
  if(is.null(lambda)) lambda <- P / N / tau_sq
  r <- P / N
  eta <- r / tau_sq
  
  X <- matrix(rnorm(N * P),nrow = N, ncol = P)
  phi <- tau * rnorm(P) / sqrt(P)
  b   <- tau * rnorm(P) / sqrt(P)
  beta <- omega * phi + b
  
  ## Matrix Multiplies
  
  Xphi <- X %*% phi
  Xbeta <- X %*% beta
  
  ## Generate Data
  
  A <- (Xphi + rnorm(N)) |> as.numeric()
  Y <- (Xbeta + gamma * A + rnorm(N)) |> as.numeric()
  
  ## Fit the A-model
  
  V_bar   <- solve(t(X) %*% X + N * lambda * diag(P))
  phi_hat <- V_bar %*% t(X) %*% A
  C_hat   <- X %*% phi_hat
  
  ## Make more matricies
  
  Psi_hat     <- cbind(C_hat, X)
  H_bar_inv   <- t(Psi_hat) %*% Psi_hat + N * lambda * rbind(0,cbind(0, diag(P)))
  H_bar       <- solve(H_bar_inv)
  
  Qbar        <- diag(N) - Psi_hat %*% H_bar %*% t(Psi_hat)
  gamma_hat <- sum(t(A) %*% Qbar %*% Y) / sum(t(A) %*% Qbar %*% A)
  
  Z <- cbind(A, C_hat, X)
  M_bar_inv <- t(Z) %*% Z + N * lambda * rbind(0, 0, cbind(0, 0, diag(P)))
  M_bar <- solve(M_bar_inv)
  theta_hat <- M_bar %*% t(Z) %*% Y
  gamma_hat_2 <- theta_hat[1]
  
  return(list(gamma_hat = gamma_hat, X = X))
  # return(c(gamma_hat = gamma_hat, gamma_hat_2 = gamma_hat_2))
  
  ## Denominator
  d <- t(A) %*% Qbar %*% A
  VA <- diag(N) + X %*% t(X) / N / eta
  VAT <- diag(P) + t(X) %*% X / N / eta
  S <- t(X) %*% X / N
  V <- N * V_bar
  K_bar   <- solve(t(X) %*% X + N * eta * diag(P))
  phi_tilde <- K_bar %*% t(X) %*% A
  C_tilde   <- X %*% phi_tilde
  numerator <- t(A) %*% Qbar %*% C_tilde
  
  eval <- eigen(S)[["values"]]
  
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
  
  ## Components of H_bar
  p <- 1 / sum(t(C_hat) %*% (diag(N) - X %*% V_bar %*% t(X)) %*% C_hat)
  q <- as.numeric(-V_bar %*% t(X) %*% C_hat * p)
  W <- V_bar + p * V_bar %*% t(X) %*% C_hat %*% t(C_hat) %*% X %*% V_bar
  
  ## Checking
  t(A) %*% A - t(A) %*% Psi_hat %*% H_bar %*% t(Psi_hat) %*% A
  
  t(A) %*% A - p * (t(A) %*% C_hat)^2 - (t(A) %*% X %*% q) * (t(A) %*% C_hat) - 
    (t(A) %*% C_hat) * (t(q) %*% t(X) %*% A) - t(A) %*% X %*% V_bar %*% t(X) %*% A - 
    p * t(A) %*% X %*% V_bar %*% t(X) %*% C_hat %*% t(C_hat) %*% X %*% V_bar %*% t(X) %*% A
  
  phi_10 <- psi_10 + psi_11 / eta
  phi_21 <- psi_21 + psi_22 / eta
  
  (1-r)/r + lambda * (phi_10 - phi_21^2 / phi_32)
  
  ## Checking numerator
  
  numerator
  t(A) %*% C_tilde - 
    p * (t(A) %*% C_hat) * (t(C_hat) %*% C_tilde) + 
    p * (t(C_hat) %*% C_hat) * (t(C_hat) %*% C_tilde) +
    p * (t(C_tilde) %*% X %*% V_bar %*% t(X) %*% C_hat) * (t(C_hat) %*% A) - 
    t(C_hat) %*% C_tilde - 
    p * (t(C_hat) %*% C_hat) * (t(C_hat) %*% X %*% V_bar %*% t(X) %*% C_tilde)
  
  psi_01 / eta - 
    phi_11 * psi_12 / lambda / eta / phi_32 + 
    phi_22 * psi_12 / lambda / eta / phi_32 + 
    phi_11 * psi_23 / lambda / eta / phi_32 -
    psi_12 / eta - 
    phi_22 * psi_23 / lambda / eta / phi_32
  
  # Notes:
  # t(A) %*% C_hat -> phi_11
  # p -> 1 / (lambda * phi_32)
  # t(A) %*% X %*% q -> -phi_22 / (lambda * phi_32)
  
  
  ## OLD STUFF
  
  C_hat_A <- mtrace(X %*% V_bar %*% t(X) %*% A %*% t(A))
  
  ## Thing I want to approximate:
  sum(A^2) - (C_hat_A^2 * p + 2 * C_hat_A * sum(t(A) %*% X %*% q) + sum(t(A) %*% X %*% W %*% t(X) %*% A))
  sum(A^2) - (C_hat_A^2 * p - 2 * C_hat_A * mtrace(diag(P) - lambda * V - lambda * (V - lambda * V %*% V) +  V %*% S %*% V %*% S %*% S / eta) * p + sum(t(A) %*% X %*% W %*% t(X) %*% A))
  
  (C_hat_A_t <- mtrace(V %*% S) + mtrace(V %*% S %*% S) / eta)
  (C_hat_A_t <- P - lambda * mtrace(V) + mtrace(V %*% S %*% S) / eta)
  (C_hat_A_t <- P - lambda * m(lambda,r) * P + mtrace(V %*% S %*% S) / eta)
  (C_hat_A_t <- P - lambda * m(lambda,r) * P + (mtrace(S) - lambda * (P - lambda * m(lambda,r) * P)) / eta)
  
  # (term_2_t <- P - 2 * lambda * P * m(lambda,r)  - lambda^2 * mm(lambda,r) * P + 
  #     (sum_eig - 2 * lambda * P + 3 * lambda^2 * P * m(lambda,r) + lambda^3 * P * mm(lambda,r))/ eta)
  
  ## Approximating p
  sum(t(C_hat) %*% (diag(N) - X %*% V_bar %*% t(X)) %*% C_hat)
  N * lambda * mtrace(X %*% V_bar %*% V_bar %*% t(X) %*% X %*% V_bar %*% t(X) %*% VA)
  lambda * mtrace(V %*% V %*% V %*% S %*% S) + lambda / eta * mtrace(V %*% V %*% V %*% S %*% S %*% S)
  
   ## New notation:
  ## C_hat_A / P -> psi(1,1) + psi(1,2)/eta
  ## p^-1 / P -> lambda * psi_32 + lambda / eta * psi_33
  ## t(A) %*% X %*% q -> (psi_22 + psi_23 / eta) / (lambda * psi_32 + lambda / eta * psi_33)
  ## t(A) %*% A / P -> 1/r + psi_01/eta
  ## t(A) %*% X %*% V_bar %*% t(X) %*% A /P -> psi_11 + psi_12 / eta
  ## Other W piece: (psi_22 + psi_23 / eta)^2 / (lambda * psi_32 + lambda / eta * psi_33)
  ## t(A) %*% X %*% W %*% t(X) %*% A / P -> psi_11 + psi_12 / eta + (psi_22 + psi_23 / eta)^2 / (lambda * psi_32 + lambda / eta * psi_33)
  
  ## Final expression:
  sum(A^2) - (C_hat_A^2 * p + 2 * C_hat_A * sum(t(A) %*% X %*% q) + sum(t(A) %*% X %*% W %*% t(X) %*% A))
  
  sum_a_sq <- 1/r + psi_01/eta
  p_t <- 1 / (lambda * psi_32 + lambda/eta * psi_33)
  C_hat_A_t <- psi_11 + psi_12/eta
  cross_term_t <- -(psi_22 + psi_23 / eta) / (lambda * psi_32 + lambda / eta * psi_33)
  w_term_t <- psi_11 + psi_12 / eta + (psi_22 + psi_23 / eta)^2 / (lambda * psi_32 + lambda / eta * psi_33)
  
  sum_a_sq - C_hat_A_t^2 * p_t - 2 * C_hat_A_t * cross_term_t - w_term_t
  
  (1-r)/r + lambda * psi_11 / eta - lambda * (psi_21 + psi_22 / eta)^2 / (psi_32 + psi_33 / eta) + lambda * psi_10
  
  ## Final form of the denominator
  (1-r)/r + lambda * (psi_10 - (psi_21 + psi_22/eta)^2/(psi_32 + psi_33/eta) + psi_11/eta)
  
  t1 <- 1/r + psi_01/eta
  t2 <- psi_11 + psi_12 / eta
  t3 <- 1/(lambda * psi_32 + lambda/eta * psi_33)
  t4 <- (psi_22 + psi_23 / eta) / (lambda * psi_32 + lambda / eta * psi_33)
  t6 <- psi_11 + psi_12 / eta + (psi_22 + psi_23 / eta)^2 / (lambda * psi_32 + lambda / eta * psi_33)
  
  t1 - t2 * t2 * t3 - 2 * t2 * t4 - t6
  
  eval <- eigen(S)[["values"]]
  
  # psi_mat <- matrix(NA, nrow = 3, ncol = 3)
  
  
  
  ## Working on the numerator
  sum(t(A) %*% Qbar %*% X %*% phi_tilde) # Correct
  sum(t(A) %*% X %*% K_bar %*% t(X) %*% A - t(A) %*% Psi_hat %*% H_bar %*% t(Psi_hat) %*% X %*% K_bar %*% t(X) %*% A) # Correct
  N * lambda * t(A) %*% Psi_hat %*% H_bar %*% c(0, phi_tilde) # Correct
  N * lambda * (sum(A * C_hat) * sum(q * phi_tilde) + t(A) %*% X %*% W %*% phi_tilde) ## Correct
  N * lambda * (P * (psi_11 + psi_12/eta) * (-p * psi_22 * P / N / eta) + t(A) %*% X %*% W %*% phi_tilde) ## Correlation?
  N * lambda * (P * (psi_11 + psi_12/eta) * (-p * psi_22 * P / N / eta) + 
                  t(A) %*% X %*% V_bar %*% phi_tilde + 
                  t(A) %*% X %*% V_bar %*% t(X) %*% C_hat %*% t(C_hat) %*% X %*% V_bar %*% phi_tilde * p) ## Correct
  N * lambda * (P * (psi_11 + psi_12/eta) * (-p * psi_22 * P / N / eta) + 
                  r * psi_11 / eta + 
                  r * (psi_22 + psi_23 / eta) * P * psi_22 / eta * p) ## Correct
  N * lambda * ((psi_11 + psi_12/eta) * (-p_t * psi_22 * r / eta) + 
                  r * psi_11 / eta + 
                  r * (psi_22 + psi_23 / eta) * psi_22 / eta * p_t) ## Correct
  
  final_numerator <- lambda * ((psi_11 + psi_12/eta) * (-p_t * psi_22 / eta) + 
                                     psi_11 / eta + 
                                     (psi_22 + psi_23 / eta) * psi_22 / eta * p_t)
  # final_numerator <- psi_11 / eta - lambda * psi_22 / eta * (psi_21 + psi_22/eta) / (lambda * psi_32 + lambda / eta * psi_33)
  final_denominator <- (1-r)/r + lambda * (psi_10 - (psi_21 + psi_22/eta)^2/(psi_32 + psi_33/eta) + psi_11/eta)
  
  omega * final_numerator / final_denominator  
  
}


ridge_test_direct_2 <- function(eta, lambda, gamma, omega, N, r) {
  
  tau_sq <- eta / r
  tau <- sqrt(tau_sq)
  P <- N * r

  ## Generate data
  
  X <- matrix(rnorm(N * P),nrow = N, ncol = P)
  phi <- tau * rnorm(P) / sqrt(P)
  b   <- tau * rnorm(P) / sqrt(P)
  beta <- omega * phi + b
  
  ## Matrix Multiplies

  Xphi <- X %*% phi
  Xbeta <- X %*% beta

  ## Generate Data

  A <- (Xphi + rnorm(N)) |> as.numeric()
  Y <- (Xbeta + gamma * A + rnorm(N)) |> as.numeric()

  # ## Fit the A-model
  
  Vbar_inv <- t(X) %*% X + N * lambda * diag(P)
  XtA <- t(X) %*% A
  phi_hat <- solve(Vbar_inv, XtA)
  C_hat <- X %*% phi_hat
  
  ## Make more matricies

  Psi_hat     <- cbind(C_hat, X)
  H_bar_inv   <- t(Psi_hat) %*% Psi_hat + N * lambda * rbind(0,cbind(0, diag(P)))
  H_bar       <- solve(H_bar_inv)
  Qbar        <- diag(N) - Psi_hat %*% H_bar %*% t(Psi_hat)
  gamma_hat <- sum(t(A) %*% Qbar %*% Y) / sum(t(A) %*% Qbar %*% A)
  
   
  # Z <- cbind(A, C_hat, X)
  # M_bar_inv <- t(Z) %*% Z + N * lambda * rbind(0, 0, cbind(0, 0, diag(P)))
  # M_bar <- solve(M_bar_inv)
  # theta_hat <- M_bar %*% t(Z) %*% Y
  # gamma_hat_2 <- theta_hat[1]
  # 
  return(gamma_hat)
  
}
