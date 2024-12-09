theoretical_numerator <- function(N,P,lambda,eta,tau_sq) {
  r <- P / N
  N - P + lambda * eta / 
  (lambda - eta) * (m(lambda,r) * (N * tau_sq * lambda - P) - 
                    m(eta,r) * (N * tau_sq * eta - P))
}

theoretical_denominator <- function(N,P,lambda,eta,tau_sq) {
  r <- P / N
  N - P + lambda * (P * m(lambda,r) + N * tau_sq) - 
  lambda^2 * tau_sq * N * m(lambda,r)
}

ridge_test_2 <- function(lambda, eta, gamma, omega, N, r) {
  P <- N * r
  tau_sq <- r / eta
  tau <- sqrt(tau_sq)
  
  X <- matrix(rnorm(N * P),nrow = N, ncol = P)
  phi <- tau * rnorm(P) / sqrt(P)
  beta <- rnorm(P, phi * omega, tau / sqrt(P))
  
  ## Matrix Multiplies
  
  Xphi <- X %*% phi
  Xbeta <- X %*% beta
  
  ## Generate Data
  
  A <- (Xphi + rnorm(N)) |> as.numeric()
  Y <- (Xbeta + gamma * A + rnorm(N)) |> as.numeric()
  
  ## Compute gamma tilde
  
  AtA <- sum(A^2)
  AtY <- sum(A * Y)
  XtA <- t(X) %*% A
  XtY <- t(X) %*% Y
  Vbar <- t(X) %*% X + diag(P) * N * lambda
  
  numerator <- AtY - t(XtA) %*% solve(Vbar, XtY)
  denominator <- AtA - t(XtA) %*% solve(Vbar, XtA)
  
  gamma_tilde <- numerator / denominator
  return(gamma_tilde)
  
  
}

ridge_test <- function(tau = 3, lambda = NULL, gamma = 3, omega = 1, N = 300, P = 300) {
  
  tau_sq <- tau^2
  if(is.null(lambda)) lambda <- P / N / tau_sq
  
  X <- matrix(rnorm(N * P),nrow = N, ncol = P)
  phi <- tau * rnorm(P) / sqrt(P)
  # phi <- rep(1 , P) * tau / sqrt(P)
  beta <- rnorm(P, phi * omega, tau / sqrt(P))
  
  ## Matrix Multiplies
  
  Xphi <- X %*% phi
  Xbeta <- X %*% beta
  
  ## Generate Data
  
  A <- (Xphi + rnorm(N)) |> as.numeric()
  Y <- (Xbeta + gamma * A + rnorm(N)) |> as.numeric()
  
  ## More matrix multiplies
  
  XtX <- t(X) %*% X
  Q <- X %*% solve(N * lambda * diag(P) + XtX, t(X))
  Qbar <- diag(N) - Q
  
  ## Computing full ridge
  
  M <- rbind(t(A) %*% X, t(X) %*% X + N * lambda * diag(P))
  M <- cbind(rbind(sum(A^2), t(X) %*% A), M)
  XtY <- t(cbind(A, X)) %*% Y
  ridge_est <- solve(M, XtY)
  # gamma_tilde_full <- ridge_est[1]
  # AX <- cbind(A,X)
  # M2 <- t(AX) %*% AX
  
  ## Other Stuff
  
  Delta_star <- sum(phi * beta) / (1 + sum(phi^2))
  
  R_hat <- Qbar %*% A
  mu_hat <- Q %*% A
  
  bias_hat <- omega * (1 - sum(R_hat * R_hat) / sum(R_hat * A))
  gamma_tilde <- as.numeric(t(A) %*% Qbar %*% Y) / as.numeric(t(A) %*% Qbar %*% A)
  bias <- gamma_tilde - gamma
  
  rho <- P/N
  m_lambda <- m(lambda, rho)
  theoretical_bias <- omega * rho * (1 - lambda * m_lambda)
  # theo_denom <- N * tau_sq * (lambda - lambda^2 * m_lambda) + P * lambda * m_lambda
  theo_denom <- N * tau_sq + N - N * tau_sq * (1 - lambda + lambda^2 * m_lambda) - P * (1 - lambda * m_lambda)
  
  # var_Y <- tau_sq * omega^2 / P * X %*% t(X) + diag(N)
  K_eta <- t(X) %*% X + P/tau_sq * diag(P)
  var_Y <- tau_sq / P * X %*% t(X) + diag(N) + omega^2 * X %*% solve(K_eta) %*% t(X)
  var_A <- tau_sq / P * X %*% t(X) + diag(N)
  test_var <- sum(diag(Qbar %*% var_Y %*% Qbar %*% var_A)) / theoretical_denominator(N,P,lambda,eta,tau_sq)^2
  test_var_2 <- sum(diag(Qbar %*% var_Y %*% Qbar %*% var_A)) / sum(diag((Qbar %*% (tau_sq / P * X %*% t(X) + diag(N)))))^2
  missing_piece <- omega * t(A) %*% Qbar %*% X %*% solve(K_eta, t(X)) %*% A / (t(A) %*% Qbar %*% A)
  # mse = mean(abs(Y - mu_hat)^2)
  mse = mean(abs(X %*% beta + gamma_tilde * A - X %*% ridge_est[-1] - gamma * A)^2)
  
  NN <- t(A) %*% Qbar %*% Y
  DD <- t(A) %*% Qbar %*% A
  
  out <- tibble(
    gamma_tilde = gamma_tilde,
    error = gamma_tilde - gamma,
    rb_error = bias_hat,
    abias = theoretical_bias,
    denom = as.numeric(t(A) %*% Qbar %*% A),
    tr_V = mtrace(solve(t(X) %*% X / N + lambda * diag(P))),
    NN = NN,
    DD = DD,
    # num = as.numeric(t(A) %*% Qbar %*% A),
    theo_denom = theo_denom,
    mse = mse,
    # T1 = sum(diag(tau_sq / P * X %*% t(X))), 
    # T2 = N,
    # T3 = sum(diag(tau_sq / P * Q %*% X %*% t(X))),
    # T4 = sum(diag(Q)),
    test_var = test_var,
    test_var_2 = test_var_2,
    missing_piece = missing_piece)
  # ) |> mutate(ham = T1 + T2 - T3 - T4)
  
  return(out)
  
  # return(c(bias_hat, bias, omega * P / N, risk = sum(R_hat^2), 
  #          inner = sum(mu_hat * R_hat),
  #          denom = sum(t(A) %*% Qbar %*% A),
  #          tinner = N,
  #          tr_inv = tr_inv, abias = theoretical_bias,
  #          Delta_star = Delta_star,
  #          gamma_tilde = gamma_tilde,
  #          Delta_approx = Delta_star * min(P/N,1)))
  
}
