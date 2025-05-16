library(tidyverse)
library(BLR)
library(bayesm) # for multivariate Inverse Wishart model

# Loading source for data generation
source("generate_data_source.R")

## Load all R scripts in the src/ directory ----
load_src_files <- function(src_path = "../src/") {
  source_files <- list.files(src_path, "*.R$")
  map(paste0(src_path, source_files), source)
}
load_src_files()

# Function to fit Hahn model ----
fit_BLRs <- function(data) {
  # Note: since BLR doesn't allow control of low-level cat statements (typical bs R package SE), I suppress this manually

  invisible(capture.output({ # Suppress output from BLR calls
    # 1. Hahn & Linero
    naive <- BLR(y = data$Y, XR = data$X, XF = as.matrix(data$D))
    fitted_ps <- BLR(y = data$D, XR = data$X)
    hcph <- BLR(y = data$Y, XF = cbind(data$D - fitted_ps$yHat), XR = data$X)
    linero <- BLR(y = data$Y, XF = cbind(data$D, fitted_ps$yHat), XR = data$X)

    # 2. FDML
    fitted_d_step1 <- fitted_ps
    fitted_y_step1 <- BLR(y = data$Y, XR = data$X)
      
    d_res_step2 <- data$D - fitted_d_step1$mu - data$X %*% fitted_d_step1$bR
    y_res_step2 <- data$Y - fitted_y_step1$mu - data$X %*% fitted_y_step1$bR

    fitted_fdml_full <- lm(y_res_step2 ~ d_res_step2)

    n <- dim(data$X)[1] # number of observations
    ix_step1 <- 1:(floor(n/2))
    fitted_d_step1 <- BLR(y = data$D[ix_step1], XR = data$X[ix_step1,])
    fitted_y_step1 <- BLR(y = data$Y[ix_step1], XR = data$X[ix_step1,])

    ix_step2 <- (floor(n/2) + 1):n
    d_res_step2 <- data$D[ix_step2] - fitted_d_step1$mu - data$X[ix_step2,] %*% fitted_d_step1$bR
    y_res_step2 <- data$Y[ix_step2] - fitted_y_step1$mu - data$X[ix_step2,] %*% fitted_y_step1$bR
    
    # 2.3 Cross‐fitted FDML
    K <- 10
    n <- nrow(data$X)
    folds <- sample(rep(1:K, length.out = n))
    d_residuals <- numeric(n)
    y_residuals <- numeric(n)
    
    for (k in seq_len(K)) {
      train_ix <- which(folds != k)
      test_ix  <- which(folds == k)
      
      # first‐stage fits on training data
      ps_fit <- BLR(y = data$D[train_ix], XR = data$X[train_ix, ])
      y_fit  <- BLR(y = data$Y[train_ix], XR = data$X[train_ix, ])
      
      # predict on hold‐out fold
      d_pred <- as.numeric(data$X[test_ix, ] %*% ps_fit$bR + ps_fit$mu)
      y_pred <- as.numeric(data$X[test_ix, ] %*% y_fit$bR  + y_fit$mu)
      
      # residualize
      d_residuals[test_ix] <- data$D[test_ix] - d_pred
      y_residuals[test_ix] <- data$Y[test_ix] - y_pred
    }
    
  })) # end of silenced BLR calls

  fitted_fdml_cf <- lm(y_residuals ~ d_residuals)
  fitted_fdml_split <- lm(y_res_step2 ~ d_res_step2)
  
  # 2.4 Alternative FDML
  y_res_alt <- data$Y - naive$mu - data$X %*% naive$bR
  d_res_alt <- data$D - fitted_ps$yHat
  fitted_fdml_alt <- lm(y_res_alt ~ d_res_alt)
  
  # 3. OLS
  ols   <- lm(data$Y ~ data$D + data$X)
  
  # 4. Oracle estimator (closed‐form posterior of alpha)
  n <- nrow(data$X)
  p <- ncol(data$X)
  sigma_V2 <- 1 - data$R_D2
  lambda   <- ((1 - data$R_D2)/data$R_D2) * p
  
  # 1st‐stage posterior for gamma
  V_gamma   <- sigma_V2 * solve(crossprod(data$X) + lambda * diag(p))
  mu_gamma <- V_gamma %*% crossprod(data$X, data$D) / sigma_V2
  
  # compute the implied prior on beta|D,X
  a        <- data$rho * sqrt(R_Y2/data$R_D2)
  mu_beta  <- a * mu_gamma
  Sigma_bg <- (data$R_Y2*(1 - data$rho^2)/p) * diag(p)
  V_beta   <- Sigma_bg + a^2 * V_gamma
  
  # marginal covariance of Y
  sigma_eps2 <- 1 - data$R_Y2
  V_Y <- sigma_eps2 * diag(n) + data$X %*% V_beta %*% t(data$X)
  
  # posterior for alpha
  V_alpha   <- 1 / as.numeric(crossprod(data$D, solve(V_Y, data$D)))
  alpha_hat <- V_alpha * as.numeric(crossprod(data$D, solve(V_Y, data$Y - data$X %*% mu_beta)))
  
  oracle <- list(
    alpha_hat = alpha_hat,       # so extract_results_blr() can grab it as fit$bF[1]
    V_alpha   = V_alpha         # you can write a tiny get_interval_oracle() if you want CIs
  )
  
  # Ensure the return object is not printed
  invisible(list(Naive = naive, HCPH = hcph, Linero = linero, 
                 "FDML-Full" = fitted_fdml_full, "FDML-Split" = fitted_fdml_split, 
                 "FDML-XFit" = fitted_fdml_cf, "FDML-Alt" = fitted_fdml_alt,
                 OLS = ols, Oracle = oracle))
}

fit_mvn_iw_model <- function(data) {

  # reg_data holds data for each regression as separate list
  reg_data <- NULL
  reg_data[[1]] <- list(y = data$Y, X = data$X)
  reg_data[[2]] <- list(y = data$D, X = data$X)

  # Get 1000 draws
  invisible(capture.output({
    draws <- rsurGibbs(
    Data = list(regdata = reg_data),
    Prior = list(
      betabar = rep(0, ncol(data$X)*2), # prior mean (2*p)x1
      A = diag(1/ncol(data$X), ncol(data$X)*2), # prior precision (2*p)x(2*p)
      nu = 4, # IW prior degrees of freedom
      V = diag(1, 2, 2) # IW prior scale matrix 2x2
      ),
    Mcmc = list(R=12000, keep=1, nprint=0))
    }))
  # Return the transformed draws for alpha
  Sigma_draws <- draws$Sigmadraw[-(1:2000), , drop = FALSE] # R x (2*2), remove burn-in
  alpha_draws <- Sigma_draws[, 2] / Sigma_draws[, 4]
}

# Function to extract results for IW model
extract_results_IW <- function(draws, alpha, method_name, additional_results_info) {
  alpha_hat <- mean(draws)
  interval <- unname(quantile(draws, c(0.025, 0.975)))
  squared_error <- (alpha_hat - alpha)^2
  catch <- check_interval(alpha, interval)
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]
  
  table <- data.frame(
    alpha_hat = alpha_hat, 
    squared_error = squared_error,
    LCL = LCL,
    UCL = UCL,
    catch = catch,
    interval_width = interval_width,
    Method = method_name
  )

  # Append additional_results_info (setting parameters) to pass through for downstream analysis
  for (name in names(additional_results_info)) {
    table[[name]] <- additional_results_info[[name]]
  }
  table
}

# Function to extract results for Hahn and Linero models from BLR object ----
extract_results_blr <- function(fit, alpha, method_name, additional_results_info) {
  alpha_hat <- fit$bF[1]
  interval <- get_interval(fit)
  squared_error <- (alpha_hat - alpha)^2
  catch <- check_interval(alpha, interval)
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]
  
  table <- data.frame(
    alpha_hat = alpha_hat, 
    squared_error = squared_error,
    LCL = LCL,
    UCL = UCL,
    catch = catch,
    interval_width = interval_width,
    Method = method_name
  )

  # Append additional_results_info (setting parameters) to pass through for downstream analysis
  for (name in names(additional_results_info)) {
    table[[name]] <- additional_results_info[[name]]
  }
  table
}

# Function to extract results for FDML model from lm object ----
extract_results_lm <- function(fit, alpha, method_name, additional_results_info) {
  alpha_hat <- summary(fit)$coefficients[2, 1]
  interval <- unname(confint(fit)[2,])
  squared_error <- (alpha_hat - alpha)^2
  catch <- check_interval(alpha, interval)
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]
  
  table <- data.frame(
    alpha_hat = alpha_hat, 
    squared_error = squared_error,
    LCL = LCL,
    UCL = UCL,
    catch = catch,
    interval_width = interval_width,
    Method = method_name
  )

  # Append additional_results_info (setting parameters) to pass through for downstream analysis
  for (name in names(additional_results_info)) {
    table[[name]] <- additional_results_info[[name]]
  }
  table
}

# Function to extract results for oracle ----
extract_results_oracle <- function(fit, alpha, method_name, info) {
  alpha_hat  <- fit$alpha_hat
  se_alpha   <- sqrt(fit$V_alpha)
  interval   <- alpha_hat + c(-1,1)*1.96*se_alpha
  squared_error <- (alpha_hat - alpha)^2
  catch <- check_interval(alpha, interval)
  interval_width <- interval[2] - interval[1]
  LCL <- interval[1]
  UCL <- interval[2]
  
  table <- data.frame(
    alpha_hat    = alpha_hat,
    squared_error = squared_error,
    LCL          = LCL,
    UCL          = UCL,
    catch        = catch,
    interval_width = interval_width,
    Method       = method_name
  )
  
  # Append additional_results_info (setting parameters) to pass through for downstream analysis
  for (name in names(additional_results_info)) {
    table[[name]] <- additional_results_info[[name]]
  }
  table
}

# Main simulation function for BRLs for a given setting ----
sim_iter_BLRs <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit_BRL <- fit_BLRs(data)
  
  # extract results for each BRL model
  BRLs_extraction <- lapply(names(fit_BRL), function(model_name) {
    if (model_name %in% c("FDML-Full", "FDML-Split", "FDML-XFit", "FDML-Alt", "OLS")) {
      extract_results_lm(fit_BRL[[model_name]], data$alpha, model_name, additional_results_info = list(R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, n = n, p = p))
    } else if (model_name == "Oracle") {
      extract_results_oracle(fit_BRL[[model_name]], data$alpha, model_name, additional_results_info = list(R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, n = n, p = p))
    } else {
      extract_results_blr(fit_BRL[[model_name]], data$alpha, model_name, additional_results_info = list(R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, n = n, p = p))
    }
  })

  # combine results
  do.call(rbind, BRLs_extraction)
}

# Main simulation function for BDML-IW for a given setting ----
sim_iter_BDML_iw <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit_IW <- fit_mvn_iw_model(data)
  
  # extract results for IW model
  IW_extraction <- extract_results_IW(fit_IW, data$alpha, "BDML-IW", additional_results_info = list(R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, n = n, p = p))

  # combine results
  IW_extraction
}
