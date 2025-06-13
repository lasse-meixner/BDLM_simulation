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

# ------------------------------------------------------------------------------- #
# Split fitting into three separate functions:
#   1) fit_BLRs_baseline: Naive, HCPH, Linero
#   2) fit_BLRs_FDML: all FDML variants
#   3) fit_BLRs_OLS_oracle: OLS and Oracle
# ------------------------------------------------------------------------------- #

# 1. Fit baseline BLR methods (Naive, HCPH, Linero) ----
fit_BLRs_baseline <- function(data) {
  invisible(capture.output({
    naive     <- BLR(y = data$Y,  XR = data$X, XF = as.matrix(data$D))
    fitted_ps <- BLR(y = data$D,  XR = data$X)
    hcph      <- BLR(y = data$Y,  XF = cbind(data$D - fitted_ps$yHat), XR = data$X)
    linero    <- BLR(y = data$Y,  XF = cbind(data$D, fitted_ps$yHat), XR = data$X)
    fit_list  <- list(
      Naive  = naive,
      HCPH   = hcph,
      Linero = linero
    )
  }))
  invisible(fit_list)
}

# 2. Fit FDML variants (FDML-Full, FDML-Split, FDML-XFit, FDML-Alt) ----
fit_BLRs_FDML <- function(data) {
  invisible(capture.output({
    # Need naive and first-stage fits
    naive     <- BLR(y = data$Y, XR = data$X, XF = as.matrix(data$D))
    fitted_ps <- BLR(y = data$D, XR = data$X)
    # 2.1 Full-sample FDML
    fitted_y_step1   <- BLR(y = data$Y, XR = data$X)
    d_res_full       <- data$D - fitted_ps$mu - data$X %*% fitted_ps$bR
    y_res_full       <- data$Y - fitted_y_step1$mu - data$X %*% fitted_y_step1$bR
    fitted_fdml_full <- lm(y_res_full ~ d_res_full)

    # 2.2 Split-sample FDML
    n               <- nrow(data$X)
    ix_step1        <- 1:floor(n/2)
    fit_d_step1     <- BLR(y = data$D[ix_step1], XR = data$X[ix_step1, ])
    fit_y_step1     <- BLR(y = data$Y[ix_step1], XR = data$X[ix_step1, ])
    ix_step2        <- (floor(n/2) + 1):n
    d_res_split     <- data$D[ix_step2] - fit_d_step1$mu - data$X[ix_step2, ] %*% fit_d_step1$bR
    y_res_split     <- data$Y[ix_step2] - fit_y_step1$mu - data$X[ix_step2, ] %*% fit_y_step1$bR
    fitted_fdml_split <- lm(y_res_split ~ d_res_split)

    # 2.3 Cross‐fitted FDML
    K              <- 10
    folds          <- sample(rep(1:K, length.out = n))
    d_residuals    <- numeric(n)
    y_residuals    <- numeric(n)
    for (k in seq_len(K)) {
      train_ix <- which(folds != k)
      test_ix  <- which(folds == k)
      ps_fit   <- BLR(y = data$D[train_ix], XR = data$X[train_ix, ])
      y_fit    <- BLR(y = data$Y[train_ix], XR = data$X[train_ix, ])
      d_pred   <- as.numeric(data$X[test_ix, ] %*% ps_fit$bR + ps_fit$mu)
      y_pred   <- as.numeric(data$X[test_ix, ] %*% y_fit$bR  + y_fit$mu)
      d_residuals[test_ix] <- data$D[test_ix] - d_pred
      y_residuals[test_ix] <- data$Y[test_ix] - y_pred
    }
    fitted_fdml_cf <- lm(y_residuals ~ d_residuals)

    # 2.4 Alternative FDML (using naive as first stage)
    y_res_alt       <- data$Y - naive$mu - data$X %*% naive$bR
    d_res_alt       <- data$D - fitted_ps$yHat
    fitted_fdml_alt <- lm(y_res_alt ~ d_res_alt)

    fit_list <- list(
      `FDML-Full`  = fitted_fdml_full,
      `FDML-Split` = fitted_fdml_split,
      `FDML-XFit`  = fitted_fdml_cf,
      `FDML-Alt`   = fitted_fdml_alt
    )
  }))
  invisible(fit_list)
}

# 3. Fit OLS and Oracle methods ----
fit_BLRs_OLS_oracle <- function(data) {
  # 3.1 OLS
  ols <- lm(data$Y ~ data$D + data$X)

  # 3.2 Oracle estimator (closed‐form posterior for alpha)
  n        <- nrow(data$X)
  p        <- ncol(data$X)
  sigma_V2 <- 1 - data$R_D2

  if (data$R_D2 == 0) {
    V_beta <- (data$R_Y2 / p) * diag(p)
    Ybar   <- data$Y
  } else {
    lambda   <- ((1 - data$R_D2) / data$R_D2) * p

    # 1st-stage posterior for gamma
    V_gamma  <- sigma_V2 * solve(crossprod(data$X) + lambda * diag(p))
    mu_gamma <- V_gamma %*% crossprod(data$X, data$D) / sigma_V2

    # compute the implied prior on beta | D, X
    a        <- data$rho * sqrt(data$R_Y2 / data$R_D2)
    mu_beta  <- a * mu_gamma
    Sigma_bg <- (data$R_Y2 * (1 - data$rho^2) / p) * diag(p)
    V_beta   <- Sigma_bg + a^2 * V_gamma
    Ybar     <- data$Y - data$X %*% mu_beta
  }

  # marginal covariance of Y
  sigma_eps2 <- 1 - data$R_Y2
  V_Y        <- sigma_eps2 * diag(n) + data$X %*% V_beta %*% t(data$X)

  # posterior for alpha
  V_alpha   <- 1 / as.numeric(crossprod(data$D, solve(V_Y, data$D)))
  alpha_hat <- V_alpha * as.numeric(crossprod(data$D, solve(V_Y, Ybar)))

  oracle <- list(
    alpha_hat = alpha_hat,
    V_alpha   = V_alpha
  )

  fit_list <- list(
    OLS    = ols,
    Oracle = oracle
  )
  invisible(fit_list)
}

fit_mvn_iw_model <- function(data) {
  reg_data <- list(
    list(y = data$Y, X = data$X),
    list(y = data$D, X = data$X)
  )
  invisible(capture.output({
    draws <- rsurGibbs(
    Data = list(regdata = reg_data),
    Prior = list(
      betabar = rep(0, ncol(data$X)*2), # prior mean (2*P)x1
      A = diag(ncol(data$X), ncol(data$X)*2), # prior precision (2*P)x(2*P)
      nu = 4, # IW prior degrees of freedom
      V = diag(1, 2, 2) # IW prior scale matrix 2x2
      ),
    Mcmc = list(R=3000, keep=1, nprint=0, burn=1000))
    }))
  # Return the transformed draws for alpha
  Sigma_draws <- draws$Sigmadraw # 1000 x (2*2)
  alpha_draws <- Sigma_draws[, 2] / Sigma_draws[, 4]
}

## Function to fit exact James-Stein IW shrinkage model ----
fit_mvn_iw_js_mat_model <- function(data) {
  # --- 0. setup & checks ---
  X    <- as.matrix(data$X)
  Tn   <- nrow(X)
  p <- ncol(X)
  if (Tn <= p) stop("Need T > p for James–Stein precision.")

  # --- 1. ordinary LS fits for gamma_hat & delta_hat + residual variances ---
  ols_fs       <- lm(data$A ~ data$X - 1)
  gamma_hat    <- coef(ols_fs)
  ols_fs_res_var <- sum(resid(ols_fs)^2)/(Tn - p)

  ols_ss       <- lm(data$Y ~ data$X - 1)
  delta_hat    <- coef(ols_ss)
  ols_ss_res_var <- sum(resid(ols_ss)^2)/(Tn - p)

  # --- 2. Ledoit–Wolf (2022) shrinkage covariance of X ---
  # Center X (subtract column means) but don’t scale by sd
  Xc <- X # scale(X, center = TRUE, scale = FALSE)
          # In the previous OLS regression, we ignored the intercept, so we don't want to center here.
          # Also, in many empirical work with ML, people usually normalize the data beforehand.

  # Compute the sample covariance matrix S_T = (1/T) Xc' Xc
  ST <- crossprod(Xc) / Tn

  # Compute the target scale ℓ_T = average variance = (1/p) tr(S_T)
  ellT <- sum(diag(ST)) / p

  # Form the deviation matrix D = S_T – ℓ_T I
  D <- ST
  diag(D) <- diag(ST) - ellT

  # Compute squared Frobenius “distance” d^2 = ||S_T – ℓ_T I||_F^2
  d2 <- sum(D^2)

  # Estimate b^2_tilde = (1/T) ∑_t || x_t x_t' – S_T ||_F^2
  b2_tilde <- 0
  for (t in 1:Tn) {
    xt       <- Xc[t, ]             # the t-th centered observation (1×p)
    outer_xt <- tcrossprod(xt)      # x_t x_t'  (p×p)
    # accumulate ||x_t x_t' – S_T||^2
    b2_tilde <- b2_tilde + sum((outer_xt - ST)^2)
  }
  b2_tilde <- b2_tilde / Tn         # average over T

  # Shrinkage intensity c_T
  cT <- min(b2_tilde / d2, 1) # always nonnegative by construction

  # Build the shrunk covariance:
  # Σ_shrunk = c_T · ℓ_T·I + (1 - c_T) · S_T
  Sigma_shrunk <- cT * ellT * diag(p) + (1 - cT) * ST

  # Recover a shrunk cross-product analogous to X'X:
  # X'X_shrunk = T · Σ_shrunk
  XtX_shrunk <- Sigma_shrunk * Tn

  # --- 3. exact James–Stein prior precision with shrunk X'X ---
  nom_mat      <- (p-2) * XtX_shrunk # same numerator for both delta and gamma
  denom_gamma        <- max(
    as.numeric(t(gamma_hat) %*% XtX_shrunk %*% gamma_hat) - ((p-2)*ols_fs_res_var),
    1e-6
  )
  gamma_precision_mat<- nom_mat / denom_gamma # denominator is always positive by construction

  denom_delta        <- max(
    as.numeric(t(delta_hat) %*% XtX_shrunk %*% delta_hat) - ((p-2)*ols_ss_res_var),
    1e-6
  )
  delta_precision_mat<- nom_mat / denom_delta

  A_shrinkage <- rbind(
    cbind(delta_precision_mat, matrix(0, p, p)),
    cbind(matrix(0, p, p), gamma_precision_mat)
  )

  # --- 4. run the same Gibbs sampler ---
  reg_data <- list(
    list(y = data$Y, X = data$X),
    list(y = data$A, X = data$X)
  )
  
  invisible(capture.output({
    draws <- rsurGibbs(
      Data  = list(regdata = reg_data),
      Prior = list(
        betabar = rep(0, 2*p),
        A       = A_shrinkage,
        nu      = 4,
        V       = diag(1, 2, 2)
      ),
      Mcmc = list(R = 3000, keep = 1, burnin = 1000))
    }))

    Sigma_draws  <- draws$Sigmadraw   # 3000 x 4
    alpha_draws  <- Sigma_draws[,2] / Sigma_draws[,4]
    return(alpha_draws)
  }



# Function to fit approximate James-Stein IW shrinkage model ----
fit_mvn_iw_js_I_model <- function(data) {

  # assert that OLS is well-defined, i.e. that n > p
  if (nrow(data$X) <= ncol(data$X)) {
    message("Number of observations must be greater than number of predictors to fit James-Stein prior precision.")
    return(NULL)
  }

  # reg_data holds data for each regression as separate list
  reg_data <- NULL
  reg_data[[1]] <- list(y = data$Y, X = data$X)
  reg_data[[2]] <- list(y = data$A, X = data$X)

  # build shrinkage prior precision matrix
  # FS OLS
  ols_fs <- lm(data$A ~ data$X - 1)
  ols_fs_res_var <- sum(ols_fs$residuals^2) / (nrow(data$X) - ncol(data$X))
  gamma_hat <- ols_fs$coefficients
  # SS OLS
  ols_ss <- lm(data$Y ~ data$X - 1)
  ols_ss_res_var <- sum(ols_ss$residuals^2) / (nrow(data$X) - ncol(data$X))
  delta_hat <- ols_ss$coefficients

  # shrinkage prior precision matrix
  nom_gamma <- (ncol(data$X) - 2)
  denom_gamma  <- max(as.numeric((t(gamma_hat) %*% gamma_hat) - ((ncol(data$X) - 2) * (ols_fs_res_var/nrow(data$X)))), 1e-6)
  gamma_precision <- if (denom_gamma <= 0) {10} else {nom_gamma / denom_gamma}

  nom_delta <- (ncol(data$X) - 2)
  denom_delta  <- max(as.numeric((t(delta_hat) %*% delta_hat) - ((ncol(data$X) - 2) * (ols_ss_res_var/nrow(data$X)))), 1e-6)
  delta_precision <- if (denom_delta <= 0) {10} else {nom_delta / denom_delta}

  # build exact JS prior shrinkage precision matrix by placing the precision values on the diagonal
  A_shrinkage <- diag(
    c(rep(delta_precision, ncol(data$X)), rep(gamma_precision, ncol(data$X)))
  )

  # Get 1000 draws
  invisible(capture.output({
    draws <- rsurGibbs(
    Data = list(regdata = reg_data),
    Prior = list(
      betabar = rep(0, ncol(data$X)*2), # prior mean (2*P)x1
      A = A_shrinkage, # prior precision (2*P)x(2*P)
      nu = 4, # IW prior degrees of freedom
      V = diag(1, 2, 2) # IW prior scale matrix 2x2
      ),
    Mcmc = list(R=1000, keep=1, nprint=0))
    }))
  # Return the transformed draws for alpha
  Sigma_draws <- draws$Sigmadraw # 1000 x (2*2)
  alpha_draws <- Sigma_draws[, 2] / Sigma_draws[, 4]
}

# ------------------------------------------------------------------------------- #
# Extraction functions for different model objects 
# ------------------------------------------------------------------------------- #
extract_results_IW <- function(draws, alpha, method_name, additional_results_info) {
  alpha_hat     <- mean(draws)
  interval      <- unname(quantile(draws, c(0.025, 0.975)))
  squared_error <- (alpha_hat - alpha)^2
  catch         <- check_interval(alpha, interval)
  interval_width<- interval[2] - interval[1]
  LCL <- interval[1]; UCL <- interval[2]
  table <- data.frame(
    alpha_hat      = alpha_hat,
    squared_error  = squared_error,
    LCL            = LCL,
    UCL            = UCL,
    catch          = catch,
    interval_width = interval_width,
    Method         = method_name
  )
  for (nm in names(additional_results_info)) {
    table[[nm]] <- additional_results_info[[nm]]
  }
  table
}

extract_results_blr <- function(fit, alpha, method_name, additional_results_info) {
  alpha_hat     <- fit$bF[1]
  interval      <- get_interval(fit)
  squared_error <- (alpha_hat - alpha)^2
  catch         <- check_interval(alpha, interval)
  interval_width<- interval[2] - interval[1]
  LCL <- interval[1]; UCL <- interval[2]
  table <- data.frame(
    alpha_hat      = alpha_hat,
    squared_error  = squared_error,
    LCL            = LCL,
    UCL            = UCL,
    catch          = catch,
    interval_width = interval_width,
    Method         = method_name
  )
  for (nm in names(additional_results_info)) {
    table[[nm]] <- additional_results_info[[nm]]
  }
  table
}

extract_results_lm <- function(fit, alpha, method_name, additional_results_info) {
  alpha_hat     <- summary(fit)$coefficients[2, 1]
  interval      <- unname(confint(fit)[2, ])
  squared_error <- (alpha_hat - alpha)^2
  catch         <- check_interval(alpha, interval)
  interval_width<- interval[2] - interval[1]
  LCL <- interval[1]; UCL <- interval[2]
  table <- data.frame(
    alpha_hat      = alpha_hat,
    squared_error  = squared_error,
    LCL            = LCL,
    UCL            = UCL,
    catch          = catch,
    interval_width = interval_width,
    Method         = method_name
  )
  for (nm in names(additional_results_info)) {
    table[[nm]] <- additional_results_info[[nm]]
  }
  table
}

extract_results_oracle <- function(fit, alpha, method_name, additional_results_info) {
  alpha_hat     <- fit$alpha_hat
  se_alpha      <- sqrt(fit$V_alpha)
  interval      <- alpha_hat + c(-1, 1) * 1.96 * se_alpha
  squared_error <- (alpha_hat - alpha)^2
  catch         <- check_interval(alpha, interval)
  interval_width<- interval[2] - interval[1]
  LCL <- interval[1]; UCL <- interval[2]
  table <- data.frame(
    alpha_hat      = alpha_hat,
    squared_error  = squared_error,
    LCL            = LCL,
    UCL            = UCL,
    catch          = catch,
    interval_width = interval_width,
    Method         = method_name
  )
  for (nm in names(additional_results_info)) {
    table[[nm]] <- additional_results_info[[nm]]
  }
  table
}

# ------------------------------------------------------------------------------- #
# sim_iter wrapper functions for different BLR methods: generate data, fit model, extract results
# ------------------------------------------------------------------------------- #

# 1. Simulate baseline BLR methods (Naive, HCPH, Linero) ----
sim_iter_BLRs_baseline <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data       <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit_list   <- fit_BLRs_baseline(data)
  selected   <- c("Naive", "HCPH", "Linero")
  results    <- lapply(selected, function(model_name) {
    fit_obj <- fit_list[[model_name]]
    extract_results_blr(
      fit = fit_obj,
      alpha = data$alpha,
      method_name = model_name,
      additional_results_info = list(
        R_Y2  = R_Y2,
        R_D2  = R_D2,
        rho   = rho,
        alpha = alpha,
        n     = n,
        p     = p
      )
    )
  })
  do.call(rbind, results)
}

# 2. Simulate FDML variants (FDML-Full, FDML-Split, FDML-XFit, FDML-Alt) ----
sim_iter_BLRs_FDML <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data     <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit_list <- fit_BLRs_FDML(data)
  selected <- c("FDML-Full", "FDML-Split", "FDML-XFit", "FDML-Alt")
  results  <- lapply(selected, function(model_name) {
    fit_obj <- fit_list[[model_name]]
    extract_results_lm(
      fit = fit_obj,
      alpha = data$alpha,
      method_name = model_name,
      additional_results_info = list(
        R_Y2  = R_Y2,
        R_D2  = R_D2,
        rho   = rho,
        alpha = alpha,
        n     = n,
        p     = p
      )
    )
  })
  do.call(rbind, results)
}

# 3. Simulate OLS and Oracle methods ----
sim_iter_BLRs_OLS_oracle <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data     <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  fit_list <- fit_BLRs_OLS_oracle(data)

  # OLS result
  ols_res <- extract_results_lm(
    fit = fit_list[["OLS"]],
    alpha = data$alpha,
    method_name = "OLS",
    additional_results_info = list(
      R_Y2  = R_Y2,
      R_D2  = R_D2,
      rho   = rho,
      alpha = alpha,
      n     = n,
      p     = p
    )
  )

  # Oracle result
  oracle_res <- extract_results_oracle(
    fit = fit_list[["Oracle"]],
    alpha = data$alpha,
    method_name = "Oracle",
    additional_results_info = list(
      R_Y2  = R_Y2,
      R_D2  = R_D2,
      rho   = rho,
      alpha = alpha,
      n     = n,
      p     = p
    )
  )

  rbind(ols_res, oracle_res)
}

# Main simulation function for BDML-IW for a given setting (unchanged) ----
sim_iter_bdml_iw <- function(n, p, R_Y2, R_D2, rho, alpha, seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  data   <- generate_data(n, p, R_Y2, R_D2, rho, alpha)
  draws  <- fit_mvn_iw_model(data)
  extract_results_IW(
    draws = draws,
    alpha = data$alpha,
    method_name = "BDML-IW",
    additional_results_info = list(
      R_Y2  = R_Y2,
      R_D2  = R_D2,
      rho   = rho,
      alpha = alpha,
      n     = n,
      p     = p
    )
  )
}
