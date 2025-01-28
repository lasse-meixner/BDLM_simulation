functions {
  /* Efficient computation of the R2D2 prior
   * Args:
   *   z: standardized population-level coefficients
   *   phi: local weight parameters
   *   tau2: global scale parameter
   * Returns:
   *   population-level coefficients following the R2D2 prior
   */
  vector R2D2(vector z, vector phi, real tau2) {
    return z .* sqrt(phi * tau2);
  }
}
data {
  int<lower=1> J;  // number of predictors
  int<lower=0> N;  // number of observations
  array[N] vector[J] x;  // predictors
  array[N] vector[2] y;  // 2-dimensional response
  
  // R2D2 prior parameter
  real<lower=0> b;
}
parameters {
  // Coefficients for each dimension
  vector[J] zb1;  // standardized coefficients for dim 1
  vector[J] zb2;  // standardized coefficients for dim 2
  
  // R2D2 parameters for each dimension
  simplex[J] R2D2_phi1;  // local weights for dim 1
  simplex[J] R2D2_phi2;  // local weights for dim 2
  real<lower=0,upper=1> R2D2_R2_1;  // R2 parameter for dim 1
  real<lower=0,upper=1> R2D2_R2_2;  // R2 parameter for dim 2
  
  // Error covariance matrix
  cholesky_factor_corr[2] L_Omega;
  vector<lower=0>[2] L_sigma;
}
transformed parameters {
  real<lower=0> a_pi;  
  real<lower=0> a;  
  vector<lower=0>[J] R2D2_alpha1;  // concentration vector for dim 1
  vector<lower=0>[J] R2D2_alpha2;  // concentration vector for dim 2
  real<lower=0> R2D2_mean_R2_1;    // mean of R2 prior for dim 1
  real<lower=0> R2D2_mean_R2_2;    // mean of R2 prior for dim 2
  real<lower=0> R2D2_prec_R2_1;    // precision of R2 prior for dim 1
  real<lower=0> R2D2_prec_R2_2;    // precision of R2 prior for dim 2
  
  vector[J] b1;  // coefficients for dim 1
  vector[J] b2;  // coefficients for dim 2
  real R2D2_tau2_1;  // global scale for dim 1
  real R2D2_tau2_2;  // global scale for dim 2
  matrix[2, 2] L_Sigma;
  
  // Initialize R2D2 hyperparameters
  a_pi = 1/(J^(b/2) * N^(b/2) * log(N));
  a = a_pi*J;
  R2D2_alpha1 = rep_vector(a_pi, J);
  R2D2_alpha2 = rep_vector(a_pi, J);
  R2D2_mean_R2_1 = a/(a+b);
  R2D2_mean_R2_2 = a/(a+b);
  R2D2_prec_R2_1 = a + b;
  R2D2_prec_R2_2 = a + b;
  
  // Compute global scale parameters
  R2D2_tau2_1 = L_sigma[1]^2 * R2D2_R2_1 / (1 - R2D2_R2_1);
  R2D2_tau2_2 = L_sigma[2]^2 * R2D2_R2_2 / (1 - R2D2_R2_2);
  
  // Compute actual regression coefficients
  b1 = R2D2(zb1, R2D2_phi1, R2D2_tau2_1);
  b2 = R2D2(zb2, R2D2_phi2, R2D2_tau2_2);
  
  // Construct L_Sigma
  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
}
model {
  array[N] vector[2] mu;
  // Compute mean predictions
  for (n in 1:N) {
    mu[n][1] = dot_product(x[n], b1);
    mu[n][2] = dot_product(x[n], b2);
  }
  
  // Priors for R2D2 parameters
  R2D2_R2_1 ~ beta(R2D2_mean_R2_1 * R2D2_prec_R2_1, 
                   (1 - R2D2_mean_R2_1) * R2D2_prec_R2_1);
  R2D2_R2_2 ~ beta(R2D2_mean_R2_2 * R2D2_prec_R2_2, 
                   (1 - R2D2_mean_R2_2) * R2D2_prec_R2_2);
  
  R2D2_phi1 ~ dirichlet(R2D2_alpha1);
  R2D2_phi2 ~ dirichlet(R2D2_alpha2);
  
  zb1 ~ std_normal();
  zb2 ~ std_normal();
  
  // Priors for covariance parameters
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);
  
  // Likelihood
  y ~ multi_normal_cholesky(mu, L_Sigma);
}
generated quantities {
  real alpha;
  alpha = L_Omega[2,1] * L_sigma[1] / L_sigma[2];
}
