data {
  int<lower=1> J;
  int<lower=0> N;
  array[N] vector[J] x;
  array[N] vector[2] y;
}
parameters {
  matrix[2, J] beta;
  cholesky_factor_corr[2] L_Omega;
  vector<lower=0>[2] L_sigma;
  vector<lower=0>[2] sigma_beta;  // Different standard deviations for beta rows
}
model {
  array[N] vector[2] mu;
  matrix[2, 2] L_Sigma;

  for (n in 1:N) {
    mu[n] = beta * x[n];
  }

  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  // Separate priors for each row of beta
  beta[1] ~ normal(0, sigma_beta[1]);
  beta[2] ~ normal(0, sigma_beta[2]);
  sigma_beta ~ inv_gamma(2, 2);

  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);

  y ~ multi_normal_cholesky(mu, L_Sigma);
}
generated quantities {
  real alpha;
  alpha = L_Omega[2,1] * L_sigma[1] / L_sigma[2];
}


