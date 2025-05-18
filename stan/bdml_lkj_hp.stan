data {
  int<lower=1> p;
  int<lower=0> n;
  array[n] vector[p] X;
  array[n] vector[2] Y;
}
parameters {
  matrix[2, p] beta;
  cholesky_factor_corr[2] L_Omega;
  vector<lower=0>[2] L_sigma;
  vector<lower=0>[2] sigma_beta;  // Different standard deviations for beta rows
}
model {
  array[n] vector[2] mu;
  matrix[2, 2] L_Sigma;

  for (i in 1:n) {
    mu[n] = beta * X[n];
  }

  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  // Separate priors for each row of beta
  beta[1] ~ normal(0, sigma_beta[1]);
  beta[2] ~ normal(0, sigma_beta[2]);

  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);

  Y ~ multi_normal_cholesky(mu, L_Sigma);

  // Inv-Gamma(2,2) for sigma_beta^2
  target += inv_gamma_lpdf( square(sigma_beta) | 2, 2 )
          + log(2 * sigma_beta);
}
generated quantities {
  real alpha;
  alpha = L_Omega[2,1] * L_sigma[1] / L_sigma[2];
}


