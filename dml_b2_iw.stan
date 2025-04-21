data {
  int<lower=1> J;
  int<lower=0> N;
  array[N] vector[J] x;
  array[N] vector[2] y;
}
parameters {
  matrix[2, J] beta;
  cov_matrix[2] Sigma;               // Covariance matrix for y
  vector<lower=0>[2] sigma_beta;     // Different standard deviations for beta rows
}
model {
  array[N] vector[2] mu;
  
  // Linear predictor
  for (n in 1:N)
    mu[n] = beta * x[n];
  
  // Priors for beta coefficients
  beta[1] ~ normal(0, sigma_beta[1]);
  beta[2] ~ normal(0, sigma_beta[2]);
  sigma_beta ~ inv_gamma(2, 2);
  
  // Inverse Wishart prior for the covariance matrix
  // Here, nu = 4 and the scale matrix is diagonal with 1 on the diagonal.
  Sigma ~ inv_wishart(4, diag_matrix(rep_vector(1, 2)));
  
  // Likelihood using the full covariance matrix
  y ~ multi_normal(mu, Sigma);
}
generated quantities {
  // Compute the correlation from the covariance matrix
  real alpha;
  alpha = Sigma[1,2] / Sigma[2,2];
}
