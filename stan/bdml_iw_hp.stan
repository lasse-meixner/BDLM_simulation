data {
  int<lower=1> p;
  int<lower=0> n;
  array[n] vector[p] X;
  array[n] vector[2] Y;
}
parameters {
  matrix[2, p] beta;
  cov_matrix[2] Sigma;               // Covariance matrix for Y
  vector<lower=0>[2] sigma_beta;     // Different standard deviations for beta rows
}
model {
  array[n] vector[2] mu;
  
  // Linear predictor
  for (i in 1:n)
    mu[i] = beta * X[i];
  
  // Priors for beta coefficients
  beta[1] ~ normal(0, sigma_beta[1]);
  beta[2] ~ normal(0, sigma_beta[2]);
  
  // Inverse Wishart prior for the covariance matrix
  // Here, nu = 4 and the scale matrix is diagonal with 1 on the diagonal.
  Sigma ~ inv_wishart(4, diag_matrix(rep_vector(1, 2)));
  
  // Likelihood using the full covariance matrix
  Y ~ multi_normal(mu, Sigma);

  // Inv-Gamma(2,2) for sigma_beta^2
  target += inv_gamma_lpdf( square(sigma_beta) | 2, 2 )
          + log(2 * sigma_beta);
}
generated quantities {
  // Compute the correlation from the covariance matrix
  real alpha;
  alpha = Sigma[1,2] / Sigma[2,2];
}
