data {
  int<lower=1> k;
  int<lower=1> p;
  int<lower=0> n;
  array[n] vector[p] X;
  array[n] vector[k] Y;
}
parameters {
  matrix[k, p] beta;
  cholesky_factor_corr[k] L_Omega;
  vector<lower=0>[k] L_sigma;
}
model {
  array[n] vector[k] mu;
  matrix[k, k] L_Sigma;

  for (i in 1:n) {
    mu[i] = beta * X[i];

  }

  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  to_vector(beta) ~ normal(0, 5);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);

  Y ~ multi_normal_cholesky(mu, L_Sigma);
}
generated quantities {
  real alpha;
  alpha = L_Omega[2,1] * L_sigma[1] / L_sigma[2];
}


