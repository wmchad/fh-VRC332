data {
  int N; //Number of observations
  int<lower=1> N_subgroups;
  int<lower=1> obs_to_subgroup[N];
  real y[N];
}

parameters {
  vector[N_subgroups] beta;
  real mu_beta;
  real<lower=0> sigma_beta;
  real<lower=0> sigma;
}

model {
  real beta_obs[N];
  for (n in 1:N)
    beta_obs[n] = beta[obs_to_subgroup[n]];

  beta ~ normal(mu_beta, sigma_beta);
  mu_beta ~ normal(0, 10);
  sigma_beta ~ cauchy(0, 10);

  sigma ~ cauchy(0, 10);

  y ~ normal(beta_obs, sigma);
}
