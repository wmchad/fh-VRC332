data {
  int<lower=1> N; //Number of observations
  int<lower=1> N_subgroups;
  int<lower=1> obs_to_subgroup[N];
  int<lower=1> N_t;
  int<lower=1> obs_to_t[N];
  real y[N];
}

parameters {
  // real<lower=0, upper=1> omega;
  vector<lower=0, upper=1>[N_subgroups] omega;
  vector[N_subgroups] beta;
  real<lower=0> sigma_beta;
  vector[N_t] beta_t;
  real<lower=0> sigma_beta_t;
  real<lower=0> sigma;
}

model {
  real beta_obs[N];
  for (n in 1:N) {
    beta_obs[n] = omega[obs_to_subgroup[n]] * (beta[obs_to_subgroup[n]] + beta_t[obs_to_t[n]]);
  }
  
  beta ~ normal(0, sigma_beta);
  beta_t ~ normal(0, sigma_beta_t);
  sigma_beta ~ cauchy(0, 10);
  sigma_beta_t ~ cauchy(0, 1);

  sigma ~ cauchy(0, 10);
  omega ~ beta(0.1, 0.1);

  // Trying to get a mixture of 0 and a normal
  // for (n in 1:N) {
  //   target += log_sum_exp(log(omega[n]) + normal_lpdf(y[n] | 0, sigma),
  //                         log(1-omega[n]) + normal_lpdf(y[n] | beta_obs[n], sigma));
  // }
  // target += log_sum_exp(log(omega) + normal_lpdf(y | 0, 0.0001),
  //                       log(1-omega) + normal_lpdf(y | beta_obs, sigma));
  y ~ normal(beta_obs, sigma);
}
