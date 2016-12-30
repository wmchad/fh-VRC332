data {
  int<lower=1> N; // number of observations
  int<lower=1> D; // dimension (number of timepoints)
  int<lower=1> K; // number of clusters
  int<lower=1> N_subgroups;
  int<lower=1> obs_to_subgroup[N];
  vector[D] y[N]; // observations
}

transformed data {
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}

parameters {
  vector[D] mu[K]; // cluster means
}

transformed parameters {
  real<upper=0> soft_z[N_subgroups, K]; // log unnormalized clusters
  for (n in 1:N_subgroups)
    for (k in 1:K)
      soft_z[n, k] = neg_log_K;

  for (n in 1:N)
    for (k in 1:K)
      soft_z[obs_to_subgroup[n], k] = soft_z[obs_to_subgroup[n], k] -
        0.5 * 0.04 * dot_self(mu[k] - y[n]);
}

model {
  // prior
  for (k in 1:K)
    mu[k] ~ normal(0, 5);

  // likelihood
  for (n in 1:N_subgroups)
    target += log_sum_exp(soft_z[n]);
}
