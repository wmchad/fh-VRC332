data {
  int<lower=1> N; // number of observations
  int<lower=1> D; // dimension (number of timepoints)
  int<lower=1> K; // number of clusters
  int<lower=1> N_subgroups;
}

transformed data {
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}

parameters {
  vector[D] mu[K]; // cluster means
}

transformed parameters {
  real<upper=0> soft_z[N, K]; // log unnormalized clusters
  for (n in 1:N)
    for (k in 1:K)
      soft_z[n, k] = neg_log_K - 0.5 * dot_self(mu[k] - y[n]);
}

model {
  // prior
  for (k in 1:K)
    mu[k] ~ normal(0, 1);

  // likelihood
  for (n in 1:N)
    target += log_sum_exp(soft_z[n]);
}
