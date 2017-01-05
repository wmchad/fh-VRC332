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
  vector<lower=0>[D] sigma[K]; // cluster variances (sdev)
}

transformed parameters {
  real<upper=0> soft_z[N_subgroups, K]; // log unnormalized clusters
  for (n in 1:N_subgroups)
    for (k in 1:K)
      soft_z[n, k] = neg_log_K;

  for (n in 1:N)
    for (k in 1:K)
      soft_z[obs_to_subgroup[n], k] = soft_z[obs_to_subgroup[n], k] +
        normal_lpdf(y[n] | mu[k], sigma[k]);
}

model {
  // prior
  for (k in 1:K) {
    mu[k] ~ normal(0, 2);
    sigma[k] ~ cauchy(0, 10);
  }
  
  // likelihood
  for (n in 1:N_subgroups)
    target += log_sum_exp(soft_z[n]);
}
