//Using the median estimates from perf-predict.stan and a vector of inputs predict
//out of sample observations for use in sensitivity analysis

data {
  int<lower=1> N1;
  int<lower=1> D1;
  vector[D1] x1[N1];
  vector[N1] y1;
  int<lower=1> N2;
  int<lower=1> D2;
  vector[D2] x2[N2];
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  vector[D1] x[N];
  matrix[N, N] K;
  matrix[N, N] L_K;
  for (n1 in 1:N1) x[n1] = x1[n1];
  for (n2 in 1:N2) x[N1 + n2] = x2[n2];
  
    K = cov_exp_quad(x, alpha, rho);
  
    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;
    
    L_K = cholesky_decompose(K);
}

parameters {
  vector[N] eta;
}

transformed parameters {
  vector[N] f;
    f = L_K * eta;
}

model {
  eta ~ normal(0, 1);
  y1 ~ normal(f[1:N1], sigma);
}

generated quantities {
  vector[N2] y2;
  for (n2 in 1:N2)
    y2[n2] = normal_rng(f[N1 + n2], sigma);
}

