data {
  int<lower=1> N;
  vector[N] delta;
  matrix[N, N] Sigma;
}
transformed data {
  real m = mean(delta);
  real s = sd(delta);
  vector[N] z = (delta - m) / s;
  matrix[N, N] Sigma_z = Sigma / square(s);
}
parameters {
  real mu_z;
  real<lower=0> tau_z;
}
model {
  mu_z ~ normal(0, 1);
  tau_z ~ cauchy(0, 1);
  z ~ multi_normal(rep_vector(mu_z, N), Sigma_z + diag_matrix(rep_vector(square(tau_z), N)));
}
generated quantities {
  real mu = mu_z * s + m;
  real tau = tau_z * s;
}
