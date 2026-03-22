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
  vector[N] zz;
}
transformed parameters{
  vector[N] theta = mu_z + tau_z*zz;
}
model {
  mu_z ~ normal(0, 1);
  zz ~ std_normal();
  tau_z ~ cauchy(0, 1);
  z ~ multi_normal(theta, Sigma_z);
}
generated quantities {
  real mu = mu_z * s + m;
  real tau = tau_z * s;
}
