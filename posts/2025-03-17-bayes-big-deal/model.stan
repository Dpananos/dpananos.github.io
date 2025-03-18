data{
  int n;
  vector[n] log_rr;
  vector[n] stderr;
  real beta;
}
parameters{
  real mu;
  real<lower=0> tau;
  vector[n] z;
}
transformed parameters {
  vector[n] theta = mu + tau * z;
}
model {
  tau ~ inv_gamma(3, beta);
  z ~ std_normal();
  log_rr ~ normal(theta, stderr);
}