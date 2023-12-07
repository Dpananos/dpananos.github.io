functions {
  
real lawless_generalized_gamma_lpdf(real x, real k, real mu, real sigma) {
  real w;
  real d;
  w = (log(x)-mu)/sigma;
  d = (k-.5)*log(k) - log(sigma) - lgamma(k) +
    (sqrt(k)*w - k*exp(1/sqrt(k)*w)) - y;
  return d;
}

real lawless_generalized_gamma_cdf(real x, real k, real mu, real sigma) {
  real w;
  real d;
  w = (log(x) - mu)/sigma;
  d = gamma_p(k, k*exp(w/sqrt(k)));
  return d;
}

real lawless_generalized_gamma_lccdf(real x, real k, real mu, real sigma) {
 
 return log_diff_exp(0, lawless_generalized_gamma_lcdf(x| k, mu, sigma));
 
}

}
data {
  int n;
  array[n] int censored;
  array[n] real time;
}
parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
  real<lower=0> k;
}
model {
  mu ~ student_t(30, 0, 1);
  sigma ~ student_t(30, 1, 0.05);
  k ~ student_t(30, 1, 0.05);
  
  for(i in 1:n) {
    if(censored[i] == 1) {
      target += lawless_generalized_gamma_lccdf(time[i] | k, mu, sigma);
    } else {
      target += lawless_generalized_gamma_lpdf(time[i] | k, mu, sigma);
    }
  }
}

