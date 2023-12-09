functions {
  
  real generalized_gamma_lpdf(real x, real k, real mu, real sigma) {
    real w;
    real d;
    w = (log(x) - mu) / sigma;
    d = (k - .5) * log(k) - log(sigma) - lgamma(k) + (sqrt(k) * w - k * exp(1 / sqrt(k) * w)) - log(x);
    
    return d;
  }
  
  real generalized_gamma_cdf(real x, real k, real mu, real sigma) {
    real w;
    real d;
    w = (log(x) - mu) / sigma;
    d = gamma_p(k, k * exp(w / sqrt(k)));
    
    return d;
  }
  
  real generalized_gamma_lcdf(real x, real k, real mu, real sigma) {
    real w;
    real d;
    w = (log(x) - mu) / sigma;
    d = log(gamma_p(k, k * exp(1 / sqrt(k) * w)));
    return d;
  }
  
  real generalized_gamma_lccdf(real x, real k, real mu, real sigma) {
    return log_diff_exp(0, generalized_gamma_lcdf(x | k, mu, sigma));
  }
  
  real generalized_gamma_rng(real k, real mu, real sigma) {
    real Q = 1.0 / sqrt(k);
    real gamma = gamma_rng(Q^-2, 1);
    real w = log(Q^2 * gamma) / Q;
    return exp(mu + sigma * w);
  }
  
}
data {
  int n;
  int n_trt;
  array[n] int trt;
  array[n] real time;
  
  array[n] int censored;
  
  // Predict survival curve
  int nt;
  array[nt] real pred_time;
  
  real mu_df;
  real sigma_df;
  real k_df;
  
  real mu_loc;
  real sigma_loc;
  real k_loc;
  
  real mu_scale;
  real sigma_scale;
  real k_scale;

}
parameters {
  vector[n_trt] mu;
  vector<lower=0>[n_trt] sigma;
  vector<lower=0>[n_trt] k;
}
model {
  mu ~ student_t(mu_df, mu_loc, mu_scale);
  sigma ~ student_t(sigma_df, sigma_loc, sigma_scale);
  k ~ student_t(k_df, k_loc, k_scale);
  
  for (i in 1 : n) {
    if (censored[i] == 1) {
      target += generalized_gamma_lccdf(time[i] | k[trt[i]], mu[trt[i]], sigma[trt[i]]);
    } else {
      target += generalized_gamma_lpdf(time[i] | k[trt[i]], mu[trt[i]], sigma[trt[i]]);
    }
  }
}
generated quantities {
  array[n, n_trt] real time_ppc;
  array[nt, n_trt] real survival_curve;
  array[nt, n_trt] real hazard;
  
  for(i in 1 : nt) {
    for (j in 1 : n_trt) {
      
      survival_curve[i, j] = generalized_gamma_cdf(pred_time[i]| k[j], mu[j], sigma[j]);
      hazard[i, j] = exp(generalized_gamma_lpdf(pred_time[i]| k[j], mu[j], sigma[j])) / survival_curve[i, j];
      
    }
  }
  
  
  for(i in 1 : n) {
    for (j in 1 : n_trt) {
      
      time_ppc[i, j] = generalized_gamma_rng(k[j], mu[j], sigma[j]);
      
    }
  }
  
  
}

