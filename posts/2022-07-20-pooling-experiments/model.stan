data{
  // For model fitting 
  int n; //number of observations in data
  int n_experiment; //number of experiments performed
  vector[n] estimated_log_relative_lift; //log relative lift estimated from each experiment
  vector[n] estimated_sd_relative_lift; //standard deviation of teh log relative lift estimated from experiemnt
  int experiment[n]; //experiment numericlal identifier
  
  //For generated quantities
  int n_experiments_per_year;
  real n_group; //Sample size in each group
}
parameters{
  real mu_metric; 
  real<lower=0> sig_ex;
  vector[n_experiment] z_ex;
}
transformed parameters{
  vector[n] true_log_rr = mu_metric + z_ex[experiment] * sig_ex;
}
model{
  mu_metric ~ student_t(3, 0, 2.5);
  sig_ex ~ cauchy(0, 1);
  z_ex ~ std_normal();
  estimated_log_relative_lift ~ normal(true_log_rr, estimated_sd_relative_lift);
}
generated quantities{
  
  real possible_rr = exp(normal_rng(mu_metric, sig_ex));
  real rr_over_year[n_experiments_per_year];
  real es;
  real power[n_experiments_per_year];
  real detected_lift[n_experiments_per_year];
  real forecasted_lift[n_experiments_per_year];
  
  
  for(i in 1:n_experiments_per_year){
    rr_over_year[i] =  exp(normal_rng(mu_metric, sig_ex));
    es = 2*asin(sqrt(rr_over_year[i] * 0.01)) - 2*asin(sqrt(0.01));
    power[i] = 1 - normal_cdf( 1.644854 - es * sqrt(n_group/2), 0, 1);    
    detected_lift[i] = power[i] * log(rr_over_year[i]);
  }
  forecasted_lift = exp(cumulative_sum(detected_lift));
}
