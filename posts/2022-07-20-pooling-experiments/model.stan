data{
  // For model fitting 
  int n; //number of observations in data
  int n_experiment; //number of experiments performed
  vector[n] estimated_log_relative_lift; //log relative lift estimated from each experiment
  vector[n] estimated_sd_relative_lift; //standard deviation of teh log relative lift estimated from experiemnt
  int experiment[n]; //experiment numericlal identifier
  
  //For generated quantities
  int n_months;
  int n_experiments_per_month;
  real n_group; //Sample size in each group
  real z_alpha;
  real baseline_rate;
  real page_views;
}
parameters{
  real mu; 
  real<lower=0> sigma;
  vector[n_experiment] z;
}
transformed parameters{
  vector[n] true_log_rr = mu + z[experiment] * sigma;
}
model{
  mu ~ student_t(3, 0, 2.5);
  sigma ~ cauchy(0, 1);
  z ~ std_normal();
  estimated_log_relative_lift ~ normal(true_log_rr, estimated_sd_relative_lift);
}
generated quantities{
  
  real possible_rl = exp(normal_rng(mu, sigma)); // Possiblle future lifts
  real lifts[n_months, n_experiments_per_month]; // Lifts we generate in each month
  real power[n_months, n_experiments_per_month]; // Power for the n experiments per month
  real es; // Effect size for the relative lift
  vector[n_months] lift_generated;
  vector[n_months] forecasted_lift;
  vector[n_months] incremental_activations;
  
  // For each month
  for(i in 1:n_months){
    // for each experiment in each month
    for(j in 1:n_experiments_per_month){
      // Draw a lift from the future lift distribution
      lifts[i, j] = exp(normal_rng(mu, sigma));
      // Compute the effect size for the drawn lift against the baseline rate
      es = 2*asin(sqrt(lifts[i, j] * baseline_rate)) - 2*asin(sqrt(baseline_rate)); 
      // Determine the power to detect that lift
      power[i, j] = 1 - normal_cdf( z_alpha - es * sqrt(n_group/2)| 0, 1); 
    }
    // Compute the weighted harmonic means of the lifts generated in the month
    // This is the lift generated in the month
    lift_generated[i] = exp(sum(to_vector(power[i, ]) .* log(to_vector(lifts[i, ]))));
    
  }
  // Compute the rolling lift.  This is a cumulative product written as a sum on the log scale.
  forecasted_lift = to_vector(exp(cumulative_sum(log(lift_generated))));

  // Compute the incremental activations based on the lfit above.
  // Again, this is the incremental activations generated up to and including month i
  incremental_activations = cumulative_sum(forecasted_lift .* rep_vector(baseline_rate * page_views, n_months)) - cumulative_sum(rep_vector(baseline_rate * page_views, n_months) );

  
  
}

