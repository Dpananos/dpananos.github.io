functions{
  real sbg_lpdf(real time, real alpha, real beta){
    
    return lbeta(alpha+1, beta+time-1) - lbeta(alpha, beta);
  }
  
  real sbg_lccdf(real time, real alpha, real beta){
    
    return lbeta(alpha, beta + time) - lbeta(alpha, beta);
  }
  
}
data{
  // Observations
  int n;
  // Predictors (2 in this case)
  int p;
  array[n] real time;
  // Design matrix
  matrix[n, p] X;
  array[n] int lost;
  // These are the hacky bits.  I just don't want to work hard to 
  // Compute these in stan, so I do them in R.
  array[n] int sum_lost;
  array[n] int n_total;
}
transformed data{
  real truncation_time = max(time);
}
parameters{
  vector[p] zeta;
  vector[p] gamma;
}
transformed parameters{
  vector<lower=0, upper=1>[n] mu = inv_logit(X * zeta);
  vector<lower=0>[n] kappa = exp(X * gamma);
  
  // Kind of like beta regression.
  vector[n] alpha = mu./kappa;
  vector[n] beta = (1-mu)./kappa;
}
model{
  zeta ~ student_t(3.5, 0, 1);
  gamma ~ student_t(3.5, 0, 1);
  
  for(i in 1:n){
    if (time[i] == truncation_time){
      target += lost[i] * sbg_lpdf(time[i]| alpha[i], beta[i]);
      target += (n_total[i] - sum_lost[i]) * sbg_lccdf(truncation_time| alpha[i], beta[i]);
    }
    else{
      target += lost[i] * sbg_lpdf(time[i]| alpha[i], beta[i]);
    }
  }
}
generated quantities{
   array[n] real expected_surviving;
  
  for(i in 1:n){
    expected_surviving[i] = exp(sbg_lccdf(time[i]| alpha[i], beta[i]));
  }
}
