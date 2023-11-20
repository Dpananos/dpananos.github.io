functions{
  real sbg_lpdf(real time, real alpha, real beta){
    
    return lbeta(alpha+1, beta+time-1) - lbeta(alpha, beta);
  }
  
  real sbg_lccdf(real time, real alpha, real beta){
    
    return lbeta(alpha, beta + time) - lbeta(alpha, beta);
  }
  
}
data{
  // Fitting the model
  int N;
  int n_total;
  array[N] int lost;
  array[N] real time;
  
  // Making Predictions
  int N_pred;
  array[N_pred] real pred_times;
}
transformed data{
  real truncation_time = max(time);
}
parameters{
  real<lower=0> alpha;
  real<lower=0> beta;
}
model{
  alpha ~ cauchy(0, 1);
  beta ~ cauchy(0, 1);
  
  for(i in 1:N){
    target += lost[i] * sbg_lpdf(time[i]| alpha, beta);
  }
  target += (n_total - sum(lost)) * sbg_lccdf(truncation_time| alpha, beta);
}
generated quantities{
  
  array[N_pred] real expected_surviving;
  
  for(i in 1:N_pred){
    expected_surviving[i] = exp(sbg_lccdf(pred_times[i]| alpha, beta));
  }
  
}
