data{
  vector[3] y;
  vector[3] A;
  
  vector[3] y_se;
  vector[3] A_se;  
}

transformed data{
  
  matrix[3, 2] Y = append_col(y, A);
  matrix[3, 2] tau = append_col(y_se, A_se);
}

parameters{
  real Gamma;
  real b_0;
  real b_1;
  cholesky_factor_corr[2] Omega;
}

transformed parameters{
  vector[3] first_stage = Gamma * trt;
  vector[3] second_stage = b_0 + b_1 .* first_stage;
  matrix[3, 2] mu = append_col(second_stage, first_stage);
}

model{
  
}