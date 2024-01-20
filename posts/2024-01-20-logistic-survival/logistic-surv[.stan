data {
  int N;
  int p;
  int Npred;
  

  array[N] int s;
  array[N] int n;
  real h;
  
  matrix[N, p] X;
  matrix[Npred, p] Xpred;
  
  real student_t_dof;
  
}
parameters {
  
  vector[p] beta;
  
}
model {
  beta ~ student_t(student_t_dof, 0, 1);
  
  s ~ binomial_logit(n, X * beta);
}
generated quantities {
  vector[Npred] haz = inv_logit(Xpred*beta);
  vector[Npred] Haz = append_row(0.0, cumulative_sum( h/2 *(haz[1:(Npred-1)] + haz[2:Npred]) ) );
  vector[Npred] S = exp(-Haz);
}