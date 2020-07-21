data {
  int<lower=2> K;              // number of categories
  int<lower=0> N;              // number of observations
  int<lower=1> D;              // number of covariattes
  int<lower=1,upper=K> y[N];   // vector of categorical outcomes
  row_vector[D] x[N];          // matrix of covariates
}

parameters {
  vector[D] beta;              // covariate estimates
  ordered[K-1] c;              // cut-points for cumulative odds model
}


model {
  
  // outcome model
  
  for (n in 1:N)
    y[n] ~ ordered_logistic(x[n] * beta, c);
}

generated quantities {
  
  real<lower = 0> OR;
  
  OR = exp(-beta[1]);

}

