//---------------------------------------------//
//  rstan model for Plasma interim evaluation  //
//                                             //
//   author:                KSG/DW             //
//   last modified date:    07/16/2020         //
//                                             //
//---------------------------------------------//

data {
  int<lower=0> N;              // number of observations
  int<lower=2> L;              // number of WHO categories
  int<lower=1> K;              // number of studies
  int<lower=1,upper=L> y[N];   // vector of categorical outcomes
  int<lower=1,upper=K> kk[N];  // site for individual
  int<lower=0,upper=1> ctrl[N];  // treatment or control
  int<lower=1,upper=3> cc[N];  // specific control arm
}

parameters {
  ordered[L-1] tau[K];    // cut-points for cumulative odds model
  vector[K] b;            // random treatment effect
  vector[3] delta;        // treatment effect
}

transformed parameters{ 
  
  vector[N] yhat;
  
  for (i in 1:N)  
      yhat[i] = ctrl[i] * delta[cc[i]] + ctrl[i] * b[kk[i]];
}

model {
  
  // priors
  
  delta ~ normal(0, 10);
  
  for (k in 1:K)
    for (l in 1:(L-1))
      tau[k, l] ~ uniform(-10, 10);
  
  b ~ normal(0, 10);
  
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau[kk[i]]);
}

generated quantities {
      
  real<lower = 0> OR[3];

  for (i in 1:3)
    OR[i] = exp(-delta[i]); 
  
}

