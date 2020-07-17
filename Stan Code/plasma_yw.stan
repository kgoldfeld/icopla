//---------------------------------------------//
//  rstan model for Plasma interim evaluation  //
//                                             //
//   author:                KSG                //
//   last modified date:    05/29/2020         //
//                                             //
//---------------------------------------------//

data {
  int<lower=0> N;              // number of observations
  int<lower=2> K;              // number of WHO categories
  int<lower=1> J;              // number of sites
//  int<lower=1> D;              // number of covariates
  int<lower=1,upper=K> y[N];   // vector of categorical outcomes
  int<lower=1,upper=J> jj[N];  // site for individual
  int<lower=0,upper=1> rx[N];   // treatment arm for individual
  int<lower=0,upper=1> ss[N];   // strata for individual
//  row_vector[D] x[N];          // matrix of covariates
}

parameters {
//  vector[D] beta;              // covariate estimates
  ordered[K-1] tau;            // cut-points for cumulative odds model
  vector[J] alpha;             // site effects
  real gamma;                  // strata effect
  real delta;                  // treatment effect
}

transformed parameters{ 
  
  vector[N] yhat;
  
  for (i in 1:N)  
      yhat[i] = alpha[jj[i]] + ss[i] * gamma + rx[i] * delta; //+ x[i] * beta;
}

model {
  
  // priors
  
  alpha ~ normal(0, 10);
  gamma ~ normal(0, 10);
  delta ~ normal(0, 0.3537);
  
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau);
}

generated quantities {
  
  real<lower = 0> OR;

  OR = exp(delta); 
  
}

