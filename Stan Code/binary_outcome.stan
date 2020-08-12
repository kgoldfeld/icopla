//---------------------------------------------//
//  rstan model for Plasma interim evaluation  //
//    (fix eta)                                //
//   author:                KSG/DW             //
//   last modified date:    08/11/2020         //
//                                             //
//---------------------------------------------//

data {
  int<lower=0> N;                // number of observations
  int<lower=0> C;                 // number of control types
  int<lower=1> K;                // number of studies
  int y[N];                     // vector of categorical outcomes
  int<lower=1,upper=K> kk[N];    // site for individual
  int<lower=0,upper=1> ctrl[N];  // treatment or control
  int<lower=1,upper=C> cc[K];            // specific control for site
  
  // real eta_0;

}

parameters {
  real alpha;               // overall intercept for treatment
  vector[K] beta_0; 
  real<lower=0> sigma_b;

  vector[K] delta_k;        // site specific treatment effect
  real<lower=0>  eta_0;     // sd of delta_k (around delta)

  vector[C] delta_c;          // control-specific effect
  // real<lower=0> eta;      // sd of delta (around Delta)
  
  real Delta;               // overall control effect
}

transformed parameters{ 
  
  vector[N] yhat;

  for (i in 1:N)  
      yhat[i] = alpha +  beta_0[kk[i]] + ctrl[i] * ( delta_k[kk[i]] );

}

model {
  
  // priors
  
  alpha ~ normal(0, 10);

  beta_0 ~ normal(0, sigma_b);
  sigma_b ~ cauchy(0, 5);

  eta_0 ~ cauchy(0, 5);

  for (k in 1:K)
      delta_k[k] ~ normal(delta_c[cc[k]], eta_0);

  delta_c ~ normal(Delta, 5);
 
  Delta ~ normal(0, 5);
  
  
  // outcome model
  
  y ~ bernoulli_logit(yhat);
}


