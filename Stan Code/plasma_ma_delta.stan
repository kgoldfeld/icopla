//---------------------------------------------//
//  rstan model for Plasma interim evaluation  //
//                                             //
//   author:                KSG/DW             //
//   last modified date:    07/20/2020         //
//                                             //
//---------------------------------------------//

data {
  int<lower=0> N;                // number of observations
  int<lower=2> L;                // number of WHO categories
  int<lower=1> K;                // number of studies
  int<lower=1,upper=L> y[N];     // vector of categorical outcomes
  int<lower=1,upper=K> kk[N];    // site for individual
  int<lower=0,upper=1> ctrl[N];  // treatment or control
  int<lower=1,upper=3> cc[K];    // specific control for site
  real<lower=0> prior_tau_sd;    // prior sd
  real<lower=0> prior_Delta_sd;  // prior sd

}

parameters {
  ordered[L-1] tau[K];      // cut-points for cumulative odds model
  
  vector[K] delta_k;        // site specific treatment effect
  real<lower=0>  eta_0;     // sd of delta_k (around delta)
  
  vector[3] delta;          // control-specific effect
  real<lower=0> eta;        // sd of delta (around Delta)
  
  real Delta;               // overall control effect
}

transformed parameters{ 
  
  vector[N] yhat;
  
  for (i in 1:N)  
      // yhat[i] = ctrl[i] * delta[cc[i]] + ctrl[i] * b[kk[i]];
      yhat[i] = ctrl[i] * delta_k[kk[i]];
}

model {
  
  // priors
  
  // eta_0 ~ exponential(0.25);
  // eta ~ exponential(0.25);
  
  eta_0 ~ cauchy(0, 2.5);
  eta ~ cauchy(0, 2.5);
  
  for (k in 1:K)
    delta_k[k] ~ normal(delta[cc[k]], eta_0);
  
  for (c in 1:3)
    delta[c] ~ normal(Delta, eta);
    
  Delta ~ normal(0, prior_Delta_sd);
  
  for (k in 1:K)
    for (l in 1:(L-1))
      tau[k, l] ~ normal(0, prior_tau_sd);
  
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau[kk[i]]);
}

generated quantities {
      
  real OR;
  OR = exp(-Delta); 
  
}

