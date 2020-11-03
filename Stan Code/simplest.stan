//---------------------------------------------//
  //  rstan model for Plasma primary analysis  //
  //    (with alpha)                              //
  //   author:                KSG/DW             //
  //   last modified date:    10/05/2020         //
  //                                             //
  //---------------------------------------------//
  
  data {
    int<lower=0> N;                // number of observations
    int<lower=2> L;                // number of WHO categories
    int<lower=1,upper=L> y_1[N];     // vector of categorical outcomes
    int<lower=0,upper=1> y_2[N];
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    
    real<lower=0> prior_tau_sd;    // prior sd
    real<lower=0> prior_Delta_sd;  // prior sd
  }

parameters {
  
  real Delta[2];               // overall control effect
  
  real alpha[2];               // overall intercept for treatment
  ordered[L-1] tau;      // cut-points for cumulative odds model
}

transformed parameters{ 
  
  vector[N] yhat[2];            // (2 x N matrix)
  
 
  for (c in 1:2)
    for (i in 1:N)  
      yhat[c,i] = alpha[c] + ctrl[i] *Delta[c];
}

model {
  
  // priors
  alpha[1] ~ normal(0,0.1);
  
  for (c in 1:2){  
    Delta[c] ~ normal(0, prior_Delta_sd);
    }
    
  for (l in 1:(L-1))
    tau[l] ~ student_t(3, 0, prior_tau_sd);
  
  // outcome model
  
  for (i in 1:N){
    y_1[i] ~ ordered_logistic(yhat[1,i], tau);
    y_2[i] ~ bernoulli_logit(yhat[2,i]);
   }
}
generated quantities {
  
  real OR[2];
  
  OR[1] = exp(-Delta[1]); 
  OR[2] = exp(-Delta[2]);
}
