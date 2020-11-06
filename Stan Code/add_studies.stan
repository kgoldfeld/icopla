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
    int<lower=1> K;                // number of studies
    int<lower=1,upper=K> kk[N];    // site for individual
    int<lower=1,upper=L> y_1[N];     // vector of categorical outcomes
    int<lower=0,upper=1> y_2[N];
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=1,upper=3> cc[K];    // specific control for site
    int<lower=1> D;              // number of covariates
    row_vector[D] x[N];          // matrix of covariates  N x D matrix 
    
    real<lower=0> prior_div;    // prior dispersion parameters for intercept 
    real<lower=0> prior_Delta_sd;  // prior sd
    real<lower=0> prior_beta_sd;             // sd of covariates
    real<lower=0> eta;             // sd of delta  (around Delta)
    real<lower=0> prior_eta_0;
  }

parameters {
  
  real Delta[2];               // overall control effect
  real<lower=0>  eta_0[2];     // sd of delta_k (around delta)
  
  real alpha;               // overall intercept for treatment
  ordered[L-1] tau_1[K];      // cut-points for cumulative odds model
  real tau_2[K];       //site_specific intercept
  
  // non-central parameterization
  vector[D] z_beta[2];
  vector[3] z_delta[2];
  vector[K] z_ran_rx[2];      // site-specific effect 
}

transformed parameters{ 
  
  vector[N] yhat[2];            // (2 x N matrix)
  vector[D] beta[2];              // covariate estimates (2 x D matrix)
  vector[3] delta[2];          // control-specific effect
  vector[K] delta_k[2];
  
   
  for (c in 1:2){
     beta[c] = prior_beta_sd * z_beta[c];
     delta[c] = eta * z_delta[c] + Delta[c];
   }
   
  for (c in 1:2)  
   for (k in 1:K){
     delta_k[c,k] = eta_0[c] * z_ran_rx[c,k] + delta[c,cc[k]]; 
   }
   
  for (i in 1:N){
     yhat[1,i] = alpha +  x[i] * beta[1] + ctrl[i] * delta_k[1,kk[i]];
     yhat[2,i] = tau_2[kk[i]] +  x[i] * beta[2] + ctrl[i] * delta_k[2,kk[i]];
   }
}

model {
  
  // priors
  alpha ~ normal(0,0.1);
  
  for (c in 1:2){  
    Delta[c] ~ normal(0, prior_Delta_sd);
    z_delta[c] ~ std_normal();
    z_beta[c] ~ std_normal();
    z_ran_rx[c] ~ std_normal(); 
    eta_0[c] ~ student_t(3,0, prior_eta_0);
    }
    
   for (k in 1:K){
    tau_1[k] ~ student_t(3, 0, prior_div);
    tau_2[k] ~ student_t(3, 0, prior_div);
   }
  
  // outcome model
  
  for (i in 1:N){
    y_1[i] ~ ordered_logistic(yhat[1,i], tau_1[kk[i]]);
    y_2[i] ~ bernoulli_logit(yhat[2,i]);
   }
}
generated quantities {
  
  real OR[2];
  
  OR[1] = exp(-Delta[1]); 
  OR[2] = exp(-Delta[2]);
}
