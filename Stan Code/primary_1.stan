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
    int<lower=1,upper=L> y_1[N];     // vector of categorical outcomes
    int<lower=0,upper=1> y_2[N];
    int<lower=0,upper=1> y_3[N];
    int<lower=1,upper=K> kk[N];    // site for individual
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=1,upper=3> cc[K];    // specific control for site
    int<lower=1> D;              // number of covariates
    
    real<lower=0> prior_tau_sd;    // prior sd
    real<lower=0> prior_Delta_sd;  // prior sd
    real<lower=0> eta;             // sd of delta  (around Delta)
    real<lower=0> prior_eta_0;
    row_vector[D] x[N];          // matrix of covariates  N x D matrix
    real<lower=0> prior_beta_sd;             // sd of covariates
  }

parameters {
  
  real Delta[3];               // overall control effect
  
  real alpha[3];               // overall intercept for treatment
  ordered[L-1] tau[K];      // cut-points for cumulative odds model
  
  real<lower=0>  eta_0[3];     // sd of delta_k (around delta)
  
  
  // non-central parameterization
  
  vector[K] z_ran_rx[3];      // site-specific effect (3 x K matrix)
  vector[3] z_delta[3];
  vector[D] z_beta[3];
}

transformed parameters{ 
  
  vector[3] delta[3];          // control-specific effect (3 x 3 matrix)
  vector[K] delta_k[3];        // site specific treatment effect (3 x K matrix)
  vector[N] yhat[3];            // (3 x N matrix)
  vector[D] beta[3];              // covariate estimates (3 x D matrix)
  
  for (c in 1:3)
  delta[c] = eta * z_delta[c] + Delta[c];
  
  for (c in 1:3)
  beta[c] = prior_beta_sd * z_beta[c];
  
  for (c in 1:3)
    for (k in 1:K)
    delta_k[c,k] = eta_0[c] * z_ran_rx[c,k] + delta[c,cc[k]]; 
  
  for (c in 1:3)
    for (i in 1:N)  
      yhat[c,i] = alpha[c] + ctrl[i] * delta_k[c,kk[i]] + x[i] * beta[c];
}

model {
  
  // priors
  for (c in 1:3){  
    z_ran_rx[c] ~ std_normal(); 
    z_delta[c] ~ std_normal();
    z_beta[c] ~ std_normal();
  
    alpha[c] ~ normal(0,0.1);
    eta_0[c] ~ student_t(3,0, prior_eta_0);
  
    Delta[c] ~ normal(0, prior_Delta_sd);
    }
  for (k in 1:K)
    for (l in 1:(L-1))
      tau[k, l] ~ student_t(3, 0, prior_tau_sd);
  
  // outcome model
  
  for (i in 1:N)
    y_1[i] ~ ordered_logistic(yhat[1,i], tau[kk[i]]);

  y_2 ~ bernoulli_logit(yhat[2]);
  y_3 ~ bernoulli_logit(yhat[3]);

}
generated quantities {
  
  real OR[3];
  
  OR[1] = exp(-Delta[1]); 
  OR[2] = exp(-Delta[2]);
  OR[3] = exp(-Delta[3]);
  
}
