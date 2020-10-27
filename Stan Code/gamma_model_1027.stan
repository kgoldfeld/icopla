//---------------------------------------------//
  //  rstan model for Plasma interim evaluation  //
  //                                             //
  //   author:                KSG/DW             //
  //   last modified date:    09/10/2020         //
  //                                             //
  //---------------------------------------------//
  
  data {
    int<lower=0> N;                // number of observations
    int<lower=2> L;                // number of WHO categories
    int<lower=1> K;                // number of studies
    int<lower=1,upper=K> kk[N];    // site for individual
    int<lower=1,upper=L> y[N];     // vector of categorical outcomes
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=1,upper=3> cc[K];    // specific control for site
    int<lower=1> D;              // number of covariates
    row_vector[D] x[N];          // matrix of covariates  N x D matrix
    row_vector[3] ds[N];    // matrix of duration strata N x 3 matrix
    
    real<lower=0> eta;             // sd of delta  (around Delta)
    real<lower=0> prior_gamma_sd;             // sd of gamma  (around Gamma)
    real<lower=0> prior_Delta_sd;  // prior sd
    real<lower=0> prior_Gamma_sd;  // prior sd
    real<lower=0> prior_tau_sd;  // prior sd
    real<lower=0> prior_phi_s; 
    real<lower=0> prior_eta_0;  
  }

parameters {
  
  real Delta;               // overall control effect
  vector[3] Gamma;         // overall strata effect
  real alpha;               // overall intercept for treatment
  ordered[L-1] tau[K];      // cut-points for cumulative odds model (K X [L-1] matrix)
  real<lower=0>  eta_0;     // sd of delta_k (around delta)
  
  vector[K] z_ran_rx;      // site-specific effect
  vector[3] z_phi[K];      // K X 3 matrix 
  vector[3] z_delta;
  vector[D] z_beta;
  vector[3] z_gamma[3];    // 3 X 3 matrix
}

transformed parameters{ 
  
  vector[N] yhat;
  vector[D] beta;           // covariate estimates 
  vector[3] delta;          // control-specific effect
  vector[K] delta_k;        // site specific treatment effect
  vector[3] gamma[3];       // control-specific duration strata effect (3 X 3 matrix)
  vector[3] gamma_k[K];     // site-specific duration strata effect (K X 3 matrix)
  
  beta = 2.5 * z_beta;
  delta = eta * z_delta + Delta;
  
  for (k in 1:K)
    delta_k[k] = eta_0 * z_ran_rx[k] + delta[cc[k]]; 
    
  for (c in 1:3) 
        gamma[c] = prior_gamma_sd * z_gamma[c] + Gamma; //gamma_c[s] ~N(Gamma[s],sd=0.25)
  
  
  for (k in 1:K)
       gamma_k[k] = prior_phi_s * z_phi[k] + gamma[cc[k]];
       
  for (i in 1:N){
    yhat[i] = alpha + x[i] * beta + ctrl[i] * (delta_k[kk[i]] + ds[i]*gamma_k[kk[i]]);
  }
}

model {
  
  // priors
  
  z_beta ~ std_normal();
  z_delta ~ std_normal();
  z_ran_rx ~ std_normal(); 
  alpha ~ normal(0,0.1);
  eta_0 ~ student_t(3, 0, prior_eta_0);
  
  Delta ~ student_t(3, 0, prior_Delta_sd);
  Gamma ~ student_t(3, 0, prior_Gamma_sd);
  
  for (k in 1:K)
      z_phi[k] ~ std_normal();
      
  for (c in 1:3)
      z_gamma[c] ~ std_normal();

  for (k in 1:K)
      tau[k] ~ student_t(3, 0, prior_tau_sd);
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau[kk[i]]);
}
