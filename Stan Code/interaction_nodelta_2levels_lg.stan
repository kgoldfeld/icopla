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
    int<lower=0,upper=1> y_2[N];     // vector of binary outcomes
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=1,upper=3> cc[K];    // specific control for site
    int<lower=1> D;              // number of covariates
    row_vector[D] x[N];          // matrix of covariates  N x D matrix
    row_vector[2] ds[N];    // matrix of age N x 2 matrix
    
    real<lower=0> prior_gamma_sd;             // sd of gamma  (around Gamma)
    real<lower=0> prior_Gamma_sd;  // prior sd
    real<lower=0> prior_tau_sd;  // prior sd
    real<lower=0> prior_phi_s; 
  }

parameters {
  
  vector[2] Gamma;         // overall strata effect
  real alpha;               // overall intercept for treatment
  real tau[K];      // cut-points for cumulative odds model
 
  
  
  vector[2] z_phi[K];      // K X 2 matrix 
  vector[D] z_beta;
  vector[2] z_gamma[3];    // 3 X 2 matrix
}

transformed parameters{ 
  
  vector[N] yhat;
  vector[D] beta;           // covariate estimates 
  
  vector[2] gamma[3];       // control-specific diabetes effect (3 X 2 matrix)
  vector[2] gamma_k[K];     // site-specific diabetes effect (K X 2 matrix)
  
  beta = 2.5 * z_beta;
 
  
    
  for (c in 1:3) 
        gamma[c] = prior_gamma_sd * z_gamma[c] + Gamma; //gamma_c[s] ~N(Gamma[s],sd=0.25)
  
  
  for (k in 1:K)
       gamma_k[k] = prior_phi_s * z_phi[k] + gamma[cc[k]];
       
  for (i in 1:N){
    yhat[i] = tau[kk[i]] + x[i] * beta + ctrl[i] * (ds[i]*gamma_k[kk[i]]);
  }
}

model {
  
  // priors
  
 z_beta ~ std_normal();
 
 Gamma ~ student_t(3, 0, prior_Gamma_sd);
  
  for (k in 1:K)
      z_phi[k] ~ std_normal();
      
  for (c in 1:3)
      z_gamma[c] ~ std_normal();

  for (k in 1:K)
      tau[k] ~ student_t(3, 0, prior_tau_sd);
  // outcome model
  
  for (i in 1:N)
    y_2[i] ~ bernoulli_logit(yhat[i]);
}

generated quantities {
  
   vector[2] OR;
  
  OR = exp(-Gamma); 
  
}

