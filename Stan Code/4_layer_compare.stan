//---------------------------------------------//
  //  rstan model for logistic regression (4 layer)  //
  //    (with alpha)                              //
  //   author:                KSG/DW             //
  //   last modified date:    10/05/2020         //
  //                                             //
  //---------------------------------------------//
  
  data {
    int<lower=0> N;                // number of observations
    int<lower=1> K;                // number of RCT
    int<lower=1> S;                //number of site
    int<lower=2> L;
    int<lower=1,upper=K> kk[N];    // RCT for individual
    int<lower=1,upper=S> site[N];   //site for individual
    int<lower=1,upper=L> y_1[N];
    int<lower=0,upper=1>  y_2[N];     
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=1,upper=3> cc[K];    // specific control for RCTs
    int<lower=1,upper=K> site_rct[S]; //specific RCT for site
    int<lower=1> D;              // number of covariates
    row_vector[D] x[N];          // matrix of covariates  N x D matrix 
    
    real<lower=0> prior_div_logistic;    // prior dispersion parameters for intercept 
    real<lower=0> prior_div_co;    // prior dispersion parameters for intercept
    real<lower=0> prior_Delta_sd;  // prior sd
    real<lower=0> prior_beta_sd;             // sd of covariates
    real<lower=0> eta;             // sd of delta_c  (around Delta)
    real<lower=0> prior_eta_0;    //sd of eta
  }

parameters {
  real alpha[2];
  real<lower=0>  eta_0[4];     // sd of delta_k (around delta_c)
  ordered[L-1] tau_1[K]; //RCT_specific intercept for CO model
  real tau_2[K];       //RCT_specific intercept for logistic model
  
  ordered[L-1] tau_3[K]; //RCT_specific intercept for primary CO model
  real tau_4[K];       //RCT_specific intercept for logistic model
  // non-central parameterization
  vector[D] z_beta[4];
  vector[3] z_delta[4];
  vector[K] z_rct[4];      // for RCT-specific effect 
  vector[S] z_site[2];            //for site-specific effect 
  vector[S] z_phi[2];             //site random effect
  real z_Delta[4];               // overall control effect
}

transformed parameters{ 
  
  vector[N] yhat[4];            
  vector[D] beta[4];              // covariate estimates (2 x D matrix)
  vector[3] delta[4];          // control-specific effect
  vector[K] delta_k[4];         // RCT-specific effect
  vector[S] delta_site[2];     // Site-specific effect
  vector[S] phi[2];            // Site-specific random effect
  real Delta[4];               // overall control effect
   
  //////// With site as an extra layer   
  for (c in 1:2){
     beta[c] = prior_beta_sd * z_beta[c];
     Delta[c] = prior_Delta_sd * z_Delta[c];
     delta[c] = eta * z_delta[c] + Delta[c];
     phi[c] = 0.25 * z_phi[c]; //phi~Normal(0,0.25)
  }
  for (c in 1:2){
   for (k in 1:K){
     delta_k[c,k] = eta_0[c] * z_rct[c,k] + delta[c,cc[k]]; 
   }
  }
  
  for (c in 1:2){
   for (s in 1:S){
     delta_site[c,s] = 0.1 * z_site[c,s] + delta_k[c,site_rct[s]];
   }
  }
  for (i in 1:N){
     yhat[1,i]=alpha[1] +x[i]*beta[1]+ctrl[i]*delta_site[1,site[i]];
     yhat[2,i] = phi[2,site[i]] + tau_2[kk[i]] + x[i] * beta[2] + ctrl[i] * delta_site[2,site[i]];
   }
  
  //////// Primary model(without site as an extra layer) 
   for (c in 3:4){
     beta[c] = prior_beta_sd * z_beta[c];
     Delta[c] = prior_Delta_sd * z_Delta[c];
     delta[c] = eta * z_delta[c] + Delta[c];
  }
  for (c in 3:4){
   for (k in 1:K){
     delta_k[c,k] = eta_0[c] * z_rct[c,k] + delta[c,cc[k]]; 
   }
  }
  
  for (i in 1:N){
     yhat[3,i]= alpha[2] +x[i]*beta[3]+ctrl[i]*delta_k[3,kk[i]];
     yhat[4,i] = tau_4[kk[i]] + x[i] * beta[4] + ctrl[i] * delta_k[4,kk[i]];
   }
}

model{
  
  // priors
  
  alpha ~ normal(0,0.1);
  for (c in 1:2){  
    z_Delta[c] ~ std_normal();
    z_delta[c] ~ std_normal();
    z_beta[c] ~ std_normal();
    z_rct[c] ~ std_normal(); 
    z_site[c] ~ std_normal(); 
    z_phi[c] ~ std_normal(); 
    eta_0[c] ~ student_t(3,0, prior_eta_0);
  } 
  for (c in 3:4){  
    z_Delta[c] ~ std_normal();
    z_delta[c] ~ std_normal();
    z_beta[c] ~ std_normal();
    z_rct[c] ~ std_normal(); 
    eta_0[c] ~ student_t(3,0, prior_eta_0);
  } 
    
   for (k in 1:K){
    tau_1[k] ~ student_t(3, 0, prior_div_co);
    tau_2[k] ~ student_t(3, 0, prior_div_logistic);
    tau_3[k] ~ student_t(3, 0, prior_div_co);
    tau_4[k] ~ student_t(3, 0, prior_div_logistic);
   }
  
  // outcome model
  
  for (i in 1:N){
    y_1[i] ~ ordered_logistic(yhat[1,i],tau_1[kk[i]]);
    y_2[i] ~ bernoulli_logit(yhat[2,i]);
    //////without site as an extra layer
    y_1[i] ~ ordered_logistic(yhat[3,i],tau_3[kk[i]]);
    y_2[i] ~ bernoulli_logit(yhat[4,i]);
   }
}

generated quantities {
  
  real OR[4];
  
  OR[1] = exp(-Delta[1]); 
  OR[2] = exp(-Delta[2]); 
  OR[3] = exp(-Delta[3]); 
  OR[4] = exp(-Delta[4]); 
}
