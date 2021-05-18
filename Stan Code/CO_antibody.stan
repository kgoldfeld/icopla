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
    int<lower=0,upper=1> ctrl[N];  // treatment or control;control:CP units=0
    int<lower=1,upper=2> t[K];    // low or not low indicator for study
    int <lower=1,upper=2> ss[N];           // vector of levels of CP of the study that the individual belongs
    int<lower=1> D;              // number of covariates
    row_vector[D] x[N];          // matrix of covariates  N x D matrix
   
  }

parameters {
  
  
  real alpha;               // overall intercept for treatment
  ordered[L-1] tau[K];      // cut-points for cumulative odds model (K X [L-1] matrix)
  real<lower=0>  eta;     
  
  vector[2] z_delta;
  vector[D] z_beta;
  vector[K] z_ran_rx[2];        // (2 x K)
}

transformed parameters{ 
  
  vector[N] yhat;
  vector[D] beta;           // covariate estimates 
  vector[2] delta;               // levels effect
  vector[K] delta_k[2];        // site specific effect (2 x K)
  
  beta = 2.5 * z_beta;
  delta = 0.354 * z_delta ;
  
    for (s in 1:2)
     for (k in 1:K)
    delta_k[s,k] = eta * z_ran_rx[s,k] + delta[t[k]]; 
  
  
  for (i in 1:N){
    yhat[i] =  alpha   + x[i] * beta + ctrl[i] * (delta_k[ss[i],kk[i]]);
  }
}

model {
  
  // priors
  
  z_beta ~ std_normal();
  z_delta ~ std_normal();
 
  eta ~ student_t(3, 0, 0.25);
  alpha ~ normal(0,0.1);
  
  
  for (i in 1:2)
    z_ran_rx[i] ~ std_normal(); 
  
  for (k in 1:K)
    tau[k] ~ student_t(3, 0, 8);
  // outcome model
  
 for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau[kk[i]]);
}

