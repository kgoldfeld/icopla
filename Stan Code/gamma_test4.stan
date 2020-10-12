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
    int<lower=1,upper=L> y[N];     // vector of categorical outcomes
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=0,upper=1> ds[N];    // matrix of duration strata N x 1 matrix
    real<lower=1> D;              // number of covariates
    int<lower=0,upper=1> x[N];          // matrix of covariates  N x D matrix
  }

parameters {
  
  real Delta;               // overall control effect
  real Gamma;         // overall strata effect
  real alpha;               // overall intercept for treatment
  ordered[L-1] tau;      // cut-points for cumulative odds model (K X [L-1] matrix)
  real z_beta;
}

transformed parameters{ 
  
  vector[N] yhat;
  real beta;           // covariate estimates 
   
  beta = 2.5 * z_beta;
    
  for (i in 1:N){
    yhat[i] = alpha + x[i] * beta + ctrl[i] * (Delta + ds[i]*Gamma);
  }
}

model {
  
  // priors
  
  z_beta ~ std_normal();
  alpha ~ student_t(3, 0, 2.5);
  
  // }
  
   for (l in 1:(L-1))
      tau[l] ~ student_t(3, 0, 1);
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau);
}
