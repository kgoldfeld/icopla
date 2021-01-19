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
    row_vector[1] ds[N];    // vector of sex N x 1 
  }

parameters {
  
  real Delta;               // overall control effect
  real Gamma;         // overall strata effect
  real alpha;               // overall intercept for treatment
  ordered[L-1] tau[K];      // cut-points for cumulative odds model (K X [L-1] matrix)
  real<lower=0>  eta_0;     // sd of delta_k (around delta)
  
  vector[K] z_ran_rx;      // site-specific effect
  vector[1] z_phi[K];      // K X 1 matrix 
  vector[3] z_delta;
  vector[D] z_beta;
  vector[1] z_gamma[3];    // 3 X 1 matrix
}

transformed parameters{ 
  
  vector[N] yhat;
  vector[D] beta;           // covariate estimates 
  vector[3] delta;          // control-specific effect
  vector[K] delta_k;        // site specific treatment effect
  vector[1] gamma[3];       // control-specific sex effect (3 X 1 matrix)
  vector[1] gamma_k[K];     // site-specific sex effect (K X 1 matrix)
  
  beta = 2.5 * z_beta;
  delta = 0.1 * z_delta + Delta;
  
  for (k in 1:K)
    delta_k[k] = eta_0 * z_ran_rx[k] + delta[cc[k]]; 
    
  for (c in 1:3) 
        gamma[c] = 0.25 * z_gamma[c] + Gamma; //gamma_c[s] ~N(Gamma[s],sd=0.25)
  
  
  for (k in 1:K)
       gamma_k[k] = 1 * z_phi[k] + gamma[cc[k]];
       
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
  eta_0 ~ student_t(3, 0, 0.25);
  
  Delta ~ student_t(3, 0, 2.5);
  Gamma ~ student_t(3, 0, 1.5);
  
  for (k in 1:K)
      z_phi[k] ~ std_normal();
      
  for (c in 1:3)
      z_gamma[c] ~ std_normal();

  for (k in 1:K)
      tau[k] ~ student_t(3, 0, 8);
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau[kk[i]]);
}

generated quantities {
  vector[1] gd_c[3];       // control-specific sex effect (3 X 1 matrix)
  real GD;
  
  GD = Delta + Gamma;    // control effect when sex=1
  
  
  for (c in 1:3){
     gd_c[c] =delta[c] +gamma[c];//
  }
  
}
