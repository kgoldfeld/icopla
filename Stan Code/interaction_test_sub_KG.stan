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
    int<lower=0,upper=1> y_2[N];   // vector of binary outcomes
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=1,upper=3> cc[K];    // specific control for site
    int<lower=1> D;                // number of covariates
    row_vector[D] x[N];            // matrix of covariates  N x D matrix
    int <lower=1> ss[N];           // vector of levels
    int <lower=1> S;               // number of subgroups
  }

parameters {
  
  vector[S] Gamma;         // overall strata effect
  
  // vector[K] z_ran_rx;         // site-specific effect
  vector[S] z_phi[K];            // K X S matrix 
  vector[D] z_beta;
  vector[S] z_gamma[3];    // 3 X S matrix
}

transformed parameters{ 
  
  vector[N] yhat;
  vector[D] beta;              // covariate estimates 
  vector[S] gamma[3];          // control-specific who_enroll effect (3 X S matrix)
  vector[S] gamma_k[K];     // site-specific who_enroll effect (K X n_ds matrix)
  
  beta = 2.5 * z_beta;
  
  for (c in 1:3) 
    gamma[c] = 0.25 * z_gamma[c] + Gamma; //gamma_c[s] ~N(Gamma[s],sd=0.25)
  
  
  for (k in 1:K)
    gamma_k[k] = 1 * z_phi[k] + gamma[cc[k]];
  
  for (i in 1:N){
    yhat[i] = x[i] * beta + ctrl[i] * gamma_k[kk[i], ss[i]];
  }
}

model {
  
  // priors
  
  z_beta ~ std_normal();
  // z_ran_rx ~ std_normal(); 

  Gamma ~ student_t(3, 0, 1.5);
  
  for (k in 1:K)
    z_phi[k] ~ std_normal();
  
  for (c in 1:3)
    z_gamma[c] ~ std_normal();
  
  for (i in 1:N)
    y_2[i] ~ bernoulli_logit(yhat[i]);
}



