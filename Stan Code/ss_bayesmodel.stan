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
    int<lower=1,upper=L> y[N];     // vector of categorical outcomes
    int<lower=1,upper=K> kk[N];    // site for individual
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=1,upper=3> cc[K];    // specific control for site
    int<lower=1,upper=4> ss[N];           //duration strata
 
    real<lower=0> prior_tau_sd;    // prior sd
    real<lower=0> prior_Delta_sd;  // prior sd
    real<lower=0> prior_Gamma_sd;  // prior sd
    real<lower=0> eta;             // sd of delta  (around Delta)
    real<lower=0> prior_gamma_sd;             // sd of gamma  (around Gamma)
    real<lower=0> prior_eta_0;     
    real<lower=0> prior_phi_s;    
    real<lower=0> prior_beta_sd;             // sd of covariates
  }

parameters {
  
  real Delta;               // overall control effect
  vector[4] Gamma;         // overall strata effect
  real alpha;               // overall intercept for treatment
  ordered[L-1] tau[K];      // cut-points for cumulative odds model
  
  real<lower=0>  eta_0;     // sd of delta_k (around delta)
  vector[4]  phi;     // sd of gamma_k (around gamma_c)
  
  // non-central parameterization
  
  vector[K] z_ran_rx;      // site-specific effect 
  vector[4] z_phi[K]; 
  vector[3] z_delta;
  vector[4] z_beta;
  vector[4] z_gamma[3];
}

transformed parameters{ 
  
  vector[3] delta;          // control-specific effect
  vector[K] delta_k;        // site specific treatment effect
  vector[4] gamma[3];         // control-specific duration strata effect (3*4)
  vector[N] yhat;
  vector[4] beta;              // covariate estimates of ss
  vector[4] gamma_k[K];  // site-specific duration strata effect (k*4)
  
  delta = eta * z_delta + Delta;
  
  for (s in 1:4)
    beta[s] = prior_beta_sd * z_beta[s]; //beta~N(0,sd=10)
    
  for (s in 1:4)
    for (c in 1:3)
      gamma[c,s] = prior_gamma_sd * z_gamma[c,s] + Gamma[s]; //gamma_c[s] ~N(Gamma[s],sd=0.1)
  
  for (k in 1:K){
    delta_k[k] = eta_0 * z_ran_rx[k] + delta[cc[k]]; 
  }
  
  
  for (k in 1:K)
    for (s in 1:4)
      gamma_k[k,s]= phi[s] * z_phi[k,s] + gamma[cc[k],s]; //gamma_k[s]~N(gamma_c[s],phi[s])
      
  for (i in 1:N)  
    yhat[i] = alpha + ctrl[i] * (delta_k[kk[i]] + ss[i] * gamma_k[kk[i],ss[i]]) +  ss[i] * beta[ss[i]];
}

model {
  
  // priors
  
  z_ran_rx ~ std_normal(); 
  z_delta ~ std_normal();

  alpha ~ normal(0, 0.1);
  eta_0 ~ cauchy(0, prior_eta_0);
  
  Delta ~ normal(0, prior_Delta_sd);
  
  for (s in 1:4) {
    phi[s] ~ cauchy(0, prior_phi_s);
    Gamma[s] ~ normal(0, prior_Gamma_sd);
    z_beta[s] ~ std_normal();
  }
  
  for (c in 1:3)
    for (s in 1:4)
      z_gamma[c,s] ~ std_normal();
  
  for (k in 1:K)
    for (s in 1:4)
      z_phi[k,s] ~ std_normal();
      
  for (l in 1:(L-1))
    for (k in 1:K)
      tau[k, l] ~ student_t(3, 0, prior_tau_sd);
  
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau[kk[i]]);
}

generated quantities {
  
  real OR;
  
  OR = exp(-Delta); 
  
}
