
data {
  int<lower=0> N;                // number of observations
  int<lower=2> L;                // number of WHO categories
  int<lower=1> K;                // number of RCTs
  int<lower=1,upper=L> y[N];     // vector of categorical outcomes
  int<lower=1,upper=K> kk[N];    // RCT for individual
  int<lower=0,upper=1> ctrl[N];  // treatment or control
  int<lower=1,upper=3> cc[K];    // specific control for RCT
  int<lower=1> D;                // number of covariates
  row_vector[D] x[N];            // strata indicators  N x D matrix
  
  real<lower=0> prior_div;       // prior sd of tau
  real<lower=0> prior_Delta_sd;  // prior sd of overall control effect
  real<lower=0> eta;             // prior sd of delta
  real<lower=0> prior_beta_sd;   // prior sd of beta
  real<lower=0> prior_eta_0;     //prior sd of eta_0 (the sd of delta_k)
}

parameters {
  
 // real alpha;                    // overall intercept for treatment
  ordered[L-1] z_tau[K];           // cut-points for cumulative odds model (K X [L-1] matrix)
  real<lower=0>  z_eta_0;          // sd of delta_k (around delta)
  
  
  // non-central parameterization
  
  vector[K] z_ran_rx;            
  vector[3] z_delta;
  vector[D] z_beta;
  real z_Delta;
}

transformed parameters{ 
  
  vector[3] delta;               // control-specific effect
  vector[K] delta_k;             // RCT-specific treatment effect
  vector[D] beta;                // covariate estimates
  real Delta;                    // overall control effect 
  ordered[L-1] tau[K]; 
  real<lower=0>  eta_0; 
  vector[N] yhat;
  
  Delta = prior_Delta_sd * z_Delta;
  delta = eta * z_delta + Delta;
  beta = prior_beta_sd * z_beta;
  eta_0 = prior_eta_0 * z_eta_0;
  
  for (k in 1:K)
    delta_k[k] = eta_0 * z_ran_rx[k] + delta[cc[k]]; 
    
  for (k in 1:K)
    for (l in 1:(L-1)){
      tau[k, l] = prior_div * z_tau[k,l];
    }
  
  for (i in 1:N)  
    yhat[i] = ctrl[i] * delta_k[kk[i]] + x[i] * beta;
}

model {
  
  // priors
  //alpha ~ normal(0,0.1);
  z_ran_rx ~ std_normal(); 
  z_delta ~ std_normal();
  z_beta ~ std_normal();
  z_Delta ~ std_normal();
  z_eta_0 ~ std_normal();
  
  
  for (k in 1:K)
    for (l in 1:(L-1)){
      z_tau[k,l] ~ std_normal();
    }
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau[kk[i]]);
}

generated quantities {
  
  real OR;                 // overall CCP effect (odds ratio)
  real negDelta;           // overall CCP effect on log scale
  
  OR = exp(-Delta); 
  negDelta=-1*Delta;
  
}

