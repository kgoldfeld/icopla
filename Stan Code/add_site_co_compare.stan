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
    int<lower=1,upper=L> y[N];
    int<lower=0,upper=1> ctrl[N];  // treatment or control
    int<lower=1,upper=3> cc[K];    // specific control for RCTs
    int<lower=1,upper=K> site_rct[S]; //specific RCT for site
    int<lower=1> D;              // number of covariates
    row_vector[D] x[N];          // matrix of covariates  N x D matrix 
    
    real<lower=0> prior_div;    // prior dispersion parameters for intercept
    real<lower=0> prior_Delta_sd;  // prior sd
    real<lower=0> prior_beta_sd;             // sd of covariates
    real<lower=0> eta;             // sd of delta_c  (around Delta)
    real<lower=0> prior_eta_0;    //sd of eta
  }

parameters {
  real alpha;
  real<lower=0>  eta_0;     // sd of delta_k (around delta_c)
  ordered[L-1] tau[K]; //RCT_specific intercept for CO model
  
  // non-central parameterization
  vector[D] z_beta;
  vector[3] z_delta;
  vector[K] z_rct;      // for RCT-specific effect 
  vector[S] z_site;            //for site-specific effect 
  vector[S] z_phi;             //site random effect
  real z_Delta;               // overall control effect
}

transformed parameters{ 
  
  vector[N] yhat;            
  vector[D] beta;              // covariate estimates (2 x D matrix)
  vector[3] delta;          // control-specific effect
  vector[K] delta_k;         // RCT-specific effect
  vector[S] delta_site;     // Site-specific effect
  vector[S] phi;            // Site-specific random effect
  real Delta;               // overall control effect
   

     beta = prior_beta_sd * z_beta;
     Delta = prior_Delta_sd * z_Delta;
     delta = eta * z_delta + Delta;
     phi = 0.25 * z_phi; //phi~Normal(0,0.25)
   
   for (k in 1:K){
     delta_k[k] = eta_0 * z_rct[k] + delta[cc[k]]; 
   }
   
   for (s in 1:S){
     delta_site[s] = 0.1 * z_site[s] + delta_k[site_rct[s]];
   }
   
  for (i in 1:N){
     yhat[i] = alpha + phi[site[i]] + x[i] * beta + ctrl[i] * delta_site[site[i]];
   }
}
  
model{
  
  // priors
  
  alpha ~ normal(0,0.1);
  z_Delta ~ std_normal();
  z_delta ~ std_normal();
  z_beta ~ std_normal();
  z_rct ~ std_normal(); 
  z_site ~ std_normal(); 
  z_phi ~ std_normal(); 
  eta_0 ~ student_t(3,0, prior_eta_0);
    
    for (k in 1:K)
     for (l in 1:(L-1))
      tau[k, l] ~ student_t(3, 0, prior_div);
      
  // outcome model
  
  for (i in 1:N){
    y[i] ~ ordered_logistic(yhat[i],tau[kk[i]]);
   
   }
}

generated quantities {
  
  real OR;
  real negDelta;
  
  OR = exp(-Delta); 
  negDelta=-1*Delta;
  
}
