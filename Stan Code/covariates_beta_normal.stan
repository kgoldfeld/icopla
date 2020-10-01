//---------------------------------------------//
 //  rstan model for Plasma interim evaluation  //
  //    (add covariates beta normal dist        //
  //   author:                KSG/DW             //
  //   last modified date:    09/25/2020         //
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
    int<lower=1,upper=3> WHO_enroll[N]; //WHO baseline strata
    int<lower=1,upper=3> age[N];           //age strata
    int<lower=1,upper=2> sex[N];           //sex strata
    int<lower=1,upper=4> ss[N];           //sex strata
    
    real<lower=0> prior_tau_sd;    // prior sd
    real<lower=0> prior_Delta_sd;  // prior sd
    real<lower=0> eta;             // sd of delta  (around Delta)
    real<lower=0> prior_eta_0;
  }

parameters {
  
  real Delta;               // overall control effect
  
  real alpha;               // overall intercept for treatment
  ordered[L-1] tau[K];      // cut-points for cumulative odds model
  
  real<lower=0>  eta_0;     
  
  // non-central parameterization
  
  vector[K] z_ran_rx;      // site-specific effect 
  vector[3] z_delta;
  vector[3] z_beta_age;     
  vector[3] z_beta_enroll;  
  vector[2] z_beta_sex;     
  vector[4] z_beta_ss;    
}

transformed parameters{ 
  
  vector[3] delta;          // control-specific effect
  vector[K] delta_k;        // site specific treatment effect
  vector[3] beta_age;     //covariate estimates of age
  vector[3] beta_enroll;  //covariate estimates of WHO-enroll
  vector[2] beta_sex;     //covariate estimates of sex
  vector[4] beta_ss;     //covariate estimates of sex
  vector[N] yhat;
  
  delta = eta * z_delta + Delta;
  beta_age = 2.5 * z_beta_age;
  beta_enroll = 2.5 * z_beta_enroll;
  beta_sex = 2.5 * z_beta_sex;
  beta_ss = 2.5 *  z_beta_ss;
  
  for (k in 1:K)
    delta_k[k] = eta_0 * z_ran_rx[k] + delta[cc[k]]; 
  
  for (i in 1:N)  
    yhat[i] = alpha + ctrl[i] * delta_k[kk[i]] + beta_ss[ss[i]] + beta_age[age[i]] + beta_enroll[WHO_enroll[i]] + beta_sex[sex[i]];
}

model {
  
  // priors
  
  z_ran_rx ~ std_normal(); 
  z_delta ~ std_normal();
  z_beta_age ~ std_normal();     
  z_beta_enroll ~ std_normal();  
  z_beta_sex ~ std_normal();     
  z_beta_ss ~ std_normal();    
  
  alpha ~student_t(3, 0, 2.5);
  eta_0 ~ cauchy(0, prior_eta_0);
  
  Delta ~ normal(0, prior_Delta_sd);
  
  for (k in 1:K)
    for (l in 1:(L-1))
      tau[k, l] ~ student_t(3, 0, prior_tau_sd);
  
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau[kk[i]]);
}

generated quantities {
  
  real OR;
  
  OR = exp(-Delta); 
  
}
