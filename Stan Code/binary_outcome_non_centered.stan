//---------------------------------------------//
//  rstan model for binary outcoe              //
//    (fix eta)                                //
//   author:                KSG                //
//   last modified date:    08/11/2020         //
//                                             //
//---------------------------------------------//

data {
  int<lower=0> N;                 // number of observations
  int<lower=0> C;                 // number of control types
  int<lower=1> K;                 // number of studies
  int y[N];                       // vector of categorical outcomes
  int<lower=1,upper=K> kk[N];     // site for individual
  int<lower=0,upper=1> ctrl[N];   // treatment or control
  int<lower=1,upper=C> cc[K];     // specific control for site
  
  // real eta_0;

}

parameters {
  real alpha;               // overall intercept for treatment
   
  real<lower=0> sigma_b;
  real<lower=0>  eta_0;     // sd of delta_k (around delta)


  vector[C] delta_c;          // control-specific effect

  real Delta;               // overall control effect
  
   // non-centered paramerization
  
  vector[K] z_ran_rx;   // site level random effects (by period)
  vector[K] z_ran_int;  // individual level random effects 
}

transformed parameters{ 
  
  vector[N] yhat;
  vector[K] beta_0;
  vector[K] delta_k;        // site specific treatment effect

  beta_0 = sigma_b * z_ran_int + alpha;
  
  for (i in 1:K)
    delta_k[i] = eta_0 * z_ran_rx[i] + delta_c[cc[i]]; 
  
  for (i in 1:N)  
      yhat[i] = beta_0[kk[i]] + ctrl[i] * delta_k[kk[i]];

}

model {
  
  // priors
  
  alpha ~ student_t(3, 0, 2.5);

  z_ran_int ~ std_normal();  
  z_ran_rx ~ std_normal();  

  sigma_b ~ cauchy(0, 2.5);
  eta_0 ~ cauchy(0, 2.5);
  
  delta_c ~ normal(Delta, 5);
  Delta ~ normal(0, 5);
  
  // outcome model
  
  y ~ bernoulli_logit(yhat);
}


