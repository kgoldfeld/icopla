//---------------------------------------------//
 //  rstan model for Plasma interim evaluation  //
  //    (add covariates beta normal dist        //
  //   author:                KSG/DW             //
  //   last modified date:    09/25/2020         //
  //                                             //
  //---------------------------------------------//
  
  data {
    int<lower=0> N;                      // number of observations
    int<lower=2> L;                      // number of WHO categories
    int<lower=1,upper=L> y[N];           // vector of categorical outcomes
    matrix<lower=0,upper=1>[N, 3] ss;    //sex strata
  }

parameters {
  
  real alpha;               // overall intercept for treatment
  ordered[L-1] tau;         // cut-points for cumulative odds model
  
  // non-central parameterization
  
  vector[3] z_beta_ss;    
}

transformed parameters{ 
  
  vector[3] beta_ss;     //covariate estimates of sex
  vector[N] yhat;
  
  beta_ss = 2.5 *  z_beta_ss;
  
  yhat = ss * beta_ss; // try to iterate as well
  
  // This should work as well if you need to do this:
  
  // for (i in 1: N) {
  //   yhat[i] = ss[i] * beta_ss;
  // }
  
  
}

model {
  
  // priors
  
  z_beta_ss ~ std_normal();    
  
  alpha ~student_t(3, 0, 2.5);

    for (l in 1:(L-1))
      tau[l] ~ student_t(3, 0, 1);
  
  // outcome model
  
  for (i in 1:N)
    y[i] ~ ordered_logistic(yhat[i], tau);
}

 
