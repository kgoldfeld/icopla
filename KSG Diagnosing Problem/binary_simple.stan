
data {
  int<lower=0> N;                // number of observations
  int y[N];                      // vector of categorical outcomes
  int<lower=0,upper=1> rx[N];  // treatment or control
}

parameters {
  real alpha;               // overall intercept for treatment
  real Delta;               // overall control effect
}

transformed parameters{ 
  
  vector[N] yhat;

  for (i in 1:N)  
      yhat[i] = alpha + rx[i] * Delta;

}

model {
  
  alpha ~ normal(0, 10);
  Delta ~ normal(0, 10);
  
  y ~ bernoulli_logit(yhat);
}


