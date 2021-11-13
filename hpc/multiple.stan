
data {
  int<lower=0> N;
  int<lower=0,upper=1> rx[N];
  vector[N] y;
  real p_mu;
  real p_sigma;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {

  real yhat[N];
  
  for (i in 1:N) 
    yhat[i] = alpha + beta * rx[i];
    
}

model {
  alpha ~ student_t(3, 0, 10);
  beta ~ student_t(3, p_mu, p_sigma);
  
  y ~ normal(yhat, sigma);
}

