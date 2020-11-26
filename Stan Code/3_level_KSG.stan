data {
    int<lower=1> N;
    int<lower=1> R;
    int<lower=1> S;
    int<lower=1,upper=S> site[N];
    int<lower=1,upper=R> rct[N];
    
    int<lower=0,upper=1> rx[N];
    real y[N];
    
    int<lower=1,upper=R> site_rct[S];
   }

parameters {
    real b0;
    real b1;
    
    real z_b0_r[R];
    real z_b1_r[R];
    
    real z_b0_s[S];
    real z_b1_s[S];

    real<lower=0> sigma_s0;
    real<lower=0> sigma_s1;
    real<lower=0> sigma_r0;
    real<lower=0> sigma_r1;
    real<lower=0> sigma;

}

transformed parameters {
    real b0_rct[R];
    real b1_rct[R];
    
    real b0_site[S];
    real b1_site[S];
    
    real yhat[N];

    
    for (r in 1:R) {
        b0_rct[r] = b0 + sigma_r0 * z_b0_r[r];
        b1_rct[r] = b1 + sigma_r1 * z_b1_r[r];
    }
    
    for (s in 1:S) {
        b0_site[s] = b0_rct[site_rct[s]] + sigma_s0 * z_b0_s[s];
        b1_site[s] = b1_rct[site_rct[s]] + sigma_s1 * z_b1_s[s];
    }
    
    for (n in 1:N) {
        yhat[n] = b0_site[site[n]] + b1_site[site[n]] * rx[n];
    }
   
}

model {
    

    z_b0_s ~ std_normal();
    z_b1_s ~ std_normal();
    z_b0_r ~ std_normal();
    z_b1_r ~ std_normal();
    
    b0 ~ normal(0, 10);
    b1 ~ normal(0, 10);
    
    sigma_s0 ~ exponential(1);
    sigma_s1 ~ exponential(1);
    sigma_r0 ~ exponential(1);
    sigma_r1 ~ exponential(1);
    
    sigma ~ exponential(1);
    
    y ~ normal(yhat, sigma);

 }
 
 generated quantities {
     real var_s0 = sigma_s0^2;
     real var_s1 = sigma_s1^2;
     
     real var_r0 = sigma_r0^2;
     real var_r1 = sigma_r1^2;
     
     real var_i = sigma ^ 2;
     
 }
     