library(simstudy)
library(data.table)
library(rstan)
library(parallel)

# ---

def_i <- defDataAdd(
    varname = "y", 
    formula = "-1 + rx * 1.2", 
    dist = "binary", 
    link = "logit")

di <- genData(150)
di <- trtAssign(di, grpName = "rx")
di <- addColumns(def_i, di)
  
rt_c <- stanc("~/r/binary_simpler.stan")
sm_c <- stan_model(stanc_ret = rt_c, verbose=FALSE)
  
N <- nrow(di) ;                                 ## number of observations
y <- as.numeric(di$y)                           ## individual outcome
rx <- di$rx                                     ## treatment arm for individual

print("prepare data as list")

sampdat <- list(N=N, y=y, rx=rx)

fit <- sampling(sm_c, data = sampdat, 
                 iter = 3000, warmup = 500, 
                 cores = 4L, chains = 4)

print("done sampling")
