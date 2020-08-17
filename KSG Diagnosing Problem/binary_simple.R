library(simstudy)
library(data.table)
library(shinystan)
library(rstan)
library(parallel)

# ---

def_i <- defDataAdd(
    varname = "y", 
    formula = "-1 + rx * 1.2", 
    dist = "binary", 
    link = "logit")

set.seed(1711)

di <- genData(50)
di <- trtAssign(di, grpName = "rx")
di <- addColumns(def_i, di)
  
rt_c <- stanc("/gpfs/home/goldfk01/r/binary_simple.stan")
sm_c <- stan_model(stanc_ret = rt_c, verbose=FALSE)
  
N <- nrow(di) ;                                 ## number of observations
y <- as.numeric(di$y)                           ## individual outcome
rx <- di$rx                                     ## treatment arm for individual

sampdat <- list(N=N, y=y, rx=rx)

fit <- sampling(sm_c, data = sampdat, 
                 iter = 3000, warmup = 500, 
                 cores = 4L, chains = 4, seed = 4938)

s <- summary(fit, pars = c("alpha", "Delta"))
s$summary
