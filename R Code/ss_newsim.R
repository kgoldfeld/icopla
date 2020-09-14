library(data.table)
library(dplyr)
library(rstan)
library(simstudy)


### Function to generate base probabilities for each study

genBaseProbs <- function(n, base, similarity, digits = 2) {
  
  # n: number of studies
  n_levels <- length(base) ## WHO level
  
  #n series of random number from Dirichlet distribution, sums close to 1;
  #similarity * base: weight; similarity:similarity to base 
  
  x <- gtools::rdirichlet(n, similarity * base) 
  
  ### to ensure that each vector sums exactly to 1
  
  x <- round(floor(x*1e8)/1e8, digits) #
  xpart <- x[, 1:(n_levels-1)]      # delete the base prob of the final level
  partsum <- apply(xpart, 1, sum)
  x[, n_levels] <- 1 - partsum      # the base prob of the final level = 1- the rest
  
  return(x)
}


#### Data definitions

defC <- defDataAdd(varname = "a", formula = 0, variance = .005, 
                   dist = "normal")    
defC <- defDataAdd(defC, varname = "b", formula = 0, variance= .01, 
                   dist = "normal")                  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 

defC2 <- defDataAdd(
  varname="C_rv",
  formula="C * control",
  dist = "nonrandom"
) # C_rv=1/2/3: patient received control treatment C=1/2/3

defC2 <- defDataAdd(defC2, varname = "ss", 
                    formula = "0.25;0.25;0.25;0.25", dist = "categorical")  #symptom duration strata
defC2 <- defDataAdd(
  defC2,
  varname = "z", 
  formula = "0.1 * ss + (0.4 + b + (0.05 + a) * ss) * (C_rv==1) + (0.5 + b + (0.1 + a) * ss) * (C_rv==2) + (0.6 + b + (0.15 + a) * ss) * (C_rv==3)", 
  dist = "nonrandom")

#### Generate data

set.seed(184916)
nsites <- 21

basestudy <- genBaseProbs(n = nsites,
                          base =  c(.20, .30, .35, .15),
                          similarity = 100)

dstudy <- genData(nsites, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
# add control treatment variable, randomly assign control treatment in each study;
# rx=1:received control;rx=0:received convalescent plasma
dind <- trtAssign(dind, strata="study", grpName = "control") 
dind <- addColumns(defC2, dind)

dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})

dind <- rbindlist(dl)


## Define data for stan 

N = nrow(dind)                             ## number of observations 
L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
K <- dind[, length(unique(study))]            ## number of studies 
y <- as.numeric(dind$ordY)                    ## individual outcome 
kk <- dind$study                              ## study for individual 
ctrl <- dind$control                          ## treatment arm for individual 
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
ss <- dind$ss
eta <- .1 #sd of delta
prior_gamma_sd <- .1 #sd of gamma_c
prior_eta_0 <- 0.25  
prior_phi_s <- 0.25 
prior_tau_sd <- 2.5 
prior_Delta_sd <- .354  
prior_Gamma_sd <- 0.354
prior_beta_sd <- 10

studydata <- list(
  N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc,ss=ss,
  prior_tau_sd = prior_tau_sd, prior_Delta_sd = prior_Delta_sd,
  prior_eta_0 = prior_eta_0, eta = eta, prior_beta_sd=prior_beta_sd,
  prior_gamma_sd=prior_gamma_sd,prior_phi_s=prior_phi_s,
  prior_Gamma_sd=prior_Gamma_sd)


rt <- stanc("./ss_bayesmodel.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

a <- Sys.time()
fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500,
                 cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
Sys.time() - a

## model estimation

sparams <- get_sampler_params(fit, inc_warmup=FALSE)
div <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))

p.eff <- mean(extract(fit, pars = "OR")[[1]] < 1)
p.clinic <- mean(extract(fit, pars = "OR")[[1]] < 0.8)
Delta.perc <- quantile(extract(fit, pars = 'Delta')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) #control treatment effect

print(fit, pars = c("gamma","Gamma","delta_k", "eta_0","delta", "Delta", "OR","beta" ),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

