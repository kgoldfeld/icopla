#---------------------------------------------#
#  simulation code for meta-analysis(fix eta) #
#                                             #
#   author:                KSG/DW             #
#   last modified date:    08/12/2020         #
#                                             #
#---------------------------------------------#

library(simstudy)
library(rstan)
library(data.table)
library(bayesplot)

# Stan model

rt <- stanc("./Stan Code/explore_divergence_nc.stan");
rt <- stanc("./Stan Code/explore_divergence_c.stan");

sm <- stan_model(stanc_ret = rt, verbose=FALSE)

### Function to generate base probabilities for each study

genBaseProbs <- function(n, base, similarity, digits = 2) {
    
    # n: number of studies
    n_levels <- length(base) ## WHO level
    
    # n series of random number from Dirichlet distribution, sums close to 1;
    # similarity * base: weight; similarity:similarity to base 
    
    x <- gtools::rdirichlet(n, similarity * base) 
    
    ### to ensure that each vector sums exactly to 1
    
    x <- round(floor(x*1e8)/1e8, digits) #
    xpart <- x[, 1:(n_levels-1)]      # delete the base prob of the final level
    partsum <- apply(xpart, 1, sum)
    x[, n_levels] <- 1 - partsum      # the base prob of the final level = 1- the rest
    
    return(x)
  }

#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= .02, 
                   dist = "normal")    #each study has a random treatment effect
defC <- defDataAdd(defC, varname = "size", formula = "50+50*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                    formula = "(0.4 + b ) * (C_rv==1) + (0.5 + b ) * (C_rv==2) + (0.6 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

#### Generate data

#set.seed(184916)

nsites <- 9

dstudy <- genData(nsites, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "control") #add control treatment variable, randomly assign control treatment in each study;rx=1:received control;rx=0:received convalescent plasma
dind <- addColumns(defC2,dind)
dind <- addColumns(defA1, dind) # add z

### Study specific base probabilities

basestudy <- genBaseProbs(
  n = nsites,
  base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060),
  similarity = 100
)

#basestudy <- genBaseProbs(n = nsites,
#               base =  c(.15, .15, .20, .30, .20) ,
#               similarity = 200)

dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})

dind <- rbindlist(dl)

### Check proportions

dind[, table(study, ordY, control)]

getProp <- function(d) {
  d[, round(prop.table(table(C_rv, ordY), margin = 1), 2)]
}

getProp(dind)

lapply(1:nsites, function(x) getProp(dind[study == x]))

### prepare data for model fitting

N = nrow(dind) ;                              ## number of observations
L <- dind[, length(unique(ordY))]             ## number of levels of outcome
K <- dind[, length(unique(study))]            ## number of studies
y <- as.numeric(dind$ordY)                    ## individual outcome
kk <- dind$study                              ## study for individual
ctrl <- dind$control                          ## treatment arm for individual
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study

eta <- .2 #.354
prior_eta_0 <- 0.5
prior_tau_sd <- 2.5
prior_Delta_sd <- .354 # 0.354

studydata <- list(
  N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc,
  prior_tau_sd = prior_tau_sd, prior_Delta_sd = prior_Delta_sd,
  prior_eta_0 = prior_eta_0, eta = eta)

## model estimation

  a <- Sys.time()
  fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                        cores = 4L, chains = 4, control = list(adapt_delta = 0.9))
  Sys.time() - a

print(fit, pars = c("alpha","delta_k", "eta_0","delta", "Delta", "OR"),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

print(fit, pars = c("alpha","delta_k", "delta","eta_0"),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

print(fit, pars = c("Delta_num", "Delta", "OR"),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

print(fit, pars = c("gamma","delta_k", "eta_0","delta", "Delta", "OR" ),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

print(fit, pars = c("tau"),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

print(fit, pars = c("yhat"),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)


traceplot(fit, pars = c("alpha","delta", "Delta","eta_0"), ncol = 1)
traceplot(fit, pars = c("delta_k"))
traceplot(fit, pars = c("tau"))
traceplot(fit, pars = c("eta_0"))
traceplot(fit, pars = c("eta"))

#----


posterior_ncp <- as.array(fit)
lp_ncp <- log_posterior(fit)
np_ncp <- nuts_params(fit)

color_scheme_set("darkgray")
mcmc_parcoord(posterior_ncp, np = np_ncp, 
    pars = c("eta_0", "Delta", "delta[1]", "delta[2]", "delta[3]",
             "delta_k[1]", "delta_k[2]", "delta_k[3]", "delta_k[4]",
             "delta_k[5]", "delta_k[6]", "delta_k[7]", "delta_k[8]", "delta_k[9]", 
             "alpha"))

# "tau[1,1]", "tau[1,2]", "tau[1,3]", "tau[1,5]", "tau[1,8]", "tau[1,9]", "tau[1,10]"


mcmc_pairs(posterior_ncp, np = np_ncp, pars = c("Delta", "eta_0"), 
           off_diag_args = list(size = 0.75))


mcmc_scatter(
  posterior_ncp, 
  pars = c("Delta","delta_k[7]"), 
  # transform = list("delta_k[7]" = "log"), # can abbrev. 'transformations'
  np = np_ncp, 
  size = 1
)



color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_ncp, pars = "eta_0", np = np_ncp) + 
  xlab("Post-warmup iteration")


mcmc_trace(posterior_ncp, pars = "eta_0", np = np_ncp, window = c(1600, 1800)) + 
  xlab("Post-warmup iteration")

mcmc_trace(posterior_ncp, pars = "eta_0", np = np_ncp, window = c(2200, 2400)) + 
  xlab("Post-warmup iteration")

mcmc_trace(posterior_ncp, pars = "Delta", np = np_ncp, window = c(1600, 1800)) + 
  xlab("Post-warmup iteration")

mcmc_trace(posterior_ncp, pars = "Delta", np = np_ncp, window = c(2200, 2400)) + 
  xlab("Post-warmup iteration")

mcmc_trace(posterior_ncp, pars = "tau[3,6]", np = np_ncp, window = c(1600, 1800)) + 
  xlab("Post-warmup iteration")

mcmc_trace(posterior_ncp, pars = "tau[3,6]", np = np_ncp, window = c(2200, 2400)) + 
  xlab("Post-warmup iteration")

mcmc_trace(posterior_ncp, pars = "alpha", np = np_ncp, window = c(0, 2500)) + 
  xlab("Post-warmup iteration")

color_scheme_set("red")
mcmc_nuts_divergence(np_ncp, lp_ncp)
