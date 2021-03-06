#---------------------------------------------#
#  simulation code for meta-analysis(fix eta) #
#                                             #
#   author:                KSG/DW             #
#   last modified date:    07/20/2020         #
#                                             #
#---------------------------------------------#

library("simstudy")
library("rstan")
library("data.table")

# Stan model

# rt <- stanc("~\\plasma_ma_delta.stan")
# rt <- stanc("./Stan Code/plasma_ma_delta_induced_dir.stan");
rt <- stanc("./Stan Code/explore_divergence.stan");
# rt <- stanc("./Stan Code/explore_divergence - gamma.stan");

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

defC <- defDataAdd(varname = "b", formula = 0, variance= .16, 
                   dist = "normal")    #each study has a random treatment effect
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                    formula = "(0.4 + b ) * (C_rv==1) + (0.5 + b ) * (C_rv==2) + (0.6 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

#### Generate data

#set.seed(184916)

nsites <- 18

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
cn <- dind[, .N, keyby = .(C)]$N       ## N for specific control arm for study

eta <- 2.5 #.354
# eta_0 <- .354 # .1, 1, 2.5
prior_eta_0 <- 5
prior_tau_sd <- 5
# prior_Delta_sd <- 5 # 0.354

studydata <- list(
  N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc, cn=cn,
  prior_tau_sd = prior_tau_sd, prior_Delta_sd = prior_Delta_sd,
  eta_0 = eta_0, prior_eta_0 = prior_eta_0, eta = eta)

## model estimation

fit <-  sampling(sm, data=studydata, iter = 5000, warmup = 500, 
                      cores = 4L, chains = 4, control = list(adapt_delta = 0.8))

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
