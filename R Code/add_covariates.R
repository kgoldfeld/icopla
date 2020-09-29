library(simstudy)
library(rstan)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)


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

defC <- defDataAdd(varname = "b", formula = 0, variance= .01, 
                   dist = "normal")    #each study has a random treatment effect
defC <- defDataAdd(defC, varname = "size", formula = "50+50*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defC2 <- defDataAdd(defC2,varname="sex",formula="sex_indicator+1", 
                    dist = "nonrandom")# sex=2:male;sex=1:female
###Assume male, older people, higher who_baseline do worse
defA1 <- defDataAdd(varname = "z", 
                    formula = "0.1*sex + 0.15*age + 0.2*WHO_enroll + (0.3 + b ) * (C_rv==1) + (0.4 + b ) * (C_rv==2) + (0.5 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

set.seed(184916)

nsites <- 12

basestudy <- genBaseProbs(n = nsites,
                          base =  c(.20, .30, .35, .15),
                          similarity = 100)


## Generate the data
dstudy <- genData(nsites, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "control") #add control treatment variable, randomly assign control treatment in each study;rx=1:received control;rx=0:received convalescent plasma
dind <- trtAssign(dind,strat="study",grpName = "sex_indicator")# add sex;male:1;female:0
dind <- trtAssign(dind,nTrt = 3,strat="study",grpName = "WHO_enroll")#add who_enrollline: 3 categories (4;5;6)
dind <- trtAssign(dind,nTrt = 3,strat="study",ratio=c(1,1,2),grpName = "age")#add age: 3 categories (ratio:1:1:2)
dind <- addColumns(defC2,dind)
dind <- addColumns(defA1, dind) # add z

##check proportions
#dind[,.N,keyby=.(study,age)]
#dind[,.N,keyby=.(study,WHO_enroll)]
#dind[,.N,keyby=.(study,sex)]

dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})

dind <- rbindlist(dl)

### Data for stan model
N = nrow(dind)                             ## number of observations
L <- dind[, length(unique(ordY))]             ## number of levels of outcome
K <- dind[, length(unique(study))]            ## number of studies
y <- as.numeric(dind$ordY)                    ## individual outcome
kk <- dind$study                              ## study for individual
ctrl <- dind$control                          ## treatment arm for individual
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study
age <- dind$age
WHO_enroll <- dind$WHO_enroll
sex <- dind$sex

eta <- .1 
prior_eta_0 <- 0.25 
prior_tau_sd <- 2.5
prior_Delta_sd <- .354 


studydata <- list(
  N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc,
  age=age,WHO_enroll=WHO_enroll,sex=sex,
  prior_tau_sd = prior_tau_sd, prior_Delta_sd = prior_Delta_sd,
  prior_eta_0 = prior_eta_0, eta = eta)

### Stan model

rt <- stanc("./add_covariates.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

ptm <- proc.time()
fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                 cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
proc.time() - ptm

sparams <- get_sampler_params(fit, inc_warmup=FALSE)
div <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
p.eff <- mean(extract(fit, pars = "OR")[[1]] < 1)
p.clinic <- mean(extract(fit, pars = "OR")[[1]] < 0.8)
Delta.perc <- quantile(extract(fit, pars = 'Delta')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) #control treatment effect
beta_age.perc <- apply(extract(fit, pars = 'beta_age')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
colnames(beta_age.perc) <- c("age_1","age_2","age_3") # 3 control type
beta_sex.perc <- apply(extract(fit, pars = 'beta_sex')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
colnames(beta_sex.perc) <- c("female","male") # 3 control type
beta_enroll.perc <- apply(extract(fit, pars = 'beta_enroll')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
colnames(beta_enroll.perc) <- c("enroll=4","enroll=5","enroll=6") # 3 control type


tab1 <- data.table(p.eff, p.clinic,Delta.perc,beta_age.perc,beta_sex.perc,beta_enroll.perc,div) 
