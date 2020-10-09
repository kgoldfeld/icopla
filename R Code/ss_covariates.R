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

#### Generate data and fit model
iter <- function(iternum, defC, defC2, defS,defC3, basestudy, nsites, sm){
  
  ## Generate the data
  dstudy <- genData(nsites, id = "study")       # generate studies
  dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
  dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
  dstudy <- addColumns(defC, dstudy)
  
  dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
  dind <- trtAssign(dind, strata="study", grpName = "control") 
  dind <- addColumns(defC2, dind)
  dind <- addCondition(defS, dind, newvar = "z_ss")
  dind <- addColumns(defC3, dind)
  
  
  dl <- lapply(1:nsites, function(i) {
    b <- basestudy[i,]
    dx <- dind[study == i]
    genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
  })
  
  dind <- rbindlist(dl)  
  
  ##data for stan
  N <- nrow(dind)                                ## number of observations 
  L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
  K <- dind[, length(unique(study))]            ## number of studies 
  y <- as.numeric(dind$ordY)                    ## individual outcome 
  kk <- dind$study                              ## study for individual 
  ctrl <- dind$control                          ## treatment arm for individual 
  cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
  ss <- dind$ss
  age <- dind$age
  who_enroll <- dind$who_enroll
  sex <- dind$sex 
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss), data = dind)[, -1]
  D <- ncol(x)
  
  eta <- 1 #sd of delta 
  prior_gamma_sd <- 0.25                           #sd of gamma_c
  prior_eta_0 <- 0.25  
  prior_phi_s <- 2.5 
  prior_tau_sd <- 2.5 
  prior_Delta_sd <- 1
  prior_Gamma_sd <- 1.25 
  prior_beta_sd <- 0.25
  
  studydata <- list(
    N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc,ss=ss,D=D, x=x,
    prior_tau_sd = prior_tau_sd, prior_Delta_sd = prior_Delta_sd,
    prior_eta_0 = prior_eta_0, eta = eta, prior_beta_sd=prior_beta_sd,
    prior_gamma_sd = prior_gamma_sd, prior_phi_s = prior_phi_s,
    prior_Gamma_sd = prior_Gamma_sd)
  
  ##Fit the model
  fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                   cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  
  sparams <- get_sampler_params(fit, inc_warmup=FALSE)
  div <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  Delta.perc <- quantile(extract(fit, pars = 'Delta')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) #control treatment effect
  delta.perc <- apply(extract(fit, pars = 'delta')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  colnames(delta.perc) <- c("delta_1","delta_2","delta_3") # 3 control type
  Gamma.perc <- apply(extract(fit, pars = 'Gamma')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) )
  colnames(Gamma.perc) <- c("Gamma_1","Gamma_2","Gamma_3","Gamma_4") # 4 strata
  
  alpha.perc <- quantile(extract(fit, pars = 'alpha')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  beta.perc <- apply(extract(fit, pars = 'beta')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  colnames(beta.perc) <- c("enroll_5","enroll_6","age_2","age_3","male","ss_2","ss_3","ss_4")
  
  data.table(iternum,Delta.perc,delta.perc,Gamma.perc,alpha.perc,beta.perc,div) 
  
}
### Stan model

rt <- stanc("/gpfs/home/dw2625/r/ss_covariates.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)


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

defC2 <- defDataAdd(defC2,varname="sex", formula = 0.5, 
                    dist = "binary")   # sex=1:male;sex=0:female
defC2 <- defDataAdd(defC2,varname="who_enroll", formula = "1/3;1/3;1/3", 
                    dist = "categorical")   # 1==who 4, 2==who 5, 3 = who 6
defC2 <- defDataAdd(defC2,varname="age", formula = "0.25;0.25;0.50", 
                    dist = "categorical")   # 3 age categories

defC2 <- defDataAdd(defC2, varname = "ss", 
                    formula = "0.25;0.25;0.25;0.25", dist = "categorical")  #symptom duration strata

defS <- defCondition(
  condition = "ss==1",  
  formula = "(0.09 + a) * (C_rv==1) + (0.10 + a) * (C_rv==2) + (0.11 + a) * (C_rv==3)", 
  dist = "nonrandom")
defS <- defCondition(defS,
                     condition = "ss==2",  
                     formula = "(0.14 + a) * (C_rv==1) + (0.15 + a) * (C_rv==2) + (0.16 + a) * (C_rv==3)", 
                     dist = "nonrandom")
defS <- defCondition(defS,
                     condition = "ss==3",  
                     formula = "(0.19 + a) * (C_rv==1) + (0.20 + a) * (C_rv==2) + (0.21 + a) * (C_rv==3)", 
                     dist = "nonrandom")
defS <- defCondition(defS,
                     condition = "ss==4",  
                     formula = "(0.24 + a) * (C_rv==1) + (0.25 + a) * (C_rv==2) + (0.26 + a) * (C_rv==3)", 
                     dist = "nonrandom")

# beta_duration = c(0, .05, .1, .15), beta_age = c(0, .075, .15) and beta_who = c(0, .06, .12) 
defC3 <- defDataAdd(
  varname = "z", 
  formula = "0.05*(ss-1) + 0.1*sex + 0.075*(age-1) + 0.06*(who_enroll-1) + z_ss + (0.4 + b) * (C_rv==1) + (0.5 + b) * (C_rv==2) + (0.6 + b) * (C_rv==3)", 
  dist = "nonrandom")

#### Generate data

set.seed(184916)
nsites <- 9

basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060), 
                          similarity = 100) #c(.20, .30, .35, .15),

NJOBS <- 90
set.seed(382332)
seeds <- sample(1:1000000,NJOBS,replace=FALSE)

job <- Slurm_lapply(1:1080,
                    iter, 
                    defC=defC,
                    defC2 = defC2,
                    defS =defS,
                    defC3 =defC3,
                    basestudy=basestudy,
                    nsites=nsites,
                    sm = sm,
                    njobs = 90, 
                    mc.cores = 4,
                    seeds=seeds,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_84",
                    sbatch_opt = list(time = "8:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data
site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/", date_stamp, "/ss_cov_start0.rda"))



