library(simstudy)
library(cmdstanr)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)
library(lme4)
library(posterior)
library(bayesplot)
library(glue)


set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")


mco <- cmdstan_model("/gpfs/home/dw2625/r/SIM/non_infor_vs_0354.stan")
#priorco <- cmdstan_model("./co_prior.stan")
# mco <- cmdstan_model("./non_infor_vs_0354.stan")

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
  iter <- function(iternum, defC, defC2, basestudy, nsites,mco)  {
    set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
    ## Generate the data
    dstudy <- genData(nsites, id = "study")       # generate studies
    dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
    dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
    dstudy <- addColumns(defC, dstudy)
    
    dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
    dind <- trtAssign(dind, strata="study", grpName = "control") #add control treatment variable, randomly assign control treatment in each study;rx=1:received control;rx=0:received convalescent plasma
    dind <- addColumns(defC2,dind)
   
    setkey(dind, "id")
    
    dl <- lapply(1:nsites, function(i) {
      b <- basestudy[i,]
      dx <- dind[study == i]
      dx <- genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
      dx[, ordY := factor(ordY, levels = c(1:11))]
      dx[]
    })
    
    dind <- rbindlist(dl)
    
    ### Data for stan model
    N = nrow(dind)                                ## number of observations
    L <- dind[, length(unique(ordY))]             ## number of levels of outcome
    K <- dind[, length(unique(study))]            ## number of studies
    y <- as.numeric(dind$ordY)                    ## individual outcome
    kk <- dind$study                              ## study for individual
    ctrl <- dind$control                          ## treatment arm for individual
    cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study
    x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss), data = dind)[, -1]
    D <- ncol(x)
    
 #################non-informative prior################   
    prior_div <- 100   #SD of site-specific intercept
    prior_Delta_sd <- 100
    eta <- 100     #sd of delta_c
    prior_eta_0 <- 100
    prior_beta_sd <- 100
    
    studydata <- list(
      N=N, L= L,K=K, y=y, ctrl=ctrl,cc=cc,kk=kk,prior_div=prior_div,
      prior_Delta_sd=prior_Delta_sd,eta=eta,
      prior_eta_0 = prior_eta_0,x=x,D=D,prior_beta_sd=prior_beta_sd)
    
    fit_co <- mco$sample(
      step_size = 0.1,
      data = studydata,
      chains = 4L,
      parallel_chains = 4L,
      refresh = 500,
      iter_warmup = 500,
      iter_sampling = 2500,
      adapt_delta=0.8
    )
    
    diagnostics_df <- as_draws_df(fit_co$sampler_diagnostics())
    div_num_prior <- sum(diagnostics_df[, 'divergent__'])
    
    res_co <- data.table(fit_co$summary(variables = c("negDelta","delta","beta")))[,.(variable,median)]
  
    res_prior <- t(res_co[,"median"])
    colnames(res_prior) <- res_co$variable
    ########model with 0354#########
    
    prior_div <- 100   #SD of site-specific intercept
    prior_Delta_sd <- .354
    eta <- 100     #sd of delta_c
    prior_eta_0 <- 100
    prior_beta_sd <- 100
    
    studydata <- list(
      N=N, L= L,K=K, y=y, ctrl=ctrl,cc=cc,kk=kk,prior_div=prior_div,
      prior_Delta_sd=prior_Delta_sd,eta=eta,
      prior_eta_0 = prior_eta_0,x=x,D=D,prior_beta_sd=prior_beta_sd)
    
    fit_co <- mco$sample(
      step_size = 0.1,
      data = studydata,
      chains = 4L,
      parallel_chains = 4L,
      refresh = 500,
      iter_warmup = 500,
      iter_sampling = 2500,
      adapt_delta=0.8
    )
    
    diagnostics_df <- as_draws_df(fit_co$sampler_diagnostics())
    div_num_co <- sum(diagnostics_df[, 'divergent__'])
    
    res_co <- data.table(fit_co$summary(variables = c("negDelta","delta","beta")))[,.(variable,median)]
    res_final <- t(res_co[,"median"])
    colnames(res_final) <- res_co$variable
    
    data.table(iternum,res_prior,div_num_prior,res_final,div_num_co)
    
  }

#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= .01, 
                   dist = "normal")    #each study has a random treatment effect
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defC2 <- defDataAdd(defC2,varname="sex", formula = 0.5, 
                    dist = "binary")   # sex=1:male;sex=0:female
defC2 <- defDataAdd(defC2,varname="who_enroll", formula = "1/3;1/3;1/3", 
                    dist = "categorical")   # 1==who 4, 2==who 5, 3 = who 6
defC2 <- defDataAdd(defC2,varname="age", formula = "0.25;0.25;0.50", 
                    dist = "categorical")   # 3 age categories

defC2 <- defDataAdd(defC2, varname = "ss", 
                    formula = "0.2;0.2;0.2;0.2;0.2", dist = "categorical")  #symptom duration strata

###Assume male, older people, higher who_baseline do worse

# beta_duration = c(0, .05, .1, .15), beta_age = c(0, .075, .15) and beta_who = c(0, .06, .12) 
defC2 <- defDataAdd(defC2, varname = "z", 
                    formula = "0.05*(ss-1) + 0.1*sex + 0.075*(age-1) + 0.06*(who_enroll-1) + (0.3 + b ) * (C_rv==1) + (0.4 + b ) * (C_rv==2) + (0.5 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

nsites <- 9

basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060),
                          similarity = 100)


NJOBS <- 90
# set.seed(382332)
# seeds <- sample(1:1000000,NJOBS,replace=FALSE)

job <- Slurm_lapply(1:2520,
                    iter, 
                    defC=defC,
                    defC2 = defC2,
                    basestudy=basestudy,
                    nsites=nsites,
                    mco = mco,
                    njobs = 90, 
                    mc.cores = 4,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "sim_242",
                    sbatch_opt = list(time = "38:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data

site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/SIM/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/SIM/", date_stamp, "/noninfor_0354_coef.rda"))

###   beta order in res_co
# 1: negDelta -0.24661850
# 2: delta[1]  0.22587200
# 3: delta[2]  0.23195650
# 4: delta[3]  0.29403900
# 5:  beta[1]  0.21102850
# 6:  beta[2]  0.17723700
# 7:  beta[3]  0.09341880
# 8:  beta[4]  0.05975005
# 9:  beta[5]  0.19049700
# 10:  beta[6]  0.16150450
# 11:  beta[7]  0.17174400
# 12:  beta[8]  0.30428500
# 13:  beta[9]  0.26899050
# 14:    alpha -0.00658741

