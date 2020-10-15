library(simstudy)
library(rstan)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)



#### Generate data and fit model
iter <- function(iternum, defC, defS,defC3, defC2, nsites, basestudy, sm){
  
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
  dind <- genOrdCat(dind, adjVar = "z", basestudy, catVar = "ordY")
  
  
  ##data for stan
  N <- nrow(dind)                                ## number of observations 
  L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
  K <- dind[, length(unique(study))]            ## number of studies 
  y <- as.numeric(dind$ordY)                    ## individual outcome 
  kk <- dind$study                              ## study for individual 
  ctrl <- dind$control                          ## treatment arm for individual 
  cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
  x <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]
  D <- ncol(x)
  ds <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]
  
  eta <- 1 #sd of delta_c
  prior_gamma_sd <- 0.25   
  prior_Delta_sd <- 1
  prior_Gamma_sd <- 1 
  prior_eta_0 <- 0.25  
  prior_phi_s <- 2.5 
  
  studydata <- list(
    N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x, eta = eta,prior_gamma_sd = prior_gamma_sd, 
    prior_Gamma_sd = prior_Gamma_sd,prior_Delta_sd = prior_Delta_sd,prior_eta_0 = prior_eta_0,
    prior_phi_s = prior_phi_s)
  
  ##Fit the model
  fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                   cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  
  sparams <- get_sampler_params(fit, inc_warmup=FALSE)
  div <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  Delta.perc <- quantile(extract(fit, pars = 'Delta')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) #control treatment effect
  delta.perc <- apply(extract(fit, pars = 'delta')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  colnames(delta.perc) <- c("delta_1","delta_2","delta_3") # 3 control type

  Gamma.perc <- apply(extract(fit, pars = 'Gamma')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) )
  colnames(Gamma.perc) <- c("Gamma_2","Gamma_3","Gamma_4") # 3 strata; 1 reference
  alpha.perc <- quantile(extract(fit, pars = 'alpha')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  beta.perc <- apply(extract(fit, pars = 'beta')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  colnames(beta.perc) <- c("ss_2","ss_3","ss_4")
  
  
  data.table(iternum,Delta.perc,delta.perc,Gamma.perc,alpha.perc,beta.perc,div) 
  
}
### Stan model

rt <- stanc("/gpfs/home/dw2625/r/4_add_more_study.stan");
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

defC2 <- defDataAdd(defC2, varname = "ss", 
                    formula = "0.25;0.25;0.25;0.25", dist = "categorical")  #symptom duration strata

defS <- defCondition(
  condition = "ss==1",  
  formula = "0 * (C_rv==1) + 0 * (C_rv==2) + 0 * (C_rv==3)", 
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


defC3 <- defDataAdd(
  varname = "z", 
  formula = "0.05*(ss-1) + z_ss + (0.4 + b) * (C_rv==1) + (0.5 + b) * (C_rv==2) + (0.6 + b) * (C_rv==3)", 
  dist = "nonrandom")

nsites <- 15
basestudy <- c(.35, .25, .20, .15, .05)




NJOBS <- 90
set.seed(382332)
seeds <- sample(1:1000000,NJOBS,replace=FALSE)

job <- Slurm_lapply(1:990,
                    iter, 
                    defC = defC,
                    defS =defS,
                    defC2=defC2,
                    defC3 =defC3,
                    basestudy=basestudy,
                    nsites=nsites,
                    sm = sm,
                    njobs = 75, 
                    mc.cores = 4,
                    seeds=seeds,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_111",
                    sbatch_opt = list(time = "24:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data
site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/", date_stamp, "/add_more_study.rda"))


