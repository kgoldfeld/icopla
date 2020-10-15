library(simstudy)
library(rstan)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)



#### Generate data and fit model
iter <- function(iternum, defC, defS,defC3, basestudy, sm){
  
  ## Generate the data
  dstudy <- genData(1200)
  dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group

  dind <- trtAssign(dstudy, strata="C", grpName = "control") 
  dind <- addColumns(defC, dind)
  dind <- addCondition(defS, dind, newvar = "z_ss")
  dind <- addColumns(defC3, dind)
  dind <- genOrdCat(dind, adjVar = "z", basestudy, catVar = "ordY")
  
  
  ##data for stan
  N <- nrow(dind)                                ## number of observations 
  L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
  K <- dind[, length(unique(C))]            ## number of studies 
  y <- as.numeric(dind$ordY)                    ## individual outcome 
  ctrl <- dind$control                          ## treatment arm for individual 
  cc <- dind$C       ## specific control arm for study 
  x <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]
  D <- ncol(x)
  ds <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]
  
  eta <- 1 #sd of delta_c
  prior_gamma_sd <- 0.25   
  prior_Delta_sd <- 1
  prior_Gamma_sd <- 1 
  
  studydata <- list(
    N=N, L= L, K=K, y=y, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x, eta = eta,prior_gamma_sd = prior_gamma_sd, 
    prior_Gamma_sd = prior_Gamma_sd,prior_Delta_sd = prior_Delta_sd)
  
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

rt <- stanc("/gpfs/home/dw2625/r/3_add_more_control.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
 

defC <- defDataAdd(
  varname="C_rv",
  formula="C * control",
  dist = "nonrandom"
) # C_rv=1/2/3: patient received control treatment C=1/2/3

defC <- defDataAdd(defC, varname = "ss", 
                    formula = "0.25;0.25;0.25;0.25", dist = "categorical")  #symptom duration strata

defS <- defCondition(
  condition = "ss==1",  
  formula = "0 * (C_rv==1) + 0 * (C_rv==2) + 0 * (C_rv==3)", 
  dist = "nonrandom")

defS <- defCondition(defS,
                     condition = "ss==2",  
                     formula = "(0.09) * (C_rv==1) + (0.1) * (C_rv==2) + (0.11) * (C_rv==3)", 
                     dist = "nonrandom")
defS <- defCondition(defS,
                     condition = "ss==3",  
                     formula = "(0.14) * (C_rv==1) + (0.15) * (C_rv==2) + (0.16) * (C_rv==3)", 
                     dist = "nonrandom")
defS <- defCondition(defS,
                     condition = "ss==4",  
                     formula = "(0.21) * (C_rv==1) + (0.20) * (C_rv==2) + (0.22) * (C_rv==3)", 
                     dist = "nonrandom")


defC3 <- defDataAdd(
  varname = "z", 
  formula = "0.05*(ss-1) + z_ss + (0.4) * (C_rv==1) + (0.5) * (C_rv==2) + (0.6) * (C_rv==3)", 
  dist = "nonrandom")

basestudy <- c(.35, .25, .20, .15, .05)




NJOBS <- 75
set.seed(382332)
seeds <- sample(1:1000000,NJOBS,replace=FALSE)

job <- Slurm_lapply(1:450,
                    iter, 
                    defC = defC,
                    defS =defS,
                    defC3 =defC3,
                    basestudy=basestudy,
                     sm = sm,
                    njobs = 75, 
                    mc.cores = 4,
                    seeds=seeds,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_110",
                    sbatch_opt = list(time = "12:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data
site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/", date_stamp, "/add_more_control.rda"))


