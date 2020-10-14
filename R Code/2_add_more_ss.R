library(simstudy)
library(rstan)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)



#### Generate data and fit model
iter <- function(iternum, defC2, defS,defC3, basestudy, sm){
  
  ## Generate the data
  dstudy <- genData(1000, defC2)
  dstudy <- trtAssign(dstudy, nTrt = 2, grpName = "C") # allocate to control group
  dstudy <- addCondition(defS, dstudy, newvar = "z_ss")
  dind <- addColumns(defC3, dstudy)
  dind <- genOrdCat(dind, adjVar = "z", basestudy, catVar = "ordY")
  
  
  ##data for stan
  N <- nrow(dind)                                ## number of observations 
  L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
  y <- as.numeric(dind$ordY)                    ## individual outcome 
  ctrl <- dind$C                          ## treatment arm for individual 
  x <- model.matrix(ordY ~ factor(ss), data = dind)[, -1] # covariants --dummy variable matrix
  D <-  ncol(x)
  ds <- model.matrix(ordY ~ factor(ss), data = dind)[, -1] 
  
  
  studydata <- list(
    N=N, L= L, y=y, ctrl=ctrl, ds=ds,D=D, x=x)
  
  ##Fit the model
  fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                   cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  
  sparams <- get_sampler_params(fit, inc_warmup=FALSE)
  div <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  Delta.perc <- quantile(extract(fit, pars = 'Delta')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) #control treatment effect
  Gamma.perc <- apply(extract(fit, pars = 'Gamma')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) )
  colnames(Gamma.perc) <- c("Gamma_2","Gamma_3","Gamma_4") # 3 strata; 1 reference
  
  alpha.perc <- quantile(extract(fit, pars = 'alpha')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  beta.perc <- apply(extract(fit, pars = 'beta')[[1]],2,function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  colnames(beta.perc) <- c("ss_2","ss_3","ss_4")
  
  
  data.table(iternum,Delta.perc,Gamma.perc,alpha.perc,beta.perc,div) 
  
}
### Stan model

rt <- stanc("/gpfs/home/dw2625/r/add_more_ss.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

rt <- stanc("./add_more_ss.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

defC2 <-  defData(varname = "ss", formula = "0.25;0.25;0.25;0.25", dist = "categorical")

defS <- defCondition(
  condition = "ss==1",  
  formula = "0 * (C==1)", 
  dist = "nonrandom")
defS <- defCondition(defS,
                     condition = "ss==2",  
                     formula = "(0.10) * (C==1)", 
                     dist = "nonrandom")

defS <- defCondition(defS,
                     condition = "ss==3",  
                     formula = "(0.15) * (C==1)", 
                     dist = "nonrandom")

defS <- defCondition(defS,
                     condition = "ss==4",  
                     formula = "(0.20) * (C==1)", 
                     dist = "nonrandom")

defC3 <- defDataAdd(
  varname = "z", 
  formula = "0.05*(ss-1) + z_ss + (0.4) * (C==1)", 
  dist = "nonrandom")

#### Generate data

basestudy <- c(.35, .25, .20, .15, .05)



NJOBS <- 75
set.seed(382332)
seeds <- sample(1:1000000,NJOBS,replace=FALSE)

job <- Slurm_lapply(1:450,
                    iter, 
                    defC2 = defC2,
                    defS =defS,
                    defC3 =defC3,
                    basestudy=basestudy,
                     sm = sm,
                    njobs = 75, 
                    mc.cores = 4,
                    seeds=seeds,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_109",
                    sbatch_opt = list(time = "12:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data
site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/", date_stamp, "/add_more_ss.rda"))


