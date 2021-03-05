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
iter <- function(iternum, defC, defS,defC3, defC2, nsites, basestudy, sm,lgm){
  
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
  y <- as.numeric(dind$ordY)                    ## individual CO outcome 
  y_2 <- ifelse(y %in% seq(8,11),1,0)            ## individual logistic outcome 
  kk <- dind$study                              ## study for individual 
  ctrl <- dind$control                          ## treatment arm for individual 
  cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss), data = dind)[, -1]
  D <- ncol(x)
  ds <- model.matrix(~factor(ss) - 1, data = dind)
  
  
  prior_gamma_sd <- 0.25   
  prior_Gamma_sd <- 1.5 
  prior_tau_sd <- 8 
  prior_phi_s <- 1 #sd of gamma_k
  
  studydata <- list(
    N=N, L= L, K=K, y=y,y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x,prior_gamma_sd = prior_gamma_sd, 
    prior_Gamma_sd = prior_Gamma_sd,
    prior_phi_s = prior_phi_s,prior_tau_sd=prior_tau_sd)
  
  ##Fit the model
  fit_co <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                   cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  fit_l <-  sampling(lgm, data=studydata, iter = 3000, warmup = 500, 
                      cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  ##div
  sparams <- get_sampler_params(fit_l, inc_warmup=FALSE)
  div_l <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  sparams <- get_sampler_params(fit_co, inc_warmup=FALSE)
  div_co <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  
  ### P(OR<1 or 0.8) for diabetes 2 levels
  p.eff.l <- mean(extract(fit_l, pars = "OR")[[1]][,1] <1)
  p.clinic.l <- mean(extract(fit_l, pars = "OR")[[1]][,1] < 0.8)
  p.eff.l.d <- mean(extract(fit_l, pars = "OR")[[1]][,2] <1)
  p.clinic.l.d <- mean(extract(fit_l, pars = "OR")[[1]][,2] < 0.8)
  
  
  p.eff.co<- mean(extract(fit_co, pars = "OR")[[1]][,1] <1)
  p.clinic.co <- mean(extract(fit_co, pars = "OR")[[1]][,1] < 0.8)
  p.eff.co.d <- mean(extract(fit_co, pars = "OR")[[1]][,2] <1)
  p.clinic.co.d <- mean(extract(fit_co, pars = "OR")[[1]][,2] < 0.8)
  
  p <- data.table(p.eff.l,p.clinic.l,p.eff.l.d,p.clinic.l.d, p.eff.co,p.clinic.co,p.eff.co.d,p.clinic.co.d)
  
  ### Mean, median, sd
  sl <- summary(fit_l, pars = c("OR","Gamma","gamma","beta"))
  slt <- data.table(t(sl$summary)[c(1,4,3),]) ##output: posterior mean(1) and median(4)and sd(3)
  colnames(slt) <- paste0(colnames(slt),"l")
  
  sco  <- summary(fit_co, pars = c("OR","Gamma","gamma","beta"))
  scot <- data.table(t(sco$summary)[c(1,4,3),]) ##output: posterior mean(1) and median(4)and sd(3)
  colnames(scot) <- paste0(colnames(scot),"co")
  
  
  ### Output of each iteration
  data.table(cbind(iternum,slt,scot,p,div_l,div_co))
  
}

### Stan model

rt <- stanc("/gpfs/home/dw2625/r/interaction_nodelta/Code/interaction_nodelta_2levels_co.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

lg <- stanc("/gpfs/home/dw2625/r/interaction_nodelta/Code/interaction_nodelta_2levels_lg.stan");
lgm <- stan_model(stanc_ret = rt, verbose=FALSE)

#rt <- stanc("Stan code/interaction_nodelta_2levels_co.stan");
#sm <- stan_model(stanc_ret = rt, verbose=FALSE)
#lg <- stanc("Stan code/interaction_nodelta_2levels_lg.stan");
#lgm <- stan_model(stanc_ret = lg, verbose=FALSE)

### Data definition
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

defC2 <- defDataAdd(defC2, varname = "ss", formula = 0.25, 
                    dist = "binary")  #diabetes

defS <- defCondition(
  condition = "ss==0",  
  formula = "0 * (C_rv==1) + 0 * (C_rv==2) + 0 * (C_rv==3)", 
  dist = "nonrandom")


defS <- defCondition(defS,
                     condition = "ss==1",  
                     formula = "(0.19 + a) * (C_rv==1) + (0.20 + a) * (C_rv==2) + (0.21 + a) * (C_rv==3)", 
                     dist = "nonrandom")


defC3 <- defDataAdd(
  varname = "z", 
  formula = "0.1*ss + z_ss + 0.1*sex + 0.075*(age-1) + 0.06*(who_enroll-1) + (0.3 + b) * (C_rv==1) + (0.4 + b) * (C_rv==2) + (0.5 + b) * (C_rv==3)", 
  dist = "nonrandom")

nsites <- 9
basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060), 
                          similarity = 100) 

### Iteration

job <- Slurm_lapply(1:990,
                    iter, 
                    defC = defC,
                    defS =defS,
                    defC2=defC2,
                    defC3 =defC3,
                    basestudy=basestudy,
                    nsites=nsites,
                    sm = sm,
                    lgm=lgm,
                    njobs = 90, 
                    mc.cores = 4,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_280",
                    sbatch_opt = list(time = "18:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data
site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/interaction_nodelta/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/interaction_nodelta/", date_stamp, "/interaction_nodelta_v1.rda"))


