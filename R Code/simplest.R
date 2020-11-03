library(simstudy)
library(rstan)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)
library(MASS)

#### Generate data and fit model
iter <- function(iternum, defC, basestudy, sm){
  
  ## Generate the data
  dstudy <- genData(1000)
  dstudy <- trtAssign(dstudy, nTrt = 2, grpName = "C") # allocate to control group
  dind <- addColumns(defC, dstudy)
  dind <- genOrdCat(dind, adjVar = "z", basestudy, catVar = "ordY")
  
  
  ##data for stan
  N <- nrow(dind)                                ## number of observations 
  L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
  y_1 <- as.numeric(dind$ordY)                  ## WHO 11 points scale 
  y_2 <- ifelse(y_1 %in% seq(8,11),1,0)         ##WHO 11 points scale 7-10 
  ctrl <- dind$C                          ## treatment arm for individual 

  prior_tau_sd <- 5 
  prior_Delta_sd <- .354
  
  studydata <- list(
    N=N, L= L, y_1=y_1,y_2=y_2, ctrl=ctrl,prior_tau_sd=prior_tau_sd,prior_Delta_sd=prior_Delta_sd)
  
  ##Fit the model
  fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                   cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  
  ##Bayesian output
  sparams <- get_sampler_params(fit, inc_warmup=FALSE)
  div <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  
  p.eff.1 <- mean(extract(fit, pars = "OR[1]")[[1]] <1)
  p.clinic.1 <- mean(extract(fit, pars = "OR[1]")[[1]] < 0.8)
  p.eff.2<- mean(extract(fit, pars = "OR[2]")[[1]] <1)
  p.clinic.2 <- mean(extract(fit, pars = "OR[2]")[[1]] < 0.8)
  
  p <- data.table(p.eff.1,p.clinic.1, p.eff.2,p.clinic.2)
  
  s <- summary(fit, pars = c("alpha", "Delta"), probs = 0.5)
  st <- data.table(t(s$summary)[c(1,4),]) ##output: posterior mean(1) and median(4)
  
  ##Data for frequentist method
  studydata$y_1 <- factor(studydata$y_1)
  studydata$y_2 <- factor(studydata$y_2)
  studydata$ctrl <- factor(studydata$ctrl)
  
  ##Frequentist output
  fit_po <- clm(formula = y_1 ~ ctrl,
                data=studydata,link = "logit")#random int
  freq_po <- t(summary(fit_po)$coefficients[11,c(1,4)])
  colnames(freq_po) <- c("Estimat_PO","P_value_PO")
  
  fit_lg <- glm(y_2 ~ ctrl,data=studydata,family=binomial)
  freq_lg <- t(summary(fit_lg)$coefficients[2,c(1,4)])
  colnames(freq_lg) <- c("Estimat_LG","P_value_LG")
  
  data.table(iternum,cbind(p,st,div,freq_po,freq_lg))  
  
}
### Stan model

rt <- stanc("/gpfs/home/dw2625/r/simplest.stan");
#rt <- stanc("./simplest.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

defC <- defDataAdd(
  varname = "z", 
  formula = "(0.4) * (C==1)", 
  dist = "nonrandom")

#### Generate data

basestudy <-  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060)


job <- Slurm_lapply(1:450,
                    iter, 
                    defC =defC,
                    basestudy=basestudy,
                    sm = sm,
                    njobs = 90, 
                    mc.cores = 4,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_194",
                    sbatch_opt = list(time = "8:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data
site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/", date_stamp, "/simplest_end_point_[3].rda"))


