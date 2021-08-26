library(simstudy)
library(cmdstanr)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)
library(glue)
library(collapse)
library(posterior)
library(bayesplot) 

load("/gpfs/home/dw2625/r/SIM/Predictive/draws_dt_v2.rda")
load("/gpfs/home/dw2625/r/SIM/Predictive/dind_v2.rda")
dind$control <- ifelse(dind$bl_treatment==1,0,1)  ## treatment arm for individual (control=1;CP=0)
dind$study <- as.numeric(dind$study)
newdt <- dind[, ordY := NULL]

newdt[, calt := control * study]
newdt$id <- newdt$ID

load("./dind_v2.rda")

est_from_draw <- function(n_draw, draws_dt, dind, newdt) {
  
  draw <- as.data.frame(draws_dt[sample(.N, 1)])
  
  x <- model.matrix(~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = dind)[, -1]
  D <- ncol(x)
  K <- dind[, length(unique(study))] 
  
  
  alpha <- draw$alpha
  beta <- as.vector(x = draw[, glue("beta[{1:D}]")], mode = "numeric")
  delta_k <- as.vector(draw[, glue("delta_k[{1:K}]")], mode = "numeric")
  
  
  tau <- as.vector(draw[grep("^tau", colnames(draw))], mode = "numeric")
  tau <- matrix(tau, nrow = K)
  tau <- cbind(tau, Inf)
  cprop <- t(apply(tau, 1, stats::plogis))
  xprop <- t(apply(cprop, 1, diff))
  baseline <- cbind(cprop[,1], xprop) 
  
  
  #main effect: alpha + beta*x + control*delta_k
  main_coefs <- c(alpha, beta, delta_k)
  x <- model.matrix(~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = newdt)
  #x <- model.matrix(~factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = newdt)
  calt <- model.matrix(~factor(calt) - 1, data = newdt)[,-1]
  zmat <- cbind(x, calt)
  newdt$z <- zmat %*% main_coefs

  newdt[, calt := NULL]
  setkey(newdt, "id")
  dl <- lapply(1:8, function(i) {
    b <- baseline[i,]
    dx <- newdt[study == i]
    dx <- genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
    dx[, ordY := factor(ordY, levels = c(1:11))]
    dx[]
  })
  
  dt_new <- rbindlist(dl)  
  
  ###################Plot#################
  dsum <- dt_new[, .(N = sum(.N)), keyby = .(control,ordY)]
  
  return(data.table(n_draw,dsum))
  
}

NJOBS <- 90 


job <- Slurm_lapply(1:10000, 
                    est_from_draw, draws_dt=draws_dt, dind=dind, newdt=newdt, 
                    njobs = 90,  
                    mc.cores = 4, 
                    tmp_path = "/gpfs/scratch/dw2625", 
                    job_name = "sim_95", 
                    sbatch_opt = list(time = "34:00:00"), 
                    plan = "wait", 
                    overwrite=TRUE) 

### Save data 

site_plasma_all <- Slurm_collect(job) # data is a list 
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message 
result <- rbindlist(site_plasma) # converting list to data.table 

date_stamp <- gsub("-", "", Sys.Date()) 
dir.create(file.path("/gpfs/home/dw2625/r/SIM/Predictive/", date_stamp), showWarnings = FALSE) 
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/SIM/Predictive", date_stamp, "/predictive_compile.rda")) 

