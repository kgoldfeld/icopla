library(simstudy)
library(cmdstanr)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)
library(glue)
library(collapse)
library(posterior)


load("/gpfs/home/dw2625/r/SIM/draws_dt_np.rda")
load("/gpfs/home/dw2625/r/SIM/dind_np_local.rda")
#######draw from posterior iteration#################

newdt <- dind[, ordY := NULL]
newdt <- dind[, z := NULL]

newdt[, calt := control * study]
setkey(newdt, "id")

load("/gpfs/home/dw2625/r/SIM/dind_np_local.rda")##load dind again
#load:dind;basestudy;fit_co

est_from_draw <- function(n_draw, draws_dt, dind, newdt) {
  
  draw <- as.data.frame(draws_dt[sample(.N, 1)])
  
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ds) + factor(ss), data = dind)[, -1]
  D <- ncol(x)
  K <- dind[, length(unique(study))] 
  
  
  alpha <- draw$alpha
  beta <- as.vector(x = draw[, glue("beta[{1:D}]")], mode = "numeric")
  delta_k <- as.vector(draw[, glue("delta_k[{1:K}]")], mode = "numeric")
  gamma_k <- as.vector(draw[grep("^gamma_k", colnames(draw))], mode = "numeric")
  gamma_k <- matrix(gamma_k, nrow = 9)
  
  
  tau <- as.vector(draw[grep("^tau", colnames(draw))], mode = "numeric")
  tau <- matrix(tau, nrow = K)
  tau <- cbind(tau, Inf)
  cprop <- t(apply(tau, 1, stats::plogis))
  xprop <- t(apply(cprop, 1, diff))
  baseline <- cbind(cprop[,1], xprop) 
  
  
  #main effect: alpha + beta*x + control*delta_k
  main_coefs <- c(alpha, beta, delta_k)
  x <- model.matrix(~factor(who_enroll) + factor(age) + 
                      factor(sex) + factor(ds) +factor(ss), data = newdt)
  calt <- model.matrix(~factor(calt) - 1, data = newdt)[,-1]
  zmat <- cbind(x, calt)
  newdt$z_main <- zmat %*% main_coefs
  
  #interaction effect: control*gamma_k[s]
  x_ss <- model.matrix(~factor(ss) - 1, data = newdt)
  newdt$z_inter<- diag(calt %*%gamma_k %*% t(x_ss))
  newdt$z <- newdt$z_main + newdt$z_inter
  newdt[, calt := NULL]
  
  dl <- lapply(1:9, function(i) {
    b <- baseline[i,]
    dx <- newdt[study == i]
    dx <- genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
    dx[, ordY := factor(ordY, levels = c(1:11))]
    dx[]
  })
  
  dt_new <- rbindlist(dl)  
  
  ###################Plot#################
  dsum <- dt_new[, .(N = sum(.N)), keyby = .(ordY)]
  
  return(data.table(n_draw,dsum))
  
}


NJOBS <- 90


job <- Slurm_lapply(1:2100,
                    est_from_draw, draws_dt=draws_dt, dind=dind, newdt=newdt,
                    njobs = 90, 
                    mc.cores = 4,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "sim_84",
                    sbatch_opt = list(time = "44:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data

site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
result <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/SIM/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/SIM/", date_stamp, "/np_predictive_iter2000.rda"))


