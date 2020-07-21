library(rstan)
library(slurmR)
library(simstudy)
library(data.table)

### Stan model

print("compiling stan")

rt <- stanc("~/r/ordmodel.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

### simulations

print("running simulations")

iter <- function(iter, maxn, sm) {
  
  genDT <- function(nobs, baseprobs, defA) {
    dT <- genData(nobs)
    dT <- trtAssign(dT, grpName = "exposed")
    dT <- addColumns(defA, dT)
    
    dT <- genOrdCat(dT, adjVar = "z", baseprobs, catVar = "r")
    dT[]
  }
  
  prepData <- function(n, dxfull) {
    
    dx <- dxfull[1:n,]
    
    K <- length(dx[, levels(r)])
    N <- dx[, .N]
    x <- as.matrix(dx[, .(exposed)])
    D <- ncol(x)
    y <- as.vector(dx[, as.numeric(r)])
    
    list(K = K, N = N, D= D, y = y, x = x)
  }
  
  
  doMCMC <- function(iter, n, dxfull) {
    
    sampdat <- prepData(n, dxfull)
    
    fit <-  sampling(sm, data=sampdat, seed = 3327, 
                     iter = 3000, warmup = 500, show_messages = FALSE, 
                     cores = 4, refresh = 0)
    
    p.eff <- mean(extract(fit, pars = "OR")[[1]] < 1)
    p.clinic <- mean(extract(fit, pars = "OR")[[1]] < 0.8)
    
    data.table(iter, n, p.eff, p.clinic)
  }
  
  cat("\f", iter, format(Sys.time(), "%X"))
  
  defA <- defDataAdd(varname = "z", formula = "0.30 * exposed", dist = "nonrandom")
  dxfull <- genDT(maxn, baseprobs = c(.15, .20, .35, .15, .10, .05), defA)
  rbindlist(lapply(seq(100, maxn, 100), function(x) doMCMC(iter, x, dxfull)))
  
}

job <- Slurm_lapply(1:100, iter, maxn = 500, sm = sm,
             njobs = 90,  # 50, 100
             mc.cores = 4,  # 10, 5
             tmp_path = "/gpfs/scratch/goldfk01",
             job_name = "bayes_interim",
             sbatch_opt = list(time = "00:30:00"),
             plan = "wait")

### Save data

print("saving data")

res <- Slurm_collect(job)

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("~/r/", date_stamp), showWarnings = FALSE)
save(res, file = paste0("~/r/", date_stamp, "/bayes_interim.rda"))



