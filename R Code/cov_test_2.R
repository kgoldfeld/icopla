
library(simstudy)
library(rstan)
library(slurmR)
library(data.table)

iteration <- function(iternum, baseprobs, d1, sm) {
  
  dd <- genData(1000, d1)
  dd <- genOrdCat(dd, adjVar = "z", baseprobs, "y")
  
  N = nrow(dd)                                  ## number of observations
  L <- dd[, length(unique(y))]                  ## number of levels of outcome
  y <- as.numeric(dd$y)                         ## individual outcome
  ss <- model.matrix(y~factor(ss), data = dd)[, -1]
  
  
  studydata <- list(N=N, L= L, y=y, ss=ss)
  
  fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                   cores = 4L, chains = 4)
  
  posterior <- as.array(fit)
  
  x <- summary(
    fit, 
    pars = c("beta_ss", "tau"),
    probs = c(0.025, 0.5, 0.975)
  )
  
  dpars <- data.table(iternum = iternum, par = rownames(x$summary), x$summary)
  
  dpars[]
}

directory <- "/gpfs/data/troxellab/ksg/r/"

baseprobs <- c(.35, .25, .20, .15, .05)

d1 <- defData(varname = "ss", formula = "0.25;0.25;0.25;0.25", dist = "categorical")
d1 <- defData(d1, varname = "z", formula = "0 + .1 * ss", dist = "nonrandom") 

rt_c <- stanc("/gpfs/data/troxellab/ksg/r/cov_test_2.stan")
sm_c <- stan_model(stanc_ret = rt_c, verbose=FALSE)

job <- Slurm_lapply(as.list(1:1050),
                    iteration, 
                    baseprobs = baseprobs,
                    d1 = d1,
                    sm = sm_c,
                    njobs = 75, 
                    mc.cores = 4,
                    tmp_path = "/gpfs/data/troxellab/ksg/scratch",
                    overwrite = TRUE,
                    job_name = "iter_ct",
                    sbatch_opt = list(time = "02:00:00"),
                    plan = "wait")

job

res <- Slurm_collect(job)

x <-  Filter((function(l) length(l) == 1), res)
if (length(x) > 0) print(x[[1]])

x <- Filter((function(l) length(l) > 1), res)
print(paste("Number of records with valid estimation:", length(x)))

ests <- rbindlist(x)

save(ests, file = "/gpfs/data/troxellab/ksg/data/cov_test_2.rda")
