library(simstudy)
library(data.table)
library(lme4)
library(slurmR)

iter <- function(param) {
  
  # generate data based on definitions and parameter
  
  d1 <- defData(varname = "a", formula= 0, variance = param["bvar"], id = "grp")
  
  d2 <- defDataAdd(varname = "x", formula = 0, variance = 16)
  d2 <- defDataAdd(d2, varname = "y", formula = "a + x", variance = param["wvar"])
  
  dx <- genData(param["nclust"], d1)
  dx <- genCluster(dx, cLevelVar = "grp", param["npclust"], "id")
  dx <- addColumns(d2, dx)
  
  # estimate model and return results
  
  lmefit <- lmer(y ~ x + (1|grp), data = dx)
  data.table(t(param), est = coef(summary(lmefit))["x",1])
  
}

### Set simulation parameters

sims <- 1:1080
nclust <- c(50, 100)
npclust <- c(30, 60)
wvar <- c(4, 9, 16)
bvar <- c(1, 1.5, 2)

dparam <- expand.grid(sim = sims, nclust = nclust, npclust = npclust, wvar=wvar, bvar= bvar)
lparam <- asplit(dparam, 1)

#### Generate results for each set (row) of paramers

job <- Slurm_lapply(lparam, iter, 
                    njobs = 50, 
                    mc.cores = 1,
                    tmp_path = "/gpfs/scratch/goldfk01",
                    job_name = "lme_iter",
                    sbatch_opt = list(time = "00:05:00"),
                    plan = "wait")

res <- Slurm_collect(job)
lme_iter <- rbindlist(res)  # creates a single data.table from the list
save(lme_iter, file = "/gpfs/scratch/goldfk01/lme_iter.rda")


