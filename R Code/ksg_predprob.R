library(cmdstanr)
library(posterior)
library(bayesplot)
library(glue)
library(simstudy)
library(data.table)
library(slurmR)

#--- Simulate a single data set ---#

defC <- defDataAdd(varname = "a", formula = 0, variance = .005, 
  dist = "normal")    
defC <- defDataAdd(defC, varname = "b", formula = 0, variance= .01, 
  dist = "normal")  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
  dist = "nonrandom") 

defC2 <- defDataAdd(varname="C_rv", formula="C * control", dist = "nonrandom")
defC2 <- defDataAdd(defC2, varname="sex", formula = 0.5, 
  dist = "binary")   # sex=1:male;sex=0:female
defC2 <- defDataAdd(defC2, varname="who_enroll", formula = "1/3;1/3;1/3", 
  dist = "categorical")   # 1==who 4, 2==who 5, 3 = who 6
defC2 <- defDataAdd(defC2, varname="age", formula = "0.25;0.25;0.50", 
  dist = "categorical")   # 3 age categories
defC2 <- defDataAdd(defC2, varname = "ss", formula = "0.25;0.25;0.25;0.25", 
  dist = "categorical")  #symptom duration strata
defC2 <- defDataAdd(defC2, varname = "z", 
  formula = "0.05*(ss-1) + 0.1*sex + 0.075*(age-1) + 0.06*(who_enroll-1) + (0.3 + b) * (C_rv==1) + (0.4 + b) * (C_rv==2) + (0.5 + b) * (C_rv==3)", 
  dist = "nonrandom")

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

nsites <- 9
basestudy <- genBaseProbs(n = nsites,
  base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060), 
  similarity = 100)

set.seed(381273)

dstudy <- genData(nsites, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)
  
dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "control") 
dind <- addColumns(defC2, dind)
  
dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})
  
dind <- rbindlist(dl)  
  
dind <- dind[, .(id, study, C, control, sex, who_enroll, age, ss, z, ordY)]
  
#--- Create a data set with observed outcomes ---#

defM <- defMiss(varname = "ordY", formula = 0.30, logit.link = FALSE)
  
missMat <- genMiss(dind, defM, idvars = "id")
dind <- genObs(dind, missMat, idvars = "id")
  
dobs <- dind[!is.na(ordY)]

#--- Fit Bayes model to original data ---#

dt_to_list <- function(dx) {
  
  N <- nrow(dx)                               ## number of observations 
  L <- dx[, length(unique(ordY))]             ## number of levels of outcome 
  K <- dx[, length(unique(study))]            ## number of studies 
  y <- as.numeric(dx$ordY)                    ## individual outcome 
  kk <- dx$study                              ## study for individual 
  ctrl <- dx$control                          ## treatment arm for individual 
  cc <- dx[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss), data = dx)[, -1]
  D <- ncol(x)
  
  list(N=N, L=L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc, x=x, D=D)
}

set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
mod <- cmdstan_model("/gpfs/data/troxellab/ksg/r/predprob.stan")

fit <- mod$sample(
  step_size = 0.1,
  data = dt_to_list(dobs),
  chains = 4L,
  parallel_chains = 4L,
  refresh = 500,
  iter_warmup = 2000,
  iter_sampling = 2500
)
  
#--- Save posterior distribution ---#
  
draws_dt <- data.table(as_draws_df(fit$draws()))
mean(draws_dt[, exp(-Delta)] < 1)
# data.table(nuts_params(fit))[Parameter == "divergent__", mean(Value)]

#--- use original with missing outcome - represents "new" data ---#

newdt <- dind[is.na(ordY)]
newdt[, calt := control * study]
newdt[, ordY := NULL]
  
#--- Function to for generating outcomes and estimating model for new data ---#

est_from_draw <- function(n_draw, Draws, dt_obs, dt_new, D, K, s_model) {
  
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  
  draw <- as.data.frame(Draws[sample(.N, 1)])
  
  alpha <- draw$alpha
  beta <- as.vector(x = draw[, glue("beta[{1:D}]")], mode = "numeric")
  delta_k <- as.vector(draw[, glue("delta_k[{1:K}]")], mode = "numeric")
  
  tau <- as.vector(draw[grep("^tau", colnames(draw))], mode = "numeric")
  tau <- matrix(tau, nrow = K)
  tau <- cbind(tau, Inf)
  cprop <- t(apply(tau, 1, stats::plogis))
  xprop <- t(apply(cprop, 1, diff))
  baseline <- cbind(cprop[,1], xprop) 
  
  coefs <- c(alpha, beta, delta_k)
  x <- model.matrix(~factor(who_enroll) + factor(age) + 
                          factor(sex) + factor(ss), data = dt_new)
  calt <- model.matrix(~factor(calt) - 1, data = dt_new)[,-1]
  zmat <- cbind(x, calt)
  dt_new$z <- zmat %*% coefs
  dt_new[, calt := NULL]
  
  dl <- lapply(1:K, function(i) {
    b <- baseline[i,]
    dx <- dt_new[study == i]
    genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
  })
  
  dt_new <- rbindlist(dl)  
  
  dx <- rbind(dt_obs, dt_new)
  
  # fit model for combined data set of "complete" data and simulated new data
  
  fit_pp <- s_model$sample(
    step_size = 0.1,
    data = dt_to_list(dx),
    chains = 4L,
    parallel_chains = 4L,
    refresh = 500,
    iter_warmup = 2000,
    iter_sampling = 2500
  )
  
  draws_pp <- data.table(as_draws_df(fit_pp$draws()))
  
  prob_success <- mean(draws_pp[, exp(-Delta)] < 1)
  return(data.table(n_draw, prob_success))
  
  # dtNuts <- data.table(nuts_params(fit_pp))
  # pct_div <- dtNuts[Parameter == "divergent__", mean(Value)]
  
  # return(data.table(prob_success, pct_div))
}

#--- Multiple draws from the posterior ---#

job <- Slurm_lapply(
  X = 1L:1080L, 
  FUN = est_from_draw, 
  Draws = draws_dt,
  dt_obs = dobs,
  dt_new = newdt,
  D = dt_to_list(dobs)$D,
  K = dt_to_list(dobs)$K,
  s_model = mod,
  njobs = 90L, 
  mc.cores = 4L,
  job_name = "i_pp",
  tmp_path = "/gpfs/data/troxellab/ksg/scratch",
  plan = "wait",
  sbatch_opt = list(time = "03:00:00", partition = "cpu_short"),
  export = c("dt_to_list"),
  overwrite = TRUE
)

job
res <- Slurm_collect(job)

save(res, file = "/gpfs/data/troxellab/ksg/r/predprob.rda")
  