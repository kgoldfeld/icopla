library("simstudy")
library("rstan")
library("data.table")
library("slurmR")

### Function to generate base probabilities for each study

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

#### Generate data


iter <- function(iternum, defC, defC1, defA1, basestudy,nsites, sm) {

dstudy <- genData(nsites, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "control") #add control treatment variable, randomly assign control treatment in each study;rx=1:received control;rx=0:received convalescent plasma
dind <- addColumns(defC2,dind)
dind <- addColumns(defA1, dind) # add z


dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})

dind <- rbindlist(dl)

### Data for Stan model

N = nrow(dind) ;                              ## number of observations
L <- dind[, length(unique(ordY))]             ## number of levels of outcome
K <- dind[, length(unique(study))]            ## number of studies
y <- as.numeric(dind$ordY)                    ## individual outcome
kk <- dind$study                              ## study for individual
ctrl <- dind$control                          ## treatment arm for individual
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study

studydata <- list(N=N, L=L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc)

fit <-  sampling(sm, data=studydata, seed = 1327, iter = 3000, warmup = 500, cores = 4L)

pars <- c("delta_k", "eta_0","delta", "eta", "Delta", "tau")

##Option 1: print all parameters
#s <- summary(fit, pars = pars, probs = c(0.05, 0.5, 0.95))$summary
#data.table(iternum,rownames(s),s)

##Option 2: export statistics of interest
p.eff <- mean(extract(fit, pars = "OR")[[1]] < 1)
p.clinic <- mean(extract(fit, pars = "OR")[[1]] < 0.8)
Delta.mean <- mean(extract(fit, pars = 'Delta')[[1]]) #control treatment effect
eta.mean <- mean(extract(fit, pars = 'eta')[[1]])
eta_0.mean <- mean(extract(fit, pars = 'eta_0')[[1]])
data.table(iternum, p.eff, p.clinic,Delta.mean,eta.mean,eta_0.mean)
}

### Stan model
rt <- stanc("C:\\Users\\aaron\\Desktop\\Covid_19_project\\Plasma_Git\\icopla\\Stan Code\\plasma_ma_delta.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

set.seed(18459)

#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= .5, 
                   dist = "normal")                  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                    formula = "(0.5 + b ) * (C_rv==1) + (1 + b ) * (C_rv==2) + (2 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

nsites <- 9

###Study specific base probabilities

basestudy <- genBaseProbs(n = nsites, 
                          base =  c(.20, .30, .30, .20) ,
                          similarity = 1e3)

#test on laptop 
#job <- lapply(1:3,iter,defC=defC, defC1=defC1, defA1=defA1, basestudy=basestudy,nsites=nsites, sm)

NJOBS <- 990
set.seed(382332)
seeds <- sample(1:1000000,NJOBS,replace=FALSE)

job <- Slurm_lapply(1:NJOBS,
                    iter, 
                    defC=defC,
                    defC1 = defC1,
                    defA1 =defA1,
                    basestudy=basestudy,
                    nsites=nsites,
                    sm = sm,
                    njobs = 90, 
                    mc.cores = 4,
                    seeds=seeds,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "bayessd",
                    sbatch_opt = list(time = "12:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data

site_bayesord <- Slurm_collect(job) # data is a list
site_bayesord <- rbindlist(site_bayesord) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("~/R/", date_stamp), showWarnings = FALSE)
save(site_bayesord, file = paste0("~/R/", date_stamp, "/bayes_plasma.rda"))