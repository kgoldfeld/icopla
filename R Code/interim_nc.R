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
  iter <- function(iternum, defC, defC2, defA1, basestudy, nsites, sm)  {
   ##size of interim 
    prepData <- function(n,dxfull){
      ##choose the first n patients in each study
      data_small <- dxfull%>% filter(size==75) %>% group_by(study) %>% do(head(.,n)) %>% as.data.table()
      data_large <- dxfull%>% filter(size==150) %>% group_by(study) %>% do(head(.,2*n)) %>% as.data.table()
      rbind(data_small,data_large)
      
    }
    
   ##fit the model in Frequentist and Bayesian way
    dofit <- function(iternum,n,dxfull){
      
      sampdat <- prepData(n,dxfull)
      N = nrow(sampdat)                             ## number of observations
      L <- sampdat[, length(unique(ordY))]             ## number of levels of outcome
      K <- sampdat[, length(unique(study))]            ## number of studies
      y <- as.numeric(sampdat$ordY)                    ## individual outcome
      kk <- sampdat$study                              ## study for individual
      ctrl <- sampdat$control                          ## treatment arm for individual
      cc <- sampdat[, .N, keyby = .(study, C)]$C       ## specific control arm for study
      eta <- .1 
      prior_eta_0 <- 0.25 
      prior_tau_sd <- 2.5
      prior_Delta_sd <- .354 
      
      studydata <- list(
        N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc, 
        prior_tau_sd = prior_tau_sd, prior_Delta_sd = prior_Delta_sd, prior_eta_0 = prior_eta_0, eta = eta)
      
      fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                       cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
      
      sparams <- get_sampler_params(fit, inc_warmup=FALSE)
      div <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
      p.eff <- mean(extract(fit, pars = "OR")[[1]] < 1)
      p.clinic <- mean(extract(fit, pars = "OR")[[1]] < 0.8)
      Delta.perc <- quantile(extract(fit, pars = 'Delta')[[1]], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) #control treatment effect
      
      if(n %in% c(25,38,50,68,75)){
      #####Frequentist method: 5 looks
      sampdat$ordY <- as.factor(sampdat$ordY)
      sampdat$control <- as.factor(sampdat$control)
      sampdat$study <- as.factor(sampdat$study)
      sampdat$C <- as.factor(sampdat$C)
      
      #####Frequentist PO model fitting
      fit_f <- clmm(formula = ordY ~ 1 + control + (1+control | C/study),
                    data=sampdat,link = "logit")#random int
      
      freq <- t(summary(fit_f)$coefficients[11,])
      data.table(iternum,n, p.eff, p.clinic,Delta.perc,freq,div)
      
      }else{
        freq <- matrix("NULL",nrow=1,ncol=4)#other looks: only run Bayesian methods
        colnames(freq) <- c("Estimate",   "Std. Error", "z value","Pr(>|z|)")
        data.table(iternum,n, p.eff, p.clinic,Delta.perc,freq,div) 
            }
    }

    ## Generate the data
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
    
    #interim look: 16%;33%;40%;50%;60%;67%;80%;90%;100%
   rbindlist(lapply(c(12,25,30,38,45,50,60,68,75),function(x) dofit(iternum,x,dxfull=dind))) 
   # rbindlist(lapply(c(25,38,50,68,75),function(x) dofit(iternum,x,dxfull=dind)))
  }
### Stan model

rt <- stanc("/gpfs/home/dw2625/r/explore_divergence_nc.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)



#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= .01, 
                   dist = "normal")    #each study has a random treatment effect
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                    formula = "(0 + b ) * (C_rv==1) + (0 + b ) * (C_rv==2) + (0 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

set.seed(184916)

nsites <- 9

basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060),
                          similarity = 100)


NJOBS <- 90
set.seed(382332)
seeds <- sample(1:1000000,NJOBS,replace=FALSE)

job <- Slurm_lapply(1:1440,
                    iter, 
                    defC=defC,
                    defC2 = defC2,
                    defA1 =defA1,
                    basestudy=basestudy,
                    nsites=nsites,
                    sm = sm,
                    njobs = 90, 
                    mc.cores = 4,
                    seeds=seeds,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_28",
                    sbatch_opt = list(time = "10:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data

site_plasma <- Slurm_collect(job) # data is a list
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/", date_stamp, "/interim_nc_000_1260_025eta0.rda"))



