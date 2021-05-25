library(simstudy)
library(cmdstanr)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)
library(lme4)
library(posterior)
library(bayesplot)
library(glue)


set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
mco <- cmdstan_model("/gpfs/home/dw2625/r/SIM/primary_co.stan")
ml <- cmdstan_model("/gpfs/home/dw2625/r/SIM/primary_l.stan")

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
  iter <- function(iternum, defC, defC2, defA1, basestudy, nsites, mco,ml)  {
    set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
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
      y_2 <- ifelse(y %in% seq(8,11),1,0) 
      kk <- sampdat$study                              ## study for individual
      ctrl <- sampdat$control                          ## treatment arm for individual
      cc <- sampdat[, .N, keyby = .(study, C)]$C       ## specific control arm for study
      prior_div <- 8   #SD of site-specific intercept
      prior_Delta_sd <- .354
      #prior_beta_sd <- 2.5
      eta <- 0.1      #sd of delta_c
      prior_eta_0 <- 0.25
      
      studydata <- list(
        N=N, L= L,K=K, y=y, y_2=y_2,ctrl=ctrl,cc=cc,kk=kk,prior_div=prior_div,
        prior_Delta_sd=prior_Delta_sd,eta=eta,
        prior_eta_0 = prior_eta_0)
      
      fit_co <- mco$sample(
        step_size = 0.1,
        data = studydata,
        chains = 4L,
        parallel_chains = 4L,
        refresh = 500,
        iter_warmup = 500,
        iter_sampling = 2500,
        adapt_delta=0.8
      )
      
      draws_co <- data.table(as_draws_df(fit_co$draws()))
      diagnostics_df <- as_draws_df(fit_co$sampler_diagnostics())
      div_co <- sum(diagnostics_df[, 'divergent__'])
      
      p.eff.co <- mean(draws_co[, OR] < 1)
      p.clinic.co <- mean(draws_co[, OR] < 0.8)
      Delta.med.co <- median(draws_co[, Delta])
      
      fit_l <- ml$sample(
        step_size = 0.1,
        data = studydata,
        chains = 4L,
        parallel_chains = 4L,
        refresh = 500,
        iter_warmup = 500,
        iter_sampling = 2500,
        adapt_delta=0.8
      )
      
      draws_l <- data.table(as_draws_df(fit_l$draws()))
      
      diagnostics_df <- as_draws_df(fit_l$sampler_diagnostics())
      div_l <- sum(diagnostics_df[, 'divergent__'])
      

      p.eff.l <- mean(draws_l[, OR] < 1)
      p.clinic.l <- mean(draws_l[, OR] < 0.8)
      Delta.med.l <- median(draws_l[, Delta])
      
      if(n %in% c(15,30,45,60,75)){
      #####Frequentist method: 5 looks
      sampdat$Y <- ifelse(sampdat$ordY %in% seq(8,11),1,0) 
      sampdat$ordY <- as.factor(sampdat$ordY)
      sampdat$control <- as.factor(sampdat$control)
      sampdat$study <- as.factor(sampdat$study)
      sampdat$C <- as.factor(sampdat$C)
      
      #####Frequentist PO model fitting
      freq_co <- clmm(formula = ordY ~ 1 + control + (1+control | C/study),
                    data=sampdat,link = "logit")#random int
      
      freq_co_out <- t(summary(freq_co)$coefficients[length(unique(sampdat$ordY)),])
      colnames(freq_co_out) <- paste0(colnames(freq_co_out),"_co")
      
      freq_l <- glmer(formula = Y ~ 1 + control + (1+control | C/study),
                      data=sampdat,family="binomial")
      
      freq_l_out <- t(summary(freq_l)$coefficients[length(unique(sampdat$Y)),])
      colnames(freq_l_out) <- paste0(colnames(freq_l_out),"_l")
      
      data.table(iternum,n, p.eff.co, p.clinic.co,p.eff.l,p.clinic.l,Delta.med.co,Delta.med.l,freq_co_out,freq_l_out,
                 div_co,div_l)
      
      }else{
        freq_co_out <- matrix("NULL",nrow=1,ncol=4)#other looks: only run Bayesian methods
        colnames(freq_co_out) <- c("Estimate_co",   "Std. Error_co", "z value_co","Pr(>|z|)_co")
        freq_l_out <- matrix("NULL",nrow=1,ncol=4)#other looks: only run Bayesian methods
        colnames(freq_l_out) <- c("Estimate_l",   "Std. Error_l", "z value_l","Pr(>|z|)_l")
        
        data.table(iternum,n, p.eff.co, p.clinic.co,p.eff.l,p.clinic.l,Delta.med.co,Delta.med.l,freq_co_out,freq_l_out,
                   div_co,div_l)
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
    
    #interim look: 20%;33%;
   rbindlist(lapply(c(15,25,30,38,45,50,60,68,75),function(x) dofit(iternum,x,dxfull=dind))) 
   # rbindlist(lapply(c(25,38,50,68,75),function(x) dofit(iternum,x,dxfull=dind)))
  }
### Stan model



# rt <- stanc("./Stan/primary_co.stan");
# mco <- stan_model(stanc_ret = rt, verbose=FALSE)
# 
# rt <- stanc("./Stan/primary_l.stan");
# ml <- stan_model(stanc_ret = rt, verbose=FALSE)

#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= .01, 
                   dist = "normal")    #each study has a random treatment effect
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                    formula = "(0.4 + b ) * (C_rv==1) + (0.5 + b ) * (C_rv==2) + (0.6 + b ) * (C_rv==3)", 
                    dist = "nonrandom")


nsites <- 9

basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060),
                          similarity = 100)


NJOBS <- 90
set.seed(382332)
seeds <- sample(1:1000000,NJOBS,replace=FALSE)

job <- Slurm_lapply(1:2520,
                    iter, 
                    defC=defC,
                    defC2 = defC2,
                    defA1 =defA1,
                    basestudy=basestudy,
                    nsites=nsites,
                    mco = mco,
                    ml=ml,
                    njobs = 90, 
                    mc.cores = 4,
                    seeds=seeds,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "sim_3",
                    sbatch_opt = list(time = "30:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data

site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/SIM/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/SIM/", date_stamp, "/interim_456_2520.rda"))



