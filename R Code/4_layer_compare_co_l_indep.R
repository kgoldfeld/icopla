library(simstudy)
library(rstan)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)
library(MASS)
library(lme4)

genBaseProbs <- function(n, base, similarity, digits = 2) {
  
  # n: number of studies
  n_levels <- length(base) ## who level
  
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

iter <- function(iternum, defC, defC2,defC3,defC4, basestudy, nRCTs, mco,ml,mpl,mpco){
  
  ## Generate the data
  
  dstudy <- genData(nRCTs, id = "RCTs")       # generate RCTs
  dstudy <- trtAssign(dstudy,nTrt=3,grpName = "C") #C: 2 control treatment
  ##dstudy[, .N, keyby = .(C,RCTs)] #each control has 3 RCT 
  dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1)) #large=1:RCT with size 150; large=0:size 75
  dstudy <- addColumns(defC, dstudy) # add each RCT's random effect
  ##dstudy[, .N, keyby = .(C,large)]
  
  gen.sites <- defData(varname = "nsites",formula = 3,dist="nonrandom")#number of sites under each RCT
  dstudy <- addColumns(gen.sites, dstudy) # number of sites under each RCT
  dsite <- genCluster(dstudy,"RCTs",numIndsVar = "nsites",level1ID = "idsites")
  
  dsite <- addColumns(defC2, dsite) #Add each site's random effect
  dsite <- trtAssign(dsite, nTrt = 2, strata = "RCTs", grpName = "largesite", ratio = c(2,1)) #each RCT has 2 small sites and 1 large site
  ###dsite[, .N, keyby = .(RCTs,largesite)]
  
  dsite <- addCondition(defC3, dsite, newvar = "site_n") #the number of subject in each site
  ###dsite[, .N, keyby = .(RCTs,largesite)]
  ###dsite[, .N, keyby = .(RCTs,site_n)]
  
  dind <- genCluster(dsite, "idsites", numIndsVar = "site_n", "id")
  
  dind <- trtAssign(dind, strata="idsites", grpName = "control") #add control treatment variable, randomly assign control treatment in each site;rx=1:received control;rx=0:received convalescent plasma
  ###dind[, .N, keyby = .(idsites,control)]
  #####################################
  dind <- addColumns(defC4, dind)
  
  dl <- lapply(1:nRCTs, function(i) {
    b <- basestudy[i,]
    dx <- dind[RCTs == i]
    genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
  })
  
  
  dind <- rbindlist(dl)
  
  
  ###Data for stan
  
  N <- nrow(dind)                                ## number of observations 
  L <- dind[, length(unique(ordY))]             ## number of levels of outcome
  S <- dind[, length(unique(idsites))]         ##number of sites
  
  y_1 <- as.numeric(dind$ordY)                      ##outcome for PO model
  y_2 <- ifelse(dind$ordY %in% seq(8,11),1,0)         ##outcome 
  ctrl <- dind$control                          ## treatment arm for individual: 0/1
  cc <- dind[, .N, keyby = .(RCTs, C)][,C]      ## specific control arm for RCTs:1/2/3
  K <- dind[, length(unique(RCTs))]            ## number of RCTs
  kk <- dind[ ,RCTs]                              ## RCTs for individual
  site <- dind[,idsites]                        ##site for individual
  site_rct <- dind[, .N, keyby = .(idsites, RCTs)][, RCTs]
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss), data = dind)[, -1]
  D <- ncol(x)
  
  prior_div <- 8   #SD of site-specific intercept
  prior_Delta_sd <- .354
  prior_beta_sd <- 2.5
  eta <- 0.1 #sd of delta_c
  prior_eta_0 <- 0.25 #sd of eta_0:sd of delta_k
  
  studydata <- list(
    N=N,K=K,S=S,L=L, y_1=y_1,y_2=y_2, ctrl=ctrl,cc=cc,kk=kk, site=site, site_rct=site_rct,
    x=x,D=D,prior_div=prior_div,
    prior_Delta_sd=prior_Delta_sd,prior_beta_sd=prior_beta_sd, eta=eta,
    prior_eta_0 = prior_eta_0)
  
  ptm_l <- Sys.time()
  fit_l <-  sampling(ml, data=studydata,iter = 3000, warmup = 500, 
                   cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  paft_l <- Sys.time() 

  time_l <- difftime( paft_l,  ptm_l, units='min')
    
  ptm_co <- Sys.time()
  fit_co <-  sampling(mco, data=studydata,iter = 3000, warmup = 500, 
                     cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  paft_co <- Sys.time() 
  time_co <- difftime(paft_co,  ptm_co, units='min')
  
  ptm_pl <- Sys.time()
  fit_pl <-  sampling(mpl, data=studydata,iter = 3000, warmup = 500, 
                     cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  paft_pl <- Sys.time() 
  time_pl <- difftime( paft_pl,  ptm_pl, units='min')
  
  ptm_pco <- Sys.time()
  fit_pco <-  sampling(mpco, data=studydata,iter = 3000, warmup = 500, 
                      cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  paft_pco <- Sys.time() 
  time_pco <- difftime(paft_pco ,  ptm_pco, units='min')
  
  ##Bayesian output
  sparams <- get_sampler_params(fit_l, inc_warmup=FALSE)
  div_l <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  sparams <- get_sampler_params(fit_co, inc_warmup=FALSE)
  div_co <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  
  sparams <- get_sampler_params(fit_pl, inc_warmup=FALSE)
  div_pl <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  sparams <- get_sampler_params(fit_pco, inc_warmup=FALSE)
  div_pco <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  
  ## with site as an extra layer 
  p.eff.l <- mean(extract(fit_l, pars = "OR")[[1]] <1)
  p.clinic.l <- mean(extract(fit_l, pars = "OR")[[1]] < 0.8)
  p.eff.co<- mean(extract(fit_co, pars = "OR")[[1]] <1)
  p.clinic.co <- mean(extract(fit_co, pars = "OR")[[1]] < 0.8)
  
  p.eff.pl <- mean(extract(fit_pl, pars = "OR")[[1]] <1)
  p.clinic.pl <- mean(extract(fit_pl, pars = "OR")[[1]] < 0.8)
  p.eff.pco<- mean(extract(fit_pco, pars = "OR")[[1]] <1)
  p.clinic.pco <- mean(extract(fit_pco, pars = "OR")[[1]] < 0.8)
  

  p <- data.table(p.eff.l,p.clinic.l, p.eff.co,p.clinic.co,p.eff.pl,p.clinic.pl, p.eff.pco,p.clinic.pco)
  
  ## [1,] and [2,]: simulation results from the model with site as an extra layer:[1,]:CO;[2,]:logistic
  ## [3,] and [4,]: simulation results from the model without site as an extra layer: [3,]:CO;[4,]:logistic
  
  sl <- summary(fit_l, pars = c("OR","Delta","delta","beta"), probs = 0.5)
  slt <- data.table(t(sl$summary)[c(1,4,3),]) ##output: posterior mean(1) and median(4)and sd(3)
  colnames(slt) <- paste0(colnames(slt),"l")
  
  sco  <- summary(fit_co, pars = c("OR","Delta","delta","beta","alpha"), probs = 0.5)
  scot <- data.table(t(sco$summary)[c(1,4,3),]) ##output: posterior mean(1) and median(4)and sd(3)
  colnames(scot) <- paste0(colnames(scot),"co")
  
  spl <- summary(fit_pl, pars = c("OR","Delta","delta","beta"), probs = 0.5)
  splt <- data.table(t(spl$summary)[c(1,4,3),]) ##output: posterior mean(1) and median(4)and sd(3)
  colnames(splt) <- paste0(colnames(splt),"pl")
  
  spco  <- summary(fit_pco, pars = c("OR","Delta","delta","beta","alpha"), probs = 0.5)
  spcot <- data.table(t(spco$summary)[c(1,4,3),]) ##output: posterior mean(1) and median(4)and sd(3)
  colnames(spcot) <- paste0(colnames(spcot),"pco")
  
  data.table(cbind(iternum,slt,scot,splt,spcot,p,div_l,div_co,div_pl,div_pco,time_l,time_co,time_pl,time_pco))
}

rc <- stanc("/gpfs/home/dw2625/r/add_site_co_compare.stan");
#rc <- stanc("./add_site_co_compare.stan");
mco <- stan_model(stanc_ret = rc, verbose=FALSE)

rl <- stanc("/gpfs/home/dw2625/r/add_site_l_compare.stan");
#rl <- stanc("./add_site_l_compare.stan");
ml <- stan_model(stanc_ret = rl, verbose=FALSE)

pl <- stanc("/gpfs/home/dw2625/r/primary_l.stan");
#pl <- stanc("./primary_l.stan");
mpl <- stan_model(stanc_ret = pl, verbose=FALSE)

pco <- stanc("/gpfs/home/dw2625/r/primary_co.stan");
#pco <- stanc("./primary_co.stan");
mpco <- stan_model(stanc_ret = pco, verbose=FALSE)

nRCTs <- 9

basestudy <- genBaseProbs(n = nRCTs,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060),
                          similarity = 100) ##assume each RCTs has its specific intercept

#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= .01, 
                   dist = "normal")    #each RCT+site random treatment effect

defC2 <- defDataAdd(varname = "r", formula = 0, variance= 0.0625, 
                    dist = "normal")    #each site has a random intercept add on RCT's specific intercept: N(0,sd=0.25)
defC2 <- defDataAdd(defC2,varname = "a", formula = 0, variance= 0.01, 
                    dist = "normal")   #each site has a random slope

defC3 <- defCondition(condition = "large==1",  
                      formula = "2*(10 + 45*largesite)", 
                      dist = "nonrandom") ##RCTs with size 150: 3 sites with size 20,20 and 110
defC3 <- defCondition(defC3,condition = "large==0",  
                      formula = "(10 + 45*largesite)", 
                      dist = "nonrandom") ##RCTs with size 75: 3 sites with size 10,10 and 55

defC4 <- defDataAdd(varname="C_rv", formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3

defC4 <- defDataAdd(defC4,varname="sex", formula = 0.5, 
                    dist = "binary")   # sex=1:male;sex=0:female
defC4 <- defDataAdd(defC4,varname="who_enroll", formula = "1/3;1/3;1/3", 
                    dist = "categorical")   # 1==who 4, 2==who 5, 3 = who 6
defC4 <- defDataAdd(defC4,varname="age", formula = "0.25;0.25;0.50", 
                    dist = "categorical")   # 3 age categories

defC4 <- defDataAdd(defC4, varname = "ss", 
                    formula = "0.25;0.25;0.25;0.25", dist = "categorical")  #symptom duration strata

defC4 <- defDataAdd(defC4, varname = "z", 
                    formula = "r + 0.05*(ss-1) + 0.1*sex + 0.075*(age-1) + 0.06*(who_enroll-1) + (0.3 + b + a ) * (C_rv==1) + (0.4 + b + a ) * (C_rv==2) + (0.5 + b + a) * (C_rv==3)", 
                    dist = "nonrandom")

job <- Slurm_lapply(1:1080,
                    iter, 
                    defC=defC,
                    defC2 = defC2,
                    defC3 = defC3,
                    defC4 = defC4,
                    basestudy=basestudy,
                    nRCTs=nRCTs,
                    mco = mco,
                    ml=ml,
                    mpl=mpl,
                    mpco=mpco,
                    njobs = 90, 
                    mc.cores = 4,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_254",
                    sbatch_opt = list(time = "24:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data
site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/", date_stamp, "/4model_compare_indep_v2.rda"))

