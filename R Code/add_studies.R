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

#### Generate data and fit model
iter <- function(iternum, defC, defC2, basestudy, nsites, sm){
  
  ## Generate the data
  dstudy <- genData(nsites, id = "study")       # generate studies
  dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
  dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
  dstudy <- addColumns(defC, dstudy)
  
  dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
  dind <- trtAssign(dind, strata="study", grpName = "control") #add control treatment variable, randomly assign control treatment in each study;rx=1:received control;rx=0:received convalescent plasma
  dind <- addColumns(defC2, dind)
  
  
  dl <- lapply(1:nsites, function(i) {
    b <- basestudy[i,]
    dx <- dind[study == i]
    genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
  })
  
  
  dind <- rbindlist(dl)
  
  ##data for Stan
  N <- nrow(dind)                                ## number of observations 
  L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
  y_1 <- as.numeric(dind$ordY)                  ## WHO 11 points scale 
  y_2 <- ifelse(y_1 %in% seq(8,11),1,0)         ##WHO 11 points scale 7-10 
  ctrl <- dind$control                          ## treatment arm for individual 
  cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study
  K <- dind[, length(unique(study))]            ## number of studies
  kk <- dind$study                              ## study for individual
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss), data = dind)[, -1]
  D <- ncol(x)
  
  prior_div <- 8   #SD of site-specific intercept
  prior_Delta_sd <- .354
  prior_beta_sd <- 2.5
  eta <- 0.1 #sd of delta_c
  prior_eta_0 <- 0.25 
  
  studydata <- list(
    N=N, L= L,K=K, y_1=y_1,y_2=y_2, ctrl=ctrl,cc=cc,kk=kk, x=x,D=D,prior_div=prior_div,
    prior_Delta_sd=prior_Delta_sd, prior_beta_sd= prior_beta_sd,eta=eta,
    prior_eta_0 = prior_eta_0)
  
  ##Fit the model
  fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                   cores = 4L, chains = 4, control = list(adapt_delta = 0.8))
  
  ##Bayesian output
  sparams <- get_sampler_params(fit, inc_warmup=FALSE)
  div <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  
  p.eff.1 <- mean(extract(fit, pars = "OR[1]")[[1]] <1)
  p.clinic.1 <- mean(extract(fit, pars = "OR[1]")[[1]] < 0.8)
  p.eff.2<- mean(extract(fit, pars = "OR[2]")[[1]] <1)
  p.clinic.2 <- mean(extract(fit, pars = "OR[2]")[[1]] < 0.8)
  
  p <- data.table(p.eff.1,p.clinic.1, p.eff.2,p.clinic.2)
  
  s <- summary(fit, pars = c("alpha", "Delta","delta","beta"), probs = 0.5)
  st <- data.table(t(s$summary)[c(1,4),]) ##output: posterior mean(1) and median(4)
  
  
  ##Data for frequentist method
  dind$y_1 <- factor(studydata$y_1)
  dind$y_2 <- factor(studydata$y_2)
  dind$C <- factor(dind$C)
  dind$ctrl <- factor(dind$control)
  dind$study <- factor(dind$study)
  dind$who_enroll <- factor(dind$who_enroll)
  dind$age <- factor(dind$age)
  dind$sex <- factor(dind$sex)
  dind$ss <- factor(dind$ss)
  
  ##Frequentist output
  fit_po <- clmm(formula = y_1 ~ ctrl+ (1+ctrl|C/study) + who_enroll+age+sex+ss,
                 data=dind,link = "logit")#random int
  freq_po <- t(summary(fit_po)$coefficients[11:19,c(1,4)])##the first row: estimation; the second row: p_value
  colnames(freq_po) <- c( "ctrl1_PO", "who_enroll2_PO", "who_enroll3_PO", "age2_PO",
                          "age3_PO","sex1_PO", "ss2_PO","ss3_PO","ss4_PO"  )
  
  fit_lg <- glmer(y_2 ~ ctrl+ (1+ctrl|C/study) + who_enroll + age + sex + ss,data=dind,family=binomial)
  freq_lg <- t(summary(fit_lg)$coefficients[-1,c(1,4)])
  colnames(freq_lg) <- c( "ctrl1_LG", "who_enroll2_LG", "who_enroll3_LG", "age2_LG",
                          "age3_LG","sex1_LG", "ss2_LG","ss3_LG","ss4_LG"  )
  
  data.table(iternum,cbind(p,st,div,freq_po,freq_lg))  
  
}
### Stan model

rt <- stanc("/gpfs/home/dw2625/r/add_studies.stan");
#rt <- stanc("./add_studies.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

#### Data definitions
defC <- defDataAdd(varname = "b", formula = 0, variance= .01, 
                   dist = "normal")    #each study has a random treatment effect
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 

defC2 <- defDataAdd(varname="C_rv", formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defC2 <- defDataAdd(defC2,varname="sex", formula = 0.5, 
                    dist = "binary")   # sex=1:male;sex=0:female
defC2 <- defDataAdd(defC2,varname="who_enroll", formula = "1/3;1/3;1/3", 
                    dist = "categorical")   # 1==who 4, 2==who 5, 3 = who 6
defC2 <- defDataAdd(defC2,varname="age", formula = "0.25;0.25;0.50", 
                    dist = "categorical")   # 3 age categories

defC2 <- defDataAdd(defC2, varname = "ss", 
                    formula = "0.25;0.25;0.25;0.25", dist = "categorical")  #symptom duration strata

###Assume male, older people, higher who_baseline do worse

# beta_duration = c(0, .05, .1, .15), beta_age = c(0, .075, .15) and beta_who = c(0, .06, .12) 
defC2 <- defDataAdd(defC2, varname = "z", 
                    formula = "0.05*(ss-1) + 0.1*sex + 0.075*(age-1) + 0.06*(who_enroll-1) + (0.3 + b ) * (C_rv==1) + (0.4 + b ) * (C_rv==2) + (0.5 + b ) * (C_rv==3)", 
                    dist = "nonrandom")


nsites <- 9

basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060),
                          similarity = 100)



job <- Slurm_lapply(1:990,
                    iter, 
                    defC=defC,
                    defC2 = defC2,
                    basestudy=basestudy,
                    nsites=nsites,
                    sm = sm,
                    njobs = 90, 
                    mc.cores = 4,
                    tmp_path = "/gpfs/scratch/dw2625",
                    job_name = "nc_212",
                    sbatch_opt = list(time = "8:00:00"),
                    plan = "wait",
                    overwrite=TRUE)

### Save data
site_plasma_all <- Slurm_collect(job) # data is a list
site_plasma <- site_plasma_all[lapply(site_plasma_all, function(x) length(x))>1] #filter out the error message
site_plasma <- rbindlist(site_plasma) # converting list to data.table

date_stamp <- gsub("-", "", Sys.Date())
dir.create(file.path("/gpfs/home/dw2625/r/", date_stamp), showWarnings = FALSE)
save(site_plasma, file = paste0("/gpfs/home/dw2625/r/", date_stamp, "/[3]_04_add_studies_endpoints.rda"))


