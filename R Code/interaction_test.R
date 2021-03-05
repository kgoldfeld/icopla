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

rt <- stanc("./Stan Code/interaction_m2_who_enroll3levels.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

defC <- defDataAdd(varname = "a", formula = 0, variance = .005, 
                   dist = "normal")    
defC <- defDataAdd(defC, varname = "b", formula = 0, variance= .01, 
                   dist = "normal")                  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 

defC2 <- defDataAdd(
  varname="C_rv",
  formula="C * control",
  dist = "nonrandom"
) # C_rv=1/2/3: patient received control treatment C=1/2/3

defC2 <- defDataAdd(defC2,varname="sex", formula = 0.5, 
                    dist = "binary")   # sex=1:male;sex=0:female
defC2 <- defDataAdd(defC2,varname="who_enroll", formula = "1/3;1/3;1/3", 
                    dist = "categorical")   # 1==who 4, 2==who 5, 3 = who 6
defC2 <- defDataAdd(defC2,varname="age", formula = "0.25;0.25;0.50", 
                    dist = "categorical")   # 3 age categories

defC2 <- defDataAdd(defC2, varname = "ss", 
                    formula = "0.3;0.3;0.4", dist = "categorical")  #symptom duration strata

defS <- defCondition(
  condition = "ss==1",  
  formula = "0 * (C_rv==1) + 0 * (C_rv==2) + 0 * (C_rv==3)", 
  dist = "nonrandom")

defS <- defCondition(defS,
                     condition = "ss==2",  
                     formula = "(0.09 + a) * (C_rv==1) + (0.1 + a) * (C_rv==2) + (0.11 + a) * (C_rv==3)", 
                     dist = "nonrandom")
defS <- defCondition(defS,
                     condition = "ss==3",  
                     formula = "(0.19 + a) * (C_rv==1) + (0.20 + a) * (C_rv==2) + (0.21 + a) * (C_rv==3)", 
                     dist = "nonrandom")


defC3 <- defDataAdd(
  varname = "z", 
  formula = "0.05*(ss-1) + z_ss + 0.1*sex + 0.075*(age-1) + 0.06*(who_enroll-1) + (0.3 + b) * (C_rv==1) + (0.4 + b) * (C_rv==2) + (0.5 + b) * (C_rv==3)", 
  dist = "nonrandom")

nsites <- 9
basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060), 
                          similarity = 100) 

## Generate the data
dstudy <- genData(nsites, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "control") 
dind <- addColumns(defC2, dind)
dind <- addCondition(defS, dind, newvar = "z_ss")
dind <- addColumns(defC3, dind)

dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})

dind <- rbindlist(dl)  

##data for stan
N <- nrow(dind)                                ## number of observations 
L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
K <- dind[, length(unique(study))]            ## number of studies 
y <- as.numeric(dind$ordY)                    ## individual outcome
y_2 <- ifelse(y %in% seq(8,11),1,0)            ## individual logistic outcome 
kk <- dind$study                              ## study for individual 
ctrl <- dind$control                          ## treatment arm for individual 
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss), data = dind)[, -1]
D <- ncol(x)

##################ss=1 as the reference level###################
ds <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]


studydata <- list(
  N=N, L= L, K=K, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x)

##Fit the model
fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                 cores = 4L, chains = 4, control = list(adapt_delta = 0.8))

###############Table
s <- summary(fit, pars=c("Delta",paste0('GD[',1:2,']')), probs = c(0.025,0.5,0.975))$summary[,c(4,5,6)] #posterior median
OR_1 <- sapply(extract(fit,pars=c("Delta",paste0('GD[',1:2,']'))),function(x) mean(x>0))
OR_08 <- sapply(extract(fit,pars=c("Delta",paste0('GD[',1:2,']'))),function(x) mean(x>-log(0.8)))

results_or_lg1 <- cbind(exp(-s)[,3:1],OR_1,OR_08)

rownames(results_or_lg1) <- c('Level 1(ref)','Level 2', 'Level 3')
colnames(results_or_lg1) <- c("2.5%","50%","97.5","P(OR<1)","P(OR<0.8)")
results_or_lg1


##################ss=3 as the reference level###################
dind$ss <- relevel(as.factor(dind$ss),
                               ref = "3")

ds <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]


studydata <- list(
  N=N, L= L, K=K, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x)

##Fit the model
fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                 cores = 4L, chains = 4, control = list(adapt_delta = 0.8))

###############Table
sparams <- get_sampler_params(fit, inc_warmup=FALSE)
div_co1 <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))

s <- summary(fit, pars=c("Delta",paste0('GD[',1:2,']')), probs = c(0.025,0.5,0.975))$summary[,c(4,5,6)] #posterior median
OR_1 <- sapply(extract(fit,pars=c("Delta",paste0('GD[',1:2,']'))),function(x) mean(x>0))
OR_08 <- sapply(extract(fit,pars=c("Delta",paste0('GD[',1:2,']'))),function(x) mean(x>-log(0.8)))

results_or_lg1 <- cbind(exp(-s)[,3:1],OR_1,OR_08)

rownames(results_or_lg1) <- c('Level 3(ref)','Level 1', 'Level 2')
colnames(results_or_lg1) <- c("2.5%","50%","97.5","P(OR<1)","P(OR<0.8)")
results_or_lg1
