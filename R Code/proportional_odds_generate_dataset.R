library(simstudy)
library(cmdstanr)
library(data.table)
library(slurmR)
library(ordinal)
library(dplyr)
library(glue)
library(collapse)
library(posterior)


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


############Define data###############
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
defC2 <- defDataAdd(defC2, varname = "ds", 
                    formula = "0.2;0.2;0.2;0.2;0.2", dist = "categorical")  #symptom duration strata

defC2 <- defDataAdd(defC2, varname = "ss", 
                    formula = "0.3;0.3;0.4", dist = "categorical")  #covariates

defS <- defCondition(
  condition = "ss==1",  
  formula = "(0.04 + a) * (C_rv==1) + (0.05 + a) * (C_rv==2) + (0.06 + a) * (C_rv==3)", 
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
  formula = "0.05*(ss-1)+0.05*(ds-1) + z_ss + 0.1*sex + 0.075*(age-1) + 0.06*(who_enroll-1) + (0.3 + b) * (C_rv==1) + (0.4 + b) * (C_rv==2) + (0.5 + b) * (C_rv==3)", 
  dist = "nonrandom")

nsites <- 9
basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060), 
                          similarity = 100) 

############Data generation################
dstudy <- genData(nsites, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "control") 
dind <- addColumns(defC2, dind)
dind <- addCondition(defS, dind, newvar = "z_ss")
dind <- addColumns(defC3, dind)


############Under proportionality assumption##############
setkey(dind, "id")
dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  dx <- genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
  dx[, ordY := factor(ordY, levels = c(1:11))]
  dx[]
})

dind <- rbindlist(dl)  

save(dind,file="./data/dind_local_v2.rda")

load(file="./data/dind_local_v2.rda")

##check if proportionality holds in study 3
#Calculation of the observed cumulative odds ratio at each response level 
#doesn't provide an entirely clear picture about proportionality,
#but the sample size is relatively small given the number of categories.
data <- subset(dind,study==3)
prop <- prop.table(data[, table(sex, ordY)], 1)
cumprop <- data.table(apply(prop, 1, cumsum))
cumodds <- cumprop[, .(odds0 = `0`/(1 - `0`), odds1=`1`/(1 - `1`))]
cumodds[1:10, odds1/odds0]

#--- Create a data set ---#
dt_to_list <- function(dx) {
  
  N <- nrow(dx)                               ## number of observations 
  L <- dx[, length(unique(ordY))]             ## number of levels of outcome 
  K <- dx[, length(unique(study))]            ## number of studies 
  y <- as.numeric(dx$ordY)                    ## individual outcome 
  kk <- dx$study                              ## study for individual 
  ctrl <- dx$control                          ## treatment arm for individual 
  cc <- dx[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ds) + factor(ss), data = dx)[, -1]
  D <- ncol(x)
  S <-dind[, length(unique(ss))]
  ##################ss=1 as the reference level###################
  ss <- dind$ss
  
  list( N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc, ss=ss,D=D, x=x,S=S)
}

##############Fit the model##################
mod <- cmdstan_model("./interaction_co_interactionSlevels.stan")

fit_co <- mod$sample(
  data = dt_to_list(dind),
  refresh = 0,
  chains = 4L,
  parallel_chains = 4L,
  iter_warmup = 500,
  iter_sampling = 2500,
  step_size = 0.1,
  show_messages = FALSE,
  adapt_delta=0.98,
  max_treedepth = 15
)

save(fit_co,file="./data/fit_co.rda")

#######draw from posterior iteration#################
draws_dt <- data.table(as_draws_df(fit_co$draws()))

save(draws_dt,file="./data/draws_dt.rda")

