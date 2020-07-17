library("simstudy")
library("rstan")

#### Data definitions

defC <- defDataAdd(varname = "a", formula = 0, variance = 2, 
                   dist = "normal")    #each study has a random intercept
defC <- defDataAdd(defC, varname = "b", formula = 0, variance= .5, 
                   dist = "normal")                  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "50+50*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                    formula = "a + (0.5 + b ) * (C_rv==1) + (1 + b ) * (C_rv==2) + (2 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

baseprobs <- c(.20, .30, .35, .15) 

#### Generate data

set.seed(18459)

dstudy <- genData(9, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "control") #add control treatment variable, randomly assign control treatment in each study;rx=1:received control;rx=0:received convalescent plasma
dind <- addColumns(defC2,dind)
dind <- addColumns(defA1, dind) # add z
dind <- genOrdCat(dind, adjVar = "z", baseprobs, catVar = "ordY") #add outcome's level

###

N = nrow(dind) ;                              ## number of observations
L <- dind[, length(unique(ordY))]             ## of levels of outcome
K <- dind[, length(unique(study))]            ## number of studies
y <- as.numeric(dind$ordY)                    ## individual outcome
kk <- dind$study                              ## study for individual
ctrl <- dind$control                          ## treatment arm for individual
cc <- dind$C                                  ## specific control arm

studydata <- list(N=N, L=L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc)

##

rt <- stanc("Programs/Group/plasma_ma.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

Sys.time()
fit <-  sampling(sm, data=studydata, seed = 1327, iter = 2000, warmup = 500, cores = 4L)
Sys.time()

pars <- c("b", "delta", "tau")
print(fit, pars = pars, probs = c(0.05, 0.5, 0.95))
plot(fit, plotfun = "trace", pars = "b", 
     inc_warmup = FALSE, ncol = 3)

