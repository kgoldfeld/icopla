library("simstudy")
library("rstan")
library("data.table")

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

#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= .005, 
                   dist = "normal")                  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                    formula = "(0.8 + b ) * (C_rv==1) + (0.9 + b ) * (C_rv==2) + (1 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

#### Generate data

set.seed(184978)

nsites <- 9

dstudy <- genData(nsites, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "control") #add control treatment variable, randomly assign control treatment in each study;rx=1:received control;rx=0:received convalescent plasma
dind <- addColumns(defC2,dind)
dind <- addColumns(defA1, dind) # add z

###Study specific base probabilities

basestudy <- genBaseProbs(n = nsites,
                  base =  c(.1,.107,.095,.085,.09,.09,.108,.1,0.09,0.075, 0.06) ,
                  similarity = 100)

# basestudy <- genBaseProbs(n = nsites,
#                 base =  c(.15, .15, .2, .3, .2) ,
#                 similarity = 200)

dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})

dind <- rbindlist(dl)

### Check proportions

dind[, table(study, ordY, control)]

getProp <- function(d) {
  d[, round(prop.table(table(C_rv, ordY), margin = 1), 2)]
}

getProp(dind)

lapply(1:9, function(x) getProp(dind[study == x]))

###

N = nrow(dind) ;                              ## number of observations
L <- dind[, length(unique(ordY))]             ## number of levels of outcome
K <- dind[, length(unique(study))]            ## number of studies
y <- as.numeric(dind$ordY)                    ## individual outcome
kk <- dind$study                              ## study for individual
ctrl <- dind$control                          ## treatment arm for individual
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study

studydata <- list(N=N, L= L, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc)

##

# rt <- stanc("~\\plasma_ma_delta.stan")
# rt <- stanc("./Stan Code/plasma_ma_delta_induced_dir.stan");
rt <- stanc("./Stan Code/plasma_ma_delta.stan");

sm <- stan_model(stanc_ret = rt, verbose=FALSE)

fit <-  sampling(sm, data=studydata, seed = 1328, iter = 3000, warmup = 500, 
                 cores = 4L, chains = 4, init = 0)





