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

nRCTs <- 9

basestudy <- genBaseProbs(n = nRCTs,
                          base =c(0.175, 0.290, 0.275, 0.260),
                          similarity = 100) ##assume each RCTs has its specific intercept

#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= .01, 
                    dist = "normal")    #each RCT has a random treatment effect

defC2 <- defDataAdd(varname = "r", formula = 0, variance= 0.0625, 
                    dist = "normal")    #each site has a random treatment effect: N(0,sd=0.25)

defC3 <- defCondition(condition = "large==1",  
                      formula = "2*(10 + 45*largesite)", 
                      dist = "nonrandom") ##RCTs with size 150: 3 sites with size 20,20 and 110
defC3 <- defCondition(defC3,condition = "large==0",  
                      formula = "(10 + 45*largesite)", 
                      dist = "nonrandom") ##RCTs with size 75: 3 sites with size 10,10 and 55

defC4 <- defDataAdd(varname="C_rv", formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3

defC4 <- defDataAdd(defC4, varname = "z", 
                    formula = "r + (0.3 + b ) * (C_rv==1) + (0.4 + b ) * (C_rv==2) + (0.5 + b ) * (C_rv==3)", 
                    dist = "nonrandom")


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


