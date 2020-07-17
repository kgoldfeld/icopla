library(simstudy)
library(data.table)

### Function to generate base probabilities for each study

genBaseProbs <- function(n, base, similarity = 5, digits = 2) {
  
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

defC <- defDataAdd(varname = "b", formula = 0, variance= .5, 
                   dist = "normal")                  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                   dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                   formula = "(0.5 + b ) * (C_rv==1) + (1 + b ) * (C_rv==2) + (2 + b ) * (C_rv==3)", 
                   dist = "nonrandom")

#### Generate data

set.seed(18459)

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
    base = c(0.03, 0.03, 0.06, 0.09, 0.09, .15, .15, .13, .10, .10, .07),
    similarity = 100)

dl <- lapply(1:nsites, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})

dind <- rbindlist(dl)

