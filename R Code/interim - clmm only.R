library("simstudy")
library("rstan")
library("data.table")
library("slurmR")
library("ordinal")
library("dplyr")
library(parallel)

# Function to generate base probabilities for each study

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

iter <- function(iternum, defC, defC2, defA1, basestudy, nsites)  {
  
    # size of interim 
  
    prepData <- function(n, dxfull) {
      
      #choose the first n patients in each study
      
      data_small <- dxfull%>% filter(size==75) %>% group_by(study) %>% do(head(.,n)) %>% as.data.table()
      data_large <- dxfull%>% filter(size==150) %>% group_by(study) %>% do(head(.,2*n)) %>% as.data.table()
      rbind(data_small,data_large)
      
    }
    
    dofit <- function(iternum,n,dxfull){
      
      sampdat <- prepData(n,dxfull)
      
      if(n %in% c(25,50,75)) {
        
        # Frequentist method: 5 looks
        
        sampdat$ordY <- as.factor(sampdat$ordY)
        sampdat$control <- as.factor(sampdat$control)
        sampdat$study <- as.factor(sampdat$study)
        sampdat$C <- as.factor(sampdat$C)
        
        # Frequentist PO model fitting
        
        fit_f <- clmm(formula = ordY ~ 1 + control + (1 + control | C/study),
                      data=sampdat,link = "logit")#random int
        
        freq <- t(summary(fit_f)$coefficients[length(unique(sampdat$ordY)),])
        data.table(iternum,n, freq)

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
    
    #interim look: 16%;33%;50%;60%;67%;80%;90%;100%
    
   # rbindlist(lapply(c(12,25,38,45,50,60,68,75),
    
   rbindlist(lapply(c(25,50,75),
                    function(x) dofit(iternum, x, dxfull = dind)))
}




#### Data definitions

defC <- defDataAdd(varname = "b", formula = 0, variance= 0.005, 
                   dist = "normal")    #each study has a random treatment effect
defC <- defDataAdd(defC, varname = "size", formula = "75+75*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C * control",
                    dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
                    formula = "(0 + b ) * (C_rv==1) + (0 + b ) * (C_rv==2) + (0 + b ) * (C_rv==3)", 
                    dist = "nonrandom")

set.seed(184916)

nsites <- 12

basestudy <- genBaseProbs(
  n = nsites,
#  base =  c(0.100, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.100, 0.090, 0.075, 0.060),
  base =  c(0.2, 0.3, 0.4, 0.1),
  similarity = 50
)

res6 <- mclapply(
  1:20, 
  function(x) {
    iter(x, defC, defC2, defA1, basestudy, nsites)
  },
  mc.cores = 4
)

