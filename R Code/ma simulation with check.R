library("simstudy")
library("ordinal")


#### Data definitions

defC <- defDataAdd(varname = "a", formula = 0, variance = 2, 
                dist = "normal")    #each study has a random intercept
defC <- defDataAdd(defC, varname = "b", formula = 0, variance= .5, 
                dist = "normal")                  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "100000+0*large", 
                dist = "nonrandom") 
defC2 <- defDataAdd(varname="C_rv",formula="C*rx",
                dist = "nonrandom") # C_rv=1/2/3: patient received control treatment C=1/2/3
defA1 <- defDataAdd(varname = "z", 
  formula = "a + (0.5 + b ) * (C_rv==1) + (1 + b ) * (C_rv==2) + (2 + b ) * (C_rv==3)", 
  dist = "nonrandom")

baseprobs <- c(0.03, 0.03, 0.06, 0.09, 0.09, .15, .15, .13, .10, .10, .07) 

#### Generate data

set.seed(18459)

dstudy <- genData(9, id = "study")       # generate studies
dstudy <- trtAssign(dstudy, nTrt = 3, grpName = "C") # allocate to control group
dstudy <- trtAssign(dstudy, nTrt = 2, strata = "C", grpName = "large", ratio = c(2,1))
dstudy <- addColumns(defC, dstudy)

dind <- genCluster(dstudy, "study", numIndsVar = "size", "id")
dind <- trtAssign(dind, strata="study", grpName = "rx") #add control treatment variable, randomly assign control treatment in each study;rx=1:received control;rx=0:received convalescent plasma
dind <- addColumns(defC2,dind)
dind <- addColumns(defA1, dind) # add z
dind <- genOrdCat(dind, adjVar = "z", baseprobs, catVar = "ordY") #add outcome's level


### Check proportions

getProp <- function(dx) {
  dx[, round(prop.table(table(C_rv, ordY), margin = 1), 2)]
}

getProp(dind)
lapply(1:9, function(x) getProp(dind[study == x]))

# Check study 4

getProp(dind[study==4])

1 / (1 + exp(-(log(.03/.97) - (-.9732))))
1 / (1 + exp(-(log(.03/.97) - (-.9732 + ( -1 + 0.1843612)) )))

1 / (1 + exp(-(log(.06/.94) - (-.9732))))
1 / (1 + exp(-(log(.06/.94) - (-.9732 + ( -1 + 0.1843612)) )))

1 / (1 + exp(-(log(.12/.88) - (-.9732))))
1 / (1 + exp(-(log(.03/.97) - (-.9732 + ( -1 + 0.1843612)) )))



