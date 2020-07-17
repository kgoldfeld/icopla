library(simstudy)
library(data.table)

### Function to generate base probabilities for each study

genBaseProbs <- function(n, base, precision = 1, digits = 2) {
  
  n_levels <- length(base)
  x <- gtools::rdirichlet(n, precision * base)

  ### to ensure that each vector sums exactly to 1
  
  x <- round(floor(x*1e8)/1e8, digits)
  xpart <- x[, 1:(n_levels-1)]
  partsum <- apply(xpart, 1, sum)
  x[, n_levels] <- 1 - partsum
   
  return(x)
}

set.seed(18459)

dstudy <- genData(9, id = "study")       # generate studies
dind <- genCluster(dstudy, "study", numIndsVar = 50000, "id")
dind[, z:= 0 ]

#### Old appraoch - same cut-points across studies (without random effect)

baseprobs <- c(.2, .3, .3, .1, .1)

dind1 <- genOrdCat(dind, adjVar = "z", baseprobs, catVar = "ordY") #add outcome's level
dind1[, round(prop.table(table(study, ordY), margin = 1), 2)]

### New approach - Study specific cut-points

basestudy <- genBaseProbs(n = 9, base = c(.5, 2, 3, 3, 1))
basestudy

dl <- lapply(1:9, function(i) {
  b <- basestudy[i,]
  dx <- dind[study == i]
  genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
})

dind2 <- rbindlist(dl)
dind2[, round(prop.table(table(study, ordY), margin = 1), 2)]



    