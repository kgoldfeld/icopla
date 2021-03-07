library(simstudy)
library(cmdstanr)
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


defC <- defDataAdd(varname = "a", formula = 0, variance = .005, 
                   dist = "normal")    
defC <- defDataAdd(defC, varname = "b", formula = 0, variance= .01, 
                   dist = "normal")                  #each study has a random slope
defC <- defDataAdd(defC, varname = "size", formula = "150+150*large", 
                   dist = "nonrandom") 
defC2 <- defDataAdd(
  varname="C_rv",
  formula="C * control",
  dist = "nonrandom"
) # C_rv=1/2/3: patient received control treatment C=1/2/3
defC2 <- defDataAdd(defC2, varname = "ss", formula = "0.25;0.5;0.25",
  dist = "categorical") 

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
  formula = "0.05*(ss-1) + z_ss + (0.4 + b) * (C_rv==1) + (0.4 + b) * (C_rv==2) + (0.4 + b) * (C_rv==3)", 
  dist = "nonrandom")

nsites <- 9
basestudy <- genBaseProbs(n = nsites,
                          base =  c(0.050, 0.107, 0.095, 0.085, 0.090, 0.090, 0.108, 0.150, 0.090, 0.075, 0.060), 
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

dind[, table(ss)]
##################ss=1 as the reference level###################

dind$ss <- relevel(as.factor(dind$ss), ref = "1")
x <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]
D <- ncol(x)

sum(x[,1])
sum(x[,2])

ds <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]
n_ds <- ncol(ds)

studydata <- list(
  N=N, L= L, K=K, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds, n_ds=n_ds, D=D, x=x)


##Fit the model

sm <- cmdstan_model("KSG/interaction_test.stan")

fit_sm <- sm$sample(
  step_size = 0.1,
  data = studydata,
  chains = 4L,
  parallel_chains = 4L,
  refresh = 500,
  iter_warmup = 500,
  iter_sampling = 2000
)

##################ss=3 as the reference level###################

dind$ss <- relevel(as.factor(dind$ss), ref = "3")

x <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]
D <- ncol(x)

sum(x[,1])
sum(x[,2])

ds <- model.matrix(ordY ~ factor(ss), data = dind)[, -1]
n_ds <- ncol(ds)

studydata <- list(
  N=N, L= L, K=K, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds, n_ds=n_ds, D=D, x=x)

sm_relevel <- cmdstan_model("KSG/interaction_test_relevel.stan")

fit_sm_relevel <- sm_relevel$sample(
  step_size = 0.1,
  data = studydata,
  chains = 4L,
  parallel_chains = 4L,
  refresh = 500,
  iter_warmup = 500,
  iter_sampling = 2000
)

#---

dind$ss <- factor(dind$ss, levels = c(1, 2, 3))
S <- length(unique(dind$ss))

x <- model.matrix(ordY ~ factor(ss) - 1, data = dind)
D <- ncol(x)

ds <- model.matrix(ordY ~ factor(ss) - 1, data = dind)
n_ds <- ncol(ds)

studydata <- list(
  N=N, L= L, K=K, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ss=dind$ss, D=D, x=x, S=S)

sm_sub <- cmdstan_model("KSG/interaction_test_sub.stan")

fit_sm_sub <- sm_sub$sample(
  step_size = 0.1,
  data = studydata,
  chains = 4L,
  parallel_chains = 4L,
  refresh = 500,
  iter_warmup = 500,
  iter_sampling = 2000
)



#---

fit_sm$summary(variables = c("Delta", "Gamma[1]", "Gamma[2]","GD[1]", "GD[2]", "GD[3]"))
fit_sm_relevel$summary(variables = c("Delta", "Gamma[1]", "Gamma[2]","GD[1]", "GD[2]", "GD[3]"))
fit_sm_sub$summary(variables = c("Gamma[1]", "Gamma[2]","Gamma[3]"))

draws_array <- as_draws_array(fit_sm$draws())
dcov <- data.table(as.vector(draws_array[,,"Delta"]), as.vector(draws_array[,,"Gamma[1]"]))
ggplot(data= dcov, aes(x = V1, y = V2)) +
  geom_point()
dcov <- data.table(as.vector(draws_array[,,"Delta"]), as.vector(draws_array[,,"Gamma[2]"]))
ggplot(data= dcov, aes(x = V1, y = V2)) +
  geom_point()

draws_array <- as_draws_array(fit_sm_relevel$draws())
dcov <- data.table(as.vector(draws_array[,,"Delta"]), as.vector(draws_array[,,"Gamma[1]"]))
ggplot(data= dcov, aes(x = V1, y = V2)) +
  geom_point()
dcov <- data.table(as.vector(draws_array[,,"Delta"]), as.vector(draws_array[,,"Gamma[2]"]))
ggplot(data= dcov, aes(x = V1, y = V2)) +
  geom_point()
