#-----------------------------------#
#  Plasma Bayes Regression          #
#                                   #
#  author: KSG                      #
#  last modified: 06/30/2020        #
#                                   #
#-----------------------------------#
###########
############ Simulation

### Required libraries
library(simstudy)
library(data.table)
library(ggplot2)
library(rstan)
library(ordinal)

### Simulate data

### Define the data

baseprobs <- c(0.03, 0.03, 0.06, 0.09, 0.09, .15, .15, .13, .10, .10, .07) # need to update based on actual data

# Site level effects

defC <- defData(varname = "a", formula = 0, variance = 2, 
                dist = "normal", id = "site")
defC <- defData(defC, varname = "b", formula = 0, variance= .5, 
                dist = "normal")

# Individual level data

defA1 <- defDataAdd(varname = "risk", formula = .3, dist = "binary")
defA2 <- defDataAdd(varname = "z", 
                    formula = "a + (-1.0 + b ) * rx + 1.5 * risk", 
                    dist = "nonrandom")

### Stan model - should be outside the loop if stan code is not changing
# compiling is a very slow process

rt <- stanc("Programs/plasma_yw.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)


### Generate the data - 7 sites, 100 patients per site (can be variable)

sim.results <- vector(mode = 'list', length = 1000)

### simulation loop

for (i in 1:1000) {

set.seed(18271 + i)

dsite <- genData(7, defC) # generate sites and site effects

dind <- genCluster(dsite, "site", numIndsVar = 100, "id") 
dind <- addColumns(defA1, dind) # add risk variable
dind <- trtAssign(dind, strata=c("site", "risk") ,grpName = "rx") 
dind <- addColumns(defA2, dind) # add z

dind <- genOrdCat(dind, adjVar = "z", baseprobs, catVar = "ordY")

### Data for Stan model

N <- dind[, .N]
K <- length(dind[, levels(ordY)])
J <- length(unique(dind[, site]))
y <- as.vector(dind[, as.numeric(ordY)])
rx <- as.vector(dind[, rx])
jj <- as.vector(dind[, site])
ss <- as.vector(dind[, risk])
#x <- model.matrix(who_final ~ factor(who_enroll) + age + sex, data = dind)[, -1]
#D <- ncol(x)


testdat <- list(N = N, K = K, J = J, #D = D,
    y = y, jj = jj, rx = rx, ss = ss) #, x = NULL)

### Fit data - for actual estimation, should use iter=10000/warmup = 1000

fit <-  sampling(sm, data=testdat, seed = 3327, 
    iter = 3000, warmup = 500, show_messages = FALSE, cores = 4)

### Inspect main parameter estimates

#pname <- c("alpha", "gamma", "delta", "beta", "tau", "OR") # no beta if no other covariates
pname <- c("alpha", "gamma", "delta", "tau", "OR")
alpha.post <- extract(fit, pars = "alpha")[[1]]
gamma.post <- extract(fit, pars = 'gamma')[[1]]
delta.post <- extract(fit, pars = 'delta')[[1]]
tau.post <- extract(fit, pars = "tau")[[1]]
or.post <- extract(fit, pars = 'OR')[[1]]

### Plot diagnostics

# probably we do not this in simulation right? - definitely don't do in simulation
#stan_trace(object = fit, pars = pname) 

sim.results[[i]] <- c(
  or_l_1 = mean(extract(fit, pars = "OR")[[1]] < 1),
  or_l_0.8 = mean(extract(fit, pars = "OR")[[1]] < 0.8),
  or.mean = mean(or.post),
  or.median = median(or.post)
  #some other parameters of interest,
  #....
) # need to confirm what to export here

write.table(dind,paste0('sim_data',i,'.csv'),row.names = FALSE, sep = ',')

} # finish loop

# save sim.results
sim.results <- reduce(sim.results,rbind)
write.table(sim.results,'sim_results_priorX.csv',row.names = FALSE,sep = ',')

### Look at posterior density of OR

#plot_dens <- function(fit, pars, p = c(0.05, 0.95), 
#                      fill = "grey80", xlab = NULL) {
#  
#  qs <- quantile(extract(fit, pars = pars)[[1]], probs = p)
#  
#  x.dens <- density(extract(fit, pars = pars)[[1]])
#  df.dens <- data.frame(x = x.dens$x, y = x.dens$y)
#  
#  p <- stan_dens(fit, pars = c(pars), fill = fill, alpha = .1) +
#    geom_area(data = subset(df.dens, x >= qs[1] & x <= qs[2]), 
#              aes(x=x,y=y), fill = fill, alpha = .4)
#  
#  if (is.null(xlab)) return(p)
#  else return(p + xlab(xlab))
#}
#
#plot_dens(fit, "OR", fill = "#a1be97")




