library(collapse)
library(data.table)
library(dplyr)
library(ggplot2)
load("C:/Users/Danni/OneDrive - NYU Langone Health/R code - COMPILE/COMPILE_final/posterior_predictive_checks/data/predictive_compile.rda")
who <- unlist2d(site_plasma,idcols="replicate")
who <- as.data.table(who)
who[, cumprob := cumsum(N)/sum(N) , keyby = .(n_draw,control)]

##########95%CI######
CI <- function(dx,x){
  rep_data <- subset(dx, (control %in% 0)&(ordY %in% x))
  rep <- rep_data$cumprob
  CI_0 <- round(quantile(rep,probs=c(.025,.975))*100,2)
  
  rep_data <- subset(dx, (control %in% 1)&(ordY %in% x))
  rep <- rep_data$cumprob
  CI_1 <- round(quantile(rep,probs=c(.025,.975))*100,2)
  data.table(paste0("[",CI_0[1],", ",CI_0[2],"]"),  paste0("[",CI_1[1],", ",CI_1[2],"]"))
}

lapply(1:11, function(x) CI(dx=who,x))

#################observed data
load("C:/Users/Danni/OneDrive - NYU Langone Health/R code - COMPILE/COMPILE_final/posterior_predictive_checks/data/dind_v2.rda")
dind$control <- ifelse(dind$bl_treatment==1,0,1)  ## treatment arm for individual (control=1;CP=0)
dind$study <- as.numeric(dind$study)

dsum <- dind[, .(N = sum(.N)), keyby = .(control,ordY)]

dsum[, cumprob := cumsum(N)/sum(N) , keyby = control]
##Test statistics in observed data
dsum[, percentage := round(cumprob*100,2)]

#######p-value##########
pvalue <- function(dx,do,x){
  rep_data <- subset(dx, (control %in% 0)&(ordY %in% x))
  rep <- rep_data$cumprob
  orig_data <- subset(do, (control %in% 0)&(ordY %in% x))
  orig <- orig_data$cumprob
  p_0 <- round(mean(rep > orig),2)
  
  rep_data <- subset(dx, (control %in% 1)&(ordY %in% x))
  rep <- rep_data$cumprob
  orig_data <- subset(do, (control %in% 1)&(ordY %in% x))
  orig <- orig_data$cumprob
  p_1 <- round(mean(rep > orig),2)
  
  data.table(p_0,p_1)
}

lapply(1:11, function(x) pvalue(dx=who,do=dsum,x))

