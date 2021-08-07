#############proportional
load("./OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/by_trt/predictive_iter9900_v2.rda")
#load("./data/predictive_iter2000.rda")
who <- unlist2d(site_plasma,idcols="replicate")
who <- as.data.table(who)

who[, cumprob := cumsum(N)/sum(N) , keyby = .(n_draw,control)]

load("./OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/by_trt/dind_local_v2.rda")
#load("./data/dind_local_v2.rda")

dsum <- dind[, .(N = sum(.N)), keyby = .(control,ordY)]

dsum[, cumprob := cumsum(N)/sum(N) , keyby = control]

round(dsum$cumprob*100,2)
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

######Non-proportional###########
load("./OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/by_trt/np_predictive_iter9900_v4.rda")
#load("./data/predictive_iter2000.rda")
who <- unlist2d(site_plasma,idcols="replicate")
who <- as.data.table(who)

who[, cumprob := cumsum(N)/sum(N) , keyby = .(n_draw,control)]

load("./OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/by_trt/dind_np_local_v3.rda")
#load("./data/dind_local_v2.rda")

dsum <- dind[, .(N = sum(.N)), keyby = .(control,ordY)]

dsum[, cumprob := cumsum(N)/sum(N) , keyby = control]
round(dsum$cumprob*100,2)

lapply(1:11, function(x) pvalue(dx=who,do=dsum,x))
lapply(1:11, function(x) CI(dx=who,x))
