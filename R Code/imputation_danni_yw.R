library(data.table)
library(dplyr)
library(mice)
library(devtools)
library(naniar)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(slurmR)


set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
mod <- cmdstan_model("/gpfs/home/dw2625/r/imputation/primary_co.stan")
modlg <- cmdstan_model("/gpfs/home/dw2625/r/imputation/primary_l.stan")

dat <- read.csv('/gpfs/home/dw2625/r/imputation/COMPILE_merged_20210309.csv',na.strings = '-999')
#dat <- read.csv('data/COMPILE_merged_20210309.csv',na.strings = '-999')
#modlg <- cmdstan_model("./Mice/primary_l.stan")
dat$ID <- paste(dat$tr_rct_id,dat$site_id,dat$pt_id,sep='_')

##withdraw or screen failure before day 14
dat <- dat[!(dat$ID %in% c('AA_2_25','AA_2_72','AA_15_3','AA_15_11','AA_10_30','AA_14_88','AA_8_25','AA_9_51','AA_9_8','AA_8_41',
                         'AA_10_19','AA_8_6','BB_1_46','BB_1_47','DD_1_5156')),]
dat$out_who14 <- ifelse(dat$out_who14 %in% 999,-999,dat$out_who14)
## outcome: WHO_14
dat$WHO_14 <- ifelse(dat$tr_rct_id %in% c('DD','EE','RR'),dat$out_who15,ifelse(dat$tr_rct_id == 'BB',dat$out_who13,dat$out_who14))
dat[is.na(dat$WHO_14) & !is.na(dat$out_who14),]$WHO_14 <- dat[is.na(dat$WHO_14) & !is.na(dat$out_who14),]$out_who14 # patients who did not have WHO 14 but had WHO 14
dat[is.na(dat$WHO_14) & !is.na(dat$out_who13),]$WHO_14 <- dat[is.na(dat$WHO_14) & !is.na(dat$out_who13),]$out_who13 # patients who did not have WHO 14 but had WHO 13
dat[is.na(dat$WHO_14) & !is.na(dat$out_who15),]$WHO_14 <- dat[is.na(dat$WHO_14) & !is.na(dat$out_who15),]$out_who15 # patients who did not have WHO 14 but had WHO 15

### nearest WHO from left
dat$who_neighbor_left <- sapply(1:nrow(dat), function(i) {
   who_0_12 <- dat[i,paste0('out_who',0:12)]
   if (all(is.na(who_0_12))) {return(NA)} else {
   return(who_0_12[max(which(!is.na(who_0_12)))])
   }
}) %>% as.numeric(.)

dat$who_neighbor_left_to_day14 <- sapply(1:nrow(dat), function(i) {
  who_0_12 <- dat[i,paste0('out_who',0:12)]
  if (all(is.na(who_0_12))) {return(NA)} else {
    return(14 - max(which(!is.na(who_0_12))))
  }
}) %>% as.numeric(.)

### nearest WHO from right
dat$who_neighbor_right <- sapply(1:nrow(dat), function(i) {
  who_16_90 <- dat[i,paste0('out_who',16:90)]
   if (all(is.na(who_16_90))) {return(NA)} else {
     return(who_16_90[min(which(!is.na(who_16_90)))])
  }
 }) %>% as.numeric(.)

dat$who_neighbor_right_to_day14 <- sapply(1:nrow(dat), function(i) {
  who_16_90 <- dat[i,paste0('out_who',16:90)]
  if (all(is.na(who_16_90))) {return(NA)} else {
    return(min(which(!is.na(who_16_90))) + 1)
  }
}) %>% as.numeric(.)


### nearest 
dat$who_nearest_neighbor_to_day14 <-  sapply(1:nrow(dat), function(i) {
  who_right_day <- dat[i,"who_neighbor_right_to_day14"]
  who_left_day <- dat[i,"who_neighbor_left_to_day14"]
  if(is.na(who_left_day)&is.na(who_right_day)){return(NA)}else{
    return(min(who_left_day,who_right_day,na.rm=TRUE))
  }
})%>% as.numeric(.)


dat$who_nearest_neighbor <- sapply(1:nrow(dat), function(i) {
  who_nearest_day <- dat[i,"who_nearest_neighbor_to_day14"]
  who_right_day <- dat[i,"who_neighbor_right_to_day14"]
  who_left_day <- dat[i,"who_neighbor_left_to_day14"]
  
  if (is.na(who_nearest_day)){
    return(NA)}
  else if (!is.na(who_nearest_day)&(who_nearest_day %in% who_left_day)){
    ## tie or =who_left_day and not= who_right_day
    return(dat[i,"who_neighbor_left"])
  }
  else if(!is.na(who_nearest_day)&(who_nearest_day %in% who_right_day)&(!(who_nearest_day %in% who_left_day))){
    return(dat[i,"who_neighbor_right"])
  }
})%>% as.numeric(.)
  

dat$discharge_by_day7 <- ifelse(!is.na(dat$out_discharge) & dat$out_discharge %in% 0:7, 1, 0)
dat$discharge_by_day15 <- ifelse(!is.na(dat$out_discharge) & dat$out_discharge %in% 0:15, 1, 0)

### baseline covariates
b.cov <- c("tr_rct_id","bl_age","bl_sex","bl_symp_days","bl_enrollqtr",'bl_blood','bl_cardio','bl_pulm','bl_diabetes','out_who0')
cov <- c("ID","tr_rct_id","bl_age","bl_sex","bl_symp_days","bl_enrollqtr",'bl_blood','bl_cardio','bl_pulm','bl_diabetes','out_who0','discharge_by_day7','discharge_by_day15',
         'who_nearest_neighbor_to_day14','who_nearest_neighbor','WHO_14')
dat <- dat %>% mutate_at(vars(c("tr_rct_id","bl_sex","bl_symp_days",'bl_cardio','bl_pulm','bl_diabetes','out_who0','discharge_by_day7','discharge_by_day15')),
                         function(x) factor(x))
dat$bl_blood <- factor(dat$bl_blood,0:3,c('O','A','B','AB'))
dat$bl_enrollqtr <- factor(dat$bl_enrollqtr,2:5,c('Apr-June 2020','July-Sept 2020','Oct-Dec 2020',
                                                  'Jan-Mar 2021'))
dat$who_nearest_neighbor <- factor(dat$who_nearest_neighbor,0:10)
dat$WHO_14 <- factor(dat$WHO_14,0:10)

########################CP group#################################
dind <- dat[dat$bl_treatment ==1,cov]

# show the missing data pattern
#vis_miss(dind)

########## set pred matrix
pred.mat <- matrix(1,nrow = length(cov), ncol = length(cov))
diag(pred.mat) <- 0 
rownames(pred.mat) = colnames(pred.mat) <- cov
# Post-baseline variables will not be used to predict baseline covariates
pred.mat[c('discharge_by_day7','discharge_by_day15',
           'who_nearest_neighbor_to_day14','who_nearest_neighbor','WHO_14'),b.cov] <- 0

##ID has not missing value and will not be used to predict other covariates
pred.mat[("ID"),cov] <- 0
pred.mat[,"ID"] <- 0

data_imp_CP <- mice(
  dind,
  m = 50,
  predictorMatrix=pred.mat
)

#save(data_imp_CP,file="./data/data_imp_CP.rda")
########################Control group#################################
dind <- dat[dat$bl_treatment ==0,cov]

# show the missing data pattern
#vis_miss(dind)

########## set pred matrix
pred.mat <- matrix(1,nrow = length(cov), ncol = length(cov))
diag(pred.mat) <- 0 
rownames(pred.mat) = colnames(pred.mat) <- cov
# Post-baseline variables will not be used to predict baseline covariates
pred.mat[c('discharge_by_day7','discharge_by_day15',
           'who_nearest_neighbor_to_day14','who_nearest_neighbor','WHO_14'),b.cov] <- 0

##ID has not missing value and will not be used to predict other covariates
pred.mat[("ID"),cov] <- 0
pred.mat[,"ID"] <- 0

data_imp_Ctrl <- mice(
  dind,
  m = 50,
  predictorMatrix=pred.mat
)

#save(data_imp_Ctrl,file="./data/data_imp_Ctrl.rda")

dt_to_list <- function(dx) {
  N <- nrow(dx)                               ## number of observations 
  L <- dx[, length(unique(ordY))]             ## number of levels of outcome 
  K <- dx[, length(unique(study))]            ## number of studies 
  y <- as.numeric(dx$ordY)                    ## individual outcome 
  y_2 <- ifelse(dx$ordY %in% seq(8,11),1,0)         ##logistic outcome  
  kk <- dx$study                              ## study for individual 
  ctrl <- dx$control                          ## treatment arm for individual 
  cc <- dx[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = dx)[, -1]
  D <- ncol(x)
  
  prior_div <- 8   #SD of site-specific intercept
  prior_Delta_sd <- .354
  prior_beta_sd <- 2.5
  eta <- 0.1      #sd of delta_c
  prior_eta_0 <- 0.25 
  
  list(N=N, L=L, K=K, y=y, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, x=x, D=D,
       prior_div=prior_div,prior_Delta_sd=prior_Delta_sd,
       prior_beta_sd=prior_beta_sd,eta=eta,prior_eta_0=prior_eta_0)
}

data_extra <- dat[,c("ID","tr_cntrl","bl_treatment")]

Missing_imp <- function(iternum, data_extra,data_imp_CP, data_imp_Ctrl, lg_model,c_model){
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  
  data_CP <- complete(data_imp_CP, action = iternum)
  data_Ctrl <- complete(data_imp_Ctrl, action = iternum)
  data_imp <- rbind(data_CP,data_Ctrl)
  dind <- merge(data_extra,data_imp,by="ID")
  
  ## data for Stan
  dind$ordY <- as.numeric(dind$WHO_14) #ordY:1-11
  dind$who_enroll <- factor(dind$out_who0,4:6,5:7) # since we change the scale from 0-10 to 1-11
  dind$cat_age <- cut(dind$bl_age, c(0,50,65,120))
  dind$age <- factor(dind$cat_age,c("(0,50]","(50,65]","(65,120]"),1:3)
  dind$sex <- dind$bl_sex #0:male; 1:female
  dind$ss <- dind$bl_symp_days
  dind$qtr <- dind$bl_enrollqtr
  dind$study <- factor(dind$tr_rct_id,c("AA","BB","CC","DD","EE","GG","KK","RR"),1:8)##
  dind$study <- as.numeric(dind$study)
  dind$C <- factor(dind$tr_cntrl,c(0,1,2),1:3) #control arm for study; in SAP:0:standard of care;1:Non_CP;2:Saline
  dind$control <- ifelse(dind$bl_treatment==1,0,1)  ## treatment arm for individual (control=1;CP=0)
  dind <- dind%>% as.data.table()
  dx <- dind[, .(ID,study, C, control, sex, who_enroll, age, ss,qtr, ordY)]
  
  
  fit_lg <- lg_model$sample(
    step_size = 0.1,
    data = dt_to_list(dx),
    chains = 4L,
    parallel_chains = 4L,
    refresh = 500,
    iter_warmup = 2500,
    iter_sampling = 2500,
    adapt_delta=0.98,
    max_treedepth =15
    
  )
  
  draws_lg <- data.table(as_draws_df(fit_lg$draws()))
  Delta_lg <- draws_lg[, Delta]
  
  fit_co <- c_model$sample(
    step_size = 0.1,
    data = dt_to_list(dind),
    chains = 4L,
    parallel_chains = 4L,
    refresh = 500,
    iter_warmup = 2500,
    iter_sampling = 2500,
    adapt_delta=0.98,
    max_treedepth =15
  )
  
  draws_co <- data.table(as_draws_df(fit_co$draws()))
  Delta_co <- draws_co[, Delta]
  
  return(list(LG=Delta_lg,CO=Delta_co))
}

job <- Slurm_lapply(
  X = 1:50, 
  FUN = Missing_imp, 
  data_extra =data_extra,
  data_imp_CP = data_imp_CP,
  data_imp_Ctrl = data_imp_Ctrl,
  lg_model = modlg,
  c_model =mod,
  njobs = 50, 
  mc.cores = 4L,
  job_name = "i_lg10",
  tmp_path = "/gpfs/scratch/dw2625",
  plan = "wait",
  sbatch_opt = list(time = "36:00:00", partition = "cpu_median"),
  export = c("dt_to_list"),
  overwrite = TRUE
)

job
res <- Slurm_collect(job)

###Combine results
lg_res <- lapply(res,function(x) return(x[['LG']])) %>% Reduce(c,.)
co_res <- lapply(res,function(x) return(x[['CO']])) %>% Reduce(c,.)
save(result, file = "/gpfs/home/dw2625/r/imputation/result/result_test.rda")

prob_success_lg1 <- mean(exp(-1*lg_res) < 1)
prob_success_lg08 <- mean(exp(-1*lg_res) < 0.8)

prob_success_co1 <- mean(exp(-1*co_res) < 1)
prob_success_co08 <- mean(exp(-1*co_res) < 0.8)


data.table(prob_success_lg1,prob_success_lg08,prob_success_co1,prob_success_co08)

