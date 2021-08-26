library(posterior) 
library(data.table)
library(cmdstanr) 

dat <- read.csv('data/COMPILE_merged_20210525.csv',na.strings = '-999')
dat$ID <- paste(dat$tr_rct_id,dat$site_id,dat$pt_id,sep='_')
dat$out_who0[is.na(dat$out_who0)] <- dat$out_who1[is.na(dat$out_who0)]
## outcome: WHO_14
dat$WHO_14 <- ifelse(dat$tr_rct_id %in% c('DD','EE','RR'),dat$out_who15,ifelse(dat$tr_rct_id == 'BB',dat$out_who13,dat$out_who14))
dat[is.na(dat$WHO_14) & !is.na(dat$out_who14),]$WHO_14 <- dat[is.na(dat$WHO_14) & !is.na(dat$out_who14),]$out_who14 # patients who did not have WHO 14 but had WHO 14
dat[is.na(dat$WHO_14) & !is.na(dat$out_who13),]$WHO_14 <- dat[is.na(dat$WHO_14) & !is.na(dat$out_who13),]$out_who13 # patients who did not have WHO 14 but had WHO 13
dat[is.na(dat$WHO_14) & !is.na(dat$out_who15),]$WHO_14 <- dat[is.na(dat$WHO_14) & !is.na(dat$out_who15),]$out_who15 # patients who did not have WHO 14 but had WHO 15
#da

for (i in dat[dat$tr_rct_id == 'EE',]$ID) {
  tmp <- dat[dat$ID == i,]
  day_discharge <- tmp$out_discharge
  if (!is.na(day_discharge) & day_discharge > 0 & day_discharge < 30) {
    if (!is.na(tmp[,paste0('out_who',day_discharge + 1)])) {
      dat[dat$ID == i,paste0('out_who',day_discharge)] <- tmp[,paste0('out_who',day_discharge + 1)]
    }
  }
}

# replace missing WHO scores by the scores at the time of discharge
dat$WHO_14_imp <- sapply(1:nrow(dat), function(i) {
  day_discharge <- dat[i,]$out_discharge
  outcome <- dat[i,]$WHO_14
  outcome_day <- ifelse(dat[i,]$tr_rct_id %in% c('DD','EE','RR'),15,ifelse(dat[i,]$tr_rct_id == 'BB',13,14))
  if (!is.na(outcome))     
  {return(outcome)} else 
  {if (!is.na(day_discharge) & day_discharge < outcome_day & day_discharge >= 0) {
    return (dat[i,paste0('out_who',day_discharge)])
  } else {
    return(outcome)
  } 
  }
}) %>% as.numeric(.)
dat$WHO_14_imp[dat$ID %in% c('AA_2_25','AA_2_72','AA_15_16','AA_15_3','AA_10_33','AA_10_16','AA_14_18','AA_9_7','AA_8_25','AA_9_51','AA_8_41','AA_10_15','AA_8_6','AA_21_3','AA_21_1','BB_1_46','BB_1_47','DD_1_5156','CC_1_65',"CC_5_9","CC_8_3",'CC_29_2','EE_10_18')] <- NA

dat$WHO_14_imp[dat$out_deathday %in% 0:15] <- 10

dind <- na.omit(dat[,c("ID","WHO_14_imp","out_who0","bl_age","bl_sex","bl_symp_days","tr_rct_id","bl_treatment","tr_cntrl","bl_enrollqtr")])
dind <- dind%>% as.data.table()
dind$ordY <- dind$WHO_14_imp + 1  #ordY:1-11
dind$who_enroll <- dind$out_who0 + 1 # since we change the scale from 0-10 to 1-11
dind$cat_age <- cut(dind$bl_age, c(0,50,65,120))
dind$age <- factor(dind$cat_age,c("(0,50]","(50,65]","(65,120]"),1:3)
dind$sex <- dind$bl_sex #0:male; 1:female
dind$ss <- recode(dind$bl_symp_days,'1'=1,'2'=2,'3'=3,'4'=4,'5'=5)
dind$qtr <- dind$bl_enrollqtr
dind$study <- factor(dind$tr_rct_id,c("AA","BB","CC","DD","EE","GG","KK","RR"),1:8)##
dind$C <- factor(dind$tr_cntrl,c(0,1,2),1:3) #control arm for study; in SAP:0:standard of care;1:Non_CP;2:Saline

### Data for stan model
N = nrow(dind)                                ## number of observations
L <- length(unique(dind$ordY))             ## number of levels of outcome
K <-  length(unique(dind$study))            ## number of studies
y <- as.numeric(dind$ordY)                    ## individual outcome for PO model
#y_2 <- ifelse(y_1 %in% seq(8,11),1,0)         ##WHO 11 points scale 7-10 
kk <- as.numeric(dind$study)                              ## study for individual
ctrl <- ifelse(dind$bl_treatment==1,0,1)  ## treatment arm for individual (control=1;CP=0)
cc <- as.numeric(dind[, .N, keyby = .(study, C)]$C)       ## specific control arm for study
x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss)+factor(qtr), data = dind)[, -1]
D <- ncol(x)

#prior
prior_div <- 8   #SD of site-specific intercept
prior_Delta_sd <- .354
prior_beta_sd <- 2.5
eta <- 0.1      #sd of delta_c
prior_eta_0 <- 0.25 

studydata <- list(
  N=N, L= L,K=K, y=y, ctrl=ctrl,cc=cc,kk=kk, x=x,D=D,prior_div=prior_div,
  prior_Delta_sd=prior_Delta_sd, prior_beta_sd= prior_beta_sd,eta=eta,
  prior_eta_0 = prior_eta_0)

mod <- cmdstan_model("./Stan/primary_co.stan")

fit_co <- mod$sample( 
  data = studydata, 
  refresh = 0, 
  chains = 4L, 
  parallel_chains = 4L, 
  iter_warmup = 2500, 
  iter_sampling = 2500, 
  step_size = 0.1, 
  show_messages = FALSE, 
  adapt_delta=0.99, 
  max_treedepth = 15 
) 

save(fit_co,file = "./posterior_predictive_checks/data/ABCDEGKR_1_fit_v2.rda")


##########draw from posterior iteration################# 
draws_dt <- data.table(as_draws_df(fit_co$draws())) 
save(draws_dt ,file = "./posterior_predictive_checks/data/draws_dt_v2.rda")
save(dind,file = "./posterior_predictive_checks/data/dind_v2.rda")
