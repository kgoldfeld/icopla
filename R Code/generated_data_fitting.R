library(cmdstanr)
library(posterior)
library(bayesplot)
library(glue)
library(simstudy)
library(data.table)
library(slurmR)
library(dplyr)

set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
mod <- cmdstan_model("/gpfs/home/dw2625/r/predco.stan")
modlg <- cmdstan_model("/gpfs/home/dw2625/r/predlogistic.stan")

#set_cmdstan_path(path = "./cmdstan-2.25.0")
#mod <- cmdstan_model("./predco.stan")
#modlg <- cmdstan_model("./predlogistic.stan")

dind <- read.csv('./generated_data.csv',na.strings = '-999')
dind <- dind%>%data.table()
dind$ordY <- dind$WHO_14_imp + 1  #ordY:1-11
dind$who_enroll <- dind$out_who0 + 1 # since we change the scale from 0-10 to 1-11
dind$cat_age <- cut(dind$bl_age, c(0,50,65,120))
dind$age <- factor(dind$cat_age,c("(0,50]","(50,65]","(65,120]"),1:3)
dind$sex <- dind$bl_sex #0:male; 1:female
dind$ss <- recode(dind$bl_symp_days,'1'=1,'2'=2,'3'=3,'4'=3,'5'=4)
dind$qtr <- dind$bl_enrollqtr
dind$study <- factor(dind$tr_rct_id,c("AA","BB","CC","DD","EE","GG","KK","RR"),1:8)##
dind$study <- as.numeric(dind$study)
dind$C <- factor(dind$tr_cntrl,c(0,1,2),1:3) #control arm for study; 
dind$control <- ifelse(dind$bl_treatment==1,0,1)  ## treatment arm for individual (control=1;CP=0)

dind <- dind[, .(id,study, C, control, sex, who_enroll, age, ss,qtr, ordY)]

#--- Create a data set with observed outcomes ---#
dobs <- dind[!is.na(ordY)]
#dobs[, .N, keyby = .(ordY)]
#sum(dind[, .N, keyby = .(id)]$N ==1)

#--- Fit Bayes model to original data ---#

dt_to_list <- function(dx) {
  
  N <- nrow(dx)                               ## number of observations 
  L <- dx[, length(unique(ordY))]             ## number of levels of outcome 
  K <- dx[, length(unique(study))]            ## number of studies 
  y <- as.numeric(dx$ordY)                    ## individual outcome 
  y_2 <- ifelse(dx$ordY %in% seq(8,11),1,0)         ##logistic outcome  
  kk <- as.numeric(dx$study)                              ## study for individual 
  ctrl <- dx$control                          ## treatment arm for individual 
  cc <- dx[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
  x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = dx)[, -1]
  D <- ncol(x)
  
  list(N=N, L=L, K=K, y=y, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, x=x, D=D)
}

fit <- mod$sample(
  step_size = 0.1,
  data = dt_to_list(dobs),
  chains = 4L,
  parallel_chains = 4L,
  refresh = 500,
  iter_warmup = 2500,
  iter_sampling = 2500
)

#--- Save posterior distribution ---#

draws_dt <- data.table(as_draws_df(fit$draws()))
mean(draws_dt[, exp(-Delta)] < 1)

diagnostics_df <- as_draws_df(fit$sampler_diagnostics())
div_num <- sum(diagnostics_df[, 'divergent__'])
div_num


#--- use original with missing outcome - represents "new" data ---#

newdt <- dind[is.na(ordY)]
newdt[, calt := control * study]
newdt[, ordY := NULL]
setkey(newdt, "id")
#--- Function to for generating outcomes and estimating model for new data ---#

est_from_draw <- function(n_draw, Draws, dt_obs, dt_new, D, K, s_model,c_model) {
  
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  
  draw <- as.data.frame(Draws[sample(.N, 1)])
  
  alpha <- draw$alpha
  beta <- as.vector(x = draw[, glue("beta[{1:D}]")], mode = "numeric")
  delta_k <- as.vector(draw[, glue("delta_k[{1:K}]")], mode = "numeric")
  
  tau <- as.vector(draw[grep("^tau", colnames(draw))], mode = "numeric")
  tau <- matrix(tau, nrow = K)
  tau <- cbind(tau, Inf)
  cprop <- t(apply(tau, 1, stats::plogis))
  xprop <- t(apply(cprop, 1, diff))
  baseline <- cbind(cprop[,1], xprop) 
  #dt_new[, .N, keyby = .(study, control)]
  coefs <- c(alpha, beta, delta_k[c(1,2,5,7)])#only RCT 1,2,5,7 has missing outcome in the control group
  x <- model.matrix(~factor(who_enroll) + factor(age) + 
                      factor(sex) + factor(ss) +factor(qtr), data = dt_new)
  calt <- model.matrix(~factor(calt) - 1, data = dt_new)[,-1]
  zmat <- cbind(x, calt)
  dt_new$z <- zmat %*% coefs
  dt_new[, calt := NULL]
  
  ##Obviously, studies that have completed data collection 
  #and already provided us with full data sets will not have any new simulated data in this analysis.
  #do not need to change: if a study is complete, dt_new won't including it:dt_new[study == i]
  dl <- lapply(c(1,2,5,7), function(i) {
    b <- baseline[i,]
    dx <- dt_new[study == i]
    genOrdCat(dx, adjVar = "z", b, catVar = "ordY")
  })
  
  dt_new <- rbindlist(dl)  
  dt_new$ordY
  #dt_new[, .N, keyby = .(ordY)]
  dt_new$ordY  <- dt_new$ordY %>% as.character(.) %>% as.numeric(.)
  dt_new <- dt_new[, .(id,study, C, control, sex, who_enroll, age, ss, qtr, ordY)]
  
  dx <- rbind(dt_obs, dt_new)
  
  # fit logistic model for combined data set of "complete" data and simulated new data
  
  fit_lg <- s_model$sample(
    step_size = 0.1,
    data = dt_to_list(dx),
    chains = 4L,
    parallel_chains = 4L,
    refresh = 500,
    iter_warmup = 2000,
    iter_sampling = 2500
  )
  
  draws_lg <- data.table(as_draws_df(fit_lg$draws()))
  
  prob_success_lg1 <- mean(draws_lg[, exp(-Delta)] < 1)
  prob_success_lg08 <- mean(draws_lg[, exp(-Delta)] < 0.8)
  
  diagnostics_df <- as_draws_df(fit_lg$sampler_diagnostics())
  div_num_lg <- sum(diagnostics_df[, 'divergent__'])
  
  
  # fit CO model for combined data set of "complete" data and simulated new data
  
  fit_pp <- c_model$sample(
    step_size = 0.1,
    data = dt_to_list(dx),
    chains = 4L,
    parallel_chains = 4L,
    refresh = 500,
    iter_warmup = 2000,
    iter_sampling = 2500,
    adapt_delta=0.8
  )
  
  draws_pp <- data.table(as_draws_df(fit_pp$draws()))
  
  prob_success_pp1 <- mean(draws_pp[, exp(-Delta)] < 1)
  prob_success_pp08 <- mean(draws_pp[, exp(-Delta)] < 0.8)
  
  diagnostics_df <- as_draws_df(fit_pp$sampler_diagnostics())
  div_num_pp <- sum(diagnostics_df[, 'divergent__'])
  
  return(data.table(n_draw, prob_success_lg1,prob_success_lg08,div_num_lg,prob_success_pp1,prob_success_pp08,div_num_pp))
  
  
  # dtNuts <- data.table(nuts_params(fit_pp))
  # pct_div <- dtNuts[Parameter == "divergent__", mean(Value)]
  
  # return(data.table(prob_success, pct_div))
}

#--- Multiple draws from the posterior ---#

job <- Slurm_lapply(
  X = 1L:1080L, 
  FUN = est_from_draw, 
  Draws = draws_dt,
  dt_obs = dobs,
  dt_new = newdt,
  D = dt_to_list(dobs)$D,
  K = dt_to_list(dobs)$K,
  s_model = modlg,
  c_model =mod,
  njobs = 90L, 
  mc.cores = 4L,
  job_name = "i_lg",
  tmp_path = "/gpfs/scratch/dw2625",
  plan = "wait",
  sbatch_opt = list(time = "03:00:00", partition = "cpu_short"),
  export = c("dt_to_list"),
  overwrite = TRUE
)

job
res <- Slurm_collect(job)

save(res, file = "/gpfs/home/dw2625/r/pred_co_l.rda")
