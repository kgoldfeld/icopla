library(simstudy)
library(data.table)
library(parallel)
library(cmdstanr)
library(posterior)
library(slurmR)

freq_fit <- function(dx) {
  
  lmfit <- lm(y ~ rx, data = dx)
  coef(summary(lmfit))["rx", "Pr(>|t|)"]
  
}

bayes_fit <- function(mod, dx) {
  
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  
  data_list <- list(N=nrow(dx), y=dx$y, rx=dx$rx, p_mu=0, p_sigma=5)
  
  fit <- mod$sample(
    data = data_list,
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 500,
    iter_sampling = 2500,
    step_size = 0.1,
    show_messages = FALSE
  )
  
  df <- data.frame(as_draws_rvars(fit$draws(variables = "beta")))
  
  #cred_int <- quantile(df$beta, probs = c(0.025, 0.975))
  #!(0 %between% cred_int)
  
  ((mean(df$beta > 0) > 0.95) & (mean(df$beta > 0.2 ) > 0.5))
  
}

seq_mods <- function(dx, mod, start, end, by) {
  
  freq_ps <- sapply(seq(start, end, by = by), function(x) freq_fit(dx[1:x]))
  freq_effect <- any(freq_ps < 0.05)
  
  bayes_ci <- sapply(seq(start, end, by = by), function(x) bayes_fit(mod, dx[1:x]))
  bayes_effect <- any(bayes_ci)
  
  return(data.table(freq_effect, bayes_effect))
  
}

s_replicate <- function(count, mod, start, end, by) {
  
  def <- defData(varname = "trt", formula = "1/2;1/2", dist = "categorical")#trt=1 or 2
  def <- defData(def, varname = "rx", formula = "trt-1",  
                     dist = "nonrandom")  
  def <- defData(def, varname = "y", formula = 0, variance = 1, dist = "normal")
  
  dd <- genData(end, def)
  
  seq_mods(dd, mod, start, end, by)
  
}

set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
mod <- cmdstan_model("/gpfs/data/troxellab/danniw/stan/multiple_v2.stan")
#mod <- cmdstan_model("multiple_v2.stan")
# s_replicate(count=1,mod, start = 100, end = 1000, by = 100)
# rbindlist(mclapply(1:2, function(x) s_replicate(mod, start = 1000,end = 1000, by = 0), mc.cores = 4))

print("stating job")

job <- Slurm_lapply(
  X = 1L:2L,
  FUN = s_replicate,
  mod = mod,
  start = 100L,
  end = 1000L,
  by = 100L,
  njobs = 2L,
  mc.cores = 4L,
  job_name = "i_mult_1",
  tmp_path = "/gpfs/data/troxellab/danniw/scratch",
  plan = "wait",
  sbatch_opt = list(time = "12:00:00", partition = "cpu_short", `mem-per-cpu` = "5G"),
  export = c("seq_mods", "freq_fit", "bayes_fit"),
  overwrite = TRUE
)

print("finishing job")

res <- Slurm_collect(job)
save(res, file = "/gpfs/data/troxellab/danniw/data/mult.rda")
