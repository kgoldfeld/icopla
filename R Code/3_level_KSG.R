library(parallel)
library(lme4)
library(cmdstanr)
library(posterior)
library(bayesplot)

s_define <- function() {
  dc <- defData(varname = "b0_rct", formula = 0, variance = 4, id = "rct")
  dc <- defData(dc, varname = "b1_rct", formula = 0, variance = 2)
  dc <- defData(dc, varname = "nsites", formula = 4, dist = "noZeroPoisson")
  
  ds_add <- defDataAdd(varname = "b0_site", formula = 0, variance = 1)
  ds_add <- defDataAdd(ds_add, varname = "b1_site", formula = 0, variance = 1)
  ds_add <- defDataAdd(ds_add, varname = "n", formula = 30, variance = 1,
                       dist = "negBinomial")
  ds_add <- defDataAdd(ds_add, varname = "npats", formula = "pmax(n, 5)",
                       dist = "nonrandom")
  
  di_add <- defDataAdd(varname = "y", 
                       formula = "b0_rct + b0_site + (3 + b1_rct + b1_site) * rx", variance = 3)
  
  list(dc = dc, ds_add = ds_add, di_add = di_add)
}

s_generate <- function(deflist) {
  
  list2env(x = deflist, envir = environment())
  
  d_rct <- genData(20, dc, id = "rct")
  d_site <- genCluster(d_rct, "rct", numIndsVar = "nsites", "site")
  d_site <- addColumns(ds_add, d_site)
  dd <- genCluster(d_site, "site", numIndsVar = "npats", "id")
  dd <- trtAssign(dd, 2, strata = "site", grpName = "rx") 
  dd <- addColumns(di_add, dd)
  dd[]
  
}

s_estimate <- function(s_model, dd) {
  
  N <- nrow(dd)
  R <- dd[, length(unique(rct))]
  S <- dd[, length(unique(site))]
  site <- dd[, site]
  rct <- dd[, rct]
  rx <- dd[, rx]
  y <- dd[, y]
  site_rct <- dd[, .N, keyby = .(site, rct)][, rct]
  
  fit <- mod1$sample(
    data = list(N=N,R=R,S=S,site=site,rct=rct,rx=rx,y=y,site_rct=site_rct),
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    iter_warmup = 500,
    iter_sampling = 1000,
    step_size = 0.1,
    adapt_delta = 0.8
  )
  
  return(fit)
}

s_model <- cmdstan_model("Stan Code/3_level_KSG.stan")

sim_def <- s_define()
dd <- s_generate(sim_def)
res <- s_estimate(s_model, dd)

mcmc_trace(res$draws(variables = c("b1")))
res$summary(c("b0","b1","sigma", "sigma_r0", "sigma_s0","sigma_r1", "sigma_s1" ))

