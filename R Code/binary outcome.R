library(simstudy)
library(data.table)
library(ggplot2)
library(rstan)
library(bayesplot)

def_c <- defData(varname = "delta_c", formula = .5, variance = .05, id = "c_type")

def_s <- defDataAdd(varname = "b0", formula = 0, variance = .1)
def_s <- defDataAdd(def_s, varname = "delta_k", formula = "delta_c", variance = .05)

def_i <- defDataAdd(
  varname = "y", 
  formula = "-1 + b0 + rx * delta_k", 
  dist = "binary", 
  link = "logit")

dc <- genData(3, def_c)

ds <- genCluster(dc, "c_type", numIndsVar = 3, level1ID = "site")
ds <- addColumns(def_s, ds)
ds[, order := .SD[, .(order = .I)]$order, keyby = c_type]


di <- genCluster(ds, "site", 100, "id")
di <- trtAssign(di, 2, strata = "site", grp = "rx")
di <- addColumns(def_i, di)

dp <- di[, mean(y), keyby = .(site, order, rx, c_type)]

ggplot(data = dp, aes(x = factor(rx), y = V1, group = site)) +
  geom_line() +
  facet_grid(c_type ~ order) +
  theme(panel.grid = element_blank()) +
  ylim(0, 1)

#----

rt <- stanc("./Stan Code/binary_outcome.stan")
rt <- stanc("./Stan Code/binary_outcome_non_centered.stan")

sm <- stan_model(stanc_ret = rt, verbose=FALSE)

N <- nrow(di) ;                                ## number of observations
C <- di[, length(unique(c_type))]             ## of control types
K <- di[, length(unique(site))]               ## number of studies
y <- as.numeric(di$y)               ## individual outcome
kk <- di$site                                 ## study for individual
ctrl <- di$rx                                        ## treatment arm for individual
cc <- di[, .N, keyby = .(site, c_type)]$c_type       ## specific control arm for study

studydata <- list(
  N=N, C=C, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc)

fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                 cores = 4L, chains = 4, control = list(adapt_delta = 0.95))



print(fit, pars = c("alpha","delta_c", "Delta", "sigma_b", "eta_0"),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

print(fit, pars = c("beta_0"),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

#--- diagnostics

mcmc_trace(fit, pars = c("eta_0"))

posterior_ncp <- as.array(fit)
lp_ncp <- log_posterior(fit)
np_ncp <- nuts_params(fit)

color_scheme_set("darkgray")

mcmc_pairs(posterior_ncp, np = np_ncp, pars = c("Delta", "eta_0", "sigma_b"), 
           off_diag_args = list(size = 0.75))

scatter_eta_ncp <- mcmc_scatter(
  posterior_ncp, 
  pars = c("Delta", "eta_0"), 
  transform = list(eta_0 = "log"), # can abbrev. 'transformations'
  np = np_ncp, 
  size = 1
)
scatter_eta_ncp

scatter_s_ncp <- mcmc_scatter(
  posterior_ncp, 
  pars = c("Delta", "sigma_b"), 
  transform = list(sigma_b = "log"), # can abbrev. 'transformations'
  np = np_ncp, 
  size = 1
)
scatter_s_ncp

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_ncp, pars = "eta_0", np = np_ncp) + 
  xlab("Post-warmup iteration")

color_scheme_set("red")
mcmc_nuts_divergence(np_ncp, lp_ncp)
