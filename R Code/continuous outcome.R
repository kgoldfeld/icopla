def_c <- defData(varname = "delta_c", formula = 3, variance = 1, id = "c_type")

def_s <- defDataAdd(varname = "b0", formula = 0, variance = 1)
def_s <- defDataAdd(def_s, varname = "delta_k", formula = "delta_c", variance = 2)

def_i <- defDataAdd(varname = "y", formula = "b0 + rx * delta_k", variance = 4)

dc <- genData(3, def_c)

ds <- genCluster(dc, "c_type", numIndsVar = 3, level1ID = "site")
ds <- addColumns(def_s, ds)
ds[, order := .SD[, .(order = .I)]$order, keyby = c_type]

di <- genCluster(ds, "site", 80, "id")
di <- trtAssign(di, 2, strata = "site", grp = "rx")
di <- addColumns(def_i, di)

ggplot(data = di, aes(x = factor(rx), y = y)) +
  geom_jitter(height = 0, width = .3, size = .75) +
  facet_grid(c_type ~ order) +
  theme(panel.grid = element_blank())

#----

rt <- stanc("./Stan Code/continuous_outcome.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

N <- nrow(di) ;                                ## number of observations
C <- di[, length(unique(c_type))]             ## of control types
K <- di[, length(unique(site))]               ## number of studies
y <- as.numeric(scale(di$y, scale = FALSE))               ## individual outcome
kk <- di$site                                 ## study for individual
ctrl <- di$rx                                        ## treatment arm for individual
cc <- di[, .N, keyby = .(site, c_type)]$c_type       ## specific control arm for study

studydata <- list(
  N=N, C=C, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc)

fit <-  sampling(sm, data=studydata, iter = 3000, warmup = 500, 
                 cores = 4L, chains = 4, control = list(adapt_delta = 0.8))

print(fit, pars = c("alpha","delta_c", "Delta", "sigma", "sigma_b", "eta_0"),
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2)

traceplot(fit, pars = c("sigma", "alpha", "delta_c", "Delta","delta_k", "eta_0"))

pairs(fit, pars = c("delta_c", "Delta","eta_0", "sigma_b", "sigma"))


