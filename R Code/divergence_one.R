library(simstudy)
library(data.table)
library(rstan)

# ---

define_data <- function() {
  
  def_s <- defDataAdd(varname = "b0", formula = 0, variance = .05)
  
  def_s <- defDataAdd(
    def_s, 
    varname = "delta_k", 
    formula = "(c_type==1) * 0.4 + (c_type==2) * 0.5 + (c_type==3) * 0.6", 
    dist = "nonrandom"
  )
  
  def_i <- defDataAdd(
    varname = "y", 
    formula = "-1 + b0 + rx * delta_k", 
    dist = "binary", 
    link = "logit")
  
  list(def_s = def_s, def_i = def_i)
  
}
  
generate_data <- function(deflist) {
    
    print("generate data")
    
    dc <- genData(3, id = "c_type")
    
    ds <- genCluster(dc, "c_type", numIndsVar = 7, level1ID = "site")
    ds <- addColumns(deflist[['def_s']], ds)
    ds[, order := .SD[, .(order = .I)]$order, keyby = c_type]
    
    di <- genCluster(ds, "site", 20, "id")
    di <- trtAssign(di, 2, strata = "site", grp = "rx")
    di <- addColumns(deflist[['def_i']], di)
    
    di[]
    
}
  
estimate_pars <- function(dd, sm, ad) {
    
    print("prepare data as list")
    
    N <- nrow(dd) ;                                 ## number of observations
    C <- dd[, length(unique(c_type))]               ## of control types
    K <- dd[, length(unique(site))]                 ## number of studies
    y <- as.numeric(dd$y)                           ## individual outcome
    kk <- dd$site                                   ## study for individual
    ctrl <- dd$rx                                   ## treatment arm for individual
    cc <- dd[, .N, keyby = .(site, c_type)]$c_type  ## specific control arm for study
    
    sampdat <- list(N=N, C=C, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc)
    
    print("sample posterior")
    
    fit <-  sampling(sm, data = sampdat, seed = 332,
                     iter = 3000, warmup = 500, show_messages = FALSE,
                     cores = 1, chains = 1, refresh = 0,
                     control = list(adapt_delta = ad)
    )
    
    print("done sampling")
    
    return(1)
    
    
    
}
  
extract_data <- function(fit) {
    
    print("extract data")
    
    posterior <- as.array(fit)
    
    pDelta <- posterior[, , 'Delta']
    dsample <- data.table(sample = as.vector(pDelta))
    
    x <- summary(
      fit, 
      pars = c("alpha","delta_c", "Delta", "sigma_b", "eta_0"),
      probs = c(0.025, 0.5, 0.975)
    )
    
    dpars <- data.table(par = rownames(x$summary), x$summary)
    
    list(dpars, dsample)
}



rt_c <- stanc("~/r/binary_outcome.stan")
sm_c <- stan_model(stanc_ret = rt_c, verbose=FALSE)
  
xdef <- define_data()
dd <- generate_data(xdef)
ests <- estimate_pars(dd, sm_c, ad = 0.8)
# res <- extract_data(ests)

print(ests)
  
save(ests, file = "/gpfs/home/goldfk01/r/divergence_one.rda")

