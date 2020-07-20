
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

#### Trace plots

plot(fit, plotfun = "trace", pars = c("delta_k"), 
     inc_warmup = FALSE, ncol = 3)

plot(fit, plotfun = "trace", pars = c("delta", "Delta"), 
     inc_warmup = FALSE, ncol = 1)


plot(fit, plotfun = "trace", pars = c("eta_0", "eta"), 
     inc_warmup = FALSE, ncol = 1)

### Tau (cut-points)

params <- extract(fit)

getTau <- function(study) {
  x <- data.table(params$tau[, study, ])
  x <- melt(x, measure.vars = 1:ncol(x),
    variable.name = "level", value.name = "tau",
    variable.factor = TRUE, verbose = FALSE) 
  x[, level := as.numeric(level)]
  x[]
}

p <- list()

for (i in 1:9) {
  p[[i]] <- ggplot(data = getTau(i)) +
    geom_density(aes(x=tau, fill = level, group = level), color = c_dark_highlight) +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.title = element_blank()) +
    scale_fill_gradient2(low = c_light, high = c_dark)
}

gridExtra::grid.arrange(grobs = p, nrow = 3)

### Study effects

delta_k <- data.table(params$delta_k)
delta_k <- melt(delta_k, measure.vars = 1:ncol(delta_k),
      variable.name = "study", value.name = "delta_k",
      variable.factor = TRUE)  

delta_k[, study := as.numeric(study)]
delta_k <- merge(delta_k, dstudy, by = "study")

library(paletteer)

ggplot(data = delta_k, aes(x=delta_k, group = study)) +
  geom_density(aes(fill = factor(C)), alpha = .6) +
  facet_grid(C ~ .) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_paletteer_d("wesanderson::Moonrise2") +
  xlim(-0.5, 2.5)

### Control-type efffects

delta <- data.table(params$delta)
delta <- melt(delta, measure.vars = 1:ncol(delta),
              variable.name = "control", value.name = "delta",
              variable.factor = TRUE)  

delta[, control := as.numeric(control)]

ggplot(data = delta, aes(x=delta)) +
  geom_density(aes(fill = factor(control))) +
  facet_grid(control ~ .) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_paletteer_d("wesanderson::Moonrise2") +
  xlim(-0.5, 2.5)

### Overall effect

Delta <- data.table(Delta = params$Delta)

ggplot(data = Delta, aes(x=Delta)) +
  geom_density(fill = c_dark)+
  theme(panel.grid = element_blank()) +
  xlim(-3, 3) 

### OR 

OR <- data.table(OR = params$OR)

ggplot(data = OR, aes(x=OR)) +
  geom_density(fill = c_dark)+
  theme(panel.grid = element_blank()) +
  xlim(0,1.5) 
