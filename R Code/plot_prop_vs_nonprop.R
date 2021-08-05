library(collapse)
library(data.table)
library(ggthemes)
library(ggExtra)
library(ggpubr)

######Under proprotionality assumption#############
load("~/OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/predictive_iter2000.rda")

who <- unlist2d(lapply(site_plasma, function(x) x[[3]]), idcols = "replicate",DT = TRUE)
who <- who[,-1]
prop <- who/900
cumprop <- data.table(apply(prop,1,cumsum))

cumprop[,id:= seq(0,10,1)]
cum_prop <- melt(cumprop,id.vars = "id",measure.vars = colnames(cumprop)[-length(colnames(cumprop))])

load("~/OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/dind_local_v2.rda")

dsum <- dind[, .(N = sum(.N)), keyby = .(ordY)]
dsum[, prop := N/sum(N)]
dsum$cumprop <- cumsum(dsum$prop)
dsum[, who := seq(0,10,1)]

#################Relatex propontionality assumption
load("~/OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/dind_np_local.rda")
dsum_unprop <- dind[, .(N = sum(.N)), keyby = .(ordY)]
dsum_unprop[, prop := N/sum(N)]
dsum_unprop$cumprop <- cumsum(dsum_unprop$prop)
dsum_unprop[, who := seq(0,10,1)]

load("~/OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/np_predictive_iter2000.rda")

who <- unlist2d(lapply(site_plasma, function(x) x[[3]]), idcols = "replicate",DT = TRUE)
who <- who[,-1]
prop <- who/900
cumprop <- data.table(apply(prop,1,cumsum))

cumprop[,id:= seq(0,10,1)]
cum_nonprop <- melt(cumprop,id.vars = "id",measure.vars = colnames(cumprop)[-length(colnames(cumprop))])



# p1 <- ggplot(data=dsum,aes(x = who, y = cumprop,group=1)) +
#   geom_line(alpha = 1,color="#018571") +
#   geom_line(data = cum_prop, alpha = .01,color="#018571",
#             aes(x=id,y=value,group=variable)) + 
#   geom_line(data=dsum_unprop,alpha=1,color="#FE6900",aes(x = who, y = cumprop,group=1))+
#   geom_line(data = cum_nonprop, alpha = .01,color="#FE6900",
#             aes(x=id,y=value,group=variable)) + 
#   geom_point() +
#   ylab("cumulative probability") +
#   xlab("ordinal category") +
#   theme(panel.grid = element_blank(),
#         legend.position = "none") 

col <- c("#018571","#018571")
p1 <- ggplot(data=dsum,aes(x = who, y = cumprop,group=1)) +
  geom_line(alpha = 1,color=col[1]) +
  geom_line(data = cum_prop, alpha = .0035,color=col[1],
            aes(x=id,y=value,group=variable)) + 
  geom_point(color=col[1]) +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 8, face = "bold")) +
  ylab("cumulative probability") + labs(title="Proportional cumulative odds")+
  xlab("WHO score") +scale_x_continuous(breaks=0:10)+
  theme_minimal() +
  removeGridX()

p2 <- ggplot(data=dsum_unprop,aes(x = who, y = cumprop,group=1)) +
  geom_line(alpha = 1,color=col[2]) +
  geom_line(data = cum_nonprop, alpha = .0035,color=col[2],
            aes(x=id,y=value,group=variable)) + 
  geom_point(color=col[2]) + labs(title="Non-proportional cumulative odds")+
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 8, face = "bold")) +
  ylab("cumulative probability") +
  xlab("WHO score") +scale_x_continuous(breaks=0:10)+
  theme(panel.grid = element_blank(),
        legend.position = "none") +theme_minimal() +
  removeGridX()

ggarrange(p1,p2, nrow = 1,ncol=2)

