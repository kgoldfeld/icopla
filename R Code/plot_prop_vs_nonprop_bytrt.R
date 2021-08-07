library(collapse)
library(data.table)
library(ggthemes)
library(ggExtra)
library(ggpubr)

######Under proprotionality assumption#############
load("~/OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/by_trt/predictive_iter9900_v2.rda")
#load("./data/predictive_iter2000.rda")
who <- unlist2d(site_plasma,idcols="replicate")
who <- as.data.table(who)

subwho <- subset(who,n_draw %in% 1:2000)
subwho[, cumprob := cumsum(N)/sum(N) , keyby = .(n_draw,control)]

load("~/OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/by_trt/dind_local_v2.rda")
#load("./data/dind_local_v2.rda")

dsum <- dind[, .(N = sum(.N)), keyby = .(control,ordY)]

dsum[, cumprob := cumsum(N)/sum(N) , keyby = control]



#########Plot###########
p1 <- ggplot(data=dsum,aes(x = ordY, y = cumprob,group= control,color=factor(control))) +
  geom_line(alpha = 1) +
  geom_line(data = subwho, alpha = .002,
            aes(group=interaction(control,n_draw))) +
  geom_point() +
  labs(x= "WHO score",y="cumulative probability",color="Treatment",
       subtitle="Proportional cumulative odds")+
  scale_color_manual(labels=c("CCP","Control"),
                     values=c("red","#018571"))+
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 8, face = "bold")) +
  ylab("cumulative probability") +
  xlab("WHO score") +scale_x_discrete(labels= as.character(factor(0:10)))+
  theme(panel.grid = element_blank(),
        legend.position = "none") +theme_minimal() +
  removeGridX()+ theme(legend.position="bottom")

#################Relatex propontionality assumption
load("~/OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/by_trt/np_predictive_iter9900_v4.rda")
#load("./data/predictive_iter2000.rda")
who <- unlist2d(site_plasma,idcols="replicate")
who <- as.data.table(who)
subwho <- subset(who,n_draw %in% 1:2000)
subwho[, cumprob := cumsum(N)/sum(N) , keyby = .(n_draw,control)]

load("~/OneDrive - NYU Langone Health/Covid_19_project/Convalescent plasma project/Simulation_paper/Sim_redo/predictive_posterior_prob/data/by_trt/dind_np_local_v3.rda")
#load("./data/dind_local_v2.rda")

dsum <- dind[, .(N = sum(.N)), keyby = .(control,ordY)]

dsum[, cumprob := cumsum(N)/sum(N) , keyby = control]

p2 <- ggplot(data=dsum,aes(x = ordY, y = cumprob,group= control,color=factor(control))) +
  geom_line(alpha = 1) +
  geom_line(data = subwho, alpha = .002,
            aes(group=interaction(control,n_draw))) +
  geom_point() +
  labs(x= "WHO score",y="cumulative probability",color="Treatment",
       subtitle="Non-proportional cumulative odds")+
  scale_color_manual(labels=c("CCP","Control"),
                     values=c("red","#018571"))+
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 8, face = "bold")) +
  ylab("cumulative probability") +
  xlab("WHO score") +scale_x_discrete(labels= as.character(factor(0:10)))+
  theme(panel.grid = element_blank(),
        legend.position = "none") +theme_minimal() +
  removeGridX()+ theme(legend.position="bottom")
# 
 ggarrange(p1,p2, nrow = 1,ncol=2,common.legend = TRUE,legend = "bottom")

