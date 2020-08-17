##stop the study if satisfied the stopping rule
#got data after running interim.R
load("./R Code/20200815/interim_nc_000.rda")
load("./R Code/20200815/interim_nc_123.rda")
load("./R Code/20200815/interim_nc_234.rda")
load("./R Code/20200815/interim_nc_345.rda")
load("./R Code/20200815/interim_nc_456.rda")

#remove Delta 
data <- subset(site_plasma,select=-Delta.perc)
data <- distinct(data)

##stopping rule

##satisfied stopping rule: stop1=1
data$stop1 <- ifelse(data$p.eff >= 0.95,1,0)*ifelse(data$p.clinic >= 0.5,1,0)
data$stop2 <- ifelse(data$p.eff >= 0.9,1,0)*ifelse(data$p.clinic >= 0.6,1,0)
#Frequentist fitting satisfied p_value <0.05
data$sig_005 <- ifelse(data$`Pr(>|z|)` < 0.05,1,0)
#Frequentist fitting satisfied p_value <0.025
data$sig_0025 <- ifelse(data$`Pr(>|z|)` < 0.025,1,0)
#Frequentist fitting satisfied p_value <0.025
data$sig_001 <- ifelse(data$`Pr(>|z|)` < 0.01,1,0)

##no. of simulation sets with NA
naiter <- subset(data,is.na(sig_001))$iternum
newdata <- data[iternum %in% naiter,]

#remove simulation sets with NA
data <- anti_join(data,newdata,by="iternum")

total_size <- length(unique(data$iternum))


###if study satisfied the stopping rule at one interim look, for later interim look: stop1 =100

##stopping rule 1

for (i in 1:990){
  ##stopping rule 1
  s <- data %>% filter(iternum==i) %>% select(stop1)
  
  ##8 interim looks for Bayesian
<<<<<<< HEAD
  for (k in 1:7){
    if(is.na(s[k,])){
      break
    } else if (s[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(stop1= ifelse(row_number() %in% c((k+1):8), 100, stop1))
=======
  for (k in 1:8){
    if(is.na(s[k,])){
      break
    } else if (s[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(stop1= ifelse(row_number() %in% c((k+1):9), 100, stop1))
>>>>>>> 8b677725a4ce62d1748ced3a66b06f9700e414d7
      break
    }
  }
}

##stopping rule 2
for (i in 1:990){
  ##stopping rule 2
  s2 <- data %>% filter(iternum==i) %>% select(stop2)
  
<<<<<<< HEAD
  for (k in 1:7){
    if(is.na(s2[k,])){
      break
    } else if (s2[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(stop2= ifelse(row_number() %in% c((k+1):8), 100, stop2))
=======
  for (k in 1:8){
    if(is.na(s2[k,])){
      break
    } else if (s2[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>% mutate(stop2= ifelse(row_number() %in% c((k+1):9), 100, stop2))
>>>>>>> 8b677725a4ce62d1748ced3a66b06f9700e414d7
      break
    }
  }
}



##stopping rule: p-value<0.05
for (i in 1:990){
  
  f5 <- data %>% filter(iternum==i) %>% select(sig_005)
  
<<<<<<< HEAD
  for (k in 1:7){
    if(is.na(f5[k,])){
      break
    } else if (f5[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_005= ifelse(row_number() %in% c((k+1):8), 100, sig_005))
=======
  for (k in 1:8){
    if(is.na(f5[k,])){
      break
    } else if (f5[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_005= ifelse(row_number() %in% c((k+1):9), 100, sig_005))
>>>>>>> 8b677725a4ce62d1748ced3a66b06f9700e414d7
      break
    }
  }
}


##stopping rule: p-value<0.025
for (i in 1:990){
  
  f25 <- data %>% filter(iternum==i) %>% select(sig_0025)
  
<<<<<<< HEAD
  for (k in 1:7){
    if(is.na(f25[k,])){
      break
    } else if (f25[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_0025= ifelse(row_number() %in% c((k+1):8), 100, sig_0025))
=======
  for (k in 1:8){
    if(is.na(f25[k,])){
      break
    } else if (f25[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_0025= ifelse(row_number() %in% c((k+1):9), 100, sig_0025))
>>>>>>> 8b677725a4ce62d1748ced3a66b06f9700e414d7
      break
    }
  }
}


##stopping rule: p-value<0.01
for (i in 1:990){
  
  f1 <- data %>% filter(iternum==i) %>% select(sig_001)
  
<<<<<<< HEAD
  for (k in 1:7){
    if(is.na(f1[k,])){
      break
    } else if (f1[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_001= ifelse(row_number() %in% c((k+1):8), 100, sig_001))
=======
  for (k in 1:8){
    if(is.na(f1[k,])){
      break
    } else if (f1[k,]==1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_001= ifelse(row_number() %in% c((k+1):9), 100, sig_001))
>>>>>>> 8b677725a4ce62d1748ced3a66b06f9700e414d7
      break
    }
  }
}


#fraction of simulations satisfied stopping rule
#interim look 16%;33%;40%;50%;60%;67%;80%;90%;100% for Bayesian
#33%,50%,67%,90%,100% for frequentist method

#rule1

##only include stop =0/1
size <- data %>% filter(stop1 ==0 | stop1 == 1) %>% group_by(n)%>%summarize(stop1_size=sum(stop1)/total_size,
                                                                            .groups = 'drop')%>% as.data.frame()
##change to %
p <- round(size$stop1_size*100,2)
p
#type one error or power
sum(p)

#rule2
size <- data %>% filter(stop2 ==0 | stop2 == 1) %>% group_by(n)%>%summarize(stop2_size=sum(stop2)/total_size,
                                                                            .groups = 'drop')
##change to %
p <- round(size$stop2_size*100,2)
p
#type one error or power
sum(p)

##P<0.05
size <- data %>% filter(sig_005 ==0 | sig_005 == 1) %>% group_by(n)%>%summarize(stop3_size=sum(sig_005)/total_size,
                                                                                .groups = 'drop')
##change to %
p <- round(size$stop3_size*100,2)
p
#type one error or power
sum(p)

##P<0.025
size <- data %>% filter(sig_0025 ==0 | sig_0025 == 1) %>% group_by(n)%>%summarize(stop4_size=sum(sig_0025)/total_size,
                                                                                  .groups = 'drop')
##change to %
p <- round(size$stop4_size*100,2)
p
#type one error or power
sum(p)
##P<0.01
size <- data %>% filter(sig_001 ==0 | sig_001 == 1) %>% group_by(n)%>%summarize(stop5_size=sum(sig_001)/total_size,
                                                                                .groups = 'drop')
##change to %
p <- round(size$stop5_size*100,2)
p
#type one error or power
sum(p)

##check divergence issue
library("ggplot2")

p <- ggplot(data, aes(x=as.factor(n), y=(div/12000)*100,group=as.factor(n))) 
p + geom_boxplot(outlier.size=1, notch=FALSE)+
  scale_x_discrete(labels=c("16","33","40","50","60","67","80","90","100"))+
  labs(title="% of divergent transitions of each interim look",x="% of sample", 
       y = "% of divergent transitions")+
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) 

data %>%  group_by(n)%>%summarize(stop5_size=max(div)/12000, .groups = 'drop')
