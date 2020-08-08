##stop the study if satisfied the stopping rule

#got data after running interim.R
#load("./8_4/20200803/interim_123.rda")
load("./8_4/20200803/interim_234.rda")
#load("./8_4/20200803/interim_345.rda")
#load("./8_4/20200803/interim_456.rda")

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

for (i in 1:total_size){
  ##stopping rule 1
  s <- data %>% filter(iternum==i) %>% select(stop1)
  
  ##8 interim looks for Bayesian
  for (k in 1:7){
    
    if (s[k]==1){
      data[iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(stop1= ifelse(row_number() %in% c((k+1):8), 100, stop1))
      break
    }
  }
}

##stopping rule 2
for (i in 1:total_size){
  ##stopping rule 2
  s2 <- data %>% filter(iternum==i) %>% select(stop2)
  
  for (k in 1:7){
    
    if (s2[k]==1){
      data[iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(stop2= ifelse(row_number() %in% c((k+1):8), 100, stop2))
      break
    }
  }
}



##stopping rule: p-value<0.05
for (i in 1:total_size){
  
  f5 <- data %>% filter(iternum==i) %>% select(sig_005)
  
  for (k in 1:7){
    
    if (f5[k]==1){
      data[iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_005= ifelse(row_number() %in% c((k+1):8), 100, sig_005))
      break
    }
  }
}


##stopping rule: p-value<0.025
for (i in 1:total_size){
  
  f25 <- data %>% filter(iternum==i) %>% select(sig_0025)
  
  for (k in 1:7){
    
    if (f25[k]==1){
      data[iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_0025= ifelse(row_number() %in% c((k+1):8), 100, sig_0025))
      break
    }
  }
}


##stopping rule: p-value<0.001
for (i in 1:total_size){
  
  f1 <- data %>% filter(iternum==i) %>% select(sig_001)
  
  for (k in 1:7){
    
    if (f1[k]==1){
      data[iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(sig_001= ifelse(row_number() %in% c((k+1):8), 100, sig_001))
      break
    }
  }
}


#fraction of simulations satisfied stopping rule
#rule1

##only include stop =0/1
size <- data %>% filter(stop1 ==0 | stop1 == 1) %>% group_by(n)%>%summarize(stop1_size=sum(stop1)/total_size,
                                                                            .groups = 'drop')
##change to %
round(size*100,3)

#rule2
size <- data %>% filter(stop2 ==0 | stop2 == 1) %>% group_by(n)%>%summarize(stop2_size=sum(stop2)/total_size,
                                                                            .groups = 'drop')
##change to %
round(size*100,3)

##P<0.05
size <- data %>% filter(sig_005 ==0 | sig_005 == 1) %>% group_by(n)%>%summarize(stop3_size=sum(sig_005)/total_size,
                                                                                .groups = 'drop')
##change to %
round(size*100,3)

##P<0.025
size <- data %>% filter(sig_0025 ==0 | sig_0025 == 1) %>% group_by(n)%>%summarize(stop4_size=sum(sig_0025)/total_size,
                                                                                  .groups = 'drop')
##change to %
round(size*100,3)

##P<0.01
size <- data %>% filter(sig_001 ==0 | sig_001 == 1) %>% group_by(n)%>%summarize(stop5_size=sum(sig_001)/total_size,
                                                                                .groups = 'drop')
##change to %
round(size*100,3)