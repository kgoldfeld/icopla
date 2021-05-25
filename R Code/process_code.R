##stop the study if satisfied the stopping rule
#got data after running interim.R
library("dplyr")

load("./20210422/interim_000_2520.rda")
load("./20210424/interim_456_2520.rda")
load("./20210424/interim_123_2520.rda")
###############Find simulations sets: not converge when fitting frequentist model
data <- site_plasma %>% filter(n %in% c(15,30,45,60,75))
##stopping rule

data$`Pr(>|z|)_l` <- as.numeric(data$`Pr(>|z|)_l`)
data$`Pr(>|z|)_co` <- as.numeric(data$`Pr(>|z|)_co`)
#satisfied O'Brien-Fleming method (yes/no)

data$stop2 <- sapply(1:nrow(data), function(i) {
  n <- data[i,"n"]
  if(n==15){
    p_co <- ifelse(data[i,"Pr(>|z|)_co"] < 0.000005,1,0)
    p_l <- ifelse(data[i,"Pr(>|z|)_l"] < 0.000005,1,0)
    return(p_co*p_l)
  }
  else if(n==30){
    p_co <- ifelse(data[i,"Pr(>|z|)_co"] < 0.0013,1,0)
    p_l <- ifelse(data[i,"Pr(>|z|)_l"] < 0.0013,1,0)
    return(p_co*p_l)
  }
  else if(n==45){
    p_co <- ifelse(data[i,"Pr(>|z|)_co"] < 0.0085,1,0)
    p_l <- ifelse(data[i,"Pr(>|z|)_l"] < 0.0085,1,0)
    return(p_co*p_l)
  }
  else if(n==60){
    p_co <- ifelse(data[i,"Pr(>|z|)_co"] < 0.0228,1,0)
    p_l <- ifelse(data[i,"Pr(>|z|)_l"] < 0.0228,1,0)
    return(p_co*p_l)
  }
  else if(n==75){
    p_co <- ifelse(data[i,"Pr(>|z|)_co"] < 0.0417,1,0)
    p_l <- ifelse(data[i,"Pr(>|z|)_l"] < 0.0417,1,0)
    return(p_co*p_l)
  }
}) 


##No. of simulation sets with NA for frequentist methods
naiter <- subset(data,is.na(stop2))$iternum

#############Find simulation sets with >=1% divergence transitions when interim looks: 100% of sample
bayes_dat <- site_plasma
##stopping rule

##satisfied bayesian stopping rule and divergent transitions < 1% (yes:1/no:0)
bayes_dat$stop1 <- ifelse(bayes_dat$p.eff.co >= 0.95,1,0)*ifelse(bayes_dat$p.clinic.co >= 0.5,1,0)*ifelse(bayes_dat$p.eff.l >= 0.95,1,0)*ifelse(bayes_dat$p.clinic.l >= 0.5,1,0)*ifelse(bayes_dat$div_co<=100,1,0)*ifelse(bayes_dat$div_l<=100,1,0)

#bayes_dat$stop1 <- ifelse(bayes_dat$p.eff.co >= 0.95,1,0)*ifelse(bayes_dat$p.clinic.co >= 0.5,1,0)
##Bayesian CO model
#find simulation sets with >=1% divergence transitions when interim looks: 100% of sample
bayes_dat$div <- ifelse(bayes_dat$div_co <= (0.01*10000),1,0)*ifelse(bayes_dat$div_l <= (0.01*10000),1,0)
diviter <- subset(bayes_dat,(n %in% 75 & div %in% 0))$iternum
newdata_div <- bayes_dat[bayes_dat$iternum %in% c(diviter,naiter),]


#remove simulation sets with >= 1% divergence transitions and frequentist model does not converge
bayes_dat <- anti_join(bayes_dat,newdata_div,by="iternum")

#the number of simulation sets left
total_size <- length(unique(bayes_dat$iternum))


###if study satisfied the stopping rule at one interim look, for later interim look: stop1 =100

##stopping rule 1

for (i in 1:2520){
  ##stopping rule 1
  s <- bayes_dat %>% filter(iternum==i) %>% select(stop1)
  
  ##8 interim looks for Bayesian
  
  for (k in 1:8){
    if(is.na(s[k,])){
      break ###this iteration has been filtered out since not converge for frequentist method or divergence issue
    } else if (s[k,]%in%1){
      bayes_dat[bayes_dat$iternum==i,] <- bayes_dat %>% filter(iternum==i) %>%  mutate(stop1= ifelse(row_number() %in% c((k+1):9), 100, stop1))
      break
    }
  }
}



#rule1

##only include stop =0/1
size <- bayes_dat %>% filter(stop1 %in% c(0,1)) %>% group_by(n)%>%summarize(stop1_size=sum(stop1)/total_size,.groups = 'drop')
##change to %
p <- round(size$stop1_size*100,2)
p
#type one error or power
sum(p)


############Freq method typw I error rate or power
newdata <- data[data$iternum %in% c(diviter,naiter),]
#remove 
data <- anti_join(data,newdata,by="iternum")
length(unique(data$iternum))


#the number of simulation sets left
total_size <- length(unique(data$iternum))


##frquentist stopping rule
for (i in 1:2520){
  
  f <- data %>% filter(iternum==i) %>% select(stop2)
  for (k in 1:4){
    if(is.na(f[k,])){
      break
    } else if (f[k,]%in% 1){
      data[data$iternum==i,] <- data %>% filter(iternum==i) %>%  mutate(stop2= ifelse(row_number() %in% c((k+1):5), 100, stop2))
      break
    }
  }
}


#rule2
size <- data %>% filter(stop2 %in% c(0,1)) %>% group_by(n)%>%summarize(stop2_size=sum(stop2)/total_size,
                                                                       .groups = 'drop')
##change to %
p <- round(size$stop2_size*100,2)
p
#type one error or power
sum(p)




