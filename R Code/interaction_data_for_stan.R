######################
#Data for Stan model #
######################

##dind is a dataset with available outcome and parsimounious covariates.
#Parsimonious adjustment includes age, sex, WHO score at baseline, days since symptom onset and enrollment
#quarters.
dind <- dind%>% as.data.table()

dind$ordY <- dind$WHO_14_imp + 1  #ordY:1-11
dind$who_enroll <- dind$out_who0 + 1 # since we change the scale from 0-10 to 1-11
dind$cat_age <- cut(dind$bl_age, c(0,50,65,120))
dind$age <- factor(dind$cat_age,c("(0,50]","(50,65]","(65,120]"),1:3)
dind$sex <- dind$bl_sex #0:male; 1:female
dind$ss <- dind$bl_symp_days
dind$qtr <- dind$bl_enrollqtr
dind$study <- factor(dind$tr_rct_id,c("AA","BB","CC","DD","EE","GG","KK","RR"),1:8) %>% unclass(.)##
dind$C <- factor(dind$tr_cntrl,c(0,1,2),1:3) %>% unclass(.) #control arm for study; in SAP:0:standard of care;1:Non_CP;2:Saline


#########corresponding Stan file: interaction_m1/2_sex2levels
#####Sex: 2 levels categorical variable
N <- nrow(dind)                                ## number of observations 
L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
K <- dind[, length(unique(study))]            ## number of studies 
y <- as.numeric(dind$ordY)                    ## individual outcome 
y_2 <- ifelse(y %in% seq(8,11),1,0)
kk <- dind$study                              ## study for individual 
ctrl <- ifelse(dind$bl_treatment==1,0,1)      ## treatment arm for individual 
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = dind)[, -1]
D <- ncol(x)
ds <- model.matrix(ordY ~ factor(sex), data = dind)[, -1] ##who_enroll: 2 levels
ds <- as.matrix(ds)

studydata <- list(N=N, L= L, K=K, y=y, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x)

#########corresponding Stan file: interaction_m1/2_symp5levels
#####Symptom duration: 5 levels categorical variable
N <- nrow(dind)                                ## number of observations 
L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
K <- dind[, length(unique(study))]            ## number of studies 
y <- as.numeric(dind$ordY)                    ## individual outcome 
y_2 <- ifelse(y %in% seq(8,11),1,0)
kk <- dind$study                              ## study for individual 
ctrl <- ifelse(dind$bl_treatment==1,0,1)      ## treatment arm for individual 
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = dind)[, -1]
D <- ncol(x)
ds <- model.matrix(ordY ~ factor(ss), data = dind)[, -1] ##ss: 5 levels


studydata <- list(N=N, L= L, K=K, y=y, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x)

#########corresponding Stan file: interaction_m1/2_age3levels
#####Age: 3 levels categorical variable
N <- nrow(dind)                                ## number of observations 
L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
K <- dind[, length(unique(study))]            ## number of studies 
y <- as.numeric(dind$ordY)                    ## individual outcome 
y_2 <- ifelse(y %in% seq(8,11),1,0)
kk <- dind$study                              ## study for individual 
ctrl <- ifelse(dind$bl_treatment==1,0,1)      ## treatment arm for individual 
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = dind)[, -1]
D <- ncol(x)
ds <- model.matrix(ordY ~ factor(age), data = dind)[, -1] ##age: 3 levels


studydata <- list(N=N, L= L, K=K, y=y, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x)

#########corresponding Stan file: interaction_m1/2_who_enroll3levels
#####WHO_enroll: 3 levels categorical variable
N <- nrow(dind)                                ## number of observations 
L <- dind[, length(unique(ordY))]             ## number of levels of outcome 
K <- dind[, length(unique(study))]            ## number of studies 
y <- as.numeric(dind$ordY)                    ## individual outcome 
y_2 <- ifelse(y %in% seq(8,11),1,0)
kk <- dind$study                              ## study for individual 
ctrl <- ifelse(dind$bl_treatment==1,0,1)      ## treatment arm for individual 
cc <- dind[, .N, keyby = .(study, C)]$C       ## specific control arm for study 
x <- model.matrix(ordY ~ factor(who_enroll) + factor(age) + factor(sex) + factor(ss) + factor(qtr), data = dind)[, -1]
D <- ncol(x)
ds <- model.matrix(ordY ~ factor(who_enroll), data = dind)[, -1] ##who_enroll: 3 levels

studydata <- list(N=N, L= L, K=K, y=y, y_2=y_2, kk=kk, ctrl=ctrl, cc=cc, ds=ds,D=D, x=x)

