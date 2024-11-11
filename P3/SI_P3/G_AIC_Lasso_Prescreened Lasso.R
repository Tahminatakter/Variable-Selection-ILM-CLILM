library(EpiILM)
library(formattable)
#install.packages("adaptMCMC")
require(adaptMCMC)
library('R.utils')
#install.packages("glmnet")
library('glmnet')
library(MCMCpack)
###################Variable Screening methods (under SI framework)
#AIC, Lasso, Pre-screened Lasso

Bcoef<-numeric()
Bocoef<-numeric()
Fcoef<-matrix(0,22,21)
Fcoef[,21]<-0:21
Bic<-numeric()

#Scenaio 1 with cov1, cov2, cov3
#B1<-c(5.5,5.0,5.0,5.0,4.5,5.0,5.5,5.0,4.5,5.0,5.0,5.5,5.0,4.5,4.5,4.5,6.0,5.5,5.0,5.0)
#B2<-c(5.5,5.0,5.5,5.0,5.0,5.0,5.0,6.0,5.0,5.0,4.5,6.0,4.5,5.0,4.5,5.0,7.5,5.5,5.5,5.0)
#B3<-c(6.5,4.5,5.0,5.0,5.0,5.5,5.0,4.5,5.0,5.0,5.0,5.5,4.5,5.5,5.0,5.0,6.5,5.0,5.5,4.5)
#Scenaio 2 with cov1, cov2, cov3
#B4<-c(5.5,5.0,5.0,4.5,5.0,5.0,5.0,4.5,5.0,4.5,4.5,5.0,4.5,5.0,5.0,5.0,6.5,6.0,5.0,5.0)
#B5<-c(5.0,5.0,5.0,5.0,5.0,5.0,5.0,4.5,5.0,5.0,4.5,5.5,5.0,5.0,5.0,5.5,8.0,6.0,5.5,5.0)
#B6<-c(5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,4.5,5.0,4.5,4.5,4.5,5.0,6.5,5.0,6.0,5.0)
#Scenaio 1 with cov7 and Scenaio 2 with cov7 
B<-c(6.0,5.5,5.5,5.5,5.5,5.5,5.5,5.0,5.0,5.0,5.0,5.5,4.5,6.0,5.0,5.0,6.0,5.0,6.0,5.0)
#B72<-c(5.5,5.0,5.5,4.5,5.0,5.0,5.0,4.5,5.0,5.0,6.0,5.5,5.0,5.0,5.0,5.0,6.0,5.0,5.5,5.5)


for (gg in 1:20){

x<-runif(500,0,10)
y<-runif(500,0,10)
A1<-runif(500,0,5)
A2<-runif(500,0,5)
A3<-runif(500,0,5)
A4<-runif(500,0,5)
A5<-runif(500,0,5)
A6<-runif(500,0,5)
A7<-runif(500,0,5)
A8<-runif(500,0,5)
A9<-runif(500,0,5)
A10<-runif(500,0,5)
A11<-runif(500,0,5)  
A12<-runif(500,0,5)
A13<-runif(500,0,5)  
A14<-runif(500,0,5)  
A15<-runif(500,0,5)
A16<-runif(500,0,5)  
A17<-runif(500,0,5)
A18<-runif(500,0,5)
A19<-runif(500,0,5)
A20<-runif(500,0,5)

#Data generate
out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)
#out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)

#out_cov
t_end<-max(out_cov$inftime)

SI<-rep(0)
for (i in 1:t_end){
  SI[i]<-length(out_cov$inftime[out_cov$inftime>i])
}

#Number of new infections
newinf_true2<-rep(0)
for (i in 1:t_end){
  newinf_true2[i]<-length(out_cov$inftime[out_cov$inftime==i])
}

newinf_true2

#cumulative infections
newinf_true<-cumsum(newinf_true2)

#################Logistic

#For all t
ndist<-matrix(0,500, t_end)

for(t in 1:t_end ){
  
  xx1<-matrix(0,500, newinf_true[t])
  
  for(i in 1:500) {
    
    for(j in 1:newinf_true[t])
    {
      xx1[i,j] <- (dist(rbind(c(x[i],y[i]),c(x[which(out_cov$inftime<(t+1) & out_cov$inftime!=0)[j]],y[which(out_cov$inftime<(t+1) & out_cov$inftime!=0)[j]]))))^(-B[gg])
      
    }
    ndist[i,t]<-log(sum(xx1[i,])) 
    
  }
  
}



#####
ind<-rep(1:500, each=(t_end-1))
t<-rep(2:t_end,500)
inf<-rep(out_cov$inftime, each=(t_end-1))
X<-rep(x,each=(t_end-1))
Y<-rep(y,each=(t_end-1))
#
AA1<-rep(A1,each=(t_end-1))
AA2<-rep(A2,each=(t_end-1))
AA3<-rep(A3,each=(t_end-1))
AA4<-rep(A4,each=(t_end-1))
AA5<-rep(A5,each=(t_end-1))
AA6<-rep(A6,each=(t_end-1))
AA7<-rep(A7,each=(t_end-1))
AA8<-rep(A8,each=(t_end-1))
AA9<-rep(A9,each=(t_end-1))
AA10<-rep(A10,each=(t_end-1))

AA11<-rep(A11,each=(t_end-1))
AA12<-rep(A12,each=(t_end-1))
AA13<-rep(A13,each=(t_end-1))
AA14<-rep(A14,each=(t_end-1))
AA15<-rep(A15,each=(t_end-1))
AA16<-rep(A16,each=(t_end-1))
AA17<-rep(A17,each=(t_end-1))
AA18<-rep(A18,each=(t_end-1))
AA19<-rep(A19,each=(t_end-1))
AA20<-rep(A20,each=(t_end-1))

dat1<-data.frame(ind,X,Y,t,inf,AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,AA9,AA10,AA11,AA12,AA13,AA14,AA15,AA16,AA17,AA18,AA19,AA20)



dat1$inf2<-rep(0,500*(t_end-1))

for(i in 1:(500*(t_end-1))){
  if(t[i]==inf[i]){dat1$inf2[i]<-1}
}



dat1$dis1<-rep(0,500*(t_end-1))


for(t in 2: t_end){
  j=t
  for(i in 1:500){
    dat1$dis1[j-1]<-ndist[i,t-1]
    j=j+(t_end-1)
    
  }
}


# dat1
dat2<-data.frame(individual=dat1$ind, time=dat1$t, Y=dat1$inf2, distance=dat1$dis1,AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,AA9,AA10,AA11,AA12,AA13,AA14,AA15,AA16,AA17,AA18,AA19,AA20)
head(dat2)

###
for(i in 1:(500*(t_end-1))){
  dat2$in1[i]<-inf[i]
}

dat3<-dat2[dat2$in1==0 | dat2$time<=dat2$in1,]

dat2<-dat3
head(dat2)

##############
XX<-matrix(0,length(dat2$distance),21)
for(i in 1:length(dat2$distance)){
  for(j in 1:20){
    
    XX[i,j]<-dat2[i,j+4]
  }
}

for(i in 1:length(dat2$distance)){
  XX[i,21]<-dat2$distance[i]
}

dim(XX)
head(XX)


######################
####Forward, Backward and both step AIC
######################

r2<-glm(dat2$Y~XX[,1]+XX[,2]+XX[,3]+XX[,4]+XX[,5]+XX[,6]+XX[,7]+XX[,8]+XX[,9]+XX[,10]+XX[,11]+XX[,12]+XX[,13]+XX[,14]+XX[,15]+XX[,16]+XX[,17]+XX[,18]+XX[,19]+XX[,20]+XX[,21],family="binomial")
#summary(r2)

#To fix distance covariate
scope <- list(lower = ~ XX[,21], upper = ~ XX[,1]+XX[,2]+XX[,3]+XX[,4]+XX[,5]+XX[,6]+XX[,7]+XX[,8]+XX[,9]+XX[,10]+XX[,11]+XX[,12]+XX[,13]+XX[,14]+XX[,15]+XX[,16]+XX[,17]+XX[,18]+XX[,19]+XX[,20]+XX[,21])


#Backward stepwise AIC
Br2_step<-step(r2,scope = scope)
summary(Br2_step)
Bcoef[gg]<-summary(Br2_step)[1]

# both step AIC
Bor2_step<-step(r2, scope = scope, direction="both")
summary(Bor2_step)
Bocoef[gg]<-summary(Bor2_step)[1]

#Forward stepwise AIC
Fr2_step<-step(r2, scope = scope, direction="forward")
summary(Fr2_step)
Fcoef[,gg]<-summary(Fr2_step)$coefficient[,4]<0.1

#BIC method
Bir2_step<-step(r2,scope = scope, k=log(length(dat2$Y)))
summary(Bir2_step)
Bic[gg]<-summary(Bir2_step)[1]

#Backward
Bcoef
#Both
Bocoef
#Forward
Fcoef
Fcoef1<-Fcoef[2:22,]
Fcoef1
d1<-data.frame(Fcoef1)
formattable(d1)
#BIC
Bic

######################
##################Lasso
#####################


#Without Intercept
#r1_cv<-cv.glmnet(XX,dat2$Y,family="binomial",penalty.factor=c(rep(1,20),0),intercept=F,alpha=1)

#With Intercept
r1_cv<-cv.glmnet(XX,dat2$Y,family="binomial",penalty.factor=c(rep(1,20),0),alpha=1)
#r1_cv$lambda.min
#r1_cv$lambda.1se

#lambda.min: the lambda at which the minimal MSE is achived
#lambda.1se: the largest lambda at which the MSE is within 1se of the minimal MSE
Fsig_min[,gg]<-(coef(r1_cv, r1_cv$lambda.min)!=0)[,1]
Fsig_se[,gg]<-(coef(r1_cv, r1_cv$lambda.1se)!=0)[,1]

}

Csig_min<-Fsig_min[2:22,]
Csig_se<-Fsig_se[2:22,]


############For Cov7 (similarly for Cov1 and Cov2)
data1<-data.frame(Csig_min)
formattable(data1)
#######
sen<-numeric()
spe<-numeric()
acc<-numeric()

for(k in 1:20){
  sen[k]<-(sum(data1[1:7,k]))/7
}

for(k in 1:20){
  spe[k]<-(13-sum(data1[8:20,k]))/13
}

for(k in 1:20){
  acc[k]<-(sum(data1[1:7,k])+13-sum(data1[8:20,k]))/20
}


sdata1<-data.frame(sen,spe,acc)
formattable(sdata1)

mean(sen)
mean(spe)
mean(acc)

sd(sen)
sd(spe)
sd(acc)

#Calculate sensitivity,specificity and accuracy

data1<-data.frame(Csig_se)
formattable(data1)
#######
sen<-numeric()
spe<-numeric()
acc<-numeric()

for(k in 1:20){
  sen[k]<-(sum(data1[1:3,k]))/3
}

for(k in 1:20){
  spe[k]<-(17-sum(data1[4:20,k]))/17
}

for(k in 1:20){
  acc[k]<-(sum(data1[1:3,k])+17-sum(data1[4:20,k]))/20
}


sdata1<-data.frame(sen,spe,acc)
formattable(sdata1)

mean(sen)
mean(spe)
mean(acc)


######################
############Pre-screened Lasso
######################


#Use single epidemic and then run Pre-screened Lasso
###############Find sig. covariates in first stage

pp1<-numeric()
for(k in 1:20){
  m1<-glm(dat2$Y~XX[,21]+XX[,k], family=binomial(link = "logit"))
  pp1[k]<-summary(m1)$coefficient[3,4]
}

pp1

data2<-data.frame(Variable=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"),pp1)
formattable(data2)

####
or1<-ifelse(data2$pp1<0.05,1,0)
or2<-which(or1==1)
XX1<-XX[,c(or2,21)]
colnames(XX1)<-c(or2,21)

#r1_cv<-cv.glmnet(XX1,dat2$Y,family="binomial",penalty.factor=c(rep(1,(ncol(XX1)-1)),0),intercept=F,alpha=1)
r1_cv<-cv.glmnet(XX1,dat2$Y,family="binomial",penalty.factor=c(rep(1,(ncol(XX1)-1)),0),alpha=1)


Fsig_min<-(coef(r1_cv, r1_cv$lambda.min)!=0)[,1]
Fsig_se<-(coef(r1_cv, r1_cv$lambda.1se)!=0)[,1]  


Fsig_min
Fsig_se

