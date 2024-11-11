library(EpiILM)
library(formattable)
#install.packages("adaptMCMC")
require(adaptMCMC)
library(MCMCpack)

###################Find Beta under SI framework
Gpow<-numeric()
GAIC<-numeric()

for (gg in 1:20){

AIC1<-numeric()

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

#cumulative infections
newinf_true<-cumsum(newinf_true2)

#################

ndist<-matrix(0,500, t_end)

for(t in 1:t_end ){
  
  xx1<-matrix(0,500, newinf_true[t])
  
  for(i in 1:500) {
    
    for(j in 1:newinf_true[t])
    {
      xx1[i,j] <- (dist(rbind(c(x[i],y[i]),c(x[which(out_cov$inftime<(t+1) & out_cov$inftime!=0)[j]],y[which(out_cov$inftime<(t+1) & out_cov$inftime!=0)[j]]))))
    }
    ndist[i,t]<-sum(xx1[i,]) 
  }
  
}

#####
ind<-rep(1:500, each=(t_end-1))
t<-rep(2:t_end,500)
inf<-rep(out_cov$inftime, each=(t_end-1))
X<-rep(x,each=(t_end-1))
Y<-rep(y,each=(t_end-1))
dat1<-data.frame(ind,X,Y,t,inf)

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
dat2<-data.frame(individual=dat1$ind, time=dat1$t, Y=dat1$inf2, distance=dat1$dis1)

####
for(i in 1:(500*(t_end-1))){
  dat2$in1[i]<-inf[i]
}

dat3<-dat2[dat2$in1==0 | dat2$time<=dat2$in1,]
dat2<-dat3
dat2$ndistance=log(dat2$distance)

##############
m1<-glm(Y~ndistance, family=binomial(link = "logit"), data=dat2)
AIC1[1]<-summary(m1)$aic


#  Model 1 to 20

for(k in 1:20){
  
  #################Logistic
  
  pw<-k*0.5
  
  ndist<-matrix(0,500, t_end)
  
  for(t in 1:t_end ){
    
    xx1<-matrix(0,500, newinf_true[t])
    
    for(i in 1:500) {
      
      for(j in 1:newinf_true[t])
      {
        xx1[i,j] <- (dist(rbind(c(x[i],y[i]),c(x[which(out_cov$inftime<(t+1) & out_cov$inftime!=0)[j]],y[which(out_cov$inftime<(t+1) & out_cov$inftime!=0)[j]]))))^(-pw)
      }
      ndist[i,t]<-sum(xx1[i,]) 
    }
    
  }
  
  #####
  ind<-rep(1:500, each=(t_end-1))
  t<-rep(2:t_end,500)
  inf<-rep(out_cov$inftime, each=(t_end-1))
  X<-rep(x,each=(t_end-1))
  Y<-rep(y,each=(t_end-1))
  dat1<-data.frame(ind,X,Y,t,inf)
  
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
  dat2<-data.frame(individual=dat1$ind, time=dat1$t, Y=dat1$inf2, distance=dat1$dis1)
  
  
  ####Important change occure
  for(i in 1:(500*(t_end-1))){
    dat2$in1[i]<-inf[i]
  }
  
  dat3<-dat2[dat2$in1==0 | dat2$time<=dat2$in1,]
  
  dat2<-dat3

  dat2$ndistance=log(dat2$distance)
  ##############
  m1<-glm(Y~ndistance, family=binomial(link = "logit"), data=dat2)
  summary(m1)
  
  AIC1[k+1]<-summary(m1)$aic

}

#AIC1

data1<-data.frame(Power=seq(0,10,0.5),AIC1)
power1=data1$Power[which(data1$AIC1==min(data1$AIC1))]

Gpow[gg]<-power1
GAIC[gg]<-data1$AIC1[which(data1$AIC1==min(data1$AIC1))]


}

data2<-data.frame(Gpow,GAIC)
formattable(data2)

#######################

