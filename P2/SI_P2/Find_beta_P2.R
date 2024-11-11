library(EpiILM)
library(formattable)
#install.packages("adaptMCMC")
require(adaptMCMC)
library(MCMCpack)

###################Find Beta
Gpow<-numeric()
GAIC<-numeric()


for (gg in 1:20){

AIC1<-numeric()

x<-runif(500,0,10)
y<-runif(500,0,10)

#out_cov<-epidata(type="SI", n=500, tmax=10,x=x, y=y, sus.par=0.7, beta=4)
#out_cov<-epidata(type="SI", n=500, tmax=10,x=x, y=y, sus.par=0.5, beta=3)
#out_cov<-epidata(type="SI", n=500, tmax=20,x=x, y=y, sus.par=0.2, beta=4)
out_cov<-epidata(type="SI", n=500, tmax=20,x=x, y=y, sus.par=0.9, beta=5)

t_end<-max(out_cov$inftime)

#plot(out_cov, plottype = "spatial")
#plot(out_cov, plottype = "curve", curvetype = "newinfect")

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

#################Logistic

#For all t
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

###
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
 
  pw<-k*0.5
  
  #For all t
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
  
  
  ##
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
#power1
#formattable(data1)


Gpow[gg]<-power1
GAIC[gg]<-data1$AIC1[which(data1$AIC1==min(data1$AIC1))]

}

data2<-data.frame(Gpow,GAIC)
formattable(data2)

#######################



