library(EpiILM)
library(formattable)
#install.packages("adaptMCMC")
require(adaptMCMC)
library('R.utils')

###################Posterior predictive distribution
Gmean<-numeric()
Gsd<-numeric()
Gpci<-numeric()

for(gg in 1:3){
  x<-runif(500,0,10)
  y<-runif(500,0,10)
  
  #out_cov<-epidata(type="SI", n=500, tmax=10,x=x, y=y, sus.par=0.7, beta=4)
  #out_cov<-epidata(type="SI", n=500, tmax=10,x=x, y=y, sus.par=0.5, beta=3)
  #out_cov<-epidata(type="SI", n=500, tmax=20,x=x, y=y, sus.par=0.2, beta=4)
  out_cov<-epidata(type="SI", n=500, tmax=20,x=x, y=y, sus.par=0.9, beta=5)
  
  #out_cov
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
        xx1[i,j] <- (dist(rbind(c(x[i],y[i]),c(x[which(out_cov$inftime<(t+1) & out_cov$inftime!=0)[j]],y[which(out_cov$inftime<(t+1) & out_cov$inftime!=0)[j]]))))^(-5)
        
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
  
  library(MCMCpack)
  ##############
  logpriorfun <- function(beta){
    sum(dcauchy(beta, log=TRUE))
  }
  
  posterior <- MCMClogit(dat2$Y~dat2$distance,
                         user.prior.density=logpriorfun,
                         logfun=TRUE)
  
  
  
  nal1<-matrix(0,length(1001:10000),2) 
  nal1<-posterior[1001:10000,1:2]
  al1<-matrix(0,500,2)
  al1<-nal1[sample(nrow(nal1),500,replace=F),]
  
  ############
  re<-500
  yatg<-matrix(0, re,t_end)
  for(m in 1:re){
    yat1<-numeric()
    yat1[1]<-1
    t=2
    datt1<-dat2[dat2$time==t,]
    
    stme<-which(out_cov$inftime>(t-1) | out_cov$inftime==0)
    xx1<-matrix(0,length(stme), newinf_true[t-1])
    ndist<-numeric()
    
    for(i in 1:length(stme)) {
      
      for(j in 1:newinf_true[t-1])
      {
        xx1[i,j] <- (dist(rbind(c(x[stme[i]],y[stme[i]]),c(x[which(out_cov$inftime<t & out_cov$inftime!=0)[j]],y[which(out_cov$inftime<t & out_cov$inftime!=0)[j]]))))^(-5)
        
        
      }
      ndist[i]<-log(sum(xx1[i,], na.rm=TRUE)) 
    }
    
    
    po<-numeric()
    yhat<-numeric()
    
    for(i in 1:length(ndist)){
      po[i]<-exp(al1[m,1]+al1[m,2]*ndist[i])/(1+exp(al1[m,1]+al1[m,2]*ndist[i])) 
      yhat[i]<-rbern(1,po[i])
    }
    yat<-sum(yhat==1)
    yat1[t]<-yat
    infn<-which(out_cov$inftime==1)
    ########
    add.pos<-(infn)
    old.pos<-seq_len(length(yhat)+length(infn))
    old.pos<-old.pos[!old.pos %in% add.pos]
    
    new.vec<-rep(NA, length(old.pos))
    new.vec[add.pos]<-rep(1, length(infn))
    new.vec[old.pos]<-yhat
    nyat<-new.vec
    ####
    t=3
    datt1<-dat2[dat2$time==t,]
    
    stme<-which(nyat==0)
    itme<-which(nyat==1)
  
    xx1<-matrix(0,length(stme), length(itme))
    ndist<-numeric()
    
    for(i in 1:length(stme)) {
      
      for(j in 1:length(itme))
      {
        xx1[i,j] <- (dist(rbind(c(x[stme[i]],y[stme[i]]),c(x[itme[j]],y[itme[j]]))))^(-5)
        
        
      }
      ndist[i]<-log(sum(xx1[i,], na.rm=TRUE) )
    }
    
    
    po<-numeric()
    yhat<-numeric()
    
    for(i in 1:length(ndist)){
      po[i]<-exp(al1[m,1]+al1[m,2]*ndist[i])/(1+exp(al1[m,1]+al1[m,2]*ndist[i])) 
      yhat[i]<-rbern(1,po[i])
    }
    yat<-sum(yhat==1)
    yat1[t]<-yat
    
    infn<-which(nyat==1)
    
    ####
    add.pos<-infn
    old.pos<-seq_len(length(yhat)+length(infn))
    old.pos<-old.pos[!old.pos %in% add.pos]
    
    new.vec<-rep(NA, length(old.pos))
    new.vec[add.pos]<-rep(1, length(infn))
    new.vec[old.pos]<-yhat
    nyat<-new.vec
    inn2<-which(nyat==1)
    ####
    for(t in 4:t_end){
      datt1<-dat2[dat2$time==t,]
      
      itme<-inn2[inn2<=500]
      cou<-1:500
      stme<-setdiff(cou,itme)
      
      if(length(stme)>1){
        
        xx1<-matrix(0,length(stme), length(itme))
        ndist<-numeric()
        
        for(i in 1:length(stme)) {
          
          for(j in 1:length(itme))
          {
            xx1[i,j] <- (dist(rbind(c(x[stme[i]],y[stme[i]]),c(x[itme[j]],y[itme[j]]))))^(-5)
            
            
          }
          ndist[i]<-log(sum(xx1[i,], na.rm=TRUE) )
        }
        
        
        po<-numeric()
        yhat<-numeric()
        
        for(i in 1:length(ndist)){
          po[i]<-exp(al1[m,1]+al1[m,2]*ndist[i])/(1+exp(al1[m,1]+al1[m,2]*ndist[i])) 
          yhat[i]<-rbern(1,po[i])
        }
        yat<-sum(yhat==1)
        yat1[t]<-yat
        
        
        infn<-which(nyat==1)
        
        ####
        add.pos<-infn
        old.pos<-seq_len(length(yhat)+length(infn))
        old.pos<-old.pos[!old.pos %in% add.pos]
        
        new.vec<-rep(NA, length(old.pos))
        new.vec[add.pos]<-rep(1, length(infn))
        new.vec[old.pos]<-yhat
        nyat<-new.vec
        inn2<-which(nyat==1)
        
        
      }else{yat1[t]<-0}
    }
    yatg[m,]<-yat1 
  }
  
  
  mest<-numeric()
  for(i in 1:t_end){
    mest[i]<-mean(yatg[,i])
  }
  ##MSE
  b1<-numeric()
  for(j in 1:500){
    b1[j]<-sum((newinf_true2-yatg[j,])^2)/t_end
  }
  #b1
  Gmean[gg]<-mean(b1)
  Gsd[gg]<-sd(b1)
  
  lowerq<-rep(0)
  upperq<-rep(0)
  for(i in 1:t_end){
    lowerq[i]<-quantile(yatg[,i], 0.025)
    upperq[i]<-quantile(yatg[,i], 0.975)
  }
  
  
  ###Proportion of time CI contains true value 
  
  pci<-numeric()
  pci[1]<-1
  for(i in 2:t_end){
    pci[i]<-ifelse(newinf_true2[i]>lowerq[i] & newinf_true2[i]<upperq[i],1,0)
  }
  
  Gpci[gg]<-sum(pci)/t_end
  
}


Gmean
Gsd
Gpci 

data2<-data.frame(Gmean,Gsd,Gpci)
formattable(data2)




