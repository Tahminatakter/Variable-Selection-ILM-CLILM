#Backward stepwise AIC
#Mar 14,2022
#Use AIC for variable selection (Forward stepwise)
#install.packages("EpiILM")
library(EpiILM)
library(formattable)
#install.packages("adaptMCMC")
require(adaptMCMC)
###
xx<-c(3,6,9,1:19)
palpha<-c(-105.899, rep(0,19))

sen<-numeric()
spe<-numeric()
acc<-numeric()
Epi<-numeric()
rrep=20

for(gg1 in 1:rrep){

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
  out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)
  
  #out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)
  
  #out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)
  
  #out_cov
  t_end<-max(out_cov$inftime)
  
  SI<-rep(0)
  for (i in 1:t_end){
    SI[i]<-length(out_cov$inftime[out_cov$inftime>i])
  }
  
  
  newinf_true<-rep(0)
  for (i in 1:t_end){
    newinf_true[i]<-length(out_cov$inftime[out_cov$inftime==i])
  }
  
  newinf_true<-cumsum(newinf_true)
  
  
  AIC<-numeric()
  Co<-numeric()
  
  f1<-function(alpha){
    lam=0
    yy=-lam*abs(alpha)
    return(exp(yy))
  }
  
  aic1<-numeric()
  ####1 variable selection
  for(d1 in 1:20){
    fa<-d1
    xnam1 <- paste("A", fa, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
    
    da<-data.frame(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,A20)
    da1<-da[fa]
    
    
    p1<-function(xx){
      loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=xx[1:(length(fa)+1)], beta=xx[length(fa)+2],Sformula = fmla1)
      palpha1<-dunif(xx[1], min=0, max=5, log=TRUE)
      for(i in 2:(length(fa)+1)){
        palpha[i]<-log(f1(xx[i])*(xx[i]>=0))
      }
      palpha2<-sum(palpha)
      pbeta<-dunif(xx[length(fa)+2], min=0, max=10, log=TRUE)
      return(-(loglikeILM+palpha1+palpha2+pbeta))
    }
    
    ob<-optim(c(0.7,rep(0,length(fa)),5),p1,control=list(maxit=1000000))
    
    loglike1<-epilike(out_cov,tmax=t_end, sus.par=ob$par[1:(length(fa)+1)], beta=ob$par[length(fa)+2],Sformula = fmla1)
    aic1[d1]<-2*(length(fa)+2)-2*loglike1
    
  }
  
  
  AIC[1]<-min(aic1)
  Co[1]<-which(aic1==min(aic1, na.rm=T))
  re<-Co[1]
  
  
  
  ####
  ####2 variable selection
  aic1<-numeric()
  
  for(d1 in (1:20)[-re]){
    fa<-c(re,d1)
    xnam1 <- paste("A", fa, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
    
    da<-data.frame(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,A20)
    da1<-da[fa]
    
    
    p1<-function(xx){
      loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=xx[1:(length(fa)+1)], beta=xx[length(fa)+2],Sformula = fmla1)
      palpha1<-dunif(xx[1], min=0, max=5, log=TRUE)
      for(i in 2:(length(fa)+1)){
        palpha[i]<-log(f1(xx[i])*(xx[i]>=0))
      }
      palpha2<-sum(palpha)
      pbeta<-dunif(xx[length(fa)+2], min=0, max=10, log=TRUE)
      return(-(loglikeILM+palpha1+palpha2+pbeta))
    }
    
    ob<-optim(c(0.7,rep(0,length(fa)),5),p1,control=list(maxit=1000000))
    
    loglike1<-epilike(out_cov,tmax=t_end, sus.par=ob$par[1:(length(fa)+1)], beta=ob$par[length(fa)+2],Sformula = fmla1)
    aic1[d1]<-2*(length(fa)+2)-2*loglike1
    
    
  }
  
  AIC[2]<-min(aic1, na.rm=T)
  Co[2]<-which(aic1==min(aic1, na.rm=T))
  re1<-Co[2]
  
  
  
  ####3 variable selection
  for(e1 in 3:20){
    aic1<-numeric()
    re<-c(re,re1)
    for(d1 in (1:20)[-re]){
      fa<-c(re,d1)
      xnam1 <- paste("A", fa, sep="")
      fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
      
      
      da<-data.frame(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,A20)
      da1<-da[fa]
      
      
      p1<-function(xx){
        loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=xx[1:(length(fa)+1)], beta=xx[length(fa)+2],Sformula = fmla1)
        palpha1<-dunif(xx[1], min=0, max=5, log=TRUE)
        for(i in 2:(length(fa)+1)){
          palpha[i]<-log(f1(xx[i])*(xx[i]>=0))
        }
        palpha2<-sum(palpha)
        pbeta<-dunif(xx[length(fa)+2], min=0, max=10, log=TRUE)
        return(-(loglikeILM+palpha1+palpha2+pbeta))
      }
      
      ob<-optim(c(0.7,rep(0,length(fa)),5),p1,control=list(maxit=1000000))
      
      
      loglike1<-epilike(out_cov,tmax=t_end, sus.par=ob$par[1:(length(fa)+1)], beta=ob$par[length(fa)+2],Sformula = fmla1)
      aic1[d1]<-2*(length(fa)+2)-2*loglike1
      
    }
    
    
    AIC[e1]<-min(aic1, na.rm=T)
    #Co[e1]<-min(dat1$sq[dat1$aic1==AIC[e1]], na.rm=T)
    Co[e1]<-which(aic1==min(aic1, na.rm=T))
    re1<-Co[e1]
  }
  
  
  ######
  ####
  ###
  
  sq2<-1:20
  dat2<-data.frame(sq2, AIC,Co)
  dat2$sq2[which.min(dat2$AIC)]
  
  ######
  
  ph1<-dat2$sq[which.min(dat2$AIC)]
  covs<-numeric()
  for(i in 1:ph1){
    covs[i]<-dat2$Co[i]
  }
  
  
  Ncovs<-matrix(0,1,20)
  for(r in 1:20){
    for(i in 1:ph1){
      if (covs[i]==r) {Ncovs[1,r]=r} 
    }
  }
  
  
  ###Cov1
  
  x1=0 
  for(j in 1){
    if(Ncovs[1,j]==j){x1=x1+1}
  }
  
  
  x2=0 
  for(j in 2:20){
    if(Ncovs[1,j]==j){x2=x2+1}
  }
  
  sen[gg1]<-x1/1
  spe[gg1]<-(19-x2)/19
  acc[gg1]<-(x1+19-x2)/20
  Epi[gg1]<-gg1+10
  
  
}

#######
datT<-data.frame(Epi,sen,spe,acc)
formattable(datT)
datT
######


