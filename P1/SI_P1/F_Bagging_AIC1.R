
#install.packages("EpiILM")
library(EpiILM)
#install.packages("adaptMCMC")
require(adaptMCMC)
library(formattable)
library(tidyr)
xx<-c(3,6,9,1:19)
palpha<-c(-105.899, rep(0,19))
#Cov1
#bagging
#Use forward AIC method
sen<-numeric()
spe<-numeric()
acc<-numeric()
Epi<-numeric()

for(gg1 in 1:20){
 
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
  
  out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)
  
  #out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)
  
  #out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)
  
  t_end<-max(out_cov$inftime)
  SI<-rep(0)
  for (i in 1:t_end){
    SI[i]<-length(out_cov$inftime[out_cov$inftime>i])
  }
  #SI
  
  newinf_true<-rep(0)
  for (i in 1:t_end){
    newinf_true[i]<-length(out_cov$inftime[out_cov$inftime==i])
  }
  
  newinf_true<-cumsum(newinf_true)
  #newinf_true
  
  rt=100
  
  facov<-matrix(0,rt,4)
  esig<-matrix(0,rt,4)
  
  
  for(g in 1:rt){
    #rt=1
    fa1<-sample (c(1:20), size=4, replace =F)
    
    xnam1 <- paste("A", fa1, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
    AIC<-numeric()
    Co<-numeric()
    
    f1<-function(alpha){
      lam=0
      yy=-lam*abs(alpha)
      return(exp(yy))
    }
    
    
    
    ###Null variable
    
    
    p1<-function(xx){
      loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=xx[1], beta=xx[2],Sformula = NULL)
      palpha1<-dunif(xx[1], min=0, max=5, log=TRUE)
      pbeta<-dunif(xx[2], min=0, max=10, log=TRUE)
      return(-(loglikeILM+palpha1+pbeta))
    }
    
    ob<-optim(c(0.7,5),p1,control=list(maxit=1000000))
    
    loglike1<-epilike(out_cov,tmax=t_end, sus.par=ob$par[1], beta=ob$par[2],Sformula = NULL)
    
    aic1<-2*(2)-2*loglike1
    
    AIC[1]<-aic1
    Co[1]<-0
    
    ####1 variable selection
    aic1<-numeric()
    for(d1 in 1:4){
      fa<-fa1[d1]
      xnam1 <- paste("A", fa, sep="")
      fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
      
    
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
      AIC[d1+1]<-aic1[d1]
      Co[d1+1]<-fa1[d1]
    }
    
    
    ####2 variable selection(12,23,34)
    aic1<-numeric()
    for(d1 in 1:3){
      fa<-c(fa1[d1], fa1[d1+1])
      xnam1 <- paste("A", fa, sep="")
      fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
      
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
    
    AIC[6]<-aic1[1]
    AIC[7]<-aic1[2]
    AIC[8]<-aic1[3]
    
    
    ####2 variable selection(13,24)
    aic1<-numeric()
    
    for(d1 in 1:2){
      fa<-c(fa1[d1], fa1[d1+2])
      xnam1 <- paste("A", fa, sep="")
      fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
      
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
    
    AIC[9]<-aic1[1]
    AIC[10]<-aic1[2]
    
    
    #2 variable selection(14)
    
    fa<-c(fa1[1], fa1[4])
    xnam1 <- paste("A", fa, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
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
    AIC[11]<-2*(length(fa)+2)-2*loglike1
    
    ####3 variable selection(123,234)
    
    aic1<-numeric()
    
    for(d1 in 1:2){
      fa<-c(fa1[d1],fa1[d1+1],fa1[d1+2])
      xnam1 <- paste("A", fa, sep="")
      fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
      
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
    
    AIC[12]<-aic1[1]
    AIC[13]<-aic1[2]
    
    
    ####3 variable selection(124,134)
    
    aic1<-numeric()
    
    for(d1 in 1:2){
      fa<-c(fa1[1],fa1[d1+1],fa1[4])
      xnam1 <- paste("A", fa, sep="")
      fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
      
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
    AIC[14]<-aic1[1]
    AIC[15]<-aic1[2]
    
    
    ####4 variable selection(1234)
    
    fa<-c(fa1[1],fa1[2],fa1[3],fa1[4])
    xnam1 <- paste("A", fa, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
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
    AIC[16]<-2*(length(fa)+2)-2*loglike1
    
    ######
    ###
    
    sq2<-1:16
    seris<-c(0,1,2,3,4,12,23,34,13,24,14,123,234,124,134,1234)
    dat2<-data.frame(sq2, AIC,seris)
    seach<-dat2$seris[which.min(dat2$AIC)]
    dat2
    #g=1
    
    x1<-0
    x2<-0
    x3<-0
    x4<-0
    
    if (seach==0){
      x1=0
    }else if (seach==1){
      x1=fa1[1]
    }else if (seach==2){
      x1=fa1[2]
    }else if (seach==3){
      x1=fa1[3]
    }else if (seach==4){
      x1=fa1[4]
    }else if (seach==12){
      x1=fa1[1] 
      x2=fa1[2]
    }else if (seach==23){
      x1=fa1[2]
      x2=fa1[3] 
    }else if (seach==34){
      x1=fa1[3] 
      x2=fa1[4]
    }else if (seach==13){
      x1=fa1[1] 
      x2=fa1[3]
    }else if (seach==24){
      x1=fa1[2] 
      x2=fa1[4]
    }else if (seach==14){
      x1=fa1[1] 
      x2=fa1[4]
    }else if (seach==123){
      x1=fa1[1] 
      x2=fa1[2] 
      x3=fa1[3]
    }else if (seach==234){
      x1=fa1[2] 
      x2=fa1[3] 
      x3=fa1[4]
    }else if (seach==124){
      x1=fa1[1] 
      x2=fa1[2] 
      x3=fa1[4]
    }else if (seach==134){
      x1=fa1[1] 
      x2=fa1[3] 
      x3=fa1[4]
    }else {
      x1=fa1[1] 
      x2=fa1[2] 
      x3=fa1[3] 
      x4=fa1[4]
    }
    
    facov[g,]<-fa1
    esig[g,]<-c(x1,x2,x3,x4)
    
  }
  
  dat<-data.frame(facov,esig)
  
  ###########
  #####
  ###
  x1<-numeric()
  for(i in 1:20){
    x1[i]=0 
    for(mn in 1:100){
      if (dat$X1[mn]==i) {x1[i]=x1[i]+1} 
    }
  }
  x2<-numeric()
  for(i in 1:20){
    x2[i]=0 
    for(mn in 1:100){
      if (dat$X2[mn]==i) {x2[i]=x2[i]+1} 
    }
  }
  
  x3<-numeric()
  for(i in 1:20){
    x3[i]=0 
    for(mn in 1:100){
      if (dat$X3[mn]==i) {x3[i]=x3[i]+1} 
    }
  }
  
  x4<-numeric()
  for(i in 1:20){
    x4[i]=0 
    for(mn in 1:100){
      if (dat$X4[mn]==i) {x4[i]=x4[i]+1} 
    }
  }
  
  tocov<-numeric()
  tocov<-x1+x2+x3+x4
  ######
  
  X2.1<-dat$X1.1
  sx1<-numeric()
  for(i in 1:20){
    sx1[i]=0 
    for(mn in 1:length(X2.1)){
      if (X2.1[mn]==i) {sx1[i]=sx1[i]+1} 
    }
  }
  
  #sx1
  
  
  X3.1<-dat$X2.1
  sx2<-numeric()
  
  for(i in 1:20){
    sx2[i]=0 
    for(mn in 1:length(X3.1)){
      if (X3.1[mn]==i) {sx2[i]=sx2[i]+1} 
    }
  }
  
  #sx2
  
  X4.1<-dat$X3.1
  sx3<-numeric()
  for(i in 1:20){
    sx3[i]=0 
    for(mn in 1:length(X4.1)){
      if (X4.1[mn]==i) {sx3[i]=sx3[i]+1} 
    }
  }
  
  #sx3
  
  sx4<-numeric()
  X5.1<-dat$X4.1
  for(i in 1:20){
    sx4[i]=0 
    for(mn in 1:length(X5.1)){
      if (X5.1[mn]==i) {sx4[i]=sx4[i]+1} 
    }
  }
  
  #sx4
  
  sicov<-numeric()
  sicov<-sx1+sx2+sx3+sx4
  prop1<-sicov/tocov
  
  
  #For cov1
  sen[gg1]<-sum(prop1[1]>=0.5)/1
  spe[gg1]<-sum(prop1[2:20]<0.5)/19
  acc[gg1]<-(sum(prop1[1]>=0.5)+sum(prop1[2:20]<0.5))/20
  Epi[gg1]<-gg1+10
}

datT<-data.frame(Epi,sen,spe,acc)
library(formattable)
formattable(datT)
datT

