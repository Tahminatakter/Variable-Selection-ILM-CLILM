#Feb 28,2022
#Find optimal lambda and cutoff using AIC
##COV1
#install.packages("EpiILM")
library(EpiILM)
#install.packages("adaptMCMC")
require(adaptMCMC)
library(formattable)
library(tidyr)

xx<-c(3,6,9,1:19)
palpha<-c(-105.899, rep(0,19))
###
lambda_Avg_1<-numeric()
cut1<-numeric()
nz2<-numeric()
nk2<-numeric()
z2<-numeric()
k2<-numeric()
nac1<-numeric()

lambda_Avg_2<-numeric() 
cut2<-numeric()
nh1<-numeric()
nng2<-numeric()
nh2<-numeric()
ng2<-numeric()
nh3<-numeric()
nacu1<-numeric()

sz=20
for(g in 1:sz){

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
  
  f1<-function(alpha){
    lam=0
    yy=-lam*abs(alpha)
    return(exp(yy))
  }
  
  xnam <- paste("A", 1:20, sep="")
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  
  p<-function(xx){
    loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=c(xx[1], xx[2], xx[3], xx[4], xx[5], xx[6], xx[7], xx[8], xx[9], xx[10], xx[11], xx[12], xx[13], xx[14], xx[15], xx[16],xx[17], xx[18], xx[19], xx[20], xx[21]), beta=xx[22],Sformula = fmla)
    palpha1<-dunif(xx[1], min=0, max=5, log=TRUE)
    for(j in 2:21){
      palpha[j]<-log(f1(xx[j])*(xx[j]>=0))
    }
    palpha2<-sum(palpha)
    pbeta<-dunif(xx[22], min=0, max=10, log=TRUE)
    return(-(loglikeILM+palpha1+palpha2+pbeta))
  }
  
  ob<-optim(c(0.7,rep(0,20),5),p,control=list(maxit=1000000))
  
  #Estimate parameter 
  para<-ob$par
  
  #Optime cutoff and lambda at first stage lasso
  cu1<-c(seq(0.01,0.07, 0.01), 0.1)
  loglike1<-matrix(0,length(cu1),21)
  aic1<-matrix(0,length(cu1),21)
  pn1<-numeric()
  
  for(u1 in 1:length(cu1)){
    
    pn1[u1]<-sum(ifelse(para[2:21]>=cu1[u1],1,0))
    
    fa<-which(para[2:21]>=cu1[u1],1,0)
    xnam1 <- paste("A", fa, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
    da<-data.frame(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,A20)
    da1<-da[fa]
    
    for(l in 1:21){
      f1<-function(alpha){
        lam=(l-1)*0.5
        yy=-lam*abs(alpha)
        return(exp(yy))
      }
      ####2nd stage confusion
      
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
      loglike1[u1,l]<-epilike(out_cov,tmax=t_end, sus.par=ob$par[1:(length(fa)+1)], beta=ob$par[length(fa)+2],Sformula = fmla1)
      aic1[u1,l]<-2*(pn1[u1]+2)-2*loglike1[u1,l]
    }
  }
  
  
  
  
  ##
  
  cu1<-c(seq(0.01,0.07,0.01),0.1)
  dat1_2<-data.frame(cu1, aic1)
  
  
  library(tidyr)
  da3 <- gather(dat1_2, X, aic1, X1:X21, factor_key=TRUE)
  
  da3$Lambda_2a<-rep(seq(0,10,0.5),each=8)
  
  
  lambda_Avg_1[g]<-da3$Lambda_2a[which.min(da3$aic1)]
  #lambda_Avg_1
  
  cut1[g]<-da3$cu1[which.min(da3$aic1)]
  #cut1
  
  #####Calculate Sen, Spe, Accu based on optimal lambda and cutoff at first stage
  
  f1<-function(alpha){
    lam=lambda_Avg_1[g]
    yy=-lam*abs(alpha)
    return(exp(yy))
  }
  
  
  xnam <- paste("A", 1:20, sep="")
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  
  p<-function(xx){
    loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=c(xx[1], xx[2], xx[3], xx[4], xx[5], xx[6], xx[7], xx[8], xx[9], xx[10], xx[11], xx[12], xx[13], xx[14], xx[15], xx[16],xx[17], xx[18], xx[19], xx[20], xx[21]), beta=xx[22],Sformula = fmla)
    palpha1<-dunif(xx[1], min=0, max=5, log=TRUE)
    for(j in 2:21){
      palpha[j]<-log(f1(xx[j])*(xx[j]>=0))
    }
    palpha2<-sum(palpha)
    pbeta<-dunif(xx[22], min=0, max=10, log=TRUE)
    return(-(loglikeILM+palpha1+palpha2+pbeta))
  }
  
  ob<-optim(c(0.7,rep(0,20),5),p,control=list(maxit=1000000))
  
  para<- ob$par
  
  
  
  nz2[g]<-sum(ifelse(para[2]>=cut1[g],1,0))
  nk2[g]<-nz2[g]/1
  z2[g]<-sum(ifelse(para[3:21]<cut1[g],1,0))
  k2[g]<-z2[g]/19
  nac1[g]<-(sum(ifelse(para[2]>=cut1[g],1,0))+sum(ifelse(para[3:21]<cut1[g],1,0)))/20
  
  
  
  ####Find optimal cutoff and lambda at second Stage lasso
  
  
  cu1<-c(seq(0.01,0.07, 0.01), 0.1)
  loglike2<-matrix(0,length(cu1),21)
  aic2<-matrix(0,length(cu1),21)
  pn2<-numeric()
  
  
  for(u1 in 1:length(cu1)){
    pn2[u1]<-sum(ifelse(para[2:21]>=cu1[u1],1,0))
    
    fa<-which(para[2:21]>=cu1[u1],1,0)
    xnam1 <- paste("A", fa, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
    
    da<-data.frame(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,A20)
    da1<-da[fa]
    #head(da1)
    
    ###
    
    for(l in 1:21){
      f1<-function(alpha){
        lam=(l-1)*0.5
        yy=-lam*abs(alpha)
        return(exp(yy))
      }
      
      
      ####2nd stage 
      
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
      loglike2[u1,l]<-epilike(out_cov,tmax=t_end, sus.par=ob$par[1:(length(fa)+1)], beta=ob$par[length(fa)+2],Sformula = fmla1)
      aic2[u1,l]<-2*(pn2[u1]+2)-2*loglike2[u1,l]
    }
  }
  
  cu1<-c(seq(0.01,0.07,0.01),0.1)
  dat1_2<-data.frame(cu1, aic2)
  
  
  library(tidyr)
  da3 <- gather(dat1_2, X, aic2, X1:X21, factor_key=TRUE)
  
  da3$Lambda_2a<-rep(seq(0,10,0.5),each=8)
  
  
  lambda_Avg_2[g]<-da3$Lambda_2a[which.min(da3$aic2)]
  #lambda_Avg_2
  
  cut2[g]<-da3$cu1[which.min(da3$aic2)]
  #cut2  
  
  fa1<-which(para[2:21]>=cut1[g],1,0)
  
  fa<-fa1
  
  ###Calculate sen, spec, accu based on optimal lambda and cutoff at second stage
  if(nz2[g]==0 && z2[g]==19){
    nng2[g]<-0/1
    ng2[g]<-z2[g]/19
    nacu1[g]<-z2[g]/20
  }else{
    f1<-function(alpha){
      lam=lambda_Avg_2[g]
      yy=-lam*abs(alpha)
      return(exp(yy))
    }
    
    
    xnam1 <- paste("A", fa, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
    ####2nd stage confusion
    
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
    
    para1<- ob$par
    
    
    if(nz2[g]==1 && z2[g]<19){
      nh1[g]<-sum(ifelse(para1[2]>=cut2[g],1,0))
      nng2[g]<-sum(ifelse(para1[2]>=cut2[g],1,0))/1
      nh2[g]<-(z2[g]+sum(ifelse(para1[3:(length(fa)+1)]<cut2[g],1,0)))
      ng2[g]<-(z2[g]+sum(ifelse(para1[3:(length(fa)+1)]<cut2[g],1,0)))/19
      nh3[g]<-(sum(ifelse(para1[2]>=cut2[g],1,0))+z2[g]+sum(ifelse(para1[3:(length(fa)+1)]<cut2[g],1,0)))
      nacu1[g]<-(sum(ifelse(para1[2]>=cut2[g],1,0))+z2[g]+sum(ifelse(para1[3:(length(fa)+1)]<cut2[g],1,0)))/20
    }else if(nz2[g]==1 && z2[g]==19){
      nh1[g]<-sum(ifelse(para1[2]>=cut2[g],1,0))
      nng2[g]<-sum(ifelse(para1[2]>=cut2[g],1,0))/1
      nh2[g]<-z2[g]
      ng2[g]<-z2[g]/19
      nh3[g]<-(sum(ifelse(para1[2]>=cut2[g],1,0))+z2[g])
      nacu1[g]<-(sum(ifelse(para1[2]>=cut2[g],1,0))+z2[g])/20
      
    }else{ 
      nh1[g]<-0
      nng2[g]<-0/1
      nh2[g]<-(z2[g]+sum(ifelse(para1[2:(length(fa)+1)]<cut2[g],1,0)))
      ng2[g]<-(z2[g]+sum(ifelse(para1[2:(length(fa)+1)]<cut2[g],1,0)))/19
      nh3[g]<-(z2[g]+sum(ifelse(para1[2:(length(fa)+1)]<cut2[g],1,0)))
      nacu1[g]<-(z2[g]+sum(ifelse(para1[2:(length(fa)+1)]<cut2[g],1,0)))/20
      
    } 
    
  }
}



##
Epidemic<-11:30

dat<-data.frame(Epidemic,lambda_Avg_1,cut1, round(nk2,3), round(k2,3), nac1, lambda_Avg_2,cut2, round(nng2,3), round(ng2,3), nacu1)
library(formattable)
formattable(dat)
dat

