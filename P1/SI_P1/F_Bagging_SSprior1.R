#install.packages("EpiILM")
library(EpiILM)
#install.packages("adaptMCMC")
require(adaptMCMC)
library(formattable)
#library(tidyr)
xx<-c(3,6,9,1:19)
palpha<-c(-105.899, rep(0,19))
#Cov1
#bagging
#Use spike and slab method
#Sample 4 covariates
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
  
  re=100
  cov1<-matrix(0,re,5)
  facov<-matrix(0,re,5)
  
  for(g in 1:re){
    
    fa<-sample (c(1:20), size=4, replace =F)
    
    xnam1 <- paste("A", fa, sep="")
    fmla1 <- as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
    
    
    #da<-data.frame(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,A20)
    #da1<-da[fa]
    
    
    #xx<-c(3,6,9,1:19)
    #palpha<-c(-105.899, rep(0,19))
    
    pp1<-function(xx){
      loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=xx[1:(length(fa)+1)], beta=xx[length(fa)+2],Sformula = fmla1)
      return(loglikeILM)
    }
    
    
    startvalue =c(0.7,rep(0.45,4),5)
    iterations=500
    chain = array(dim = c(iterations+1,6))
    chain[1,] = startvalue
    skp=0.5
    
    for (i in 1:iterations){
      #i=1
      #Block 1
      #curr = startvalue[2]
      proposal=rbern(1,skp)*runif(1, min=0, max=5)
      startvalue[2]<-proposal 
      new1= startvalue
      
      #probab = (prior(proposal)*exp(pp1(new1)))/(prior(curr)*exp(pp1(chain[i,])))
      probab = exp(pp1(new1)-pp1(chain[i,]))
      probab1=min(1,probab)
      if (runif(1) < probab1){
        chain[i+1,] = new1
      }else{
        chain[i+1,] = chain[i,]
      }
      
      
      #Block 2:20
      
      for (k in 2:4){
        startvalue<-chain[i+1,]
        #curr = startvalue[k+1]
        proposal=rbern(1,skp)*runif(1, min=0, max=5)
        
        startvalue[k+1]<-proposal 
        new1= startvalue
        probab = exp(pp1(new1)-pp1(chain[i+1,]))
        probab1=min(1,probab)
        if (runif(1) < probab1){
          chain[i+1,] = new1
        }else{
          chain[i+1,] = chain[i+1,]
        }
        
      }
      
      
      
      
      #Block 0
      startvalue<-chain[i+1,]
      #curr = startvalue[k+1]
      proposal=runif(1, min=0, max=5)
      startvalue[1]<-proposal 
      new1= startvalue
      probab = exp(pp1(new1)-pp1(chain[i+1,]))
      probab1=min(1,probab)
      if (runif(1) < probab1){
        chain[i+1,] = new1
      }else{
        chain[i+1,] = chain[i+1,]
      }
      
      #Block 22
      startvalue<-chain[i+1,]
      #curr = startvalue[k+1]
      proposal=runif(1, min=0, max=6)
      startvalue[6]<-proposal 
      new1= startvalue
      probab = exp(pp1(new1)-pp1(chain[i+1,]))
      probab1=min(1,probab)
      if (runif(1) < probab1){
        chain[i+1,] = new1
      }else{
        chain[i+1,] = chain[i+1,]
      }
      
      
      startvalue<-chain[i+1,]
    }
    
    burnIn<-50
    chaink0<-matrix(0,iterations,6)
    
    for(i in 1:iterations){
      for (k in 1:6){
        if (chain[i,k]==0){
          chaink0[i,k]<-0
        }else{
          chaink0[i,k]<-1
        }
        
      }
    }
    
    
    
    chok0<-numeric()
    for(q in 1:6){
      chok0[q]<-sum(chaink0[-(1:burnIn),q])/length(chaink0[-(1:burnIn),q])
    }
    
    
    
    sig<-numeric()
    for(r in 2:5){
      if(chok0[r]<0.5){
        sig[r]<-0
      }else{
        sig[r]<-fa[(r-1)]
      }
      cov1[g,r]<-sig[r]
      facov[g,r]<-fa[(r-1)]
    }
    
    
  }
  
  #####
  ###
  ###
  ###########3:23pm 
  #####
  ###
  
  
  min(cov1[5,2:5])
  max(cov1[5,2:5])
  
  min1<-numeric()
  max1<-numeric()
  
  for (i in 1:100){
    min1[i]<-min(cov1[i,2:5])
    max1[i]<-max(cov1[i,2:5])
  }
  
  
  ###
  min12<-numeric()
  max12<-numeric()
  
  for (i in 1:100){
    min12[i]<-min(facov[i,2:5])
    max12[i]<-max(facov[i,2:5])
  }
  
  
  
  ####
  
  
  min(facov[6,2:5])
  max(facov[6,2:5])
  
  sumcov<-numeric()
  for(i in 1:100){
    sumcov[i]=sum(cov1[i,])
  }
  
  datt<-data.frame(facov[,2:5],min12,max12, cov1[,2:5], min1,max1)
  dat<-datt
  
  #####
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
  
  X2.1<-dat$X1.1[!is.na(dat$X1.1)]
  sx1<-numeric()
  if(length(X2.1)==0){
    for(i in 1:20){sx1[i]=0}
  }else{
    for(i in 1:20){
      sx1[i]=0 
      for(mn in 1:length(X2.1)){
        if (X2.1[mn]==i) {sx1[i]=sx1[i]+1} 
      }
    }
  }
  #sx1
  
  
  X3.1<-dat$X2.1[!is.na(dat$X2.1)]
  sx2<-numeric()
  if(length(X3.1)==0){
    for(i in 1:20){sx2[i]=0}
  }else{
    for(i in 1:20){
      sx2[i]=0 
      for(mn in 1:length(X3.1)){
        if (X3.1[mn]==i) {sx2[i]=sx2[i]+1} 
      }
    }
  }
  #sx2
  
  X4.1<-dat$X3.1[!is.na(dat$X3.1)]
  sx3<-numeric()
  if(length(X4.1)==0){
    for(i in 1:20){sx3[i]=0}
  }else{
    for(i in 1:20){
      sx3[i]=0 
      for(mn in 1:length(X4.1)){
        if (X4.1[mn]==i) {sx3[i]=sx3[i]+1} 
      }
    }
  }
  #sx3
  
  sx4<-numeric()
  X5.1<-dat$X4.1[!is.na(dat$X4.1)]
  if(length(X5.1)==0){
    for(i in 1:20){sx4[i]=0}
  }else{for(i in 1:20){
    sx4[i]=0 
    for(mn in 1:length(X5.1)){
      if (X5.1[mn]==i) {sx4[i]=sx4[i]+1} 
    }
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

