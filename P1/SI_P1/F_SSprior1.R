
#install.packages("EpiILM")
library(EpiILM)
#install.packages("adaptMCMC")
require(adaptMCMC)
library(formattable)
#library(tidyr)
xx<-c(3,6,9,1:19)
palpha<-c(-105.899, rep(0,19))


#####COV1
#install.packages("Rmisc")
#library(Rmisc)
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
  
  
  xnam <- paste("A", 1:20, sep="")
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  
  
  pp1<-function(xx){
    loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=c(xx[1], xx[2], xx[3], xx[4], xx[5], xx[6], xx[7], xx[8], xx[9], xx[10], xx[11], xx[12], xx[13], xx[14], xx[15], xx[16],xx[17], xx[18], xx[19], xx[20], xx[21]), beta=xx[22],Sformula = fmla)
    return(loglikeILM)
  }
  
  
  
  startvalue =c(0.7,rep(0.45,20),5)
  iterations=1000
  chain = array(dim = c(iterations+1,22))
  chain[1,] = startvalue
  
  
  chose10<-numeric()
  chosek0<-matrix(0,iterations, 20)
  
  for (i in 1:iterations){
    
    proposal=rbern(1,0.5)*runif(1, min=0, max=5)
    startvalue[2]<-proposal 
    new1= startvalue
    
    probab = exp(pp1(new1)-pp1(chain[i,]))
    probab1=min(1,probab)
    if (runif(1) < probab1){
      chain[i+1,] = new1
    }else{
      chain[i+1,] = chain[i,]
    }
    
    
    #Block 2:20
    
    for (k in 2:20){
      startvalue<-chain[i+1,]
      proposal=rbern(1,0.5)*runif(1, min=0, max=5)
      
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
    proposal=runif(1, min=0, max=6)
    startvalue[22]<-proposal 
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
  
  #############################################
  
  burnIn<-50
  chaink0<-matrix(0,iterations,22)
  
  for(i in 1:iterations){
    for (k in 1:22){
      if (chain[i,k]==0){
        chaink0[i,k]<-0
      }else{
        chaink0[i,k]<-1
      }
      
    }
  }
  
  
  
  chok0<-numeric()
  for(q in 1:22){
    chok0[q]<-sum(chaink0[-(1:burnIn),q])/length(chaink0[-(1:burnIn),q])
  }
  
  
  
  
  sen[gg1]<-sum(chok0[2]>=0.5)/1
  spe[gg1]<-sum(chok0[3:21]<0.5)/19
  acc[gg1]<-(sum(chok0[2]>=0.5)+sum(chok0[3:21]<0.5))/20
  Epi[gg1]<-gg1+10
  
}

datT<-data.frame(Epi,sen,spe,acc)
library(formattable)
formattable(datT)
datT

