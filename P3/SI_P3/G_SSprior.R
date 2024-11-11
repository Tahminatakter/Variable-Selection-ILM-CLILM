library(EpiILM)
library(formattable)
#install.packages("adaptMCMC")
require(adaptMCMC)
library('R.utils')
#install.packages("glmnet")
library('glmnet')


cut1<-c(0.5,0.2)
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


para<-matrix(0,length(B),22)
sen<-matrix(0,length(B),length(cut1))
spe<-matrix(0,length(B),length(cut1))
acc<-matrix(0,length(B),length(cut1))


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
#out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)
out_cov<-epidata(type="SI", n=500, tmax=50,x=x, y=y,Sformula = ~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+A17+A18+A19+A20, sus.par=c(0.7,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0,0,0,0,0,0,0,0,0,0,0,0,0), beta=5)

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
    #ndist[i,t]<-sum(xx1[i,]) 
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
library(MCMCpack)
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


# Fit logistic regression model
model <- glm(dat2$Y ~ XX, family = binomial)

# Extract coefficients
beta <- coef(model)
# Compute log-likelihood function
pp1 <- function(beta) {
  y<-dat2$Y
  X<-cbind(1,XX)
  eta <- X %*% beta
  p <- 1 / (1 + exp(-eta))
  log_likelihood <- sum(ifelse(y == 1, log(p), log(1 - p)))
  return(log_likelihood)
}

# Compute log-likelihood
log_likelihood <- pp1(beta)
#print(log_likelihood)


###########


startvalue =c(beta)
iterations=10000
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


burnIn<-5000


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

chok0

#ifelse(chok0<0.5,0,1)
#ifelse(chok0<0.2,0,1)

##########For cov7 (similarly calculate for Cov1 and Cov2)


para[gg,]<-chok0

#Sen, Spe, and acc:Cov7
for (g in 1: length(cut1)){
  sen[gg,g]<-sum(ifelse(para[gg,2:8]>=cut1[g],1,0))/7
  spe[gg,g]<-sum(ifelse(para[gg,9:21]<cut1[g],1,0))/13
  acc[gg,g]<-(sum(ifelse(para[gg,2:8]>=cut1[g],1,0))+sum(ifelse(para[gg,9:21]<cut1[g],1,0)))/20
  
}



}

data2<-data.frame(sen=round(sen[,1],3),spe=round(spe[,1],3),acc=round(acc[,1],3))
formattable(data2)

data3<-data.frame(sen=round(sen[,2],3),spe=round(spe[,2],3),acc=round(acc[,2],3))
formattable(data3)

data4<-data.frame(msen=mean(round(sen[,1],3)),mspe=mean(round(spe[,1],3)),macc=mean(round(acc[,1],3)),sdsen=sd(round(sen[,1],3)),sdspe=sd(round(spe[,1],3)),sdacc=sd(round(acc[,1],3)))
formattable(data4)

data5<-data.frame(msen=mean(round(sen[,2],3)),mspe=mean(round(spe[,2],3)),macc=mean(round(acc[,2],3)),sdsen=sd(round(sen[,2],3)),sdspe=sd(round(spe[,2],3)),sdacc=sd(round(acc[,2],3)))
formattable(data5)




