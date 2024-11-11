library(EpiILM)
library(formattable)
#install.packages("adaptMCMC")
require(adaptMCMC)
###
###################Posterior predictive distribution (Figure)

x<-runif(500,0,10)
y<-runif(500,0,10)


#out_cov<-epidata(type="SI", n=500, tmax=10,x=x, y=y, sus.par=0.7, beta=4)
#out_cov<-epidata(type="SI", n=500, tmax=10,x=x, y=y, sus.par=0.5, beta=3)
#out_cov<-epidata(type="SI", n=500, tmax=20,x=x, y=y, sus.par=0.2, beta=4)
out_cov<-epidata(type="SI", n=500, tmax=20,x=x, y=y, sus.par=0.9, beta=5)

#out_cov
t_end<-max(out_cov$inftime)

########
#plot(out_cov, plottype = "spatial") 
#plot(out_cov, plottype = "curve", curvetype = "newinfect")

############
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

#################ILM (Method1)

p1<-function(xx){
  loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=xx[1], beta=xx[2],Sformula = NULL)
  palpha1<-dunif(xx[1], min=0, max=5, log=TRUE)
  pbeta<-dunif(xx[2], min=0, max=10, log=TRUE)
  return(-(loglikeILM+palpha1+pbeta))
}

ob<-optim(c(0.6,6),p1,control=list(maxit=1000000))
####

p<-function(xx){
  loglikeILM<-epilike(out_cov,tmax=t_end, sus.par=xx[1], beta=xx[2],Sformula = NULL)
  palpha1<-dunif(xx[1], min=0, max=5, log=TRUE)
  pbeta<-dunif(xx[2], min=0, max=10, log=TRUE)
  return(loglikeILM+palpha1+pbeta)
}

test<-MCMC (p, n=10000, init=c(0.6,6), scale=c(0.02,0.02), adapt=TRUE, acc.rate=0.5)
#plot(as.mcmc(test$samples), density=FALSE)
mcm1<-window(test$samples, start=1001)
summary(mcm1)


mpost<-mcm1[sample(nrow(mcm1),500, replace=F),]
ninf<-replace(out_cov$inftime,out_cov$inftime!=1,0)

out1<-vector(mode="list",length=500)
for(i in 1:500){
  dat<-epidata(type="SI", n=500, tmax=t_end,Sformula = NULL, sus.par=mpost[i,1], beta=mpost[i,2],x=x, y=y, inftime=ninf)
  out1[[i]]<-dat$inftime
}

#Number of new infection at each time point
infect1<-matrix(0, ncol=t_end, nrow=500)
for(j in 1:500){
  newinf<-rep(0)
  for(i in 1:t_end){
    newinf[i]<-length(out1[[j]][out1[[j]]==i])
  }
  infect1[j,]<-newinf
}

mest<-numeric()
for(i in 1:t_end){
  mest[i]<-mean(infect1[,i])
}

##MSE
b1<-numeric()
for(j in 1:500){
  b1[j]<-sum((newinf_true2-infect1[j,])^2)/t_end
}

newinf_true2
mean(b1)
sd(b1)


#Posterioir predictive realizations
lowerq<-rep(0)
upperq<-rep(0)
for(i in 1:t_end){
  lowerq[i]<-quantile(infect1[,i], 0.025)
  upperq[i]<-quantile(infect1[,i], 0.975)
}

time<-rep(1:t_end)
plot(time, newinf_true2, xlim=c(min(time), max(time)), ylim=c(0,350),ylab="No. of new infections", xlab="Time",main="(a) ILM", type="b", pch=20, lwd=2,xaxt="n",col="red")
axis(1, at=1:max(time))
par(new=TRUE)
for(i in 1:500){
  lines(time, infect1[i,], col="green",lwd=0.2)
}



par(new=TRUE)
plot(time, newinf_true2, xlim=c(min(time), max(time)), ylim=c(0,350),ylab="No. of new infections", xlab="Time", type="b", pch=20, lwd=2,xaxt="n",col="red")
axis(1, at=1:max(time))
lines(lowerq, type="b", pch=20)
lines(upperq, type="b", pch=20)
lines(mest,col="blue",  type="b",pch=20)
legend(1, 350, legend=c("Real data", "Simulated data"),
       col=c("red", "green"), lty=1, cex=0.8)


##########












