################################################
### Illustration Recursion Mixture on R ###########
################################################

rm(list = setdiff(ls(), lsf.str())) #this command removes evertything except functions

source("functions4R.R")
ptm <- proc.time()
#set problem parameter
tau=1
n=2500 #CHECK: SOMETIMES I SAMPLE MORE TO HAVE A FINER GRID
###CHECK:ACTUAL SAMPLE IS PICKED BELOW HERE
betaA=0.5
sigma2A=0.5
y=rep(0,n)
pop=c(1,2,3,4,5)
x=sample(pop,n,replace=TRUE)
#w=1-abs(x-mean(pop))/(mean(pop))
#w=1-(1/3+x-min(x))/(2*(max(x)-min(x)))
#w=1-(1/3+x-min(x))/(1+x-min(x))



#generate y
for (i in 1:n){
    y[i]=rnorm(1,mean=betaA*x[i],sd=sqrt(sigma2A*x[i]))  
}



#Create the grid of point in which to evaluate the distribution (not needed I guess)

ytest=y[order(y)]
xtest=pop[order(pop)]
#wtest=1-abs(xtest-mean(xtest))/(mean(xtest))
#wtest=1-(1/3+xtest-min(xtest))/(1+xtest-min(xtest))
wtest=rep(0.5,length(xtest))





#prior
mup<-1
sigmaP<-2
cdn<-pnorm(ytest,mean=mup,sd=sqrt(sigmaP))
cdn<-matrix(rep(cdn,length(xtest)),ncol=length(xtest))
cdn0<- cdn #save the prior in order to use it for plot 

size<-dim(cdn)
###### Recursive procedure ######

#simple version 

#First we do the Student T copula
cdnOR<-matrix(0,size[1],size[2])

for (i in 1:length(xtest)){
  neword<-order(abs(x-xtest[i]),decreasing=TRUE)
  x0<-x[neword]
  y0<-y[neword]
  x1<-x0[x0!=xtest[i]]
  y1<-y0[x0!=xtest[i]]
  x2<-x0[x0==xtest[i]]
  y2<-y0[x0==xtest[i]]
  fitORinter<-recursive(y1,x1,length(x1),cdn,xtest,"Student-T")
  #monte carlo
  mc=10
  cdnMCOR<-matrix(0,size[1],size[2])
  for (t in 1:mc){
    perm=sample(seq(1,length(y2),1),replace=FALSE)
    y2=y2[perm]
    x2=x2[perm]
    fitORMC<-recursive(y2,x2,length(y2),fitORinter$cdn,xtest,"Student-T")
    cdnMCOR<-cdnMCOR+fitORMC$cdn
  }
  cdnOR[,i]<-cdnMCOR[,i]/mc
}


#The the Gaussian copula 
cdnOR_G<-matrix(0,size[1],size[2])

for (i in 1:length(xtest)){
  neword<-order(abs(x-xtest[i]),decreasing=TRUE)
  x0<-x[neword]
  y0<-y[neword]
  x1<-x0[x0!=xtest[i]]
  y1<-y0[x0!=xtest[i]]
  x2<-x0[x0==xtest[i]]
  y2<-y0[x0==xtest[i]]
  fitORinter<-recursive(y1,x1,length(x1),cdn,xtest,"Gaussian")
  #monte carlo
  mc=10
  cdnMCOR_G<-matrix(0,size[1],size[2])
  for (t in 1:mc){
    perm=sample(seq(1,length(y2),1),replace=FALSE)
    y2=y2[perm]
    x2=x2[perm]
    fitORMC_G<-recursive(y2,x2,length(y2),fitORinter$cdn,xtest,"Gaussian")
    cdnMCOR_G<-cdnMCOR_G+fitORMC_G$cdn
  }
  cdnOR_G[,i]<-cdnMCOR_G[,i]/mc
}
######################################
######### Plot CDF ##################
######################################

ptrue <- function(grid, x) {
  pnorm(grid,mean=betaA*x,sd=sqrt(sigma2A*x))
}
par(mfrow=c(2,2))

dd=1
xx=ecdf(y[x==xtest[dd]])
plot(ytest,cdnOR[,dd],lwd=3,type="l",lty=1,
     xlab="y",ylab="Cdf",ylim=c(0,1),main="x=1",cex.lab=1.3,cex.main=1.3)
lines(ytest,cdnOR_G[,dd],lwd=3,
      type="l",lty=1,col="red")
lines(ytest,ptrue(ytest,xtest[dd]),lwd=3,
      type="l",lty=1,col="blue")
lines(ytest,cdn0[,dd],lwd=3,type="l",lty=1,
      col="green")
legend(1.5,0.6, c("True","RR Stud-T","RR Gauss","Prior"),lty=c(1,1,1,1),
       col=c("Blue","Black","Red","Green"),bty="n",cex=1.2,y.intersp=0.2)


dd=5
xx=ecdf(y[x==xtest[dd]])
plot(ytest,cdnOR[,dd],lwd=3,type="l",lty=1,
     xlab="y",ylab="Cdf",ylim=c(0,1),main="x=5",cex.lab=1.3,cex.main=1.3)
lines(ytest,cdnOR_G[,dd],lwd=3,
      type="l",lty=1,col="red")
lines(ytest,ptrue(ytest,xtest[dd]),lwd=3,
      type="l",lty=1,col="blue")
lines(ytest,cdn0[,dd],lwd=3,type="l",lty=1,
      col="green")
legend(1.5,0.6, c("True","RR Stud-T","RR Gauss","Prior"),lty=c(1,1,1,1),
       col=c("Blue","Black","Red","Green"),bty="n",cex=1.2,y.intersp=0.2)

proc.time() - ptm







