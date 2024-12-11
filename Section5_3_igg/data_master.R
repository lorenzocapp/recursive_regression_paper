################################################
###  Recursion Mixture o R ###########
################################################

rm(list = setdiff(ls(), lsf.str())) #this command removes evertything except functions
ptm <- proc.time()
#set problem parameter
library(DPpackage)
library(np)
tau=1
###CHECK:ACTUAL SAMPLE IS PICKED BELOW HERE
dataigg=data(igg)
y=igg$igg
x=igg$age
x[x<=1]=0.5
x[x>1 & x <=2]=1.5
x[x>2 & x <=3]=2.5
x[x>3 & x <=4]=3.5
x[x>4 & x <=5]=4.5
x[x>5 ]=5.5
n=length(y)
pop=unique(x)


perm=sample(seq(1,n,1),replace=FALSE)
x=x[perm]
y=y[perm]





#Create the grid of point in which to evaluate the distribution (not needed I guess)
ygr=unique(y)
if (length(ygr)<=400){
  add=400-length(ygr)-1
  miny=0.2*min(y)
  maxy=1.1*max(y)
  int=(maxy-miny)/add
  yadd=seq(miny,maxy,int)
} else {
  yadd=numeric(0)
}


ytest=sort(c(ygr,yadd))
xtest=pop[order(pop)]


#################################
###### Recursive procedure ######
#################################


#prior
mup<-0
sigmaP<-4
cdn<-pnorm(ytest,mean=mup,sd=sqrt(sigmaP))
#cdn=pgamma(ytest,shape=3,rate=2)
cdn<-matrix(rep(cdn,length(xtest)),ncol=length(xtest))
cdn0<- cdn #save the prior in order to use it for plot 

size<-dim(cdn)

#simple version 

#monte carlo
ptm <- proc.time()
mc=50

cdnMC2<-matrix(0,size[1],size[2])

for (t in 1:mc){
  cdn<-cdn0
  perm=sample(seq(1,n,1),replace=FALSE)
  y=y[perm]
  x=x[perm]
  fitMC2<-recursive(y,x,n,cdn,xtest,"StudentT",8)
  cdnMC2<-cdnMC2+fitMC2$cdn
}
cdnMC2<-cdnMC2/mc
proc.time() - ptm

#monte carlo
ptm <- proc.time()
mc=20

cdnMC1<-matrix(0,size[1],size[2])

for (t in 1:mc){
  cdn<-cdn0
  perm=sample(seq(1,n,1),replace=FALSE)
  y=y[perm]
  x=x[perm]
  fitMC1<-recursive(y,x,n,cdn,xtest,"Gaussian",8)
  cdnMC1<-cdnMC1+fitMC1$cdn
}
cdnMC1<-cdnMC1/mc
proc.time() - ptm

print(sum(fit$saveP))

#monte carlo
ptm <- proc.time()
mc=20

cdnMC<-matrix(0,size[1],size[2])

for (t in 1:mc){
  cdn<-cdn0
  perm=sample(seq(1,n,1),replace=FALSE)
  y=y[perm]
  x=x[perm]
  fitMC<-recursive(y,x,n,cdn,xtest,"StudentT",8)
  cdnMC<-cdnMC+fitMC$cdn
}
cdnMC<-cdnMC/mc
proc.time() - ptm

#################################
###### Bayesian procedure ######
#################################


#prior
# W <- x
# S0 <- 5
# m0 <- 4
# 
# prior <- list(a0 = 8, b0 = 1, m0 = m0,
#               S0 = S0, tau1 = 6.01,
#               taus1 = 6.01, taus2 = 2.01,
#               nu = 1, psiinv = 0.2)

W <- x
S0 <- 2
m0 <- 0

prior <- list(a0 = 1, b0 = 1, m0 = m0,
              S0 = S0, tau1 = 1.01,
              taus1 = 1.01, taus2 = 2.01,
              nu = 1, psiinv = 0.2)

#grid in which to predict
#Wpred<- matrix(rep(c(1,xtest),n),nrow=n,ncol=length(c(1,xtest)),byrow=TRUE)
Wpred<-t(xtest)
Wpred<-t(Wpred)



# Initial state
state <- NULL


# MCMC parameters
nburn <- 5000
nsave <- 5000
nskip <- 3
ndisplay <- 100
mcmc <- list(nburn=nburn,
             nsave=nsave,
             nskip=nskip,
             ndisplay=ndisplay)


# Fit the model
fitBay <- LDDPdensity(formula=y~0+W,zpred=Wpred, prior=prior,mcmc=mcmc,
                    compute.band = FALSE,state=NULL,status=TRUE,
                    grid=ytest)



#MCMC diagnostics
plot(fitBay, hpd = FALSE)




#compute the CDF through trapezoidal rule. 
trap_cdf<- function(f,ytest){
  diffy<-ytest[2:length(ytest)]-ytest[1:length(ytest)-1]
  sumbase<-f[1:length(f)-1]+f[2:length(f)]
  prod<-diffy*sumbase/2
  trap_cdf<-c(0,cumsum(prod))
  trap_cdf<-trap_cdf/trap_cdf[length(trap_cdf)] #in this step I am normalizing the stuff
}

cdnBay<-matrix(0,length(ytest),length(xtest))
for (i in 1:length(xtest)){
  cdnBay[,i]<-trap_cdf(fitBay$densp.m[i,],ytest)
}


######################################
######### Kernel smoother ############
##¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶


ptm <- proc.time()
Fhat <- npcdist(y ~ x,
                 tol = 0.1,
                 ftol = 0.1)
summary(Fhat)

ndd <- data.frame(y = ytest, x = rep(4,length(ytest)))
predDD=predict(Fhat,newdata=ndd)
proc.time() - ptm
######### Plot CDF ##################
##¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶¶




par(mfrow=c(2,3))
dd=1
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC2[,dd],lwd=1,type="l",lty=2,
     xlab="y",ylab="C.d.f.",ylim=c(0,1))
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnMC1[,dd],lwd=1,type="l",lty=2,
      col="blue")
lines(ytest,cdnMC[,dd],lwd=1,type="l",lty=2,
      col="green")
lines(ytest,cdnBay[,dd],lwd=1,type="l",lty=2,
      col="orange")
lines(ytest,predDD,lwd=1,type="l",lty=2,
      col="purple")

dd=2
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC2[,dd],lwd=1,type="l",lty=2,
     xlab="y",ylab="C.d.f.",ylim=c(0,1))
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnMC1[,dd],lwd=1,type="l",lty=2,
      col="blue")
lines(ytest,cdnMC[,dd],lwd=1,type="l",lty=2,
      col="green")
lines(ytest,cdnBay[,dd],lwd=1,type="l",lty=2,
      col="orange")
lines(ytest,predDD,lwd=1,type="l",lty=2,
      col="purple")

dd=3
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC2[,dd],lwd=1,type="l",lty=2,
     xlab="y",ylab="C.d.f.",ylim=c(0,1))
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnMC1[,dd],lwd=1,type="l",lty=2,
      col="blue")
lines(ytest,cdnMC[,dd],lwd=1,type="l",lty=2,
      col="green")
lines(ytest,cdnBay[,dd],lwd=1,type="l",lty=2,
      col="orange")
lines(ytest,predDD,lwd=1,type="l",lty=2,
      col="purple")

dd=4
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC2[,dd],lwd=1,type="l",lty=2,
     xlab="y",ylab="C.d.f.",ylim=c(0,1))
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnMC1[,dd],lwd=1,type="l",lty=2,
      col="blue")
lines(ytest,cdnMC[,dd],lwd=1,type="l",lty=2,
      col="green")
lines(ytest,cdnBay[,dd],lwd=1,type="l",lty=2,
      col="orange")
lines(ytest,predDD,lwd=1,type="l",lty=2,
      col="purple")

dd=5
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC2[,dd],lwd=1,type="l",lty=2,
     xlab="y",ylab="C.d.f.",ylim=c(0,1))
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnMC1[,dd],lwd=1,type="l",lty=2,
      col="blue")
lines(ytest,cdnMC[,dd],lwd=1,type="l",lty=2,
      col="green")
lines(ytest,cdnBay[,dd],lwd=1,type="l",lty=2,
      col="orange")
lines(ytest,predDD,lwd=1,type="l",lty=2,
      col="purple")

dd=6
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC2[,dd],lwd=1,type="l",lty=2,
     xlab="y",ylab="C.d.f.",ylim=c(0,1))
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnMC1[,dd],lwd=1,type="l",lty=2,
      col="blue")
lines(ytest,cdnMC[,dd],lwd=1,type="l",lty=2,
      col="green")
lines(ytest,cdnBay[,dd],lwd=1,type="l",lty=2,
      col="orange")
lines(ytest,predDD,lwd=1,type="l",lty=2,
      col="purple")




proc.time() - ptm
# 
#Paper

par(mfrow=c(4,2),mar=c(2,2.5,1.5,1),oma=c(2,2,0,0))
dd=1
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC1[,dd],lwd=2,type="l",lty=1,
     xlab="y",ylab="Cdf",ylim=c(0,1),main="x=0.5",cex.axis=1.5,cex=1.5,cex.lab=1.5)
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnBay[,dd],lwd=2,type="l",lty=1,
      col="orange")
lines(ytest,predDD,lwd=2,type="l",lty=1,
      col="purple")
legend(11,0.8, c("Emp","Rec","Bay","Kern"),lty=c(1,1,1,1),
       col=c("Red","Black","Orange","Purple"),bty="n")

dd=2
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC1[,dd],lwd=2,type="l",lty=1,
     xlab="y",ylab="Cdf",ylim=c(0,1),main="x=1.5",cex.axis=1.5,cex=1.5,cex.lab=1.5)
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnBay[,dd],lwd=2,type="l",lty=1,
      col="orange")
lines(ytest,predDD,lwd=2,type="l",lty=1,
      col="purple")
legend(11,0.8, c("Emp","Rec","Bay","Kern"),lty=c(1,1,1,1),
       col=c("Red","Black","Orange","Purple"),bty="n")
dd=3
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC1[,dd],lwd=2,type="l",lty=1,
     xlab="y",ylab="Cdf",ylim=c(0,1),main="x=2.5",cex.axis=1.5,cex=1.5,cex.lab=1.5)
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnBay[,dd],lwd=2,type="l",lty=1,
      col="orange")
lines(ytest,predDD,lwd=2,type="l",lty=1,
      col="purple")
legend(11,0.8, c("Emp","Rec","Bay","Kern"),lty=c(1,1,1,1),
       col=c("Red","Black","Orange","Purple"),bty="n")

dd=4
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC1[,dd],lwd=2,type="l",lty=1,
     xlab="y",ylab="Cdf",ylim=c(0,1),main="x=3.5",cex.axis=1.5,cex=1.5,cex.lab=1.5)
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnBay[,dd],lwd=2,type="l",lty=1,
      col="orange")
lines(ytest,predDD,lwd=2,type="l",lty=1,
      col="purple")
legend(11,0.8, c("Emp","Rec","Bay","Kern"),lty=c(1,1,1,1),
       col=c("Red","Black","Orange","Purple"),bty="n")

dd=5
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC1[,dd],lwd=2,type="l",lty=1,
     xlab="y",ylab="Cdf",ylim=c(0,1),main="x=4.5",cex.axis=1.5,cex=1.5,cex.lab=1.5)
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnBay[,dd],lwd=2,type="l",lty=1,
      col="orange")
lines(ytest,predDD,lwd=2,type="l",lty=1,
      col="purple")
legend(11,0.8, c("Emp","Rec","Bay","Kern"),lty=c(1,1,1,1),
       col=c("Red","Black","Orange","Purple"),bty="n")

dd=6
aa=ecdf(y[x==xtest[dd]])
ndd <- data.frame(y = ytest, x = rep(xtest[dd],length(ytest)))
predDD=predict(Fhat,newdata=ndd)
plot(ytest,cdnMC1[,dd],lwd=2,type="l",lty=2,
     xlab="y",ylab="Cdf",ylim=c(0,1),main="x=5.5",cex.axis=1.5,cex=1.5,cex.lab=1.5)
lines(ytest,aa(ytest),lwd=2,
      type="l",lty=1,col="red")
lines(ytest,cdnBay[,dd],lwd=2,type="l",lty=1,
      col="orange")
lines(ytest,predDD,lwd=2,type="l",lty=1,
      col="purple")
legend(11,0.8, c("Emp","Rec","Bay","Kern"),lty=c(1,1,1,1),
       col=c("Red","Black","Orange","Purple"),bty="n")

