################################################
### Illustration Bayesian Mixture on R ###########
################################################


rm(list = setdiff(ls(), lsf.str())) #this command removes evertything except functions

library(DPpackage)
ptm <- proc.time()
#set problem parameter
tau=1
n=300 #CHECK: SOMETIMES I SAMPLE MORE TO HAVE A FINER GRID
betaA=1
betaB=2
sigma2A=1.2
sigma2B=1
y=rep(0,n)
pop=c(0.8,1,1.3,1.5)
x=sample(pop,n,replace=TRUE)
#w=1-abs(x-mean(pop))/(mean(pop))
#w=1-(1/3+x-min(x))/(2*(max(x)-min(x)))
w=1-(1/3+x-min(x))/(1+x-min(x))


u=runif(n)

#generate y
for (i in 1:n){
  if (u[i]<=w[i]){
    y[i]=rnorm(1,betaA*x[i],sigma2A)  }
  else{
    y[i]=rnorm(1,betaB*x[i],sigma2B) 
  }
}



#Create the grid of point in which to evaluate the distribution (not needed I guess)

ytest=y[order(y)]
xtest=pop[order(pop)]
#wtest=1-abs(xtest-mean(xtest))/(mean(xtest))
wtest=1-(1/3+xtest-min(xtest))/(1+xtest-min(xtest))


#reduce the sample
n=200
y=y[1:n]
x=x[1:n]

#true distribution
ptrue <- function(grid, x,wtest) {
  wtest*pnorm(grid,mean=betaA*x,sd=sqrt(sigma2A))+(1-wtest)*pnorm(grid,mean=betaB*x,sd=sqrt(sigma2B))
}


#BAYESIAN PROCEDURE 



#prior
W <- x
S0 <- 0.5
m0 <- 0

prior <- list(a0 = 5, b0 = 5, m0 = m0,
              S0 = S0, tau1 = 6.01,
              taus1 = 6.01, taus2 = 6.01,
              nu = 6, psiinv = S0)

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
fit1 <- LDDPdensity(formula=y~0+W,zpred=Wpred, prior=prior,mcmc=mcmc,
                    compute.band = FALSE,state=NULL,status=TRUE,
                    grid=ytest)


#MCMC diagnostics
#plot(fit1, hpd = FALSE)




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
cdnBay[,i]<-trap_cdf(fit1$densp.m[i,],ytest)
}

#plot CDF
dd=4
plot(ytest,fit$cdn[,2],lwd=1,type="l",lty=2,
     xlab="values",ylab="density",ylim=c(0,1))
lines(ytest,fit$cdn[,4],lwd=2,
      type="l",lty=1,col="red")

