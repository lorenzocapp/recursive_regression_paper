##########################################
######## Function needed for R #######
##########################################
rm(list=ls())

##### PARAMETERS ########
r=50
tau=0.75
n=250 #CHECK: SOMETIMES I SAMPLE MORE TO HAVE A FINER GRID
###CHECK:ACTUAL SAMPLE IS PICKED BELOW HERE
dfA=3
betaA=-0.75
betaB=1
sigma2A=0.5
sigma2B=1
y=rep(0,n)
pop=c(0.8,1,1.3,1.5)
ngrid=300

q=c(0.001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999)
#grid and weight for x
xtest=pop[order(pop)]
wtest=1-(1/3+xtest-min(xtest))/(1+xtest-min(xtest))
nMCint=50000
uMCint=runif(nMCint)


########


library(DPpackage)
#compute the checkloss function 
checkloss<- function(q,xtest,ytest,CDN,samplemonte){
  checkloss=matrix(0,length(q),length(xtest))
  for (i in 1:length(pop)){
    q_inv<-approx(CDN[,i], ytest, q, method = "linear", rule = 2)
    for (b in 1:length(q)){
      bb=as.numeric(I(samplemonte[,i]>q_inv$y[b]))
      checkloss[b,i]=mean(bb*samplemonte[,i])*q_inv$y[b]
    }}
  
  return(checkloss)
}

checklossTRUE<- function(q,xtest,ytest,samplemonte){
  checkloss=matrix(0,length(q),length(xtest))
  for (i in 1:length(pop)){
    q_inv<-  wtest[i]*(betaA*xtest[i]+sqrt(sigma2A)*qt(q,dfA))+(1-wtest[i])*qnorm(q,mean=betaB*xtest[i],sd=sqrt(sigma2B))
    for (b in 1:length(q)){
      bb=as.numeric(I(samplemonte[,i]>q_inv[b]))
      checkloss[b,i]=mean(bb*samplemonte[,i])*q_inv[b]
    }}
  
  return(checkloss)
}





#compute the Cramer-vonMises test statistics
CvM<-function(diff,pdf,ytest){
  f=diff*pdf
  CvM<-trap(f,ytest)
}


#compute the trapezoidal rule. 
trap<- function(f,ytest){
  diffy<-matrix(rep(ytest[2:length(ytest)]-ytest[1:length(ytest)-1],dim(f)[2]),ncol=dim(f)[2])
  sumbase<-f[1:dim(f)[1]-1,]+f[2:dim(f)[1],]
  prod<-diffy*sumbase/2
  trap<-apply(prod,2,sum)
}

#true distribution
ptrue <- function(grid, x,wtest) {
  wtest*pt((grid-betaA*x)/sqrt(sigma2A),dfA)+(1-wtest)*pnorm(grid,mean=betaB*x,sd=sqrt(sigma2B))
}

#density
dtrue <- function(grid, x,wtest) {
  wtest*dt((grid-betaA*x)/sqrt(sigma2A),dfA)+(1-wtest)*dnorm(grid,mean=betaB*x,sd=sqrt(sigma2B))
}

#recursive procedure
recursive<-function(y,x,n,cdn,xtest,copula){
  size=dim(cdn)
  #needed input
  rate=0.75 #fix the rate of the det sequence
  df=8 # degrees of freedom Student t
  alpha=(1+seq(1,n,1))^(-rate)
  rhoMax=0.98 #fix the max correlation picked in the procedure
  dist_x=abs(matrix(rep(x,size[2]),ncol=size[2])-matrix(rep(xtest,length(x)),ncol=size[2],byrow=TRUE))
  rho=rhoMax*(1-dist_x/(tau+dist_x))
  
  #preallocate matrix
  saveP<-rep(0,n)
  cdnSAVE<-array(0,c(size[1],size[2],n))
  
  if (identical(copula,"StudentT")){
    for (i in 1:n){
      cdnPrev<-cdn
      if (any(is.finite(cdn))){
        idx=which(x[i]==xtest)
        idy=which(y[i]==ytest)
        u2=cdnPrev[idy,idx]  
        num=qt(cdnPrev,df)-matrix(rep(qt(u2,df)*rho[i,],size[1]), ncol=size[2],byrow=TRUE)
        den=matrix(rep(sqrt((df+qt(u2,df)^2)/(df+1)*(1-rho[i,]^2)),size[1]), ncol=size[2],byrow=TRUE)
        partCop=pt(num/den,df+1)
        cdn=(1-alpha[i])*cdnPrev+alpha[i]*partCop
        cdnSAVE[,,i]=cdn
      }
      else {
        saveP[i]=1  
      }
      
    }
    
    
  } else {
    for (i in 1:n){
      cdnPrev<-cdn
      if (any(is.finite(cdn))){
        idx=which(x[i]==xtest)
        idy=which(y[i]==ytest)
        u2=cdnPrev[idy,idx]  
        num=qnorm(cdnPrev)-matrix(rep(qnorm(u2)*rho[i,],size[1]), ncol=size[2],byrow=TRUE)
        den=matrix(rep(sqrt(1-rho[i,]^2),size[1]), ncol=size[2],byrow=TRUE)
        partCop=pnorm(num/den)
        cdn=(1-alpha[i])*cdnPrev+alpha[i]*partCop
        cdnSAVE[,,i]=cdn
      }
      else {
        saveP[i]=1  
      }
      
    }
    
    
  }
  
  output=list(cdn=cdn,
              saveP=saveP,
              rho=rho,
              cdnSAVE=cdnSAVE)
  return(output)
}

##############################################################################################################################

##########################################
######## Start of the procedure #######
##########################################

#prepare sample for MC integration for check loss
samplemonte=matrix(0,nMCint,length(pop))
for (u in 1:length(pop)){
  for (l in 1:nMCint){
    if (uMCint[l]<=wtest[u]){
      samplemonte[l,u]=betaA*xtest[u]+sqrt(sigma2A)*rt(1,dfA)  }
    else{
      samplemonte[l,u]=rnorm(1,betaB*xtest[u],sigma2B) 
    }}}
samplemonte=as.matrix(samplemonte)

#preallocate matrices
dist_cdnTsave=matrix(0,nrow=ngrid,ncol=length(pop)*r)
dist_cdnMixsave=matrix(0,nrow=ngrid,ncol=length(pop)*r)
dist_cdnMCsave=matrix(0,nrow=ngrid,ncol=length(pop)*r)
dist_cdnORsave=matrix(0,nrow=ngrid,ncol=length(pop)*r)
dist_cdnBAYsave=matrix(0,nrow=ngrid,ncol=length(pop)*r)
dist_cdn2save=matrix(0,nrow=ngrid,ncol=length(pop)*r)
K_T=matrix(0,nrow=r,ncol=length(pop))
K_Mix=matrix(0,nrow=r,ncol=length(pop))
K_MC=matrix(0,nrow=r,ncol=length(pop))
K_OR=matrix(0,nrow=r,ncol=length(pop))
K_BAY=matrix(0,nrow=r,ncol=length(pop))
K_2=matrix(0,nrow=r,ncol=length(pop))
Cram_T<-matrix(0,nrow=r,ncol=length(pop))
Cram_Mix<-matrix(0,nrow=r,ncol=length(pop))
Cram_MC<-matrix(0,nrow=r,ncol=length(pop))
Cram_OR<-matrix(0,nrow=r,ncol=length(pop))
Cram_BAY<-matrix(0,nrow=r,ncol=length(pop))
Cram_2<-matrix(0,nrow=r,ncol=length(pop))
CheckT=matrix(0,nrow=length(q),ncol=length(pop)*r)
CheckMix=matrix(0,nrow=length(q),ncol=length(pop)*r)
CheckMC=matrix(0,nrow=length(q),ncol=length(pop)*r)
CheckOR=matrix(0,nrow=length(q),ncol=length(pop)*r)
CheckBAY=matrix(0,nrow=length(q),ncol=length(pop)*r)




######### START!!! #########
for (k in 1:r){
  
  print(k)
  #generate samples 
  x=sample(pop,n,replace=TRUE)
  w=1-(1/3+x-min(x))/(1+x-min(x))
  
  u=runif(n)
  #generate y
  for (i in 1:n){
    if (u[i]<=w[i]){
      y[i]=betaA*x[i]+sqrt(sigma2A)*rt(1,dfA)  }
    else{
      y[i]=rnorm(1,betaB*x[i],sigma2B) 
    }
  }
  
  
  
  #Create the grid of point in which to evaluate the distribution (not needed I guess)
  if (n<ngrid){
    extra=ngrid-n
    miny=min(y)*1.1
    maxy=max(y)*1.1
    dist=(maxy-miny)/(extra-1)
    extragr=seq(miny,maxy,dist)
  } else {
    extragr=numeric(0)
  }
  
  ytest=sort(c(extragr,y))
  xtest=pop[order(pop)]
  wtest=1-(1/3+xtest-min(xtest))/(1+xtest-min(xtest))
  
  ############ ############  ############
  ###### Recursive procedure ######
  ############ ############  ############
  
  
  #prior
  mup<-0
  sigmaP<-3
  cdn<-pnorm(ytest,mean=mup,sd=sqrt(sigmaP))
  cdn<-matrix(rep(cdn,length(xtest)),ncol=length(xtest))
  cdn0<- cdn #save the prior in order to use it for plot 
  
  size<-dim(cdn)
  
  
  #simple version 
  
  fit<-recursive(y,x,n,cdn,xtest,"StudentT")
  cdnTave=apply(fit$cdnSAVE,c(1,2),mean)
  
  
  ptm <- proc.time()
  
  #monte carlo
  mc=30
  cdnMC<-matrix(0,size[1],size[2])
  cdn<-cdn0
  for (t in 1:mc){
    perm=sample(seq(1,n,1),replace=FALSE)
    y=y[perm]
    x=x[perm]
    fitMC<-recursive(y,x,n,cdn,xtest,"StudentT")
    cdnMC<-cdnMC+fitMC$cdn
  }
  cdnMC<-cdnMC/mc
  
  proc.time() - ptm
  
  
  
  #order observation: due to the small sample
  
  cdnOR<-matrix(0,size[1],size[2])
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
    mc=30
    cdnMCOR<-matrix(0,size[1],size[2])
    for (t in 1:mc){
      perm=sample(seq(1,length(y2),1),replace=FALSE)
      y2=y2[perm]
      x2=x2[perm]
      fitORMC<-recursive(y2,x2,length(y2),fitORinter$cdn,xtest,"Gaussian")
      cdnMCOR<-cdnMCOR+fitORMC$cdn
    }
    cdnOR[,i]<-cdnMCOR[,i]/mc
  }
  
  
  
  cdnMix=(cdnOR+cdnMC)/2
  #################################
  ###### Bayesian procedure ######
  #################################
  ptm <- proc.time()
  
  #prior
  W <- x
  S0 <- 4
  m0 <- 0
  
  prior <- list(a0 = 8, b0 = 1, m0 = m0,
                S0 = S0, tau1 = 6.01,
                taus1 = 6.01, taus2 = 2.01,
                nu = 1, psiinv = 0.25)
  
  #grid in which to predict
  Wpred<- matrix(rep(c(1,xtest),n),nrow=n,ncol=length(c(1,xtest)),byrow=TRUE)
  Wpred<-t(xtest)
  Wpred<-t(Wpred)
  
  
  
  # Initial state
  state <- NULL
  
  
  # MCMC parameters
  nburn <- 500
  nsave <- 500
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
  proc.time() - ptm
  #################################
  ###### Compute the Metrics ######
  #################################
  
  
  cdnTrue=matrix(0,size[1],size[2])
  pdfTrue=matrix(0,size[1],size[2])
  for (i in 1:length(xtest)){
    cdnTrue[,i]=ptrue(ytest,xtest[i],wtest[i])
    pdfTrue[,i]=dtrue(ytest,xtest[i],wtest[i])
  }
  
  
  dist_cdnTsave[,(4*(k-1)+1):(4*(k-1)+4)]=abs(cdnTave-cdnTrue)
  dist_cdnMixsave[,(4*(k-1)+1):(4*(k-1)+4)]=abs(cdnMix-cdnTrue)
  dist_cdnMCsave[,(4*(k-1)+1):(4*(k-1)+4)]=abs(cdnMC-cdnTrue)
  dist_cdnORsave[,(4*(k-1)+1):(4*(k-1)+4)]=abs(cdnOR-cdnTrue)
  dist_cdnBAYsave[,(4*(k-1)+1):(4*(k-1)+4)]=abs(cdnBay-cdnTrue)
  dist_cdn2save[,(4*(k-1)+1):(4*(k-1)+4)]=abs(cdnBay-cdnMC)
  
  K_T[k,]=apply(abs(cdnTave-cdnTrue),2,max)
  K_Mix[k,]=apply(abs(cdnMix-cdnTrue),2,max)
  K_MC[k,]=apply(abs(cdnMC-cdnTrue),2,max)
  K_OR[k,]=apply(abs(cdnOR-cdnTrue),2,max)
  K_BAY[k,]=apply(abs(cdnBay-cdnTrue),2,max)
  K_2[k,]=apply(abs(cdnBay-cdnMC),2,max)
  
  Cram_T[k,]<-CvM(abs(cdnTave-cdnTrue),pdfTrue,ytest)
  Cram_Mix[k,]<-CvM(abs(cdnMix-cdnTrue),pdfTrue,ytest)
  Cram_MC[k,]<-CvM(abs(cdnMC-cdnTrue),pdfTrue,ytest)
  Cram_OR[k,]<-CvM(abs(cdnOR-cdnTrue),pdfTrue,ytest)
  Cram_BAY[k,]<-CvM(abs(cdnBay-cdnTrue),pdfTrue,ytest)
  Cram_2[k,]<-CvM(abs(cdnBay-cdnMC),pdfTrue,ytest)
  
  
  
  CheckT[,(4*(k-1)+1):(4*(k-1)+4)]<-checkloss(q,xtest,ytest,cdnTave,samplemonte)
  CheckMix[,(4*(k-1)+1):(4*(k-1)+4)]<-checkloss(q,xtest,ytest,cdnMix,samplemonte)
  CheckMC[,(4*(k-1)+1):(4*(k-1)+4)]<-checkloss(q,xtest,ytest,cdnMC,samplemonte)
  CheckOR[,(4*(k-1)+1):(4*(k-1)+4)]<-checkloss(q,xtest,ytest,cdnOR,samplemonte)
  CheckBAY[,(4*(k-1)+1):(4*(k-1)+4)]<-checkloss(q,xtest,ytest,cdnBay,samplemonte)
  
  
  
  
  
} #this ends the main for loop


savemat_bis<-list(distT=dist_cdnTsave,
              distMix=dist_cdnMixsave,
              distMC=dist_cdnMCsave,
              distOR=dist_cdnORsave,
              distBAY=dist_cdnBAYsave,
              dist2=dist_cdn2save,
              K_T=K_T,
              K_Mix=K_Mix,
              K_MC=K_MC,
              K_OR=K_OR,
              K_BAY=K_BAY,
              K_2=K_2,
              Cram_T=Cram_T,
              Cram_Mix=Cram_Mix,
              Cram_MC=Cram_MC,
              Cram_OR=Cram_OR,
              Cram_BAY=Cram_BAY,
              Cram_2=Cram_2,
              CheckT=CheckT,
              CheckMix=CheckMix,
              CheckMC=CheckMC,
              CheckOR=CheckOR,
              CheckBAY=CheckBAY,
              n=n)
setwd("/Users/lorenzocappello/Google Drive/phd/Walker PhD/nonparametrci categorical regression/illustration section 5 T- Gauss copia")
#save(savemat_bis,file="illv3_n=250.RData")