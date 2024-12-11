##########################################
######## Recursive functions for R #######
##########################################

#true distribution
dtrue <- function(grid, x) {
  w1*pnorm(grid,mean=betaA*x,sd=sqrt(sigma2A))+(1-w1)*pnorm(grid,mean=betaB*x,sd=sqrt(sigma2B))
}

#recursive procedure
recursive<-function(y,x,n,cdn,xtest,copula){
  #needed input
  rate=0.8 #fix the rate of the det sequence
  df=8 # degrees of freedom Student t
  alpha=(1+seq(1,n,1))^(-rate)
  rhoMax=0.97 #fix the max correlation picked in the procedure
  dist_x=abs(matrix(rep(x,size[2]),ncol=size[2])-matrix(rep(xtest,length(x)),ncol=size[2],byrow=TRUE))
  rho=rhoMax*(1-dist_x/(tau+dist_x))
  
  #preallocate matrix
  saveP<-rep(0,n)
  cdnSAVE<-array(0, dim=c(size[1],size[2],(n+1)))
  cdnSAVE[,,1]<-cdn
  
  if (identical(copula,"StudentT")){
  for (i in 1:n){
    cdnPrev<-cdn
    if (any(is.finite(cdn))){
      idx=which(x[i]==xtest)
      idy=which(y[i]==ytest)
      u2=cdnPrev[idy,idx]  
      num=qt(cdnPrev,df)-matrix(rep(qt(u2,df)*rho[i,],size[1]), ncol=size[2],byrow=TRUE)
      den=matrix(rep(sqrt((tau+qt(u2,df)^2)/(tau+1)*(1-rho[i,]^2)),size[1]), ncol=size[2],byrow=TRUE)
      partCop=pt(num/den,df)
      cdn=(1-alpha[i])*cdnPrev+alpha[i]*partCop
      cdnSAVE[,,i+1]=cdn
    }
    else {
      saveP[i]=1  
    }
    
  }
  
  output=list(cdn=cdn,
                 saveP=saveP,
              cdnSAVE=cdnSAVE,
              rho=rho)
  return(output)
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
        cdnSAVE[,,i+1]=partCop-cdnPrev
      }
      else {
        saveP[i]=1  
      }
      
    }
    
    output=list(cdn=cdn,
                saveP=saveP,
                cdnSAVE=cdnSAVE,
                rho=rho)
    return(output)
    
  }

}

#######
#recursive POLYAK  procedure
#########
recursivePOLYAK<-function(y,x,n,cdn,xtest,copula){
  #needed input
  rate=0.55 #fix the rate of the det sequence
  df=10 # degrees of freedom Student t
  alpha=(1+seq(1,n,1))^(-rate)
  rhoMax=0.99 #fix the max correlation picked in the procedure
  dist_x=abs(matrix(rep(x,size[2]),ncol=size[2])-matrix(rep(xtest,length(x)),ncol=size[2],byrow=TRUE))
  rho=rhoMax*(1-dist_x/(tau+dist_x))
  
  #preallocate matrix
  saveP<-rep(0,n)
  cdnSAVE<-array(0, dim=c(size[1],size[2],(n+1)))
  cdnSAVE[,,1]<-cdn
  


    for (i in 1:n){
      cdnPrev<-cdn
      if (any(is.finite(cdn))){
        idx=which(x[i]==xtest)
        idy=which(y[i]==ytest)
        u2=cdnPrev[idy,idx]  
        num=qt(cdnPrev,df)-matrix(rep(qt(u2,df)*rho[i,],size[1]), ncol=size[2],byrow=TRUE)
        den=matrix(rep(sqrt((tau+qt(u2,df)^2)/(tau+1)*(1-rho[i,]^2)),size[1]), ncol=size[2],byrow=TRUE)
        partCop=pt(num/den,df)
        index=max(1,i-10) #make a rolling window
        cdnAv=apply(cdnSAVE[,,index:i],c(1,2),mean)
        cdn=(1-alpha[i])*cdnAv+alpha[i]*partCop
        cdnSAVE[,,i+1]=cdn
      }
      else {
        saveP[i]=1  
      }
      
    }
    
    output=list(cdn=cdn,
                saveP=saveP)
    return(output)

}



#recursive Capped
recursiveCAPPED<-function(y,x,n,cdn,xtest,copula){
  #needed input
  rate=0.8 #fix the rate of the det sequence
  df=1 # degrees of freedom Student t
  alpha=(1+seq(1,n,1))^(-rate)
  rhoMax=0.99 #fix the max correlation picked in the procedure
  dist_x=abs(matrix(rep(x,size[2]),ncol=size[2])-matrix(rep(xtest,size[1]),ncol=size[2],byrow=TRUE))
  rho=rhoMax*(1-dist_x/(tau+dist_x))
  
  #preallocate matrix
  saveP<-rep(0,n)
 par=0.4
  

    for (i in 1:n){
      cdnPrev<-cdn
      if (any(is.finite(cdn))){
        idx=which(x[i]==xtest)
        idy=which(y[i]==ytest)
        u2=cdnPrev[idy,idx]  
        num=qt(cdnPrev,df)-matrix(rep(qt(u2,df)*rho[i,],size[1]), ncol=size[2],byrow=TRUE)
        den=matrix(rep(sqrt((tau+qt(u2,df)^2)/(tau+1)*(1-rho[i,]^2)),size[1]), ncol=size[2],byrow=TRUE)
        partCop=pt(num/den,df)
        diffC=partCop-cdnPrev
        sig=sign(diffC)
        diffC=abs(diffC)
        diffC[diffC>par]=par
        cdn=cdnPrev+alpha[i]*(sig*diffC)
      }
      else {
        saveP[i]=1  
      }
      
    }
    
    output=list(cdn=cdn,
                saveP=saveP)
    return(output)

}
