##########################################
######## Recursive functions for R #######
##########################################


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
    q_inv<-  wtest[i]*qnorm(q,mean=betaA*xtest[i],sd=sqrt(sigma2A))+(1-wtest[i])*qnorm(q,mean=betaB*xtest[i],sd=sqrt(sigma2B))
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
  wtest*pnorm(grid,mean=betaA*x,sd=sqrt(sigma2A))+(1-wtest)*pnorm(grid,mean=betaB*x,sd=sqrt(sigma2B))
}

#density
dtrue <- function(grid, x,wtest) {
  wtest*dnorm(grid,mean=betaA*x,sd=sqrt(sigma2A))+(1-wtest)*dnorm(grid,mean=betaB*x,sd=sqrt(sigma2B))
}

#recursive procedure
recursive<-function(y,x,n,cdn,xtest,copula,df){
  
  #needed input
  rate=0.8 #fix the rate of the det sequence
  alpha=(1+seq(1,n,1))^(-rate)
  rhoMax=0.99 #fix the max correlation picked in the procedure
  dist_x=abs(matrix(rep(x,size[2]),ncol=size[2])-matrix(rep(xtest,length(x)),ncol=size[2],byrow=TRUE))
  rho=rhoMax*(1-dist_x/(tau+dist_x))
  
  #preallocate matrix
  saveP<-rep(0,n)
  
  if (identical(copula,"StudentT")){
  for (i in 1:n){
    cdnPrev<-cdn
    if (any(is.finite(cdn))){
      idx=which(x[i]==xtest)
      idy=min(which(y[i]==ytest))
      u2=cdnPrev[idy,idx]  
      num=qt(cdnPrev,df)-matrix(rep(qt(u2,df)*rho[i,],size[1]), ncol=size[2],byrow=TRUE)
      #den=matrix(rep(sqrt((tau+qt(u2,df)^2)/(tau+1)*(1-rho[i,]^2)),size[1]), ncol=size[2],byrow=TRUE)
      den=matrix(rep((1-rho[i,]^2),size[1]), ncol=size[2],byrow=TRUE)
      partCop=pt(num/den,df)
      cdn=(1-alpha[i])*cdnPrev+alpha[i]*partCop
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
        idy=min(which(y[i]==ytest))
        u2=cdnPrev[idy,idx]  
        num=qnorm(cdnPrev)-matrix(rep(qnorm(u2)*rho[i,],size[1]), ncol=size[2],byrow=TRUE)
        den=matrix(rep(sqrt(1-rho[i,]^2),size[1]), ncol=size[2],byrow=TRUE)
        partCop=pnorm(num/den)
        cdn=(1-alpha[i])*cdnPrev+alpha[i]*partCop
      }
      else {
        saveP[i]=1  
      }
      
    }
  
    
  }

  output=list(cdn=cdn,
              saveP=saveP,
              rho=rho)
  return(output)
}