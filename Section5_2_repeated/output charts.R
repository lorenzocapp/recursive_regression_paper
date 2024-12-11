
#various plot
hist(K_T[,1],breaks=20)
hist(K_T[,2],add=T,col="red",breaks=20)
hist(K_T[,3],add=T,col="blue",breaks=20)
hist(K_T[,4],add=T,col="green",breaks=20)


hist(K_T[,1],breaks=20)
hist(K_MC[,1],add=T,col="red",breaks=20)
hist(K_G[,1],add=T,col="blue",breaks=20)
hist(K_BAY[,1],add=T,col="green",breaks=20)
 
idx=3
boxplot(savemat_bis$Cram_MC10[,idx], savemat_bis$Cram_Ker[,idx], savemat_bis$Cram_MC30[,idx],
        savemat_bis$Cram_MC50[,idx], savemat_bis$Cram_BAY[,idx], savemat_bis$Cram_G30[,idx],
        names=c("MC10","Ker","MC30","MC50","Bayes","G30"),main="Cramer von Mises x4")
idx=4
boxplot(savemat_bis$K_MC10[,idx], savemat_bis$K_Ker[,idx], savemat_bis$K_MC30[,idx],
        savemat_bis$K_MC50[,idx], savemat_bis$K_BAY[,idx], savemat_bis$K_G30[,idx],
        names=c("MC10","Ker","MC30","MC50","Bayes","G30"),main="Kolm Smirn x1")

idx=4
boxplot(savemat_bis$Cram_T[,idx], savemat_bis$Cram_Ker[,idx], savemat_bis$Cram_MC[,idx],
        savemat_bis$Cram_OR[,idx], savemat_bis$Cram_BAY[,idx], savemat_bis$Cram_2[,idx],
        names=c("StudT","Gauss","MCT","MCG","Bayes","2"),main="Cramer von Mises x4")
idx=3
boxplot(Cram_T[,idx], Cram_Ker[,idx], Cram_MC[,idx],
        Cram_OR[,idx], Cram_BAY[,idx], Cram_2[,idx],
        names=c("StudT","Gauss","MCT","MCG","Bayes","2"),main="Cramer von Mises x4")



idx=4
boxplot(savemat_bis$K_T[,idx], savemat_bis$K_Ker[,idx], savemat_bis$K_MC[,idx], 
        savemat_bis$K_OR[,idx], savemat_bis$K_BAY[,idx], savemat_bis$K_2[,idx],
        names=c("StudT","Kern","MCT","OR","Bayes","2"),main="Kolm Smirnov x4")

par(mfrow=c(4,2))


plot(ytest,cdnTrue[,1],lwd=1,type="l",lty=2,
     xlab="values",ylab="density",ylim=c(0,1))
lines(ytest,cdnTrue[,2],lwd=1,
      type="l",lty=2,col="red")
lines(ytest,cdnTrue[,3],lwd=1,type="l",lty=2,
      col="blue")
lines(ytest,cdnTrue[,4],lwd=1,type="l",lty=2,
      col="green")


####manipulation of the L1 distance-> obtain the average at each grid point
## and plot it 
idx=4
meanL1_Bay<-apply(savemat_bis$distBAY[,seq(idx,(r*4),4)],1,mean)
meanL1_MC<-apply(savemat$distMC[,seq(idx,(r*4),4)],1,mean)
meanL1_MCG<-apply(savemat$distMCG[,seq(idx,(r*4),4)],1,mean)
meanL1_T<-apply(savemat$distT[,seq(idx,(r*4),4)],1,mean)
meanL1_G<-apply(savemat$distG[,seq(idx,(r*4),4)],1,mean)
meanL1_2<-apply(savemat$dist2[,seq(idx,(r*4),4)],1,mean)

plot(ytest,meanL1_T,lwd=1,type="l",lty=2,col="black",
     xlab="values",ylab="L1dist",main="L1 dist x4")
lines(ytest,meanL1_G,lwd=1,type="l",lty=2,col="red")
lines(ytest,meanL1_MC,lwd=1,type="l",lty=2,col="blue")
lines(ytest,meanL1_MCG,lwd=1,type="l",lty=2,col="orange")
lines(ytest,meanL1_Bay,lwd=1,type="l",lty=2,col="green")
legend(3.5,0.075, c("StudT","Gauss","MC_T","MCG","Bayes"),lty=c(1,1,1,1),
        col=c("black","blue","red","orange","green"))



####check loss function analysis
CheckTrue=checklossTRUE(q,xtest,ytest,samplemonte)
CheckTrueall=matrix(rep(CheckTrue,r),nrow=9,ncol=(4*r))
ratioCheck=(savemat_bis$CheckMC10-savemat_bis$CheckBAY)/CheckTrueall

idx=4
boxplot(ratioCheck[,seq(idx,(r*4),4)],use.cols=FALSE,
        main="CheckLoss x4",
        names=c("0.001","0.01","0.1","0.25","0.5","0.75","0.9","0.99","0.999"))
abline(h=0)
 

  
idx=4
boxplot(ratioCheck[-2,seq(idx,(r*4),4)],use.cols=FALSE,
        main="CheckLoss x4 -0.01",
        names=c("0.001","0.1","0.25","0.5","0.75","0.9","0.99","0.999"))
abline(h=0)
