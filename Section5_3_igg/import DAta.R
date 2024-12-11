######################################
####### data load and manipulaton ####
######################################


setwd("/Users/lorenzocappello/Google Drive/phd/Walker PhD/nonparametrci categorical regression/dataset section 5")

#schoolgirls x group
group1<-schoolgirls[schoolgirls$group==1,]
group2<-schoolgirls[schoolgirls$group==2,]
group3<-schoolgirls[schoolgirls$group==3,]

hist(group1$height,breaks=30,freq=FALSE, xlim=c(min(schoolgirls$height),max(schoolgirls$height)))
hist(group2$height,col=rgb(0,0,1,0.5),breaks=30, add=T,freq=FALSE)
hist(group3$height,col=rgb(1,0,0,0.5),breaks=30,add=T,freq=FALSE)

#schoolgirls x age
group1<-schoolgirls[schoolgirls$age==6,]
group2<-schoolgirls[schoolgirls$age==7,]
group3<-schoolgirls[schoolgirls$age==8,]
group4<-schoolgirls[schoolgirls$age==9,]
group5<-schoolgirls[schoolgirls$age==10,]

hist(group1$height,breaks=30,freq=FALSE, xlim=c(min(schoolgirls$height),max(schoolgirls$height)))
hist(group2$height,col=rgb(0,0,1,0.5),breaks=30, add=T,freq=FALSE)
hist(group3$height,col=rgb(1,0,0,0.5),breaks=30,add=T,freq=FALSE)
hist(group4$height,col=rgb(0,1,0,0.5),breaks=30,add=T,freq=FALSE)
hist(group5$height,col=rgb(1,1,1,0.5),breaks=30,add=T,freq=FALSE)



#fleabeatles
#COMMENT: 1) small sample size 2) three obvious categories, difficult to associate a
#number 3) fjft , mwafp the 2 best, maybe awfs


group1<-fleabeetles[fleabeetles$species==1,]
group2<-fleabeetles[fleabeetles$species==2,]
group3<-fleabeetles[fleabeetles$species==3,]

hist(group1$fjft,breaks=30,freq=FALSE, xlim=c(min(fleabeetles$fjft),max(fleabeetles$fjft)))
hist(group2$fjft,col=rgb(0,0,1,0.5),breaks=30, add=T,freq=FALSE)
hist(group3$fjft,col=rgb(1,0,0,0.5),breaks=30,add=T,freq=FALSE)





#igg data
#COMMENT: GOOD BUT NEEDS TO BE CATEGORIZED
group1<-igg$igg[igg$age<=1]
group2<-igg$igg[igg$age>1 & igg$age <=3]
group3<-igg$igg[igg$age>3 & igg$age <=5]
group4<-igg$igg[igg$age>5 ]

hist(group1,breaks=30,freq=FALSE, xlim=c(min(igg$igg),max(igg$igg)))
hist(group2,col=rgb(0,0,1,0.5),breaks=30, add=T,freq=FALSE)
hist(group3,col=rgb(1,0,0,0.5),breaks=30,add=T,freq=FALSE)
hist(group4,col=rgb(0,1,0,0.5),breaks=30,add=T,freq=FALSE)


#sports data
length(unique(sports$sport))
 group1<-sports$rcc[sports$sport=="Swim"]
 group2<-sports$rcc[sports$sport=="W_Polo"]
 group3<-sports$rcc[sports$sport=="B_Ball"]
 group4<-sports$rcc[sports$sport=="Row"]
 
 hist(group3,breaks=30,xlim=c(min(sports$rcc),max(sports$rcc)))
 hist(group1,col=rgb(0,0,1,0.5),breaks=30, add=T)
hist(group2,col=rgb(1,0,0,0.5),breaks=30,add=T)
hist(group4,col=rgb(0,1,0,0.5),breaks=30,add=T)

hist(sports$bmi,breaks=30)


#Yacht data
# mydata = read.table("yacht_hydrodynamics.txt")
# data1=data.frame(mydata[3],mydata[7])
# unique(data1[1])
# y=log(data1[2])
# x=data1[1]
# 
# 
# 
# group1<-c(y[x==4.78],y[x==4.77],y[x==4.76])
# group2<-c(y[x==4.36],y[x==4.34])
# group3<-c(y[x==5.11],y[x==5.14],y[x==5.10])
# 
# hist(group3,breaks=30,freq=FALSE)
# hist(group1,col=rgb(0,0,1,0.5),breaks=30,freq=FALSE, add=T)
# hist(group2,col=rgb(1,0,0,0.5),breaks=30,freq=FALSE, add=T)
