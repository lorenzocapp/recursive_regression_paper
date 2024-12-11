##Plot results



# par(mfrow=c(2,2))
# for (i in 1:4){
# toplot<-cbind(rr$sK_MC7[,i],rr$sK_MC8[,i],rr$sK_MC9[,i],rr$sK_MC10[,i],rr$sK_Ker[,i],rr$sK_BAY[,i])
# boxplot(toplot)
# }


# par(mfrow=c(2,4))
# for (i in 1:4){
#   toplot<-c(rr$sCram_MC7[,i],rr$sCram_MC8[,i],rr$sCram_MC9[,i],rr$sCram_MC10[,i],rr$sCram_Ker[,i],rr$sCram_Bay[,i])
#   boxplot(toplot)
# }
# 
# library(ggplot2)
# ggplot( data.frame(toplot))+geom_boxplot()

library(ggplot2)
library(ggpubr)

i=1
KolSmir<-c(rr$sK_MC7[,i],rr$sK_MC8[,i],rr$sK_MC9[,i],rr$sK_MC10[,i],rr$sK_Ker[,i],rr$sK_BAY[,i])
model<-c(rep("RR_0.7",50),rep("RR_0.8",50),rep("RR_0.9",50),rep("RR_1.0",50),rep("Ker",50),
         rep("Bay",50))
data<-data.frame(model,KolSmir)
p1<-qplot( x=model, y=KolSmir , data=data, geom=c("boxplot","jitter"),show.legend=FALSE)+theme_bw()+ theme(axis.title = element_blank(),
                                                                                                           axis.text  = element_text(size=12,angle=90))


i=2
KolSmir<-c(rr$sK_MC7[,i],rr$sK_MC8[,i],rr$sK_MC9[,i],rr$sK_MC10[,i],rr$sK_Ker[,i],rr$sK_BAY[,i])
model<-c(rep("RR_0.7",50),rep("RR_0.8",50),rep("RR_0.9",50),rep("RR_1.0",50),rep("Ker",50),
         rep("Bay",50))
data<-data.frame(model,KolSmir)
p2<-qplot( x=model, y=KolSmir , data=data, geom=c("boxplot","jitter"),show.legend=FALSE)+theme_bw()+ theme(axis.title = element_blank(),
                                                                                                           axis.text  = element_text(size=12,angle=90))

i=3
KolSmir<-c(rr$sK_MC7[,i],rr$sK_MC8[,i],rr$sK_MC9[,i],rr$sK_MC10[,i],rr$sK_Ker[,i],rr$sK_BAY[,i])
model<-c(rep("RR_0.7",50),rep("RR_0.8",50),rep("RR_0.9",50),rep("RR_1.0",50),rep("Ker",50),
         rep("Bay",50))
data<-data.frame(model,KolSmir)
p3<-qplot( x=model, y=KolSmir , data=data, geom=c("boxplot","jitter"),show.legend=FALSE)+theme_bw()+ theme(axis.title = element_blank(),
                                                                                                           axis.text  = element_text(size=12,angle=90))

i=4
KolSmir<-c(rr$sK_MC7[,i],rr$sK_MC8[,i],rr$sK_MC9[,i],rr$sK_MC10[,i],rr$sK_Ker[,i],rr$sK_BAY[,i])
model<-c(rep("RR_0.7",50),rep("RR_0.8",50),rep("RR_0.9",50),rep("RR_1.0",50),rep("Ker",50),
         rep("Bay",50))
data<-data.frame(model,KolSmir)
p4<-qplot( x=model, y=KolSmir , data=data, geom=c("boxplot","jitter"),show.legend=FALSE)+theme_bw()+ theme(axis.title = element_blank(),
                                                                                                           axis.text  = element_text(size=12,angle=90))




i=1
KolSmir<-c(rr$sCram_MC7[,i],rr$sCram_MC8[,i],rr$sCram_MC9[,i],rr$sCram_MC10[,i],rr$sCram_Ker[,i],rr$sCram_Bay[,i])
model<-c(rep("RR_0.7",50),rep("RR_0.8",50),rep("RR_0.9",50),rep("RR_1.0",50),rep("Ker",50),
                                                                rep("Bay",50))
data<-data.frame(model,KolSmir)
p5<-qplot( x=model, y=KolSmir , data=data, geom=c("boxplot","jitter"),show.legend=FALSE)+theme_bw()+ theme(axis.title = element_blank(),
                                                                                                      axis.text  = element_text(size=12,angle=90))


i=2
KolSmir<-c(rr$sCram_MC7[,i],rr$sCram_MC8[,i],rr$sCram_MC9[,i],rr$sCram_MC10[,i],rr$sCram_Ker[,i],rr$sCram_Bay[,i])
model<-c(rep("RR_0.7",50),rep("RR_0.8",50),rep("RR_0.9",50),rep("RR_1.0",50),rep("Ker",50),
         rep("Bay",50))
data<-data.frame(model,KolSmir)
p6<-qplot( x=model, y=KolSmir , data=data, geom=c("boxplot","jitter"),show.legend=FALSE)+theme_bw()+ theme(axis.title = element_blank(),
                                                                                                           axis.text  = element_text(size=12,angle=90))
i=3
KolSmir<-c(rr$sCram_MC7[,i],rr$sCram_MC8[,i],rr$sCram_MC9[,i],rr$sCram_MC10[,i],rr$sCram_Ker[,i],rr$sCram_Bay[,i])
model<-c(rep("RR_0.7",50),rep("RR_0.8",50),rep("RR_0.9",50),rep("RR_1.0",50),rep("Ker",50),
         rep("Bay",50))
data<-data.frame(model,KolSmir)
p7<-qplot( x=model, y=KolSmir , data=data, geom=c("boxplot","jitter"),show.legend=FALSE)+theme_bw()+ theme(axis.title = element_blank(),
                                                                                                           axis.text  = element_text(size=12,angle=90))
i=4
KolSmir<-c(rr$sCram_MC7[,i],rr$sCram_MC8[,i],rr$sCram_MC9[,i],rr$sCram_MC10[,i],rr$sCram_Ker[,i],rr$sCram_Bay[,i])
model<-c(rep("RR_0.7",50),rep("RR_0.8",50),rep("RR_0.9",50),rep("RR_1.0",50),rep("Ker",50),
         rep("Bay",50))
data<-data.frame(model,KolSmir)
p8<-qplot( x=model, y=KolSmir , data=data, geom=c("boxplot","jitter"),show.legend=FALSE)+theme_bw()+ theme(axis.title = element_blank(),
                                                                                                           axis.text  = element_text(size=12,angle=90))

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=3,ncol=4)

                                                                                                       axis.text  = element_text(size=12))

