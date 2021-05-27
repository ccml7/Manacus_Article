# Supplementary figure S2
data<-read.csv("~/Dropbox/Tesis Camilo/Data/Data_Manacus.csv")
dat4fig	<-	data.frame(Subspecies=data$Specie,Sex=data$Sex,Source=data$Source
						,Wing=scale(data$Wing),Tail=scale(data$Tail),Culmen=scale(data$Exposed_Culmen)
						,Tarsus=scale(data$Tarsus),Mass=scale(data$Body_mass))
dat4fig$Sex[dat4fig$Sex==""]	<- NA
dat4fig$Sex	<-	as.factor(as.character(dat4fig$Sex))

#Subspecies differences
summary(aov(Wing~Subspecies,data=dat4fig))
summary(aov(Tail~Subspecies,data=dat4fig))
summary(aov(Culmen~Subspecies,data=dat4fig))
summary(aov(Tarsus~Subspecies,data=dat4fig))
summary(aov(Mass~Subspecies,data=dat4fig))

# Sex differences
summary(aov(Wing~Sex,data=dat4fig))
summary(aov(Tail~Sex,data=dat4fig))
summary(aov(Culmen~Sex,data=dat4fig))
summary(aov(Tarsus~Sex,data=dat4fig))
summary(aov(Mass~Sex,data=dat4fig))

# Source differences
summary(aov(Wing~Source,data=dat4fig))
summary(aov(Tail~Source,data=dat4fig))
summary(aov(Culmen~Source,data=dat4fig))
summary(aov(Tarsus~Source,data=dat4fig))
summary(aov(Mass~Source,data=dat4fig))

tiff("~/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v2/FigureS2.tiff",width=7,height=7,units="in",res=300
		,compression="lzw",type="cairo",family='serif')
mat=matrix(c(1,1,2,2,0,0,0,0,0,3,3,0),ncol=4,byrow=TRUE)
layout(mat,heights=c(1,0.05,1))
par(mar=c(3,3,1,1),oma=c(1,1,1.5,1),mgp=c(1.75,0.5,0),tcl=-0.25)
# Figure A
boxplot(Wing~Subspecies,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20"
		,at=1:2,xlim=c(0,15),ylim=c(-5,5),xaxt="n",ylab="Scaled Trait Value")
boxplot(Tail~Subspecies,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=4:5,add=TRUE,xaxt="n",yaxt="n")
boxplot(Culmen~Subspecies,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=7:8,add=TRUE,xaxt="n",yaxt="n")
boxplot(Tarsus~Subspecies,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=10:11,add=TRUE,xaxt="n",yaxt="n")
boxplot(Mass~Subspecies,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=13:14,add=TRUE,xaxt="n",yaxt="n")
axis(1,at=c(1.5,4.5,7.5,10.5,13.5),labels=c("Wing","Tail","Culmen","Tarsus","Body Mass"))
mtext("*",side=1,at=c(1.5,10.5),cex=2,col="red",line=-1)
abline(v=c(3,6,9,12),col="grey",lty=2)
legend("topright",legend=c("manacus","vitellinus"),fill=c('#1f78b4','#b2df8a'))
mtext("A.",at=-1.5,cex=1.5)
# Figure B
boxplot(Wing~Sex,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20"
		,at=1:2,xlim=c(0,15),ylim=c(-5,5),xaxt="n",ylab="Scaled Trait Value")
boxplot(Tail~Sex,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=4:5,add=TRUE,xaxt="n",yaxt="n")
boxplot(Culmen~Sex,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=7:8,add=TRUE,xaxt="n",yaxt="n")
boxplot(Tarsus~Sex,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=10:11,add=TRUE,xaxt="n",yaxt="n")
boxplot(Mass~Sex,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=13:14,add=TRUE,xaxt="n",yaxt="n")
axis(1,at=c(1.5,4.5,7.5,10.5,13.5),labels=c("Wing","Tail","Culmen","Tarsus","Body Mass"))
mtext("*",side=1,at=c(1.5,4.5,10.5,13.5),cex=2,col="red",line=-1)
abline(v=c(3,6,9,12),col="grey",lty=2)
legend("topright",legend=c("Female","Male"),fill=c('#1f78b4','#b2df8a'))
mtext("B.",at=-1.5,cex=1.5)
# Figure C
boxplot(Wing~Source,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20"
		,at=1:2,xlim=c(0,15),ylim=c(-5,5),xaxt="n",ylab="Scaled Trait Value")
boxplot(Tail~Source,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=4:5,add=TRUE,xaxt="n",yaxt="n")
boxplot(Culmen~Source,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=7:8,add=TRUE,xaxt="n",yaxt="n")
boxplot(Tarsus~Source,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=10:11,add=TRUE,xaxt="n",yaxt="n")
boxplot(Mass~Source,data=dat4fig,col=c('#1f78b4','#b2df8a'),border="grey20",at=13:14,add=TRUE,xaxt="n",yaxt="n")
axis(1,at=c(1.5,4.5,7.5,10.5,13.5),labels=c("Wing","Tail","Culmen","Tarsus","Body Mass"))
mtext("*",side=1,at=c(1.5,4.5,7.5,10.5),cex=2,col="red",line=-1)
abline(v=c(3,6,9,12),col="grey",lty=2)
legend("topright",legend=c("Museum","Net"),fill=c('#1f78b4','#b2df8a'))
mtext("C.",at=-1.5,cex=1.5)
dev.off()


