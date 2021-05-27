library(lme4)
library(ape)
library(geosphere)
library(nlme)
library(mgcv)

dat	<-	read.csv("~/Dropbox/Tesis Camilo/Data/Data_Manacus.txt")
geo.dist	<-	distm(dat[,c("Long","Lat")])/1000
geo.pcoa	<-	cmdscale(as.dist(geo.dist),k=1)
sc.pcoa		<-	scale(geo.pcoa)

new.dat		<-	data.frame(Wing=dat$Wing,Tail=dat$Tail,Culmen=dat$Exposed_Culmen
							,Tarsus=dat$Tarsus,Mass=dat$Body_mass
							,bio1=scale(dat$bio1),bio3=scale(dat$bio3)
							,bio7=scale(dat$bio7),bio12=scale(dat$bio12)
							,bio15=scale(dat$bio15),bio18=scale(dat$bio18)
							,Altitude=scale(dat$Altitude),geog=sc.pcoa
							,Sex=factor(dat$Sex,levels=c("Female","Male")),Source=dat$Source)
model1<-lmer(Mass~bio1+bio7+bio12+bio15+bio18+sc.pcoa+(1|Sex)+(1|Source),data=new.dat,REML=T)
model2<-lmer(Mass~bio1+bio7+bio12+bio15+bio18+sc.pcoa+(1|Source),data=new.dat,REML=T)	
model3<-lmer(Mass~bio1+bio7+bio12+bio15+bio18+sc.pcoa+(1|Sex),data=new.dat,REML=T)	
model4<-lmer(Mass~bio1+bio7+bio12+bio15+bio18+(1|Sex),data=new.dat,REML=T)

AIC(model2)
anova(model1,model3,model4)

tmp<-lme(Mass~bio1+bio7+bio12+bio15+bio18+sc.pcoa,random=~1|Sex,data=new.dat,na.action=na.omit)	

tmp1<-lme(Mass~bio1+bio7+bio12+bio18,random=~1|Sex,data=new.dat,na.action=na.omit)	

pred.dat	<-	data.frame(bio1=rep(seq(-2.2,1.5,length=100),2),bio7=rep(0,200)
							,bio12=rep(0,200),bio15=rep(0,200),bio18=rep(0,200)
							,Sex=c(rep("Male",100),rep("Female",100)))

fitted.vals	<-	predict(tmp1,newdata=pred.dat)

plot(pred.dat$bio1[1:100],fitted.vals[1:100],type="l",ylim=range(new.dat$Mass,na.rm=TRUE))
points(pred.dat$bio1[101:200],fitted.vals[101:200],type="l",col="red")
points(new.dat$bio1[new.dat$Sex=="Male"],new.dat$Mass[new.dat$Sex=="Male"],pch=19)
points(new.dat$bio1[new.dat$Sex=="Female"],new.dat$Mass[new.dat$Sex=="Female"],pch=19,col="red")

tmp2	<-	gamm(Mass~s(bio1)+s(bio7)+s(bio12)+s(bio18)+s(sc.pcoa),random=list(Sex=~1,Source=~1),data=new.dat)
tmp3	<-	gamm(Mass~s(bio1)+s(bio7)+s(bio12)+s(bio18)+s(sc.pcoa),random=list(Sex=~1),data=new.dat)
tmp4	<-	gamm(Mass~s(bio1)+s(bio12)+s(bio18),random=list(Sex=~1),data=new.dat)
tmp5	<-	gamm(Mass~s(bio1)+s(bio12),random=list(Sex=~1),data=new.dat,na.action=na.omit)
tmp6	<-	gamm(Mass~s(bio12),random=list(Sex=~1),data=new.dat)

anova(tmp2$lme,tmp3$lme,tmp4$lme,tmp5$lme,tmp6$lme)

pred.dat	<-	data.frame(Xbio1=rep(seq(-2.2,1.5,length=100),2),Xbio12=rep(1,200)
							,Sex=factor(c(rep("Male",100),rep("Female",100))))

tmp5.fit	<-	predict(tmp5,newdata=pred.dat,na.action=na.omit,level=0)

plot(pred.dat$bio1[1:100],tmp4.fit[1:100],type="l",ylim=range(new.dat$Mass,na.rm=TRUE))
points(pred.dat$bio1[101:200],tmp4.fit[101:200],type="l",col="red")
points(new.dat$bio1[new.dat$Sex=="Male"],new.dat$Mass[new.dat$Sex=="Male"],pch=19)
points(new.dat$bio1[new.dat$Sex=="Female"],new.dat$Mass[new.dat$Sex=="Female"],pch=19,col="red")

plot(tmp4.fit)
