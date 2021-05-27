## Manacus GAM Models ##
library(lme4)
library(ape)
library(geosphere)
library(nlme)
library(mgcv)
library(MuMIn)

dat	<-	read.csv("~/Dropbox/Tesis Camilo/Data/Data_Manacus.txt")
geo.dist	<-	distm(dat[,c("Long","Lat")])/1000
geo.pcoa	<-	cmdscale(as.dist(geo.dist),k=1)
sc.pcoa		<-	scale(geo.pcoa)

dat$Mass[dat$Sex == ''] <- NA

dat$Body_mass

dat$PC2

new.dat		<-	data.frame(Wing=dat$Wing, Tail=dat$Tail, Culmen=dat$Exposed_Culmen
							,Pc1=dat$PC1, Mass=dat$Body_mass
							,bio1=scale(dat$bio1), bio3=scale(dat$bio3)
							,bio7=scale(dat$bio7), bio12=scale(dat$bio12)
							,bio15=scale(dat$bio15), bio18=scale(dat$bio18)
							,Altitude=scale(dat$Altitude), geog=sc.pcoa[,1]
							,PC1=dat$PC1, PC2=dat$PC2, PC3=dat$PC3
							,Sex=factor(dat$Sex,levels=c("Female","Male")),Source=dat$Source
							,Subspecies=dat$Specie)

## Body Mass ##
#mass.model1	<-	gamm(Mass~s(bio1)+s(bio7)+s(bio12)+s(bio18)+s(sc.pcoa),random=list(Sex=~1,Source=~1),data=new.dat)
#mass.model2	<-	gamm(Mass~s(bio1)+s(bio7)+s(bio12)+s(bio18)+s(sc.pcoa),random=list(Sex=~1),data=new.dat)
#mass.model3	<-	gamm(Mass~s(bio1)+s(bio12)+s(bio18),random=list(Sex=~1),data=new.dat)
#mass.model4	<-	gamm(Mass~s(bio1)+s(bio12),random=list(Sex=~1),data=new.dat,na.action=na.omit)
#mass.model5	<-	gamm(Mass~s(bio12),random=list(Sex=~1),data=new.dat)

mass.model	 <-		lme(Mass~geog+Subspecies,random=list(Sex=~1,Source=~1))
mass.model1	 <- 	lme(Mass~bio1+bio7+bio12+bio15+bio18+geog+Subspecies
						,random=list(Sex=~1, Source=~1),data=new.dat, na.action=na.omit)	
mass.model2  <- 	lme(Mass~bio1+bio7+bio12+bio15+bio18,random=list(Sex=~1),data=new.dat, na.action=na.omit)	
mass.model3  <- 	lme(Mass~bio1+bio7+bio12+bio15,random=list(Sex=~1),data=new.dat,na.action=na.omit)	

AIC(mass.model3)

pred.dat	<-	data.frame(bio1=rep(seq(-2.20,1.50,length=100),2),bio7=rep(0,200), bio12=rep(1,200)
							,bio15=rep(0,200), Sex=factor(c(rep("Male",100),rep("Female",100))))

wing.fit	<-	predict(mass.model3,newdata=pred.dat,na.action=na.omit,level=0, se=T)

wing.fit

plot(pred.dat$bio1,wing.fit$fit,type="l",ylim=range(new.dat$Mass,na.rm=TRUE), 
	xlab='Scaled bio1', ylab='Body Mass Fitted Values')
lines(pred.dat$bio1, wing.fit$fit + wing.fit$se.fit, lty=1, lwd=.5, col='red')
lines(pred.dat$bio1, wing.fit$fit - wing.fit$se.fit, lty=1, lwd=.5, col='red')
points(new.dat$bio1[new.dat$Sex=="Male"],new.dat$Mass[new.dat$Sex=="Male"],pch=19)
points(new.dat$bio1[new.dat$Sex=="Female"],new.dat$Mass[new.dat$Sex=="Female"],pch=19,col="red")
legend('bottomleft', legend=c('Male', 'Female'), col=c('black','red'), pch=19)

wing.fit$se.fit

## Wing ##
wing.model1	 <- 	lme(Wing~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
wing.model2  <- 	lme(Wing~bio1+bio7+bio12+bio15+bio18,random=list(Source=~1),data=new.dat,na.action=na.omit)	
wing.model3  <- 	lme(Wing~bio7+geog,random=list(Source=~1),data=new.dat,na.action=na.omit)	

AIC(wing.model1)
## mass.model1	<-	gamm(Wing~s(bio1)+s(bio7)+s(bio12)+s(bio18)+s(sc.pcoa),random=list(Sex=~1,Source=~1),data=new.dat)
## mass.model2	<-	gamm(Wing~s(bio1)+s(bio7)+s(bio12)+s(bio18)+s(sc.pcoa),random=list(Sex=~1),data=new.dat)
## mass.model3	<-	gamm(Wing~s(bio1)+s(bio12)+s(bio18),random=list(Sex=~1),data=new.dat)
## mass.model4	<-	gamm(Wing~s(bio1)+s(bio12),random=list(Sex=~1),data=new.dat,na.action=na.omit)
## mass.model5	<-	gamm(Wing~s(bio12),random=list(Sex=~1),data=new.dat)

anova(mass.model1$lme, mass.model2$lme, mass.model3$lme,
	mass.model4$lme, mass.model5$lme)

summary(mass.model1)
r.squaredGLMM(mass.model1)

## Tail ##
tail.model1	 <- 	lme(Tail~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
tail.model2	 <- 	lme(Tail~bio1+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
tail.model3	 <- 	lme(Tail~bio15,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	

tail.model3.gamm	<-	gamm(Tail~s(bio1)+s(bio12)+s(bio15)+s(bio18)+s(geog),random=list(Sex=~1,Source=~1),data=new.dat, family=quasipoisson)
tail.model4.gamm	<-	gamm(Tail~s(bio1)+s(bio12)+s(bio15),random=list(Sex=~1,Source=~1),data=new.dat, family=quasipoisson)

AIC(tail.model4.gamm)
summary(tail.model4.gamm$gam)

AIC(tail.model3)
summary(tail.model3)
anova(tail.model3)

range(new.dat$bio1)

pred.dat	<-	data.frame(bio15=rep(seq(-3.4,4.2,length=100),2),
	Source=factor(c(rep("Museum",100),rep("Net",100))),
	Sex=factor(c(rep("Female",100),rep("Male",100))))

pred.dat	<-	data.frame(bio1=rep(seq(-2.2,1.5,length=100),2),
	bio12=rep(1,200), bio15=rep(1,200),
	Source=factor(c(rep("Museum",100),rep("Net",100))),
	Sex=factor(c(rep("Female",100),rep("Male",100))))

tail.fit	<-	predict(tail.model3,newdata=pred.dat,na.action=na.omit,level=0, se=T)
tail.fit	<-	predict(tail.model4.gamm$gam,newdata=pred.dat,na.action=na.omit, se=T, type='response')
tail.fit$fit

plot(pred.dat$bio15[1:100],tail.fit[1:100],type="l",ylim=range(new.dat$Tail,na.rm=TRUE), 
	xlab='Scaled bio15', ylab='Observed Tail Length', lwd=2)
lines(pred.dat$bio15[101:200],tail.fit[101:200],col='blue', lwd=2)
points(new.dat$bio15[new.dat$Sex=="Male"],new.dat$Tail[new.dat$Sex=="Male"],pch=19)
points(new.dat$bio15[new.dat$Sex=="Female"],new.dat$Tail[new.dat$Sex=="Female"],pch=19, col='blue')
legend('bottomright', legend=c('Male','Female'), col=c('black', 'blue'), pch=19)
r.squaredGLMM(tail.model3)

## plot(pred.dat$bio1[1:100],tail.fit$fit[1:100],type="l",ylim=range(new.dat$Tail,na.rm=TRUE), 
## 	xlab='Scaled bio15', ylab='Observed Tail Length', lwd=2)
## lines(pred.dat$bio1[101:200],tail.fit$fit[101:200],col='blue', lwd=2)
## points(new.dat$bio1[new.dat$Sex=="Male"],new.dat$Tail[new.dat$Sex=="Male"],pch=19)
## points(new.dat$bio1[new.dat$Sex=="Female"],new.dat$Tail[new.dat$Sex=="Female"],pch=19, col='blue')
## legend('bottomright', legend=c('Male','Female'), col=c('black', 'blue'), pch=19)


summary(tail.model3)
##### Culmen #####
culmen.model1	 <- 	lme(Culmen~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
culmen.model2	 <- 	lme(Culmen~bio12,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit) ## AIC reduction but non significant term	


summary(culmen.model1)
anova(culmen.model2)
AIC(culmen.model1)
r.squaredGLMM(culmen.model2)



##### Pc1 ####
pc1.model1	 <- 	lme(Pc1~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
pc1.model2	 <- 	lme(Pc1~bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
pc1.model3	 <- 	lme(Pc1~bio7+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
pc1.model4	 <- 	lme(Pc1~bio7,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	

## Model Selection ##
anova(pc1.model4)
summary(pc1.model4)
AIC(pc1.model4)
r.squaredGLMM(pc1.model4)

pred.dat	<-	data.frame(bio7=rep(seq(-4.5, 2.7, length=100),2),
	Source=factor(c(rep("Net",100),rep("Museum",100))),
	Sex=factor(c(rep("Female",100),rep("Male",100))))

pc1.fit	<-	predict(pc1.model4,newdata=pred.dat,na.action=na.omit,level=0, se=T)

plot(pred.dat$bio7[1:100],pc1.fit$fit[1:100],type="l",ylim=range(new.dat$Pc1,na.rm=TRUE), 
	xlab='Scaled bio7', ylab='Observed Pc1 Length')
lines(pred.dat$bio7[101:200],pc1.fit$fit[101:200],col='red')
points(new.dat$bio7[new.dat$Source=="Net"],new.dat$Pc1[new.dat$Source=="Net"],pch=19)
points(new.dat$bio7[new.dat$Source=="Museum"],new.dat$Pc1[new.dat$Source=="Museum"],pch=19, col='blue')
lines(pred.dat$bio7, pc1.fit$fit + pc1.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio7, pc1.fit$fit - pc1.fit$se.fit, lty=2, lwd=2)
legend('bottomright', legend=c('Net','Museum'), col=c('black', 'blue'), pch=19)

r.squaredGLMM(pc1.model4)

## PC2 ##
PC2.model1	 <- 	lme(PC2~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
PC2.model2	 <- 	lme(PC2~bio1+bio7+bio15+bio18,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
PC2.model3	 <- 	lme(PC2~bio1+bio7+bio15,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	

anova(PC2.model3)
summary(PC2.model)
AIC(PC2.model3)
r.squaredGLMM(PC2.model3)

range(new.dat$bio1)
pred.dat	<-	data.frame(bio1=rep(seq(-2.3, 1.5, length=100),2),
						bio7=rep(0,200), bio15=rep(0,200),
						Source=factor(c(rep("Net",100),rep("Museum",100))),
						Sex=factor(c(rep("Female",100),rep("Male",100))))

PC2.fit	<-	predict(PC2.model3,newdata=pred.dat,na.action=na.omit,level=0, se=T)

plot(pred.dat$bio1[1:100],PC2.fit$fit[1:100],type="l",ylim=range(new.dat$PC2,na.rm=TRUE), 
	xlab='Scaled bio1', ylab='PC2 Value')
lines(pred.dat$bio1[101:200],PC2.fit$fit[101:200],col='red')
points(new.dat$bio1[new.dat$Source=="Net"],new.dat$PC2[new.dat$Source=="Net"],pch=19)
points(new.dat$bio1[new.dat$Source=="Museum"],new.dat$PC2[new.dat$Source=="Museum"],pch=19, col='blue')
lines(pred.dat$bio1, PC2.fit$fit + PC2.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio1, PC2.fit$fit - PC2.fit$se.fit, lty=2, lwd=2)
legend('bottomright', legend=c('Net','Museum'), col=c('black', 'blue'), pch=19)

#################### PC2 ######################

PC2.model1	 <- 	lme(PC2~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
PC2.model2	 <- 	lme(PC2~bio7+bio12+geog,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	
PC2.model3	 <- 	lme(PC2~bio7,random=list(Sex=~1, Source=~1),data=new.dat,na.action=na.omit)	

anova(PC2.model3)
summary(PC2.model)
AIC(PC2.model3)
r.squaredGLMM(PC2.model3)

range(new.dat$bio7)
pred.dat	<-	data.frame(bio7=rep(seq(-4.5, 2.65, length=100),2),
						Source=factor(c(rep("Net",100),rep("Museum",100))),
						Sex=factor(c(rep("Female",100),rep("Male",100))))

PC2.fit	<-	predict(PC2.model3,newdata=pred.dat,na.action=na.omit,level=0, se=T)

plot(pred.dat$bio7[1:100],PC2.fit$fit[1:100],type="l",ylim=range(new.dat$PC2,na.rm=TRUE), 
	xlab='Scaled bio7', ylab='PC2 Value')
lines(pred.dat$bio7[101:200],PC2.fit$fit[101:200],col='red')
points(new.dat$bio7[new.dat$Source=="Net"],new.dat$PC2[new.dat$Source=="Net"],pch=19)
points(new.dat$bio7[new.dat$Source=="Museum"],new.dat$PC2[new.dat$Source=="Museum"],pch=19, col='blue')
lines(pred.dat$bio7, PC2.fit$fit + PC2.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio7, PC2.fit$fit - PC2.fit$se.fit, lty=2, lwd=2)
legend('bottomright', legend=c('Net','Museum'), col=c('black', 'blue'), pch=19)


################################## MODELS BY SUBSPECIES #############################################
############# M. manacus manacus Models ############################

## Body Mass mm ##
new.dat.mm <- new.dat[new.dat$Subspecies=='Manacus manacus',]

mass.model1	 <- 	lme(Mass~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.mm, na.action=na.omit)	
mass.model2  <- 	lme(Mass~bio1+bio7+bio12+bio15+bio18,random=list(Sex=~1),data=new.dat.mm, na.action=na.omit)	
mass.model3  <- 	lme(Mass~bio1+bio7+bio12,random=list(Sex=~1),data=new.dat.mm,na.action=na.omit)	

anova(mass.model3)
AIC(mass.model3)

pred.dat	<-	data.frame(bio1=rep(seq(-2.20,1.50,length=100),2),bio7=rep(1,200), bio12=rep(1,200)
							,Sex=factor(c(rep("Male",100),rep("Female",100))))

mass.fit	<-	predict(mass.model3,newdata=pred.dat,na.action=na.omit,level=0, se=T)

plot(pred.dat$bio1,mass.fit$fit,type="l",ylim=range(new.dat.mm$Mass,na.rm=TRUE), 
	xlab='Scaled bio1', ylab='Body Mass Fitted Values')
lines(pred.dat$bio1, mass.fit$fit + mass.fit$se.fit, lty=1, lwd=.5, col='red')
lines(pred.dat$bio1, mass.fit$fit - mass.fit$se.fit, lty=1, lwd=.5, col='red')
points(new.dat$bio1[new.dat$Sex=="Male"],new.dat$Mass[new.dat$Sex=="Male"],pch=19)
points(new.dat$bio1[new.dat$Sex=="Female"],new.dat$Mass[new.dat$Sex=="Female"],pch=19,col="red")
legend('bottomleft', legend=c('Male', 'Female'), col=c('black','red'), pch=19)

## Wing mm ##
wing.model1	 <- 	lme(Wing~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)	
wing.model2  <- 	lme(Wing~bio1+bio12+bio15+bio18,random=list(Sex=~1,Source=~1),data=new.dat.mm,na.action=na.omit)	
wing.model3  <- 	lme(Wing~bio15+bio18,random=list(Sex=~1,Source=~1),data=new.dat.mm,na.action=na.omit)	


AIC(wing.model3)
anova(wing.model3)
summary(wing.model3)

## mass.model1	<-	gamm(Wing~s(bio1)+s(bio7)+s(bio12)+s(bio18)+s(sc.pcoa),random=list(Sex=~1,Source=~1),data=new.dat)
## mass.model2	<-	gamm(Wing~s(bio1)+s(bio7)+s(bio12)+s(bio18)+s(sc.pcoa),random=list(Sex=~1),data=new.dat)
## mass.model3	<-	gamm(Wing~s(bio1)+s(bio12)+s(bio18),random=list(Sex=~1),data=new.dat)
## mass.model4	<-	gamm(Wing~s(bio1)+s(bio12),random=list(Sex=~1),data=new.dat,na.action=na.omit)
## mass.model5	<-	gamm(Wing~s(bio12),random=list(Sex=~1),data=new.dat)


###### Tail mm ###########
tail.model1	 <- 	lme(Tail~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)	
tail.model2	 <- 	lme(Tail~bio1+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)	
tail.model3	 <- 	lme(Tail~bio12+bio15,random=list(Sex=~1,Source=~1),data=new.dat.mm,na.action=na.omit)	

summary(tail.model3)
AIC(tail.model3)
range(new.dat.mm$bio15)

pred.dat <- data.frame(bio12=rep(seq(-1.5,1.7,length=100),2),
	bio15=rep(seq(-3.2,4.2,length=100),2),
	Source=factor(c(rep("Net",100),rep("Museum",100))),
	Sex=factor(c(rep("Male",100),rep("Female",100))))

tail.fit	<-	predict(tail.model3,se=T, newdata=pred.dat, na.action=na.omit, level=0)


plot(pred.dat$bio15,tail.fit$fit,type="l",ylim=range(new.dat.mm$Tail,na.rm=TRUE), 
	xlab='Scaled bio15', ylab='Observed Tail values')
points(pred.dat$bio15[101:200],tail.fit$fit[101:200],type="l",col="red")
lines(pred.dat$bio15, tail.fit$fit + tail.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio15[1:100], tail.fit$fit[1:100] - tail.fit$se.fit[1:100], lty=2, lwd=2)
points(new.dat$bio15[new.dat$Source=="Net"],new.dat$Tail[new.dat$Source=="Net"],pch=19)
points(new.dat$bio15[new.dat$Source=="Museum"],new.dat$Tail[new.dat$Source=="Museum"],pch=19,col="red")
legend('bottomleft', legend=c('Net', 'Museum'), col=c('black','red'), pch=19)

##### Culmen mm #####
culmen.model1	 <- 	lme(Culmen~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)
culmen.model2	 <- 	lme(Culmen~bio18,random=list(Sex=~1,Source=~1),data=new.dat.mm,na.action=na.omit)

anova(culmen.model2)
AIC(culmen.model2)

## Pc1 ##
pc1.model1	 <- 	lme(Pc1~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)
summary(pc1.model1)
AIC(pc1.model1)

## PC2 ##
PC2.model1	 <- 	lme(PC2~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)
PC2.model2	 <- 	lme(PC2~bio1+bio12+bio15+bio18,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)

anova(PC2.model2)
AIC(PC2.model2)

range(new.dat.mm$bio18)
pred.dat <- data.frame(bio1=rep(seq(-2.2,1.5,length=100),2),
	bio12=rep(1,200),
	bio15=rep(1,200),
	bio18=rep(1,200),
	Source=factor(c(rep("Net",100),rep("Museum",100))),
	Sex=factor(c(rep("Male",100),rep("Female",100))))

PC2.fit	<-	predict(PC2.model2, se=T, newdata=pred.dat, na.action=na.omit, level=0)


plot(pred.dat$bio1,PC2.fit$fit,type="l",ylim=range(new.dat.mm$PC2,na.rm=TRUE), 
	xlab='Scaled bio1', ylab='Observed PC2 values')
points(pred.dat$bio1[101:200],PC2.fit$fit[101:200],type="l",col="red")
lines(pred.dat$bio1, PC2.fit$fit + PC2.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio1[1:100], PC2.fit$fit[1:100] - PC2.fit$se.fit[1:100], lty=2, lwd=2)
	points(new.dat.mm$bio1[new.dat$Source=="Net"],new.dat$PC2[new.dat$Source=="Net"],pch=19)
	points(new.dat.mm$bio1[new.dat$Source=="Museum"],new.dat$PC2[new.dat$Source=="Museum"],pch=19,col="red")
legend('bottomleft', legend=c('Net', 'Museum'), col=c('black','red'), pch=19)

## PC2 ####
PC2.model1	 <- 	lme(PC2~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)
PC2.model2	 <- 	lme(PC2~bio1+bio12+bio15+bio18,random=list(Sex=~1, Source=~1),data=new.dat.mm,na.action=na.omit)

anova(PC2.model1)
AIC(PC2.model2)


### Models for Manacus Vitellinus ##
## Body Mass #
new.dat.vt <- new.dat[new.dat$Subspecies=='Manacus vitellinus',]

mass.model1	 <- 	lme(Mass~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.vt, na.action=na.omit)	
mass.model1 <- lm(Mass~bio1+bio7+bio12+bio15+bio18+geog,data=new.dat.vt, na.action=na.omit)
## No significant Terms ###
drop1(mass.model1)
anova(mass.model1)


## Wing ##
wing.model1	 <- 	lme(Wing~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)	
wing.model2	 <- 	lme(Wing~bio7,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)	

anova(wing.model2)
AIC(wing.model2)

range(new.dat.vt$bio7)
pred.dat <- data.frame(bio7=rep(seq(-2.35,1.3,length=100),2),
	Source=factor(c(rep("Net",100),rep("Museum",100))),
	Sex=factor(c(rep("Male",100),rep("Female",100))))

wing.fit	<-	predict(wing.model2, se=T, newdata=pred.dat, na.action=na.omit, level=0)


plot(pred.dat$bio7,wing.fit$fit,type="l",ylim=range(new.dat.vt$Wing,na.rm=TRUE), 
	xlab='Scaled bio7', ylab='Observed wing values')
points(pred.dat$bio7[101:200],wing.fit$fit[101:200],type="l",col="red")
lines(pred.dat$bio7, wing.fit$fit + wing.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio7[1:100], wing.fit$fit[1:100] - wing.fit$se.fit[1:100], lty=2, lwd=2)
points(new.dat.vt$bio7[new.dat.vt$Source=="Net"],new.dat.vt$Wing[new.dat.vt$Source=="Net"],pch=19)
points(new.dat.vt$bio7[new.dat.vt$Source=="Museum"],new.dat.vt$Wing[new.dat.vt$Source=="Museum"],pch=19,col="red")
legend('bottomleft', legend=c('Net', 'Museum'), col=c('black','red'), pch=19)


## Tail ##
tail.model1	 <- 	lme(Tail~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)	
summary(tail.model1)

# Culmen #
culmen.model1	 <- 	lme(Culmen~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)
culmen.model2	 <- 	lme(Culmen~bio12,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)

anova(culmen.model1)
AIC(culmen.model2)


range(new.dat.vt$bio12)
pred.dat <- data.frame(bio12=rep(seq(-1.4,4.6,length=100),2),
	Source=factor(c(rep("Net",100),rep("Museum",100))),
	Sex=factor(c(rep("Male",100),rep("Female",100))))

culmen.fit	<-	predict(culmen.model2, se=T, newdata=pred.dat, na.action=na.omit, level=0)

plot(pred.dat$bio12,culmen.fit$fit,type="l",ylim=range(new.dat.vt$Culmen,na.rm=TRUE), 
	xlab='Scaled Geog', ylab='Observed culmen values')
points(pred.dat$bio12[101:200],culmen.fit$fit[101:200],type="l",col="red")
lines(pred.dat$bio12, culmen.fit$fit + culmen.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio12[1:100], culmen.fit$fit[1:100] - culmen.fit$se.fit[1:100], lty=2, lwd=2)
points(new.dat.vt$bio12[new.dat.vt$Source=="Net"],new.dat.vt$Culmen[new.dat.vt$Source=="Net"],pch=19)
points(new.dat.vt$bio12[new.dat.vt$Source=="Museum"],new.dat.vt$Culmen[new.dat.vt$Source=="Museum"],pch=19,col="red")
legend('bottomleft', legend=c('Net', 'Museum'), col=c('black','red'), pch=19)

## Pc1 ##
pc1.model1	 <- 	lme(Pc1~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)
pc1.model2	 <- 	lme(Pc1~bio12,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)

anova(pc1.model2)
AIC(pc1.model2)

pred.dat <- data.frame(bio12=rep(seq(-1.4,4.6,length=100),2),
	Source=factor(c(rep("Net",100),rep("Museum",100))),
	Sex=factor(c(rep("Male",100),rep("Female",100))))

pc1.fit	<-	predict(pc1.model2, se=T, newdata=pred.dat, na.action=na.omit, level=0)

plot(pred.dat$bio12,pc1.fit$fit,type="l",ylim=range(new.dat.vt$Pc1,na.rm=TRUE), 
	xlab='Scaled Geog', ylab='Observed pc1 values')
points(pred.dat$bio12[101:200],pc1.fit$fit[101:200],type="l",col="red")
lines(pred.dat$bio12, pc1.fit$fit + pc1.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio12[1:100], pc1.fit$fit[1:100] - pc1.fit$se.fit[1:100], lty=2, lwd=2)
points(new.dat.vt$bio12[new.dat.vt$Source=="Net"],new.dat.vt$Pc1[new.dat.vt$Source=="Net"],pch=19)
points(new.dat.vt$bio12[new.dat.vt$Source=="Museum"],new.dat.vt$Pc1[new.dat.vt$Source=="Museum"],pch=19,col="red")
legend('bottomleft', legend=c('Net', 'Museum'), col=c('black','red'), pch=19)


## PC2 ##
PC2.model1	 <- 	lme(PC2~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)
PC2.model2	 <- 	lme(PC2~bio1,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)

anova(PC2.model2)
AIC(PC2.model2)
range(new.dat.vt$bio1)

pred.dat <- data.frame(bio1=rep(seq(-2.2,1,length=100),2),
	Source=factor(c(rep("Net",100),rep("Museum",100))),
	Sex=factor(c(rep("Male",100),rep("Female",100))))

PC2.fit	<-	predict(PC2.model2, se=T, newdata=pred.dat, na.action=na.omit, level=0)

plot(pred.dat$bio1,PC2.fit$fit,type="l",ylim=range(new.dat.vt$PC2,na.rm=TRUE), 
	xlab='Scaled Bio1', ylab='Observed PC2 values')
points(pred.dat$bio1[101:200],PC2.fit$fit[101:200],type="l",col="red")
lines(pred.dat$bio1, PC2.fit$fit + PC2.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio1[1:100], PC2.fit$fit[1:100] - PC2.fit$se.fit[1:100], lty=2, lwd=2)
points(new.dat.vt$bio1[new.dat.vt$Source=="Net"],new.dat.vt$PC2[new.dat.vt$Source=="Net"],pch=19)
points(new.dat.vt$bio1[new.dat.vt$Source=="Museum"],new.dat.vt$PC2[new.dat.vt$Source=="Museum"],pch=19,col="red")
legend('bottomleft', legend=c('Net', 'Museum'), col=c('black','red'), pch=19)

## PC2 Vitellinus ##
PC2.model1	 <- 	lme(PC2~bio1+bio7+bio12+bio15+bio18+geog,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)
PC2.model2	 <- 	lme(PC2~bio12,random=list(Sex=~1, Source=~1),data=new.dat.vt,na.action=na.omit)

anova(PC2.model2)
AIC(PC2.model2)

range(new.dat.vt$bio12)

pred.dat <- data.frame(bio12=rep(seq(-1.4,4.6,length=100),2),
	Source=factor(c(rep("Net",100),rep("Museum",100))),
	Sex=factor(c(rep("Male",100),rep("Female",100))))

PC2.fit	<-	predict(PC2.model2, se=T, newdata=pred.dat, na.action=na.omit, level=0)

plot(pred.dat$bio12,PC2.fit$fit,type="l",ylim=range(new.dat.vt$PC2,na.rm=TRUE), 
	xlab='Scaled Bio1', ylab='Observed PC2 values')
points(pred.dat$bio12[101:200],PC2.fit$fit[101:200],type="l",col="red")
lines(pred.dat$bio12, PC2.fit$fit + PC2.fit$se.fit, lty=2, lwd=2)
lines(pred.dat$bio12[1:100], PC2.fit$fit[1:100] - PC2.fit$se.fit[1:100], lty=2, lwd=2)
points(new.dat.vt$bio12[new.dat.vt$Source=="Net"],new.dat.vt$PC2[new.dat.vt$Source=="Net"],pch=19)
points(new.dat.vt$bio12[new.dat.vt$Source=="Museum"],new.dat.vt$PC2[new.dat.vt$Source=="Museum"],pch=19,col="red")
legend('bottomleft', legend=c('Net', 'Museum'), col=c('black','red'), pch=19)


