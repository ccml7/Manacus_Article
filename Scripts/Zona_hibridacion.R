logitm <- function(b,c,d,e,x){(c+(d-c)/(1+exp(b*(log(x)-log(e)))))}

#Funcion para calcular la probabilidad del modelo
logitm.lik	<-	function(pars,x.var,y.var){
	b = pars[1]
	c = pars[2]
	d = pars[3]
	e = exp(pars[4])
	x = x.var
	y =	y.var
	y.hat	<-	logitm(b,c,d,e,x)
	sum.sq	<-	sum((y-y.hat)^2,na.rm=TRUE)
	return(sum.sq)
}

### AIC for least squares model
AIC.lssq <- function(n,p,trss){
	
	aic <- n*(log((2*pi*trss)/n)+1) + 2*(p+1)
	return(aic)
}

logitm.ls	<-	function(pars,x.var,y.var,CI=TRUE,reps=1000,method="BFGS",control=NULL,prediction=NULL){
	if(missing(pars)){pars=rep(1,4)}
	mle	<-	optim(pars,logitm.lik,method=method,x.var=x.var,y.var=y.var,control=control)
	b	<-	mle$par[1]
	c	<-	mle$par[2]
	d	<-	mle$par[3]
	e	<-	exp(mle$par[4])
	n	<-	length(y.var)
	p	<-	length(mle$par)
	trss<-	mle$value
	aic	<-	AIC.lssq(n,p,trss)
	SSq	<-	sum((y.var-mean(y.var,na.rm=TRUE))^2,na.rm=TRUE)
	r2	<-	1-trss/SSq
	if(CI==TRUE){
		params	<-	matrix(NA,nrow=reps,ncol=length(pars))
		preds	<-	matrix(NA,nrow=length(x.var),ncol=reps)
		preds1	<-	matrix(NA,nrow=length(prediction),ncol=reps)	
		predicted		<-	logitm(b,c,d,e,x.var)
		predicted1		<-	logitm(b,c,d,e,prediction)
		for(i in 1:reps){
			
			new.dat			<-	predicted	+	rnorm(length(x.var),0,sqrt(SSq/length(x.var)))
			params[i,]		<-	optim(par=c(b=1,c=1,d=1,e=1),logitm.lik,x.var=x.var,y.var=new.dat
									,control=control)$par
			preds[,i]		<-	logitm(params[i,][1],params[i,][2],params[i,][3],exp(params[i,][4]),x.var)
			if(!is.null(prediction)){preds1[,i]	<-	logitm(params[i,][1],params[i,][2],params[i,][3],exp(params[i,][4]),prediction)}
			
		}

		ci.b	<-	quantile(params[,1],probs=c(0.025,0.975))
		ci.c	<-	quantile(params[,2],probs=c(0.025,0.975))
		ci.d	<-	quantile(params[,3],probs=c(0.025,0.975))
		ci.e	<-	quantile(params[,4],probs=c(0.025,0.975))
		
		fit.ci		<-	apply(preds,1,quantile,probs=c(0.025,0.975))
		fit			<-	cbind(predicted,t(fit.ci))
		colnames(fit)	<-	c("fitted","2.5%","97.5%")
		if(!is.null(prediction)){
		pred.ci		<-	apply(preds1,1,quantile,probs=c(0.025,0.975))
		Predicted	<-	cbind(predicted1,t(pred.ci))
		colnames(Predicted)	<-	c("fitted","2.5%","97.5%")}else{Predicted<-NA}			
		parameters	<-	matrix(c(b,ci.b,c,ci.c,d,ci.d,exp(e),exp(ci.e)),ncol=3,byrow=T
						,dimnames=list(c("b","c","d","e"),c("MLE","2.5%","97.5%")))
		results	<-	list(par=parameters,AIC=aic,R2=r2,Residual.ssq=trss,Total.ssq=SSq,Fitted=fit,Predicted=Predicted)
	}
	
	else{
		parameters	<-	c(b,c,d,exp(e))
		names(parameters)	<-	letters[2:5]
		
		results	<-	list(par=parameters,AIC=aic,R2=r2,Residual.ssq=trss,Total.ssq=SSq)
		
	}
	
	return(results)
}

drm.boot	<-	function(object,newdata,reps=1000){
	
	x.var	<-	object$data[,1]
	if(missing(newdata)){pred.x	<- x.var}else{pred.x <- newdata}
	
	drm.coefs	<-	coef(object)
	b	<-	drm.coefs[1]
	c	<-	drm.coefs[2]
	d	<-	drm.coefs[3]
	e	<-	drm.coefs[4]
	
	mod.SSq	<-	sum(object$predres[,2]^2)	
	preds	<-	matrix(NA,nrow=length(pred.x),ncol=reps)
	params	<-	matrix(NA,nrow=reps,ncol=4)
	predicted		<-	logitm(b,c,d,e,x.var)
	
	for(i in 1:reps){
	
		new.dat			<-	predicted	+	rnorm(length(x.var),0,sqrt(mod.SSq/length(x.var)))
		mod				<-	tryCatch(drm(new.dat~x.var,fct=LL.4()),error=function(e)NULL)
		if(is.null(mod)){
			params[i,]	<-	rep(NA,4)
			preds[,i]	<-	NA
			}else{
				params[i,]	<-	coef(mod)
				preds[,i]	<-	logitm(params[i,1],params[i,2],params[i,3],params[i,4],pred.x)						
			}

	}

	prediction		<-	t(apply(preds,1,quantile,probs=c(0.5,0.025,0.975),na.rm=TRUE))
	parameters		<-	t(apply(params,2,quantile,probs=c(0.025,0.975),na.rm=TRUE))
	parameters		<-	cbind(drm.coefs,parameters)
	eff.boot		<-	reps-length(which(is.na(params[,1])))
	colnames(prediction)	<-colnames(parameters)	<-	c("Fitted","2.5%","97.5%")
	results	<-	list(Parameters=parameters,Prediction=prediction, Effective.boot=eff.boot)
	return(results)
	
}


ML.logitm.lik	<-	optim(par=c(a=1,b=1,cent=1),logitm.lik
							,x.var=Distance2West,y.var=BodyMass)

# Exploracion de la zona de hibridacion
library(lme4)
library(ape)
library(geosphere)
library(nlme)
library(mgcv)
library(MuMIn)
library(raster)

dat				<-	read.csv("~/Dropbox/Tesis Camilo/Data/Data_Manacus.csv")
clim.vars.path	<-	list.files("~/Documents/GIS/bioclim",full.names=TRUE,pattern=".bil")
elev.path		<-	"/Users/JuanPabloGomez/Documents/GIS/alt_30s_bil/alt.bil"
clim.vars		<-	stack(clim.vars.path)
elev			<-	raster(elev.path)
colombia		<-	readOGR("/Users/JuanPabloGomez/Documents/GIS/COL_adm/","COL_adm0")

clim.dat		<-	extract(clim.vars,dat[,c("Long","Lat")])
elev.dat		<-	extract(elev,dat[,c("Long","Lat")])

clim.dat		<-	cbind(clim.dat,elev.dat)

clim.pca		<-	prcomp(~.,data=data.frame(clim.dat),scale.=TRUE,center=TRUE)
morph.pca		<-	prcomp(~Wing+Tail+Exposed_Culmen+Tarsus,scale.=TRUE,center=TRUE,data=dat)

pops		<-	read.delim("~/Dropbox/Tesis Camilo/Data/Populations_CZ.txt")

ConZone.pop1<-	pops[pops$Population=='Acandi', c('Longitude', 'Latitude')]


#ConZone.dist<-	distm(rbind(ConZone.pop1,pops[,c("Longitude","Latitude")]))
#ConZone.dist<-	ConZone.dist[2:21,1]/1000
#names(ConZone.dist)	<- pops$Population

ConZone.dist<-	distm(pops[,c("Longitude","Latitude")])
ConZone.dist<-	ConZone.dist[1,]/1000
names(ConZone.dist)	<- pops$Population

ConZone.inds	<-	which(!is.na(match(dat$Population,pops$Population)))	
ConZone.dist
morph.pca.x		<-	morph.pca$x[ConZone.inds,]
clim.pca.x		<-	clim.pca$x[ConZone.inds,]
BodyMass		<-	dat$Body_mass[ConZone.inds]
dat.Zone		<-	dat[ConZone.inds,]

Distance2West	<-	ConZone.dist[match(dat.Zone$Population,names(ConZone.dist))]


dat.Zone$Population	<-	as.factor(as.character(dat.Zone$Population))

summary(lm(BodyMass~clim.pca.x[,1]+clim.pca.x[,2]+clim.pca.x[,3]+clim.pca.x[,4]+Distance2West))
summary(lm(morph.pca.x[,1]~clim.pca.x[,1]+clim.pca.x[,2]+clim.pca.x[,3]+clim.pca.x[,4]+Distance2West))
summary(lm(morph.pca.x[,2]~clim.pca.x[,1]+clim.pca.x[,2]+clim.pca.x[,3]+clim.pca.x[,4]+Distance2West))
summary(lm(morph.pca.x[,3]~clim.pca.x[,1]+clim.pca.x[,2]+clim.pca.x[,3]+clim.pca.x[,4]+Distance2West))

par(mfrow=c(2,2))
plot(Distance2West,BodyMass,pch=19)
abline(v=193.4,col="red")
plot(Distance2West,morph.pca.x[,1],pch=19)
abline(v=193.4,col="red")
plot(Distance2West,morph.pca.x[,2],pch=19)
abline(v=193.4,col="red")
plot(Distance2West,morph.pca.x[,3],pch=19)
abline(v=193.4,col="red")

plot(clim.pca.x[,1],BodyMass,pch=19)

plot(clim.pca.x[,2],morph.pca.x[,1],pch=19)

plot(Distance2West,morph.pca.x[,2],pch=19)
abline(v=193.4,col="red")
plot(dat.Zone$Distance2West,dat.Zone$Body_mass,pch=19)
abline(v=193.4,col="red")

sampling.size	<-	table(dat.Zone$Population,dat.Zone$Specie)[order(ConZone.dist),]
sampling.size

all.locs	<-	table(dat.Zone$Population,dat.Zone$Specie)

pops$ConZone.dist <- ConZone.dist
pops <- data.frame(pops, ConZone.dist)
names(pops)
require(hzar)

input <- pops[,c(1,28,12:17,2)]

manaPC1 <- hzar.doNormalData1DPops(distance=input$ConZone.dist,
								   muObs=input$PC1,
								   nEff=input$n,
								   varObs=input$PC1_SD,
								   siteID=input$Population)

manaPC2 <- hzar.doNormalData1DPops(distance=input$ConZone.dist,
								   muObs=input$PC2,
								   nEff=input$n,
								   varObs=input$PC2_SD,
								   siteID=input$Population)

manaPC3 <- hzar.doNormalData1DPops(distance=input$ConZone.dist,
								   muObs=input$PC3,
								   nEff=input$n,
								   varObs=input$PC3_SD,
								   siteID=input$Population)


cline.PC1 <- hzar.makeCline1DNormal(data=manaPC1)
PC1.model <- hzar.first.fitRequest.gC(gModel=cline.PC1, obsData=manaPC1, verbose=T)
PC1.fit <- hzar.doFit(PC1.model)
PC1.ModelData <- hzar.dataGroup.add(PC1.fit)
PC1.cline <- hzar.make.obsDataGroup(PC1.fit)
hzar.plot.cline(PC1.cline)

cline.PC2 <- hzar.makeCline1DNormal(data=manaPC2)
PC2.model <- hzar.first.fitRequest.gC(gModel=cline.PC2, obsData=manaPC2, verbose=T)
PC2.fit <- hzar.doFit(PC2.model)
PC2.cline <- hzar.make.obsDataGroup(PC2.fit)
hzar.plot.cline(PC2.cline)

cline.PC3 <- hzar.makeCline1DNormal(data=manaPC3)
PC3.model <- hzar.first.fitRequest.gC(gModel=cline.PC3, obsData=manaPC3, verbose=T)
PC3.fit <- hzar.doFit(PC3.model)
PC3.cline <- hzar.make.obsDataGroup(PC3.fit)
hzar.plot.cline(PC3.cline)


names(pops)

ML.logitm.lik	<-	optim(par=c(a=1,b=1,cent=1),logitm.lik
							,x.var=Distance2West,y.var=BodyMass
							,control=list(maxit=1500,parscale=c(0.1,100,1000)))

#####################################################################################
require(drc)

ConZone.dist<-	distm(pops[,c("Longitude","Latitude")])
ConZone.dist<-	ConZone.dist[1,]/1000
names(ConZone.dist)	<- pops$Population

pops <- data.frame(pops, ConZone.dist)

pops <- pops[order(pops$ConZone.dist),]
PC1.log <- drm(PC1~ConZone.dist, data=pops, fct=LL.4())
PC2.log <- drm(PC2~ConZone.dist, data=pops, fct=LL.4())
PC3.log <- drm(PC3~ConZone.dist, data=pops, fct=LL.4())

summary(PC1.log)
summary(PC2.log)
summary(PC3.log)

fit.PC1 <- predict(PC1.log, pred.data=seq(0,range(pops$ConZone.dist)[2]), se.fit=TRUE)
fit.PC2 <- predict(PC2.log, pred.data=seq(0,range(pops$ConZone.dist)[2]), se.fit=TRUE)
fit.PC3 <- predict(PC3.log, pred.data=seq(0,range(pops$ConZone.dist)[2]), se.fit=TRUE)

grey <- rgb(150,150,150, alpha=200, maxColorValue=255)

png('/home/camilo/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v2/Cz_PC.png', units='cm', height = 10, width=25, res=150)
par(mfrow=c(1,3), family='serif')
plot(PC1~ConZone.dist, data=pops, col=grey, xlab='')
points(pops$ConZone.dist[pops$Specie == 'Manacus manacus'], pops$PC1[pops$Specie == 'Manacus manacus'], col='black', pch=19)
points(pops$ConZone.dist[pops$Specie == 'Manacus vitellinus'], pops$PC1[pops$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(pops$ConZone.dist, fit.PC1[,1], col='black')
lines(pops$ConZone.dist, fit.PC1[,1]+fit.PC1[,2], col='blue')
lines(pops$ConZone.dist, fit.PC1[,1]-fit.PC1[,2], col='blue')
abline(v=193.4,col="red", lty='dashed')


plot(PC2~ConZone.dist, data=pops, col=grey, xlab='')
points(pops$ConZone.dist[pops$Specie == 'Manacus manacus'], pops$PC2[pops$Specie == 'Manacus manacus'], col='black', pch=19)
points(pops$ConZone.dist[pops$Specie == 'Manacus vitellinus'], pops$PC2[pops$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(pops$ConZone.dist, fit.PC2[,1], col='black')
lines(pops$ConZone.dist, fit.PC2[,1]+fit.PC2[,2], col='blue')
lines(pops$ConZone.dist, fit.PC2[,1]-fit.PC2[,2], col='blue')
abline(v=193.4,col="red", lty='dashed')

plot(PC3~ConZone.dist, data=pops, col=grey, xlab='')
points(pops$ConZone.dist[pops$Specie == 'Manacus manacus'], pops$PC3[pops$Specie == 'Manacus manacus'], col='black', pch=19)
points(pops$ConZone.dist[pops$Specie == 'Manacus vitellinus'], pops$PC3[pops$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(pops$ConZone.dist, fit.PC3[,1], col='black')
lines(pops$ConZone.dist, fit.PC3[,1]+fit.PC3[,2], col='blue')
lines(pops$ConZone.dist, fit.PC3[,1]-fit.PC3[,2], col='blue')
abline(v=193.4,col="red", lty='dashed')

par(mar=c(0,0,0,0),oma=c(0,0,0,0),fig=c(0,1,0,1),new=TRUE,cex=1.25)
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("top",legend=c("Manacus manacus","Manacus vitellinus"),pch=19,lty=1,col=c("black","goldenrod1"),bty="n",ncol=2)
dev.off()


#####################################################################################
require(drc)

# Mi optimización
BM.log	<- 	logitm.ls(pars=runif(4,-10,10),x.var=Distance2West,y.var=BodyMass,CI=TRUE,prediction=seq(0,400,length=1000),reps=2000,method="Nelder-Mead")
PC1.log	<-	logitm.ls(pars=PC1.inits,x.var=Distance2West,y.var=morph.pca.x[,1],CI=FALSE,prediction=seq(0,400,length=1000),reps=1000,method="Nelder-Mead")
PC2.log <- 	logitm.ls(pars=runif(4,-10,10),x.var=Distance2West,y.var=morph.pca.x[,2],CI=TRUE,prediction=seq(0,400,length=1000),reps=2000,method="Nelder-Mead")
PC3.log <- 	logitm.ls(pars=runif(4,-10,10),x.var=Distance2West,y.var=morph.pca.x[,3],CI=TRUE,prediction=seq(0,400,length=1000),reps=2000,method="Nelder-Mead")

#BM.log$Fitted
fit.BM	<- logitm.pred(BM.log,x.var=Distance2West,y.var=BodyMass,pred.x=seq(0,400,length=1000))
fit.PC1 <- PC1.log$Fitted
#logitm.pred(PC1.log$par,x.var=Distance2West,y.var=morph.pca.x[,1],pred.x=seq(0,400,length=1000))
fit.PC2 <- PC2.log$Fitted
#logitm.pred(PC2.log$par,x.var=Distance2West,y.var=morph.pca.x[,2],pred.x=seq(0,400,length=1000))
fit.PC3 <- PC3.log$Fitted
#logitm.pred(PC3.log$par,x.var=Distance2West,y.var=morph.pca.x[,3],pred.x=seq(0,400,length=1000))


grey <- rgb(150,150,150, alpha=200, maxColorValue=255)
x.pred	<-	seq(0,400,length=1000)

png('/home/camilo/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v2/Cz_PC.png', units='cm', height = 10, width=25, res=150)
par(mfrow=c(2,2), family='serif')
plot(BodyMass~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], BodyMass[dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], BodyMass[dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, BM.log$Predicted[,1], col='black')
lines(x.pred, BM.log$Predicted[,2], col='blue')
lines(x.pred, BM.log$Predicted[,3], col='blue')
#abline(v=193.4,col="red", lty='dashed')
#abline(v=coef(BM.log)[4],col="red", lty='dashed')

plot(morph.pca.x[,1]~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,1][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,1][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, PC1.log$Predicted[,1], col='black')
lines(x.pred, PC1.log$Predicted[,2], col='blue')
lines(x.pred, PC1.log$Predicted[,3], col='blue')
#abline(v=193.4,col="red", lty='dashed')
#abline(v=coef(PC1.log)[4],col="red", lty='dashed')

plot(morph.pca.x[,2]~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,2][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,2][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, PC2.log$Predicted[,1], col='black')
lines(x.pred, PC2.log$Predicted[,2], col='blue')
lines(x.pred, PC2.log$Predicted[,3], col='blue')
#abline(v=193.4,col="red", lty='dashed')
#abline(v=coef(PC2.log)[4],col="red", lty='dashed')


plot(morph.pca.x[,3]~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,3][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,3][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, PC3.log$Predicted[,1], col='black')
lines(x.pred, PC3.log$Predicted[,2], col='blue')
lines(x.pred, PC3.log$Predicted[,3], col='blue')
#abline(v=193.4,col="red", lty='dashed')
#abline(v=coef(PC2.log)[4],col="red", lty='dashed')
par(mar=c(0,0,0,0),oma=c(0,0,0,0),fig=c(0,1,0,1),new=TRUE,cex=1.25)
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("top",legend=c("Manacus manacus","Manacus vitellinus"),pch=19,lty=1,col=c("black","goldenrod1"),bty="n",ncol=2)
dev.off()

# Con el paquete drc
BM.log		<-	drm(BodyMass~Distance2West,fct=LL.4())
PC1.log		<-	drm(morph.pca.x[,1]~Distance2West,fct=LL.4())

PC1.log		<-	logitm.ls(PC1.inits,x.var=Distance2West,y.var=morph.pca.x[,1],CI=FALSE)
PC2.log		<-	drm(morph.pca.x[,2]~Distance2West,fct=LL.4(),control=drmc(method="Nelder-Mead"),start=rep(1,4))
PC3.log		<-	drm(morph.pca.x[,3]~Distance2West,fct=LL.4(),control=drmc(method="Nelder-Mead"),start=rep(1,4))

fit.PC1	<- drm.boot(PC1.log,newdata=seq(0,400,length=1000),reps=1100)
fit.PC1 <- predict(PC1.log,newdata=data.frame(Distance2West=seq(0,400,length=1000)),interval="confidence")
fit.PC2 <- predict(PC2.log,newdata=data.frame(Distance2West=seq(0,400,length=1000)),se.fit=TRUE)
fit.PC3 <- predict(PC3.log,newdata=data.frame(Distance2West=seq(0,400,length=1000)),se.fit=TRUE)

x.pred	<-	seq(0,400,length=1000)

png('/home/camilo/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v2/Cz_PC.png', units='cm', height = 10, width=25, res=150)
par(mfrow=c(2,2), family='serif')
plot(BodyMass~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], BodyMass[dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], BodyMass[dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, fit.BM[,1], col='black')
lines(x.pred, fit.BM[,1]+1.96*fit.BM[,2], col='blue')
lines(x.pred, fit.BM[,1]-1.96*fit.BM[,2], col='blue')
abline(v=coef(BM.log)[4],col="red", lty='dashed')

plot(morph.pca.x[,1]~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,1][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,1][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, fit.PC1$Prediction[,1], col='black')
lines(x.pred, fit.PC1$Prediction[,2], col='blue')
lines(x.pred, fit.PC1$Prediction[,3], col='blue')
abline(v=coef(PC1.log)[4],col="red", lty='dashed')

plot(morph.pca.x[,2]~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,2][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,2][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, fit.PC2[,1], col='black')
lines(x.pred, fit.PC2[,1]+1.96*fit.PC1[,2], col='blue')
lines(x.pred, fit.PC2[,1]-1.96*fit.PC1[,2], col='blue')
abline(v=coef(PC2.log)[4],col="red", lty='dashed')


plot(morph.pca.x[,3]~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,3][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,3][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, fit.PC3[,1], col='black')
lines(x.pred, fit.PC3[,1]+1.96*fit.PC3[,2], col='blue')
lines(x.pred, fit.PC3[,1]-1.96*fit.PC3[,2], col='blue')
abline(v=coef(PC3.log)[4],col="red", lty='dashed')

par(mar=c(0,0,0,0),oma=c(0,0,0,0),fig=c(0,1,0,1),new=TRUE,cex=1.25)
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("top",legend=c("Manacus manacus","Manacus vitellinus"),pch=19,lty=1,col=c("black","goldenrod1"),bty="n",ncol=2)
dev.off()

wing.log	<-	drm(dat.Zone$Wing~Distance2West,fct=LL.4())
Tail.log	<-	drm(dat.Zone$Tail~Distance2West,fct=LL.4())
bill.log	<-	drm(dat.Zone$Exposed_Culmen~Distance2West,fct=LL.4(),control=drmc(method="Nelder-Mead"))
Tarsus.log	<-	drm(dat.Zone$Tarsus~Distance2West,fct=LL.4())


fit.Tars	<-	predict(Tarsus.log,newdata=data.frame(Distance2West=seq(0,400,length=1000)),se.fit=TRUE)
fit.bill	<-	predict(bill.log,newdata=data.frame(Distance2West=seq(0,400,length=1000)),se.fit=TRUE)
fit.tail	<-	predict(Tail.log,newdata=data.frame(Distance2West=seq(0,400,length=1000)),se.fit=TRUE)


plot(dat.Zone$Tarsus~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], dat.Zone$Tarsus[dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], dat.Zone$Tarsus[dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, fit.Tars[,1], col='black')
lines(x.pred, fit.Tars[,1]+1.96*fit.Tars[,2], col='blue')
lines(x.pred, fit.Tars[,1]-1.96*fit.Tars[,2], col='blue')
abline(v=coef(Tarsus.log)[4],col="red", lty='dashed')

plot(dat.Zone$Exposed_Culmen~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], dat.Zone$Exposed_Culmen[dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], dat.Zone$Exposed_Culmen[dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, fit.bill[,1], col='black')
lines(x.pred, fit.bill[,1]+1.96*fit.bill[,2], col='blue')
lines(x.pred, fit.bill[,1]-1.96*fit.bill[,2], col='blue')
abline(v=coef(bill.log)[4],col="red", lty='dashed')

plot(dat.Zone$Tail~Distance2West, xlab='',type="n")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], dat.Zone$Tail[dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], dat.Zone$Tail[dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(x.pred, fit.tail[,1], col='black')
lines(x.pred, fit.tail[,1]+1.96*fit.tail[,2], col='blue')
lines(x.pred, fit.tail[,1]-1.96*fit.tail[,2], col='blue')
abline(v=coef(Tail.log)[4],col="red", lty='dashed')
