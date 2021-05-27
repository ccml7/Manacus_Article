# Funciones para hacer bootstrap parametrico
# Funcion logistica
logitm <- function(b,c,d,e,x){(c+(d-c)/(1+exp(b*(log(x)-log(e)))))}

# Bootstrap parametrico
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
# Paquetes
# Exploracion de la zona de hibridacion
library(geosphere)
library(drc)

# Preparacion de los datos
dat				<-	read.csv("~/Dropbox/Tesis Camilo/Data/Data_Manacus.csv")
morph.pca		<-	prcomp(~Wing+Tail+Exposed_Culmen+Tarsus,scale.=TRUE,center=TRUE,data=dat)
unique.locs		<-	unique(dat$Population)
unique.locs.pos	<-	match(unique.locs,dat$Population)

all.locs	 	<-	dat[unique.locs.pos,c("Population","Long","Lat")]

all.locs.dist	<-	distm(all.locs[,c("Long","Lat")])/1000
PV.pos			<-	which(all.locs$Population=="Puerto Valdivia")

locs.include.pos		<-	which(all.locs.dist[PV.pos,]<=300)

pops		<-	all.locs[locs.include.pos,]

locs.include	<-	pops$Population

ConZone.pop1<-	which(pops$Population=="Acandi")

ConZone.dist<-	distm(pops[,c("Long","Lat")])
ConZone.dist<-	ConZone.dist[ConZone.pop1,]/1000
names(ConZone.dist)	<- pops$Population

ConZone.inds	<-	which(!is.na(match(dat$Population,locs.include)))
morph.pca.x		<-	morph.pca$x[ConZone.inds,]
BodyMass		<-	dat$Body_mass[ConZone.inds]
dat.Zone		<-	dat[ConZone.inds,]

Distance2West	<-	ConZone.dist[match(dat.Zone$Population,names(ConZone.dist))]


# Con el paquete drc
BM.log		<-	drm(BodyMass~Distance2West,fct=LL.4())
BM.lm		<-	lm(BodyMass~1)
BM.lm1		<-	lm(BodyMass~Distance2West)
PC1.log		<-	drm(morph.pca.x[,1]~Distance2West,fct=LL.4())
PC1.lm		<-	lm(morph.pca.x[,1]~1)
PC1.lm1		<-	lm(morph.pca.x[,1]~Distance2West)
PC2.log		<-	drm(morph.pca.x[,2]~Distance2West,fct=LL.4())
PC2.lm		<-	lm(morph.pca.x[,2]~1)
PC2.lm1		<-	lm(morph.pca.x[,2]~Distance2West)
PC3.log		<-	drm(morph.pca.x[,3]~Distance2West,fct=LL.4())
PC3.lm		<-	lm(morph.pca.x[,3]~1)
PC3.lm1		<-	lm(morph.pca.x[,3]~Distance2West)

BM.BIC.vec	<-	c(BIC(BM.log),BIC(BM.lm),BIC(BM.lm1))
BM.best		<-	which.min(BM.BIC.vec)
BM.delta.bic<-	BM.BIC.vec-BM.BIC.vec[BM.best]
PC1.BIC.vec	<-	c(BIC(PC1.log),BIC(PC1.lm),BIC(PC1.lm1))
PC1.best	<-	which.min(PC1.BIC.vec)
PC1.delta.bic<-	PC1.BIC.vec - PC1.BIC.vec[PC1.best]
PC2.BIC.vec	<-	c(BIC(PC2.log),BIC(PC2.lm),BIC(PC2.lm1))
PC2.best	<-	which.min(PC2.BIC.vec)
PC2.delta.bic<-	PC2.BIC.vec - PC2.BIC.vec[PC2.best]
PC3.BIC.vec	<-	c(BIC(PC3.log),BIC(PC3.lm),BIC(PC3.lm1))
PC3.best	<-	which.min(PC3.BIC.vec)
PC3.delta.bic<-	PC3.BIC.vec - PC3.BIC.vec[PC3.best]

fit.BM	<-	predict(BM.lm,newdata=data.frame(Distance2West=seq(0,500,length=1000)),se.fit=TRUE)
fit.PC1	<- 	drm.boot(PC1.log,newdata=seq(0,500,length=1000),reps=1100)
fit.PC2	<-	predict(PC2.lm1,newdata=data.frame(Distance2West=seq(0,500,length=1000)),se.fit=TRUE)
fit.PC3	<- 	predict(PC3.lm,newdata=data.frame(Distance2West=seq(0,500,length=1000)),se.fit=TRUE)

x.pred	<-	seq(0,500,length=1000)

newDist2west	<-seq(0,500,length=1000)	

tiff('~/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v3/Figure3.tiff', units='in', height = 6, width=6, res=300
	,compression="lzw",type="cairo",family="times")
par(mfrow=c(2,2),mar=c(2,3,1,1),oma=c(2,1,2.5,1),mgp=c(2,0.5,0),tcl=-0.25)
plot(BodyMass~Distance2West, xlab='',type="n",ylab="Body Mass (g)")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], BodyMass[dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], BodyMass[dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(newDist2west, fit.BM$fit, col='red')
segments(x0=0,y0=fit.BM$fit+1.96*fit.BM$se.fit,x1=500,col="red",lty=2)
segments(x0=0,y0=fit.BM$fit-1.96*fit.BM$se.fit,x1=500,col="red",lty=2)
mtext("A.)",adj=-0.15)

plot(morph.pca.x[,1]~Distance2West, xlab='',type="n",ylab="Morphological PC1")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,1][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,1][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(newDist2west, fit.PC1$Prediction[,1], col='red')
lines(newDist2west, fit.PC1$Prediction[,2], col='red',lty=2)
lines(newDist2west, fit.PC1$Prediction[,3], col='red',lty=2)
mtext("B.)",adj=-0.15)

plot(morph.pca.x[,2]~Distance2West, xlab='',type="n",ylab="Morphological PC2")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,2][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,2][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(newDist2west, fit.PC2$fit, col='red')
lines(newDist2west, fit.PC2$fit + 1.96*fit.PC2$se.fit, col='red',lty=2)
lines(newDist2west, fit.PC2$fit - 1.96*fit.PC2$se.fit, col='red',lty=2)
mtext("C.)",adj=-0.15)

plot(morph.pca.x[,3]~Distance2West, xlab='',type="n",ylab="Morphological PC3")
points(Distance2West[dat.Zone$Specie == 'Manacus manacus'], morph.pca.x[,3][dat.Zone$Specie == 'Manacus manacus'], col='black', pch=19)
points(Distance2West[dat.Zone$Specie == 'Manacus vitellinus'], morph.pca.x[,3][dat.Zone$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)
lines(newDist2west, fit.PC3$fit, col='red')
segments(x0=0,y0=fit.PC3$fit+1.96*fit.PC3$se.fit,x1=500,col="red",lty=2)
segments(x0=0,y0=fit.PC3$fit-1.96*fit.PC3$se.fit,x1=500,col="red",lty=2)
mtext("D.)",adj=-0.15)
mtext("Distance from Westernmost population (km)",side=1,outer=TRUE)

par(mar=c(0,0,0,0),oma=c(0,0,0,0),fig=c(0,1,0,1),new=TRUE,cex=1.25)
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("top",legend=c("white phenotype","yellow phenotype"),pch=19,lty=1,col=c("black","goldenrod1"),bty="n",ncol=2)
dev.off()

