library(lme4)
library(ape)
library(geosphere)
library(nlme)
library(mgcv)
library(MuMIn)
library(raster)

dat	<-	read.csv("~/Dropbox/Tesis Camilo/Data/Data_Manacus.csv")
geo.dist	<-	read.delim("~/Dropbox/Tesis Camilo/Data/LCP_Matrix.txt")
geo.dist	<-	as.dist(geo.dist)
geo.dist.sc	<-	scale(geo.dist)
geo.pcoa	<-	cmdscale(as.dist(geo.dist),k=1)
sc.pcoa		<-	scale(geo.pcoa)

clim.vars.path	<-	list.files("/Volumes/JUANP/GIS/bioclim",full.names=TRUE,pattern=".bil")
elev.path		<-	"/Volumes/JUANP/GIS/Elevation_1km/elevation_1KMmd_SRTM.tif"
clim.vars		<-	stack(clim.vars.path)
elev			<-	raster(elev.path)

clim.dat		<-	extract(clim.vars,dat[,c("Long","Lat")])
elev.dat		<-	extract(elev,dat[,c("Long","Lat")])

clim.dat		<-	cbind(clim.dat,elev.dat)

clim.pca		<-	prcomp(~.,data=data.frame(clim.dat),scale.=TRUE,center=TRUE)
morph.pca		<-	prcomp(~Wing+Tail+Exposed_Culmen+Tarsus,scale.=TRUE,center=TRUE,data=dat)

new.dat		<-	data.frame(Wing=dat$Wing, Tail=dat$Tail, Culmen=dat$Exposed_Culmen
							,PC1=dat$PC1, Mass=dat$Body_mass
							,bio1=scale(dat$bio1), bio3=scale(dat$bio3)
							,bio7=scale(dat$bio7), bio12=scale(dat$bio12)
							,bio15=scale(dat$bio15), bio18=scale(dat$bio18)
							,Altitude=scale(dat$Altitude), geog=sc.pcoa[,1]
							,Sex=factor(dat$Sex,levels=c("Female","Male")),Source=dat$Source
							,Subspecies=dat$Specie, Population=dat$Population)

new.dat		<-	data.frame(Mass=dat$Body_mass,mPC1=morph.pca$x[,1],mPC2=morph.pca$x[,2]
							,mPC3=morph.pca$x[,3],cPC1=clim.pca$x[,1],cPC2=clim.pca$x[,2]
							,cPC3=clim.pca$x[,3],cPC4=clim.pca$x[,4],geog=sc.pcoa[,1]
							,Sex=dat$Sex,Source=dat$Source
							,Subspecies=dat$Specie)
							
full.dat		<-	data.frame(Subspecies=dat$Specie,Sex=dat$Sex,Source=dat$Source,VoucherID=dat$Voucher_Code
							,Longitude=dat$Long,Latitude=dat$Lat
							,Mass=dat$Body_mass,Wing=dat$Wing,Tail=dat$Tail,Culmen=dat$Exposed_Culmen
							,mPC1=morph.pca$x[,1],mPC2=morph.pca$x[,2],mPC3=morph.pca$x[,3]
							,bio1=clim.dat[,"bio_1"], bio3=clim.dat[,"bio_3"],bio7=clim.dat[,"bio_7"], bio12=clim.dat[,"bio_12"]
							,bio15=clim.dat[,"bio_15"], bio18=clim.dat[,"bio_18"],Altitude=clim.dat[,"elev.dat"]
							,cPC1=clim.pca$x[,1],cPC2=clim.pca$x[,2],cPC3=clim.pca$x[,3]
							,cPC4=clim.pca$x[,4])
full.dat$Subspecies	<-	as.character(full.dat$Subspecies)
full.dat$Subspecies[full.dat$Subspecies=="Manacus manacus"]	<-	"White"
full.dat$Subspecies[full.dat$Subspecies=="Manacus vitellinus"]	<-	"Yellow"
full.dat$Sex	<-	as.character(full.dat$Sex)
full.dat$Sex[full.dat$Sex==""]	<-	"Unknown"

write.table(full.dat,"~/Dropbox/Tesis Camilo/Articulo_Manacus/Supplementary_data_museum_nets.txt",sep="\t",row.names=FALSE
			,quote=FALSE)

mass.model.geo	<-	lmer(Mass~geog+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit)	
summary(mass.model.geo)

mass.model.env <- lmer(Mass~cPC1+cPC2+cPC3+cPC4+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit)
summary(mass.model.env)

mass.model.add <- lmer(Mass~geog+cPC1+cPC2+cPC3+cPC4+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit)
summary(mass.model.add)

mass.model.null		<- lm(Mass~1,dat=new.dat)
mass.model.null.re	<- lmer(Mass~1+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,control=lmerControl(optimizer="Nelder_Mead"))
	

anova(mass.model.null.re,mass.model.geo,mass.model.env,mass.model.add)

mass.model.env.Sp		<-	lmer(Mass~cPC1+cPC2+cPC3+cPC4+(1|Subspecies),data=new.dat,na.action=na.omit)
mass.model.env.Source	<-	lmer(Mass~cPC1+cPC2+cPC3+cPC4+(1|Source),data=new.dat,na.action=na.omit)
mass.model.env.Sex		<-	lmer(Mass~cPC1+cPC2+cPC3+cPC4+(1|Sex),data=new.dat,na.action=na.omit)

anova(mass.model.env,mass.model.env.Sp,mass.model.env.Source,mass.model.env.Sex)
#mass.fit	<-	predict(mass.model.env,newdata=pred.dat.mass,na.action=na.omit,level=0, se=T)


############# PC1 ################
pc1.model.geo	 <- 	lmer(mPC1~geog+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat, na.action=na.omit,REML=FALSE,control=lmerControl(optimizer="Nelder_Mead"))	
summary(pc1.model.geo)

pc1.model.env	 <- 	lmer(mPC1~cPC1+cPC2+cPC3+cPC4+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit,REML=FALSE)	
summary(pc1.model.env)

pc1.model.add	 <- 	lmer(mPC1~cPC1+cPC2+cPC3+cPC4+geog+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))	
summary(pc1.model.add)

pc1.model.null	<-	lmer(mPC1~(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))

anova(pc1.model.null,pc1.model.geo,pc1.model.env,pc1.model.add)

pc1.model.env.Sp	 <- 	lmer(mPC1~cPC1+cPC2+cPC3+cPC4+(1|Subspecies),data=new.dat,na.action=na.omit)	
pc1.model.env.Source <-	lmer(mPC1~cPC1+cPC2+cPC3+cPC4+(1|Source),data=new.dat,na.action=na.omit)	
pc1.model.env.Sex	<-	lmer(mPC1~cPC1+cPC2+cPC3+cPC4+(1|Sex),data=new.dat,na.action=na.omit)	

pc1.model.env.nor	<-	lm(mPC1~cPC1+cPC2+cPC3+cPC4,data=new.dat,na.action=na.omit)

anova(pc1.model.env,pc1.model.env.Sp,pc1.model.env.Source,pc1.model.env.Sex)

#range(new.dat$bio15)

#pred.dat.pc1	<-	data.frame(bio1=rep(0,200), bio7=rep(0,200), bio12=rep(0,200)
							,bio15=seq(-3.4,4.2, length=200),  bio18=rep(0,200), Sex=rep("Male",200)
							,Subspecies=rep('Manacus manacus',200))

#pc1.fit	<-	predict(pc1.model.env,newdata=pred.dat.pc1,na.action=na.omit,level=0, se=T)

#range(new.dat$pc3)



############ PC2 #################
pc2.model.geo	 <- 	lmer(mPC2~geog+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat, na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))	
summary(pc2.model.geo)

pc2.model.env	 <- 	lmer(mPC2~cPC1+cPC2+cPC3+cPC4+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit)	
summary(pc2.model.env)
anova(pc2.model.env)

pc2.model.add	 <- 	lmer(mPC2~cPC1+cPC2+cPC3+cPC4+geog+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))	
summary(pc2.model.add)

pc2.model.null	<-	lmer(mPC2~(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat,na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))

anova(pc2.model.null,pc2.model.geo,pc2.model.env,pc2.model.add)

pc2.model.geo.Sp	<-	lmer(mPC2~(1|Subspecies),data=new.dat
								,na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))
pc2.model.geo.Source	<-	lmer(mPC2~(1|Source),data=new.dat
								,na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))
pc2.model.geo.Sex	<-	lmer(mPC2~(1|Sex),data=new.dat
								,na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))

anova(pc2.model.null,pc2.model.geo.Sp,pc2.model.geo.Source,pc2.model.geo.Sex)
							

#range(new.dat$bio7)
#pred.dat.pc2	<-	data.frame(bio1=rep(0,200), bio7=rep(seq(-4.5, 2.7, length=100),2), bio12=rep(0,200)
							,bio15=rep(0,200),  bio18=rep(0,200), Sex=rep("Male",200)
							,Subspecies=factor(c(rep('Manacus manacus',100), rep('Manacus vitellinus',100))))

#pc2.fit	<-	predict(pc2.model.env,newdata=pred.dat.pc2,na.action=na.omit,level=0, se=T)



############ PC3 ################
pc3.model.geo	 <- 	lmer(mPC3~geog+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat, na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))
summary(pc3.model.geo)

pc3.model.env	 <- 	lmer(mPC3~cPC1+cPC2+cPC3+cPC4+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat, na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))
summary(pc3.model.env)
anova(pc3.model.env)


pc3.model.add	 <- 	lmer(mPC3~geog+cPC1+cPC2+cPC3+cPC4+(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat, na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))
summary(pc3.model.add)

pc3.model.null	<-	lmer(mPC3~(1|Subspecies)+(1|Sex)+(1|Source),data=new.dat, na.action=na.omit,control=lmerControl(optimizer="Nelder_Mead"))
pc3.model.null1	<-	lm(mPC3~1,data=new.dat)

anova(pc3.model.null,pc3.model.geo,pc3.model.env,pc3.model.add)

pc3.model.env.Sp	<-	lmer(mPC3~cPC1+cPC2+cPC3+cPC4+(1|Subspecies),data=new.dat, na.action=na.omit)
pc3.model.env.Source<-	lmer(mPC3~cPC1+cPC2+cPC3+cPC4+(1|Source),data=new.dat, na.action=na.omit)
pc3.model.env.Sex	<-	lmer(mPC3~cPC1+cPC2+cPC3+cPC4+(1|Sex),data=new.dat, na.action=na.omit)
pc3.model.env.Sex.Source	<-	lmer(mPC3~cPC1+cPC2+cPC3+cPC4+(1|Sex)+(1|Source),data=new.dat, na.action=na.omit)


anova(pc3.model.env,pc3.model.env.Sp,pc3.model.env.Source,pc3.model.env.Sex)


#pred.dat.pc3	<-	data.frame(bio1=rep(0,200), bio3=rep(0,200), bio12=seq(-1.5,4.6,length=200)
							,bio15=rep(0,200),  bio18=rep(0,200), Sex=rep("Male",200)
							,Subspecies=rep('Manacus manacus', 200), Source=rep('Museum',200))

#pc3.fit	<-	predict(pc3.model.env,newdata=pred.dat.pc3,na.action=na.omit,level=0, se=T)


# Figure 
# Body mass
bod.mass.dat.PC1	<-	data.frame(cPC1=rep(seq(min(new.dat$cPC1),max(new.dat$cPC1),length=500),3),cPC2=rep(0,1500)
									,cPC3=rep(0,500),cPC4=rep(0,500),Sex=c(rep("Male",500),rep("Female",500),rep("Unknown",500)))
bod.mass.PC1	<-	predict(mass.model.env.Sex,newdata=bod.mass.dat.PC1,se.fit=TRUE,re.form=NA)
bod.mass.PC1.re		<-	predict(mass.model.env.Sex,newdata=bod.mass.dat.PC1,re.form=~(1|Sex))				

bod.mass.dat.PC4	<-	data.frame(cPC1=rep(0,1500),cPC2=rep(0,1500)
									,cPC3=rep(0,500),cPC4=rep(seq(min(new.dat$cPC4),max(new.dat$cPC4),length=500),3)
									,Sex=c(rep("Male",500),rep("Female",500),rep("Unknown",500)))
bod.mass.PC4		<-	predict(mass.model.env.Sex,newdata=bod.mass.dat.PC4,re.form=NA,se.fit=TRUE)
bod.mass.PC4.re		<-	predict(mass.model.env.Sex,newdata=bod.mass.dat.PC4,re.form=~(1|Sex))

	
# PC1
PC1.newdat	<-	data.frame(cPC1=rep(0,1500),cPC2=rep(0,1500),cPC3=rep(0,1500)
							,cPC4=rep(seq(min(new.dat$cPC4),max(new.dat$cPC4),length=500),3)
							,Sex=c(rep("Male",500),rep("Female",500),rep("Unknown",500))
							,Source=c(rep("Museum",500),rep("Net",500),rep(NA,500))
							,Subspecies=c(rep("Manacus manacus",500),rep("Manacus vitellinus",500),rep(NA,500)))
							
PC1.pred	<-	predict(pc1.model.env,newdata=PC1.newdat,se.fit=TRUE,re.form=NA)
PC1.pred.Sex	<-	predict(pc1.model.env,newdata=PC1.newdat,re.form=~(1|Sex))
PC1.pred.Source	<-	predict(pc1.model.env,newdata=PC1.newdat,re.form=~(1|Source),allow.new.levels=TRUE)
PC1.pred.Sp	<-	predict(pc1.model.env,newdata=PC1.newdat,re.form=~(1|Subspecies),allow.new.levels=TRUE)

#PC2
PC2.newdat	<-	data.frame(geog=rep(seq(min(new.dat$geog),max(new.dat$geog),length=500),3)
							,Sex=c(rep("Male",500),rep("Female",500),rep("Unknown",500))
							,Source=c(rep("Museum",500),rep("Net",500),rep(NA,500))
							,Subspecies=c(rep("Manacus manacus",500),rep("Manacus vitellinus",500),rep(NA,500)))
PC2.pred		<-	predict(pc2.model.geo,newdata=PC2.newdat,se.fit=TRUE,re.form=NA)
PC2.pred.Sex	<-	predict(pc2.model.geo,newdata=PC2.newdat,re.form=~(1|Sex))
PC2.pred.Source	<-	predict(pc2.model.geo,newdata=PC2.newdat,re.form=~(1|Source),allow.new.levels=TRUE)
PC2.pred.Sp		<-	predict(pc2.model.geo,newdata=PC2.newdat,re.form=~(1|Subspecies),allow.new.levels=TRUE)
	
#PC3
PC3.newdat	<-	data.frame(cPC1=rep(0,1500),cPC2=rep(0,1500)
							,cPC3=rep(0,1500),cPC4=rep(seq(min(new.dat$cPC4),max(new.dat$cPC4),length=500),3)
							,Sex=c(rep("Male",500),rep("Female",500),rep("Unknown",500))
							,Source=c(rep("Museum",500),rep("Net",500),rep(NA,500))
							,Subspecies=c(rep("Manacus manacus",500),rep("Manacus vitellinus",500),rep(NA,500)))
PC3.pred		<-	predict(pc3.model.env,newdata=PC3.newdat,se.fit=TRUE,re.form=NA)
PC3.pred.Source	<-	predict(pc3.model.env,newdata=PC3.newdat,re.form=~(1|Source),allow.new.levels=TRUE)
PC3.pred.Sex	<-	predict(pc3.model.env,newdata=PC3.newdat,re.form=~(1|Sex))
PC3.pred.Sp	<-	predict(pc3.model.env,newdata=PC3.newdat,re.form=~(1|Subspecies),allow.new.levels=TRUE)

# Figure 2
tiff("~/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v3/Figure2.tiff",width=7,height=7.7,units="in"
	,res=300,compression="lzw",type="cairo",family="serif")
mat	<-	matrix(c(1,1,2,2,0,0,0,0,3,3,4,4,0,0,0,0,0,5,5,0),ncol=4,byrow=TRUE)
layout(mat,heights=c(1,0.05,1,0.05,1))	
par(mar=c(3,3.5,1.5,1.5),oma=c(1,1,3,1),mgp=c(1.75,0.5,0),tcl=-0.25,cex.axis=1.5,cex.lab=1.5)						

# Body mass
plot(new.dat$cPC1,new.dat$Mass,pch=19,xlab="Climatic PC1",ylab="Body Mass (g)",type="n")
points(new.dat$cPC1[new.dat$Subspecies=="Manacus manacus"],new.dat$Mass[new.dat$Subspecies=="Manacus manacus"],pch=19)
points(new.dat$cPC1[new.dat$Subspecies=="Manacus vitellinus"],new.dat$Mass[new.dat$Subspecies=="Manacus vitellinus"],pch=19,col="goldenrod1")
	points(bod.mass.dat.PC1$cPC1[1:500],bod.mass.PC1$fit[1:500],pch=19, type="l",col="red")
	lines(bod.mass.dat.PC1$cPC1[1:500], bod.mass.PC1$fit[1:500] + bod.mass.PC1$se.fit[1:500], lty=2, lwd=.7, col='red')
	lines(bod.mass.dat.PC1$cPC1[1:500], bod.mass.PC1$fit[1:500] - bod.mass.PC1$se.fit[1:500], lty=2, lwd=.7, col='red')
mtext("A.",adj=-0.2,cex=1.5)


plot(new.dat$cPC4,new.dat$Mass,pch=19,xlab="Climatic PC4",ylab="Body Mass (g)",type="n")
points(new.dat$cPC4[new.dat$Subspecies=="Manacus manacus"],new.dat$Mass[new.dat$Subspecies=="Manacus manacus"],pch=19)
points(new.dat$cPC4[new.dat$Subspecies=="Manacus vitellinus"],new.dat$Mass[new.dat$Subspecies=="Manacus vitellinus"],pch=19,col="goldenrod1")
	points(bod.mass.dat.PC4$cPC4[1:500],bod.mass.PC4$fit[1:500],pch=19, type="l",col="red")
	lines(bod.mass.dat.PC4$cPC4[1:500], bod.mass.PC4$fit[1:500] + bod.mass.PC4$se.fit[1:500], lty=2, lwd=.7, col='red')
	lines(bod.mass.dat.PC4$cPC4[1:500], bod.mass.PC4$fit[1:500] - bod.mass.PC4$se.fit[1:500], lty=2, lwd=.7, col='red')
mtext("B.",adj=-0.2,cex=1.5)	
# PC1

plot(new.dat$cPC4,new.dat$mPC1,pch=19,xlab="Climatic PC4 ",ylab="Morphological PC1",type="n")
points(new.dat$cPC4[new.dat$Subspecies=="Manacus manacus"],new.dat$mPC1[new.dat$Subspecies=="Manacus manacus"],pch=19)
points(new.dat$cPC4[new.dat$Subspecies=="Manacus vitellinus"],new.dat$mPC1[new.dat$Subspecies=="Manacus vitellinus"],pch=19,col="goldenrod1")
	points(PC1.newdat$cPC4[1:500],PC1.pred$fit[1:500],pch=19, type="l",col="red")
	lines(PC1.newdat$cPC4[1:500], PC1.pred$fit[1:500] + PC1.pred$se.fit[1:500], lty=2, lwd=.7, col='red')
	lines(PC1.newdat$cPC4[1:500], PC1.pred$fit[1:500] - PC1.pred$se.fit[1:500], lty=2, lwd=.7, col='red')
mtext("C.",adj=-0.2,cex=1.5)	
# PC3

plot(new.dat$cPC4,new.dat$mPC3,pch=19,xlab="Climatic PC4",ylab="Morphological PC3")
points(new.dat$cPC4[new.dat$Subspecies=="Manacus manacus"],new.dat$mPC3[new.dat$Subspecies=="Manacus manacus"],pch=19)
points(new.dat$cPC4[new.dat$Subspecies=="Manacus vitellinus"],new.dat$mPC3[new.dat$Subspecies=="Manacus vitellinus"],pch=19,col="goldenrod1")
	points(PC3.newdat$cPC4[1:500],PC3.pred$fit[1:500],pch=19, type="l",col="red")
	lines(PC3.newdat$cPC4[1:500], PC3.pred$fit[1:500] + PC3.pred$se.fit[1:500], lty=2, lwd=.7, col='red')
	lines(PC3.newdat$cPC4[1:500], PC3.pred$fit[1:500] - PC3.pred$se.fit[1:500], lty=2, lwd=.7, col='red')
mtext("D.",adj=-0.2,cex=1.5)

# PC2

plot(new.dat$geog,new.dat$mPC2,pch=19,xlab="Geographical PCoA (LCP)",ylab="Morphological PC2",type="n")
points(new.dat$geog[new.dat$Subspecies=="Manacus manacus"],new.dat$mPC2[new.dat$Subspecies=="Manacus manacus"],pch=19)
points(new.dat$geog[new.dat$Subspecies=="Manacus vitellinus"],new.dat$mPC2[new.dat$Subspecies=="Manacus vitellinus"],pch=19,col="goldenrod1")
#	points(PC2.newdat$geog[1:500],PC2.pred$fit[1:500],pch=19, type="l",col="red")
#	lines(PC2.newdat$geog[1:500], PC2.pred$fit[1:500] + PC2.pred$se.fit[1:500], lty=2, lwd=.7, col='red')
#	lines(PC2.newdat$geog[1:500], PC2.pred$fit[1:500] - PC2.pred$se.fit[1:500], lty=2, lwd=.7, col='red')
mtext("E.",adj=-0.2,cex=1.5)

par(mar=c(0,0,0,0),oma=c(0,0,0,0),fig=c(0,1,0,1),new=TRUE,cex=1.25)
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("top",legend=c("white phenotype","yellow phenotype"),pch=19,lty=1,col=c("black","goldenrod1"),bty="n",ncol=2)

dev.off()


# Figure S2
tiff("~/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v3/FigureS2.tiff",width=7.5,height=7.7,units="in"
	,res=300,compression="lzw",type="cairo",family="serif")
mat	<-	matrix(c(12,12,12,12,12,12
				,1,1,1,2,2,2
				,3,3,4,4,5,5
				,6,6,7,7,8,8
				,9,9,10,10,11,11),ncol=6,byrow=TRUE)
layout(mat,heights=c(0.25,1,1,1,1),widths=c(1,1,1,1,1,1))	
par(mar=c(3,2,1.5,1.5),oma=c(1,3,1,1),mgp=c(1.75,0.5,0),tcl=-0.25,cex.axis=1.1,cex.lab=1.1)						

# Body mass
plot(new.dat$cPC2,new.dat$Mass,pch=19,xlab="Climatic PC2",ylab="",type="n")
points(new.dat$cPC2[new.dat$Sex=="Male"],new.dat$Mass[new.dat$Sex=="Male"])
points(new.dat$cPC2[new.dat$Sex=="Female"],new.dat$Mass[new.dat$Sex=="Female"],pch=19,col="darkgreen")
points(new.dat$cPC2[new.dat$Sex=="Unknown"],new.dat$Mass[new.dat$Sex=="Unknown"],pch=19,col="grey")
	points(bod.mass.dat.PC2$cPC2[1:500],bod.mass.PC2$fit[1:500],pch=19, type="l",col="black",lty=2,lwd=3)
	points(bod.mass.dat.PC2$cPC2[1:500],bod.mass.PC2.re[1:500],pch=19, type="l",col="black")
	lines(bod.mass.dat.PC2$cPC2[501:1000], bod.mass.PC2.re[501:1000], col='darkgreen')
	lines(bod.mass.dat.PC2$cPC2[1001:1500], bod.mass.PC2.re[1001:1500], col='grey')
mtext("A.",adj=-0.1,cex=1.25)


plot(new.dat$cPC4,new.dat$Mass,pch=19,xlab="Climatic PC4",ylab="",type="n")
points(new.dat$cPC4[new.dat$Sex=="Male"],new.dat$Mass[new.dat$Sex=="Male"])
points(new.dat$cPC4[new.dat$Sex=="Female"],new.dat$Mass[new.dat$Sex=="Female"],pch=19,col="darkgreen")
points(new.dat$cPC4[new.dat$Sex=="Unknown"],new.dat$Mass[new.dat$Sex=="Unknown"],pch=19,col="grey")
	points(bod.mass.dat.PC4$cPC4[1:500],bod.mass.PC4$fit[1:500],pch=19, type="l",col="black",lwd=3,lty=2)
	points(bod.mass.dat.PC4$cPC4[1:500],bod.mass.PC4.re[1:500],pch=19, type="l",col="black")
	lines(bod.mass.dat.PC4$cPC4[501:1000], bod.mass.PC4.re[501:1000], col='darkgreen')
	lines(bod.mass.dat.PC4$cPC4[1001:1500], bod.mass.PC4.re[1001:1500], col='grey')

mtext("B.",adj=-0.1,cex=1.25)

# PC1
plot(new.dat$cPC4,new.dat$mPC1,pch=19,xlab="Climatic PC4",ylab="",type="n")
points(new.dat$cPC4[new.dat$Subspecies=="Manacus manacus"],new.dat$mPC1[new.dat$Subspecies=="Manacus manacus"],pch=19)
points(new.dat$cPC4[new.dat$Subspecies=="Manacus vitellinus"],new.dat$mPC1[new.dat$Subspecies=="Manacus vitellinus"],pch=19,col="goldenrod1")
	points(PC1.newdat$cPC4[1:500],PC1.pred$fit[1:500],pch=19, type="l",col="black",lty=2,lwd=3)
	lines(PC1.newdat$cPC4[1:500], PC1.pred.Sp[1:500], col='black')
	lines(PC1.newdat$cPC4[501:1000], PC1.pred.Sp[501:1000], col='goldenrod1')
mtext("C.",adj=-0.2,cex=1.25)	

plot(new.dat$cPC4,new.dat$mPC1,pch=19,xlab="Climatic PC4",ylab="",type="n")
points(new.dat$cPC4[new.dat$Sex=="Male"],new.dat$mPC1[new.dat$Sex=="Male"])
points(new.dat$cPC4[new.dat$Sex=="Female"],new.dat$mPC1[new.dat$Sex=="Female"],pch=19,col="darkgreen")
points(new.dat$cPC4[new.dat$Sex=="Unknown"],new.dat$mPC1[new.dat$Sex=="Unknown"],pch=19,col="grey")
points(PC1.newdat$cPC4[1:500],PC1.pred$fit[1:500],pch=19, type="l",col="black",lty=2,lwd=3)
lines(PC1.newdat$cPC4[1:500], PC1.pred.Sex[1:500], col='black')
lines(PC1.newdat$cPC4[501:1000], PC1.pred.Sex[501:1000], col='darkgreen')
lines(PC1.newdat$cPC4[1001:1500], PC1.pred.Sex[1001:1500], col='grey')

mtext("D.",adj=-0.2,cex=1.25)

plot(new.dat$cPC4,new.dat$mPC1,pch=19,xlab="Climatic PC4",ylab="",type="n")
points(new.dat$cPC4[new.dat$Source=="Museum"],new.dat$mPC1[new.dat$Source=="Museum"],pch=19,col='#1f78b4')
points(new.dat$cPC4[new.dat$Source=="Net"],new.dat$mPC1[new.dat$Source=="Net"],pch=19,col='#b2df8a')
points(PC1.newdat$cPC4[1:500],PC1.pred$fit[1:500],pch=19, type="l",col="black",lty=2,lwd=3)
lines(PC1.newdat$cPC4[1:500], PC1.pred.Source[1:500],col='#1f78b4')
lines(PC1.newdat$cPC4[501:1000], PC1.pred.Source[501:1000],col='#b2df8a')

mtext("E.",adj=-0.2,cex=1.25)

# PC2
plot(new.dat$geog,new.dat$mPC2,pch=19,xlab="Geographical PCoA (LCP)",ylab="",type="n")
points(new.dat$geog[new.dat$Subspecies=="Manacus manacus"],new.dat$mPC2[new.dat$Subspecies=="Manacus manacus"],pch=19)
points(new.dat$geog[new.dat$Subspecies=="Manacus vitellinus"],new.dat$mPC2[new.dat$Subspecies=="Manacus vitellinus"],pch=19,col="goldenrod1")
	points(PC2.newdat$geog[1:500],rep(fixef(pc2.model.null),500), type="l",col="black",lty=2,lwd=3)
	lines(PC2.newdat$geog[1:500], rep(ranef(pc2.model.null)$Subspecies[1,1],500), col='black')
	lines(PC2.newdat$geog[501:1000], rep(ranef(pc2.model.null)$Subspecies[2,1],500), col='goldenrod1')
mtext("F.",adj=-0.2,cex=1.25)	

plot(new.dat$geog,new.dat$mPC2,pch=19,xlab="Geographical PCoA (LCP)",ylab="",type="n")
points(new.dat$geog[new.dat$Sex=="Male"],new.dat$mPC2[new.dat$Sex=="Male"])
points(new.dat$geog[new.dat$Sex=="Female"],new.dat$mPC2[new.dat$Sex=="Female"],pch=19,col="darkgreen")
points(new.dat$geog[new.dat$Sex=="Unknown"],new.dat$mPC2[new.dat$Sex=="Unknown"],pch=19,col="grey")
points(PC2.newdat$geog[1:500],rep(fixef(pc2.model.null),500),pch=19, type="l",col="black",lty=2,lwd=3)
lines(PC2.newdat$geog[1:500], rep(ranef(pc2.model.null)$Sex[1,1],500), col='black')
lines(PC2.newdat$geog[501:1000], rep(ranef(pc2.model.null)$Sex[2,1],500), col='darkgreen')
lines(PC2.newdat$geog[1001:1500], rep(ranef(pc2.model.null)$Sex[3,1],500), col='grey')

mtext("G.",adj=-0.2,cex=1.25)

plot(new.dat$geog,new.dat$mPC2,pch=19,xlab="Geographical PCoA (LCP)",ylab="",type="n")
points(new.dat$geog[new.dat$Source=="Museum"],new.dat$mPC2[new.dat$Source=="Museum"],pch=19,col='#1f78b4')
points(new.dat$geog[new.dat$Source=="Net"],new.dat$mPC2[new.dat$Source=="Net"],pch=19,col='#b2df8a')
points(PC2.newdat$geog[1:500],rep(fixef(pc2.model.null),500),pch=19, type="l",col="black",lty=2,lwd=3)
lines(PC2.newdat$geog[1:500], rep(ranef(pc2.model.null)$Source[1,1],500),col='#1f78b4')
lines(PC2.newdat$geog[501:1000], rep(ranef(pc2.model.null)$Source[2,1],500),col='#b2df8a')

mtext("H.",adj=-0.2,cex=1.25)

# PC3
plot(new.dat$cPC4,new.dat$mPC3,pch=19,xlab="Climatic PC4",ylab="",type="n")
points(new.dat$cPC4[new.dat$Subspecies=="Manacus manacus"],new.dat$mPC3[new.dat$Subspecies=="Manacus manacus"],pch=19)
points(new.dat$cPC4[new.dat$Subspecies=="Manacus vitellinus"],new.dat$mPC3[new.dat$Subspecies=="Manacus vitellinus"],pch=19,col="goldenrod1")
	points(PC3.newdat$cPC4[1:500],PC3.pred$fit[1:500],pch=19, type="l",col="black",lty=2,lwd=3)
	lines(PC3.newdat$cPC4[1:500], PC3.pred.Sp[1:500], col='black')
	lines(PC3.newdat$cPC4[501:1000], PC3.pred.Sp[501:1000], col='goldenrod1')
mtext("I.",adj=-0.2,cex=1.25)	

plot(new.dat$cPC4,new.dat$mPC3,pch=19,xlab="Climatic PC4",ylab="",type="n")
points(new.dat$cPC4[new.dat$Sex=="Male"],new.dat$mPC3[new.dat$Sex=="Male"])
points(new.dat$cPC4[new.dat$Sex=="Female"],new.dat$mPC3[new.dat$Sex=="Female"],pch=19,col="darkgreen")
points(new.dat$cPC4[new.dat$Sex=="Unknown"],new.dat$mPC3[new.dat$Sex=="Unknown"],pch=19,col="grey")
points(PC3.newdat$cPC4[1:500],PC3.pred$fit[1:500],pch=19, type="l",col="black",lty=2,lwd=3)
lines(PC3.newdat$cPC4[1:500], PC3.pred.Sex[1:500], col='black')
lines(PC3.newdat$cPC4[501:1000], PC3.pred.Sex[501:1000], col='darkgreen')
lines(PC3.newdat$cPC4[1001:1500], PC3.pred.Sex[1001:1500], col='grey')

mtext("J.",adj=-0.2,cex=1.25)

plot(new.dat$cPC4,new.dat$mPC3,pch=19,xlab="Climatic PC4",ylab="",type="n")
points(new.dat$cPC4[new.dat$Source=="Museum"],new.dat$mPC3[new.dat$Source=="Museum"],pch=19,col='#1f78b4')
points(new.dat$cPC4[new.dat$Source=="Net"],new.dat$mPC3[new.dat$Source=="Net"],pch=19,col='#b2df8a')
points(PC3.newdat$cPC4[1:500],PC3.pred$fit[1:500],pch=19, type="l",col="black",lty=2,lwd=3)
lines(PC3.newdat$cPC4[1:500], PC3.pred.Source[1:500],col='#1f78b4')
lines(PC3.newdat$cPC4[501:1000], PC3.pred.Source[501:1000],col='#b2df8a')

mtext("K.",adj=-0.2,cex=1.25)


par(mar=c(0,0,0,0))
plot(0,0,axes=FALSE,xlab="",ylab="",type="n")
legend("top",c("Male","Female","Unknown","white phenotype","yellow phenotype","Museum","Net","Prediction")
		,pch=c(1,rep(19,6),-1),lty=c(rep(1,7),2),lwd=c(rep(1,7),2)
		,col=c("black","darkgreen","grey","black","goldenrod1",'#1f78b4','#b2df8a')
		,bty="n",ncol=4,cex=1.2,pt.cex=1.1)
#legend("top","Prediction",lty=2,bty="n",cex=1.25,lwd=2)

mtext("Body Mass (g)",side=2,outer=TRUE,adj=0.9,line=0.75,cex=0.9)
mtext("Morphological PC1",side=2,outer=TRUE,adj=0.625,line=0.75,cex=0.9)
mtext("Morphological PC2",side=2,outer=TRUE,adj=0.36,line=0.75,cex=0.9)
mtext("Morphological PC3",side=2,outer=TRUE,adj=0.06,line=0.75,cex=0.9)
dev.off()
##### Contact Zones Analysis ##
pop_cz <- c('Acandi', 'Amalfi', 'Anori',
		 'Chigorodo', 'Choco', 'Humaga',
		 'Remedios', 'Riosucio', 'Tacui-Cuni',
		 'Tulenapa', 'Unguia', 'Cimitarra',
		 'El Cucui', 'Gomez Plata', 'Isagen',
		 'Isagen2', 'Maceo','Playas',
		 'Barbacoas', 'San juan')

crs_utm <- "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

which(data$Population == 'Acandi')
pop1 <- dat[142, c('Long', 'Lat')]
dat_sp <- dat
coordinates(pop1) <- ~Long+Lat
coordinates(dat_sp) <- ~Long+Lat
crs(pop1) <- wgs ; crs(dat_sp) <- wgs
pop1 <- spTransform(pop1, crs_utm) ; dat <- spTransform(dat_sp, crs_utm)

dist_pop1 <- gDistance(pop1, dat, byid = T) / 1000

names(data)
names(data)[21] <- 'pop1.D'
data <- cbind(dat, dist_pop1)
write.table(data, "/home/camilo/Dropbox/Tesis Camilo/Data/Data_Manacus.txt", sep='\t', dec='.', row.names=F)

pop_cz <- c('Acandi', 'Amalfi', 'Anori', 'Chigorodo', 'Choco', 'Humaga','Remedios', 'Riosucio', 'Tacui-Cuni',
	'Tulenapa', 'Unguia', 'Cimitarra', 'El Cucui', 'Gomez Plata', 'Isagen', 'Isagen2', 'Maceo',
	'Puerto Salgar', 'Playas')

ind_cz <- data[which((match(data$Population, pop_cz)) > 0),]
grey <- rgb(150,150,150, alpha=200, maxColorValue=255)


png('/home/camilo/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v2/Cz_PC.png', units='cm', height = 10, width=25, res=150)
par(mfrow=c(1,3), family='serif')
plot(PC1~pop1.D, data=ind_cz, col=grey, xlab='')
points(ind_cz$pop1.D[ind_cz$Specie == 'Manacus manacus'], ind_cz$PC1[ind_cz$Specie == 'Manacus manacus'], col='black', pch=19)
points(ind_cz$pop1.D[ind_cz$Specie == 'Manacus vitellinus'], ind_cz$PC1[ind_cz$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)

plot(PC2~pop1.D, data=ind_cz, col=grey, xlab='Distance to Population 1 (m)')
points(ind_cz$pop1.D[ind_cz$Specie == 'Manacus manacus'], ind_cz$PC2[ind_cz$Specie == 'Manacus manacus'], col='black', pch=19)
points(ind_cz$pop1.D[ind_cz$Specie == 'Manacus vitellinus'], ind_cz$PC2[ind_cz$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)

plot(PC3~pop1.D, data=ind_cz, col=grey, xlab='')
points(ind_cz$pop1.D[ind_cz$Specie == 'Manacus manacus'], ind_cz$PC3[ind_cz$Specie == 'Manacus manacus'], col='black', pch=19)
points(ind_cz$pop1.D[ind_cz$Specie == 'Manacus vitellinus'], ind_cz$PC3[ind_cz$Specie == 'Manacus vitellinus'], col='goldenrod1', pch=19)

par(mar=c(0,0,0,0),oma=c(0,0,0,0),fig=c(0,1,0,1),new=TRUE,cex=1.25)
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("top",legend=c("Manacus manacus","Manacus vitellinus"),pch=19,lty=1,col=c("black","goldenrod1"),bty="n",ncol=2)
dev.off()
points

cz <- ()