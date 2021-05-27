# Figure to show differences in morphology between sexes and with green individuals of unknwon sex.
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

new.dat2		<-	new.dat
new.dat2$Sex	<-	as.character(new.dat2$Sex)
new.dat2$Sex[new.dat2$Sex==""]	<-	"Unknown"

tiff("~/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v3/FigureS4.tiff",width=7,height=7.7,units="in"
	,res=300,compression="lzw",type="cairo",family="serif")
par(mfrow=c(2,2), oma=c(1,1,1,1),mar=c(3,3,1,1),mgp=c(1.5,0.25,0),tcl=-0.25)
boxplot(mPC1~Sex,data=new.dat2,ylab="Mrophological Principal Component 1")
mtext("A.",side=3,adj=-0.1)
boxplot(mPC2~Sex,data=new.dat2,ylab="Mrophological Principal Component 2")
mtext("B.",side=3,adj=-0.1)
boxplot(mPC3~Sex,data=new.dat2,ylab="Mrophological Principal Component 3")
mtext("C.",side=3,adj=-0.1)
boxplot(Mass~Sex,data=new.dat2,ylab="Body Mass (g)")
mtext("D.",side=3,adj=-0.1)
dev.off()