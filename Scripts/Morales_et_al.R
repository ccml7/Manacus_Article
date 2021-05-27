library(raster)
library(gdm)
library(geosphere)
library(FD)
library(dismo)
manacus		<-	read.delim("~/Dropbox/Tesis Camilo/Data/Data_Manacus.txt")
manacus2	<-	read.delim("~/Dropbox/Tesis Camilo/Data/Pop_Data.txt")
wc.files	<-	list.files("~/Documents/GIS/wc2_10km",full.names=TRUE)[c(1,3,12,15,18)]
env.files	<-	list.files("~/Documents/GIS/earthenv",full.names=TRUE)
dem.file	<-	list.files("~/Documents/GIS/Elevation_1km",full.names=TRUE)

clim.vars	<-	stack(wc.files)

clim.locs	<-	extract(clim.vars,manacus2[,c("Longitude","Latitude")])
clim.locs	<-	scale(clim.locs)
colnames(clim.locs)	<-	c("bio1","bio3","bio12","bio15","bio18")

# pc1.means	<-	tapply(manacus[,"PC1"],manacus$Population,FUN="mean",na.rm=TRUE)
# pc2.means	<-	tapply(manacus[,"PC2"],manacus$Population,FUN="mean",na.rm=TRUE)
# pc3.means	<-	tapply(manacus[,"PC3"],manacus$Population,FUN="mean",na.rm=TRUE)
# BM.means	<-	tapply(manacus[,"Body_mass"],manacus$Population,FUN="mean",na.rm=TRUE)
# BM.means	<-	BM.means[is.finite(BM.means)]

# BM.males	<-	tapply(manacus[manacus$Sex=="Male","Body_mass"],manacus$Population[manacus$Sex=="Male"],FUN="mean",na.rm=TRUE)
# BM.males	<-	BM.males[is.finite(BM.males)]

pc1.df		<-	data.frame(manacus2[,c("Locality_ID","Longitude","Latitude")],clim.locs)
pc1.dist	<-	as.matrix(gowdis(matrix(manacus2$PC1,ncol=1)))
pc1.dist	<-	data.frame(Locality_ID=manacus2$Locality_ID,pc1.dist)

# BM.dist		<-	as.matrix(gowdis(matrix(BM.means,ncol=1)))
# BM.dist		<-	data.frame(Population=names(BM.means),BM.dist)
# BM.df		<-	pc1.df[match(names(BM.means),pc1.df$Population),]
# BMM.dist	<-	as.matrix(gowdis(matrix(BM.males,ncol=1)))
# BMM.dist	<-	data.frame(Population=names(BM.males),BMM.dist)
# BMM.df		<-	pc1.df[match(names(BM.males),pc1.df$Population),]

pc1.gdm.dat	<-	formatsitepair(bioData=pc1.dist,bioFormat=3,siteColumn="Locality_ID",XColumn="Longitude",YColumn="Latitude",predData=pc1.df)
BM.gdm.dat	<-	formatsitepair(bioData=BM.dist,bioFormat=3,siteColumn="Population",XColumn="Long",YColumn="Lat",predData=BM.df)
BMM.gdm.dat	<-	formatsitepair(bioData=BMM.dist,bioFormat=3,siteColumn="Population",XColumn="Long",YColumn="Lat",predData=BMM.df[,1:8])

mod.pc1<-gdm(pc1.gdm.dat,geo=TRUE)
mod.BM	<-	gdm(BM.gdm.dat,geo=F)
mod.BMM	<-	gdm(BMM.gdm.dat,geo=TRUE)

anovas	<-	apply(morpho,2,function(x)aov(x~Sex+Species))
anovas2	<-	apply(morpho.pc,2,function(x)aov(x~sex))
lapply(anovas,summary)
lapply(anovas2,summary)

plot(0,0,xlim=c(0,14),ylim=c(1,5))
boxplot(morpho$Wing~sex,at=1:2,add=TRUE)
boxplot(morpho$Tail~sex,at=4:5,add=TRUE)
boxplot(morpho$Exposed_Culmen~sex,at=6:7,add=TRUE)
boxplot(morpho$Tarsus~sex,at=9:10,add=TRUE)
boxplot(morpho$Body_mass~sex,at=12:13,add=TRUE)

library(vegan)
library(ape)
library(geosphere)

dist.mat	<-	distm(manacus[,c("Long","Lat")])
xy.vec		<-	pcoa(dist.mat)

mod.pc1	<-	aov(PC1~Sex+Species,data=manacus)
mod.pc2	<-	aov(PC2~Sex+Species,data=manacus)



summary(lm(morpho.pc$PC1~xy.vec$values[,1]+sex-1))
summary(lm(morpho.pc$PC2~xy.vec$values[,1]+sex-1))
summary(lm(morpho.pc$PC3~xy.vec$values[,1]+sex-1))

library(lme4)
mod	<-	lmer(PC1~Sex-1+1|Population,data=manacus)
mod1<-	lm(PC1~Sex-1,data=manacus)

boxplot(PC1~Population*Sex,data=manacus)
?pcoa



range(manacus$PC2)

plot(0,0,xlim=c(-4,4),ylim=c(-4,4),type="n")
points(manacus$PC1[manacus$Sex=="Male"],manacus$PC2[manacus$Sex=="Male"])
points(manacus$PC1[manacus$Sex=="Female"],manacus$PC2[manacus$Sex=="Female"],pch=19,col="forestgreen")
points(manacus$PC1[which(is.na(manacus$Sex))],manacus$PC2[which(is.na(manacus$Sex))],pch=19)





table(manacus$Population,manacus$Sex)