library(raster)
library(foreign)

lat.elev	<-	raster("/Volumes/JUANP/GIS/LatinAmerica/elev_1x1_latam.tif")
x.lim	<-	range(pops$Long)
y.lim	<-	range(pops$Lat)

manacus		<-	unique(dat.Zone[dat.Zone$Specie=="Manacus manacus","Population"])
vitellinus	<-	unique(dat.Zone[dat.Zone$Specie=="Manacus vitellinus","Population"])

pops$Manacus	<-	pops$Population%in%manacus
pops$Vitellinus	<-	pops$Population%in%vitellinus

plot(lat.elev,xlim=x.lim,ylim=y.lim)
points(pops$Long[pops$Vitellinus],pops$Lat[pops$Vitellinus]
		,pch=19,col="red")
points(pops$Long[pops$Manacus],pops$Lat[pops$Manacus]
		,pch=19,col="black")
		
plot(morph.pca$x[,1],morph.pca$x[,2],type="n")
points(morph.pca$x[dat$Specie=="Manacus manacus",1],morph.pca$x[dat$Specie=="Manacus manacus",2],pch=19)
points(morph.pca$x[dat$Specie=="Manacus vitellinus",1],morph.pca$x[dat$Specie=="Manacus vitellinus",2],pch=19,col="goldenrod")
