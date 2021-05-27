library(geosphere)
dat	<-	read.csv("~/Dropbox/Tesis Camilo/Data/Data_Manacus.csv")
obs.pcoa	<-	matrix(NA,nrow=1000,ncol=3)

sd.vec	<-	c(1,5,10)	

for(i in 1:1000){
	for(j in 1:3){
		lat	<-	runif(311,min=min(dat$Lat),max=max(dat$Lat))
		long<-	runif(311,min=min(dat$Long),max=max(dat$Long))

		pos	<-	cbind(long,lat)
		geo.dist.sim	<-	distm(pos)/1000

	
		morph.dist.sim	<-	matrix(abs(rnorm(length(geo.dist.sim),mean=0.02*geo.dist.sim,sd=sd.vec[j])),ncol=311)
		diag(morph.dist.sim)	<-	0
		geo.dist.sim	<-	as.dist(geo.dist.sim,diag=TRUE)
		morph.dist.sim	<-	as.dist(morph.dist.sim,diag=TRUE)

		morph.pcoa	<-	cmdscale(morph.dist.sim,k=1)
		dist.pcoa	<-	cmdscale(geo.dist.sim,k=1)

		obs.pcoa[i,j]	<-	abs(coef(lm(morph.pcoa~dist.pcoa)))[2]
	}
}

new.slopes	<-	rbind(rep(0.02,3),obs.pcoa)
p.vals		<-	apply(new.slopes,2,function(x)rank(x)[1])/1001

tiff("~/Dropbox/Tesis Camilo/Articulo_Manacus/Figures/v2/Hist_Simulations.tiff"
	,width=8,height=4,units="in",res=300,compression="lzw",type="cairo",family="serif")
par(mfrow=c(1,3),mar=c(3,3,1,1),oma=c(2,2,1,1),mgp=c(0,0.75,0),tcl=-0.25)
hist(obs.pcoa[,1],col="black",border="white",xlab="",ylab="",main="",xaxp=c(0.01994,0.02004,2))
box(bty="l")
mtext("SD = 1",3,adj=0,font=2)
hist(obs.pcoa[,2],col="black",border="white",xlab="",ylab="",main="",xaxp=c(0.0198,0.0203,2))
box(bty="l")
mtext("SD = 5",3,adj=0,font=2)
hist(obs.pcoa[,3],col="black",border="white",xlab="",ylab="",main="",xaxp=c(0.0194,0.0206,2))
box(bty="l")
mtext("SD = 10",3,adj=0,font=2)
mtext("Slope of Morphology vs. Geographic PCoA",1,outer=TRUE)
mtext("Frequency",2,outer=TRUE)
dev.off()