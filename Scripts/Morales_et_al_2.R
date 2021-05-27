library(raster)
library(gdm)
library(geosphere)
library(FD)
library(dismo)
library(vegan)
library(FD)

setwd('')

### Data sets ####
data <- read.csv('Pop_Data.txt', header=T, sep='\t')
data <- read.csv('Pop_Data_Male.txt', header=T, sep='\t')
data <- read.csv('Pop_Data_BM.txt', header=T, sep='\t')
data <- read.table('Pop_Data_BM_Male.txt', header=T, sep='\t')

#### 52 Populations | Morphological Distances ####

## PC1 ##
PC1 <- as.matrix(data$PC1)
gPC1 <- as.matrix(vegdist(PC1, 'gower'))
Locality_ID <- 1:dim(data)[1]
gPC1 <- cbind(Locality_ID, gPC1)

## culmen ##
culmen <- as.matrix(data$culmen)
gculmen <- as.matrix(vegdist(culmen, 'gower'))
Locality_ID <- 1:dim(data)[1]
gculmen <- cbind(Locality_ID, gculmen)

## Wing ##
Wing <- as.matrix(data$Wing_mean)
gWing <- as.matrix(vegdist(Wing, 'gower'))
Locality_ID <- 1:dim(data)[1]
gWing <- cbind(Locality_ID, gWing)

## Tail ##
Tail <- as.matrix(data$Tail_mean)
gTail <- as.matrix(vegdist(Tail, 'gower'))
Locality_ID <- 1:dim(data)[1]
gTail <- cbind(Locality_ID, gTail)

## Culmen ##
Culmen <- as.matrix(data$Exposed_Culmen_mean)
gCulmen <- as.matrix(vegdist(Culmen, 'gower'))
Locality_ID <- 1:dim(data)[1]
gCulmen <- cbind(Locality_ID, gCulmen)

## Tarsus ##
Tarsus <- as.matrix(data$Tarsus_mean)
gTarsus <- as.matrix(vegdist(Tarsus, 'gower'))
Locality_ID <- 1:dim(data)[1]
gTarsus <- cbind(Locality_ID, gTarsus)


### 32 Populations | Body mass Matrix ###
## Body_mass ##

data <- read.csv('/home/camilo/Documentos/Manacus_Articulo/Data/Pop_Data_BM.txt', header=T, sep='\t')

BM <- as.matrix(data$Body_mass_mean)
gBM <- as.matrix(vegdist(BM, 'gower'))
Locality_ID <- 1:dim(data)[1]
gBM <- cbind(Locality_ID, gBM)


############ GDM Models #################

## PC1 ##
gdm.dis.pc1 <- formatsitepair(gPC1, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=data[,18:26], siteColumn='Locality_ID')
gdm.pc1 <- gdm(gdm.dis.pc1, geo=T)
gdm.pc1$explained

## PC2 ##
gdm.dis.pc2 <- formatsitepair(gPC2, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=data[,18:26], siteColumn='Locality_ID')
gdm.pc2 <- gdm(gdm.dis.pc2, geo=T)
gdm.pc2$explained

## Culmen ##
gdm.dis.culmen <- formatsitepair(gculmen, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=data[,18:26], siteColumn='Locality_ID')
gdm.culmen <- gdm(gdm.dis.culmen, geo=T)
gdm.culmen$explained

## Wing ##
gdm.dis.wing <- formatsitepair(gWing, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=data[,18:26], siteColumn='Locality_ID')
gdm.wing <- gdm(gdm.dis.wing, geo=T)
gdm.wing$explained

## Tail ##
gdm.dis.tail <- formatsitepair(gTail, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=data[,18:26], siteColumn='Locality_ID')
gdm.tail <- gdm(gdm.dis.tail, geo=T)
gdm.tail$explained

## Tarsus ##
gdm.dis.tarsus <- formatsitepair(gTarsus, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=data[,18:26], siteColumn='Locality_ID')
gdm.tarsus <- gdm(gdm.dis.tarsus, geo=T)
gdm.tarsus$explained

## Culmen ##
gdm.dis.culmen <- formatsitepair(gCulmen, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=data[,18:26], siteColumn='Locality_ID')
gdm.culmen <- gdm(gdm.dis.culmen, geo=T)
gdm.culmen$explained

## Body mass ##
gdm.dis.BM <- formatsitepair(gBM, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=data[,20:28], siteColumn='Locality_ID')
gdm.BM <- gdm(gdm.dis.BM, geo=T)
gdm.BM$explained


########################################################################################################
#########################################################################################################
#### Manacus Morphological Matrices #######

manacus <- which(data$Specie == 'Manacus manacus')
pop_mana <- data[manacus,]
pop_mana[,'Locality_ID'] <- 1:dim(pop_mana)[1]

## PC1 ##
PC1_m <- as.matrix(pop_mana$PC1)
gPC1_m <- as.matrix(vegdist(PC1_m, 'gower'))
Locality_ID <- 1:dim(pop_mana)[1]
gPC1_m <- cbind(Locality_ID, gPC1_m)

## PC2 ##
PC2_m <- as.matrix(pop_mana$PC2)
gPC2_m <- as.matrix(vegdist(PC2_m, 'gower'))
Locality_ID <- 1:dim(pop_mana)[1]
gPC2_m <- cbind(Locality_ID, gPC2_m)

## Wing ##
Wing_m <- as.matrix(pop_mana$Wing_mean)
gWing_m <- as.matrix(vegdist(Wing_m, 'gower'))
Locality_ID <- 1:dim(pop_mana)[1]
gWing_m <- cbind(Locality_ID, gWing_m)

## Tail ##
Tail_m <- as.matrix(pop_mana$Tail_mean)
gTail_m <- as.matrix(vegdist(Tail_m, 'gower'))
Locality_ID <- 1:dim(pop_mana)[1]
gTail_m <- cbind(Locality_ID, gTail_m)

## Culmen ##
Culmen_m <- as.matrix(pop_mana$Exposed_Culmen_mean)
gCulmen_m <- as.matrix(vegdist(Culmen_m, 'gower'))
Locality_ID <- 1:dim(pop_mana)[1]
gCulmen_m <- cbind(Locality_ID, gCulmen_m)

## Tarsus ##
Tarsus_m <- as.matrix(pop_mana$Tarsus_mean)
gTarsus_m <- as.matrix(vegdist(Tarsus_m, 'gower'))
Locality_ID <- 1:dim(pop_mana)[1]
gTarsus_m <- cbind(Locality_ID, gTarsus_m)

## Body mass ##
data <- read.csv('Pop_Data_BM.txt', header=T, sep='\t')

BM_m <- as.matrix(pop_mana$Body_mass_mean)
gBM_m <- as.matrix(vegdist(BM_m, 'gower'))
Locality_ID <- 1:dim(pop_mana)[1]
gBM_m <- cbind(Locality_ID, gBM_m)


###################### manacus GDM models #######################################################
names(pop_mana)
## PC1 ##
gdm.dis.pc1.m <- formatsitepair(gPC1_m, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_mana[,18:26], siteColumn='Locality_ID')
gdm.pc1.m <- gdm(gdm.dis.pc1.m, geo=T)
gdm.pc1.m$explained

## PC2 ##
gdm.dis.pc2.m <- formatsitepair(gPC2_m, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_mana[,18:26], siteColumn='Locality_ID')
gdm.pc2.m <- gdm(gdm.dis.pc2.m, geo=T)
gdm.pc2.m$explained

## Culmen ##
gdm.dis.culmen.m <- formatsitepair(gCulmen_m, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_mana[,18:26], siteColumn='Locality_ID')
gdm.culmen.m <- gdm(gdm.dis.culmen.m, geo=T)
gdm.culmen.m$explained

## Wing ##
gdm.dis.wing.m <- formatsitepair(gWing_m, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_mana[,18:26], siteColumn='Locality_ID')
gdm.wing.m <- gdm(gdm.dis.wing.m, geo=T)
gdm.wing.m$explained

## Tail ##
gdm.dis.tail.m <- formatsitepair(gTail_m, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_mana[,18:26], siteColumn='Locality_ID')
gdm.tail.m <- gdm(gdm.dis.tail.m, geo=T)
gdm.tail.m$explained

## Tarsus ##
gdm.dis.tarsus.m <- formatsitepair(gTarsus_m, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_mana[,18:26], siteColumn='Locality_ID')
gdm.tarsus.m <- gdm(gdm.dis.tarsus.m, geo=T)
gdm.tarsus.m$explained

## Body mass ##
gdm.dis.BM.m <- formatsitepair(gBM_m, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_mana[,20:28], siteColumn='Locality_ID')
gdm.BM.m <- gdm(gdm.dis.BM.m, geo=T)
gdm.BM.m$explained

################################################################################################
#### vitellinus morphological matrices ########

vitellinus <- which(data$Specie == 'Manacus vitellinus')
pop_vite <- data[vitellinus,]
pop_vite[,'Locality_ID'] <- 1:dim(pop_vite)[1]

## PC1 ##
PC1_v <- as.matrix(pop_vite$PC1)
gPC1_v <- as.matrix(vegdist(PC1_v, 'gower'))
Locality_ID <- 1:dim(pop_vite)[1]
gPC1_v <- cbind(Locality_ID, gPC1_v)

## PC2 ##
PC2_v <- as.matrix(pop_vite$PC2)
gPC2_v <- as.matrix(vegdist(PC2_v, 'gower'))
Locality_ID <- 1:dim(pop_vite)[1]
gPC2_v <- cbind(Locality_ID, gPC2_v)


## Wing ##
Wing_v <- as.matrix(pop_vite$Wing_mean)
gWing_v <- as.matrix(vegdist(Wing_v, 'gower'))
Locality_ID <- 1:dim(pop_vite)[1]
gWing_v <- cbind(Locality_ID, gWing_v)

## Tail ##
Tail_v <- as.matrix(pop_vite$Tail_mean)
gTail_v <- as.matrix(vegdist(Tail_v, 'gower'))
Locality_ID <- 1:dim(pop_vite)[1]
gTail_v <- cbind(Locality_ID, gTail_v)

## Culmen ##
Culmen_v <- as.matrix(pop_vite$Exposed_Culmen_mean)
gCulmen_v <- as.matrix(vegdist(Culmen_v, 'gower'))
Locality_ID <- 1:dim(pop_vite)[1]
gCulmen_v <- cbind(Locality_ID, gCulmen_v)

## Tarsus ##
Tarsus_v <- as.matrix(pop_vite$Tarsus_mean)
gTarsus_v <- as.matrix(vegdist(Tarsus_v, 'gower'))
Locality_ID <- 1:dim(pop_vite)[1]
gTarsus_v <- cbind(Locality_ID, gTarsus_v)

## Body mass ##
data <- read.csv('Pop_Data_BM.txt', header=T, sep='\t')

BM_v <- as.matrix(pop_vite$Body_mass_mean)
gBM_v <- as.matrix(vegdist(BM_v, 'gower'))
Locality_ID <- 1:dim(pop_vite)[1]
gBM_v <- cbind(Locality_ID, gBM_v)


#### vitellinus GDM Models ###

gdm.dis.pc1.v <- formatsitepair(gPC1_v, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_vite[,18:26], siteColumn='Locality_ID')
gdm.pc1.v <- gdm(gdm.dis.pc1.v, geo=T)
gdm.pc1.v$explained

## PC2 ##
gdm.dis.pc2.v <- formatsitepair(gPC2_v, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_vite[,18:26], siteColumn='Locality_ID')
gdm.pc2.v <- gdm(gdm.dis.pc2.v, geo=T)
gdm.pc2.v$explained

## Culmen ##
gdm.dis.culmen.v <- formatsitepair(gCulmen_v, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_vite[,18:26], siteColumn='Locality_ID')
gdm.culmen.v <- gdm(gdm.dis.culmen.v, geo=T)
gdm.culmen.v$explained

## Wing ##
gdm.dis.wing.v <- formatsitepair(gWing_v, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_vite[,18:26], siteColumn='Locality_ID')
gdm.wing.v <- gdm(gdm.dis.wing.v, geo=T)
gdm.wing.v$explained

## Tail ##
gdm.dis.tail.v <- formatsitepair(gTail_v, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_vite[,18:26], siteColumn='Locality_ID')
gdm.tail.v <- gdm(gdm.dis.tail.v, geo=T)
gdm.tail.v$explained

## Tarsus ##
gdm.dis.tarsus.v <- formatsitepair(gTarsus_v, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_vite[,18:26], siteColumn='Locality_ID')
gdm.tarsus.v <- gdm(gdm.dis.tarsus.v, geo=T)
gdm.tarsus.v$explained

## Body mass ##
gdm.dis.BM.v <- formatsitepair(gBM_v, bioFormat=3, XColumn='Longitude', YColumn='Latitude', 
	predData=pop_vite[,20:28], siteColumn='Locality_ID')
gdm.BM.v <- gdm(gdm.dis.BM.v, geo=T)
gdm.BM.v$explained