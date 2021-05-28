##############################################################################/
##############################################################################/
#Code for plotting the phenotype maps
##############################################################################/
##############################################################################/

#loading the necessary packages
library(rgdal)
library(rgeos)
library(plotrix)
library(mapplots)
library(raster)
library(RColorBrewer)


#load geographical data
load("data/regions.RData")
load("data/departe.RData")
load("data/regionsLight.RData")
load("data/departeLight.RData")

#changing the projection of the map to WGS84
departe.wgs<-spTransform(departe,
                         CRS("+proj=longlat +datum=WGS84"))
#changing the projection of the map
regions.wgs<-spTransform(regions,
                         CRS("+proj=longlat +datum=WGS84"))
departeL.wgs<-spTransform(departeLight,
                          CRS("+proj=longlat +datum=WGS84"))
#changing the projection of the map
regionsL.wgs<-spTransform(regionsLight,
                          CRS("+proj=longlat +datum=WGS84"))

#load the resistance results for the 2020 campaign
databruteTOT<-read.delim(
  "data/data_carto_pheno.txt",
  header=TRUE,
  sep="\t",
  colClasses=c("character","numeric","numeric",
               "factor","factor","factor","factor",
               "factor")
)
levels(databruteTOT$SeqMeth)<-c(19,17)

PheOri<-databruteTOT[databruteTOT$SampRound=="origin",]
PheNew<-databruteTOT[databruteTOT$SampRound=="newLoc",]

#turning this dataframe into a spatial dataframe (wgs84)
PheOri.wgs<-SpatialPointsDataFrame(coords=PheOri[,c(3,2)],
                                  data=PheOri,
                                  proj4string=
                                    CRS("+proj=longlat +datum=WGS84")
)
PheOri<-spTransform(PheOri.wgs,CRS("+init=epsg:2154"))


#turning this dataframe into a spatial dataframe (wgs84)
PheNew.wgs<-SpatialPointsDataFrame(coords=PheNew[,c(3,2)],
                                   data=PheNew,
                                   proj4string=
                                     CRS("+proj=longlat +datum=WGS84")
)
PheNew<-spTransform(PheNew.wgs,CRS("+init=epsg:2154"))

#to catch the coordinates
#locator()


##############################################################################/
#defining additional function for the mapping####
##############################################################################/


#function for a scale, found in "Auxiliary Cartographic Functions in R: 
#North Arrow, Scale Bar, and Label with a Leader Arrow", Tanimura et al 2007, 
#J of Statistical software
#The code has been slightly modified in order to convert the meter in km
scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])/1000
  labels <- append(labels,paste((x[5]-loc[1])/1000,unit))
  text(x[c(1,3,5)],y[4],labels=labels,adj=c(0.5,0),cex=division.cex)
}


##############################################################################/
#code for the maps by mutation in WGS####
##############################################################################/

#one map for the 4 possibilities
op<-par(mar=c(0,1,0,0))
#Imazamox and/or Tribenuron phenotypes
plot(departe,lwd=0.8,border=grey(0.9))
plot(regions,lwd=1.8,add=TRUE)
segments(PheNew$Longitude,PheNew$Latitude,
         PheOri$Longitude,PheOri$Latitude,
         col=grey(0.4,1),lwd=0.7)
plot(PheOri,pch=21,cex=0.3,bg=grey(0.4,1),col="black",add=TRUE)
plot(PheNew[(PheNew$Imz=="S" & PheNew$Tbn=="S")|
             (PheNew$Imz=="S" & PheNew$Tbn=="Fail")|
             (PheNew$Imz=="Fail" & PheNew$Tbn=="S"),],
     pch=21,col="transparent",
     bg=rgb(80,140,20,255,maxColorValue=255),cex=1.2,
     add=TRUE)
plot(PheNew[PheNew$Imz=="R" & (PheNew$Tbn=="S"|PheNew$Tbn=="Fail"),],
     pch=21,col="transparent",
     bg=rgb(244,114,208,255,maxColorValue=255),cex=1.2,
     add=TRUE)
plot(PheNew[PheNew$Tbn=="R" & (PheNew$Imz=="S"|PheNew$Imz=="Fail"),],
     pch=21,col="transparent",
     bg=rgb(250,104,0,255,maxColorValue=255),cex=1.2,
     add=TRUE)
plot(PheNew[PheNew$Imz=="R" & PheNew$Tbn=="R",],
     pch=21,col="transparent",
     bg=rgb(220,31,31,255,maxColorValue=255),cex=1.2,
     add=TRUE)
plot(PheNew[PheNew$RS=="Mut",],
     pch=21,
     col=rgb(10,60,200,255,maxColorValue=255),
     bg="transparent",cex=1.2,
     add=TRUE)
plot(PheNew[PheNew$RS=="Mut",],
     pch=21,
     col=rgb(10,60,200,255,maxColorValue=255),
     bg="transparent",cex=1.0,
     add=TRUE)
legend(10000,6750000,
       legend=c("No resistant\nplants detected",
                "Detection of resistant\nplants with imazamox",
                "Detection of resistant\nplants with tribenuron",
                "Detection of resistant\nplants with both herbicides",
                "Mutation of the ALS\ngene detected"),
       cex=0.6,pt.cex=1,y.intersp=2.3,x.intersp=1,
       pch=c(19,19,19,19,1),
       col=c(rgb(80,140,20,255,maxColorValue=255),
             rgb(244,114,208,255,maxColorValue=255),
             rgb(250,104,0,255,maxColorValue=255),
             rgb(220,31,31,255,maxColorValue=255),
             rgb(10,60,200,255,maxColorValue=255)),
       bg="transparent",bty="n",xpd=TRUE)
scalebar(c(191260,6060000),300000,"km",division.cex=0.8)
text(594045.1,6502086,labels="NAQ",cex=0.9,font=2,col=grey(0.3))
text(594045.1,6751660,labels="CVL",cex=0.9,font=2,col=grey(0.3))
text(742229.5,6515085,labels="ARA",cex=0.9,font=2,col=grey(0.3))
text(708433,6333104,labels="OCC",cex=0.9,font=2,col=grey(0.3))
text(412064,6720000,labels="PDL",cex=0.9,font=2,col=grey(0.3))
text(870000,6698000,labels="BFC",cex=0.9,font=2,col=grey(0.3))
# #if you want to add the name of the populations
# text(x=PheNew$Longitude,y=PheNew$Latitude,
#      labels=PheNew@data$Code_ID,cex=0.2)
par(op)

#export to .pdf 5 x 5 inches


##############################################################################/
#END
##############################################################################/