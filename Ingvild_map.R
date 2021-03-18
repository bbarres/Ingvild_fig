##############################################################################/
##############################################################################/
#Code for plotting the maps
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
load("data/regionsLight.RData")
load("data/departeLight.RData")
#changing the projection of the map
departe.wgs <- spTransform(departe,
                           CRS("+proj=longlat +datum=WGS84"))
#changing the projection of the map
regions.wgs <- spTransform(regions,
                           CRS("+proj=longlat +datum=WGS84"))

#crop a subpart of map
departe.wgs.1<-crop(departe.wgs,extent(-1.25,1.57,44.9,46.6))
departe.wgs.2<-crop(departe.wgs,extent(2.13,4.95,46.5,48.2))
departe.wgs.3<-crop(departe.wgs,extent(3.59,6.41,44.3,46.0))
departe.wgs.4<-crop(departe.wgs,extent(0.15,2.97,42.9,44.6))
regions.wgs.1<-crop(regions.wgs,extent(-1.25,1.57,44.9,46.6))
regions.wgs.2<-crop(regions.wgs,extent(2.13,4.95,46.5,48.2))
regions.wgs.3<-crop(regions.wgs,extent(3.59,6.41,44.3,46.0))
regions.wgs.4<-crop(regions.wgs,extent(0.15,2.97,42.9,44.6))

#load the resistance results for the 2020 campaign
databruteTOT<-read.delim(
  "data/data_carto_france.txt",
  header = TRUE,
  sep = "\t",
  colClasses = c("character","numeric","numeric",
                 "numeric","numeric","numeric","numeric",
                 "numeric","numeric","numeric","numeric",
                 "numeric","numeric","factor")
)
databruteTOT$RS<-rowSums(databruteTOT[,4:13])
levels(databruteTOT$SeqMeth)<-c(19,17)

#turning this dataframe into a spatial dataframe
ambroFiel<-SpatialPointsDataFrame(coords=databruteTOT[,c(3,2)],
                                  data=databruteTOT,
                                  proj4string=CRS("+proj=longlat +datum=WGS84")
                                  )

#to catch the coordinates
#locator()

#defining a vector of 5 colors, one for each population
colopop<-brewer.pal(5,"Dark2")


##############################################################################/
#code for the maps by mutation####
##############################################################################/

#all mutations on one map
op<-par(mar=c(0,0,0,0))
plot(departe.wgs,lwd=0.8,border=grey(0.7))
plot(regions.wgs,lwd=1.8,add=TRUE)
plot(ambroFiel[ambroFiel$RS==0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=0.8,
     add=TRUE)
plot(ambroFiel[ambroFiel$RS!=0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=0.8,
     add=TRUE)
polygon(x=c(-1.25,1.57,1.57,-1.25,-1.25),
        y=c(46.6,46.6,44.9,44.9,46.6),
        col="transparent",lwd=3,border="lightblue")
polygon(x=c(2.13,4.95,4.95,2.13,2.13),
        y=c(48.2,48.2,46.5,46.5,48.2),
        col="transparent",lwd=3,border="lightblue")
polygon(x=c(3.59,6.41,6.41,3.59,3.59),
        y=c(46.0,46.0,44.3,44.3,46.0),
        col="transparent",lwd=3,border="lightblue")
polygon(x=c(0.15,2.97,2.97,0.15,0.15),
        y=c(44.6,44.6,42.9,42.9,44.6),
        col="transparent",lwd=3,border="lightblue")
par(op)

#export to .pdf 7 x 7 inches

#all mutations on one map
op<-par(mar=c(0,0,1,0),mfrow=c(2,5))
for (i in 4:13) {
  plot(departe.wgs,lwd=0.8,border=grey(0.7),
       main=colnames(ambroFiel@data)[i])
  plot(regions.wgs,lwd=1.8,add=TRUE)
  plot(ambroFiel[ambroFiel@data[,i]==0,],
       pch=as.numeric(as.character(ambroFiel[ambroFiel@data[,i]==0,]$SeqMeth)),
       col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
       add=TRUE)
  plot(ambroFiel[ambroFiel@data[,i]!=0,],
       pch=as.numeric(as.character(ambroFiel[ambroFiel@data[,i]!=0,]$SeqMeth)),
       col=rgb(255,0,0,255,maxColorValue=255),cex=1.5,
       add=TRUE)
}
par(op)

#export to .pdf 16 x 7 inches


##############################################################################/
#code with cropped maps####
##############################################################################/

op<-par(mar=c(0,0,0,0))
plot(departe.wgs.1,lwd=0.8,border=grey(0.7))
plot(regions.wgs.1,lwd=2.5,add=TRUE)
plot(ambroFiel[ambroFiel$RS==0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambroFiel[ambroFiel$RS!=0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
par(op)

op<-par(mar=c(0,0,0,0))
plot(departe.wgs.2,lwd=0.8,border=grey(0.7))
plot(regions.wgs.2,lwd=2.5,add=TRUE)
plot(ambroFiel[ambroFiel$RS==0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambroFiel[ambroFiel$RS!=0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
par(op)

op<-par(mar=c(0,0,0,0))
plot(departe.wgs.3,lwd=0.8,border=grey(0.7))
plot(regions.wgs.3,lwd=2.5,add=TRUE)
plot(ambroFiel[ambroFiel$RS==0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambroFiel[ambroFiel$RS!=0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
par(op)

op<-par(mar=c(0,0,0,0))
plot(departe.wgs.4,lwd=0.8,border=grey(0.7))
plot(regions.wgs.4,lwd=2.5,add=TRUE)
plot(ambroFiel[ambroFiel$RS==0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambroFiel[ambroFiel$RS!=0,],
     pch=as.numeric(as.character(ambroFiel[ambroFiel$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
par(op)



##############################################################################/
#END
##############################################################################/