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
load("data/regions.RData")
load("data/departe.RData")
load("data/regionsLight.RData")
load("data/departeLight.RData")

#load the resistance results for the 2020 campaign
databruteTOT<-read.delim(
  "data/data_carto_france.txt",
  header=TRUE,
  sep="\t",
  colClasses=c("character","numeric","numeric",
               "numeric","numeric","numeric","numeric",
               "numeric","numeric","numeric","numeric",
               "numeric","numeric","numeric","factor",
               "factor","factor","factor")
)
databruteTOT$RS<-rowSums(databruteTOT[,4:14])
levels(databruteTOT$SeqMeth)<-c(19,17)

#turning this dataframe into a spatial dataframe (wgs84)
ambro.wgs<-SpatialPointsDataFrame(coords=databruteTOT[,c(3,2)],
                                  data=databruteTOT,
                                  proj4string=CRS("+proj=longlat +datum=WGS84")
)
ambro<-spTransform(ambro.wgs,CRS("+init=epsg:2154"))


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

#one single map without R/S information
op<-par(mar=c(0,0,1,0))
#original screening
plot(departe,lwd=0.8,border=grey(0.7))
plot(regions,lwd=1.8,add=TRUE)
plot(ambro[ambro$SampRound!="Screening",],
     pch=0,
     col=rgb(31,133,235,175,maxColorValue=255),cex=1,
     add=TRUE)
plot(ambro[ambro$SampRound=="Screening",],
     pch=1,
     col=rgb(213,179,0,175,maxColorValue=255),cex=1,
     add=TRUE)
legend(100000,7120000,
       legend=c("Targeted sampling (N=43)","Random sampling (N=212)"),
       cex=0.8,pt.cex=1,y.intersp=1.2,x.intersp=0.8,
       pch=c(1,0),
       col=c(rgb(213,179,0,175,maxColorValue=255),
             rgb(31,133,235,175,maxColorValue=255)),
       bg="transparent",bty="n")
scalebar(c(191260,6060000),300000,"km",division.cex=1)
text(594045.1,6502086,labels="NAQ",cex=0.9,font=2,col=grey(0.3))
text(594045.1,6751660,labels="CVL",cex=0.9,font=2,col=grey(0.3))
text(742229.5,6515085,labels="ARA",cex=0.9,font=2,col=grey(0.3))
text(708433,6333104,labels="OCC",cex=0.9,font=2,col=grey(0.3))
text(412064,6720000,labels="PDL",cex=0.9,font=2,col=grey(0.3))
par(op)

#export to .pdf 5 x 5 inches


##############################################################################/
#END
##############################################################################/