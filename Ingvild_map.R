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

#crop a subparts of the map
departe.1<-crop(departe,extent(364666.8,590540.8,6431338.8,6612096.8))
departe.2<-crop(departe,extent(633281.3,844866.6,6600367.6,6790638.4))
departe.3<-crop(departe,extent(747072.2,963827.9,6355840.0,6550171.5))
departe.4<-crop(departe,extent(467049.9,697619.1,6204268.4,6388990.7))
departe.5<-crop(departe,extent(487049.9,602334.5,6306629.5,6398990.7))
regions.1<-crop(regions,extent(364666.8,590540.8,6431338.8,6612096.8))
regions.2<-crop(regions,extent(633281.3,844866.6,6600367.6,6790638.4))
regions.3<-crop(regions,extent(747072.2,963827.9,6355840.0,6550171.5))
regions.4<-crop(regions,extent(467049.9,697619.1,6204268.4,6388990.7))
regions.5<-crop(regions,extent(487049.9,602334.5,6306629.5,6398990.7))

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

#all mutations on one map
op<-par(mar=c(0,0,0,0))
plot(departe.wgs,lwd=0.8,border=grey(0.7))
plot(regions.wgs,lwd=1.8,add=TRUE)
plot(ambro.wgs[ambro.wgs$RS==0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=0.8,
     add=TRUE)
plot(ambro.wgs[ambro.wgs$RS!=0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS!=0,]$SeqMeth)),
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
op<-par(mar=c(0,0,1,0),mfrow=c(2,6))
for (i in 4:14) {
  plot(departeL.wgs,lwd=0.8,border=grey(0.7),
       main=colnames(ambro.wgs@data)[i])
  plot(regionsL.wgs,lwd=1.8,add=TRUE)
  plot(ambro.wgs[ambro.wgs@data[,i]==0,],
       pch=as.numeric(as.character(ambro.wgs[ambro.wgs@data[,i]==0,]$SeqMeth)),
       col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
       add=TRUE)
  plot(ambro.wgs[ambro.wgs@data[,i]!=0,],
       pch=as.numeric(as.character(ambro.wgs[ambro.wgs@data[,i]!=0,]$SeqMeth)),
       col=rgb(255,0,0,255,maxColorValue=255),cex=1.5,
       add=TRUE)
}
par(op)

#export to .pdf 18 x 7 inches


##############################################################################/
#code with cropped maps wgs84####
##############################################################################/

op<-par(mar=c(0,0,0,0))
plot(departe.wgs.1,lwd=0.8,border=grey(0.7))
plot(regions.wgs.1,lwd=2.5,add=TRUE)
plot(ambro.wgs[ambro.wgs$RS==0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambro.wgs[ambro.wgs$RS!=0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
par(op)

op<-par(mar=c(0,0,0,0))
plot(departe.wgs.2,lwd=0.8,border=grey(0.7))
plot(regions.wgs.2,lwd=2.5,add=TRUE)
plot(ambro.wgs[ambro.wgs$RS==0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambro.wgs[ambro.wgs$RS!=0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
par(op)

op<-par(mar=c(0,0,0,0))
plot(departe.wgs.3,lwd=0.8,border=grey(0.7))
plot(regions.wgs.3,lwd=2.5,add=TRUE)
plot(ambro.wgs[ambro.wgs$RS==0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambro.wgs[ambro.wgs$RS!=0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
par(op)

op<-par(mar=c(0,0,0,0))
plot(departe.wgs.4,lwd=0.8,border=grey(0.7))
plot(regions.wgs.4,lwd=2.5,add=TRUE)
plot(ambro.wgs[ambro.wgs$RS==0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambro.wgs[ambro.wgs$RS!=0,],
     pch=as.numeric(as.character(ambro.wgs[ambro.wgs$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
par(op)


##############################################################################/
#code for the maps by mutation in Lambert 93####
##############################################################################/

#all mutations on one map
op<-par(mar=c(0,0,0,0))
plot(departe,lwd=0.8,border=grey(0.7))
plot(regions,lwd=1.8,add=TRUE)
plot(ambro[ambro$RS==0,],
     pch=as.numeric(as.character(ambro[ambro$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=0.8,
     add=TRUE)
plot(ambro[ambro$RS!=0,],
     pch=as.numeric(as.character(ambro[ambro$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=0.8,
     add=TRUE)
polygon(x=c(364666.8,590540.8,590540.8,364666.8,364666.8),
        y=c(6612096.8,6612096.8,6431338.8,6431338.8,6612096.8),
        col="transparent",lwd=3,border="darkblue")
polygon(x=c(633281.3,844866.6,844866.6,633281.3,633281.3),
        y=c(6790638.4,6790638.4,6600367.6,6600367.6,6790638.4),
        col="transparent",lwd=3,border="darkblue")
polygon(x=c(747072.2,963827.9,963827.9,747072.2,747072.2),
        y=c(6550171.5,6550171.5,6355840.0,6355840.0,6550171.5),
        col="transparent",lwd=3,border="darkblue")
# polygon(x=c(467049.9,697619.1,697619.1,467049.9,467049.9),
#         y=c(6388990.7,6388990.7,6204268.4,6204268.4,6388990.7),
#         col="transparent",lwd=3,border="darkblue")
polygon(x=c(487049.9,602334.5,602334.5,487049.9,487049.9),
        y=c(6398990.7,6398990.7,6306629.5,6306629.5,6398990.7),
        col="transparent",lwd=3,border="darkblue")
scalebar(c(191260,6060000),300000,"km",division.cex=1)
par(op)
#export to .pdf 7 x 7 inches

#each mutations on separate map
op<-par(mar=c(0,0,1,0),mfrow=c(2,6))
for (i in 4:14) {
  plot(departeLight,lwd=0.8,border=grey(0.7),
       main=colnames(ambro@data)[i])
  plot(regionsLight,lwd=1.8,add=TRUE)
  plot(ambro[ambro@data[,i]==0,],
       pch=as.numeric(as.character(ambro[ambro@data[,i]==0,]$SeqMeth)),
       col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
       add=TRUE)
  plot(ambro[ambro@data[,i]!=0,],
       pch=as.numeric(as.character(ambro[ambro@data[,i]!=0,]$SeqMeth)),
       col=rgb(255,0,0,255,maxColorValue=255),cex=1.5,
       add=TRUE)
}
scalebar(c(191260,6060000),300000,"km",division.cex=1)
par(op)
#export to .pdf 18 x 7 inches

#each position on separate map
op<-par(mar=c(0,0,1,0),mfrow=c(2,3))

#position 197
plot(departeLight,lwd=0.8,border=grey(0.7),
       main="Position 197")
plot(regionsLight,lwd=1.8,add=TRUE)
plot(ambro,pch=as.numeric(as.character(ambro$SeqMeth)),
     col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
     add=TRUE)
plot(ambro[ambro@data[,6]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,6]!=0,]$SeqMeth)),
     col="red",cex=1.5,
     add=TRUE)
plot(ambro[ambro@data[,10]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,10]!=0,]$SeqMeth)),
     col="blue",cex=1.5,
     add=TRUE)
plot(ambro[ambro@data[,11]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,11]!=0,]$SeqMeth)),
     col="purple",cex=1.5,
     add=TRUE)
legend(57000,7160000,legend=c("Gln197","Ser197","Thr197"),cex=1,pt.cex=1.8,
       y.intersp=0.5,x.intersp=0.5,
       pch=19,col=c("red","blue","purple"),bg="transparent",bty="n")
scalebar(c(191260,6060000),300000,"km",division.cex=1)

#position 205
plot(departeLight,lwd=0.8,border=grey(0.7),
     main="Position 205")
plot(regionsLight,lwd=1.8,add=TRUE)
plot(ambro,pch=as.numeric(as.character(ambro$SeqMeth)),
     col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
     add=TRUE)
plot(ambro[ambro@data[,12]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,12]!=0,]$SeqMeth)),
     col="red",cex=1.5,
     add=TRUE)
plot(ambro[ambro@data[,13]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,13]!=0,]$SeqMeth)),
     col="blue",cex=1.5,
     add=TRUE)
plot(ambro[ambro@data[,14]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,13]!=0,]$SeqMeth)),
     col="purple",cex=1.5,
     add=TRUE)
legend(57000,7160000,legend=c("Thr205.1","Thr205.2","Val205"),cex=1,
       pt.cex=1.8,y.intersp=0.5,x.intersp=0.5,
       pch=19,col=c("red","blue","purple"),bg="transparent",bty="n")

#position 376
plot(departeLight,lwd=0.8,border=grey(0.7),
     main="Position 376")
plot(regionsLight,lwd=1.8,add=TRUE)
plot(ambro,pch=as.numeric(as.character(ambro$SeqMeth)),
     col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
     add=TRUE)
plot(ambro[ambro@data[,7]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,7]!=0,]$SeqMeth)),
     col="red",cex=1.5,
     add=TRUE)
legend(57000,7160000,legend=c("Glu376"),cex=1,pt.cex=1.8,
       y.intersp=0.5,x.intersp=0.5,
       pch=19,col=c("red"),bg="transparent",bty="n")

#position 574
plot(departeLight,lwd=0.8,border=grey(0.7),
     main="Position 574")
plot(regionsLight,lwd=1.8,add=TRUE)
plot(ambro,pch=as.numeric(as.character(ambro$SeqMeth)),
     col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
     add=TRUE)
plot(ambro[ambro@data[,4]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,4]!=0,]$SeqMeth)),
     col="red",cex=1.5,
     add=TRUE)
plot(ambro[ambro@data[,9]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,9]!=0,]$SeqMeth)),
     col="blue",cex=1.5,
     add=TRUE)
legend(57000,7160000,legend=c("Arg574","Leu574"),cex=1,pt.cex=1.8,
       y.intersp=0.5,x.intersp=0.5,
       pch=19,col=c("red","blue"),bg="transparent",bty="n")

#position 578
plot(departeLight,lwd=0.8,border=grey(0.7),
     main="Position 578")
plot(regionsLight,lwd=1.8,add=TRUE)
plot(ambro,pch=as.numeric(as.character(ambro$SeqMeth)),
     col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
     add=TRUE)
plot(ambro[ambro@data[,8]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,8]!=0,]$SeqMeth)),
     col="red",cex=1.5,
     add=TRUE)
legend(57000,7160000,legend=c("Ile578"),cex=1,pt.cex=1.8,
       y.intersp=0.5,x.intersp=0.5,
       pch=19,col=c("red"),bg="transparent",bty="n")

#position 654
plot(departeLight,lwd=0.8,border=grey(0.7),
     main="Position 654")
plot(regionsLight,lwd=1.8,add=TRUE)
plot(ambro,pch=as.numeric(as.character(ambro$SeqMeth)),
     col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
     add=TRUE)
plot(ambro[ambro@data[,5]!=0,],
     pch=as.numeric(as.character(ambro[ambro@data[,5]!=0,]$SeqMeth)),
     col="red",cex=1.5,
     add=TRUE)
legend(57000,7160000,legend=c("Arg654"),cex=1,pt.cex=1.8,
       y.intersp=0.5,x.intersp=0.5,
       pch=19,col=c("red"),bg="transparent",bty="n")
scalebar(c(191260,6060000),300000,"km",division.cex=1)

par(op)

#export to .pdf 12 x 7 inches


##############################################################################/
#code with cropped maps Lambert 93####
##############################################################################/

op<-par(mar=c(0,0,0,0))
plot(departe.1,lwd=0.8,border=grey(0.7))
plot(regions.1,lwd=2.5,add=TRUE)
plot(ambro[ambro$RS==0,],
     pch=as.numeric(as.character(ambro[ambro$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambro[ambro$RS!=0,],
     pch=as.numeric(as.character(ambro[ambro$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
scalebar(c(380000,6440000),50000,"km",division.cex=1)
par(op)
#export to .pdf 7 x 5.5 inches

op<-par(mar=c(0,0,0,0))
plot(departe.2,lwd=0.8,border=grey(0.7))
plot(regions.2,lwd=2.5,add=TRUE)
plot(ambro[ambro$RS==0,],
     pch=as.numeric(as.character(ambro[ambro$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambro[ambro$RS!=0,],
     pch=as.numeric(as.character(ambro[ambro$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
scalebar(c(780000,6610000),50000,"km",division.cex=1)
par(op)
#export to .pdf 7 x 5.5 inches

op<-par(mar=c(0,0,0,0))
plot(departe.3,lwd=0.8,border=grey(0.7))
plot(regions.3,lwd=2.5,add=TRUE)
plot(ambro[ambro$RS==0,],
     pch=as.numeric(as.character(ambro[ambro$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambro[ambro$RS!=0,],
     pch=as.numeric(as.character(ambro[ambro$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
scalebar(c(762000,6530000),50000,"km",division.cex=1)
par(op)
#export to .pdf 7 x 5.5 inches

# op<-par(mar=c(0,0,0,0))
# plot(departe.4,lwd=0.8,border=grey(0.7))
# plot(regions.4,lwd=2.5,add=TRUE)
# plot(ambro[ambro$RS==0,],
#      pch=as.numeric(as.character(ambro[ambro$RS==0,]$SeqMeth)),
#      col=rgb(0,0,0,150,maxColorValue=255),cex=2,
#      add=TRUE)
# plot(ambro[ambro$RS!=0,],
#      pch=as.numeric(as.character(ambro[ambro$RS!=0,]$SeqMeth)),
#      col=rgb(255,0,0,150,maxColorValue=255),cex=2,
#      add=TRUE)
# scalebar(c(630000,6370000),50000,"km",division.cex=1)
# par(op)
# #export to .pdf 7 x 5.5 inches

op<-par(mar=c(0,0,0,0))
plot(departe.5,lwd=0.8,border=grey(0.7))
plot(regions.5,lwd=2.5,add=TRUE)
plot(ambro[ambro$RS==0,],
     pch=as.numeric(as.character(ambro[ambro$RS==0,]$SeqMeth)),
     col=rgb(0,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
plot(ambro[ambro$RS!=0,],
     pch=as.numeric(as.character(ambro[ambro$RS!=0,]$SeqMeth)),
     col=rgb(255,0,0,150,maxColorValue=255),cex=2,
     add=TRUE)
scalebar(c(493000,6388000),25000,"km",division.cex=1)
par(op)
#export to .pdf 7 x 5.5 inches


##############################################################################/
#maps screening vs monitoring####
##############################################################################/

op<-par(mar=c(0,0,1,0),mfrow=c(1,2))
#original screening
plot(departe,lwd=0.8,border=grey(0.7),
     main="Screening (N=45)")
plot(regions,lwd=1.8,add=TRUE)
plot(ambro[ambro$RS==0 & ambro$SampRound=="Screening",],
     pch=1,
     col=rgb(50,100,0,150,maxColorValue=255),cex=1.1,
     add=TRUE)
plot(ambro[ambro$RS!=0 & ambro$SampRound=="Screening",],
     pch=1,
     col=rgb(255,0,0,250,maxColorValue=255),cex=1.1,
     add=TRUE)
#NGS monitoring
plot(departe,lwd=0.8,border=grey(0.7),
     main="Monitoring (N=210)")
plot(regions,lwd=1.8,add=TRUE)
plot(ambro[ambro$RS==0 & ambro$SampRound!="Screening",],
     pch=3,
     col=rgb(50,100,0,150,maxColorValue=255),cex=0.4,
     add=TRUE)
plot(ambro[ambro$RS!=0 & ambro$SampRound!="Screening",],
     pch=3,
     col=rgb(255,0,0,250,maxColorValue=255),cex=0.4,
     add=TRUE)
scalebar(c(191260,6060000),300000,"km",division.cex=1)
par(op)
#export to .pdf 12 x 5 inches


##############################################################################/
#maps phenotype####
##############################################################################/

#2 maps with only independent information
op<-par(mar=c(0,0,1,0),mfrow=c(1,2))
#Imazamox phenotype
plot(departe,lwd=0.8,border=grey(0.7),
     main="Imazamox (N=39)")
plot(regions,lwd=1.8,add=TRUE)
plot(ambro[ambro$Imz=="S",],
     pch=19,
     col=rgb(50,100,0,150,maxColorValue=255),cex=1,
     add=TRUE)
plot(ambro[ambro$Imz=="R",],
     pch=19,
     col=rgb(255,0,0,250,maxColorValue=255),cex=1,
     add=TRUE)
#Tribenuron phenotype
plot(departe,lwd=0.8,border=grey(0.7),
     main="Tribenuron (N=36)")
plot(regions,lwd=1.8,add=TRUE)
plot(ambro[ambro$Tbn=="S",],
     pch=19,
     col=rgb(50,100,0,150,maxColorValue=255),cex=1,
     add=TRUE)
plot(ambro[ambro$Tbn=="R",],
     pch=19,
     col=rgb(255,0,0,250,maxColorValue=255),cex=1,
     add=TRUE)
scalebar(c(191260,6060000),300000,"km",division.cex=1)
par(op)
#export to .pdf 12 x 5 inches

#2 maps with both information
op<-par(mar=c(0,0,1,0),mfrow=c(1,2))
#Imazamox phenotype
plot(departe,lwd=0.8,border=grey(0.7),
     main="Imazamox (N=39)")
plot(regions,lwd=1.8,add=TRUE)
plot(ambro[ambro$Imz=="S",],
     pch=19,
     col=rgb(50,100,0,150,maxColorValue=255),cex=1,
     add=TRUE)
plot(ambro[ambro$Imz=="R" & ambro$Tbn=="R",],
     pch=19,
     col=rgb(255,0,0,250,maxColorValue=255),cex=1,
     add=TRUE)
plot(ambro[ambro$Imz=="R" & ambro$Tbn=="R",],
     pch=19,
     col=rgb(255,0,0,250,maxColorValue=255),cex=1,
     add=TRUE)
#Tribenuron phenotype
plot(departe,lwd=0.8,border=grey(0.7),
     main="Tribenuron (N=36)")
plot(regions,lwd=1.8,add=TRUE)
plot(ambro[ambro$Tbn=="S",],
     pch=19,
     col=rgb(50,100,0,150,maxColorValue=255),cex=1,
     add=TRUE)
plot(ambro[ambro$Tbn=="R",],
     pch=19,
     col=rgb(255,0,0,250,maxColorValue=255),cex=1,
     add=TRUE)
scalebar(c(191260,6060000),300000,"km",division.cex=1)
par(op)
#export to .pdf 12 x 5 inches


##############################################################################/
#END
##############################################################################/