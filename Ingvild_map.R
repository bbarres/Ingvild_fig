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
library(RColorBrewer)

#load geographical data
load("data/regions.RData")
load("data/departe.RData")
#changing the projection of the map
departe.wgs <- spTransform(departe,
                           CRS("+proj=longlat +datum=WGS84"))
#changing the projection of the map
regions.wgs <- spTransform(regions,
                           CRS("+proj=longlat +datum=WGS84"))

#load the resistance results for the 2020 campaign
databruteTOT<-read.delim(
  "data/data_DC_classes_2020.txt",
  header = TRUE,
  sep = "\t",
  colClasses = c("character", "character", "character",
                 "character", "factor")
)

#defining a vector of 5 colors, one for each population
colopop<-brewer.pal(5,"Dark2")


##############################################################################/
#code for the maps by mutation 
##############################################################################/



op<-par(mar=c(0,0,0,0))
plot(departe.wgs,lwd=0.8,border=grey(0.7))
plot(regions.wgs,lwd=2.5,add=TRUE)

par(op)



##############################################################################/
#END
##############################################################################/