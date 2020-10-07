##############################################################################/
##############################################################################/
#RNAseq data analyses
##############################################################################/
##############################################################################/

#loading the necessary packages
library(DESeq2)
library(RColorBrewer)

#loading and preparing the dataset
logstattout<-read.csv("data/log2.csv",sep=";",header=TRUE,
                      stringsAsFactors=FALSE)
#creating the "popula" variable to tag each individual
popula<-as.factor(
  sapply(strsplit(as.character(logstattout$sorteddatatoutpheno),
                  split="-"),"[[",1))
logstattout<-cbind(logstattout,popula)
#defining a vector of 5 colors, one for each population
colopop<-brewer.pal(5,"Dark2")


##############################################################################/
#independant files, ordered by population then by ascending log2 
##############################################################################/






##############################################################################/
#END
##############################################################################/