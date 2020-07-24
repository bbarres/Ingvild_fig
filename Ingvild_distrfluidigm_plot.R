##############################################################################/
##############################################################################/
#Plotting distribution of 
##############################################################################/
##############################################################################/

#loading the necessary packages
library(dplyr)
library(ggplot2)
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

for (i in 3:(dim(logstattout)[2]-1)) {
  
  datadd<-logstattout[,c(1,2,i,90)]
  #replacing missing data by "0"
  datadd[is.na(datadd[,3]),3]<-0
  #here you determine the order of the data
  datadd<-arrange(datadd,popula,get(colnames(logstattout)[i]))
  datadd$ordina<-seq.int(nrow(datadd))
  p<-ggplot(data=datadd,aes(x=reorder(sorteddatatoutpheno,ordina),
                            y=get(colnames(logstattout)[i]),
                            fill=sorteddatatoutpheno_1),
            position=position_dodge()) + 
    geom_bar(stat="identity") +
    theme_minimal()
  
  p + theme(axis.text.x = element_text(angle=90,size=12,hjust=0,vjust=0.5,
                                       colour=colopop[datadd$popula])) + 
    scale_fill_manual(values=c("#FF0000","#33CC33")) + 
    theme(legend.position="bottom") + 
    xlab("individus") + ylab("log2(ddct)") + 
    labs(fill = "Phenotype") + 
    ggtitle(colnames(logstattout)[i],
            subtitle="pval : ajouter ligne tableau") + #here something to adjust
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(plot.subtitle = element_text(hjust = 0.5)) + 
    ylim(-10,10)
  
  ggsave(file=paste("./output/log(",colnames(logstattout)[i],").png",sep=""),
         width=11,height=8)
  
}


##############################################################################/
#one pdf file, ordered by population then by ascending log2 
##############################################################################/



##############################################################################/
#independant files, ordered by ascending log2 
##############################################################################/

for (i in 3:(dim(logstattout)[2]-1)) {
  
  datadd<-logstattout[,c(1,2,i,90)]
  #replacing missing data by "0"
  datadd[is.na(datadd[,3]),3]<-0
  #here you determine the order of the data
  datadd<-arrange(datadd,get(colnames(logstattout)[i]))
  datadd$ordina<-seq.int(nrow(datadd))
  p<-ggplot(data=datadd,aes(x=reorder(sorteddatatoutpheno,ordina),
                            y=get(colnames(logstattout)[i]),
                            fill=sorteddatatoutpheno_1),
            position=position_dodge()) + 
    geom_bar(stat="identity") +
    theme_minimal()
  
  p + theme(axis.text.x = element_text(angle=90,size=12,hjust=0,vjust=0.5,
                                       colour=colopop[datadd$popula])) + 
    scale_fill_manual(values=c("#FF0000","#33CC33")) + 
    theme(legend.position="bottom") + 
    xlab("individus") + ylab("log2(ddct)") + 
    labs(fill = "Phenotype") + 
    ggtitle(colnames(logstattout)[i],
            subtitle="pval : ajouter ligne tableau") + #here something to adjust
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(plot.subtitle = element_text(hjust = 0.5)) + 
    ylim(-10,10)
  
  ggsave(file=paste("./output/log(",colnames(logstattout)[i],").png",sep=""),
         width=11,height=8)
  
}


##############################################################################/
#END
##############################################################################/