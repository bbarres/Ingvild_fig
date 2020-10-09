##############################################################################/
##############################################################################/
#RNAseq data analyses
##############################################################################/
##############################################################################/

#loading the necessary packages
library("DESeq2")
library("RColorBrewer")
library("pheatmap")

#loading and preparing the dataset
countdata=read.table("data/BWA-counts-120H-148H-173H-2019x2020.txt",header=T)
head(countdata)
summary(countdata)
coldata=read.table("data/coldata-120H-148H-173H-UT-2019x2020.txt",header=T)
head(coldata)
colnames(countdata)=rownames(coldata)

coldata$Status<-as.factor(coldata$Status)
coldata$Annee<-as.factor(coldata$Annee)
coldata$Population<-as.factor(coldata$Population)

dds=DESeqDataSetFromMatrix(countdata, coldata, ~Status)

nrow(dds)
head(counts(dds))
colData(dds)

#ordering factor levels (default: alphabetical order)
dds$Status=relevel(dds$Status, ref = "s")
ddsTxi$Treat=relevel(ddsTxi$Treat, ref = "AT")


##############################################################################/
#independant files, ordered by population then by ascending log2 
##############################################################################/






##############################################################################/
#END
##############################################################################/