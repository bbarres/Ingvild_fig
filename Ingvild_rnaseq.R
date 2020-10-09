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
countdata<-countdata[,order(colnames(countdata))]
head(countdata)
coldata=read.table("data/coldata-120H-148H-173H-UT-2019x2020.txt",header=T)
head(coldata)
coldata<-coldata[order(coldata$Run),]
colnames(countdata)=rownames(coldata)

coldata$Status<-as.factor(coldata$Status)
coldata$Annee<-as.factor(coldata$Annee)
coldata$Population<-as.factor(coldata$Population)

coldata19<-coldata[coldata$Annee=="2019",]
coldata20<-coldata[coldata$Annee=="2020",]
count2019<-countdata[,coldata$Annee=="2019"]
count2020<-countdata[,coldata$Annee=="2020"]

dds=DESeqDataSetFromMatrix(countdata,coldata, ~Status)
dds1=DESeqDataSetFromMatrix(countdata,coldata, ~Annee + Population + Status)
dds2=DESeqDataSetFromMatrix(countdata,coldata, ~Annee*Population + Status)
dds19=DESeqDataSetFromMatrix(count2019,coldata19, ~Population + Status + 
                               Population:Status)
dds20=DESeqDataSetFromMatrix(count2020,coldata20, ~Population + Status + 
                               Population:Status)

nrow(dds1)
head(counts(dds1))
colData(dds1)

#ordering factor levels (default: alphabetical order)
dds$Status=relevel(dds$Status,ref="s")
dds1$Status=relevel(dds1$Status,ref="s")
dds2$Status=relevel(dds2$Status,ref="s")
dds19$Status=relevel(dds19$Status,ref="s")
dds20$Status=relevel(dds20$Status,ref="s")


##############################################################################/
#DESeq analysis####
##############################################################################/

#analysis with only status as a variable####
dds<-DESeq(dds)
res<-results(dds)
head(res[order(res$padj),])
summary(res)

resultsNames(dds)
res<-results(dds,contrast=c("Status","r","s"))
summary(res)
head(res)

resOrdered=res[order(res$padj),]
head(resOrdered)
write.table(resOrdered,"output/BWA-Res-Deseq-120H-148H-173H-19x20.csv",
            sep=";",dec=",")

#MA plot and dispersion graphic
plotMA(res,ylim=c(-2, 2) )
plotDispEsts(dds,ylim=c(1e-6, 1e1) )


#analysis taking into account year and population####
dds1<-DESeq(dds1)
res1<-results(dds1)
head(res1[order(res1$padj),])
summary(res1)

resultsNames(dds1)
res1<-results(dds1,contrast=c("Status","r","s"))
summary(res1)
head(res1)
#there is more than a quarter of the genes that are up regulated and 
#another quarter that are down regulated between years
summary(results(dds1,contrast=c("Annee","2020","2019")))
#there is more than a 15% of the genes that are up regulated and 
#more than 13% that are down regulated between pop 120H and 173H
summary(results(dds1,contrast=c("Population","173H","120H")))

resOrdered1=res1[order(res1$padj),]
head(resOrdered1)
write.table(resOrdered1,
            "output/BWA-Res-Deseq-120H-148H-173H-19x20_Year+Pop.csv",
            sep=";",dec=",")

#MA plot and dispersion graphic
plotMA(res1,ylim=c(-2, 2) )
plotDispEsts(dds1,ylim=c(1e-6, 1e1) )


#2019analysis taking into account population and interac status:population####
dds19<-DESeq(dds19)
res19<-results(dds19)
head(res19[order(res19$padj),])
summary(res19)
resultsNames(dds19)
#results for population 120H
res19_120H<-results(dds19,contrast=c("Status","r","s"))
summary(res19_120H)
head(res19_120H[order(res19_120H$padj),])
plot(sort(log(res19_120H$padj))[1:500])
abline(h=c(log(0.05),log(0.01),log(0.001)),
       col=c("green","blue","red"))
#results for population 148H
res19_148H<-results(dds19,contrast=list(c("Status_r_vs_s",
                                          "Population148H.Statusr")))
summary(res19_148H)
head(res19_148H[order(res19_148H$padj),])
plot(sort(log(res19_148H$padj))[1:500])
abline(h=c(log(0.05),log(0.01),log(0.001)),
       col=c("green","blue","red"))
#results for population 173H
res19_173H<-results(dds19,contrast=list(c("Status_r_vs_s",
                                          "Population173H.Statusr")))
summary(res19_173H)
head(res19_173H[order(res19_173H$padj),])
plot(sort(log(res19_173H$padj))[1:500])
abline(h=c(log(0.05),log(0.01),log(0.001)),
       col=c("green","blue","red"))
#MA plot and dispersion graphic
plotMA(res19_120H,ylim=c(-2, 2))
plotMA(res19_148H,ylim=c(-2, 2))
plotMA(res19_173H,ylim=c(-2, 2))
plotDispEsts(dds19,ylim=c(1e-6, 1e1))


#2020analysis taking into account population and interac status:population####
dds20<-DESeq(dds20)
res20<-results(dds20)
head(res20[order(res20$padj),])
summary(res20)
#results for population 120H
resultsNames(dds20)
res20_120H<-results(dds20,contrast=c("Status","r","s"))
summary(res20_120H)
head(res20_120H)
plot(sort(log(res20_120H$padj))[1:500])
abline(h=c(log(0.05),log(0.01),log(0.001)),
       col=c("green","blue","red"))
#results for population 148H
res20_148H<-results(dds20,contrast=list(c("Status_r_vs_s",
                                          "Population148H.Statusr")))
summary(res20_148H)
head(res20_148H)
plot(sort(log(res20_148H$padj))[1:500])
abline(h=c(log(0.05),log(0.01),log(0.001)),
       col=c("green","blue","red"))
#results for population 173H
res20_173H<-results(dds20,contrast=list(c("Status_r_vs_s",
                                          "Population173H.Statusr")))
summary(res20_173H)
head(res20_173H)
plot(sort(log(res20_173H$padj))[1:500])
abline(h=c(log(0.05),log(0.01),log(0.001)),
       col=c("green","blue","red"))
#MA plot and dispersion graphic
plotMA(res20_148H,ylim=c(-2, 2) )
plotDispEsts(dds20,ylim=c(1e-6, 1e1) )

#exporting all the results
combiRez<-cbind(res19_120H,res19_148H,res19_173H,
                res20_120H,res20_148H,res20_173H)
write.table(combiRez,file="output/combirez.csv",sep=";")


##############################################################################/
#END
##############################################################################/