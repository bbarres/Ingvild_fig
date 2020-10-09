library("DESeq2")

#Count matrix input

countdata=read.table("data/BWA-counts-120H-148H-173H-2019x2020.txt", header=T)
head(countdata)
summary(countdata)
coldata=read.table("coldata-UT-120H-148H-173H-19x20.txt", header=T)
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

#DATA VisUALISATION
###################

#variance stabilizing transformation
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

#heatmap


library("RColorBrewer")

#distance euclidienne
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Population,rld$Status,rld$Annee,sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#PCA plot
library(ggplot2)

a<-plotPCA(rld, intgroup = c("Status","Annee"))
a+scale_color_manual(values=c("red","orange","green","darkgreen"))
a+geom_label(aes(label=rld$Population))+scale_color_manual(values=c("red", "orange","green","darkgreen"))


#Running the differential expression pipeline
#############################################

dds <- DESeq(dds)

res <- results(dds)
head(res)
summary(res)

resultsNames(dds)
res<-results(dds, contrast=c("Status","r","s"))
summary(res)
head(res)

resOrdered=res[order(res$padj),]
head(resOrdered)
write.table(resOrdered, "BWA-Res-Deseq-120H-148H-173H-19x20.csv", sep=";",dec=",")


#MA plot and dispersion graphic
plotMA(res, ylim = c(-2, 2) )

plotDispEsts( dds, ylim = c(1e-6, 1e1) )