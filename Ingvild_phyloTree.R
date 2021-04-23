##############################################################################/
##############################################################################/
#Code for phylogenetic tree
##############################################################################/
##############################################################################/

#loading the necessary packages
library(seqinr)
library(ape)
library(phangorn)
library(RColorBrewer)

#loading the alignment
AmbAlign<-read.alignment(file="data/Alignt-haplotypes-R-Sanger-NGS-v7.fas",
                         format="fasta")
AmbAlign<-as.DNAbin(AmbAlign)
AmbAlign<-phyDat(AmbAlign,"DNA")


##############################################################################/
#Computation of the distance matrix####
##############################################################################/

AmbDist<-dist.dna(as.DNAbin(AmbAlign))


AmbTree<-upgma(AmbDist)
AmbTreePars<-optim.parsimony(AmbTree,AmbAlign)


##############################################################################/
#building the tree####
##############################################################################/

op<-par(mfrow=c(2,2))
plot(AmbTree,"cladogram")
plot(AmbTree,"fan")
tiplabels(pch=c(19,19,21,21,22,23),bg=c(1,2,3),cex=2)

bottTree<-bootstrap.phyDat(AmbAlign,pratchet,bs=100)
bottTree2<-bootstrap.phyDat(AmbAlign,
                            FUN=function(x)upgma(dist.dna(as.DNAbin(x))),bs=100)
plotBS(AmbTree,bottTree,"phylogram")
plotBS(AmbTree,bottTree2,"cladogram")
tiplabels(pch=c(21,21,22,23,20,24),bg=c(1,2,3),cex=3,offset=0.002)
plotBS(AmbTree,bottTree2,"fan",show.tip.label=FALSE,edge.width=5)
tiplabels(pch=c(21,21,22,23,25,24),bg=c(1,2,3),cex=3,offset=0.001)
par(op)
#export to .pdf 8 x 8 inches

AmbTree$tip.label

##############################################################################/
#END
##############################################################################/