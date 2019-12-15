if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylumi")
suppressPackageStartupMessages(require('methylumi'))
suppressPackageStartupMessages(require('TCGAMethylation450k'))
suppressPackageStartupMessages(require('FDb.InfiniumMethylation.hg19'))
#################################################################################
######EJEMPLO##########
library(MethylMix)
data("METcancer")
data("METnormal")
data("GEcancer")
METnormal[,1:4]
head(METcancer)
head(GEcancer)
MethylMixResults<-MethylMix(METcancer, GEcancer, METnormal)
MethylMixResults$MethylationDrivers
plots<-MethylMix_PlotModel("MGMT", MethylMixResults, METcancer)
plots$CorrelationPlot
plots$MixtureModelPlot
plots$MixtureModelPlot
MethylMixResults$Classifications
MethylMixResults$MethylationStates
############################################################################
library(readr)
meth_brca <- readLines('randomized_meth_brca.txt')
meth_BRCA <- read.table(textConnection(meth_brca[-2]), header = TRUE, sep = "\t")
colnames(meth_BRCA)[1:4] <- c('probe', 'gene_symbol', 'chr', 'start')
dim(meth_BRCA)
############################################################################
install.packages("stringr")
library(stringr)
Samples<-colnames(meth_BRCA[,5:225])
Samples<-gsub(Samples, pattern = ".", replacement = "-", fixed = TRUE)
idx2<-sapply(strsplit(Samples, split = '-'), function(x) return(x[4]))
table(idx2)
posiciones11b<-which(idx2=="11B")
posiciones11a<-which(idx2=="11A")
Samples1<-Samples
Samples1<-Samples1[posiciones11a]
Samples2<-Samples[posiciones11b]
Samples_Normal<-c(Samples1, Samples2)
idx3<-match(Samples_Normal, Samples)
Samples_Tumor<-Samples[-idx3]
##########################################################################
colnames(meth_BRCA)<-gsub(colnames(meth_BRCA), pattern = ".", replacement = "-", fixed = TRUE)
idx2<-match( Samples_Tumor,colnames(meth_BRCA))
MetCancer<-meth_BRCA[,idx2]
idx2<-match( Samples_Normal,colnames(meth_BRCA))
METnormal<-meth_BRCA[,idx2]
load("C:/Users/Lucio/Desktop/Bioinformatica/TFM/Bases de datos/Breast_rna.RData")
colnames(Breast_rna)<-gsub(colnames(Breast_rna), pattern = ".", replacement = "-", fixed = TRUE)
colnames(MetCancer)<-substr(colnames(MetCancer),1,12)
idx2<-match(colnames(METcancer[,2:196]),colnames(Breast_rna))
Breast_rna<-Breast_rna[,idx2]
#########################################################################
Samples<-colnames(meth_BRCA[,5:225])
Samples<-substr(Samples,1,12)
idx<-colnames(Breast_rna)
idx2<-match(Samples, idx)
Breast_rna<-Breast_rna[,idx2]
GEcancer<-Breast_rna
METcancer<-METCancer
###########################################
#Nombre de filas por simbolo del gen
METcancer<-cbind(meth_BRCA[,1:4],METcancer)
METnormal<-cbind(meth_BRCA[,1:4],METnormal)
idx<-METcancer$gene_symbol
idx<-which(is.na(METcancer$gene_symbol))
METcancer<-METcancer[-idx,]
idx<-which(is.na(METnormal$gene_symbol))
#####################################
METcancer<-METcancer[,c(2,5:199)]
METnormal<-as.data.frame(METnormal[,c(2,5:26)])
GEcancer<-as.matrix(GEcancer)
METnormal<-as.matrix(METnormal)
METcancer<-as.matrix(METcancer)
############################################################################
METnormal[,1]
rownames(METnormal)<-sapply(strsplit(METnormal[,1], split = ';'), function(x) return(x[1]))
rownames(METcancer)<-sapply(strsplit(METcancer[,1], split = ';'), function(x) return(x[1]))
rownames<-rownames(GEcancer)
rownames(GEcancer)<-gsub(rownames(GEcancer), pattern = "|", replacement = ";", fixed = TRUE)
GEcancer<-as.matrix(GEcancer)
rownames(GEcancer)<-sapply(strsplit(rownames(GEcancer), split = ';'), function(x) return(x[1]))
idx<-match(row.names(METnormal),row.names(GEcancer))
GEcancer<-as.data.frame(GEcancer)
rownames(GEcancer)<-rownames
GEcancer<-GEcancer[idx,]
GEcancer<-as.matrix(GEcancer)
METcancer<-METcancer[,2:195]
METnormal<-METnormal[,2:23]
METnormal<-as.matrix(METnormal)
METcancer<-as.matrix(METcancer)
rownames(GEcancer)<-gsub(rownames(GEcancer), pattern = "|", replacement = ";", fixed = TRUE)
rownames(GEcancer)<-sapply(strsplit(rownames(GEcancer), split = ';'), function(x) return(x[1]))
GEcancer<-as.data.frame(GEcancer)
GEcancer<-as.matrix(GEcancer)
METnormal<-as.data.frame(METnormal)



genes<-rownames(METnormal)
METnormal<-apply(METnormal,2,as.numeric)
rownames(METnormal)<-genes

genes<-rownames(METcancer)
METnormal<-apply(METcancer,2,as.numeric)
rownames(METcancer)<-genes


METnormal<-as.matrix(METnormal)
METcancer<-as.data.frame(METcancer)
METcancer<-as.matrix(METcancer)
colnames(METnormal)<-substr(colnames(METnormal),1,12)
levels(METnormal)
idx<-which(is.na(METnormal))

row.has.na <- apply(METcancer, 1, function(x){any(is.na(x))})
sum(row.has.na)
METcancer <- METcancer[!row.has.na,]

row.has.na <- apply(METnormal, 1, function(x){any(is.na(x))})
sum(row.has.na)
METnormal <- METnormal[!row.has.na,]
idx<-match(rownames(METcancer), rownames(METnormal))
METnormal<-METnormal[idx,]
idx<-match(rownames(METcancer), rownames(GEcancer))
GEcancer<-GEcancer[idx,]
####################
######EJEMPLO##########
library(MethylMix)
saveRDS(METnormal, file = "METnormal.rds")
saveRDS(METcancer,file = "METcancer.rds")
saveRDS(GEcancer, file = "GEcancer.rds")
METcancer<-as.matrix(METcancer)
METnormal<-METnormal[1:5,]
GEcancer<-GEcancer[1:5,]
METnormal[,1:4]
head(METcancer)
head(GEcancer)
match(colnames(METcancer), colnames(GEcancer))
METcancer<-apply(METcancer,2, as.numeric)
class(METcancer)

MethylMixResults<-MethylMix(METcancer[1:997,], GEcancer[1:997,], METnormal[1:997,], filter = TRUE, NoNormalMode = FALSE, listOfGenes = rownames(METcancer), OutputRoot = "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2")

Drivers<-MethylMixResults$MethylationDrivers
Drivers
plots<-MethylMix_PlotModel("PAX8", MethylMixResults, METcancer)
plots$CorrelationPlot
plots$MixtureModelPlot
plots$MixtureModelPlot
MethylMixResults$Classifications
MethylMixResults$MethylationStates
# Plot MGMT also with its normal methylation variation
plots <- MethylMix_PlotModel("PAX8", MethylMixResults, METcancer, METnormal = METnormal)
plots$MixtureModelPlot
# Also plot the inverse correlation with gene expression (creates two
# separate plots)
plots <- MethylMix_PlotModel("TAP1", MethylMixResults, METcancer, GEcancer, 
                             METnormal)
plots$MixtureModelPlot
plots$CorrelationPlot
# Plot all functional and differential genes
gene<-Drivers
for (gene in MethylMixResults$MethylationDrivers) {
  MethylMix_PlotModel(gene, MethylMixResults, METcancer, METnormal = METnormal)
}

############################################################################
write.xlsx2(Drivers, file = "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/Drivers_breast.xlsx")
library(xlsx)      
############################################################################
library(ConsensusClusterPlus)
MethylMixResults<- MethylMix(METcancer, GEcancer, METnormal)

ComplexHeatmap::Heatmap(Metilation)
Metilation<-cbind(as.data.frame(METcancer), as.data.frame(METnormal))

Samplestype<-c(1:216)
Samplestype[1:194]<-"Tumor"
Samplestype[195:216]<-"Normal"
library(gplots)
DMvalues<- MethylMixResults$MethylationStates
cons_cluster<- ConsensusClusterPlus(d=DMvalues, maxK=10, reps=1000, pItem=0.8, distance='euclidean', clusterAlg="km")
heatmap.2(as.matrix(DMvalues),scale='row',col=bluered(149),trace='none',
          main = "DM Values Centered CpG sites",margins = c(9, 22),density.inf="none",symkey=TRUE)
BiocManager::install("pvcluster")

pvclust(as.matrix(Metilation), method.hclust="average",
        method.dist="correlation", use.cor="pairwise.complete.obs",
        nboot=1000, parallel=FALSE, r=seq(.5,1.4,by=.1),
        store=FALSE, weight=FALSE, iseed=NULL, quiet=FALSE)

library(pvclust)
Boxplot<-as.data.frame(t(DMvalues))
boxplot(Boxplot, las =2, col = rainbow(39), main = "Predictive Genes DM Values")


idx1<-match( Drivers, rownames(Metilation))
Metilation<-Metilation[idx1,]
BetaValues<-Metilation
boxplot(t(BetaValues), las = 2, col = rainbow(39), main = "Predictive Genes Beta Values")
