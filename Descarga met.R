if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylumi")
suppressPackageStartupMessages(require('methylumi'))
suppressPackageStartupMessages(require('TCGAMethylation450k'))
suppressPackageStartupMessages(require('FDb.InfiniumMethylation.hg19'))
#################################################################################

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
METcancer<-MetCancer
idx2<-match(colnames(METcancer[,1:195]),colnames(Breast_rna))
Breast_rna<-Breast_rna[,idx2]
#########################################################################
Samples<-colnames(meth_BRCA[,5:225])
Samples<-substr(Samples,1,12)
idx<-colnames(Breast_rna)
idx2<-match(Samples, idx)
#Breast_rna<-Breast_rna[,idx2]
GEcancer<-Breast_rna
#METcancer<-METCancer
###########################################
#Nombre de filas por simbolo del gen
METcancer<-cbind(meth_BRCA[,1:4],METcancer)
METnormal<-cbind(meth_BRCA[,1:4],METnormal)
idx<-METcancer$gene_symbol
idx<-which(is.na(METcancer$gene_symbol))
METcancer<-METcancer[-idx,]
METnormal<-METnormal[-idx,]
#####################################
#METcancer<-METcancer[,c(2,5:199)]
#METnormal<-as.data.frame(METnormal[,c(2,5:26)])
#GEcancer<-as.matrix(GEcancer)
#METnormal<-as.matrix(METnormal)
##METcancer<-as.matrix(METcancer)
############################################################################
#rownames<-METcancer
METcancer<-as.data.frame(METcancer)
METnormal<-as.data.frame(METnormal)
idx<-which(duplicated(METcancer$gene_symbol))
METcancer<-METcancer[-idx,]
METnormal<-METnormal[-idx,]
rownames(METcancer)<-METcancer$gene_symbol
rownames(METnormal)<-METcancer$gene_symbol
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
METcancer<-METcancer[,5:199]
METnormal<-METnormal[,5:30]
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
rownames(METnormal)<-genes

METnormal<-as.matrix(METnormal)
METcancer<-as.data.frame(METcancer)
METcancer<-as.matrix(METcancer)
#colnames(METnormal)<-substr(colnames(METnormal),1,12)
levels(METnormal)
#idx<-which(is.na(METnormal))

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
saveRDS(METnormal, file = "METnormal_BRCA.rds")
saveRDS(METcancer,file = "METcancer_BRCA.rds")
saveRDS(GEcancer, file = "GEcancer_BRCA.rds")
#METcancer<-as.matrix(METcancer)
#METnormal<-METnormal[1:5,]
#GEcancer<-GEcancer[1:5,]
#METnormal[,1:4]
#head(METcancer)
#head(GEcancer)
#match(colnames(METcancer), colnames(GEcancer))
METcancer<-apply(METcancer,2, as.numeric)
class(METcancer)

MethylMixResults<-MethylMix(METcancer, GEcancer, METnormal, filter = TRUE, listOfGenes = rownames(METcancer), OutputRoot = "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2")

Drivers<-MethylMixResults$MethylationDrivers
Drivers
plots<-MethylMix_PlotModel("ZNF667", MethylMixResults, METcancer, GEcancer)
plots$CorrelationPlot
plots$MixtureModelPlot
plots$MixtureModelPlot
MethylMixResults$Classifications
MethylMixResults$MethylationStates
# Plot MGMT also with its normal methylation variation
plots <- MethylMix_PlotModel("ZNF274", MethylMixResults, METcancer, METnormal = METnormal, GEcancer)
plots$CorrelationPlot
plots$MixtureModelPlot
# Also plot the inverse correlation with gene expression (creates two
# separate plots)
par(mfrow=c(2,2))
plots <- MethylMix_PlotModel("NUDT12", MethylMixResults, METcancer, GEcancer, 
                             METnormal)
plots$MixtureModelPlot
plots$CorrelationPlot
# Plot all functional and differential genes
gene<-Drivers
for (i in gene) {
 plots<-MethylMix_PlotModel(i, MethylMixResults, METcancer, METnormal = METnormal)
 plots$MixtureModelPlot 
}
plots$MixtureModelPlot
gene
############################################################################
write.xlsx2(Drivers, file = "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/Drivers_Breast.xlsx")
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
heatmap.2(as.matrix(METcancer),scale='row',col=bluered(149),trace='none',
          main = "  Beta Values",margins = c(9, 22),density.inf="none",symkey=TRUE)
heatmap.2(as.matrix(DMvalues),scale='row',col=bluered(149),trace='none',
          main = "  DM Values",margins = c(9, 22),density.inf="none",symkey=TRUE)
par(mfrow=c(1,2))
heatmap.2(as.matrix(Metilation),scale='row',col=bluered(149),trace='none',
          main = "  Beta Values",margins = c(9, 22),density.inf="none",symkey=TRUE)
BiocManager::install("pvcluster")

pvclust(as.matrix(Metilation), method.hclust="average",
        method.dist="correlation", use.cor="pairwise.complete.obs",
        nboot=1000, parallel=FALSE, r=seq(.5,1.4,by=.1),
        store=FALSE, weight=FALSE, iseed=NULL, quiet=FALSE)

library(pvclust)
Boxplot<-as.data.frame(t(DMvalues))
boxplot(Boxplot[,1:40], las =2, col = rainbow(39), main = "Predictive Genes DM Values", cex.names=0.2)


idx1<-match( Drivers, rownames(Metilation))
Metilation<-Metilation[idx1,]
BetaValues<-Metilation
BetaValues<-as.data.frame(t(BetaValues))
boxplot(BetaValues[,1:40], las = 2, col = rainbow(40), main = "Predictive Genes Beta Values")

par(mfrow=c(3,3))
genes<-Drivers
for (i in genes){
boxplot(METcancer[i,], METnormal[i,], names = c("Tumor", "Normal"), col = rainbow(2), main = i)
}

#############################################################################################
#Análisis expresión ARN y CNVmut
library(readr)
cna <- read_delim("cna.txt", "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
rownames<-cna$SAMPLE_ID
cna<-as.data.frame(t(cna))
colnames(cna)<-rownames
cna<-cna[-1,]
cna<-cna[-1,]
colnames(cna)<-substr(colnames(cna),1,12)
idx<-match(colnames(GEcancer),colnames(cna))
cna<-cna[,idx]
cna<-as.numeric(cna)
cna<-apply(cna, 2, as.numeric)
cna<-replace(cna, cna==2, "Ampl")
cna<-replace(cna, cna==1 | cna==0 | cna==-1 | cna==-2, "No_ampl")
rownames(cna)<-Drivers
cna<-as.data.frame(t(cna))
GEcancer<-as.data.frame(t(GEcancer))

idx<-match(colnames(cna), colnames(GEcancer))
GEcancer<-GEcancer[,idx]
#########################################################################################
mut_exp<-merge(GEcancer, cna, by = "row.names")
rownames(mut_exp)<-mut_exp$Row.names
mut_exp<-mut_exp[,-1]
library(ggplot2)
#mut_exp<-as.data.frame(t(mut_exp))

########################################
for (i in genes) {
  
}
mut_exp<-mut_exp[order(as.matrix(mut_exp$ZNF667.x), decreasing = TRUE),]
col<-which(mut_exp$ZNF667.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:39]<-apply(mut_exp[,1:39],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$ZNF667.x,  col = mut_exp$col, main= "ZNF667")

mut_exp[,1:39]<-log(mut_exp[,1:39])
boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
###########################################################


