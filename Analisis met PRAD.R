library(readr)
meth_PRAD <- readLines('randomized_meth_PRAD.txt')
meth_PRAD <- read.table(textConnection(meth_PRAD[-2]), header = TRUE, sep = "\t")
colnames(meth_PRAD)[1:4] <- c('probe', 'gene_symbol', 'chr', 'start') 
dim(meth_PRAD)
############################################################################

install.packages("stringr")
library(stringr)
Samples<-colnames(meth_PRAD[,5:141])
Samples<-gsub(Samples, pattern = ".", replacement = "-", fixed = TRUE)
idx2<-sapply(strsplit(Samples, split = '-'), function(x) return(x[4]))
table(idx2)
#posiciones11b<-which(idx2=="11B")
posiciones11a<-which(idx2=="11A")
Samples1<-Samples
Samples1<-Samples1[posiciones11a]
#Samples2<-Samples[posiciones11b]
Samples_Normal<-Samples1
idx3<-match(Samples_Normal, Samples)
Samples_Tumor<-Samples[-idx3]
##########################################################################

colnames(meth_PRAD)<-gsub(colnames(meth_PRAD), pattern = ".", replacement = "-", fixed = TRUE)
idx2<-match( Samples_Tumor,colnames(meth_PRAD))
MetCancer<-meth_PRAD[,idx2]
idx2<-match( Samples_Normal,colnames(meth_PRAD))
METnormal<-meth_PRAD[,idx2]
load("C:/Users/Lucio/Desktop/Bioinformatica/TFM/Bases de datos/PRAD_rna.RData")
#colnames(PRAD_rna)<-gsub(colnames(PRAD_rna), pattern = ".", replacement = "-", fixed = TRUE)
colnames(MetCancer)<-substr(colnames(MetCancer),1,12)
METcancer<-MetCancer
remove(MetCancer, PRAD_Samples, rna_expression)
METcancer<-as.data.frame(t(METcancer))
PRAD_rna<-as.data.frame(t(PRAD_rna))
idx2<-match(rownames(METcancer[,1:127]),rownames(PRAD_rna))
PRAD_rna<-PRAD_rna[idx2,]
METcancer<-as.data.frame(t(METcancer))
PRAD_rna<-as.data.frame(t(PRAD_rna))

#########################################################################
Samples<-colnames(meth_PRAD[,5:141])
Samples<-substr(Samples,1,12)
#idx<-colnames(PRAD_rna)
#idx2<-match(Samples, idx)
#PRAD_rna<-PRAD_rna[,idx2]
GEcancer<-PRAD_rna
remove(PRAD_rna)
#METcancer<-METCancer
###########################################
#Nombre de filas por simbolo del gen
METcancer$gene_symbol<-meth_PRAD$gene_symbol
#METcancer<-cbind(meth_PRAD[,1:4],METcancer)
METnormal$gene_symbol<-meth_PRAD$gene_symbol
#METnormal<-cbind(meth_PRAD[,1:4],METnormal)
#Miramos cuales son NA en Genes
idx<-which(is.na(METcancer$gene_symbol))
METcancer<-METcancer[-idx,]
idx<-which(is.na(METnormal$gene_symbol))
METnormal<-METnormal[-idx,]
#Miramos so hay duplicados
idx<-which(duplicated(METnormal$gene_symbol))
METnormal<-METnormal[-idx,]
METcancer<-METcancer[-idx,]
rownames(METcancer)<-METcancer$gene_symbol
rownames(METnormal)<-METnormal$gene_symbol
idx<-match(rownames(METcancer), rownames(GEcancer))
#####################################
#METcancer<-METcancer[,c(2,5:73)]
#METnormal<-as.data.frame(METnormal[,c(2,5:15)])
#GEcancer<-as.matrix(GEcancer)
#METnormal<-as.matrix(METnormal)
#METcancer<-as.matrix(METcancer)
############################################################################
#METnormal[,1]
#rownames(METnormal)<-sapply(strsplit(METnormal[,1], split = ';'), function(x) return(x[1]))
#rownames(METcancer)<-sapply(strsplit(METcancer[,1], split = ';'), function(x) return(x[1]))
rownames<-rownames(GEcancer)
GEcancer1<-GEcancer
rownames(GEcancer)<-gsub(rownames(GEcancer), pattern = "|", replacement = ";", fixed = TRUE)
GEcancer<-as.matrix(GEcancer)
rownames(GEcancer)<-sapply(strsplit(rownames(GEcancer), split = ';'), function(x) return(x[1]))
idx<-match(row.names(METnormal),row.names(GEcancer))
GEcancer<-as.data.frame(GEcancer)
#rownames(GEcancer)<-rownames
GEcancer<-GEcancer[idx,]
METcancer<-METcancer[,-128]
METnormal<-METnormal[,-11]
remove(GEcancer1, meth_PRAD)
GEcancer<-as.matrix(GEcancer)
METcancer<-as.matrix(METcancer)
METnormal<-as.matrix(METnormal)
#METcancer<-METcancer[,2:195]
#METnormal<-METnormal[,2:23]
#METnormal<-as.matrix(METnormal)
#METcancer<-as.matrix(METcancer)
#rownames(GEcancer)<-gsub(rownames(GEcancer), pattern = "|", replacement = ";", fixed = TRUE)
#rownames(GEcancer)<-sapply(strsplit(rownames(GEcancer), split = ';'), function(x) return(x[1]))
#GEcancer<-as.data.frame(GEcancer)
#GEcancer<-as.matrix(GEcancer)
#METnormal<-as.data.frame(METnormal)



genes<-rownames(METnormal)
#METnormal<-apply(METnormal,2,as.numeric)
#rownames(METnormal)<-genes

#genes<-rownames(METcancer)
#METnormal<-apply(METcancer,2,as.numeric)
#rownames(METcancer)<-genes


#METnormal<-as.matrix(METnormal)
#METcancer<-as.data.frame(METcancer)
#METcancer<-as.matrix(METcancer)
#colnames(METnormal)<-substr(colnames(METnormal),1,12)
#levels(METnormal)
#idx<-which(is.na(METnormal))
#QUITAMOS NAs
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
#########################################################
####################
######EJEMPLO##########
library(MethylMix)
saveRDS(METnormal, file = "METnormal_PRAD.rds")
saveRDS(METcancer,file = "METcancer_PRAD.rds")
saveRDS(GEcancer, file = "GEcancer_PRAD.rds")


MethylMixResults<-MethylMix(METcancer, GEcancer, METnormal, filter = TRUE, NoNormalMode = FALSE, listOfGenes = rownames(METcancer), OutputRoot = "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2")

Drivers<-MethylMixResults$MethylationDrivers
Drivers
plots<-MethylMix_PlotModel("NEU1", MethylMixResults, METcancer, METnormal = METnormal, GEcancer = GEcancer)
plots$CorrelationPlot
plots$MixtureModelPlot
plots$MixtureModelPlot
MethylMixResults$Classifications
MethylMixResults$MethylationStates
# Plot MGMT also with its normal methylation variation
plots <- MethylMix_PlotModel("CASZ1", MethylMixResults, METcancer, METnormal = METnormal)
plots$MixtureModelPlot
plots <- MethylMix_PlotModel("KLK1", MethylMixResults, METcancer, METnormal = METnormal)
plots$MixtureModelPlot
# Also plot the inverse correlation with gene expression (creates two
# separate plots)
plots <- MethylMix_PlotModel("KLK1", MethylMixResults, METcancer, GEcancer, 
                             METnormal)
plots$MixtureModelPlot
plots$CorrelationPlot
# Plot all functional and differential genes
gene<-Drivers
for (gene in MethylMixResults$MethylationDrivers) {
  plots<-MethylMix_PlotModel(gene, MethylMixResults, METcancer, METnormal = METnormal)
  plots$MixtureModelPlot
}

############################################################################
write.xlsx2(Drivers, file = "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/Drivers_PRAD.xlsx")
library(xlsx)      
############################################################################
library(ConsensusClusterPlus)
MethylMixResults<- MethylMix(METcancer, GEcancer, METnormal)

ComplexHeatmap::Heatmap(Metilation)
Metilation<-cbind(as.data.frame(METcancer), as.data.frame(METnormal))

Samplestype<-c(1:123)
Samplestype[1:118]<-"Tumor"
Samplestype[119:123]<-"Normal"
library(gplots)
DMvalues<- MethylMixResults$MethylationStates
cons_cluster<- ConsensusClusterPlus(d=DMvalues, maxK=10, reps=1000, pItem=0.8, distance='euclidean', clusterAlg="km")
heatmap.2(as.matrix(METcancer),scale='row',col=bluered(149),trace='none',
          main = "DM Values Centered CpG sites",margins = c(9, 22),density.inf="none",symkey=TRUE)
BiocManager::install("pvcluster")

pvclust(as.matrix(Metilation), method.hclust="average",
        method.dist="correlation", use.cor="pairwise.complete.obs",
        nboot=1000, parallel=FALSE, r=seq(.5,1.4,by=.1),
        store=FALSE, weight=FALSE, iseed=NULL, quiet=FALSE)

library(pvclust)
Boxplot<-as.data.frame(t(DMvalues))
par(mfrow=c(1,1))
boxplot(Boxplot, las =2, col = rainbow(39), main = "Predictive Genes DM Values")


idx1<-match( Drivers, rownames(Metilation))
Metilation<-Metilation[idx1,]
BetaValues<-Metilation
boxplot(t(BetaValues), las = 2, col = rainbow(39), main = "Predictive Genes Beta Values")

METcancer_PRAD<-as.data.frame(t(METcancer_PRAD))
METnormal_PRAD<-as.data.frame(t(METnormal_PRAD))


#############################################################################################
#Análisis expresión ARN y CNVmut
library(readr)
cna <- read_delim("cna_PRAD.txt", "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
rownames<-cna$SAMPLE_ID
cna<-as.data.frame(t(cna))
colnames(cna)<-rownames
cna<-cna[-1,]
cna<-cna[-1,]
colnames(cna)<-substr(colnames(cna),1,12)
GEcancer<-as.data.frame(t(GEcancer))
cna<-as.data.frame(t(cna))
idx<-match(rownames(GEcancer),rownames(cna))
cna<-cna[idx,]

#cna<-as.numeric(cna)
genes<-colnames(cna)
samples<-rownames(cna)
cna<-apply(cna, 2, as.numeric)
cna<-replace(cna, cna==2, "Ampl")
cna<-replace(cna, cna==1 | cna==0 | cna==-1 | cna==-2, "No_ampl")
rownames(cna)<-samples

GEcancer<-as.data.frame(t(GEcancer))
cna<-as.data.frame(t(cna))
idx<-match(rownames(cna), rownames(GEcancer))
GEcancer<-GEcancer[idx,]
cna<-as.data.frame(t(cna))
GEcancer<-as.data.frame(t(GEcancer))
#########################################################################################
mut_exp<-merge(GEcancer, cna, by = "row.names")
rownames(mut_exp)<-mut_exp$Row.names
mut_exp<-mut_exp[,-1]
library(ggplot2)
#mut_exp<-as.data.frame(t(mut_exp))

########################################

mut_exp<-mut_exp[order(as.matrix(mut_exp$CASZ1.x), decreasing = TRUE),]
genes

col<-which(mut_exp$CASZ1.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$CASZ1.x,  col = mut_exp$col, main= "CACYBP")

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
########################################

########################################

mut_exp<-mut_exp[order(as.matrix(mut_exp$ATP10D.x), decreasing = TRUE),]
genes

col<-which(mut_exp$ATP10D.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$ATP10D.x,  col = mut_exp$col, main= "ATP10D")

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
#######################################

mut_exp<-mut_exp[order(as.matrix(mut_exp$ARPC1B.x), decreasing = TRUE),]
genes

col<-which(mut_exp$ARPC1B.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$ARPC1B.x,  col = mut_exp$col, main= "ARPC1B")

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))

############################################
mut_exp<-mut_exp[order(as.matrix(mut_exp$ST6GAL1.x), decreasing = TRUE),]
genes

col<-which(mut_exp$ST6GAL1.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$ST6GAL1.x,  col = mut_exp$col, main= "ST6GAL1")

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))

mut_exp<-mut_exp[order(as.matrix(mut_exp$RHOD.x), decreasing = TRUE),]
genes

col<-which(mut_exp$RHOD.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$RHOD.x,  col = mut_exp$col, main= "RHOD")

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
mut_exp<-mut_exp[order(as.matrix(mut_exp$CTF1.x), decreasing = TRUE),]
genes

col<-which(mut_exp$CTF1.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$CTF1.x,  col = mut_exp$col, main= "CTF1")

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
mut_exp<-mut_exp[order(as.matrix(mut_exp$SOX9.x), decreasing = TRUE),]
genes

col<-which(mut_exp$SOX9.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$SOX9.x,  col = mut_exp$col, main= "SOX9")

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
mut_exp<-mut_exp[order(as.matrix(mut_exp$GABRG3.x), decreasing = TRUE),]
genes

col<-which(mut_exp$GABRG3.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$GABRG3.x,  col = mut_exp$col, main= "GABRG3")

#mut_exp[,1:73]<-log(mut_exp[,1:73])
#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
mut_exp<-mut_exp[order(as.matrix(mut_exp$DNAJA2.x), decreasing = TRUE),]
genes

col<-which(mut_exp$DNAJA2.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$DNAJA2.x,  col = mut_exp$col, main= "DNAJA2")


#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
mut_exp<-mut_exp[order(as.matrix(mut_exp$CLCN4.x), decreasing = TRUE),]
genes

col<-which(mut_exp$CLCN4.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$CLCN4.x,  col = mut_exp$col, main= "CLCN4")


#boxplot(mut_exp$TAP1.x, mut_exp$AKT1.y, col = rainbow(2))
mut_exp<-mut_exp[order(as.matrix(mut_exp$MPP6.x), decreasing = TRUE),]
genes

col<-which(mut_exp$MPP6.y=="Ampl")
mut_exp$col<-"blue"
mut_exp$col[c(col)]<-"red"
rownames<-rownames(mut_exp)
mut_exp[,1:127]<-apply(mut_exp[,1:127],2,as.numeric)
#mut_exp<-as.data.frame(mut_exp)
#rownames(mut_exp)<-rownames
barplot(mut_exp$MPP6.x,  col = mut_exp$col, main= "MPP6")



#########################################################################
#Comparamos muestras de tejido tumoral vs Normal

genes<-Drivers
par(mfrow=c(3,3))
for (i in genes)
{
  boxplot(METcancer[i,],METnormal[i,], las = 2, col = rainbow(2), main = i, names = c("tumor", "normal"))
}
