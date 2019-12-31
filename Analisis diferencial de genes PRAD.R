#############################################################################################
#Preparamos la base de datos arn para PRAD
############################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
library(TCGAbiolinks)
query.exp.hg19 <- GDCquery(project = "TCGA-PRAD",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq",
                           file.type  = "results",
                           experimental.strategy = "RNA-Seq",
                           legacy = TRUE)
GDCdownload(query.exp.hg19)


data <- GDCprepare(query.exp.hg19)
library(maftools)


#data.table::data.table(as.data.frame(colData(data)), 
#                       options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#                       rownames = FALSE)

library(SummarizedExperiment)
PRAD_RNA_nonorm<-SummarizedExperiment::assay(data)
PRAD_RNA_nonorm<-as.data.frame(PRAD_RNA_nonorm)
Samples_PRAD<-colnames(PRAD_RNA_nonorm)
Samples_PRAD<-substr(Samples_PRAD,1,12)
idx<-which(duplicated(Samples_PRAD))
Samples_PRAD<-Samples_PRAD[-idx]
PRAD_RNA_nonorm<-PRAD_RNA_nonorm[,-idx]
colnames(PRAD_RNA_nonorm)<-Samples_PRAD
#############################################################################################
#Clinical data for PRAD
#################################################################
Clinical_PRAD<-SummarizedExperiment::colData(data)
Clinical_PRAD<-as.data.frame(Clinical_PRAD)
sampleInfo<-Clinical_PRAD$definition
sampleInfo<-as.data.frame(sampleInfo)
sampleInfo$sample<-Clinical_PRAD$sample
colnames(sampleInfo)<-c("Tumor_Type", "Sample")

sampleinfo<-sampleInfo
remove(sampleInfo)

rownames<-sampleinfo$Sample
sampleinfo<-sampleinfo[,-2]
sampleinfo<-as.data.frame(sampleinfo)
rownames(sampleinfo)<-rownames

rownames<-rownames(sampleinfo)
which(duplicated(rownames))
rownames<-substr(rownames(sampleinfo),1,12)
idx<-which(duplicated(rownames))
sampleinfo<-sampleinfo[-idx,]
rownames<-rownames[-idx]
sampleinfo<-as.data.frame(sampleinfo)
rownames(sampleinfo)<-rownames
levels(sampleinfo$sampleinfo)
table(sampleinfo$sampleinfo)


#table(sampleinfo$sampleinfo)
#idx<-which(sampleinfo$sampleinfo=="Recurrent Solid Tumor")
#sampleinfo<-sampleinfo[-idx,]
#sampleinfo<-as.data.frame(sampleinfo)
#rownames<-rownames[-idx]
#rownames(sampleinfo)<-rownames
sampleinfo$sampleinfo<-droplevels(sampleinfo$sampleinfo)
table(sampleinfo$sampleinfo)

colnames<-colnames(meth_PRAD[,5:88])
colnames<-substr(colnames,1,12)
colnames<-gsub(colnames, pattern = ".", replacement = "-", fixed = TRUE)
idx<-match(colnames, rownames(sampleinfo))
sampleinfo<-sampleinfo[idx,]
sampleinfo<-as.data.frame(sampleinfo)
rownames(sampleinfo)<-colnames
idx<-match(colnames, colnames(PRAD_RNA_nonorm))
PRAD_RNA_nonorm<-as.data.frame(t(PRAD_RNA_nonorm))
PRAD_RNA_nonorm<-PRAD_RNA_nonorm[idx,]
PRAD_RNA_nonorm<-as.data.frame(t(PRAD_RNA_nonorm))
idx<-which(is.na(sampleinfo$sampleinfo))
sampleinfo<-as.data.frame(sampleinfo[-idx,])
PRAD_RNA_nonorm<-as.data.frame(t(PRAD_RNA_nonorm))
PRAD_RNA_nonorm<-PRAD_RNA_nonorm[-idx,]
PRAD_RNA_nonorm<-as.data.frame(t(PRAD_RNA_nonorm))
rownames(sampleinfo)<-colnames(PRAD_RNA_nonorm)
genes<-meth_PRAD$gene_symbol
rownames(PRAD_RNA_nonorm)<-gsub(rownames(PRAD_RNA_nonorm), pattern = "|", replacement = ";", fixed = TRUE)
rownames<-rownames(PRAD_RNA_nonorm)

rownames<-sapply(strsplit(rownames(PRAD_RNA_nonorm), split = ';'), function(x) return(x[1]))
idx<-which(duplicated(rownames))
rownames<-rownames[-idx]
PRAD_RNA_nonorm<-PRAD_RNA_nonorm[-idx,]
rownames(PRAD_RNA_nonorm)<-rownames

#idx<-which(duplicated(rownames))
rownames<-rownames[-idx]
#PRAD_RNA_nonorm<-PRAD_RNA_nonorm[-idx,]
#rownames(PRAD_RNA_nonorm)<-rownames
idx<-match(genes, rownames(PRAD_RNA_nonorm))
PRAD_RNA_nonorm<-PRAD_RNA_nonorm[idx,]
idx<-which(is.na(PRAD_RNA_nonorm[,1]))
PRAD_RNA_nonorm<-PRAD_RNA_nonorm[-idx,]
#####################################################################

design <- model.matrix(~ 0 + sampleinfo$sampleinfo)
colnames(design)<-c("PT","NT")
y <- DGEList(PRAD_RNA_nonorm)
y$counts[1:5,1:5]
v <- voom(y,design,plot = TRUE)

#idx6<-which(is.na(sampleinfo))
#sampleinfo<-sampleinfo[-idx6,]
#sampleinfo<-as.data.frame(sampleinfo)
#idx7<-match(rownames(sampleinfo), colnames(Breast_rna))
#Breast_rna1<-Breast_rna[,idx7]
#Breast_rna<-Breast_rna1
#rm(Breast_rna1)
#y <- DGEList(Breast_rna)
#y$counts[1:5,1:5]
#v <- voom(y,design,plot = TRUE)
colnames(sampleinfo)<-"sampleinfo"
levels(sampleinfo$sampleinfo) <- c(levels(sampleinfo), "PT", "NT") 

cont.matrix <- makeContrasts(PTvsNT=PT - NT,levels=sampleinfo$sampleinfo)
cont.matrix
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
colnames(vfit$coefficients)

design <- model.matrix(~ 0 + sampleinfo$sampleinfo)
colnames(design)<-c("PT","NT")
y <- DGEList(COAD_RNA_nonorm)
y$counts[1:5,1:5]
v <- voom(y,design,plot = TRUE)

efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))
res3<-topTable(efit, adjust="BH",number=nrow(y), sort.by = "P")
topTable(efit, adjust="BH", sort.by = "P")
####################################
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
################
res4<-topTreat(tfit, coef=1, n=Inf)
library(xlsx)
write.xlsx(res4, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/res4_PRAD_methylationgenes.xlsx")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
library(Glimma)
Glimma::glXYPlot(x=efit$coefficients[,1], y=efit$lods[,1],
                 xlab="logFC", ylab="B", main="Top 100 DE genes",
                 counts=y$counts, groups=sampleinfo$sampleinfo, status=dt, side.main="ENTREZID", folder="volcano")
length(which(res4$adj.P.Val < 0.05))

#Gen set enrichment analysis

source("https://bioconductor.org/biocLite.R")
library(CAMERA)
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c6_v5p2.rdata"))##ONCOGENIC
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))##CURATOR
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c7_v5p2.rdata"))###Inmunogenic
library(org.Hs.eg.db)
library(AnnotationDbi)
#rownames(v)<-gsub(rownames(v), pattern = "|", replacement = ";", fixed = TRUE)
#rownames(v)<-sapply(strsplit(rownames(v), split = ';'), function(x) return(x[1]))

ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(v),columns=c("ENTREZID","SYMBOL","GENENAME"), keytype = "SYMBOL")
idx1 <- ids2indices(Hs.c2,id=ann$ENTREZID)
idx2 <- ids2indices(Hs.c6,id=ann$ENTREZID)
idx3 <- ids2indices(Hs.c7,id=ann$ENTREZID)
cam.upvslow_1<- camera(v,idx1,design,contrast=cont.matrix)
cam.upvslow_2 <- camera(v,idx2,design,contrast=cont.matrix)
cam.upvslow_3 <- camera(v,idx3,design,contrast=cont.matrix)
head(cam.upvslow_1,10)
head(cam.upvslow_2,10)
head(cam.upvslow_3,10)
write.xlsx(cam.upvslow_1, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_PRAD_oncogenic_methylationgenes.xlsx")
write.xlsx(cam.upvslow_2, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_PRAD_curator_methylatongenes.xlsx")
write.xlsx(cam.upvslow_3, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_PRAD_inmunogenic_methylationgenes.xlsx")
################################################################################################################################
idx<-match(Drivers,rownames(res4))
res4_drivers<-res4[idx,]
write.xlsx(res4_drivers, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/res4_PRAD_drivers.xlsx")
