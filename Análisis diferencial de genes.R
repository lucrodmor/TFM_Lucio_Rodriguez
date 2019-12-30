#<<<<<<< HEAD
###BRCA#######

##########Asociamos a cada muestra su correspondiente información clínica.
sampleinfo<-as.data.frame(clinical_BRCA[,c(1,41)])
levels(sampleinfo$breast_carcinoma_progesterone_receptor_status)
idx<-which(sampleinfo$breast_carcinoma_progesterone_receptor_status=="Indeterminate")
sampleinfo<-sampleinfo[-idx,]
idx1<-which(sampleinfo$breast_carcinoma_progesterone_receptor_status=="")
sampleinfo<-sampleinfo[-idx1,]

id<-sampleinfo$bcr_patient_barcode
sampleinfo<-as.data.frame(t(sampleinfo))
colnames(sampleinfo)<-id

colnames(sampleinfo) <- gsub(x = colnames(sampleinfo), pattern = "-", replacement = ".")  
idx2<-match(colnames(Breast_rna), colnames(sampleinfo))
which(is.na(idx2))

#rbind(Breast_rna, sampleinfo, by = colnames, na.rm=TRUE)
install.packages("gtools")
library(gtools)
Data<-smartbind(Breast_rna, sampleinfo)
idx4<-rownames(Breast_rna)
Data<-as.data.frame(Data)
rownames(Data)[1:20531]<-idx4
#################################
Data<-Data[-20532,]
rownames(Data)[20532]<-"Progesterona"
idx<-which(is.na(Data[20532,]))
Data<-Data[,-idx]
sampleinfo<-Data[20532,]
Breast_rna1<-Data[1:20531,]
Breast_rna<-as.matrix(Breast_rna)
genes<-rownames(Breast_rna)
Breast_rna<-apply(Breast_rna,2,as.numeric)
rownames(Breast_rna)<-genes
col.has.na <- apply(Breast_rna, 2, function(x){any(is.na(x))})
sum(col.has.na)
Breast_rna <- Breast_rna[,!col.has.na]
sampleinfo<-sampleinfo[,!col.has.na]
###############################
#idx5<-match(colnames(Breast_rna),colnames(sampleinforder))
#sampleinforder<-sampleinforder[,idx5]
#sampleinforder<-sampleinforder[-1,]
#sampleinfo<-sampleinforder
#######################################################################
sampleinfo<-as.data.frame(t(sampleinfo))
colnames(sampleinfo)<-"Progesterona"
levels(sampleinfo$Progesterona)
sampleinfo$Progesterona<-relevel(sampleinfo$Progesterona, ref = "Positive")
design <- model.matrix(~ 0 + sampleinfo$Progesterona)
########################################################################3
#Creamos objeto DGElist
####################################################################
Breast_rna1<-Breast_rna
library(limma)
library(edgeR)

y <- DGEList(Breast_rna)
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
cont.matrix <- makeContrasts(PositivevsNegative=Positive - Negative,levels=sampleinfo$Progesterona)
cont.matrix
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
colnames(vfit$coefficients)


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
write.xlsx(res4, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/res4_brca.xlsx")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
library(Glimma)
Glimma::glXYPlot(x=efit$coefficients[,1], y=efit$lods[,1],
         xlab="logFC", ylab="B", main="Top 100 DE genes",
         counts=y$counts, groups=sampleinfo$Progesterona, status=dt, side.main="ENTREZID", folder="volcano")


#Gen set enrichment analysis

source("https://bioconductor.org/biocLite.R")
library(CAMERA)
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c6_v5p2.rdata"))##ONCOGENIC
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))##CURATOR
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c7_v5p2.rdata"))###Inmunogenic
library(org.Hs.eg.db)
library(AnnotationDbi)
rownames(v)<-gsub(rownames(v), pattern = "|", replacement = ";", fixed = TRUE)
rownames(v)<-sapply(strsplit(rownames(v), split = ';'), function(x) return(x[1]))

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
write.xlsx(cam.upvslow_1, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_oncogenic.xlsx")
write.xlsx(cam.upvslow_2, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_curator.xlsx")
write.xlsx(cam.upvslow_3, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_inmunogenic.xlsx")
################################################################################################################################
#Gen set enrichment with driver genes from methylation data
rownames(Breast_rna)<-gsub(rownames(Breast_rna), pattern = "|", replacement = ";", fixed = TRUE)

rownames<-sapply(strsplit(rownames(Breast_rna), split = ';'), function(x) return(x[1]))
idx<-which(duplicated(rownames))
Breast_rna<-Breast_rna[-idx,]
rownames(Breast_rna)<-sapply(strsplit(rownames(Breast_rna), split = ';'), function(x) return(x[1]))
idx<-match(rownames(METcancer), rownames(Breast_rna))
Breast_rna<-Breast_rna[idx,]
rownames(sampleinfo)<-gsub(rownames(sampleinfo), pattern = ".", replacement = "-", fixed = TRUE)
idx<-match(colnames(Breast_rna), rownames(sampleinfo))
sampleinfo<-sampleinfo[idx,]
sampleinfo<-as.data.frame(sampleinfo)
rownames(sampleinfo)<-colnames(Breast_rna)
idx<-match(Drivers, rownames(Breast_rna))
Breast_rna<-Breast_rna[idx,]



sampleinfo<-as.data.frame(t(sampleinfo))
colnames(sampleinfo)<-"Progesterona"
levels(sampleinfo$Progesterona)
sampleinfo$Progesterona<-relevel(sampleinfo$Progesterona, ref = "Positive")
design <- model.matrix(~ + sampleinfo$Progesterona)


ncol(y$counts)
nrow(design)
nrow(sampleinfo)
y <- DGEList(Breast_rna)
y$counts[1:5,1:5]
v <- voom(y,design,plot = TRUE)
idx<-c(102, 111, 113, 150, 151, 152, 154, 155, 156, 157, 190)
Breast_rna<-Breast_rna[,-idx]
sampleinfo<-sampleinfo[-idx,]

sampleinfo<-as.data.frame(sampleinfo)
colnames(sampleinfo)<-"Progesterona"
levels(sampleinfo$Progesterona)
#sampleinfo$Progesterona<-relevel(sampleinfo$Progesterona, ref = "Positive")
design <- model.matrix(~ + sampleinfo$Progesterona)
y <- DGEList(Breast_rna)
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
cont.matrix <- makeContrasts(PositivevsNegative=Positive - Negative,levels=sampleinfo$Progesterona)
cont.matrix
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
colnames(vfit$coefficients)


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
write.xlsx(res4, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/res4_brca_predictivegenes.xlsx")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
library(Glimma)
Glimma::glXYPlot(x=efit$coefficients[,1], y=efit$lods[,1],
                 xlab="logFC", ylab="B", main="Top 100 DE genes",
                 counts=y$counts, groups=sampleinfo$Progesterona, status=dt, side.main="ENTREZID", folder="volcano")


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
write.xlsx(cam.upvslow_1, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_oncogenic_drivergenes.xlsx")
write.xlsx(cam.upvslow_2, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_curator_drivergenes.xlsx")
write.xlsx(cam.upvslow_3, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_inmunogenic_drivergenes.xlsx")
################################################################################################################################
################################################################################################################################
#Gen set enrichment with 10000 genes from methylation data
rownames(Breast_rna)<-gsub(rownames(Breast_rna), pattern = "|", replacement = ";", fixed = TRUE)

rownames<-sapply(strsplit(rownames(Breast_rna), split = ';'), function(x) return(x[1]))
idx<-which(duplicated(rownames))
Breast_rna<-Breast_rna[-idx,]
rownames(Breast_rna)<-sapply(strsplit(rownames(Breast_rna), split = ';'), function(x) return(x[1]))
rownames(METcancer)<-gsub(rownames(METcancer), pattern = "|", replacement = ";", fixed = TRUE)
rownames(METcancer)<-sapply(strsplit(rownames(METcancer), split = ';'), function(x) return(x[1]))
idx<-match(rownames(METcancer), rownames(Breast_rna))
Breast_rna<-Breast_rna[idx,]
rownames(sampleinfo)<-gsub(rownames(sampleinfo), pattern = ".", replacement = "-", fixed = TRUE)
idx<-match(colnames(Breast_rna), rownames(sampleinfo))

##################################
##########Asociamos a cada muestra su correspondiente información clínica.
sampleinfo<-as.data.frame(clinical_BRCA[,c(1,41)])
levels(sampleinfo$breast_carcinoma_progesterone_receptor_status)
idx<-which(sampleinfo$breast_carcinoma_progesterone_receptor_status=="Indeterminate")
sampleinfo<-sampleinfo[-idx,]
idx1<-which(sampleinfo$breast_carcinoma_progesterone_receptor_status=="")
sampleinfo<-sampleinfo[-idx1,]

id<-sampleinfo$bcr_patient_barcode
sampleinfo<-as.data.frame(t(sampleinfo))
colnames(sampleinfo)<-id

colnames(sampleinfo) <- gsub(x = colnames(sampleinfo), pattern = "-", replacement = ".")  
idx2<-match(colnames(Breast_rna), colnames(sampleinfo))
which(is.na(idx2))
sampleinfo<-as.data.frame(t(sampleinfo))
rownames<-rownames(sampleinfo)
sampleinfo<-sampleinfo[,-1]
sampleinfo<-as.data.frame(sampleinfo)
rownames(sampleinfo)<-rownames
idx<-match(colnames(Breast_rna), rownames(sampleinfo))
sampleinfo<-sampleinfo[idx,]
rownames(sampleinfo)<-colnames(Breast_rna)
idx<-which(is.na(sampleinfo))
rownames<-rownames(sampleinfo)
sampleinfo<-sampleinfo[-idx,]
rownames<-rownames[-idx]
rownames(sampleinfo)<-rownames
Breast_rna<-Breast_rna[-idx]
####################################



#sampleinfo<-as.data.frame(t(sampleinfo))
colnames(sampleinfo)<-"Progesterona"
levels(sampleinfo$Progesterona)
sampleinfo$Progesterona<-relevel(sampleinfo$Progesterona, ref = "Positive")
design <- model.matrix(~0 + sampleinfo$Progesterona)


ncol(y$counts)
nrow(design)
nrow(sampleinfo)
row.has.na <- apply(Breast_rna, 1, function(x){any(is.na(x))})
sum(row.has.na)
Breast_rna <- Breast_rna[!row.has.na,]


#sampleinfo$Progesterona<-relevel(sampleinfo$Progesterona, ref = "Positive")
design <- model.matrix(~0 + sampleinfo$Progesterona)
y <- DGEList(Breast_rna)
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
cont.matrix <- makeContrasts(PositivevsNegative=Positive - Negative,levels=sampleinfo$Progesterona)
cont.matrix
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
colnames(vfit$coefficients)


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
write.xlsx(res4, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/res4_brca_methylationgenes.xlsx")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
library(Glimma)
Glimma::glXYPlot(x=efit$coefficients[,1], y=efit$lods[,1],
                 xlab="logFC", ylab="B", main="Top 100 DE genes",
                 counts=y$counts, groups=sampleinfo$Progesterona, status=dt, side.main="ENTREZID", folder="volcano")


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
write.xlsx(cam.upvslow_1, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_oncogenic_methylationgenes.xlsx")
write.xlsx(cam.upvslow_2, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_curator_methylatongenes.xlsx")
write.xlsx(cam.upvslow_3, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_inmunogenic_methylationgenes.xlsx")
################################################################################################################################
idx<-match(Drivers,rownames(v))
v<-v[idx,]
nrow(v)
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
write.xlsx(cam.upvslow_1, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_oncogenic_methylationgenes.xlsx")
write.xlsx(cam.upvslow_2, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_curator_methylatongenes.xlsx")
write.xlsx(cam.upvslow_3, "C:/Users/Lucio/Desktop/Bioinformatica/TFM/PEC2/gsea_brca_inmunogenic_methylationgenes.xlsx")
######################