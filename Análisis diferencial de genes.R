<<<<<<< HEAD
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

rbind(Breast_rna, sampleinfo, by = colnames, na.rm=TRUE)
install.packages("gtools")
library(gtools)
Data<-smartbind(Breast_rna, sampleinfo)
idx4<-rownames(Breast_rna)
Data<-as.data.frame(Data)
rownames(Data[1:20531,])<-idx4
idx5<-match(colnames(Breast_rna),colnames(sampleinforder))
sampleinforder<-sampleinforder[,idx5]
sampleinforder<-sampleinforder[-1,]
sampleinfo<-sampleinforder
#######################################################################
sampleinfo<-as.data.frame(t(sampleinfo))
colnames(sampleinfo)<-"Progesterona"
levels(sampleinfo$Progesterona)
design <- model.matrix(~ 0 + sampleinfo$Progesterona)
########################################################################3
#Creamos objeto DGElist
####################################################################

library(limma)
library(edgeR)

y <- DGEList(Breast_rna)
y$counts[1:5,1:5]
v <- voom(y,design,plot = TRUE)
idx6<-which(is.na(sampleinfo))
sampleinfo<-sampleinfo[-idx6,]
sampleinfo<-as.data.frame(sampleinfo)
idx7<-match(rownames(sampleinfo), colnames(Breast_rna))
Breast_rna1<-Breast_rna[,idx7]
Breast_rna<-Breast_rna1
rm(Breast_rna1)
y <- DGEList(Breast_rna)
y$counts[1:5,1:5]
v <- voom(y,design,plot = TRUE)
cont.matrix <- makeContrasts(NegativevsPositive=Negative - Positive,levels=sampleinfo$Progesterona)
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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
library(Glimma)
glXYPlot(x=efit$coefficients[,1], y=efit$lods[,1],
         xlab="logFC", ylab="B", main="Top 100 DE genes",
         counts=y$counts, groups=sampleinfo$Progesterona, status=dt, side.main="ENTREZID", folder="volcano")


=======
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

rbind(Breast_rna, sampleinfo, by = colnames, na.rm=TRUE)
install.packages("gtools")
library(gtools)
Data<-smartbind(Breast_rna, sampleinfo)
idx4<-rownames(Breast_rna)
Data<-as.data.frame(Data)
rownames(Data[1:20531,])<-idx4
idx5<-match(colnames(Breast_rna),colnames(sampleinforder))
sampleinforder<-sampleinforder[,idx5]
sampleinforder<-sampleinforder[-1,]
sampleinfo<-sampleinforder
#######################################################################
sampleinfo<-as.data.frame(t(sampleinfo))
colnames(sampleinfo)<-"Progesterona"
levels(sampleinfo$Progesterona)
design <- model.matrix(~ 0 + sampleinfo$Progesterona)
########################################################################3
#Creamos objeto DGElist
####################################################################

library(limma)
library(edgeR)

y <- DGEList(Breast_rna)
y$counts[1:5,1:5]
v <- voom(y,design,plot = TRUE)
idx6<-which(is.na(sampleinfo))
sampleinfo<-sampleinfo[-idx6,]
sampleinfo<-as.data.frame(sampleinfo)
idx7<-match(rownames(sampleinfo), colnames(Breast_rna))
Breast_rna1<-Breast_rna[,idx7]
Breast_rna<-Breast_rna1
rm(Breast_rna1)
y <- DGEList(Breast_rna)
y$counts[1:5,1:5]
v <- voom(y,design,plot = TRUE)
cont.matrix <- makeContrasts(NegativevsPositive=Negative - Positive,levels=sampleinfo$Progesterona)
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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
library(Glimma)
glXYPlot(x=efit$coefficients[,1], y=efit$lods[,1],
         xlab="logFC", ylab="B", main="Top 100 DE genes",
         counts=y$counts, groups=sampleinfo$Progesterona, status=dt, side.main="ENTREZID", folder="volcano")


>>>>>>> f84912fa63e203ac39d02adbc1278f2b97769360
