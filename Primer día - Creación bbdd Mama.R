if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
#########################################################
#Crearemos 5 bases de datos de expresi?n de ARN.  Cada una de ellas referente a cada tipo de tumor estudiado: 
#Hemos descargado los archivos de expresi?n de arn a nuestro escritorio. Procederemos a importarlos. 

rna_expression<-read.table(file = 'pancan_normalized_rnaseq.tsv', sep = '\t', header = TRUE)
colnames(rna_expression)<-substr(colnames(rna_expression),1,nchar(colnames(rna_expression))-13)
#rna_expression<-as.data.frame(t(rna_expression))
rownames<-rna_expression[,1]
rownames(rna_expression)<-rownames
rna_expression<-rna_expression[,-1]
save(rna_expression, file = "rna_expression.RData")
#############################################################
head(rna_expression)
#Breast

#Importamos Breast_Samples.txt
library(readr)
Breast_Samples <- read_delim("Breast_Samples.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

Breast_Samples<-Breast_Samples[,2:3]
Breast_Samples$SAMPLE_ID<-substr(Breast_Samples$SAMPLE_ID,1,nchar(Breast_Samples$SAMPLE_ID)-3)
rownames<-Breast_Samples$SAMPLE_ID
Breast_Samples<-as.data.frame(t(Breast_Samples))
colnames(Breast_Samples)<-rownames
rownames<-Breast_Samples[,1]
rownames(Breast_Samples)<-rownames
Breast_Samples<-Breast_Samples[,-1]
save(Breast_Samples, file = "Breast_Samples.RData")
############################################################
#Reducimos la base de dato a solo los pacientes que se definan por tener Cancer de Mama.
Breast_Samples<-as.data.frame(t(Breast_Samples))
idx<-match(row.names(Breast_Samples), colnames(rna_expression))
Breast_rna<-as.data.frame(rna_expression[,idx])
save(Breast_rna, file = "Breast_rna.RData")
#####################################################################################
#Creación bbdd Metilación
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MethylMix")

install.packages(c( "foreach", "doParallel") )
source("http://bioconductor.org/biocLite.R")
library(doParallel)
cl<-makeCluster(5)
registerDoParallel(cl)

##################################
library(MethylMix)
cancerSite<-"COAD"
targetDirectory<-paste0(getwd(), "/")
#####################################

# Downloading methylation data
METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory)
GetData(cancerSite, targetDirectory)

# Processing methylation data
METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)



object.size(METdirectories)
gc()

# Saving methylation processed data
saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, 
                                        "_Processed.rds"))

# Downloading gene expression data
GEdirectories <- Download_GeneExpression(cancerSite, targetDirectory)
# Processing gene expression data
GEProcessedData <- Preprocess_GeneExpression(cancerSite, GEdirectories)
# Saving gene expression processed data
saveRDS(GEProcessedData, file = paste0(targetDirectory, "GE_", cancerSite, "_Processed.rds"))

# Clustering probes to genes methylation data
METProcessedData <- readRDS(paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])

# Putting everything together in one file
toSave <- list(METcancer = res[[1]], METnormal = res[[2]], GEcancer = GEProcessedData[[1]], 
               GEnormal = GEProcessedData[[2]], ProbeMapping = res$ProbeMapping)
saveRDS(toSave, file = paste0(targetDirectory, "data_", cancerSite, ".rds"))
######################################################################################################

######################################################################################################


########################################################################################################
cancerSite<-"COAD"
targetDirectory<-paste0(getwd(), "/")
GetData(cancerSite, targetDirectory)
#####################################################################################################
#####################################################################################################
install.packages("Github")

install.packages("TCGA2STAT_1.0.tar.gz", repos = NULL, type = "source")









