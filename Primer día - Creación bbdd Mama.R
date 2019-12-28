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
########################################################################################
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
#LUAD

#Importamos LUAD_Samples.txt
library(readr)
LUAD_Samples <- read_delim("LUAD_Samples.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

LUAD_Samples<-LUAD_Samples[,2:3]
LUAD_Samples$SAMPLE_ID<-substr(LUAD_Samples$SAMPLE_ID,1,nchar(LUAD_Samples$SAMPLE_ID)-3)
rownames<-LUAD_Samples$SAMPLE_ID
LUAD_Samples<-as.data.frame(t(LUAD_Samples))
colnames(LUAD_Samples)<-rownames
rownames<-LUAD_Samples[,1]
rownames(LUAD_Samples)<-rownames
LUAD_Samples<-LUAD_Samples[,-1]
save(LUAD_Samples, file = "LUAD_Samples.RData")
############################################################
#Reducimos la base de dato a solo los pacientes que se definan por tener Cancer de Pulmon.
LUAD_Samples<-as.data.frame(t(LUAD_Samples))
idx<-match(row.names(LUAD_Samples), colnames(rna_expression))
LUAD_rna<-as.data.frame(rna_expression[,idx])
save(LUAD_rna, file = "LUAD_rna.RData")
########################################################################################

#PRAD

#Importamos PRAD_Samples.txt
library(readr)
PRAD_Samples <- read_delim("PRAD_Samples.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

PRAD_Samples<-PRAD_Samples[,2:3]
PRAD_Samples$SAMPLE_ID<-substr(PRAD_Samples$SAMPLE_ID,1,nchar(PRAD_Samples$SAMPLE_ID)-3)
rownames<-PRAD_Samples$SAMPLE_ID
PRAD_Samples<-as.data.frame(t(PRAD_Samples))
colnames(PRAD_Samples)<-rownames
rownames<-PRAD_Samples[,1]
rownames(PRAD_Samples)<-rownames
PRAD_Samples<-PRAD_Samples[,-1]
save(PRAD_Samples, file = "PRAD_Samples.RData")
############################################################
#Reducimos la base de dato a solo los pacientes que se definan por tener Cancer de Pulmon.
PRAD_Samples<-as.data.frame(t(PRAD_Samples))
idx<-match(row.names(PRAD_Samples), colnames(rna_expression))
PRAD_rna<-as.data.frame(rna_expression[,idx])
save(PRAD_rna, file = "PRAD_rna.RData")
########################################################################################







