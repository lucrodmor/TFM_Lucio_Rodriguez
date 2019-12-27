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
colnames(rna_expression)<-gsub(colnames(rna_expression), pattern = ".", replacement = "-", fixed = TRUE)
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
#####COAD#######
###################################################################################


#########################################################
#Crearemos 5 bases de datos de expresi?n de ARN.  Cada una de ellas referente a cada tipo de tumor estudiado: 
#Hemos descargado los archivos de expresi?n de arn a nuestro escritorio. Procederemos a importarlos. 

rna_expression<-read.table(file = 'pancan_normalized_rnaseq.tsv', sep = '\t', header = TRUE)
colnames(rna_expression)<-substr(colnames(rna_expression),1,nchar(colnames(rna_expression))-13)
#rna_expression<-as.data.frame(t(rna_expression))
rownames<-rna_expression[,1]
rownames(rna_expression)<-rownames
rna_expression<-rna_expression[,-1]
colnames(rna_expression)<-gsub(colnames(rna_expression), pattern = ".", replacement = "-", fixed = TRUE)
save(rna_expression, file = "rna_expression.RData")
#############################################################
head(rna_expression)
#COAD

#Importamos COAD_Samples.txt
library(readr)
COAD_Samples <- read_delim("COAD_Samples.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

COAD_Samples<-COAD_Samples[,2:3]
COAD_Samples$SAMPLE_ID<-substr(COAD_Samples$SAMPLE_ID,1,nchar(COAD_Samples$SAMPLE_ID)-3)
rownames<-COAD_Samples$SAMPLE_ID
COAD_Samples<-as.data.frame(t(COAD_Samples))
colnames(COAD_Samples)<-rownames
rownames<-COAD_Samples[,1]
rownames(COAD_Samples)<-rownames
COAD_Samples<-COAD_Samples[,-1]
save(COAD_Samples, file = "COAD_Samples.RData")
############################################################
#Reducimos la base de dato a solo los pacientes que se definan por tener Cancer de Mama.
COAD_Samples<-as.data.frame(t(COAD_Samples))
idx<-match(row.names(COAD_Samples), colnames(rna_expression))
COAD_rna<-as.data.frame(rna_expression[,idx])
save(COAD_rna, file = "COAD_rna.RData")






