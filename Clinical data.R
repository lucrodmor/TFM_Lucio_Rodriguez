
source("http://bioconductor.org/biocLite.R")
useDevel()
BiocManager::install("TCGAbiolinks")
###LUAD######
# Step 1.2 download expression data
# Fetch clinical data directly from the clinical XML files.
# if barcode is not set, it will consider all samples.
# We only set it to make the example faster
query <- GDCquery(project = "TCGA-LUAD",
                  file.type = "xml",
                  data.category = "Clinical") class(query)

class(query$results)
x<-which(duplicated(query$results[[1]]$cases))
query$results[[1]]<-query$results[[1]][-x]

#########################

GDCdownload(query)
##########################
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

clinical_LuAD<-clinical

data.table::data.table(clinical_LuAD, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

##################################
#BREAST#######
# Step 1.2 download expression data
# Fetch clinical data directly from the clinical XML files.
# if barcode is not set, it will consider all samples.
# We only set it to make the example faster
query <-TCGAbiolinks::GDCquery(project = "TCGA-BRCA",
                  file.type = "xml",
                  data.category = "Clinical") 

class(query$results)
x<-which(duplicated(query$results[[1]]$cases))
query$results[[1]]<-query$results[[1]][-x]

#########################

TCGAbiolinks::GDCdownload(query)
##########################
clinical_BRCA <-TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "patient")
######################################################
saveRDS(clinical_BRCA, file = "clinical_BRCA.RData")
saveRDS(clinical_LuAD, file = "clinical_LuAD.RData")
#####################################################3
###COAD######
# Step 1.2 download expression data
# Fetch clinical data directly from the clinical XML files.
# if barcode is not set, it will consider all samples.
# We only set it to make the example faster
query <-TCGAbiolinks::GDCquery(project = "TCGA-COAD",
                  file.type = "xml",
                  data.category = "Clinical")

class(query$results)
x<-which(duplicated(query$results[[1]]$cases))
query$results[[1]]<-query$results[[1]][-x]

#########################

TCGAbiolinks::GDCdownload(query)
##########################
clinical <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "patient")

clinical_COAD<-clinical

data.table::data.table(clinical_LuAD, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

##################################
#PROAD#######
# Step 1.2 download expression data
# Fetch clinical data directly from the clinical XML files.
# if barcode is not set, it will consider all samples.
# We only set it to make the example faster
query <-TCGAbiolinks::GDCquery(project = "TCGA-PRAD",
                               file.type = "xml",
                               data.category = "Clinical") 

class(query$results)
x<-which(duplicated(query$results[[1]]$cases))
query$results[[1]]<-query$results[[1]][-x]

#########################

TCGAbiolinks::GDCdownload(query)
##########################
clinical_PRAD <-TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "patient")
######################################################
saveRDS(clinical_PRAD, file = "clinical_PRAD.RData")
#####################################################3
query <-TCGAbiolinks::GDCquery(project = "TCGA-PRAD",
                               file.type = "bam",
                               data.category = "Sequencing Reads") 

class(query$results)
x<-which(duplicated(query$results[[1]]$cases))
query$results[[1]]<-query$results[[1]][-x]

#########################

TCGAbiolinks::GDCdownload(query)
##########################

######################################################
saveRDS(clinical_PRAD, file = "clinical_PRAD.RData")
#####################################################3
#Download ARNSeq data for PRAD
query.exp<- TCGAbiolinks::GDCquery(project = "TCGA-PRAD",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           workflow.type="HTSeq-Counts",
                           legacy = TRUE)









