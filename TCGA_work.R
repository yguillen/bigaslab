## Yolanda Guilén ##
## October 2018 ##

##################### TCGA dataset RNA-seq pipeline ####################

#required packages
install.packages("DT")
install.packages("devtools")

#basic libraries
library(data.table)
library(ggplot2)
library(nlme)
library(readxl)
library(devtools)

devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')


library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)


## For acute myeloid leukemia
query.exp <- GDCquery(project = "TCGA-AML", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Blood Derived Cancer - Bone Marrow","Bone Marrow Normal"))

GDCdownload(query.exp)

query.exp <- GDCquery(project = "TCGA-BRCA", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
GDCdownload(query.exp)
brca.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "brcaExp.rda")


brca.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "brcaExp.rda")


