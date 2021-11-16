### GSEA analyses ###
setwd("/Users/instalar/Desktop/RNASeq_CPorch/")

source("https://bioconductor.org/biocLite.R")

biocLite("gage")
biocLite("gageData")

library(gage)
library(gageData)
library(readxl)


## Example manual gage
data(gse16873)
data(go.sets.hs)
data(go.subs.hs)

data(go.sets.mm)

#numeric vector for control samples columns
hn=(1:6)*2-1
#numeric vector for experiment/condition columns
dcis=(1:6)*2

data(kegg.gs)
data(go.gs)

gse16873.go.p <- gage(gse16873, gsets = go.gs, ref = hn, samp = dcis,same.dir = TRUE)
head(gse16873.go.p$greater[, 1:5],4)

# Significant genes
gse16873.go.sig<-sigGeneSet(gse16873.go.p, outname="gse16873.go")
str(gse16873.go.sig, strict.width='wrap')

gse16873.go.esg.up <- esset.grp(gse16873.go.p$greater,
                                gse16873, gsets = go.gs, ref = hn, samp = dcis,
                                test4up = T, output = T, outname = "gse16873.go.up", make.plot = T,heatmap=T)






