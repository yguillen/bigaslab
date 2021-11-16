## Expression profile of T-ALL human samples extracted from GEO ##
## modification of the script using GEOquery ##

source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("GEOquery")

#BiocManager::install("GEOquery")
BiocManager::install("hgu133plus2.db")
BiocManager::install("hthgu133pluspmprobe")
BiocManager::install("sva")

library(biomaRt)
library(Biobase)
library(affy)
library(ggrepel)
library(dplyr)
library(data.table)
library(GEOquery)
library(hgu133plus2.db)
library(hthgu133pluspmcdf)
library(hthgu133pluspmprobe)
library(genefilter)
library(sva)

# GSE62156
# GSE28703
# GSE14618
# GSE26713
# GSE8879
# E-MTAB-604
# GSE37389
# GSE56488
# GSE7615

# GSE33469
# GSE33470

# Download raw data CEL files for each project
#setwd("/Volumes/grcmc/YGUILLEN/beta_catenin_project/TALL_bcatenin_data/GEO_Expression/RAW_DATA") 
setwd("/Users/yguillen/Desktop/temp/beta_catenin_project/TALL_bcatenin_data/GEO_Expression/RAW_DATA") 
set.seed(100)

#### GSE62156 ###
#####
# Get phenotype data
GSE62156 <- getGEO("GSE62156", GSEMatrix = TRUE)
show(GSE62156)
pheno_GSE62156<-GSE62156$GSE62156_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE62156_raw = getGEOSuppFiles("GSE62156")
GSE62156_raw

## normalize CEL files (method old 3' microarrays) and output data
GSE62156_affy<-ReadAffy(celfile.path = "GSE62156/GSE62156_RAW/")
esetGSE62156 <- affy::rma(GSE62156_affy)

write.exprs(esetGSE62156,file="GSE62156/rma_GSE62156.txt")
exp_GSE62156<-data.frame(affy::exprs(esetGSE62156))
colnames(exp_GSE62156)<-gsub('_.*','',colnames(exp_GSE62156))

# Merge pheno data and expression
row.names(pheno_GSE62156)
dim(pheno_GSE62156)
table(pheno_GSE62156$`molecular subgroup:ch1`)

## ANNOTATE GENES
hgu133plus2()

Annot <- toTable(hgu133plus2SYMBOL)
colnames(Annot)
row.names(Annot)<-Annot$probe_id

expannot_GSE62156<-merge(Annot,exp_GSE62156,by.x=0,by.y=0,all.y=T)
row.names(expannot_GSE62156)<-expannot_GSE62156$Row.names
expannot_GSE62156$Row.names<-NULL


## PCA to remove outliers
pca_GSE62156<-expannot_GSE62156[,4:ncol(expannot_GSE62156)]
pca_GSE62156<-t(pca_GSE62156)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE62156, 2, var)!=0)
pca_GSE62156<-pca_GSE62156[, gsecond, drop = FALSE]


res.pca_GSE62156<-prcomp(pca_GSE62156,scale=TRUE)
(res.pca_GSE62156$sdev^2/sum(res.pca_GSE62156$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE62156<-as.data.frame(res.pca_GSE62156$x)
pcatab_GSE62156$names<-row.names(res.pca_GSE62156)
rm(res.pca_GSE62156)
row.names(pcatab_GSE62156)

pcatab_GSE62156<-pcatab_GSE62156[,c(1:3)]
row.names(pcatab_GSE62156)

row.names(pheno_GSE62156)
setdiff(row.names(pheno_GSE62156),row.names(pcatab_GSE62156))

pcatab_GSE62156<-merge(pheno_GSE62156,pcatab_GSE62156,by.x=0,by.y=0)
row.names(pcatab_GSE62156)<-pcatab_GSE62156$Row.names
pcatab_GSE62156$Row.names<-NULL

pcatab_GSE62156$source_name_ch1
pcatab_GSE62156$characteristics_ch1<-gsub('molecular subgroup: ','',pcatab_GSE62156$characteristics_ch1)
pcatab_GSE62156$characteristics_ch1<-gsub(' subgroup','',pcatab_GSE62156$characteristics_ch1)
pcatab_GSE62156$characteristics_ch1<-gsub(' subroup','',pcatab_GSE62156$characteristics_ch1)
pcatab_GSE62156$characteristics_ch1<-gsub('unknown.*','unknown',pcatab_GSE62156$characteristics_ch1)
table(pcatab_GSE62156$characteristics_ch1)

library(ggnewscale)

#Outliers converting pc1 to z-scores
out62146<-data.frame(geo_accession=pcatab_GSE62156$geo_accession,
           zscore_PC1=abs(scale(pcatab_GSE62156$PC1, center = TRUE, scale = TRUE)))

pcatab_GSE62156<-merge(pcatab_GSE62156,out62146,by="geo_accession")

gse62156_plot<-ggplot(pcatab_GSE62156, aes(PC1, PC2,label=geo_accession)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=characteristics_ch1)) +
  scale_color_brewer(palette = "Paired")+
  new_scale_color()+
#  stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
 # scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE62156[pcatab_GSE62156$zscore_PC1>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (6.5% variance)")+
  xlab("PC1 (10.0% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse62156_plot






#####

#### GSE28703 ###
#####
# Get phenotype data
GSE28703 <- getGEO("GSE28703", GSEMatrix = TRUE)
show(GSE28703)
pheno_GSE28703<-GSE28703$GSE28703_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE28703_raw = getGEOSuppFiles("GSE28703")
GSE28703_raw

## normalize CEL files (method old 3' microarrays) and output data
GSE28703_affy<-ReadAffy(celfile.path = "GSE28703/GSE28703_RAW/")
esetGSE28703 <- affy::rma(GSE28703_affy)

write.exprs(esetGSE28703,file="GSE28703/rma_GSE28703.txt")
exp_GSE28703<-data.frame(affy::exprs(esetGSE28703))
colnames(exp_GSE28703)<-gsub('_.*','',colnames(exp_GSE28703))

# Merge pheno data and expression
row.names(pheno_GSE28703)
dim(pheno_GSE28703)
table(pheno_GSE28703$source_name_ch1)

## ANNOTATE GENES
# from downloaded annotation GEO GPL13158
GPL13158<-read.delim("GSE28703/GPL13158-5065.txt")


Annot <- subset(GPL13158,select=c("ID","Gene.Symbol"))
colnames(Annot)<-c("probe_id","symbol")
row.names(Annot)<-Annot$probe_id

expannot_GSE28703<-merge(Annot,exp_GSE28703,by.x=0,by.y=0,all.y=T)
row.names(expannot_GSE28703)<-expannot_GSE28703$Row.names


## PCA to remove outliers
pca_GSE28703<-expannot_GSE28703[,4:ncol(expannot_GSE28703)]
pca_GSE28703<-t(pca_GSE28703)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE28703, 2, var)!=0)
pca_GSE28703<-pca_GSE28703[, gsecond, drop = FALSE]


res.pca_GSE28703<-prcomp(pca_GSE28703,scale=TRUE)
(res.pca_GSE28703$sdev^2/sum(res.pca_GSE28703$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE28703<-as.data.frame(res.pca_GSE28703$x)
pcatab_GSE28703$names<-row.names(res.pca_GSE28703)
rm(res.pca_GSE28703)
row.names(pcatab_GSE28703)

pcatab_GSE28703<-pcatab_GSE28703[,c(1:3)]
row.names(pcatab_GSE28703)

row.names(pheno_GSE28703)
setdiff(row.names(pheno_GSE28703),row.names(pcatab_GSE28703))

pcatab_GSE28703<-merge(pheno_GSE28703,pcatab_GSE28703,by.x=0,by.y=0)
row.names(pcatab_GSE28703)<-pcatab_GSE28703$Row.names
pcatab_GSE28703$Row.names<-NULL

pcatab_GSE28703$source_name_ch1<-gsub('.*early.*','ETP',pcatab_GSE28703$source_name_ch1)
pcatab_GSE28703$source_name_ch1<-gsub('.*non-ETP.*','non-ETP',pcatab_GSE28703$source_name_ch1)

table(pcatab_GSE28703$source_name_ch1)

#Outliers converting pc1 to z-scores
out28703<-data.frame(geo_accession=pcatab_GSE28703$geo_accession,
                     zscore_PC1=abs(scale(pcatab_GSE28703$PC1, center = TRUE, scale = TRUE)))

pcatab_GSE28703<-merge(pcatab_GSE28703,out28703,by="geo_accession")

library(ggnewscale)

gse28703_plot<-ggplot(pcatab_GSE28703, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=source_name_ch1)) +
  scale_color_brewer(palette = "Paired")+
  new_scale_color()+
  stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE28703[pcatab_GSE28703$zscore_PC1>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  ylab("PC2 (6.9% variance)")+
  xlab("PC1 (8.6% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse28703_plot



#####

#### GSE14618
#####
# Get phenotype data
GSE14618 <- getGEO("GSE14618", GSEMatrix = TRUE)
show(GSE14618)
pheno_GSE14618_gpl570<-GSE14618$`GSE14618-GPL570_series_matrix.txt.gz`@phenoData@data
pheno_GSE14618_gpl96<-GSE14618$`GSE14618-GPL96_series_matrix.txt.gz`@phenoData@data


# Download raw data CEL files for each project
GSE14618_raw = getGEOSuppFiles("GSE14618")
GSE14618_raw

# Divide manually in set 1 and set 2 within raw directory

## normalize CEL files (method old 3' microarrays) and output data
GSE14618_affy_1<-ReadAffy(celfile.path = "GSE14618/GSE14618_RAW/set1/")
esetGSE14618_1 <- affy::rma(GSE14618_affy_1)

GSE14618_affy_2<-ReadAffy(celfile.path = "GSE14618/GSE14618_RAW/set2/")
esetGSE14618_2 <- affy::rma(GSE14618_affy_2)


write.exprs(esetGSE14618_1,file="GSE62156/rma_GSE62156_1.txt")
exp_GSE14618_1<-data.frame(affy::exprs(esetGSE14618_1))
colnames(exp_GSE14618_1)<-gsub('\\..*','',colnames(exp_GSE14618_1))

write.exprs(esetGSE14618_2,file="GSE62156/rma_GSE62156_2.txt")
exp_GSE14618_2<-data.frame(affy::exprs(esetGSE14618_2))
colnames(exp_GSE14618_2)<-gsub('\\..*','',colnames(exp_GSE14618_2))


# Merge pheno data and expression
row.names(pheno_GSE14618_gpl570)
colnames(exp_GSE14618_2)

row.names(pheno_GSE14618_gpl96)
colnames(exp_GSE14618_1)


## ANNOTATE GENES
hgu133plus2()

Annot <- toTable(hgu133plus2SYMBOL)
colnames(Annot)
row.names(Annot)<-Annot$probe_id

expannot_GSE14618_2<-merge(Annot,exp_GSE14618_2,by.x=0,by.y=0,all.y=T)


# from downloaded annotation GEO GPL96
GPL96<-read.delim("../GPL96-57554.txt")


Annot <- subset(GPL96,select=c("ID","Gene.Symbol"))
colnames(Annot)<-c("probe_id","symbol")
row.names(Annot)<-Annot$probe_id

expannot_GSE14618_1<-merge(Annot,exp_GSE14618_1,by.x=0,by.y=0,all.y=T)


## PCA to remove outliers

row.names(expannot_GSE14618_1)<-expannot_GSE14618_1$Row.names

pca_GSE14618_1<-expannot_GSE14618_1[,4:ncol(expannot_GSE14618_1)]
pca_GSE14618_1<-t(pca_GSE14618_1)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE14618_1, 2, var)!=0)
pca_GSE14618_1<-pca_GSE14618_1[, gsecond, drop = FALSE]


res.pca_GSE14618_1<-prcomp(pca_GSE14618_1,scale=TRUE)
(res.pca_GSE14618_1$sdev^2/sum(res.pca_GSE14618_1$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE14618_1<-as.data.frame(res.pca_GSE14618_1$x)
pcatab_GSE14618_1$names<-row.names(res.pca_GSE14618_1)
rm(res.pca_GSE14618_1)
row.names(pcatab_GSE14618_1)

pcatab_GSE14618_1<-pcatab_GSE14618_1[,c(1:3)]
row.names(pcatab_GSE14618_1)

row.names(pheno_GSE14618_gpl96)
setdiff(row.names(pheno_GSE14618_gpl96),row.names(pcatab_GSE14618_1))

pcatab_GSE14618_1<-merge(pheno_GSE14618_gpl96,pcatab_GSE14618_1,by.x=0,by.y=0)
row.names(pcatab_GSE14618_1)<-pcatab_GSE14618_1$Row.names
pcatab_GSE14618_1$Row.names<-NULL

pcatab_GSE14618_1$source_name_ch1<-gsub('Primary T-ALL cells; newly diagnosed patient; ','',pcatab_GSE14618_1$source_name_ch1)

table(pcatab_GSE14618_1$source_name_ch1)


#Outliers converting pc1 to z-scores
out14618_1<-data.frame(geo_accession=pcatab_GSE14618_1$geo_accession,
                     zscore_PC1=abs(scale(pcatab_GSE14618_1$PC1, center = TRUE, scale = TRUE)))

pcatab_GSE14618_1<-merge(pcatab_GSE14618_1,out14618_1,by="geo_accession")

library(ggnewscale)

gse14618_1plot<-ggplot(pcatab_GSE14618_1,aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=source_name_ch1)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE14618_1[pcatab_GSE14618_1$zscore_PC1>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (9.5% variance)")+
  xlab("PC1 (25.6% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse14618_1plot


row.names(expannot_GSE14618_2)<-expannot_GSE14618_2$Row.names

pca_GSE14618_2<-expannot_GSE14618_2[,4:ncol(expannot_GSE14618_2)]
pca_GSE14618_2<-t(pca_GSE14618_2)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE14618_2, 2, var)!=0)
pca_GSE14618_2<-pca_GSE14618_2[, gsecond, drop = FALSE]


res.pca_GSE14618_2<-prcomp(pca_GSE14618_2,scale=TRUE)
(res.pca_GSE14618_2$sdev^2/sum(res.pca_GSE14618_2$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE14618_2<-as.data.frame(res.pca_GSE14618_2$x)
pcatab_GSE14618_2$names<-row.names(res.pca_GSE14618_2)
rm(res.pca_GSE14618_2)
row.names(pcatab_GSE14618_2)

pcatab_GSE14618_2<-pcatab_GSE14618_2[,c(1:3)]
row.names(pcatab_GSE14618_2)

row.names(pheno_GSE14618_gpl570)
setdiff(row.names(pheno_GSE14618_gpl570),row.names(pcatab_GSE14618_2))

pcatab_GSE14618_2<-merge(pheno_GSE14618_gpl570,pcatab_GSE14618_2,by.x=0,by.y=0)
row.names(pcatab_GSE14618_2)<-pcatab_GSE14618_2$Row.names
pcatab_GSE14618_2$Row.names<-NULL

pcatab_GSE14618_2$source_name_ch1<-gsub('Primary T-ALL cells; newly diagnosed patient; ','',pcatab_GSE14618_2$source_name_ch1)

table(pcatab_GSE14618_2$source_name_ch1)

#Outliers converting pc1 to z-scores
out14618_2<-data.frame(geo_accession=pcatab_GSE14618_2$geo_accession,
                       zscore_PC1=abs(scale(pcatab_GSE14618_2$PC1, center = TRUE, scale = TRUE)))

pcatab_GSE14618_2<-merge(pcatab_GSE14618_2,out14618_2,by="geo_accession")

library(ggnewscale)

gse14618_2plot<-ggplot(pcatab_GSE14618_2, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=source_name_ch1)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE14618_2[pcatab_GSE14618_2$zscore_PC1>=3,],aes(label=geo_accession),size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (7.1% variance)")+
  xlab("PC1 (22.4% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse14618_2plot



#####

#### GSE26713 ####
#####
# Get phenotype data
GSE26713 <- getGEO("GSE26713", GSEMatrix = TRUE)
show(GSE26713)
pheno_GSE26713<-GSE26713$GSE26713_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE26713_raw = getGEOSuppFiles("GSE26713")
GSE26713_raw

## normalize CEL files (method old 3' microarrays) and output data
GSE26713_affy<-ReadAffy(celfile.path = "GSE26713/GSE26713_RAW/")
esetGSE26713 <- affy::rma(GSE26713_affy)

write.exprs(esetGSE26713,file="GSE26713/rma_GSE26713.txt")
exp_GSE26713<-data.frame(affy::exprs(esetGSE26713))
colnames(exp_GSE26713)<-gsub('\\..*','',colnames(exp_GSE26713))

# Merge pheno data and expression
row.names(pheno_GSE26713)
dim(pheno_GSE26713)

## ANNOTATE GENES
hgu133plus2()

Annot <- toTable(hgu133plus2SYMBOL)
colnames(Annot)
row.names(Annot)<-Annot$probe_id

expannot_GSE26713<-merge(Annot,exp_GSE26713,by.x=0,by.y=0,all.y=T)

## PCA to remove outliers

row.names(expannot_GSE26713)<-expannot_GSE26713$Row.names

pca_GSE26713<-expannot_GSE26713[,4:ncol(expannot_GSE26713)]
pca_GSE26713<-t(pca_GSE26713)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE26713, 2, var)!=0)
pca_GSE26713<-pca_GSE26713[, gsecond, drop = FALSE]


res.pca_GSE26713<-prcomp(pca_GSE26713,scale=TRUE)
(res.pca_GSE26713$sdev^2/sum(res.pca_GSE26713$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE26713<-as.data.frame(res.pca_GSE26713$x)
pcatab_GSE26713$names<-row.names(res.pca_GSE26713)
rm(res.pca_GSE26713)
row.names(pcatab_GSE26713)

pcatab_GSE26713<-pcatab_GSE26713[,c(1:3)]
row.names(pcatab_GSE26713)

row.names(pheno_GSE26713)
setdiff(row.names(pheno_GSE26713),row.names(pcatab_GSE26713))

pcatab_GSE26713<-merge(pheno_GSE26713,pcatab_GSE26713,by.x=0,by.y=0)
row.names(pcatab_GSE26713)<-pcatab_GSE26713$Row.names
pcatab_GSE26713$Row.names<-NULL

pcatab_GSE26713$source_name_ch1<-gsub(' patient .*','',pcatab_GSE26713$source_name_ch1)
pcatab_GSE26713$source_name_ch1<-gsub('normal bone marrow control','normal BM',pcatab_GSE26713$source_name_ch1)

table(pcatab_GSE26713$`cytogenetics:ch1`)

#Outliers converting pc1 to z-scores
out26713<-data.frame(geo_accession=pcatab_GSE26713$geo_accession,
                       zscore_PC1=abs(scale(pcatab_GSE26713$PC1, center = TRUE, scale = TRUE)))

pcatab_GSE26713<-merge(pcatab_GSE26713,out26713,by="geo_accession")

library(ggnewscale)

gse26713_plot<-ggplot(pcatab_GSE26713, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=`cytogenetics:ch1`)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE26713[pcatab_GSE26713$zscore_PC1>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (10.4% variance)")+
  xlab("PC1 (16.7% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse26713_plot

# Remove outliers
outlier<-pcatab_GSE26713[pcatab_GSE26713$zscore_PC1>=3,]$geo_accession
expannot_GSE26713<-expannot_GSE26713[,!(colnames(expannot_GSE26713) %in% outlier)]

#####

#### GSE32215 ####
#####
# Get phenotype data
GSE32215 <- getGEO("GSE32215", GSEMatrix = TRUE)
show(GSE32215)
pheno_GSE32215<-GSE32215$GSE32215_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE32215_raw = getGEOSuppFiles("GSE32215")
GSE32215_raw

## normalize CEL files (method old 3' microarrays) and output data
GSE32215_affy<-ReadAffy(celfile.path = "GSE32215/GSE32215_RAW/")
esetGSE32215 <- affy::rma(GSE32215_affy)

write.exprs(esetGSE32215,file="GSE32215/rma_GSE32215.txt")
exp_GSE32215<-data.frame(affy::exprs(esetGSE32215))
colnames(exp_GSE32215)<-gsub('\\..*','',colnames(exp_GSE32215))

# Merge pheno data and expression
row.names(pheno_GSE32215)
dim(pheno_GSE32215)

## ANNOTATE GENES
hgu133plus2()

Annot <- toTable(hgu133plus2SYMBOL)
colnames(Annot)
row.names(Annot)<-Annot$probe_id

expannot_GSE32215<-merge(Annot,exp_GSE32215,by.x=0,by.y=0,all.y=T)

## PCA to remove outliers

row.names(expannot_GSE32215)<-expannot_GSE32215$Row.names
row.names(pheno_GSE32215)

# Remove reanalized samples from GSE26713
disc_rean<-pheno_GSE32215[!grepl('reanalysis',pheno_GSE32215$description),]$geo_accession
expannot_GSE32215<-expannot_GSE32215[,colnames(expannot_GSE32215) %in% disc_rean]

pca_GSE32215<-expannot_GSE32215[,1:ncol(expannot_GSE32215)]
pca_GSE32215<-t(pca_GSE32215)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE32215, 2, var)!=0)
pca_GSE32215<-pca_GSE32215[, gsecond, drop = FALSE]


res.pca_GSE32215<-prcomp(pca_GSE32215,scale=TRUE)
(res.pca_GSE32215$sdev^2/sum(res.pca_GSE32215$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE32215<-as.data.frame(res.pca_GSE32215$x)
pcatab_GSE32215$names<-row.names(res.pca_GSE32215)
rm(res.pca_GSE32215)
row.names(pcatab_GSE32215)

pcatab_GSE32215<-pcatab_GSE32215[,c(1:3)]
row.names(pcatab_GSE32215)

row.names(pheno_GSE32215)
setdiff(row.names(pheno_GSE32215),row.names(pcatab_GSE32215))

pcatab_GSE32215<-merge(pheno_GSE32215,pcatab_GSE32215,by.x=0,by.y=0)
row.names(pcatab_GSE32215)<-pcatab_GSE32215$Row.names
pcatab_GSE32215$Row.names<-NULL

pcatab_GSE32215$characteristics_ch1<-gsub('cell type: ','',pcatab_GSE32215$characteristics_ch1)

table(pcatab_GSE32215$characteristics_ch1)

#Outliers converting pc1 to z-scores
out32215<-data.frame(geo_accession=pcatab_GSE32215$geo_accession,
                     zscore_PC1=abs(scale(pcatab_GSE32215$PC1, center = TRUE, scale = TRUE)),
                     zscore_PC2=abs(scale(pcatab_GSE32215$PC2, center = TRUE, scale = TRUE)))

pcatab_GSE32215<-merge(pcatab_GSE32215,out32215,by="geo_accession")



library(ggnewscale)

gse32215_plot<-ggplot(pcatab_GSE32215, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=characteristics_ch1)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE32215[pcatab_GSE32215$zscore_PC1>=3 | pcatab_GSE32215$zscore_PC2>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (8.0% variance)")+
  xlab("PC1 (15.7% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse32215_plot

# Remove outliers
outlier<-pcatab_GSE32215[pcatab_GSE32215$zscore_PC1>=3,]$geo_accession
expannot_GSE322153<-expannot_GSE32215[,!(colnames(expannot_GSE32215) %in% outlier)]


#####


#### GSE8879 
#####
# Get phenotype data
GSE8879 <- getGEO("GSE8879", GSEMatrix = TRUE)
show(GSE8879)
pheno_GSE8879<-GSE8879$GSE8879_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE8879_raw = getGEOSuppFiles("GSE8879")
GSE8879_raw

## normalize CEL files (method old 3' microarrays) and output data
GSE8879_affy<-ReadAffy(celfile.path = "GSE8879/GSE8879_RAW/")
esetGSE8879 <- affy::rma(GSE8879_affy)

write.exprs(esetGSE8879,file="GSE8879/rma_GSE8879.txt")

exp_GSE8879<-data.frame(affy::exprs(esetGSE8879))
colnames(exp_GSE8879)<-gsub('\\..*','',colnames(exp_GSE8879))

# Merge pheno data and expression
row.names(pheno_GSE8879)
dim(pheno_GSE8879)
table(pheno_GSE8879$source_name_ch1)

## ANNOTATE GENES
# from downloaded annotation GEO GPL96
#GPL96

Annot <- subset(GPL96,select=c("ID","Gene.Symbol"))
colnames(Annot)<-c("probe_id","symbol")
row.names(Annot)<-Annot$probe_id

expannot_GSE8879<-merge(Annot,exp_GSE8879,by.x=0,by.y=0,all.y=T)

## PCA to remove outliers

row.names(expannot_GSE8879)<-expannot_GSE8879$Row.names

pca_GSE8879<-expannot_GSE8879[,4:ncol(expannot_GSE8879)]
pca_GSE8879<-t(pca_GSE8879)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE8879, 2, var)!=0)
pca_GSE8879<-pca_GSE8879[, gsecond, drop = FALSE]


res.pca_GSE8879<-prcomp(pca_GSE8879,scale=TRUE)
(res.pca_GSE8879$sdev^2/sum(res.pca_GSE8879$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE8879<-as.data.frame(res.pca_GSE8879$x)
pcatab_GSE8879$names<-row.names(res.pca_GSE8879)
rm(res.pca_GSE8879)
row.names(pcatab_GSE8879)

pcatab_GSE8879<-pcatab_GSE8879[,c(1:3)]
row.names(pcatab_GSE8879)

row.names(pheno_GSE8879)
setdiff(row.names(pheno_GSE8879),row.names(pcatab_GSE8879))

pcatab_GSE8879<-merge(pheno_GSE8879,pcatab_GSE8879,by.x=0,by.y=0)
row.names(pcatab_GSE8879)<-pcatab_GSE8879$Row.names
pcatab_GSE8879$Row.names<-NULL

pcatab_GSE8879$characteristics_ch1<-gsub('cell type: diagnostic leukemic blasts of early T-cell precursor acute lymphoblastic leukemia ','',pcatab_GSE8879$characteristics_ch1)
pcatab_GSE8879$characteristics_ch1<-gsub('cell type: diagnostic leukemic blasts of T-cell precursor acute lymphoblastic leukemia','non-ETP',pcatab_GSE8879$characteristics_ch1)
pcatab_GSE8879$characteristics_ch1<-gsub('\\(ETP\\)','ETP',pcatab_GSE8879$characteristics_ch1)

table(pcatab_GSE8879$characteristics_ch1)

#Outliers converting pc1 to z-scores
out8879<-data.frame(geo_accession=pcatab_GSE8879$geo_accession,
                     zscore_PC1=abs(scale(pcatab_GSE8879$PC1, center = TRUE, scale = TRUE)))

pcatab_GSE8879<-merge(pcatab_GSE8879,out8879,by="geo_accession")

library(ggnewscale)

gse8879_plot<-ggplot(pcatab_GSE8879, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=characteristics_ch1)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE8879[pcatab_GSE8879$zscore_PC1>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (8.9% variance)")+
  xlab("PC1 (15.0% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse8879_plot



#####


#### E-MTAB
#####

load("../E-MTAB-604/E-MTAB-604.eSet.r")


# Read phenotype data
pheno_emtab504<-read.delim("../E-MTAB-604/samples.txt",sep="\t")
row.names(pheno_emtab504)<-gsub('ALEK_','',pheno_emtab504$Sample.Name)
row.names(pheno_emtab504)<-gsub('_HS','',row.names(pheno_emtab504))

#Normalize expression data
esetemtab604 <- affy::rma(study)
write.exprs(esetemtab604,file="../E-MTAB-604/data.txt")
colnames(esetemtab604)<-gsub('_U133_2','',colnames(esetemtab604))

class(esetemtab604)

write.exprs(esetemtab604,file="../E-MTAB-604/rma_E-MTAB-604.txt")
exp_emtab604<-data.frame(affy::exprs(esetemtab604))
colnames(exp_emtab604)

# Merge pheno data and expression
row.names(pheno_emtab504)
dim(pheno_emtab504)
table(pheno_emtab504$Characteristics..DiseaseState.)

## ANNOTATE GENES
hgu133plus2()

Annot <- toTable(hgu133plus2SYMBOL)
colnames(Annot)
row.names(Annot)<-Annot$probe_id

expannot_emtab604<-merge(Annot,exp_emtab604,by.x=0,by.y=0,all.y=T)

## PCA to remove outliers

row.names(expannot_emtab604)<-expannot_emtab604$Row.names

pca_emtab604<-expannot_emtab604[,4:ncol(expannot_emtab604)]
pca_emtab604<-t(pca_emtab604)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_emtab604, 2, var)!=0)
pca_emtab604<-pca_emtab604[, gsecond, drop = FALSE]


res.pca_emtab604<-prcomp(pca_emtab604,scale=TRUE)
(res.pca_emtab604$sdev^2/sum(res.pca_emtab604$sdev^2))[1:5]


#Proportion of variance

pcatab_emtab604<-as.data.frame(res.pca_emtab604$x)
pcatab_emtab604$names<-row.names(res.pca_emtab604)
rm(res.pca_emtab604)
row.names(pcatab_emtab604)

pcatab_emtab604<-pcatab_emtab604[,c(1:3)]
row.names(pcatab_emtab604)

row.names(pheno_emtab504)
setdiff(row.names(pheno_emtab504),row.names(pcatab_emtab604))

pcatab_emtab604<-merge(pheno_emtab504,pcatab_emtab604,by.x=0,by.y=0)
row.names(pcatab_emtab604)<-pcatab_emtab604$Row.names
pcatab_emtab604$Row.names<-NULL

pcatab_emtab604$Characteristics..DiseaseState.

table(pcatab_emtab604$Characteristics..DiseaseState.)
table(pcatab_emtab604$Characteristics..OrganismPart.)
table(pcatab_emtab604$Characteristics..ClinicalInformation.)
table(pcatab_emtab604$Characteristics..CellLine.)


pcatab_emtab604$Characteristics..ClinicalInformation.<-gsub(' ALEK_.*','',pcatab_emtab604$Characteristics..ClinicalInformation.)
pcatab_emtab604$Characteristics..ClinicalInformation.<-gsub('T-ALL *','T-ALL',pcatab_emtab604$Characteristics..ClinicalInformation.)


#Outliers converting pc1 to z-scores
outemtab604<-data.frame(Sample.Name=pcatab_emtab604$Sample.Name,
                    zscore_PC1=abs(scale(pcatab_emtab604$PC1, center = TRUE, scale = TRUE)))

pcatab_emtab604<-merge(pcatab_emtab604,outemtab604,by="Sample.Name")

library(ggnewscale)

emtab604_plot<-ggplot(pcatab_emtab604, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=Characteristics..ClinicalInformation.)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_emtab604[pcatab_emtab604$zscore_PC1>=3,],aes(label=Sample.Name), size=2,fontface = "italic",alpha=0.7)+
  geom_label_repel(data=pcatab_emtab604[pcatab_emtab604$Factor.Value.DiseaseState.!="T-cell Acute Lymphoblastic Leukemia",],aes(label=Characteristics..CellLine.), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (8.5% variance)")+
  xlab("PC1 (9.9% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

emtab604_plot

#####

#### GSE33469 
#####
# Get phenotype data
GSE33469 <- getGEO("GSE33469", GSEMatrix = TRUE)
show(GSE33469)
pheno_GSE33469<-GSE33469$GSE33469_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE33469_raw = getGEOSuppFiles("GSE33469")
GSE33469_raw

## normalize (method old 3' microarrays) and output data
GSE33469_affy<-read.delim("GSE33469/GSE33469_non-normalized.txt")
row.names(GSE33469_affy)<-GSE33469_affy$ID_REF
GSE33469_affy$ID_REF<-NULL

esetGSE33469<-ExpressionSet(as.matrix(GSE33469_affy))
affy::rma(esetGSE33469)


frm#Normalize expression data
exp_GSE33469<-data.frame(affy::exprs(esetGSE33469))
colnames(exp_GSE33469)<-row.names(pheno_GSE33469)

# Merge pheno data and expression
row.names(pheno_GSE33469)
colnames(exp_GSE33469)


dim(pheno_GSE33469)


## ANNOTATE GENES
# from downloaded annotation GEO GPL10558
#GPL10558
GPL10558<-read.delim("../GPL10558-50081.txt")

Annot<-GPL10558
Annot<-subset(Annot,select=c("ID","Symbol"))
colnames(Annot)<-c("probe_id","symbol")
Annot<-Annot[!(is.na(Annot$probe_id) | Annot$probe_id==""), ]

Annot<-Annot[!duplicated(Annot$probe_id),]

row.names(Annot)<-Annot$probe_id

row.names(exp_GSE33469)
expannot_GSE33469<-merge(Annot,exp_GSE33469,by.x=0,by.y=0,all.y=T)

## PCA to remove outliers

row.names(expannot_GSE33469)<-expannot_GSE33469$Row.names

pca_GSE33469<-expannot_GSE33469[,4:ncol(expannot_GSE33469)]
pca_GSE33469<-t(pca_GSE33469)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE33469, 2, var)!=0)
pca_GSE33469<-pca_GSE33469[, gsecond, drop = FALSE]


res.pca_GSE33469<-prcomp(pca_GSE33469,scale=TRUE)
(res.pca_GSE33469$sdev^2/sum(res.pca_GSE33469$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE33469<-as.data.frame(res.pca_GSE33469$x)
pcatab_GSE33469$names<-row.names(res.pca_GSE33469)
rm(res.pca_GSE33469)
row.names(pcatab_GSE33469)

pcatab_GS33469<-pcatab_GSE33469[,c(1:3)]
row.names(pcatab_GSE33469)

row.names(pheno_GSE33469)
setdiff(row.names(pheno_GSE33469),row.names(pcatab_GSE33469))

pcatab_GSE33469<-merge(pheno_GSE33469,pcatab_GSE33469,by.x=0,by.y=0)
row.names(pcatab_GSE33469)<-pcatab_GSE33469$Row.names
pcatab_GSE33469$Row.names<-NULL



#Outliers converting pc1 to z-scores
out33469<-data.frame(geo_accession=pcatab_GSE33469$geo_accession,
                    zscore_PC1=abs(scale(pcatab_GSE33469$PC1, center = TRUE, scale = TRUE)))

pcatab_GSE33469<-merge(pcatab_GSE33469,out33469,by="geo_accession")

library(ggnewscale)

gse33469_plot<-ggplot(pcatab_GSE33469, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=characteristics_ch1)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE33469[pcatab_GSE33469$zscore_PC1>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (4.8% variance)")+
  xlab("PC1 (27.5% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse33469_plot

#####

#### GSE33470 
#####
# Get phenotype data
GSE33470 <- getGEO("GSE33470", GSEMatrix = TRUE)
show(GSE33470)
pheno_GSE33470<-GSE33470$GSE33470_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE33470_raw = getGEOSuppFiles("GSE33470")
GSE33470_raw

## normalize (method old 3' microarrays) and output data
GSE33470_affy<-read.delim("GSE33470/GSE33470_non-normalized.txt")
row.names(GSE33470_affy)<-GSE33470_affy$ID_REF

GSE33470_affy$ID_REF<-NULL
esetGSE33470<-ExpressionSet(as.matrix(GSE33470_affy))

#Normalize expression data
exp_GSE33470<-data.frame(affy::exprs(esetGSE33470))
colnames(exp_GSE33470)<-row.names(pheno_GSE33470)

# Merge pheno data and expression
row.names(pheno_GSE33470)
colnames(exp_GSE33470)


dim(pheno_GSE33470)
table(pheno_GSE33470$source_name_ch1)

## ANNOTATE GENES
# from downloaded annotation GEO GPL10558
#GPL10558
GPL10558<-read.delim("../GPL10558-50081.txt")

Annot<-GPL10558
Annot<-subset(Annot,select=c("ID","Symbol"))
colnames(Annot)<-c("probe_id","symbol")
Annot<-Annot[!(is.na(Annot$probe_id) | Annot$probe_id==""), ]

Annot<-Annot[!duplicated(Annot$probe_id),]

row.names(Annot)<-Annot$probe_id

row.names(exp_GSE33470)
expannot_GSE33470<-merge(Annot,exp_GSE33470,by.x=0,by.y=0,all.y=T)

## PCA to remove outliers

row.names(expannot_GSE33470)<-expannot_GSE33470$Row.names
row.names(pheno_GSE33470)

pca_GSE33470<-expannot_GSE33470[,4:ncol(expannot_GSE33470)]
pca_GSE33470<-t(pca_GSE33470)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE33470, 2, var)!=0)
pca_GSE33470<-pca_GSE33470[, gsecond, drop = FALSE]


res.pca_GSE33470<-prcomp(pca_GSE33470,scale=TRUE)
(res.pca_GSE33470$sdev^2/sum(res.pca_GSE33470$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE33470<-as.data.frame(res.pca_GSE33470$x)
pcatab_GSE33470$names<-row.names(res.pca_GSE33470)
rm(res.pca_GSE33470)
row.names(pcatab_GSE33470)

pcatab_GSE33470<-pcatab_GSE33470[,c(1:3)]
row.names(pcatab_GSE33470)

row.names(pheno_GSE33470)
setdiff(row.names(pheno_GSE33470),row.names(pcatab_GSE33470))

pcatab_GSE33470<-merge(pheno_GSE33470,pcatab_GSE33470,by.x=0,by.y=0)
row.names(pcatab_GSE33470)<-pcatab_GSE33470$Row.names
pcatab_GSE33470$Row.names<-NULL

pcatab_GSE33470$source_name_ch1<-gsub('T-Cell population ','',pcatab_GSE33470$source_name_ch1)

table(pcatab_GSE33470$source_name_ch1)

#Outliers converting pc1 to z-scores
out33470<-data.frame(geo_accession=pcatab_GSE33470$geo_accession,
                     zscore_PC1=abs(scale(pcatab_GSE33470$PC1, center = TRUE, scale = TRUE)),
                     zscore_PC2=abs(scale(pcatab_GSE33470$PC2, center = TRUE, scale = TRUE)))

pcatab_GSE33470<-merge(pcatab_GSE33470,out33470,by="geo_accession")



library(ggnewscale)

gse33470_plot<-ggplot(pcatab_GSE33470, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=source_name_ch1)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE33470[pcatab_GSE33470$zscore_PC1>=3 | pcatab_GSE33470$zscore_PC2>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (16.8% variance)")+
  xlab("PC1 (54.77% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse33470_plot

# Remove outliers
outlier<-pcatab_GSE33470[pcatab_GSE33470$zscore_PC1>=3,]$geo_accession
expannot_GSE33470<-expannot_GSE33470[,!(colnames(expannot_GSE33470) %in% outlier)]


#####

#### GSE10609 ####
#####
# Get phenotype data
GSE10609 <- getGEO("GSE10609", GSEMatrix = TRUE)
show(GSE10609)
pheno_GSE10609<-GSE10609$GSE10609_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE10609_raw = getGEOSuppFiles("GSE10609")
GSE10609_raw

## normalize CEL files (method old 3' microarrays) and output data
GSE10609_affy<-ReadAffy(celfile.path = "GSE10609/GSE10609_RAW/")
esetGSE10609 <- affy::rma(GSE10609_affy)

write.exprs(esetGSE10609,file="GSE10609/rma_GSE10609.txt")
exp_GSE10609<-data.frame(affy::exprs(esetGSE10609))
colnames(exp_GSE10609)<-gsub('\\..*','',colnames(exp_GSE10609))

# Merge pheno data and expression
row.names(pheno_GSE10609)
dim(pheno_GSE10609)

## ANNOTATE GENES
hgu133plus2()

Annot <- toTable(hgu133plus2SYMBOL)
colnames(Annot)
row.names(Annot)<-Annot$probe_id

expannot_GSE10609<-merge(Annot,exp_GSE10609,by.x=0,by.y=0,all.y=T)

## PCA to remove outliers

row.names(expannot_GSE10609)<-expannot_GSE10609$Row.names

pca_GSE10609<-expannot_GSE10609[,4:ncol(expannot_GSE10609)]
pca_GSE10609<-t(pca_GSE10609)

#select those with variance != 0 (error in prcomp otherwise)
gsecond<-(apply(pca_GSE10609, 2, var)!=0)
pca_GSE10609<-pca_GSE10609[, gsecond, drop = FALSE]


res.pca_GSE10609<-prcomp(pca_GSE10609,scale=TRUE)
(res.pca_GSE10609$sdev^2/sum(res.pca_GSE10609$sdev^2))[1:5]


#Proportion of variance

pcatab_GSE10609<-as.data.frame(res.pca_GSE10609$x)
pcatab_GSE10609$names<-row.names(res.pca_GSE10609)
rm(res.pca_GSE10609)
row.names(pcatab_GSE10609)

pcatab_GSE10609<-pcatab_GSE10609[,c(1:3)]
row.names(pcatab_GSE10609)

row.names(pheno_GSE10609)
setdiff(row.names(pheno_GSE10609),row.names(pcatab_GSE10609))

pcatab_GSE10609<-merge(pheno_GSE10609,pcatab_GSE10609,by.x=0,by.y=0)
row.names(pcatab_GSE10609)<-pcatab_GSE10609$Row.names
pcatab_GSE10609$Row.names<-NULL

pcatab_GSE10609$source_name_ch1<-gsub('pediatric T-ALL patient .*, ','',pcatab_GSE10609$source_name_ch1)

table(pcatab_GSE10609$source_name_ch1)

#Outliers converting pc1 to z-scores
out10609<-data.frame(geo_accession=pcatab_GSE10609$geo_accession,
                     zscore_PC1=abs(scale(pcatab_GSE10609$PC1, center = TRUE, scale = TRUE)))

pcatab_GSE10609<-merge(pcatab_GSE10609,out10609,by="geo_accession")

library(ggnewscale)

gse10609_plot<-ggplot(pcatab_GSE10609, aes(PC1, PC2)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=source_name_ch1)) +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=source_name_ch1,linetype=source_name_ch1),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_GSE10609[pcatab_GSE10609$zscore_PC1>=3,],aes(label=geo_accession), size=2,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  ylab("PC2 (11.7% variance)")+
  xlab("PC1 (18.1% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

gse10609_plot



#####


#### GSE37389 
#####
# Get phenotype data
GSE37389 <- getGEO("GSE37389", GSEMatrix = TRUE)
show(GSE37389)
pheno_GSE37389<-GSE37389$GSE37389_series_matrix.txt.gz@phenoData@data


# Download raw data CEL files for each project
GSE37389_raw = getGEOSuppFiles("GSE37389")
GSE37389_raw

## Agilent 4 data
files=list.files(path="GSE37389/GSE37389_RAW/")
GSE37389_affy<-read.maimages(files=files,path="GSE37389/GSE37389_RAW",source="agilent",green.only=TRUE,other.columns="gIsWellAboveBG")

GSE37389_affy$targets

# Correct and normalize
GSE37389_norm <- backgroundCorrect(GSE37389_affy,method="normexp")

GSE37389_norm <- normalizeBetweenArrays(GSE37389_norm,method="quantile")

head(GSE37389_norm$genes$ProbeName)
rownames(GSE37389_norm) <-GSE37389_norm$genes$ProbeName

Control <- GSE37389_norm$genes$ControlType==1L
IsExpr <- rowSums(GSE37389_norm$other$gIsWellAboveBG > 0) >= 59

GSE37389_norm <- GSE37389_norm[!Control & IsExpr, ]


#####

#### GSE56488 
#####
# Get phenotype data
GSE56488 <- getGEO("GSE56488", GSEMatrix = TRUE)
show(GSE56488)
pheno_GSE56488<-GSE56488$GSE56488_series_matrix.txt.gz@phenoData@data

# Download raw data CEL files for each project
GSE56488_raw = getGEOSuppFiles("GSE56488")
GSE56488_raw

## normalize with oligo package, affy does not support this platform
celFiles <- list.celfiles("GSE56488/GSE56488_RAW/", full.name=TRUE)
GSE56488_affy<-oligo::read.celfiles(celFiles)



#####



### CORRELATION PROBES GENES ####

### Step 1 Select a reference data set Dref
head(expannot_GSE26713)
dim(expannot_GSE26713)
# 54675 probes (including NAs)
dim(expannot_GSE26713[!is.na(expannot_GSE26713$symbol),])
# 41935 probes (excluding NAs)

## Select common gene set: genes represented in all datasets

commongenes<-Reduce(intersect,list(expannot_GSE62156[!is.na(expannot_GSE62156$symbol),]$symbol,
                      expannot_GSE28703[!is.na(expannot_GSE28703$symbol),]$symbol,
                      expannot_GSE14618_1[!is.na(expannot_GSE14618_1$symbol),]$symbol,
                      expannot_GSE14618_2[!is.na(expannot_GSE14618_2$symbol),]$symbol,
                      expannot_GSE26713[!is.na(expannot_GSE26713$symbol),]$symbol,
                      expannot_GSE8879[!is.na(expannot_GSE8879$symbol),]$symbol,
                      expannot_emtab604[!is.na(expannot_emtab604$symbol),]$symbol))

commongenes
length(commongenes)

# Select common genes from gene set reference
mad_GSE62156<-expannot_GSE62156[!is.na(expannot_GSE62156$symbol) & expannot_GSE62156$symbol %in% commongenes,]
mad_GSE28703<-expannot_GSE28703[!is.na(expannot_GSE28703$symbol) & expannot_GSE28703$symbol %in% commongenes,]
mad_GSE14618_1<-expannot_GSE14618_1[!is.na(expannot_GSE14618_1$symbol) & expannot_GSE14618_1$symbol %in% commongenes,]
mad_GSE14618_2<-expannot_GSE14618_2[!is.na(expannot_GSE14618_2$symbol) & expannot_GSE14618_2$symbol %in% commongenes,]
mad_GSE26713<-expannot_GSE26713[!is.na(expannot_GSE26713$symbol) & expannot_GSE26713$symbol %in% commongenes,]
mad_GSE8879<-expannot_GSE8879[!is.na(expannot_GSE8879$symbol) & expannot_GSE8879$symbol %in% commongenes,]
mad_emtab604<-expannot_emtab604[!is.na(expannot_emtab604$symbol) & expannot_emtab604$symbol %in% commongenes,]


### Reference dataset 5000 genes with probeset id largest MAD

mad_GSE26713$MAD <- apply(mad_GSE26713[,4:ncol(mad_GSE26713)], 1, mad, na.rm = TRUE)
max_GSE26713<-aggregate(MAD~symbol, mad_GSE26713, max)
max_GSE26713<-merge(mad_GSE26713, max_GSE26713)
grefDref<-(max_GSE26713[rev(order(max_GSE26713$MAD)),][1:5000,])
grefDref<-grefDref[,c(1,5:ncol(grefDref))]
grefDref$symbol
#####


# Estimate MAD per each dataset and select Genes reference for each dataset
### GSE62156 ####

# 1. select expression values for all probes corresponding to reference genes in grefDref and extract probe ids for these genes
# with highest MAD values
mad_GSE62156$MAD <- apply(mad_GSE62156[,4:ncol(mad_GSE62156)], 1, mad, na.rm = TRUE)
grefGSE62156<-mad_GSE62156[mad_GSE62156$symbol %in% grefDref$symbol,]
grefGSE62156<-aggregate(MAD~symbol, grefGSE62156, max)
grefGSE62156<-merge(mad_GSE62156, grefGSE62156)

# 2. select expression values for the rest of genes, all probe ids
mad_GSE62156<-mad_GSE62156[!(mad_GSE62156$symbol %in% grefGSE62156$symbol),]

# Do correlation for each of the probe in mod_GSE62156 with each of the genes in grefGSE62156
dim(mad_GSE62156)
tmad_GSE62156<-t(mad_GSE62156[,4:(ncol(mad_GSE62156)-1)])
dim(tmad_GSE62156)

dim(grefGSE62156)
row.names(grefGSE62156)<-grefGSE62156$symbol
tgrefGSE62156<-t(grefGSE62156[,5:ncol(grefGSE62156)])
dim(tgrefGSE62156)

matGSE62156<-merge(tgrefGSE62156,tmad_GSE62156,by=0)
dim(matGSE62156)
row.names(matGSE62156)<-matGSE62156$Row.names
matGSE62156$Row.names<-NULL

# Create matrix of correlations, and select correlation of each probe (,5001:ncol(corGSE62156)) with reference gene (1:5000,)
corGSE62156<-cor(matGSE62156)
corGSE62156<-corGSE62156[1:5000,5001:ncol(matGSE62156)]
corGSE62156[1:5,1:5]
dim(corGSE62156)

#####

### GSE28703 ####

# 1. select expression values for all probes corresponding to reference genes in grefDref and extract probe ids for these genes
# with highest MAD values
mad_GSE28703$MAD <- apply(mad_GSE28703[,4:ncol(mad_GSE28703)], 1, mad, na.rm = TRUE)
grefGSE28703<-mad_GSE28703[mad_GSE28703$symbol %in% grefDref$symbol,]
grefGSE28703<-aggregate(MAD~symbol, grefGSE28703, max)
grefGSE28703<-merge(mad_GSE28703, grefGSE28703)

# 2. select expression values for the rest of genes, all probe ids
mad_GSE28703<-mad_GSE28703[!(mad_GSE28703$symbol %in% grefGSE28703$symbol),]

# Do correlation for each of the probe in mod_GSE62156 with each of the genes in grefGSE28703
dim(mad_GSE28703)
tmad_GSE28703<-t(mad_GSE28703[,4:(ncol(mad_GSE28703)-1)])
dim(grefGSE28703)
row.names(grefGSE28703)<-grefGSE28703$symbol
tgrefGSE28703<-t(grefGSE28703[,5:ncol(grefGSE28703)])

matGSE28703<-merge(tgrefGSE28703,tmad_GSE28703,by=0)
row.names(matGSE28703)<-matGSE28703$Row.names
matGSE28703$Row.names<-NULL

# Create matrix of correlations, and select correlation of each probe (,5001:ncol(corGSE62156)) with reference gene (1:5000,)
corGSE28703<-cor(matGSE28703)
corGSE28703<-corGSE28703[1:5000,5001:ncol(corGSE28703)]

#####

### GSE14618_1 ####

# 1. select expression values for all probes corresponding to reference genes in grefDref and extract probe ids for these genes
# with highest MAD values
mad_GSE14618_1$MAD <- apply(mad_GSE14618_1[,4:ncol(mad_GSE14618_1)], 1, mad, na.rm = TRUE)
grefGSE14618_1<-mad_GSE14618_1[mad_GSE14618_1$symbol %in% grefDref$symbol,]
grefGSE14618_1<-aggregate(MAD~symbol, grefGSE14618_1, max)
grefGSE14618_1<-merge(mad_GSE14618_1, grefGSE14618_1)

# 2. select expression values for the rest of genes, all probe ids
mad_GSE14618_1<-mad_GSE14618_1[!(mad_GSE14618_1$symbol %in% grefGSE14618_1$symbol),]

# Do correlation for each of the probe in mod_GSE62156 with each of the genes in grefGSE28703
dim(mad_GSE14618_1)
tmad_GSE14618_1<-t(mad_GSE14618_1[,4:(ncol(mad_GSE14618_1)-1)])
dim(grefGSE14618_1)
row.names(grefGSE14618_1)<-grefGSE14618_1$symbol
tgrefGSE14618_1<-t(grefGSE14618_1[,5:ncol(grefGSE14618_1)])

matGSE14618_1<-merge(tgrefGSE14618_1,tmad_GSE14618_1,by=0)
row.names(matGSE14618_1)<-matGSE14618_1$Row.names
matGSE14618_1$Row.names<-NULL

# Create matrix of correlations, and select correlation of each probe (,5001:ncol(corGSE62156)) with reference gene (1:5000,)
corGSE14618_1<-cor(matGSE14618_1)
corGSE14618_1<-corGSE14618_1[1:5000,5001:ncol(corGSE14618_1)]

#####

### GSE14618_2 ####

# 1. select expression values for all probes corresponding to reference genes in grefDref and extract probe ids for these genes
# with highest MAD values
mad_GSE14618_2$MAD <- apply(mad_GSE14618_2[,4:ncol(mad_GSE14618_2)], 1, mad, na.rm = TRUE)
grefGSE14618_2<-mad_GSE14618_2[mad_GSE14618_2$symbol %in% grefDref$symbol,]
grefGSE14618_2<-aggregate(MAD~symbol, grefGSE14618_2, max)
grefGSE14618_2<-merge(mad_GSE14618_2, grefGSE14618_2)

# 2. select expression values for the rest of genes, all probe ids
mad_GSE14618_2<-mad_GSE14618_2[!(mad_GSE14618_2$symbol %in% grefGSE14618_2$symbol),]

# Do correlation for each of the probe in mod_GSE62156 with each of the genes in grefGSE28703
dim(mad_GSE14618_2)
tmad_GSE14618_2<-t(mad_GSE14618_2[,4:(ncol(mad_GSE14618_2)-1)])
dim(grefGSE14618_2)
row.names(grefGSE14618_2)<-grefGSE14618_2$symbol
tgrefGSE14618_2<-t(grefGSE14618_2[,5:ncol(grefGSE14618_2)])

matGSE14618_2<-merge(tgrefGSE14618_2,tmad_GSE14618_2,by=0)
row.names(matGSE14618_2)<-matGSE14618_2$Row.names
matGSE14618_2$Row.names<-NULL

# Create matrix of correlations, and select correlation of each probe (,5001:ncol(corGSE62156)) with reference gene (1:5000,)
corGSE14618_2<-cor(matGSE14618_2)
corGSE14618_2<-corGSE14618_2[1:5000,5001:ncol(corGSE14618_2)]

#####

### GSE26713 Reference dataset ####

# 1. select expression values for all probes corresponding to reference genes in grefDref and extract probe ids for these genes
# with highest MAD values
mad_GSE26713$MAD <- apply(mad_GSE26713[,4:ncol(mad_GSE26713)], 1, mad, na.rm = TRUE)
grefGSE26713<-mad_GSE26713[mad_GSE26713$symbol %in% grefDref$symbol,]
grefGSE26713<-aggregate(MAD~symbol, grefGSE26713, max)
grefGSE26713<-merge(mad_GSE26713, grefGSE26713)

# 2. select expression values for the rest of genes, all probe ids
mad_GSE26713<-mad_GSE26713[!(mad_GSE26713$symbol %in% grefGSE26713$symbol),]

# Do correlation for each of the probe in mod_GSE62156 with each of the genes in grefGSE28703
dim(mad_GSE26713)
tmad_GSE26713<-t(mad_GSE26713[,4:(ncol(mad_GSE26713)-1)])
dim(grefGSE26713)
row.names(grefGSE26713)<-grefGSE26713$symbol
tgrefGSE26713<-t(grefGSE26713[,5:ncol(grefGSE26713)])

matGSE26713<-merge(tgrefGSE26713,tmad_GSE26713,by=0)
row.names(matGSE26713)<-matGSE26713$Row.names
matGSE26713$Row.names<-NULL

# Create matrix of correlations, and select correlation of each probe (,5001:ncol(corGSE62156)) with reference gene (1:5000,)
corGSE26713<-cor(matGSE26713)
corGSE26713<-corGSE26713[1:5000,5001:ncol(corGSE26713)]

#####

### GSE8879 ####

# 1. select expression values for all probes corresponding to reference genes in grefDref and extract probe ids for these genes
# with highest MAD values
mad_GSE8879$MAD <- apply(mad_GSE8879[,4:ncol(mad_GSE8879)], 1, mad, na.rm = TRUE)
grefGSE8879<-mad_GSE8879[mad_GSE8879$symbol %in% grefDref$symbol,]
grefGSE8879<-aggregate(MAD~symbol, grefGSE8879, max)
grefGSE8879<-merge(mad_GSE8879, grefGSE8879)

# 2. select expression values for the rest of genes, all probe ids
mad_GSE8879<-mad_GSE8879[!(mad_GSE8879$symbol %in% grefGSE8879$symbol),]

# Do correlation for each of the probe in mod_GSE62156 with each of the genes in grefGSE28703
dim(mad_GSE8879)
tmad_GSE8879<-t(mad_GSE8879[,4:(ncol(mad_GSE8879)-1)])
dim(grefGSE8879)
row.names(grefGSE8879)<-grefGSE8879$symbol
tgrefGSE8879<-t(grefGSE8879[,5:ncol(grefGSE8879)])

matGSE8879<-merge(tgrefGSE8879,tmad_GSE8879,by=0)
row.names(matGSE8879)<-matGSE8879$Row.names
matGSE8879$Row.names<-NULL

# Create matrix of correlations, and select correlation of each probe (,5001:ncol(corGSE62156)) with reference gene (1:5000,)
corGSE8879<-cor(matGSE8879)
corGSE8879<-corGSE8879[1:5000,5001:ncol(corGSE8879)]

#####

### emtab604 ####

# 1. select expression values for all probes corresponding to reference genes in grefDref and extract probe ids for these genes
# with highest MAD values
mad_emtab604$MAD <- apply(mad_emtab604[,4:ncol(mad_emtab604)], 1, mad, na.rm = TRUE)
grefemtab604<-mad_emtab604[mad_emtab604$symbol %in% grefDref$symbol,]
grefemtab604<-aggregate(MAD~symbol, grefemtab604, max)
grefemtab604<-merge(mad_emtab604, grefemtab604)

# 2. select expression values for the rest of genes, all probe ids
mad_emtab604<-mad_emtab604[!(mad_emtab604$symbol %in% grefemtab604$symbol),]

# Do correlation for each of the probe in mod_GSE62156 with each of the genes in grefGSE28703
dim(mad_emtab604)
tmad_emtab604<-t(mad_emtab604[,4:(ncol(mad_emtab604)-1)])
dim(grefemtab604)
row.names(grefemtab604)<-grefemtab604$symbol
tgrefemtab604<-t(grefemtab604[,5:ncol(grefemtab604)])

matemtab604<-merge(tgrefemtab604,tmad_emtab604,by=0)
row.names(matemtab604)<-matemtab604$Row.names
matemtab604$Row.names<-NULL

# Create matrix of correlations, and select correlation of each probe (,5001:ncol(corGSE62156)) with reference gene (1:5000,)
coremtab604<-cor(matemtab604)
coremtab604<-coremtab604[1:5000,5001:ncol(coremtab604)]

#####

### Correlation ####

corGSE62156_dref<-as.matrix(merge(corGSE26713,corGSE62156,by=0))
row.names(corGSE62156_dref)<-corGSE62156_dref$Row.names
corGSE62156_dref$Row.names<-NULL
colnames(corGSE62156_dref)<-gsub('\\.x','_ref',colnames(corGSE62156_dref))
colnames(corGSE62156_dref)<-gsub('\\.y','_GSE62156',colnames(corGSE62156_dref))

#only corresponding positions in dataframes
testcor<-outer(1:ncol(corGSE26713),(ncol(corGSE26713)+1):ncol(corGSE62156_dref), function(x,y) cor.test(corGSE62156_dref[,x], corGSE62156_dref[,y]))

testcor<-outer(1:10,11, function(x,y) cor.test(corGSE62156_dref[,x],corGSE62156_dref[,y]))


## Example with one gene, example PAX8 (no ref genes)
# probeset ids for PAX8 in GSE26713 reference dataset
row.names(mad_GSE26713[mad_GSE26713$symbol=="PAX8",])
test1<-corGSE26713[,colnames(corGSE26713) %in% row.names(mad_GSE26713[mad_GSE26713$symbol=="PAX8",])]
colnames(test1)<-gsub('^','ref',colnames(test1))

# probeset ids for PAX8 in the rest of datasets
row.names(mad_GSE62156[mad_GSE62156$symbol=="PAX8",])
test2<-corGSE62156[,colnames(corGSE62156) %in% row.names(mad_GSE62156[mad_GSE62156$symbol=="PAX8",])]
dim(test2)

row.names(mad_GSE28703[mad_GSE28703$symbol=="PAX8",])
test3<-corGSE28703[,colnames(corGSE28703) %in% row.names(mad_GSE28703[mad_GSE28703$symbol=="PAX8",])]
dim(test3)

row.names(mad_GSE14618_1[mad_GSE14618_1$symbol=="PAX8",])
test4<-corGSE14618_1[,colnames(corGSE14618_1) %in% row.names(mad_GSE14618_1[mad_GSE14618_1$symbol=="PAX8",])]
dim(test4)

row.names(mad_GSE14618_2[mad_GSE14618_2$symbol=="PAX8",])
test5<-corGSE14618_2[,colnames(corGSE14618_2) %in% row.names(mad_GSE14618_2[mad_GSE14618_2$symbol=="PAX8",])]
dim(test5)

row.names(mad_GSE26713[mad_GSE26713$symbol=="PAX8",])
test6<-corGSE26713[,colnames(corGSE26713) %in% row.names(mad_GSE26713[mad_GSE26713$symbol=="PAX8",])]
dim(test6)

row.names(mad_GSE8879[mad_GSE8879$symbol=="PAX8",])
test7<-corGSE8879[,colnames(corGSE8879) %in% row.names(mad_GSE8879[mad_GSE8879$symbol=="PAX8",])]
dim(test7)

row.names(mad_emtab604[mad_emtab604$symbol=="PAX8",])
test8<-coremtab604[,colnames(coremtab604) %in% row.names(mad_emtab604[mad_emtab604$symbol=="PAX8",])]
dim(test8)


testall<-merge(test2,test3,by=0)
row.names(testall)<-testall$Row.names
testall$Row.names<-NULL

testall<-merge(testall,test4,by=0)
row.names(testall)<-testall$Row.names
testall$Row.names<-NULL

testall<-merge(testall,test5,by=0)
row.names(testall)<-testall$Row.names
testall$Row.names<-NULL

testall<-merge(testall,test6,by=0)
row.names(testall)<-testall$Row.names
testall$Row.names<-NULL

testall<-merge(testall,test7,by=0)
row.names(testall)<-testall$Row.names
testall$Row.names<-NULL

testall<-merge(testall,test8,by=0)
row.names(testall)<-testall$Row.names
testall$Row.names<-NULL

testi<-as.data.frame(cor(test1,testall))
testi

# Select reference probeset = that with more times the highest correlation
df<-data.frame(times=apply(testi, 2, function(x) which.max(x)))
which.max(table(df$times))[1]
row.names(testi[which.max(table(df$times))[1],])

# Select dataframe probeset
test1

#####

##### EXPRESSION B-CAT TARGETS IN MICROARRAYS ####
#B-CAT TARGETS
sort(unique(rpmi_al[rpmi_al$peak=="peak",]$hgnc_symbol))

# Microarrays dataset
bcat_pattern<-sort(unique(rpmi_al[rpmi_al$peak=="peak",]$hgnc_symbol))
bcat_pattern<-bcat_pattern[-1]

# EXTRACTING FROM BIOMART ANNOTATIONS OF PROBES
#hg <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
listgenes=sort(unique(rpmi_al[rpmi_al$peak=="peak",]$hgnc_symbol))
bcat_att<-getBM(attributes = c("external_gene_name",'hgnc_symbol', 'chromosome_name',
                     'start_position', 'end_position','affy_hg_u133_plus_2','ensembl_gene_id','ensembl_transcript_id_version'),
      filters = 'external_gene_name', 
      values = listgenes, 
      mart = hg)


                                                                                                                                     bcat_sig1<-filter(expannot_GSE14618_1, grepl(paste0(bcat_pattern,'\\>',collapse = '|'),symbol))
bcat_sig2<-filter(expannot_GSE14618_1, grepl(paste0(setdiff(bcat_pattern,bcat_sig1$symbol),' ',collapse = '|'),symbol))
bcat_sig3<-filter(expannot_GSE14618_1, grepl(paste0(' ',setdiff(bcat_pattern,bcat_sig1$symbol),collapse = '|'),symbol))

bcat_sig4<-filter(expannot_GSE14618_2, grepl(paste0(bcat_pattern,'\\>',collapse = '|'),symbol))
#bcat_sig5<-filter(expannot_GSE14618_2, grepl(paste0(setdiff(bcat_pattern,bcat_sig4$symbol),' ',collapse = '|'),symbol))
#bcat_sig6<-filter(expannot_GSE14618_2, grepl(paste0(' ',setdiff(bcat_pattern,bcat_sig4$symbol),collapse = '|'),symbol))


bcat_sig<-rbind(bcat_sig1,bcat_sig2)
bcat_sig<-rbind(bcat_sig,bcat_sig3)

bcat_sig_2<-bcat_sig4

bcat_sig$Row.names<-NULL
bcat_sig$probe_id<-NULL

bcat_sig_2$Row.names<-NULL
bcat_sig_2$probe_id<-NULL

# Select unique probeid per gene
bcat_sig<-bcat_sig[order(bcat_sig$GSM365046,decreasing = TRUE),]
bcat_sig<-bcat_sig[!duplicated(bcat_sig$symbol),]

row.names(bcat_sig)<-bcat_sig$symbol
bcat_sig$symbol<-NULL


bcat_sig_2<-bcat_sig_2[order(bcat_sig_2$GSM365113,decreasing = TRUE),]
bcat_sig_2<-bcat_sig_2[!duplicated(bcat_sig_2$symbol),]

row.names(bcat_sig_2)<-bcat_sig_2$symbol
bcat_sig_2$symbol<-NULL


# BRCA2, bcat and KAISO levels in _1 and AXIN2
ctnnb1_exp<-grepl('CTNNB1',expannot_GSE14618_1$symbol)
kaiso_exp<-grepl('ZBTB33',expannot_GSE14618_1$symbol)
brca1_exp<-grepl('BRCA1',expannot_GSE14618_1$symbol)

df1<-expannot_GSE14618_1[ctnnb1_exp,]
df2<-expannot_GSE14618_1[kaiso_exp,]
df3<-expannot_GSE14618_1[brca1_exp,]

df_1<-as.data.frame(t(Reduce(function(x, y) merge(x, y, all=TRUE), list(df1, df2, df3))))
df_1<-df_1[-c(1:2),]
colnames(df_1)<-df_1[1,]
df_1<-df_1[-1,]
colnames(df_1)

# BRCA2, bcat and KAISO levels in _2
ctnnb1_exp2<-grepl('CTNNB1',expannot_GSE14618_2$symbol)
kaiso_exp2<-grepl('ZBTB33',expannot_GSE14618_2$symbol)
brca1_exp2<-grepl('BRCA1',expannot_GSE14618_2$symbol)
tcf_exp2<-grepl('TCF7',expannot_GSE14618_2$symbol)
lef_exp2<-grepl('LEF1',expannot_GSE14618_2$symbol)
                  
df1_2<-expannot_GSE14618_2[ctnnb1_exp2,]
df2_2<-expannot_GSE14618_2[kaiso_exp2,]
df3_2<-expannot_GSE14618_2[brca1_exp2,]
df4_2<-expannot_GSE14618_2[tcf_exp2,][c(1:2),]
df5_2<-expannot_GSE14618_2[lef_exp2,][c(1:3),]

df_2<-as.data.frame(t(Reduce(function(x, y) merge(x, y, all=TRUE), list(df1_2, df2_2, df3_2, df4_2, df5_2))))
df_2<-df_2[-c(1:2),]
colnames(df_2)<-df_2[1,]
df_2<-df_2[-1,]
colnames(df_2)<-c("bcat_exp_prob1","bcat_exp_prob2","brca1_exp_prob1",
                  "tcf_exp_prob1","tcf_exp_prob2","lef_exp_prob1",
                    "brca1_exp_prob2","kaiso_exp_prob1",
                  "lef_exp_prob2","lef_exp_prob3","bcat_exp_prob3","kaiso_exp_prob2")
  
library(RColorBrewer)
cols <- colorRampPalette(rev(brewer.pal(10,name="RdBu")))(50)
paletteLength<-length(cols)

myBreaks <- c(seq(min(scale(t(bcat_sig))), 0, length.out=round(ceiling(paletteLength*0.45))),
              seq(0.1,2.5,length.out = round(ceiling(paletteLength*0.45))),
              seq(2.6, max(scale(t(bcat_sig))), length.out=round(floor(paletteLength*0.10))))

phet<-pheatmap(as.data.frame(t(scale(t(bcat_sig)))),
               annotation_col = data.frame(row.names=row.names(pheno_GSE14618_gpl96),
                                           First_event=gsub(' \\(.*','',pheno_GSE14618_gpl96$description),
                                           bcat_exp=as.numeric(df_1$CTNNB1),
                                           brca1_exp_prob1=as.numeric(df_1[,3]),
                                           brca1_exp_prob2=as.numeric(df_1[,2]),
                                           kaiso_exp=as.numeric(df_1$ZBTB33)),
             #  annotation_row = colcol[,c(6,11,12,13,14)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 3,
               cutree_rows = 5,
               color=cols,
             breaks=myBreaks)

clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=5)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=3)))


table(patients$cluster)
table(clusters$cluster)
clusters$Gene<-row.names(clusters)
clusters[clusters$cluster==3,]$Gene

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(t(bcat_sig))))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(t(bcat_sig))))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


colnames(df_1)<-c("bcat_exp","brca1_exp_prob2","brca1_exp_prob1","kaiso_exp")
patinfo<-merge(df_1,patients,by="row.names")
patinfo$bcat_exp<-as.numeric(patinfo$bcat_exp)
patinfo$brca1_exp_prob1<-as.numeric(patinfo$brca1_exp_prob1)
patinfo$brca1_exp_prob2<-as.numeric(patinfo$brca1_exp_prob2)

ggplot(patinfo,aes(x=reorder(as.factor(cluster), bcat_exp, FUN = median),y=scale(bcat_exp)))+
  geom_point(aes(size=order))+
  geom_boxplot(alpha=0.2,lwd=1)+
  geom_hline(yintercept = median(scale(patinfo$bcat_exp)),color="red",lwd=2)+
  geom_hline(yintercept = quantile(scale(patinfo$bcat_exp))[4],color="grey",lwd=2)+
  geom_hline(yintercept = quantile(scale(patinfo$bcat_exp))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))


### pheno data survival for _2
survival_GSE14618<-data.frame(read_excel("../GSE14618/survival_GSE14618.xlsx"))
survival_GSE14618<-survival_GSE14618[survival_GSE14618$GEO_Array_ID!="ND",]
row.names(survival_GSE14618)<-survival_GSE14618$GEO_Array_ID

#pheatmap only with survival data
t_bcat_sig_2<-t(bcat_sig_2)
t_bcat_sig_2<-t_bcat_sig_2[row.names(t_bcat_sig_2) %in% row.names(survival_GSE14618), ]

pheno_2<-pheno_GSE14618_gpl570[row.names(pheno_GSE14618_gpl570) %in% row.names(t_bcat_sig_2), ]
df_2_surv<-df_2[row.names(df_2) %in% row.names(t_bcat_sig_2), ]

pheno_2<-merge(pheno_2,df_2_surv,by="row.names")
row.names(pheno_2)<-pheno_2$Row.names
pheno_2$Row.names<-NULL
pheno_2<-merge(pheno_2,survival_GSE14618,by="row.names")
row.names(pheno_2)<-pheno_2$Row.names
pheno_2$Row.names<-NULL

pheno_2$bcat_exp_prob1<-as.numeric(pheno_2$bcat_exp_prob1)
pheno_2$bcat_exp_prob2<-as.numeric(pheno_2$bcat_exp_prob2)
pheno_2$bcat_exp_prob3<-as.numeric(pheno_2$bcat_exp_prob3)
pheno_2$brca1_exp_prob1<-as.numeric(pheno_2$brca1_exp_prob1)
pheno_2$brca1_exp_prob2<-as.numeric(pheno_2$brca1_exp_prob2)
pheno_2$kaiso_exp_prob1<-as.numeric(pheno_2$kaiso_exp_prob1)
pheno_2$kaiso_exp_prob2<-as.numeric(pheno_2$kaiso_exp_prob2)
pheno_2$lef_exp_prob1<-as.numeric(pheno_2$lef_exp_prob1)
pheno_2$lef_exp_prob2<-as.numeric(pheno_2$lef_exp_prob2)
pheno_2$lef_exp_prob3<-as.numeric(pheno_2$lef_exp_prob3)
pheno_2$tcf_exp_prob1<-as.numeric(pheno_2$tcf_exp_prob1)
pheno_2$tcf_exp_prob2<-as.numeric(pheno_2$tcf_exp_prob2)

pheno_2$Pheno<-gsub(', non-ABD','',pheno_2$T.ALL_Subset)
pheno_2$Pheno<-gsub(';.*','',pheno_2$Pheno)

myBreaks <- c(seq(min(scale(t_bcat_sig_2)), 0, length.out=round(ceiling(paletteLength*0.45))),
              seq(0.1,2.5,length.out = round(ceiling(paletteLength*0.45))),
              seq(2.6, max(scale(t_bcat_sig_2)), length.out=round(floor(paletteLength*0.10))))

phet<-pheatmap(as.data.frame(t(scale(t_bcat_sig_2))),
         annotation_col = pheno_2[,c(53,56,61,63)],
         #  annotation_row = colcol[,c(6,11,12,13,14)],
         show_colnames=F,
         show_rownames = F,
         fontsize = 6,
         cutree_cols = 4,
         cutree_rows = 3,
         color=cols,
         breaks=myBreaks)


# Heatmap with b-cat targets within DEGs TARGET remission relapse
deg_bcattarg_sym<-target_prog_degs_sig[row.names(target_prog_degs_sig) %in% colnames(bcatexp_survival[,21:(ncol(bcatexp_survival)-4)]),]$NAME

myBreaks <- c(seq(min(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% deg_bcattarg_sym])), 0, length.out=round(ceiling(paletteLength*0.45))),
              seq(0.1,2.5,length.out = round(ceiling(paletteLength*0.45))),
              seq(2.6, max(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% deg_bcattarg_sym])), length.out=round(floor(paletteLength*0.10))))

phet<-pheatmap(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% deg_bcattarg_sym]))),
               annotation_col = pheno_2[,c(61,56,37:43)],
               #  annotation_row = colcol[,c(6,11,12,13,14)],
               show_colnames=F,
               show_rownames = T,
               fontsize = 6,
               cutree_cols = 3,
               cutree_rows = 3,
               color=cols,
               breaks=myBreaks)



clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=3)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=4)))


table(patients$cluster)
table(clusters$cluster)
clusters$Gene<-row.names(clusters)
clusters[clusters$cluster==1,]$Gene

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(t_bcat_sig_2)))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% unique(deg_clust_GSE14618[deg_clust_GSE14618$hgnc_symbol !="NA",]$hgnc_symbol)])))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% deg_bcattarg_sym])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(t_bcat_sig_2)))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% unique(deg_clust_GSE14618[deg_clust_GSE14618$hgnc_symbol !="NA",]$hgnc_symbol)])))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% deg_bcattarg_sym])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(pheno_2,patients,by="row.names")
patinfo$bcat_exp_prob1<-as.numeric(patinfo$bcat_exp_prob1)
patinfo$bcat_exp_prob2<-as.numeric(patinfo$bcat_exp_prob2)
patinfo$bcat_exp_prob3<-as.numeric(patinfo$bcat_exp_prob3)

patinfo$brca1_exp_prob1<-as.numeric(patinfo$brca1_exp_prob1)
patinfo$brca1_exp_prob2<-as.numeric(patinfo$brca1_exp_prob2)

patinfo$kaiso_exp_prob1<-as.numeric(patinfo$kaiso_exp_prob1)
patinfo$kaiso_exp_prob2<-as.numeric(patinfo$kaiso_exp_prob2)

patinfo$tcf_exp_prob1<-as.numeric(patinfo$tcf_exp_prob1)
patinfo$tcf_exp_prob2<-as.numeric(patinfo$tcf_exp_prob2)

patinfo$lef_exp_prob1<-as.numeric(patinfo$lef_exp_prob1)
patinfo$lef_exp_prob2<-as.numeric(patinfo$lef_exp_prob2)

patinfo<-patinfo[,-2]

ggplot(patinfo,aes(x=reorder(as.factor(cluster), bcat_exp_prob1, FUN = median),y=scale(kaiso_exp_prob2)))+
#ggplot(patinfo,aes(x=Outcome,y=scale(bcat_exp_prob1)))+
  geom_point(aes(color=Outcome),size=5)+
  geom_boxplot(alpha=0.2,lwd=1)+
#  geom_hline(yintercept = median(scale(patinfo$bcat_exp_prob1)),color="red",lwd=2)+
#  geom_hline(yintercept = median(scale(patinfo$bcat_exp_prob1))-0.2,color="grey",lwd=2)+
#  geom_hline(yintercept = median(scale(patinfo$kaiso_exp_prob2))+0.25,color="grey",lwd=2)+
#  geom_hline(yintercept = 0,color="grey",lwd=1)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))

patinfo$scale_bcat_prob1<-scale(patinfo$bcat_exp_prob1)
mean(patinfo$scale_bcat_prob1)
median(patinfo[patinfo$cluster==6,]$scale_bcat_prob1)

clust1<-getBM(attributes = c("external_gene_name",'hgnc_symbol', 'chromosome_name',
                             'start_position', 'end_position','ensembl_gene_id'),
              filters = 'external_gene_name', 
              values = clusters$Gene, 
              mart = hg)
colnames(clust1)[1]<-"Gene"

clust_microarray<-merge(clusters,clust1,by="Gene")
clust1_microarray<-unique(clust_microarray[clust_microarray$cluster==1,]$ensembl_gene_id)



### Survival curves in bcat_sig_2 (GSE14618_2)

# Create variable brca1 expression
pheno_2$brcagroup_prob1<-NA
pheno_2$brcagroup_prob1<-ifelse(pheno_2$brca1_exp_prob1>median(pheno_2$brca1_exp_prob1),
                                   "high_brca1",
                                   "low_brca1")

pheno_2$brcagroup_prob2<-ifelse(pheno_2$brca1_exp_prob2>median(pheno_2$brca1_exp_prob2),
                                          "high_brca1",
                                          "low_brca1")



# Create variable bcat expression
pheno_2$bcat_prob1<-NA
pheno_2$bcat_prob1<-ifelse(pheno_2$bcat_exp_prob1>median(pheno_2$bcat_exp_prob1),
                                "high_bcat",
                                "low_bcat")

pheno_2$bcat_prob2<-ifelse(pheno_2$bcat_exp_prob2>median(pheno_2$bcat_exp_prob2),
                           "high_bcat",
                           "low_bcat")

pheno_2$tcf_prob1<-ifelse(pheno_2$tcf_exp_prob1>median(pheno_2$tcf_exp_prob1),
                           "high_tcf",
                           "low_tcf")

pheno_2$lef_prob1<-ifelse(pheno_2$lef_exp_prob1>median(pheno_2$lef_exp_prob1),
                          "high_lef",
                          "low_lef")

pheno_2$bcat_extr_l<-ifelse(pheno_2$bcat_exp_prob1<quantile(pheno_2$bcat_exp_prob1)[2],
                           "extr_low_bcat",
                           "no_extr_low_bcat")

pheno_2$bcat_extr_h<-ifelse(pheno_2$bcat_exp_prob1>quantile(pheno_2$bcat_exp_prob1)[4],
                          "extr_high_bcat",
                          "no_extr_high_bcat")

pheno_2$bcat_extr<-paste(pheno_2$bcat_extr_l,pheno_2$bcat_extr_h,sep="_")

# Create variable kaiso expression
pheno_2$kaiso_prob1<-NA
pheno_2$kaiso_prob1<-ifelse(pheno_2$kaiso_exp_prob1>median(pheno_2$kaiso_exp_prob1),
                           "high_kaiso",
                           "low_kaiso")

pheno_2$kaiso_prob2<-ifelse(pheno_2$kaiso_exp_prob2>median(pheno_2$kaiso_exp_prob2),
                           "high_kaiso",
                           "low_kaiso")


pheno_2_pat<-merge(patients,pheno_2,by="row.names")
row.names(pheno_2_pat)<-pheno_2_pat$Row.names
pheno_2_pat$Row.names<-NULL

pheno_2_pat$newcluster<-NA
pheno_2_pat<-pheno_2_pat %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 1,
                                cluster == 2 ~ 23,
                                cluster == 3 ~ 23,
                                cluster == 4 ~ 4
  ))


pheno_2_pat$Event<-as.integer(pheno_2_pat$Event)
pheno_2_pat$Alive<-as.integer(pheno_2_pat$Alive)

surv_C_fit_TALL <- survfit(Surv(EFS_yrs, Event) ~ bcat_extr, data=pheno_2_pat)
ggsurvplot(surv_C_fit_TALL, data = pheno_2_pat,
           pval = TRUE,
           size=2,
           palette = c("red", "gray60","gray87"),
           ggtheme = theme_bw()+theme(legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      text = element_text(size=20),
                                      axis.text.x = element_text(angle=45, size=20,hjust=1),
                                      axis.text.y = element_text(size=20,hjust=1)),
           font.x = c(20),
           font.y = c(20))



## DEGs between clusters or relapse vs remission
dt<-data.frame(row.names=GSE14618[[1]]@phenoData[[2]],
           patient=GSE14618[[1]]@phenoData[[2]])

patients

dt<-merge(dt,patients,by=0,all.x=TRUE)


dt$newcluster<-NA
dt<-dt %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 145,
                                cluster == 2 ~ 2367,
                                cluster == 3 ~ 2367,
                                cluster == 4 ~ 145,
                                cluster == 5 ~ 145,
                                cluster == 6 ~ 2367,
                                cluster == 7 ~ 2367
  ))

dt<-dt[complete.cases(dt), ]
row.names(dt)<-dt$Row.names
dt$Row.names<-NULL


dt$newcluster<-as.character(dt$newcluster)



# by cluster
group<-dt$newcluster
design<-model.matrix(~group)

degs_GSE14618<-exprs(GSE14618[[1]])[,colnames(exprs(GSE14618[[1]])) %in% row.names(dt)]

fit <-lmFit(degs_GSE14618,design)
fit<-eBayes(fit)

top.all <- topTable(fit,n=nrow(degs_GSE14618),adjust="BH",coef=2) #This will give you all the pvalue when comparing the samples (for 2 condition)

# by relapse/remission or indfailure/remission or indfailure/relapse

dt<-data.frame(row.names=GSE14618[[1]]@phenoData[[2]],
               patient=GSE14618[[1]]@phenoData[[2]])
dt<-merge(dt,pheno_2,by=0)

row.names(dt)<-dt$patient

#induction failure and relapse
dt_rem_rel<-dt[dt$Outcome=="Long-term event-free survivor" |  dt$Outcome=="Relapse",]
dt_rem_ind_fail<-dt[dt$Outcome=="Induction failure" |  dt$Outcome=="Long-term event-free survivor",]
dt_rel_ind_fail<-dt[dt$Outcome=="Induction failure" |  dt$Outcome=="Relapse",]

group<-dt_rem_rel$Outcome

group<-dt_rem_ind_fail$Outcome

group<-dt_rel_ind_fail$Outcome

#change variable design
design<-model.matrix(~group)

degs_GSE14618<-exprs(GSE14618[[1]])[,colnames(exprs(GSE14618[[1]])) %in% row.names(dt_rem_rel)]
degs_GSE14618<-exprs(GSE14618[[1]])[,colnames(exprs(GSE14618[[1]])) %in% row.names(dt_rem_ind_fail)]
degs_GSE14618<-exprs(GSE14618[[1]])[,colnames(exprs(GSE14618[[1]])) %in% row.names(dt_rel_ind_fail)]

fit <-lmFit(degs_GSE14618,design)
fit<-eBayes(fit)

top.all_rem_rel <- topTable(fit,n=nrow(degs_GSE14618),adjust="BH",coef=2) #This will give you all the pvalue when comparing the samples (for 2 condition)
top.all_rem_rel$comp<-"Rem_Rel"

top.all_rem_ind_fail <- topTable(fit,n=nrow(degs_GSE14618),adjust="BH",coef=2) #This will give you all the pvalue when comparing the samples (for 2 condition)
top.all_rem_ind_fail$comp<-"Rem_indfail"

top.all_rel_ind_fail <- topTable(fit,n=nrow(degs_GSE14618),adjust="BH",coef=2) #This will give you all the pvalue when comparing the samples (for 2 condition)
top.all_rel_ind_fail$comp<-"Rel_indfail"

# Annotate
Annot <- toTable(hgu133plus2SYMBOL)
colnames(Annot)
row.names(Annot)<-Annot$probe_id

deg_GSE14618surv<-merge(top.all_rem_rel,Annot,by=0,all.x=TRUE)
deg_GSE14618surv<-merge(top.all_rem_ind_fail,Annot,by=0,all.x=TRUE)
deg_GSE14618surv<-merge(top.all_rel_ind_fail,Annot,by=0,all.x=TRUE)

degall<-getBM(attributes = c("external_gene_name",'hgnc_symbol', 'chromosome_name',
                               'start_position', 'end_position','ensembl_gene_id'),
                filters = 'external_gene_name', 
                values = unique(deg_GSE14618surv$symbol), 
                mart = hg)

deg_GSE14618surv<-merge(degall,deg_GSE14618surv,by.x="external_gene_name",by.y="symbol",all.y=TRUE)
deg_GSE14618surv$reg<-ifelse((deg_GSE14618surv$logFC<0),"up","down")


deg_clust_GSE14618<-subset(deg_GSE14618surv,deg_GSE14618surv$P.Value<0.05)
genedeg<-deg_clust_GSE14618[,c(3,6,14,15)]
genedeg<-genedeg[!duplicated(genedeg$ensembl_gene_id),]
genedeg<-genedeg[!is.na(genedeg$ensembl_gene_id),]
row.names(genedeg)<-genedeg$ensembl_gene_id
colnames(genedeg)[1]<-"chr"

# b-cat targets DEGs
targdeg<-colnames(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% unique(deg_clust_GSE14618[deg_clust_GSE14618$hgnc_symbol !="NA",]$hgnc_symbol)])
unique(deg_clust_GSE14618[deg_clust_GSE14618$external_gene_name %in% targdeg,]$ensembl_gene_id)
deg_bcatmicro_sym<-unique(deg_clust_GSE14618[deg_clust_GSE14618$external_gene_name %in% targdeg,]$hgnc_symbol)

#annot genes per comparison
rem_rel<-deg_clust_GSE14618[deg_clust_GSE14618$external_gene_name %in% targdeg,]
rem_rel<-rem_rel[!duplicated(rem_rel$external_gene_name),]

rem_indfail<-deg_clust_GSE14618[deg_clust_GSE14618$external_gene_name %in% targdeg,]
rem_indfail<-rem_indfail[!duplicated(rem_indfail$external_gene_name),]

rel_indfail<-deg_clust_GSE14618[deg_clust_GSE14618$external_gene_name %in% targdeg,]
rel_indfail<-rel_indfail[!duplicated(rel_indfail$external_gene_name),]

comp_all<-rbind(rem_rel,rem_indfail)
#comp_all<-rbind(comp_all,rel_indfail)
comp_all<-comp_all[!duplicated(comp_all$hgnc_symbol),]
row.names(comp_all)<-comp_all$hgnc_symbol

# Heatmap with b-cat targets DE between relapse and remission in microarrays
# with no induction failure
#noind<-row.names(pheno_2[pheno_2$Outcome!="Induction failure",])

myBreaks <- c(seq(min(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% comp_all$hgnc_symbol])), 0, length.out=round(ceiling(paletteLength*0.45))),
              seq(0.1,2.5,length.out = round(ceiling(paletteLength*0.45))),
              seq(2.6, max(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% comp_all$hgnc_symbol])), length.out=round(floor(paletteLength*0.10))))

pheno_2$bcat_scale<-scale(pheno_2$bcat_exp_prob1)
pheno_2$kaiso_scale<-scale(pheno_2$kaiso_exp_prob1)
pheno_2$tcf_scale<-scale(pheno_2$tcf_exp_prob1)
pheno_2$lef_scale<-scale(pheno_2$lef_exp_prob1)

phet<-pheatmap(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% comp_all$hgnc_symbol]))),
               annotation_col = pheno_2[,c(54,58,62,64)],
               #  annotation_row = colcol[,c(6,11,12,13,14)],
               annotation_row = comp_all[,c(14,16)],
               show_colnames=F,
               show_rownames = T,
               fontsize = 4,
               cutree_cols = 3,
               cutree_rows = 3,
               color=cols,
               breaks=myBreaks)



clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=3)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=3)))


table(patients$cluster)
table(clusters$cluster)
clusters$Gene<-row.names(clusters)
clusters[clusters$cluster==1,]$Gene

# Order genes (clusters)
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(t_bcat_sig_2)))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% unique(deg_clust_GSE14618[deg_clust_GSE14618$hgnc_symbol !="NA",]$hgnc_symbol)])))[phet$tree_row[["order"]],]))
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% comp_all$hgnc_symbol])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(t_bcat_sig_2)))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% unique(deg_clust_GSE14618[deg_clust_GSE14618$hgnc_symbol !="NA",]$hgnc_symbol)])))[,phet$tree_col[["order"]]]))
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(t_bcat_sig_2[,colnames(t_bcat_sig_2) %in% comp_all$hgnc_symbol])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(pheno_2,patients,by="row.names")
patinfo<-patinfo[,-2]
patinfo$bcat_exp_prob1<-as.numeric(patinfo$bcat_exp_prob1)
patinfo$bcat_exp_prob2<-as.numeric(patinfo$bcat_exp_prob2)
patinfo$bcat_exp_prob3<-as.numeric(patinfo$bcat_exp_prob3)

patinfo$brca1_exp_prob1<-as.numeric(patinfo$brca1_exp_prob1)
patinfo$brca1_exp_prob2<-as.numeric(patinfo$brca1_exp_prob2)

patinfo$kaiso_exp_prob1<-as.numeric(patinfo$kaiso_exp_prob1)
patinfo$kaiso_exp_prob2<-as.numeric(patinfo$kaiso_exp_prob2)

patinfo$tcf_exp_prob1<-as.numeric(patinfo$tcf_exp_prob1)
patinfo$tcf_exp_prob2<-as.numeric(patinfo$tcf_exp_prob2)

patinfo$lef_exp_prob1<-as.numeric(patinfo$lef_exp_prob1)
patinfo$lef_exp_prob2<-as.numeric(patinfo$lef_exp_prob2)
patinfo$lef_exp_prob3<-as.numeric(patinfo$lef_exp_prob3)

ggplot(patinfo,aes(x=reorder(as.factor(cluster), bcat_exp_prob1, FUN = median),y=scale(kaiso_exp_prob2)))+
  geom_point(aes(color=Outcome),size=5)+
  geom_boxplot(alpha=0.2,lwd=1)+
  scale_color_manual(values=c("chartreuse3","hotpink","deepskyblue"))+
  #  geom_hline(yintercept = median(scale(patinfo$bcat_exp_prob1)),color="red",lwd=2)+
  #  geom_hline(yintercept = median(scale(patinfo$bcat_exp_prob1))-0.2,color="grey",lwd=2)+
  #  geom_hline(yintercept = median(scale(patinfo$kaiso_exp_prob2))+0.25,color="grey",lwd=2)+
  # geom_hline(yintercept = 0,color="grey",lwd=1)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_text(angle=45, size=16,hjust=1),
        axis.title.x = element_blank())

patinfo_melt<-subset(patinfo,select=c("Row.names","Outcome","cluster","bcat_exp_prob1","kaiso_exp_prob2","tcf_exp_prob1","lef_exp_prob1"))
patinfo_melt<-melt(patinfo_melt,id.vars=c("Row.names","Outcome","cluster"))

patinfo_melt$cluster <-as.factor(patinfo_melt$cluster)
patinfo_melt$cluster <- ordered(patinfo_melt$cluster, levels = c("2", "3", "1"))
ggplot(patinfo_melt,aes(x=as.factor(cluster),y=value))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=Outcome),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_color_manual(values=c("chartreuse3","hotpink","deepskyblue"))+
  ylab("scaled expression")+
  xlab('')+
  facet_wrap(~variable,scales = "free",ncol=2)+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=12),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_text(angle=45, size=16,hjust=1),
        axis.title.x = element_blank())

my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3") )
ggboxplot(patinfo, x = "cluster", y = "kaiso_exp_prob2",color="cluster", palette = "jco")+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(method="anova")


prop_clust<-data.frame(prop.table(table(patinfo$cluster,patinfo$Outcome),1))
colnames(prop_clust)<-c("cluster","Outcome","Proportion")
prop_clust$cluster<-as.factor(prop_clust$cluster)

prop_clust$cluster<- factor(prop_clust$cluster, levels=c("2","3","1"))

ggplot(prop_clust,aes(x=cluster,y=Proportion,fill=Outcome))+
  geom_col(color="black")+
  scale_fill_manual(values=c("chartreuse3","hotpink","deepskyblue"))+
  theme_bw()+
  theme(legend.position="bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=16))


patinfo$scale_bcat_prob1<-scale(patinfo$bcat_exp_prob1)
mean(patinfo$scale_bcat_prob1)
median(patinfo[patinfo$cluster==6,]$scale_bcat_prob1)



### Survival curves in bcat_sig_2 (GSE14618_2)

pheno_2_pat<-merge(patients,pheno_2,by="row.names")
row.names(pheno_2_pat)<-pheno_2_pat$Row.names
pheno_2_pat$Row.names<-NULL

pheno_2_pat$newcluster<-NA
pheno_2_pat<-pheno_2_pat %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 12,
                                cluster == 2 ~ 12,
                                cluster == 3 ~ 3,
                                cluster == 4 ~ 4,
                                cluster == 5 ~ 5
  ))


pheno_2_pat$Event<-as.integer(pheno_2_pat$Event)
pheno_2_pat$Alive<-as.integer(pheno_2_pat$Alive)

pheno_2_pat$Alive_bin<-gsub('1','zero',pheno_2_pat$Alive)
pheno_2_pat$Alive_bin<-gsub('0','one',pheno_2_pat$Alive_bin)

pheno_2_pat$Alive_bin<-gsub('zero',0,pheno_2_pat$Alive_bin)
pheno_2_pat$Alive_bin<-gsub('one',1,pheno_2_pat$Alive_bin)

pheno_2_pat$Alive_bin<-as.integer(pheno_2_pat$Alive_bin)

surv_C_fit_TALL <- survfit(Surv(Surv_yrs, Alive_bin) ~ cluster, data=pheno_2_pat)
ggsurvplot(surv_C_fit_TALL, data = pheno_2_pat,
           pval = TRUE,
           size=2,
           palette = c("red", "gray87","gray60"),
           ggtheme = theme_bw()+theme(legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      text = element_text(size=20),
                                      axis.text.x = element_text(angle=45, size=20,hjust=1),
                                      axis.text.y = element_text(size=20,hjust=1)),
           font.x = c(20),
           font.y = c(20))

pheno_2_pat$Age_hr<-ifelse(pheno_2_pat$Age<10,"<10",">10")
table(pheno_2_pat$CNS_Status)

pheno_2_pat$cluster<-as.factor(pheno_2_pat$cluster)
pheno_2_pat$cluster <- factor(pheno_2_pat$cluster, levels = c("3","2","1"))
#pheno_2_pat$CNS_Status <- factor(pheno_2_pat$CNS_Status, levels = c("CNS 1","CNS 2","CNS 3","Bloody tap, blasts","No data"))
pheno_2_pat$Pheno_ETP<-ifelse(pheno_2_pat$Pheno=="ETP","ETP","nonETP")
pheno_2_pat$Pheno_ETP <- factor(pheno_2_pat$Pheno_ETP, levels = c("nonETP","ETP"))

fit.coxph<- coxph(Surv(EFS_yrs, Event) ~ cluster+Gender+Age_hr+CNS_Status+Pheno_ETP,
                  data = pheno_2_pat[pheno_2_pat$CNS_Status!="No data" & !is.na(pheno_2_pat$Pheno),],
                  ties = 'exact')
summary(fit.coxph)

ggforest(fit.coxph)



# ENRICHR
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "InterPro_Domains_2019",
         "GO_Cellular_Component_2017",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2016","DisGeNET"
         )


clust1bcat<-enrichr(as.character(clusters[clusters$cluster==1,]$Gene),dbs)
clust2bcat<-enrichr(as.character(clusters[clusters$cluster==2,]$Gene),dbs)
clust3bcat<-enrichr(as.character(clusters[clusters$cluster==3,]$Gene),dbs)
clust4bcat<-enrichr(as.character(clusters[clusters$cluster==4,]$Gene),dbs)
clust5bcat<-enrichr(as.character(clusters[clusters$cluster==5,]$Gene),dbs)

degprog_clust1<-enrichr(targdeg[targdeg %in% clusters[clusters$cluster==1,]$Gene],dbs)
degprog_clust2<-enrichr(targdeg[targdeg %in% clusters[clusters$cluster==2,]$Gene],dbs)
degprog_clust3<-enrichr(targdeg[targdeg %in% clusters[clusters$cluster==3,]$Gene],dbs)

degprog<-enrichr(comp_all$hgnc_symbol,dbs)

#degprog_down<-enrichr(unique(deg_GSE14618surv[deg_GSE14618surv$adj.P.Val<0.01 & deg_GSE14618surv$logFC>0,]$hgnc_symbol),dbs)
#degprog_up<-enrichr(unique(deg_GSE14618surv[deg_GSE14618surv$adj.P.Val<0.01 & deg_GSE14618surv$logFC<0,]$hgnc_symbol),dbs)

clust1bcatrich <-clust1bcat[["GO_Biological_Process_2018"]]
clust2bcatrich <-clust2bcat[["GO_Biological_Process_2018"]]
clust3bcatrich <-clust3bcat[["GO_Biological_Process_2018"]]
clust4bcatrich <-clust4bcat[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
clust5bcatrich <-clust5bcat[["GO_Biological_Process_2018"]]

degprogrich_up <-degprog_up[["DisGeNET"]]
degprogrich_down <-degprog_down[["DisGeNET"]]

degprogrich <-degprog[["GO_Biological_Process_2018"]]

degprog_clust1 <-degprog_clust1[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
degprog_clust2 <-degprog_clust2[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
degprog_clust3 <-degprog_clust3[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]

clust1bcatrich$group<-c("cluster1")
clust2bcatrich$group<-c("cluster2")
clust3bcatrich$group<-c("cluster3")
clust4bcatrich$group<-c("cluster4")
clust5bcatrich$group<-c("cluster5")


degprog_clust1$group<-"clust1"
degprog_clust2$group<-"clust2"
degprog_clust3$group<-"clust3"

degprogrich$group<-"degprog"

degprogrich_down$group<-"down"
degprogrich_up$group<-"up"


allGO<-rbind(clust1bcatrich,clust2bcatrich,clust3bcatrich)
#,clust4bcatrich)
#             ,clust5bcatrich)

allGO<-rbind(degprogrich_up,degprogrich_down)

allGO<-degprogrich

bpsub<-subset(allGO,allGO$Adjusted.P.value<0.01)


bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=unique((bpsub$Term))[rev(order(bpsub$P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))


# for multiple groups
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$group, rev(abs(bpsub$Combined.Score)))), ]$Term))

bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$group, abs(bpsub$P.value))), ]$Term))

### for facets multiple groups
ggplot(bpsub,aes(y=-(log(P.value)),x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group))+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~group,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.position = "bottom",
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=8, hjust = -1),
        axis.text.y = element_text(size=14, hjust = 1))

nondup_cat<-as.data.frame(table(bpsub$Term))[as.data.frame(table(bpsub$Term))[,2]<2,]$Var1

# heatmap
ggplot(bpsub[bpsub$Term %in% nondup_cat,], aes(y=group, x=Term)) + 
  geom_tile(aes(fill = -(log(P.value))),colour = "black",size=0.5)+
  scale_fill_gradient2(low = "white",
                       high = "red",
                       midpoint = 7) +
  #scale_fill_gradient(low = "blue", high = "red")+
  theme_cleantable()+
  theme(axis.text.x = element_text(size=16, hjust = 1,vjust=0.5,angle=90),
        axis.text.y = element_text(size=12, hjust = 1,face="italic"))+
  ggtitle("Biological Processes")



