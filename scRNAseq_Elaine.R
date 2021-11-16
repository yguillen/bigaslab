## Load libraries

library(devtools)
library(Rtsne)
library(DESeq2)
library(scran)
require(knitr)
library(scRNAseq)
#library(bglab)
library(ggrepel)

#Set seed to reproduce the data every time
set.seed(100)

# Set working directory
#setwd("/Volumes/grcmc/YGUILLEN/scRNA_cam/Elaine_CD27/scRNASeq_processed/")
setwd("/Users/yguillen/Desktop/temp/scRNA_cam/Elaine_CD27/scRNASeq_processed/")

#counts<-read.table("counts.txt",header=TRUE,row.names=1,sep=" ",stringsAsFactors = FALSE)
counts_el<-read.table("htseq_counts.txt",header=TRUE,row.names=1,stringsAsFactors = FALSE)
counts_el<-as.matrix(counts_el)
colnames(counts_el)
rownames(counts_el)
nrow(counts_el)
dim(counts_el)



#meta before his quality control
# Do this previous in batch: sed 's/,/       /g' All_plates.tsv | sed 's/\.//g' > All_plates_comma.tsv
# CAREFULL First six columns FSC and SSC, have been multiplied x1000
metael<-read.delim("/Users/yguillen/Desktop/temp/scRNA_cam/Elaine_CD27/scRNASeq_processed/meta_all.txt",sep="\t",stringsAsFactors = FALSE)
dim(metael)
row.names(metael)<-metael$ID

colnames(counts_el)
row.names(metael)

dim(counts_el)
dim(metael)


metael<- metael[order(metael$ID),]


genetabel<-read.delim("scanpy_data/genes.tsv",header=FALSE)
colnames(genetabel)<-c("Ensembl_Gene_ID","Associated_Gene_Name","biotype")
nrow(genetabel)



# Check dimensions
dim(genetabel)
dim(counts_el)

genetabel<-genetabel[!duplicated(genetabel$Ensembl_Gene_ID),]


# Loading data into the single cell dataset (SCD) object which is required for the bglab package
colnames(metael)
head(rownames(metael))
head(colnames(counts_el))

#Probably cells lost in quality control (21 in total). It should be 0!
setdiff(colnames(counts_el),rownames(metael))


# Cell cycle assignment using cyclone.
# Use of pre-defined classifier to assign cells into their cell cycle
# phases. This classifier was constructed from a training data set
# by identifying pairs of genes where the difference in expression
# within each pair changed sign accross phases.

# classifiers for human and mouse are provided. For other systems they can be constructed using sandbag package
hs.pairs <- readRDS(system.file("exdata","human_cycle_markers.rds",package = "scran"))
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

scel <- SingleCellExperiment(assays = list(counts = counts_el))
counts_el <- assay(scel, "counts")
libsizes <- colSums(counts_el)
size.factors <- libsizes/mean(libsizes)
logcounts(scel) <- log2(t(t(counts_el)/size.factors) + 1)
assayNames(scel)

expall_el <- assay(scel, "logcounts")



#for doing cyclone we need to change the gene names ENS in this case
#Do cyclone before filterGene(scd)
# This is because the lack of expression of particular genes can provide 
# some information about the cell cycle.
cycel<-cyclone(expall_el,mm.pairs)
cyctabel<-as.data.frame(cycel$phases)
row.names(cyctabel)<-colnames(expall_el)
colnames(cyctabel)<-"Phase"

cycelG1<-as.data.frame(cycel$scores$G1)
row.names(cycelG1)<-colnames(expall_el)
colnames(cycG1)<-"G1"

cycelG2M<-as.data.frame(cycel$scores$G2M)
row.names(cycelG2M)<-colnames(expall_el)
colnames(cycelG2M)<-"G2M"

cycelS<-as.data.frame(cycel$scores$S)
row.names(cycelS)<-colnames(expall_el)
colnames(cycelS)<-"S"

#write.table(cyctab,"cyclone_scRNAseq_R.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
#write.table(cycG1,"../G1_score_cyclone.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
#write.table(cycG2M,"../G2M_score_cyclone.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
#write.table(cycS,"../S_score_cyclone.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)


######### Import scanpy results


### SCANPY RESULTS 
#setwd("/Volumes/grcmc/YGUILLEN/scRNA_cam/Elaine_CD27/scRNASeq_processed/scanpy_data/")
setwd("/Users/yguillen/Desktop/temp/scRNA_cam/Elaine_CD27/scRNASeq_processed/scanpy_data/")

UMAPEl<-read.delim("umap_coords.tsv")
colnames(UMAPEl)<-c("sampleID","UMAP_C1","UMAP_C2")
row.names(UMAPEl)<-UMAPEl$sampleID

exprmatEl<-read.delim("exprMatrix.tsv")
row.names(exprmatEl)<-exprmatEl$gene

metascEl<-read.delim("meta.tsv")
colnames(metascEl)<-c("sampleID","Louvain_cluster","Mito_perc","expr_genes","UMICount")
row.names(metascEl)<-metascEl$sampleID
metascEl$sampleID<-NULL


## Results UMAP Louvain scanpy

scanpyEl<-merge(UMAPEl,metascEl,by="row.names")
colnames(meta)[1]<-"sampleID"
scanpyEl<-merge(meta,scanpyEl,by="sampleID")
row.names(scanpyEl)<-scanpyEl$sampleID
scanpyEl$Row.names=NULL


ggplot(scanpyEl,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=4)+
  geom_point(aes(color=as.factor(Louvain_cluster)),size=3)+
  #scale_color_brewer(palette="Spectral")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")


## Results preprocessed Vink and Xiaonan
setwd("/Users/yguillen/Desktop/temp/scRNA_cam/Elaine_CD27/scRNASeq_processed/")

metapre<-read.delim("meta_all.txt")
colnames(metapre)
row.names(metapre)<-metapre$ID


ggplot(metapre,aes(x=DR1,y=DR2,labels=ID))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=as.factor(Cluster)),size=2)+
  #scale_color_brewer(palette="Spectral")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

# merge metadata from preprocessed Elaine and our scanpy process
Eldata<-merge(scanpyEl,metapre,by.x="sampleID",by.y="ID")

# Add expression data from our scanpy expression results from htseq
dim(exprmatEl)
dim(Eldata)

t_exprmatEl<-data.frame(t(exprmatEl))
t_exprmatEl<-t_exprmatEl[-1,]
class(t_exprmatEl)

ix <- 1:ncol(t_exprmatEl)
t_exprmatEl[ix] <- lapply(t_exprmatEl[ix], as.character)
t_exprmatEl[ix] <- lapply(t_exprmatEl[ix], as.numeric)

row.names(t_exprmatEl)
colnames(t_exprmatEl)<-gsub('ENS.*\\.','',colnames(t_exprmatEl))

Eldata$sampleID<-as.character(Eldata$sampleID)

#merge
Eldataexp<-merge(Eldata,t_exprmatEl,by.x="sampleID",by.y="row.names")

Gfi_e<-ggplot(Eldataexp,aes(x=DR1,y=DR2,labels=sampleID))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=Gfi1),size=2)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

Jag1_e<-ggplot(Eldataexp,aes(x=DR1,y=DR2,labels=sampleID))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=Jag1),size=2)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")


Nfkbia_e<-ggplot(Eldataexp,aes(x=DR1,y=DR2,labels=sampleID))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=Nfkbia),size=2)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

Gata2_e<-ggplot(Eldataexp,aes(x=DR1,y=DR2,labels=sampleID))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=Gata2),size=2)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

Kit_e<-ggplot(Eldataexp,aes(x=DR1,y=DR2,labels=sampleID))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=Kit),size=2)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")


Hscs_e<-ggplot(Eldataexp,aes(x=DR1,y=DR2,labels=sampleID))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=HSCscore),size=2)+
  scale_color_gradient2(low="white",high = "purple")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

library(gridExtra)
grid.arrange(Gfi_e,Jag1_e)
grid.arrange(Gfi_e,Nfkbia_e,Gata2_e,Kit_e,ncol=2)


## IMPORT MODULES ENRICHED IN JAG1+ HE CELLS FROM CEMITOOLS
cemimood<-read.delim("/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/CEMItools/cemitool_results_orig_vs_noorig/module.tsv")
geneTable

cemimood<-merge(cemimood,geneTable,by.x="genes",by.y="Ensembl_Gene_ID")
cemimood<- cemimood[order(cemimood$modules),]

genecemi<-Eldataexp[,colnames(Eldataexp) %in% cemimood$Associated_Gene_Name]
expcemi<-colnames(genecemi[, which(numcolwise(sum)(genecemi) !=0)])
genecemi<-genecemi[,colnames(genecemi) %in% expcemi]

genecemi<-cbind(Eldataexp[,c(1,11)],genecemi)

row.names(genecemi)<-genecemi$sampleID
genecemi$sampleID<-NULL

colnames(genecemi)

name_cemi<-as.data.frame(colnames(Eldataexp[,colnames(Eldataexp) %in% expcemi]))
colnames(name_cemi)<-"gene"
row.names(name_cemi)<-name_cemi$gene


name_cemi<-merge(name_cemi,cemimood,by.x="gene",by.y="Associated_Gene_Name", alL.x=TRUE,sort = FALSE)

nrow(name_cemi)

library(plyr)

library(heatmaply)

heatmaply(genecemi,
          plot_method = "plotly",
          ColSideColors = name_cemi$modules,
          k_row = 4,
          k_col = 4,
          fontsize_row = 5,
          fontsize_col = 5)

row.names(name_cemi)<-colnames(genecemi[,2:ncol(genecemi)])

# pheatmap

pheatcemi<-t(genecemi[,-1])
row.names(name_cemi)<-name_cemi$gene
name_cemi$annot<-"gene"

annotcol<-data.frame(cluster=genecemi[,c(1)])
names(annotcol)
row.names(annotcol)<-row.names(genecemi)

library(pheatmap)

annotrow<-as.data.frame(name_cemi[,c(3)])
row.names(annotrow)<-name_cemi$gene

pheatmap(pheatcemi,
        annotation_row = annotrow,
        annotation_col = annotcol,
         fontsize =5,
        cutree_rows = 5,
        cutree_cols = 2)




