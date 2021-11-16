## Load libraries

library(devtools)
library(Rtsne)
library(DESeq2)
library(scran)
require(knitr)
library(scRNAseq)
#library(bglab)
library(ggrepel)
library(gridExtra)
library(pheatmap)

#Set seed to reproduce the data every time
set.seed(100)

# Set working directory
setwd("/Users/yguillen/Desktop/temp/scRNA_cam/DLL4_Porcheri/")

counts_dl<-read.table("GSE124477_HTSeq_counts.txt",header=TRUE,row.names=1,stringsAsFactors = FALSE)
counts_dl<-as.matrix(counts_dl)
colnames(counts_dl)
rownames(counts_dl)
nrow(counts_dl)
dim(counts_dl)

#metadata
metadl<-read.delim("/Users/yguillen/Desktop/temp/scRNA_cam/DLL4_Porcheri/metadata.txt",sep="\t",stringsAsFactors = FALSE,header = FALSE)
colnames(metadl)<-c("ID","pop","treatment","Cluster_ICG")
dim(metadl)
row.names(metadl)<-metadl$ID

colnames(counts_dl)
row.names(metadl)

dim(counts_dl)
dim(metadl)


metadl<- metadl[order(metadl$ID),]


genetabdl<-read.delim("gTable2.txt",header=FALSE)
colnames(genetabdl)<-c("Ensembl_Gene_ID","Associated_Gene_Name","biotype")
nrow(genetabdl)



# Check dimensions
dim(genetabdl)
dim(counts_dl)

genetabdl<-genetabdl[!duplicated(genetabdl$Ensembl_Gene_ID),]


# Loading data into the single cell dataset (SCD) object which is required for the bglab package
colnames(metadl)
head(rownames(metadl))
head(colnames(counts_dl))

#Probably cells lost in quality control (21 in total). It should be 0!
setdiff(colnames(counts_dl),rownames(metadl))


# Cell cycle assignment using cyclone.
# Use of pre-defined classifier to assign cells into their cell cycle
# phases. This classifier was constructed from a training data set
# by identifying pairs of genes where the difference in expression
# within each pair changed sign accross phases.

# classifiers for human and mouse are provided. For other systems they can be constructed using sandbag package
hs.pairs <- readRDS(system.file("exdata","human_cycle_markers.rds",package = "scran"))
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

scdl <- SingleCellExperiment(assays = list(counts = counts_dl))
counts_dl <- assay(scdl, "counts")
libsizes <- colSums(counts_dl)
size.factors <- libsizes/mean(libsizes)
logcounts(scdl) <- log2(t(t(counts_dl)/size.factors) + 1)
assayNames(scdl)

expall_dl <- assay(scdl, "logcounts")



#for doing cyclone we need to change the gene names ENS in this case
#Do cyclone before filterGene(scd)
# This is because the lack of expression of particular genes can provide 
# some information about the cell cycle.
cycdl<-cyclone(expall_dl,mm.pairs)
cyctabdl<-as.data.frame(cycdl$phases)
row.names(cyctabdl)<-colnames(expall_dl)
colnames(cyctabdl)<-"Phase"

cycdlG1<-as.data.frame(cycdl$scores$G1)
row.names(cycdlG1)<-colnames(expall_dl)
colnames(cycdlG1)<-"G1"

cycdlG2M<-as.data.frame(cycdl$scores$G2M)
row.names(cycdlG2M)<-colnames(expall_dl)
colnames(cycdlG2M)<-"G2M"

cycdlS<-as.data.frame(cycdl$scores$S)
row.names(cycdlS)<-colnames(expall_dl)
colnames(cycdlS)<-"S"

write.table(cyctabdl,"cyclone_dll4_scRNAseq_R.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
write.table(cycdlG1,"G1_score_cyclone_dll4.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
write.table(cycdlG2M,"G2M_score_cyclone_dll4.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
write.table(cycdlS,"S_score_cyclone_dll4.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)


### IMPORT SCANPY RESULTS

setwd("/Users/yguillen/Desktop/temp/scRNA_cam/DLL4_Porcheri/cellbrowser_dll4_output/")

UMAPDl<-read.delim("umap_coords.tsv")
colnames(UMAPDl)<-c("sampleID","UMAP_C1","UMAP_C2")
row.names(UMAPDl)<-UMAPDl$sampleID

exprmatDl<-read.delim("exprMatrix.tsv")
row.names(exprmatDl)<-exprmatDl$gene

metascDl<-read.delim("meta.tsv")
colnames(metascDl)<-c("sampleID","Louvain_cluster","Mito_perc","expr_genes","UMICount")
row.names(metascDl)<-metascDl$sampleID
metascDl$sampleID<-NULL


## Results UMAP Louvain scanpy

scanpyDl<-merge(UMAPDl,metascDl,by="row.names")
colnames(metadl)[1]<-"sampleID"
scanpyDl<-merge(metadl,scanpyDl,by="sampleID")
row.names(scanpyDl)<-scanpyDl$sampleID
scanpyDl$Row.names=NULL


scanpyp<-ggplot(scanpyDl,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=4,alpha=0.5)+
  geom_point(aes(color=as.factor(Louvain_cluster)),size=3)+
  #geom_point(aes(color=treatment),size=3)+
  scale_color_brewer(palette="Spectral")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

scanpypigc<-ggplot(scanpyDl,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=4,alpha=0.5)+
  geom_point(aes(color=Cluster_ICG),size=3)+
  #geom_point(aes(color=treatment),size=3)+
  #scale_color_brewer(palette="Spectral")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

picg<-ggplot(scanpyDl,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=4,alpha=0.5)+
  geom_point(aes(color=Cluster_ICG),size=3)+
  #geom_point(aes(color=treatment),size=3)+
  scale_color_brewer(palette="Spectral")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

treatp<-ggplot(scanpyDl,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=5,alpha=0.5)+
  geom_point(aes(color=treatment),size=4)+
  #geom_point(aes(color=treatment),size=3)+
  scale_color_brewer(palette="Paired")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")



# Add expression data from our scanpy expression results from htseq
dim(exprmatDl)

t_exprmatDl<-data.frame(t(exprmatDl))
t_exprmatDl<-t_exprmatDl[-1,]
class(t_exprmatDl)

ix <- 1:ncol(t_exprmatDl)
t_exprmatDl[ix] <- lapply(t_exprmatDl[ix], as.character)
t_exprmatDl[ix] <- lapply(t_exprmatDl[ix], as.numeric)

row.names(t_exprmatDl)
colnames(t_exprmatDl)<-gsub('ENS.*\\.','',colnames(t_exprmatDl))

scanpyDl$sampleID<-as.character(scanpyDl$sampleID)

#merge
Dldataexp<-merge(scanpyDl,t_exprmatDl,by.x="sampleID",by.y="row.names")

Gfi_e<-ggplot(Dldataexp,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=5)+
  geom_point(aes(color=Gfi1),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

Jag1_e<-ggplot(Dldataexp,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=5)+
  geom_point(aes(color=Jag1),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

Notch1_e<-ggplot(Dldataexp,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=5)+
  geom_point(aes(color=Notch1),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

Hes1_e<-ggplot(Dldataexp,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=5)+
  geom_point(aes(color=Hes1),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

grid.arrange(scanpyp,treatp,
             scanpypigc,Gfi_e,
             Jag1_e,Notch1_e,ncol=2)

grid.arrange(scanpyp,treatp,
             Hes1_e,ncol=3)

## Number of control vs Dll4 treatment per cluster
chit<-as.matrix(table(Dldataexp$treatment,Dldataexp$Louvain_cluster))
chisq.test(chit)

df_treat<-melt(t(prop.table(table(Dldataexp$treatment,Dldataexp$Louvain_cluster),2)))
colnames(df_treat)<-c("Louvain_cluster","treatment","propcells")

ptreat<-ggplot(df_treat,aes(x=as.factor(Louvain_cluster),y=propcells))+
  geom_bar(stat="identity",aes(fill=treatment),color="black")+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  theme(legend.position = "bottom")


grid.arrange(scanpyp,treatp,ptreat,ncol=2)
  

ggplot(Dldataexp,aes(x=treatment,y=Jag1))+
  geom_violin(aes(fill=treatment),alpha=0.2)+
  geom_jitter(aes(color=treatment))+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  facet_wrap(~as.factor(Louvain_cluster),scales="free_y",ncol=5)+
  theme_bw()


### Genes that correlate with Jag1 or Notch1 protein display in Roshana's data, all clusters
Jag1cor<-corhey2_sig
Jag1cor<-Jag1cor[with(Jag1cor,order(Jag1cor$Pval, Jag1cor$Rho)),]


Jag1_sig<-Dldataexp[,colnames(Dldataexp) %in% c("sampleID",Jag1cor$Gene) ]
row.names(Jag1_sig)<-Jag1_sig$sampleID
Jag1_sig$sampleID<-NULL
Jag1_sig<-t(Jag1_sig)

metaJag1<-subset(Dldataexp,select=c(sampleID,Louvain_cluster,treatment))
row.names(metaJag1)<-metaJag1$sampleID
metaJag1$sampleID<-NULL

metaJag1$Louvain_cluster<-as.factor(metaJag1$Louvain_cluster)

Jag1cor_pos<-subset(Jag1cor,Jag1cor$Rho>0)
Jag1cor_pos$rho<-"positive"
Jag1cor_pos<-Jag1cor_pos[,c("Gene","rho")]

Jag1cor_neg<-subset(Jag1cor,Jag1cor$Rho<0)
Jag1cor_neg$rho<-"negative"
Jag1cor_neg<-Jag1cor_neg[,c("Gene","rho")]

Jag1coranot<-rbind(Jag1cor_pos,Jag1cor_neg)
Jag1coranot$Gene<-NULL


# Scale by each element
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

Jag1_sig_norm <- t(apply(Jag1_sig, 1, cal_z_score))


#order by treatment (no patterns..) cluster_cols=FALSE

Jag1_sig_norm<-Jag1_sig_norm[, order(metaJag1$Louvain_cluster)]
Jag1_sig_norm<-Jag1_sig_norm[,row.names(metaJag1[order(metaJag1$Louvain_cluster),])]


library(pheatmap)
pheatmap(Jag1_sig_norm,
         annotation_col = metaJag1,
         annotation_row = Jag1coranot,
         cutree_rows = 2,
         cutree_cols = 3,
         fontsize =5,
         Colv=FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         dendrogram="row")

