
#### Install packages (do only once)

#install.packages("BiocManager")

BiocManager::install(c("DESeq2",
                       "devtools"))
BiocManager::install(c("Rtsne",
                       "ChIPseeker"))
BiocManager::install(c("scran"))
BiocManager::install("scRNAseq")
BiocManager::install("umap")
BiocManager::install("rgl")
install_github("wjawaid/bglab")

## Load libraries
#library(rgl)
library(devtools)
library(Rtsne)
library(DESeq2)
library(scran)
require(knitr)
library(scRNAseq)
library(bglab)
library(ggrepel)
library(gridExtra)
library(reshape2)
library(heatmaply)

#Set seed to reproduce the data every time
set.seed(100)

setwd("/Volumes/cancer/Rosh_scRNA/htseq/")

#counts<-read.table("counts.txt",header=TRUE,row.names=1,sep=" ",stringsAsFactors = FALSE)
counts_agm<-read.table("htseq_counts_genes.txt",header=TRUE,row.names=1,stringsAsFactors = FALSE)
counts_agm<-as.matrix(counts_agm)
colnames(counts_agm)
rownames(counts_agm)
nrow(counts_agm)
dim(counts_agm)


sce_agm <- SingleCellExperiment(assays = list(counts = counts_agm))
counts_agm <- assay(sce_agm, "counts")
libsizes <- colSums(counts_agm)
size.factors <- libsizes/mean(libsizes)
logcounts(sce_agm) <- log2(t(t(counts_agm)/size.factors) + 1)
assayNames(sce_agm)

expall_agm <- assay(sce_agm, "logcounts")



# Cell cycle assignment using cyclone.
# Use of pre-defined classifier to assign cells into their cell cycle
# phases. This classifier was constructed from a training data set
# by identifying pairs of genes where the difference in expression
# within each pair changed sign accross phases.

# classifiers for human and mouse are provided. For other systems they can be constructed using sandbag package
mm.pairs <- readRDS("/Volumes/cancer/Rosh_scRNA/scanpy/mouse_cycle_markers.rds")

#for doing cyclone we need to change the gene names ENS in this case
#Do cyclone before filterGene(scd)
# This is because the lack of expression of particular genes can provide 
# some information about the cell cycle.
cyc_agm<-cyclone(expall_agm,mm.pairs)

cyctab_agm<-as.data.frame(cyc_agm$phases)
row.names(cyctab_agm)<-colnames(expall_agm)
colnames(cyctab_agm)<-"Phase"

cycG1<-as.data.frame(cyc_agm$scores$G1)
row.names(cycG1)<-colnames(expall_agm)
colnames(cycG1)<-"G1"

cycG2M<-as.data.frame(cyc_agm$scores$G2M)
row.names(cycG2M)<-colnames(expall_agm)
colnames(cycG2M)<-"G2M"

cycS<-as.data.frame(cyc_agm$scores$S)
row.names(cycS)<-colnames(expall_agm)
colnames(cycS)<-"S"

write.table(cyctab_agm,"/Volumes/cancer/Rosh_scRNA/scanpy/cyclone_scRNAseq.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
write.table(cycG1,"/Volumes/cancer/Rosh_scRNA/scanpy/G1_score_cyclone.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
write.table(cycG2M,"/Volumes/cancer/Rosh_scRNA/scanpy/G2M_score_cyclone.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)
write.table(cycS,"/Volumes/cancer/Rosh_scRNA/scanpy/S_score_cyclone.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)


## Generate FACS index data ordered
#metadata all cells, all batches
cell_metadata<-read.delim('/Volumes/cancer/Rosh_scRNA/htseq/cell_metadata.txt',header=FALSE)
colnames(cell_metadata)<-c("ID","Timepoint","Population","Batch","Batch_name","Position")
row.names(cell_metadata)<-cell_metadata$ID

cell_metadata$Position<-toupper(cell_metadata$Position)
cell_metadata$Position<-gsub('PLATE','Plate',cell_metadata$Position)

# FACS
restGL<-read.delim('/Volumes/cancer/Rosh_scRNA/scanpy/FACS_data/all_plates.csv',header = TRUE)
colnames(restGL)[1]<-"Position"
restGL<-restGL[,c(1,4:14)]
restGL<-restGL[restGL$Position!="",]


restGL<-restGL[,c(2:6,1,7:ncol(restGL))]
row.names(restGL)<-restGL$ID

# read CD45 marker from Roshana files
facsindex<-read.delim('/Volumes/cancer/Rosh_scRNA/scanpy/FACS_data/all_plates.csv',header = TRUE,sep=",",dec = ".")
facsindex_g<-as.data.frame(apply(facsindex, 2, function(y) gsub(",","", y)))
facsindex_g<-facsindex_g[,c(1,10:18)]
colnames(facsindex_g)<-c("Position","DLL4","NOTCH2","cKIT","GF1","CD41","NOTCH1","CD31","JAG1","CD45")


facsindex_g<-merge(cell_metadata,facsindex_g,all.x=TRUE,by="Position",sort=FALSE)
row.names(facsindex_g)<-facsindex_g$ID
write.table(facsindex_g,'/Volumes/cancer/Rosh_scRNA/scanpy/FACS_data/metadata_facs.txt',quote = FALSE,row.names = TRUE,sep="\t")

#############

######### Import scanpy results

### SCANPY RESULTS WITH ALL CELLS
#setwd("/Volumes/grcmc/YGUILLEN/scRNA_cam/scanpy_run/cellbrowser_output/")
#setwd("/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/cellbrowser_output/")

### SCANPY RESULTS EXCLUDING PRIMORDIAL CELLS
#setwd("/Volumes/grcmc/YGUILLEN/scRNA_cam/scanpy_run/cellbrowser_noclust6_output/")
setwd("/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/cellbrowser_noclust6_output/")

UMAPsc<-read.delim("umap_coords.tsv")
colnames(UMAPsc)<-c("sampleID","UMAP_C1","UMAP_C2")
row.names(UMAPsc)<-UMAPsc$sampleID

exprmat<-read.delim("exprMatrix.tsv")
row.names(exprmat)<-exprmat$gene

metasc<-read.delim("meta.tsv")
colnames(metasc)<-c("sampleID","Louvain_cluster","Mito_perc","expr_genes","UMICount")
row.names(metasc)<-metasc$sampleID
metasc$sampleID<-NULL

#pseudotime
pdt<-read.delim("/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/pseudotime.csv",sep=",")
row.names(pdt)<-pdt$X
pdt
dim(pdt)

row.names(metasc)
row.names(metasc)
dim(metasc)

metasc<-merge(metasc,pdt,by="row.names")
row.names(metasc)<-metasc$Row.names
metasc$Row.names<-NULL
metasc$X<-NULL

# metadata
meta

#matching clusters Louvain scanpy and UMAP Zaki
metasc_cluster<-subset(metasc,select=Louvain_cluster)
metasc_cluster$exp<-"Louvain"
colnames(metasc_cluster)<-c("cluster","exp")
metasc_cluster$cluster<-as.factor(metasc_cluster$cluster)
metasc_cluster$cell<-row.names(metasc_cluster)

metazaki_clusters$exp<-"Zaki"
colnames(metazaki_clusters)<-c("cluster","exp")
metazaki_clusters$cell<-row.names(metazaki_clusters)

#mito<-subset(metasc,select=c(Louvain_cluster,Mito_perc))
#colnames(mito)<-c("cluster","mito")
#mito$cluster<-as.factor(metasc_cluster$cluster)

metasc_cluster$cell<-row.names(metasc_cluster)


# merge info clusters both metadata
allcluster<-rbind(metazaki_clusters,metasc_cluster)
dim(metazaki_clusters)
dim(metasc_cluster)

require(easyalluvial)
require(tidyverse)

col_vector = c('coral','gold','dodgerblue','darkolivegreen','orange','darkcyan','blue','darksalmon','grey')

p<-alluvial_long(allcluster
                 , key = exp
                 , value = cluster
                 , id = cell
                 , fill_by = 'value'
                 , col_vector_flow = col_vector
                 ,NA_label = 'no_match'
                 , col_vector_value = col_vector
)


p<-p+
  theme_minimal()+
  theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12))
p

ggsave(p,filename ="../../Roshana_scRNA/Plots/scanpy/compare_clusters_Zaki.pdf" )

## Results UMAP Louvain scanpy

scanpy<-merge(UMAPsc,metasc,by="row.names")
row.names(scanpy)<-scanpy$Row.names
scanpy$Row.names=NULL

scanpy<-merge(scanpy,meta,by="row.names")
row.names(scanpy)=scanpy$Row.names
scanpy$Row.names=NULL
colnames(scanpy)[1]<-"sampleID"


#Plot UMAP
# add Zaki clusters info
scanpy_ad<-merge(scanpy,metazaki_clusters,by="row.names",all=TRUE)
row.names(scanpy_ad)<-scanpy_ad$Row.names

library(ggnewscale)

UMAP_comp<-ggplot(scanpy_ad,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=as.factor(cluster)),size=5)+
  scale_color_brewer(palette="Spectral")+
  new_scale("color")+
  geom_point(color="white",size=3)+
  geom_point(aes(color=as.factor(Louvain_cluster)),size=2)+
  scale_color_brewer(palette="Spectral")+
  theme_bw()+
  theme(legend.position = "bottom")

UMAP_comp

ggsave(UMAP_comp,filename = "../../Roshana_scRNA/Plots/scanpy/comp_clusters_UMAP.pdf")

## UMAP with only scanpy info, we can label cells
UMAP_lab<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=4)+
  geom_point(aes(color=as.factor(Louvain_cluster)),size=3)+
  scale_color_brewer(palette="Spectral")+
  geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

UMAP_lab

ggsave(UMAP_lab,filename = "/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/UMAP_lab_filtercells.pdf")

# pseudotime
UMAP_pst<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=4)+
  geom_point(aes(color=dpt_pseudotime),size=3)+
  scale_color_gradient(low="white",high="blue")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

UMAP_pst


UMAPpop<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=5,alpha=0.5)+
  geom_point(aes(color=population),size=4,alpha=0.5)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_manual(values=c("blue","red","green"))+
  geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

UMAPpop

ggsave(UMAPpop,filename = "/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/UMAP_pop.pdf")



#3D plot
#p <- plot_ly(scanpy, x = ~UMAP_C1, y = ~UMAP_C2, z = ~UMAP_C3,color=~as.factor(Louvain_cluster)) %>%
#  add_markers() %>%
#  layout(scene = list(xaxis = list(title = 'UMAP_C1'),
#                      yaxis = list(title = 'UMAP_C2'),
#                      zaxis = list(title = 'UMAP_C3'))
#  )

#p

#p <- plot_ly(scanpy, x = ~UMAP_C1, y = ~UMAP_C2, z = ~UMAP_C3,
#             marker = list(color = ~Notch2,colorscale=c("'#FFE1A2'", '#3565C6'), showscale = TRUE)) %>%
#  add_markers() %>%
#  layout(scene = list(xaxis = list(title = 'UMAP_C1'),
#                      yaxis = list(title = 'UMAP_C2'),
#                      zaxis = list(title = 'UMAP_C3'))
#  )

#p

#Select cells for each of clusters
# Number of clusters
table(scanpy$Louvain_cluster)


###### generating output clusters tables ####

clust0<-row.names(scanpy[scanpy$Louvain_cluster=="0",])
clust0mat<-exprmat[,colnames(exprmat) %in% clust0]
#Remove genes not expressed in any cell
clust0mat<-clust0mat[rowSums(clust0mat)>0, ]
#write output csv
write.table(clust0mat,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/R_tables/cluster0.tsv",sep="\t",quote = FALSE,row.names = TRUE)

clust1<-row.names(scanpy[scanpy$Louvain_cluster=="1",])
clust1mat<-exprmat[,colnames(exprmat) %in% clust1]
#Remove genes not expressed in any cell
clust1mat<-clust1mat[rowSums(clust1mat)>0, ]
#write output csv
write.table(clust1mat,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/R_tables/cluster1.tsv",sep="\t",quote = FALSE,row.names = TRUE)

clust2<-row.names(scanpy[scanpy$Louvain_cluster=="2",])
clust2mat<-exprmat[,colnames(exprmat) %in% clust2]
#Remove genes not expressed in any cell
clust2mat<-clust2mat[rowSums(clust2mat)>0, ]
#write output csv
write.table(clust2mat,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/R_tables/cluster2.tsv",sep="\t",quote = FALSE,row.names = TRUE)


clust3<-row.names(scanpy[scanpy$Louvain_cluster=="3",])
clust3mat<-exprmat[,colnames(exprmat) %in% clust3]
#Remove genes not expressed in any cell
clust3mat<-clust3mat[rowSums(clust3mat)>0, ]
#write output csv
write.table(clust3mat,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/R_tables/cluster3.tsv",sep="\t",quote = FALSE,row.names = TRUE)

clust4<-row.names(scanpy[scanpy$Louvain_cluster=="4",])
clust4mat<-exprmat[,colnames(exprmat) %in% clust4]
#Remove genes not expressed in any cell
clust4mat<-clust4mat[rowSums(clust4mat)>0, ]
#write output csv
write.table(clust4mat,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/R_tables/cluster4.tsv",sep="\t",quote = FALSE,row.names = TRUE)


clust5<-row.names(scanpy[scanpy$Louvain_cluster=="5",])
clust5mat<-exprmat[,colnames(exprmat) %in% clust5]
#Remove genes not expressed in any cell
clust5mat<-clust5mat[rowSums(clust5mat)>0, ]
#write output csv
write.table(clust5mat,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/R_tables/cluster5.tsv",sep="\t",quote = FALSE,row.names = TRUE)


#clust6<-row.names(scanpy[scanpy$Louvain_cluster=="6",])
#clust6mat<-exprmat[,colnames(exprmat) %in% clust6]
#Remove genes not expressed in any cell
#clust6mat<-clust6mat[rowSums(clust6mat)>0, ]
#write output csv
#write.table(clust6mat,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/R_tables/cluster6.tsv",sep="\t",quote = FALSE,row.names = TRUE)


# Cells from cluster 2 that resemble HSC (check limits UMAP coordinates)
#clust2_HSC<-row.names(subset(scanpy,(scanpy$UMAP_C1<(-5) & scanpy$UMAP_C2<(-3))))
#clust2_HSCmat<-exprmat[,colnames(exprmat) %in% clust2_HSC]
#Remove genes not expressed in any cell
#clust2_HSCmat<-clust2_HSCmat[rowSums(clust2_HSCmat)>0, ]
#write output csv
#write.table(clust2_HSCmat,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/R_tables/clust2_HSCmat.tsv",sep="\t",quote = FALSE,row.names = TRUE)

#####

## Plotting receptors / thresholds

# Threshold for Notch1 and Notch2 markers according to gating strategy

scanpy$th_NOTCH1<-ifelse(scanpy$Notch1<500,
                        c("neg_NOTCH1"),c("pos_NOTCH1"))

scanpy$th_NOTCH2<-ifelse(scanpy$Notch2<200,
                         c("neg_NOTCH2"),c("pos_NOTCH2"))

scanpy$th_JAG1<-ifelse(scanpy$Jag1<11000,
                         c("neg_JAG1"),c("pos_JAG1"))

scanpy$th_DLL4<-ifelse(scanpy$Dll4<120,
                         c("neg_DLL4"),c("pos_DLL4"))

scanpy$th_GFI1<-ifelse(scanpy$Gfi1<270,
                       c("neg_GFI1"),c("pos_GFI1"))


scanpy$nothall<-paste(scanpy$th_NOTCH1,scanpy$th_NOTCH2,sep="_")
scanpy$recall<-paste(scanpy$th_JAG1,scanpy$th_DLL4,sep="_")
table(scanpy$nothall)
table(scanpy$recall)

## First, scatter plots FACS index

# melt dataset
scanpymelt<-scanpy[,c(4,20:28)]

scanpymelt<-melt(scanpymelt,id.vars = "Louvain_cluster")

scanpymelt$value<-as.numeric(scanpymelt$value)

# histogram distribution markers
pdf(file = "/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/FACS_markers/marker_distribution.pdf")
ggplot(scanpymelt)+
  geom_density(aes(x=value))+
  scale_color_brewer(palette = "Spectral")+
  geom_vline(data=filter(scanpymelt, variable=="Dll4"), aes(xintercept=120),linetype='dashed',color="red") + 
  geom_vline(data=filter(scanpymelt, variable=="Notch2"), aes(xintercept=200),linetype='dashed',color="red" ) + 
  geom_vline(data=filter(scanpymelt, variable=="cKIT"), aes(xintercept=300),linetype='dashed',color="red" ) + 
  geom_vline(data=filter(scanpymelt, variable=="Gfi1"), aes(xintercept=550),linetype='dashed',color="red" ) + 
  geom_vline(data=filter(scanpymelt, variable=="CD41"), aes(xintercept=median((meta[!is.na(meta$CD41),])$CD41)),linetype='dashed',color="white" ) + 
  geom_vline(data=filter(scanpymelt, variable=="Notch1"), aes(xintercept=500),linetype='dashed',color="red" ) + 
  geom_vline(data=filter(scanpymelt, variable=="CD31"), aes(xintercept=1000),linetype='dashed',color="red" ) + 
  geom_vline(data=filter(scanpymelt, variable=="Jag1"), aes(xintercept=11000),linetype='dashed',color="red" ) + 
  geom_vline(data=filter(scanpymelt, variable=="CD45"), aes(xintercept=median((meta[!is.na(meta$CD45),])$CD45)),linetype='dashed',color="red" ) + 
  facet_wrap(~variable,scales="free")+
  theme_minimal()+
  theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())

dev.off()

pdf(file = "/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/FACS_markers/marker_distribution_per_cluster.pdf")
ggplot(scanpymelt)+
  geom_density(aes(x=value,color=as.factor(Louvain_cluster)))+
  scale_color_brewer(palette = "Spectral")+
  geom_vline(data=filter(scanpymelt, variable=="Dll4"), aes(xintercept=120),linetype='dashed',color="black") + 
  geom_vline(data=filter(scanpymelt, variable=="Notch2"), aes(xintercept=200),linetype='dashed',color="black" ) + 
  geom_vline(data=filter(scanpymelt, variable=="cKIT"), aes(xintercept=300),linetype='dashed',color="black" ) + 
  geom_vline(data=filter(scanpymelt, variable=="Gfi1"), aes(xintercept=550),linetype='dashed',color="black" ) + 
  #geom_vline(data=filter(scanpymelt, variable=="CD41"), aes(xintercept=median((meta[!is.na(meta$CD41),])$CD41)),linetype='dashed',color="white" ) + 
  geom_vline(data=filter(scanpymelt, variable=="Notch1"), aes(xintercept=500),linetype='dashed',color="black" ) + 
  geom_vline(data=filter(scanpymelt, variable=="CD31"), aes(xintercept=1000),linetype='dashed',color="black" ) + 
  geom_vline(data=filter(scanpymelt, variable=="Jag1"), aes(xintercept=11000),linetype='dashed',color="black" ) + 
  #geom_vline(data=filter(scanpymelt, variable=="CD45"), aes(xintercept=median((meta[!is.na(meta$CD45),])$CD45)),linetype='dashed',color="white" ) + 
  facet_wrap(~variable,scales="free")+
  theme_minimal()+
  theme(legend.position = "bottom",axis.text.x = element_blank(),axis.ticks.x = element_blank())

dev.off()

library(ggnewscale)

facs<-ggplot(scanpy,aes(x=cKIT,y=Gfi1))+
  geom_point(data=scanpy[scanpy$sampleID %in% c(clust4),],color="black",size=5)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust6),],color="black",size=3)+
  geom_point(color="black",size=4)+
  geom_point(aes(color=as.factor(Louvain_cluster)),size=3)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust4,clust6),],aes(color=as.factor(Louvain_cluster)),alpha=1,size=1)+
  scale_color_brewer(palette = "Spectral")+
  geom_vline(xintercept = 200,linetype="dashed")+
  geom_hline(yintercept = 270,linetype="dashed")+
  theme_bw()+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')+
  new_scale("color")+
  geom_label_repel(aes(label=as.factor(Louvain_cluster)),size=5,
                   data=(scanpy[scanpy$th_NOTCH2=="pos_NOTCH2",]))+
  scale_color_brewer(palette = "Spectral")+
  ggtitle("Notch2 positive cells")
  #ylim(c(-1800,10000))
  #xlim(c(-1100,10000))
facs


ggsave(facs,filename="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/FACS_clusters.pdf")


cd31lim<-ggplot(scanpy,aes(x=FSC.A/1000,y=CD31))+
  geom_point(data=scanpy[scanpy$sampleID %in% c(clust3),],color="red",size=3)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust3) & scanpy$UMAP_C2 < 0,],color="red",size=3)+
  geom_point(color="black",alpha=0.5,size=2)+
  geom_point(color="white",alpha=1,size=1)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust4,clust6),],aes(color=as.factor(Louvain_cluster)),alpha=1,size=1)+
  geom_rug(sides="r")+
  geom_rug(data=scanpy[scanpy$sampleID %in% c(clust3),],color="red",sides="r")+
  scale_color_brewer(palette = "Spectral")+
  geom_hline(yintercept = 1000,linetype="dashed")+
  scale_y_continuous(trans='log10')+
  theme_bw()+
  theme(legend.position = "bottom",legend.title = element_blank())+
  new_scale("color")+
  #geom_label_repel(aes(label=population,color=nothall),
  #                 data=(scanpy[(scanpy$Louvain_cluster=="4" | scanpy$Louvain_cluster=="6") &
  #                                (scanpy$th_NOTCH1=="pos_NOTCH1" | scanpy$th_NOTCH2=="pos_NOTCH2"),]),
  #                 alpha=0.8)+
  scale_color_brewer(palette = "Set2")+
  theme(legend.title = element_blank())

cd31lim


ggsave(cd31lim,filename="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/FACS_CD31_cluster3.pdf")

gfi1clust1<-ggplot(scanpy,aes(x=cKIT,y=Gfi1))+
  geom_point(data=scanpy[scanpy$sampleID %in% c(clust3),],color="red",size=3)+
  geom_point(color="black",alpha=0.5,size=2)+
  geom_point(color="white",alpha=1,size=1)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust4,clust6),],aes(color=as.factor(Louvain_cluster)),alpha=1,size=1)+
  scale_color_brewer(palette = "Spectral")+
  geom_vline(xintercept = 200,linetype="dashed")+
  geom_hline(yintercept = 270,linetype="dashed")+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  theme_bw()+
  theme(legend.position = "bottom",legend.title = element_blank())+
  new_scale("color")+
  #geom_label_repel(aes(label=population,color=nothall),
  #                 data=(scanpy[(scanpy$Louvain_cluster=="4" | scanpy$Louvain_cluster=="6") &
  #                                (scanpy$th_NOTCH1=="pos_NOTCH1" | scanpy$th_NOTCH2=="pos_NOTCH2"),]),
  #                 alpha=0.8)+
  scale_color_brewer(palette = "Set2")+
  theme(legend.title = element_blank())
gfi1clust1

ggsave(gfi1clust1,filename="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/FACS_kitgfi_cluster3.pdf")


facsnotchpop<-ggplot(scanpy,aes(x=Notch2,y=Notch1))+
  geom_point(data=scanpy[scanpy$sampleID %in% c(clust4),],color="black",size=5)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust6),],color="black",size=3)+
  geom_point(color="black",size=4)+
  geom_point(aes(color=population),alpha=1,size=3)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust4,clust6),],aes(color=as.factor(Louvain_cluster)),alpha=1,size=1)+
  scale_color_brewer(palette = "Spectral")+
  geom_vline(xintercept = 200,linetype="dashed")+
  geom_hline(yintercept = 500,linetype="dashed")+
  theme_bw()+
  theme(legend.position = "left",legend.title = element_blank())
  #xlim(c(NA,10000))+
  #scale_x_continuous(trans='log10') +
  #ylim(c(NA,10000))+
  #scale_y_continuous(trans='log10')

facsnotchpop


facsnotch<-ggplot(scanpy,aes(x=Notch2,y=Notch1))+
  geom_point(data=scanpy[scanpy$sampleID %in% c(clust4),],color="black",size=5)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust6),],color="black",size=3)+
  geom_point(color="black",size=4)+
  geom_point(aes(color=as.factor(Louvain_cluster)),alpha=1,size=3)+
  #geom_point(data=scanpy[scanpy$sampleID %in% c(clust4,clust6),],aes(color=as.factor(Louvain_cluster)),alpha=1,size=1)+
  scale_color_brewer(palette = "Spectral")+
  geom_vline(xintercept = 200,linetype="dashed")+
  geom_hline(yintercept = 500,linetype="dashed")+
  theme_bw()+
  theme(legend.position = "left",legend.title = element_blank())
#xlim(c(NA,10000))+
#scale_x_continuous(trans='log10') +
#ylim(c(NA,10000))+
#scale_y_continuous(trans='log10')

facsnotch

gridfacs<-grid.arrange(facsnotch,facsnotchpop,ncol=1)

ggsave(gridfacs,filename="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/FACS_NOTCH_clusters.pdf")



comp_thres<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=th_NOTCH1),alpha=1,size=5)+
  #scale_colour_gradient(low = "white", high = "red")+
  scale_colour_manual(values=c("coral","dodgerblue"))+
  new_scale("color")+
  geom_point(alpha=1,size=3,color="black")+
  new_scale("color")+
  geom_point(aes(color=th_NOTCH2),alpha=1,size=2)+
  #scale_colour_gradient(low = "white", high = "blue")+
  scale_colour_manual(values=c("coral","dodgerblue"))+
  theme_bw()+
  theme(legend.position = "bottom")

comp_thres

comp_thres1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,label=sampleID))+
  #geom_point(data=(subset(scanpy,scanpy$Louvain_cluster=="4" & (scanpy$nothall=="pos_NOTCH1_pos_NOTCH2" | scanpy$nothall=="neg_NOTCH1_pos_NOTCH2"))),alpha=1,size=3,color="red")+
  geom_point(alpha=1,size=2,color="black")+
  geom_point(aes(color=nothall),size=1)+
  scale_colour_manual(values=c("grey","purple","dodgerblue","coral","darkolivegreen"))+
  theme_bw()+
  theme(legend.position = "bottom")+
  new_scale("color")+
  stat_ellipse(aes(color=as.factor(Louvain_cluster)),level=0.95,size=1)+
  scale_color_manual(values=c("black","black","black","black","black","black","black"))+
  new_scale("color")+
  stat_ellipse(aes(color=as.factor(Louvain_cluster)),level=0.95,size=0.5)+
  scale_color_brewer(palette="Spectral")+
  geom_label_repel(aes(label=Plate_Well),color="black",data=subset(scanpy,scanpy$Louvain_cluster=="4" & (scanpy$nothall=="pos_NOTCH1_pos_NOTCH2" | scanpy$nothall=="neg_NOTCH1_pos_NOTCH2")), fontface = "italic",size=3)+
  ggtitle("Receptors")

comp_thres1

comp_thres2<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(alpha=1,size=2,color="black")+
  geom_point(aes(color=recall),size=1)+
  scale_colour_manual(values=c("grey","purple","dodgerblue","coral","darkolivegreen"))+
  theme_bw()+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  theme()+
  theme(legend.position = "bottom")+
  new_scale("color")+
  stat_ellipse(aes(color=as.factor(Louvain_cluster)),level=0.95,size=1)+
  scale_color_manual(values=c("black","black","black","black","black","black","black"))+
  new_scale("color")+
  stat_ellipse(aes(color=as.factor(Louvain_cluster)),level=0.95,size=0.5)+
  scale_color_brewer(palette="Spectral")+
  ggtitle("Ligands")

comp_thres2

library(gridExtra)

comp_thres<-grid.arrange(comp_thres1,comp_thres2,ncol=1)

ggsave(comp_thres,file="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/Plotscomp_thres_notch1_2_rec.pdf")


# CYCLONE data cell cycle phase score
#cyclondata<-read.delim("/Volumes/grcmc/YGUILLEN/scRNA_cam/scanpy_run/cyclone_R/cyclone_scRNAseq_R.txt",header = FALSE)
cyclondata<-read.delim("/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/cyclone_R/cyclone_scRNAseq_R.txt",header = FALSE)
colnames(cyclondata)<-c("sampleID","cycle")

scanpy_cyc<-merge(scanpy,cyclondata,by="sampleID")

cycplot<-ggplot(scanpy_cyc,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=5)+
  geom_point(aes(color=cycle),size=4)+
  scale_color_brewer(palette="Pastel1")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

cycplot

ggsave(cycplot,file="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/cyclone_scanpy.pdf")


notch1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=th_NOTCH1),alpha=1,size=3)+
  #scale_colour_gradient(low = "white", high = "red")+
  scale_colour_manual(values=c("coral","dodgerblue"))+
  theme_bw()+
  theme(legend.position = "bottom")

notch1

notch2<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=th_NOTCH2),alpha=1,size=3)+
  #scale_colour_gradient(low = "white", high = "red")+
  scale_colour_manual(values=c("coral","dodgerblue"))+
  theme_bw()+
  theme(legend.position = "bottom")

notch2

dll4<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=th_DLL4),alpha=1,size=3)+
  #scale_colour_gradient(low = "white", high = "red")+
  scale_colour_manual(values=c("coral","dodgerblue"))+
  theme_bw()+
  theme(legend.position = "bottom")

dll4

Jag1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=th_JAG1),alpha=1,size=3)+
  #scale_colour_gradient(low = "white", high = "red")+
  scale_colour_manual(values=c("coral","dodgerblue"))+
  theme_bw()+
  theme(legend.position = "bottom")

Jag1


notchmix<-grid.arrange(notch1,notch2,ncol=2)

repmix<-grid.arrange(dll4,Jag1,ncol=2)


ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(data=scanpy[scanpy$th_NOTCH1=="pos_NOTCH1" & scanpy$th_JAG1=="pos_JAG1",],color="red",alpha=1,size=4)+
  geom_point(data=scanpy[scanpy$th_NOTCH1=="pos_NOTCH1" & scanpy$th_DLL4=="pos_DLL4",],color="dodgerblue",alpha=1,size=3)+
  #scale_colour_gradient(low = "white", high = "red")+
 # scale_colour_manual(values=c("coral","dodgerblue"))+
  theme_bw()+
  theme(legend.position = "bottom")



### Flow diagrams receptors - ligands

rec<-subset(scanpy,select=c(recall,sampleID))
colnames(rec)<-c("par","sampleID")
rec$exp<-"receptor"

lig<-subset(scanpy,select=c(nothall,sampleID))
colnames(lig)<-c("par","sampleID")
lig$exp<-"ligand"

cellc<-subset(scanpy,select=c(Louvain_cluster,sampleID))
colnames(cellc)<-c("par","sampleID")
cellc$exp<-"cluster"

reclig<-rbind(cellc,rec,lig)


col_vector = c('coral','gold','dodgerblue','darkolivegreen','orange','darkcyan','blue','darksalmon','grey')

p<-alluvial_long(reclig
                 , key = exp
                 , value = par
                 , id = sampleID
                 , fill_by = 'value'
                 , col_vector_flow = col_vector
                 ,NA_label = 'no_match'
                 , col_vector_value = col_vector,
                 stratum_label_size = 2
)


p<-p+
  theme_minimal()+
  theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12))
p

ggsave(p,filename ="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/recep_ligand_match.pdf" )

#####
###### MARKERS CORRELATION ###


### Correlation matrix for markers (protein levels), independent for each cluster
# Check number columns
corscan_0<-scanpy[scanpy$Louvain_cluster==0,]
corscan_0<-corscan_0[,c(19:27)]

corscan_1<-scanpy[scanpy$Louvain_cluster==1,]
corscan_1<-corscan_1[,c(19:27)]

corscan_2<-scanpy[scanpy$Louvain_cluster==2,]
corscan_2<-corscan_2[,c(19:27)]

corscan_3<-scanpy[scanpy$Louvain_cluster==3,]
corscan_3<-corscan_3[,c(19:27)]

corscan_4<-scanpy[scanpy$Louvain_cluster==4,]
corscan_4<-corscan_4[,c(19:27)]

corscan_5<-scanpy[scanpy$Louvain_cluster==5,]
corscan_5<-corscan_5[,c(19:27)]

#corscan_6<-scanpy[scanpy$Louvain_cluster==6,]
#corscan_6<-corscan_6[,c(19:27)]


#correlate using Hmisc, with p-values
library(Hmisc)
library(corrplot)
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


# All cells, no clusters
clustall_cor<-rcorr(as.matrix(scanpy[,c(19:27)]))
#use flattenCorrMatrix to rearrange tables pvalues and R coefficients
flattenCorrMatrix(clustall_cor$r, clustall_cor$P)

#setwd("/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/FACS_markers/")

pdf(file = "corrplot_allcells.pdf")
corrplot(clustall_cor$r, order="hclust", 
         p.mat = clustall_cor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("All cells",line = 0,adj=0)
dev.off()

# Cluster 0
clust0_cor<-rcorr(as.matrix(corscan_0))
#use flattenCorrMatrix to rearrange tables pvalues and R coefficients
flattenCorrMatrix(clust0_cor$r, clust0_cor$P)

pdf(file = "corrplot_cluster0.pdf")
corrplot(clust0_cor$r, order="hclust", 
         p.mat = clust0_cor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 0",line = 0,adj=0)
dev.off()

# Cluster 1
pdf(file = "corrplot_cluster1.pdf")
clust1_cor<-rcorr(as.matrix(corscan_1))
#use flattenCorrMatrix to rearrange tables pvalues and R coefficients
flattenCorrMatrix(clust1_cor$r, clust1_cor$P)

corrplot(clust1_cor$r, order="hclust", 
         p.mat = clust1_cor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 1",line = 0,adj=0)
dev.off()

# Cluster 2
pdf(file = "corrplot_cluster2.pdf")
clust2_cor<-rcorr(as.matrix(corscan_2))
#use flattenCorrMatrix to rearrange tables pvalues and R coefficients
flattenCorrMatrix(clust2_cor$r, clust2_cor$P)

corrplot(clust2_cor$r, order="hclust", 
         p.mat = clust2_cor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 2",line = 0,adj=0)
dev.off()

# Cluster 3
pdf(file = "corrplot_cluster3.pdf")
clust3_cor<-rcorr(as.matrix(corscan_3))
#use flattenCorrMatrix to rearrange tables pvalues and R coefficients
flattenCorrMatrix(clust3_cor$r, clust3_cor$P)

corrplot(clust3_cor$r, order="hclust", 
         p.mat = clust3_cor$P, sig.level = 0.05, insig = "blank",addrect=3)

title("Cluster 3",line = 0,adj=0)
dev.off()

# Cluster 4
pdf(file = "corrplot_cluster4.pdf")
clust4_cor<-rcorr(as.matrix(corscan_4))
#use flattenCorrMatrix to rearrange tables pvalues and R coefficients
flattenCorrMatrix(clust4_cor$r, clust4_cor$P)

corrplot(clust4_cor$r, order="hclust", 
         p.mat = clust4_cor$P, sig.level = 0.05, insig = "blank",addrect=3)
title("Cluster 4",line = 0,adj=0)
dev.off()

# Cluster 5
pdf(file = "corrplot_cluster5.pdf")

clust5_cor<-rcorr(as.matrix(corscan_5))
#use flattenCorrMatrix to rearrange tables pvalues and R coefficients
flattenCorrMatrix(clust5_cor$r, clust5_cor$P)

corrplot(clust5_cor$r, order="hclust", 
         p.mat = clust5_cor$P, sig.level = 0.05, insig = "blank",addrect=3)
title("Cluster 5",line = 0,adj=0)
dev.off()

# Cluster 6
#pdf(file = "corrplot_cluster6.pdf")

#clust6_cor<-rcorr(as.matrix(corscan_6))
#use flattenCorrMatrix to rearrange tables pvalues and R coefficients
#flattenCorrMatrix(clust6_cor$r, clust6_cor$P)

#corrplot(clust6_cor$r, order="hclust", 
#         p.mat = clust6_cor$P, sig.level = 0.05, insig = "blank",addrect = 3)
#title("Cluster 6",line = 0,adj=0)
#dev.off()


#correlate using corrr package
#BiocManager::install("corrr")

library(corrr)
tabcor_0<-correlate(corscan_0)
tabcor_1<-correlate(corscan_1)
tabcor_2<-correlate(corscan_2)
tabcor_3<-correlate(corscan_3)
tabcor_4<-correlate(corscan_4)
tabcor_5<-correlate(corscan_5)
tabcor_6<-correlate(corscan_6)

library(dplyr)
# we can filter
#tabcor %>% filter(Dll4 > .5)

# fashion table
fashion(tabcor_0)
fashion(tabcor_1)
fashion(tabcor_2)
fashion(tabcor_3)
fashion(tabcor_4)
fashion(tabcor_5)
fashion(tabcor_6)

# plot correlation with dots
p0<-tabcor_0 %>%
  rearrange(method = "MDS", absolute = FALSE) %>%
  shave() %>%
  rplot()

p1<-tabcor_1 %>%
  rearrange(method = "MDS", absolute = FALSE) %>%
  shave() %>%
  rplot()

p2<-tabcor_2 %>%
  rearrange(method = "MDS", absolute = FALSE) %>%
  shave() %>%
  rplot()

p3<-tabcor_3 %>%
  rearrange(method = "MDS", absolute = FALSE) %>%
  shave() %>%
  rplot()

p4<-tabcor_4 %>%
  rearrange(method = "MDS", absolute = FALSE) %>%
  shave() %>%
  rplot()

p5<-tabcor_5 %>%
  rearrange(method = "MDS", absolute = FALSE) %>%
  shave() %>%
  rplot()

#p6<-tabcor_6 %>%
#  rearrange(method = "MDS", absolute = FALSE) %>%
#  shave() %>%
#  rplot()

grid.arrange(p0,p1,p2,p3,p4,p5)

# Plot network
p0_net<-network_plot(tabcor_0,min_cor = 0.3)
p1_net<-network_plot(tabcor_1,min_cor = 0.3,legend = FALSE)
p2_net<-network_plot(tabcor_2,min_cor = 0.3,legend = FALSE)
p3_net<-network_plot(tabcor_3,min_cor = 0.3)
p4_net<-network_plot(tabcor_4,min_cor = 0.3,legend = FALSE)
p5_net<-network_plot(tabcor_5,min_cor = 0.3,legend = FALSE)
#p6_net<-network_plot(tabcor_6,min_cor = 0.3,legend = FALSE)


grid.arrange(p0_net,p1_net,p2_net,p3_net,p4_net,p5_net)
 
# Correlation markers and genes across clusters
#####


### GENE EXPRESSION ##

#setwd("/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/")
setwd("/Users/yguillen/Desktop/temp/scRNA_cam/Roshana_scRNA/Plots/scanpy/")

### Correlation matrix for genes 

exprmat$gene<-gsub('ENS.*\\|','',exprmat$gene)

genelist<-c("Dll4","Notch2","Kit","Gfi1","Itga2b","Notch1","Pecam1","Jag1","Ptprc")
t_geneexp<-t(exprmat[exprmat$gene %in% genelist,])
colnames(t_geneexp)<-t_geneexp[1,]
t_geneexp<-t_geneexp[-1,]

colnames(t_geneexp)<-gsub('^','sc_',colnames(t_geneexp))

scanpy_genes<-merge(scanpy,t_geneexp,by="row.names")
row.names(scanpy_genes)<-scanpy_genes$Row.names
scanpy_genes$Row.names=NULL



# Correlation for each cluster, only genes

# ALL CELLS
gene_corscan<-scanpy_genes[,c(38:46)]

ix <- 1:9
gene_corscan[ix] <- lapply(gene_corscan[ix], as.character)
gene_corscan[ix] <- lapply(gene_corscan[ix], as.numeric)
# Remove genes with 0 counts
gene_corscan<- gene_corscan[,colSums(gene_corscan)!=0]

library(corrr)
genetabcor<-correlate(gene_corscan)
fashion(genetabcor)
genetabcor <- genetabcor[,colSums(is.na(genetabcor))<nrow(genetabcor)]

rplot(genetabcor)
network_plot(genetabcor,min_cor = 0.3)


all_genecor<-rcorr(as.matrix(gene_corscan))
flattenCorrMatrix(all_genecor$r, all_genecor$P)

pdf(file = "genes_corrplot_allcells.pdf")
corrplot(all_genecor$r, order="hclust", 
         p.mat = all_genecor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("All cells",line = 0,adj=0)
dev.off()


# CLUSTER 0
gene_corscan_0<-scanpy_genes[scanpy_genes$Louvain_cluster==0,]
gene_corscan_0<-gene_corscan_0[,c(38:46)]

ix <- 1:9
gene_corscan_0[ix] <- lapply(gene_corscan_0[ix], as.character)
gene_corscan_0[ix] <- lapply(gene_corscan_0[ix], as.numeric)
# Remove genes with 0 counts
gene_corscan_0<- gene_corscan_0[,colSums(gene_corscan_0)!=0]


genetabcor_0<-correlate(gene_corscan_0)
fashion(genetabcor_0)
genetabcor_0 <- genetabcor_0[,colSums(is.na(genetabcor_0))<nrow(genetabcor_0)]

rplot(genetabcor_0)
network_plot(genetabcor_0,min_cor = 0.3)


clust0_genecor<-rcorr(as.matrix(gene_corscan_0))
flattenCorrMatrix(clust0_genecor$r, clust0_genecor$P)

pdf(file = "genes_corrplot_cluster0.pdf")
corrplot(clust0_genecor$r, order="hclust", 
         p.mat = clust0_genecor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 0",line = 0,adj=0)
dev.off()

# CLUSTER 1
gene_corscan_1<-scanpy_genes[scanpy_genes$Louvain_cluster==1,]
gene_corscan_1<-gene_corscan_1[,c(38:46)]

ix <- 1:9
gene_corscan_1[ix] <- lapply(gene_corscan_1[ix], as.character)
gene_corscan_1[ix] <- lapply(gene_corscan_1[ix], as.numeric)
# Remove genes with 0 counts
gene_corscan_1<- gene_corscan_1[,colSums(gene_corscan_1)!=0]


genetabcor_1<-correlate(gene_corscan_1)
fashion(genetabcor_1)
genetabcor_1 <- genetabcor_1[,colSums(is.na(genetabcor_1))<nrow(genetabcor_1)]

rplot(genetabcor_1)
network_plot(genetabcor_1,min_cor = 0.3)

clust1_genecor<-rcorr(as.matrix(gene_corscan_1))
flattenCorrMatrix(clust1_genecor$r, clust1_genecor$P)

pdf(file = "genes_corrplot_cluster1.pdf")
corrplot(clust1_genecor$r, order="hclust", 
         p.mat = clust1_genecor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 1",line = 0,adj=0)
dev.off()

# CLUSTER 2
gene_corscan_2<-scanpy_genes[scanpy_genes$Louvain_cluster==2,]
gene_corscan_2<-gene_corscan_2[,c(38:46)]

ix <- 1:9
gene_corscan_2[ix] <- lapply(gene_corscan_2[ix], as.character)
gene_corscan_2[ix] <- lapply(gene_corscan_2[ix], as.numeric)
# Remove genes with 0 counts
gene_corscan_2<- gene_corscan_2[,colSums(gene_corscan_2)!=0]


genetabcor_2<-correlate(gene_corscan_2)
fashion(genetabcor_2)
genetabcor_2 <- genetabcor_2[,colSums(is.na(genetabcor_2))<nrow(genetabcor_2)]

rplot(genetabcor_2)
network_plot(genetabcor_2,min_cor = 0.3)

clust2_genecor<-rcorr(as.matrix(gene_corscan_2))
flattenCorrMatrix(clust2_genecor$r, clust2_genecor$P)

pdf(file = "genes_corrplot_cluster2.pdf")
corrplot(clust2_genecor$r, order="hclust", 
         p.mat = clust2_genecor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 2",line = 0,adj=0)
dev.off()

# CLUSTER 3
gene_corscan_3<-scanpy_genes[scanpy_genes$Louvain_cluster==3,]
gene_corscan_3<-gene_corscan_3[,c(38:46)]

ix <- 1:9
gene_corscan_3[ix] <- lapply(gene_corscan_3[ix], as.character)
gene_corscan_3[ix] <- lapply(gene_corscan_3[ix], as.numeric)
# Remove genes with 0 counts
gene_corscan_3<- gene_corscan_3[,colSums(gene_corscan_3)!=0]


genetabcor_3<-correlate(gene_corscan_3)
fashion(genetabcor_3)
genetabcor_3 <- genetabcor_3[,colSums(is.na(genetabcor_3))<nrow(genetabcor_3)]

rplot(genetabcor_3)
network_plot(genetabcor_3,min_cor = 0.3)

clust3_genecor<-rcorr(as.matrix(gene_corscan_3))
flattenCorrMatrix(clust3_genecor$r, clust3_genecor$P)

pdf(file = "genes_corrplot_cluster3.pdf")
corrplot(clust3_genecor$r, order="hclust", 
         p.mat = clust3_genecor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 3",line = 0,adj=0)
dev.off()

# CLUSTER 4
gene_corscan_4<-scanpy_genes[scanpy_genes$Louvain_cluster==4,]
gene_corscan_4<-gene_corscan_4[,c(38:46)]

ix <- 1:9
gene_corscan_4[ix] <- lapply(gene_corscan_4[ix], as.character)
gene_corscan_4[ix] <- lapply(gene_corscan_4[ix], as.numeric)
# Remove genes with 0 counts
gene_corscan_4<- gene_corscan_4[,colSums(gene_corscan_4)!=0]


genetabcor_4<-correlate(gene_corscan_4)
fashion(genetabcor_4)
genetabcor_4 <- genetabcor_4[,colSums(is.na(genetabcor_4))<nrow(genetabcor_4)]

rplot(genetabcor_4)
network_plot(genetabcor_4,min_cor = 0.3)

clust4_genecor<-rcorr(as.matrix(gene_corscan_4))
flattenCorrMatrix(clust4_genecor$r, clust4_genecor$P)

pdf(file = "genes_corrplot_cluster4.pdf")
corrplot(clust4_genecor$r, order="hclust", 
         p.mat = clust4_genecor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 4",line = 0,adj=0)
dev.off()

# CLUSTER 5
gene_corscan_5<-scanpy_genes[scanpy_genes$Louvain_cluster==5,]
gene_corscan_5<-gene_corscan_5[,c(38:46)]

ix <- 1:9
gene_corscan_5[ix] <- lapply(gene_corscan_5[ix], as.character)
gene_corscan_5[ix] <- lapply(gene_corscan_5[ix], as.numeric)
# Remove genes with 0 counts
gene_corscan_5<- gene_corscan_5[,colSums(gene_corscan_5)!=0]


genetabcor_5<-correlate(gene_corscan_5)
fashion(genetabcor_5)
genetabcor_5 <- genetabcor_5[,colSums(is.na(genetabcor_5))<nrow(genetabcor_5)]

rplot(genetabcor_5)
network_plot(genetabcor_5,min_cor = 0.3)

clust5_genecor<-rcorr(as.matrix(gene_corscan_5))
flattenCorrMatrix(clust5_genecor$r, clust5_genecor$P)

pdf(file = "genes_corrplot_cluster5.pdf")
corrplot(clust5_genecor$r, order="hclust", 
         p.mat = clust5_genecor$P, sig.level = 0.05, insig = "blank",addrect = 3)
title("Cluster 5",line = 0,adj=0)
dev.off()

# CLUSTER 6
#gene_corscan_6<-scanpy_genes[scanpy_genes$Louvain_cluster==6,]
#gene_corscan_6<-gene_corscan_6[,c(36:44)]

#ix <- 1:9
#gene_corscan_6[ix] <- lapply(gene_corscan_6[ix], as.character)
#gene_corscan_6[ix] <- lapply(gene_corscan_6[ix], as.numeric)
# Remove genes with 0 counts
#gene_corscan_6<- gene_corscan_6[,colSums(gene_corscan_6)!=0]


#genetabcor_6<-correlate(gene_corscan_6)
#fashion(genetabcor_6)
#genetabcor_6 <- genetabcor_6[,colSums(is.na(genetabcor_6))<nrow(genetabcor_6)]

#rplot(genetabcor_6)
#network_plot(genetabcor_6,min_cor = 0.3)

#clust6_genecor<-rcorr(as.matrix(gene_corscan_6))
#flattenCorrMatrix(clust6_genecor$r, clust6_genecor$P)

#pdf(file = "genes_corrplot_cluster6.pdf")
#corrplot(clust6_genecor$r, order="hclust", 
#         p.mat = clust6_genecor$P, sig.level = 0.05, insig = "blank",addrect = 3)
#title("Cluster 6",line = 0,adj=0)
#dev.off()

### Heatmap and PCA based on surface markers ###

markall<-as.matrix(scanpy[,c(20:28)])
complete.cases(markall)
is.na(markall)
markall<-markall[complete.cases(markall),]

# PCA
#BiocManager::install("ggfortify")

library(ggfortify)
pcamark<-prcomp(markall,scale=TRUE)
summary(pcamark)

pcamark<-as.data.frame(pcamark$x)
row.names(pcamark)

scanpy_mark_genes<-merge(scanpy_genes,pcamark[,1:3],by="row.names",all=TRUE)
row.names(scanpy_mark_genes)<-scanpy_mark_genes$Row.names
scanpy_mark_genes$Row.names<-NULL

auto<-autoplot(prcomp(markall,scale=TRUE), data = markall,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)
auto+theme_bw()

# 2D PLOT

pcamark<-ggplot(scanpy_mark_genes, aes(x=PC1, y=PC2)) +
  geom_point(color="black",size=4) +
  geom_point(aes(color=as.factor(Louvain_cluster)),size=3,alpha=0.8) +
  scale_color_brewer(palette = "Spectral")+
  new_scale("color")+
  stat_ellipse(aes(color=as.factor(Louvain_cluster)),level=0.8)+
  scale_color_brewer(palette = "Spectral")+
  #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #facet_grid(~gender)+
  #scale_color_manual(values=c("dodgerblue2","coral","darkolivegreen","grey","black","purple"))+
  theme_bw()+
  theme(legend.position="bottom",
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1))
pcamark

cord<-ggplot(scanpy_mark_genes, aes(y=PC2, x=as.factor(Louvain_cluster))) +
  geom_point(aes(color=as.factor(Louvain_cluster)),size=3) +
  geom_violin(alpha=0)+
  xlab("")+
  scale_color_brewer(palette = "Spectral")+
  theme_bw()+
  theme(legend.position="none",
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1))
cord



pdf("/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/FACS_markers/PCA_markers_scale.pdf")
grid.arrange(pcamark,cord,ncol=2)
dev.off()

# 3D plot
p <- plot_ly(scanpy_mark_genes, x = ~PC1, y = ~PC2, z = ~PC3,color=~as.factor(Louvain_cluster)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3'))
  )

p


## Heatmap

#install.packages("heatmaply")
library(heatmaply)

#https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
# https://rdrr.io/cran/heatmaply/man/heatmaply.html

# markers matrix with no NAs
markall<-as.matrix(scanpy[,c(20:28)])
complete.cases(markall)
is.na(markall)
markall<-markall[complete.cases(markall),]

cold<-merge(scanpy[,c(4,8,11)],markall,by="row.names")
row.names(cold)<-cold$Row.names
cold$Row.names<-NULL
cold$Louvain_cluster<-as.factor(cold$Louvain_cluster)

# RColorBrewer::brewer.pal(7, "Spectral")
# "#D53E4F" "#FC8D59" "#FEE08B" "#FFFFBF" "#E6F598" "#99D594" "#3288BD"

row.names(cold)<-gsub('_.*','',row.names(cold))

#normalize and scale are different. normalize substracts the minimum, and divide by
# the maximum of all observations. This preserves the shape of each variable's distribution
# while making them easily comparable on the same "scale".

#Install phantomjs for saving file
webshot::install_phantomjs()

#colors for clusters
colclust = c("0" = "#D53E4F", 
             "1" = "#FC8D59",  
             "2" = "#FEE08B", 
             "3" = "#FFFFBF", 
             "4" = "#E6F598",
             "5" = "#99D594",
             "6" = "#3288BD")

# row_side_palette = colclust

heatmaply(cold,
          plot_method = "plotly",
          scale="column")


# gene expression of markers
ix <- 38:ncol(scanpy_mark_genes)
scanpy_mark_genes[ix] <- lapply(scanpy_mark_genes[ix], as.character)
scanpy_mark_genes[ix] <- lapply(scanpy_mark_genes[ix], as.numeric)

genemarkall<-as.matrix(scanpy_mark_genes[,c(38:46)])

complete.cases(genemarkall)
is.na(genemarkall)
genemarkall<-genemarkall[complete.cases(genemarkall),]

genemarkall<-merge(scanpy[,c(4,8,11)],genemarkall,by="row.names")
row.names(genemarkall)<-genemarkall$Row.names
genemarkall$Row.names<-NULL
genemarkall$Louvain_cluster<-as.factor(genemarkall$Louvain_cluster)

row.names(genemarkall)<-gsub('_.*','',row.names(genemarkall))

heatmaply(genemarkall,
               plot_method = "plotly",
          scale="column")



## Gene expression differences between clusters
pdf("/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/Notch.pdf")

n1<-ggplot(scanpy_mark_genes,aes(x=as.factor(Louvain_cluster),y=sc_Notch1))+
  geom_jitter()+
  #scale_color_gradient(low="white",high = "steelblue")+
  #new_scale("color")+
  geom_violin(aes(color=as.factor(Louvain_cluster)),scale="width",size=2)+
  scale_color_brewer(palette = "Spectral")+
  theme_bw()+
  theme(legend.position = "none")+
  geom_boxplot(width=0.1)

n2<-ggplot(scanpy_mark_genes,aes(x=as.factor(Louvain_cluster),y=sc_Notch2))+
  geom_jitter()+
  #scale_color_gradient(low="white",high = "steelblue")+
  #new_scale("color")+
  geom_violin(aes(color=as.factor(Louvain_cluster)),scale="width",size=2)+
  scale_color_brewer(palette = "Spectral")+
  theme_bw()+
  theme(legend.position = "none")+
  geom_boxplot(width=0.1)

grid.arrange(n1,n2,ncol=2)

dev.off()


#### SUBSELECTION FOR SCANPY ####

### Select counts excluding clusters 4 and 6 (or desired) for running scanpy and see
## the clusters reducing the variance from some clusters. Include cyclone data (before filtering QC)
clust_HSC_counts<-counts[,colnames(counts) %in% c(clust0,clust1,clust2,clust3,clust4,clust5)]
# Remove ^_ en el output para input scanpy (grep -v '^_')
write.table(clust_HSC_counts,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/scanpy_run/clustsub_counts.txt",sep="\t",quote = FALSE,row.names = TRUE)

cycsub<-cyctab
cycsub$cell<-row.names(cycsub)
cycsub<-cycsub[cycsub$cell %in% c(clust0,clust1,clust2,clust3,clust4,clust5),]
cycsub$cell=NULL

write.table(cycsub,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/scanpy_run/cyclone_R/sub_cyclone_scRNAseq_R.txt",sep="\t",quote = FALSE,row.names = TRUE)


cycG1sub<-cycG1
cycG1sub$cell<-row.names(cycG1sub)
cycG1sub<-cycG1sub[cycG1sub$cell %in% c(clust0,clust1,clust2,clust3,clust4,clust5),]
cycG1sub$cell=NULL

write.table(cycG1sub,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/scanpy_run/cyclone_R/sub_G1_score_cyclone.txt",sep="\t",quote = FALSE,row.names = TRUE)

cycG2Msub<-cycG2M
cycG2Msub$cell<-row.names(cycG2Msub)
cycG2Msub<-cycG2Msub[cycG2Msub$cell %in% c(clust0,clust1,clust2,clust3,clust4,clust5),]
cycG2Msub$cell=NULL

write.table(cycG2Msub,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/scanpy_run/cyclone_R/sub_G2M_score_cyclone.txt",sep="\t",quote = FALSE,row.names = TRUE)

cycSsub<-cycS
cycSsub$cell<-row.names(cycSsub)
cycSsub<-cycSsub[cycSsub$cell %in% c(clust0,clust1,clust2,clust3,clust4,clust5),]
cycSsub$cell=NULL

write.table(cycSsub,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/scanpy_run/cyclone_R/sub_S_score_cyclone.txt",sep="\t",quote = FALSE,row.names = TRUE)

#sub metadata
metasub<-meta[meta$sampleID %in% c(clust0,clust1,clust2,clust3,clust4,clust5),]

# For scanpy input, should remove header
write.table(metasub,"/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/scanpy_run/submetadata_scRNAseq_R.txt",sep="\t",quote = FALSE,row.names = TRUE)

# Identification of PRE HSC from Zhou et al
#####
#### IDENTIFICATION OF PRE HSC FROM Zhou et al

## Zhou genes

#zhou_degs<-read.table("/Volumes/grcmc/YGUILLEN/scRNA_cam/nature17997-s1/list_DEGs_Zhou.txt",header=TRUE,sep="\t")
zhou_degs<-read.table("/Users/yguillen/Desktop/temp/scRNA_cam/nature17997-s1/list_DEGs_Zhou.txt",header=TRUE,sep="\t")
table(zhou_degs$Compare)

# Patterns (sheet 4 and 5)
#patEC_Pre<-read.table("/Volumes/grcmc/YGUILLEN/scRNA_cam/nature17997-s1/list_EC_PreHSC_pattern.txt",header=TRUE,sep="\t")
patEC_Pre<-read.table("/Users/yguillen/Desktop/temp/scRNA_cam/nature17997-s1/list_EC_PreHSC_pattern.txt",header=TRUE,sep="\t")

#patPre_E14<-read.table("/Volumes/grcmc/YGUILLEN/scRNA_cam/nature17997-s1/list_PreHSC_E14.txt",header=TRUE,sep="\t")
patPre_E14<-read.table("/Users/yguillen/Desktop/temp/scRNA_cam/nature17997-s1/list_PreHSC_E14.txt",header=TRUE,sep="\t")

## EC to T1 pre HSC
EC_vs_T1pre<-subset(zhou_degs,zhou_degs$Compare=="EC.vs.T1 pre-HSC")
EC_vs_T1pre<-EC_vs_T1pre[order(EC_vs_T1pre$Up_or_Down,EC_vs_T1pre$FDR),]

listECT1<-subset(EC_vs_T1pre,select=c("Gene","Up_or_Down"))
listECT1

#subset first up and down 200 or 50 elements (50 for all together)
listECT1<-listECT1[c(1:200,246:446),]
#listECT1<-listECT1[c(1:50,246:296)]

## T1 to T2 pre HSC
T1pre_vs_T2pre<-subset(zhou_degs,zhou_degs$Compare=="T1 pre-HSC.vs.T2 pre-HSC")
T1pre_vs_T2pre<-T1pre_vs_T2pre[order(T1pre_vs_T2pre$FDR),]

listT1T2<-subset(T1pre_vs_T2pre,select=c("Gene","Up_or_Down"))
listT1T2

#subset first 100 elements
listT1T2<-listT1T2[c(1:100),]
dim(listT1T2)

## T2 to E12
T2pre_vs_E12<-subset(zhou_degs,zhou_degs$Compare=="T2 pre-HSC.vs.E12 HSC")
T2pre_vs_E12<-T2pre_vs_E12[order(T2pre_vs_E12$Up_or_Down,T2pre_vs_E12$FDR),]

listT2E12<-subset(T2pre_vs_E12,select=c("Gene","Up_or_Down"))
length(listT2E12)

#subset first up and down 200 elements
listT2E12<-listT2E12[c(1:200,351:nrow(T2pre_vs_E12)),]
#listT2E12<-listT2E12[c(1:50,351:401)]

## E12 to E14
E12_vs_E14<-subset(zhou_degs,zhou_degs$Compare=="E12 HSC.vs.E14 HSC")
E12_vs_E14<-E12_vs_E14[order(E12_vs_E14$Up_or_Down,E12_vs_E14$FDR),]

listE12E14<-subset(E12_vs_E14,select=c("Gene","Up_or_Down"))
listE12E14
# there are only 20

## list with all genes form different transitions
ECT1<-as.data.frame(listECT1$Gene)
ECT1$stage<-"EC_to_T1"
colnames(ECT1)[1]<-"gene"
row.names(ECT1)<-ECT1$gene

T1T2<-as.data.frame(listT1T2$Gene)
T1T2$stage<-"T1_to_T2"
colnames(T1T2)[1]<-"gene"
row.names(T1T2)<-T1T2$gene

T2E12<-as.data.frame(listT2E12$Gene)
T2E12$stage<-"T2_to_E12"
colnames(T2E12)[1]<-"gene"
row.names(T2E12)<-T2E12$gene

E12E14<-as.data.frame(listE12E14$Gene)
E12E14$stage<-"E12_to_E14"
colnames(E12E14)[1]<-"gene"
row.names(E12E14)<-E12E14$gene

allgenes<-rbind(ECT1,T1T2,T2E12,E12E14)
nrow(allgenes)
allgenes<-allgenes[!duplicated(allgenes$gene),]

# Select list of DEGs from gene expression matrix (depending on the subset, change %in%)
# for all: allgenes$gene (be carefull with the number of genes)
# for EC to T1: listECT1
# for T1 to T2: listT1T2
# for T2 to E12: listT2E12
# for E12 to E14: listE12E14
# for pattern EC to Pre (T1 and T2)

#all cells
t_degs<-t(exprmat[exprmat$gene %in% as.character(listT1T2$Gene),])

colnames(t_degs)<-t_degs[1,]
t_degs<-t_degs[-1,]
t_degs<-as.data.frame(t_degs)

#for only cluster 4
#t_degs<-t_degs[row.names(t_degs) %in% clust4,]


ix <- 1:ncol(t_degs)
t_degs[ix] <- lapply(t_degs[ix], as.character)
t_degs[ix] <- lapply(t_degs[ix], as.numeric)

t_degs<-as.matrix(t_degs)

complete.cases(t_degs)
is.na(t_degs)
t_degs<-t_degs[complete.cases(t_degs),]

cold<-merge(scanpy[,c(4,11)],t_degs,by="row.names")
row.names(cold)<-cold$Row.names
cold$Row.names<-NULL
cold$Louvain_cluster<-as.factor(cold$Louvain_cluster)

# RColorBrewer::brewer.pal(7, "Spectral")
# "#D53E4F" "#FC8D59" "#FEE08B" "#FFFFBF" "#E6F598" "#99D594" "#3288BD"

row.names(cold)<-gsub('_.*','',row.names(cold))


colnames(t_degs)
name_degs<-as.data.frame(colnames(t_degs))
colnames(name_degs)<-"gene"
row.names(name_degs)<-name_degs$gene

# only for all genes together
#name_degs<-merge(name_degs,allgenes,by="gene", alL.x=TRUE,sort = FALSE)
#name_degs<- name_degs[order(name_degs$stage),]

#for each transition change dataframe (listECT1,listT1T2,listT2E12,listE12E14)
name_degs<-merge(name_degs,listT1T2,by.x="gene",by.y="Gene", alL.x=TRUE,sort = FALSE)
#name_degs<- name_degs[order(name_degs$Up_or_Down),]
nrow(name_degs)

heatmaply(cold,
          plot_method = "plotly",
          ColSideColors = name_degs$Up_or_Down,
          scale="column",
          k_row = 3)


browseURL("/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/heatmap_T1_T2_Zhou.html")





#### DEGs wilcoxon test, scanpy ####
#DEGs_clusters_pval<-read.table("/Volumes/grcmc/YGUILLEN/scRNA_cam/scanpy_run/output_DEGs_adatared_wilcox_pval_adj.csv",header=TRUE,sep=",")
DEGs_clusters_pval<-read.table("/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/output_DEGs_adatared_wilcox_pval_adj.csv",header=TRUE,sep=",")

# pvalue < 0.01
#cluster 6 DEGs
#DEGsclust6<-(DEGs_clusters_pval_adj[DEGs_clusters_pval$X6_p<0.05,])$X6_n
#DEGsclust6<-as.character(DEGsclust6)

#cluster 5 DEGs
DEGsclust5<-(DEGs_clusters_pval[DEGs_clusters_pval$X5_p<0.01,])$X5_n
DEGsclust5<-as.character(DEGsclust5)

#cluster 4 DEGs
DEGsclust4<-(DEGs_clusters_pval[DEGs_clusters_pval$X4_p<0.01,])$X4_n
DEGsclust4<-as.character(DEGsclust4)

#cluster 3 DEGs
DEGsclust3<-(DEGs_clusters_pval[DEGs_clusters_pval$X3_p<0.01,])$X3_n
DEGsclust3<-as.character(DEGsclust3)

#cluster 2 DEGs
DEGsclust2<-(DEGs_clusters_pval[DEGs_clusters_pval$X2_p<0.01,])$X2_n
DEGsclust2<-as.character(DEGsclust2)

#cluster 1 DEGs
DEGsclust1<-(DEGs_clusters_pval[DEGs_clusters_pval$X1_p<0.01,])$X1_n
DEGsclust1<-as.character(DEGsclust1)

#cluster 0 DEGs
DEGsclust0<-(DEGs_clusters_pval[DEGs_clusters_pval$X0_p<0.01,])$X0_n
DEGsclust0<-as.character(DEGsclust0)

# Select list of DEGs from gene expression matrix
t_degs<-t(exprmat[exprmat$gene %in% DEGsclust4,])
colnames(t_degs)<-t_degs[1,]
t_degs<-t_degs[-1,]
t_degs<-as.data.frame(t_degs)


ix <- 1:ncol(t_degs)
t_degs[ix] <- lapply(t_degs[ix], as.character)
t_degs[ix] <- lapply(t_degs[ix], as.numeric)

t_degs<-as.matrix(t_degs)

complete.cases(t_degs)
is.na(t_degs)
t_degs<-t_degs[complete.cases(t_degs),]

cold<-merge(t_degs,scanpy,by="row.names",sort=FALSE)

cold<-subset(cold,select=c("sampleID","Louvain_cluster","population"))
row.names(cold)<-cold$sampleID


row.names(t_degs)<-gsub('_.*','',row.names(t_degs))

tmp<-heatmaply(t_degs,
               plot_method = "plotly",
               row_side_palette = colclust,
               RowSideColors = cold$Louvain_cluster)
#               file="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/heatmap_cluster3_DEGs_scanpy.html")
tmp

browseURL("/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/heatmap_cluster3_DEGs_scanpy.html")







library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "GO_Cellular_Component_2018",
         "GO_Molecular_Function_2018",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2019_Mouse",
         "KEGG_2019_Mouse",
         "Mouse_Gene_Atlas")


clustenrich <- enrichr(DEGsclust5, dbs)

clustenrich <- clustenrich[["GO_Biological_Process_2018"]]


bpsub<-subset(clustenrich,clustenrich$Adjusted.P.value<=0.1)
bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))
bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$Combined.Score)), ]$Term)


ggplot(bpsub,aes(y=Combined.Score,x=Term),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",color="black")+
  #geom_text(aes(label=tolower(Genes)), size=4,hjust=0.5,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Blues")+
  #facet_wrap(~group,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=6, hjust = 1),
        axis.text.y = element_text(size=9, hjust = 1))


# Histogram of markers distribution


ggplot(scanpy_mark_genes,aes(x=as.factor(Louvain_cluster),y=CD31))+
  geom_violin(aes(color=as.factor(Louvain_cluster)))+
  geom_boxplot(width=0.1,alpha=0.5)+
  scale_color_brewer(palette="Spectral")+
  #geom_vline(data=filter(metaUMAPmelt, variable=="CD31"), aes(xintercept=1000),linetype='dashed',color="red" ) + 
  theme_minimal()+
  theme(legend.position = "bottom",axis.text.x = element_blank(),axis.ticks.x = element_blank())

ggplot(scanpy_mark_genes,aes(x=CD31))+
  geom_density(aes(color=as.factor(Louvain_cluster)))+
  scale_color_brewer(palette="Spectral")+
  geom_vline(aes(xintercept=1000),linetype='dashed',color="red" ) + 
  geom_vline(aes(xintercept=median((scanpy_mark_genes[!is.na(scanpy_mark_genes$CD31) & scanpy_mark_genes$Louvain_cluster=="1",])$CD31)),linetype='dashed',color="blue" ) + 
  theme_minimal()+
  theme(legend.position = "bottom",axis.text.x = element_blank(),axis.ticks.x = element_blank())


genelist2<-c("Cdh5","Pecam1")
t_geneexp<-t(exprmat[exprmat$gene %in% genelist2,])
colnames(t_geneexp)<-t_geneexp[1,]
t_geneexp<-t_geneexp[-1,]

scanpy_genes<-merge(scanpy,t_geneexp,by="row.names")
row.names(scanpy_genes)<-scanpy_genes$Row.names
scanpy_genes$Row.names=NULL


scanpy_genes$Cdh5<-as.numeric(as.character(scanpy_genes$Cdh5))
scanpy_genes$Pecam1<-as.numeric(as.character(scanpy_genes$Pecam1))

endpl1<-ggplot(scanpy_genes,aes(x=Pecam1))+
  geom_density(aes(color=as.factor(Louvain_cluster)),size=1.5)+
  geom_rug(data=(scanpy_genes[scanpy_genes$Louvain_cluster=="1" | scanpy_genes$Louvain_cluster=="4",]),color="black",size=1)+
  geom_rug(data=(scanpy_genes[scanpy_genes$Louvain_cluster=="1" | scanpy_genes$Louvain_cluster=="4",]),aes(color=as.factor(Louvain_cluster)))+
  scale_color_brewer(palette="Spectral")+
  geom_vline(aes(xintercept=median((scanpy_genes[!is.na(scanpy_genes$Pecam1) & scanpy_genes$Louvain_cluster=="4",])$Cdh5)),linetype='dashed' ) + 
  theme_minimal()+
  theme(legend.position = "bottom",axis.text.x = element_blank(),axis.ticks.x = element_blank())

ggsave(endpl1,file="/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/Plots/scanpy/Pecam1_hist.pdf")

endpl2<-ggplot(scanpy_genes,aes(x=as.factor(Louvain_cluster),y=Pecam1))+
  geom_violin(aes(color=as.factor(Louvain_cluster)),width=2,size=1.5)+
  geom_boxplot(width=0.1,alpha=0.5)+
  scale_color_brewer(palette="Spectral")+
  theme_minimal()+
  theme(legend.position = "bottom",axis.ticks.x = element_blank())
ggsave(endpl2,file="/Users/yolanda_guillen/Desktop/IMIM/scRNA_cam/Roshana_scRNA/Plots/scanpy/Pecam1_box.pdf")


endpl1
endpl2


## Select T1 - T2 preHSC from cluster 4
# Lab notebook, page 16 back. Based on DEGs T1 - T2 Zhou et al

T2sub<-c("C206","C292","C174","C235","C150","C308","C102","C344",
         "C343","C317","C042","C348","C242","C120")

T1sub<-c("C328","C011","C056","C093","C087","C342","C147",
         "C192","C148","C221","C022","C268","C110","C180")

TpreNA<-c("C039","C148","C110","C216")

scanpy$sampleID<-as.character(scanpy$sampleID)


T1preclust<-scanpy_mark_genes[grepl(paste(T1sub, collapse="|"),scanpy_mark_genes$sampleID),]
T1preclust$preHSC<-"T1"
T2preclust<-scanpy_mark_genes[grepl(paste(T2sub, collapse="|"),scanpy_mark_genes$sampleID),]
T2preclust$preHSC<-"T2"
TNApreclust<-scanpy_mark_genes[grepl(paste(TpreNA, collapse="|"),scanpy_mark_genes$sampleID),]
TNApreclust$preHSC<-"TpreNA"
nopreclust<-scanpy_mark_genes[!grepl(paste(c(T1sub,T2sub,TpreNA), collapse="|"),scanpy_mark_genes$sampleID),]
nopreclust$preHSC<-"No_preHSC"


preHSCclas<-rbind(T1preclust,T2preclust,TNApreclust,nopreclust)

# Based on CD45 levels
summary((scanpy[scanpy$Louvain_cluster=="4",])$CD45)
# setting min CD45 as 6000

T1preclust<-scanpy_mark_genes[scanpy_mark_genes$CD45<6000 & scanpy_mark_genes$Louvain_cluster=="4",]
T1preclust$preHSC<-"T1"
T1preclust<-T1preclust[row.names(T1preclust)!="NA",]
T2preclust<-scanpy_mark_genes[scanpy_mark_genes$CD45>6000 & scanpy_mark_genes$Louvain_cluster=="4",]
T2preclust$preHSC<-"T2"
T2preclust<-T2preclust[row.names(T2preclust)!="NA",]
nopreclust<-scanpy_mark_genes[scanpy_mark_genes$Louvain_cluster!="4",]
nopreclust$preHSC<-"No_preHSC"

preHSCclas<-rbind(T1preclust,T2preclust,nopreclust)

#JAG1
Jag1min=min((scanpy[!is.na(scanpy$Jag1) ,])$Jag1)
Jag1max=max((scanpy[!is.na(scanpy$Jag1) ,])$Jag1)
Jag1med=median((scanpy[!is.na(scanpy$Jag1),])$Jag1)

#NOTCH1
Notch1min=min((scanpy[!is.na(scanpy$Notch1),])$Notch1)
Notch1max=max((scanpy[!is.na(scanpy$Notch1),])$Notch1)
Notch1med=median((scanpy[!is.na(scanpy$Notch1),])$Notch1)

#NOTCH2
Notch2min=min((scanpy[!is.na(scanpy$Notch2),])$Notch2)
Notch2max=max((scanpy[!is.na(scanpy$Notch2),])$Notch2)
Notch2med=median((scanpy[!is.na(scanpy$Notch2),])$Notch2)

#DLL4
Dll4min=min((scanpy[!is.na(scanpy$Dll4) ,])$Dll4)
Dll4max=max((scanpy[!is.na(scanpy$Dll4) ,])$Dll4)
Dll4med=median((scanpy[!is.na(scanpy$Dll4),])$Dll4)

# protein FACS
ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="black",size=7)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=7)+
  #geom_point(data = preHSCclas[preHSCclas$preHSC=="TpreNA",],color="blue",size=7)+
  geom_point(color="black",size=6)+
  geom_point(aes(color=as.numeric(Notch2)),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="red2", 
                        midpoint = Notch2med,
                        breaks=seq(Notch2min,Notch2max,(Notch2max-Notch2min)/4))+
  #stat_ellipse(data=preHSCclas[preHSCclas$Louvain_cluster=="4",],aes(color=preHSC),level=0.95,size=1)+
  theme_bw()+
  theme(legend.position = "right")

# Gene
ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="black",size=7)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=7)+
  #geom_point(data = preHSCclas[preHSCclas$preHSC=="TpreNA",],color="blue",size=7)+
  geom_point(color="black",size=6)+
  geom_point(aes(color=sc_Dll4),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="red2", 
                        midpoint = median(preHSCclas[!is.na(preHSCclas$sc_Dll4),]$sc_Dll4))+
  #stat_ellipse(data=preHSCclas[preHSCclas$Louvain_cluster=="4",],aes(color=preHSC),level=0.95,size=1)+
  theme_bw()+
  theme(legend.position = "right")


preT2mark<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="black",size=7)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=7)+
  #geom_point(data = preHSCclas[preHSCclas$preHSC=="TpreNA",],color="blue",size=7)+
  geom_point(color="black",size=5)+
  geom_point(aes(color=as.numeric(CD45)),size=4)+
  scale_color_gradient(low="white",high = "coral")+
  #stat_ellipse(data=preHSCclas[preHSCclas$Louvain_cluster=="4",],aes(color=preHSC),level=0.95,size=1)+
  theme_bw()+
  theme(legend.position = "none")
preT2mark

preT2gene<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="black",size=7)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=7)+
  #geom_point(data = preHSCclas[preHSCclas$preHSC=="TpreNA",],color="blue",size=7)+
  geom_point(color="black",size=5)+
  geom_point(aes(color=sc_Nnat),size=4)+
  scale_color_gradient(low="white",high = "coral")+
  #stat_ellipse(data=preHSCclas[preHSCclas$Louvain_cluster=="4",],aes(color=preHSC),level=0.95,size=1)+
  theme_bw()+
  theme(legend.position = "none")
preT2gene

preT2<-grid.arrange(preT2mark,preT2gene)

ggsave(preT2,file="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/preT1T2_Zhou_UMAP.pdf")


## Notch ligands and receptors accross stages

preHSCclas$group<-paste(preHSCclas$Louvain_cluster,preHSCclas$preHSC,sep="_")
table(preHSCclas$group)
preHSCclas$group<-factor(preHSCclas$group,levels = c("3_No_preHSC","5_No_preHSC","0_No_preHSC","1_No_preHSC","2_No_preHSC","4_T1","4_T2"))


# Only markers in cluster 2, and 4 T1 and 4 T2

preHSCclas_sub<-preHSCclas[preHSCclas$group=="2_No_preHSC" | preHSCclas$group=="4_T1" | preHSCclas$group=="4_T2",]


gNOTCH1<-ggplot(preHSCclas_sub,aes(x=group,y=sc_Notch1))+
  #geom_point(color="black",size=3)+
  #geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(aes(fill=as.factor(Louvain_cluster)),alpha=0.5)+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=5,shape=23)  + 
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=median, geom="point",size=5,shape=23)  + 
  scale_fill_brewer(palette="Spectral")+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("Notch1")


pNOTCH1<-ggplot(preHSCclas_sub,aes(x=group,y=Notch1))+
  #geom_point(color="black",size=3)+
  #geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(aes(fill=as.factor(Louvain_cluster)),alpha=0.5)+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=5,shape=23)  + 
  geom_hline(yintercept = 500)+
  #scale_fill_brewer(palette="Spectral")+
  scale_fill_manual(values=c("orange","green","green"))+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("NOTCH1")

  
notch1clust<-ggplot(preHSCclas,aes(x=as.factor(Louvain_cluster),y=Notch1))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(alpha=0)+
  scale_color_brewer(palette="Spectral")+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=8,shape=23)  + 
  scale_fill_brewer(palette="Spectral")+
  geom_hline(yintercept = 500,linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("NOTCH1")

notch1clust

gNOTCH2<-ggplot(preHSCclas,aes(x=group,y=sc_Notch2))+
  #geom_point(color="black",size=3)+
  #geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(aes(fill=as.factor(Louvain_cluster)),alpha=0.5)+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=5,shape=23)  + 
  scale_fill_brewer(palette="Spectral")+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("Notch2")


pNOTCH2<-ggplot(preHSCclas_sub,aes(x=group,y=Notch2))+
  #geom_point(color="black",size=3)+
  #geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(aes(fill=as.factor(Louvain_cluster)),alpha=0.5)+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=5,shape=23)  + 
  geom_hline(yintercept = 200)+
  #scale_fill_brewer(palette="Spectral")+
  scale_fill_manual(values=c("orange","green","green"))+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("NOTCH2")


notch2clust<-ggplot(preHSCclas,aes(x=as.factor(Louvain_cluster),y=Notch2))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(alpha=0)+
  scale_color_brewer(palette="Spectral")+
  stat_summary(fun.y=median, geom="point",size=5,shape=23)  + 
  geom_hline(yintercept = 200,linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("NOTCH2")

notch2clust

gDLL4<-ggplot(preHSCclas,aes(x=group,y=sc_Dll4))+
  #geom_point(color="black",size=3)+
  geom_violin(aes(fill=as.factor(Louvain_cluster)),alpha=0.5)+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=5,shape=23)  + 
  scale_fill_brewer(palette="Spectral")+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("Dll4")

pDLL4<-ggplot(preHSCclas_sub,aes(x=group,y=Dll4))+
  #geom_point(color="black",size=3)+
  #geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(aes(fill=as.factor(Louvain_cluster)),alpha=0.5)+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=5,shape=23)  + 
  geom_hline(yintercept = 120)+
  #scale_fill_brewer(palette="Spectral")+
  scale_fill_manual(values=c("orange","green","green"))+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("DLL4")


dll4clust<-ggplot(preHSCclas,aes(x=as.factor(Louvain_cluster),y=Dll4))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(alpha=0)+
  scale_color_brewer(palette="Spectral")+
  stat_summary(fun.y=median, geom="point",size=5,shape=23)  + 
  geom_hline(yintercept = 120,linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("DLL4")

dll4clust


gJAG1<-ggplot(preHSCclas,aes(x=group,y=sc_Jag1))+
  #geom_point(color="black",size=3)+
  #geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(aes(fill=as.factor(Louvain_cluster)),alpha=0.5)+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=5,shape=23)  + 
  scale_fill_brewer(palette="Spectral")+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("Jag1")

pJAG1<-ggplot(preHSCclas_sub,aes(x=group,y=Jag1))+
  #geom_point(color="black",size=3)+
  #geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(aes(fill=as.factor(Louvain_cluster)),alpha=0.5)+
  stat_summary(aes(fill=as.factor(Louvain_cluster)),fun.y=median, geom="point",size=5,shape=23)  + 
  geom_hline(yintercept = 11000)+
  #scale_fill_brewer(palette="Spectral")+
  scale_fill_manual(values=c("orange","green","green"))+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("JAG1")


jag1clust<-ggplot(preHSCclas,aes(x=as.factor(Louvain_cluster),y=Jag1))+
  geom_point(color="black",size=3)+
  geom_point(aes(color=as.factor(Louvain_cluster)))+
  geom_violin(alpha=0)+
  scale_color_brewer(palette="Spectral")+
  stat_summary(fun.y=median, geom="point",size=5,shape=23)  + 
  geom_hline(yintercept = 11000,linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("JAG1")

grid.arrange(notch1clust,notch2clust,dll4clust,jag1clust,ncol=2)
grid.arrange(pNOTCH1,pNOTCH2,pDLL4,pJAG1,ncol=2)
grid.arrange(gNOTCH1,gNOTCH2,gDLL4,gJAG1,ncol=2)



### Gradient colors in UMAPs

Notch2_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=as.numeric(Notch2)),size=4)+
  scale_color_gradient2(low="white",high = "red",midpoint = 200)+
  theme_bw()+
  theme(legend.position = "bottom")
Notch2_grad

gNotch2_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=sc_Notch2),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")
gNotch2_grad


mid=500

Notch1_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=as.numeric(Notch1)),size=4)+
  scale_color_gradient2(low = "white", mid="white",high = "red", midpoint = mid)+
  theme_bw()+
  theme(legend.position = "bottom")

Notch1_grad

gNotch1_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=sc_Notch1),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")
gNotch1_grad



Gfi1_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=as.numeric(Gfi1)),size=4)+
  scale_color_gradient2(low="white",high = "red",midpoint = 550)+
  theme_bw()+
  theme(legend.position = "bottom")

Gfi1_grad


gGfi1<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=as.numeric(sc_Gfi1)),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

gGfi1


Dll4_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=as.numeric(Dll4)),size=4)+
  scale_color_gradient2(low="white",high = "red",midpoint = 120)+
  theme_bw()+
  theme(legend.position = "bottom")

Dll4_grad

gDll4_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=sc_Dll4),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")
gDll4_grad


Jag1_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=as.numeric(Jag1)),size=4)+
  scale_color_gradient2(low="white",high = "red",midpoint = 11000)+
  theme_bw()+
  theme(legend.position = "bottom")

Jag1_grad

gJag1_grad<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="black",size=5)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="blue",size=6)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=6)+
  geom_point(aes(color=sc_Jag1),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")
gJag1_grad

grid.arrange(Notch1_grad,gNotch1_grad,Notch2_grad,gNotch2_grad,
             Dll4_grad,gDll4_grad,Jag1_grad,gJag1_grad,ncol=4)




### FACS ligands - receptors ###

preHSCclas$receptor<-paste(preHSCclas$th_NOTCH1,preHSCclas$th_NOTCH2,sep="_")
table(preHSCclas$receptor,preHSCclas$Louvain_cluster)
preHSCclas$ligand<-paste(preHSCclas$th_DLL4,preHSCclas$th_JAG1,sep="_")
table(preHSCclas$ligand,preHSCclas$Louvain_cluster)

preHSCclas$comb<-paste(preHSCclas$receptor,preHSCclas$ligand,sep="_")
table(preHSCclas$comb,preHSCclas$Louvain_cluster)


reclig<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  geom_point(color="black",size=6)+
  geom_point(aes(color=comb),size=5)+
  scale_color_brewer(palette="Set1")+
  #geom_label_repel(data=subset(scanpy,(scanpy$UMAP_C1<(-10) & scanpy$UMAP_C2>7.5)),aes(label=sampleID),size=3,alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

preT2mark

Gfimark<-ggplot(preHSCclas,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T2",],color="black",size=7)+
  geom_point(data = preHSCclas[preHSCclas$preHSC=="T1",],color="green",size=7)+
  #geom_point(data = preHSCclas[preHSCclas$preHSC=="TpreNA",],color="blue",size=7)+
  geom_point(color="black",size=5)+
  geom_point(aes(color=as.numeric(Gfi1)),size=4)+
  scale_color_gradient(low="white",high = "blue")+
  #stat_ellipse(data=preHSCclas[preHSCclas$Louvain_cluster=="4",],aes(color=preHSC),level=0.95,size=1)+
  theme_bw()+
  theme(legend.position = "none")
Gfimark



grid.arrange(preT2mark,Gfimark,reclig,ncol=2)


#### PLOT GENE EXPRESSION AND METADATA

metext<-exprmat
row.names(metext)<-gsub('ENS.*\\|','',row.names(metext))
metext$gene=NULL

metext<-as.data.frame(t(metext))
dim(metext)


metext<-merge(scanpy,metext,by="row.names")
row.names(metext)<-metext$Row.names
metext$Row.names=NULL

#PLOT Hey2
Hey1<-ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Hey1),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

Dll4_p<-ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Dll4.y),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

Hey2<-ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Hey2),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

Jag1_p<-ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Jag1.y),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")


notch_targets<-grid.arrange(Hey1,Dll4_p,Hey2,Jag1_p,ncol=2)

ggsave(notch_targets,file="/Volumes/grcmc/YGUILLEN/scRNA_cam/Roshana_scRNA/Plots/scanpy/notch_targets_exp.pdf")

#Retinoic acid receptors
Rara<-ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Rara),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

Rarb<-ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Rarb),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

Rarg<-ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Rarg),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

Rxrg<-ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Rxrg),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

grid.arrange(Rara,Rarb,Rarg,Rxrg,ncol=2)

ggplot(metext,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Irx2),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")

# nfkb members
metext_nfkb<-subset(metext,select=c("UMAP_C1","UMAP_C2","Nfkbia","Nfkb1","Nfkb2","Rela","Relb","Rel","Nfkbib","Nfkbie","Nfkbiz","Nfkbid"))

ggplot(metext_nfkb,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=Nfkbid),size=4)+
  scale_color_gradient2(low="white",high = "red")+
  theme_bw()+
  theme(legend.position = "bottom")


################## COEXPRESSED GENES MATRIX #################

comat<-exprmat
row.names(comat)<-gsub('ENS.*\\|','',row.names(comat))
comat$gene=NULL

comat_t<-as.data.frame(t(comat))
dim(comat_t)

# With Jag1 PROTEIN or other genes
# Matrix with cells in HSC cluster 4
hsc_comat_t<-comat_t[row.names(comat_t) %in% row.names(scanpy[scanpy$Louvain_cluster=="4",]),]
scanpy_clust4<-subset(scanpy,scanpy$Louvain_cluster==4)

#corhey2<-apply(comat_t, 2, cor.test,scanpy$Notch1,method="spearman")
corhey2<-apply(hsc_comat_t, 2, cor.test,as.numeric(scanpy_clust4$Jag1),method="spearman")
corhey2_p<-as.data.frame(sapply(corhey2, "[[", "p.value"))
corhey2_e<-as.data.frame(sapply(corhey2, "[[", "estimate"))
corhey2<-cbind(corhey2_p,corhey2_e)
names(corhey2)<-c("Pval","Rho")
corhey2$Gene<-row.names(corhey2)

corhey2_sig<-subset(corhey2,corhey2$Pval<0.005)

corhey2_sig<-corhey2_sig[with(corhey2_sig,order(corhey2_sig$Pval, corhey2_sig$Rho)),]


### Merge gene expression with metadata
corhey2_list<-corhey2_sig$Gene
t_corexp_hey2<-t(exprmat[exprmat$gene %in% st,])
colnames(t_corexp_hey2)<-t_corexp_hey2[1,]
t_corexp_hey2<-t_corexp_hey2[-1,]

scanpy_cor<-merge(t_corexp_hey2,scanpy,by="row.names",sort=FALSE)
row.names(scanpy_cor)<-scanpy_cor$Row.names
scanpy_cor$Row.names=NULL

ncol(t_corexp_hey2)
ncol(scanpy)

# Top 10 correlated
top10_hey2<-(corhey2_sig[c(2:20),])$Gene



scanpy_melt1<-scanpy_cor[,colnames(scanpy_cor) %in% top10_hey2]
scanpy_melt2<-scanpy_cor[,c((ncol(t_corexp_hey2)+2):(ncol(t_corexp_hey2)+4))]
scanpy_melt<-cbind(scanpy_melt1,scanpy_melt2)


scanpy_melt<-melt(scanpy_melt,id.vars = c("UMAP_C1","UMAP_C2","Louvain_cluster"))

#PLOT 
ggplot(scanpy_melt,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(color="grey",size=5)+
  geom_point(aes(color=(as.numeric(as.character(value)))),size=4)+
  scale_color_gradient2(low="grey",high = "red")+
  facet_wrap(~variable,scales = "free",ncol=4)+
  theme_bw()+
  theme(legend.position = "bottom")





## Enrichment
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "GO_Cellular_Component_2018",
         "GO_Molecular_Function_2018",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2019_Mouse",
         "KEGG_2019_Mouse",
         "Mouse_Gene_Atlas")


clustenrich <- enrichr(corhey2_sig[corhey2_sig$Rho<0,]$Gene, dbs)

clustenrich <- clustenrich[["GO_Biological_Process_2018"]]


bpsub<-subset(clustenrich,clustenrich$P.value<=0.01)
bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))
bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$Combined.Score)), ]$Term)


ggplot(bpsub,aes(y=Combined.Score,x=Term),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",color="black")+
  #geom_text(aes(label=tolower(Genes)), size=4,hjust=0.5,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Blues")+
  #facet_wrap(~group,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=6, hjust = 1),
        axis.text.y = element_text(size=9, hjust = 1))







### PREPARE for GSEA
# expression removing cluster 6

normexpr<-read.delim("/Volumes/grcmc/YGUILLEN/scRNA_cam/scanpy_run/cellbrowser_noclust6_output/exprMatrix.tsv")
row.names(normexpr)<-normexpr$gene
normexpr$gene<-NULL
row.names(normexpr)<-gsub('ENS.*\\|','',row.names(normexpr))
colnames(normexpr)

# HE cluster
clust2<-row.names(scanpy[scanpy$Louvain_cluster=="2",])
clust2mat<-normexpr[,colnames(normexpr) %in% clust2]

# HSC cluster
clust4<-row.names(scanpy[scanpy$Louvain_cluster=="4",])
clust4mat<-normexpr[,colnames(normexpr) %in% clust4]

clustGSEA<-cbind(clust2mat,clust4mat)
dim(clustGSEA)

row.names(clustGSEA)<-toupper(row.names(clustGSEA))

write.table(clustGSEA,"/Volumes/grcmc/YGUILLEN/scRNA_cam/scanpy_run/cellbrowser_noclust6_output/GSEA/clustGSEA_normcounts.txt",quote = FALSE,sep="\t")

# num samples for labels file
dim(clust2mat)[2]
dim(clust4mat)[2]

name2<-"clust2"
times <- dim(clust2mat)[2]
names2<-rep(name2, times = times, length.out = NA, each = 1)

name4<-"clust4"
times <- dim(clust4mat)[2]
names4<-rep(name4, times = times, length.out = NA, each = 1)

namesfin<-c(names2,names4)
length(namesfin)
namesfin<-as.data.frame(t(namesfin))
dim(namesfin)

# output and then remove line with sample ids, and add class groups first line and second line
write.table(namesfin,"/Volumes/grcmc/YGUILLEN/scRNA_cam/scanpy_run/cellbrowser_noclust6_output/GSEA/labels.cls",row.names = FALSE,col.names = FALSE,quote = FALSE)



### April 2020 from Roshana's characterization of HE cells3 

# Hypothesis --> Jagged1 block Notch1 in HE to allow IAHC formation
# In HE (GFI1+ cKIT-), Jag1 and Notch are increased

# Do plots by cluster

Gfi1min=min((scanpy[!is.na(scanpy$Gfi1),])$Gfi1)
Gfi1max=max((scanpy[!is.na(scanpy$Gfi1),])$Gfi1)
Gfi1med=median((scanpy[!is.na(scanpy$Gfi1) & scanpy$Gfi1>270,])$Gfi1)


Notch1min=min((scanpy[!is.na(scanpy$Notch1),])$Notch1)
Notch1max=max((scanpy[!is.na(scanpy$Notch1),])$Notch1)
Notch1med=median((scanpy[!is.na(scanpy$Notch1) & scanpy$Notch1>500,])$Notch1)

Clust5_Notch<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="5",],aes(color=Notch1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                       midpoint = Notch1med,
                       breaks=seq(Notch1min,Notch1max,(Notch1max-Notch1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 5",
       subtitle = "NOTCH1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust2_Notch<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="2",],aes(color=Notch1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Notch1med,
                        breaks=seq(Notch1min,Notch1max,(Notch1max-Notch1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 2",
       subtitle = "NOTCH1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust3_Notch<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="3",],aes(color=Notch1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Notch1med,
                        breaks=seq(Notch1min,Notch1max,(Notch1max-Notch1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 3",
       subtitle = "NOTCH1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust1_Notch<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="1",],aes(color=Notch1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Notch1med,
                        breaks=seq(Notch1min,Notch1max,(Notch1max-Notch1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 1",
       subtitle = "NOTCH1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust0_Notch<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="0",],aes(color=Notch1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Notch1med,
                        breaks=seq(Notch1min,Notch1max,(Notch1max-Notch1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 0",
       subtitle = "NOTCH1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Clust4_Notch<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="4",],aes(color=Notch1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Notch1med,
                        breaks=seq(Notch1min,Notch1max,(Notch1max-Notch1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 4",
       subtitle = "NOTCH1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


grid.arrange(Clust5_Notch,Clust2_Notch,Clust3_Notch,
             Clust1_Notch,Clust0_Notch,Clust4_Notch,ncol=6)


CD45min=min((scanpy[!is.na(scanpy$CD45),])$CD45)
CD45max=max((scanpy[!is.na(scanpy$CD45),])$CD45)
CD45med=median((scanpy[!is.na(scanpy$CD45),])$CD45)

Clust5_CD45<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="5",],aes(color=CD45),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = CD45med,
                        breaks=seq(CD45min,CD45max,(CD45max-CD45min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 5",
       subtitle = "CD45 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust2_CD45<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="2",],aes(color=CD45),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = CD45med,
                        breaks=seq(CD45min,CD45max,(CD45max-CD45min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 2",
       subtitle = "CD45 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust3_CD45<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="3",],aes(color=CD45),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = CD45med,
                        breaks=seq(CD45min,CD45max,(CD45max-CD45min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 3",
       subtitle = "CD45 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust1_CD45<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="1",],aes(color=CD45),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = CD45med,
                        breaks=seq(CD45min,CD45max,(CD45max-CD45min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 1",
       subtitle = "CD45 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust0_CD45<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="0",],aes(color=CD45),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = CD45med,
                        breaks=seq(CD45min,CD45max,(CD45max-CD45min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 0",
       subtitle = "CD45 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Clust4_CD45<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="4",],aes(color=CD45),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = CD45med,
                        breaks=seq(CD45min,CD45max,(CD45max-CD45min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 4",
       subtitle = "CD45 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


grid.arrange(Clust5_CD45,Clust2_CD45,Clust3_CD45,
             Clust1_CD45,Clust0_CD45,Clust4_CD45,ncol=6)

Jag1min=min((scanpy[!is.na(scanpy$Jag1) ,])$Jag1)
Jag1max=max((scanpy[!is.na(scanpy$Jag1) ,])$Jag1)
Jag1med=median((scanpy[!is.na(scanpy$Jag1),])$Jag1)

Clust5_Jag1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="5",],aes(color=Jag1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Jag1med,
                        breaks=seq(Jag1min,Jag1max,(Jag1max-Jag1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 5",
       subtitle = "JAG1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust2_Jag1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="2",],aes(color=Jag1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Jag1med,
                        breaks=seq(Jag1min,Jag1max,(Jag1max-Jag1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 2",
       subtitle = "JAG1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Clust3_Jag1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="3",],aes(color=Jag1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Jag1med,
                        breaks=seq(Jag1min,Jag1max,(Jag1max-Jag1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 3",
       subtitle = "JAG1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Clust1_Jag1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="1",],aes(color=Jag1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Jag1med,
                        breaks=seq(Jag1min,Jag1max,(Jag1max-Jag1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 1",
       subtitle = "JAG1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))



Clust0_Jag1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="0",],aes(color=Jag1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Jag1med,
                        breaks=seq(Jag1min,Jag1max,(Jag1max-Jag1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 0",
       subtitle = "JAG1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Clust4_Jag1<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="4",],aes(color=Jag1),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Jag1med,
                        breaks=seq(Jag1min,Jag1max,(Jag1max-Jag1min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 4",
       subtitle = "JAG1 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Dll4min=min((scanpy[!is.na(scanpy$Dll4) ,])$Dll4)
Dll4max=max((scanpy[!is.na(scanpy$Dll4) ,])$Dll4)
Dll4med=median((scanpy[!is.na(scanpy$Dll4),])$Dll4)

Clust5_Dll4<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="5",],aes(color=Dll4),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Dll4med,
                        breaks=seq(Dll4min,Dll4max,(Dll4max-Dll4min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 5",
       subtitle = "DLL4 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))

Clust2_Dll4<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="2",],aes(color=Dll4),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Dll4med,
                        breaks=seq(Dll4min,Dll4max,(Dll4max-Dll4min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 2",
       subtitle = "DLL4 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Clust3_Dll4<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="3",],aes(color=Dll4),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Dll4med,
                        breaks=seq(Dll4min,Dll4max,(Dll4max-Dll4min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 3",
       subtitle = "DLL4 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Clust1_Dll4<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="1",],aes(color=Dll4),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Dll4med,
                        breaks=seq(Dll4min,Dll4max,(Dll4max-Dll4min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 1",
       subtitle = "DLL4 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))



Clust0_Dll4<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="0",],aes(color=Dll4),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Dll4med,
                        breaks=seq(Dll4min,Dll4max,(Dll4max-Dll4min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 0",
       subtitle = "DLL4 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


Clust4_Dll4<-ggplot(scanpy,aes(x=UMAP_C1,y=UMAP_C2,labels=sampleID))+
  #geom_point(color="grey",size=5.5)+
  geom_point(color="white",size=4,alpha=0.0)+
  geom_point(data=scanpy[scanpy$Louvain_cluster=="4",],aes(color=Dll4),size=4)+
  #scale_color_brewer(palette="Spectral")+
  scale_color_gradient2(low = "lightyellow", mid = "lightyellow", high ="red2", 
                        midpoint = Dll4med,
                        breaks=seq(Dll4min,Dll4max,(Dll4max-Dll4min)/4),na.value = "gray97")+
  theme_dark()+
  labs(title = "Cluster 4",
       subtitle = "DLL4 protein")+
  theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold"))


grid.arrange(Clust5_Dll4,Clust2_Dll4,Clust3_Dll4,
             Clust1_Dll4,Clust0_Dll4,Clust4_Dll4,ncol=6)

grid.arrange(Clust5_Notch,Clust3_Notch,Clust2_Notch,
             Clust1_Notch,Clust0_Notch,Clust4_Notch,
             Clust5_Dll4,Clust3_Dll4,Clust2_Dll4,
             Clust1_Dll4,Clust0_Dll4,Clust4_Dll4,
             ncol=6)

grid.arrange(Clust5_Jag1,Clust3_Jag1,Clust2_Jag1,
             Clust1_Jag1,Clust0_Jag1,Clust4_Jag1,
             Clust5_CD45,Clust3_CD45,Clust2_CD45,
             Clust1_CD45,Clust0_CD45,Clust4_CD45,
             ncol=6)



## ANIMATION
#BiocManager::install('gganimate')
library(gganimate)
library(magick)

ttime<-subset(scanpy,select=c(UMAP_C1,UMAP_C2,Jag1,Notch1,Gfi1,population,dpt_pseudotime))
ttime$time<-ttime$population
ttime$time<-gsub('Gfi1_HE','1',ttime$time)
ttime$time<-gsub('Gfi1_pos_IAHC','2',ttime$time)
ttime$time<-gsub('Gfi1_neg_IAHC','3',ttime$time)

ttime_ad<-subset(ttime,ttime$time==3)
ttime_ad$time<-gsub('3','4',ttime_ad$time)
ttime_ad$population<-"Gfi1_neg_IAHC_rep"

ttime<-rbind(ttime,ttime_ad)

ttime$time<-as.numeric(ttime$time)


a<-ggplot(ttime,aes(x=UMAP_C1,y=UMAP_C2))+
  #geom_point(color="grey",size=5.5)+
  #geom_point(color="white",size=5)+
  geom_point(aes(color=Jag1,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="red2", na.value = "grey",
                       midpoint = Jag1med,
                       breaks=seq(Jag1min,Jag1max,(Jag1max-Jag1min)/4))+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,5)}"))

Jag_gif <- animate(a, width = 200, height = 200,fps = 2)

Jag_gif

b<-ggplot(ttime,aes(x=UMAP_C1,y=UMAP_C2))+
  #geom_point(color="grey",size=5.5)+
  #geom_point(color="white",size=5)+
  geom_point(aes(color=Notch1,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="red2", na.value = "grey",
                        midpoint = Notch1med,
                        breaks=seq(Notch1min,Notch1max,(Notch1max-Notch1min)/4))+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,5)}"))

Notch_gif <- animate(b, width = 200, height = 200,fps = 2)

c<-ggplot(ttime,aes(x=UMAP_C1,y=UMAP_C2))+
  #geom_point(color="grey",size=5.5)+
  #geom_point(color="white",size=5)+
  geom_point(aes(color=Gfi1,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="red2", na.value = "grey",
                        midpoint = Gfi1med,
                        breaks=seq(Gfi1min,Gfi1max,(Gfi1max-Gfi1min)/4))+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,5)}"))

Gfi1_gif <- animate(c, width = 200, height = 200,fps = 2)


new_gif <- image_append(c(Gfi1_gif[1],Notch_gif[1],Jag_gif[1]))

for(i in 2:100){
  combined <- image_append(c(Gfi1_gif[i],Notch_gif[i], Jag_gif[i]))
  new_gif <- c(new_gif, combined)
}

new_gif

image_write(new_gif, "/Users/yguillen/Desktop/temp/scRNA_cam/Roshana_scRNA/Plots/protein_population.gif")


### Create dataframe as new variable for pseudotime based on population and jagged1
newpst<-subset(scanpy,select=c(population,Jag1,Notch1))
summary(newpst$Jag1)

newpst$orig<-ifelse(newpst$Jag1>11000 & newpst$population=="Gfi1_HE",
                         c("orig"),c("unknown"))


#merge with original number of cells
row.names(meta)
row.names(newpst)
newpst<-merge(newpst,meta,by="row.names",all=TRUE)

row.names(newpst)<-newpst$Row.names
newpst$Row.names<-NULL

newpst<-subset(newpst,select=c(orig,population.y))

write.table(newpst,"/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/metadata_origin_psdt.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)


# Or reverse, using Cluser 4 CD45 positive cells 
preHSCclas$preHSC
newpstCD45<-subset(preHSCclas,select=c(sampleID,preHSC))
newpstCD45$orig45<-ifelse(preHSCclas$preHSC=="T2",
                    c("orig45"),c("unknown"))

table(newpstCD45$orig45)

#merge with original number of cells
row.names(meta)
row.names(newpstCD45)
newpstCD45<-merge(newpstCD45,meta,by="row.names",all=TRUE)

row.names(newpstCD45)<-newpstCD45$Row.names
newpstCD45$Row.names<-NULL

newpstCD45<-subset(newpstCD45,select=c(orig45,population))


newpstCD45$orig45[is.na(newpstCD45$orig45)] <- "unknown"

write.table(newpstCD45,"/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/metadata_originCD45_psdt.txt",sep="\t",col.names = FALSE,row.names = TRUE,quote = FALSE)




# scatter plot of x and y variables
# color by groups
scatterPlot <- ggplot(newpst,aes(Notch1, Jag1, color=population)) + 
  geom_point() + 
  scale_color_manual(values = c('red','#E69F00','blue')) +
  geom_smooth(aes(color=population),method="lm",se=FALSE)+
  theme_minimal()+
  theme(legend.position=c(0,1), legend.justification=c(0,1))
scatterPlot

# Marginal density plot of x (top panel)
xdensity <- ggplot(newpst, aes(Notch1, fill=population)) + 
  geom_density(alpha=.3) + 
  scale_fill_manual(values = c('red','#E69F00','blue')) +
  theme_minimal()+
  theme(legend.position = "none")
xdensity

# Marginal density plot of y (right panel)
ydensity <- ggplot(newpst, aes(Jag1, fill=population)) + 
  geom_density(alpha=.3) + 
  scale_fill_manual(values = c('red','#E69F00','blue')) + 
  theme_minimal()+
  theme(legend.position = "none")+
  coord_flip()
ydensity

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))



# DENSITY PLOT FOR PSEUDOTIME
# orig is the variable for orig as HE and JAG1 POSITIVE
# orig45 is the variable for orig as cluster4 CD45 positive T2
origcd45<-read.delim('/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/metadata_originCD45_psdt.txt',header=FALSE)
psdt45<-read.delim('/Users/yguillen/Desktop/temp/scRNA_cam/scanpy_run/pseudotime_orig45.csv',sep=",",header=TRUE)

row.names(origcd45)<-origcd45$V1
origcd45$V1<-NULL
names(origcd45)<-c("orig45","population_rr")

row.names(psdt45)<-psdt45$X
psdt45$X<-NULL
names(psdt45)<-c("dpt_CD45_T2")

dentime<-subset(scanpy,select=c(th_NOTCH1,th_NOTCH2,th_DLL4,th_JAG1,th_GFI1,Louvain_cluster,population,dpt_pseudotime,recall))

# merge witg orig time
dentime<-merge(dentime,newpst,by="row.names")
row.names(dentime)<-dentime$Row.names
dentime$Row.names<-NULL
#merge with origcd45
dentime<-merge(dentime,psdt45,by="row.names")
row.names(dentime)<-dentime$Row.names
dentime$Row.names<-NULL
dentime<-merge(dentime,origcd45,by="row.names")
row.names(dentime)<-dentime$Row.names
dentime$Row.names<-NULL

dentime$Louvain_cluster <- factor(dentime$Louvain_cluster, levels = c("5","3","2","1","0","4"))

ggplot(dentime[!is.na(dentime$th_JAG1),],aes(x=1-dpt_CD45_T2))+
  geom_density(aes(fill=th_JAG1),size=1,adjust=1/2,alpha=0.3)+
  #geom_density(aes(color=as.factor(Louvain_cluster)),size=1)+
  #geom_rug(aes(color=th_JAG1),length = unit(0.1, "npc"))+
  scale_fill_brewer(palette="Set2")+
  facet_wrap(~as.factor(Louvain_cluster),ncol=1)+
  theme_light()+
  theme(legend.position = "bottom")


# bind dataframes
thNOTCH1<-dentime[!is.na(dentime$th_NOTCH1) & dentime$th_NOTCH1=="pos_NOTCH1",colnames(dentime) == "Louvain_cluster" | colnames(dentime) == "th_NOTCH1" | colnames(dentime) == "dpt_pseudotime" | colnames(dentime) =="dpt_CD45_T2"]
colnames(thNOTCH1)[1]<-"variable"

thNOTCH2<-dentime[!is.na(dentime$th_NOTCH2) & dentime$th_NOTCH2=="pos_NOTCH2",colnames(dentime) == "Louvain_cluster" | colnames(dentime) == "th_NOTCH2" | colnames(dentime) == "dpt_pseudotime" | colnames(dentime) =="dpt_CD45_T2"]
colnames(thNOTCH2)[1]<-"variable"

thDLL4<-dentime[!is.na(dentime$th_DLL4) & dentime$th_DLL4=="pos_DLL4",colnames(dentime) == "Louvain_cluster" | colnames(dentime) == "th_DLL4" | colnames(dentime) == "dpt_pseudotime" | colnames(dentime) =="dpt_CD45_T2"]
colnames(thDLL4)[1]<-"variable"

thJAG1<-dentime[!is.na(dentime$th_JAG1) & dentime$th_JAG1=="pos_JAG1",colnames(dentime) == "Louvain_cluster" | colnames(dentime) == "th_JAG1" | colnames(dentime) == "dpt_pseudotime" | colnames(dentime) =="dpt_CD45_T2"]
colnames(thJAG1)[1]<-"variable"

thGFI1<-dentime[!is.na(dentime$th_GFI1) & dentime$th_GFI1=="pos_GFI1",colnames(dentime) == "Louvain_cluster" | colnames(dentime) == "th_GFI1" | colnames(dentime) == "dpt_pseudotime" | colnames(dentime) =="dpt_CD45_T2"]
colnames(thGFI1)[1]<-"variable"

thall<-rbind(thNOTCH1,thNOTCH2)
thall<-rbind(thall,thDLL4)
thall<-rbind(thall,thJAG1)
thall<-rbind(thall,thGFI1)

library(ggnewscale)

ggplot(thall,aes(x=dpt_CD45_T2))+
  geom_density(color="black",size=1,adjust=0.5)+
  #geom_density(aes(color=as.factor(Louvain_cluster)),size=1)+
  geom_rug(aes(color=as.factor(Louvain_cluster)),length = unit(0.1, "npc"))+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(~variable,ncol=1)+
  theme_void()+
  theme(legend.position = "bottom")

## add gene expression
dentime<-merge(dentime,metext[,c(12,20:28,38:ncol(metext))],by=0)
row.names(dentime)<-dentime$Row.names
dentime$Row.names<-NULL

dentime<-merge(dentime,cyclondata,by=0)
row.names(dentime)<-dentime$Row.names

# order pseudotime factors clusters
dentime$Louvain_cluster<-as.factor(dentime$Louvain_cluster)
levels(dentime$Louvain_cluster)
dentime$Louvain_cluster <- factor(dentime$Louvain_cluster, levels = c("5","3","2","1","0","4"))

N1rib<-ggplot(dentime,aes(x=dpt_CD45_T2,y=Notch1.y))+
  geom_point(aes(shape=Plate),size=2,alpha=0.5)+
  geom_ribbon(aes(ymin=0, ymax=Notch1.y,group=as.factor(Louvain_cluster),color=as.factor(Louvain_cluster),fill=as.factor(Louvain_cluster)),alpha=0.2)+
  #geom_density(aes(color=as.factor(Louvain_cluster)),size=1)+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  facet_wrap(~as.factor(Louvain_cluster),ncol=1)+
  theme_bw()+
  theme(legend.position = "bottom")

Jag1rib<-ggplot(dentime,aes(x=dpt_CD45_T2,y=Jag1.y))+
  geom_point(aes(shape=Plate),size=2,alpha=0.5)+
  geom_ribbon(aes(ymin=0, ymax=Jag1.y,group=as.factor(Louvain_cluster),color=as.factor(Louvain_cluster),fill=as.factor(Louvain_cluster)),alpha=0.2)+
  #geom_density(aes(color=as.factor(Louvain_cluster)),size=1)+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  facet_wrap(~as.factor(Louvain_cluster),ncol=1)+
  theme_bw()+
  theme(legend.position = "bottom")

N2rib<-ggplot(dentime,aes(x=dpt_CD45_T2,y=Notch2.x))+
  geom_point(aes(shape=Plate),size=2,alpha=0.5)+
  geom_ribbon(aes(ymin=min((dentime[!is.na(dentime$Notch2.x),])$Notch2.x), ymax=Notch2.x,group=as.factor(Louvain_cluster),color=as.factor(Louvain_cluster),fill=as.factor(Louvain_cluster)),alpha=0.2)+
  #geom_density(aes(color=as.factor(Louvain_cluster)),size=1)+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  facet_wrap(~as.factor(Louvain_cluster),ncol=1)+
  theme_bw()+
  theme(legend.position = "bottom")

D4rib<-ggplot(dentime,aes(x=dpt_CD45_T2,y=Dll4.x))+
  geom_point(aes(shape=Plate),size=2,alpha=0.5)+
  geom_ribbon(aes(ymin=min((dentime[!is.na(dentime$Dll4.x),])$Dll4.x), ymax=Dll4.x,group=as.factor(Louvain_cluster),color=as.factor(Louvain_cluster),fill=as.factor(Louvain_cluster)),alpha=0.2)+
  #geom_density(aes(color=as.factor(Louvain_cluster)),size=1)+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  facet_wrap(~as.factor(Louvain_cluster),ncol=1)+
  theme_bw()+
  theme(legend.position = "bottom")

dentime_phase<-subset(dentime,select=c(dpt_CD45_T2,Louvain_cluster,cycle,
#                                       Jag1.x,Notch1.x,
#                                       Notch2.x,Dll4.x,
#                                       CD45,
                                       Jag1.y,Notch1.y,
#                                       Notch2.y,Dll4.y,
#                                      Ptprc,
                                       Gata2,Hes1,Hey1,Hey2,Hist1h1b))
dentime_phase_melt<-melt(dentime_phase,id.vars = c("dpt_CD45_T2","Louvain_cluster","cycle"))

genesel<-c("Jag1.x","Jag1.y",
"Notch1.x","Notch1.y",
"Notch2.x","Notch1.y",
"Dll4.x","Dll4.y",
"Gata2","Hes1","Hey1","Hey2","Ptprc","Hist1h1b")

# No bins
ggplot(dentime_phase_melt[!is.na(dentime_phase_melt$value),],aes(x=1-dpt_CD45_T2,y=value,group=as.factor(Louvain_cluster)))+
  geom_line(aes(color=as.factor(Louvain_cluster)),size=1)+
  scale_color_brewer(palette="Spectral")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=as.factor(Louvain_cluster),shape=cycle),size=2,alpha=0.8)+
  scale_color_brewer(palette="Spectral")+
#  geom_hline(data=dentime_phase_melt[!is.na(dentime_phase_melt$value) & dentime_phase_melt$variable=="Jag1.x",],aes(yintercept = 7500),color="red",size=1.5,linetype="dotted")+
  geom_hline(data=dentime_phase_melt[!is.na(dentime_phase_melt$value) & dentime_phase_melt$variable!="Jag1.x",],aes(yintercept = 1),color="red",size=1.5,linetype="dotted")+
  #geom_ribbon(aes(ymin=0, ymax=value,group=variable,color=variable,fill=variable),alpha=0)+
  #geom_density(aes(color=as.factor(Louvain_cluster)),size=1)+
  facet_wrap(variable~cycle,ncol=3,scales="free")+
#  geom_hline(aes(yintercept = 10000))+
  theme_bw()+
  theme(legend.position = "bottom")


#grid.arrange(N1rib,Jag1rib,N2rib,D4rib,ncol=4)


### Binning pseudotime
library(Hmisc)
library(data.table)
library(dplyr)

# Estimating current bins
dt <- as.data.table(dentime_phase)
setkey(dt,Louvain_cluster,dpt_CD45_T2)
dt[, diff := dpt_CD45_T2 - shift(dpt_CD45_T2, fill = first(dpt_CD45_T2)), by = Louvain_cluster]
summary(dt$diff)

# binning time windows with median
addin = 0.01

bintime <- dentime_phase_melt %>% 
  mutate(window = .$dpt_CD45_T2 %/% addin) %>% 
  group_by(window,Louvain_cluster,variable) %>% summarise(value_median = median(value)) %>%
  mutate(dpt_CD45_T2 = window*addin, stop=(window+1)*addin)
#  select(-window)

library(scales)

ggplot(bintime[!is.na(bintime$value_median),],aes(x=1-dpt_CD45_T2,y=value_median,group=as.factor(Louvain_cluster)))+
  geom_line(aes(color=as.factor(Louvain_cluster)),na.rm = TRUE,size=1)+
#  scale_color_brewer(palette="Set2")+
#  ggnewscale::new_scale_color()+
  geom_point(na.rm = TRUE,color="black",size=3)+
  geom_point(na.rm = TRUE,aes(color=as.factor(Louvain_cluster)),size=2)+
  #geom_ribbon(aes(ymin=0, ymax=value,group=variable,color=variable,fill=variable),alpha=0)+
  #geom_density(aes(color=as.factor(Louvain_cluster)),size=1)+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(~variable,ncol=1,scales="free_y")+
  geom_hline(data=bintime[!is.na(bintime$value_median) & bintime$variable=="Jag1.x",],aes(yintercept = 7500),color="red",size=1.5,linetype="dotted")+
  geom_hline(data=bintime[!is.na(bintime$value_median) & bintime$variable!="Jag1.x",],aes(yintercept = 1),color="red",size=1.5,linetype="dotted")+
  #geom_hline(data=bintime[!is.na(bintime$value_median) & bintime$variable=="Dll4.x",],aes(yintercept = 120),color="red",size=1.5,linetype="dotted")+
  theme_bw()+
  scale_x_continuous(breaks=pretty_breaks(n=40))+
  theme(legend.position = "bottom",axis.text.x = element_blank())


### EXPRESSION BY CELL CYCLE PHASE
# no bin windows


phase_exp<-metext[,c(38:ncol(metext))]
phase_exp<-merge(cyclondata,phase_exp,by="row.names")
row.names(phase_exp)<-phase_exp$Row.names
phase_exp$Row.names<-NULL
phase_exp$sampleID<-NULL

phase_exp_melt<-melt(phase_exp,id.vars = "cycle")


phase_exp_melt_sub<-phase_exp_melt[phase_exp_melt$variable %in% sample(colnames(phase_exp)[2:length(colnames(phase_exp))],100),]

#or
dentime_phase_melt

ggplot(dentime_phase_melt,aes(x=cycle,y=value))+
  #geom_point(aes(color=as.factor(Louvain_cluster)),size=2)+
  geom_boxplot(aes(color=as.factor(Louvain_cluster)),alpha=0,size=1)+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(~variable,scales="free",nrow=3)+
  theme_bw()+
  theme(legend.position = "bottom")




# FOR GENE EXPRESSION animation, pseudotime orig Jag1+ HE
ttime_gene<-subset(metext,select=c(UMAP_C1,UMAP_C2,Notch3,Gata2,Runx1,population,dpt_pseudotime,Hey1,Hey2,Hes1,Nfkbia))
# merge with cyclone data
row.names(cyclondata)<-cyclondata$sampleID
ttime_gene<-merge(cyclondata,ttime_gene,by="row.names")
row.names(ttime_gene)<-ttime_gene$Row.names
ttime_gene<-ttime_gene[,-c(1:2)]

ttime_gene$time<-ttime_gene$population
ttime_gene$time<-gsub('Gfi1_HE','1',ttime_gene$time)
ttime_gene$time<-gsub('Gfi1_pos_IAHC','2',ttime_gene$time)
ttime_gene$time<-gsub('Gfi1_neg_IAHC','3',ttime_gene$time)

ttime_gene_ad<-subset(ttime_gene,ttime_gene$time==3)
ttime_gene_ad$time<-gsub('3','4',ttime_gene_ad$time)
ttime_gene_ad$population<-"Gfi1_neg_IAHC_rep"

ttime_gene<-rbind(ttime_gene,ttime_gene_ad)

ttime_gene$time<-as.numeric(ttime_gene$time)


cyg<-ggplot(ttime_gene,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=cycle,group = population),size=4)+
  scale_color_brewer(palette="Pastel1")+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,1)}"))

cyclone_gif <- animate(cyg, width = 200, height = 200,fps = 2)

popg<-ggplot(ttime_gene,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=population,group = population),size=4)+
  scale_color_brewer(palette="Pastel1")+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,1)}"))

population_gif <- animate(popg, width = 200, height = 200,fps = 2)



tN1<-ggplot(ttime_gene,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=Hey1,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="deeppink", na.value = "grey")+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,1)}"))

tN1_gif <- animate(tN1, width = 200, height = 200,fps = 2)

tN2<-ggplot(ttime_gene,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=Hes1,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="deeppink", na.value = "grey")+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,1)}"))

tN2_gif <- animate(tN2, width = 200, height = 200,fps = 2)

tN3<-ggplot(ttime_gene,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=Hey2,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="deeppink", na.value = "grey")+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,1)}"))

tN3_gif <- animate(tN3, width = 200, height = 200,fps = 2)


gnotch_gif <- image_append(c(tN1_gif[1],tN3_gif[1],tN2_gif[1]))

for(i in 2:100){
  combined <- image_append(c(tN1_gif[i],tN3_gif[i],tN2_gif[i]))
  gnotch_gif <- c(gnotch_gif, combined)
}

gnotch_gif


image_write(gnotch_gif, "/Users/yguillen/Desktop/temp/scRNA_cam/Roshana_scRNA/Plots/target_genes_population.gif")


image_write(cyclone_gif, "/Users/yguillen/Desktop/temp/scRNA_cam/Roshana_scRNA/Plots/cycle_population.gif")
image_write(population_gif, "/Users/yguillen/Desktop/temp/scRNA_cam/Roshana_scRNA/Plots/pst_population.gif")




ag<-ggplot(ttime_gene,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=Gata2,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="deeppink", na.value = "grey")+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,1)}"))

Gata2_gif <- animate(ag, width = 200, height = 200, fps=2)


bg<-ggplot(ttime_gene,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=Notch3,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="deeppink", na.value = "grey")+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,1)}"))

Notch3_gif <- animate(bg, width = 200, height = 200,fps=2)

cg<-ggplot(ttime_gene,aes(x=UMAP_C1,y=UMAP_C2))+
  geom_point(aes(color=Runx1,group = population),size=5)+
  scale_color_gradient2(low = "white", mid = "lightyellow", high ="deeppink", na.value = "grey")+
  theme_dark()+
  theme(legend.position = "bottom",plot.title = element_text(size = 12, face = "bold"))+
  transition_time(dpt_pseudotime)+
  shadow_mark()+
  labs(title = paste("Cells:", "{round(frame_time,1)}"))

Runx1_gif <- animate(cg, width = 200, height = 200,fps=2)


gnew_gif <- image_append(c(Notch3_gif[1],Gata2_gif[1],Runx1_gif[1]))

for(i in 2:100){
  combined <- image_append(c(Notch3_gif[i],Gata2_gif[i], Runx1_gif[i]))
  gnew_gif <- c(gnew_gif, combined)
}

gnew_gif


image_write(gnew_gif, "/Users/yguillen/Desktop/temp/scRNA_cam/Roshana_scRNA/Plots/gene_population.gif")


library(pheatmap)

### HEATMAP WITH SORTED CELLS BASED ON JAG1 OR GFI1?

namesel<-c("population","Louvain_cluster","Jag1.x","Notch1.x","Dll4.x","Notch2.x","Gfi1.x","CD41","CD31","CD45")

namesel<-c("population","Louvain_cluster",
           "Jag1.y","Notch1.y","Dll4.y","Notch2.y","Itga2b","Pecam1","Ptprc",
           "Gata2","Runx1","Sox7","Hes1","Hey1","Hey2")


heatmark<-metext[metext$population=="Gfi1_HE",colnames(metext) %in% namesel]
heatmark<-heatmark[,c(3:ncol(heatmark))]
heatmark_t<-t(heatmark)

# Scale by each element
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

heatmark_norm <- t(apply(heatmark_t, 1, cal_z_score))
colnames(heatmark_norm)<-colnames(heatmark_t)
colnames(heatmark_norm)<-gsub('^.*_','',colnames(heatmark_norm))

my_sample_col <- data.frame(c(metext[metext$population=="Gfi1_HE",]$population),(metext[metext$population=="Gfi1_HE",]$Louvain_cluster))
names(my_sample_col)<-c("population","cluster")

my_sample_col$cluster<-paste(my_sample_col$population,my_sample_col$cluster)

row.names(my_sample_col) <- colnames(heatmark_norm)

heatmark_norm <- heatmark_norm[ , rev(order( heatmark_norm[row.names(heatmark_norm) == "Jag1.y"])) ]

pheatmap(heatmark_norm,
         annotation_col = my_sample_col,
         fontsize =10,
         Colv=FALSE,dendrogram="row")

######### PHEATMAP WITH GENES DOWNREGULATED IN IKB KO #####
# From RNASeq_HSC_IKB.R script
UPgenes
DOWNgenes

namesel<-c("population","Louvain_cluster")
namesel<-c(namesel,UPgenes)


heatmark<-metext[,colnames(metext) %in% namesel]
heatmark<-heatmark[,c(3:ncol(heatmark))]


#for expressed genes
exp<-colnames(heatmark[, which(numcolwise(sum)(heatmark) !=0)])
exp

heatmark<-heatmark[,colnames(heatmark) %in% exp]

heatmark_t<-t(heatmark)

# Scale by each element
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

heatmark_norm <- t(apply(heatmark_t, 1, cal_z_score))
colnames(heatmark_norm)<-colnames(heatmark_t)
colnames(heatmark_norm)<-gsub('^.*_','',colnames(heatmark_norm))

my_sample_col <- data.frame(population=metext$population,
                            Louvain_cluster=metext$Louvain_cluster)

my_sample_col$Louvain_cluster<-as.factor(my_sample_col$Louvain_cluster)

row.names(my_sample_col) <- colnames(heatmark_norm)

#heatmark_norm <- heatmark_norm[ , rev(order( heatmark_norm[row.names(heatmark_norm) == "Jag1.y"])) ]

#heatmark_norm[, order(my_sample_col$Louvain_cluster)]
#heatmark_norm<-heatmark_norm[,row.names(my_sample_col[order(my_sample_col$Louvain_cluster),])]


pheatmap(heatmark_norm,
         annotation_col = my_sample_col,
         cutree_rows = 4,
       #  cutree_cols = 6,
         fontsize =4,
         Colv=FALSE,
         dendrogram="row",
         show_colnames = FALSE)


