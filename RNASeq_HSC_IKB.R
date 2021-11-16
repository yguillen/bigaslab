
### Analyze RNASeq data from fetal liver RBP KO, NICD and WT in HSCs

library(DESeq2)
library(biomaRt)
library(geneplotter)
library(MDplot)
library(lattice)
library(genefilter)
library(plotly)
library(limma)
library(ggrepel)
#library(pathfindR)
library(devtools)
#library(bc3net)

setwd("/Users/yguillen/Desktop/temp/HSC_IKB/HTSEQ/")

output.dir="/Users/yguillen/Desktop/temp/HSC_IKB/HTSEQ"

# Import metadata files
metadata_all<-read.delim("metadata_all.txt",header = FALSE)
colnames(metadata_all)<-c("sampleID","countFile","condition","gender")


metadata_sub<-metadata_all[metadata_all$sampleID!="KO3",]
metadata_sub_wt<-metadata_sub[metadata_sub$condition=="WT",]
metadata_sub_he<-metadata_sub[metadata_sub$condition=="HE",]
metadata_sub_noko<-metadata_sub[metadata_sub$condition!="KO",]
metadata_sub_fem<-metadata_sub[metadata_sub$gender=="F",]
metadata_sub_mal<-metadata_sub[metadata_sub$gender=="M",]

meta_wt_ko<-metadata_sub[metadata_sub$condition=="WT" | metadata_sub$condition=="KO",]

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata_sub_fem,
                                          directory = output.dir,
                                          design = ~ condition)

rowData(DESeq2Table)
mcols(rowData(DESeq2Table))

# Add additional annotation information using biomaRt
ensembl76 <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

bmm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","mgi_symbol","gene_biotype","chromosome_name"),
            values= rownames(DESeq2Table),
            #            filter="external_gene_name",
            mart=ensembl76) 

head(bmm)

#Add description data to gene counts
DESeq2Features <- data.frame(id = rownames(DESeq2Table))
DESeq2Features$id <- as.character(DESeq2Features$id)


### join them together

#if htseqgene symbol
rowData <- merge(DESeq2Features, bmm, by.x="id",by.y = "external_gene_name")

#if not
#rowData <- merge(DESeq2Features, bm, by.x="id",by.y = "ensembl_gene_id")


rowData <- as(rowData, "DataFrame")

### add the annotation to the DESeq2 table
#mcols(rowData(DESeq2Table)) <- c(mcols(rowData(DESeq2Table)),rowData)

DESeq2Table
colnames(DESeq2Table)
rownames(DESeq2Table)

## Quality control and normalization of counts

#how many genes we capture, counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(DESeq2Table)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

### random sample from the count matrix
nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]


### NORMALIZATION
# Remove genes with low expression levels ()
keep <- rowSums(counts(DESeq2Table)) >= 10
DESeq2Table <- DESeq2Table[keep,]

GeneCounts <- counts(DESeq2Table)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

### make sure to get fold change WT-deletion
# for all
colData(DESeq2Table)$gender <- factor(colData(DESeq2Table)$gender, levels = c("F","M"))
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","HE","KO"))

colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","KO"))
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","HE"))

#### estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)


### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Table, blind=TRUE) ## create a distance matrix between the samples


## Create PCA (check if for all or for RA sign)
pcaData <- plotPCA(rld, intgroup=c("condition","gender"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(aes(shape=condition),color="black",size=10,alpha=0.3) +
  geom_point(aes(color=condition, fill=condition, shape=condition),size=8,alpha=0.8) +
  stat_ellipse(aes(color=condition),level=0.8)+
  #scale_color_manual(values=c("dodgerblue","coral"))+
  #scale_color_manual(values=c("grey","dodgerblue","coral","darkolivegreen2","black","dodgerblue4","coral4"))+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #facet_grid(~gender)+
  scale_shape_manual(values=c(16,15,25))+
  scale_color_manual(values=c("black","black","red"))+
  scale_fill_manual(values=c("black","black","red"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1))


ggplot(pcaData, aes(y=PC1, x=condition,label=name)) +
  geom_point(aes(color=condition),size=8) +
  geom_boxplot(alpha=0)+
  #stat_ellipse(aes(color=condition),level=0.8)+
  scale_color_manual(values=c("black","black","red"))+
  scale_shape_manual(values=c(15,17,16))+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
  xlab("")+
  #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #scale_color_manual(values=c("purple","red"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1))


###### Differential EXPRESSION ANALYSIS ###

DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)


# Statistical testing of DE genes
levels(DESeq2Table$condition)

design(DESeq2Table)

dds<-DESeq(DESeq2Table)
resultsNames(dds)

rescondition<-results(dds)


# Example plot counts
plotCounts(dds, gene=which(rownames(rescondition)=="Gata2"), intgroup="condition")
head(countab)


head(rescondition)
summary(rescondition)
sum(rescondition$pvalue < 0.05,na.rm=TRUE)
sum(rescondition$padj < 0.1,na.rm=TRUE)
plotMA(rescondition)

# Selecting specific DEGs
#To dataframe (ckeck females, males or gender)
res_df_rna<-as.data.frame(rescondition)
class(rowData)
rowData<-as.data.frame(rowData)

# CHECK gene  names
#res_df_rna<-merge(rowData,res_df_rna,by.x="id",by.y="row.names",all.y=TRUE)
res_df_rna<-merge(rowData,res_df_rna,by.x="id",by.y="row.names",all.y=TRUE)

res_df_sig_rna<-subset(res_df_rna,res_df_rna$padj<0.1)
res_df_sig_rna <- res_df_sig_rna[order(res_df_sig_rna$padj),]


write.table(res_df_sig_rna,"/Users/yguillen/Desktop/temp/HSC_IKB/DEGs_sig.txt",quote = FALSE,sep="\t")

# UP and DOWN genes
UPgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange>0))$id)
length(UPgenes)
updf<-data.frame(fc="UP",
                 Gene=UPgenes)
DOWNgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange<0))$id)
length(DOWNgenes)
dodf<-data.frame(fc="DOWN",
                 Gene=DOWNgenes)
fcdf<-rbind(updf,dodf)

# Heatmap
#Manual heatmap from log trasformed data, using all samples
rld_df<-as.data.frame(assay(rld))

#only degs
deggenes<-unique(res_df_sig_rna$id)

deg_rld<-which(rownames(rld_df) %in% deggenes)
deg_rld_df<-rld_df[deg_rld,]

deg_rld<-merge(deg_rld_df,res_df_rna[,c(1,7:12)],by.x="row.names",by.y="id")

deg_rld<- deg_rld[order(deg_rld$log2FoldChange),]
deg_rld$Row.names

library(dplyr)
library(tidyverse)

deg_rld<-deg_rld[!duplicated(deg_rld), ]

row.names(deg_rld)<-deg_rld$Row.names

heatmap.2( as.matrix(deg_rld[,2:6]), scale="row", dendrogram="column",trace="none", Rowv=FALSE,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           cexCol=1)


# pheatmap
expfem<-as.matrix(deg_rld[,2:6])

anikb<-subset(metadata_sub_fem,select=c("sampleID","condition"))
row.names(anikb)<-anikb$sampleID

row.names(fcdf)<-fcdf$Gene
fcdf$Gene<-NULL

ann_colors = list(
  condition = c("white", "firebrick"),
  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)

library(pheatmap)
pheatmap(expfem,
         annotation_col = anikb,
         annotation_row = fcdf,
         cutree_rows = 2,
         cutree_cols = 2,
         fontsize =5,
       scale="row",
       #  Colv=FALSE,
       #  cluster_cols = FALSE,
       #  show_colnames = FALSE,
         dendrogram="row")


###### Save DEGs list

#all genes
deggenes<-unique(res_df_rna$id)
deg_rld<-which(rownames(rld_df) %in% deggenes)
deg_rld_df<-rld_df[deg_rld,]

deg_rld<-merge(deg_rld_df,res_df_rna,by.x="row.names",by.y="id")


#######
# EXTRACT AND PLOT NORMALIZED COUNTS
countab<-as.data.frame(counts(dds, normalized=TRUE))

gseacount<-countab[,c(5,1,2,3,4)]
row.names(gseacount)<-toupper(row.names(gseacount))
write.table(gseacount,"/Users/yguillen/Desktop/temp/HSC_IKB/GSEA/countab_norm.txt",quote = FALSE,sep="\t")

countab<-t(countab)
countab<-melt(countab)
colnames(countab)<-c("Sample","Gene","value")
countab$group<-gsub('[0-9]','',countab$Sample)

countab$clas<-gsub('WT','IKB',countab$group)
countab$clas<-gsub('HE','IKB',countab$clas)

countab$group <- factor(countab$group, levels = c("WT","HE","KO"))

countab<-merge(fcdf,countab,by="Gene",all.y="TRUE")

listgenes<-c("Nfkbia","Gata2","Mfng","Neurl3","Sox18","Hhex","Nanog","Smed1","Ctnnb1","Gfi1b","Adgrg1","Ptprc",
             "Dnmt3a","Dnmt3b","Meis1","Hes1","Jag1")

listgenes<-c("Gata2","Neurl3","Sox18","Adgrg1","Ptprc","Smad6")

ggplot(countab[countab$Gene %in% listgenes,],aes(x=clas,y=value))+
  geom_point(aes(shape=group),color="black",size=5)+
  geom_point(aes(shape=group,color=group,fill=group),size=4)+
  scale_shape_manual(values=c(16,15,25))+
  scale_color_manual(values=c("black","black","red"))+
  scale_fill_manual(values=c("black","black","red"))+
  ylab("Normalized counts")+
  xlab('')+
  facet_wrap(fc~Gene,scales = "free")+
  theme_classic()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1))

write.table(countab,"/Users/yguillen/Desktop/temp/HSC_IKB/DEGs_fc.tab",quote = FALSE)





library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2019_Mouse",
         "KEGG_2019_Mouse")


enrichedU <- enrichr(UPgenes, dbs)
enrichedD <- enrichr(DOWNgenes, dbs)

Unrich <- enrichedU[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
Dnrich <- enrichedD[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]

Unrich$group<-c("Up")
Dnrich$group<-c("Down")

allGO<-rbind(Unrich,Dnrich)

bpsub<-subset(allGO,allGO$Adjusted.P.value<=0.05)
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
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[order(bpsub$Combined.Score), ]$Term))

ggplot(bpsub,aes(y=Combined.Score,x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group),color="black")+
  #geom_text(aes(label=tolower(Genes)), size=4,hjust=0.5,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Blues")+
  facet_wrap(~group,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=6, hjust = 1),
        axis.text.y = element_text(size=9, hjust = 1))


### RNASeq including fetal liver samples from CRuiz


setwd("/Users/yguillen/Desktop/temp/HSC_IKB/HTSEQ/HTSEQ_IKB_RBPJ_FL/")

aldir="/Users/yguillen/Desktop/temp/HSC_IKB/HTSEQ/HTSEQ_IKB_RBPJ_FL/"

# Import metadata files
metadata_fl<-read.delim("metadata_merge_fetal.txt",header = FALSE)
colnames(metadata_fl)<-c("sampleID","countFile","condition","gender","batch")


metadata_fl<-metadata_fl[metadata_fl$sampleID!="KO3",]

DESeq2Tableall <- DESeqDataSetFromHTSeqCount(sampleTable = metadata_fl,
                                          directory = aldir,
                                          design = ~ condition)

rowData(DESeq2Tableall)
mcols(rowData(DESeq2Tableall))

#Add description data to gene counts
DESeq2Featuresall <- data.frame(id = rownames(DESeq2Tableall))
DESeq2Featuresall$id <- as.character(DESeq2Featuresall$id)


### join them together

#if htseqgene symbol
rowData <- merge(DESeq2Featuresall, bm, by.x="id",by.y = "external_gene_name")
rowData <- as(rowData, "DataFrame")
DESeq2Tableall
DESeq2Tableall<-DESeq2Tableall[rownames(DESeq2Tableall) %in% res_df_sig_rna$id,]

colnames(DESeq2Tableall)
rownames(DESeq2Tableall)

## Quality control and normalization of counts

#how many genes we capture, counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(DESeq2Tableall)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

### random sample from the count matrix
nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]


### NORMALIZATION
# Remove genes with low expression levels ()
keep <- rowSums(counts(DESeq2Tableall)) >= 10
DESeq2Tableall <- DESeq2Tableall[keep,]

GeneCounts <- counts(DESeq2Tableall)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

### make sure to get fold change WT-deletion
# for all
colData(DESeq2Tableall)$condition <- factor(colData(DESeq2Tableall)$condition, levels = c("WT","HE","IKB_KO","RBP_KO","NICD"))

#### estimate size factors
DESeq2Tableall <- estimateSizeFactors(DESeq2Tableall)
sizeFactors(DESeq2Tableall)




### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Tableall, blind=TRUE) ## create a distance matrix between the samples



## Create PCA (check if for all or for RA sign)
pcaData <- plotPCA(rld, intgroup=c("condition","gender","batch"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2,label=name)) +
  geom_point(aes(shape=batch),color="black",size=10,alpha=0.3) +
  geom_point(aes(color=condition, shape=batch),size=8,alpha=0.8) +
  #stat_ellipse(aes(color=condition),level=0.8)+
  #scale_color_manual(values=c("dodgerblue","coral"))+
  #scale_color_manual(values=c("grey","dodgerblue","coral","darkolivegreen2","black","dodgerblue4","coral4"))+
  geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #facet_grid(~gender)+
  scale_color_manual(values=c("dodgerblue2","coral","darkolivegreen","grey","aquamarine4"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed()





