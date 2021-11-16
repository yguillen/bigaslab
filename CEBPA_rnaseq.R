
install.packages(c("FactoMineR", "factoextra"))

### Analyze RNASeq data from mieloid patients with different mutations

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
library(readxl)
library(FactoMineR)
library(factoextra)
library(enrichR)
library(reshape2)

set.seed(100)

#setwd("/Users/yolanda_guillen/Desktop/IMIM/CEBPA/htseq/")
setwd("/Volumes/cancer/CEBPA/htseq/")

input.dir="/Volumes/cancer/CEBPA/htseq/"
output.dir="/Volumes/cancer/CEBPA/results/"

# Import metadata files
metadata<-read.delim("../metadata/ID_sample.txt",header = TRUE)
colnames(metadata)
row.names(metadata)<-metadata$Sample
metadata$Sample<-NULL

countf<-read.delim("../metadata/met.txt",header = TRUE)
row.names(countf)<-countf$Sample

metadata<-merge(metadata,countf,by=0)
row.names(metadata)<-metadata$Row.names
metadata$Row.names<-NULL

mutdata<-as.data.frame(read_excel("../metadata/CEBPA_panel_trans.xlsx"))
class(mutdata)
row.names(mutdata)<-mutdata$Sample
mutdata$Sample<-NULL
mutdata[is.na(mutdata)]<-0


metadata<-merge(metadata,mutdata,by=0,all.x=TRUE)
row.names(metadata)<-metadata$Row.names
metadata$Row.names<-NULL


inmuno<-as.data.frame(read_excel("../metadata/RNAseq_immunophenotype_data.xlsx"))
rownames(inmuno)<-inmuno$ID
inmuno$ID<-NULL

metadata<-merge(metadata,inmuno,by="row.names",all.x=TRUE)
row.names(metadata)<-metadata$Row.names
metadata$Row.names<-NULL

# Input CIBERSORT WT1
cibersort_wt1<-read.delim("../results/CIBERSORT.Output_WT1_norm.txt")
row.names(cibersort_wt1)<-cibersort_wt1$Input.Sample
colnames(cibersort_wt1)[1]<-"Sample"

# Input CIBERSORT CEBPA
cibersort_cebpa<-read.delim("../results/CIBERSORT.Output_CEBPA_norm.txt")
row.names(cibersort_cebpa)<-cibersort_cebpa$Input.Sample
colnames(cibersort_cebpa)[1]<-"Sample"

cibersort<-rbind(cibersort_cebpa,cibersort_wt1)

metadata<-merge(metadata,cibersort,by="Sample")

metadata<-metadata[,c(1,4,2,3,23:26,5:21,27:ncol(metadata))]


metadata$exp<-ifelse(metadata$Pheno=="CEBPA","Dx_CEBPA","WT1_Rem")

metadata$Positivity<-gsub('-','minus',metadata$Positivity)
metadata$Positivity<-gsub('\\+','plus',metadata$Positivity)

metadata$Pheno_WT1<-ifelse(metadata$Pheno!=">100copias WT1",
                           "less100copias",
                           "more100copias")
metadata$Pheno_WT1<-as.factor(metadata$Pheno_WT1)


metadata$WT1mut<-ifelse(metadata$WT1>0,"WT1_mut","no_WT1_mut")

# DESeq tables for all samples

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata,
                                          directory = input.dir,
                                          design = ~ Pheno_WT1)

# DESeq tables for Dx samples

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata[metadata$exp=="Dx_CEBPA",],
                                          directory = input.dir,
                                          design = ~ WT1mut)
# DESeq tables for treated patients

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata[metadata$exp=="WT1_Rem",],
                                          directory = input.dir,
                                          design = ~ Pheno_WT1)

rowData(DESeq2Table)
mcols(rowData(DESeq2Table))

# Add additional annotation information using biomaRt
ensembl<-useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", 
                 dataset="hsapiens_gene_ensembl",
                 mirror="asia")

bm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","gene_biotype","chromosome_name"),
            values= rownames(DESeq2Table),
            filter="ensembl_gene_id",
            mart=ensembl) 

head(bm)

#Add description data to gene counts
DESeq2Features <- data.frame(id = rownames(DESeq2Table))
DESeq2Features$id <- as.character(DESeq2Features$id)


### join them together
rowData <- merge(DESeq2Features, bm, by.x="id",by.y = "ensembl_gene_id")
rowData <- as(rowData, "DataFrame")

DESeq2Table
rownames(DESeq2Table)
colnames(DESeq2Table)

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

# reorder factors for WT1 data
colData(DESeq2Table)$Pheno <- factor(colData(DESeq2Table)$Pheno, levels = c("CEBPA","WT1","10 copias WT1","10-100copias WT1",">100copias WT1"))

colData(DESeq2Table)$WT1mut <- factor(colData(DESeq2Table)$WT1mut, levels = c("WT1_mut","no_WT1_mut"))

#### estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)


### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Table, blind=TRUE) ## create a distance matrix between the samples


#OUTPUT for cibersort WT1 or CEBPA
#####

genetab<-as.data.frame(counts(DESeq2Table, normalized=TRUE))
row.names(genetab)
row.names(bm)<-bm$ensembl_gene_id
genetab<-merge(bm,genetab,by="row.names")
row.names(genetab)<-genetab$Row.names
genetab<-genetab[,c(4,7:ncol(genetab))]

genetab<-genetab[!duplicated(genetab$external_gene_name),]
genetab<-genetab[genetab$external_gene_name != "",]
# FOR LM22 gene signature use gene symbols
row.names(genetab)<-genetab$external_gene_name
genetab$external_gene_name<-NULL

# FOR custom T cell populations gene signature use ensembl ids. don't change row.names to external gene names.
#genetab$external_gene_name<-NULL

#reduce size
keep <- rowSums(genetab) >= 50
genetab<- genetab[keep,]

#write.table(genetab,"../results/norm_exp_symbol_WT1.tab",sep="\t",quote = FALSE)
write.table(genetab,"../results/norm_exp_symbol_CEBPA.tab",sep="\t",quote = FALSE)
#####

## Create PCA
pcaData <- plotPCA(rld, intgroup=c("Pheno_WT1"), returnData=TRUE)

pcaData <- plotPCA(rld, intgroup=c("WT1mut","Positivity"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
#pcaData[is.na(pcaData$WT1),]$WT1<-"no_data"


ggplot(pcaData, aes(x=PC1, y=PC2,label=name)) +
  geom_point(color="black",size=10,alpha=0.3) +
  geom_point(aes(color=WT1mut),size=8,alpha=0.8) +
  geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
#  stat_ellipse(aes(color=Pheno),level=0.8)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("dodgerblue2","darkolivegreen","coral","red"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed()


ggplot(pcaData, aes(y=PC1, x=Positivity,label=name)) +
  geom_point(aes(color=Positivity),size=8) +
  geom_boxplot(alpha=0)+
  #stat_ellipse(aes(color=condition),level=0.8)+
  scale_color_manual(values=c("dodgerblue2","darkolivegreen","coral","red"))+
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

# HEATMAP WITH METADATA INFO

#### THIS STEP ONLY TO ESTIMATE VARIANCE
#df with log transformed values
take<-rowVars( assay(rld)) >= 0.5
genevar<-row.names(assay(rld)[take,])
dim(assay(rld)[take,])
rld_df<-as.data.frame(assay(rld)[take,])



library(RColorBrewer)
library(pheatmap)
cols <- colorRampPalette(rev(brewer.pal(11,name="RdBu")))(120)

paletteLength<-length(cols)

df<-scale(t(rld_df))

#hist(df)

row.names(metadata)<-metadata$Sample
df<-merge(metadata,df,by="row.names")
row.names(df)<-df$Row.names
df$Row.names<-NULL

df[is.na(df)]<-0

# Convert to numeric

sapply(df[,c(8:50)], class)
cols.num<-colnames(df[,c(8:50)])
df[cols.num] <- sapply(df[cols.num],as.numeric)
sapply(df[,c(8:50)], class)


colnames(df)
geneannot<-data.frame(id=colnames(df[,c(54:ncol(df))]))
geneannot<-merge(geneannot,rowData,by="id",all.x=TRUE)
row.names(geneannot)<-geneannot$id


ggplot(df,aes(x=Pheno,y=scale(ENSG00000184937)))+
  geom_jitter(aes(color=Pheno))+
  scale_color_brewer(palette="Spectral")+
  geom_boxplot(aes(fill=Pheno),alpha=0.5,outlier.shape = NA)+
  scale_fill_brewer(palette="Spectral")+
  theme_classic()

phet<-pheatmap(t(df[,c(54:ncol(df))]),
         annotation_col = as.data.frame(df[,c(which( colnames(df)==c("ENSG00000184937")),8,7,6)]),show_colnames= T,
         show_rownames = F,
         fontsize = 6,
         cutree_cols = 4,
         cutree_rows = 8,
         color=cols)

clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=8)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=4)))

clusters$ensembl_gene_id<-row.names(clusters)
clusters<-merge(clusters,bm,by='ensembl_gene_id',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==4,]$external_gene_name

# Order genes (clusters)
geneorder<-data.frame(ensembl_gene_id=rownames(as.data.frame(t(df[,c(44:ncol(df))]))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$ensembl_gene_id
                                                                                                                                   
clusters<-merge(clusters,geneorder,by="ensembl_gene_id")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(df[,c(44:ncol(df))]))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL
patientorder$order<-as.numeric(patientorder$order)
                                                                                                                                        
ciber_df<-df[,c(26:47)]
ciber_df$Sample<-row.names(ciber_df)
ciber_df<-melt(ciber_df,id.vars = "Sample")

ciber_df$Sample <- factor(ciber_df$Sample, levels = row.names(patientorder))
coli <- colorRampPalette(rev(brewer.pal(11,name="Dark2")))(30)

ggplot(ciber_df,aes(y=value,x=Sample,fill=variable),alpha = 0.5)+
  geom_bar(stat = "identity",aes(fill=variable))+
  #geom_text(aes(label=tolower(Genes)), size=4,hjust=0.5,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  scale_fill_manual(values=coli)+
  #facet_wrap(~group,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.position = "top",
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=6, hjust = 1),
        axis.text.y = element_text(size=9, hjust = 1))


#Per mutations panel 
#Expression of mutated genes
genes<-bm[bm$external_gene_name %in% colnames(df[,c(9:25)]),]$ensembl_gene_id
genetab<-merge(assay(rld[row.names(rld) %in% genes,]),bm,by.x="row.names",by.y="ensembl_gene_id")
row.names(genetab)<-genetab$external_gene_name
genetab$Row.names<-NULL
row.names(genetab)<-paste(row.names(genetab),"exp",sep="_")

colcol<-as.data.frame(df[,c(9:25)])
colcol<-merge(as.data.frame(scale(t(genetab[,1:13]))),colcol,by="row.names")
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL

sapply(colcol[,c(1:ncol(colcol))], class)

p1<-pheatmap(t(df[,c(54:ncol(df))]),
         annotation_col = colcol[,c(1:17)],
         show_colnames= T,
         show_rownames = F,
         fontsize = 6,
         cutree_cols = 4,
         cutree_rows = 8,
         color=cols)

p2<-pheatmap(t(df[,c(54:ncol(df))]),
             annotation_col = colcol[,c(18:ncol(colcol))],
             show_colnames= T,
             show_rownames = F,
             fontsize = 6,
             cutree_cols = 4,
             cutree_rows = 8,
             color=cols)


# Extract genes and patients clusters
clusters<-data.frame(cluster=sort(cutree(p1$tree_row, k=8)))
patients<-data.frame(cluster=sort(cutree(p1$tree_col, k=4)))

# order genes in cluster
geneorder<-data.frame(id=rownames(t(df[,c(27:ncol(df))])[p1$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
geneorder$order<-as.numeric(geneorder$order)

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneannot,by="row.names",all.x=T)
row.names(clusters)<-clusters$Row.names
clusters$Row.names<-NULL
geneorder<-merge(geneorder,clusters,by="id")


table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==8,]$Gene

patientscol<-merge(colcol,patients,by="row.names")
row.names(patientscol)<-patientscol$Row.names
patientscol$Row.names<-NULL
patientscol<-melt(patientscol[,c(1:17,ncol(patientscol))],id.vars=c("cluster"))

ggplot(patientscol,aes(x=as.factor(cluster),y=value))+
  geom_boxplot(alpha=0.4)+
  facet_wrap(~variable)+
theme_bw()


### DESeq DEGs in WT1 phenotypes
###### Differential EXPRESSION ANALYSIS ###
DESeq2Table <- estimateDispersions(DESeq2Table)
#plotDispEsts(DESeq2Table)


# Statistical testing of DE genes
levels(DESeq2Table)
design(DESeq2Table)

dds<-DESeq(DESeq2Table)
resultsNames(dds)
rescondition<-results(dds)

# Example plot counts
plotCounts(dds, gene=which(rownames(rescondition)=="ENSG00000184557"), intgroup="WT1mut")

head(rescondition)
summary(rescondition)
sum(rescondition$pvalue < 0.01,na.rm=TRUE)
sum(rescondition$padj < 0.1,na.rm=TRUE)
#plotMA(rescondition)

# Selecting specific DEGs
table_degs<-as.data.frame(rescondition)
table_degs<-merge(table_degs,bm,by.x=0,by.y="ensembl_gene_id",all.x=TRUE)
row.names(table_degs)<-table_degs$Row.names
table_degs$Row.names<-NULL

#up regulated genes
table_upreg_padj01<-subset(table_degs,table_degs$padj<0.1 & table_degs$log2FoldChange>0)
table_downreg_padj01<-subset(table_degs,table_degs$padj<0.1 & table_degs$log2FoldChange<0)


#write.table(table_degs,"/Users/yolanda_guillen/Desktop/IMIM/CEBPA/results/DEG_WTpheno.txt",quote = FALSE,sep="\t")
#write.table(table_upreg_padj01,"/Users/yolanda_guillen/Desktop/IMIM/CEBPA/results/upregulated_WTpheno.txt",quote = FALSE,sep="\t")
#write.table(table_downreg_padj01,"/Users/yolanda_guillen/Desktop/IMIM/CEBPA/results/downregulated_WTpheno.txt",quote = FALSE,sep="\t")


countab<-as.data.frame(counts(dds, normalized=TRUE))
countab<-t(countab)
countab<-melt(countab)
colnames(countab)<-c("Sample","Gene","value")
countab<-merge(countab,bm,by.x="Gene",by.y="ensembl_gene_id")

genesdown<-subset(table_degs,table_degs$pvalue<=0.01 & table_degs$log2FoldChange<0)
genesup<-subset(table_degs,table_degs$pvalue<=0.01 & table_degs$log2FoldChange>0)

countab<-merge(metadata,countab,by="Sample")

ggplot(countab[countab$Gene %in% c(row.names(genesup[genesup$padj<0.01,])),],aes(x=as.factor(WT1mut),y=log(value)))+
  geom_point(color="black",size=5,alpha=0.5)+
  geom_point(aes(color=as.factor(WT1mut)),size=4,alpha=0.5)+
  geom_boxplot(alpha=0.2)+
  scale_shape_manual(values=c(15,16,16,17,18,19,19,20,20,0))+
  scale_color_brewer(palette = "Set3")+
  ylab("Normalized counts")+
  xlab('')+
  facet_wrap(~external_gene_name,scales = "free",ncol=5)+
  theme_classic()+
  theme(axis.text.x = element_text(size=8,angle=45,hjust=1))

## Heatmap with differentially expressed genes

matrix<-as.data.frame(t(assay(rld)))

matrix<-t(matrix[,colnames(matrix) %in% c(row.names(genesup),row.names(genesdown))])

matrix<-merge(matrix,table_degs,by="row.names")
row.names(matrix)<-matrix$Row.names
matrix$Row.names<-NULL
matrix$deg<-ifelse(matrix$log2FoldChange<0,
                   'down','up')


pheatmap(t(scale(t(matrix[,1:13]))),
         annotation_col = as.data.frame(df[,c(which( colnames(df)==c("ENSG00000184937")),3,4,5,6,7:20,25,26)]),
         annotation_row = matrix[,c(22,24)],
         show_colnames= T,
         show_rownames = F,
         fontsize = 6,
         cutree_cols = 2,
         cutree_rows = 2,
         color=cols)



# Enrichment genes

dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "InterPro_Domains_2019",
         "GO_Cellular_Component_2017",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2016")

enrichedC8 <- enrichr(as.character(unique(clusters[clusters$cluster==8,]$external_gene_name)), dbs)
C8enr <- enrichedC8[["GO_Biological_Process_2018"]]

write.table(clusters,"/Users/yolanda_guillen/Desktop/IMIM/CEBPA/results/clusters.txt",sep="\t",quote = FALSE,row.names = FALSE)

enrichedup <- enrichr(as.character(unique(table_degs[table_degs$pvalue<0.01 & table_degs$log2FoldChange>0,]$external_gene_name)),dbs)
enricheddown <- enrichr(as.character(unique(table_degs[table_degs$pvalue<0.01 & table_degs$log2FoldChange<0,]$external_gene_name)),dbs)

enrichup <- enrichedup[["GO_Biological_Process_2018"]]
enrichdown <- enricheddown[["GO_Biological_Process_2018"]]

enrichup$group<-"UP"
enrichdown$group<-"DOWN"

allGO<-rbind(enrichup,enrichdown)
bpsub<-subset(allGO,allGO$Adjusted.P.value<0.05)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=unique((bpsub$Term))[rev(order(bpsub$P.value))])


library(tidyr)

bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))

# for multiple groups
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$group, rev(abs(bpsub$Combined.Score)))), ]$Term))

bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$group, abs(bpsub$Combined.Score))), ]$Term))

# heatmap
ggplot(bpsub, aes(x=group, y=Term)) + 
  geom_tile(aes(fill = -(log(P.value))),colour = "black")+
  scale_fill_gradient2(low = "white",
                       high = "red") +
  #scale_fill_gradient(low = "blue", high = "red")+
  theme_dark()+
  theme(axis.text.x = element_text(size=10, hjust = 1,angle=45),
        axis.text.y = element_text(size=8, hjust = 1,face="italic"))+
  ggtitle("Biological Processes")


### GSEA gene lists enrichment

GSEA_curated_down<-read.delim("../results/GSEA_genesets_curated_hallmark_position_downregulated.tsv")
GSEA_curated_down$group<-"Down_WThigh"
GSEA_curated_down$Term<-row.names(GSEA_curated_down)


plot1<-ggplot(GSEA_curated_down, aes(x=group, y=reorder(Term,rev(FDR)))) + 
  geom_tile(aes(fill = -(log(FDR))),colour = "black")+
  #scale_fill_gradient2(low = "white",high = "red") +
  scale_fill_gradient(low = "white", high = "red")+
  theme_dark()+
  theme(axis.text.x = element_text(size=4, hjust = 1,angle=45),
        axis.text.y = element_text(size=7, hjust = 1,face="italic"))+
  ggtitle("Gene sets enrichment")

ggsave(plot1, file="../results/up_DEGs_WT_genesets.pdf", scale=2)



## WT1 coexpression all samples Sant Pau
genetab<-as.data.frame(counts(DESeq2Table, normalized=TRUE))
row.names(genetab)
row.names(bm)<-bm$ensembl_gene_id
genetab<-merge(bm,genetab,by="row.names")
row.names(genetab)<-genetab$Row.names
genetab<-genetab[,c(4,7:ncol(genetab))]

gene_name<-data.frame(row.names=row.names(genetab),
                      Gene=row.names(genetab),
                      Symbol=genetab$external_gene_name)

genetab<-as.data.frame(t(genetab[,-1]))
row.names(genetab)

genetab_Dx<-genetab[row.names(genetab) %in% row.names(metadata[metadata$exp=="Dx_CEBPA",]),]

## WT1 coexpression public samples expression and clinical metadata
#ohsu<-read.delim("/Users/yolanda_guillen/Desktop/IMIM/CEBPA/AML_public_data/aml_ohsu_2018/data_RNA_Seq_expression_median.txt",sep="\t",header = TRUE)
ohsu<-read.delim("/Volumes/cancer/CEBPA/cbioportal/aml_ohsu_2018/data_RNA_Seq_expression_median.txt",sep="\t",header = TRUE)
ohsu$Entrez_Gene_Id<-NULL
ohsu<-ohsu[!duplicated(ohsu$Hugo_Symbol),]
row.names(ohsu)<-ohsu$Hugo_Symbol
ohsu<-as.data.frame(t(ohsu[,-1]))
row.names(ohsu)<-gsub('\\.','-',row.names(ohsu))

#ohsu_meta<-read.delim("/Users/yolanda_guillen/Desktop/IMIM/CEBPA/AML_public_data/aml_ohsu_2018/data_clinical_sample.txt",sep="\t",header = TRUE)
ohsu_meta<-read.delim("/Volumes/cancer/CEBPA/cbioportal/aml_ohsu_2018/data_clinical_sample.txt",sep="\t",header = TRUE)
ohsu_meta<-ohsu_meta[-c(1:4),]
row.names(ohsu_meta)<-ohsu_meta$Sample.Identifier

#ohsu_mut<-read.delim("/Users/yolanda_guillen/Desktop/IMIM/CEBPA/AML_public_data/aml_ohsu_2018/data_mutations_extended.txt",sep="\t",header = TRUE)
ohsu_mut<-read.delim("/Volumes/cancer/CEBPA/cbioportal/aml_ohsu_2018/data_mutations_extended.txt",sep="\t",header = TRUE)

# Change ohsu or genetab datasets!!!
corout<-apply(ohsu, 2, cor.test,ohsu$WT1,method="spearman")

corout<-apply(genetab_Dx, 2, cor.test,genetab_Dx$ENSG00000184937,method="spearman")

corout1<-as.data.frame(sapply(corout, "[[", "p.value"))
corout2<-as.data.frame(sapply(corout, "[[", "estimate"))
corout<-cbind(corout1,corout2)
names(corout)<-c("Pval","Rho")
corout$Gene<-row.names(corout)
#FOR GENETAB DATASET
corout<-merge(gene_name,corout,by="Gene",all.y=TRUE)

corout_sig_ohsu<-subset(corout,corout$Pval<0.000000001)
corout_sig_ohsu<-corout_sig_ohsu[with(corout_sig_ohsu,order(corout_sig_ohsu$Pval, corout_sig_ohsu$Rho)),]

corout_sig_Dx<-subset(corout,corout$Pval<0.01)
corout_sig_Dx<-corout_sig_Dx[with(corout_sig_Dx,order(corout_sig_Dx$Pval, corout_sig_Dx$Rho)),]


# plots
# CEBPA Dx
genecor<-merge(metadata,genetab,by=0)
row.names(genecor)<-genecor$Row.names
genecor$Row.names<-NULL

genecor<-genecor[,colnames(genecor) %in% corout_sig_Dx$Gene | colnames(genecor) %in% colnames(metadata)]

genecor_melt<-melt(genecor,id.vars = c(colnames(metadata),"ENSG00000184937"))
genecor_melt<-merge(genecor_melt,gene_name,by.x="variable",by.y="Gene",all.x=TRUE)

ggplot(genecor_melt[genecor_melt$variable %in% corout_sig_Dx$Gene[1:10],],aes(x=ENSG00000184937,y=value))+
  geom_point(size=5)+
  geom_point(data=genecor_melt[genecor_melt$variable %in% corout_sig_Dx$Gene[1:10],],aes(color=Pheno),size=7)+
  scale_color_brewer("Spectral")+
#  new_scale_color()+
#  geom_point(aes(color=Type),size=3)+
#  scale_color_brewer(palette="Paired")+
  geom_smooth(data=genecor_melt[genecor_melt$variable %in% corout_sig_Dx$Gene[1:10],],method="lm",se = FALSE,linetype="dotted")+
  #geom_label_repel(data=DSmat[DSmat$GSEA=="line" ,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
 facet_wrap(~Symbol+exp,ncol=4,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = "bottom")


pheatmap(t(scale(genecor[genecor$exp=="WT1_Rem",28:(ncol(genecor)-1)])),
         annotation_col = genecor[genecor$exp=="WT1_Rem",c(which( colnames(genecor)==c("ENSG00000184937")),3,22,24,25)],
#         annotation_row = matrix[,c(22,24)],
         show_colnames= F,
         show_rownames = F,
         fontsize = 6,
         cutree_cols = 2,
         cutree_rows = 2,
         color=cols)


ggplot(genecor,aes(x=exp,y=ENSG00000184937)) +
  geom_point(aes(color=Pheno),size=8) +
  geom_boxplot(alpha=0)+
  #stat_ellipse(aes(color=condition),level=0.8)+
  scale_color_brewer("Spectral")+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  ylab("WT1 expression") +
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1))


# Ohsu
ohsu_data<-merge(ohsu_meta,ohsu,by=0)
row.names(ohsu_data)<-ohsu_data$Row.names
ohsu_data$Row.names<-NULL

ohsu_data<-ohsu_data[,colnames(ohsu_data) %in% corout_sig_ohsu$Gene | colnames(ohsu_data) %in% colnames(ohsu_meta)]
ohsu_data_melt<-melt(ohsu_data,id.vars = c(colnames(ohsu_meta),"WT1"))

ggplot(ohsu_data_melt[ohsu_data_melt$variable %in% corout_sig_ohsu$Gene[1:10],],aes(x=WT1,y=value))+
  geom_point(size=1)+
  geom_point(data=ohsu_data_melt[ohsu_data_melt$variable %in% corout_sig_ohsu$Gene[1:10],],aes(color=Group),size=1)+
  scale_color_manual(values=coli)+
  #  new_scale_color()+
  #  geom_point(aes(color=Type),size=3)+
  #  scale_color_brewer(palette="Paired")+
#  geom_smooth(data=genecor_melt[genecor_melt$variable %in% corout_sig_Dx$Gene[1:10],],method="lm",se = FALSE,linetype="dotted")+
  #geom_label_repel(data=DSmat[DSmat$GSEA=="line" ,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  facet_wrap(~variable,ncol=4,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = "bottom")


pheatmap(t(scale(ohsu)),
 #        annotation_col = genecor[genecor$exp=="WT1_Rem",c(which( colnames(genecor)==c("ENSG00000184937")),3,22,24,25)],
         #         annotation_row = matrix[,c(22,24)],
         show_colnames= F,
         show_rownames = F,
         fontsize = 6,
         cutree_cols = 2,
         cutree_rows = 2,
         color=cols)


ggplot(genecor,aes(x=exp,y=ENSG00000184937)) +
  geom_point(aes(color=Pheno),size=8) +
  geom_boxplot(alpha=0)+
  #stat_ellipse(aes(color=condition),level=0.8)+
  scale_color_brewer("Spectral")+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  ylab("WT1 expression") +
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1))
