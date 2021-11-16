
## RScript to analyze RNAseq data using two options 1) HTSEQ and 2) Cufflinks cuffdiff

library(DESeq2)
library(biomaRt)
library(geneplotter)
library(MDplot)
library(lattice)
library(genefilter)

#Set wd
setwd("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/HTSEQ/")

output.dir="/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/HTSEQ/"

# Import metadata files
#metadata does not include bcat KO
metadata<-read.delim("metadata.txt",header = FALSE)
colnames(metadata)<-c("sampleID","countFile","condition","population")

metadata_pos<-subset(metadata,metadata$population=="pos")
metadata_neg<-subset(metadata,metadata$population=="neg")

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata_pos,
                                          directory = output.dir,
                                          design = ~ condition)

rowData(DESeq2Table)
mcols(rowData(DESeq2Table))

# Add additional annotation information using biomaRt
ensembl76 <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name"),
            values= rownames(DESeq2Table),
            filter="mgi_symbol",
            mart=ensembl76) 
head(bm)

#Add description data to gene counts
DESeq2Features <- data.frame(external_gene_name = rownames(DESeq2Table))
DESeq2Features$external_gene_name <- as.character(DESeq2Features$external_gene_name)

### join them together
rowData <- merge(DESeq2Features, bm, by = "external_gene_name")
rowData <- as(rowData, "DataFrame")

### add the annotation to the DESeq2 table
mcols(rowData(DESeq2Table)) <- c(mcols(rowData(DESeq2Table)),rowData)

DESeq2Table
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

### make sure to get fold change WT-deletion
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT", "KO"))

#### estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

# plot densities of counts for the different samples to assess their distributions
multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))

multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))


### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Table, blind=TRUE) ## create a distance matrix between the samples

## Create PCA
DESeq2::plotPCA(rld, intgroup=c("condition","population"))
dev.off()

## Create heatmap
pdf("HeatmapPlots.pdf")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

dev.off()

## Create PCA 
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

#pcaData$type<-paste(pcaData$condition,pcaData$population,sep="_")
ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=7,color="black",alpha=0.5)+
  geom_point(size=5,aes(color=condition)) +
  #geom_text(aes(label=name),hjust=1, vjust=0)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("darkred","darkolivegreen4"))+
  theme_bw()+
  coord_fixed()

dev.off()

### Differential EXPRESSION ANALYSIS ###
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

# Statistical testing of DE genes
levels(DESeq2Table$condition)
#levels(DESeq2Table$population)

design(DESeq2Table)
#design(DESeq2Table) <- formula(~ population + condition)

dds<-DESeq(DESeq2Table)
res<-results(dds,alpha = 0.1)
head(res)

#output results for GSEA
gseatable<-(as.data.frame(counts(dds),normalized=T))
row.names(gseatable)<-toupper(row.names(gseatable))
write.table(gseatable, file = "/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/GSEA/WT_vs_KO_EPhb2neg_normalized_counts.txt", sep = '\t', quote = FALSE)

resultsNames(dds)

resOdered <- res[order(res$pvalue),]
summary(res)

#genes with < adj 0.1
sum(res$padj < 0.1, na.rm=TRUE)

table(res[,"pvalue"] < 0.05)
plotMA(res)

#Heatmap most differentially expressed genes
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 500 )

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", dendrogram="column",trace = "none", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( WT="green", KO="red" )[
             colData(rld)$condition ],cexCol=0.3)

#To dataframe
res_df<-as.data.frame(res)
#preranked file for GSEA
res_df$SYMBOL<-row.names(res_df)
res_df<- res_df[order(res_df$log2FoldChange,res_df$pvalue),]

rank<-subset(res_df,select=c("SYMBOL","log2FoldChange"))
rank<- rank[order(rank$log2FoldChange),]
write.table(rank, file = "/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/GSEA/WT_vs_KO_EPhb2pos.rnk", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)


res_df_sig<-subset(res_df,res_df$padj<=0.1)
res_df_sig_IKB <- res_df_sig[order(res_df_sig$pvalue),]

#merge stem cell signature and DEGs
DEGs_IKB<-row.names(res_df_sig)

#Signature from ChIPhist.R
SCI<-scsignature$Gene
SCI<-unique(SCIc)

#Batlle signature
BatSCI<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/ISC_Batlle.txt",header = FALSE)
BatSCI<-paste0(toupper(substr(BatSCI$V1, 1, 1)), tolower(substr(BatSCI$V1, 2, nchar(BatSCI$V1))))
BatSCI<-unique(BatSCI)

BatLgr5<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/ISC_Lgr5_Batlle.txt",header = FALSE)
BatLgr5<-paste0(toupper(substr(BatLgr5$V1, 1, 1)), tolower(substr(BatLgr5$V1, 2, nchar(BatLgr5$V1))))
BatLgr5<-unique(BatLgr5)

BatPro<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/ISC_Pro_Batlle.txt",header = FALSE)
BatPro<-paste0(toupper(substr(BatPro$V1, 1, 1)), tolower(substr(BatPro$V1, 2, nchar(BatPro$V1))))
BatPro<-unique(BatPro)

BatLTA<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/ISC_LateTA_Batlle.txt",header = FALSE)
BatLTA<-paste0(toupper(substr(BatLTA$V1, 1, 1)), tolower(substr(BatLTA$V1, 2, nchar(BatLTA$V1))))
BatLTA<-unique(BatLTA)

### venn diagram ##
geneLists <- list(DEGs = DEGs, Bat=SCI)
geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off()
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue"), alpha=c(0.3,0.3), cex = 2, cat.fontface=4, category.names=c("DEGs", "SCI"), main="Clevers intestinal SC signature")

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

# To get the list of gene present in each Venn compartment we can use the gplots package
require("gplots")

a <- venn(VENN.LIST, show.plot=FALSE)

# You can inspect the contents of this object with the str() function
str(a)

# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a,"intersections")
length(inters[[1]])
length(inters[[2]])
length(inters[[3]])

common<-inters[[3]]
write.table(common,'/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/Clevers_ISC_DEGs.txt',quote = FALSE,row.names = FALSE)

D_KO<-subset(res_df_sig,res_df_sig$log2FoldChange<0)
D_KO_genes<-rownames(D_KO)
U_KO<-subset(res_df_sig,res_df_sig$log2FoldChange>0)
U_KO_genes<-rownames(U_KO)

D_KO_mice<-D_KO_genes
U_KO_mice<-U_KO_genes

#Manual heatmap from log trasformed data, using all samples
rld_df<-as.data.frame(assay(rld))

up_rld<-rld_df[with(rld_df, row.names(rld_df) %in% U_KO_genes),]
up_rld$gene<-row.names(up_rld)
up_rld<-melt(up_rld,id.vars =c("gene"))
up_rld$state<-c("UP_KO")

down_rld<-rld_df[with(rld_df, row.names(rld_df) %in% D_KO_genes),]
down_rld$gene<-row.names(down_rld)
down_rld<-melt(down_rld,id.vars =c("gene"))
down_rld$state<-c("DOWN_KO")

updown_rld<-rbind(up_rld,down_rld)

updown_rld$condition<-updown_rld$variable
updown_rld$condition<-gsub('KO.*','KO',updown_rld$condition)
updown_rld$condition<-gsub('WT.*','WT',updown_rld$condition)

updown_rld$value<-scale(updown_rld$value)

updown_rld$gene<-as.factor(updown_rld$gene)
updown_rld<-updown_rld[order(updown_rld$state), ]

#plot
updown_rld$gene <-reorder(updown_rld$gene, rev(updown_rld$value))

ggplot(updown_rld,aes(x=variable,y=gene,fill=value))+
  geom_tile(color="gray")+
  scale_fill_gradient2(high="red", low="blue",mid="white", midpoint=(-0.1), limits=range(updown_rld$value),na.value="white")+
  facet_wrap(~state,scales="free")+
  ggtitle("DEG p-value<0.05")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(size=5))




### Crosslink RNASeq and ChIP peaks 3me Histone ###
res_df_sig_IKB

#From ChIP unique genes peaks in KO and WT
Kpeaks<-unique(uniquegenes_KH3_sites$SYMBOL)
Wpeaks<-unique(uniquegenes_WTH3_sites$SYMBOL)


cross_rna_chip_KO<-res_df[with(res_df,res_df$SYMBOL %in% Kpeaks),]
cross_rna_chip_KO$condition<-"KO"


cross_rna_chip_WT<-res_df[with(res_df,res_df$SYMBOL %in% Wpeaks),]
cross_rna_chip_WT$condition<-"WT"


#Common wt and ko
cross_rna_chip_common<-res_df[with(res_df,res_df$SYMBOL %in% common_KO_WT),]
cross_rna_chip_common$condition<-"common"

cross_rna_chip_KO<-subset(cross_rna_chip_KO,select=c(log2FoldChange,padj,SYMBOL,condition))
cross_rna_chip_WT<-subset(cross_rna_chip_WT,select=c(log2FoldChange,padj,SYMBOL,condition))
cross_rna_chip_common<-subset(cross_rna_chip_common,select=c(log2FoldChange,padj,SYMBOL,condition))

cross_rna_chips<-rbind(cross_rna_chip_KO,cross_rna_chip_WT,cross_rna_chip_common)
cross_rna_chips$sig<-ifelse(cross_rna_chips$padj<0.1,c("sig"),c("nosig"))

cross_rna_chips$condition <- factor(cross_rna_chips$condition, levels = c("KO","common","WT"))

#ISC signature
scsignature<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/stem_cell_signature_list.txt")
scsignature$source <- "ISC"
scsignature$SYMBOL <- scsignature$Gene
scsignature$Gene<-NULL

cross_rna_chips_ISC<-merge(scsignature,cross_rna_chips,by="SYMBOL",all=TRUE)
cross_rna_chips_ISC$source[is.na(cross_rna_chips_ISC$source)] <- "noISC"
cross_rna_chips_ISC$var<-paste(cross_rna_chips_ISC$source,cross_rna_chips_ISC$sig)
cross_rna_chips_ISC<-subset(cross_rna_chips_ISC,cross_rna_chips_ISC$var!="ISC NA")

ggplot(cross_rna_chips_ISC,aes(x=condition,y=log2FoldChange,label=SYMBOL))+
  geom_jitter(aes(color=var))+
  geom_violin(data = subset(cross_rna_chips_ISC,cross_rna_chips_ISC$sig=="sig"),alpha=0.2)+
  geom_label_repel(data = subset(cross_rna_chips_ISC,cross_rna_chips_ISC$var=="ISC sig"),colour = "black", fontface = "italic",size=4)+
  #geom_violin(alpha=0.2)+
  scale_color_manual(values=c("grey","blue","grey","red"))+
  ylab("log2FC KO vs WT")+
  theme_bw()



#Batlle signature
BatSCI<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/ISC_Batlle.txt",header = FALSE)
BatSCI$SYMBOL<-paste0(toupper(substr(BatSCI$V1, 1, 1)), tolower(substr(BatSCI$V1, 2, nchar(BatSCI$V1))))
BatSCI$V1<-NULL
BatSCI$source<-"ISC_Batlle"

cross_rna_chips_ISC<-merge(BatSCI,cross_rna_chips,by="SYMBOL",all=TRUE)
cross_rna_chips_ISC$source[is.na(cross_rna_chips_ISC$source)] <- "noISC"
cross_rna_chips_ISC$var<-paste(cross_rna_chips_ISC$source,cross_rna_chips_ISC$sig)
cross_rna_chips_ISC<-subset(cross_rna_chips_ISC,cross_rna_chips_ISC$var!="ISC_Batlle NA")


ggplot(cross_rna_chips_ISC,aes(x=condition,y=log2FoldChange,label=SYMBOL))+
  geom_jitter(aes(color=var))+
  geom_violin(data = subset(cross_rna_chips_ISC,cross_rna_chips_ISC$sig=="sig"),alpha=0.2)+
  #geom_label_repel(data = subset(cross_rna_chips_ISC,cross_rna_chips_ISC$var=="ISC_Batlle sig"),colour = "black", fontface = "italic",size=2.5)+
  #geom_violin(alpha=0.2)+
  scale_color_manual(values=c("grey","blue","grey","red"))+
  ylab("log2FC KO vs WT")+
  theme_bw()


#Fetal signature (Helminths)
UpGAC<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/UP_GAC.txt",header = FALSE)
UpGAC$SYMBOL<-UpGAC$V1
UpGAC$V1<-NULL
UpGAC$source<-"UpGAC"

Down_GAC<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/DOWN_GAC.txt",header = FALSE)
Down_GAC$SYMBOL<-Down_GAC$V1
Down_GAC$V1<-NULL
Down_GAC$source<-"DownGAC"


GAC<-rbind(UpGAC,Down_GAC)

cross_rna_chips_ISC<-merge(GAC,cross_rna_chips,by="SYMBOL",all=TRUE)
cross_rna_chips_ISC$source[is.na(cross_rna_chips_ISC$source)] <- "noISC"
cross_rna_chips_ISC$var<-paste(cross_rna_chips_ISC$source,cross_rna_chips_ISC$sig)
cross_rna_chips_ISC<-subset(cross_rna_chips_ISC,cross_rna_chips_ISC$var!="UpGAC NA")
cross_rna_chips_ISC<-subset(cross_rna_chips_ISC,cross_rna_chips_ISC$var!="DownGAC NA")

cross_rna_chips_ISC$points<-ifelse(cross_rna_chips_ISC$var=="DownGAC sig",c("4"),c("2"))
cross_rna_chips_ISC$points<-ifelse(cross_rna_chips_ISC$var=="UpGAC sig",c("4"),c("2"))

ggplot(cross_rna_chips_ISC,aes(x=condition,y=log2FoldChange,label=SYMBOL))+
  geom_jitter(aes(color=var))+
  geom_violin(data = subset(cross_rna_chips_ISC,cross_rna_chips_ISC$sig=="sig"),alpha=0.0)+
  #geom_label_repel(data = subset(cross_rna_chips_ISC,cross_rna_chips_ISC$var=="UpGAC sig" | cross_rna_chips_ISC$var=="DownGAC sig"),colour = "black", fontface = "italic",size=2.5)+
  #geom_violin(alpha=0.2)+
  scale_color_manual(values=c("grey","purple","grey","chocolate2","grey","cyan4"))+
  ylab("log2FC KO vs WT")+
  theme_bw()




sigDEGs_all<-subset(cross_rna_chips,cross_rna_chips$sig=="sig")
kruskal.test(sigDEGs_all$condition,sigDEGs_all$log2FoldChange)

