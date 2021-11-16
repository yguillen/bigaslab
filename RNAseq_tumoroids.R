
## RScript to analyze RNAseq data using two options 1) HTSEQ and 2) Cufflinks cuffdiff

library(DESeq2)
library(biomaRt)
library(geneplotter)
library(MDplot)
library(lattice)
library(genefilter)


# signatures Batlle
BatSCI<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/ISC_Batlle.txt",header = FALSE)
BatSCI<-BatSCI$V1

#Fetal signature
fetalsignature_up<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/gene_signature_fetal_UP.txt")
fetalsignature_up$dir<-c("up")
fetalsignature_down<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/gene_signature_fetal_DOWN.txt")
fetalsignature_down$dir<-c("down")

fetalsignature_up<-toupper(fetalsignature_up$Gene)
fetalsignature_down<-toupper(fetalsignature_down$Gene)

#Set wd
setwd("/Volumes/grcmc/YGUILLEN/LSole_tumoroids/HTSEQ/")

output.dir="/Volumes/grcmc/YGUILLEN/LSole_tumoroids/HTSEQ/"

# Import metadata files

metadata<-read.delim("metadata.txt",header = FALSE)
colnames(metadata)<-c("sampleID","countFile","condition")

# For NFKBIA exons
metadata_exon<-read.delim("exon_metadata.txt",header = FALSE)
colnames(metadata_exon)<-c("sampleID","countFile","condition")


DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata,
                                          directory = output.dir,
                                          design = ~ condition)

rowData(DESeq2Table)
mcols(rowData(DESeq2Table))

# Add additional annotation information using biomaRt
ensembl78 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name"),
            values= rownames(DESeq2Table),
            #filter="hgnc_symbol",
            mart=ensembl78) 
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
DESeq2::plotPCA(rld, intgroup=c("condition"))
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
write.table(as.data.frame(counts(dds),normalized=T), file = "/Volumes/grcmc/YGUILLEN/LSole_tumoroids/GSEA/WT_vs_KO_EPhb2pos_normalized_counts.txt", sep = '\t', quote = FALSE)

resultsNames(dds)

resOdered <- res[order(res$pvalue),]
summary(res)

#genes with < 0.05
sum(res$pvalue < 0.05, na.rm=TRUE)

table(res[,"pvalue"] < 0.05)
plotMA(res)

#Heatmap most differentially expressed genes
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", dendrogram="column",trace = "none", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( WT="green", KO="red" )[
             colData(rld)$condition ],cexCol=0.3)

#To dataframe
res_df<-as.data.frame(res)
#preranked file for GSEA
res_df$SYMBOL<-row.names(res_df)
rank<-subset(res_df,select=c("SYMBOL","log2FoldChange"))
rank<- rank[order(rank$log2FoldChange),]
write.table(rank, file = "/Volumes/grcmc/YGUILLEN/LMarrueChIP/RNASeq/GSEA/WT_vs_KO_EPhb2pos.rnk", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)


res_df_sig<-subset(res_df,res_df$pvalue<=0.05)
res_df_sig <- res_df_sig[order(res_df_sig$pvalue),]

write.table(res_df,file="/Volumes/grcmc/YGUILLEN/LSole_tumoroids/DESEq2_results.csv",sep='\t',row.names = TRUE,quote=FALSE)

## GO enrichment

D_KO<-subset(res_df_sig,res_df_sig$log2FoldChange<0)
D_KO_genes<-rownames(D_KO)
U_KO<-subset(res_df_sig,res_df_sig$log2FoldChange>0)
U_KO_genes<-rownames(U_KO)


## GO terms and pathways
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018","TF_Perturbations_Followed_by_Expression")

enrichedD_KO <- enrichr(D_KO_genes, dbs)
enrichedU_KO <- enrichr(U_KO_genes,dbs)

D_KO_enrich <- enrichedD_KO[["GO_Biological_Process_2018"]]
D_KO_enrich$group<-c("Down_KO")

U_KO_enrich <- enrichedU_KO[["GO_Biological_Process_2018"]]
U_KO_enrich$group<-c("Up_KO")

allGO<-rbind(U_KO_enrich ,D_KO_enrich)

bpsub<-subset(allGO,allGO$Adjusted.P.value<0.1)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$Adjusted.P.value))])

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))

ggplot(bpsub,aes(x=group,y=Term,fill=P.value),alpha = 0.5)+
  geom_tile(color="black")+
  geom_text(aes(label=tolower(Genes)), size=2,hjust=0.5,color="gray")+
  scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  facet_wrap(~group,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))

write.table(bpsub,file="/Volumes/grcmc/YGUILLEN/LSole_tumoroids/DESEq2/GO_BP_pval005_enrichpadj01_results.csv",sep='\t',row.names = TRUE,quote=FALSE)


## Compare UP KO tumoroids vs UP KO intestinal mice cells (Ephb2 negative)
#From mice
D_KO_mice
U_KO_mice
#From tumoroids, transform capital letters
D_KO_genes<-paste0(toupper(substr(D_KO_genes, 1, 1)), tolower(substr(D_KO_genes, 2, nchar(D_KO_genes))))
U_KO_genes<-paste0(toupper(substr(U_KO_genes, 1, 1)), tolower(substr(U_KO_genes, 2, nchar(U_KO_genes))))


### venn diagram ##
#down
#geneLists <- list(D_KO_mice = D_KO_mice, D_KO_tumor = D_KO_genes)
#up
#geneLists <- list(U_KO_mice = U_KO_mice, U_KO_tumor = U_KO_genes)

#Stem cell signature
geneLists <- list(ISC = fetalsignature_down, U_KO_tumor = U_KO_genes, D_KO_tumor = D_KO_genes)

geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off()
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("purple","darkolivegreen", "darkblue"), alpha=c(0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("ISC", "Up KO tumoroids","Down KO tumoroids"), main="Genes")

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
length(inters[[4]])

