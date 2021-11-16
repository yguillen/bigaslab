
## RScript to analyze RNAseq data using two options 1) HTSEQ and 2) Cufflinks cuffdiff
install.packages("plotly")

library(DESeq2)
library(biomaRt)
library(geneplotter)
library(MDplot)
library(lattice)
library(genefilter)
library(plotly)
library(limma)
library(ggrepel)

#Set wd
setwd("/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/HTSEQ/")

output.dir="/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/HTSEQ"


#Sarah Bonin selection based on old PCA
#B1_DEGs<-as.data.frame(read_xlsx("/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/B1_SBonin/2018-04-20_CRuiz_RNAseq_olddata_selection1_RNAseq_DESeq2.xlsx"))

#Polycomb target genes
EED<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/gene_lists/EED_benporath.txt",header=FALSE)
EED<-EED[-c(1,2),]
EED<-as.data.frame(EED)
EED$EED<-paste0(toupper(substr(as.character(EED$EED), 1, 1)), tolower(substr(as.character(EED$EED), 2, nchar(as.character(EED$EED)))))

#Suz12 target genes
Suz12<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/gene_lists/Suz12_Benporath.txt",header=FALSE)
Suz12<-Suz12[-c(1,2),]
Suz12<-as.data.frame(Suz12)
Suz12$Suz12<-paste0(toupper(substr(as.character(Suz12$Suz12), 1, 1)), tolower(substr(as.character(Suz12$Suz12), 2, nchar(as.character(Suz12$Suz12)))))

#M3H target genes
H3M<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/gene_lists/H3K27ME3_benporath.txt",header=FALSE)
H3M<-H3M[-c(1,2),]
H3M<-as.data.frame(H3M)
H3M$H3M<-paste0(toupper(substr(as.character(H3M$H3M), 1, 1)), tolower(substr(as.character(H3M$H3M), 2, nchar(as.character(H3M$H3M)))))



# Import metadata files
#metadata does not include bcat KO, includes B3 but not B1 and not Trumpp
metadata<-read.delim("metadata_Trumpp.txt",header = FALSE)

metadata_A<-read.delim("metadata_A.txt",header=FALSE)

colnames(metadata)<-c("sampleID","countFile","batch","condition","population","gender")

colnames(metadata_A)<-c("sampleID","countFile","batch","condition","population","gender")

#excluding B3_KO_4, possible false KO, and two other discarded samples B3_KO_6 based on dendrogram all samples and B1_14002
#metadata<-subset(metadata,metadata$sampleID!="B3_KO_4_18821_GTCCGC")
#metadata<-subset(metadata,metadata$sampleID!="B3_KO_6_18823_ATCACG")
#metadata<-subset(metadata,metadata$sampleID!="B1_14002")

#metadata$group <-factor(paste0(metadata$population, metadata$condition))
metadata_A$group <-factor(paste0(metadata_A$population, metadata_A$condition))

#metadata_bcat includes bcat KO samples
metadat_bcat<-read.delim("metadata_bcat.txt",header = FALSE)
colnames(metadat_bcat)<-c("sampleID","countFile","batch","condition","population")
metadata_bcat<-metadat_bcat[c(3:5,8:10),]


metadata_NP<-subset(metadata,metadata$population=="NP")
metadata_NP$group <- factor(paste0(metadata_NP$gender, metadata_NP$condition))

metadata_NP_noTrumpp<-subset(metadata_NP,metadata_NP$batch!="NPTR")

#metadata_A<-subset(metadata,metadata$population=="A")
#metadata_A$group <- factor(paste0(metadata_A$gender, metadata_A$condition))

class(metadata_A$condition)

metadata_A_B1<-subset(metadata_A,metadata_A$batch=="B1")
metadata_A_B2<-subset(metadata_A,metadata_A$batch=="B2")
metadata_A_B3<-subset(metadata_A,metadata_A$batch=="B3")
metadata_A_B2B3<-subset(metadata_A,metadata_A$batch!="B1")


#variable group as merging gender and condition
metadata_A_WT<-subset(metadata_A,metadata_A$condition=="WT")
#only HSCs
metadata_adHSC<-subset(metadata,metadata$batch=="rawHSC")

#excluding normalized preprocessed data a and dHSCs
metadata_all_HSC<-subset(metadata,metadata$batch=="proHSC")

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata_NP_noTrumpp,
                                          directory = output.dir,
                                          design = ~ condition)

rowData(DESeq2Table)
mcols(rowData(DESeq2Table))

# Add additional annotation information using biomaRt
ensembl76 <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","mgi_symbol","gene_biotype","chromosome_name"),
            values= rownames(DESeq2Table),
#            filter="external_gene_name",
            mart=ensembl76) 

head(bm)

#Add description data to gene counts
DESeq2Features <- data.frame(id = rownames(DESeq2Table))
DESeq2Features$id <- as.character(DESeq2Features$id)


### join them together
#if htseq ensemblid
rowData <- merge(DESeq2Features, bm, by.x="id",by.y = "ensembl_gene_id")

#if htseqgene symbol
rowData <- merge(DESeq2Features, bm, by.x="id",by.y = "external_gene_name")

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
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","KO"))
colData(DESeq2Table)$group <- factor(colData(DESeq2Table)$group, levels = c("FKO","FWT","MKO","MWT"))
#colData(DESeq2Table)$group <- factor(colData(DESeq2Table)$group, levels = c("AWT", "AKO","NPWT","NPKO"))
colData(DESeq2Table)$gender <- factor(colData(DESeq2Table)$gender, levels = c("F", "M"))
#colData(DESeq2Table)$batch <- factor(colData(DESeq2Table)$batch, levels = c("NP","NPTR","B2", "B3","B1"))
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("aHSC", "dHSC"))

#### estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

# plot densities of counts for the different samples to assess their distributions
multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))


### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Table, blind=TRUE) ## create a distance matrix between the samples

## Create PCA
DESeq2::plotPCA(rld, intgroup=c("population"))
dev.off()

## Create heatmap
pdf("HeatmapPlots.pdf")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

dev.off()


## Create PCA 
pcaData <- plotPCA(rld, intgroup=c("condition","batch","population","gender"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$group<-gsub(':NP:.*','',pcaData$group)

ggplot(pcaData, aes(PC1, PC2, color=condition,label=name)) +
  geom_point(size=8) +
  stat_ellipse(aes(color=condition),level=0.8)+
  #scale_color_manual(values=c("dodgerblue","coral"))+
  #scale_color_manual(values=c("grey","dodgerblue","coral","darkolivegreen2","black","dodgerblue4","coral4"))+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #facet_grid(~gender)+
  scale_color_manual(values=c("dodgerblue2","red"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed()


ggplot(pcaData, aes(y=PC1, x=population,label=name)) +
  geom_point(aes(color=population),size=8) +
  geom_boxplot(alpha=0)+
  #stat_ellipse(aes(color=condition),level=0.8)+
  scale_color_manual(values=c("dodgerblue2","red"))+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
  xlab("")+
  #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #facet_grid(~gender)+
  #scale_color_manual(values=c("purple","red"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1))


## adonis test
library(vegan)

rld_df<-as.data.frame(assay(rld))

rownames(rld_df)
colnames(rld_df)
t_rld_df<-as.data.frame(t(rld_df))
t_rld_df$sampleID<-row.names(t_rld_df)
#t_rld_df<-merge(metadata_A,t_rld_df,by="sampleID")
#t_rld_df<-merge(metadata_NP_noTrumpp,t_rld_df,by="sampleID")
t_rld_df<-merge(metadata_all_HSC,t_rld_df,by="sampleID")

mat_gene<-t_rld_df[,8:ncol(t_rld_df)]

adonis(mat_gene ~ condition,method="euclidean", data = t_rld_df)

### Differential EXPRESSION ANALYSIS ###
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

# Statistical testing of DE genes
levels(DESeq2Table$condition)
levels(DESeq2Table$gender)
levels(DESeq2Table$group)
levels(DESeq2Table$batch)
levels(DESeq2Table$population)

design(DESeq2Table)

dds<-DESeq(DESeq2Table)
resultsNames(dds)


#for batch, 3 categories
ddsLRT <- DESeq(DESeq2Table,test="LRT",reduced=~1)
resultsNames(ddsLRT)


#output results for GSEA
#ranked file
rank<-subset(res_df_rna,select=c("log2FoldChange"))
rank$NAME<-row.names(rank)
rank<- rank[order(rank$log2FoldChange),]
rank<-subset(rank,select=c("NAME","log2FoldChange"))
rank$NAME<-toupper(rank$NAME)
write.table(rank, file = "/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/GSEA/genes.rnk", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#normalized table
gseatable<-(as.data.frame(counts(dds),normalized=T))
row.names(gseatable)<-toupper(row.names(gseatable))

write.table(gseatable, file = "/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/GSEA/KO_vs_WT_normalized_counts.txt", sep = '\t', quote = FALSE)


#results considering only males or females, DEGs KO vs WT
resmales<-results(dds, contrast=c("group", "MKO", "MWT"),alpha=0.1)
resfemales<-results(dds,contrast=c("group","FWT","FKO"),alpha=0.1)

#results DEGs between males and females
resgender<-results(dds,contrast=c("group","FWT","MWT"),alpha=0.1)

#results DEGs between condition
rescondition<-results(dds,contrast=c("condition","KO","WT"),alpha=0.1)

#for POPULATION! check population or condition
rescondition<-results(dds,contrast=c("condition","dHSC","aHSC"),alpha=0.1)

#rescondition<-results(dds)

#for batch
resbatch<-results(ddsLRT)
head(resbatch)

resbatch<-results(ddsLRT, name="batch_B1_vs_B2" , test="Wald")
summary(resbatch)
sum(resbatch$padj < 0.1,na.rm=TRUE)
plotMA(resbatch)


#Example differences genes between groups
plotCounts(dds, gene=which(rownames(rescondition)=="ENSMUSG00000000167"), intgroup="condition")

plotCounts(dds, gene=which(rownames(rescondition)=="Psmd11"), intgroup="condition")

head(resmales)
head(resfemales)
head(resgender)
head(rescondition)


summary(resgender)
summary(resmales)
summary(resfemales)
summary(rescondition)

#genes with < 0.05
sum(resgender$padj < 0.1, na.rm=TRUE)
sum(resmales$padj < 0.1, na.rm=TRUE)
sum(resfemales$padj < 0.1, na.rm=TRUE)
sum(rescondition$pvalue < 0.05,na.rm=TRUE)

plotMA(resgender)
plotMA(resfemales)
plotMA(resmales)
plotMA(rescondition)

#Heatmap most differentially expressed genes (only ranked by DEGs)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 300 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", dendrogram="column",trace = "none", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( WT="grey", KO="black" )[
             colData(rld)$condition ],cexCol=0.3)

# Selecting specific DEGs
#To dataframe (ckeck females, males or gender)
res_df_rna<-as.data.frame(rescondition)
class(rowData)
rowData<-as.data.frame(rowData)
res_df_rna<-merge(rowData,res_df_rna,by.x="id",by.y="row.names",all.y=TRUE)
res_df_sig_rna<-subset(res_df_rna,res_df_rna$pvalue<=0.05)
res_df_sig_rna <- res_df_sig_rna[order(res_df_sig_rna$pvalue),]


# UP and DOWN genes
UPgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange>0))$id)
DOWNgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange<0))$id)
length(UPgenes)
length(DOWNgenes)

#genderNP<-row.names(res_df_sig_rna)

#Manual heatmap from log trasformed data, using all samples
rld_df<-as.data.frame(assay(rld))

# Selecting only WT samples
WTsamples<-(subset(metadata_NP,metadata_NP$condition=="WT"))$sampleID

# Selecting only male samples
Malesamples<-(subset(metadata_A,metadata_A$gender=="M"))$sampleID

# Selecting only female samples
Femalesamples<-(subset(metadata_A,metadata_A$gender=="F"))$sampleID

wt_rld<-which(colnames(rld_df) %in% WTsamples)
wt_rld<-rld_df[,wt_rld]
colnames(wt_rld)<-c("F","M","M")
#colnames(wt_rld)<-c("KO","KO","WT")
#colnames(wt_rld)<-c("KO","KO","KO","WT","WT","WT")

#Selecting DEGs gender from wt_rld
deggenes<-row.names(res_df_rna)
deg_rld<-which(rownames(wt_rld) %in% deggenes)
deg_rld<-wt_rld[deg_rld,]

# for condition
#all genes
deggenes<-res_df_rna$id

#only degs
deggenes<-res_df_sig_rna$id


deg_rld<-which(rownames(rld_df) %in% deggenes)
deg_rld_df<-rld_df[deg_rld,]

deg_rld<-merge(deg_rld_df,res_df_rna,by.x="row.names",by.y="id")

deg_rld<- deg_rld[order(deg_rld$log2FoldChange),]
deg_rld$Row.names

heatmap.2( as.matrix(deg_rld[,2:7]), scale="row", dendrogram="column",trace="none", Rowv=FALSE,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           cexCol=0.3)

#for condition
heatmap.2( as.matrix(deg_rld[,2:9]), scale="row", dendrogram="column",trace="none", Rowv=FALSE,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cexCol=0.3)


#correlation plot gene expression

#B1B2B3 list
B1B2B3_DEGs<-deg_rld
B1B2B3_DEGs_comp<-subset(B1B2B3_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
colnames(B1B2B3_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
B1B2B3_DEGs_comp<-B1B2B3_DEGs_comp[!duplicated(B1B2B3_DEGs_comp$Row.names),]
B1B2B3_DEGs_comp$population<-"All"


#NP list
NP_DEGs<-deg_rld
NP_DEGs_comp<-subset(NP_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
colnames(NP_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
NP_DEGs_comp<-NP_DEGs_comp[!duplicated(NP_DEGs_comp$Row.names),]
NP_DEGs_comp$population<-"NP"


#aHSC and dHSC list
HSC_DEGs<-deg_rld
HSC_DEGs_comp<-subset(HSC_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
colnames(HSC_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
HSC_DEGs_comp<-HSC_DEGs_comp[!duplicated(HSC_DEGs_comp$mgi_symbol),]
HSC_DEGs_comp$population<-"HSCs"



write.table(toupper(subset(HSC_DEGs,HSC_DEGs$log2FoldChange<0)$external_gene_name),"/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/gene_lists/down_dormantHSCS.txt",quote = FALSE,row.names = FALSE)
write.table(toupper(subset(HSC_DEGs,HSC_DEGs$log2FoldChange>0)$external_gene_name),"/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/gene_lists/up_dormantHSCS.txt",quote = FALSE,row.names = FALSE)


#colnames(Ball_comp)<-c("Gene","log2FC_B1","pval_B1","padj_B1","log2FC_B2","pval_B2","padj_B2","log2FC_B3","pval_B3","padj_B3")

#Merge all batches
Ball_comp<-rbind(B1B2B3_DEGs_comp,NP_DEGs_comp,HSC_DEGs_comp)

#Select DEgs padj < 0.1 from here
B1B2B3_DEGs_sig<-as.character((subset(B1B2B3_DEGs_comp,B1B2B3_DEGs_comp$pvalue<=0.05))$Row.names)
NP_DEGs_sig<-(subset(NP_DEGs_comp,NP_DEGs_comp$pvalue<=0.05))$Row.names

HSC_DEGs_comp$Row.names<-HSC_DEGs_comp$mgi_symbol
HSC_DEGs_sig<-unique((subset(HSC_DEGs_comp,HSC_DEGs_comp$pvalue<=0.05 & (!is.na(HSC_DEGs_comp[,1]))))$Row.names)

alldegs<-unique(c(B1B2B3_DEGs_sig,NP_DEGs_sig,HSC_DEGs_sig))
length(alldegs)

Ball_comp<-Ball_comp[(Ball_comp$Row.names %in% alldegs),]
Ball_comp$state<-ifelse(Ball_comp$log2FoldChange< 0, 
                      c("KOdown"), c("KOup")) 


Ball_comp$state<-as.factor(Ball_comp$state)
Ball_comp<- Ball_comp[order(Ball_comp$population,Ball_comp$pvalue),]

Ball_comp$Row.names<-as.character(Ball_comp$Row.names)

ggplot(Ball_comp, aes(Row.names,population)) +
  geom_tile(aes(fill = state), color = "white") +
  geom_point(aes(size=(-log(pvalue))),data = subset(Ball_comp,Ball_comp$padj<0.1))+
  scale_fill_manual(values = c("darkolivegreen","coral3")) +
  xlab("Genes") +
  ylab("Batch") +
  theme_minimal()+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=8,face="bold"),
        axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=4, hjust = 1)) +
  labs(fill = "Expression level fold change")+
  #facet_wrap(~state)+
  coord_flip()



#Merge different batches or populations to melt
B1B2B3_ball_comp<-subset(Ball_comp,Ball_comp$population=="All")
NP_ball_comp<-subset(Ball_comp,Ball_comp$population=="NP")
HSCs_ball_comp<-subset(Ball_comp,Ball_comp$population=="HSCs")

HSCs_ball_comp<-subset(HSCs_ball_comp, (!is.na(HSCs_ball_comp[,1])))

B1B2B3_all_NP_comp<-merge(B1B2B3_ball_comp,NP_ball_comp,by="Row.names",all=TRUE)
#HSCs_ball_comp$Row.names<-as.character(HSCs_ball_comp$Row.names)
#HSCs_ball_comp$Row.names<-HSCs_ball_comp$Row.names
B1B2B3_all_NP_comp<-merge(B1B2B3_all_NP_comp,HSCs_ball_comp,by="Row.names",all=TRUE)

#change names
names(B1B2B3_all_NP_comp)[2]<-"log2FoldChange_all"
names(B1B2B3_all_NP_comp)[10]<-"log2FoldChange_NP"
names(B1B2B3_all_NP_comp)[18]<-"log2FoldChange_HSCs"

B1B2B3_all_NP_comp<- B1B2B3_all_NP_comp[order(B1B2B3_all_NP_comp$pvalue),]

B1B2B3_all_NP_comp$significativeall<-ifelse(B1B2B3_all_NP_comp$pvalue.x < 0.05, 
                               c("sig_all"), c("nosig_all")) 

B1B2B3_all_NP_comp$significativeNP<-ifelse(B1B2B3_all_NP_comp$pvalue.y < 0.05, 
                               c("sig_NP"), c("nosig_NP"))  

B1B2B3_all_NP_comp$significativeHSC<-ifelse(B1B2B3_all_NP_comp$pvalue < 0.05, 
                                           c("sig_HSC"), c("nosig_HSC"))  

B1B2B3_all_NP_comp$significative<-paste(B1B2B3_all_NP_comp$significativeall,B1B2B3_all_NP_comp$significativeNP,B1B2B3_all_NP_comp$significativeHSC)


#Add signatures polycomb
#EED$signature<-"EED"
#Suz12$signature<-"Suz12"
#H3M$signature<-"H3M"

#B12_comp<-merge(B12_comp,EED,by.x="Row.names",by.y="EED",all=TRUE)
#B12_comp<-merge(B12_comp,Suz12,by.x="Row.names",by.y="Suz12",all=TRUE)
#B12_comp<-merge(B12_comp,H3M,by.x="Row.names",by.y="H3M",all=TRUE)
#B12_comp$signature<-paste(B12_comp$signature,B12_comp$signature.x,B12_comp$signature.y,sep="_")

#signature MolOUP
MolOUP<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/gene_lists/MolOUP.txt",header=FALSE)
MolODOWN<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/gene_lists/MolODOWN.txt",header=FALSE)

MolOUP$sig<-"MolUP"
MolODOWN$sig<-"MolDOWN"

MolO<-rbind(MolODOWN,MolOUP)
names(MolO)[1]<-"Row.names"

#B1B2B3_all_NP_comp<-merge(B1B2B3_all_NP_comp,MolO,by="Row.names",all=TRUE)

B1B2B3_all_NP_comp$sigtype<-paste(B1B2B3_all_NP_comp$significativeall,B1B2B3_all_NP_comp$significativeHSC)

B1B2B3_all_NP_comp$sigtype<-ifelse(B1B2B3_all_NP_comp$sigtype == "sig_all sig_HSC", 
                                            c("sig_both"), c("sig_one"))  

#NDU genes genes
B1B2B3_all_NP_comp$names<-B1B2B3_all_NP_comp$Row.names
B1B2B3_all_NP_comp$names<-gsub('Ndu.*','NDU',B1B2B3_all_NP_comp$names)
B1B2B3_all_NP_comp$names<-ifelse(B1B2B3_all_NP_comp$names == "NDU", 
                                   c("NDU"), c("no_ndu"))  


#PSM genes genes
B1B2B3_all_NP_comp$names<-B1B2B3_all_NP_comp$Row.names
B1B2B3_all_NP_comp$names<-gsub('Psm.*','PSM',B1B2B3_all_NP_comp$names)
B1B2B3_all_NP_comp$names<-ifelse(B1B2B3_all_NP_comp$names == "PSM", 
                                 c("PSM"), c("no_PSM"))  


#OLD PLOT
ggplot(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$pvalue.y<=0.05 | B1B2B3_all_NP_comp$pvalue<=0.05),aes(x=log2FoldChange_NP,y=log2FoldChange_HSCs))+
  #annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0)+
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0)+
  geom_point(color="black",size=3)+
  geom_point(color="grey",size=1)+
  #geom_density_2d(color="black",bins=15)+
  geom_density_2d(data=(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$sigtype=="sig_both")),color="purple")+
  scale_color_manual(values=c("black","grey"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = Row.names),color="darkolivegreen", data = subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$names=="NDU" & B1B2B3_all_NP_comp$state.x=="KOdown" & B1B2B3_all_NP_comp$state=="KOdown"), fontface = "italic",size=2.5)+
  #geom_label_repel(aes(label = mgi_symbol.x), color="red",data = subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$significativeall=="sig_all" & B1B2B3_all_NP_comp$state.x=="KOdown" & B1B2B3_all_NP_comp$significativeHSC=="sig_HSC" & B1B2B3_all_NP_comp$state=="KOdown"), fontface = "italic",size=2.5)+
  xlim(-10,10)+
  ylim(-10,10)+
  xlab("log2FC CD48- CD150+ KO vs. WT")+
  ylab("log2FC dHSCs vs. aHSC")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1)) 


#new plot as fetal
ggplot(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$pvalue.y<=0.05 | B1B2B3_all_NP_comp$pvalue<=0.05),aes(x=log2FoldChange_NP,y=log2FoldChange_HSCs))+
  #annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "blue",alpha=0.1)+
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "red",alpha=0.1)+
  geom_point(color="black",size=3,alpha=0.5)+
  #geom_point(data=subset(merge_WNK,merge_WNK$sigfet=="sig_fet"),color="red",size=1)+
  #geom_point(data=subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$pvalue.x<=0.05 | B1B2B3_all_NP_comp$pvalue<=0.05 & B1B2B3_all_NP_comp$sigtype=="sig_both"),color="red",size=1)+
  geom_density_2d(color="green")+
  #geom_density_2d(data=(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$sigtype=="sig_both")),color="green")+
  #geom_density_2d(data=(subset(merge_WNK,merge_WNK$sigfet=="sig_fet")),color="red")+
  #geom_smooth(data=subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$pvalue.x<=0.05 | B1B2B3_all_NP_comp$pvalue<=0.05),method="lm",color="red")+
  #scale_color_manual(values=c("black","grey"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = Row.names),color="black", data = subset(merge_WNK,merge_WNK$sigtype=="sig_all"), fontface = "italic",size=5)+
  #geom_label_repel(aes(label = Row.names),color="purple", data = subset(B1B2B3_all_NP_comp,(B1B2B3_all_NP_comp$significativeall=="sig_all" | B1B2B3_all_NP_comp$significativeHSC=="sig_HSC") & B1B2B3_all_NP_comp$names=="PSM"), fontface = "italic",size=4)+
  xlim(-10,10)+
  ylim(-10,10)+
  xlab("log2FC KO vs. WT CD34- CD135-+")+
  #ylab("log2FC NICD vs. WT fetal")+
  ylab("log2FC dHSCs vs aHSC")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1)) 


cor.test((subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$pvalue.y<=0.05 | B1B2B3_all_NP_comp$pvalue<=0.05))$log2FoldChange_NP,(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$pvalue.y<=0.05 | B1B2B3_all_NP_comp$pvalue<=0.05))$log2FoldChange_HSCs)


Upconc<-(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$pvalue.y<=0.05 & 
                        B1B2B3_all_NP_comp$state.y=="KOup" & 
                        B1B2B3_all_NP_comp$pvalue<=0.05 & 
                        B1B2B3_all_NP_comp$state=="KOdown")$Row.names)
Downconc<-(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$pvalue.y<=0.05 & 
                    B1B2B3_all_NP_comp$state.x=="KOdown" & 
                    B1B2B3_all_NP_comp$pvalue<0.05 & 
                    B1B2B3_all_NP_comp$state=="KOup")$Row.names)

# UP and DOWN subset genes
Up_NP<-unique((subset(NP_DEGs_comp,NP_DEGs_comp$pvalue<0.05 & NP_DEGs_comp$log2FoldChange>0))$Row.names)
Down_NP<-unique((subset(NP_DEGs_comp,NP_DEGs_comp$pvalue<0.05 & NP_DEGs_comp$log2FoldChange<0))$Row.names)


#3D plots
p <- plot_ly(B12_comp, x = ~log2FoldChange_B3, y = ~log2FoldChange_B1, z = ~log2FoldChange_B2, color = ~significative) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'B3'),
                      yaxis = list(title = 'B1'),
                      zaxis = list(title = 'B2')))
p

B1B2B3_UP<-subset(B12_comp,B12_comp$state=="KOup KOup KOup")
B1B2B3_UP$state_all<-"UP"
B1B2B3_DOWN<-subset(B12_comp,B12_comp$state=="KOdown KOdown KOdown")
B1B2B3_DOWN$state_all<-"DOWN"

B1B2B3_UPDOWN<-rbind(B1B2B3_UP,B1B2B3_DOWN)

write.table(B12_comp,"/Volumes/grcmc/CRISTINA/Yolanda_RNAseq/B1_B2_padj01.csv", sep = "\t",quote = FALSE,row.names = FALSE)

#geneLists <- list(B2 = (subset(B2_DEGs,B2_DEGs$pvalue<0.05))$Row.names, B3 = (subset(B3_DEGs,B3_DEGs$pvalue<0.05))$Row.names, B1 = (subset(B1_DEGs,B1_DEGs$pvalue_KO_vs_WT<0.05))$gene_name)

geneLists <- list(B2 = (subset(B2_DEGs,B2_DEGs$padj<0.1))$Row.names,Suz12 = Suz12$Suz12, H3M = H3M$H3M, EED = EED$EED)

geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off()
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue","red","grey"), alpha=c(0.3,0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("B1", "Suz12","H3M","EED"), main="Genes")

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

require("gplots")

a <- venn(VENN.LIST, show.plot=FALSE)

# You can inspect the contents of this object with the str() function
str(a)

inters <- attr(a,"intersections")
length(inters[[1]])
length(inters[[2]])
length(inters[[3]])
length(inters[[4]])
length(inters[[5]])
length(inters[[6]])
length(inters[[7]])
length(inters[[8]])
length(inters[[9]])
length(inters[[10]])
length(inters[[11]])
length(inters[[12]])
length(inters[[13]])



commonB2B3<-inters[[3]]



B2B3gene_cor<-res_df_sig_rna
B2B3genes<-row.names(B2B3gene_cor)

NPgene_cor<-res_df_sig_rna
NPgenes<-row.names(NPgene_cor)

condgenes<-c(B2B3genes,NPgenes)


NPtab<-res_df_rna[(rownames(res_df_rna) %in% condgenes),]
ABtab<-res_df_rna[(rownames(res_df_rna) %in% condgenes),]

names(NPtab)[2]<-"log2FoldChange_NP"
names(ABtab)[2]<-"log2FoldChange_AB"

alltab<-merge(NPtab,ABtab,by="row.names",all=TRUE)

alltab<- alltab[order(alltab$pvalue.x),]

alltab$significativeNP<-ifelse(alltab$pvalue.x < 0.05, 
                    c("sigNP"), c("nosigNP")) 

alltab$significativeAB<-ifelse(alltab$pvalue.y < 0.05, 
                             c("sigAB"), c("nosigAB")) 
alltab$adjust<-ifelse(alltab$padj.x < 0.1, 
                               c("adj"), c("noadjust")) 

alltab$significative<-paste(alltab$significativeAB,alltab$significativeNP)

ggplot(alltab,aes(x=log2FoldChange_NP,y=log2FoldChange_AB,label=Row.names))+
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "darkolivegreen",alpha=0.2)+
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "dodgerblue",alpha=0.2)+
  geom_point(aes(color=significative))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_label_repel(aes(label = Row.names,color=adjust), data = subset(alltab,alltab$significative=="sigAB sigNP"), fontface = "italic",size=2.5)+
  #geom_label_repel(aes(label = Row.names,color=adjust), data = subset(alltab,alltab$log2FoldChange_NP>=5), fontface = "italic",size=2.5)+  
  theme_bw()

DEGsgenderNP<-deggenes
DEGgenderAB<-deggenes

DEGsmalesNPWTKO<-deggenes
DEGsfemalesNPWTKO<-deggenes

setdiff(DEGsgenderNP,DEGgenderAB)
setdiff(DEGgenderAB,DEGsgenderNP)


DEGsfemales<-deggenes
DEGsmales<-deggenes

## GO terms and pathways
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018","Chromosome_Location","GO_Molecular_Function_2018","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X","WikiPathways_2016")


enrichedU <- enrichr(Upconc, dbs)
enrichedD <- enrichr(Downconc, dbs)


enrichedU <- enrichr(UPgenes, dbs)
enrichedD <- enrichr(DOWNgenes, dbs)

Unrich <- enrichedU[["GO_Biological_Process_2018"]]
Dnrich <- enrichedD[["GO_Biological_Process_2018"]]

#npgender<-enrichr(DEGsgenderNP,dbs)
#npgender<-npgender[["GO_Biological_Process_2018"]]

#write.table(enrichedD_KO, "/Volumes/grcmc/CRISTINA/Yolanda_RNAseq/DOWN_GO_B2B3_pval001.csv" , sep = "\t",quote = FALSE,row.names = FALSE)
#write.table(enrichedU_KO, "/Volumes/grcmc/CRISTINA/Yolanda_RNAseq/UP_GO_B2B3_pval001.csv" , sep = "\t",quote = FALSE,row.names = FALSE)

Unrich$group<-c("Up")
Dnrich$group<-c("Down")
#npgender$group<-c("Gender")

allGO<-rbind(Unrich,Dnrich)
#allGO<-npgender

bpsub<-subset(allGO,allGO$Adjusted.P.value<=0.01)
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
bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$group, bpsub$Z.score)), ]$Term)

ggplot(bpsub,aes(y=Z.score,x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group),color="black")+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0,color="black")+
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
        axis.text.y = element_text(size=6, hjust = 1))


## Subset functions
bpsub_sel<-bpsub[which(bpsub$GO %in% c("0042776","0070129","0043248","1901532","0090263",
                                       "0060071","0010972","0009060","2000736","0035567","0035722","0070498",
                                       "0010389","0043488","0000398","0030178","1901991","0007077","0051436",
                                       "0002220","0071349","0000387","0006336","1902750","0038061","0090090",
                                       "0071456","0060828","0030177","0019221","0045333","0050852","0007005",
                                       "0044772","0000165","0046777","0016310","1902533","0043123","0006357",
                                       "0043122","0006468")),]
bpsub_sel$class<-c("mit","mit","prot","stem","wnt","wnt","cyc","mit","stem","wnt","il","il","cyc","splice",
                   "splice","wnt","cyc","cyc","cyc","sig","il","splice","cyc","cyc","sig","wnt","hyp",
                   "wnt","wnt","sig","mit","sig","mit","cyc","sig","pho","pho","sig","sig","rna","sig","pho")

bpsub_sel$colorclass<-c("green","grey","coral","grey","grey","green","green","grey","pink","darkolivegreen","grey","darkolivegreen","grey","blue",
                        "blue","purple","blue","grey","pink","pink","yellow","grey","grey","pink","pink","pink","blue",
                        "yellow","yellow","pink","grey","grey","blue","orange","darkolivegreen","pink","blue","blue","orange","red","darkolivegreen","darkolivegreen")


bpsub_sel$Term <-reorder(bpsub_sel$Term, rev(bpsub_sel$P.value))
bpsub_sel$Term <- factor(bpsub_sel$Term, levels=bpsub_sel[rev(order(bpsub_sel$group, bpsub_sel$Z.score)), ]$Term)

bpsub_sel$class<-as.factor(bpsub_sel$class)


ggplot(bpsub_sel,aes(y=Z.score,x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group),color="black")+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0.5,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Blues")+
  #facet_wrap(~group,scales="free",ncol=1)+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=9, hjust = 1,colour=bpsub_sel$colorclass),
        axis.text.x = element_text(size=8, hjust = 1))


## CUFFDIFF ##
cuff_NP<-read.delim("../cufflinks/NP_gene_exp.diff")
cuff_B2B3<-read.delim("../cufflinks/B2B3_gene_exp.diff")
cuff_bcat<-read.delim("../cufflinks/Bcat_gene_exp.diff")

cuff_NP_sig<-subset(cuff_NP,cuff_NP$q_value<0.1)
sigNP_D_cuff<-subset(cuff_NP_sig,cuff_NP_sig$log2.fold_change.<0)
NP_Down_cuff<-sigNP_D_cuff$gene
sigNP_U_cuff<-subset(cuff_NP_sig,cuff_NP_sig$log2.fold_change.>0)
NP_Up_cuff<-sigNP_U_cuff$gene

cuff_B2B3_sig<-subset(cuff_B2B3,cuff_B2B3$p_value<=0.05)
sigB2B3_D_cuff<-subset(cuff_B2B3_sig,cuff_B2B3_sig$log2.fold_change.<0)
B2B3_Down_cuff<-sigB2B3_D_cuff$gene
sigB2B3_U_cuff<-subset(cuff_B2B3_sig,cuff_B2B3_sig$log2.fold_change.>0)
B2B3_Up_cuff<-sigB2B3_U_cuff$gene

cuff_bcat_sig<-subset(cuff_bcat,cuff_bcat$q_value<0.1)
sigbcat_D_cuff<-subset(cuff_bcat_sig,cuff_bcat_sig$log2.fold_change.<0)
bcat_Down_cuff<-sigbcat_D_cuff$gene
sigbcat_U_cuff<-subset(cuff_bcat_sig,cuff_bcat_sig$log2.fold_change.>0)
bcat_Up_cuff<-sigbcat_U_cuff$gene

## number of genes expressed in each population
maxNP<-apply(subset(cuff_NP,select=c(value_2)), 1, FUN=max)
expNP <- maxNP[maxNP>1 & !is.nan(maxNP)]
length(expNP)

maxB2B3<-apply(subset(cuff_B2B3,select=c(value_2)), 1, FUN=max)
expB2B3 <- maxB2B3[maxB2B3>1 & !is.nan(maxB2B3)]
length(expB2B3)

maxbcat<-apply(subset(cuff_bcat,select=c(value_2)), 1, FUN=max)
expbcat <- maxbcat[maxbcat>1 & !is.nan(maxbcat)]
length(expbcat)



###### Hallmark functions C Ruiz #####
# B1B2B3
listfun<-read_xlsx("../R_plots/B1_B2_B3/hallmarks_Ruiz_excel.xlsx")
listfun$GeneSetname<-gsub('_',' ',listfun$GeneSetname)

listfun$GeneSetname <-reorder(listfun$GeneSetname, rev(listfun$P_value))


ggplot(listfun,aes(y=-log(P_value),x=GeneSetname),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",color="black",fill="dodgerblue",alpha=0.5)+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0.5,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  #scale_fill_brewer(palette = "Blues")+
  #facet_wrap(~group,scales="free",ncol=1)+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=14, hjust = 1),
        axis.text.x = element_text(size=14, hjust = 1))+
  xlab("Downregulated patwhays in Rbpj KO")



# NP
listfun<-read_xlsx("../R_plots/B1_B2_B3/hallmarks_Ruiz_NP_excel.xlsx")
listfun$GeneSetname<-gsub('_',' ',listfun$GeneSetname)

listfun$GeneSetname<-as.character(listfun$GeneSetname)
listfun$P_value<-as.numeric(listfun$P_value)
listfun$GeneSetname <-reorder(listfun$GeneSetname, rev(listfun$P_value))





ggplot(listfun,aes(y=-log(P_value),x=GeneSetname),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",color="black",fill="dodgerblue",alpha=0.5)+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0.5,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  #scale_fill_brewer(palette = "Blues")+
  #facet_wrap(~group,scales="free",ncol=1)+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=14, hjust = 1),
        axis.text.x = element_text(size=14, hjust = 1))+
  xlab("Downregulated patwhays in Rbpj KO")


## heatmap from functions GENEsetname C. Ruiz
# Here I included pathways from hallmarks excel Cristina, and two pathways from enrichR BP, Wnt and hematopoietic

#for B1B2B3
pathsel<-read.delim("../R_plots/B1_B2_B3/hallmarks_Ruiz.txt",header = FALSE)

#for NP (enriched pathways hallmark Cristina common with NP Trumpp)
pathsel<-read.delim("../R_plots/B1_B2_B3/hallmarks_Ruiz_NP.txt",header = FALSE)

colnames(pathsel)<-c("gene","pathway")
pathsel$gene<-paste0(toupper(substr(as.character(pathsel$gene), 1, 1)), tolower(substr(as.character(pathsel$gene), 2, nchar(as.character(pathsel$gene)))))

tab_path<-merge(deg_rld_df,pathsel,by.x="row.names",by.y="gene")

names(tab_path)[1]<-"gene"
#B1B2B3
#tab_path<-tab_path[,c(1,10,5,6,7,8,2,3,4,9)]
#NP
tab_path<-tab_path[,c(1,8,5,6,7,2,3,4)]

#sc_tab_path<-as.data.frame(scale(tab_path[, -c(1,2)]))
sc_tab_path<-as.data.frame(t(apply(tab_path[, -c(1,2)], 1, rescale)))

sc_tab_path$pathway<-tab_path$pathway
sc_tab_path$gene<-tab_path$gene

sc_tab_path<-melt(sc_tab_path,id.vars=c("gene","pathway"))
sc_tab_path$gene<-as.character(sc_tab_path$gene)

ggplot(sc_tab_path, aes(x=gene,y=variable)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="grey60", high="red",mid="white",midpoint = 0.5)+
  xlab("Genes") +
  #ylab("Batch") +
  theme_minimal()+
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=8,face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, hjust = 1,face="italic")) +
  labs(fill = "Expression level fold change")+
  facet_wrap(~pathway,scales="free",ncol=2)+
  xlab("")+
  ylab("")+
  coord_flip()



