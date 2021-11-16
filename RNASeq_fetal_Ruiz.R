
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

setwd("/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/HTSEQ/")

output.dir="/Volumes/grcmc/YGUILLEN/RNASeq_CRuiz/HTSEQ"

# Import metadata files
fmetadata_all<-read.delim("metadata_all.txt",header = FALSE)
colnames(fmetadata_all)<-c("sampleID","countFile","batch","condition","population","gender","stage")

# metadata for aHSCs, dHSCs and MPP
metadata_Trumpp<-read.delim("metadata_Trumpp.txt",header=FALSE)
colnames(metadata_Trumpp)<-c("sampleID","countFile","batch","condition","population","gender")
metadata_Trumpp<-subset(metadata_Trumpp,(metadata_Trumpp$batch=="proHSC" | metadata_Trumpp$batch=="proMPP"))
metadata_Trumpp<-subset(metadata_Trumpp,(metadata_Trumpp$condition!="aHSC"))


#Discard samples KO_4 and KO_6
fmetadata_all<-subset(fmetadata_all,(fmetadata_all$sampleID!="B3_KO_4_18821_GTCCGC" & 
                                     fmetadata_all$sampleID!="B3_KO_6_18823_ATCACG"))

#exclude "false" fetal KO
fmetadata_all<-subset(fmetadata_all,fmetadata_all$sampleID!="S16")

## For fetal KO vs. WT
metadata_fetal_KO<-subset(fmetadata_all,fmetadata_all$population=="FE" & fmetadata_all$condition!="NICD")
metadata_fetal_KO<-fmetadata_fetal_KO[fmetadata_fetal_KO$sampleID!="S16",]

## For fetal Notch vs. WT
metadata_fetal_NICD<-subset(metadata_all,metadata_all$population=="FE" & (metadata_all$condition=="WT" | metadata_all$condition=="NICD"))

## For adult KO vs. WT
metadata_adult<-subset(metadata_all,metadata_all$population=="A")

# For adult WT vs fetal WT
metadata_WT<-subset(metadata_all,metadata_all$condition=="WT" & (metadata_all$population=="FE" | metadata_all$population=="A" ))

# For HSCs vs MPP1
metadata_Trumpp

# For WT B1B2B3 vs NP
metadata_pop<-metadata_all[c(4,5,6,7,12,13,14),]




DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata_all,
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

#if htseqgene symbol
rowData <- merge(DESeq2Features, bm, by.x="id",by.y = "external_gene_name")

#if not
rowData <- merge(DESeq2Features, bm, by.x="id",by.y = "ensembl_gene_id")


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
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","KO","NICD"))

#for fetal KO vs WT and adult
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","KO"))

#for fetal NICD vs WT
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","NICD"))

#for fetal vs adult WT
#colData(DESeq2Table)$stage <- factor(colData(DESeq2Table)$condition, levels = c("AD","FE"))

#for dHSC vs MPP1
colData(DESeq2Table)$population <- factor(colData(DESeq2Table)$population, levels = c("MPP1","dHSC"))


#for WT B1B2B3 vs NP
colData(DESeq2Table)$population <- factor(colData(DESeq2Table)$population, levels = c("A","NP"))


#### estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)


### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Table, blind=TRUE) ## create a distance matrix between the samples




#sub rld for retinoic acid
subRA<-rld[row.names(rld) %in% RAsign,]
subRA
subRA_df<-as.data.frame(assay(subRA))


## Create PCA (check if for all or for RA sign)
pcaData <- plotPCA(rld, intgroup=c("condition","population"), returnData=TRUE)

#pcaData <- plotPCA(rld, intgroup=c("stage","gender"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition,label=name)) +
  geom_point(color="black",size=10,alpha=0.3) +
  geom_point(size=8,alpha=0.8) +
  stat_ellipse(aes(color=condition),level=0.8)+
  #scale_color_manual(values=c("dodgerblue","coral"))+
  #scale_color_manual(values=c("grey","dodgerblue","coral","darkolivegreen2","black","dodgerblue4","coral4"))+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #facet_grid(~gender)+
  scale_color_manual(values=c("dodgerblue2","coral","darkolivegreen"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed()


ggplot(pcaData, aes(y=PC1, x=condition,label=name)) +
  geom_point(aes(color=condition,shape=population),size=8) +
  geom_boxplot(alpha=0)+
  #stat_ellipse(aes(color=condition),level=0.8)+
  scale_color_manual(values=c("dodgerblue2","red","darkolivegreen"))+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
  xlab("")+
  #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  facet_wrap(~population,scales = "free")+
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
levels(DESeq2Table$stage)
levels(DESeq2Table$gender)
levels(DESeq2Table$population)

design(DESeq2Table)

dds<-DESeq(DESeq2Table)
resultsNames(dds)

rescondition<-results(dds)

# EXTRACT COUNTS
countab<-counts(dds, normalized=TRUE)
countabRA<-countab[row.names(countab) %in% RAsign,]
rownames<-row.names(countabRA)
cols<-colnames(countabRA)

countabRA<-as.data.frame(t(countabRA))
countabRA<-merge(countabRA,metadata_all,by.x="row.names",by.y="sampleID")
countabRA<-melt(countabRA,id.vars=c("Row.names","countFile","batch","condition","population","gender","stage"))
countabRA$group<-paste(countabRA$population,countabRA$condition,sep="_")
  
countabRA$group<-as.factor(countabRA$group)
countabRA$group <- factor(countabRA$group, levels = c("A_WT", "A_KO", "FE_WT","FE_KO","FE_NICD","NP_WT","NP_KO"))

ggplot(countabRA,aes(x=group,y=value))+
  geom_boxplot(aes(color=condition,linetype=population))+
  geom_point(aes(color=condition))+
  scale_color_manual(values=c("darkolivegreen","dodgerblue","coral"))+
  facet_wrap(~variable,scales="free")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),legend.position="bottom")


ggplot(data=countabRA[countabRA$population=="A",],aes(x=group,y=value))+
  geom_boxplot()+
  geom_point(aes(color=condition))+
  scale_color_manual(values=c("coral","darkolivegreen"))+
  facet_wrap(~variable,scales="free")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),legend.position="bottom")


ggplot(data=countabRA[countabRA$population=="FE",],aes(x=group,y=value))+
  geom_boxplot()+
  geom_point(aes(color=condition))+
  scale_color_manual(values=c("coral","dodgerblue","darkolivegreen"))+
  facet_wrap(~variable,scales="free")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),legend.position="bottom")


ggplot(data=countabRA[countabRA$population=="NP",],aes(x=group,y=value))+
  geom_boxplot()+
  geom_point(aes(color=condition))+
  scale_color_manual(values=c("coral","darkolivegreen"))+
  facet_wrap(~variable,scales="free")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),legend.position="bottom")


# Example plot counts
plotCounts(dds, gene=which(rownames(rescondition)=="Hoxa5"), intgroup="group")
plotCounts(dds, gene=which(rownames(rescondition)=="ENSMUSG00000039191"), intgroup="population")
plotCounts(dds, gene=which(rownames(rescondition)=="Lgals1"), intgroup="stage")

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

res_df_sig_rna<-subset(res_df_rna,res_df_rna$pvalue<0.05)
res_df_sig_rna <- res_df_sig_rna[order(res_df_sig_rna$pvalue),]


# UP and DOWN genes
UPgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange>0))$id)
DOWNgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange<0))$id)
length(UPgenes)
length(DOWNgenes)

# Heatmap
#Manual heatmap from log trasformed data, using all samples
rld_df<-as.data.frame(assay(rld))

#only degs
deggenes<-unique(res_df_sig_rna$id)

deg_rld<-which(rownames(rld_df) %in% deggenes)
deg_rld_df<-rld_df[deg_rld,]

deg_rld<-merge(deg_rld_df,res_df_rna,by.x="row.names",by.y="id")

deg_rld<- deg_rld[order(deg_rld$log2FoldChange),]
deg_rld$Row.names

library(dplyr)
library(tidyverse)

deg_rld<-deg_rld[!duplicated(deg_rld), ]

heatmap.2( as.matrix(deg_rld[,2:6]), scale="row", dendrogram="column",trace="none", Rowv=FALSE,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           cexCol=0.3)



###### Save DEGs list

#all genes
deggenes<-unique(res_df_rna$id)
deg_rld<-which(rownames(rld_df) %in% deggenes)
deg_rld_df<-rld_df[deg_rld,]

deg_rld<-merge(deg_rld_df,res_df_rna,by.x="row.names",by.y="id")


# WT fetal vs WT adult
Wadfet_DEGs<-deg_rld
Wadfet_DEGs_comp<-subset(Wadfet_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
#colnames(KWfet_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
Wadfet_DEGs_comp<-Wadfet_DEGs_comp[!duplicated(Wadfet_DEGs_comp$Row.names),]
Wadfet_DEGs_comp$population<-"WTfetad"

UpWTfet<-unique((Wadfet_DEGs_comp[Wadfet_DEGs_comp$padj<0.1 & Wadfet_DEGs_comp$log2FoldChange>0,])$Row.names)
DownWTfet<-unique((Wadfet_DEGs_comp[Wadfet_DEGs_comp$padj<0.1 & Wadfet_DEGs_comp$log2FoldChange<0,])$Row.names)

# WT vs KO adult
KWad_DEGs<-deg_rld
KWad_DEGs_comp<-subset(KWad_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
#colnames(KWfet_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
KWad_DEGs_comp<-KWad_DEGs_comp[!duplicated(KWad_DEGs_comp$Row.names),]
KWad_DEGs_comp$population<-"KOWTad"


# WT vs KO fetal
KWfet_DEGs<-deg_rld
KWfet_DEGs_comp<-subset(KWfet_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
KWfet_DEGs_comp<-KWfet_DEGs_comp[!duplicated(KWfet_DEGs_comp$Row.names),]
#colnames(KWfet_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
KWfet_DEGs_comp$population<-"KOWTfet"

#WT vs Notch
NWfet_DEGs<-deg_rld
NWfet_DEGs_comp<-subset(NWfet_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
NWfet_DEGs_comp<-NWfet_DEGs_comp[!duplicated(NWfet_DEGs_comp$Row.names),]
#colnames(KWfet_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
NWfet_DEGs_comp$population<-"NWTfet"


#dHSCs vs MPP1
dHMPP_DEGs<-deg_rld
dHMPP_DEGs_comp<-subset(dHMPP_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
dHMPP_DEGs_comp<-dHMPP_DEGs_comp[!duplicated(dHMPP_DEGs_comp$mgi_symbol),]
#colnames(KWfet_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
dHMPP_DEGs_comp$population<-"dHMPP"


#WT B1B2B3 vs WT NP
WTpop_DEGs<-deg_rld
WTpop_DEGs_comp<-subset(WTpop_DEGs,select=c(Row.names,log2FoldChange,pvalue,padj,gene_biotype,description,mgi_symbol))
WTpop_DEGs_comp<-WTpop_DEGs_comp[!duplicated(WTpop_DEGs_comp$Row.names),]
#colnames(KWfet_DEGs_comp)<-c("Row.names","log2FoldChange","pvalue","padj","gene_biotype","description","mgi_symbol")
WTpop_DEGs_comp$population<-"WTpop"



#Merge all batches for WT vs KO adult, fetal and Notch vs KO
NKW_comp<-rbind(KWfet_DEGs_comp,NWfet_DEGs_comp,KWad_DEGs_comp)



#Select DEgs padj < 0.1 from here
KW_sig<-unique(as.character((subset(KWfet_DEGs_comp,KWfet_DEGs_comp$pvalue<=0.05))$Row.names))
KWad_sig<-unique(as.character((subset(KWad_DEGs_comp,KWad_DEGs_comp$pvalue<=0.05))$Row.names))
NW_sig<-unique(as.character((subset(NWfet_DEGs_comp,NWfet_DEGs_comp$pvalue<=0.05))$Row.names))


## Heatmap

#for fetal
alldegs<-unique(c(KW_sig,NW_sig))
length(alldegs)

#for fetal and adult
alldegs<-unique(c(KW_sig,NW_sig,KWad_sig))
length(alldegs)


NKW_comp<-NKW_comp[(NKW_comp$Row.names %in% alldegs),]
NKW_comp$state<-ifelse(NKW_comp$log2FoldChange< 0, 
                        c("Conddown"), c("Condup")) 


NKW_comp$state<-as.factor(NKW_comp$state)
NKW_comp<- NKW_comp[order(NKW_comp$population,NKW_comp$pvalue),]

NKW_comp$Row.names<-as.character(NKW_comp$Row.names)

ggplot(NKW_comp, aes(Row.names,population)) +
  geom_tile(aes(fill = state), color = "white") +
  geom_point(aes(size=(-log(pvalue))),data = subset(NKW_comp,NKW_comp$padj<0.1))+
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


## Correlation FC WTKO and NotchWT and others

KWfet_comp<-subset(NKW_comp,NKW_comp$population=="KOWTfet")
NWfet_comp<-subset(NKW_comp,NKW_comp$population=="NWTfet")
KWad_comp<-subset(NKW_comp,NKW_comp$population=="KOWTad")


merge_WNK<-merge(KWfet_comp,NWfet_comp,by="Row.names",all=TRUE)
merge_WNK<-merge(merge_WNK,KWad_comp,by="Row.names",all=TRUE)


# CHECK POSITIONS OF LOG2FOLDCHANGE
names(merge_WNK)[2]<-"log2FoldChangeKWfet"
names(merge_WNK)[10]<-"log2FoldChangeNWfet"
names(merge_WNK)[18]<-"log2FoldChangeKWad"

merge_WNK<-merge_WNK[order(merge_WNK$pvalue.x),]


merge_WNK$sigWK<-ifelse(merge_WNK$pvalue.x < 0.05 & !is.na(merge_WNK$pvalue.x), 
                        c("sig_WKfet"), c("nosig_WKfet")) 

merge_WNK$sigWN<-ifelse(merge_WNK$pvalue.y < 0.05 & !is.na(merge_WNK$pvalue.y), 
                        c("sig_WNfet"), c("nosig_WNfet")) 

merge_WNK$sigWKad<-ifelse(merge_WNK$pvalue < 0.05 & !is.na(merge_WNK$pvalue), 
                        c("sig_WKad"), c("nosig_WKad")) 

merge_WNK$sigKO<-ifelse(merge_WNK$sigWKad=="sig_WKad" & merge_WNK$sigWK =="sig_WKfet", 
                          c("sig_WKall"), c("nosig_WKall"))

merge_WNK$sigfet<-ifelse(merge_WNK$sigWK=="sig_WKfet" & merge_WNK$sigWN =="sig_WNfet", 
                        c("sig_fet"), c("nosig_fet"))



merge_WNK$sig<-paste(merge_WNK$sigWK,merge_WNK$sigWN,sep="_")
merge_WNK$sig<-paste(merge_WNK$sig,merge_WNK$sigWKad,sep="_")

merge_WNK$sigtype<-ifelse(merge_WNK$sig == "sig_WKfet_sig_WNfet_sig_WKad", 
                                   c("sig_all"), c("nosig_all"))  

merge_WNK<-merge_WNK[!duplicated(merge_WNK),]

ggplot(subset(merge_WNK,merge_WNK$pvalue.x<=0.05 | merge_WNK$pvalue.y<=0.05),aes(x=log2FoldChangeKWfet,y=log2FoldChangeNWfet))+
  #annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "blue",alpha=0.1)+
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "red",alpha=0.1)+
  geom_point(color="black",size=3,alpha=0.5)+
  #geom_point(data=subset(merge_WNK,merge_WNK$sigfet=="sig_fet"),color="red",size=1)+
  geom_point(data=subset(merge_WNK,merge_WNK$sigfet=="sig_fet"),color="red",size=1)+
  geom_density_2d(color="green")+
  #geom_density_2d(data=(subset(merge_WNK,merge_WNK$sigKO=="sig_WKall")),color="red")+
  #geom_density_2d(data=(subset(merge_WNK,merge_WNK$sigfet=="sig_fet")),color="red")+
  geom_smooth(data=(subset(merge_WNK,merge_WNK$sigfet=="sig_fet")),method="lm",color="red")+
  #scale_color_manual(values=c("black","grey"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = Row.names),color="black", data = subset(merge_WNK,merge_WNK$sigtype=="sig_all"), fontface = "italic",size=5)+
  geom_label_repel(aes(label = Row.names),color="black", data = merge_WNK[merge_WNK$Row.names %in% degscond,], fontface = "italic",size=5)+
  xlim(-5,5)+
  ylim(-5,5)+
  xlab("log2FC KO vs. WT fetal")+
  ylab("log2FC NICD vs. WT fetal")+
  #ylab("log2FC KO vs. WT adult")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1)) 


cor.test((subset(merge_WNK,merge_WNK$sigfet=="sig_fet"))$log2FoldChangeKWfet,(subset(merge_WNK,merge_WNK$sigfet=="sig_fet"))$log2FoldChangeNWfet)

# 3D plot
p <- plot_ly(subset(merge_WNK,merge_WNK$pvalue.x<=0.05 | merge_WNK$pvalue.y<=0.05), x = ~log2FoldChangeKWfet, y = ~log2FoldChangeNWfet, z = ~log2FoldChangeKWad, color = ~sigfet) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'KO vs WT fetal'),
                      yaxis = list(title = 'NICD vs WT fetal'),
                      zaxis = list(title = 'KO vs WT adult')))
p



# Group of genes with no expression in one condition, and significant FC in other.
unique((merge_WNK[merge_WNK$sig=="NA_sig_WNfet_NA" | merge_WNK$sig=="NA_sig_WNfet_nosig_WKad",])$Row.names)

## Functional enrichment of UP-UP significant and DOWN-DOWN significant
Upconc<-unique((subset(merge_WNK,merge_WNK$pvalue.x<=0.05 & 
                  merge_WNK$log2FoldChangeKWfet>0 &
                    merge_WNK$pvalue.y<=0.05 &
                    merge_WNK$log2FoldChangeNWfet>0 &
                    merge_WNK$pvalue<=0.05 &
                  merge_WNK$log2FoldChangeKWad>0)$Row.names))
Downconc<-unique((subset(merge_WNK,merge_WNK$pvalue.x<=0.05 & 
                         merge_WNK$log2FoldChangeKWfet<0 &
                         merge_WNK$pvalue.y<=0.05 &
                         merge_WNK$log2FoldChangeNWfet<0 &
                         merge_WNK$pvalue<=0.05 &
                         merge_WNK$log2FoldChangeKWad<0)$Row.names))

degscond<-c(Upconc,Downconc)

# Only for UP and DOWN in KO vs WT fetal
UpKOfet<-unique(na.omit(KWfet_DEGs_comp[KWfet_DEGs_comp$pvalue<0.05 & KWfet_DEGs_comp$log2FoldChange>0,])$Row.names)
DownKOfet<-unique(na.omit(KWfet_DEGs_comp[KWfet_DEGs_comp$pvalue<0.05 & KWfet_DEGs_comp$log2FoldChange<0,])$Row.names)

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2019_Mouse",
         "KEGG_2019_Mouse")


enrichedU <- enrichr(Upconc, dbs)
enrichedD <- enrichr(Downconc, dbs)

Unrich <- enrichedU[["KEGG_2019_Mouse"]]
Dnrich <- enrichedD[["KEGG_2019_Mouse"]]

Unrich$group<-c("Up")
Dnrich$group<-c("Down")

allGO<-rbind(Unrich,Dnrich)

bpsub<-subset(allGO,allGO$P.value<=0.01)
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
bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$group, bpsub$Combined.Score)), ]$Term)


ggplot(bpsub,aes(y=Combined.Score,x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group),color="black")+
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



#### Compare MPP1 vs dHSC to infer MPP1 gene signature

dHMPP_sig<-unique(as.character((subset(dHMPP_DEGs_comp,dHMPP_DEGs_comp$pvalue<=0.05))$mgi_symbol))
WTpop_sig<-unique(as.character((subset(WTpop_DEGs_comp,WTpop_DEGs_comp$pvalue<=0.05))$Row.names))

length(dHMPP_sig)
length(WTpop_sig)


alldegs<-unique(c(dHMPP_sig,WTpop_sig))
length(alldegs)

dHMPP_DEGs_comp$Row.names=dHMPP_DEGs_comp$mgi_symbol

MPWT_comp<-rbind(dHMPP_DEGs_comp,WTpop_DEGs_comp)

MPWT_comp<-MPWT_comp[MPWT_comp$Row.names %in% alldegs,]

MPWT_comp$state<-ifelse(MPWT_comp$log2FoldChange< 0, 
                       c("Down"), c("Up")) 


## Correlation FC WTKO and NotchWT and others

dHMPP_comp<-subset(MPWT_comp,MPWT_comp$population=="dHMPP")
WTpop_comp<-subset(MPWT_comp,MPWT_comp$population=="WTpop")

# For WTs population
merge_WTH<-merge(WTpop_comp,dHMPP_comp,by="Row.names",all=TRUE)

# CHECK POSITIONS OF LOG2FOLDCHANGE
names(merge_WTH)[2]<-"log2FoldChangeWTANP"
names(merge_WTH)[10]<-"log2FoldChangeHSCM"

merge_WTH<-merge_WTH[order(merge_WTH$pvalue.x),]


merge_WTH$sigWTpop<-ifelse(merge_WTH$padj.x < 0.1 & !is.na(merge_WTH$padj.x), 
                        c("sig_WTpopfet"), c("nosig_WTpopfet")) 

merge_WTH$sigMPPH<-ifelse(merge_WTH$padj.y < 0.1 & !is.na(merge_WTH$padj.y), 
                        c("sig_MPPH"), c("nosig_MPPH")) 


merge_WTH$sig<-paste(merge_WTH$sigWTpop,merge_WTH$sigMPPH,sep="_")

merge_WTH$sigtype<-ifelse(merge_WTH$sig == "sig_WTpopfet_nosig_MPPH", 
                          c("sig_all"), c("nosig_all"))  

ggplot(subset(merge_WTH,merge_WTH$padj.x<=0.1 | merge_WTH$padj.y<=0.1),aes(x=log2FoldChangeWTANP,y=log2FoldChangeHSCM))+
  #annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "blue",alpha=0.1)+
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "red",alpha=0.1)+
  geom_point(color="black",size=3,alpha=0.5)+
  geom_point(data=subset(merge_WTH,merge_WTH$sigtype=="sig_all"),color="red",size=1)+
  #geom_point(data=subset(merge_WNK,merge_WNK$sigKO=="sig_WKall"),color="red",size=1)+
  geom_density_2d(data=subset(merge_WTH,merge_WTH$sigtype=="sig_all"),color="green")+
  #geom_density_2d(data=(subset(merge_WNK,merge_WNK$sigKO=="sig_WKall")),color="red")+
  #geom_density_2d(data=(subset(merge_WNK,merge_WNK$sigfet=="sig_fet")),color="red")+
  geom_smooth(data=(subset(merge_WTH,merge_WTH$sigtype=="sig_all")),method="lm",color="red")+
  #scale_color_manual(values=c("black","grey"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = mgi_symbol),color="black", data = subset(merge_WTH,merge_WTH$sigtype=="sig_all"), fontface = "italic",size=5)+
  xlim(-10,5)+
  ylim(-10,10)+
  xlab("log2FC WT all vs. WT NP")+
  ylab("log2FC MPP1 vs. dHSC")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1))

cor.test((subset(merge_WTH,merge_WTH$sigtype=="sig_all"))$log2FoldChangeWTANP,(subset(merge_WTH,merge_WTH$sigtype =="sig_all"))$log2FoldChangeHSCM)

# 3D plot
p <- plot_ly(subset(merge_WNK,merge_WNK$pvalue.x<=0.05 | merge_WNK$pvalue.y<=0.05), x = ~log2FoldChangeKWfet, y = ~log2FoldChangeNWfet, z = ~log2FoldChangeKWad, color = ~sigfet) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'KO vs WT fetal'),
                      yaxis = list(title = 'NICD vs WT fetal'),
                      zaxis = list(title = 'KO vs WT adult')))
p



#Combination
length(unique((merge_WTH[merge_WTH$sig!="nosig_WTpopfet_nosig_MPPH" & merge_WTH$log2FoldChangeWTANP<0 & merge_WTH$log2FoldChangeHSCM<0,])$Row.names))
length(unique((merge_WTH[merge_WTH$sig!="nosig_WTpopfet_nosig_MPPH" & merge_WTH$log2FoldChangeWTANP<0 & merge_WTH$log2FoldChangeHSCM>0,])$Row.names))
length(unique((merge_WTH[merge_WTH$sig!="nosig_WTpopfet_nosig_MPPH" & merge_WTH$log2FoldChangeWTANP>0 & merge_WTH$log2FoldChangeHSCM>0,])$Row.names))
length(unique((merge_WTH[merge_WTH$sig!="nosig_WTpopfet_nosig_MPPH" & merge_WTH$log2FoldChangeWTANP>0 & merge_WTH$log2FoldChangeHSCM<0,])$Row.names))


#Individual
length(unique((merge_WTH[merge_WTH$sig!="nosig_WTpopfet_nosig_MPPH" & merge_WTH$log2FoldChangeWTANP>0,])$Row.names))
length(unique((merge_WTH[merge_WTH$sig!="nosig_WTpopfet_nosig_MPPH" & merge_WTH$log2FoldChangeWTANP<0,])$Row.names))
length(unique((merge_WTH[merge_WTH$sig!="nosig_WTpopfet_nosig_MPPH" & merge_WTH$log2FoldChangeHSCM>0,])$Row.names))
length(unique((merge_WTH[merge_WTH$sig!="nosig_WTpopfet_nosig_MPPH" & merge_WTH$log2FoldChangeHSCM<0,])$Row.names))


## Manually creation of dataframe for alluvial plots


# Create dataframe to generate flow diagram
flowdig <- data.frame("NP" = c("up","up","down","down"), 
                      "dHSC" = c("up","down","up","down"), 
                      "Freq" = c(401,519,750,430))
flowdig$state<-c("conc","disc","disc","conc")

#flowdig <- data.frame("Population" = c("NP","NP","dHSC","dHSC"), 
#                      "State" = c("up","down","up","down"), 
#                      "Freq" = c(920,1180,1151,949))

library(ggalluvial)

ggplot(flowdig,aes(y = Freq, axis1=NP, axis2=dHSC,fill=state)) +
  geom_flow(aes(fill=state),aes.bind = TRUE, reverse = FALSE) +
  scale_fill_manual(values=c("dodgerblue","red"))+
  geom_stratum(reverse = FALSE) +
  geom_text(stat = "stratum", label.strata = TRUE, reverse = FALSE) +
  scale_x_discrete(limits = c("NP", "dHSC"))+
  theme_minimal()


### Create dataframe for alluvial from each gene

## Chech concordant and discordant expression flow between conditions
flowtab1<-merge_WNK[,c(1,3,8,9)]
#flowtab1<-flowtab1[flowtab1$pvalue.x<=0.05,]
flowtab1<-flowtab1[,c(1,3,4)]
names(flowtab1)<-c("Row.names","population","state")
flowtab1<-flowtab1[!is.na(flowtab1$Row.names),]

flowtab2<-merge_WNK[,c(1,11,16,17)]
#flowtab2<-flowtab2[flowtab2$pvalue.y<=0.05,]
flowtab2<-flowtab2[,c(1,3,4)]
names(flowtab2)<-c("Row.names","population","state")
flowtab2<-flowtab2[!is.na(flowtab2$Row.names),]

flowtab3<-merge_WNK[,c(1,19,24,25)]
#flowtab3<-flowtab3[flowtab3$pvalue<=0.05,]
flowtab3<-flowtab3[,c(1,3,4)]
names(flowtab3)<-c("Row.names","population","state")
flowtab3<-flowtab3[!is.na(flowtab3$Row.names),]


## Add dHSC vs aHSC data (using B1B2B3_all_NP_comp dataframe)
flowtab4<-B1B2B3_all_NP_comp[,c(1,19,24,25)]
#flowtab4<-flowtab4[flowtab4$pvalue<=0.05,]
flowtab4<-flowtab4[,c(1,3,4)]
names(flowtab4)<-c("Row.names","population","state")
flowtab4<-flowtab4[!is.na(flowtab3$Row.names),]
flowtab4$state<-revalue(flowtab4$state, c("KOdown"="Conddown", "KOup"="Condup"))


## For MPP1/dHSC vs WT CD150+ / CD135-
flowtab5<-merge_WTH[,c(1,8,9)]
names(flowtab5)<-c("Row.names","population","state")
flowtab5<-flowtab5[!is.na(flowtab5$Row.names),]

flowtab6<-merge_WTH[,c(1,16,17)]
names(flowtab6)<-c("Row.names","population","state")
flowtab6<-flowtab6[!is.na(flowtab6$Row.names),]


# For WT AND KOs and others
#flowtab<-rbind(flowtab1,flowtab2,flowtab3,flowtab4)
#flowtab$state<-revalue(flowtab$state, c("Conddown"="Down", "Condup"="Up"))

# for MPP1 vs dHSCs
flowtab<-rbind(flowtab5,flowtab6)


flowtab<-flowtab[!duplicated(flowtab),]


flowtab$state <- as.factor(flowtab$state)
flowtab$Row.names<-as.factor(flowtab$Row.names)
flowtab$population<-as.factor(flowtab$population)

flowtab<-flowtab[!is.na(flowtab$state),]

ggplot(flowtab,
       aes(x = population, stratum = state, alluvium = Row.names,
           fill = state, label = state)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray",alpha=0.5) +
  geom_stratum() +
  theme(legend.position = "bottom")+
  theme_bw()

devtools::install_github("erblast/easyalluvial")
require(easyalluvial)
require(tidyverse)


# MANUAL
# https://www.datisticsblog.com/2018/10/intro_easyalluvial/#easyalluvial

col_vector = c('coral','grey','dodgerblue')

p<-alluvial_long(flowtab
              , key = population
              , value = state
              , id = Row.names
              , fill_by = 'value'
              , col_vector_flow = col_vector
              ,NA_label = 'no_match'
              , col_vector_value = col_vector
)


p+
  theme_minimal()+
  theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12))


## Retinoic acid signaling
RAsign<-c("Adam11","Ahnak","Aim2","Bambi","Btg2","Casc3","Cd38","Cdk4","Cebpe","Chn2",
          "Col1a1","Cyp26b1","Egr1","Ets1","Gbx2","Hgs","Hoxa5","Hoxb2","Hoxb4","Il6st","Itga4",
          "Mdk","Meis1","Meis2","Mfng","Mmp11","Nptn","Nrgn","Nrip1","Nrp1","Pik3c3","Pkp2","Prkch",
          "Ptch1","Pycr2","Rarb","Rbp1","Serpinb8","Serpinb9","Serping1","Slc5a3","Spp1","Stard10",
          "Tgm2","Tubgcp2","Zeb2")

RA_met<-c("Aldh1a1","Aldh1a2","Aldh1a3","Ces1g","Crabp2","Cyp26a1","Cyp26b1","Cyp26c1","Dhrs9","Fabp5","Lrat",
          "Rbp1","Rbp4","Rdh1","Rdh10","Rdh7","Stra6","Ttr")

RAdat<-merge_WNK[merge_WNK$Row.names %in% RAsign,]
RAmet<-merge_WNK[merge_WNK$Row.names %in% RA_met,]
