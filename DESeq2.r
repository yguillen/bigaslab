
### Pipeline DESeq2 RNAseq data ###

biocLite("pheatmap")
biocLite("vegan")
biocLite("tidyr")

library(pheatmap)
library(vegan)
library(tidyr)


setwd("/Volumes/grcmc/YGUILLEN/RNASeq_CPorch/")

# Example Sarah C. Porch
#metadata
samples<-read_xlsx("samples_summary.xlsx")
samples<-samples[-2,]

#HTSEq counts raw, NOT normalized
htseqraw<-read.delim("raw_htseq_counts.txt")
colnames(htseqraw)<-gsub('^X','',colnames(htseqraw))

#HTSEQ counts normalized
rnaex<-read.delim("normalized_log2_htseq_counts_annot.txt2",sep=" ")
rnaex$X<-NULL
colnames(rnaex)<-gsub('^X','',colnames(rnaex))

#DESEQ data
deseq<-read_xlsx("20171221_ABigas_RNAseq_DESeq2.xlsx")



## PCA using all gene expression data, normalized htseq
row.names(rnaex)<-rnaex$ensembl_id
sub<-rnaex[,8:ncol(rnaex)]
rnaex_t<-t(sub)

pca <- prcomp(rnaex_t, scale=F)
summary(pca)

pca_htseq<-as.data.frame(pca$x[,c(1:2)])
pca_htseq$Sample<-row.names(pca_htseq)
pca_data<-merge(samples,pca_htseq,by="Sample")

pca_data$type<-paste(pca_data$Condition,pca_data$Population,sep="_")

ggplot(pca_data,aes(x=PC1,y=PC2))+
  geom_point(aes(color=type),size=5)+
  stat_ellipse(aes(color=Population))+
  scale_color_manual(values=c("coral4","darkolivegreen","dodgerblue4","coral3","darkolivegreen2","dodgerblue","coral4","darkolivegreen","dodgerblue4"))+
  ggtitle("PCA all genes")+
  theme_bw()


## adonis test
rownames(rnaex_t)
colnames(rnaex_t)
rnaex_meta<-merge(samples,rnaex_t,by.x="Sample",by.y="row.names")

# All populations
rnaex_gene<-rnaex_meta[,6:ncol(rnaex_meta)]
adonis(rnaex_gene ~ Condition+Population, data = rnaex_meta)


#per population
rnaex_meta_pop1<-subset(rnaex_meta,rnaex_meta$Population=="45neg_kitneg_41pos")
rnaex_gene<-rnaex_meta_pop1[,6:ncol(rnaex_meta_pop1)]
adonis(rnaex_gene ~ Condition, data = rnaex_meta_pop1)

rnaex_meta_pop2<-subset(rnaex_meta,rnaex_meta$Population=="45neg_kitpos")
rnaex_gene<-rnaex_meta_pop2[,6:ncol(rnaex_meta_pop2)]
adonis(rnaex_gene ~ Condition, data = rnaex_meta_pop2)

rnaex_meta_pop3<-subset(rnaex_meta,rnaex_meta$Population=="45pos_kitpos")
rnaex_gene<-rnaex_meta_pop3[,6:ncol(rnaex_meta_pop3)]
adonis(rnaex_gene ~ Condition, data = rnaex_meta_pop3)


## select subset genes by population
table(samples$Population)

#population classificaton subsets
#pop1 = 45- kit- 41+
pop1<-(subset(samples,samples$Population=="45neg_kitneg_41pos")$Sample)
pop1control<-(subset(samples,samples$Population=="45neg_kitneg_41pos" & samples$Condition=="Control")$Sample)
pop1ab<-(subset(samples,samples$Population=="45neg_kitneg_41pos" & samples$Condition=="Ab")$Sample)

#pop2 45- kit+
pop2<-(subset(samples,samples$Population=="45neg_kitpos")$Sample)
pop2control<-(subset(samples,samples$Population=="45neg_kitpos" & samples$Condition=="Control")$Sample)
pop2ab<-(subset(samples,samples$Population=="45neg_kitpos" & samples$Condition=="Ab")$Sample)

#pop3 45+ kit+
pop3<-(subset(samples,samples$Population=="45pos_kitpos")$Sample)
pop3control<-(subset(samples,samples$Population=="45pos_kitpos" & samples$Condition=="Control")$Sample)
pop3ab<-(subset(samples,samples$Population=="45pos_kitpos" & samples$Condition=="Ab")$Sample)


pop1rna<-subset(rnaex,select=c(pop1))
pop1cont<-subset(rnaex,select=c(pop1control))
pop1treat<-subset(rnaex,select=c(pop1ab))

pop2rna<-subset(rnaex,select=c(pop2))
pop2cont<-subset(rnaex,select=c(pop2control))
pop2treat<-subset(rnaex,select=c(pop2ab))

pop3rna<-subset(rnaex,select=c(pop3))
pop3cont<-subset(rnaex,select=c(pop3control))
pop3treat<-subset(rnaex,select=c(pop3ab))






# HEATMAP expression
#ranked gene list pop2
rankedgenes<-read.delim("ranked_gene_list_pop45negkitpos.txt",header = FALSE)
rankedgenes<-rankedgenes[,c(1,5)]
colnames(rankedgenes)<-c("gene_name","score")

hist(rankedgenes$score,breaks = 60)
ggplot(rankedgenes,aes(x=score))+
  geom_density()+
  theme_bw()
rankedgenessig<-subset(rankedgenes,rankedgenes$score<(-4.5) | rankedgenes$score>4.5)
rankedgenessig$pos<-row.names(rankedgenessig)
rankedgenessig<-rankedgenessig[,c(1,3)]
rankedgenessig$pos<-as.numeric(rankedgenessig$pos)

#core notch genes GSEA population kit+ 45- (pop2)
corenotch<-read.delim("core_notch1GSEA_genes.txt",header=FALSE)
colnames(corenotch)<-c("gene_name")
corenotch$pos<-row.names(corenotch)
corenotch$pos<-as.numeric(corenotch$pos)

#pop2
pop2corecont<-subset(rnaex,select=c("gene_name",pop2control))
pop2corecont$type<-c("CONTROL")
pop2corecont<-melt(pop2corecont,id.vars=c("gene_name","type"))

pop2coretreat<-subset(rnaex,select=c("gene_name",pop2ab))
pop2coretreat$type<-c("AB")
pop2coretreat<-melt(pop2coretreat,id.vars=c("gene_name","type"))

pop2core<-rbind(pop2corecont,pop2coretreat)
pop2core$gene_name<-as.character(pop2core$gene_name)
corenotch$gene_name<-as.character(corenotch$gene_name)

#for ranked gene list
pop2corenotch<-merge(pop2core,rankedgenessig,by="gene_name",all.y=TRUE)

#for core notch
pop2corenotch<-merge(pop2core,corenotch,by="gene_name",all.y=TRUE)
#position in ranked gene list of core notch list
rankedgenes$pos<-row.names(rankedgenes)
rankednotch<-merge(corenotch,rankedgenes,by="gene_name",all.x=TRUE)
rankednotch$pos.y<-as.numeric(rankednotch$pos.y)
rankednotch<-rankednotch[order(rankednotch$pos.y), ]
hist(rankednotch$score,breaks = 60)
ggplot(rankednotch,aes(x=score))+
  geom_density()+
  theme_bw()

pop2corenotch$gene_name<-as.factor(pop2corenotch$gene_name)
pop2corenotch<-pop2corenotch[order(pop2corenotch$pos), ]

#plot
pop2corenotch$gene_name <-reorder(pop2corenotch$gene_name, rev(pop2corenotch$pos))

ggplot(pop2corenotch,aes(x=variable,y=gene_name,fill=value))+
  geom_tile(color="gray")+
  scale_fill_gradient2(high="red", low="white",mid="yellow", midpoint=5, limits=range(pop2corenotch$value),na.value="white")+
  facet_wrap(~type,scales="free")+
  ggtitle("CD45- kit+")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(size=5))



## Heatmap for genes notch1 core, with log2foldchange values of deseq2
ncol(deseq)
class(deseq)
deseq<-as.data.frame(deseq)

#log2fold
pop2deseq<-deseq[,c(2,33:35)]

pop2deseqnotch<-merge(pop2deseq,corenotch,by="gene_name",all.y=TRUE)
#create variable pval significative por foldchange
pop2deseqnotch$sig<-NA
pop2deseqnotch$`pvalue_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`<-as.numeric(pop2deseqnotch$`pvalue_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`)
pop2deseqnotch$sig <- ifelse(pop2deseqnotch$`pvalue_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-` <= 0.05, 
                        c("sig"), c("nosig")) 

pop2deseqnotch$gene_name<-as.factor(pop2deseqnotch$gene_name)
pop2deseqnotch<-pop2deseqnotch[order(pop2deseqnotch$`log2FoldChange_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`), ]
pop2deseqnotch<-melt(pop2deseqnotch,id.vars=c("gene_name","pos","sig"))
pop2deseqnotch$value<-as.numeric(pop2deseqnotch$value)

#plot
pop2deseqnotch$gene_name <-reorder(pop2deseqnotch$gene_name, rev(pop2deseqnotch$pos))

ggplot(subset(pop2deseqnotch,pop2deseqnotch$variable=="log2FoldChange_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-"),aes(x=variable,y=gene_name,fill=value))+
  geom_tile(color="gray")+
  scale_fill_gradient2(high="red", low="white",mid="yellow", midpoint=3, limits=range(pop2deseqnotch$value),na.value="white")+
  facet_wrap(~sig,scales="free")+
  ggtitle("CD45- kit+")+
  theme_minimal()+
  theme(axis.ticks.x = element_blank(),axis.text.y = element_text(size=8))


### Excel data Pierre
pathpierre<-read_xlsx('pierre/abigas_Tables 2.xlsx',sheet = 3)
pathpierre<-pathpierre[-34,]
names(pathpierre)[5]<-"gene_name"

#Extract expression levels of genes in pathpierre
pop2pierre<-merge(pop2core,pathpierre,by="gene_name")

#Extract log2foldchange from Deseq
pop2pierre<-merge(pop2deseq,pop2pierre,by="gene_name")
pop2pierre$gene_name <-reorder(pop2pierre$gene_name, pop2pierre$`log2FoldChange_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`)

ggplot(pop2pierre,aes(x=category,y=gene_name,fill=`log2FoldChange_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`))+
  geom_tile(color="gray")+
  scale_fill_gradient2(high="red", low="blue",mid="white", midpoint=0, limits=range(pop2pierre$`log2FoldChange_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`),na.value="white")+
  facet_wrap(~GO,scales="free",nrow=2)+
  ggtitle("CD45- kit+")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_text(size=12))





## GO terms and pathways
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "Chromosome_Location","GO_Biological_Process_2018" , "TF_Perturbations_Followed_by_Expression" ,"KEGG_2018")


#### For pop2, log2foldchange significative pe treatment (p no adjusted < 0.05)
pop2deseq$`pvalue_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`<-as.numeric(pop2deseq$`pvalue_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`)
pop2deseq_sig<-subset(pop2deseq,pop2deseq$`pvalue_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-`<=0.05)
#create variable direction expression (up or down in control)
pop2deseq_sig$exp<-NA
pop2deseq_sig$exp<-ifelse(pop2deseq_sig$`log2FoldChange_Ctrl_cKit+_CD45-_vs_DII4_cKit+_CD45-` > 0, 
       c("DOWN_DII4"), c("UP_DII4")) 

enrich_uptreat <- enrichr((subset(pop2deseq_sig,pop2deseq_sig$exp=="UP_DII4"))$gene_name, dbs)
enrichdowntreat <- enrichr((subset(pop2deseq_sig,pop2deseq_sig$exp=="DOWN_DII4"))$gene_name, dbs)

uptreat <- enrich_uptreat[["GO_Biological_Process_2018"]]
uptreat$exp<-c("UP_DII4")
downtreat <- enrichdowntreat[["GO_Biological_Process_2018"]]
downtreat$exp<-c("DOWN_DII4")

allupdown<-rbind(uptreat,downtreat)
allupdown<-subset(allupdown,allupdown$Adjusted.P.value<=0.05)
allupdown<- allupdown[order(allupdown$Adjusted.P.value),]

allupdown$Term <- factor(allupdown$Term, levels=(allupdown$Term)[rev(order(allupdown$Adjusted.P.value))])

ggplot(allupdown,aes(x=exp,y=Term,fill=Adjusted.P.value),alpha = 0.5)+
  geom_tile(color="black",width=1,height=0.8)+
  scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(allupdown$Adjusted.P.value),na.value="white")+
  #scale_fill_viridis()+
  facet_wrap(~exp,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=6),axis.title=element_text(size=14,face="bold"))



#PCA with only significative up and down by phenotype


#number of genes expressed
pop3cont$total<-rowSums(pop3cont)
pop3cont<-subset(pop3cont,pop3cont$total>1)
pop3treat$total<-rowSums(pop3treat)
pop3treat<-subset(pop3treat,pop3treat$total>1)

pop2rna$total<-rowSums(pop2rna)
pop2rna<-subset(pop2rna,pop2rna$total>1)
pop3rna$total<-rowSums(pop3rna)
pop3rna<-subset(pop3rna,pop3rna$total>1)

##subset treatment control each population
table(samples$Condition,samples$Population)

### PCA for each population
rnaex_t<-t(pop3rna)
rownames(rnaex_t)

pca <- prcomp(rnaex_t, scale=F)
summary(pca)

pca_htseq<-as.data.frame(pca$x[,c(1:2)])
pca_htseq$Sample<-row.names(pca_htseq)
pca_data<-merge(samples,pca_htseq,by="Sample")

pca_data$Population

pcaplot<-ggplot(pca_data,aes(x=PC1,y=PC2))+
  geom_point(aes(color=Condition),size=5)+
  scale_color_manual(values=c("dodgerblue","dodgerblue4"))+
  stat_ellipse(aes(color=Condition),level=0.8)+
  xlab(paste0("PC1: ",round(((summary(pca))$importance)[2]*100),"% variance")) +
  ylab(paste0("PC2: ",round(((summary(pca))$importance)[5]*100),"% variance"))+
  ggtitle(unique(pca_data$Population))+
  theme_bw()

boxp<-ggplot(pca_data, aes(y=PC1, x=Condition,label=Sample)) +
  geom_point(aes(color=Condition),size=6) +
  geom_boxplot(alpha=0)+
  scale_color_manual(values=c("dodgerblue","dodgerblue4"))+
  ylab(paste0("PC1: ",round(((summary(pca))$importance)[2]*100),"% variance")) +
  xlab("")+
  ggtitle(unique(pca_data$Population))+
  theme_bw()

grid.arrange(pcaplot,boxp)

t.test(pca_data$PC1~pca_data$Condition)

# adonis separately per population
rnaex_meta<-merge(samples,rnaex_t,by.x="Sample",by.y="row.names")
rnaex_gene<-rnaex_meta[,6:ncol(rnaex_meta)]
adonis(rnaex_gene ~ Condition, data = rnaex_meta)




##### Differentially expressed genes DESeq2 ###

#prepare metadata
sampleTable <- data.frame(sampleName = samples$Sample,
                          condition = samples$Condition,
                          population = samples$Population)
#remove samples not included in RNAseq output reads
sampleTable<-sampleTable[-2,]
#add row.names = colnames htseq table
row.names(sampleTable)<-sampleTable$sampleName

#make characters
sampleTable$condition<-as.character(sampleTable$condition)
sampleTable$sampleName<-as.character(sampleTable$sampleName)
sampleTable$population<-as.character(sampleTable$population)

countData<-as.matrix(htseqraw)
colnames(countData)

# Object DESeq from matrix
DESeq2Table <- DESeqDataSetFromMatrix(countData =  countData,
                                      colData = sampleTable,
                                      design = ~ condition+population)

rowData(DESeq2Table)
mcols(DESeq2Table)


## Filter data, keeping only rows that have at least 10 reads total
keep<-rowSums(counts(DESeq2Table))>=10
DESeq2Table<-DESeq2Table[keep,]

DESeq2Table$condition <- factor(DESeq2Table$condition, levels = c("Control","Ab"))

## Differential expression analysis
DESeq2Table<-DESeq(DESeq2Table)

res<-results(DESeq2Table)
res

#heatmap
select <- order(rowMeans(counts(DESeq2Table,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(DESeq2Table)[,c("condition","population")])
ntd<-normTransform(DESeq2Table)

pheatmap(assay(ntd)[select,],cluster_rows=FALSE,show_rownames=FALSE,cluster_cols=FALSE,annotation_col=df)


#add biomart Data
ensembl76 <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "description"),
            filter="ensembl_gene_id",
            values= rownames(DESeq2Table),
            mart=ensembl76 )

#quality control and normalization
GeneCounts <- counts(DESeq2Table)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)



######### scRNAseq C. Procheri ########

genes_up<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CPorch/scRNAseq_up.txt",header = FALSE)
colnames(genes_up)<-c("gene")
genes_up$state<-"UP"

genes_down<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CPorch/scRNAseq_down.txt",header = FALSE)
colnames(genes_down)<-c("gene")
genes_down$state<-"DOWN"


## genes RNAseq Pierre treatment vs. control CD45
DEG_DII4<-read_xlsx("/Volumes/grcmc/YGUILLEN/RNASeq_CPorch/DEG_genes_cd45.xlsx")
DEG_DII4<- DEG_DII4[order(DEG_DII4$pvalue),]
max(DEG_DII4$pvalue)

DEG_DII4_down<-subset(DEG_DII4,DEG_DII4$FoldChange<0)
DEG_DII4_up<-subset(DEG_DII4,DEG_DII4$FoldChange>0)

## GO terms and pathways
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018","GO_Molecular_Function_2018","GO_Biological_Process_2017")
#"NCI-Nature_2016"

enrichedUP <- enrichr(genes_up$gene, dbs)
enrichedDOWN <- enrichr(genes_down$gene,dbs)

UPenrich <- enrichedUP[["GO_Biological_Process_2017"]]
UPenrich$group<-c("up")

DOWNenrich <- enrichedDOWN[["GO_Biological_Process_2017"]]
DOWNenrich$group<-c("down")

allGO<-rbind(UPenrich ,DOWNenrich)
require(reshape2)
require(tidyr)
allGO<-separate(data = allGO, col = Overlap, into = c("counts", "pathway"), sep = "/")
allGO<-subset(allGO,allGO$counts>=5)

bpsub<-subset(allGO,allGO$Adjusted.P.value<=0.05)
bpsub<- bpsub[order(bpsub$P.value),]
#bpsub$Term<-gsub('_Homo.*','',bpsub$Term)

bpsub$Term<-as.factor(bpsub$Term)
class(bpsub$Combined.Score)

bpsub$Genes<-gsub(';','/',bpsub$Genes)
bpsub$Term<-gsub('.GO.*','',bpsub$Term)

bpsub$Term <- factor(bpsub$Term, levels=bpsub[order(bpsub$group, bpsub$Combined.Score), ]$Term)
ggplot(bpsub,aes(x=Term,y=Combined.Score,fill=group))+
  geom_bar(position="dodge",stat = "identity",color="black")+
  #geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Blues")+
  theme(axis.text.y =element_text(size=9),axis.title=element_text(size=14,face="bold"),axis.text.x = element_text(size=12))+
  coord_flip()


### DAVID enrichment 45neg kit positive
DAVID_RNAseq<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CPorch/pierre/DAVID_RNAseq_BP.txt",header = TRUE,sep="\t",dec = ",")
DAVID_RNAseq$X<-NULL

DAVID_scRNAseq<-read.delim("/Volumes/grcmc/YGUILLEN/RNASeq_CPorch/pierre/DAVID_scRNAseq_BP_MF.txt",header = TRUE,sep="\t",dec = ",")
DAVID_scRNAseq$X<-NULL

DAVID_RNAseq$Term<-as.factor(DAVID_RNAseq$Term)
DAVID_RNAseq$Term<-gsub('^GO.*:.*_','',DAVID_RNAseq$Term)
DAVID_RNAseq$Term <- factor(DAVID_RNAseq$Term, levels=DAVID_RNAseq[order(DAVID_RNAseq$State, DAVID_RNAseq$Fold.Enrichment), ]$Term)

ggplot(DAVID_RNAseq,aes(x=Term,y=Fold.Enrichment,fill=State))+
  geom_bar(position="dodge",stat = "identity",color="black")+
  #geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Blues")+
  ggtitle("Functional Enrichment of DEGs IgG vs. a-DII4")+
  theme(axis.text.y =element_text(size=12),axis.title=element_text(size=14,face="bold"),axis.text.x = element_text(size=12))+
  coord_flip()


DAVID_scRNAseq$Term<-as.factor(DAVID_scRNAseq$Term)
DAVID_scRNAseq$Term<-gsub('^GO.*:.*~','',DAVID_scRNAseq$Term)
DAVID_scRNAseq$Term <- factor(DAVID_scRNAseq$Term, levels=DAVID_scRNAseq[order(DAVID_scRNAseq$State, DAVID_scRNAseq$Fold.Enrichment), ]$Term)

ggplot(DAVID_scRNAseq,aes(x=Term,y=Fold.Enrichment,fill=State))+
  geom_bar(position="dodge",stat = "identity",color="black")+
  #geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Blues")+
  ggtitle("Functional Enrichment of DEGs IgG vs. a-DII4")+
  theme(axis.text.y =element_text(size=12),axis.title=element_text(size=14,face="bold"),axis.text.x = element_text(size=12))+
  coord_flip()
