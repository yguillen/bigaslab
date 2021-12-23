### YGUILLEN SCRIPT TO MANAGE CLINICAL DATA FROM TCell project (EGA AND SRA toolkit data)
# 05 JUNE 2019
# UPDATED MAY 2021

## iNstall extra packages
BiocManager::install("topGO")
BiocManager::install("heatmaply")
BiocManager::install("readxl")
BiocManager::install("ggrepel")
BiocManager::install(c("GO.db","ggnewscale","topGO","gridExtra","Hmisc","corrplot","gRbase","pheatmap","remotes"))
BiocManager::install("biomaRt")
BiocManager::install("scatterpie")
BiocManager::install("ggpubr")
remotes::install_github("grst/immunedeconv")

#basic libraries
library(data.table)
library(ggplot2)
library(nlme)
library(readxl)
library(ggrepel)
library(ggnewscale)
library(GO.db)
library(topGO)
library(gridExtra)
library(Hmisc)
library(corrplot)
library(DESeq2)
library(gplots)
library(genefilter)
library(gRbase)
library(plyr)
library(pheatmap)
library(heatmaply)
library(limma)
library(remotes)
library(dplyr)
library(scatterpie)
library(biomaRt)

set.seed(100)

# Import datasets
##### 
setwd("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/") 
#setwd("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Public_data/") 

# IntoGen cancer driver genes
candriver<-read.delim("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Drivers_intoGens/Compendium_Cancer_Genes.tsv",header = TRUE)
#candriver<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Drivers_intoGens/Compendium_Cancer_Genes.tsv",header = TRUE)
unique(candriver$SYMBOL)

# metadata for EGA
#progega<-read_xlsx("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Public_data/EGA/EGAD00001000849/metadata_EGA.xlsx")
#row.names(progega)<-progega$Sample
#colnames(progega)<-c("SampID","Prog")

## Structural and mutational data from TARGET
somatic_target<-read.delim("/Users/yguillen/Documents/Projects/TARGET/stjude.org_TARGET_TALL_WXS_Diagnosis_IlluminaHiSeq_somatic.maf.txt")
structural_target<-read.delim("/Users/yguillen/Documents/Projects/TARGET/TARGET_TALL_mRNA-seq_sv.tsv")

# Selecting RBP genes form different sources
Huang_RBP<-as.data.frame(read_xlsx('/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/RBP_gene_lists/Huang_et_2018_PNAS.xlsx',sheet = "Sheet1"))
Wang_RBP<-as.data.frame(read_xlsx('/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/RBP_gene_lists/Wang_Aifantis_et_2019_CancerCell.xlsx',sheet = "Sheet1"))
Sebest_RBP<-as.data.frame(read_xlsx('/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/RBP_gene_lists/Sebestyen_et_2016_GenRes.xlsx',sheet = "Sheet1"))
rMATs_RBP<-read.delim('/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/RBP_gene_lists/rMAPS2_list.txt')


# Genes id and ENS id
geneid<-read.delim('/Volumes/cancer/TCell/vast_output/genes_ENS.txt',header = TRUE,sep=' ')
#geneid<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/genes_ENS.txt',header = TRUE,sep=' ')

## For all samples
metadata<-read.delim("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/metadata_all.txt",sep="\t")
#metadata<-read_xlsx("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Public_data/metadata_all.xlsx")

metadata_target<-as.data.frame(read_xlsx("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/metadata_gender_TARGET.xlsx"))
row.names(metadata_target)<-metadata_target$GSM
metadata_target$GSM<-NULL

#####

# Metadata summary
#####

metadata_sub<-subset(metadata,select=c(SampID,total_reads))
metadata_sub<-aggregate(x = metadata_sub$total_reads, by = list(metadata_sub$SampID), FUN = "sum")
metadata_sub$thres<-ifelse(metadata_sub$x>=65000000,c("validvast"),c("novalidvast"))

metadata<-merge(metadata,metadata_sub,by.x="SampID",by.y="Group.1")

metadata$Event_free_surv_days<-as.numeric(metadata$Event_free_surv_days)

metadata$Reads_pair1<-NULL
metadata$total_reads<-NULL
metadata<-metadata[!duplicated(metadata$SampID),]
names(metadata)[16]<-"Reads"

ggplot(metadata,aes(x=SampID,y=Reads))+
  geom_segment(aes(yend = Reads, x = SampID, xend = SampID, y = 0),color="grey",size=3)+
  geom_segment(data=metadata[metadata$First_event!="NA",],aes(yend = Reads, x = SampID, xend = SampID, y = 0,color=First_event),size=2)+
  geom_point(aes(y=Reads+2000000,shape=Stage),alpha=0.5,size=4)+
  scale_shape_manual(values=c(13, 1,9,8))+
  #geom_jitter(aes(color=Pheno,shape=thres),size=3)+
  #geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_hline(yintercept=65000000)+
  scale_color_brewer(palette = "Set3")+
  facet_grid(~GSEA,scales="free",space="free")+
  theme_light()+
  #ylim(c(0,300000000))+
  xlab("")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, hjust = 1),
        plot.title=element_text(family="Tahoma",angle=45,size=5),
        text=element_text(family="Tahoma"),
        legend.position = "bottom")

table(metadata$First_event,metadata$Vital_status)

#####

# CIBERSORT adding data
#####

## Expression profile CIBERSORT from above vastout
#cibersort<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CIBERSORT/CIBERSORT.Output_Job5.txt",header = TRUE)
cibersort<-read.delim("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/CIBERSORT/CIBERSORT.Output_Job5.txt",header = TRUE)
row.names(cibersort)<-cibersort$Input.Sample
colnames(cibersort)[1]<-"SampID"

ciberbar<-cibersort[,c(1:(ncol(cibersort)-3))]
corciber<-merge(ciberbar,metadata,by="SampID")

corciber$Blasts<-as.numeric(corciber$Blasts)
cor.test(corciber$Blasts,corciber$Thy1)
plot(corciber$Blasts,corciber$Thy1)

ciberbar<-melt(ciberbar,id.vars="SampID")

ciberbar<-merge(metadata,ciberbar,by="SampID")

ciberbar$SampID <- factor(ciberbar$SampID, levels=unique(ciberbar[rev(order(ciberbar$value)), ]$SampID))

sampord<-ciberbar[ciberbar$variable=="Thy2",][rev(order(ciberbar[ciberbar$variable=="Thy2",]$value)),]$SampID
sampord
ciberbar$SampID <- factor(ciberbar$SampID, levels=sampord)

ciberbar$Blasts<-as.numeric(ciberbar$Blasts)

ggplot(ciberbar, aes(x=as.factor(SampID), y=value,fill=variable))+
  geom_bar(stat="identity", fill="black",width = 1) +
  geom_bar(stat="identity", aes(fill=variable),width = 0.8) +
  scale_fill_brewer(palette="Paired")+
  facet_grid(~First_event,scales="free",space="free")+
  #geom_point(aes(y=genes+500,color=Reads),size=5)+
  new_scale_color()+
  geom_point(data=ciberbar[!is.na(ciberbar$First_event) & ciberbar$variable=="Thy2",],aes(y=1.3),color="black",size=3,shape=17)+
  geom_point(data=ciberbar[!is.na(ciberbar$First_event) & ciberbar$variable=="Thy2",],aes(y=1.3,color=First_event),size=2,shape=17)+
  scale_color_brewer(palette="Set1")+
  new_scale_color()+
  geom_point(data=ciberbar[ciberbar$variable=="Thy2",],aes(y=1.1),color="black",size=3)+
  geom_point(data=ciberbar[ciberbar$GSEA!="line" & ciberbar$variable=="Thy2",],aes(y=1.1,color=Pheno),size=2)+
  scale_color_brewer(palette="Set3")+
  new_scale_color()+
  geom_point(data=ciberbar[ciberbar$variable=="Thy2",],aes(y=1.2),color="black",size=3)+
  geom_point(data=ciberbar[ciberbar$GSEA!="line" & ciberbar$variable=="Thy2",],aes(y=1.2,color=GSEA),size=2)+
  scale_color_brewer(palette="Dark2")+
  geom_point(data=ciberbar[ciberbar$GSEA!="line" & ciberbar$variable=="Thy2" & !is.na(ciberbar$Blasts),],aes(y=1.4,size=Blasts/10),color="grey")+
#  scale_colour_gradient(low = "white", high = "black")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom")


# Summary of variable by group
tapply(ciberbar[!is.na(ciberbar$Blasts),]$Blasts, ciberbar[!is.na(ciberbar$Blasts),]$First_event, summary)

#####


## Read expressed genes VAST-TOOLS OUTPUT
#####
#setwd("/Volumes/cancer/TCell/vast_output/") 
setwd("/Users/yguillen/Desktop/temp/TCell/vast_output/") 


#Total expressed genes in primary samples
exp_genes<-read.delim('expr_out/total_expressed_genes.txt',header = FALSE)

#Total expressed genes in cell lines
#exp_genes_lines<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/Gekas_RNAseq/vast_output/expr_out/total_expressed_genes.txt',header = FALSE)
exp_genes_lines<-read.delim('/Users/yguillen/Desktop/temp/Gekas_RNASeq/vast_output/expr_out/total_expressed_genes.txt',header = FALSE)

exp_genes<-rbind(exp_genes,exp_genes_lines)

colnames(exp_genes)<-c("SampID","genes")

# matching metadata and expression matrix, check names
metadata$SampID
exp_genes$SampID<-gsub('EGAR[0-9]*_','',exp_genes$SampID)
exp_genes$SampID<-gsub('.sra_1','',exp_genes$SampID)
exp_genes$SampID<-gsub('_1','',exp_genes$SampID)
setdiff(metadata$SampID,exp_genes$SampID)

# merge
exp_genes<-merge(exp_genes,metadata,by="SampID")
exp_genes$Source<-as.factor(exp_genes$Source)
levels(exp_genes$Source)

# ORDER samples by genes
exp_genes$SampID <- factor(exp_genes$SampID, levels = exp_genes$SampID[order(exp_genes$genes)])



ggplot(exp_genes, aes(x=as.factor(SampID), y=genes,fill=GSEA)) +# Note that id is a factor. If x is numeric, there is some space between the first bar
  # This add the bars with a blue color
  geom_bar(stat="identity", fill="black",width = 1) +
  geom_bar(stat="identity", aes(fill=GSEA),width = 0.8) +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  #ylim(1000,20000) +
  scale_fill_brewer(palette="Paired")+
  new_scale_color()+
  geom_point(aes(y=genes+500,color=Reads/1000000),size=2)+
  geom_point(data=subset(exp_genes,exp_genes$genes>18000),aes(y=genes+1000),size=3,shape=8,color="red")+
  scale_colour_gradient(low = "white", high = "black")+
  theme_bw() +
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0)+
  theme(axis.text.x = element_text(angle=0, hjust = 1),
        strip.text = element_text(size=6),
        legend.position = "bottom")+
  theme(axis.text.x = element_blank())


exp_genes$SampID <- factor(exp_genes$SampID, levels = exp_genes$SampID[order(exp_genes$genes)])
exp_genes$Type<-factor(exp_genes$Type,levels = c("Thymus","TALL","TLBL","BM"))


r1<-ggplot(exp_genes, aes(x=Type, y=Reads)) +# Note that id is a factor. If x is numeric, there is some space between the first bar
  # This add the bars with a blue color
  geom_point(aes(color=Type),size=4)+
  scale_color_manual(values=c("grey","coral","coral","grey"))+
  new_scale_color()+
  geom_violin(aes(color=Type),alpha=0.3,width = 0.8) +
  scale_color_manual(values=c("grey","coral","coral","grey"))+
  new_scale_color()+
  stat_summary(fun.y=median, geom="point",size=7)  + 
  stat_summary(aes(color=Type),fun.y=median, geom="point",size=4)  + 
  scale_color_manual(values=c("grey","coral","coral","grey"))+
  facet_grid(~GSEA,scales="free_x",space = "free")+
  theme_bw()

r2<-ggplot(exp_genes, aes(x=Type, y=genes)) +# Note that id is a factor. If x is numeric, there is some space between the first bar
  # This add the bars with a blue color
  geom_point(aes(color=Type),size=4)+
  scale_color_manual(values=c("grey","coral","coral","grey"))+
  new_scale_color()+
  geom_violin(aes(color=Type),alpha=0.3,width = 0.8) +
  scale_color_manual(values=c("grey","coral","coral","grey"))+
  new_scale_color()+
  stat_summary(fun.y=median, geom="point",size=7)  + 
  stat_summary(aes(color=Type),fun.y=median, geom="point",size=4)  + 
  scale_color_manual(values=c("grey","coral","coral","grey"))+
  facet_grid(~GSEA,scales="free_x",space = "free")+
  theme_bw()


## Number of expressed genes and reads per sample
grid.arrange(r1,r2,ncol=1)


ggplot(exp_genes, aes(x=Reads, y=genes)) +# Note that id is a factor. If x is numeric, there is some space between the first bar
  # This add the bars with a blue color
  geom_point(aes(color=Type),size=5)+
  geom_smooth(aes(color=GSEA),color="black",method="lm")+
  scale_color_manual(values=c("grey","coral","coral","grey"))+
  facet_wrap(~GSEA,scales="free_x")+
  theme_bw()


## Correlation # Reads and # total genes expressed? per cohort

wilcox.test((exp_genes[exp_genes$GSEA=="GSE57982",])$Reads~(exp_genes[exp_genes$GSEA=="GSE57982",])$Type)
wilcox.test((exp_genes[exp_genes$GSEA=="GSE57982",])$genes~(exp_genes[exp_genes$GSEA=="GSE57982",])$Type)

wilcox.test((exp_genes[exp_genes$GSEA=="GSE109231",])$Reads~(exp_genes[exp_genes$GSEA=="GSE109231",])$Type)
wilcox.test((exp_genes[exp_genes$GSEA=="GSE109231",])$genes~(exp_genes[exp_genes$GSEA=="GSE109231",])$Type)

wilcox.test((exp_genes[exp_genes$GSEA=="GSE57982",])$Reads~(exp_genes[exp_genes$GSEA=="GSE57982",])$Type)
wilcox.test((exp_genes[exp_genes$GSEA=="GSE57982",])$genes~(exp_genes[exp_genes$GSEA=="GSE57982",])$Type)

wilcox.test((exp_genes[exp_genes$GSEA=="EGA_TLE",])$Reads~(exp_genes[exp_genes$GSEA=="EGA_TLE",])$Type)
wilcox.test((exp_genes[exp_genes$GSEA=="EGA_TLE",])$genes~(exp_genes[exp_genes$GSEA=="EGA_TLE",])$Type)


# correlation expressed genes and cohort
kruskal.test(genes~as.factor(GSEA),data=exp_genes)

# cor.test expressed genes and reads per cohort

cor.test((exp_genes[exp_genes$GSEA=="EGA_TLE",])$Reads,(exp_genes[exp_genes$GSEA=="EGA_TLE",])$genes)
cor.test((exp_genes[exp_genes$GSEA=="EGA_R",])$Reads,(exp_genes[exp_genes$GSEA=="EGA_R",])$genes)
cor.test((exp_genes[exp_genes$GSEA=="GSE109231",])$Reads,(exp_genes[exp_genes$GSEA=="GSE109231",])$genes)
cor.test((exp_genes[exp_genes$GSEA=="GSE110633",])$Reads,(exp_genes[exp_genes$GSEA=="GSE110633",])$genes)
cor.test((exp_genes[exp_genes$GSEA=="GSE57982",])$Reads,(exp_genes[exp_genes$GSEA=="GSE57982",])$genes)
cor.test((exp_genes[exp_genes$GSEA=="GSE69239",])$Reads,(exp_genes[exp_genes$GSEA=="GSE69239",])$genes)
cor.test((exp_genes[exp_genes$GSEA=="line",])$Reads,(exp_genes[exp_genes$GSEA=="line",])$genes)
cor.test((exp_genes[exp_genes$GSEA=="TARGET",])$Reads,(exp_genes[exp_genes$GSEA=="TARGET",])$genes)

#####


# Upload cRPKM and counts table from VASTOOLS and annotate genes
#####
###### UPLOAD COUNTS TABLE HERE cRPKM_AND_COUNTS-HsaXXX.tab
# Input counts (seq(4..)) or RPKM (seq(3..)) (normalized) T-ALL RNASeq
tallcount<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/cRPKM_AND_COUNTS-Hsa389.tab")
#tallcount<-read.delim("/Volumes/cancer/TCell/vast_output/cRPKM_AND_COUNTS-Hsa264.tab")

tallcount<-tallcount[,c(1,2,seq(4, ncol(tallcount), by=2))]

colnames(tallcount)<-gsub('EGAR[0-9]*_','',colnames(tallcount))
colnames(tallcount)<-gsub('.sra_1','',colnames(tallcount))
colnames(tallcount)<-gsub('_1','',colnames(tallcount))
colnames(tallcount)<-gsub('.Counts','',colnames(tallcount))
#colnames(tallcount)<-gsub('.cRPKM','',colnames(tallcount))
colnames(tallcount)


# Input counts cell lines shbcat
linecount<-read.delim("/Users/yguillen/Desktop/temp/Gekas_RNASeq/vast_output/cRPKM_AND_COUNTS-Hsa4.tab")
linecount<-linecount[,c(1,2,seq(4, ncol(linecount), by=2))]
colnames(linecount)<-gsub('.Counts','',colnames(linecount))
#colnames(linecount)<-gsub('.cRPKM','',colnames(linecount))
colnames(linecount)

dim(linecount)
dim(tallcount)

## Distribution of counts, from vasttools counts table line 1106
colnames(tallcount)

#subset random 1000 genes
tallcount_melt<-sample_n(tallcount[,3:ncol(tallcount)],1000)

tallcount_melt<-melt(tallcount_melt)
colnames(tallcount_melt)<-c("SampID","value")
tallcount_melt<-merge(metadata,tallcount_melt,by="SampID")


ggplot(tallcount_melt[!is.na(tallcount_melt$value),],aes(x=value))+
  geom_density(aes(color=GSEA),adjust=1)+
  scale_color_brewer(palette="Paired")+
  #  ylim(c(0,5))+
  theme_light()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")+
  xlim(c(0,1000))



# only mitochondrial genes
tallcount_melt<-tallcount[grepl("MT-",tallcount$NAME),]
tallcount_melt<-tallcount_melt[,3:ncol(tallcount)]

tallcount_melt<-melt(tallcount_melt)
colnames(tallcount_melt)<-c("SampID","value")
tallcount_melt<-merge(metadata,tallcount_melt,by="SampID")


ggplot(tallcount_melt,aes(x=value))+
  geom_density(aes(color=GSEA),adjust=1)+
  scale_color_brewer(palette="Paired")+
  #  ylim(c(0,5))+
  theme_light()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")+
  xlim(c(0,500000))


## Annotation of genes
#genome version hg19
library(biomaRt)

#version 37
hg19<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

#latest version
hg<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

# Using other mirrors if slow...
#hg <- useEnsembl(biomart = "ensembl", 
#                   dataset = "hsapiens_gene_ensembl", 
#                   mirror = "asia")

bm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
            values= geneid$Gene,
            filter="ensembl_gene_id",
            mart=hg) 

length(unique(bm$ensembl_gene_id))
length(unique(geneid$Gene))

## genes not found
length(setdiff(unique(geneid$Gene),unique(bm$ensembl_gene_id)))
#689

## GO TERM table all genes 
vas<-toTable(GOTERM)
colnames(vas)[2]<-"repgo_id"
vascc<-vas[vas$Ontology=="CC",]
vasbp<-vas[vas$Ontology=="BP",]
vasmf<-vas[vas$Ontology=="MF",]

geneGO<-merge(geneid,bm,by.x="NAME",by.y="external_gene_name",all.x=TRUE)

#empty values as NA
geneGO$go_id<-gsub("^$", NA, geneGO$go_id)

# all terms
#geneGO<-merge(geneGO,vas,by="go_id",all.x=TRUE)

#only subcategories terms
geneGOCC<-merge(geneGO,vascc,by="go_id",all.x=TRUE)
geneGOBP<-merge(geneGO,vasbp,by="go_id",all.x=TRUE)
geneGOMF<-merge(geneGO,vasmf,by="go_id",all.x=TRUE)


## Principal Component Analysis ####### (caution with  NAs)

#vastout<-read.delim("INCLUSION_filter.tab")

# cRPKM with merged replicataes already
#vastout<-read.delim("/Volumes/cancer/TCell/vast_output/cRPKM-Hsa264.tab")
vastout<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/cRPKM-Hsa389.tab")


# cRPKM from RPMI RNASeq Christos
#cellines<-read.delim("/Volumes/cancer/Gekas_RNAseq/vast_output/cRPKM-Hsa4.tab")
cellines<-read.delim("/Users/yguillen/Desktop/temp/Gekas_RNAseq/vast_output/cRPKM-Hsa4.tab")

vastout<-merge(cellines,vastout,by="ID")
vastout$NAME.y<-NULL
colnames(vastout)[2]<-"NAME"


colnames(vastout)<-gsub('_1','',colnames(vastout))
colnames(vastout)<-gsub('EGAR.*_','',colnames(vastout))
colnames(vastout)<-gsub('.sra','',colnames(vastout))

row.names(vastout)<-vastout$ID
vastout$ID=NULL
vastout$NAME=NULL

rownames<-colnames(vastout)

#####

## OUTPUT FOR EPIC AND CIBERSORT
#####
### Export vastout with gene names for EPIC and CIBERSORT deconvolution analysis
row.names(vastout)
epic<-vastout
epic$geneName<-row.names(epic)

dim(epic)

epic<-merge(epic,geneid,by.x="geneName",by.y="Gene")
epic$geneName<-NULL
epic<-epic[,c(ncol(epic),1:(ncol(epic))-1)]
dim(epic)

write.table(epic,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/EPIC/TALL_EPIC_exp.tsv",row.names = FALSE,quote = FALSE,sep="\t")
#write.table(epic,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/EPIC/TALL_EPIC_exp.tsv",row.names = FALSE,quote = FALSE,sep="\t")

#only T-Cell types for cibersort
epic<-vastout
epic$geneName<-row.names(epic)

epic$geneName<-NULL
epic<-epic[,-c(5:10,55:56,392:393)]
dim(epic)

# CAREFULL WITH DIMENSIONS! OLD AN NEW
write.table(epic,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CIBERSORT/mixture_tall.tab",row.names = TRUE,quote = FALSE,sep="\t")
#write.table(epic,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/CIBERSORT/mixture_tall.tab",row.names = TRUE,quote = FALSE,sep="\t")

#####

# PCA with expression data and 3 approaches for batch correction

#####

#fOR PCA complete.cases
dim(vastout)
vastout<-vastout[complete.cases(vastout), ]

vastout_t<-as.data.frame(t(vastout))
colnames(vastout_t)
row.names(vastout_t)

DSmat <- data.frame(sapply(vastout_t, function(x) as.numeric(as.character(x))))
row.names(DSmat)<-row.names(vastout_t)


#select those with variance != 0 (error in prcomp otherwise)
cond<-(apply(DSmat, 2, var)!=0)
DSmatpca<-DSmat[, cond, drop = FALSE]


#PCA all data
res.pca<-prcomp(DSmatpca,scale=TRUE)
(res.pca$sdev^2/sum(res.pca$sdev^2))[1:5]


#Proportion of variance

pcatab<-as.data.frame(res.pca$x)
pcatab$names<-row.names(pcatab)
row.names(pcatab)
colnames(pcatab)

pcatab<-pcatab[,c(1:3)]
row.names(pcatab)

metadata$SampID

setdiff(metadata$SampID,row.names(pcatab))

# FOR ALL DATA
pcatab<-merge(metadata,pcatab,by.x="SampID",by.y="row.names",all=TRUE)

pcatab$var<-paste(pcatab$Source,pcatab$Pheno,sep="_")

library(ggnewscale)

pcatab$GSEA<-as.factor(pcatab$GSEA)


# Detecting outliers, estimate z-score
outDS<-data.frame(SampID=pcatab$SampID,
                  zscore_PC1=abs(scale(pcatab$PC1, center = TRUE, scale = TRUE)),
                  zscore_PC2=abs(scale(pcatab$PC2, center = TRUE, scale = TRUE)))

pcatab<-merge(pcatab,outDS,by="SampID")


ggplot(pcatab, aes(PC1, PC2,label=SampID)) +
  #geom_point(data=pcatab[pcatab$source=="BM",],size=10,color="black") +
  #geom_point(data=pcatab[pcatab$source=="Thymus",],size=10,color="grey") +
  #geom_point(data=pcatab[pcatab$source=="TALL",],size=10,color="red") +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=GSEA)) +
  #geom_point(data=pcatab[pcatab$GSEA=="GSE57982",],size=6,color="black") +
  scale_color_brewer(palette = "Paired")+
  new_scale_color()+
  stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
  scale_color_brewer(palette = "Paired")+
  #geom_label_repel(data=pcatab[pcatab$GSEA=="GSE69239",],aes(label=Source), size=4,fontface = "italic",alpha=0.7)+
  #geom_label_repel(data=pcatab[pcatab$Source=="Thymus" ,],aes(label=Source), size=4,fontface = "italic",alpha=0.7,color="grey")+
  #geom_label_repel(data=pcatab[pcatab$GSEA=="EGA_TLE" ,],aes(label=SampID), size=4,fontface = "italic",alpha=0.7,color="grey")+
  geom_label_repel(data=pcatab[pcatab$zscore_PC1>=2.5 | pcatab$zscore_PC2>=2.5,],aes(label=SampID), size=3,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  #ylim(-50,50)+
  ylab("PC2 (10.86% variance)")+
  xlab("PC1 (17.64% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()


ggplot(pcatab, aes(y=PC1, x=GSEA)) +
  geom_point(aes(color=GSEA),size=8) +
  geom_boxplot(alpha=0)+
  #stat_ellipse(aes(color=condition),level=0.8)+
  scale_color_brewer(palette="Paired")+
  #geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  #ylab(paste0("PC1: ",percentVar[1],"% variance")) +
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



# PCA with no EGA_TLE

#DSmat_noTLE <- data.frame(sapply(vastout_t, function(x) as.numeric(as.character(x))))
#row.names(DSmat_noTLE)<-row.names(vastout_t)
DSmat_noTLE<-vastout_t[!grepl('TLE',row.names(vastout_t)) & !grepl('Thymus',row.names(vastout_t)) & !grepl('R7',row.names(vastout_t)) & !grepl('R8',row.names(vastout_t)) & !grepl('R9',row.names(vastout_t)) & !grepl('R10',row.names(vastout_t)),]
row.names(DSmat_noTLE)

#select those with variance != 0 (error in prcomp otherwise)
cond<-(apply(DSmat_noTLE, 2, var)!=0)
DSmatpca_noTLE<-DSmat_noTLE[, cond, drop = FALSE]


res.pcaTLE<-prcomp(DSmatpca_noTLE,scale=TRUE)
(res.pcaTLE$sdev^2/sum(res.pcaTLE$sdev^2))[1:5]


#Proportion of variance

pcatab_TLE<-as.data.frame(res.pcaTLE$x)
pcatab_TLE$names<-row.names(pcatab_TLE)
row.names(pcatab_TLE)

pcatab_TLE<-pcatab_TLE[,c(1:3)]
row.names(pcatab_TLE)

metadata$SampID

setdiff(metadata$SampID,row.names(pcatab_TLE))

pcatab_TLE<-merge(metadata,pcatab_TLE,by.x="SampID",by.y="row.names",all.y=TRUE)
pcatab_TLE$var<-paste(pcatab_TLE$Source,pcatab_TLE$Pheno,sep="_")

# Detecting outliers, estimate z-score
outDSTLE<-data.frame(SampID=pcatab_TLE$SampID,
                  zscore_PC1=abs(scale(pcatab_TLE$PC1, center = TRUE, scale = TRUE)),
                  zscore_PC2=abs(scale(pcatab_TLE$PC2, center = TRUE, scale = TRUE)))

pcatab_TLE<-merge(pcatab_TLE,outDSTLE,by="SampID")
dim(pcatab_TLE)

library(ggnewscale)

pcatab_TLE$GSEA<-as.factor(pcatab_TLE$GSEA)

ggplot(pcatab_TLE, aes(PC1, PC2,label=SampID)) +
  #geom_point(data=pcatab[pcatab$source=="BM",],size=10,color="black") +
  #geom_point(data=pcatab[pcatab$source=="Thymus",],size=10,color="grey") +
  #geom_point(data=pcatab[pcatab$source=="TALL",],size=10,color="red") +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=GSEA)) +
  #geom_point(data=pcatab[pcatab$GSEA=="GSE57982",],size=6,color="black") +
  scale_color_brewer(palette = "Paired")+
  new_scale_color()+
  stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
  scale_color_brewer(palette = "Paired")+
  #scale_linetype_manual(values=c("solid","dotted"))+
  #geom_label_repel(data=pcatab_TLE[pcatab_TLE$GSEA=="GSE69239",],aes(label=Source), size=4,fontface = "italic",alpha=0.7)+
  #geom_label_repel(data=pcatab_TLE[pcatab_TLE$Source=="Thymus" ,],aes(label=Source), size=4,fontface = "italic",alpha=0.7,color="grey")+
  geom_label_repel(data=pcatab_TLE[pcatab_TLE$zscore_PC1>=2.5 | pcatab_TLE$zscore_PC2>=2.5,],aes(label=SampID), size=3,fontface = "italic",alpha=0.7,color="black")+
  #scale_color_manual(values=c("purple","red"))+
  #ylim(-50,50)+
  ylab("PC2 (6.86% variance)")+
  xlab("PC1 (12.07% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()


library(plotly)

#3D plots
p <- plot_ly(x = pcatab$PC1, y = pcatab$PC2, z = pcatab$PC3, color = ~pcatab$GSEA,marker=list(size=5)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))


## BATCH EFFECT
# TEST EXCLUDING TLE

### APPROACH 1
# dataset all samples
vastout_t

# dataset with no TLE
vastout_t_TLE<-vastout_t[!grepl('TLE',row.names(vastout_t)) & !grepl('Thymus',row.names(vastout_t)) & !grepl('R7',row.names(vastout_t)) & !grepl('R8',row.names(vastout_t)) & !grepl('R9',row.names(vastout_t)) & !grepl('R10',row.names(vastout_t)),]

mandist<-dist(vastout_t_TLE)
heatmap(as.matrix(mandist))


require(MASS)

## Non-metric multidimensional scaling: sammon's non-linear mapping
sam1<-sammon(mandist,trace=FALSE)
sampos<-as.data.frame(sam1$points)
colnames(sampos)<-c("sam1","sam2")
sampos$SampID<-row.names(sampos)
sampos$SampID<-sub('EGAR.*_','',sampos$SampID)
rm(mandist)

pcatab_TLE<-merge(pcatab_TLE,sampos,by="SampID")


ggplot(pcatab_TLE, aes(sam1,sam2,label=SampID)) +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=GSEA)) +
  scale_color_brewer(palette = "Paired")+
  new_scale_color()+
  stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
  scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_TLE[pcatab_TLE$GSEA=="GSE69239",],aes(label=Source), size=4,fontface = "italic",alpha=0.7)+
  geom_label_repel(data=pcatab_TLE[pcatab_TLE$Source=="Thymus" ,],aes(label=Source), size=4,fontface = "italic",alpha=0.7,color="grey")+
 # geom_label_repel(data=pcatab[pcatab$GSEA=="EGA",],aes(label=SampID), size=4,fontface = "italic",alpha=0.7,color="grey")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()


## 2nd approach Quantification of batch effect
# all datasets
DSmat<-merge(metadata,DSmat,by.x="SampID",by.y="row.names",sort=FALSE)

row.names(DSmat)<-DSmat$SampID

#without TLE
DSmat_noTLE<-merge(metadata,DSmat_noTLE,by.x="SampID",by.y="row.names",sort=FALSE)

row.names(DSmat_noTLE)<-DSmat_noTLE$SampID

# ANOVA test
#table with pvalues and coefficients of GSEA variable and phenotype for each gene

# CHECK COLUMN 11, IN OLD SCRIPT IS 18, it might be other
# noTLE
dataaov<-sapply(DSmat_noTLE[,18:ncol(DSmat_noTLE)], function(y) summary(aov(y~DSmat_noTLE$GSEA+DSmat_noTLE$Source)))

#all datasets
dataaov<-sapply(DSmat[,18:ncol(DSmat)], function(y) summary(aov(y~DSmat$GSEA+DSmat$Source)))

pbatch<-as.data.frame(do.call(rbind,lapply(dataaov, function(v){v$'Pr(>F)'[1:2]})))
colnames(pbatch)<-c("P_GSEA","P_Source")
Fbatch<-as.data.frame(do.call(rbind,lapply(dataaov, function(v){v$'F value'[1:2]})))
colnames(Fbatch)<-c("F_GSEA","F_Source")

dim(pbatch)
dim(Fbatch)
pf_batch<-merge(pbatch,Fbatch,by="row.names")
pf_batch<-merge(geneid,pf_batch,by.x="Gene",by.y="Row.names")

pf_batch$imp<-ifelse(pf_batch$P_GSEA < pf_batch$P_Source, 
                        c("GSEA"), c("ASource"))  

pf_batch<-pf_batch[with(pf_batch,order(pf_batch$P_Source,pf_batch$imp)),]

## linear model with dummy variables, only for GSEA variable
DSmat_noTLE$GSEA<-as.factor(DSmat_noTLE$GSEA)

X <- model.matrix(~DSmat_noTLE$GSEA + DSmat_noTLE$Pheno)
dim(X)

geneexp<-DSmat_noTLE[,18:ncol(DSmat_noTLE)]

res <- t( sapply(1:ncol(geneexp),function(j){
  y <- geneexp[,j]
  fit <- lm(y~X-1)
  summary(fit)$coef[2,c(1,4)]
} ) )

res<-as.data.frame(res)
row.names(res)<-colnames(geneexp)
names(res) <- c("Estimate","p.value")

ggplot(res,aes(x=Estimate,y=p.value))+
  geom_point()+
  theme_bw()


res<-merge(res,geneid,by.x="row.names",by.y="Gene")
res<-res[with(res,order(res$p.value)),]

rm(dataaov)
rm(geneexp)

# 3rd approach BATCH effect with limma

# https://rdrr.io/bioc/limma/man/removeBatchEffect.html
# Extracting the batch effect (GSEA variable) and generating a dataframe with effect corrected.
library(limma)

# Excluding TLE from vastout dataset
vastout_TLE<-vastout[,!grepl('TLE',colnames(vastout)) & !grepl('Thymus',colnames(vastout)) & !grepl('R7',colnames(vastout)) & !grepl('R8',colnames(vastout)) & !grepl('R9',colnames(vastout)) & !grepl('R10',colnames(vastout))]
y<-vastout_TLE

# Not excluding EGA_TLE
#y<-vastout

row.names(y)
colnames(y)

namey<-as.data.frame(colnames(y))
colnames(namey)<-"SampID"
namey<-merge(namey,metadata,by="SampID",sort=FALSE)


y2 <- removeBatchEffect(y, namey$GSEA)
rm(y)
par(mfrow=c(1,2))
boxplot(as.data.frame(y[,c(8:10,20:22,47:49)]),main="Original")
boxplot(as.data.frame(y2[,c(8:10,20:22,47:49)]),main="Batch corrected")

# PCA with corrected batch effect
ty2<-t(y2)
rm(y2)
cond<-(apply(ty2, 2, var)!=0)
ty2<-ty2[, cond, drop = FALSE]

#PCA
res.pca<-prcomp(ty2,scale=TRUE)
summary(res.pca)

(res.pca$sdev^2/sum(res.pca$sdev^2))[1:5]

#Proportion of variance

pcatab_batch<-as.data.frame(res.pca$x)
pcatab_batch$names<-row.names(pcatab_batch)
row.names(pcatab_batch)

pcatab_batch<-pcatab_batch[,c(1:3)]

metadata$SampID

setdiff(metadata$SampID,row.names(pcatab_batch))

pcatab_batch<-merge(metadata,pcatab_batch,by.x="SampID",by.y="row.names")
pcatab_batch$var<-paste(pcatab_batch$Source,pcatab_batch$Pheno,sep="_")

library(ggnewscale)

pcatab_batch$GSEA<-as.factor(pcatab_batch$GSEA)

# Detecting outliers, estimate z-score
outbatchTLE<-data.frame(SampID=pcatab_batch$SampID,
                     zscore_PC1=abs(scale(pcatab_batch$PC1, center = TRUE, scale = TRUE)),
                     zscore_PC2=abs(scale(pcatab_batch$PC2, center = TRUE, scale = TRUE)))

pcatab_batch<-merge(pcatab_batch,outbatchTLE,by="SampID")


ggplot(pcatab_batch, aes(PC1, PC2,label=SampID)) +
  #geom_point(data=pcatab[pcatab$source=="BM",],size=10,color="black") +
  #geom_point(data=pcatab[pcatab$source=="Thymus",],size=10,color="grey") +
  #geom_point(data=pcatab[pcatab$source=="TALL",],size=10,color="red") +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=GSEA)) +
  #geom_point(data=pcatab[pcatab$GSEA=="GSE57982",],size=6,color="black") +
  scale_color_brewer(palette = "Paired")+
  new_scale_color()+
  stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
  scale_color_brewer(palette = "Paired")+
  #scale_linetype_manual(values=c("solid","dotted"))+
  #geom_label_repel(data=pcatab_batch[pcatab_batch$GSEA=="GSE69239",],aes(label=Source), size=4,fontface = "italic",alpha=0.7)+
  #geom_label_repel(data=pcatab_batch[pcatab_batch$Source=="Thymus" ,],aes(label=Source), size=4,fontface = "italic",alpha=0.7,color="grey")+
  geom_label_repel(data=pcatab_batch[pcatab_batch$zscore_PC1>=3 | pcatab_batch$zscore_PC2>=3 ,],aes(label=SampID), size=3,fontface = "italic",alpha=0.7,color="black")+
  ylab("PC2 (6.08% variance)")+
  xlab("PC1 (11.80% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()


## Adding cibersort info to PCA scatterpies
pcatab_batch<-merge(pcatab_batch,cibersort,by="SampID",all.x=TRUE)

ggplot(pcatab_batch, aes(PC1, PC2,label=SampID))+
  geom_scatterpie(aes(x=PC1,y=PC2,group=SampID,r=3.5),
                  data=pcatab_batch, cols=colnames(cibersort)[c(2:11)],alpha=0.8)+
  scale_fill_brewer(palette = "Paired")+
  new_scale_color()+
  geom_point(data=pcatab_batch[pcatab_batch$GSEA=="GSE69239",],color="black",size=4,shape=17) +
  geom_point(data=pcatab_batch[pcatab_batch$GSEA=="GSE69239",],aes(color=Source),size=3,shape=17) +
  scale_color_brewer(palette = "Paired")+
  geom_label_repel(data=pcatab_batch[pcatab_batch$GSEA=="GSE69239",],aes(label=Source), size=4,fontface = "italic",alpha=0.7)+
  #geom_point(data=pcatab[pcatab$GSEA=="GSE57982",],size=6,color="black") +
#  new_scale_color()+
#  stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
#  scale_color_brewer(palette = "Paired")+
  #scale_linetype_manual(values=c("solid","dotted"))+
#  geom_label_repel(data=pcatab_batch[pcatab_batch$Source=="Thymus" ,],aes(label=Source), size=4,fontface = "italic",alpha=0.7,color="grey")+
  #scale_color_manual(values=c("purple","red"))+
  #ylim(-50,50)+
  ylab("PC2 (6.08% variance)")+
  xlab("PC1 (11.80% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()
#####

## Representation of genes and events
#####

# corrected batch
melt_vastout<-as.data.frame(ty2)

melt_vastout$sample<-row.names(melt_vastout)
melt_vastout<-melt(melt_vastout,id.vars="sample")

table(melt_vastout$sample)

melt_vastout<-merge(melt_vastout,metadata,by.x="sample",by.y="SampID")

table(melt_vastout$Pheno)

melt_vastout$Pheno<-as.factor(melt_vastout$Pheno)

melt_vastout$Pheno<-factor(melt_vastout$Pheno,levels = c("BM","Thymus","pre-T-ALL","IMM","c-T-ALL","Cortical","c-T-LBL","Pre-cortical","Post-cortical",
                                                         "mature","NA","control","shbcat"))

table(melt_vastout$Lesion)
melt_vastout$Lesion<-factor(melt_vastout$Lesion,levels = c("BM","Thymus","ETV6-NCOA2","HOXA","LMO1/2","LMO2_LYL1","NKX2_1","TAL","TAL1",
                                                         "TAL2","TLX","TLX1","TLX3","Unknown","NA","control","shbcat"))


table(melt_vastout$Source)
melt_vastout$Source<-factor(melt_vastout$Source,levels = c("HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8","Thymus","TALL","RPMI","Jurkat"))

melt_vastout<-merge(melt_vastout,geneid,by.x="variable",by.y="Gene")

# without line
#melt_vastout_noline<-melt_vastout[melt_vastout$GSEA!="line",]


# add variable type, leukemic or not
melt_vastout$state<-ifelse(melt_vastout$Type=="TALL" | melt_vastout$Type=="TLBL", 
                          c("leukemia"), c("normal")) 

# add variable for cell line

melt_vastout$state<-ifelse(melt_vastout$Type=="TALL" | melt_vastout$Type=="TLBL", 
                           c("leukemia"), c("normal")) 


# for batch corrected, free_x. If not corrected, free.
c1<-ggplot(data=melt_vastout[melt_vastout$NAME=="NOTCH1",],aes(x=state,y=value))+
  geom_point(color="black",size=6)+
  geom_point(data=melt_vastout[melt_vastout$GSEA!="line" & melt_vastout$NAME=="NOTCH1",],aes(color=Type),size=5)+
  geom_boxplot(data=melt_vastout[melt_vastout$NAME=="NOTCH1",],alpha=0)+
  #geom_point(data=melt_vastout[melt_vastout$GSEA=="GSE69239" & melt_vastout$NAME=="RBM39",],size=5)+
  geom_point(data=melt_vastout[melt_vastout$GSEA=="line" & melt_vastout$NAME=="NOTCH1",],aes(color=Pheno),size=5)+
  #scale_shape_manual(values=c(0,8,2,5,6,7,10,11,12,9))+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(~GSEA,scales="free_x",nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        strip.text = element_text(size=9),
        legend.position = "bottom")

c1
rm(c1)

## Non corrected

#melt non corrected batch
#melt_vastout_noncor<-vastout_t

#melt_vastout_noncor$sample<-row.names(melt_vastout_noncor)
#melt_vastout_noncor<-melt(melt_vastout_noncor,id.vars="sample")

#melt_vastout_noncor$value<-as.numeric(melt_vastout_noncor$value)

#table(melt_vastout_noncor$sample)

#melt_vastout_noncor<-merge(melt_vastout_noncor,metadata,by.x="sample",by.y="SampID")

#table(melt_vastout_noncor$Pheno)

#melt_vastout_noncor$Pheno<-as.factor(melt_vastout_noncor$Pheno)

#melt_vastout_noncor$Pheno<-factor(melt_vastout_noncor$Pheno,levels = c("BM","Thymus","pre-T-ALL","IMM","c-T-ALL","Cortical","c-T-LBL","Pre-cortical","Post-cortical",
#                                                         "mature","NA","control","shbcat"))

#table(melt_vastout_noncor$Source)
#melt_vastout_noncor$Source<-factor(melt_vastout_noncor$Source,levels = c("HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8","Thymus","TALL","RPMI","Jurkat"))


#melt_vastout_noncor$Type<-factor(melt_vastout_noncor$Type,levels = c("BM","Thymus","TALL","TLBL"))

#melt_vastout_noncor<-merge(melt_vastout_noncor,geneid,by.x="variable",by.y="Gene")

# without line
#melt_vastout_noline_noncor<-melt_vastout_noncor[melt_vastout_noncor$GSEA!="line",]


# add variable type, leukemic or not
#melt_vastout_noncor$state<-ifelse(melt_vastout_noncor$Type=="TALL" | melt_vastout_noncor$Type=="TLBL", 
#                           c("leukemia"), c("normal")) 


# If not corrected, free.
#nc<-ggplot(data=melt_vastout_noncor[melt_vastout_noncor$NAME=="NOTCH1",],aes(x=state,y=value))+
#  geom_point(color="black",size=6)+
#  geom_point(data=melt_vastout_noncor[melt_vastout_noncor$GSEA!="line" & melt_vastout_noncor$NAME=="NOTCH1",],aes(color=Type),size=5)+
#  geom_boxplot(data=melt_vastout_noncor[melt_vastout_noncor$NAME=="NOTCH1",],alpha=0)+
#  #geom_point(data=melt_vastout[melt_vastout$GSEA=="GSE69239" & melt_vastout$NAME=="RBM39",],size=5)+
#  geom_point(data=melt_vastout_noncor[melt_vastout_noncor$GSEA=="line" & melt_vastout_noncor$NAME=="NOTCH1",],aes(color=Pheno),size=5)+
#  #scale_shape_manual(values=c(0,8,2,5,6,7,10,11,12,9))+
#  scale_color_brewer(palette="Spectral")+
#  facet_wrap(~GSEA,scales="free",nrow=1)+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle=45, hjust = 1),
#        strip.text = element_text(size=9),
#        legend.position = "bottom")
#nc




#grid.arrange(c1,nc,ncol=2)

### Batch effect corrected, grouping T-ALL

melt_vastout$Type<-factor(melt_vastout$Type,levels = c("BM","Thymus","TALL","TLBL"))

ggplot(data=melt_vastout[melt_vastout$NAME=="NOTCH1",],aes(x=Type,y=value))+
  geom_point(color="black",size=6)+
  geom_point(data=melt_vastout[melt_vastout$GSEA!="line" & melt_vastout$NAME=="NOTCH1",],aes(color=Type),size=5)+
  geom_boxplot(data=melt_vastout[melt_vastout$NAME=="NOTCH1",],alpha=0)+
  #geom_point(data=melt_vastout[melt_vastout$GSEA=="GSE69239" & melt_vastout$NAME=="RBM39",],size=5)+
  geom_point(data=melt_vastout[melt_vastout$GSEA=="line" & melt_vastout$NAME=="NOTCH1",],aes(color=Pheno),size=5)+
  #scale_shape_manual(values=c(0,8,2,5,6,7,10,11,12,9))+
  scale_color_brewer(palette="Spectral")+
  #facet_wrap(~GSEA,scales="free",nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1,size=15),
        legend.position = "bottom")

#####


## Wnt signaling signature
#####

## Check expression Wnt signaling genes (https://web.stanford.edu/group/nusselab/cgi-bin/wnt/main)

wntsig<-read.delim("/Volumes/grcmc/YGUILLEN/beta_catenin_project/gene_lists/Wnt_targets.txt",header = FALSE)
#wntsig<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/gene_lists/Wnt_targets.txt",header = FALSE)
colnames(wntsig)<-"NAME"
dim(wntsig)

wntsig<-merge(wntsig,geneid,by="NAME")
dim(wntsig)


# NO TLE Samples
wntexp<-vastout_t_TLE[,colnames(vastout_t_TLE) %in% wntsig$Gene,]
dim(wntexp)

# select matrix for cell lines
linescell<-c("control_Jurkat","control_RPMI","shbcat_Jurkat","shbcat_RPMI")
#wntexp<-wntexp[row.names(wntexp) %in% linescell,]

# select matrix for GSE57982
GSE57982<-row.names(metadata[metadata$GSEA=="GSE57982",])

#wntexp<-wntexp[row.names(wntexp) %in% GSE57982,]

# select matrix for cell lines and cohort GSE57982
#wntexp<-wntexp[row.names(wntexp) %in% c(linescell,GSE57982),]


library(heatmaply)

#https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
# https://rdrr.io/cran/heatmaply/man/heatmaply.html

# markers matrix with no NAs
complete.cases(wntexp)
is.na(wntexp)
wntexp<-wntexp[complete.cases(wntexp),]
dim(wntexp)

#wntexp<-wntexp[, -which(numcolwise(sum)(wntexp) ==0)]
length(colnames(wntexp))
length((geneid[geneid$Gene %in% colnames(wntexp),])$NAME)

colnames(wntexp)<-(geneid[geneid$Gene %in% colnames(wntexp),])$NAME


# Wnt targets with no expression
noexp<-colnames(wntexp[, which(numcolwise(sum)(wntexp) ==0)])
noexp

# Wnt targets with expressed,Remove genes with sum 0 in all samples
exp<-colnames(wntexp[, which(numcolwise(sum)(wntexp) !=0)])
exp

wntexp<-wntexp[,colnames(wntexp) %in% exp]


row.names(metadata)<-metadata$SampID
wntexp<-merge(subset(metadata,select=c(Source,Pheno,Type,GSEA)),wntexp,by="row.names",sort=FALSE)
dim(metadata)
dim(wntexp)

row.names(wntexp)<-wntexp$Row.names
wntexp$Row.names<-NULL

# RColorBrewer::brewer.pal(7, "Spectral")
# "#D53E4F" "#FC8D59" "#FEE08B" "#FFFFBF" "#E6F598" "#99D594" "#3288BD"

#normalize and scale are different. normalize substracts the minimum, and divide by
# the maximum of all observations. This preserves the shape of each variable's distribution
# while making them easily comparable on the same "scale".


heatmaply(scale(wntexp[,5:ncol(wntexp)]),
          plot_method = "plotly",
          #row_side_palette = colclust,
          RowSideColors = wntexp[,c(2:4)],
         scale="row",
          k_row=2,
          k_col=2)

library(RColorBrewer)
cols <- colorRampPalette(rev(brewer.pal(10,name="RdBu")))(50)
paletteLength<-length(cols)
#cols <- colorRampPalette(c("red", "green", "blue"))(50)
myBreaks <- c(seq(min(t(scale(wntexp[,5:ncol(wntexp)]))), 2.5, length.out=round(ceiling(paletteLength*0.90)) + 1), 
              seq(2.6, max(t(scale(wntexp[,5:ncol(wntexp)]))), length.out=round(floor(paletteLength*0.10)-1))) 

pheatmap(as.data.frame(t(scale(wntexp[,5:ncol(wntexp)]))),
         annotation_col = as.data.frame(wntexp[,c(3,4)]),
         show_colnames=T,
         fontsize = 6,
         cutree_cols = 4,
         cutree_rows = 4,
         color=cols, breaks=myBreaks)


# Empty space for subsequent data
rm(ty2)
rm(res.pca)
rm(res.pcaTLE)
rm(c1)

#####



#### PROFILING of gene expression to identify MOLECULAR SUBGROUPS 
#####

# merge both datasets and remove names columns
talline<-merge(linecount,tallcount,by="ID")
talline$NAME.y<-NULL
names(talline)[2]<-"NAME"

dim(talline)
colnames(talline)

## Remove TLE data
metadata[metadata$GSEA=="EGA_TLE",]$SampID
talline<-talline[,!colnames(talline) %in% metadata[metadata$GSEA=="EGA_TLE",]$SampID]

# set row.names and colnames
row.names(talline)<-talline$ID
talline$NAME<-NULL
talline$ID<-NULL
colnames(talline)

# To matrix
talline<-as.matrix(talline)
talline<-talline[complete.cases(talline),]

metadata$SampID
row.names(metadata)
colnames(metadata)

# Separate by cohort
table(metadata$GSEA)

## DESeq package

# Need to be in the same order. As in metadata the order is alphabetical, we reorder matrix counts alphabetically
colnames(talline)
row.names(metadata)

talline<-talline[ , order(colnames(talline))]


## All samples
metTCel<-metadata
TCelcounts<-talline


## All T-ALL and T-LBL samples not thymus, no TLE
metTCel<-metadata[((metadata$Source=="TALL" | metadata$Source=="TLBL") & metadata$SampID %in% colnames(talline)),]
TCelcounts<-talline[,colnames(talline) %in% metTCel$SampID]


# For estimating gene variance each cohort independently, exlude thymus!

# T-cell development
metTCel<-metadata[metadata$GSEA=="GSE69239",]
TCelcounts<-talline[,colnames(talline) %in% metTCel$SampID]

# T-ALL EGA_R
metTCel<-metadata[metadata$GSEA=="EGA_R",]
TCelcounts<-talline[,colnames(talline) %in% metTCel$SampID]
dim(TCelcounts)


# T-ALL GSE110633
metTCel<-metadata[metadata$GSEA=="GSE110633",]
TCelcounts<-talline[,colnames(talline) %in% metTCel$SampID]

# T-ALL GSE109231
metTCel<-metadata[metadata$GSEA=="GSE109231" & metadata$Source!="Thymus",]
TCelcounts<-talline[,colnames(talline) %in% metTCel$SampID]

# T-ALL GSE57982 (use ~1 in deseq object instead of Pheno)
metTCel<-metadata[metadata$GSEA=="GSE57982" & metadata$Source!="Thymus",]
TCelcounts<-talline[,colnames(talline) %in% metTCel$SampID]

# T-ALL TARGET
metTCel<-metadata[metadata$GSE=="TARGET",]
TCelcounts<-talline[,colnames(talline) %in% metTCel$SampID]


# OUTPUT TO CEMITOOL (R version not compatible!)
write.table(metTCel[,c(1,6)],"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TCelldev_cemitool_TCell.tsv",quote = FALSE,row.names=FALSE,sep="\t")
write.table(TCelcounts,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TCelldev_cemitool_TCell_exp.tsv",quote = FALSE,sep="\t")


#### Cluster genes with DESEq

# make deseq object
TALLCountTable <- DESeqDataSetFromMatrix(
  countData=TCelcounts,
  colData = metTCel,
  ~ Pheno)

colnames(TALLCountTable)

#how many genes we capture, counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(TALLCountTable)
class(GeneCounts)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)


### random sample from the count matrix
nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]


### NORMALIZATION
# Remove genes with low expression levels ()
# 15 counts in minimun three samples
keep <- rowSums(counts(TALLCountTable) >= 15) >= 3
TALLCountTable <- TALLCountTable[keep,]
row.names(TALLCountTable)

#### estimate size factors
TALLCountTable <- estimateSizeFactors(TALLCountTable)
sizeFactors(TALLCountTable)
dim(TALLCountTable)

#rld<-rlogTransformation(TALLCountTable,blind = FALSE)
rld<-vst(TALLCountTable,blind = FALSE)


#### THIS STEP ONLY TO ESTIMATE VARIANCE FOR EACH COHORT INDEPENDENTLY (no batch effect applied yet)
take<-rowVars( assay(rld)) >= 0.5
genevar<-row.names(assay(rld)[take,])
dim(assay(rld)[take,])

geneid[geneid$Gene %in% genevar,]

# CAREfUL!
genevar_GSE69239<-genevar
genevar_R<-genevar
genevar_GSE110633<-genevar
genevar_GSE109231<-genevar
genevar_GSE57982<-genevar
genevar_TARGET<-genevar

length(genevar_GSE69239)
length(genevar_R)
length(genevar_GSE110633)
length(genevar_GSE109231)
length(genevar_GSE57982)
length(genevar_TARGET)

# Sum of genes that have a variance higher than 0.5 in each cohort independently TO REST LATER (excluding TLE)
vartake<-unique(c(genevar_R,genevar_GSE110633,genevar_GSE109231,genevar_GSE57982,genevar_TARGET))
length(vartake)

## Select rld from all T-ALL and T-LBL datasets
hvg<-assay(rld)[rownames(assay(rld)) %in% vartake,]

#for TARGET ONLY
hvg<-assay(rld)[rownames(assay(rld)) %in% genevar_TARGET,]

dim(hvg)

# Two approaches to scale
#hvgplot  <- hvg - rowMeans(hvg)
hvgplot  <- t(apply(hvg, 1, scale))
colnames(hvgplot)<-colnames(hvg)

row.names(metTCel)<-metTCel$SampID

anno <- subset(metTCel,select=c(Pheno,Source,GSEA,First_event))


#random subset
nsamp<-sample(row.names(hvgplot), 2000)

pheatmap(hvgplot[row.names(hvgplot) %in% nsamp,],annotation_col = anno,fontsize = 8,
                    cutree_cols = 7,
                    cutree_rows = 6,
                    show_rownames = FALSE, show_colnames = FALSE)

geneid[geneid$Gene %in% nsamp,]
names(pheatplot)



### apply removebatch effect in all TALL and TLBL datasets (rld from all data, including only genes with var 0.5 unique in each cohort)
rld_batch<-removeBatchEffect(as.data.frame(assay(rld)),(colData(rld))$GSEA)

# For cemitool analysis later (include all genes, not only highly varirable, with batch effect removal)
write.table(rld_batch,'/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TALL_cemitool_clusters_uns_exp.tsv',quote = FALSE,sep="\t")
#write.table(rld_batch,'/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/CEMITools/TALL_cemitool_clusters_uns_exp.tsv',quote = FALSE,sep="\t")

# Select highly variable for further analysis
rld_batch<-rld_batch[rownames(assay(rld_batch)) %in% vartake,]
dim(rld_batch)

#for TARGET
rld_batch<-assay(rld)[rownames(assay(rld)) %in% genevar_TARGET,]
dim(rld_batch)


# Select highly variable genes from remove batch effect dataset
takebatch<-rowVars( rld_batch) >= 0.5
length(takebatch)

genevarbatch<-row.names(rld_batch[takebatch,])
length(genevarbatch)
dim(rld_batch[takebatch,])

#sort by variance
sortvar<-data.frame(variance=rowVars(rld_batch),Gene=row.names(rld_batch))
sortvar<-merge(sortvar,geneid,by="Gene")
sortvar<-sortvar[rev(order(sortvar$variance)),]

sortvar <- transform(sortvar, NAME=reorder(NAME, -variance) ) 

#Plot variance
ggplot(sortvar[1:50,],aes(x=NAME,y=variance))+
  geom_point(size=5)+
  theme_classic()+
  theme(axis.text.y=element_text(size=12,face="bold.italic"),
        axis.title=element_text(size=12,face="bold.italic"))+
  coord_flip()



# FOR VISUALIZING, CHECK NUMBER OF GENES IN TOPVARGENES
#heatmap
#heatmap.2( assay(rld)[take,], scale="row",
#           trace="none", dendrogram="column",
#           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


#pheatmap




# From highly variable genes, how are they co-expressed?
#hvg<-t(hvg)
#chvg<-rcorr(hvg)

#corrplot(chvg$r[1:200,1:200], order="hclust", 
#         p.mat = chvg$P[1:200,1:200], sig.level = 0.05, insig = "blank",addrect =6,tl.cex=0.3,
#         tl.col = "black")

#chvg$r

## Detecting modules of correlated genes from CEMITools output





#####



## SUPERVISED CLUSTERING
#####
#Using dataframe with scaled expression 
# supervised clustering, needed two conditions, not multiple labels
# all genes, not selecting only highly variable rowVar > 0.5. (line 1102)
assay(rld)
supwil<-t(assay(rld) - rowMeans(assay(rld)))
supwil[1:10,1:10]

row.names(supwil)
#Exclude samples that have no phenotype assigned ()
supwil<-supwil[row.names(supwil) %in% (metadata[metadata$Pheno!="NA",])$SampID,]

#vector with phenotype
vecmet<-metadata[metadata$SampID %in% row.names(supwil),]$Pheno
table(vecmet)

# only two possible categories
vecmet<-gsub("c-T-ALL","noIMM",vecmet)
vecmet<-gsub("c-T-LBL","noIMM",vecmet)
vecmet<-gsub("HOXA","noIMM",vecmet)
vecmet<-gsub("mature","noIMM",vecmet)
vecmet<-gsub("pre-T-ALL","IMM",vecmet)
vecmet<-gsub("TAL","noIMM",vecmet)
vecmet<-gsub("TLX","noIMM",vecmet)

table(vecmet)

vecmet<-gsub("IMM","0",vecmet)
vecmet<-gsub("no0","1",vecmet)


## wilma module
library(supclust)

set.seed(724)

## Fitting Wilma # CAREGULL WITH THE NUMBER OF GENES
fit  <- wilma(supwil, vecmet, noc = 3, trace = 1)

## Working with the output
fit
summary(fit)
plot(fit)
fitted(fit)

######

# T-CELL GENE SIGNATURE FROM CEMITOOLS
#####

tcell_cemi<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TCelldev_cemitool_results/module.tsv")
mean_tcell<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TCelldev_cemitool_results/summary_mean.tsv")

mean_tcell<-melt(mean_tcell)
mean_tcell$variable<-factor(mean_tcell$variable,levels = c("HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8"))
mean_tcell$group<-"group"

ggplot(mean_tcell,aes(x=variable,y=value))+
  geom_point()+
  geom_line(aes(group=group,color=modules),size=2)+
  scale_color_brewer(palette="Set2")+
  facet_wrap(~modules,scales="free_y",ncol=1)+
  theme_bw()

# MODULE 2 gene names
tcell_cemi[tcell_cemi$modules=="M2",]
sort(geneid[geneid$Gene %in% tcell_cemi[tcell_cemi$modules=="M2",]$genes,]$NAME)


# Select genes from pheatmap matrix to represent using all samples (go back to line 1033!!!!! and run until rld) or select rld_batch
#rld
#tmat<-assay(rld)[rownames(assay(rld)) %in% tcell_cemi$genes,]

#rld batch
tmat<-rld_batch[rownames(rld_batch) %in% tcell_cemi$genes,]


#tmatplot  <- tmat - rowMeans(tmat)
tmatplot  <- t(apply(tmat, 1, scale))
colnames(tmatplot)<-colnames(tmat)


annotcel <- subset(metadata,select=c(SampID,Type,Pheno))
row.names(annotcel)<-annotcel$SampID
annotcel$SampID<-NULL
annotcel <- annotcel[row.names(annotcel) %in% colnames(rld_batch),]
annotcel$Pheno <- gsub('immature','IMM',annotcel$Pheno)

row.names(tcell_cemi)<-tcell_cemi$genes
tcell_cemi_map<-tcell_cemi
tcell_cemi_map$genes<-NULL

# Create color scale
library(RColorBrewer)
s1<-length(levels(as.factor(annotcel$Type)))
s2<-length(unique(tcell_cemi$modules))
s3<-length(levels(as.factor(annotcel$Pheno)))

myColors = list(
  Type = c(TALL=brewer.pal(s1,"BrBG")[1],
           TLBL=brewer.pal(s1,"BrBG")[2]),
  modules = c(M1=brewer.pal(s2,"Set2")[1],
              M2=brewer.pal(s2,"Set2")[2],
              M3=brewer.pal(s2,"Set2")[3],
              M4=brewer.pal(s2,"Set2")[4],
              Not.Correlated=brewer.pal(s2,"Set2")[5]),
  Pheno = c("NA"=brewer.pal(s3,"Paired")[1],
            "c-T-ALL"=brewer.pal(s3,"Paired")[2],
            "c-T-LBL"=brewer.pal(s3,"Paired")[9],
            mature=brewer.pal(s3,"Paired")[3],
            IMM=brewer.pal(s3,"Paired")[4],
            TAL=brewer.pal(s3,"Paired")[5],
            TLX=brewer.pal(s3,"Paired")[6],
            HOXA=brewer.pal(s3,"Paired")[7],
            "pre-T-ALL"=brewer.pal(s3,"Paired")[8]))

pheatmap(tmatplot,
         annotation_col = annotcel,
         annotation_row = tcell_cemi_map,
         show_rownames = FALSE,
         fontsize = 14,
         cutree_cols = 9,
         cutree_rows = 7,
         fontsize_col = 5,
         annotation_colors = myColors)

#####


### CONSENSUS CLUSTER PLUS (UNSUPERVISED)
#####
#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(Biobase)

## Using filtered data from Deseq2, using normalized counts, no batch effect removed
# rld will be all data to analize T-ALL

#ccp<-assay(rld)[vartake,]
#ccp<-assay(rld)[rownames(assay(rld)) %in% vartake,]

#Using removebatch effect data
ccp<-rld_batch[takebatch,]

# only TARGET cohort
ccp<-rld_batch

dim(ccp)


# Remove control thymus and T-cell dev data for signatures
colnames(ccp)
length(rownames(ccp))

#median center genes
dc = sweep(ccp,1, apply(ccp,1,median))

# run consensus cluster, with standard options
rcc = ConsensusClusterPlus(dc,
                           maxK=4,
                           reps=100,
                           pItem=0.8,
                           pFeature=1,
                           title="con_clust_TALL",
                           distance="pearson",
                           clusterAlg="hc",plot='pdf')


rm(dc)

rcc[[4]]$consensusClass

#For .example, the top five rows and columns of results for k=2:
rcc[[4]][["consensusMatrix"]][1:10,1:10]

rcc[[4]][["consensusTree"]]

rcc[[4]][["consensusClass"]]

#consensus matrix results
rcc[[4]][["ml"]]


icl = calcICL(rcc,title="itcon_clust_TALL_noTLE",plot='pdf')
icl[["itemConsensus"]][1:4,]



## add to metadata consensus class

cuns<-data.frame(rcc[[4]][["consensusClass"]])
colnames(cuns)<-"uns_clust"
cuns$SampID<-row.names(cuns)

#write.table(cuns,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/CEMITools/TALL_cemitool_clusters_uns.tsv",row.names = FALSE,sep="\t",quote = FALSE)
write.table(cuns,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TALL_cemitool_clusters_uns.tsv",row.names = FALSE,sep="\t",quote = FALSE)

hvgplot  <- t(apply(rld_batch[takebatch,], 1, scale))
colnames(hvgplot)<-colnames(rld_batch[takebatch,])

#only TARGET cohort
hvgplot  <- t(apply(rld_batch, 1, scale))
colnames(hvgplot)<-colnames(rld_batch)
dim(hvgplot)

announ <- subset(metadata[metadata$SampID %in% colnames(hvgplot),],select=c(SampID,Pheno,Source,GSEA,First_event))
row.names(announ)<-announ$SampID
announ$SampID<-NULL
table(announ$First_event)
dim(announ)

ccp_met<-merge(cuns,announ,by="row.names")
row.names(ccp_met)<-row.names(announ)
ccp_met<-subset(ccp_met,select=c("Source","uns_clust","Pheno","GSEA","First_event"))
ccp_met$uns_clust<-gsub("^","C",ccp_met$uns_clust)


row.names(ccp_met)

library(RColorBrewer)

s1<-length(unique(ccp_met$Source))
s2<-length(unique(ccp_met$uns_clust))
s3<-length(unique(ccp_met$Pheno))
s4<-length(unique(ccp_met$First_event))

myColors2 = list(
  Source = c(TALL=brewer.pal(s1,"Set2")[1],
             TLBL=brewer.pal(s1,"Set2")[2]),
  uns_clust = c(C1=brewer.pal(s2,"BrBG")[1],
                C2=brewer.pal(s2,"BrBG")[2],
                C3=brewer.pal(s2,"BrBG")[3],
                C4=brewer.pal(s2,"BrBG")[4]),
#                C5=brewer.pal(s2,"BrBG")[5]),
  Pheno = c("NA"="grey",
            "c-T-ALL"="chocolate",
            "c-T-LBL"="burlywood2",
            mature="dodgerblue",
            IMM="darkseagreen1",
            Cortical="chocolate",
            "Post-cortical"="dodgerblue4",
            "Pre-cortical"="darkseagreen3",
            "pre-T-ALL"="darkseagreen2"),
  First_event = c("NA"="white",
           "Censored"=brewer.pal(s4,"Dark2")[1],
           "Death"=brewer.pal(s4,"Dark2")[2],
           "Disease"=brewer.pal(s4,"Dark2")[3],
           "Progression"=brewer.pal(s4,"Dark2")[4],
           "Refractory"=brewer.pal(s4,"Dark2")[5],
           "Relapse"=brewer.pal(s4,"Dark2")[6],
           "Remission"=brewer.pal(s4,"Dark2")[7],
           "Second Malignant Neoplasm"=brewer.pal(s4,"Dark2")[8],
           "Secondary leukemia"=brewer.pal(s4,"Dark2")[8]))

myColors2
#order samples by clust category

hvgplot[, order(ccp_met$uns_clust)]
hvgplot<-hvgplot[,row.names(ccp_met[order(ccp_met$uns_clust),])]

#hvgplot[, order(ccp_met$First_event)]
#hvgplot<-hvgplot[,row.names(ccp_met[order(ccp_met$First_event),])]

pheatmap(hvgplot,
         annotation_col = ccp_met,
                    fontsize = 8,
                    cutree_rows = 6,
                    show_rownames = FALSE, 
                    cluster_cols=FALSE,
                    show_colnames = FALSE,
                    annotation_colors = myColors2)




table(ccp_met$uns_clust,ccp_met$First_event)

# Run CEMITOOls using uns cluster metadata and all expression matrix (rld_batch)
uns_cemi<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TALL_allexp_unclust_cemitools_results_all/module.tsv")
mean_uns_cemi<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TALL_allexp_unclust_cemitools_results_all/summary_mean.tsv")

#uns_cemi<-read.delim("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/CEMITools/TALL_allexp_unclust_cemitool_results/module.tsv")
#mean_uns_cemi<-read.delim("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/CEMITools/TALL_allexp_unclust_cemitool_results/summary_mean.tsv")

mean_uns_TALL<-melt(mean_uns_cemi)
mean_uns_TALL<-merge(mean_uns_TALL,ccp_met,by.x="variable",by.y="row.names")
mean_uns_TALL$uns_clust<-as.factor(mean_uns_TALL$uns_clust)


mean_uns_TALL$variable <- factor(mean_uns_TALL$variable, levels=unique(mean_uns_TALL[order(mean_uns_TALL$uns_clust), ]$variable))

mean_uns_TALL$group="group"

## intercepts per group
C1_num<-((length(mean_uns_TALL[mean_uns_TALL$uns_clust=="C1",]$uns_clust))/6)+0.5
C2_num<-(length(mean_uns_TALL[mean_uns_TALL$uns_clust=="C2",]$uns_clust)/6)+C1_num
C3_num<-(length(mean_uns_TALL[mean_uns_TALL$uns_clust=="C3",]$uns_clust)/6)+C2_num
C4_num<-(length(mean_uns_TALL[mean_uns_TALL$uns_clust=="C4",]$uns_clust)/6)+C3_num
#C5_num<-(length(mean_uns_TALL[mean_uns_TALL$uns_clust=="C5",]$uns_clust)/6)+C4_num

ggplot(mean_uns_TALL,aes(x=variable,y=value))+
  geom_point()+
  geom_line(aes(group=group,color=modules),size=2)+
  geom_vline(xintercept = c(C1_num,C2_num,C3_num,C4_num))+
  scale_color_brewer(palette="Set2")+
  facet_wrap(~modules,scales="free_y",ncol=1)+
  theme_bw()+
  theme(axis.text.x = element_blank(),legend.position = "bottom")


# MODULE 2 gene names
sort(geneid[geneid$Gene %in% uns_cemi[uns_cemi$modules=="M1",]$genes,]$NAME)

#####


## Functional enrichment for different lists, including b-cat data and unsupervised, cemitools modules
#####

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "InterPro_Domains_2019",
         "GO_Cellular_Component_2017",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2016",
         "DisGeNET")

enrichedM1 <- enrichr(as.character(geneid[geneid$Gene %in% uns_cemi[uns_cemi$modules=="M1",]$genes,]$NAME), dbs)
enrichedM2 <- enrichr(as.character(geneid[geneid$Gene %in% uns_cemi[uns_cemi$modules=="M2",]$genes,]$NAME), dbs)
enrichedM3 <- enrichr(as.character(geneid[geneid$Gene %in% uns_cemi[uns_cemi$modules=="M3",]$genes,]$NAME), dbs)
enrichedM4 <- enrichr(as.character(geneid[geneid$Gene %in% uns_cemi[uns_cemi$modules=="M4",]$genes,]$NAME), dbs)
enrichedM5 <- enrichr(as.character(geneid[geneid$Gene %in% uns_cemi[uns_cemi$modules=="M5",]$genes,]$NAME), dbs)
enrichedM6 <- enrichr(as.character(geneid[geneid$Gene %in% uns_cemi[uns_cemi$modules=="Not.Correlated",]$genes,]$NAME), dbs)

enrichedpos <- enrichr(as.character((corout_sig[corout_sig$Rho>0.4,])$NAME), dbs)
enrichedneg <- enrichr(as.character((corout_sig[corout_sig$Rho<(-0.4),])$NAME), dbs)

enrichcomb <- enrichr(as.character(rpmi_al[rpmi_al$peak=="peak" & rpmi_al$comb!="peak_noDEG_noncor_nofit" & !is.na(rpmi_al$hgnc_symbol),]$hgnc_symbol),dbs)
enrichonlypeak <- enrichr(as.character(rpmi_al[rpmi_al$peak=="peak" & rpmi_al$comb=="peak_noDEG_noncor_nofit" & !is.na(rpmi_al$external_gene_name.x),]$external_gene_name.x),dbs)
enrichallpeak <- enrichr(as.character(rpmi_al[rpmi_al$peak=="peak" & !is.na(rpmi_al$hgnc_symbol),]$hgnc_symbol),dbs)
enrichnopeak <- enrichr(as.character(rpmi_al[rpmi_al$peak!="peak" & !is.na(rpmi_al$hgnc_symbol) & (rpmi_al$comb=="nopeak_DEG_cor_fit" | rpmi_al$comb=="nopeak_DEG_cor_nofit" | rpmi_al$comb=="nopeak_DEG_nocor_nofit"),]$hgnc_symbol),dbs)

heatbcat<-enrichr(as.character(geneid[geneid$Gene %in% rownames(as.data.frame(t(scale(bcatexp[bcatexp$GSEA=="TARGET",][,6:ncol(bcatexp[bcatexp$GSEA=="TARGET",])])))[phet$tree_row[["order"]],])[1:200],]$NAME),dbs)

clustup<-enrichr(as.character(genesup$NAME),dbs)
clustdown<-enrichr(as.character(genesdown$NAME),dbs)

clust1bcat<-enrichr(as.character(clusters[clusters$cluster==1,]$NAME),dbs)
clust2bcat<-enrichr(as.character(clusters[clusters$cluster==2,]$NAME),dbs)
clust3bcat<-enrichr(as.character(clusters[clusters$cluster==3,]$NAME),dbs)
clust4bcat<-enrichr(as.character(clusters[clusters$cluster==4,]$NAME),dbs)
clust5bcat<-enrichr(as.character(clusters[clusters$cluster==5,]$NAME),dbs)
clust6bcat<-enrichr(as.character(clusters[clusters$cluster==6,]$NAME),dbs)
clust7bcat<-enrichr(as.character(clusters[clusters$cluster==7,]$NAME),dbs)
clust8bcat<-enrichr(as.character(clusters[clusters$cluster==8,]$NAME),dbs)
clust9bcat<-enrichr(as.character(clusters[clusters$cluster==9,]$NAME),dbs)

clust1bcat<-enrichr(as.character(clusters[clusters$cluster==5 | clusters$cluster==2,]$NAME),dbs)
clust2bcat<-enrichr(as.character(clusters[clusters$cluster==1 | clusters$cluster==3 | clusters$cluster==4,]$NAME),dbs)

enrichalltarget<-enrichr(as.character(unique(rpmi_al[rpmi_al$peak=="peak",]$external_gene_name.y)),dbs)


posnrich <- enrichedpos[["InterPro_Domains_2019"]]
negrich <- enrichedneg[["InterPro_Domains_2019"]]

heatbcatrich<-heatbcat[["GO_Biological_Process_2018"]]

combrich <- enrichcomb[["GO_Biological_Process_2018"]]
onpeakrich <- enrichonlypeak[["GO_Biological_Process_2018"]]
allpeakrich <- enrichallpeak[["GO_Biological_Process_2018"]]
nopeakrich <- enrichnopeak[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]

clustuprich <- clustup[["GO_Biological_Process_2018"]]
clustdownrich <- clustdown[["GO_Biological_Process_2018"]]

M1rich <- enrichedM1[["GO_Biological_Process_2018"]]
M2rich <- enrichedM2[["GO_Biological_Process_2018"]]
M3rich <- enrichedM3[["GO_Biological_Process_2018"]]
M4rich <- enrichedM4[["GO_Biological_Process_2018"]]
M5rich <- enrichedM5[["GO_Biological_Process_2018"]]
M6rich <- enrichedM6[["GO_Biological_Process_2018"]]

enrichalltarget_g<-enrichalltarget[["GO_Biological_Process_2018"]]
enrichalltarget_g$group<-c("all_bcat_targets")

clustuprich$group<-c("UP")
clustdownrich$group<-c("DOWN")

clust1bcatrich <-clust1bcat[["GO_Biological_Process_2018"]]
clust2bcatrich <-clust2bcat[["GO_Biological_Process_2018"]]
clust3bcatrich <-clust3bcat[["GO_Biological_Process_2018"]]
clust4bcatrich <-clust4bcat[["GO_Biological_Process_2018"]]
clust5bcatrich <-clust5bcat[["GO_Biological_Process_2018"]]
clust6bcatrich <-clust6bcat[["GO_Biological_Process_2018"]]
clust7bcatrich <-clust7bcat[["GO_Biological_Process_2018"]]
clust8bcatrich <-clust8bcat[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
clust9bcatrich <-clust9bcat[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]



clust1bcatrich$group<-c("cluster1")
clust2bcatrich$group<-c("cluster2")
clust3bcatrich$group<-c("cluster3")
clust4bcatrich$group<-c("cluster4")
clust5bcatrich$group<-c("cluster5")
clust6bcatrich$group<-c("cluster6")
clust7bcatrich$group<-c("cluster7")
clust8bcatrich$group<-c("cluster8")
clust9bcatrich$group<-c("cluster9")

posnrich$group<-c("positive")
negrich$group<-c("negative")

combrich$group<-c("peak_combination")
onpeakrich$group<-c("only_peak")
allpeakrich$group<-c("all_peak")
nopeakrich$group<-c("no_peak")

heatbcatrich$group<-c("peak")

M1rich$group="M1"
M2rich$group="M2"
M3rich$group="M3"
M4rich$group="M4"
M5rich$group="M5"
M6rich$group="Not_cor"


allGO<-rbind(M1rich,M2rich,M3rich,M4rich,M5rich,M6rich)

allGO<-rbind(combrich,onpeakrich,allpeakrich,nopeakrich)

allGO<-rbind(clust1bcatrich,clust2bcatrich,clust3bcatrich,clust4bcatrich,clust5bcatrich,clust6bcatrich,clust7bcatrich)
#,clust8bcatrich,clust9bcatrich)

allGO<-enrichalltarget_g

allGO<-rbind(clustuprich,clustdownrich)

bpsub<-subset(allGO,allGO$Adjusted.P.value<0.05)
bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=unique((bpsub$Term))[rev(order(bpsub$P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))



#bpsub$Score<-ifelse(bpsub$group=="negative",bpsub$Combined.Score*(-1),bpsub$Combined.Score)
bpsub$Score<-ifelse(bpsub$group=="all_peak",bpsub$Combined.Score*(-1),bpsub$Combined.Score)

# for multiple groups
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$group, rev(abs(bpsub$P.value)))), ]$Term))

# for one group
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$Combined.Score)), ]$Term))


ggplot(bpsub,aes(y=Combined.Score,x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group))+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~group,scales="free",ncol=4)+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=8, hjust = 1),
        axis.text.y = element_text(size=10, hjust = 1))



bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$group, bpsub$Adjusted.P.value)), ]$Term))

### for facets multiple groups
ggplot(bpsub,aes(y=Combined.Score,x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group))+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~group,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.position = "none",
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=8, hjust = -1),
        axis.text.y = element_text(size=14, hjust = 1))

# heatmap
ggplot(bpsub, aes(x=group, y=Term)) + 
  geom_tile(aes(fill = -(log(P.value))),colour = "black")+
  scale_fill_gradient2(low = "white",
                       high = "red") +
  #scale_fill_gradient(low = "blue", high = "red")+
  theme_dark()+
  theme(axis.text.x = element_text(size=10, hjust = 1,angle=45),
        axis.text.y = element_text(size=9, hjust = 1,face="italic"))+
  coord_flip()


#  ggtitle("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")


unlist(strsplit(bpsub[bpsub$Term=="BRCA1 ENCODE",]$Genes,split=";"))
unlist(strsplit(bpsub[bpsub$Term=="ZBTB33 ENCODE",]$Genes,split=";"))

dt_bp<-data.frame(table(bpsub$Term))
dt_bp<-dt_bp[dt_bp$Freq>=4,]

# heatmap
ggplot(bpsub[bpsub$Term %in% unique(dt_bp$Var1),], aes(x=group, y=Term)) + 
  geom_tile(aes(fill = -(log(P.value))),colour = "black")+
  scale_fill_gradient2(low = "white",
                       high = "red") +
  #scale_fill_gradient(low = "blue", high = "red")+
  theme_minimal()+
  theme(axis.text.x = element_text(size=10, hjust = 1,angle=45),
        axis.text.y = element_text(size=9, hjust = 1,face="italic"),
        panel.grid.major = element_blank())+
  ggtitle("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")


enrich_NES<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CEMITools/TALL_allexp_unclust_cemitools_results_all/enrichment_nes.tsv")
row.names(enrich_NES)<-enrich_NES$pathway
enrich_NES$pathway<-NULL

# Heatmap with combined score of each bp

M1heat<-M1rich[M1rich$Adjusted.P.value<=0.01,c(1,8)]
colnames(M1heat)[2]<-"M1"
M2heat<-M2rich[M1rich$Adjusted.P.value<=0.01,c(1,8)]
colnames(M2heat)[2]<-"M2"
M3heat<-M3rich[M3rich$Adjusted.P.value<=0.01,c(1,8)]
colnames(M3heat)[2]<-"M3"
M4heat<-M4rich[M4rich$Adjusted.P.value<=0.01,c(1,8)]
colnames(M4heat)[2]<-"M4"
M5heat<-M5rich[M5rich$Adjusted.P.value<=0.01,c(1,8)]
colnames(M5heat)[2]<-"M5"

Moduleheat<-merge(M1heat,M2heat,by="Term",all=TRUE)
Moduleheat<-merge(Moduleheat,M3heat,by="Term",all=TRUE)
Moduleheat<-merge(Moduleheat,M4heat,by="Term",all=TRUE)
Moduleheat<-merge(Moduleheat,M5heat,by="Term",all=TRUE)
Moduleheat$Term<-gsub(' \\(GO.*','',Moduleheat$Term)
Moduleheat[is.na(Moduleheat)]<-'0'
row.names(Moduleheat)<-Moduleheat$Term
savename<-row.names(Moduleheat)
Moduleheat$Term<-NULL


Moduleheat<-apply(Moduleheat, 2 , function(x) as.numeric(x))
row.names(Moduleheat)<-savename

library(RColorBrewer)
cols <- colorRampPalette(c("white", "paleturquoise1", "paleturquoise2","paleturquoise3","paleturquoise4"))(50)
paletteLength<-length(cols)
#cols <- colorRampPalette(c("red", "green", "blue"))(50)

myBreaks <- c(seq(min(Moduleheat),20, length.out=round(ceiling(paletteLength*0.05))), 
              seq(20.1, 500, length.out=round(floor(paletteLength*0.90))),
              seq(500.1, max(Moduleheat), length.out = round(floor(paletteLength*0.05)))) 

pheatmap(Moduleheat,
         show_colnames=T,
         fontsize = 6,
        # cutree_cols = 4,
         cutree_rows = 5,
        treeheight_col = 0,
         color=cols,
        breaks=myBreaks)


cols <- colorRampPalette(c("blue", "white", "red"))(10)
paletteLength<-length(cols)

# order
enrich_NES<-enrich_NES[c("M5","M3","M4","M1","M2"),]

pheatmap(t(enrich_NES),
         show_colnames=T,
         fontsize = 6,
         cluster_cols=FALSE,
         # cutree_cols = 4,
         treeheight_col = 0,
         color = cols)
#####

## DEGs REMISSION - RELAPSE in TARGET cohort using DESeq
#####

# T-ALL TARGET (preliminar)
pretarg<-exp_genes[exp_genes$GSEA=="TARGET",]$SampID
targmet<-metadata[metadata$SampID %in% pretarg,]
# Remove censored and NAs
targmet<-targmet[targmet$First_event!="NA" & targmet$First_event!="Censored",]

TCeltarg<-talline[,colnames(talline) %in% targmet$SampID]
dim(TCeltarg)
dim(targmet)

# Create variable remission vs. non-remission
targmet$prognosis<-ifelse(targmet$First_event == "Remission", 
                               c("Remission"), c("non_Remission")) 
table(targmet$First_event,targmet$prognosis)

#### Cluster genes with DESEq

# make deseq object
TARGETCountTable <- DESeqDataSetFromMatrix(
  countData=TCeltarg,
  colData = targmet,
  ~ prognosis)

colnames(TARGETCountTable)

#how many genes we capture, counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(TARGETCountTable)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

### random sample from the count matrix
nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]


### NORMALIZATION
# Remove genes with low expression levels ()
keep <- rowSums(counts(TARGETCountTable)) >= 10
TARGETCountTable <- TARGETCountTable[keep,]

GeneCounts <- counts(TARGETCountTable)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

colData(TARGETCountTable)$prognosis <- factor(colData(TARGETCountTable)$prognosis, levels = c("Remission","non_Remission"))

#### estimate size factors
TARGETCountTable <- estimateSizeFactors(TARGETCountTable)
sizeFactors(TARGETCountTable)


### produce rlog-transformed data
rldTARGET <- vst(TARGETCountTable, blind=TRUE) ## create a distance matrix between the samples


###### Differential EXPRESSION ANALYSIS ###
TARGETCountTable <- estimateDispersions(TARGETCountTable)
plotDispEsts(TARGETCountTable)


# Statistical testing of DE genes
levels(TARGETCountTable$prognosis)
design(TARGETCountTable)

dds<-DESeq(TARGETCountTable)
resultsNames(dds)
rescondition<-results(dds)

# Example plot counts
plotCounts(dds, gene=which(rownames(rescondition)=="ENSG00000145592"), intgroup="prognosis")

head(rescondition)
summary(rescondition)
sum(rescondition$pvalue < 0.05,na.rm=TRUE)
sum(rescondition$padj < 0.1,na.rm=TRUE)
plotMA(rescondition)

# Selecting specific DEGs
#To dataframe (ckeck females, males or gender)
target_prog_degs<-as.data.frame(rescondition)
target_prog_degs<-merge(target_prog_degs,geneid,by.x=0,by.y="Gene",all.x=TRUE)

countab<-as.data.frame(counts(dds, normalized=TRUE))
countab<-t(countab)
countab<-melt(countab)
colnames(countab)<-c("SampID","Gene","value")
countab<-merge(countab,geneid,by="Gene")

countab<-merge(countab,metadata,by="SampID")
genesdown<-subset(target_prog_degs,target_prog_degs$padj<=0.1 & target_prog_degs$log2FoldChange<0)$Row.names
genesup<-subset(target_prog_degs,target_prog_degs$padj<=0.1 & target_prog_degs$log2FoldChange>0)$Row.names

countab$prognosis<-ifelse(countab$First_event == "Remission", 
                          c("Remission"), c("non_Remission")) 


countab$Age_Dx_days<-as.numeric(countab$Age_Dx_days)

ggplot(countab[countab$Gene %in% genesdown[1:10],],aes(x=prognosis,y=log(value)))+
  geom_point(color="black",size=5,aes(shape=Lesion),alpha=0.5)+
  geom_point(aes(color=First_event,shape=Lesion),size=4,alpha=0.5)+
  geom_boxplot(alpha=0.2)+
  scale_shape_manual(values=c(15,16,16,17,18,19,19,20,20,0))+
  scale_color_brewer(palette = "Paired")+
  ylab("Normalized counts")+
  xlab('')+
  facet_wrap(~NAME,scales = "free",ncol=5)+
  theme_classic()+
  theme(axis.text.x = element_text(size=8))

rm(TCeltarg)
rm(nz.counts)

## Heatmap for DEGs remission vs. nonremission
listgenes<-subset(target_prog_degs,target_prog_degs$pvalue<=0.01)$Row.names

# Using norm_counts TARGET
deg_rem<-counts(dds, normalized=TRUE)[row.names(counts(dds, normalized=TRUE)) %in% genesdown,]

row.names(ccp_met)
colnames(deg_rem)

myColors2 = list(
  Source = c(TALL=brewer.pal(s1,"Set2")[1],
             TLBL=brewer.pal(s1,"Set2")[2]),
  uns_clust = c(C1=brewer.pal(s2,"BrBG")[1],
                C2=brewer.pal(s2,"BrBG")[2],
                C3=brewer.pal(s2,"BrBG")[3],
                C4=brewer.pal(s2,"BrBG")[4]),
  #                C5=brewer.pal(s2,"BrBG")[5]),
  Pheno = c("NA"=brewer.pal(s3,"Paired")[1],
            "c-T-ALL"=brewer.pal(s3,"Paired")[2],
            "c-T-LBL"=brewer.pal(s3,"Paired")[9],
            mature=brewer.pal(s3,"Paired")[3],
            IMM=brewer.pal(s3,"Paired")[4],
            Cortical=brewer.pal(s3,"Paired")[5],
            "Post-cortical"=brewer.pal(s3,"Paired")[6],
            "Pre-cortical"=brewer.pal(s3,"Paired")[7],
            "pre-T-ALL"=brewer.pal(s3,"Paired")[8]),
  First_event = c("NA"="white",
                  "Censored"=brewer.pal(s4,"Dark2")[1],
                  "Death"=brewer.pal(s4,"Dark2")[2],
                  "Disease"=brewer.pal(s4,"Dark2")[3],
                  "Progression"=brewer.pal(s4,"Dark2")[4],
                  "Refractory"=brewer.pal(s4,"Dark2")[5],
                  "Relapse"=brewer.pal(s4,"Dark2")[6],
                  "Remission"=brewer.pal(s4,"Dark2")[7],
                  "Second Malignant Neoplasm"=brewer.pal(s4,"Dark2")[8],
                  "Secondary leukemia"=brewer.pal(s4,"Dark2")[8]))

#order samples by clust category

deg_rem[, order(ccp_met$uns_clust)]
deg_rem<-deg_rem[,row.names(ccp_met[order(ccp_met$uns_clust),])]

deg_rem[, order(ccp_met[ccp_met$GSEA=="TARGET",]$First_event)]
deg_rem<-deg_rem[,row.names(ccp_met[ccp_met$GSEA=="TARGET",][order(ccp_met[ccp_met$GSEA=="TARGET",]$First_event),])]

pheatmap(deg_rem,
         annotation_col = ccp_met[ccp_met$GSEA=="TARGET",],
         fontsize = 8,
         scale="column",
         cutree_rows = 6,
         show_rownames = FALSE,
         cluster_cols=FALSE,
         show_colnames = FALSE,
         annotation_colors = myColors2)
#####

## Survival curve with DEGs genes
#####

library(survival)
library(survminer)
library(dplyr)

sup_target<-countab[countab$NAME=="FAT3",]
sup_target<-subset(sup_target,select=c(First_event,Event_free_surv_days,Age_Dx_days,value))
sup_target$value<-log(sup_target$value)
sup_target<-sup_target[sup_target$First_event!="Censored" & sup_target$First_event!="NA",]


sup_target<-sup_target %>% mutate(FAT3group =  ifelse(value >=5, "FAT3_high", "FAT3_low"))
sup_target<-sup_target %>% mutate(fustat =  ifelse(First_event == "Remission", 1, 0))

sup_target$FAT3group<- factor(sup_target$FAT3group)

surv_object <- Surv(time = sup_target$Event_free_surv_days, event = sup_target$fustat)
fit1 <- survfit(surv_object ~ FAT3group, data = sup_target)
summary(fit1)
ggsurvplot(fit1, data = sup_target, pval = TRUE)


# DEGs from vastoutput

#####

### DEGs with vast-tools selecting RBPs
#####

DEGs_TLE<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/EGA_TLE_expr_out/gene_expression/DiffGE-Hsa124-fold2_EGAR00001193349-vs-EGAR00001193362-minReads50-minRPKM2-strict-range1.5.tab")
DEGs_TLE$cohort<-"EGA_TLE"
DEGs_GSE109231<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/GSE109231_expr_out/gene_expression/DiffGE-Hsa124-fold2_SRR6495835-vs-SRR6495838-minReads50-minRPKM2-strict-range1.5.tab")
DEGs_GSE109231$cohort<-"GSE109231"
DEGs_GSE57982<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/GSE57982_expr_out/gene_expression/DiffGE-Hsa124-fold2_GSM1399182-vs-GSM1399180-minReads50-minRPKM2-strict-range1.5.tab")
DEGs_GSE57982$cohort<-"GSE57982"
DEGs_TALL_DP<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/TALL_DP_expr_out/gene_expression/DiffGE-Hsa389-fold2_EGAR00001193349-vs-DP-minReads50-minRPKM2-strict-range1.5.tab")
DEGs_TALL_DP$cohort<-"TALL_DP"

DEGs_cohorts<-rbind(DEGs_TLE[,c(1,2,ncol(DEGs_TLE)-1,ncol(DEGs_TLE))],DEGs_GSE109231[,c(1,2,ncol(DEGs_GSE109231)-1,ncol(DEGs_GSE109231))])
DEGs_cohorts<-rbind(DEGs_cohorts,DEGs_GSE57982[,c(1,2,ncol(DEGs_GSE57982)-1,ncol(DEGs_GSE57982))])
DEGs_cohorts<-rbind(DEGs_cohorts,DEGs_TALL_DP[,c(1,2,ncol(DEGs_TALL_DP)-1,ncol(DEGs_TALL_DP))])

# cross with RBP data
Huang_RBP<-merge(geneid,Huang_RBP,by="NAME",all.y=TRUE)
Huang_RBP$ref<-"Huang"

Wang_RBP<-merge(geneid,Wang_RBP,by="NAME",all.y=TRUE)
Wang_RBP$ref<-"Wang"

Sebest_RBP<-merge(geneid,Sebest_RBP,by="NAME",all.y=TRUE)
Sebest_RBP$ref<-"Sebest"

rMATs_RBP<-rMATs_RBP[,c(1,2)]
rMATs_RBP<-rMATs_RBP[!duplicated(rMATs_RBP$Ensembl.Gene.ID), ]
rMATs_RBP<-merge(geneid,rMATs_RBP,by.x="Gene",by.y="Ensembl.Gene.ID",all.y=TRUE)
rMATs_RBP<-subset(rMATs_RBP,select=c(NAME,Gene))
rMATs_RBP$ref<-"rMATs"

RBP_refs<-merge(Huang_RBP,Wang_RBP,by="NAME",all=TRUE)
RBP_refs<-merge(RBP_refs,Sebest_RBP,by="NAME",all=TRUE)
RBP_refs<-merge(RBP_refs,rMATs_RBP[!is.na(rMATs_RBP$NAME),],by="NAME",all=TRUE)
colnames(RBP_refs)<-c("NAME","Gene.x","ref.x","Gene.y","ref.y","Gene.z","ref.z","Gene.p","ref.p")

RBP_refs$refs<-paste(RBP_refs$ref.x,RBP_refs$ref.y,RBP_refs$ref.z,RBP_refs$ref.p,sep="_")
RBP_refs<-RBP_refs[,c(1,2,10)]
RBP_refs$refs<-gsub('_NA','',RBP_refs$refs)
RBP_refs$refs<-gsub('NA_','',RBP_refs$refs)
table(RBP_refs$refs)
colnames(RBP_refs)[2]<-"GENE"

DEGs_cohorts<-merge(DEGs_cohorts,RBP_refs,by="NAME",all.x=TRUE)
DEGs_cohorts$state<-ifelse(DEGs_cohorts$Log2_Fold_Ch<0,c("UP_TALL"),c("DOWN_TALL"))
table(DEGs_cohorts$cohort,DEGs_cohorts$refs,DEGs_cohorts$state)

DEGs_cohorts_RBPs<-DEGs_cohorts[!is.na(DEGs_cohorts$refs),]
freqRBPs<-as.data.frame(table(DEGs_cohorts_RBPs$NAME))
freqRBPs<-freqRBPs[order(freqRBPs$Freq,decreasing = TRUE),]
names(freqRBPs)[1]<-"NAME"
DEGs_cohorts_RBPs<-merge(DEGs_cohorts_RBPs,freqRBPs,by="NAME")

# heatmap

DEGs_cohorts_RBPs_unique<-DEGs_cohorts_RBPs[!duplicated(DEGs_cohorts_RBPs$GENE.x),]
row.names(DEGs_cohorts_RBPs_unique)<-DEGs_cohorts_RBPs_unique$GENE.x

#merge with bcat info
DEGs_cohorts_RBPs_unique<-merge(DEGs_cohorts_RBPs_unique,colcol,by="row.names",all.x=TRUE)
row.names(DEGs_cohorts_RBPs_unique)<-DEGs_cohorts_RBPs_unique$Row.names
DEGs_cohorts_RBPs_unique$Row.names<-NULL

myBreaks <- c(seq(min(t(scale(vastout_t[,colnames(vastout_t_TLE) %in% c(DEGs_cohorts_RBPs_unique$GENE.x,"ENSG00000168036")]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(vastout_t[,colnames(vastout_t_TLE) %in% c(DEGs_cohorts_RBPs_unique$GENE.x,"ENSG00000168036")]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(vastout_t_TLE[,colnames(vastout_t_TLE) %in% c(DEGs_cohorts_RBPs_unique$GENE.x,"ENSG00000168036")]))),
               annotation_col = as.data.frame(metadata[metadata$SampID %in% row.names(vastout_t_TLE),c(2,6,13)]),
                              annotation_row = DEGs_cohorts_RBPs_unique[,c(4,7,8,58,78:81)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
#               cutree_cols = 15,
#               cutree_rows = 3,
               color=cols, breaks=myBreaks)


clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=4)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=4)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==2,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% deg_bcattarg])))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% row.names(target_prog_degs_sig)])))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% deg_bcattarg])))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% row.names(target_prog_degs_sig)])))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)



### MUTATIONS IN RBPs in TARGET cohort
somatic_target_RBP<-merge(somatic_target,RBP_refs,by.x="Hugo_Symbol",by.y="NAME",all.x=TRUE)
length(unique(somatic_target_RBP[!is.na(somatic_target_RBP$refs),]$Hugo_Symbol))
somatic_target_RBP$Tumor_Sample_Barcode<-gsub('-0.*','',somatic_target_RBP$Tumor_Sample_Barcode)
somatic_target_RBP<-merge(somatic_target_RBP,metadata,by.x="Tumor_Sample_Barcode",by.y="Sample")

#merge with DEGs RBPs
test<-somatic_target_RBP[somatic_target_RBP$Hugo_Symbol %in% unique(DEGs_cohorts_RBPs$NAME),]
test<-data.frame(table(test$Hugo_Symbol))
test<-test[test$Freq>0,]

ggplot(test,aes(x=Freq))+
  geom_histogram(bins=30,fill="coral")+
  geom_density(color="coral")+
  theme_bw()

mutfreq<-data.frame(table(somatic_target_RBP$Tumor_Sample_Barcode))
colnames(mutfreq)<-c("Tumor_Sample_Barcode","Freq")
somatic_target_RBP<-merge(somatic_target_RBP,mutfreq,by="Tumor_Sample_Barcode")

mutRBPfreq<-data.frame(table(somatic_target_RBP[!is.na(somatic_target_RBP$refs),]$Tumor_Sample_Barcode))
colnames(mutRBPfreq)<-c("Tumor_Sample_Barcode","FreqRBPs")
somatic_target_RBP<-merge(somatic_target_RBP,mutRBPfreq,by="Tumor_Sample_Barcode",all.x="TRUE")

mutDEGRBPfreq<-data.frame(table(somatic_target_RBP$Tumor_Sample_Barcode))
colnames(mutDEGRBPfreq)<-c("Tumor_Sample_Barcode","FreqDEGRBPs")
somatic_target_RBP<-merge(somatic_target_RBP,mutDEGRBPfreq,by="Tumor_Sample_Barcode",all.x="TRUE")


melt_target_RBP<-subset(somatic_target_RBP,select=c("Tumor_Sample_Barcode","Freq","FreqRBPs","FreqDEGRBPs","First_event","Pheno"))
melt_target_RBP<-melt_target_RBP[!duplicated(melt_target_RBP$Tumor_Sample_Barcode), ]
melt_target_RBP<-melt(melt_target_RBP,id.vars=c("Tumor_Sample_Barcode","First_event","Pheno"))


# ORDER samples by mutations
melt_target_RBP$value<-as.numeric(melt_target_RBP$value)
melt_target_RBP$Tumor_Sample_Barcode <- factor(melt_target_RBP$Tumor_Sample_Barcode, levels = melt_target_RBP$Tumor_Sample_Barcode[rev(order(melt_target_RBP[melt_target_RBP$variable=="Freq",]$value))])

ggplot(melt_target_RBP,aes(x=Tumor_Sample_Barcode,y=value,fill=variable))+
  geom_bar(data=(melt_target_RBP[melt_target_RBP$variable=="Freq",]),stat="identity",fill="lavender",color="black")+
  geom_bar(data=(melt_target_RBP[melt_target_RBP$variable=="FreqRBPs",]),stat="identity",fill="dodgerblue",color="black")+
  geom_bar(data=(melt_target_RBP[melt_target_RBP$variable=="FreqDEGRBPs",]),stat="identity",fill="pink",color="black")+
  facet_grid(~First_event,space="free",scales = "free")+
  theme_light()+
  theme(axis.text.x = element_blank(),strip.text.x = element_text(size = 14,angle = 90))
  



#####

## alternative splicing analysis INCLUSION TABLES
#####

# To create select file: awk -F '\t' '{for (i=7; i<=NF; i+=2) printf "%s\t", $i; printf "\n"}' INCLUSION_LEVELS_FULL-Hsa264-hg38.tab
difAS<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/INCLUSION_LEVELS_select_gene_event_Hsa389-hg38.tab",sep="\t")
#difAS<-read.delim("/Volumes/cancer/TCell/vast_output/INCLUSION_LEVELS_select_Hsa254-hg38.tab",sep=" ")
difAS$X<-NULL
dim(difAS)

#gene_event<-read.delim("/Volumes/cancer/TCell/vast_output/gene_event_info_Hsa264.tab")
#gene_event<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/gene_event_info_Hsa.tab")
#dim(gene_event)

#unannot<-read.delim("/Volumes/cancer/TCell/vast_output/unannot_vast.txt")
#unannot<-read.delim("/Users/yguillen/Desktop/temp/TCell/vast_output/unannot_vast.txt",sep="\t")

#gene_event<-merge(gene_event,unannot,by="EVENT",all.x=TRUE)
#gene_event$GENE<-paste(gene_event$GENE.x,gene_event$GENE.y)
#gene_event$GENE<-gsub(' NA','',gene_event$GENE)
#gene_event<-gene_event[,c(1,8,3,4,5,6)]

#dim(gene_event)

#difAS<-cbind(gene_event,difAS)

#rm(gene_event)
#rm(unannot)

#difAS_bcat<-read.delim("/Volumes/cancer/Gekas_RNAseq/vast_output/INCLUSION_LEVELS_FULL-Hsa4-hg38.tab")
difAS_bcat<-read.delim("/Users/yguillen/Desktop/temp/Gekas_RNASeq/vast_output/INCLUSION_LEVELS_FULL-Hsa4-hg38.tab")

# Select only values with ratios of inclusion
#difAS_prim<-difAS[,c(1,2,6,seq(7, ncol(difAS), by=2))]

difAS_bcat<-difAS_bcat[,c(1:6,seq(7, ncol(difAS_bcat), by=2))]

difAS_all<-cbind(difAS,difAS_bcat[,7:ncol(difAS_bcat)])
rm(difAS)

## Discard events that min and max are <10 and >90
selev<-data.frame("GENE"=difAS_all$GENE,
           "EVENT"=difAS_all$EVENT,
           "maxv"=apply(difAS_all[,7:ncol(difAS_all)], 1 , function(x) as.numeric(max(x,na.rm = TRUE))),
           "minv"=apply(difAS_all[,7:ncol(difAS_all)], 1 , function(x) as.numeric(min(x,na.rm = TRUE))))

#selev_s<-selev[is.finite(selev$maxv) & is.finite(selev$minv) & (selev$maxv-selev$minv>10 ),]
discard<-row.names(selev[selev$maxv<10 | selev$minv>90 | (is.na(selev$minv) & is.na(selev$maxv)),])
selev_s<-difAS_all[!row.names(difAS_all) %in% discard,]
rm(discard)

#MELT
selev_s<-selev_s[,-c(3,4,5)]

difASm <- melt(selev_s,id.vars=c("EVENT","COMPLEX","GENE"))

dim(difASm)
class(difASm$value)
class(difASm$variable)

difASm$variable<-gsub('EGAR[0-9]*_','',difASm$variable)
difASm$variable<-gsub('.sra_1','',difASm$variable)
difASm$variable<-gsub('_1','',difASm$variable)

table(difASm$variable)

difASm<-merge(metadata,difASm,by.x="SampID",by.y="variable")


difASm$value<-as.numeric(difASm$value)

dim(difASm)

difASm$Pheno<-factor(difASm$Pheno,levels = c("BM","Thymus","pre-T-ALL","IMM","c-T-ALL","Cortical","c-T-LBL","Pre-cortical","Post-cortical",
                                             "mature","NA","control","shbcat"))


difASm$Source<-as.factor(difASm$Source)

difASm$Source<-factor(difASm$Source,levels = c("HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8","Thymus","TALL","TLBL","Jurkat","RPMI"))

difASm$Type<-factor(difASm$Type,levels = c("BM","Thymus","TALL","TLBL"))

ggplot(data=difASm[difASm$EVENT=="HsaINT0050180" & difASm$Type!="TLBL",],aes(x=Type,y=value))+
  geom_point(color="black",size=6)+
  geom_point(data=difASm[difASm$EVENT=="HsaINT0050180" & difASm$GSEA=="line" & difASm$Type!="TLBL",],aes(color=Pheno),size=5)+
  scale_color_brewer(palette="Set2")+
  new_scale_color()+
  geom_point(aes(color=Type),size=4)+
  geom_boxplot(data=difASm[difASm$EVENT=="HsaINT0050180" & difASm$Type!="TLBL",],alpha=0)+
  #geom_point(data=difASm[difASm$GSEA=="GSE69239" & difASm$EVENT=="HsaINT0030865",],aes(shape=Source),size=5)+
  scale_shape_manual(values=c(0,8,2,5,6,7,10,11,12,9))+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(~GSEA,scales="free_x",nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))



## For Tcell development

ggplot(data=difASm[difASm$EVENT=="HsaINT0159248" & difASm$GSEA=="GSE69239",],aes(x=Source,y=value))+
  geom_point(color="black",size=6)+
  geom_point(data=difASm[difASm$EVENT=="HsaINT0159248" & difASm$GSEA=="GSE69239",],aes(color=Type),size=5)+
  geom_ribbon(aes(group=GSEA,ymin=0,ymax=value),color="black",alpha=0.1)+
  scale_color_brewer(palette="Set2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))



## To test differences

ggplot(data=difASm[difASm$EVENT=="HsaINT0159248" & difASm$GSE!="line",],aes(x=GSEA,y=value))+
  geom_point(color="black",size=6)+
  geom_point(data=difASm[difASm$EVENT=="HsaINT0159248" & difASm$GSEA!="line",],aes(color=Type),size=5)+
  scale_color_brewer(palette="Set2")+
  geom_boxplot(data=difASm[difASm$EVENT=="HsaINT0159248" & difASm$GSEA!="line",],alpha=0)+
  scale_color_brewer(palette="Spectral")+
  #facet_wrap(~GSEA,scales="free_x",nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

##

library(ggpubr)

my_comparisons1 <- list(c("TALL", "Thymus"))
my_comparisons2 <- list(c("TALL", "BM"))

ggboxplot(difASm[difASm$EVENT=="HsaINT0159248",], x = "Type", y = "value",
          color = "Type", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons1)+
  stat_compare_means(label.y = 100)   


table(difASm$COMPLEX)

difASm$comp_type<-difASm$COMPLEX
difASm$comp_type<-gsub('C1','Exon_skipping',difASm$comp_type)
difASm$comp_type<-gsub('C2','Exon_skipping',difASm$comp_type)
difASm$comp_type<-gsub('C3','Exon_skipping',difASm$comp_type)
difASm$comp_type<-gsub('ANN','Exon_skipping',difASm$comp_type)
difASm$comp_type<-gsub('IR-S','Intron_retention',difASm$comp_type)
difASm$comp_type<-gsub('IR-C','Intron_retention',difASm$comp_type)
difASm$comp_type<-gsub('S','Exon_skipping',difASm$comp_type)


table(difASm$comp_type)


# --------- #
## Changing na by NA
#difAsm_NA<-difASm
#difAsm_NA$comp_type[is.na(difAsm_NA$value)] <- "NA"

# Remove NA
difASm_nona_thr<-difASm[!is.na(difASm$value),]

#threvent<-as.character(unique(difASm_nona_thr$EVENT))

# ADD info cancer driver genes
unique(candriver$SYMBOL)

difASm_nona_thr$driver<-ifelse(difASm_nona_thr$GENE %in% unique(candriver$SYMBOL), 
                               c("driver"), c("non-driver")) 


# selecting or not

samptab<-as.data.frame(table(difASm_nona_thr$comp_type,difASm_nona_thr$SampID))
samptab<-merge(samptab,metadata,by.x="Var2",by.y="SampID")

ggplot(samptab[samptab$GSE!="EGA_TLE",],aes(x=Reads,y=Freq))+
  #geom_bar(stat="identity", aes(fill=Var1),width = 0.7) +
  geom_point(aes(color=Var1))+
  geom_label_repel(data=samptab[samptab$GSE!="EGA_TLE" & samptab$Reads<75000000,],aes(label=Sample),size=2)+
  scale_color_brewer(palette="Paired")+
  geom_vline(xintercept=70000000,linetype ="dashed")+
  facet_wrap(~Var1,scales="free_y",ncol=3)+
  theme_bw()

ggplot(samptab[samptab$GSE!="EGA_TLE",],aes(x=GSEA,y=Freq))+
  #geom_bar(stat="identity", aes(fill=Var1),width = 0.7) +
  geom_point(aes(color=Var1))+
  geom_boxplot(alpha=0)+
  #geom_smooth(method="lm",aes(color=Var1))+
  scale_color_brewer(palette="Paired")+
  facet_wrap(~Var1,scales="free_y",ncol=3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

# REMOVING EGA_TLE cohort
samptab<-samptab[samptab$GSEA!="EGA_TLE",]

cor.test(samptab[samptab$Var1=="Alt3" & samptab$Reads>70000000,]$Freq,samptab[samptab$Var1=="Alt3" & samptab$Reads>70000000,]$Reads)
kruskal.test(samptab[samptab$Var1=="Alt3",]$GSEA,samptab[samptab$Var1=="Alt3" & samptab$GSE!="EGA_TLE",]$Freq)

cor.test(samptab[samptab$Var1=="Alt5" & samptab$Reads>70000000,]$Freq,samptab[samptab$Var1=="Alt5" & samptab$Reads>70000000,]$Reads)
kruskal.test(samptab[samptab$Var1=="Alt5",]$GSEA,samptab[samptab$Var1=="Alt5",]$Freq)

cor.test(samptab[samptab$Var1=="Exon_skipping" & samptab$Reads>70000000,]$Freq,samptab[samptab$Var1=="Exon_skipping" & samptab$Reads>70000000,]$Reads)
kruskal.test(samptab[samptab$Var1=="Exon_skipping",]$GSEA,samptab[samptab$Var1=="Exon_skipping",]$Freq)

cor.test(samptab[samptab$Var1=="Intron_retention" & samptab$Reads>70000000,]$Freq,samptab[samptab$Var1=="Intron_retention" & samptab$Reads>70000000,]$Reads)
kruskal.test(samptab[samptab$Var1=="Intron_retention",]$GSEA,samptab[samptab$Var1=="Intron_retention",]$Freq)

cor.test(samptab[samptab$Var1=="MIC" & samptab$Reads>70000000,]$Freq,samptab[samptab$Var1=="MIC" & samptab$Reads>70000000,]$Reads)
kruskal.test(samptab[samptab$Var1=="MIC",]$GSEA,samptab[samptab$Var1=="MIC",]$Freq)

tab_eve<-as.data.frame(table(difASm_nona_thr$comp_type,difASm_nona_thr$driver,difASm_nona_thr$SampID))
tab_eve_met<-merge(metadata,tab_eve,by.x="SampID",by.y="Var3")
tab_eve_met$Type<-factor(tab_eve_met$Type,levels = c("Thymus","TALL","TLBL","BM"))

ggplot(data=tab_eve_met[tab_eve_met$Reads>70000000,],aes(x=Type,y=Freq))+
  geom_point(color="black",size=6)+
  geom_point(aes(color=Type),size=5)+
  geom_boxplot(alpha=0)+
  scale_color_brewer(palette="Spectral")+
  facet_grid(Var1~GSEA,scales="free",space="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

# REMOVE EGA_TLE
tab_eve_met<-tab_eve_met[tab_eve_met$GSEA!="EGA_TLE",]

tab_eve_met_total<-as.data.frame(tab_eve_met %>% 
  group_by(SampID) %>% 
  summarise(amt = sum(Freq)))

sumtab<-merge(tab_eve_met_total,metadata,by="SampID")
sumtab$amt<-as.numeric(sumtab$amt)

sumtab$Type<-factor(sumtab$Type,levels = c("Thymus","BM","TALL","TLBL"))
sumtab$SampID<- factor(sumtab$SampID, levels = unique(sumtab$SampID)[order(sumtab$amt, decreasing = TRUE)])


ggplot(data=sumtab,aes(x=Type,y=amt))+
  geom_point(color="black",size=6)+
  geom_point(aes(color=Type),size=5)+
  geom_boxplot(alpha=0)+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(~GSEA,scales="free_x",nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

wilcox.test(sumtab[sumtab$GSE=="GSE109231",]$amt~sumtab[sumtab$GSE=="GSE109231",]$Type)
wilcox.test(sumtab[sumtab$GSE=="GSE57982",]$amt~sumtab[sumtab$GSE=="GSE57982",]$Type)
wilcox.test(sumtab[sumtab$GSE=="GSE69239",]$amt~sumtab[sumtab$GSE=="GSE69239",]$Type)


tab_eve_met<-merge(tab_eve_met,tab_eve_met_total,by="SampID")
tab_eve_met$prop<-(tab_eve_met$Freq/tab_eve_met$amt)*100

## barplot

tab_eve_met$SampID <- factor(tab_eve_met$SampID, levels=unique(tab_eve_met[rev(order(tab_eve_met$prop)), ]$SampID))

ggplot(tab_eve_met,aes(x=SampID,y=prop))+
  #geom_bar(data=(tab_eve_met[tab_eve_met$Type=="Thymus" | tab_eve_met$Type=="BM",]),stat="identity", fill="black",width = 1.5) +
  geom_bar(stat="identity", aes(fill=Var1),width = 0.7) +
  scale_fill_brewer(palette="Paired")+
  facet_grid(~GSEA+Type,scales = "free",space = "free")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 7, angle = 90))


#Info cancer driver genes
subtab<-subset(tab_eve_met,select=c(SampID,Var1,Freq))
subtab<-as.data.frame (subtab %>%
                        group_by(SampID,Var1) %>%
                        summarise(tot = sum(Freq)))

tab_eve_met_driv_tot<-dcast(subtab, SampID ~ Var1)
tab_eve_met_driv<-merge(tab_eve_met,tab_eve_met_driv_tot,by="SampID")


## Number of alternative splicing events per group across Tcell development

tab_eve_met_driv$Source<-factor(tab_eve_met_driv$Source,levels = c("HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8","Thymus","TALL","RPMI","Jurkat"))
tab_eve_met_driv$Var2<-factor(tab_eve_met_driv$Var2,levels=c("non-driver","driver"))

ggplot(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" | tab_eve_met_driv$Source=="Thymus" & tab_eve_met_driv$Type=="Thymus",],aes(x=Source,y=prop))+
  geom_point(color="black",size=4)+
  geom_ribbon(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="BM",],aes(group=Var1,ymin=0,ymax=prop),color="darkolivegreen",size=1,alpha=0)+
  geom_ribbon(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & (tab_eve_met_driv$Type=="Thymus" | tab_eve_met_driv$Source=="CLP"),],aes(group=Var1,ymin=0,ymax=prop),color="orange",size=1.5,alpha=0)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="BM" & tab_eve_met_driv$Source!="CLP",],color="darkolivegreen",size=3)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="Thymus",],color="orange",size=3)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Source=="CLP",],color="white",size=3)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$Source=="Thymus" & tab_eve_met_driv$Type=="Thymus",],aes(color=GSEA),size=3)+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(Var1~Var2,scales="free",ncol=2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        strip.text = element_text(size=9),
        legend.position = "bottom")


ggplot(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" | tab_eve_met_driv$Source=="Thymus" & tab_eve_met_driv$Type=="Thymus",],aes(x=Source,y=Freq))+
  geom_point(color="black",size=3)+
  geom_ribbon(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="BM",],aes(group=Var1,ymin=0,ymax=Freq),color="darkolivegreen",size=1,alpha=0)+
  geom_ribbon(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & (tab_eve_met_driv$Type=="Thymus" | tab_eve_met_driv$Source=="CLP"),],aes(group=Var1,ymin=0,ymax=Freq),color="orange",size=1.5,alpha=0)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="BM" & tab_eve_met_driv$Source!="CLP",],color="darkolivegreen",size=2)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="Thymus",],color="orange",size=2)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Source=="CLP",],color="white",size=2)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$Source=="Thymus" & tab_eve_met_driv$Type=="Thymus",],aes(color=GSEA),size=2)+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(Var1~Var2,scales="free",ncol=2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        strip.text = element_text(size=9),
        legend.position = "bottom")+
  coord_flip()

ggplot(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" | tab_eve_met_driv$Source=="Thymus" & tab_eve_met_driv$Type=="Thymus",],aes(x=Source,y=prop))+
  geom_point(color="black",size=3)+
  geom_ribbon(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="BM",],aes(group=Var1,ymin=0,ymax=prop),color="darkolivegreen",size=1,alpha=0)+
  geom_ribbon(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & (tab_eve_met_driv$Type=="Thymus" | tab_eve_met_driv$Source=="CLP"),],aes(group=Var1,ymin=0,ymax=prop),color="orange",size=1.5,alpha=0)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="BM" & tab_eve_met_driv$Source!="CLP",],color="darkolivegreen",size=2)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Type=="Thymus",],color="orange",size=2)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$GSEA=="GSE69239" & tab_eve_met_driv$Source=="CLP",],color="white",size=2)+
  geom_point(data=tab_eve_met_driv[tab_eve_met_driv$Source=="Thymus" & tab_eve_met_driv$Type=="Thymus",],aes(color=GSEA),size=2)+
  scale_color_brewer(palette="Spectral")+
  facet_wrap(Var1~Var2,scales="free",ncol=2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust = 1,size=7),
        axis.text.y = element_text(size=7),
        strip.text = element_text(size=7),
        legend.position = "bottom")+
  coord_flip()



#For Alt3
tab_eve_Alt3<-subset(tab_eve_met_driv,tab_eve_met_driv$Var1=="Alt3")
tab_eve_Alt3$amount<-(tab_eve_Alt3$Freq/tab_eve_Alt3$Alt3)*100
#For Alt5
tab_eve_Alt5<-subset(tab_eve_met_driv,tab_eve_met_driv$Var1=="Alt5")
tab_eve_Alt5$amount<-(tab_eve_Alt5$Freq/tab_eve_Alt5$Alt5)*100
#For Exon_skipping
tab_eve_Exon_skipping<-subset(tab_eve_met_driv,tab_eve_met_driv$Var1=="Exon_skipping")
tab_eve_Exon_skipping$amount<-(tab_eve_Exon_skipping$Freq/tab_eve_Alt5$Exon_skipping)*100
#For Intron_retention
tab_eve_Intron_retention<-subset(tab_eve_met_driv,tab_eve_met_driv$Var1=="Intron_retention")
tab_eve_Intron_retention$amount<-(tab_eve_Intron_retention$Freq/tab_eve_Intron_retention$Intron_retention)*100
#For MIC
tab_eve_MIC<-subset(tab_eve_met_driv,tab_eve_met_driv$Var1=="MIC")
tab_eve_MIC$amount<-(tab_eve_MIC$Freq/tab_eve_MIC$MIC)*100

eventstab<-rbind(tab_eve_Alt3,tab_eve_Alt5,tab_eve_Exon_skipping,tab_eve_Intron_retention,tab_eve_MIC)

ggplot(eventstab[eventstab$Var1=="Intron_retention",],aes(x=SampID,y=amount))+
  #geom_bar(data=(tab_eve_met[tab_eve_met$Type=="Thymus" | tab_eve_met$Type=="BM",]),stat="identity", fill="black",width = 1.5) +
  geom_bar(stat="identity", aes(fill=Var2),width = 0.7) +
  scale_fill_manual(values=c("coral","darkolivegreen"))+
  facet_grid(~GSEA+Type,scales = "free",space = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 7, angle = 90))


eventstab$SampID <- factor(eventstab$SampID, levels=unique(eventstab[rev(order(eventstab$Var1,eventstab$amount)), ]$SampID))
eventstab$SampID<- factor(eventstab$SampID, levels = unique(eventstab$SampID)[order(eventstab$Var1,eventstab$amount, decreasing = TRUE)])

ggplot(eventstab[eventstab$Var1=="Exon_skipping",],aes(x=SampID,y=amount))+
  #geom_bar(data=(tab_eve_met[tab_eve_met$Type=="Thymus" | tab_eve_met$Type=="BM",]),stat="identity", fill="black",width = 1.5) +
  geom_bar(stat="identity", aes(fill=Var2),width = 0.7) +
  scale_fill_manual(values=c("coral","darkolivegreen"))+
  facet_grid(~GSEA+Type,scales = "free",space = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 7, angle = 90))

ggplot(eventstab[eventstab$Var1=="MIC",],aes(x=SampID,y=amount))+
  #geom_bar(data=(tab_eve_met[tab_eve_met$Type=="Thymus" | tab_eve_met$Type=="BM",]),stat="identity", fill="black",width = 1.5) +
  geom_bar(stat="identity", aes(fill=Var2),width = 0.7) +
  scale_fill_manual(values=c("coral","darkolivegreen"))+
  facet_grid(~GSEA+Type,scales = "free",space = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 7, angle = 90))


tab_eve_met$Type<-factor(tab_eve_met$Type,levels = c("BM","Thymus","TALL","TLBL"))
tab_eve_met$Var2<-factor(tab_eve_met$Var2,levels=c("non-driver","driver"))


ggplot(tab_eve_met,aes(x=Type,y=prop))+
  geom_point(color="black",size=4)+
  geom_point(aes(color=Type),size=3)+
  geom_boxplot(aes(fill=Type),alpha=0.5) +
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  facet_wrap(Var2~Var1,scales = "free",ncol=5)+
  theme_bw()


tab_eve_met$state<-ifelse(tab_eve_met$Type=="TALL" | tab_eve_met$Type=="TLBL", 
                       c("leukemia"), c("normal")) 


ggplot(tab_eve_met,aes(x=state,y=prop))+
  geom_point(color="black",size=4)+
  geom_point(aes(color=Var1),size=3)+
  geom_violin(aes(fill=Var1),alpha=0.2) +
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  facet_wrap(Var2~Var1,scales = "free",ncol=5)+
  theme_bw()


ggplot(tab_eve_met,aes(x=Type,y=prop))+
  geom_point(color="black",size=4)+
  geom_point(aes(color=Var1),size=3)+
  geom_violin(aes(fill=Var1),alpha=0.2) +
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  facet_wrap(Var2~Var1,scales = "free",ncol=5)+
  theme_bw()

my_comparisons <- list( c("TALL", "Thymus"), c("TLBL", "Thymus"), c("TALL", "BM"),c("TLBL", "BM") )

ggboxplot(tab_eve_met[tab_eve_met$Var1=="Exon_skipping" & tab_eve_met$Var2=="driver",], x = "Type", y = "prop",
          color = "Type", palette = "jco")+ 
  #stat_compare_means(method="wilcox.test",comparisons = my_comparisons, label.y = c(45,46,47,48))+
  stat_compare_means(method="wilcox.test",comparisons = my_comparisons, label.y = c(1.9,2.1,2.3,2.4))+
  ylim(c(1.8,2.4))
  #ylim(c(42,52))



#BiocManager::install('pgirmess')
library(pgirmess)

kruskalmc(tab_eve_met[tab_eve_met$Var1=="Alt3" & tab_eve_met$Var2=="non-driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="Alt3" & tab_eve_met$Var2=="non-driver",]$state,probs = 0.05)
kruskalmc(tab_eve_met[tab_eve_met$Var1=="Alt5" & tab_eve_met$Var2=="non-driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="Alt5" & tab_eve_met$Var2=="non-driver",]$state,probs = 0.05)
kruskalmc(tab_eve_met[tab_eve_met$Var1=="Exon_skipping" & tab_eve_met$Var2=="non-driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="Exon_skipping" & tab_eve_met$Var2=="non-driver",]$state,probs = 0.05)
kruskalmc(tab_eve_met[tab_eve_met$Var1=="Intron_retention" & tab_eve_met$Var2=="non-driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="Intron_retention" & tab_eve_met$Var2=="non-driver",]$state,probs = 0.05)
kruskalmc(tab_eve_met[tab_eve_met$Var1=="MIC" & tab_eve_met$Var2=="non-driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="MIC" & tab_eve_met$Var2=="non-driver",]$state,probs = 0.05)


kruskalmc(tab_eve_met[tab_eve_met$Var1=="Alt3" & tab_eve_met$Var2=="driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="Alt3" & tab_eve_met$Var2=="driver",]$state,probs = 0.05)
kruskalmc(tab_eve_met[tab_eve_met$Var1=="Alt5" & tab_eve_met$Var2=="driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="Alt5" & tab_eve_met$Var2=="driver",]$state,probs = 0.05)
kruskalmc(tab_eve_met[tab_eve_met$Var1=="Exon_skipping" & tab_eve_met$Var2=="driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="Exon_skipping" & tab_eve_met$Var2=="driver",]$state,probs = 0.05)
kruskalmc(tab_eve_met[tab_eve_met$Var1=="Intron_retention" & tab_eve_met$Var2=="driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="Intron_retention" & tab_eve_met$Var2=="driver",]$state,probs = 0.05)
kruskalmc(tab_eve_met[tab_eve_met$Var1=="MIC" & tab_eve_met$Var2=="driver",]$prop ~ tab_eve_met[tab_eve_met$Var1=="MIC" & tab_eve_met$Var2=="driver",]$state,probs = 0.05)

# Classification per first event
ggplot(tab_eve_met[tab_eve_met$First_event=="Remission" | tab_eve_met$First_event=="Relapse",],aes(x=First_event,y=prop))+
  #geom_point(color="black",size=4)+
  geom_point(aes(color=Var1,size=Event_free_surv_days))+
  geom_violin(aes(fill=Var1),alpha=0.2) +
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  facet_wrap(Var2~Var1,scales = "free",ncol=5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,size=7))

my_comparisons <- list( c("Relapse", "Remission"))

ggboxplot(tab_eve_met[tab_eve_met$Var1=="Intron_retention" & tab_eve_met$Var2=="non-driver",], x = "First_event", y = "prop",
          color = "First_event", palette = "jco")+ 
  stat_compare_means(method="wilcox.test",comparisons = my_comparisons, label.y = c(40))+
  #stat_compare_means(method="wilcox.test",comparisons = my_comparisons, label.y = c(1,1.2,1.3,1.4))+
  #ylim(c(1,2))
  ylim(c(40,50))

# Number of events per gene
# Selecting gene
ggplot(difASm_nona_thr[difASm_nona_thr$GENE=="AFF1",],aes(x=SampID,y=value))+
  geom_point(aes(color=EVENT,shape=Type))+
  facet_grid(~GSEA,scales="free",space="free")+
  theme_classic()

# Distribution of PSI events and samples
distevent<-as.data.frame(table(difASm_nona_thr$EVENT))
distevent<-distevent[distevent$Freq!=0,]
distevent$prop<-distevent$Freq/max(distevent$Freq)
hist(distevent$prop)

distevent<-merge(distevent,difAS_all[,c(2,7:ncol(difAS_all))],by.x="Var1",by.y="EVENT",all.x="TRUE")

# Number of samples for which more than 65% of AS events had null PSI values
na_count <-data.frame(nullPSI=sapply(distevent[,4:ncol(distevent)], function(x) sum(length(which(is.na(x))))),
                      TotalAS=nrow(distevent))
na_count$freqNULLPSI<-na_count$nullPSI/na_count$TotalAS*100
thresh_AS<-row.names(na_count[na_count$freqNULLPSI<65,])


#### PCA with alternative splicing events

#fOR PCA complete.cases
row.names(selev_s)<-selev_s$EVENT
selev_s<-selev_s[,4:ncol(selev_s)]

# Exclude samples with the higher number of unreported splicing events = R3
splpca<-selev_s[,(colnames(selev_s) %in% thresh_AS) ]

# EXCLUDE TLE COHORT
splpca<-splpca[,!grepl('TLE',colnames(splpca)) & !grepl('Thymus',colnames(splpca)) & !grepl('R7',colnames(splpca)) & !grepl('R8',colnames(splpca)) & !grepl('R9',colnames(splpca)) & !grepl('R10',colnames(splpca))]

splpca<-splpca[complete.cases(splpca),]
which(is.na(splpca))



splpca_t<-t(splpca)

#DSmat_splpca <- sapply(splpca_t, function(x) as.numeric(x))
#dim(DSmat_splpca)
#class(DSmat_splpca)

DSmat_splpca<-splpca_t
row.names(DSmat_splpca)<-row.names(splpca_t)
rm(splpca_t)
rm(splpca)

row.names(DSmat_splpca)<-gsub('_1','',row.names(DSmat_splpca))
row.names(DSmat_splpca)<-gsub('.sra','',row.names(DSmat_splpca))
row.names(DSmat_splpca)<-gsub('EGAR.*_','',row.names(DSmat_splpca))


#select those with variance != 0 (error in prcomp otherwise)
cond<-(apply(DSmat_splpca, 2, var)>=20)
DSmat_splpca<-DSmat_splpca[, cond, drop = FALSE]

### Following vast-tools paper, select only PSI with SD > 20
#DSmat_splpca<-DSmat_splpca[, sapply(DSmat_splpca, function(x) { sd(x) > 20} )]


res.splTLE<-prcomp(DSmat_splpca,scale=TRUE)
(res.splTLE$sdev^2/sum(res.splTLE$sdev^2))[1:5]


#Proportion of variance

pcatab_spl<-as.data.frame(res.splTLE$x)
pcatab_spl$names<-row.names(res.splTLE)
row.names(pcatab_spl)

pcatab_spl<-pcatab_spl[,c(1:3)]
row.names(pcatab_spl)

metadata$SampID

setdiff(metadata$SampID,row.names(pcatab_spl))

pcatab_spl<-merge(metadata,pcatab_spl,by.x="SampID",by.y="row.names",all.y=TRUE)
pcatab_spl$var<-paste(pcatab_spl$Source,pcatab_spl$Pheno,sep="_")

library(ggnewscale)

pcatab_spl$GSEA<-as.factor(pcatab_spl$GSEA)
pcatab_spl$Age_Dx_days<-as.numeric(pcatab_spl$Age_Dx_days)

ggplot(pcatab_spl, aes(PC1, PC2,label=SampID)) +
  #geom_point(data=pcatab[pcatab$source=="BM",],size=10,color="black") +
  #geom_point(data=pcatab[pcatab$source=="Thymus",],size=10,color="grey") +
  #geom_point(data=pcatab[pcatab$source=="TALL",],size=10,color="red") +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=GSEA)) +
  #geom_point(data=pcatab[pcatab$GSEA=="GSE57982",],size=6,color="black") +
  scale_color_brewer(palette = "Blues")+
  #new_scale_color()+
  #stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
  #scale_color_brewer(palette = "Blues")+
  new_scale_color()+
  geom_point(size=5,color="white")+
  geom_point(data=pcatab_spl[pcatab_spl$First_event!="NA",],size=5,aes(color=First_event))+
  scale_color_manual(values=c("grey","black","brown1","coral4","coral","darkolivegreen4","goldenrod"))+
  #scale_linetype_manual(values=c("solid","dotted"))+
  geom_label_repel(data=pcatab_spl[pcatab_spl$GSEA=="GSE69239",],aes(label=Source), size=4,fontface = "italic",alpha=0.7)+
  geom_label_repel(data=pcatab_spl[pcatab_spl$Source=="Thymus" ,],aes(label=Source), size=4,fontface = "italic",alpha=0.7,color="grey")+
  #  geom_label_repel(data=pcatab_TLE[pcatab_TLE$GSEA=="EGA_TLE" ,],aes(label=SampID), size=4,fontface = "italic",alpha=0.7,color="grey")+
  #scale_color_manual(values=c("purple","red"))+
  #ylim(-50,50)+
  ylab("PC2 (6.207% variance)")+
  xlab("PC1 (27.50% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()



ggplot(pcatab_spl, aes(PC1, PC2,label=SampID)) +
  #geom_point(data=pcatab[pcatab$source=="BM",],size=10,color="black") +
  #geom_point(data=pcatab[pcatab$source=="Thymus",],size=10,color="grey") +
  #geom_point(data=pcatab[pcatab$source=="TALL",],size=10,color="red") +
  geom_point(size=8,color="black")+
 # geom_point(size=7,aes(color=GSEA)) +
  #geom_point(data=pcatab[pcatab$GSEA=="GSE57982",],size=6,color="black") +
  #scale_color_brewer(palette = "Blues")+
  #new_scale_color()+
  #stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
  #scale_color_brewer(palette = "Blues")+
  #new_scale_color()+
  geom_point(size=5,color="white")+
  geom_point(data=pcatab_spl[pcatab_spl$Pheno!="shbcat" & pcatab_spl$Pheno!="NA" & pcatab_spl$Pheno!="control" & pcatab_spl$Pheno!="BM",],size=5,aes(color=Pheno))+
  scale_color_manual(values=c("darkolivegreen2","darkolivegreen3","darkolivegreen1",
                              "brown1",
                              "dodgerblue","dodgerblue4",
                              "brown4","brown3",
                              "grey"))+
  #scale_linetype_manual(values=c("solid","dotted"))+
  geom_label_repel(data=pcatab_spl[pcatab_spl$GSEA=="GSE69239",],aes(label=Source), size=4,fontface = "italic",alpha=0.7)+
  geom_label_repel(data=pcatab_spl[pcatab_spl$Source=="Thymus" ,],aes(label=Source), size=4,fontface = "italic",alpha=0.7,color="grey")+
  #  geom_label_repel(data=pcatab_TLE[pcatab_TLE$GSEA=="EGA_TLE" ,],aes(label=SampID), size=4,fontface = "italic",alpha=0.7,color="grey")+
  #scale_color_manual(values=c("purple","red"))+
  #ylim(-50,50)+
  ylab("PC2 (6.20% variance)")+
  xlab("PC1 (27.50% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()

#####

## Surival curves with AS events
#####
# Survival with CD47 total levels and HsaEX0013876
cd47_gene<-melt_vastout[melt_vastout$variable=="ENSG00000196776",]
cd47_event<-as.data.frame(t(difAS_all[difAS_all$GENE=="CD47",]))

surv_cd47_gene<-subset(cd47_gene,select=c(GSM,First_event,Event_free_surv_days,Age_Dx_days,value))
surv_cd47_gene<-surv_cd47_gene[surv_cd47_gene$First_event!="Censored" &
                                 surv_cd47_gene$First_event!="Second Malignant Neoplasm" &
                                 !is.na(surv_cd47_gene$First_event),]

boxplot(surv_cd47_gene$value)

surv_cd47_gene<-surv_cd47_gene %>% mutate(CD47_group =  ifelse(value >=200, "CD47_high", "CD47_low"))
surv_cd47_gene<-surv_cd47_gene %>% mutate(disease_free =  ifelse(First_event == "Remission", 0, 1))

surv_cd47_gene$CD47_group<- factor(surv_cd47_gene$CD47_group)

row.names(surv_cd47_gene)<-surv_cd47_gene$GSM

surv_object <- Surv(time = surv_cd47_gene$Event_free_surv_days/365, event = surv_cd47_gene$disease_free)
fit1 <- survfit(surv_object ~ CD47_group, data = surv_cd47_gene)
summary(fit1)
ggsurvplot(fit1, data = surv_cd47_gene, pval = TRUE)


fit.coxph<- coxph(Surv(Event_free_surv_days/365, disease_free) ~ CD47_group, 
                  data = surv_cd47_gene,
                  ties = 'exact')
summary(fit.coxph)

ggforest(fit.coxph)


## skipping event

cd47_skip<-difASm[difASm$EVENT=="HsaEX0013876",]
row.names(cd47_skip)<-cd47_skip$SampID

surv_cd47_skip<-subset(cd47_skip,select=c(SampID,First_event,Event_free_surv_days,Age_Dx_days,value))
surv_cd47_skip<-surv_cd47_skip[surv_cd47_skip$First_event!="Censored" &
                                 surv_cd47_skip$First_event!="Second Malignant Neoplasm" &
                                 !is.na(surv_cd47_skip$First_event) &
                                 surv_cd47_skip$First_event!="Secondary leukemia" &
                                 surv_cd47_skip$First_event!="Disease",]


boxplot(surv_cd47_skip$value)

surv_cd47_skip<-surv_cd47_skip %>% mutate(CD47skip_group =  ifelse(value <=90, "CD47skip_low", "CD47skip_high"))
surv_cd47_skip<-surv_cd47_skip %>% mutate(disease_free =  ifelse(First_event == "Remission", 0, 1))

surv_cd47_skip$CD47skip_group<- factor(surv_cd47_skip$CD47skip_group)

surv_object <- Surv(time = surv_cd47_skip$Event_free_surv_days/365, event = surv_cd47_skip$disease_free)
fit1 <- survfit(surv_object ~ CD47skip_group, data = surv_cd47_skip)
summary(fit1)
ggsurvplot(fit1, data = surv_cd47_skip, pval = TRUE)


fit.coxph<- coxph(Surv(Event_free_surv_days/365, disease_free) ~ CD47skip_group, 
                  data = surv_cd47_skip,
                  ties = 'exact')
summary(fit.coxph)

ggforest(fit.coxph)



# Merging
row.names(surv_cd47_skip)<-surv_cd47_skip$SampID

surv_cd47<-merge(surv_cd47_gene,surv_cd47_skip,by=0)

plot(surv_cd47$value.x,surv_cd47$value.y)


fit.coxph<- coxph(Surv(Event_free_surv_days.x/365, disease_free.x) ~ value.x+value.y, 
                  data = surv_cd47,
                  ties = 'exact')
summary(fit.coxph)

ggforest(fit.coxph)

#####

## Somatic and structural mutations associated to certain alternative splicing events
#####
#Notch mutations and splicing events
Notch_mut<-somatic_target
Notch_mut$Sample<-Notch_mut$Tumor_Sample_Barcode
Notch_mut$Sample<-gsub('-0.*','',Notch_mut$Sample)
Notch_mut<-merge(Notch_mut,metadata,by="Sample")
Notch_mut<-as.data.frame(table(Notch_mut[grepl('NOTCH1',Notch_mut$Hugo_Symbol),]$SampID))
row.names(Notch_mut)<-Notch_mut$Var1
colnames(Notch_mut)<-c('Sample','NOTCH1_freq')

difASm_mut<-merge(difASm,Notch_mut,by="Sample",all.x=TRUE)

ggplot(data=difASm[difASm$EVENT=="HsaINT0046187",],aes(x=Type,y=value))+
  geom_jitter(data=difASm[difASm$EVENT=="HsaINT0046187",],aes(color=Source),size=2)+
#  geom_jitter(data=difASm[difASm$EVENT=="HsaEX0013876" & !(difASm$SampID %in% unique(Notch_mut$Sample)),],color="black",size=5)+
#  geom_jitter(data=difASm[difASm$EVENT=="HsaEX0013876" & difASm$SampID %in% unique(Notch_mut$Sample),],color="red",size=5)+
#  scale_color_brewer(palette="Set2")+
  geom_boxplot(data=difASm[difASm$EVENT=="HsaINT0046187",],alpha=0)+
#  scale_color_brewer(palette="Spectral")+
  facet_wrap(~GSEA,scales="free_x",nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ylim(0,100)



#####



### UNSUPERVISED CLUSTERING OF SAMPLES USING PSI VALUES
#####
# Exclude samples with the higher number of unreported splicing events
splpca<-selev_s[,(colnames(selev_s) %in% thresh_AS) ]

# EXCLUDE TLE COHORT
splpca<-splpca[,!grepl('TLE',colnames(splpca)) & !grepl('Thymus',colnames(splpca)) & !grepl('R7',colnames(splpca)) & !grepl('R8',colnames(splpca)) & !grepl('R9',colnames(splpca)) & !grepl('R10',colnames(splpca))]

splpca<-splpca[complete.cases(splpca),]
which(is.na(splpca))

indx <- sapply(splpca, is.factor)
splpca[indx] <- lapply(splpca[indx], function(x) as.numeric(as.character(x)))

colnames(splpca)<-gsub('_1','',colnames(splpca))
colnames(splpca)<-gsub('.sra','',colnames(splpca))
colnames(splpca)<-gsub('EGAR.*_','',colnames(splpca))

row.names(DSmat_splpca)<-gsub('_1','',row.names(DSmat_splpca))
row.names(DSmat_splpca)<-gsub('.sra','',row.names(DSmat_splpca))
row.names(DSmat_splpca)<-gsub('EGAR.*_','',row.names(DSmat_splpca))


#select those with variance != 0 (error in prcomp otherwise)
cond<-(apply(DSmat_splpca, 2, var)>=20)
DSmat_splpca<-DSmat_splpca[, cond, drop = FALSE]





#median center genes
dc = sweep(ccp,1, apply(ccp,1,median))

# run consensus cluster, with standard options
rcc = ConsensusClusterPlus(dc,
                           maxK=4,
                           reps=100,
                           pItem=0.8,
                           pFeature=1,
                           title="con_clust_TALL",
                           distance="pearson",
                           clusterAlg="hc",plot='pdf')


rcc[[4]]$consensusClass

#For .example, the top five rows and columns of results for k=2:
rcc[[4]][["consensusMatrix"]][1:10,1:10]

rcc[[4]][["consensusTree"]]

rcc[[4]][["consensusClass"]]

#consensus matrix results
rcc[[4]][["ml"]]


icl = calcICL(rcc,title="itcon_clust_TALL_noTLE",plot='pdf')
icl[["itemConsensus"]][1:4,]



## add to metadata consensus class

cuns<-data.frame(rcc[[4]][["consensusClass"]])
colnames(cuns)<-"uns_clust"
cuns$SampID<-row.names(cuns)

write.table(cuns,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/CEMITools/TALL_cemitool_clusters_uns.tsv",row.names = FALSE,sep="\t",quote = FALSE)


hvgplot  <- t(apply(rld_batch[takebatch,], 1, scale))
colnames(hvgplot)<-colnames(rld_batch[takebatch,])

#only TARGET cohort
hvgplot  <- t(apply(rld_batch, 1, scale))
colnames(hvgplot)<-colnames(rld_batch)


announ <- subset(metadata[metadata$SampID %in% colnames(hvgplot),],select=c(SampID,Pheno,Source,GSEA,First_event))
row.names(announ)<-announ$SampID
announ$SampID<-NULL
table(announ$First_event)

ccp_met<-merge(cuns,announ,by="row.names")
row.names(ccp_met)<-row.names(announ)
ccp_met<-subset(ccp_met,select=c("Source","uns_clust","Pheno","GSEA","First_event"))
ccp_met$uns_clust<-gsub("^","C",ccp_met$uns_clust)


row.names(ccp_met)

library(RColorBrewer)

s1<-length(unique(ccp_met$Source))
s2<-length(unique(ccp_met$uns_clust))
s3<-length(unique(ccp_met$Pheno))
s4<-length(unique(ccp_met$First_event))

myColors2 = list(
  Source = c(TALL=brewer.pal(s1,"Set2")[1],
             TLBL=brewer.pal(s1,"Set2")[2]),
  uns_clust = c(C1=brewer.pal(s2,"BrBG")[1],
                C2=brewer.pal(s2,"BrBG")[2],
                C3=brewer.pal(s2,"BrBG")[3],
                C4=brewer.pal(s2,"BrBG")[4]),
  #                C5=brewer.pal(s2,"BrBG")[5]),
  Pheno = c("NA"=brewer.pal(s3,"Paired")[1],
            "c-T-ALL"=brewer.pal(s3,"Paired")[2],
            "c-T-LBL"=brewer.pal(s3,"Paired")[9],
            mature=brewer.pal(s3,"Paired")[3],
            IMM=brewer.pal(s3,"Paired")[4],
            Cortical=brewer.pal(s3,"Paired")[5],
            "Post-cortical"=brewer.pal(s3,"Paired")[6],
            "Pre-cortical"=brewer.pal(s3,"Paired")[7],
            "pre-T-ALL"=brewer.pal(s3,"Paired")[8]),
  First_event = c("NA"="white",
                  "Censored"=brewer.pal(s4,"Dark2")[1],
                  "Death"=brewer.pal(s4,"Dark2")[2],
                  "Disease"=brewer.pal(s4,"Dark2")[3],
                  "Progression"=brewer.pal(s4,"Dark2")[4],
                  "Refractory"=brewer.pal(s4,"Dark2")[5],
                  "Relapse"=brewer.pal(s4,"Dark2")[6],
                  "Remission"=brewer.pal(s4,"Dark2")[7],
                  "Second Malignant Neoplasm"=brewer.pal(s4,"Dark2")[8],
                  "Secondary leukemia"=brewer.pal(s4,"Dark2")[8]))

#order samples by clust category

hvgplot[, order(ccp_met$uns_clust)]
hvgplot<-hvgplot[,row.names(ccp_met[order(ccp_met$uns_clust),])]

#hvgplot[, order(ccp_met$First_event)]
#hvgplot<-hvgplot[,row.names(ccp_met[order(ccp_met$First_event),])]

pheatplot<-pheatmap(hvgplot,
                    annotation_col = ccp_met,
                    fontsize = 8,
                    cutree_rows = 6,
                    show_rownames = FALSE, 
                    cluster_cols=FALSE,
                    show_colnames = FALSE,
                    annotation_colors = myColors2)

pheatplot

#####

### GO classification ###
#####
library(GO.db)

# GO TERM table all genes
vas<-toTable(GOTERM)
colnames(vas)[2]<-"repgo_id"
vasbp<-vas[vas$Ontology=="BP",]
#####

#### RBP genes profile
#####

# Select from corrected batch cRPKM expression values, RBPs
matDEGs_RBPs<-y2[row.names(y2) %in% unique(DEGs_cohorts_RBPs$GENE.x),]

# 1st do PCA
# PCA with corrected batch effect
tmatDEGs_RBPs<-t(matDEGs_RBPs)
rm(matDEGs_RBPs)

cond<-(apply(tmatDEGs_RBPs, 2, var)!=0)
tmatDEGs_RBPs<-tmatDEGs_RBPs[, cond, drop = FALSE]

#PCA
res.pca<-prcomp(tmatDEGs_RBPs,scale=TRUE)
summary(res.pca)

(res.pca$sdev^2/sum(res.pca$sdev^2))[1:5]

#Proportion of variance

pcatab_batch<-as.data.frame(res.pca$x)
pcatab_batch$names<-row.names(pcatab_batch)
row.names(pcatab_batch)

pcatab_batch<-pcatab_batch[,c(1:3)]

metadata$SampID

setdiff(metadata$SampID,row.names(pcatab_batch))

pcatab_batch<-merge(metadata,pcatab_batch,by.x="SampID",by.y="row.names")
pcatab_batch$var<-paste(pcatab_batch$Source,pcatab_batch$Pheno,sep="_")

library(ggnewscale)

pcatab_batch$GSEA<-as.factor(pcatab_batch$GSEA)

# Detecting outliers, estimate z-score
outbatchTLE<-data.frame(SampID=pcatab_batch$SampID,
                        zscore_PC1=abs(scale(pcatab_batch$PC1, center = TRUE, scale = TRUE)),
                        zscore_PC2=abs(scale(pcatab_batch$PC2, center = TRUE, scale = TRUE)))

pcatab_batch<-merge(pcatab_batch,outbatchTLE,by="SampID")


ggplot(pcatab_batch, aes(PC1, PC2,label=SampID)) +
  #geom_point(data=pcatab[pcatab$source=="BM",],size=10,color="black") +
  #geom_point(data=pcatab[pcatab$source=="Thymus",],size=10,color="grey") +
  #geom_point(data=pcatab[pcatab$source=="TALL",],size=10,color="red") +
  geom_point(size=8,color="black")+
  geom_point(size=7,aes(color=Type)) +
  #geom_point(data=pcatab[pcatab$GSEA=="GSE57982",],size=6,color="black") +
  scale_color_brewer(palette = "Paired")+
  #new_scale_color()+
  #stat_ellipse(aes(color=GSEA,linetype=GSEA),level=0.9)+
  #scale_color_brewer(palette = "Paired")+
  #scale_linetype_manual(values=c("solid","dotted"))+
  #geom_label_repel(data=pcatab_batch[pcatab_batch$GSEA=="GSE69239",],aes(label=Source), size=4,fontface = "italic",alpha=0.7)+
  #geom_label_repel(data=pcatab_batch[pcatab_batch$Source=="Thymus" ,],aes(label=Source), size=4,fontface = "italic",alpha=0.7,color="grey")+
  geom_label_repel(data=pcatab_batch[pcatab_batch$zscore_PC1>=3 | pcatab_batch$zscore_PC2>=3 ,],aes(label=SampID), size=3,fontface = "italic",alpha=0.7,color="black")+
  ylab("PC2 (7.58% variance)")+
  xlab("PC1 (28.36% variance)")+
  theme_minimal()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()


#####

#### ALTERNATIVE SPLICING EVENTS ACROSS T-CELL DEVELOPMENT AND THYMUS comparison
#####
## Barplot with up-and down-regulated AS events
#tev<-read_xlsx("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/table_ASevents_TCell.xlsx")
tev<-read_xlsx("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/table_ASevents_TCell.xlsx")
tev<-melt(tev,id.vars=c("Type","Class"))
tev$variable<-factor(tev$variable,levels = c("CD8","CD4","DP","DN","Thy2","Thy1","BCP","CLP","LMPP","HSC"))

ggplot(tev,aes(x = variable, y = value, group = Type, fill = Class)) +
  geom_bar(stat = "identity", width = 0.75,color="black") +
  coord_flip()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill =  "grey90")) +
  scale_fill_brewer(palette="Paired")+
  xlab("Stage")+
  ylab("AS Events")+
  theme(legend.position = "bottom",
        axis.text.y = element_text(size=9),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  theme_bw()

# T-Cell samples
Tsamp<- c("GENE","EVENT","HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8")
Tord<-c("HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8")

# vastools diff alt splicing each cell type vs all
#HSCdev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/HSC_rest_minPSI10.txt',header = TRUE)
HSCdev<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/Tdev_expr_out/isoforms_diff/HSC_rest_minPSI10.txt',header = TRUE)
HSCdev<-HSCdev[,c(1,2,seq(7, ncol(HSCdev), by=2))]
HSCdev<-HSCdev[,colnames(HSCdev) %in% Tsamp]
HSCdev$type<-"HSC"

HSCdev_med<-HSCdev[,c(1:2,9,10,6,3,11,12,7,8,4,5)]
HSCdev_med$med<-apply(HSCdev_med[-c(1:3)], 1, FUN=median, na.rm=TRUE)
HSCdev_med$dif<-HSCdev_med$HSC-HSCdev_med$med
HSCdev_med$state<-ifelse(HSCdev_med$dif<0,c("HSC_down"),c("HSC_UP"))
table(HSCdev_med[grepl("EX",HSCdev_med$EVENT),]$state)
table(HSCdev_med[grepl("INT",HSCdev_med$EVENT),]$state)
table(HSCdev_med[grepl("ALT",HSCdev_med$EVENT),]$state)

  
LMPPdev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/LMPP_rest_minPSI10.txt',header = TRUE)
LMPPdev<-LMPPdev[,c(1,2,seq(7, ncol(LMPPdev), by=2))]
LMPPdev<-LMPPdev[,colnames(LMPPdev) %in% Tsamp]
LMPPdev$type<-"LMPP"

CLPdev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/CLP_rest_minPSI10.txt',header = TRUE)
CLPdev<-CLPdev[,c(1,2,seq(7, ncol(CLPdev), by=2))]
CLPdev<-CLPdev[,colnames(CLPdev) %in% Tsamp]
CLPdev$type<-"CLP"

BCPdev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/BCP_rest_minPSI10.txt',header = TRUE)
BCPdev<-BCPdev[,c(1,2,seq(7, ncol(BCPdev), by=2))]
BCPdev<-BCPdev[,colnames(BCPdev) %in% Tsamp]
BCPdev$type<-"BCP"

Thy1dev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/Thy1_rest_minPSI10.txt',header = TRUE)
Thy1dev<-Thy1dev[,c(1,2,seq(7, ncol(Thy1dev), by=2))]
Thy1dev<-Thy1dev[,colnames(Thy1dev) %in% Tsamp]
Thy1dev$type<-"Thy1"

Thy2dev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/Thy2_rest_minPSI10.txt',header = TRUE)
Thy2dev<-Thy2dev[,c(1,2,seq(7, ncol(Thy2dev), by=2))]
Thy2dev<-Thy2dev[,colnames(Thy2dev) %in% Tsamp]
Thy2dev$type<-"Thy2"

DNdev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/DN_rest_minPSI10.txt',header = TRUE)
DNdev<-DNdev[,c(1,2,seq(7, ncol(DNdev), by=2))]
DNdev<-DNdev[,colnames(DNdev) %in% Tsamp]
DNdev$type<-"DN"

DPdev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/DP_rest_minPSI10.txt',header = TRUE)
DPdev<-DPdev[,c(1,2,seq(7, ncol(DPdev), by=2))]
DPdev<-DPdev[,colnames(DPdev) %in% Tsamp]
DPdev$type<-"DP"

CD4dev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/CD4_rest_minPSI10.txt',header = TRUE)
CD4dev<-CD4dev[,c(1,2,seq(7, ncol(CD4dev), by=2))]
CD4dev<-CD4dev[,colnames(CD4dev) %in% Tsamp]
CD4dev$type<-"CD4"

CD8dev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/CD8_rest_minPSI10.txt',header = TRUE)
CD8dev<-CD8dev[,c(1,2,seq(7, ncol(CD8dev), by=2))]
CD8dev<-CD8dev[,colnames(CD8dev) %in% Tsamp]
CD8dev$type<-"CD8"


#CD4CD8dev<-read.delim('/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/CD4CD8_rest_minPSI10.txt',header = TRUE)
CD4CD8dev<-CD4CD8dev[,c(1,2,seq(7, ncol(CD4CD8dev), by=2))]
CD4CD8dev<-CD4CD8dev[,colnames(CD4CD8dev) %in% Tsamp]
CD4CD8dev$type<-"CD4CD8"


evtype<-merge(HSCdev[,c(2,ncol(HSCdev))],LMPPdev[,c(2,ncol(LMPPdev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,CLPdev[,c(2,ncol(CLPdev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,BCPdev[,c(2,ncol(BCPdev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,Thy1dev[,c(2,ncol(Thy1dev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,Thy2dev[,c(2,ncol(Thy2dev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,DNdev[,c(2,ncol(DNdev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,DPdev[,c(2,ncol(DPdev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,CD4dev[,c(2,ncol(CD4dev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,CD8dev[,c(2,ncol(CD8dev))],by="EVENT",all=TRUE)
evtype<-merge(evtype,CD4CD8dev[,c(2,ncol(CD4CD8dev))],by="EVENT",all=TRUE)

colnames(evtype)<-c("EVENT","HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8","CD4CD8")
evtype$class<-paste(evtype$HSC,evtype$LMPP,evtype$CLP,evtype$BCP,evtype$Thy1,evtype$Thy2,evtype$DN,evtype$DP,evtype$CD4,evtype$CD8,evtype$CD4CD8,sep="_")
evtype$class<-gsub('_NA','',evtype$class)
evtype$class<-gsub('NA_','',evtype$class)

evtype<-evtype[,c(1,ncol(evtype))]
row.names(evtype)<-evtype$EVENT


# Events exclusive for one cell type
exc_cell<-evtype[!grepl('_',evtype$class),]
exc_cell$EVENT<-NULL
evtype$EVENT<-NULL

# Merge PSI for all datasets
mdev<-rbind(HSCdev,LMPPdev,CLPdev,BCPdev,Thy1dev,Thy2dev,DNdev,DPdev,CD4dev,CD8dev,CD4CD8dev)

mdev$type<-NULL

#select unique events
mdev<-mdev[!duplicated(mdev), ]
row.names(mdev)<-mdev$EVENT
mdev$EVENT<-NULL
mdev$GENE<-NULL

# Select exclusive events for one cell type
mdev<-mdev[row.names(mdev) %in% row.names(exc_cell),]

mdev[, order(Tord)]
mdev<-mdev[,Tord]

s1<-length(unique(exc_cell$class))

dColors = list(
  class = c(HSC=brewer.pal(s1,"BrBG")[1],
           LMPP=brewer.pal(s1,"BrBG")[2],
           CLP=brewer.pal(s1,"BrBG")[3],
           BCP=brewer.pal(s1,"BrBG")[4],
           Thy1=brewer.pal(s1,"BrBG")[5],
           Thy2=brewer.pal(s1,"BrBG")[6],
           DN=brewer.pal(s1,"BrBG")[7],
           DP=brewer.pal(s1,"BrBG")[8],
           CD4CD8=brewer.pal(s1,"BrBG")[9],
           CD4=brewer.pal(s1,"Paired")[1],
           CD8=brewer.pal(s1,"Paired")[2]))

pheatmap(mdev,
#         annotation_col = annotcel,
         annotation_row = exc_cell,
         show_rownames = FALSE,
          cluster_cols=FALSE,
          cluster_rows = FALSE,
         fontsize = 10,
#         cutree_cols = 9,
#         cutree_rows = 7,
         fontsize_col = 10,
      annotation_colors = dColors)

# Gene identifiers
mdevid<-merge(selev_s[,c(1:2)],mdev,by.y=0,by.x="EVENT")
mdevid<-merge(evtype,mdevid,by.x=0,by.y="EVENT")
row.names(mdevid)<-mdevid$Row.names
mdevid$Row.names<-NULL

mdevid$EVENT<-row.names(mdevid)
mdevid_melt<-melt(mdevid,id.vars=c("EVENT","class","GENE"))


ggplot(data=difASm[difASm$EVENT=="HsaINT0046168" & difASm$GSEA=="GSE69239",],aes(x=Source,y=value))+
  geom_point(color="black",size=6)+
  geom_point(data=difASm[difASm$EVENT=="HsaINT0046168" & difASm$GSEA=="GSE69239",],aes(color=Type),size=5)+
  geom_ribbon(data=difASm[difASm$EVENT=="HsaINT0046168" & difASm$GSEA=="GSE69239" & difASm$Type=="BM",],aes(group=GSEA,ymin=0,ymax=value),color="darkolivegreen",size=1,alpha=0)+
  geom_ribbon(data=difASm[difASm$EVENT=="HsaINT0046168" & difASm$GSEA=="GSE69239" & (difASm$Type=="Thymus" | difASm$Source=="CLP"),],aes(group=GSEA,ymin=0,ymax=value),color="orange",size=1,alpha=0)+
  scale_color_brewer(palette="Set2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))


## Functional enrichment of genes with AS exclusive for each cell

## Functional enrichment

hsc_go_altex<-read.delim("/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/AltEx-HSC_rest_minPSI10.txt",header = FALSE)
hsc_go_IRUP<-read.delim("/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/IR_UP-HSC_rest_minPSI10.txt",header = FALSE)
hsc_go_IRDOWN<-read.delim("/Volumes/cancer/TCell/vast_output/Tdev_expr_out/isoforms_diff/IR_DOWN-HSC_rest_minPSI10.txt",header = FALSE)


library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "InterPro_Domains_2019",
         "GO_Cellular_Component_2017",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2016")

enrCELL <- enrichr(as.character(hsc_go_IRDOWN$V1), dbs)

enrtype <- enrCELL[["GO_Biological_Process_2018"]]

bpsub<-subset(enrtype,enrtype$Adjusted.P.value<0.01)
bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=unique((bpsub$Term))[rev(order(bpsub$P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))

# for one group
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$Adjusted.P.value)), ]$Term))


ggplot(bpsub,aes(y=-log10(Adjusted.P.value),x=Term),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity")+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=8, hjust = 1),
        axis.text.y = element_text(size=10, hjust = 1))

#####


### read differential AS events ###
#####
# EGA_TLE
EGA_AS<-read.delim('/Volumes/cancer/TCell/vast_output/EGA_TLE_expr_out/isoforms_diff/EGA_TLE_diffisoforms.txt',header = TRUE)
EGA_AS<-EGA_AS[,c(1,2,seq(7, ncol(EGA_AS), by=2))]
EGA_AS$type<-"EGA_TLE"

# GSE109231
GSE109231_AS<-read.delim('/Volumes/cancer/TCell/vast_output/GSE109231_expr_out/isoforms_diff/GSE109231_minPSI10.txt',header = TRUE)
GSE109231_AS<-GSE109231_AS[,c(1,2,seq(7, ncol(GSE109231_AS), by=2))]
GSE109231_AS$type<-"GSE109231"

# GSE57982
GSE57982_AS<-read.delim('/Volumes/cancer/TCell/vast_output/GSE57982_expr_out/isoforms_diff/GSE57982_minPSI10.txt',header = TRUE)
GSE57982_AS<-GSE57982_AS[,c(1,2,seq(7, ncol(GSE57982_AS), by=2))]
GSE57982_AS$type<-"GSE57982"

## T-ALL samples vs DP
TALLDP_AS<-read.delim('/Volumes/cancer/TCell/vast_output/TALL_DP_expr_out/isoforms_diff/TALL_DP_minPSI10.txt',header = TRUE)
TALLDP_AS<-TALLDP_AS[,c(1,2,seq(7, ncol(TALLDP_AS), by=2))]
TALLDP_AS$type<-"TALL_DP"


astype<-merge(EGA_AS[,c(2,ncol(EGA_AS))],GSE109231_AS[,c(2,ncol(GSE109231_AS))],by="EVENT",all=TRUE)
astype<-merge(astype,GSE57982_AS[,c(2,ncol(GSE57982_AS))],by="EVENT",all=TRUE)
astype<-merge(astype,TALLDP_AS[,c(2,ncol(TALLDP_AS))],by="EVENT",all=TRUE)


colnames(astype)<-c("EVENT","EGA_TLE","GSE109231","GSE57982","TALL_DP")
astype$class<-paste(astype$EGA_TLE,astype$GSE109231,astype$GSE57982,astype$TALL_DP,sep="_")
astype$class<-gsub('_NA','',astype$class)
astype$class<-gsub('NA_','',astype$class)

astype<-astype[,c(1,ncol(astype))]
row.names(astype)<-astype$EVENT


# Events shared by two cohorts
exc_leu<-astype[grepl('EGA_TLE_',astype$class) | grepl('_GSE',astype$class),]
exc_leu$EVENT<-NULL

astype$EVENT<-NULL

# Merge PSI for all datasets
mleu<-rbind(EGA_AS,GSE109231_AS,GSE57982_AS,TALLDP_AS)
mleu$type<-NULL

#select unique events
mleu<-mleu[!duplicated(mleu), ]
row.names(mleu)<-mleu$EVENT
mleu$EVENT<-NULL
mleu$GENE<-NULL

# Select common events for two cohorts
mleu<-mleu[row.names(mleu) %in% row.names(exc_leu),]


s1<-length(unique(astype$class))

colnames(mleu)<-gsub('_1','',colnames(mleu))
colnames(mleu)<-gsub('EGAR.*_','',colnames(mleu))  
  
samgsea<-subset(metadata,select=c(SampID,GSEA,Type,Event_free_surv_days,Pheno))
row.names(samgsea)<-samgsea$SampID
samgsea$SampID<-NULL
#samgsea$Prog<-gsub('^D.*','D',samgsea$Prog)

samgsea<-samgsea[row.names(samgsea)!="R3",]


lColors = list(
  class = c(EGA_TLE=brewer.pal(s1,"Set3")[1],
            EGA_TLE_GSE109231=brewer.pal(s1,"Set3")[2],
            EGA_TLE_GSE109231_TALL_DP=brewer.pal(s1,"Set3")[3],
            EGA_TLE_TALL_DP=brewer.pal(s1,"Set3")[4],
            GSE109231=brewer.pal(s1,"Set3")[5],
            GSE109231_GSE57982=brewer.pal(s1,"Set3")[6],
            GSE109231_TALL_DP=brewer.pal(s1,"Set3")[7],
            GSE57982=brewer.pal(s1,"Set3")[8],
            GSE57982_TALL_DP=brewer.pal(s1,"Set3")[9],
            TALL_DP=brewer.pal(s1,"Set3")[1]))

pheatmap(mleu[,colnames(mleu)!="R3"],
         annotation_col = samgsea,
         annotation_row = astype,
         show_rownames = FALSE,
#         cluster_cols=FALSE,
#         cluster_rows = FALSE,
#scale="row",
         fontsize = 6,
                  cutree_cols = 5,
                  cutree_rows = 4,
         fontsize_col = 7,
         annotation_colors = lColors)


# Gene identifiers
idmleu<-merge(selev_s[,c(1:2)],mleu,by.y=0,by.x="EVENT")
idmleu<-merge(astype,idmleu,by.x=0,by.y="EVENT")
row.names(idmleu)<-idmleu$Row.names
idmleu$Row.names<-NULL

idmleu$EVENT<-row.names(idmleu)
idmleu_melt<-melt(idmleu,id.vars=c("EVENT","class","GENE"))
idmleu_melt<-merge(idmleu_melt,candriver,by.x="GENE",by.y="SYMBOL",all.x="TRUE")


# Genes common in at least two comparisons
idmleu_com<-idmleu_melt[idmleu_melt$EVENT %in% row.names(exc_leu),]
idmleu_com<-merge(metadata,idmleu_com,by.x="SampID",by.y="variable")



idmleu_com$Type<-factor(idmleu_com$Type,levels = c("Thymus","TALL","TLBL","BM"))


ggplot(idmleu_com[idmleu_com$class=="EGA_TLE_TALL_DP",],aes(x=Type,y=value))+
  geom_point(color="black",size=5,aes(shape=class))+
  geom_point(aes(color=Type,shape=class),size=4)+
  geom_boxplot(alpha=0)+
  scale_color_brewer(palette="Paired")+
  facet_wrap(EVENT+GENE~GSEA,scales="free")+
  theme_light()


## Multivenn diagram
ven1<-data.frame(table(astype$class))
ven1$Freq<-as.numeric(ven1$Freq)
ven1<-ven1[order(ven1$Freq)]

ven1$Var1<-as.factor(ven1$Var1)
ven1$Var1 <- factor(ven1$Var1, levels=unique((ven1$Var1))[order(ven1$Freq)])

ggplot(ven1, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color="black",fill="orange",alpha=0.5)+
  geom_text(aes(label=Freq), vjust=-0.1, color="black", size=10)+
  theme_minimal()+
  theme(axis.text.x=element_text(size=14,angle=90,face="bold.italic"),
        axis.text.y=element_text(size=14,face="bold.italic"),
        axis.title=element_blank())

dim(EGA_AS)
dim(GSE109231_AS)
dim(GSE57982_AS)
dim(TALLDP_AS)

ven2<-data.frame(Comp=c("EGA_TLE","GSE109231","GSE57982","TALL_DP"),
                Frequency=c(dim(EGA_AS)[1],dim(GSE109231_AS)[1],dim(GSE57982_AS)[1],dim(TALLDP_AS)[1]))

ven2$Comp<-as.factor(ven2$Comp)
ven2$Comp <- factor(ven2$Comp, levels=unique((ven2$Comp))[rev(order(ven2$Frequency))])

ggplot(ven2, aes(x=Comp, y=Frequency)) +
  geom_bar(stat="identity", color="black",fill="blue",alpha=0.5)+
  geom_text(aes(label=Frequency), vjust=-0.1, color="black", size=10)+
  theme_minimal()+
  theme(axis.text.x=element_text(size=14,angle=90,face="bold.italic"),
        axis.text.y=element_text(size=14,face="bold.italic"),
        axis.title=element_blank())


## Functional enrichment

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "InterPro_Domains_2019",
         "GO_Cellular_Component_2017",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2016")

enrAS <- enrichr(as.character(unique(idmleu[idmleu$class=="EGA_TLE_GSE109231",]$GENE)), dbs)

enrleu <- enrAS[["GO_Biological_Process_2018"]]

bpsub<-subset(enrleu,enrleu$P.value<0.01)
bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=unique((bpsub$Term))[rev(order(bpsub$P.value))])

library(tidyr)

#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))

# for one group
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$Adjusted.P.value)), ]$Term))


ggplot(bpsub[1:5,],aes(y=-log10(Adjusted.P.value),x=Term),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",color="black",fill="darkolivegreen")+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=10, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1))





ensembl76 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

anoEGA <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                  values=unique(EGA_AS$GENE),
                  filter="external_gene_name",
                  mart=hg) 




EGACC<-merge(EGA_AS,anoEGA,by.x="GENE",by.y="external_gene_name",all.x=TRUE)
#empty values as NA
EGACC$go_id<-gsub("^$", NA, EGACC$go_id)

EGACC<-merge(vascc,EGACC,by="go_id",all.x=TRUE)



EGA_GO_MEMBRA<-EGACC[grepl('plasma membrane',EGACC$Term),]
unique(EGA_GO_MEMBRA$GENE)

EGA_GO_MEMBRA<-EGACC[grepl('apical plasma membrane',EGACC$Term),]

######

## MOVED TO BCAT PROJECT SESSION:

## B-CAT CORRELATION NETWORWK MODELS
#####
###### CORRELATION MATRIX 

#Bcat gene = ENSG00000168036

## Correlation between one target genes and others!

## Dataframes categorized by T-ALL, or normal T-cell development (DSmat is not corrected)
# DO NOT INCLUDE TLE samples
Tdev<-subset(DSmat_noTLE,DSmat_noTLE$GSEA=="GSE69239")
Tall<-subset(DSmat_noTLE,DSmat_noTLE$Source=="TALL")
Tlbl<-subset(DSmat_noTLE,DSmat_noTLE$Source=="TLBL")
Tleuk<-subset(DSmat_noTLE,DSmat_noTLE$Source=="TALL" | DSmat_noTLE$Source=="TLBL")
Ttarget<-subset(DSmat_noTLE,DSmat_noTLE$GSEA=="TARGET")

# Cell lines to predict later
Tline<-subset(DSmat_noTLE,DSmat_noTLE$GSEA=="line")

# In TALL and T-LBL databases, bcat coexpressed genes or kaiso or brca1
corout<-apply(Tleuk[,18:ncol(Tleuk)], 2, cor.test,Tleuk$ENSG00000168036,method="spearman")

#For only target
corout<-apply(Ttarget[,18:ncol(Ttarget)], 2, cor.test,Ttarget$ENSG00000168036,method="spearman")


corout1<-as.data.frame(sapply(corout, "[[", "p.value"))
corout2<-as.data.frame(sapply(corout, "[[", "estimate"))
corout<-cbind(corout1,corout2)
names(corout)<-c("Pval","Rho")
corout$Gene<-row.names(corout)
corout<-merge(geneid,corout,by="Gene",all.y=TRUE)


corout_sig<-subset(corout,corout$Rho>=0.4)
corout_sig<-corout_sig[with(corout_sig,order(corout_sig$Pval, corout_sig$Rho)),]

negcorctnnb1<-(corout_sig[corout_sig$Rho<0,])$NAME
length(negcorctnnb1)
poscorctnnb1<-(corout_sig[corout_sig$Rho>0,])$NAME
length(poscorctnnb1)

# vector for functional enrichment
(corout_sig[corout_sig$Rho>0.4,])$NAME

# KAISO for all
corout_kai<-apply(Tleuk[,18:ncol(Tleuk)], 2, cor.test,Tleuk$ENSG00000177485,method="spearman")

# KAISO for TARGET
corout_kai<-apply(Ttarget[,18:ncol(Ttarget)], 2, cor.test,Ttarget$ENSG00000177485,method="spearman")

corout1_kai<-as.data.frame(sapply(corout_kai, "[[", "p.value"))
corout2_kai<-as.data.frame(sapply(corout_kai, "[[", "estimate"))
corout_kai<-cbind(corout1_kai,corout2_kai)
names(corout_kai)<-c("Pval","Rho")
corout_kai$Gene<-row.names(corout_kai)
corout_kai<-merge(geneid,corout_kai,by="Gene",all.y=TRUE)


# BRCA1 for all
corout_brca<-apply(Tleuk[,18:ncol(Tleuk)], 2, cor.test,Tleuk$ENSG00000012048,method="spearman")

# BRCA1 for TARGET
corout_brca<-apply(Ttarget[,18:ncol(Ttarget)], 2, cor.test,Ttarget$ENSG00000012048,method="spearman")

corout1_brca<-as.data.frame(sapply(corout_brca, "[[", "p.value"))
corout2_brca<-as.data.frame(sapply(corout_brca, "[[", "estimate"))
corout_brca<-cbind(corout1_brca,corout2_brca)
names(corout_brca)<-c("Pval","Rho")
corout_brca$Gene<-row.names(corout_brca)
corout_brca<-merge(geneid,corout_brca,by="Gene",all.y=TRUE)


# Model correlation example RBM39
as.character((geneid[geneid$NAME=="RBM39",])$Gene)
as.character((geneid[geneid$NAME=="CTNNB1",])$Gene)

#model in Tall + tlbl
model_RBM39<-lm(ENSG00000131051~ENSG00000168036,data=Tleuk)

#observed CTNNB1 (constant)
obs_CTNNB1<-data.frame(
  ENSG00000168036=c(Tline[,colnames(Tline) == "ENSG00000168036"]))

row.names(obs_CTNNB1)<-row.names(Tline[,colnames(Tline) == "ENSG00000168036" | colnames(Tline) == "ENSG00000131051" ,])

# Observed x gene (specific for each gene)
obs_RBM39<-data.frame(
  ENSG00000131051=c(Tline[,colnames(Tline) == "ENSG00000131051" ,]))

row.names(obs_RBM39)<-row.names(Tline[,colnames(Tline) == "ENSG00000168036" | colnames(Tline) == "ENSG00000131051" ,])

# prediction
pre_RBM39<-predict(lm(ENSG00000131051~ENSG00000168036,data=Tleuk),
                   newdata = obs_CTNNB1,
                   interval="prediction")
pre_RBM39


# Save data models for each correlated gene and prediction
#Select matrix with expression of correlated genes with bcat
modcor<-Tleuk[,colnames(Tleuk) %in% c(as.character(corout_sig$Gene),"ENSG00000168036")]
dim(modcor)

# Output matrix with expected expression based on model
exp_fit<-data.frame(t(apply(modcor, 2, function(x) predict(lm(x~ENSG00000168036,data=modcor),
                                                           newdata = obs_CTNNB1,
                                                           interval="prediction")[,1])))

colnames(exp_fit)<-gsub('^',"expected_",colnames(exp_fit))

exp_lwr<-data.frame(t(apply(modcor, 2, function(x) predict(lm(x~ENSG00000168036,data=modcor),
                                                           newdata = obs_CTNNB1,
                                                           interval="prediction")[,2])))
colnames(exp_lwr)<-gsub('^',"lwr_",colnames(exp_lwr))

exp_upr<-data.frame(t(apply(modcor, 2, function(x) predict(lm(x~ENSG00000168036,data=modcor),
                                                           newdata = obs_CTNNB1,
                                                           interval="prediction")[,3])))
colnames(exp_upr)<-gsub('^',"upr_",colnames(exp_upr))

# Observed x gene (specific for each gene)
obs_gene<-data.frame(t(Tline[,colnames(Tline) %in% c(as.character(corout_sig$Gene),"ENSG00000168036")]))
colnames(obs_gene)<-gsub('^','observed_',colnames(obs_gene))

# bind observed and expected
est_mod<-cbind(exp_fit,exp_lwr,exp_upr,obs_gene)

est_mod$as_cont_jurk<-ifelse(est_mod$lwr_control_Jurkat < est_mod$observed_control_Jurkat &  est_mod$observed_control_Jurkat < est_mod$upr_control_Jurkat, 
                             c("YES"),c("NO")) 

est_mod$as_cont_RPMI<-ifelse(est_mod$lwr_control_RPMI < est_mod$observed_control_RPMI &  est_mod$observed_control_RPMI < est_mod$upr_control_RPMI, 
                             c("YES"),c("NO")) 

est_mod$as_sh_jurk<-ifelse(est_mod$lwr_shbcat_Jurkat < est_mod$observed_shbcat_Jurkat &  est_mod$observed_shbcat_Jurkat < est_mod$upr_shbcat_Jurkat, 
                           c("YES"),c("NO")) 

est_mod$as_sh_RPMI<-ifelse(est_mod$lwr_shbcat_RPMI < est_mod$observed_shbcat_RPMI &  est_mod$observed_shbcat_RPMI < est_mod$upr_shbcat_RPMI, 
                           c("YES"),c("NO")) 


est_mod<-merge(geneid,est_mod,by.x="Gene",by.y="row.names")

coreg_bcat<-est_mod[est_mod$as_sh_jurk =="YES" & est_mod$as_sh_RPMI=="YES" & est_mod$as_cont_jurk=="YES" & est_mod$as_cont_RPMI=="YES",]
nocont<-est_mod[est_mod$as_cont_jurk =="NO" & est_mod$as_cont_RPMI=="NO" & est_mod$as_sh_jurk =="NO" & est_mod$as_sh_RPMI=="NO",]

coreg_bcat$NAME

gene<-as.character((geneid[geneid$NAME=="ZNF76",])$Gene)
bcat<-as.character((geneid[geneid$NAME=="CTNNB1",])$Gene)

ggplot(rbind(Tleuk,Tline,Tdev),aes_string(x=bcat,y=gene))+
  geom_point(size=5)+
  geom_point(data=Tline,aes(color=Pheno,shape=Sample),size=7)+
  scale_color_brewer("Spectral")+
  new_scale_color()+
  geom_point(aes(color=Type),size=3)+
  scale_color_brewer(palette="Paired")+
  geom_smooth(data=rbind(Tleuk,Tdev),method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,aes(color=GSEA))+
  #geom_smooth(data=DSmat[DSmat$Source!="TALL" ,],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="purple")+
  #geom_smooth(data=gene_sam[gene_sam$RBPJ<75 & gene_sam$Source!="TALL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="red")+
  #geom_label_repel(data=DSmat[DSmat$GSEA=="line" ,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  facet_wrap(~GSEA,ncol=2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))
#coord_fixed()


# Plots
gene<-as.character((geneid[geneid$NAME=="BRCA1",])$Gene)
bcat<-as.character((geneid[geneid$NAME=="CTNNB1",])$Gene)
kaiso<-as.character((geneid[geneid$NAME=="ZBTB3",])$Gene)

DSmat_noTLE$Source<-factor(DSmat_noTLE$Source,levels = c("HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8","Thymus","TALL","TLBL","Jurkat","RPMI"))

ggplot(DSmat_noTLE[DSmat_noTLE$GSEA=="TARGET",],aes_string(x=bcat,y=kaiso))+
  geom_point(size=5)+
  geom_point(aes(color=First_event),size=4)+
  scale_color_brewer(palette="Paired")+
  new_scale_color()+
  #  geom_point(data=DSmat_noTLE[DSmat_noTLE$GSEA=="GSE69239",],size=3.5)+
  #  geom_point(data=DSmat_noTLE[DSmat_noTLE$GSEA=="GSE69239",],size=2)+
  #  geom_smooth(data=DSmat_noTLE[DSmat_noTLE$Source=="TALL" | DSmat_noTLE$Source=="TLBL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,aes(color=GSEA))+
  #geom_smooth(data=DSmat[DSmat$Source!="TALL" ,],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="purple")+
  #geom_smooth(data=gene_sam[gene_sam$RBPJ<75 & gene_sam$Source!="TALL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="red")+
  #  geom_label_repel(data=DSmat_noTLE[DSmat_noTLE$GSEA=="line" ,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  geom_smooth(method="lm")+
  facet_wrap(~GSEA,ncol=2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))
#coord_fixed()

cor.test(DSmat_noTLE[DSmat_noTLE$GSEA=="TARGET",]$ENSG00000185670,DSmat_noTLE[DSmat_noTLE$GSEA=="TARGET",]$ENSG00000168036)

#Plot only normal T-cell cohort
ggplot(DSmat_noTLE[DSmat_noTLE$Source!="TALL" & DSmat_noTLE$Source!="TLBL" & DSmat_noTLE$GSEA!="line",],aes_string(x=gene,y=kaiso))+
  geom_point(size=7)+
  geom_point(aes(color=Type),size=6)+
  scale_color_brewer(palette="Paired")+
  geom_smooth(data=DSmat_noTLE[DSmat_noTLE$Source!="TALL" & DSmat_noTLE$Source!="TLBL" & DSmat_noTLE$GSEA!="line",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="blue")+
  #geom_smooth(data=gene_sam[gene_sam$RBPJ<75 & gene_sam$Source!="TALL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="red")+
  #geom_label_repel(data=gene_sam[gene_sam$DDX3X>350,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

gene
bcat

#Plot only T-ALL cohort
ggplot(DSmat_noTLE[DSmat_noTLE$Source=="TALL" | DSmat_noTLE$Source=="TLBL" | DSmat_noTLE$GSEA=="line",],aes(x=log(ENSG00000012048),y=log(ENSG00000185670)))+
  geom_point(size=7)+
  geom_point(aes(color=Type),size=6)+
  scale_color_brewer(palette="Paired")+
  geom_smooth(data=DSmat_noTLE[DSmat_noTLE$Source=="TALL" | DSmat_noTLE$Source=="TLBL" | DSmat_noTLE$GSEA=="line",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="purple")+
  #geom_smooth(data=gene_sam[gene_sam$RBPJ<75 & gene_sam$Source!="TALL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="red")+
  #geom_label_repel(data=gene_sam[gene_sam$DDX3X>350,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

#####

# list of DEGs shbcat in RPMI from chip.R script cross with b-cat chips and correlation models (INCLUDES VIOLIN PLOT FOR RNASEQ AND CHIPS)
#####

res_df_rna_rpmi<-read.delim('/Users/yguillen/Desktop/temp/Gekas_RNASeq/DESeq2/RPMI/res_df_rna_rpmi.txt')
res_df_rna_rpmi

res_df_rna_jurkat<-read.delim('/Users/yguillen/Desktop/temp/Gekas_RNASeq/DESeq2/Jurkat/res_df_rna_jurkat.txt')
res_df_rna_jurkat

# Kaiso chips
kaiso_peaks<-read_xlsx("/Users/yguillen/Desktop/temp/beta_catenin_project/KAISO_ChIPs/kaiso_in4.xlsx")
kaiso_peaks_tab<-kaiso_peaks[,colnames(kaiso_peaks) %in% c("annotation","ENSEMBL")]
kaiso_peaks_tab<-kaiso_peaks_tab[!duplicated(kaiso_peaks_tab),]
kaiso_peaks_tab$annotation<-gsub('Intron\\s(.*)','Intron',kaiso_peaks_tab$annotation)
kaiso_peaks_tab$annotation<-gsub('Exon\\s(.*)','Exon',kaiso_peaks_tab$annotation)
kaiso_target<-data.frame(table(kaiso_peaks_tab$ENSEMBL))
row.names(kaiso_target)<-kaiso_target$Var1
kaiso_target$Var1<-NULL
colnames(kaiso_target)<-("kaiso")



# b-cat chips
# genes with peaks in control 
# RD antibody
bcat_cpeaks<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_DBmerge_nearest.genes.txt",header = TRUE)
dim(bcat_cpeaks)
length(unique(bcat_cpeaks$gene))
length(unique(bcat_cpeaks$Name))

# SC antibody
#bcat_cpeaks<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_RBCSmerge_nearest.genes.txt",header = TRUE)
#dim(bcat_cpeaks)
#length(unique(bcat_cpeaks$Name))

# if a peak has two annotated genes, they are all taken into account
bcat_cpeaks_ens<-getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                       values= unique(bcat_cpeaks$gene),
                       filter="external_gene_name",
                       mart=hg) 

length(unique(bcat_cpeaks_ens$external_gene_name))

bcat_cpeaks<-merge(bcat_cpeaks_ens,bcat_cpeaks,by.x="external_gene_name",by.y="gene",all.y=TRUE)
dim(bcat_cpeaks)
length(unique(bcat_cpeaks$Name))
length(unique(bcat_cpeaks$external_gene_name))

# genes with peaks in lithium
#RD Antibody
bcat_lpeaks<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_DBLmerge_nearest.genes.txt",header = TRUE)
dim(bcat_lpeaks)
length(unique(bcat_lpeaks$Name))

#SC Antibody
#bcat_lpeaks<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_RBLSmerge_nearest.genes.txt",header = TRUE)


bcat_lpeaks_ens<-getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                       values= unique(bcat_lpeaks$gene),
                       filter="external_gene_name",
                       mart=hg) 

bcat_lpeaks<-merge(bcat_lpeaks_ens,bcat_lpeaks,by.x="external_gene_name",by.y="gene",all.y=TRUE)
dim(bcat_lpeaks)
length(unique(bcat_lpeaks$external_gene_name))

# genes with peaks in TCF 
tcf_peaks_dt<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_TCFmerge_nearest.genes.txt",header = TRUE)

tcf_peaks<-getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                 values= unique(tcf_peaks_dt$gene),
                 filter="external_gene_name",
                 mart=hg)

tcf_peaks<-merge(tcf_peaks_dt,tcf_peaks,by.y="external_gene_name",by.x="gene",all.x = TRUE)

length(unique(tcf_peaks$gene))

tcf_peaks_tab<-data.frame(TCF=unique(tcf_peaks[,colnames(tcf_peaks) %in% c("ensembl_gene_id")]))
tcf_target<-data.frame(table(tcf_peaks_tab$TCF))
row.names(tcf_target)<-tcf_target$Var1
tcf_target$Var1<-NULL
colnames(tcf_target)<-("TCF")


# genes with peaks in LEF
lef_peaks_dt<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_LEFmerge_nearest.genes.txt",header = TRUE)
lef_peaks<-getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                 values= unique(lef_peaks_dt$gene),
                 filter="external_gene_name",
                 mart=hg)

lef_peaks<-merge(lef_peaks_dt,lef_peaks,by.y="external_gene_name",by.x="gene",all.x = TRUE)

lef_peaks_tab<-data.frame(LEF=unique(lef_peaks[,colnames(lef_peaks) %in% c("ensembl_gene_id")]))
lef_target<-data.frame(table(lef_peaks_tab$LEF))
row.names(lef_target)<-lef_target$Var1
lef_target$Var1<-NULL
colnames(lef_target)<-("LEF")

# RNASeq with chipseq with correlation values and model correlated genes
row.names(corout)<-corout$Gene

# RPMI
bcat_allpeaks<-rbind(bcat_cpeaks,bcat_lpeaks)
dim(bcat_allpeaks)
length(unique(bcat_allpeaks$Name))
length(unique(bcat_allpeaks$external_gene_name))

# filter genes which peak is found in > replicates
rep_bcat<-read.delim("/Volumes/cancer/bcat_Project/peakcall/bed_inter/unique_repl_RD.txt",header = FALSE)
colnames(rep_bcat)<-c("chr","Start","End","Name","Score")

bcat_allpeaks<-bcat_allpeaks[!(bcat_allpeaks$Name %in% rep_bcat$Name),]

rpmi_al<-merge(bcat_allpeaks[!duplicated(bcat_allpeaks$ensembl_gene_id),],res_df_rna_rpmi,by.x="ensembl_gene_id",by.y="id",all=TRUE)
rpmi_al<-merge(rpmi_al,est_mod,by.x="ensembl_gene_id",by.y="Gene",all=TRUE)
rpmi_al<-merge(rpmi_al,corout,by.x="ensembl_gene_id",by.y="Gene",all=TRUE)

#variables
rpmi_al$peak<-ifelse((!is.na(rpmi_al$external_gene_name.x)),c("peak"),c("nopeak"))
rpmi_al$deg<-ifelse(!is.na(rpmi_al$pvalue) & rpmi_al$pvalue<0.05,c("DEG"),c("noDEG"))
rpmi_al$cor<-ifelse(!is.na(rpmi_al$Rho) & abs(rpmi_al$Rho)>=0.4 & rpmi_al$Pval<0.01, c("cor"),c("noncor"))
rpmi_al$fit<-ifelse(!is.na(rpmi_al$as_cont_RPMI) & rpmi_al$as_cont_RPMI=="YES" & rpmi_al$as_sh_RPMI=="YES",c("fit"),c("nofit"))

rpmi_al$comb<-paste(rpmi_al$peak,rpmi_al$deg,sep="_")
rpmi_al$comb<-paste(rpmi_al$comb,rpmi_al$cor,sep="_")
rpmi_al$comb<-paste(rpmi_al$comb,rpmi_al$fit,sep="_")

row.names(rpmi_al)<-rpmi_al$ensembl_gene_id

table(rpmi_al$comb)

## number of genes per each category
clas_gene<-data.frame(class=c("Peaks",
                              "DEGs",
                              "coexpression",
                              "fit_model"),
                      frequency=c(dim(rpmi_al[rpmi_al$peak=="peak",])[1],
                                  dim(rpmi_al[rpmi_al$deg=="DEG",])[1],
                                  dim(rpmi_al[rpmi_al$cor=="cor",])[1],
                                  dim(rpmi_al[rpmi_al$fit=="fit",])[1]))

clas_gene$class<- factor(clas_gene$class, levels = clas_gene$class[order(clas_gene$frequency, decreasing = TRUE)])


ggplot(clas_gene, aes(x=class, y=frequency)) +
  geom_bar(stat="identity", color="black",fill="orange",alpha=0.5)+
  geom_text(aes(label=frequency), vjust=-0.1, color="black", size=10)+
  theme_minimal()+
  theme(axis.text.x=element_text(size=14,angle=90),
        axis.text.y=element_text(size=14),
        axis.title=element_blank())



bardata<-data.frame(table(rpmi_al$comb))
colnames(bardata)<-c("combination","Freq")

bardata$combination<- factor(bardata$combination, levels = bardata$combination[order(bardata$Freq, decreasing = FALSE)])

ggplot(bardata[bardata$combination!="nopeak_noDEG_noncor_nofit",], aes(x=combination, y=Freq)) +
  geom_bar(stat="identity", color="black",fill="#56B4E9",alpha=0.5)+
  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=6)+
  theme_minimal()+
  theme(axis.text.x=element_text(size=10,angle=90,face="bold.italic"),
        axis.text.y=element_text(size=14,face="bold.italic"),
        axis.title=element_blank())

rpmi_al[rpmi_al$comb=="peak_DEG_cor_fit",]$external_gene_name.x
rpmi_al[rpmi_al$comb=="peak_DEG_cor_nofit",]$external_gene_name.x


## VIOLIN PLOT WITH BCAT,TCF AND LEF PEAKS AND RNASEQ FC shBCAT
unique(tcf_peaks$ensembl_gene_id)

tcf_peaksdup<-tcf_peaks[!duplicated(tcf_peaks$ensembl_gene_id), ]
#tcf_peaksdup<-tcf_peaks[!duplicated(tcf_peaks$gene), ]

tcf_peaksdup<-tcf_peaksdup[,c(1,6,12)]
tcf_peaksdup$TCF<-"TCF"

lef_peaksdup<-lef_peaks[!duplicated(lef_peaks$ensembl_gene_id), ]
#lef_peaksdup<-lef_peaks[!duplicated(lef_peaks$gene), ]

lef_peaksdup<-lef_peaksdup[,c(1,6,12)]
lef_peaksdup$LEF<-"LEF"


#rpmi_al<-merge(rpmi_al,tcf_peaksdup,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.x="TRUE")
#rpmi_al<-merge(rpmi_al,lef_peaksdup,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.x="TRUE")

#rpmi_al$peaks<-paste(rpmi_al$peak,rpmi_al$TCF,sep="_")
#rpmi_al$peaks<-paste(rpmi_al$peaks,rpmi_al$LEF,sep="_")

#ggplot(rpmi_al,aes(x=peaks,y=log2FoldChange,label=gene.x))+
#  geom_jitter(data=rpmi_al[rpmi_al$padj>0.1 & !is.na(rpmi_al$padj),],color="grey")+
#  geom_jitter(data=rpmi_al[rpmi_al$padj<=0.1 & !is.na(rpmi_al$padj),],color="red")+
#  #geom_label_repel(data=rpmi_al[rpmi_al$padj<=0.05 & !is.na(rpmi_al$padj),],aes(label=gene.y),size=2)+
#  facet_grid(~peak,scales= "free_x")+
#  geom_violin(data=rpmi_al[rpmi_al$padj>0.1 & !is.na(rpmi_al$padj),],color="grey",alpha=0.3)+
#  geom_violin(data=rpmi_al[rpmi_al$padj<=0.1 & !is.na(rpmi_al$padj),],color="red",alpha=0.3)+
#  theme_bw()


#####

## b-cat signature
#####

## Genes with peak b-cat

#rpmi_al has all the data peaks, DEGs, cor, fit
rpmi_al
length(unique(rpmi_al$external_gene_name.y))
length(unique(rpmi_al$ensembl_gene_id))

# output for ssGSEA all data
ssGSEA<-as.data.frame(t(merge(metadata,vastout_t,by="row.names")))
colnames(ssGSEA)<-ssGSEA[1,]
ssGSEA<-ssGSEA[19:nrow(ssGSEA),]
ssGSEA<-merge(geneid,ssGSEA,by.x="Gene",by.y="row.names")
ssGSEA$Gene<-NULL
ssGSEA$Description<-'NA'
ssGSEA<-ssGSEA[,c(1,ncol(ssGSEA),2:(ncol(ssGSEA)-1))]
dim(ssGSEA)
write.table(ssGSEA,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/ggsea_bcatexp_allgenes_incTLE.gct",sep="\t",row.names=FALSE,quote = FALSE)

# All Samples
# excluding genes with no conditions
bcatexp<-vastout_t[,colnames(vastout_t) %in% rpmi_al[rpmi_al$comb!="nopeak_noDEG_noncor_nofit",]$ensembl_gene_id,]

#including genes with condition peak
bcatexp<-vastout_t[,colnames(vastout_t) %in% rpmi_al[rpmi_al$peak=="peak",]$ensembl_gene_id,]


dim(bcatexp)
colcol<-(rpmi_al[rpmi_al$ensembl_gene_id %in% colnames(bcatexp),c(6,11,9,16,49,50,52:54)])
row.names(colcol)<-colnames(bcatexp)
dim(colcol)
colcol$dss<-colcol$Start-colcol$genestart

#merge with TARGET correlation networwk, not all
row.names(corout)<-corout$Gene
colnames(corout)<-c("Gene","NAME","Pval_target","Rho_target")
colcol<-merge(colcol,corout,by="row.names")
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL
colcol$cortarget<-ifelse(!is.na(colcol$Pval_target) & abs(colcol$Rho_target)>=0.4 & colcol$Pval_target<0.01, c("cortar"),c("noncortar"))

colcol$corsen<-ifelse(colcol$Rho_target>0,'pos','neg')
colcol$corsen<-paste(colcol$cortarget,colcol$corsen,sep="_")


# KAISO, TCF and LEF targets

row.names(colcol)
row.names(kaiso_target)
dim(colcol)

colcol<-merge(colcol,kaiso_target,by="row.names",all.x=TRUE)
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL
colcol<-merge(colcol,tcf_target,by="row.names",all.x=TRUE)
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL
colcol<-merge(colcol,lef_target,by="row.names",all.x=TRUE)
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL

colcol$kaiso[is.na(colcol$kaiso)]<-0
colcol$kaiso[colcol$kaiso>0]<-1
colcol$kaiso<-as.factor(colcol$kaiso)

colcol$TCF[is.na(colcol$TCF)]<-0
colcol$TCF[colcol$TCF>0]<-1
colcol$TCF<-as.factor(colcol$TCF)

colcol$LEF[is.na(colcol$LEF)]<-0
colcol$LEF[colcol$LEF>0]<-1
colcol$LEF<-as.factor(colcol$LEF)

# markers matrix with no NAs
complete.cases(bcatexp)
is.na(bcatexp)
bcatexp<-bcatexp[complete.cases(bcatexp),]
dim(bcatexp)

length(colnames(bcatexp))
length((geneid[geneid$Gene %in% colnames(bcatexp),])$NAME)

# targets with no expression
noexp<-colnames(bcatexp[, which(numcolwise(sum)(bcatexp) ==0)])
noexp

# targets with expressed,Remove genes with sum 0 in all samples
exp<-colnames(bcatexp[, which(numcolwise(sum)(bcatexp) !=0)])
exp

bcatexp<-bcatexp[,colnames(bcatexp) %in% exp]

row.names(metadata)<-metadata$SampID
bcatexp<-merge(subset(metadata,select=c(Source,Pheno,Type,GSEA,First_event,Event_free_surv_days,Vital_status,Lesion,Age_Dx_days)),bcatexp,by="row.names",sort=FALSE)
dim(metadata)
dim(bcatexp)

## Add mutations metadata info
somatic_target_bcat<-somatic_target
somatic_target_bcat$Tumor_Sample_Barcode<-gsub('-0.*','',somatic_target_bcat$Tumor_Sample_Barcode)
somatic_target_bcat<-merge(somatic_target_bcat,metadata,by.x="Tumor_Sample_Barcode",by.y="Sample")

somatic_target_bcat_notch<-as.data.frame(table(somatic_target_bcat[grepl('NOTCH1',somatic_target_bcat$Hugo_Symbol),]$SampID))
row.names(somatic_target_bcat_notch)<-somatic_target_bcat_notch$Var1
colnames(somatic_target_bcat_notch)<-c('Samplemut','NOTCH1_freq')

somatic_target_bcat_fbxw7<-as.data.frame(table(somatic_target_bcat[grepl('FBXW7',somatic_target_bcat$Hugo_Symbol),]$SampID))
row.names(somatic_target_bcat_fbxw7)<-somatic_target_bcat_fbxw7$Var1
colnames(somatic_target_bcat_fbxw7)<-c('Samplemut','FBXW7_freq')

somatic_target_bcat_lef<-as.data.frame(table(somatic_target_bcat[grepl('LEF1',somatic_target_bcat$Hugo_Symbol),]$SampID))
row.names(somatic_target_bcat_lef)<-somatic_target_bcat_lef$Var1
colnames(somatic_target_bcat_lef)<-c('Samplemut','LEF1_freq')

somatic_target_mut<-merge(somatic_target_bcat_notch,somatic_target_bcat_fbxw7,by="Samplemut",all=TRUE)
somatic_target_mut<-merge(somatic_target_mut,somatic_target_bcat_lef,by="Samplemut",all=TRUE)
row.names(somatic_target_mut)<-somatic_target_mut$Samplemut

row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL


somatic_target_mut$NOTCH1_freq[is.na(somatic_target_mut$NOTCH1_freq)] = "0"
somatic_target_mut$FBXW7_freq[is.na(somatic_target_mut$FBXW7_freq)] = "0"
somatic_target_mut$LEF1_freq[is.na(somatic_target_mut$LEF1_freq)] = "0"


somatic_target_mut$NOTCH1_freq<-as.numeric(somatic_target_mut$NOTCH1_freq)
somatic_target_mut$LEF1_freq<-as.numeric(somatic_target_mut$LEF1_freq)
somatic_target_mut$FBXW7_freq<-as.numeric(somatic_target_mut$FBXW7_freq)

bcatexp<-merge(somatic_target_mut,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)

row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

bcatexp$LEF1_freq<-as.factor(bcatexp$LEF1_freq)

### Add bcat expression variable
bcat_lev<-data.frame(row.names = row.names(vastout_t_TLE),
                     bcat = vastout_t_TLE[,colnames(vastout_t_TLE) == "ENSG00000168036"])

bcatexp<-merge(bcat_lev,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

### Add kaiso expression variable
kaiso_lev<-data.frame(row.names = row.names(vastout_t_TLE),
                      kaiso_exp = vastout_t_TLE[,colnames(vastout_t_TLE) == "ENSG00000177485"])

bcatexp<-merge(kaiso_lev,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

### Add MYC expression variable
myc_lev<-data.frame(row.names = row.names(vastout_t_TLE),
                    myc_exp = vastout_t_TLE[,colnames(vastout_t_TLE) == "ENSG00000136997"])

bcatexp<-merge(myc_lev,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL


### Add BRCA1 expression variable
brca1_lev<-data.frame(row.names = row.names(vastout_t_TLE),
                      brca1_exp = vastout_t_TLE[,colnames(vastout_t_TLE) == "ENSG00000012048"])

bcatexp<-merge(brca1_lev,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

# add gender info and overall survival
bcatexp<-merge(metadata_target[,c(2:5)],bcatexp,by="row.names",all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

#write.table(bcatexp,"/Users/yguillen/Downloads/bcat_RD_target_genes_expression_TALL.txt",quote = FALSE,row.names = FALSE,sep="\t")

### Using all patients, no matter the first event status
library(RColorBrewer)
cols <- colorRampPalette(rev(brewer.pal(10,name="RdBu")))(50)

paletteLength<-length(cols)

#cols <- colorRampPalette(c("red", "green", "blue"))(50)

## Add or not GSE69239 ( | bcatexp$GSEA=="GSE69239")

myBreaks <- c(seq(min(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])]))),
               annotation_col = as.data.frame(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,c(1,2,3,5,7,8,10,11,12,14,17,19)]),
               annotation_row = colcol[,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 3,
               cutree_rows = 10,
               color=cols, breaks=myBreaks)

clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=10)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=3)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==4,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(bcatexp,patients,by="row.names")
ggplot(patinfo,aes(x=as.factor(cluster),y=bcat))+
  geom_boxplot(aes(color=as.factor(cluster)))+
  theme_bw()

# boxplot expression bcat, tcf, lef and kaiso and brca1 clusters patients
clust_bcat<-data.frame(cluster="cluster1",
                       exp_bcat=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000168036,
                       exp_kaiso=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000177485,
                       exp_tcf=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000081059,
                       exp_lef=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000138795,
                       exp_brca1=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000012048)

clust_bcat<-rbind(clust_bcat,
                  data.frame(cluster="cluster2",
                             exp_bcat=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000168036,
                             exp_kaiso=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000177485,
                             exp_tcf=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000081059,
                             exp_lef=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000138795,
                             exp_brca1=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000012048))

clust_bcat<-rbind(clust_bcat,
                  data.frame(cluster="cluster3",
                             exp_bcat=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000168036,
                             exp_kaiso=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000177485,
                             exp_tcf=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000081059,
                             exp_lef=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000138795,
                             exp_brca1=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000012048))

clust_bcat<-melt(clust_bcat)

ggplot(clust_bcat,aes(x=cluster,y=value,fill=value))+
  geom_boxplot(aes(fill=cluster),alpha=0.3)+
  facet_wrap(~variable,scales="free",nrow = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle=45,hjust = 1),
        axis.text.y = element_text(size=14),
        strip.text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set3")


write.table(clusters,"/Users/yguillen/Downloads/clusters_TARGET_bcatpeaks.txt",sep="\t",quote = FALSE,row.names = FALSE)


## Create correlation across genes using only expression data in patients with no NAs
bcatexp_survival <-bcatexp[bcatexp$GSEA=="TARGET",]
table(bcatexp_survival$First_event)

#Discard NAs and Censored
bcatexp_survival<-bcatexp_survival[(bcatexp_survival$First_event!="NA"),]
bcatexp_survival<-bcatexp_survival[(bcatexp_survival$First_event!="Censored"),]
names<-row.names(bcatexp_survival)
length(names)

# Remission = 0, no remission = 1
bcatexp_survival<-bcatexp_survival %>% mutate(disease_free =  ifelse(First_event == "Remission", 0, 1))
bcatexp_survival$disease_free<-as.factor(bcatexp_survival$disease_free)
dim(bcatexp_survival)
str(bcatexp_survival)

table(bcatexp_survival$disease_free)
row.names(bcatexp_survival)<-names

# Create variable brca1 expression
bcatexp_survival$brcagroup<-NA
bcatexp_survival$brcagroup<-ifelse(bcatexp_survival$brca1_exp>median(bcatexp_survival$brca1_exp),
                                   "high_brca1",
                                   "low_brca1")


# Create variable kaiso expression
bcatexp_survival$kaisogroup<-NA
bcatexp_survival$kaisogroup<-ifelse(bcatexp_survival$kaiso_exp>median(bcatexp_survival$kaiso_exp),
                                    "high_kaiso",
                                    "low_kaiso")

# Create variable bcat expression
bcatexp_survival$bcatgroup<-NA
bcatexp_survival$bcatgroup<-ifelse(bcatexp_survival$bcat>median(bcatexp_survival$bcat),
                                   "high_bcat",
                                   "low_bcat")


library(RColorBrewer)
cols <- colorRampPalette(rev(brewer.pal(10,name="RdBu")))(50)

paletteLength<-length(cols)

bcatexp_survival$Age_Dx_days<-as.numeric(bcatexp_survival$Age_Dx_days)
bcatexp_survival$Age_Dx_days<-bcatexp_survival$Age_Dx_days/365
colnames(bcatexp_survival)[21]<-c("Age_Dx_years")

#cols <- colorRampPalette(c("red", "green", "blue"))(50)
myBreaks <- c(seq(min(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)]))),
               annotation_col = as.data.frame(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",c(7,8,14,17,19)]),
               annotation_row = colcol[,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 5,
               cutree_rows = 7,
               color=cols, breaks=myBreaks)

clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=7)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=5)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==2,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(bcatexp_survival,patients,by="row.names")

patinfo$group<-patinfo$cluster
patinfo$group<-gsub('4','1',patinfo$group)
patinfo$group<-gsub('5','3',patinfo$group)

ggplot(patinfo[patinfo$First_event!="Second Malignant Neoplasm",],aes(x=reorder(as.factor(group),bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3)+
  geom_boxplot(alpha=0.2,lwd=1)+
  #  geom_hline(yintercept = median(scale(bcatexp_survival$brca1_exp))+0.1,color="grey",lwd=1,linetype = "dotted")+
  #  geom_hline(yintercept = median(scale(bcatexp_survival$brca1_exp))-0.1,color="grey",lwd=1,linetype = "dotted")+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[4],color="grey",lwd=2)+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))

patinfo$scaled_brca1<-as.numeric(scale(patinfo$brca1_exp))
median(patinfo$scaled_brca1)
median(patinfo[patinfo$cluster=="1",]$scaled_brca1)

## survival curves with all patients, REDO PATIENTS WITH HEATMAP

bcatexp_survival_cluster<-merge(patients,bcatexp_survival,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Row.names<-NULL


library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)
bcatexp_survival_cluster$Vital_status<-ifelse(bcatexp_survival_cluster$Vital_status=="Alive",0,1)
bcatexp_survival_cluster$Vital_status<-as.integer(bcatexp_survival_cluster$Vital_status)

# Create variable group of clusters
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 1,
                                cluster == 2 ~ 2,
                                cluster == 3 ~ 3,
                                cluster == 4 ~ 1,
                                cluster == 5 ~ 3
                                #                                cluster == 6 ~ 23679,
                                #                                cluster == 7 ~ 23679,
                                #                                cluster == 8 ~ 1854,
                                #                                cluster == 9 ~ 23679
  ))

# Use vital_status or diseas_free and overall_survival_time_days or Event_Free_survival_days
surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])
surv_pvalue(surv_C_fit_TALL)
ggsurvplot(surv_C_fit_TALL, data = bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",], 
           pval = TRUE,
           size=2,
           palette = c("red", "gray87","gray60"),
           ggtheme = theme_bw()+theme(legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      text = element_text(size=20),
                                      axis.text.x = element_text(angle=45, size=20,hjust=1),
                                      axis.text.y = element_text(size=20,hjust=1)),
           font.x = c(20),
           font.y = c(20))

pairwise_survdiff(Surv(Event_free_surv_days, disease_free) ~ cluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])


library(peperr)
bcatexp_survival_cluster$Gender_surv<-ifelse(bcatexp_survival_cluster$Gender=="Female",1,2)
bcatexp_survival_cluster$newcluster_surv<-gsub('235','high_brca1_clusters',bcatexp_survival_cluster$newcluster)
bcatexp_survival_cluster$newcluster_surv<-gsub('1749','medium_brca1_clusters',bcatexp_survival_cluster$newcluster_surv)
bcatexp_survival_cluster$newcluster_surv<-gsub('6810','low_brca1_clusters',bcatexp_survival_cluster$newcluster_surv)

bcatexp_survival_cluster$newcluster_surv<-ifelse(bcatexp_survival_cluster$newcluster=="235","high_brca1_cl","medium/low_brca1_cl")


fit.coxph<- coxph(Surv(Event_free_surv_days, disease_free) ~ newcluster_surv+bcatgroup+kaisogroup+brcagroup+Pheno, 
                  data = bcatexp_survival_cluster[,bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm"],
                  ties = 'exact')
summary(fit.coxph)

ggforest(fit.coxph)


# by brcagroup or kaiso group, remove genes with zero counts

#brca high or low, removing or not second malignant neoplasm (& bcatexp_survival$First_event != "Second Malignant Neoplasm")
i <- (colSums(bcatexp_survival[bcatexp_survival$bcatgroup=="low_bcat",][,20:(ncol(bcatexp_survival)-4)], na.rm=T) != 0)
matnonzero <- bcatexp_survival[bcatexp_survival$bcatgroup=="low_bcat",][,20:(ncol(bcatexp_survival)-4),][, i]

myBreaks1<- c(seq(min(t(scale(matnonzero))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(matnonzero))), length.out=round(floor(paletteLength*0.10)))) 

matnonzero<-cbind(bcatexp_survival[bcatexp_survival$bcatgroup=="low_bcat",][,c(1,3,4,6,7,8,10,13,14,15,(ncol(bcatexp_survival)-3):ncol(bcatexp_survival))],matnonzero)

phet<-pheatmap(as.data.frame(t(scale(matnonzero[,20:ncol(matnonzero)]))),
               annotation_col = as.data.frame(matnonzero[,c(1:14)]),
               annotation_row = colcol[row.names(colcol) %in% colnames(matnonzero[,14:ncol(matnonzero)]),c(6,11,12,13,14)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 7,
               cutree_rows = 4,
               color=cols,
               breaks = myBreaks1)

clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=4)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=7)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
sort(clusters[clusters$cluster==2,]$NAME)

patinfo<-merge(matnonzero,patients,by="row.names")
patinfo$cluster<-as.factor(patinfo$cluster)
table(patinfo$cluster,patinfo$Vital_status)

ggplot(patinfo,aes(x=reorder(as.factor(cluster), brca1_exp, FUN = median),y=scale(brca1_exp)))+
  geom_point(aes(color=as.factor(First_event)),size=3)+
  geom_boxplot(alpha=0.2,lwd=1)+
  geom_hline(yintercept = median(scale(matnonzero$brca1_exp)),color="red",lwd=2)+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[4],color="grey",lwd=2)+
  geom_hline(yintercept = quantile(scale(matnonzero$brca1_exp))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))



## survival curves with selection of groups of patients

bcatexp_survival_cluster<-merge(patients,matnonzero,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Row.names<-NULL

table(bcatexp_survival_cluster$cluster,bcatexp_survival_cluster$disease_free)
table(bcatexp_survival_cluster$cluster)

library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)

# Create variable group of clusters brca high
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 1,
                                cluster == 2 ~ 23457,
                                cluster == 3 ~ 23457,
                                cluster == 4 ~ 23457,
                                cluster == 5 ~ 23457,
                                cluster == 6 ~ 6,
                                cluster == 7 ~ 23457
  ))

surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ brcagroup, data=bcatexp_survival_cluster)
ggsurvplot(surv_C_fit_TALL, data = bcatexp_survival_cluster, pval = TRUE,
           font.x = c(20),
           font.y = c(20))
#           palette = c("green", "indianred","lightseagreen","purple"))



### DEGs clusters within high brca or
## DEGs remission vs relapse

#Groups
# T-ALL TARGET (preliminar)
countdata<-tallcount
row.names(countdata)<-countdata$ID
countdata$NAME<-NULL
countdata$ID<-NULL

countdata<-as.matrix(countdata[, colnames(countdata) %in% rownames(patients)])
countdata[is.na(countdata)] <-0
codata<-as.data.frame(bcatexp_survival[rownames(bcatexp_survival) %in% rownames(patients),][,c(1,5,7,8,10:12,14,17,19,(ncol(bcatexp_survival)-3):ncol(bcatexp_survival))])

# using clusters patients brca1 group heatmap
codata<-merge(patients,codata,by="row.names")
row.names(codata)<-codata$Row.names
codata$Row.names<-NULL
#highest BRCA1
codata[codata$cluster==3,]$cluster<-235
codata[codata$cluster==2,]$cluster<-235
codata[codata$cluster==5,]$cluster<-235

#lowest BRCA1
codata[codata$cluster!=235,]$cluster<-14678910
#codata[codata$cluster==6,]$cluster<-6810
#codata[codata$cluster==8,]$cluster<-6810
#codata[codata$cluster==10,]$cluster<-6810

#codata[codata$cluster!=235,]$cluster<-6810

# or using REMISSION vs RELAPSE without second malignant neoplasm
codata<-as.data.frame(bcatexp_survival[rownames(bcatexp_survival) %in% rownames(patients),][,c(1,5,7,8,10:12,14,17,19,(ncol(bcatexp_survival)-3):ncol(bcatexp_survival))])
codata<-codata[codata$First_event!="Second Malignant Neoplasm",]
dim(codata)
#row.names(codata)<-codata$Row.names
#codata$Row.names<-NULL
codata$disease_free<-as.factor(codata$disease_free)

table(codata$disease_free,codata$First_event)

# make deseq object for clusters
TARGETCountTable <- DESeqDataSetFromMatrix(
  countData=countdata[,colnames(countdata) %in% row.names(codata[codata$cluster=="235" | codata$cluster=="14678910",])],
  colData = codata[codata$cluster=="235" | codata$cluster=="14678910",],
  ~ cluster)

#make deseq object for disease free REMISSION vs RELAPSE TARGET
TARGETCountTable <- DESeqDataSetFromMatrix(
  countData=countdata[,colnames(countdata) %in% row.names(codata)],
  colData = codata,
  ~ disease_free)


colnames(TARGETCountTable)

#how many genes we capture, counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(TARGETCountTable)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

### random sample from the count matrix
nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]


### NORMALIZATION
# Remove genes with low expression levels ()
keep <- rowSums(counts(TARGETCountTable)) >= 10
TARGETCountTable <- TARGETCountTable[keep,]

GeneCounts <- counts(TARGETCountTable)

#colData(TARGETCountTable)$cluster <- factor(colData(TARGETCountTable)$cluster, levels = c("1","2"))

#### estimate size factors
TARGETCountTable <- estimateSizeFactors(TARGETCountTable)
sizeFactors(TARGETCountTable)


### produce rlog-transformed data
rldTARGET <- varianceStabilizingTransformation(TARGETCountTable, blind=TRUE) ## create a distance matrix between the samples


###### Differential EXPRESSION ANALYSIS ###
TARGETCountTable <- estimateDispersions(TARGETCountTable)
plotDispEsts(TARGETCountTable)


# Statistical testing of DE genes
levels(TARGETCountTable$disease_free)
design(TARGETCountTable)

dds<-DESeq(TARGETCountTable)
resultsNames(dds)
rescondition<-results(dds)

# Example plot counts
dds$cluster<-as.factor(dds$cluster)
plotCounts(dds, gene=which(rownames(rescondition)=="ENSG00000165323"), intgroup="disease_free")

head(rescondition)
summary(rescondition)
sum(rescondition$pvalue < 0.05,na.rm=TRUE)
sum(rescondition$padj < 0.1,na.rm=TRUE)
plotMA(rescondition)

# Selecting specific DEGs
target_prog_degs<-as.data.frame(rescondition)
target_prog_degs<-merge(target_prog_degs,geneid,by.x=0,by.y="Gene",all.x=TRUE)

countab<-as.data.frame(counts(dds, normalized=TRUE))
countab<-t(countab)
countab<-melt(countab)
colnames(countab)<-c("SampID","Gene","value")
countab<-merge(countab,geneid,by="Gene")

countab<-merge(countab,codata,by.x="SampID",by.y="row.names")
genesdown<-subset(target_prog_degs,target_prog_degs$pvalue<=0.05 & target_prog_degs$log2FoldChange<0)
genesdown<-genesdown[order(genesdown$pvalue),]
genesup<-subset(target_prog_degs,target_prog_degs$pvalue<=0.05 & target_prog_degs$log2FoldChange>0)
genesup<-genesup[order(genesup$pvalue),]

genesdown[genesdown$Row.names %in% colnames(bcatexp_survival[,c(21:(ncol(bcatexp_survival)-4))]),]
genesup[genesup$Row.names %in% colnames(bcatexp_survival[,c(21:(ncol(bcatexp_survival)-4))]),]

countab$group<-ifelse(countab$First_event=="Progression" | countab$First_event=="Death" | countab$First_event=="Relapse",
                      "Relapse","Remission")

# Gene signature
ggplot(countab[countab$Gene %in% c(deg_bcattarg),],aes(x=group,y=log(value)))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=First_event),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_shape_manual(values=c(15,16,16,17,18,19,19,20,20,0))+
  scale_color_brewer(palette = "Spectral")+
  ylab("Normalized counts")+
  xlab('')+
  facet_wrap(~NAME,scales = "free",ncol=4)+
  theme_classic()+
  theme(axis.text.x = element_blank())

# Add b-cat, kaiso, lef and tcf
ggplot(countab[countab$Gene %in% c("ENSG00000168036","ENSG00000177485","ENSG00000081059","ENSG00000138795"),],aes(x=group,y=log(value)))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=First_event),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_shape_manual(values=c(15,16,16,17,18,19,19,20,20,0))+
  scale_color_brewer(palette = "Spectral")+
  ylab("Normalized counts")+
  xlab('')+
  facet_wrap(~NAME,scales = "free",ncol=4)+
  theme_classic()


rm(TCeltarg)
rm(nz.counts)

## Differentially expressed genes between 235 vs. 10869417 or relapse vs remission

target_prog_degs_sig<-(subset(target_prog_degs,target_prog_degs$pvalue<0.05))
row.names(target_prog_degs_sig)<-target_prog_degs_sig$Row.names
target_prog_degs_sig$Row.names<-NULL

deg_bcattarg<-row.names(target_prog_degs_sig[row.names(target_prog_degs_sig) %in% colnames(bcatexp_survival[,21:(ncol(bcatexp_survival)-4)]),])

## Common targets bcat and brca1
brca_bcat<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/gene_lists/common_brca1_bcat.txt",header = FALSE)
colnames(brca_bcat)<-"GENE"

brca_bcat_genes<-rpmi_al[rpmi_al$external_gene_name.y %in% brca_bcat$V1 & !is.na(rpmi_al$external_gene_name.y),]$ensembl_gene_id

myBreaks <- c(seq(min(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcattarg]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcattarg]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% deg_bcattarg]))),
               annotation_col = as.data.frame(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",c(7,8,14,17,19,21)]),
               annotation_row = colcol[row.names(colcol) %in% deg_bcattarg,c(17,18,19)],
               show_colnames=F,
               show_rownames = T,
               fontsize = 6,
               cutree_cols = 4,
               cutree_rows = 4,
               color=cols, breaks=myBreaks)


listbcat<-c("ENSG00000008300","ENSG00000064703","ENSG00000065029","ENSG00000093183","ENSG00000096401",
            "ENSG00000109805","ENSG00000124574","ENSG00000125871","ENSG00000133119","ENSG00000139629",
            "ENSG00000164068","ENSG00000175463","ENSG00000105810")


myBreaks <- c(seq(min(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% listbcat]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% listbcat]))), length.out=round(floor(paletteLength*0.10)))) 


phet<-pheatmap(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% listbcat]))),
               annotation_col = as.data.frame(bcatexp_survival[,c(7,8,14,17,19,21)]),
               #               annotation_row = colcol[,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = T,
               fontsize = 6,
               cutree_cols = 6,
               cutree_rows = 2,
               color=cols, breaks=myBreaks)



myBreaks <- c(seq(min(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes]))),
               annotation_col = as.data.frame(bcatexp_survival[,c(5,7,8,10:12,14,17,19)]),
               annotation_row = colcol[row.names(colcol) %in% brca_bcat_genes,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 8,
               cutree_rows = 5,
               color=cols, breaks=myBreaks)

# random selection within DEGs relapse vs remission
listbcat<-sample(row.names(target_prog_degs_sig),23)

myBreaks <- c(seq(min(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% listbcat]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% listbcat]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% listbcat]))),
               annotation_col = as.data.frame(bcatexp_survival[,c(5,7,8,10:12,14,17,19)]),
               #               annotation_row = colcol[,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 15,
               cutree_rows = 3,
               color=cols, breaks=myBreaks)


clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=4)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=4)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==2,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% deg_bcattarg])))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% row.names(target_prog_degs_sig)])))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% deg_bcattarg])))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% row.names(target_prog_degs_sig)])))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(bcatexp_survival,patients,by="row.names")

patinfo$groups<-ifelse(patinfo$cluster=="1" | patinfo$cluster=="4", "P1","P2")


### Add TCF7 and LEF1 expression variable
tcf_lev<-data.frame(row.names = row.names(vastout_t),
                    tcf_exp = vastout_t[,colnames(vastout_t) == "ENSG00000081059"])

lef_lev<-data.frame(row.names = row.names(vastout_t),
                    lef_exp = vastout_t[,colnames(vastout_t) == "ENSG00000138795"])

tcflefexp<-merge(tcf_lev,lef_lev,by="row.names",sort=FALSE,all=TRUE)
row.names(tcflefexp)<-tcflefexp$Row.names
tcflefexp$Row.names<-NULL

row.names(patinfo)<-patinfo$Row.names
patinfo<-merge(patinfo,tcflefexp,by="row.names")
patinfo<-patinfo[,-2]

patinfo$outcome_group<-ifelse(patinfo$First_event=="Progression" | patinfo$First_event=="Death" | patinfo$First_event=="Relapse",
                              "Relapse","Remission")


ggplot(patinfo[patinfo$First_event!="Second Malignant Neoplasm",],aes(x=reorder(outcome_group,bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3,alpha=0.8)+
  geom_boxplot(alpha=0.2,lwd=1)+
  #  geom_hline(yintercept = median(scale(bcatexp_survival$bcat))+0.1,color="grey",lwd=1,linetype = "dotted")+
  # geom_hline(yintercept = median(scale(bcatexp_survival$brca1_exp))-0.1,color="grey",lwd=1,linetype = "dotted")+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[4],color="grey",lwd=2)+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))

patinfo_melt<-subset(patinfo,select=c("Row.names","First_event","outcome_group","groups","bcat","kaiso_exp","tcf_exp","lef_exp"))
patinfo_melt<-melt(patinfo_melt,id.vars=c("Row.names","First_event","outcome_group","groups"))

ggplot(patinfo_melt,aes(x=outcome_group,y=log(value)))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=First_event),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_color_brewer(palette = "Spectral")+
  ylab("log normalized expression")+
  xlab('')+
  facet_wrap(~variable,scales = "free",ncol=4)+
  theme_classic()

wilcox.test(log(patinfo[patinfo$First_event!="Second Malignant Neoplasm",]$lef_exp)~patinfo[patinfo$First_event!="Second Malignant Neoplasm",]$outcome_group)

patinfo$scaled_brca1<-as.numeric(scale(patinfo$brca1_exp))
median(patinfo$scaled_brca1)
median(patinfo[patinfo$cluster=="1",]$scaled_brca1)

## survival curves with all patients, REDO PATIENTS WITH HEATMAP

bcatexp_survival_cluster<-merge(patients,bcatexp_survival,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Sample<-row.names(bcatexp_survival_cluster)
bcatexp_survival_cluster$Row.names<-NULL

library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)
bcatexp_survival_cluster$Vital_status<-ifelse(bcatexp_survival_cluster$Vital_status=="Alive",0,1)
bcatexp_survival_cluster$Vital_status<-as.integer(bcatexp_survival_cluster$Vital_status)

# Create variable group of clusters
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 14,
                                cluster == 2 ~ 23,
                                cluster == 3 ~ 23,
                                cluster == 4 ~ 14
  ))

# Use vital_status or diseas_free and overall_survival_time_days or Event_Free_survival_days
surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])

surv_C_fit_TALL <- survfit(Surv(Overall_survival_time_days, Vital_status) ~ newcluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])

surv_pvalue(surv_C_fit_TALL)
ggsurvplot(surv_C_fit_TALL, 
           pval = TRUE,
           size=2,
           palette = c("red", "gray60"),
           ggtheme = theme_bw()+theme(legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      text = element_text(size=20),
                                      axis.text.x = element_text(angle=45, size=20,hjust=1),
                                      axis.text.y = element_text(size=20,hjust=1)),
           font.x = c(20),
           font.y = c(20))


# Add metadata (Age, WBC count, CNS, phenotype, lession)
metadata_target_add<-read_excel("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/metadata_all.xlsx",sheet="TARGET")
row.names(metadata_target_add)<-metadata_target_add$Sample

datat<-merge(bcatexp_survival_cluster,metadata_target_add,by.x="Sample",by.y="GSM")
row.names(datat)<-datat$Sample
datat$Row.names<-NULL  

datat<-subset(datat,datat$First_event!="Second Malignant Neoplasm")
datat$newcluster_surv<-ifelse(datat$newcluster=="14","C_14","C_23")

datat$Age_Dx<-ifelse(datat$`Age at Diagnosis in years`<10,"<10",">10")
datat$WBC_Dx<-ifelse(datat$`WBC at Diagnosis`<200,"<200",">200")
datat$MRD_29<-ifelse(datat$`MRD Day 29`>0,">0","0")

datat$newcluster_surv <- factor(datat$newcluster_surv, levels = c("C_23","C_14"))
datat$Age_Dx <- factor(datat$Age_Dx, levels = c(">10","<10"))
datat$MRD_29 <- factor(datat$MRD_29, levels = c("0",">0"))

datat$CNS_hr <- gsub('2.*','2',datat$`CNS Status at Diagnosis`)
datat$CNS_hr <- gsub('3.*','3',datat$CNS_hr)

fit.coxph<- coxph(Surv(Event_free_surv_days, disease_free) ~ newcluster_surv+Gender.x+Age_Dx+WBC_Dx+MRD_29+CNS_hr, 
                  data = datat[datat$CNS_hr!=".",],
                  ties = 'exact')
summary(fit.coxph)

ggforest(fit.coxph)
##

write.table(clusters,"/Users/yguillen/Downloads/clusters_TARGET_bcatpeaks.txt",sep="\t",quote = FALSE,row.names = FALSE)

#brcadegs degs between clusters high brca1 target 235 and low brca1 target 6810
brcadegs<-subset(target_prog_degs,target_prog_degs$padj<=0.01)
#targets bcat from brcadets
brcadegs<-brcadegs[brcadegs$Row.names %in% colnames(bcatexp_survival),]

# Display genes in each cluster and order from the heatmap and heatmap only these
# clust1_microarray from script T_ALL_GEOquery.R includes all genes in cluster1 separating patients 145 vs 2367

# targdeg are b-cat targets differentially expressed between 145 and 2367 (deg_clust_GSE14618)
deg_bcat_prog<-unique(deg_clust_GSE14618[deg_clust_GSE14618$external_gene_name %in% targdeg,]$ensembl_gene_id)

# corout_brca matrix with p and R for brca1
corout_brca[corout_brca$Pval<0.05,]$Gene

#brca and kaiso
corout_brca[corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,]$Gene
corout_kai[corout_kai$Pval<0.01 & abs(corout_kai$Rho)>0.2,]$Gene

kaibrcacor<-unique(c(corout_brca[corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,]$Gene,
                     corout_kai[corout_kai$Pval<0.01 & abs(corout_kai$Rho)>0.2,]$Gene))

myBreaks <- c(seq(min(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcat_prog]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcat_prog]))), length.out=round(floor(paletteLength*0.10)))) 

row.names(corout_brca)<-corout_brca$Gene

phetsel<-pheatmap(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcat_prog]))),
                  annotation_col = as.data.frame(bcatexp_survival[,c(1,5,7,8,10:12,14,17,19,(ncol(bcatexp_survival)-3):ncol(bcatexp_survival))]),
                  #         annotation_row = corout_brca[corout_brca$Gene %in% colnames(bcatexp_survival) & corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,][,3:4],
                  show_colnames=F,
                  show_rownames = F,
                  fontsize = 6,
                  cutree_cols =8,
                  cutree_rows = 6,
                  color=cols, breaks=myBreaks)


# deg_clust_GSE14618 has all genes DEGs between the previous two groups
deg_micro<-vastout_t_TLE[,colnames(vastout_t_TLE) %in% unique(genedeg$ensembl_gene_id),]
deg_micro<-deg_micro[row.names(deg_micro) %in% row.names(bcatexp_survival),]

myBreaks <- c(seq(min(t(scale(deg_micro))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(deg_micro))), length.out=round(floor(paletteLength*0.10)))) 


phetsel<-pheatmap(as.data.frame(t(scale(deg_micro))),
                  annotation_col = as.data.frame(bcatexp_survival[,c(5,7,8,10:12,14,17,19)]),
                  annotation_row = genedeg[,c(1,3)],
                  show_colnames=F,
                  show_rownames = F,
                  fontsize = 6,
                  cutree_cols =8,
                  cutree_rows = 6,
                  color=cols, breaks=myBreaks)


clusters<-data.frame(cluster=sort(cutree(phetsel$tree_row, k=6)))
patients<-data.frame(cluster=sort(cutree(phetsel$tree_col, k=8)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==2,]$NAME

# Order genes (clusters)
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(deg_micro)))[phetsel$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% clusters[clusters$cluster==6 | clusters$cluster==3,]$Gene])))[phetsel$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcat_prog])))[phetsel$tree_row[["order"]],]))
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% corout_brca[corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,]$Gene])))[phetsel$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(deg_micro)))[,phetsel$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% clusters[clusters$cluster==6 | clusters$cluster==3,]$Gene])))[,phetsel$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcat_prog])))[,phetsel$tree_col[["order"]]]))
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% corout_brca[corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,]$Gene])))[,phetsel$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(bcatexp_survival,patients,by="row.names")

ggplot(patinfo[patinfo$First_event!="Second Malignant Neoplasm",],aes(x=reorder(as.factor(cluster), bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3)+
  geom_boxplot(alpha=0.2,lwd=1)+
  geom_hline(yintercept = median(scale(bcatexp_survival$bcat)),color="red",lwd=2)+
  geom_hline(yintercept = quantile(scale(bcatexp_survival$bcat))[4],color="grey",lwd=2)+
  geom_hline(yintercept = quantile(scale(bcatexp_survival$bcat))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))



## survival curves with all patients, REDO PATIENTS WITH HEATMAP

bcatexp_survival_cluster<-merge(patients,bcatexp_survival,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Row.names<-NULL


library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)
bcatexp_survival_cluster$Vital_status<-ifelse(bcatexp_survival_cluster$Vital_status=="Alive",0,1)
bcatexp_survival_cluster$Vital_status<-as.integer(bcatexp_survival_cluster$Vital_status)

# Create variable group of clusters
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 13658,
                                cluster == 2 ~ 247,
                                cluster == 3 ~ 13658,
                                cluster == 4 ~ 247,
                                cluster == 5 ~ 13658,
                                cluster == 6 ~ 13658,
                                cluster == 7 ~ 247,
                                cluster == 8 ~ 13658
                                #                                cluster == 9 ~ 49731112,
                                #                                cluster == 10 ~ 2810,
                                #                                cluster == 11 ~ 49731112,
                                #                                cluster == 12 ~ 49731112
  ))


surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster)

surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])

ggsurvplot(surv_C_fit_TALL, pval = TRUE,
           font.x = c(20),
           font.y = c(20))
#           palette = c("green", "indianred","lightseagreen","purple"))

pairwise_survdiff(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster)

surv_C_fit_TALL$time




#### Reducing b-cat targets signature

flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

#correlation matrix with all targets
rcorTALL<-rcorr(as.matrix(bcatexp_survival[,20:(ncol(bcatexp_survival)-4)]), type = "spearman")

#correlation matrix with targets cluster 1 microarray data
rcorTALL<-rcorr(as.matrix(bcatexp_survival[,colnames(bcatexp_survival) %in% clust1_microarray]), type = "spearman")

# correlation matrix
corrplot(rcorTALL$r, order="hclust",tl.col="black", tl.cex = 0.1, type = "lower",
         p.mat = rcorTALL$P, sig.level = 0.05, insig = "blank", method = "color")


ordercorTALL <- flat_cor_mat(rcorTALL$r, rcorTALL$P)
ordercorTALL<-ordercorTALL[!is.na(ordercorTALL$p),]

ordercorTALL<-merge(ordercorTALL,geneid,by.x="row",by.y="Gene",all.x=TRUE)

# selection by p-value (0.85 FOR ALL, 0.8 FOR CLUST1_MICROARRAY SUBSET)
pair1<-as.vector(unique(ordercorTALL[ordercorTALL$p<0.01 & abs(ordercorTALL$cor)>0.8,]$row))
pair2<-as.vector(unique(ordercorTALL[ordercorTALL$p<0.01 & abs(ordercorTALL$cor)>0.8,]$column))

genes<-unique(c(pair1,pair2))
length(genes)


clusterTALL <-bcatexp_survival[,genes]
dim(clusterTALL)

clusterTALL<-t(clusterTALL)
clusterTALL<-merge(geneid,clusterTALL,by.y="row.names",by.x="Gene",all.y=TRUE)
row.names(clusterTALL)<-clusterTALL$NAME
clusterTALL$Gene<-NULL
clusterTALL$NAME<-NULL
clusterTALL<-t(clusterTALL)

rclusterTALL<-rcorr(as.matrix(clusterTALL), type = "pearson")
corrplot(rclusterTALL$r, order="hclust",tl.col="black", tl.cex = 0.5,
         p.mat = rclusterTALL$P, sig.level = 0.01, insig = "blank",addrect = 2)



genes

# remove two elements from genes list x[x != "b"], THESE ARE TWO HEMOGLOBIN SUBUNITS, THEY HAVE A HIGH CORRELATION, BUT ONLY AMONG THEM
genes<-genes[genes != "ENSG00000244734"]
genes<-genes[genes != "ENSG00000086506"]

#clusterTALL <-bcatexp_survival[,genes]
clusterTALL <-bcatexp_survival[,deg_bcattarg]


dim(clusterTALL)

rclusterTALL<-rcorr(as.matrix(clusterTALL), type = "pearson")

cor(as.matrix(clusterTALL))
corrplot(rclusterTALL$r, order="hclust",tl.col="black", tl.cex = 0.5,
         p.mat = rclusterTALL$P, sig.level = 0.05, insig = "blank",addrect = 4,tl.pos = "n")


# Using the reduced gene list from correlation networks
bcat_signature<-genes
sort(geneid[geneid$Gene %in% bcat_signature,]$NAME)

bcat_matrix<-bcatexp_survival[,colnames(bcatexp_survival) %in% bcat_signature,]

# or using gene liste crossing DEGs leukemia vs thymus in scRNASeq data
scsig<-c("ROMO1","SERINC5","TCF7","FAM40A","BCLAF1",
         "RBM39","RPS28","DDX3X","NDUFA11","SLC38A2","UBE2D3")

scsig<-geneid[geneid$NAME %in% scsig,]$Gene

bcat_matrix<-bcatexp_survival[,colnames(bcatexp_survival) %in% scsig,]


bcat_matrix<-merge(subset(bcatexp_survival,select=c(bcat,kaiso_exp,brca1_exp,
                                                    bcatgroup,kaisogroup,brcagroup,disease_free,Vital_status,First_event,Pheno,
                                                    NOTCH1_freq,FBXW7_freq,LEF1_freq,Event_free_surv_days)),bcat_matrix,by="row.names",sort=FALSE)
dim(metadata)
dim(bcat_matrix)

row.names(bcat_matrix)<-bcat_matrix$Row.names
bcat_matrix$Row.names<-NULL


myBreaks1 <- c(seq(min(t(scale(bcat_matrix[,15:ncol(bcat_matrix)]))), 2.5, length.out=round(ceiling(paletteLength*0.90))),
               seq(2.6,max(t(scale(bcat_matrix[,15:ncol(bcat_matrix)]))),length.out = round(ceiling(paletteLength*0.10)))) 

colcolmat<-colcol[row.names(colcol) %in% colnames(as.data.frame(bcat_matrix[,15:ncol(bcat_matrix)])),]

phetsel<-pheatmap(as.data.frame(t(scale(bcat_matrix[,15:ncol(bcat_matrix)]))),
                  annotation_col = as.data.frame(bcat_matrix[,c(1:13)]),
                  annotation_row = colcolmat[,c(14,16,17,18,19)],
                  show_colnames=F,
                  show_rownames = F,
                  fontsize = 6,
                  cutree_cols = 6,
                  cutree_rows = 2,
                  color=cols, breaks=myBreaks1)

# order by bcat group
matplot<-as.data.frame(t(scale(bcat_matrix[,15:ncol(bcat_matrix)])))
matplot[, order(bcat_matrix$bcatgroup)]
matplot<-matplot[,row.names(bcat_matrix[order(bcat_matrix$bcatgroup),])]

pheatmap(matplot,
         annotation_col = as.data.frame(bcat_matrix[,c(1:13)]),
         annotation_row = colcolmat[,c(14,16,17,18,19)],
         show_colnames=F,
         show_rownames = F,
         fontsize = 6,
         cutree_cols = FALSE,
         cluster_cols = FALSE,
         cutree_rows = 5,
         color=cols, breaks=myBreaks1)

##

clusters<-data.frame(cluster=sort(cutree(phetsel$tree_row, k=2)))
patients<-data.frame(cluster=sort(cutree(phetsel$tree_col, k=6)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==1,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcat_matrix[,15:ncol(bcat_matrix)])))[phetsel$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcat_matrix[,15:ncol(bcat_matrix)])))[,phetsel$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)



## survival curves

bcatexp_survival_cluster<-merge(patients,bcat_matrix,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Row.names<-NULL

table(bcatexp_survival_cluster$cluster,bcatexp_survival_cluster$Vital_status)

# Create variable group of clusters
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 13,
                                cluster == 2 ~ 2456,
                                cluster == 3 ~ 13,
                                cluster == 4 ~ 2456,
                                cluster == 5 ~ 2456,
                                cluster == 6 ~ 2456))


library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)
surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster)
ggsurvplot(surv_C_fit_TALL, data = bcatexp_survival_cluster, pval = TRUE,
           font.x = c(20),
           font.y = c(20))
#           palette = c("green", "indianred","lightseagreen","purple"))

pairwise_survdiff(Surv(Event_free_surv_days, disease_free) ~ cluster, data=bcatexp_survival_cluster)

surv_C_fit_TALL$time


# Display genes in each cluster and order from the heatmap and heatmap only these

myBreaks <- c(seq(min(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene]))), length.out=round(floor(paletteLength*0.10)))) 


phetsel<-pheatmap(as.data.frame(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene]))),
                  annotation_col = as.data.frame(bcat_matrix[,c(1:13)]),
                  annotation_row = colcolmat[,c(14,16,17,18,19)],
                  show_colnames=F,
                  show_rownames = T,
                  fontsize = 6,
                  cutree_cols = 4,
                  #  cutree_rows = 3,
                  color=cols, breaks=myBreaks)


patients<-data.frame(cluster=sort(cutree(phetsel$tree_col, k=4)))
matsel<-as.data.frame(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene])))

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(matsel[,phetsel$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)



patinfosel<-merge(t(as.data.frame(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene])))),patients,by="row.names")
row.names(patinfosel)<-patinfosel$Row.names
patinfosel$Row.names<-NULL
patinfosel<-merge(bcat_matrix[,c(1:13)],patinfosel,by="row.names")
row.names(patinfosel)<-patinfosel$Row.names
patinfosel$Row.names<-NULL

patinfosel$cluster<-as.factor(patinfosel$cluster)
table(patinfosel$cluster,patinfosel$Vital_status)

ggplot(patinfosel,aes(x=reorder(as.factor(cluster), bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3)+
  geom_boxplot(alpha=0.2,lwd=1)+
  geom_hline(yintercept = median(scale(bcat_matrix$bcat)),color="red",lwd=2)+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[4],color="grey",lwd=2)+
  geom_hline(yintercept = quantile(scale(bcat_matrix$bcat))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))

# output for ssGSEA
write.table(clusters[clusters$cluster==4,]$NAME,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/cluster4.txt",row.names=FALSE,quote = FALSE)
write.table(clusters[clusters$cluster==5 | clusters$cluster==2,]$NAME,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/cluster2_5.txt",row.names=FALSE,quote = FALSE)
write.table(clusters[clusters$cluster==1 | clusters$cluster==3,]$NAME,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/cluster1_3.txt",row.names=FALSE,quote = FALSE)

# Results ssGSEA
ssGSEAout<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/output/TALL_clusters_mod.gct")
ssGSEAout<-as.data.frame(t(ssGSEAout))
colnames(ssGSEAout)<-ssGSEAout[1,]
ssGSEAout<-ssGSEAout[-c(1:2),]

## For other cohorts, first remove genes with sum zero values
i <- (colSums(bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",22:ncol(bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",])], na.rm=T) != 0)
matnonzero <- bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",22:ncol(bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",])][, i]

matnonzero<-cbind(bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",1:21],matnonzero)

myBreaks1 <- c(seq(min(t(scale(matnonzero[row.names(matnonzero)!="Thymus",22:ncol(matnonzero)]))), 0, length.out=round(ceiling(paletteLength*0.40)) + 1),
               seq(0.1,2,length.out = round(ceiling(paletteLength*0.40))),
               seq(2.1, max(t(scale(matnonzero[row.names(matnonzero)!="Thymus",22:ncol(matnonzero)]))), length.out=round(floor(paletteLength*0.20)-1))) 




phet1<-pheatmap(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus",22:ncol(matnonzero)]))),
                annotation_col = as.data.frame(matnonzero[row.names(matnonzero)!="Thymus",c(5,17)]),
                #         annotation_row = colcol[row.names(colcol) %in% colnames(matnonzero[,21:ncol(matnonzero)]),c(3,6,8,7,9)],
                show_colnames=T,
                show_rownames = F,
                fontsize = 6,
                cutree_cols = 3,
                cutree_rows = 3,
                color=cols,
                breaks = myBreaks1)


myBreaks1 <- c(seq(min(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA"  & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% deg_bcattarg]))), 0, length.out=round(ceiling(paletteLength*0.40)) + 1),
               seq(0.1,2,length.out = round(ceiling(paletteLength*0.40))),
               seq(2.1, max(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA" & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% deg_bcattarg]))), length.out=round(floor(paletteLength*0.20)-1))) 


matnonzero$Age_Dx_days<-(as.numeric(matnonzero$Age_Dx_days))/365
colnames(matnonzero)[21]<-"Age_Dx"

phet1<-pheatmap(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA" & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% deg_bcattarg]))),
                annotation_col = as.data.frame(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA"  & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",c(14,17,19,21)]),
                #                annotation_row = colcol[row.names(colcol) %in% colnames(matnonzero[,colnames(matnonzero) %in% deg_bcattarg]),c(6,16)],
                show_colnames=F,
                show_rownames = T,
                fontsize = 6,
                cutree_cols = 3,
                cutree_rows = 4,
                color=cols,
                breaks = myBreaks1)

phet1<-pheatmap(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA" & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% colnames(bcatexp[,22:ncol(bcatexp)])]))),
                annotation_col = as.data.frame(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA"  & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",c(14,17,19,21)]),
                #                annotation_row = colcol[row.names(colcol) %in% colnames(matnonzero[,colnames(matnonzero) %in% deg_bcattarg]),c(6,16)],
                show_colnames=F,
                show_rownames = F,
                fontsize = 6,
                cutree_cols = 3,
                cutree_rows = 4,
                color=cols,
                breaks = myBreaks1)

phet1<-pheatmap(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA" & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% geneid[geneid$NAME %in% targdeg,]$Gene]))),
                annotation_col = as.data.frame(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA"  & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",c(14,17,19,21)]),
                #                annotation_row = colcol[row.names(colcol) %in% colnames(matnonzero[,colnames(matnonzero) %in% deg_bcattarg]),c(6,16)],
                show_colnames=F,
                show_rownames = F,
                fontsize = 6,
                cutree_cols = 3,
                cutree_rows = 4,
                color=cols,
                breaks = myBreaks1)


clusters1<-data.frame(cluster=sort(cutree(phet1$tree_row, k=4))) 
patients1<-data.frame(cluster=sort(cutree(phet1$tree_col, k=3)))




table(patients1$cluster)
table(clusters1$cluster)
clusters1$Gene<-row.names(clusters1)
clusters1[clusters1$cluster==1,]$Gene

# Order genes (clusters)
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% deg_bcattarg])))[phet1$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA"  & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% deg_bcattarg])))[phet1$tree_row[["order"]],]))
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA"  & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% colnames(bcatexp[,22:ncol(bcatexp)])])))[phet1$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters1<-merge(clusters1,geneorder,by="Gene")
clusters1$order<-as.numeric(clusters1$order)

clusters1<-merge(clusters1,geneid,by="Gene",all.x=TRUE)

# Order patients (patients)
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% deg_bcattarg])))[,phet1$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA" & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% deg_bcattarg])))[,phet1$tree_col[["order"]]]))
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA" & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% colnames(bcatexp[,22:ncol(bcatexp)])])))[,phet1$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients1<-merge(patients1,patientorder,by=0)
row.names(patients1)<-patients1$Row.names
patients1$Row.names<-NULL
patients1$order<-as.numeric(patients1$order)


patinfo1<-merge(matnonzero,patients1,by="row.names")


### Add TCF7 and LEF1 and bcat and kaiso in EGA_TLE expression variable
tcf_lev<-data.frame(row.names = row.names(vastout_t),
                    tcf_exp_EGA = vastout_t[,colnames(vastout_t) == "ENSG00000081059"])

lef_lev<-data.frame(row.names = row.names(vastout_t),
                    lef_exp_EGA = vastout_t[,colnames(vastout_t) == "ENSG00000138795"])

bcat_lev<-data.frame(row.names = row.names(vastout_t),
                     bcat_exp_EGA = vastout_t[,colnames(vastout_t) == "ENSG00000168036"])

kaiso_lev<-data.frame(row.names = row.names(vastout_t),
                      kaiso_exp_EGA = vastout_t[,colnames(vastout_t) == "ENSG00000177485"])

tf_expTLE<-merge(tcf_lev,lef_lev,by="row.names",all=TRUE)
row.names(tf_expTLE)<-tf_expTLE$Row.names
tf_expTLE$Row.names<-NULL
tf_expTLE<-merge(tf_expTLE,bcat_lev,by="row.names",all=TRUE)
row.names(tf_expTLE)<-tf_expTLE$Row.names
tf_expTLE$Row.names<-NULL
tf_expTLE<-merge(tf_expTLE,kaiso_lev,by="row.names",all=TRUE)
row.names(tf_expTLE)<-tf_expTLE$Row.names
tf_expTLE$Row.names<-NULL

row.names(patinfo1)<-patinfo1$Row.names
patinfo1<-merge(patinfo1,tf_expTLE,by="row.names")
patinfo1<-patinfo1[,-2]

patinfo1$outcome_group<-ifelse(patinfo1$First_event=="Progression" | patinfo1$First_event=="Death" | patinfo1$First_event=="Relapse",
                               "Relapse","Remission")


ggplot(patinfo1,aes(x=reorder(as.factor(cluster), bcat, FUN = median),y=scale(tcf_exp_EGA)))+
  geom_point(aes(color=First_event),size=5)+
  geom_boxplot(alpha=0.2,lwd=1)+
  #  scale_color_manual(values=c("chartreuse3","hotpink","deepskyblue"))+
  #  geom_hline(yintercept = median(scale(patinfo$bcat_exp_prob1)),color="red",lwd=2)+
  #  geom_hline(yintercept = median(scale(patinfo$bcat_exp_prob1))-0.2,color="grey",lwd=2)+
  #  geom_hline(yintercept = median(scale(patinfo$kaiso_exp_prob2))+0.25,color="grey",lwd=2)+
  # geom_hline(yintercept = 0,color="grey",lwd=1)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_text(angle=45, size=16,hjust=1),
        axis.title.x = element_blank(),
        panel.border = element_rect(size = 2))

patinfo1$bcat_scale<-as.numeric(scale(patinfo1$bcat_exp_EGA))
patinfo1$kaiso_scale<-as.numeric(scale(patinfo1$kaiso_exp_EGA))
patinfo1$tcf_scale<-as.numeric(scale(patinfo1$tcf_exp_EGA))
patinfo1$lef_scale<-as.numeric(scale(patinfo1$lef_exp_EGA))

wilcox.test(patinfo1[patinfo1$cluster!=3,]$lef_exp_EGA~as.factor(patinfo1[patinfo1$cluster!=3,]$cluster))

kruskal.test(patinfo1$lef_exp_EGA~patinfo1$outcome_group)

patinfo1_melt<-subset(patinfo1,select=c("Row.names","First_event","cluster","outcome_group","bcat_scale","kaiso_scale","tcf_scale","lef_scale"))
patinfo1_melt<-melt(patinfo1_melt,id.vars=c("Row.names","First_event","cluster","outcome_group"))

ggplot(patinfo1_melt,aes(x=outcome_group,y=value))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=First_event),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_color_brewer(palette = "Spectral")+
  ylab("scaled expression")+
  xlab('')+
  facet_wrap(~variable,scales = "free",ncol=4)+
  theme_classic()


prop_clust<-data.frame(prop.table(table(patinfo1[patinfo1$First_event!="NA" | patinfo1$First_event!="Disease" | patinfo1$First_event!="Death" | patinfo1$First_event!="Secondary leukemia",]$cluster,patinfo1[patinfo1$First_event!="NA"  | patinfo1$First_event!="Secondary leukemia" | patinfo1$First_event!="Death",]$First_event),1))
colnames(prop_clust)<-c("cluster","First_event","Proportion")
prop_clust$cluster<-as.factor(prop_clust$cluster)

prop_clust$cluster<- factor(prop_clust$cluster, levels=c("2","3","1"))

ggplot(prop_clust,aes(x=cluster,y=Proportion,fill=First_event))+
  geom_col(color="black")+
  scale_fill_manual(values=c("dodgerblue","green","darkolivegreen","cyan"))+
  theme_bw()+
  theme(legend.position="right",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=16))



## Read ssGSEA from b-cat signature
#UP AND DOWN SEPARATE all genes
ssbcat<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/ssGSEA_output/TALL_ALL.gct")

#UP AND DOWN SEPARATE excluding SAP18,HLA-DOA and KRT5
ssbcat<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/ssGSEA_output/TALL_reduced.gct")

ssbcat<-as.data.frame(t(ssbcat))
colnames(ssbcat)<-ssbcat[1,]
ssbcat<-ssbcat[-c(1,2),]

metadata_ssbcat<-merge(ssbcat,metadata,by=0,all.y=TRUE)
row.names(metadata_ssbcat)<-metadata_ssbcat$Row.names
metadata_ssbcat$Row.names<-NULL

metadata_ssbcat$deg_up_relapse<-as.numeric(metadata_ssbcat$deg_up_relapse)
metadata_ssbcat$deg_down_relapse<-as.numeric(metadata_ssbcat$deg_down_relapse)

cohorts<-c("TARGET","EGA_R","EGA_TLE")

metadata_ssbcat<-subset(metadata_ssbcat,(metadata_ssbcat$First_event!="Second Malignant Neoplasm" & 
                                           metadata_ssbcat$First_event!="Secondary leukemia" &
                                           metadata_ssbcat$First_event!="Disease"))

ggplot(metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts , ],aes(x=deg_up_relapse,y=deg_down_relapse))+
  geom_point(data=metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),color="black",size=4)+
  geom_point(data=metadata_ssbcat[metadata_ssbcat$First_event!="NA" & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),size=3)+
  scale_color_brewer(palette="Set2")+
  geom_density_2d(data=metadata_ssbcat[metadata_ssbcat$First_event!="Remission",],color="darkolivegreen",bins=20)+
  geom_vline(xintercept = 6500)+
  geom_hline(yintercept = 4000)+
  theme_bw()+
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.title.y = element_text(size=16))+
  coord_fixed(ratio=1)


#UNIQUE
ssbcat_unique<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/ssGSEA_output/TALL_ALL_unique.gct")
ssbcat_unique<-as.data.frame(t(ssbcat_unique))
colnames(ssbcat_unique)<-ssbcat_unique[1,]
ssbcat_unique$method<-"ssGSEA"
ssbcat_unique<-ssbcat_unique[-c(1,2),]

metadata_ssbcat<-merge(ssbcat_unique,metadata,by=0,all.y=TRUE)
row.names(metadata_ssbcat)<-metadata_ssbcat$Row.names
metadata_ssbcat$Row.names<-NULL

metadata_ssbcat$bcat_relapse<-as.numeric(metadata_ssbcat$bcat_relapse)

cohorts<-c("TARGET","EGA_R","EGA_TLE")

ggplot(metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts, ],aes(x=First_event,y=bcat_relapse))+
  geom_point(data=metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),color="black",size=4)+
  geom_point(data=metadata_ssbcat[!is.na(metadata_ssbcat$First_event) & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),size=3)+
  geom_boxplot(data=metadata_ssbcat[!is.na(metadata_ssbcat$First_event) & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],alpha=0.1)+
  scale_color_brewer(palette="Set2")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=16))


# Empty space for subsequent data
rm(ty2)
rm(res.pca)
rm(res.pcaTLE)
rm(c1)

#####
