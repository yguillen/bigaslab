#Load libraries
library(data.table)
library(GenomicAlignments)
library(GO.db)
library(DiffBind)
library(ChIPQC)
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(BayesPeak)
library(parallel)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ReactomePA)
library(trackViewer)
library(Gviz)
library(rtracklayer)
library(Sushi)
library(ChIPpeakAnno)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(DO.db)
library(VennDiagram)
library(readxl)
library(enrichR)

library(rJava)
library(venneuler)

setwd("/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/peakcall/")

#RNAseq expression analysis IKK Colomer 2018
rnaCK<-read.delim("/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/RNAseq_BJC_Colomer/list_IKK_RNAseq_genes.txt",dec=",")
#rnaCK<-subset(rnaCK,rnaCK$padj_KO_vs_WT<0.01)
rnaCK$state<-ifelse(rnaCK$log2FoldChange_KO_vs_WT < 0, 
                         c("Down_KO"), c("Up_KO")) 
rnaCKgenesUP<-(subset(rnaCK,rnaCK$state=="Up_KO"))$external_gene_name
rnaCKgenesDOWN<-(subset(rnaCK,rnaCK$state=="Down_KO"))$external_gene_name
allrnaCKgenes<-unique(rnaCK$external_gene_name)

length(unique(rnaCKgenesUP))
length(unique(rnaCKgenesDOWN))

#Top Brd4 targets from RNAseq KO BrD4 and ChIPSEQ 
brd4<-read.delim("/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/papers/brd4_targets_Hussog.txt",header = FALSE)
brd4$V1<-paste0(toupper(substr(brd4$V1, 1, 1)), substr(brd4$V1, 2, nchar(brd4$V1)))


## ChIP for Brd4

## CKC KO IKKa
peakCK1 <- readPeakFile("CKC1_peaks_peaks.narrowPeak")
peakCK2 <- readPeakFile("CKC2_peaks_peaks.narrowPeak")
peakCK3 <- readPeakFile("CKC3_peaks_peaks.narrowPeak")

#Merged CK1 and CK2
peakCK <- readPeakFile("CKC_peaks_peaks.narrowPeak")

# KO with TNF
#peakCKT1 <- readPeakFile("CKT1_peaks_peaks.narrowPeak")
#peakCKT2 <- readPeakFile("CKT2_peaks_peaks.narrowPeak")

#KO IR
peakCKIR <- readPeakFile("CKIR_peaks_peaks.narrowPeak")
#2N run
peakKO1IR <- readPeakFile("KOI1_peaks_peaks.narrowPeak")
peakKO2IR <- readPeakFile("KOI2_peaks_peaks.narrowPeak")
peakKOIR <- readPeakFile("KOI_peaks_peaks.narrowPeak")

## WT
peakCW1 <- readPeakFile("CWC1_peaks_peaks.narrowPeak")
peakCW2 <- readPeakFile("CWC2_peaks_peaks.narrowPeak")
peakCW3 <- readPeakFile("CWC3_peaks_peaks.narrowPeak")

#merged CW1 and CW2
peakCW <- readPeakFile("CWC_peaks_peaks.narrowPeak")

# WT With TNF
#peakCWT1 <- readPeakFile("CWT1_peaks_peaks.narrowPeak")
#peakCWT2 <- readPeakFile("CWT2_peaks_peaks.narrowPeak")

#WT with IR
peakCWIR <- readPeakFile("CWIR_peaks_peaks.narrowPeak")
#2nd run
peakWT1IR <- readPeakFile("WTI1_peaks_peaks.narrowPeak")
peakWT2IR <- readPeakFile("WTI2_peaks_peaks.narrowPeak")
peakWTIR <- readPeakFile("WTI_peaks_peaks.narrowPeak")


## Histone marks
#WT
peakWT3H<- readPeakFile("WTH3_peaks_peaks.broadPeak")
peakWT4H<- readPeakFile("WTH4_peaks_peaks.broadPeak")
peakWTH<- readPeakFile("WTH_peaks_peaks.broadPeak")

#KO merged
peakKOH3 <- readPeakFile("KOH3_peaks_peaks.broadPeak")
peakKOH4 <- readPeakFile("KOH4_peaks_peaks.broadPeak")
peakKOH <- readPeakFile("KOH_peaks_peaks.broadPeak")

#Overall all chromosmes peak coverage
peakbed<-readPeakFile("CWC3_peaks_summits.bed")
covplot(peakbed, weightCol="V5")


#genes database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#hsdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

### PEAK ANNOTATION
#KO
peakAnnoCK1 <- annotatePeak(peakCK1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoCK2 <- annotatePeak(peakCK2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoCK3 <- annotatePeak(peakCK3,tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db") 

# merged CK1 and CK2
peakAnnoCK <- annotatePeak(peakCK, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# KO With TNF
#peakAnnoCKT1 <- annotatePeak(peakCKT1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db" )
#peakAnnoCKT2 <- annotatePeak(peakCKT2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db" )

# KO irradiated
peakAnnoCKIR <- annotatePeak(peakCKIR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db" )
peakAnnoKO1IR <- annotatePeak(peakKO1IR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db" )
peakAnnoKO2IR <- annotatePeak(peakKO2IR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db" )
peakAnnoKOIR <- annotatePeak(peakKOIR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db" )

#WT
peakAnnoCW1 <- annotatePeak(peakCW1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoCW2 <- annotatePeak(peakCW2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoCW3 <- annotatePeak(peakCW3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

#merged CW1 and CW2
peakAnnoCW <- annotatePeak(peakCW, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# WT With TNF
#peakAnnoCWT1 <- annotatePeak(peakCWT1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
#peakAnnoCWT2 <- annotatePeak(peakCWT2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# WT irradiated
peakAnnoCWIR <- annotatePeak(peakCWIR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWT1IR <- annotatePeak(peakWT1IR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWT2IR <- annotatePeak(peakWT2IR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWTIR <- annotatePeak(peakWTIR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")


### Histone marks ###
#WT
peakAnnoWTH3 <- annotatePeak(peakWT3H, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWTH4 <- annotatePeak(peakWT4H, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWTH <- annotatePeak(peakWTH, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

#KO
peakAnnoKOH3 <- annotatePeak(peakKOH3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoKOH4 <- annotatePeak(peakKOH4, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoKOH <- annotatePeak(peakKOH, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

####### Patient-derived tumoroids (human genome txdb) ######
# WT
peakPWT1<-readPeakFile("PWT1_peaks_peaks.narrowPeak")
peakPWT2<-readPeakFile("PWT2_peaks_peaks.narrowPeak")

# KO
peakPKO1<-readPeakFile("PKO1_peaks_peaks.narrowPeak")
peakPKO2<-readPeakFile("PKO2_peaks_peaks.narrowPeak")

#Overall all chromosmes peak coverage
peakbed<-readPeakFile("PWT1_peaks_summits.bed")
covplot(peakbed, weightCol="V5")

#Peak annotation

peakAnnoPWT1 <- annotatePeak(peakPWT1, tssRegion=c(-3000, 3000), TxDb=hsdb, annoDb="org.Hs.eg.db")
peakAnnoPWT2 <- annotatePeak(peakPWT2, tssRegion=c(-3000, 3000), TxDb=hsdb, annoDb="org.Hs.eg.db")

peakAnnoPKO1 <- annotatePeak(peakPKO1, tssRegion=c(-3000, 3000), TxDb=hsdb, annoDb="org.Hs.eg.db")
peakAnnoPKO2 <- annotatePeak(peakPKO2, tssRegion=c(-3000, 3000), TxDb=hsdb, annoDb="org.Hs.eg.db")


##DATAFRAME WITH genes id linked to peaks

#CK1 KO
dataf_peakannoCK1<-as.data.frame(peakAnnoCK1)
dataf_peakannoCK1$annotation<-gsub('Intron.*','Intron',dataf_peakannoCK1$annotation)
dataf_peakannoCK1$annotation<-gsub('Exon.*','Exon',dataf_peakannoCK1$annotation)
dataf_peakannoCK1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCK1$annotation)
dataf_peakannoCK1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCK1$annotation)
dataf_peakannoCK1$peaksite<-paste(dataf_peakannoCK1$SYMBOL,dataf_peakannoCK1$annotation,sep='_')
table(dataf_peakannoCK1$annotation)

#Remove distal intergenic
#dataf_peakannoCK1<-subset(dataf_peakannoCK1,dataf_peakannoCK1$annotation!="Distal Intergenic")
#table(dataf_peakannoCK1$annotation)

CK1_genes<-unique(dataf_peakannoCK1$SYMBOL)
CK1_peaks<-unique(dataf_peakannoCK1$peaksite)

#CK2 KO
dataf_peakannoCK2<-as.data.frame(peakAnnoCK2)
dataf_peakannoCK2$annotation<-gsub('Intron.*','Intron',dataf_peakannoCK2$annotation)
dataf_peakannoCK2$annotation<-gsub('Exon.*','Exon',dataf_peakannoCK2$annotation)
dataf_peakannoCK2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCK2$annotation)
dataf_peakannoCK2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCK2$annotation)
dataf_peakannoCK2$peaksite<-paste(dataf_peakannoCK2$SYMBOL,dataf_peakannoCK2$annotation,sep='_')
table(dataf_peakannoCK2$annotation)

#Remove distal intergenic
#dataf_peakannoCK2<-subset(dataf_peakannoCK2,dataf_peakannoCK2$annotation!="Distal Intergenic")
#table(dataf_peakannoCK2$annotation)

CK2_genes<-unique(dataf_peakannoCK2$SYMBOL)
CK2_peaks<-unique(dataf_peakannoCK2$peaksite)

#CK3 
dataf_peakannoCK3<-as.data.frame(peakAnnoCK3)
dataf_peakannoCK3$annotation<-gsub('Intron.*','Intron',dataf_peakannoCK3$annotation)
dataf_peakannoCK3$annotation<-gsub('Exon.*','Exon',dataf_peakannoCK3$annotation)
dataf_peakannoCK3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCK3$annotation)
dataf_peakannoCK3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCK3$annotation)
dataf_peakannoCK3$peaksite<-paste(dataf_peakannoCK3$SYMBOL,dataf_peakannoCK3$annotation,sep='_')
table(dataf_peakannoCK3$annotation)

#Remove distal intergenic
#dataf_peakannoCK3<-subset(dataf_peakannoCK3,dataf_peakannoCK3$annotation!="Distal Intergenic")
#table(dataf_peakannoCK3$annotation)

CK3_genes<-unique(dataf_peakannoCK3$SYMBOL)
CK3_peaks<-unique(dataf_peakannoCK3$peaksite)


# KO Common peaks genes in at least two samples 
CKcomgenes<-as.data.frame(table(sort(c(CK1_genes,CK2_genes,CK3_genes))))
CKcomgenes<-CKcomgenes[order(CKcomgenes$Freq,decreasing = TRUE),]
CKcomgenes_all<-subset(CKcomgenes,CKcomgenes$Freq>=2)
length(unique(CKcomgenes_all$Var1))

#Merge 2, 3 ,3 dataframes
colnames(dataf_peakannoCK1)[6]<-c("Name")
colnames(dataf_peakannoCK1)[7]<-c("Score")
colnames(dataf_peakannoCK1)[9]<-c("FE")
colnames(dataf_peakannoCK1)[10]<-c("log10pval")
colnames(dataf_peakannoCK1)[11]<-c("log10qval")
colnames(dataf_peakannoCK1)[12]<-c("param")

colnames(dataf_peakannoCK2)[6]<-c("Name")
colnames(dataf_peakannoCK2)[7]<-c("Score")
colnames(dataf_peakannoCK2)[9]<-c("FE")
colnames(dataf_peakannoCK2)[10]<-c("log10pval")
colnames(dataf_peakannoCK2)[11]<-c("log10qval")
colnames(dataf_peakannoCK2)[12]<-c("param")

colnames(dataf_peakannoCK3)[6]<-c("Name")
colnames(dataf_peakannoCK3)[7]<-c("Score")
colnames(dataf_peakannoCK3)[9]<-c("FE")
colnames(dataf_peakannoCK3)[10]<-c("log10pval")
colnames(dataf_peakannoCK3)[11]<-c("log10qval")
colnames(dataf_peakannoCK3)[12]<-c("param")

dataf_peakanno_allCK<-rbind(dataf_peakannoCK1,dataf_peakannoCK2,dataf_peakannoCK3)

#Select common genes from all databases
dataf_peakanno_allCK<-dataf_peakanno_allCK[with(dataf_peakanno_allCK, dataf_peakanno_allCK$SYMBOL %in% CKcomgenes_all$Var1),]

#sort by FE
dataf_peakanno_allCK<-dataf_peakanno_allCK[order(dataf_peakanno_allCK$FE,decreasing = TRUE),]


#CKT1 KO + TNF
dataf_peakannoCKT1<-as.data.frame(peakAnnoCKT1)
dataf_peakannoCKT1$annotation<-gsub('Intron.*','Intron',dataf_peakannoCKT1$annotation)
dataf_peakannoCKT1$annotation<-gsub('Exon.*','Exon',dataf_peakannoCKT1$annotation)
dataf_peakannoCKT1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCKT1$annotation)
dataf_peakannoCKT1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCKT1$annotation)
dataf_peakannoCKT1$peaksite<-paste(dataf_peakannoCKT1$SYMBOL,dataf_peakannoCKT1$annotation,sep='_')
table(dataf_peakannoCKT1$annotation)

CKT1_genes<-unique(dataf_peakannoCKT1$SYMBOL)
CKT1_peaks<-unique(dataf_peakannoCKT1$peaksite)

#CKT2 KO + TNF
dataf_peakannoCKT2<-as.data.frame(peakAnnoCKT2)
dataf_peakannoCKT2$annotation<-gsub('Intron.*','Intron',dataf_peakannoCKT2$annotation)
dataf_peakannoCKT2$annotation<-gsub('Exon.*','Exon',dataf_peakannoCKT2$annotation)
dataf_peakannoCKT2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCKT2$annotation)
dataf_peakannoCKT2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCKT2$annotation)
dataf_peakannoCKT2$peaksite<-paste(dataf_peakannoCKT2$SYMBOL,dataf_peakannoCKT2$annotation,sep='_')
table(dataf_peakannoCKT2$annotation)

CKT2_genes<-unique(dataf_peakannoCKT2$SYMBOL)
CKT2_peaks<-unique(dataf_peakannoCKT2$peaksite)


#CKIR KO + Irradiated
dataf_peakannoCKIR<-as.data.frame(peakAnnoCKIR)
dataf_peakannoCKIR$annotation<-gsub('Intron.*','Intron',dataf_peakannoCKIR$annotation)
dataf_peakannoCKIR$annotation<-gsub('Exon.*','Exon',dataf_peakannoCKIR$annotation)
dataf_peakannoCKIR$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCKIR$annotation)
dataf_peakannoCKIR$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCKIR$annotation)
dataf_peakannoCKIR$peaksite<-paste(dataf_peakannoCKIR$SYMBOL,dataf_peakannoCKIR$annotation,sep='_')
table(dataf_peakannoCKIR$annotation)

#Remove distal intergenic
#dataf_peakannoCKIR<-subset(dataf_peakannoCKIR,dataf_peakannoCKIR$annotation!="Distal Intergenic")
#table(dataf_peakannoCKIR$annotation)

CKIR_genes<-unique(dataf_peakannoCKIR$SYMBOL)
CKIR_peaks<-unique(dataf_peakannoCKIR$peaksite)

#KOIR 2ND RUN SAMP 1
dataf_peakannoKO1IR<-as.data.frame(peakAnnoKO1IR)
dataf_peakannoKO1IR$annotation<-gsub('Intron.*','Intron',dataf_peakannoKO1IR$annotation)
dataf_peakannoKO1IR$annotation<-gsub('Exon.*','Exon',dataf_peakannoKO1IR$annotation)
dataf_peakannoKO1IR$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoKO1IR$annotation)
dataf_peakannoKO1IR$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoKO1IR$annotation)
dataf_peakannoKO1IR$peaksite<-paste(dataf_peakannoKO1IR$SYMBOL,dataf_peakannoKO1IR$annotation,sep='_')
table(dataf_peakannoKO1IR$annotation)

#Remove distal intergenic
#dataf_peakannoKO1IR<-subset(dataf_peakannoKO1IR,dataf_peakannoKO1IR$annotation!="Distal Intergenic")
#table(dataf_peakannoKO1IR$annotation)

KO1IR_genes<-unique(dataf_peakannoKO1IR$SYMBOL)
KO1IR_peaks<-unique(dataf_peakannoKO1IR$peaksite)

#KOIR 2ND RUN SAMP 2
dataf_peakannoKO2IR<-as.data.frame(peakAnnoKO2IR)
dataf_peakannoKO2IR$annotation<-gsub('Intron.*','Intron',dataf_peakannoKO2IR$annotation)
dataf_peakannoKO2IR$annotation<-gsub('Exon.*','Exon',dataf_peakannoKO2IR$annotation)
dataf_peakannoKO2IR$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoKO2IR$annotation)
dataf_peakannoKO2IR$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoKO2IR$annotation)
dataf_peakannoKO2IR$peaksite<-paste(dataf_peakannoKO2IR$SYMBOL,dataf_peakannoKO2IR$annotation,sep='_')
table(dataf_peakannoKO2IR$annotation)

#Remove distal intergenic
#dataf_peakannoKO2IR<-subset(dataf_peakannoKO2IR,dataf_peakannoKO1IR$annotation!="Distal Intergenic")
#table(dataf_peakannoKO2IR$annotation)

KO2IR_genes<-unique(dataf_peakannoKO2IR$SYMBOL)
KO2IR_peaks<-unique(dataf_peakannoKO2IR$peaksite)


#KOIR 2ND RUN merged
dataf_peakannoKOIR<-as.data.frame(peakAnnoKOIR)
dataf_peakannoKOIR$annotation<-gsub('Intron.*','Intron',dataf_peakannoKOIR$annotation)
dataf_peakannoKOIR$annotation<-gsub('Exon.*','Exon',dataf_peakannoKOIR$annotation)
dataf_peakannoKOIR$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoKOIR$annotation)
dataf_peakannoKOIR$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoKOIR$annotation)
dataf_peakannoKOIR$peaksite<-paste(dataf_peakannoKOIR$SYMBOL,dataf_peakannoKOIR$annotation,sep='_')
table(dataf_peakannoKOIR$annotation)

#Remove distal intergenic
#dataf_peakannoKOIR<-subset(dataf_peakannoKOIR,dataf_peakannoKOIR$annotation!="Distal Intergenic")
#table(dataf_peakannoKOIR$annotation)

KOIR_genes<-unique(dataf_peakannoKOIR$SYMBOL)
KOIR_peaks<-unique(dataf_peakannoKOIR$peaksite)


length(CK1_genes)
length(CK2_genes)
length(CK3_genes)

#length(CKT1_genes)
#length(CKT2_genes)

length(CKIR_genes)
length(KO1IR_genes)
length(KO2IR_genes)
length(KOIR_genes)

KOIR_merge_genes<-c(CKIR_genes,KO1IR_genes,KO2IR_genes,KOIR_genes)
KOIR_merge_genes<-unique(KOIR_merge_genes)
length(KOIR_merge_genes)

# KO Common peaks genes in at least two samples 
KOIRcomgenes<-as.data.frame(table(sort(c(CKIR_genes,KO1IR_genes,KO2IR_genes,KOIR_genes))))
KOIRcomgenes<-KOIRcomgenes[order(KOIRcomgenes$Freq,decreasing = TRUE),]
#Select minimun 1 replicate for KOIR
KOIRcomgenes_all<-subset(KOIRcomgenes,KOIRcomgenes$Freq>=1)
length(unique(KOIRcomgenes_all$Var1))

#Merge 2, 3 ,3 dataframes
colnames(dataf_peakannoCKIR)[6]<-c("Name")
colnames(dataf_peakannoCKIR)[7]<-c("Score")
colnames(dataf_peakannoCKIR)[9]<-c("FE")
colnames(dataf_peakannoCKIR)[10]<-c("log10pval")
colnames(dataf_peakannoCKIR)[11]<-c("log10qval")
colnames(dataf_peakannoCKIR)[12]<-c("param")

colnames(dataf_peakannoKO1IR)[6]<-c("Name")
colnames(dataf_peakannoKO1IR)[7]<-c("Score")
colnames(dataf_peakannoKO1IR)[9]<-c("FE")
colnames(dataf_peakannoKO1IR)[10]<-c("log10pval")
colnames(dataf_peakannoKO1IR)[11]<-c("log10qval")
colnames(dataf_peakannoKO1IR)[12]<-c("param")

colnames(dataf_peakannoKO2IR)[6]<-c("Name")
colnames(dataf_peakannoKO2IR)[7]<-c("Score")
colnames(dataf_peakannoKO2IR)[9]<-c("FE")
colnames(dataf_peakannoKO2IR)[10]<-c("log10pval")
colnames(dataf_peakannoKO2IR)[11]<-c("log10qval")
colnames(dataf_peakannoKO2IR)[12]<-c("param")

colnames(dataf_peakannoKOIR)[6]<-c("Name")
colnames(dataf_peakannoKOIR)[7]<-c("Score")
colnames(dataf_peakannoKOIR)[9]<-c("FE")
colnames(dataf_peakannoKOIR)[10]<-c("log10pval")
colnames(dataf_peakannoKOIR)[11]<-c("log10qval")
colnames(dataf_peakannoKOIR)[12]<-c("param")


dataf_peakanno_allKOIR<-rbind(dataf_peakannoCKIR,dataf_peakannoKO1IR,dataf_peakannoKO2IR,dataf_peakannoKOIR)

#Select common genes from all databases
dataf_peakanno_allKOIR<-dataf_peakanno_allKOIR[with(dataf_peakanno_allKOIR, dataf_peakanno_allKOIR$SYMBOL %in% KOIRcomgenes_all$Var1),]

#sort by FE
dataf_peakanno_allKOIR<-dataf_peakanno_allKOIR[order(dataf_peakanno_allKOIR$FE,decreasing = TRUE),]






##DATAFRAME WITH genes id linked to peaks

#WT 1
dataf_peakannoCW1<-as.data.frame(peakAnnoCW1)
dataf_peakannoCW1$annotation<-gsub('Intron.*','Intron',dataf_peakannoCW1$annotation)
dataf_peakannoCW1$annotation<-gsub('Exon.*','Exon',dataf_peakannoCW1$annotation)
dataf_peakannoCW1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCW1$annotation)
dataf_peakannoCW1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCW1$annotation)
dataf_peakannoCW1$peaksite<-paste(dataf_peakannoCW1$SYMBOL,dataf_peakannoCW1$annotation,sep='_')
table(dataf_peakannoCW1$annotation)


#Remove distal intergenic
#dataf_peakannoCW1<-subset(dataf_peakannoCW1,dataf_peakannoCW1$annotation!="Distal Intergenic")
#table(dataf_peakannoCW1$annotation)

CW1_genes<-unique(dataf_peakannoCW1$SYMBOL)
CW1_peaks<-unique(dataf_peakannoCW1$peaksite)

#WT 2
dataf_peakannoCW2<-as.data.frame(peakAnnoCW2)
dataf_peakannoCW2$annotation<-gsub('Intron.*','Intron',dataf_peakannoCW2$annotation)
dataf_peakannoCW2$annotation<-gsub('Exon.*','Exon',dataf_peakannoCW2$annotation)
dataf_peakannoCW2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCW2$annotation)
dataf_peakannoCW2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCW2$annotation)
dataf_peakannoCW2$peaksite<-paste(dataf_peakannoCW2$SYMBOL,dataf_peakannoCW2$annotation,sep='_')
table(dataf_peakannoCW2$annotation)

#Remove distal intergenic
#dataf_peakannoCW2<-subset(dataf_peakannoCW2,dataf_peakannoCW2$annotation!="Distal Intergenic")
#table(dataf_peakannoCW2$annotation)

CW2_genes<-unique(dataf_peakannoCW2$SYMBOL)
CW2_peaks<-unique(dataf_peakannoCW2$peaksite)


#WT 3
dataf_peakannoCW3<-as.data.frame(peakAnnoCW3)
dataf_peakannoCW3$annotation<-gsub('Intron.*','Intron',dataf_peakannoCW3$annotation)
dataf_peakannoCW3$annotation<-gsub('Exon.*','Exon',dataf_peakannoCW3$annotation)
dataf_peakannoCW3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCW3$annotation)
dataf_peakannoCW3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCW3$annotation)
dataf_peakannoCW3$peaksite<-paste(dataf_peakannoCW3$SYMBOL,dataf_peakannoCW3$annotation,sep='_')
table(dataf_peakannoCW3$annotation)

#Remove distal intergenic
#dataf_peakannoCW3<-subset(dataf_peakannoCW3,dataf_peakannoCW3$annotation!="Distal Intergenic")
#table(dataf_peakannoCW3$annotation)

CW3_genes<-unique(dataf_peakannoCW3$SYMBOL)
CW3_peaks<-unique(dataf_peakannoCW3$peaksite)


# WT Common peaks genes in at least two samples 
CWcomgenes<-as.data.frame(table(sort(c(CW1_genes,CW2_genes,CW3_genes))))
CWcomgenes<-CWcomgenes[order(CWcomgenes$Freq,decreasing = TRUE),]
CWcomgenes_all<-subset(CWcomgenes,CWcomgenes$Freq>=2)
length(unique(CWcomgenes_all$Var1))


#Merge 2, 3 ,3 dataframes
colnames(dataf_peakannoCW1)[6]<-c("Name")
colnames(dataf_peakannoCW1)[7]<-c("Score")
colnames(dataf_peakannoCW1)[9]<-c("FE")
colnames(dataf_peakannoCW1)[10]<-c("log10pval")
colnames(dataf_peakannoCW1)[11]<-c("log10qval")
colnames(dataf_peakannoCW1)[12]<-c("param")

colnames(dataf_peakannoCW2)[6]<-c("Name")
colnames(dataf_peakannoCW2)[7]<-c("Score")
colnames(dataf_peakannoCW2)[9]<-c("FE")
colnames(dataf_peakannoCW2)[10]<-c("log10pval")
colnames(dataf_peakannoCW2)[11]<-c("log10qval")
colnames(dataf_peakannoCW2)[12]<-c("param")

colnames(dataf_peakannoCW3)[6]<-c("Name")
colnames(dataf_peakannoCW3)[7]<-c("Score")
colnames(dataf_peakannoCW3)[9]<-c("FE")
colnames(dataf_peakannoCW3)[10]<-c("log10pval")
colnames(dataf_peakannoCW3)[11]<-c("log10qval")
colnames(dataf_peakannoCW3)[12]<-c("param")

dataf_peakanno_allCW<-rbind(dataf_peakannoCW1,dataf_peakannoCW2,dataf_peakannoCW3)
dataf_peakanno_allCW<-dataf_peakanno_allCW[with(dataf_peakanno_allCW, dataf_peakanno_allCW$SYMBOL %in% CWcomgenes_all$Var1),]

#sort by FE
dataf_peakanno_allCW<-dataf_peakanno_allCW[order(dataf_peakanno_allCW$FE,decreasing = TRUE),]

CW_merge_genes<-c(CW1_genes,CW2_genes,CW3_genes)
CW_merge_genes<-unique(CW_merge_genes)
length(CW_merge_genes)



# WT + TNF Rep 1
dataf_peakannoCWT1<-as.data.frame(peakAnnoCWT1)
dataf_peakannoCWT1$annotation<-gsub('Intron.*','Intron',dataf_peakannoCWT1$annotation)
dataf_peakannoCWT1$annotation<-gsub('Exon.*','Exon',dataf_peakannoCWT1$annotation)
dataf_peakannoCWT1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCWT1$annotation)
dataf_peakannoCWT1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCWT1$annotation)
dataf_peakannoCWT1$peaksite<-paste(dataf_peakannoCWT1$SYMBOL,dataf_peakannoCWT1$annotation,sep='_')
table(dataf_peakannoCWT1$annotation)

CWT1_genes<-unique(dataf_peakannoCWT1$SYMBOL)
CWT1_peaks<-unique(dataf_peakannoCWT1$peaksite)

# WT + TNF Rep 2
dataf_peakannoCWT2<-as.data.frame(peakAnnoCWT2)
dataf_peakannoCWT2$annotation<-gsub('Intron.*','Intron',dataf_peakannoCWT2$annotation)
dataf_peakannoCWT2$annotation<-gsub('Exon.*','Exon',dataf_peakannoCWT2$annotation)
dataf_peakannoCWT2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCWT2$annotation)
dataf_peakannoCWT2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCWT2$annotation)
dataf_peakannoCWT2$peaksite<-paste(dataf_peakannoCWT2$SYMBOL,dataf_peakannoCWT2$annotation,sep='_')
table(dataf_peakannoCWT2$annotation)

CWT2_genes<-unique(dataf_peakannoCWT2$SYMBOL)
CWT2_peaks<-unique(dataf_peakannoCWT2$peaksite)



#WT + IRRADIATED
dataf_peakannoCWIR<-as.data.frame(peakAnnoCWIR)
dataf_peakannoCWIR$annotation<-gsub('Intron.*','Intron',dataf_peakannoCWIR$annotation)
dataf_peakannoCWIR$annotation<-gsub('Exon.*','Exon',dataf_peakannoCWIR$annotation)
dataf_peakannoCWIR$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoCWIR$annotation)
dataf_peakannoCWIR$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoCWIR$annotation)
dataf_peakannoCWIR$peaksite<-paste(dataf_peakannoCWIR$SYMBOL,dataf_peakannoCWIR$annotation,sep='_')
table(dataf_peakannoCWIR$annotation)

#Remove distal intergenic
#dataf_peakannoCWIR<-subset(dataf_peakannoCWIR,dataf_peakannoCWIR$annotation!="Distal Intergenic")
#table(dataf_peakannoCWIR$annotation)

#sort by FE
dataf_peakannoCWIR<-dataf_peakannoCWIR[order(dataf_peakannoCWIR$X9.02985,decreasing = TRUE),]


CWIR_genes<-unique(dataf_peakannoCWIR$SYMBOL)
CWIR_peaks<-unique(dataf_peakannoCWIR$peaksite)

length(CWIR_genes)

length(CW1_genes)
length(CW2_genes)
length(CW3_genes)


#WT irradiated 2n run SAMP 1
dataf_peakannoWT1IR<-as.data.frame(peakAnnoWT1IR)
dataf_peakannoWT1IR$annotation<-gsub('Intron.*','Intron',dataf_peakannoWT1IR$annotation)
dataf_peakannoWT1IR$annotation<-gsub('Exon.*','Exon',dataf_peakannoWT1IR$annotation)
dataf_peakannoWT1IR$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoWT1IR$annotation)
dataf_peakannoWT1IR$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoWT1IR$annotation)
dataf_peakannoWT1IR$peaksite<-paste(dataf_peakannoWT1IR$SYMBOL,dataf_peakannoWT1IR$annotation,sep='_')
table(dataf_peakannoWT1IR$annotation)

#Remove distal intergenic
#dataf_peakannoWT1IR<-subset(dataf_peakannoWT1IR,dataf_peakannoWT1IR$annotation!="Distal Intergenic")
#table(dataf_peakannoWT1IR$annotation)

#sort by FE
dataf_peakannoWT1IR<-dataf_peakannoWT1IR[order(dataf_peakannoWT1IR$X9.31548,decreasing = TRUE),]


WT1IR_genes<-unique(dataf_peakannoWT1IR$SYMBOL)
WT1IR_peaks<-unique(dataf_peakannoWT1IR$peaksite)


#WT irradiated 2n run SAMP 2
dataf_peakannoWT2IR<-as.data.frame(peakAnnoWT2IR)
dataf_peakannoWT2IR$annotation<-gsub('Intron.*','Intron',dataf_peakannoWT2IR$annotation)
dataf_peakannoWT2IR$annotation<-gsub('Exon.*','Exon',dataf_peakannoWT2IR$annotation)
dataf_peakannoWT2IR$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoWT2IR$annotation)
dataf_peakannoWT2IR$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoWT2IR$annotation)
dataf_peakannoWT2IR$peaksite<-paste(dataf_peakannoWT2IR$SYMBOL,dataf_peakannoWT2IR$annotation,sep='_')
table(dataf_peakannoWT2IR$annotation)

#Remove distal intergenic
#dataf_peakannoWT2IR<-subset(dataf_peakannoWT2IR,dataf_peakannoWT2IR$annotation!="Distal Intergenic")
#table(dataf_peakannoWT2IR$annotation)

#sort by FE
dataf_peakannoWT2IR<-dataf_peakannoWT2IR[order(dataf_peakannoWT2IR$X13.34014,decreasing = TRUE),]


WT2IR_genes<-unique(dataf_peakannoWT2IR$SYMBOL)
WT2IR_peaks<-unique(dataf_peakannoWT2IR$peaksite)


#WT irradiated merged 2nd run
dataf_peakannoWTIR<-as.data.frame(peakAnnoWTIR)
dataf_peakannoWTIR$annotation<-gsub('Intron.*','Intron',dataf_peakannoWTIR$annotation)
dataf_peakannoWTIR$annotation<-gsub('Exon.*','Exon',dataf_peakannoWTIR$annotation)
dataf_peakannoWTIR$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoWTIR$annotation)
dataf_peakannoWTIR$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoWTIR$annotation)
dataf_peakannoWTIR$peaksite<-paste(dataf_peakannoWTIR$SYMBOL,dataf_peakannoWTIR$annotation,sep='_')
table(dataf_peakannoWTIR$annotation)

#Remove distal intergenic
#dataf_peakannoWTIR<-subset(dataf_peakannoWTIR,dataf_peakannoWTIR$annotation!="Distal Intergenic")
#table(dataf_peakannoWTIR$annotation)

#sort by FE
dataf_peakannoWTIR<-dataf_peakannoWT2IR[order(dataf_peakannoWT2IR$X13.34014,decreasing = TRUE),]


WTIR_genes<-unique(dataf_peakannoWTIR$SYMBOL)
WTIR_peaks<-unique(dataf_peakannoWTIR$peaksite)


# WT Common peaks genes in at least two samples 
WTIRcomgenes<-as.data.frame(table(sort(c(CWIR_genes,WT1IR_genes,WT2IR_genes,WTIR_genes))))
WTIRcomgenes<-WTIRcomgenes[order(WTIRcomgenes$Freq,decreasing = TRUE),]
#select minimum 1 replicate
WTIRcomgenes_all<-subset(WTIRcomgenes,WTIRcomgenes$Freq>=2)
length(unique(WTIRcomgenes_all$Var1))


#Merge 2, 3 ,3 dataframes
colnames(dataf_peakannoCWIR)[6]<-c("Name")
colnames(dataf_peakannoCWIR)[7]<-c("Score")
colnames(dataf_peakannoCWIR)[9]<-c("FE")
colnames(dataf_peakannoCWIR)[10]<-c("log10pval")
colnames(dataf_peakannoCWIR)[11]<-c("log10qval")
colnames(dataf_peakannoCWIR)[12]<-c("param")

colnames(dataf_peakannoWT1IR)[6]<-c("Name")
colnames(dataf_peakannoWT1IR)[7]<-c("Score")
colnames(dataf_peakannoWT1IR)[9]<-c("FE")
colnames(dataf_peakannoWT1IR)[10]<-c("log10pval")
colnames(dataf_peakannoWT1IR)[11]<-c("log10qval")
colnames(dataf_peakannoWT1IR)[12]<-c("param")

colnames(dataf_peakannoWT2IR)[6]<-c("Name")
colnames(dataf_peakannoWT2IR)[7]<-c("Score")
colnames(dataf_peakannoWT2IR)[9]<-c("FE")
colnames(dataf_peakannoWT2IR)[10]<-c("log10pval")
colnames(dataf_peakannoWT2IR)[11]<-c("log10qval")
colnames(dataf_peakannoWT2IR)[12]<-c("param")

colnames(dataf_peakannoWTIR)[6]<-c("Name")
colnames(dataf_peakannoWTIR)[7]<-c("Score")
colnames(dataf_peakannoWTIR)[9]<-c("FE")
colnames(dataf_peakannoWTIR)[10]<-c("log10pval")
colnames(dataf_peakannoWTIR)[11]<-c("log10qval")
colnames(dataf_peakannoWTIR)[12]<-c("param")


dataf_peakanno_allWTIR<-rbind(dataf_peakannoCWIR,dataf_peakannoWT1IR,dataf_peakannoWT2IR,dataf_peakannoWTIR)
dataf_peakanno_allWTIR<-dataf_peakanno_allWTIR[with(dataf_peakanno_allWTIR, dataf_peakanno_allWTIR$SYMBOL %in% WTIRcomgenes_all$Var1),]

#sort by FE
dataf_peakanno_allWTIR<-dataf_peakanno_allWTIR[order(dataf_peakanno_allWTIR$FE,decreasing = TRUE),]

WTIR_merge_genes<-c(CWIR_genes,WT1IR_genes,WT2IR_genes,WTIR_genes)
WTIR_merge_genes<-unique(WTIR_merge_genes)
length(WTIR_merge_genes)







## HUMAN TUMOROIDS P
# WT1
dataf_peakannoPWT1<-as.data.frame(peakAnnoPWT1)
dataf_peakannoPWT1$annotation<-gsub('Intron.*','Intron',dataf_peakannoPWT1$annotation)
dataf_peakannoPWT1$annotation<-gsub('Exon.*','Exon',dataf_peakannoPWT1$annotation)
dataf_peakannoPWT1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoPWT1$annotation)
dataf_peakannoPWT1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoPWT1$annotation)
dataf_peakannoPWT1$peaksite<-paste(dataf_peakannoPWT1$SYMBOL,dataf_peakannoPWT1$annotation,sep='_')
table(dataf_peakannoPWT1$annotation)

#Remove distal intergenic
#dataf_peakannoPWT1<-subset(dataf_peakannoPWT1,dataf_peakannoPWT1$annotation!="Distal Intergenic")
#table(dataf_peakannoPWT1$annotation)

PWT1_genes<-unique(dataf_peakannoPWT1$SYMBOL)
PWT1_peaks<-unique(dataf_peakannoPWT1$peaksite)

length(PWT1_genes)


# WT2
dataf_peakannoPWT2<-as.data.frame(peakAnnoPWT2)
dataf_peakannoPWT2$annotation<-gsub('Intron.*','Intron',dataf_peakannoPWT2$annotation)
dataf_peakannoPWT2$annotation<-gsub('Exon.*','Exon',dataf_peakannoPWT2$annotation)
dataf_peakannoPWT2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoPWT2$annotation)
dataf_peakannoPWT2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoPWT2$annotation)
dataf_peakannoPWT2$peaksite<-paste(dataf_peakannoPWT2$SYMBOL,dataf_peakannoPWT2$annotation,sep='_')
table(dataf_peakannoPWT2$annotation)

#Remove distal intergenic
#dataf_peakannoPWT2<-subset(dataf_peakannoPWT2,dataf_peakannoPWT2$annotation!="Distal Intergenic")
#table(dataf_peakannoPWT2$annotation)


PWT2_genes<-unique(dataf_peakannoPWT2$SYMBOL)
PWT2_peaks<-unique(dataf_peakannoPWT2$peaksite)

length(PWT2_genes)


# KO1
dataf_peakannoPKO1<-as.data.frame(peakAnnoPKO1)
dataf_peakannoPKO1$annotation<-gsub('Intron.*','Intron',dataf_peakannoPKO1$annotation)
dataf_peakannoPKO1$annotation<-gsub('Exon.*','Exon',dataf_peakannoPKO1$annotation)
dataf_peakannoPKO1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoPKO1$annotation)
dataf_peakannoPKO1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoPKO1$annotation)
dataf_peakannoPKO1$peaksite<-paste(dataf_peakannoPKO1$SYMBOL,dataf_peakannoPKO1$annotation,sep='_')
table(dataf_peakannoPKO1$annotation)

#Remove distal intergenic
#dataf_peakannoPKO1<-subset(dataf_peakannoPKO1,dataf_peakannoPKO1$annotation!="Distal Intergenic")
#table(dataf_peakannoPKO1$annotation)

PKO1_genes<-unique(dataf_peakannoPKO1$SYMBOL)
PKO1_peaks<-unique(dataf_peakannoPKO1$peaksite)

length(PKO1_genes)

# KO2
dataf_peakannoPKO2<-as.data.frame(peakAnnoPKO2)
dataf_peakannoPKO2$annotation<-gsub('Intron.*','Intron',dataf_peakannoPKO2$annotation)
dataf_peakannoPKO2$annotation<-gsub('Exon.*','Exon',dataf_peakannoPKO2$annotation)
dataf_peakannoPKO2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoPKO2$annotation)
dataf_peakannoPKO2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoPKO2$annotation)
dataf_peakannoPKO2$peaksite<-paste(dataf_peakannoPKO2$SYMBOL,dataf_peakannoPKO2$annotation,sep='_')
table(dataf_peakannoPKO2$annotation)

#Remove distal intergenic
#dataf_peakannoPKO2<-subset(dataf_peakannoPKO2,dataf_peakannoPKO2$annotation!="Distal Intergenic")
#table(dataf_peakannoPKO2$annotation)

PKO2_genes<-unique(dataf_peakannoPKO2$SYMBOL)
PKO2_peaks<-unique(dataf_peakannoPKO2$peaksite)

length(PKO2_genes)


# WT Common peaks genes in at least two samples 
PWTcomgenes<-as.data.frame(table(sort(c(PWT1_genes,PWT2_genes))))
PWTcomgenes<-PWTcomgenes[order(PWTcomgenes$Freq,decreasing = TRUE),]
PWTcomgenes_all<-subset(PWTcomgenes,PWTcomgenes$Freq>=2)
length(unique(PWTcomgenes_all$Var1))

#Merge 1 AND 2 WT dataframes
colnames(dataf_peakannoPWT1)[6]<-c("Name")
colnames(dataf_peakannoPWT1)[7]<-c("Score")
colnames(dataf_peakannoPWT1)[9]<-c("FE")
colnames(dataf_peakannoPWT1)[10]<-c("log10pval")
colnames(dataf_peakannoPWT1)[11]<-c("log10qval")
colnames(dataf_peakannoPWT1)[12]<-c("param")

colnames(dataf_peakannoPWT2)[6]<-c("Name")
colnames(dataf_peakannoPWT2)[7]<-c("Score")
colnames(dataf_peakannoPWT2)[9]<-c("FE")
colnames(dataf_peakannoPWT2)[10]<-c("log10pval")
colnames(dataf_peakannoPWT2)[11]<-c("log10qval")
colnames(dataf_peakannoPWT2)[12]<-c("param")


dataf_peakanno_allPWT<-rbind(dataf_peakannoPWT1,dataf_peakannoPWT2)

dataf_peakanno_allPWT<-dataf_peakanno_allPWT[with(dataf_peakanno_allPWT, dataf_peakanno_allPWT$SYMBOL %in% PWTcomgenes_all$Var1),]



#KO Common peaks genes in at least two samples 
PKOcomgenes<-as.data.frame(table(sort(c(PKO1_genes,PKO2_genes))))
PKOcomgenes<-PKOcomgenes[order(PKOcomgenes$Freq,decreasing = TRUE),]
PKOcomgenes_all<-subset(PKOcomgenes,PKOcomgenes$Freq>=2)
length(unique(PKOcomgenes_all$Var1))

#Merge 1 AND 2 WT dataframes
colnames(dataf_peakannoPKO1)[6]<-c("Name")
colnames(dataf_peakannoPKO1)[7]<-c("Score")
colnames(dataf_peakannoPKO1)[9]<-c("FE")
colnames(dataf_peakannoPKO1)[10]<-c("log10pval")
colnames(dataf_peakannoPKO1)[11]<-c("log10qval")
colnames(dataf_peakannoPKO1)[12]<-c("param")

colnames(dataf_peakannoPKO2)[6]<-c("Name")
colnames(dataf_peakannoPKO2)[7]<-c("Score")
colnames(dataf_peakannoPKO2)[9]<-c("FE")
colnames(dataf_peakannoPKO2)[10]<-c("log10pval")
colnames(dataf_peakannoPKO2)[11]<-c("log10qval")
colnames(dataf_peakannoPKO2)[12]<-c("param")


dataf_peakanno_allPKO<-rbind(dataf_peakannoPKO1,dataf_peakannoPKO2)

dataf_peakanno_allPKO<-dataf_peakanno_allPKO[with(dataf_peakanno_allPKO, dataf_peakanno_allPKO$SYMBOL %in% PKOcomgenes_all$Var1),]







### Histone marks mm10

# WTH3
dataf_peakannoWTH3<-as.data.frame(peakAnnoWTH3)
dataf_peakannoWTH3$annotation<-gsub('Intron.*','Intron',dataf_peakannoWTH3$annotation)
dataf_peakannoWTH3$annotation<-gsub('Exon.*','Exon',dataf_peakannoWTH3$annotation)
dataf_peakannoWTH3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoWTH3$annotation)
dataf_peakannoWTH3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoWTH3$annotation)
dataf_peakannoWTH3$peaksite<-paste(dataf_peakannoWTH3$SYMBOL,dataf_peakannoWTH3$annotation,sep='_')
table(dataf_peakannoWTH3$annotation)

#Remove distal intergenic
#dataf_peakannoWTH3<-subset(dataf_peakannoWTH3,dataf_peakannoWTH3$annotation!="Distal Intergenic")
#table(dataf_peakannoWTH3$annotation)

#sort by FE
dataf_peakannoWTH3<-dataf_peakannoWTH3[order(dataf_peakannoWTH3$X31,decreasing = TRUE),]


WTH3_genes<-unique(dataf_peakannoWTH3$SYMBOL)
WTH3_peaks<-unique(dataf_peakannoWTH3$peaksite)



# WTH4
dataf_peakannoWTH4<-as.data.frame(peakAnnoWTH4)
dataf_peakannoWTH4$annotation<-gsub('Intron.*','Intron',dataf_peakannoWTH4$annotation)
dataf_peakannoWTH4$annotation<-gsub('Exon.*','Exon',dataf_peakannoWTH4$annotation)
dataf_peakannoWTH4$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoWTH4$annotation)
dataf_peakannoWTH4$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoWTH4$annotation)
dataf_peakannoWTH4$peaksite<-paste(dataf_peakannoWTH4$SYMBOL,dataf_peakannoWTH4$annotation,sep='_')
table(dataf_peakannoWTH4$annotation)

#Remove distal intergenic
#dataf_peakannoWTH3<-subset(dataf_peakannoWTH3,dataf_peakannoWTH3$annotation!="Distal Intergenic")
#table(dataf_peakannoWTH3$annotation)

#sort by FE
dataf_peakannoWTH4<-dataf_peakannoWTH4[order(dataf_peakannoWTH4$X32==TRUE),]


WTH4_genes<-unique(dataf_peakannoWTH4$SYMBOL)
WTH4_peaks<-unique(dataf_peakannoWTH4$peaksite)



colnames(dataf_peakannoWTH3)[6]<-c("Name")
colnames(dataf_peakannoWTH3)[7]<-c("Score")
colnames(dataf_peakannoWTH3)[9]<-c("FE")
colnames(dataf_peakannoWTH3)[10]<-c("log10pval")
colnames(dataf_peakannoWTH3)[11]<-c("log10qval")

colnames(dataf_peakannoWTH4)[6]<-c("Name")
colnames(dataf_peakannoWTH4)[7]<-c("Score")
colnames(dataf_peakannoWTH4)[9]<-c("FE")
colnames(dataf_peakannoWTH4)[10]<-c("log10pval")
colnames(dataf_peakannoWTH4)[11]<-c("log10qval")




# KO histone marks
dataf_peakannoKOH3<-as.data.frame(peakAnnoKOH3)
dataf_peakannoKOH3$annotation<-gsub('Intron.*','Intron',dataf_peakannoKOH3$annotation)
dataf_peakannoKOH3$annotation<-gsub('Exon.*','Exon',dataf_peakannoKOH3$annotation)
dataf_peakannoKOH3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoKOH3$annotation)
dataf_peakannoKOH3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoKOH3$annotation)
dataf_peakannoKOH3$peaksite<-paste(dataf_peakannoKOH3$SYMBOL,dataf_peakannoKOH3$annotation,sep='_')
table(dataf_peakannoKOH3$annotation)

#Remove distal intergenic
#dataf_peakannoKOH3<-subset(dataf_peakannoKOH3,dataf_peakannoKOH3$annotation!="Distal Intergenic")
#table(dataf_peakannoKOH3$annotation)

#sort by FE
dataf_peakannoKOH3<-dataf_peakannoKOH3[order(dataf_peakannoKOH3$X48,decreasing = TRUE),]


KOH3_genes<-unique(dataf_peakannoKOH3$SYMBOL)
KOH3_peaks<-unique(dataf_peakannoKOH3$peaksite)

# KOH4
dataf_peakannoKOH4<-as.data.frame(peakAnnoKOH4)
dataf_peakannoKOH4$annotation<-gsub('Intron.*','Intron',dataf_peakannoKOH4$annotation)
dataf_peakannoKOH4$annotation<-gsub('Exon.*','Exon',dataf_peakannoKOH4$annotation)
dataf_peakannoKOH4$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoKOH4$annotation)
dataf_peakannoKOH4$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoKOH4$annotation)
dataf_peakannoKOH4$peaksite<-paste(dataf_peakannoKOH4$SYMBOL,dataf_peakannoKOH4$annotation,sep='_')
table(dataf_peakannoKOH4$annotation)

#Remove distal intergenic
#dataf_peakannoKOH4<-subset(dataf_peakannoKOH4,dataf_peakannoKOH4$annotation!="Distal Intergenic")
#table(dataf_peakannoKOH4$annotation)

#sort by FE
dataf_peakannoKOH4<-dataf_peakannoKOH4[order(dataf_peakannoKOH4$X56,decreasing = TRUE),]


KOH4_genes<-unique(dataf_peakannoKOH4$SYMBOL)
KOH4_peaks<-unique(dataf_peakannoKOH4$peaksite)




colnames(dataf_peakannoKOH3)[6]<-c("Name")
colnames(dataf_peakannoKOH3)[7]<-c("Score")
colnames(dataf_peakannoKOH3)[9]<-c("FE")
colnames(dataf_peakannoKOH3)[10]<-c("log10pval")
colnames(dataf_peakannoKOH3)[11]<-c("log10qval")

colnames(dataf_peakannoKOH4)[6]<-c("Name")
colnames(dataf_peakannoKOH4)[7]<-c("Score")
colnames(dataf_peakannoKOH4)[9]<-c("FE")
colnames(dataf_peakannoKOH4)[10]<-c("log10pval")
colnames(dataf_peakannoKOH4)[11]<-c("log10qval")



## Peak distribution
#CKC1
distCK1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoCK1$annotation)))))
colnames(distCK1) <- as.character(unlist(distCK1[1,]))
distCK1<-distCK1[-1,]
distCK1$Sample<-"CK1"

#CKC2
distCK2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoCK2$annotation)))))
colnames(distCK2) <- as.character(unlist(distCK2[1,]))
distCK2<-distCK2[-1,]
distCK2$Sample<-"CK2"

#CKC3
distCK3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoCK3$annotation)))))
colnames(distCK3) <- as.character(unlist(distCK3[1,]))
distCK3<-distCK3[-1,]
distCK3$Sample<-"CK3"

#CWC1
distCW1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoCW1$annotation)))))
colnames(distCW1) <- as.character(unlist(distCW1[1,]))
distCW1<-distCW1[-1,]
distCW1$Sample<-"CW1"

#CWC2
distCW2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoCW2$annotation)))))
colnames(distCW2) <- as.character(unlist(distCW2[1,]))
distCW2<-distCW2[-1,]
distCW2$Sample<-"CW2"

#CWC3
distCW3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoCW3$annotation)))))
colnames(distCW3) <- as.character(unlist(distCW3[1,]))
distCW3<-distCW3[-1,]
distCW3$Sample<-"CW3"

#CWIR
distCWIR<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakanno_allWTIR$annotation)))))
colnames(distCWIR) <- as.character(unlist(distCWIR[1,]))
distCWIR<-distCWIR[-1,]
distCWIR$Sample<-"CWIR"

#WTI1
distWTI1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoWT1IR$annotation)))))
colnames(distWTI1) <- as.character(unlist(distWTI1[1,]))
distWTI1<-distWTI1[-1,]
distWTI1$Sample<-"WTI1R"

#WTI2
distWTI2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoWT2IR$annotation)))))
colnames(distWTI2) <- as.character(unlist(distWTI2[1,]))
distWTI2<-distWTI2[-1,]
distWTI2$Sample<-"WTI2R"


#WTI2
distWTI<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoWTIR$annotation)))))
colnames(distWTI) <- as.character(unlist(distWTI[1,]))
distWTI<-distWTI[-1,]
distWTI$Sample<-"WTIR"


#CKIR
distCKIR<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoCKIR$annotation)))))
colnames(distCKIR) <- as.character(unlist(distCKIR[1,]))
distCKIR<-distCKIR[-1,]
distCKIR$Sample<-"CKIR"

#KO1IR
distKO1IR<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoKO1IR$annotation)))))
colnames(distKO1IR) <- as.character(unlist(distKO1IR[1,]))
distKO1IR<-distKO1IR[-1,]
distKO1IR$Sample<-"KO1IR"


#KO2IR
distKO2IR<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoKO2IR$annotation)))))
colnames(distKO2IR) <- as.character(unlist(distKO2IR[1,]))
distKO2IR<-distKO2IR[-1,]
distKO2IR$Sample<-"KO2IR"

#KOIR
distKOIR<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoKOIR$annotation)))))
colnames(distKOIR) <- as.character(unlist(distKOIR[1,]))
distKOIR<-distKOIR[-1,]
distKOIR$Sample<-"KOIR"


#HISTONES MARKS
#WTH3
distWTH3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoWTH3$annotation)))))
colnames(distWTH3) <- as.character(unlist(distWTH3[1,]))
distWTH3<-distWTH3[-1,]
distWTH3$Sample<-"WTH3"

#WTH4
distWTH4<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoWTH4$annotation)))))
colnames(distWTH4) <- as.character(unlist(distWTH4[1,]))
distWTH4<-distWTH4[-1,]
distWTH4$Sample<-"WTH4"

#KOH3
distKOH3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoKOH3$annotation)))))
colnames(distKOH3) <- as.character(unlist(distKOH3[1,]))
distKOH3<-distKOH3[-1,]
distKOH3$Sample<-"KOH3"

#KOH4
distKOH4<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoKOH4$annotation)))))
colnames(distKOH4) <- as.character(unlist(distKOH4[1,]))
distKOH4<-distKOH4[-1,]
distKOH4$Sample<-"KOH4"


#tumoroids human
#PW1
distPWT1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoPWT1$annotation)))))
colnames(distPWT1) <- as.character(unlist(distPWT1[1,]))
distPWT1<-distPWT1[-1,]
distPWT1$Sample<-"PWT1"

#PW2
distPWT2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoPWT2$annotation)))))
colnames(distPWT2) <- as.character(unlist(distPWT2[1,]))
distPWT2<-distPWT2[-1,]
distPWT2$Sample<-"PWT2"


#PKO1
distPKO1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoPKO1$annotation)))))
colnames(distPKO1) <- as.character(unlist(distPKO1[1,]))
distPKO1<-distPKO1[-1,]
distPKO1$Sample<-"PKO1"


#PKO2
distPKO2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoPKO2$annotation)))))
colnames(distPKO2) <- as.character(unlist(distPKO2[1,]))
distPKO2<-distPKO2[-1,]
distPKO2$Sample<-"PKO2"


#iGUALAR COLUMNAS ENTRE DISTS!

distCWIR$`5' UTR`<-0
distbrd4<-rbind(distCK1,distCK2,
                distCK3,distCW1,distCW2,
                distCW3,distCWIR)

#Human tumoroids
distbrd4<-rbind(distPWT1,distPWT2,
                distPKO1,distPKO2)

#Histones
distbrd4<-rbind(distWTH3,distWTH4,
                distKOH3,distKOH4)


colnames(distbrd4)[1]<-c("Three_UTR")
colnames(distbrd4)[2]<-c("Five_UTR")
colnames(distbrd4)[3]<-c("Distal_intergenic")

#distbrd4$condition<-c("KO","KO","KO","WT","WT","WT","WTIR")
distbrd4$condition<-c("H3","H4","H3","H4")

distall_melt<-melt(distbrd4,id.vars = c("Sample","condition"))

distall_melt$value<-as.numeric(distall_melt$value)

ggplot(distall_melt,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat="identity",color="black") +
  facet_grid(~condition,scale="free",space = "free_x")+
  scale_fill_brewer(palette="Pastel1")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")



### Bed file gene annotation, selecting from bed genes with one peak in at least two replicates for each condition
bedgenes<-read.delim("/Volumes/cancer/db_files/mm10/mm10_genes.bed",header = FALSE)
colnames(bedgenes)<-c("Chr","Start","End","Name","Source","Strand")

#all genes
#for BRd4 peaks
selbed<-unique(c(as.character(CKcomgenes_all$Var1),as.character(CWcomgenes_all$Var1CKIR_genes),as.character(WTIRcomgenes_all$Var1)))
selbed<-bedgenes[with(bedgenes, bedgenes$Name %in% selbed),]

write.table(selbed,"/Volumes/cancer/CKC/peakcall/browser_files/genes_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")


#for histones H3
selbed<-unique(WTH3_genes,KOH3_genes)
selbed<-bedgenes[with(bedgenes, bedgenes$Name %in% selbed),]

write.table(selbed,"/Volumes/cancer/CKC/peakcall/browser_files/genes_peaks_histoneH3.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")


#for histones H4
selbed<-unique(WTH4_genes,KOH4_genes)
selbed<-bedgenes[with(bedgenes, bedgenes$Name %in% selbed),]

write.table(selbed,"/Volumes/cancer/CKC/peakcall/browser_files/genes_peaks_histoneH4.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

# for unique WTH3 and KOH3, select exclusive genes WTH3 and KOH3
setdiff(WTH3_genes,KOH3_genes)
setdiff(KOH3_genes,WTH3_genes)

#One bed per each group of genes

KOgenes<-bedgenes[with(bedgenes, bedgenes$Name %in% unique(CKcomgenes_all$Var1)),]
write.table(KOgenes,"/Volumes/cancer/CKC/peakcall/browser_files/KO_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

WTgenes<-bedgenes[with(bedgenes, bedgenes$Name %in% unique(CWcomgenes_all$Var1)),]
write.table(WTgenes,"/Volumes/cancer/CKC/peakcall/browser_files/WT_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

WTIRgenes<-bedgenes[with(bedgenes, bedgenes$Name %in% unique(WTIRcomgenes_all$Var1)),]
write.table(WTIRgenes,"/Volumes/cancer/CKC/peakcall/browser_files/WTIR_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

#peaks only in WTIR
#onlyWTIRgenes<-setdiff(WTIRcomgenes_all$Var1,CWcomgenes_all$Var1)

#onlyWTIRgenes<-bedgenes[with(bedgenes, bedgenes$Name %in% onlyWTIRgenes),]
#write.table(onlyWTIRgenes,"/Volumes/cancer/CKC/peakcall/browser_files/onlyWTIR_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

WTH3bedgenes<-bedgenes[with(bedgenes, bedgenes$Name %in% WTH3_genes),]
write.table(WTH3bedgenes,"/Volumes/cancer/CKC/peakcall/browser_files/WTH3_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

WTH4bedgenes<-bedgenes[with(bedgenes, bedgenes$Name %in% WTH4_genes),]
write.table(WTH4bedgenes,"/Volumes/cancer/CKC/peakcall/browser_files/WTH4_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

KOH3bedgenes<-bedgenes[with(bedgenes, bedgenes$Name %in% KOH3_genes),]
write.table(KOH3bedgenes,"/Volumes/cancer/CKC/peakcall/browser_files/KOH3_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

KOH4bedgenes<-bedgenes[with(bedgenes, bedgenes$Name %in% KOH4_genes),]
write.table(KOH4bedgenes,"/Volumes/cancer/CKC/peakcall/browser_files/KOH4_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")



## Cross RNASeq and peaks
#common peaks
rnaCK$external_gene_name
common_genes_peaks<-as.character(comgenes_all$Var1)
#unique genes CK
uniquegenes_CK_all
#unique genes CW
uniquegenes_CW_all


#Check peaks with brd4 targets paper (using human cell lines)
### venn diagram ##

# MICE WT genes vs KO genes peaks, minimum common 2 samples within each group
geneLists<-list(WT = CWcomgenes_all$Var1,KO = CKcomgenes_all$Var1,WTIR = WTIRcomgenes_all$Var1)
geneLists<-list(WT = CWcomgenes_all$Var1,KO = CKcomgenes_all$Var1)

# HUMAN WT genes vs KO in PD tumoroids
geneLists<-list(PWT = PWTcomgenes_all$Var1,PKO = PKOcomgenes_all$Var1)


## Histone WT vs KO
geneLists<-list(WTH3=WTH3_genes,KOH3=KOH3_genes)
geneLists<-list(WTH4=WTH4_genes,KOH4=KOH4_genes)


#geneLists<-list(WTH3=WTH3_genes,CWcomgenes_all$Var1)


geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off() 

# FOR THREE GROUPS
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue","pink"), alpha=c(0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("WT","KO","WTIR"), main="Genes IKK KO mice")

# FOR TWO GROUPS
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue"), alpha=c(0.3,0.3), cex = 2, cat.fontface=4, category.names=c("WTH3","KOH3"), main="Genes H3 WT KO")


# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

# To get the list of gene present in each Venn compartment we can use the gplots package
require("gplots")

a <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(a,"intersections")
length(inters[[1]])
length(inters[[2]])
length(inters[[3]])

length(inters[[4]])
length(inters[[5]])
length(inters[[6]])
length(inters[[7]])

# Gene lists venn
unique_WT<-inters[[3]]
unique_KO<-inters[[2]]
unique_WTIR<-inters[[4]]

common_all<-inters[[1]]
common_WTIR_KO<-inters[[6]]
common_WT_WTIR<-inters[[7]]
common_WT_KO<-inters[[5]]

write.table(unique_WT,'/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/unique_WT.txt',quote = FALSE,row.names = FALSE)
write.table(unique_KO,'/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/unique_KO.txt',quote = FALSE,row.names = FALSE)
write.table(unique_WTIR,'/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/unique_WTIR.txt',quote = FALSE,row.names = FALSE)
write.table(common_all,'/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/common_all.txt',quote = FALSE,row.names = FALSE)
write.table(common_WTIR_KO,'/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/common_WTIR_KO.txt',quote = FALSE,row.names = FALSE)
write.table(common_WT_WTIR,'/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/common_WT_WTIR',quote = FALSE,row.names = FALSE)
write.table(common_WT_KO,'/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/common_WT_KO',quote = FALSE,row.names = FALSE)


## Select genes in venn categories from datframe histones H3
KO_H3_unique<-inters[[1]]
uniquepeaks_KOH3<-dataf_peakannoKOH3[with(dataf_peakannoKOH3, dataf_peakannoKOH3$SYMBOL %in% KO_H3_unique),]
uniquepeaks_KOH3<-uniquepeaks_KOH3[order(uniquepeaks_KOH3$Score==TRUE),]


WT_H3_unique<-inters[[2]]
uniquepeaks_WTH3<-dataf_peakannoWTH3[with(dataf_peakannoWTH3, dataf_peakannoWTH3$SYMBOL %in% WT_H3_unique),]
uniquepeaks_WTH3<-uniquepeaks_WTH3[order(uniquepeaks_WTH3$Score==TRUE),]


## Select genes in venn categories from datframe histones H4
KO_H4_unique<-inters[[1]]
uniquepeaks_KOH4<-dataf_peakannoKOH4[with(dataf_peakannoKOH4, dataf_peakannoKOH4$SYMBOL %in% KO_H4_unique),]
uniquepeaks_KOH4<-uniquepeaks_KOH4[order(uniquepeaks_KOH4$Score==TRUE),]


WT_H4_unique<-inters[[2]]
uniquepeaks_WTH4<-dataf_peakannoWTH4[with(dataf_peakannoWTH4, dataf_peakannoWTH4$SYMBOL %in% WT_H4_unique),]
uniquepeaks_WTH4<-uniquepeaks_WTH4[order(uniquepeaks_WTH4$Score==TRUE),]






# PD human
commonPD<-inters[[3]]
KO_unique_PD<-inters[[1]]
WT_unique_PD<-inters[[2]]

common_hum<-data.frame(commonPD)
common_hum$class<-"commonPD"
colnames(common_hum)[1]<-"gene"

KO_unique_PD<-data.frame(KO_unique_PD)
KO_unique_PD$class<-"KO_unique_PD"
colnames(KO_unique_PD)[1]<-"gene"

WT_unique_PD<-data.frame(WT_unique_PD)
WT_unique_PD$class<-"WT_unique_PD"
colnames(WT_unique_PD)[1]<-"gene"

all_hum_genes_peaks<-rbind(common_hum,KO_unique_PD,WT_unique_PD)


# Mice
commonmice<-c(inters[[1]],inters[[5]])
KO_unique_mice<-inters[[2]]
WT_unique_mice<-inters[[3]]
WTIR_unique_mice<-inters[[4]]
common_WT_WTIR_mice<-inters[[7]]
common_KO_WTIR_mice<-inters[[6]]

micegenes<-data.frame(commonmice)
micegenes$class<-"commonmice"
colnames(micegenes)[1]<-"gene"

KO_unique_mice<-data.frame(KO_unique_mice)
KO_unique_mice$class<-"KO_unique_mice"
colnames(KO_unique_mice)[1]<-"gene"

WT_unique_mice<-data.frame(WT_unique_mice)
WT_unique_mice$class<-"WT_unique_mice"
colnames(WT_unique_mice)[1]<-"gene"

WTIR_unique_mice<-data.frame(WTIR_unique_mice)
WTIR_unique_mice$class<-"WTIR_unique_mice"
colnames(WTIR_unique_mice)[1]<-"gene"

common_WT_WTIR_mice<-data.frame(common_WT_WTIR_mice)
common_WT_WTIR_mice$class<-"common_WT_WTIR_mice"
colnames(common_WT_WTIR_mice)[1]<-"gene"

common_KO_WTIR_mice<-data.frame(common_KO_WTIR_mice)
common_KO_WTIR_mice$class<-"common_KO_WTIR_mice"
colnames(common_KO_WTIR_mice)[1]<-"gene"

all_mice_genes_peaks<-rbind(micegenes,KO_unique_mice,WT_unique_mice,WTIR_unique_mice,common_WT_WTIR_mice,common_KO_WTIR_mice)
all_mice_genes_peaks$gene<-toupper(all_mice_genes_peaks$gene)



## CROSS MICE AND human PD tumoroids
hum_mice_cross<-merge(all_hum_genes_peaks,all_mice_genes_peaks,by="gene")







# IRRADIATED ANALYSIS
common<-inters[[3]]
WTIR_unique<-inters[[4]]
WT_unique<-inters[[1]]

KOWTIR<-inters[[7]]


uniquepeaks_CWIR<-dataf_peakannoCWIR[with(dataf_peakannoCWIR, dataf_peakannoCWIR$SYMBOL %in% WTIR_unique),]
colnames(uniquepeaks_CWIR)[6]<-c("Name")
colnames(uniquepeaks_CWIR)[7]<-c("Score")
colnames(uniquepeaks_CWIR)[9]<-c("FE")
colnames(uniquepeaks_CWIR)[10]<-c("log10pval")
colnames(uniquepeaks_CWIR)[11]<-c("log10qval")


uniquepeaks_CWIR<-subset(uniquepeaks_CWIR,uniquepeaks_CWIR$log10qval>=2)
uniquepeaks_CWIR<-uniquepeaks_CWIR[order(uniquepeaks_CWIR$log10qval,decreasing = TRUE),]
length(unique(uniquepeaks_CWIR$SYMBOL))


write.table(uniquepeaks_CWIR,"/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/2ndrun/uniquepeaks_CWIR.txt",quote = FALSE,row.names = FALSE,sep="\t")




#Cross with RNA-Seq
rnaCK_common<-rnaCK[with(rnaCK, rnaCK$external_gene_name %in% common),]
rnaCK_common$peak<-c("common")
rnaCK_KO<-rnaCK[with(rnaCK, rnaCK$external_gene_name %in% KO_unique),]
rnaCK_KO$peak<-c("KO")
rnaCK_WT<-rnaCK[with(rnaCK, rnaCK$external_gene_name %in% WT_unique),]
rnaCK_WT$peak<-c("WT")

rnaCK_peaks<-rbind(rnaCK_common,rnaCK_KO,rnaCK_WT)
rnaCK_peaks$sig<-ifelse(rnaCK_peaks$padj<0.01,c("sig"),c("nosig"))
rnaCK_peaks$var<-paste(rnaCK_peaks$state,rnaCK_peaks$sig)

rnaCK_peaks$peak<-as.factor(rnaCK_peaks$peak)
rnaCK_peaks$peak <- factor(x = rnaCK_peaks$peak, levels = c("WT","common","KO"))

levels(rnaCK_peaks$peak)

ggplot(rnaCK_peaks,aes(x=peak,y=log2FoldChange_KO_vs_WT,label=external_gene_name))+
  geom_jitter(aes(color=sig))+
  geom_violin(data = subset(rnaCK_peaks,cross_rna_chips_ISC$sig=="sig"),alpha=0.2)+
  #geom_label_repel(data = subset(cross_rna_chips_ISC,cross_rna_chips_ISC$var=="ISC sig"),colour = "black", fontface = "italic",size=4)+
  #geom_violin(alpha=0.2)+
  scale_color_manual(values=c("grey","red"))+
  ylab("log2FC KO vs WT")+
  theme_bw()

write.table(rnaCK_peaks,'/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/cross_RNAseq_peaks.txt',quote = FALSE,row.names = FALSE,sep='\t')


genes_venns<-data.frame(matrix(NA_character_, nrow = 392, ncol = 4))
genes_venns<-data.frame(Genes_DEGs_WT,Genes_DEGs_KO,Genes_KO_peak,Genes_WT_peak)

write(Genes_DEGs_WT,file="/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/Genes_DEGs_WT.txt")
write(Genes_DEGs_KO,file="/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/Genes_DEGs_KO.txt")
write(Genes_KO_peak,file="/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/Genes_KO_peak.txt")
write(Genes_WT_peak,file="/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/Genes_WT_peak.txt")

#Compare Genes peaks KO, DEGs RNA_seq
geneLists<-list(DRNA=rnaCKgenesDOWN,URNA=rnaCKgenesUP,KO_peaks=Genes_KO_peak)

geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off() 
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue","pink"), alpha=c(0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("DownRNA","UpRNA","KO_peak"), main="Genes IKK KO mice")

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

# To get the list of gene present in each Venn compartment we can use the gplots package
require("gplots")
a <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(a,"intersections")
length(inters[[1]])
length(inters[[2]])
length(inters[[3]])
length(inters[[4]])
length(inters[[5]])

Up_peakKO<-inters[[2]]
Down_peakKO<-inters[[1]]
NoDEG_peakKO<-inters[[5]]

#Select genes from dataframe and distribution
Up_peakKO_df<-dataf_peakannoCK[with(dataf_peakannoCK, dataf_peakannoCK$SYMBOL %in% Up_peakKO),]
table(Up_peakKO_df$annotation)
Down_peakKO_df<-dataf_peakannoCK[with(dataf_peakannoCK, dataf_peakannoCK$SYMBOL %in% Down_peakKO),]
table(Down_peakKO_df$annotation)
NoDEG_peakKO_df<-dataf_peakannoCK[with(dataf_peakannoCK, dataf_peakannoCK$SYMBOL %in% NoDEG_peakKO),]
table(NoDEG_peakKO_df$annotation)







############ Functional Enrichment

## GO terms and pathways
library(devtools)

library(enrichR)

dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "GO_Molecular_Function_2018",
         "Chromosome_Location",
         "Epigenomics_Roadmap_HM_ChIP-seq",
         "ENCODE_TF_ChIP-seq_2018")

enriched_unique_WT <- enrichr(unique_WT, dbs)
enriched_unique_KO <- enrichr(unique_KO,dbs)
enriched_unique_WTIR <- enrichr(unique_WTIR,dbs)
enriched_common_all <- enrichr(common_all,dbs)
enriched_common_WTIR_KO<-enrichr(common_WTIR_KO,dbs)
enriched_common_WT_WTIR<-enrichr(common_WT_WTIR,dbs)
enriched_common_WT_KO<-enrichr(common_WT_KO,dbs)

#downrnako_peak<- enrichr(Genes_KODOWNRNA_downpeak,dbs)

#printEnrich(enrichedWT, "output_WT.csv" , sep = "\t", columns = c(1:9))
#printEnrich(enrichedK, "output_K.csv" , sep = "\t", columns = c(1:9))

unique_WTGO <- enriched_unique_WT[["GO_Biological_Process_2018"]]
unique_WTGO$cat<-"WTGO"
unique_KOGO <- enriched_unique_KO[["GO_Biological_Process_2018"]]
unique_KOGO$cat<-"KOGO"
unique_WTIRGO <- enriched_unique_WTIR[["GO_Biological_Process_2018"]]
unique_WTIRGO$cat<-"WTIRGO"
common_allGO <- enriched_common_all[["GO_Biological_Process_2018"]]
common_allGO$cat<-"common_allGO"
common_WTIR_KOGO <- enriched_common_WTIR_KO[["GO_Biological_Process_2018"]]
common_WTIR_KOGO$cat <- "WTIR_KOGO"
common_WT_WTIRGO <- enriched_common_WT_WTIR[["GO_Biological_Process_2018"]]
common_WT_WTIRGO$cat<-"WT_WTIRGO"
common_WT_KOGO <- enriched_common_WT_KO[["GO_Biological_Process_2018"]]
common_WT_KOGO$cat<-"WT_KOGO"


allGO<-rbind(unique_WTGO,unique_KOGO,unique_WTIRGO,common_allGO,common_WTIR_KOGO,common_WT_WTIRGO,common_WT_KOGO)
allGO$db<-c("GO_BP")

#allGO<-downrnako_peak[["GO_Biological_Process_2018"]]

#allGO<-separate(data = allGO, col = Overlap, into = c("counts", "pathway"), sep = "/")
#allGO<-subset(allGO,allGO$counts>=3)

bpsub<-subset(allGO,allGO$P.value<0.001)


bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)
class(bpsub$Combined.Score)
bpsub$cat<-as.factor(bpsub$cat)

bpsub$Genes<-gsub(';','/',bpsub$Genes)
bpsub$Term<-gsub('.GO.*','',bpsub$Term)

bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$cat, bpsub$P.value)), ]$Term)

bpsub$Term <- factor(bpsub$Term, levels=rev(unique(bpsub$Term)), ordered = TRUE)



#bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$P.value)), ]$Term)

ggplot(bpsub,aes(x=Term,y=-(log(P.value))))+
  geom_bar(position="dodge",stat = "identity",aes(fill=cat),color="black")+
  #geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Set1")+
  facet_grid(~cat,space = "free",scales="free")+
  theme(axis.text.y =element_text(size=8),axis.title=element_text(size=11,face="bold"),
        axis.text.x = element_text(size=8,angle=45),legend.position = "bottom")+
  coord_flip()


write.table(bpsub,"/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/bpsub_GO_BP_KOWT_WTIR_intergenic.csv",sep="\t",row.names = FALSE,quote = FALSE)
write.table(allGO,"/Volumes/grcmc/YGUILLEN/Irene_CKC_ChIP/gene_lists/bpsub_ALLnopval_BP_KOWT_WTIR_intergenic.csv",sep="\t",row.names = FALSE,quote = FALSE)


subgo<-bpsub[c(5,1,46,3,10,56,62,58,9,49,81),]

subgo<- subgo[order(subgo$P.value),]
subgo$Term<-as.factor(subgo$Term)
subgo$Term <- factor(subgo$Term, levels=subgo[rev(order(subgo$peak, subgo$P.value)), ]$Term)

ggplot(subgo,aes(x=Term,y=-(log(P.value)),fill=peak))+
  geom_bar(position="dodge",stat = "identity",color="black")+
  #geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Blues")+
  theme(axis.text.y =element_text(size=9),axis.title=element_text(size=14,face="bold"),axis.text.x = element_text(size=12))+
  coord_flip()


### Bed files for each category


unique_WT_bed<-bedgenes[with(bedgenes, bedgenes$Name %in% unique_WT),]
write.table(unique_WT_bed,"/Volumes/cancer/CKC/peakcall/browser_files/unique_WT_bed.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

unique_KO_bed<-bedgenes[with(bedgenes, bedgenes$Name %in% unique_KO),]
write.table(unique_KO_bed,"/Volumes/cancer/CKC/peakcall/browser_files/unique_KO_bed.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

unique_WTIR_bed<-bedgenes[with(bedgenes, bedgenes$Name %in% unique_WTIR),]
write.table(unique_WTIR_bed,"/Volumes/cancer/CKC/peakcall/browser_files/unique_WTIR_bed.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

common_all_bed<-bedgenes[with(bedgenes, bedgenes$Name %in% common_all),]
write.table(common_all_bed,"/Volumes/cancer/CKC/peakcall/browser_files/common_all_bed.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

common_WTIR_KO_bed<-bedgenes[with(bedgenes, bedgenes$Name %in% common_WTIR_KO),]
write.table(common_WTIR_KO_bed,"/Volumes/cancer/CKC/peakcall/browser_files/common_WTIR_KO_bed.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

common_WT_WTIR_bed<-bedgenes[with(bedgenes, bedgenes$Name %in% common_WT_WTIR),]
write.table(common_WT_WTIR_bed,"/Volumes/cancer/CKC/peakcall/browser_files/common_WT_WTIR_bed.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

common_WT_KO_bed<-bedgenes[with(bedgenes, bedgenes$Name %in% common_WT_KO),]
write.table(common_WT_KO_bed,"/Volumes/cancer/CKC/peakcall/browser_files/common_WT_KO_bed.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")


