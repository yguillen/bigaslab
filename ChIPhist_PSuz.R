
#Mouse genome
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
biocLite("EnsDb.Mmusculus.v79")
biocLite("org.Mm.eg.db")


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
library(ggrepel)

library(rJava)
library(venneuler)

setwd("/Volumes/grcmc/YGUILLEN/LMarrueChIP/")


## List of genes intestinal stem cell signature
scsignature<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/stem_cell_signature_list.txt")
fetalsignature_up<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/genes_lists/gene_signature_fetal_UP.txt")
fetalsignature_up$dir<-c("up")
fetalsignature_down<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/genes_lists/gene_signature_fetal_DOWN.txt")
fetalsignature_down$dir<-c("down")
fetalsignature<-rbind(fetalsignature_down,fetalsignature_up)

#NFKB gene targets
nfkb_targets<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/genes_lists/NFKB_targets.txt",header = FALSE)
nfkb_targets$V1<-paste0(toupper(substr(nfkb_targets$V1, 1, 1)), tolower(substr(nfkb_targets$V1, 2, nchar(nfkb_targets$V1))))
nfkb_targets<-nfkb_targets$V1

################## MACS2 PEAK CALLING with ChIPseeker ###############

## p65 WT
peakW1P6 <- readPeakFile("peakcalling_PSuzIKB_qval02/W1P6_peaks_peaks.narrowPeak")
peakW2P6 <- readPeakFile("peakcalling_PSuzIKB_qval02/W2P6_peaks_peaks.narrowPeak")

## p65 KO
peakK8P6 <- readPeakFile("peakcalling_PSuzIKB_qval02/K8P6_peaks_peaks.narrowPeak")
peakK9P6 <- readPeakFile("peakcalling_PSuzIKB_qval02/K9P6_peaks_peaks.narrowPeak")
peakKP6 <- readPeakFile("peakcalling_PSuzIKB_qval02/KP6_peaks_peaks.narrowPeak")

## Suz12 WT Broad peaks
peakW1SZ <- readPeakFile("peakcalling_PSuzIKB_qval02/W1SZ_peaks_peaks.broadPeak")
peakW2SZ <- readPeakFile("peakcalling_PSuzIKB_qval02/W2SZ_peaks_peaks.broadPeak")
peakWSZ <- readPeakFile("peakcalling_PSuzIKB_qval02/WSZ_peaks_peaks.broadPeak")

## Suz12 KO Broad peaks
peakK8SZ <- readPeakFile("peakcalling_PSuzIKB_qval02/K8SZ_peaks_peaks.broadPeak")
peakK9SZ <- readPeakFile("peakcalling_PSuzIKB_qval02/K9SZ_peaks_peaks.broadPeak")

## IKB WT Broad peaks
peakW4IKB <- readPeakFile("peakcalling_PSuzIKB_qval02/W4CI_peaks_peaks.broadPeak")
peakW5IKB <- readPeakFile("peakcalling_PSuzIKB_qval02/W5CI_peaks_peaks.broadPeak")

## IKB WT DSP Crosslink Broad peaks
peakW4DSP <- readPeakFile("peakcalling_PSuzIKB_qval02/W4DI_peaks_peaks.broadPeak")
peakW5DSP <- readPeakFile("peakcalling_PSuzIKB_qval02/W5DI_peaks_peaks.broadPeak")
peakWDSP <- readPeakFile("peakcalling_PSuzIKB_qval02/WDI_peaks_peaks.broadPeak")

# IKB WT DSP + noDSP
peakW_all<- readPeakFile("peakcalling_PSuzIKB_qval02/W_all_peaks_peaks.broadPeak")


#Overall all chromosmes peak coverage
covplot(peakW1P6,weightCol = 5)
covplot(peakW2P6, weightCol = 5)

covplot(peakK8P6, weightCol = 5)
covplot(peakK9P6, weightCol = 5)
covplot(peakKP6, weightCol = 5)

covplot(peakW1SZ, weightCol = 5)
covplot(peakW2SZ, weightCol = 5)
covplot(peakWSZ, weightCol = 5)

covplot(peakK8SZ, weightCol = 5)
covplot(peakK9SZ, weightCol = 5)

covplot(peakW4IKB, weightCol = 5)
covplot(peakW5IKB, weightCol = 5)

covplot(peakW4DSP, weightCol = 5)
covplot(peakW5DSP, weightCol = 5)
covplot(peakWDSP, weightCol = 5)

covplot(peakW_all, weightCol = 5)

#genes database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

### PEAK ANNOTATION WT P65
peakAnnoW1P6 <- annotatePeak(peakW1P6, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoW2P6 <- annotatePeak(peakW2P6, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

### PEAK ANNOTATION KO P65
peakAnnoK8P6 <- annotatePeak(peakK8P6, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoK9P6 <- annotatePeak(peakK9P6, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoKP6 <- annotatePeak(peakKP6, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

### PEAK ANNOTATION WT SZ
peakAnnoW1SZ <- annotatePeak(peakW1SZ, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoW2SZ <- annotatePeak(peakW2SZ, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWSZ <- annotatePeak(peakWSZ, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

### PEAK ANNOTATION KO SZ
peakAnnoK8SZ <- annotatePeak(peakK8SZ, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoK9SZ <- annotatePeak(peakK9SZ, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

### PEAK ANNOTATION WT IKB
peakAnnoW4IKB <- annotatePeak(peakW4IKB, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoW5IKB <- annotatePeak(peakW5IKB, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

### PEAK ANNOTATION WT DSP IKB
peakAnnoW4DSP <- annotatePeak(peakW4DSP, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoW5DSP <- annotatePeak(peakW5DSP, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWDSP <- annotatePeak(peakWDSP, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")


### PEAK ANNOTATION WT DSP IKB and no DSP
peakAnnoW_all <- annotatePeak(peakW_all, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")


# Visualize genomic annotation
plotAnnoPie(peakAnnoW5IKB)
dev.off()
plotAnnoBar(peakAnnoW5IKB,title = "Peak distribution WT IKB")
dev.off()
vennpie(peakAnnoW5IKB)
dev.off()

upsetplot(peakAnnoW5IKB, vennpie=TRUE)
dev.off()
upsetplot(peakAnnoW5IKB, vennpie=TRUE)
dev.off()

plotDistToTSS(peakAnnoW5IKB,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")



##DATAFRAME WITH genes id linked to peaks
dataf_peakAnnoW4IKB<-as.data.frame(peakAnnoW4IKB)
colnames(dataf_peakAnnoW4IKB)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoW4IKB$peaksite<-paste(dataf_peakAnnoW4IKB$SYMBOL,dataf_peakAnnoW4IKB$annotation,sep='_')
W4IKB_genes<-unique(dataf_peakAnnoW4IKB$SYMBOL)
W4IKB_peaks<-unique(dataf_peakAnnoW4IKB$peaksite)
sort(W4IKB_genes)


##DATAFRAME WITH genes id linked to peaks
dataf_peakAnnoW5IKB<-as.data.frame(peakAnnoW5IKB)
colnames(dataf_peakAnnoW5IKB)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoW5IKB$peaksite<-paste(dataf_peakAnnoW5IKB$SYMBOL,dataf_peakAnnoW5IKB$annotation,sep='_')
W5IKB_genes<-unique(dataf_peakAnnoW5IKB$SYMBOL)
W5IKB_peaks<-unique(dataf_peakAnnoW5IKB$peaksite)
sort(W5IKB_genes)


##DATAFRAME WITH genes id linked to peaks
dataf_peakAnnoW4DSP<-as.data.frame(peakAnnoW4DSP)
colnames(dataf_peakAnnoW4DSP)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoW4DSP$peaksite<-paste(dataf_peakAnnoW4DSP$SYMBOL,dataf_peakAnnoW4DSP$annotation,sep='_')
W4DSP_genes<-unique(dataf_peakAnnoW4DSP$SYMBOL)
W4DSP_peaks<-unique(dataf_peakAnnoW4DSP$peaksite)
sort(W4DSP_genes)

##DATAFRAME WITH genes id linked to peaks
dataf_peakAnnoW5DSP<-as.data.frame(peakAnnoW5DSP)
colnames(dataf_peakAnnoW5DSP)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoW5DSP$peaksite<-paste(dataf_peakAnnoW5DSP$SYMBOL,dataf_peakAnnoW5DSP$annotation,sep='_')
W5DSP_genes<-unique(dataf_peakAnnoW5DSP$SYMBOL)
W5DSP_peaks<-unique(dataf_peakAnnoW5DSP$peaksite)
sort(W5DSP_genes)

##DATAFRAME WITH genes id linked to peaks
dataf_peakAnnoWDSP<-as.data.frame(peakAnnoWDSP)
colnames(dataf_peakAnnoWDSP)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoWDSP$peaksite<-paste(dataf_peakAnnoWDSP$SYMBOL,dataf_peakAnnoWDSP$annotation,sep='_')
WDSP_genes<-unique(dataf_peakAnnoWDSP$SYMBOL)
WDSP_peaks<-unique(dataf_peakAnnoWDSP$peaksite)
sort(WDSP_genes)


IKBpeaks<-rbind(dataf_peakAnnoW4DSP,dataf_peakAnnoW5DSP,dataf_peakAnnoWDSP,dataf_peakAnnoW4IKB,dataf_peakAnnoW5IKB)

IKBpeaks<- IKBpeaks[rev(order(IKBpeaks$FE)),]
IKBpeaks$annotation<-gsub('Intron.*','Intron',IKBpeaks$annotation)
IKBpeaks$annotation<-gsub('Exon.*','Exon',IKBpeaks$annotation)
IKBpeaks$annotation<-gsub('Downstream.*','Downstream',IKBpeaks$annotation)
IKBpeaks$annotation<-gsub('Promoter.*','Promoter',IKBpeaks$annotation)
IKB_genes<-unique(IKBpeaks$SYMBOL)
IKB_peaks<-unique(IKBpeaks$peaksite)

distIKBpeaks<-as.data.frame(t(as.data.frame(prop.table(table(IKBpeaks$annotation)))))
colnames(distIKBpeaks) <- as.character(unlist(distIKBpeaks[1,]))
distIKBpeaks<-distIKBpeaks[-1,]
distIKBpeaks$Sample<-"IKB"

colnames(distIKBpeaks)[1]<-c("Three_UTR")
colnames(distIKBpeaks)[2]<-c("Five_UTR")

distall_IKBpeaks <- data.frame(sapply(distIKBpeaks[,1:6], function(x) as.numeric(as.character(x))*100))
colnames(distall_IKBpeaks)[1]<-"Distribution"
distall_IKBpeaks$condition<-"IKB"
distall_IKBpeaks$annotation<-row.names(distall_IKBpeaks)

ggplot(distall_IKBpeaks,aes(x=condition,y=Distribution,fill=annotation))+
  geom_bar(stat="identity") +
  #facet_wrap(~Condition,scale="free")+
  scale_fill_brewer(palette="Accent")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  xlab("")+ylab("Fraction of peaks")



write.table(IKBpeaks,"peakcalling_PSuzIKB_qval02/peaks_WIKB.txt",quote = FALSE,sep="\t")



##DATAFRAME P65 PEAKS WT
dataf_peakAnnoW1P6<-as.data.frame(peakAnnoW1P6)
colnames(dataf_peakAnnoW1P6)[c(6,7,9,10,11,12)]<-c("Name","Score","FE","log10pval","log10qval","summit")
dataf_peakAnnoW1P6$peaksite<-paste(dataf_peakAnnoW1P6$SYMBOL,dataf_peakAnnoW1P6$annotation,sep='_')
W1P6_genes<-unique(dataf_peakAnnoW1P6$SYMBOL)
W1P6_peaks<-unique(dataf_peakAnnoW1P6$peaksite)
sort(W1P6_genes)

dataf_peakAnnoW2P6<-as.data.frame(peakAnnoW2P6)
colnames(dataf_peakAnnoW2P6)[c(6,7,9,10,11,12)]<-c("Name","Score","FE","log10pval","log10qval","summit")
dataf_peakAnnoW2P6$peaksite<-paste(dataf_peakAnnoW2P6$SYMBOL,dataf_peakAnnoW2P6$annotation,sep='_')
W2P6_genes<-unique(dataf_peakAnnoW2P6$SYMBOL)
W2P6_peaks<-unique(dataf_peakAnnoW2P6$peaksite)
sort(W2P6_genes)

#dataframe p65 peaks KO

dataf_peakAnnoK8P6<-as.data.frame(peakAnnoK8P6)
colnames(dataf_peakAnnoK8P6)[c(6,7,9,10,11,12)]<-c("Name","Score","FE","log10pval","log10qval","summit")
dataf_peakAnnoK8P6$peaksite<-paste(dataf_peakAnnoK8P6$SYMBOL,dataf_peakAnnoK8P6$annotation,sep='_')
K8P6_genes<-unique(dataf_peakAnnoK8P6$SYMBOL)
K8P6_peaks<-unique(dataf_peakAnnoK8P6$peaksite)
sort(K8P6_genes)

dataf_peakAnnoK9P6<-as.data.frame(peakAnnoK9P6)
colnames(dataf_peakAnnoK9P6)[c(6,7,9,10,11,12)]<-c("Name","Score","FE","log10pval","log10qval","summit")
dataf_peakAnnoK9P6$peaksite<-paste(dataf_peakAnnoK9P6$SYMBOL,dataf_peakAnnoK9P6$annotation,sep='_')
K9P6_genes<-unique(dataf_peakAnnoK9P6$SYMBOL)
K9P6_peaks<-unique(dataf_peakAnnoK9P6$peaksite)
sort(K9P6_genes)

dataf_peakAnnoKP6<-as.data.frame(peakAnnoKP6)
colnames(dataf_peakAnnoKP6)[c(6,7,9,10,11,12)]<-c("Name","Score","FE","log10pval","log10qval","summit")
dataf_peakAnnoKP6$peaksite<-paste(dataf_peakAnnoKP6$SYMBOL,dataf_peakAnnoKP6$annotation,sep='_')
KP6_genes<-unique(dataf_peakAnnoKP6$SYMBOL)
KP6_peaks<-unique(dataf_peakAnnoKP6$peaksite)
sort(KP6_genes)


P65peaks_WT<-rbind(dataf_peakAnnoW1P6,dataf_peakAnnoW2P6)
P65peaks_WT<-P65peaks_WT[order(P65peaks_WT$Score,decreasing=TRUE),]

P65peaks_KO<-rbind(dataf_peakAnnoK8P6,dataf_peakAnnoK9P6,dataf_peakAnnoKP6)
P65peaks_KO<-P65peaks_KO[order(P65peaks_KO$FE),]

P65genes_WT<-unique(P65peaks_WT$SYMBOL)
P65genes_KO<-unique(P65peaks_KO$SYMBOL)

write.table(P65peaks,"peakcalling_PSuzIKB_qval02/peaks_P65.txt",quote = FALSE,sep="\t")



#Distribution and number of peaks
dataf_peakAnnoW1P6$annotation<-gsub('Exon.*','Exon',dataf_peakAnnoW1P6$annotation)
dataf_peakAnnoW1P6$annotation<-gsub('Intron.*','Intron',dataf_peakAnnoW1P6$annotation)
dataf_peakAnnoW1P6$annotation<-gsub('Downstream.*','Downstream',dataf_peakAnnoW1P6$annotation)
dataf_peakAnnoW1P6$annotation<-gsub('Promoter.*','Promoter',dataf_peakAnnoW1P6$annotation)
dataf_peakAnnoW1P6$peaksite<-paste(dataf_peakAnnoK8P6$SYMBOL,dataf_peakAnnoW1P6$annotation,sep='_')
W1P6table<-as.data.frame(table(dataf_peakAnnoW1P6$annotation))
W1P6table$sample<-"W1P6"
W1P6table$condition<-"WT"

dataf_peakAnnoW2P6$annotation<-gsub('Exon.*','Exon',dataf_peakAnnoW2P6$annotation)
dataf_peakAnnoW2P6$annotation<-gsub('Intron.*','Intron',dataf_peakAnnoW2P6$annotation)
dataf_peakAnnoW2P6$annotation<-gsub('Downstream.*','Downstream',dataf_peakAnnoW2P6$annotation)
dataf_peakAnnoW2P6$annotation<-gsub('Promoter.*','Promoter',dataf_peakAnnoW2P6$annotation)
dataf_peakAnnoW2P6$peaksite<-paste(dataf_peakAnnoW2P6$SYMBOL,dataf_peakAnnoW2P6$annotation,sep='_')
W2P6table<-as.data.frame(table(dataf_peakAnnoW2P6$annotation))
W2P6table$sample<-"W2P6"
W2P6table$condition<-"WT"

dataf_peakAnnoK8P6$annotation<-gsub('Exon.*','Exon',dataf_peakAnnoK8P6$annotation)
dataf_peakAnnoK8P6$annotation<-gsub('Intron.*','Intron',dataf_peakAnnoK8P6$annotation)
dataf_peakAnnoK8P6$annotation<-gsub('Downstream.*','Downstream',dataf_peakAnnoK8P6$annotation)
dataf_peakAnnoK8P6$annotation<-gsub('Promoter.*','Promoter',dataf_peakAnnoK8P6$annotation)
dataf_peakAnnoK8P6$peaksite<-paste(dataf_peakAnnoK8P6$SYMBOL,dataf_peakAnnoK8P6$annotation,sep='_')
K8P6table<-as.data.frame(table(dataf_peakAnnoK8P6$annotation))
K8P6table$sample<-"K8P6"
K8P6table$condition<-"KO"

dataf_peakAnnoK9P6$annotation<-gsub('Exon.*','Exon',dataf_peakAnnoK9P6$annotation)
dataf_peakAnnoK9P6$annotation<-gsub('Intron.*','Intron',dataf_peakAnnoK9P6$annotation)
dataf_peakAnnoK9P6$annotation<-gsub('Downstream.*','Downstream',dataf_peakAnnoK9P6$annotation)
dataf_peakAnnoK9P6$annotation<-gsub('Promoter.*','Promoter',dataf_peakAnnoK9P6$annotation)
dataf_peakAnnoK9P6$peaksite<-paste(dataf_peakAnnoK9P6$SYMBOL,dataf_peakAnnoK9P6$annotation,sep='_')
K9P6table<-as.data.frame(table(dataf_peakAnnoK9P6$annotation))
K9P6table$sample<-"K9P6"
K9P6table$condition<-"KO"

P65table<-rbind(W1P6table,W2P6table,K8P6table,K9P6table)
colnames(P65table)[1]<-c("Annotation")

P65table<-subset(P65table,P65table$sample!="W1P6")
P65table<-subset(P65table,P65table$sample!="K8P6")

ggplot(P65table,aes(x=Annotation,y=Freq,fill=Annotation))+
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
  facet_wrap(~sample,ncol=4)+
  theme_bw()+
  scale_fill_brewer()+
  theme(axis.text.x = element_blank())


##DATAFRAME SZ PEAKS WT
dataf_peakAnnoW1SZ<-as.data.frame(peakAnnoW1SZ)
colnames(dataf_peakAnnoW1SZ)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoW1SZ$peaksite<-paste(dataf_peakAnnoW1SZ$SYMBOL,dataf_peakAnnoW1SZ$annotation,sep='_')
W1SZ_genes<-unique(dataf_peakAnnoW1SZ$SYMBOL)
W1SZ_peaks<-unique(dataf_peakAnnoW1SZ$peaksite)
sort(W1SZ_genes)

dataf_peakAnnoW2SZ<-as.data.frame(peakAnnoW2SZ)
colnames(dataf_peakAnnoW2SZ)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoW2SZ$peaksite<-paste(dataf_peakAnnoW2SZ$SYMBOL,dataf_peakAnnoW2SZ$annotation,sep='_')
W2SZ_genes<-unique(dataf_peakAnnoW2SZ$SYMBOL)
W2SZ_peaks<-unique(dataf_peakAnnoW2SZ$peaksite)
sort(W2SZ_genes)

dataf_peakAnnoWSZ<-as.data.frame(peakAnnoWSZ)
colnames(dataf_peakAnnoWSZ)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoWSZ$peaksite<-paste(dataf_peakAnnoWSZ$SYMBOL,dataf_peakAnnoWSZ$annotation,sep='_')
WSZ_genes<-unique(dataf_peakAnnoWSZ$SYMBOL)
WSZ_peaks<-unique(dataf_peakAnnoWSZ$peaksite)
sort(WSZ_genes)

#dataframe SZ peaks KO

dataf_peakAnnoK8SZ<-as.data.frame(peakAnnoK8SZ)
colnames(dataf_peakAnnoK8SZ)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoK8SZ$peaksite<-paste(dataf_peakAnnoK8SZ$SYMBOL,dataf_peakAnnoK8SZ$annotation,sep='_')
K8SZ_genes<-unique(dataf_peakAnnoK8SZ$SYMBOL)
K8SZ_peaks<-unique(dataf_peakAnnoK8SZ$peaksite)
sort(K8SZ_genes)

dataf_peakAnnoK9SZ<-as.data.frame(peakAnnoK9SZ)
colnames(dataf_peakAnnoK9SZ)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakAnnoK9SZ$peaksite<-paste(dataf_peakAnnoK9SZ$SYMBOL,dataf_peakAnnoK9SZ$annotation,sep='_')
K9SZ_genes<-unique(dataf_peakAnnoK9SZ$SYMBOL)
K9SZ_peaks<-unique(dataf_peakAnnoK9SZ$peaksite)
sort(K9SZ_genes)

SZpeaks<-rbind(dataf_peakAnnoK8SZ,dataf_peakAnnoK9SZ,dataf_peakAnnoW1SZ,dataf_peakAnnoW2SZ,dataf_peakAnnoWSZ)

write.table(SZpeaks,"peakcalling_PSuzIKB_qval02/peaks_SZ.txt",quote = FALSE,sep="\t")


#genes
comgenes<-as.data.frame(table(sort(c(K8SZ_genes,WSZ_genes))))

#peaks
#comgenes<-as.data.frame(table(sort(c(WTH3_peaks,KH3_peaks))))

comgenes<-comgenes[order(comgenes$Freq,decreasing = TRUE),]
comgenes_all<-subset(comgenes,comgenes$Freq==2)
uniqueenes_all<-subset(comgenes,comgenes$Freq==1)

#genes
uniquegenes_KO<-setdiff(K8SZ_genes,comgenes_all$Var1)
uniquegenes_WT<-setdiff(WSZ_genes,comgenes_all$Var1)


#genes
uniquegenes_KO_sites<-dataf_peakAnnoK8SZ[with(dataf_peakAnnoK8SZ, dataf_peakAnnoK8SZ$SYMBOL %in% uniquegenes_KO),]
uniquegenes_WT_sites<-dataf_peakAnnoWSZ[with(dataf_peakAnnoWSZ, dataf_peakAnnoWSZ$SYMBOL %in% uniquegenes_WT),]


uniquegenes_KO_sites$type=c("KO")
uniquegenes_KO_sites_sub<-subset(uniquegenes_KO_sites,select=c(start,seqnames,width,type,Score,SYMBOL))
colnames(uniquegenes_KO_sites_sub)<-c("start","seqnames","width","type","cov","Gene")
uniquegenes_WT_sites$type=c("WT")
uniquegenes_WT_sites_sub<-subset(uniquegenes_WT_sites,select=c(start,seqnames,width,type,Score,SYMBOL))
colnames(uniquegenes_WT_sites_sub)<-c("start","seqnames","width","type","cov","Gene")

uniquegenes_KO_WT<-rbind(uniquegenes_KO_sites_sub,uniquegenes_WT_sites_sub)
uniquegenes_KO_WT$seqnames<-gsub('chr','',uniquegenes_KO_WT$seqnames)
uniquegenes_KO_WT$seqnames<-as.factor(uniquegenes_KO_WT$seqnames)

uniquegenes_KO_WT$seqnames <- ordered(uniquegenes_KO_WT$seqnames,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X"))

#merge with ISC signature
#scsignature$source<-"SC"

#merge with fetal signature
#fetalsignature$source<-"Fetal"

#uniquegenes_KH3_WTH3<-merge(scsignature,uniquegenes_KH3_WTH3,all.y="TRUE")
#uniquegenes_KH3_WTH3<-merge(fetalsignature,uniquegenes_KH3_WTH3,all.y="TRUE")

#uniquegenes_KH3_WTH3$source[is.na(uniquegenes_KH3_WTH3$source)] <- "noSC"
#uniquegenes_KH3_WTH3$dir[is.na(uniquegenes_KH3_WTH3$dir)] <- "NULL"

#uniquegenes_KH3_WTH3$class <- paste(uniquegenes_KH3_WTH3$dir,uniquegenes_KH3_WTH3$type,sep="_")

uniquegenes_KO_WT_noX<-subset(uniquegenes_KO_WT,uniquegenes_KO_WT$seqnames!="X")

ggplot(uniquegenes_KO_WT_noX,aes(x=start,y=cov,label=Gene))+
  geom_point(aes(color=type))+
  geom_linerange(aes(x=start, ymax=cov, ymin=0,color=type),position = position_jitter(height = 0L, seed = 1L))+
  #geom_text(aes(start, cov, label = Gene), data = subset(uniquegenes_KH3_WTH3_noX,uniquegenes_KH3_WTH3_noX$source=="SC"),check_overlap = TRUE,vjust=-0.5)+
  geom_label_repel(aes(fill = type),colour = "black", fontface = "italic",size=2.5)+
  #scale_color_manual(values=c("coral3","cornflowerblue","firebrick1","cyan"))+
  scale_color_manual(values=c("coral3","cornflowerblue"))+
  #scale_fill_manual(values=c("darksalmon","firebrick1","darkslateblue"))+
  facet_wrap(~seqnames,scales="free",ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_blank())



#venn diagram common peaks or genes KO and WT
### venn diagram ##
IKBgenes<-sort(unique(IKBpeaks$SYMBOL))

SZgenesWT<-unique(c(W1SZ_genes,W2SZ_genes,WSZ_genes))
SZgenesKO<-unique(c(K8SZ_genes,K9SZ_genes))

SZgenes<-sort(unique(SZpeaks$SYMBOL))



geneLists <- list(WT = P65genes_WT, KO = P65genes_KO)
geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off()
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue"), alpha=c(0.3,0.3), cex = 2, cat.fontface=4, category.names=c("WT","KO"), main="Genes")

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


## Cross RNASeq with IKB genes
IKBgenes
IKB_rnaseq<-res_df[with(res_df, res_df$SYMBOL %in% IKBgenes),]
IKB_rnaseq<-subset(IKB_rnaseq,select=c(log2FoldChange,padj,SYMBOL))
IKB_rnaseq$sig<-ifelse(IKB_rnaseq$padj<0.1,c("sig"),c("nosig"))

IKB_rnaseq$condition<-"IKB"

ggplot(IKB_rnaseq,aes(x=condition,y=log2FoldChange))+
  geom_jitter(aes(color=sig))+
  geom_violin(alpha=0.2)+
  scale_color_manual(values=c("grey","red"))+
  geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  geom_hline(yintercept = 0)+
  ylab("log2FC KO vs WT")+
  theme_bw()


# We can summarize the contents of each venn compartment, as follows:
# in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
lapply(inters, head) 


## Venn euler proportional diagram
vd <- venneuler(c(KO=9164, WT=8141, ISC=510, 
                  "KO&WT"=7623, "KO&ISC"=219, "WT&ISC"=182 ,"KO&WT&ISC"=170))

vd <- venneuler(c(KO=KH3_genes, WT=WTH3_genes, ISC=fetalsignature$Gene, 
                  "KO&WT", "KO&ISC", "WT&ISC","KO&WT&ISC"))
plot(vd)

vd$centers


### PATHWAYS ENRICHMENT ##
geneWTH3 <- seq2gene(peakWTH3, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathwayWTH3 <- enrichPathway(geneWTH3,organism = "mouse")
head(pathwayWTH3, 2)
dotplot(pathwayWTH3)

dev.off()

geneKH3 <- seq2gene(peakKH3, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathwayKH3 <- enrichPathway(geneKH3,organism = "mouse")
head(pathwayKH3, 2)
dotplot(pathwayKH3)


####### Compare peaks with peakAnno ###
path <- c("/Users/instalar/Desktop/LMarrueChIP/peakcalling/")

## Compare overlapping peaks WT and KO

WTH3peaks<-toGRanges("peakcalling/WH3_peaks_peaks.narrowPeak",format="narrowPeak",header=FALSE)
KH3peaks<-toGRanges("peakcalling/KH3_peaks_peaks.narrowPeak",format="narrowPeak",header=FALSE)

ol <- findOverlapsOfPeaks(WTH3peaks,KH3peaks,connectedPeaks = "keepAll")

## add metadata (mean of score) to the overlapping peaks
#ol <- addMetadata(ol, colNames="score", FUN=mean) 

# Venn diagram overlapping peaks
averagePeakWidth <- mean(width(unlist(GRangesList(ol$peaklist))))
tot <- ceiling(1.87e+09 * .03 / averagePeakWidth)
makeVennDiagram(ol, totalTest=tot, connectedPeaks="keepAll")
dev.off()

## Annotate
annoData <- toGRanges(EnsDb.Mmusculus.v79, feature="gene")
annoData[1:2]

# binding site distribution relative to features
overlaps <- ol$peaklist[["WTH3peaks///KH3peaks"]]
binOverFeature(overlaps, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")

aCR<-assignChromosomeRegion(overlaps, nucleotideLevel=FALSE, 
                            precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                            TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene)
barplot(aCR$percentage, las=3)

# Annotate overlapping peaks
overlaps.anno <- annotatePeakInBatch(overlaps, 
                                     AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))
overlaps.anno <- addGeneIDs(overlaps.anno,
                            "org.Mm.eg.db",
                            IDs2Add = "entrez_id")
head(overlaps.anno)


## GO terms and pathways
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018")
enrichedWT <- enrichr(uniquegenes_WTH3, dbs)
enrichedKH3 <- enrichr(uniquegenes_KH3,dbs)

printEnrich(enrichedWT, "output_WT.csv" , sep = "\t", columns = c(1:9))
printEnrich(enrichedKH3, "output_KH3.csv" , sep = "\t", columns = c(1:9))

WTGO <- enrichedWT[["GO_Biological_Process_2018"]]
WTGO$peak<-c("WT")

KOGO <- enrichedKH3[["GO_Biological_Process_2018"]]
KOGO$peak<-c("KO")

allGO<-rbind(WTGO,KOGO)
allGO$db<-c("GO")

allGO<-separate(data = allGO, col = Overlap, into = c("counts", "pathway"), sep = "/")
allGO<-subset(allGO,allGO$counts>=5)

bpsub<-subset(allGO,allGO$P.value<=0.05)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)
class(bpsub$Combined.Score)

bpsub$Genes<-gsub(';','/',bpsub$Genes)
bpsub$Term<-gsub('.GO.*','',bpsub$Term)

bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$peak, bpsub$P.value)), ]$Term)
ggplot(bpsub,aes(x=Term,y=-(log(P.value)),fill=peak))+
  geom_bar(position="dodge",stat = "identity",color="black")+
  geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Blues")+
  theme(axis.text.y =element_text(size=9),axis.title=element_text(size=14,face="bold"),axis.text.x = element_text(size=12))+
  coord_flip()

write.table(bpsub,"/Volumes/grcmc/YGUILLEN/LMarrueChIP/GO_BP_filter5genes_KOWT.csv",sep="\t",row.names = FALSE,quote = FALSE)


WTTF <- enrichedWT[["TF_Perturbations_Followed_by_Expression"]]
WTTF$peak<-c("WT")

KOTF <- enrichedKH3[["TF_Perturbations_Followed_by_Expression"]]
KOTF$peak<-c("KO")

allTF<-rbind(WTTF,KOTF)
allTF$db<-c("TF")


allGOTF<-rbind(allGO,allTF)
write.table(allGOTF,"GO2015TF.csv",row.names = FALSE,quote = FALSE,sep="\t")

bpsub<-subset(allGOTF,allGOTF$P.value<=0.05)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$P.value))])

ggplot(bpsub,aes(x=peak,y=Term,fill=P.value),alpha = 0.5)+
  geom_tile(color="black",width=1,height=0.8)+
  scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  facet_wrap(~db,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=9),axis.title=element_text(size=14,face="bold"))
  

dev.off()
