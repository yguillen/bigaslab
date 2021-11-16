#Trainning ChIPSeq data analyses

# Install packages
install.packages("data.table")
# AnnoPeak (peak annotations)
install.packages(c("shiny", "ggplot2", "VennDiagram", "RColorBrewer", "reshape2", "xtable", "gplots"))
install.packages("RMySQL")
install.packages("rJava")
install.packages("xlsx")

BiocManager::install("ChIPseeker")


source("https://bioconductor.org/biocLite.R")

biocLite("GenomicAlignments")
biocLite("ChIPQC")
biocLite("GenomicRanges")
biocLite("DiffBind")
biocLite("GO.db")
biocLite(c("ChIPpeakAnno", 
           "biomaRt",
           "GSEABase",
           "GOstats",
           "TxDb.Mmusculus.UCSC.mm9.knownGene",
           "TxDb.Mmusculus.UCSC.mm10.knownGene",
           "TxDb.Hsapiens.UCSC.hg18.knownGene",
           "TxDb.Hsapiens.UCSC.hg19.knownGene",
           "TxDb.Hsapiens.UCSC.hg38.knownGene",
           "org.Hs.eg.db",
           "org.Mm.eg.db"))
biocLite("BayesPeak")
biocLite("ChIPseeker")
biocLite("clusterProfiler")
biocLite("ReactomePA")
biocLite("trackViewer")
biocLite("Gviz")
biocLite("rtracklayer")
biocLite("Sushi")
devtools::install_github("jolars/eulerr")
biocLite("venneuler")
biocLite("DO.db")
biocLite("VennDiagram")
biocLite("readxl")
biocLite("ddplyr")
biocLite("geneplotter")
biocLite("MDplot")
biocLite("enrichR")
biocLite("viridis")
biocLite("vegan")
biocLite("ggrepel")
biocLite("scatterpie")

#Load libraries
library(data.table)
library(GenomicAlignments)
library(GO.db)
library(DO.db)
#library(DiffBind)
#library(ChIPQC)
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
#library(BayesPeak)
library(parallel)
library(ChIPseeker)
#library(clusterProfiler)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
#library(ReactomePA)
#library(trackViewer)
#library(Gviz)
#library(rtracklayer)
#library(Sushi)
#library(venneuler)
library(VennDiagram)
library(readxl)
library(biomaRt)
library(viridis)
library(vegan)
library(ggrepel)



#setwd("/Volumes/cancer/bcat_Project/peakcall/")
#setwd("/Users/yolanda_guillen/Desktop/IMIM/bcat_Project/peakcall_mod/")
setwd("/Users/yguillen/Desktop/temp/beta_catenin_project/peakcall_mod/")

# load blacklisted Regions from modEncode
#blacklist_HUM<-read.delim('/Users/instalar/Desktop/ChIPseq_train/consensusBlacklist.bed',header = FALSE)
#names(blacklist_HUM)<-c("Chr","Start","End","Class","Col5","Col6")

################## MACS2 PEAK CALLING ###############
# After MACS2 peakcalling:
# Six files are created: mypeakds_peaks.narrowPeak, mypeaks_summits.bed, mypeaks_peaks.xls and mypeaks_model.r (and two bedgraphs)

################## MACS2 PEAK CALLING with ChIPseeker ###############

#Bcatenin Control
peakDB <- readPeakFile("DB_peaks_peaks.narrowPeak")
peakDB1 <- readPeakFile("DB1_peaks_peaks.narrowPeak")
peakDB2 <- readPeakFile("DB2_peaks_peaks.narrowPeak")
peakDB3 <- readPeakFile("DB3_peaks_peaks.narrowPeak")
peakBRC1 <- readPeakFile("BRC1_peaks_peaks.narrowPeak")
peakBRC2 <- readPeakFile("BRC2_peaks_peaks.narrowPeak")
peakBRC <- readPeakFile("BRC_peaks_peaks.narrowPeak")

peakRCB1<- readPeakFile("RCB1_peaks_peaks.narrowPeak")



#Bcatenin Lithium
peakBRL1 <- readPeakFile("BRL1_peaks_peaks.narrowPeak")
peakBRL2 <- readPeakFile("BRL2_peaks_peaks.narrowPeak")
#peakBRL <- readPeakFile("BRL_peaks_summits.bed")
peakRBL1 <- readPeakFile("RBL1_peaks_peaks.narrowPeak")
peakRBL2 <- readPeakFile("RBL2_peaks_peaks.narrowPeak")
peakRBL <- readPeakFile("RBL_peaks_peaks.narrowPeak")
#peakRBL_merg <- readPeakFile("RBL_merge_peaks_summits.bed")


# bcatenin control and lithium in Santa Cruz antibody
peakRBCS1 <- readPeakFile("RBCS1_peaks_peaks.narrowPeak")
peakRBCS2 <- readPeakFile("RBCS2_peaks_peaks.narrowPeak")

peakRBLS1 <- readPeakFile("RBLS1_peaks_peaks.narrowPeak")
peakRBLS2 <- readPeakFile("RBLS2_peaks_peaks.narrowPeak")

## DKO bcat
peakRKB1 <-readPeakFile("RKB1_peaks_peaks.narrowPeak")
peakRKB2 <-readPeakFile("RKB2_peaks_peaks.narrowPeak")
peakRKB3 <-readPeakFile("RKB3_peaks_peaks.narrowPeak")

#LEF
peakDL <- readPeakFile("DL_peaks_peaks.narrowPeak")
peakDL1 <- readPeakFile("DL1_peaks_peaks.narrowPeak")
peakDL2 <- readPeakFile("DL2_peaks_peaks.narrowPeak")
peakDL3 <- readPeakFile("DL3_peaks_peaks.narrowPeak")

#TCF
peakDT1 <- readPeakFile("T1R1_peaks_peaks.narrowPeak")
peakDT2 <- readPeakFile("T1R2_peaks_peaks.narrowPeak")
peakDT3 <- readPeakFile("T1R3_peaks_peaks.narrowPeak")
peakDT <- readPeakFile("T1R_peaks_peaks.narrowPeak")

peakRCT1 <- readPeakFile("RCT1_peaks_peaks.narrowPeak")

## DKO TCF
peakRKT1 <-readPeakFile("RKT1_peaks_peaks.narrowPeak")
peakRKT2 <-readPeakFile("RKT2_peaks_peaks.narrowPeak")
peakRKT3 <-readPeakFile("RKT3_peaks_peaks.narrowPeak")

# Histone acet control and lithium
peakRHC1<-readPeakFile("RHC1_peaks_peaks.broadPeak")
peakRHC2<-readPeakFile("RHC2_peaks_peaks.broadPeak")
peakRHC<-readPeakFile("RHC_peaks_peaks.broadPeak")

peakRHL1<-readPeakFile("RHL1_peaks_peaks.broadPeak")
peakRHL2<-readPeakFile("RHL2_peaks_peaks.broadPeak")
peakRHL<-readPeakFile("RHL_peaks_peaks.broadPeak")



#Overall all chromosmes peak coverage
covplot(peakDB, weightCol="FE")
covplot(peakBRC,weightCol = "FE")
covplot(peakRBL,weightCol = "FE")
covplot(peakBRL2,weightCol = "FE")

covplot(peakRBCS1,weightCol = "FE")
covplot(peakRBCS2,weightCol = "FE")

covplot(peakRBLS1,weightCol = "FE")
covplot(peakRBLS2,weightCol = "FE")

covplot(peakDL, weightCol="FE")
covplot(peakDT, weightCol="FE")

covplot(peakRCB1, weightCol="FE")

peak=GenomicRanges::GRangesList(DB=peakDB,RCB=peakRCB1)
covplot(peak)


#common genes DB,DL and DT

#human genes database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

### PEAK ANNOTATION
#DB Bcat control 
peakAnnoDB <- annotatePeak(peakDB, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDB1 <- annotatePeak(peakDB1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDB2 <- annotatePeak(peakDB2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDB3 <- annotatePeak(peakDB3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnnoBRC <- annotatePeak(peakBRC, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoBRC1 <- annotatePeak(peakBRC1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoBRC2 <- annotatePeak(peakBRC2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnnoRCB1 <- annotatePeak(peakRCB1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

# bcat Lithium
#peakAnnoBRL <- annotatePeak(peakBRL, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoBRL1 <- annotatePeak(peakBRL1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoBRL2 <- annotatePeak(peakBRL2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRBL1 <- annotatePeak(peakRBL1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRBL2 <- annotatePeak(peakRBL2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRBL <- annotatePeak(peakRBL, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
#peakAnnoRBL_merg <- annotatePeak(peakRBL_merg, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

#RBC and RBL Bcat Santa Cruz 
peakAnnoRBCS1 <- annotatePeak(peakRBCS1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRBCS2 <- annotatePeak(peakRBCS2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRBLS1 <- annotatePeak(peakRBLS1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRBLS2 <- annotatePeak(peakRBLS2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")



# DKO bcat
peakAnnoRKB1 <- annotatePeak(peakRKB1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRKB2 <- annotatePeak(peakRKB2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRKB3 <- annotatePeak(peakRKB3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

# Histone control and lithium
peakAnnoRHC1 <- annotatePeak(peakRHC1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRHC2 <- annotatePeak(peakRHC2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRHC <- annotatePeak(peakRHC, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnnoRHL <- annotatePeak(peakRHL, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRHL1 <- annotatePeak(peakRHL1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRHL2 <- annotatePeak(peakRHL2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

#DL LEF
peakAnnoDL <- annotatePeak(peakDL, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDL1 <- annotatePeak(peakDL1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDL2 <- annotatePeak(peakDL2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDL3 <- annotatePeak(peakDL3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

#DT TCF
peakAnnoDT <- annotatePeak(peakDT, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDT1 <- annotatePeak(peakDT1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDT2 <- annotatePeak(peakDT2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDT3 <- annotatePeak(peakDT3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnnoRCT1 <- annotatePeak(peakRCT1,tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

# DKO TCF
peakAnnoRKT1 <- annotatePeak(peakRKT1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRKT2 <- annotatePeak(peakRKT2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoRKT3 <- annotatePeak(peakRKT3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

# Visualize genomic annotation
plotAnnoPie(peakAnnoBRL2)
dev.off()
plotAnnoBar(peakAnnoDB)
dev.off()
vennpie(peakAnnoDB)
dev.off()

plotAnnoPie(peakAnnoDL)
dev.off()
plotAnnoBar(peakAnnoDL)
dev.off()
vennpie(peakAnnoDL)
dev.off()

plotAnnoPie(peakAnnoDT)
dev.off()
plotAnnoBar(peakAnnoDT)
dev.off()
vennpie(peakAnnoDT)
dev.off()


upsetplot(peakAnnoDB, vennpie=TRUE)
dev.off()
upsetplot(peakAnnoDL, vennpie=TRUE)
dev.off()
upsetplot(peakAnnoDT, vennpie=TRUE)
dev.off()

plotDistToTSS(peakAnnoDB,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

plotDistToTSS(peakAnnoBRC2,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

plotDistToTSS(peakAnnoBRL,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

plotDistToTSS(peakAnnoDL,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

plotDistToTSS(peakAnnoDT,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

plotDistToTSS(peakAnnoRHC,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


##DATAFRAME WITH genes id linked to peaks

#bcatenin control BRC merge
dataf_peakannoBRC<-as.data.frame(peakAnnoBRC)
dataf_peakannoBRC$annotation<-gsub('Intron.*','Intron',dataf_peakannoBRC$annotation)
dataf_peakannoBRC$annotation<-gsub('Exon.*','Exon',dataf_peakannoBRC$annotation)
dataf_peakannoBRC$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBRC$annotation)
dataf_peakannoBRC$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBRC$annotation)
dataf_peakannoBRC$peaksite<-paste(dataf_peakannoBRC$SYMBOL,dataf_peakannoBRC$annotation,sep='_')
table(dataf_peakannoBRC$annotation)


#bcatenin control BRC1
dataf_peakannoBRC1<-as.data.frame(peakAnnoBRC1)
dataf_peakannoBRC1$annotation<-gsub('Intron.*','Intron',dataf_peakannoBRC1$annotation)
dataf_peakannoBRC1$annotation<-gsub('Exon.*','Exon',dataf_peakannoBRC1$annotation)
dataf_peakannoBRC1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBRC1$annotation)
dataf_peakannoBRC1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBRC1$annotation)
dataf_peakannoBRC1$peaksite<-paste(dataf_peakannoBRC1$SYMBOL,dataf_peakannoBRC1$annotation,sep='_')
table(dataf_peakannoBRC1$annotation)

#bcatenin control BRC2
dataf_peakannoBRC2<-as.data.frame(peakAnnoBRC2)
dataf_peakannoBRC2$annotation<-gsub('Intron.*','Intron',dataf_peakannoBRC2$annotation)
dataf_peakannoBRC2$annotation<-gsub('Exon.*','Exon',dataf_peakannoBRC2$annotation)
dataf_peakannoBRC2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBRC2$annotation)
dataf_peakannoBRC2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBRC2$annotation)
dataf_peakannoBRC2$peaksite<-paste(dataf_peakannoBRC2$SYMBOL,dataf_peakannoBRC2$annotation,sep='_')
table(dataf_peakannoBRC2$annotation)


#bcatenin control RCB1
dataf_peakannoRCB1<-as.data.frame(peakAnnoRCB1)
dataf_peakannoRCB1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRCB1$annotation)
dataf_peakannoRCB1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRCB1$annotation)
dataf_peakannoRCB1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRCB1$annotation)
dataf_peakannoRCB1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRCB1$annotation)
dataf_peakannoRCB1$peaksite<-paste(dataf_peakannoRCB1$SYMBOL,dataf_peakannoRCB1$annotation,sep='_')
table(dataf_peakannoRCB1$annotation)

#BRL1
dataf_peakannoBRL1<-as.data.frame(peakAnnoBRL1)
dataf_peakannoBRL1$annotation<-gsub('Intron.*','Intron',dataf_peakannoBRL1$annotation)
dataf_peakannoBRL1$annotation<-gsub('Exon.*','Exon',dataf_peakannoBRL1$annotation)
dataf_peakannoBRL1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBRL1$annotation)
dataf_peakannoBRL1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBRL1$annotation)
dataf_peakannoBRL1$peaksite<-paste(dataf_peakannoBRL1$SYMBOL,dataf_peakannoBRL1$annotation,sep='_')
table(dataf_peakannoBRL1$annotation)


#BRL2
dataf_peakannoBRL2<-as.data.frame(peakAnnoBRL2)
dataf_peakannoBRL2$annotation<-gsub('Intron.*','Intron',dataf_peakannoBRL2$annotation)
dataf_peakannoBRL2$annotation<-gsub('Exon.*','Exon',dataf_peakannoBRL2$annotation)
dataf_peakannoBRL2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBRL2$annotation)
dataf_peakannoBRL2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBRL2$annotation)
dataf_peakannoBRL2$peaksite<-paste(dataf_peakannoBRL2$SYMBOL,dataf_peakannoBRL2$annotation,sep='_')
table(dataf_peakannoBRL2$annotation)

#RBL1
dataf_peakannoRBL1<-as.data.frame(peakAnnoRBL1)
dataf_peakannoRBL1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRBL1$annotation)
dataf_peakannoRBL1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRBL1$annotation)
dataf_peakannoRBL1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRBL1$annotation)
dataf_peakannoRBL1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRBL1$annotation)
dataf_peakannoRBL1$peaksite<-paste(dataf_peakannoRBL1$SYMBOL,dataf_peakannoRBL1$annotation,sep='_')
table(dataf_peakannoRBL1$annotation)

#RBL2
dataf_peakannoRBL2<-as.data.frame(peakAnnoRBL2)
dataf_peakannoRBL2$annotation<-gsub('Intron.*','Intron',dataf_peakannoRBL2$annotation)
dataf_peakannoRBL2$annotation<-gsub('Exon.*','Exon',dataf_peakannoRBL2$annotation)
dataf_peakannoRBL2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRBL2$annotation)
dataf_peakannoRBL2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRBL2$annotation)
dataf_peakannoRBL2$peaksite<-paste(dataf_peakannoRBL2$SYMBOL,dataf_peakannoRBL2$annotation,sep='_')
table(dataf_peakannoRBL2$annotation)

#RBL
dataf_peakannoRBL<-as.data.frame(peakAnnoRBL)
dataf_peakannoRBL$annotation<-gsub('Intron.*','Intron',dataf_peakannoRBL$annotation)
dataf_peakannoRBL$annotation<-gsub('Exon.*','Exon',dataf_peakannoRBL$annotation)
dataf_peakannoRBL$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRBL$annotation)
dataf_peakannoRBL$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRBL$annotation)
dataf_peakannoRBL$peaksite<-paste(dataf_peakannoRBL$SYMBOL,dataf_peakannoRBL$annotation,sep='_')
table(dataf_peakannoRBL$annotation)

#RBL_MERGE
#dataf_peakannoRBL_merg<-as.data.frame(peakAnnoRBL_merg)
#dataf_peakannoRBL_merg$annotation<-gsub('Intron.*','Intron',dataf_peakannoRBL_merg$annotation)
#dataf_peakannoRBL_merg$annotation<-gsub('Exon.*','Exon',dataf_peakannoRBL_merg$annotation)
#dataf_peakannoRBL_merg$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRBL_merg$annotation)
#dataf_peakannoRBL_merg$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRBL_merg$annotation)
#dataf_peakannoRBL_merg$peaksite<-paste(dataf_peakannoRBL_merg$SYMBOL,dataf_peakannoRBL_merg$annotation,sep='_')
#table(dataf_peakannoRBL_merg$annotation)


#DB bcat
# DB1
dataf_peakannoDB1<-as.data.frame(peakAnnoDB1)
dataf_peakannoDB1$annotation<-gsub('Intron.*','Intron',dataf_peakannoDB1$annotation)
dataf_peakannoDB1$annotation<-gsub('Exon.*','Exon',dataf_peakannoDB1$annotation)
dataf_peakannoDB1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDB1$annotation)
dataf_peakannoDB1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDB1$annotation)
dataf_peakannoDB1$peaksite<-paste(dataf_peakannoDB1$SYMBOL,dataf_peakannoDB1$annotation,sep='_')
table(dataf_peakannoDB1$annotation)

#DB2
dataf_peakannoDB2<-as.data.frame(peakAnnoDB2)
dataf_peakannoDB2$annotation<-gsub('Intron.*','Intron',dataf_peakannoDB2$annotation)
dataf_peakannoDB2$annotation<-gsub('Exon.*','Exon',dataf_peakannoDB2$annotation)
dataf_peakannoDB2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDB2$annotation)
dataf_peakannoDB2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDB2$annotation)
dataf_peakannoDB2$peaksite<-paste(dataf_peakannoDB2$SYMBOL,dataf_peakannoDB2$annotation,sep='_')
table(dataf_peakannoDB2$annotation)

#DB3
dataf_peakannoDB3<-as.data.frame(peakAnnoDB3)
dataf_peakannoDB3$annotation<-gsub('Intron.*','Intron',dataf_peakannoDB3$annotation)
dataf_peakannoDB3$annotation<-gsub('Exon.*','Exon',dataf_peakannoDB3$annotation)
dataf_peakannoDB3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDB3$annotation)
dataf_peakannoDB3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDB3$annotation)
dataf_peakannoDB3$peaksite<-paste(dataf_peakannoDB3$SYMBOL,dataf_peakannoDB3$annotation,sep='_')
table(dataf_peakannoDB3$annotation)

#Merge DB
dataf_peakannoDB<-as.data.frame(peakAnnoDB)
dataf_peakannoDB$annotation<-gsub('Intron.*','Intron',dataf_peakannoDB$annotation)
dataf_peakannoDB$annotation<-gsub('Exon.*','Exon',dataf_peakannoDB$annotation)
dataf_peakannoDB$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDB$annotation)
dataf_peakannoDB$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDB$annotation)
dataf_peakannoDB$peaksite<-paste(dataf_peakannoDB2$SYMBOL,dataf_peakannoDB$annotation,sep='_')
table(dataf_peakannoDB$annotation)

#bcatenin SANTA CRUZ antibodies 
# control RBCS1
dataf_peakannoRBCS1<-as.data.frame(peakAnnoRBCS1)
dataf_peakannoRBCS1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRBCS1$annotation)
dataf_peakannoRBCS1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRBCS1$annotation)
dataf_peakannoRBCS1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRBCS1$annotation)
dataf_peakannoRBCS1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRBCS1$annotation)
dataf_peakannoRBCS1$peaksite<-paste(dataf_peakannoRBCS1$SYMBOL,dataf_peakannoRBCS1$annotation,sep='_')
table(dataf_peakannoRBCS1$annotation)

# control RBCS2
dataf_peakannoRBCS2<-as.data.frame(peakAnnoRBCS2)
dataf_peakannoRBCS2$annotation<-gsub('Intron.*','Intron',dataf_peakannoRBCS2$annotation)
dataf_peakannoRBCS2$annotation<-gsub('Exon.*','Exon',dataf_peakannoRBCS2$annotation)
dataf_peakannoRBCS2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRBCS2$annotation)
dataf_peakannoRBCS2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRBCS2$annotation)
dataf_peakannoRBCS2$peaksite<-paste(dataf_peakannoRBCS2$SYMBOL,dataf_peakannoRBCS2$annotation,sep='_')
table(dataf_peakannoRBCS2$annotation)

# lithium RBLS1
dataf_peakannoRBLS1<-as.data.frame(peakAnnoRBLS1)
dataf_peakannoRBLS1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRBLS1$annotation)
dataf_peakannoRBLS1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRBLS1$annotation)
dataf_peakannoRBLS1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRBLS1$annotation)
dataf_peakannoRBLS1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRBLS1$annotation)
dataf_peakannoRBLS1$peaksite<-paste(dataf_peakannoRBLS1$SYMBOL,dataf_peakannoRBLS1$annotation,sep='_')
table(dataf_peakannoRBLS1$annotation)

# lithium RBLS2
dataf_peakannoRBLS2<-as.data.frame(peakAnnoRBLS2)
dataf_peakannoRBLS2$annotation<-gsub('Intron.*','Intron',dataf_peakannoRBLS2$annotation)
dataf_peakannoRBLS2$annotation<-gsub('Exon.*','Exon',dataf_peakannoRBLS2$annotation)
dataf_peakannoRBLS2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRBLS2$annotation)
dataf_peakannoRBLS2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRBLS2$annotation)
dataf_peakannoRBLS2$peaksite<-paste(dataf_peakannoRBLS2$SYMBOL,dataf_peakannoRBLS2$annotation,sep='_')
table(dataf_peakannoRBLS2$annotation)



## bcat DKO
# RKB1
dataf_peakannoRKB1<-as.data.frame(peakAnnoRKB1)
dataf_peakannoRKB1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRKB1$annotation)
dataf_peakannoRKB1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRKB1$annotation)
dataf_peakannoRKB1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRKB1$annotation)
dataf_peakannoRKB1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRKB1$annotation)
dataf_peakannoRKB1$peaksite<-paste(dataf_peakannoRKB1$SYMBOL,dataf_peakannoRKB1$annotation,sep='_')
table(dataf_peakannoRKB1$annotation)

# RKB2
dataf_peakannoRKB2<-as.data.frame(peakAnnoRKB2)
dataf_peakannoRKB2$annotation<-gsub('Intron.*','Intron',dataf_peakannoRKB2$annotation)
dataf_peakannoRKB2$annotation<-gsub('Exon.*','Exon',dataf_peakannoRKB2$annotation)
dataf_peakannoRKB2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRKB2$annotation)
dataf_peakannoRKB2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRKB2$annotation)
dataf_peakannoRKB2$peaksite<-paste(dataf_peakannoRKB2$SYMBOL,dataf_peakannoRKB2$annotation,sep='_')
table(dataf_peakannoRKB2$annotation)

# RKB3
dataf_peakannoRKB3<-as.data.frame(peakAnnoRKB3)
dataf_peakannoRKB3$annotation<-gsub('Intron.*','Intron',dataf_peakannoRKB3$annotation)
dataf_peakannoRKB3$annotation<-gsub('Exon.*','Exon',dataf_peakannoRKB3$annotation)
dataf_peakannoRKB3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRKB3$annotation)
dataf_peakannoRKB3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRKB3$annotation)
dataf_peakannoRKB3$peaksite<-paste(dataf_peakannoRKB3$SYMBOL,dataf_peakannoRKB3$annotation,sep='_')
table(dataf_peakannoRKB3$annotation)


## bcat histones control 

dataf_peakannoRHC1<-as.data.frame(peakAnnoRHC1)
dataf_peakannoRHC1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRHC1$annotation)
dataf_peakannoRHC1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRHC1$annotation)
dataf_peakannoRHC1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRHC1$annotation)
dataf_peakannoRHC1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRHC1$annotation)
dataf_peakannoRHC1$peaksite<-paste(dataf_peakannoRHC1$SYMBOL,dataf_peakannoRHC1$annotation,sep='_')
table(dataf_peakannoRHC1$annotation)

dataf_peakannoRHC2<-as.data.frame(peakAnnoRHC2)
dataf_peakannoRHC2$annotation<-gsub('Intron.*','Intron',dataf_peakannoRHC2$annotation)
dataf_peakannoRHC2$annotation<-gsub('Exon.*','Exon',dataf_peakannoRHC2$annotation)
dataf_peakannoRHC2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRHC2$annotation)
dataf_peakannoRHC2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRHC2$annotation)
dataf_peakannoRHC2$peaksite<-paste(dataf_peakannoRHC2$SYMBOL,dataf_peakannoRHC2$annotation,sep='_')
table(dataf_peakannoRHC2$annotation)

dataf_peakannoRHC<-as.data.frame(peakAnnoRHC)
dataf_peakannoRHC$annotation<-gsub('Intron.*','Intron',dataf_peakannoRHC$annotation)
dataf_peakannoRHC$annotation<-gsub('Exon.*','Exon',dataf_peakannoRHC$annotation)
dataf_peakannoRHC$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRHC$annotation)
dataf_peakannoRHC$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRHC$annotation)
dataf_peakannoRHC$peaksite<-paste(dataf_peakannoRHC$SYMBOL,dataf_peakannoRHC$annotation,sep='_')
table(dataf_peakannoRHC$annotation)



## bcat histones lithium

dataf_peakannoRHL1<-as.data.frame(peakAnnoRHL1)
dataf_peakannoRHL1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRHL1$annotation)
dataf_peakannoRHL1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRHL1$annotation)
dataf_peakannoRHL1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRHL1$annotation)
dataf_peakannoRHL1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRHL1$annotation)
dataf_peakannoRHL1$peaksite<-paste(dataf_peakannoRHL1$SYMBOL,dataf_peakannoRHL1$annotation,sep='_')
table(dataf_peakannoRHL1$annotation)


dataf_peakannoRHL2<-as.data.frame(peakAnnoRHL2)
dataf_peakannoRHL2$annotation<-gsub('Intron.*','Intron',dataf_peakannoRHL2$annotation)
dataf_peakannoRHL2$annotation<-gsub('Exon.*','Exon',dataf_peakannoRHL2$annotation)
dataf_peakannoRHL2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRHL2$annotation)
dataf_peakannoRHL2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRHL2$annotation)
dataf_peakannoRHL2$peaksite<-paste(dataf_peakannoRHL2$SYMBOL,dataf_peakannoRHL2$annotation,sep='_')
table(dataf_peakannoRHL2$annotation)

dataf_peakannoRHL<-as.data.frame(peakAnnoRHL)
dataf_peakannoRHL$annotation<-gsub('Intron.*','Intron',dataf_peakannoRHL$annotation)
dataf_peakannoRHL$annotation<-gsub('Exon.*','Exon',dataf_peakannoRHL$annotation)
dataf_peakannoRHL$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRHL$annotation)
dataf_peakannoRHL$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRHL$annotation)
dataf_peakannoRHL$peaksite<-paste(dataf_peakannoRHL$SYMBOL,dataf_peakannoRHL$annotation,sep='_')
table(dataf_peakannoRHL$annotation)

# Write histone peaks
write.table(dataf_peakannoRHC1,"/Volumes/grcmc/DAVID/Yolanda_data/RHC1_histone_bcat.csv",quote = FALSE,row.names = FALSE,sep="\t")
write.table(dataf_peakannoRHC2,"/Volumes/grcmc/DAVID/Yolanda_data/RHC2_histone_bcat.csv",quote = FALSE,row.names = FALSE,sep="\t")
write.table(dataf_peakannoRHL1,"/Volumes/grcmc/DAVID/Yolanda_data/RHL1_histone_bcat.csv",quote = FALSE,row.names = FALSE,sep="\t")
write.table(dataf_peakannoRHL2,"/Volumes/grcmc/DAVID/Yolanda_data/RHL2_histone_bcat.csv",quote = FALSE,row.names = FALSE,sep="\t")

#LEF1

# DL1
dataf_peakannoDL1<-as.data.frame(peakAnnoDL1)
dataf_peakannoDL1$annotation<-gsub('Intron.*','Intron',dataf_peakannoDL1$annotation)
dataf_peakannoDL1$annotation<-gsub('Exon.*','Exon',dataf_peakannoDL1$annotation)
dataf_peakannoDL1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDL1$annotation)
dataf_peakannoDL1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDL1$annotation)
dataf_peakannoDL1$peaksite<-paste(dataf_peakannoDL1$SYMBOL,dataf_peakannoDL1$annotation,sep='_')
table(dataf_peakannoDL1$annotation)


# DL2
dataf_peakannoDL2<-as.data.frame(peakAnnoDL2)
dataf_peakannoDL2$annotation<-gsub('Intron.*','Intron',dataf_peakannoDL2$annotation)
dataf_peakannoDL2$annotation<-gsub('Exon.*','Exon',dataf_peakannoDL2$annotation)
dataf_peakannoDL2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDL2$annotation)
dataf_peakannoDL2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDL2$annotation)
dataf_peakannoDL2$peaksite<-paste(dataf_peakannoDL2$SYMBOL,dataf_peakannoDL2$annotation,sep='_')
table(dataf_peakannoDL2$annotation)

#DL3
dataf_peakannoDL3<-as.data.frame(peakAnnoDL3)
dataf_peakannoDL3$annotation<-gsub('Intron.*','Intron',dataf_peakannoDL3$annotation)
dataf_peakannoDL3$annotation<-gsub('Exon.*','Exon',dataf_peakannoDL3$annotation)
dataf_peakannoDL3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDL3$annotation)
dataf_peakannoDL3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDL3$annotation)
dataf_peakannoDL3$peaksite<-paste(dataf_peakannoDL3$SYMBOL,dataf_peakannoDL3$annotation,sep='_')
table(dataf_peakannoDL3$annotation)


# Merge DL
dataf_peakannoDL<-as.data.frame(peakAnnoDL)
dataf_peakannoDL$annotation<-gsub('Intron.*','Intron',dataf_peakannoDL$annotation)
dataf_peakannoDL$annotation<-gsub('Exon.*','Exon',dataf_peakannoDL$annotation)
dataf_peakannoDL$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDL$annotation)
dataf_peakannoDL$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDL$annotation)
dataf_peakannoDL$peaksite<-paste(dataf_peakannoDL$SYMBOL,dataf_peakannoDL$annotation,sep='_')
table(dataf_peakannoDL$annotation)


# Merge DT
dataf_peakannoDT<-as.data.frame(peakAnnoDT)
dataf_peakannoDT$annotation<-gsub('Intron.*','Intron',dataf_peakannoDT$annotation)
dataf_peakannoDT$annotation<-gsub('Exon.*','Exon',dataf_peakannoDT$annotation)
dataf_peakannoDT$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDT$annotation)
dataf_peakannoDT$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDT$annotation)
dataf_peakannoDT$peaksite<-paste(dataf_peakannoDT$SYMBOL,dataf_peakannoDT$annotation,sep='_')
table(dataf_peakannoDT$annotation)

#DT1
dataf_peakannoDT1<-as.data.frame(peakAnnoDT1)
dataf_peakannoDT1$annotation<-gsub('Intron.*','Intron',dataf_peakannoDT1$annotation)
dataf_peakannoDT1$annotation<-gsub('Exon.*','Exon',dataf_peakannoDT1$annotation)
dataf_peakannoDT1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDT1$annotation)
dataf_peakannoDT1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDT1$annotation)
dataf_peakannoDT1$peaksite<-paste(dataf_peakannoDT1$SYMBOL,dataf_peakannoDT1$annotation,sep='_')
table(dataf_peakannoDT1$annotation)

#DT2
dataf_peakannoDT2<-as.data.frame(peakAnnoDT2)
dataf_peakannoDT2$annotation<-gsub('Intron.*','Intron',dataf_peakannoDT2$annotation)
dataf_peakannoDT2$annotation<-gsub('Exon.*','Exon',dataf_peakannoDT2$annotation)
dataf_peakannoDT2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDT2$annotation)
dataf_peakannoDT2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDT2$annotation)
dataf_peakannoDT2$peaksite<-paste(dataf_peakannoDT2$SYMBOL,dataf_peakannoDT2$annotation,sep='_')
table(dataf_peakannoDT2$annotation)

#DT3
dataf_peakannoDT3<-as.data.frame(peakAnnoDT3)
dataf_peakannoDT3$annotation<-gsub('Intron.*','Intron',dataf_peakannoDT3$annotation)
dataf_peakannoDT3$annotation<-gsub('Exon.*','Exon',dataf_peakannoDT3$annotation)
dataf_peakannoDT3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDT3$annotation)
dataf_peakannoDT3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDT3$annotation)
dataf_peakannoDT3$peaksite<-paste(dataf_peakannoDT3$SYMBOL,dataf_peakannoDT3$annotation,sep='_')
table(dataf_peakannoDT3$annotation)

#DT
dataf_peakannoDT<-as.data.frame(peakAnnoDT)
dataf_peakannoDT$annotation<-gsub('Intron.*','Intron',dataf_peakannoDT$annotation)
dataf_peakannoDT$annotation<-gsub('Exon.*','Exon',dataf_peakannoDT$annotation)
dataf_peakannoDT$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoDT$annotation)
dataf_peakannoDT$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoDT$annotation)
dataf_peakannoDT$peaksite<-paste(dataf_peakannoDT$SYMBOL,dataf_peakannoDT$annotation,sep='_')
table(dataf_peakannoDT$annotation)

#RCT1
dataf_peakannoRCT1<-as.data.frame(peakAnnoRCT1)
dataf_peakannoRCT1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRCT1$annotation)
dataf_peakannoRCT1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRCT1$annotation)
dataf_peakannoRCT1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRCT1$annotation)
dataf_peakannoRCT1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRCT1$annotation)
dataf_peakannoRCT1$peaksite<-paste(dataf_peakannoRCT1$SYMBOL,dataf_peakannoRCT1$annotation,sep='_')
table(dataf_peakannoRCT1$annotation)


# DKO TCF
#RKT1
dataf_peakannoRKT1<-as.data.frame(peakAnnoRKT1)
dataf_peakannoRKT1$annotation<-gsub('Intron.*','Intron',dataf_peakannoRKT1$annotation)
dataf_peakannoRKT1$annotation<-gsub('Exon.*','Exon',dataf_peakannoRKT1$annotation)
dataf_peakannoRKT1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRKT1$annotation)
dataf_peakannoRKT1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRKT1$annotation)
dataf_peakannoRKT1$peaksite<-paste(dataf_peakannoRKT1$SYMBOL,dataf_peakannoRKT1$annotation,sep='_')
table(dataf_peakannoRKT1$annotation)

#RKT2
dataf_peakannoRKT2<-as.data.frame(peakAnnoRKT2)
dataf_peakannoRKT2$annotation<-gsub('Intron.*','Intron',dataf_peakannoRKT2$annotation)
dataf_peakannoRKT2$annotation<-gsub('Exon.*','Exon',dataf_peakannoRKT2$annotation)
dataf_peakannoRKT2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRKT2$annotation)
dataf_peakannoRKT2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRKT2$annotation)
dataf_peakannoRKT2$peaksite<-paste(dataf_peakannoRKT2$SYMBOL,dataf_peakannoRKT2$annotation,sep='_')
table(dataf_peakannoRKT2$annotation)

#RKT3
dataf_peakannoRKT3<-as.data.frame(peakAnnoRKT3)
dataf_peakannoRKT3$annotation<-gsub('Intron.*','Intron',dataf_peakannoRKT3$annotation)
dataf_peakannoRKT3$annotation<-gsub('Exon.*','Exon',dataf_peakannoRKT3$annotation)
dataf_peakannoRKT3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoRKT3$annotation)
dataf_peakannoRKT3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoRKT3$annotation)
dataf_peakannoRKT3$peaksite<-paste(dataf_peakannoRKT3$SYMBOL,dataf_peakannoRKT3$annotation,sep='_')
table(dataf_peakannoRKT3$annotation)


# Peak distribution
peakAnnoBRC1@annoStat


#BRC
distBRC1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBRC1$annotation)))))
colnames(distBRC1) <- as.character(unlist(distBRC1[1,]))
distBRC1<-distBRC1[-1,]
distBRC1$Sample<-"BRC1"

distBRC2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBRC2$annotation)))))
colnames(distBRC2) <- as.character(unlist(distBRC2[1,]))
distBRC2<-distBRC2[-1,]
distBRC2$Sample<-"BRC2"

distBRC<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBRC$annotation)))))
colnames(distBRC) <- as.character(unlist(distBRC[1,]))
distBRC<-distBRC[-1,]
distBRC$Sample<-"BRC"


#BRL

distBRL1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBRL1$annotation)))))
colnames(distBRL1) <- as.character(unlist(distBRL1[1,]))
distBRL1<-distBRL1[-1,]
distBRL1$Sample<-"BRL1"

distBRL2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBRL2$annotation)))))
colnames(distBRL2) <- as.character(unlist(distBRL2[1,]))
distBRL2<-distBRL2[-1,]
distBRL2$Sample<-"BRL2"

distRBL1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRBL1$annotation)))))
colnames(distRBL1) <- as.character(unlist(distRBL1[1,]))
distRBL1<-distRBL1[-1,]
distRBL1$Sample<-"RBL1"

distRBL2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRBL2$annotation)))))
colnames(distRBL2) <- as.character(unlist(distRBL2[1,]))
distRBL2<-distRBL2[-1,]
distRBL2$Sample<-"RBL2"

distRBL<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRBL$annotation)))))
colnames(distRBL) <- as.character(unlist(distRBL[1,]))
distRBL<-distRBL[-1,]
distRBL$Sample<-"RBL"

#distRBL_merg<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRBL_merg$annotation)))))
#colnames(distRBL_merg) <- as.character(unlist(distRBL_merg[1,]))
#distRBL_merg<-distRBL_merg[-1,]
#distRBL_merg$Sample<-"RBL_merg"

#DB

distDB1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDB1$annotation)))))
colnames(distDB1) <- as.character(unlist(distDB1[1,]))
distDB1<-distDB1[-1,]
distDB1$Sample<-"DB1"

distDB2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDB2$annotation)))))
colnames(distDB2) <- as.character(unlist(distDB2[1,]))
distDB2<-distDB2[-1,]
distDB2$Sample<-"DB2"

distDB3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDB3$annotation)))))
colnames(distDB3) <- as.character(unlist(distDB3[1,]))
distDB3<-distDB3[-1,]
distDB3$Sample<-"DB3"

distDB<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDB$annotation)))))
colnames(distDB) <- as.character(unlist(distDB[1,]))
distDB<-distDB[-1,]
distDB$Sample<-"DB"

# RCB1
distRCB1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRCB1$annotation)))))
colnames(distRCB1) <- as.character(unlist(distRCB1[1,]))
distRCB1<-distRCB1[-1,]
distRCB1$Sample<-"RCB1"

# RBCS1
distRBCS1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRBCS1$annotation)))))
colnames(distRBCS1) <- as.character(unlist(distRBCS1[1,]))
distRBCS1<-distRBCS1[-1,]
distRBCS1$Sample<-"RBCS1"

# RBCS2
distRBCS2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRBCS2$annotation)))))
colnames(distRBCS2) <- as.character(unlist(distRBCS2[1,]))
distRBCS2<-distRBCS2[-1,]
distRBCS2$Sample<-"RBCS2"


# RBLS1
distRBLS1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRBLS1$annotation)))))
colnames(distRBLS1) <- as.character(unlist(distRBLS1[1,]))
distRBLS1<-distRBLS1[-1,]
distRBLS1$Sample<-"RBLS1"


# RBLS2
distRBLS2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRBLS2$annotation)))))
colnames(distRBLS2) <- as.character(unlist(distRBLS2[1,]))
distRBLS2<-distRBLS2[-1,]
distRBLS2$Sample<-"RBLS2"

# DKO RKB
distRKB1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRKB1$annotation)))))
colnames(distRKB1) <- as.character(unlist(distRKB1[1,]))
distRKB1<-distRKB1[-1,]
distRKB1$Sample<-"RKB1"

colnames(distRKB1)
distRKB1$Three_UTR<-"0"
distRKB1$Five_UTR<-"0"
distRKB1$Downstream<-"0"
distRKB1<-distRKB1[,c(6,7,1,8,2,3,4,5)]
colnames(distRKB1)[3]<-c("Distal_intergenic")

distRKB2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRKB2$annotation)))))
colnames(distRKB2) <- as.character(unlist(distRKB2[1,]))
distRKB2<-distRKB2[-1,]
distRKB2$Sample<-"RKB2"

colnames(distRKB2)
distRKB2$Five_UTR<-"0"
distRKB2$Downstream<-"0"
colnames(distRKB2)[1]<-c("Three_UTR")
colnames(distRKB2)[2]<-c("Distal_intergenic")
distRKB2<-distRKB2[,c(1,7,2,8,3,4,5,6)]


distRKB3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRKB3$annotation)))))
colnames(distRKB3) <- as.character(unlist(distRKB3[1,]))
distRKB3<-distRKB3[-1,]
distRKB3$Sample<-"RKB3"

colnames(distRKB3)
distRKB3$Five_UTR<-"0"
distRKB3$Three_UTR<-"0"
colnames(distRKB3)[1]<-c("Distal_intergenic")

distRKB3<-distRKB3[,c(8,7,1,2,3,4,5,6)]

## Histones CONTROL AND LITHIUM

distRHC1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRHC1$annotation)))))
colnames(distRHC1) <- as.character(unlist(distRHC1[1,]))
distRHC1<-distRHC1[-1,]
distRHC1$Sample<-"RHC1"

distRHC2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRHC2$annotation)))))
colnames(distRHC2) <- as.character(unlist(distRHC2[1,]))
distRHC2<-distRHC2[-1,]
distRHC2$Sample<-"RHC2"

distRHC<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRHC$annotation)))))
colnames(distRHC) <- as.character(unlist(distRHC[1,]))
distRHC<-distRHC[-1,]
distRHC$Sample<-"RHC"

distRHL1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRHL1$annotation)))))
colnames(distRHL1) <- as.character(unlist(distRHL1[1,]))
distRHL1<-distRHL1[-1,]
distRHL1$Sample<-"RHL1"

distRHL2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRHL2$annotation)))))
colnames(distRHL2) <- as.character(unlist(distRHL2[1,]))
distRHL2<-distRHL2[-1,]
distRHL2$Sample<-"RHL2"

distRHL<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRHL$annotation)))))
colnames(distRHL) <- as.character(unlist(distRHL[1,]))
distRHL<-distRHL[-1,]
distRHL$Sample<-"RHL"

#DL1
distDL1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDL1$annotation)))))
colnames(distDL1) <- as.character(unlist(distDL1[1,]))
distDL1<-distDL1[-1,]
distDL1$Sample<-"DL1"

#DL2
distDL2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDL2$annotation)))))
colnames(distDL2) <- as.character(unlist(distDL2[1,]))
distDL2<-distDL2[-1,]
distDL2$Sample<-"DL2"

#DL3
distDL3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDL3$annotation)))))
colnames(distDL3) <- as.character(unlist(distDL3[1,]))
distDL3<-distDL3[-1,]
distDL3$Sample<-"DL3"

#DL
distDL<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDL$annotation)))))
colnames(distDL) <- as.character(unlist(distDL[1,]))
distDL<-distDL[-1,]
distDL$Sample<-"DL"

#DT1
distDT1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDT1$annotation)))))
colnames(distDT1) <- as.character(unlist(distDT1[1,]))
distDT1<-distDT1[-1,]
distDT1$Sample<-"DT1"

#DT2
distDT2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDT2$annotation)))))
colnames(distDT2) <- as.character(unlist(distDT2[1,]))
distDT2<-distDT2[-1,]
distDT2$Sample<-"DT2"

#DT3
distDT3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDT3$annotation)))))
colnames(distDT3) <- as.character(unlist(distDT3[1,]))
distDT3<-distDT3[-1,]
distDT3$Sample<-"DT3"

#DT merge
distDT<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoDT$annotation)))))
colnames(distDT) <- as.character(unlist(distDT[1,]))
distDT<-distDT[-1,]
distDT$Sample<-"DT"

# RCT1
distRCT1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRCT1$annotation)))))
colnames(distRCT1) <- as.character(unlist(distRCT1[1,]))
distRCT1<-distRCT1[-1,]
distRCT1$Sample<-"RCT1"

colnames(distRCT1)

# DKO TCF
distRKT1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRKT1$annotation)))))
colnames(distRKT1) <- as.character(unlist(distRKT1[1,]))
distRKT1<-distRKT1[-1,]
distRKT1$Sample<-"RKT1"

colnames(distRKT1)[1]<-c("Three_UTR")
colnames(distRKT1)[2]<-c("Five_UTR")
colnames(distRKT1)[3]<-c("Distal_intergenic")

distRKT2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRKT2$annotation)))))
colnames(distRKT2) <- as.character(unlist(distRKT2[1,]))
distRKT2<-distRKT2[-1,]
distRKT2$Sample<-"RKT2"

colnames(distRKT2)
colnames(distRKT2)[1]<-c("Three_UTR")
colnames(distRKT2)[2]<-c("Distal_intergenic")
distRKT2$Five_UTR<-0

colnames(distRKT2)
distRKT2<-distRKT2[,c(1,8,2,3,4,5,6,7)]


distRKT3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoRKT3$annotation)))))
colnames(distRKT3) <- as.character(unlist(distRKT3[1,]))
distRKT3<-distRKT3[-1,]
distRKT3$Sample<-"RKT3"

colnames(distRKT3)[1]<-c("Three_UTR")
colnames(distRKT3)[2]<-c("Five_UTR")
colnames(distRKT3)[3]<-c("Distal_intergenic")



distRBL1$Downstream<-"0"
distRBL1<-distRBL1[,c(1,2,7,3,4,5,6)]
distRBL2$Downstream<-"0"
distRBL2<-distRBL2[,c(1,2,7,3,4,5,6)]
distRBL$Downstream<-"0"
distRBL<-distRBL[,c(1,2,7,3,4,5,6)]
#distRBL_merg$Downstream<-"0"
#distRBL_merg<-distRBL_merg[,c(1,2,7,3,4,5,6)]

## DO NOT INCLUDE RCB1, BRL1 nor RBL1

distBcat<-rbind(distBRC1,distBRC2,
                distDB1,distDB2,distDB3,
                distBRL2,distRBL2,distRBL)

colnames(distBcat)[1]<-c("Three_UTR")
colnames(distBcat)[2]<-c("Distal_intergenic")
distBcat$Five_UTR<-"0"
distBcat<-distBcat[,c(1,8,2,3,4,5,6,7)]

colnames(distRCB1)[1]<-c("Three_UTR")
colnames(distRCB1)[2]<-c("Five_UTR")
colnames(distRCB1)[3]<-c("Distal_intergenic")

#distBcat<-rbind(distBcat,distRCB1)

distLef<-rbind(distDL1,distDL2,distDL3)

colnames(distLef)[1]<-c("Three_UTR")
colnames(distLef)[2]<-c("Five_UTR")
colnames(distLef)[3]<-c("Distal_intergenic")

# DO NOT INCLUDE RCT1
distTcf<-rbind(distDT1,distDT2,distDT3)

colnames(distTcf)[1]<-c("Three_UTR")
colnames(distTcf)[2]<-c("Five_UTR")
colnames(distTcf)[3]<-c("Distal_intergenic")
#distTcf$Five_UTR<-c(0,0,0)
#distTcf<-distTcf[,c(1,8,2,3,4,5,6,7)]

distHB<-rbind(distRHC1,distRHC2,distRHC,distRHL1,distRHL2,distRHL)
colnames(distHB)[1]<-c("Three_UTR")
colnames(distHB)[2]<-c("Five_UTR")
colnames(distHB)[3]<-c("Distal_intergenic")

distBcatDKO<-rbind(distRKB1,distRKB2,distRKB3)

distTcfDKO<-rbind(distRKT1,distRKT2,distRKT3)
                   

distBcatSC<-rbind(distRBCS1,distRBCS2,
                  distRBLS1,distRBLS2)

colnames(distBcatSC)[1]<-"Three_UTR"
colnames(distBcatSC)[2]<-"Five_UTR"
colnames(distBcatSC)[3]<-"Distal_intergenic"

distchip<-rbind(distBcat,distBcatSC,distBcatDKO,distTcf,distTcfDKO,distLef,distHB)

distall<- data.frame(sapply(distchip[,1:7], function(x) as.numeric(as.character(x))*100))
distall<-cbind(distchip$Sample,distall)
colnames(distall)[1]<-"Sample"

distall$condition<-c("Bcat_control","Bcat_control","Bcat_control","Bcat_control","Bcat_control",
                     "Bcat_lit","Bcat_lit","Bcat_lit",
                     "Bcat_control_SC","Bcat_control_SC",
                     "Bcat_lit_SC","Bcat_lit_SC",
                     "Bcat_DKO","Bcat_DKO","Bcat_DKO",
                     "Tcf","Tcf","Tcf",
                     "Tcf_DKO","Tcf_DKO","Tcf_DKO",
                     "Lef1","Lef1","Lef1",
                     "HC","HC","HC","HL","HL","HL")

distall_melt<-melt(distall,id.vars = c("condition","Sample"))

pdf("/Users/yolanda_guillen/Desktop/IMIM/bcat_Project/plots_chips/dist_peaks.pdf")

#not KOs
distpeaks<-ggplot(distall_melt[distall_melt$condition!="Tcf_DKO" & distall_melt$condition!="Bcat_DKO",],aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat="identity",color="black") +
  facet_grid(~condition,scale="free",space="free")+
  scale_fill_brewer(palette="Pastel1")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")

print(distpeaks)


dev.off()

# Distribution according to number fo sequences aligned?? NEED REDO ALIGNMENT ON OF THE SAMPLES, LACKING %ALIGNMENT DATA
#distall$seqalign<-c(46.99,36.62,54.72,48.57,46.30,8.16,35.14,9.08,37.59,46.67,46.67,54.58,63.43,59.78,52.25,65.38,60,33.77,40.89,35,47.63,95.51,90)

#distall_melt<-melt(distall,id.vars = c("condition","Sample","seqalign"))

#ggplot(distall_melt,aes(y=seqalign,x=value))+
#  geom_point(size=4,color="black")+
#  geom_point(aes(color=condition),size=3) +
#  geom_smooth(method=lm,linetype="dashed",
#              color="darkred", fill="blue",alpha=0.1)+
#  facet_wrap(~variable,scale="free")+
  #scale_fill_brewer(palette="Pastel1")+
#  theme_light()+
#  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
#        axis.text.y = element_text(size=12),
#        axis.title.y = element_text(size=14,face = "bold"))+
#  ylab("% seqs aligned")+xlab("Fraction of peaks")


#Correlation?
#intron<-c((cor.test(distall$Intron,distall$seqalign))$p.value,(cor.test(distall$Intron,distall$seqalign))$estimate)
#exon<-c((cor.test(distall$Exon,distall$seqalign))$p.value,(cor.test(distall$Exon,distall$seqalign))$estimate)
#promoter<-c((cor.test(distall$Promoter,distall$seqalign))$p.value,(cor.test(distall$Promoter,distall$seqalign))$estimate)
#inter<-c((cor.test(distall$Distal_intergenic,distall$seqalign))$p.value,(cor.test(distall$Distal_intergenic,distall$seqalign))$estimate)
#downs<-c((cor.test(distall$Downstream,distall$seqalign))$p.value,(cor.test(distall$Downstream,distall$seqalign))$estimate)
#five<-c((cor.test(distall$Five_UTR,distall$seqalign))$p.value,(cor.test(distall$Five_UTR,distall$seqalign))$estimate)
#three<-c((cor.test(distall$Three_UTR,distall$seqalign))$p.value,(cor.test(distall$Three_UTR,distall$seqalign))$estimate)

#cortab<-as.data.frame(rbind(three,five,inter,downs,exon,intron,promoter))
#colnames(cortab)<-c("pvalue","cor")

#cortab

## Number of peaks
BRC1peaks<-peakAnnoBRC1@peakNum
BRC2peaks<-peakAnnoBRC2@peakNum
DB1peaks<-peakAnnoDB1@peakNum
DB2peaks<-peakAnnoDB2@peakNum
DB3peaks<-peakAnnoDB3@peakNum
RCB1peaks<-peakAnnoRCB1@peakNum

RBCS1peaks<-peakAnnoRBCS1@peakNum
RBCS2peaks<-peakAnnoRBCS2@peakNum
RBLS1peaks<-peakAnnoRBLS1@peakNum
RBLS2peaks<-peakAnnoRBLS2@peakNum

BRL2peaks<-peakAnnoBRL2@peakNum
BRL1peaks<-peakAnnoBRL1@peakNum
RBL1peaks<-peakAnnoRBL1@peakNum
RBL2peaks<-peakAnnoRBL2@peakNum
RBLpeaks<-peakAnnoRBL@peakNum
#RBL2_mergpeaks<-nrow(dataf_peakannoRBL_merg)

RHC1peaks<-peakAnnoRHC1@peakNum
RHC2peaks<-peakAnnoRHC2@peakNum
RHCpeaks<-peakAnnoRHC@peakNum

RHL1peaks<-peakAnnoRHL1@peakNum
RHL2peaks<-peakAnnoRHL2@peakNum
RHLpeaks<-peakAnnoRHL@peakNum


DL1peaks<-peakAnnoDL1@peakNum
DL2peaks<-peakAnnoDL2@peakNum
DL3peaks<-peakAnnoDL3@peakNum

DT1peaks<-peakAnnoDT1@peakNum
DT2peaks<-peakAnnoDT2@peakNum
DT3peaks<-peakAnnoDT3@peakNum
RCT1peaks<-peakAnnoRCT1@peakNum

RKB1peaks<-peakAnnoRKB1@peakNum
RKB2peaks<-peakAnnoRKB2@peakNum
RKB3peaks<-peakAnnoRKB3@peakNum

RKT1peaks<-peakAnnoRKT1@peakNum
RKT2peaks<-peakAnnoRKT2@peakNum
RKT3peaks<-peakAnnoRKT3@peakNum

peaknumb<-as.data.frame(rbind(BRC1peaks,BRC2peaks,
                              DB1peaks,DB2peaks,DB3peaks,RCB1peaks,
                              RBCS1peaks,RBCS2peaks,
                              RBLS1peaks,RBLS2peaks,
                              BRL1peaks,BRL2peaks,
                              RBL1peaks,RBL2peaks,RBLpeaks,
                              RKB1peaks,RKB2peaks,RKB3peaks,
                              DL1peaks,DL2peaks,DL3peaks,
                              DT1peaks,DT2peaks,DT3peaks,RCT1peaks,
                              RKT1peaks,RKT2peaks,RKT3peaks,
                              RHC1peaks,RHC2peaks,
                              RHL1peaks,RHL2peaks))

colnames(peaknumb)<-"num_peaks"

peaknumb


## Histogram peaks Score
BRC1<-subset(dataf_peakannoBRC1,select=c(FE))
BRC1$sample<-"BRC1"
BRC1$type<-"bcat_control"

BRC2<-subset(dataf_peakannoBRC2,select=c(FE))
BRC2$sample<-"BRC2"
BRC2$type<-"bcat_control"

DB1<-subset(dataf_peakannoDB1,select=c(FE))
DB1$sample<-"DB1"
DB1$type<-"bcat_control"

DB2<-subset(dataf_peakannoDB2,select=c(FE))
DB2$sample<-"DB2"
DB2$type<-"bcat_control"

DB3<-subset(dataf_peakannoDB3,select=c(FE))
DB3$sample<-"DB3"
DB3$type<-"bcat_control"

RCB1<-subset(dataf_peakannoRCB1,select=c(FE))
RCB1$sample<-"RCB1"
RCB1$type<-"bcat_control"

RKB1<-subset(dataf_peakannoRKB1,select=c(FE))
RKB1$sample<-"RKB1"
RKB1$type<-"bcat_DKO"

RKB2<-subset(dataf_peakannoRKB2,select=c(FE))
RKB2$sample<-"RKB2"
RKB2$type<-"bcat_DKO"

RKB3<-subset(dataf_peakannoRKB3,select=c(FE))
RKB3$sample<-"RKB3"
RKB3$type<-"bcat_DKO"

BRL2<-subset(dataf_peakannoBRL2,select=c(FE))
BRL2$sample<-"BRL2"
BRL2$type<-"bcat_lithium"

BRL1<-subset(dataf_peakannoBRL1,select=c(FE))
BRL1$sample<-"BRL1"
BRL1$type<-"bcat_lithium"

RBL1<-subset(dataf_peakannoRBL1,select=c(FE))
RBL1$sample<-"RBL1"
RBL1$type<-"bcat_lithium"

RBL2<-subset(dataf_peakannoRBL2,select=c(FE))
RBL2$sample<-"RBL2"
RBL2$type<-"bcat_lithium"

RBL<-subset(dataf_peakannoRBL,select=c(FE))
RBL$sample<-"RBL"
RBL$type<-"bcat_lithium"

#RBL_merge<-subset(dataf_peakannoRBL_merg,select=c(V5))
#RBL_merge$sample<-"RBL_merge"
#RBL_merge$type<-"lithium"

RBCS1<-subset(dataf_peakannoRBCS1,select=c(FE))
RBCS1$sample<-"RBCS1"
RBCS1$type<-"bcat_control_SC"

RBCS2<-subset(dataf_peakannoRBCS2,select=c(FE))
RBCS2$sample<-"RBCS2"
RBCS2$type<-"bcat_control_SC"

RBLS1<-subset(dataf_peakannoRBLS1,select=c(FE))
RBLS1$sample<-"RBLS1"
RBLS1$type<-"bcat_lithium_SC"

RBLS2<-subset(dataf_peakannoRBLS2,select=c(FE))
RBLS2$sample<-"RBLS2"
RBLS2$type<-"bcat_lithium_SC"

RHC1<-subset(dataf_peakannoRHC1,select=c(FE))
RHC1$sample<-"RHC1"
RHC1$type<-"Histone_control"

RHC2<-subset(dataf_peakannoRHC2,select=c(FE))
RHC2$sample<-"RHC2"
RHC2$type<-"Histone_control"

RHL1<-subset(dataf_peakannoRHL1,select=c(FE))
RHL1$sample<-"RHL1"
RHL1$type<-"Histone_Lithium"

RHL2<-subset(dataf_peakannoRHL2,select=c(FE))
RHL2$sample<-"RHL2"
RHL2$type<-"Histone_Lithium"

DL1<-subset(dataf_peakannoDL1,select=c(FE))
DL1$sample<-"DL1"
DL1$type<-"lef_control"

DL2<-subset(dataf_peakannoDL2,select=c(FE))
DL2$sample<-"DL2"
DL2$type<-"lef_control"

DL3<-subset(dataf_peakannoDL3,select=c(FE))
DL3$sample<-"DL3"
DL3$type<-"lef_control"

DT1<-subset(dataf_peakannoDT1,select=c(FE))
DT1$sample<-"DT1"
DT1$type<-"tcf_control"

DT2<-subset(dataf_peakannoDT2,select=c(FE))
DT2$sample<-"DT2"
DT2$type<-"tcf_control"

DT3<-subset(dataf_peakannoDT3,select=c(FE))
DT3$sample<-"DT3"
DT3$type<-"tcf_control"

RCT1<-subset(dataf_peakannoRCT1,select=c(FE))
RCT1$sample<-"RCT1"
RCT1$type<-"tcf_control"

RKT1<-subset(dataf_peakannoRKT1,select=c(FE))
RKT1$sample<-"RKT1"
RKT1$type<-"tcf_DKO"

RKT2<-subset(dataf_peakannoRKT2,select=c(FE))
RKT2$sample<-"RKT2"
RKT2$type<-"tcf_DKO"

RKT3<-subset(dataf_peakannoRKT3,select=c(FE))
RKT3$sample<-"RKT3"
RKT3$type<-"tcf_DKO"

#Do not include RBL1
dataall<-rbind(BRC1,BRC2,DB1,DB2,DB3,RCB1,
               BRL2,RBL2,RBL,
               RBCS1,RBCS2,
               RBLS1,RBLS2,
               RKB1,RKB2,RKB3,
               DT1,DT2,DT3,RCT1,
               RKT1,RKT2,RKT3,
               RHC1,RHC2,RHL1,RHL2)

library(plyr)
mu <- ddply(dataall, "sample", summarise, grp.median=median(FE))
mu$type<-c("bcat_control","bcat_control",
           "bcat_lithium",
           "bcat_control","bcat_control","bcat_control",
           "tcf_control","tcf_control","tcf_control",
           "bcat_control_SC","bcat_control_SC",
           "bcat_lithium","bcat_lithium",
           "bcat_lithium_SC","bcat_lithium_SC",
           "bcat_control",
           "tcf_control",
           "histcont","histcont","histlit","histlit",
           "bcat_DKO","bcat_DKO","bcat_DKO",
           "tcf_DKO","tcf_DKO","tcf_DKO")

p1<-ggplot(dataall,aes(x=FE))+
  #geom_histogram(aes(y=..density..), color="black",alpha=0.1, 
  #               position="identity")+
  geom_density(aes(color=sample),size=1,alpha=0.1)+
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("black","black",
                              "grey",
                              "black","black","black",
                              "black","black","black",
                              "black","black",
                              "grey","grey",
                              "grey","grey",
                              "black",
                              "black",
                              "darkolivegreen","darkolivegreen",
                              "dodgerblue","dodgerblue",
                              "purple","purple","purple",
                              "coral","coral","coral"))+
  #geom_vline(data=mu, aes(xintercept=grp.median, color=sample),
           #  linetype="dashed",size=1)+
  xlim(0,25)+
  xlab("Peak FE")+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"),legend.position = "bottom")

p1



### Distribution score peaks to detect contamination
# DKO bcat

dataf_peakannoRKB1$sample<-"RKB1"
dataf_peakannoRKB2$sample<-"RKB2"
dataf_peakannoRKB3$sample<-"RKB3"

dataf_peakannoDB$sample<-"DB"
dataf_peakannoRBL$sample<-"RBL"
dataf_peakannoRCB1$sample<-"RCB1"

dataf_peakannoRBCS1$sample<-"RBCS1"
dataf_peakannoRBCS2$sample<-"RBCS2"
dataf_peakannoRBLS1$sample<-"RBLS1"
dataf_peakannoRBLS2$sample<-"RBLS2"

DKOpeaks<-rbind(dataf_peakannoDB,dataf_peakannoRBL,dataf_peakannoRCB1,
                dataf_peakannoRBCS1,dataf_peakannoRBCS2,
                dataf_peakannoRBLS1,dataf_peakannoRBLS2,
                dataf_peakannoRKB1,dataf_peakannoRKB2,dataf_peakannoRKB3)

ggplot(DKOpeaks,aes(x=sample,y=FE,label=SYMBOL))+
  geom_boxplot()+
  geom_point()+
  geom_label_repel(data=DKOpeaks[DKOpeaks$FE>60,],aes(label=SYMBOL), size=3,fontface = "italic",alpha=0.5)+
  facet_wrap(~sample,scales = "free_x",ncol=5)+
  theme_bw()


### Distribution score peaks to detect contamination
# DKO TCF

dataf_peakannoRKT1$sample<-"RKT1"
dataf_peakannoRKT2$sample<-"RKT2"
dataf_peakannoRKT3$sample<-"RKT3"

dataf_peakannoDT$sample<-"DT"
dataf_peakannoRCT1$sample<-"RCT1"


DKOtcfpeaks<-rbind(dataf_peakannoRKT1,dataf_peakannoRKT2,dataf_peakannoRKT3,
                dataf_peakannoDT,dataf_peakannoRCT1)

ggplot(DKOtcfpeaks,aes(x=sample,y=FE,label=SYMBOL))+
  geom_boxplot()+
  geom_point()+
  geom_label_repel(data=DKOtcfpeaks[DKOtcfpeaks$FE>50,],aes(label=SYMBOL), size=3,fontface = "italic",alpha=0.5)+
  facet_grid(~sample,scales = "free_x")+
  theme_bw()


BiocManager::install("ggpubr")
library(ggpubr)

#my_comparisons <- list( c("control", "lithium") )

my_comparisons <- list( c("bcat_control", "lithium"),c("Histone_control","Histone_Lithium"))


compare_means(FE ~ type,  data = dataall)

p2<-ggplot(dataall,aes(type,FE))+
  geom_boxplot(aes(color=type),size=2)+
  #scale_color_manual(values=c("grey","black","darkolivegreen","red"))+
  ylim(0,50)+
  stat_compare_means(method = "anova")+
  stat_compare_means(comparisons = my_comparisons)+
  #stat_compare_means(aes(label=..p.signif..),label.x=1.5,label.y=40,size=10)+
  stat_compare_means(label.y=40)+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "bottom")+
  ylab("Peak score")+xlab("")


grid.arrange(p1,p2,ncol=2)

wilcox.test(dataall$FE~dataall$type)

## Venn diagram peaks lithium vs. control bcat and TCF and LEF
# FOR GENES, NOT PEAKS

#bcat Control
BRC1_genes<-unique(dataf_peakannoBRC1$SYMBOL)
BRC2_genes<-unique(dataf_peakannoBRC2$SYMBOL)
DB1_genes<-unique(dataf_peakannoDB1$SYMBOL)
DB2_genes<-unique(dataf_peakannoDB2$SYMBOL)
DB3_genes<-unique(dataf_peakannoDB3$SYMBOL)
RCB1_genes<-unique(dataf_peakannoRCB1$SYMBOL)

bcatgenes<-c(BRC1_genes,BRC2_genes,DB1_genes,DB2_genes,DB3_genes,RCB1_genes)
bcatcomgenes<-as.data.frame(table(sort(bcatgenes)))
bcatcomgenes<-bcatcomgenes[order(bcatcomgenes$Freq,decreasing = TRUE),]

#Lithium
#BRL1_genes<-unique(dataf_peakannoBRL1$SYMBOL)
BRL2_genes<-unique(dataf_peakannoBRL2$SYMBOL)
#RBL1_genes<-unique(dataf_peakannoRBL1$SYMBOL)
RBL2_genes<-unique(dataf_peakannoRBL2$SYMBOL)
RBL_genes<-unique(dataf_peakannoRBL$SYMBOL)
#RBL_merge_genes<-unique(dataf_peakannoRBL_merg$SYMBOL)

bcatlit<-c(BRL2_genes,RBL2_genes,RBL_genes)
bcatlitcomgenes<-as.data.frame(table(sort(bcatlit)))
bcatlitcomgenes<-bcatlitcomgenes[order(bcatlitcomgenes$Freq,decreasing = TRUE),]

## BCAT control SANTA CRUZ
RBCS1_genes<-unique(dataf_peakannoRBCS1$SYMBOL)
RBCS2_genes<-unique(dataf_peakannoRBCS2$SYMBOL)

bcatgenesSC<-c(RBCS1_genes,RBCS2_genes)
bcatSCcomgenes<-as.data.frame(table(sort(bcatgenesSC)))
bcatSCcomgenes<-bcatSCcomgenes[order(bcatSCcomgenes$Freq,decreasing = TRUE),]

## BCAT lithium SANTA CRUZ
RBLS1_genes<-unique(dataf_peakannoRBLS1$SYMBOL)
RBLS2_genes<-unique(dataf_peakannoRBLS2$SYMBOL)

bcatlitgenesSC<-c(RBLS1_genes,RBLS2_genes)
bcatSClitcomgenes<-as.data.frame(table(sort(bcatlitgenesSC)))
bcatSClitcomgenes<-bcatSClitcomgenes[order(bcatSClitcomgenes$Freq,decreasing = TRUE),]


# DKO bcat
RKB1_genes<-unique(dataf_peakannoRKB1$SYMBOL)
RKB2_genes<-unique(dataf_peakannoRKB2$SYMBOL)
RKB3_genes<-unique(dataf_peakannoRKB3$SYMBOL)

bcatDKO<-c(RKB1_genes,RKB2_genes,RKB3_genes)
bcatDKOcomgenes<-as.data.frame(table(sort(bcatDKO)))
bcatDKOcomgenes<-bcatDKOcomgenes[order(bcatDKOcomgenes$Freq,decreasing = TRUE),]

# Hist control genes
RHC1_genes<-unique(dataf_peakannoRHC1$SYMBOL)
RHC2_genes<-unique(dataf_peakannoRHC2$SYMBOL)
RHC_genes<-unique(dataf_peakannoRHC$SYMBOLS)

RHCcomgenes<-c(RHC1_genes,RHC2_genes,RHC_genes)
RHCcomgenes<-as.data.frame(table(sort(RHCcomgenes)))
RHCcomgenes<-RHCcomgenes[order(RHCcomgenes$Freq,decreasing = TRUE),]

# Hist lithium genes
RHL1_genes<-unique(dataf_peakannoRHL1$SYMBOL)
RHL2_genes<-unique(dataf_peakannoRHL2$SYMBOL)
RHL_genes<-unique(dataf_peakannoRHL$SYMBOLS)

RHLcomgenes<-c(RHL1_genes,RHL2_genes,RHL_genes)
RHLcomgenes<-as.data.frame(table(sort(RHLcomgenes)))
RHLcomgenes<-RHLcomgenes[order(RHLcomgenes$Freq,decreasing = TRUE),]




#TCF
DT1_genes<-unique(dataf_peakannoDT1$SYMBOL)
DT2_genes<-unique(dataf_peakannoDT2$SYMBOL)
DT3_genes<-unique(dataf_peakannoDT3$SYMBOL)
RCT1_genes<-unique(dataf_peakannoRCT1$SYMBOL)

tcfgenes<-c(DT1_genes,DT2_genes,DT2_genes,RCT1_genes)
tcfcomgenes<-as.data.frame(table(sort(tcfgenes)))
tcfcomgenes<-tcfcomgenes[order(tcfcomgenes$Freq,decreasing = TRUE),]


# DKO TCF
RKT1_genes<-unique(dataf_peakannoRKT1$SYMBOL)
RKT2_genes<-unique(dataf_peakannoRKT2$SYMBOL)
RKT3_genes<-unique(dataf_peakannoRKT3$SYMBOLS)

tcfDKO<-c(RKT1_genes,RKT2_genes,RKT3_genes)
tcfDKOcomgenes<-as.data.frame(table(sort(tcfDKO)))
tcfDKOcomgenes<-tcfDKOcomgenes[order(tcfDKOcomgenes$Freq,decreasing = TRUE),]

#LEF
DL1_genes<-unique(dataf_peakannoDL1$SYMBOL)
DL2_genes<-unique(dataf_peakannoDL2$SYMBOL)
DL3_genes<-unique(dataf_peakannoDL3$SYMBOL)

lefgenes<-c(DL1_genes,DL2_genes,DL2_genes)
lefcomgenes<-as.data.frame(table(sort(lefgenes)))
lefcomgenes<-lefcomgenes[order(lefcomgenes$Freq,decreasing = TRUE),]

#Tables to plot
tabgenesbcat<-as.data.frame(table(bcatcomgenes$Freq))
tabgenesbcat$type<-c("bcat")
tabgenesbcatlit<-as.data.frame(table(bcatlitcomgenes$Freq))
tabgenesbcatlit$type<-c("bcat_lit")
tabgenesbcatSC<-as.data.frame(table(bcatSCcomgenes$Freq))
tabgenesbcatSC$type<-c("bcat_SC")
tabgenesbcatSClit<-as.data.frame(table(bcatSClitcomgenes$Freq))
tabgenesbcatSClit$type<-c("bcat_lit_SC")
tabgenesbcatDKO<-as.data.frame(table(bcatDKOcomgenes$Freq))
tabgenesbcatDKO$type<-c("bcatDKO")
tabgenestcf<-as.data.frame(table(tcfcomgenes$Freq))
tabgenestcf$type<-c("tcf")
tabgenestcfDKO<-as.data.frame(table(tcfDKOcomgenes$Freq))
tabgenestcfDKO$type<-c("tcfDKO")
tabgeneslef<-as.data.frame(table(lefcomgenes$Freq))
tabgeneslef$type<-c("lef")
tabgenesRHC<-as.data.frame(table(RHCcomgenes$Freq))
tabgenesRHC$type<-c("RHC")
tabgenesRHL<-as.data.frame(table(RHLcomgenes$Freq))
tabgenesRHL$type<-c("RHL")

tabgenesreps<-rbind(tabgenesbcat,tabgenesbcatlit,tabgenesbcatSC,tabgenesbcatSClit,
                    tabgenesbcatDKO,
                    tabgenestcf,tabgenestcfDKO,tabgeneslef,tabgenesRHC,tabgenesRHL)

ggbarplot(tabgenesreps, y = "Freq", x = "Var1",
          fill="steelblue",alpha=0.4,
          facet.by = "type",scales="free",
          label = TRUE, lab.pos = "in", lab.col = "white",lab.size = 3,
          xlab = "Number of replicates", ylab="Number of genes")



# Taking genes that appear in at least two different replicates
#WHEN NEW REPLICATES OF BCAT LITHIUM, BETTER!

bcatcomgenes<-subset(bcatcomgenes,bcatcomgenes$Freq>=2)
#for lithium, only one replicate minimum
bcatlitcomgenes<-subset(bcatlitcomgenes,bcatlitcomgenes$Freq>=1)

# Control histone control and lithium
#geneLists<-list(Control = as.character(RHCcomgenes$Var1),Lithium = as.character(RHLcomgenes$Var1))

# Control bcat genes vs bcat genes SC
geneLists<-list(SC = as.character(bcatgenesSC),
                RD = as.character(bcatgenes))


# Replicates lithium
#geneLists<-list(Lithium1 = BRL2_genes,Lithium2 = RBL2_genes,Lithium3=RBL_genes,Lithium4=RBL_merge_genes)

geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off() 

# FOR THREE GROUPS
#venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "blue","darkolivegreen2", "pink"), alpha=c(0.3,0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("Lithium1","Lithium2","Lithium_m","Lithium_merge"), main="Genes ChIP bcat min 2 replicates")

# FOR TWO GROUPS
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "blue"), alpha=c(0.3,0.3), cex = 2, cat.fontface=4, category.names=c("SC","RD"), main="Genes ChIP bcat 2 Abs")


# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

# To get the list of gene present in each Venn compartment we can use the gplots package
require("gplots")

a <- venn(VENN.LIST, show.plot=FALSE)

inters <- attr(a,"intersections")

length(unique(inters[[1]]))
length(unique(inters[[2]]))
length(unique(inters[[3]]))

congenes<-inters[[2]]
commongenes<-inters[[1]]
litgenes<-inters[[3]]

#bcatexpeak<-inters[[3]]
#contexpeak<-inters[[2]]

#setdiff(bcatexpeak,litgenes)


comDB1<-dataf_peakannoDB1[with(dataf_peakannoDB1, dataf_peakannoDB1$SYMBOL %in% commongenes),]
comDB1<-comDB1[!duplicated(comDB1$SYMBOL), ] 
comDB2<-dataf_peakannoDB2[with(dataf_peakannoDB2, dataf_peakannoDB2$SYMBOL %in% commongenes),]
comDB2<-comDB2[!duplicated(comDB2$SYMBOL), ] 
comDB3<-dataf_peakannoDB3[with(dataf_peakannoDB3, dataf_peakannoDB3$SYMBOL %in% commongenes),]
comDB3<-comDB3[!duplicated(comDB3$SYMBOL), ] 
comBRC1<-dataf_peakannoBRC1[with(dataf_peakannoBRC1, dataf_peakannoBRC1$SYMBOL %in% commongenes),]
comBRC1<-comBRC1[!duplicated(comBRC1$SYMBOL), ] 
comBRC2<-dataf_peakannoBRC2[with(dataf_peakannoBRC2, dataf_peakannoBRC2$SYMBOL %in% commongenes),]
comBRC2<-comBRC2[!duplicated(comBRC2$SYMBOL), ] 
dataf_peakannoRCB1$sample<-NULL
comRCB1<-dataf_peakannoRCB1[with(dataf_peakannoRCB1, dataf_peakannoRCB1$SYMBOL %in% commongenes),]
comRCB1<-comRCB1[!duplicated(comRCB1$SYMBOL), ] 

comBRL2<-dataf_peakannoBRL2[with(dataf_peakannoBRL2, dataf_peakannoBRL2$SYMBOL %in% commongenes),]
comBRL2<-comBRL2[!duplicated(comBRL2$SYMBOL), ] 

comRBL2<-dataf_peakannoRBL2[with(dataf_peakannoRBL2, dataf_peakannoRBL2$SYMBOL %in% commongenes),]
comRBL2<-comRBL2[!duplicated(comRBL2$SYMBOL), ] 

#comRBL_merg<-dataf_peakannoRBL_merg[with(dataf_peakannoRBL_merg, dataf_peakannoRBL_merg$SYMBOL %in% commongenes),]
#comRBL_merg<-comRBL_merg[!duplicated(comRBL_merg$SYMBOL), ]

#Subset scores per replicate
comball<-rbind(comBRL2,comDB1,comDB2,comDB3,comBRC1,comBRC2,comRBL2,comRCB1)
comball$name<-gsub('/.*/','',comball$name)
comball$name<-gsub('_peaks.*','',comball$name)

#subset from dataframe
comball<-subset(comball,select=c(FE,name,SYMBOL))

#rescale score values
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
comball$rescale<-range01(comball$FE)

#variable lithium vs control
comball$samp<-ifelse(comball$name== "BRL2" | comball$name== "RBL_merge" | comball$name== "RBL2", 
                        c("Lithium"), c("Control")) 

comball$SYMBOL  = factor(comball$SYMBOL , levels=unique(comball$SYMBOL [order(comball$FE, comball$name)]), ordered=TRUE)


levels(comball$SYMBOL)


ggplot(comball, aes(x=name, y=SYMBOL)) + 
  geom_tile(aes(fill = rescale),colour = "white")+
  scale_fill_gradient2(low = "#FFFFFF",
                       high = "red",
                       name="Peak Score") +
  #scale_fill_gradient(low = "blue", high = "red")+
  facet_grid(~ samp, switch = "x", scales = "free_x", space = "free_x")+
  theme_dark()+
  theme(axis.text.x = element_text(size=6, hjust = 1),
        axis.text.y = element_text(size=5, hjust = 1,face="italic"))+
  ggtitle("Common peaks Control and Lithium")


### heatmap chip with bcat peaks
### Bed file gene annotation, selecting from bed genes with one peak in at least two replicates for each condition
bedhuman<-read.delim("/Volumes/cancer/db_files/Human_GRCh38/Human_GRCh38_genes.bed",header = FALSE)
colnames(bedhuman)<-c("Chr","Start","End","Name","Source","Strand")

#all genes
#for bcat peaks
selcat<-unique(c(bcatgenes,bcatlit))
selcat<-bedhuman[with(bedhuman, bedhuman$Name %in% selcat),]
selcat$Chr<-gsub('^','chr',selcat$Chr)
 
write.table(selcat,"/Volumes/cancer/bcat_Project/peakcall/browser_files/genes_peaks_bcatall.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

#Selecting genes with no peak in bcat, but peak in histones controls and lithium
selhist<-unique(c(RHC1_genes,RHC2_genes,RHL1_genes,RHL2_genes))
selhist<-bedhuman[(bedhuman$Name %in% selhist),]

write.table(selhist,"/Volumes/cancer/bcat_Project/peakcall/browser_files/genes_histpeaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")



#Functional enrichment
## GO terms and pathways
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018","InterPro_Domains_2019",
         "Pfam_Domains_2019","Chromosome_Location","GO_Molecular_Function_2018",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X","WikiPathways_2016")

#Separately for each venn category (histones!!!)

enrichedcom <- enrichr(commongenes, dbs)
enrichedlit <- enrichr(litgenes,dbs)
enrichedcont <- enrichr(congenes,dbs)

# DKO bcat
enrichedDKObcat <- enrichr(as.character(bcatDKOcomgenes$Var1), dbs)

Unrich <- enrichedcom[["GO_Biological_Process_2018"]]
Drich <- enrichedlit[["GO_Biological_Process_2018"]]
Crich <- enrichedcont[["GO_Biological_Process_2018"]]

DKObcat <- enrichedDKObcat[["GO_Biological_Process_2018"]]

Unrich$group<-c("Common")
Drich$group<-c("Lithium")
Crich$group<-c("Control")
DKObcat$group<-c("DKObcat")


# exclusive lithium (see below)
exc_lit$Var1
enrichedexclit <- enrichr(as.character(exc_lit$Var1),dbs)
exclitrich <- enrichedexclit[["GO_Biological_Process_2018"]]
exclitrich$group<-c("exclusive_lit")

allGO<-rbind(Unrich,Drich,Crich,DKObcat,exclitrich)



#All targets, lithium control and common
bcat<-unique(c(commongenes,litgenes,congenes))
enrichedcomlit<-enrichr(unique(bcat),dbs)

COMLITrich<-enrichedcomlit[["GO_Biological_Process_2018"]]

COMLITrich$group<-"comlitcont"
allGO<-COMLITrich

##################

bpsub<-subset(allGO,allGO$Adjusted.P.value<=0.5)
bpsub<- bpsub[order(bpsub$P.value),]

#for negative and positive values, pyramid, two groups comparison
#bpsub$sig<-with(bpsub, ifelse(group == "Common", log(P.value), (log(P.value))*(-1)))


bpsub$Term<-as.factor(bpsub$Term)

bpsub$Term <- factor(bpsub$Term, levels= unique(bpsub$Term)[order(bpsub$group,bpsub$Combined.Score)])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <- factor(bpsub$Term, levels= unique(bpsub$Term)[order(bpsub$group,bpsub$Combined.Score)])

ggplot(bpsub,aes(y=Combined.Score,x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group),color="black")+
  #geom_text(aes(label=tolower(Genes)), size=4,hjust=1,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Blues")+
  #facet_wrap(~group,scales="free")+
  theme_minimal()+
  coord_flip()+
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=8, face="bold",hjust = 1),
        axis.text.y = element_text(size=10, hjust = 1),
        legend.position = "bottom")+
  facet_grid(~group,scales="free",space = "free")

RNA_bind<-(bpsub[1,][9])
write.table(RNA_bind,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/gene_lists/RNA_binding_bcat_peaks.txt",quote = FALSE,row.names = FALSE)

mRNA_bind<-(bpsub[3,][9])
write.table(mRNA_bind,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/gene_lists/mRNA_binding_bcat_peaks.txt",quote = FALSE,row.names = FALSE)


## GO ANNOTATION
library(GO.db)
library(topGO)

# genes annotated in chips
#latest version
hg<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")


# common genes bcat chips, at least in two replicates
bcatcomgenes$Var1
bcatlitcomgenes$Var1
rbind(bcatcomgenes,bcatlitcomgenes)

bmbcat <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
            values= (rbind(bcatcomgenes,bcatlitcomgenes))$Var1,
            filter="external_gene_name",
            mart=hg) 
getBM(attributes=c("ensembl_gene_id", "external_gene_name","mmusculus_homolog_ensembl_gene","mmusculus_homolog_associated_gene_name"),
      values= (rbind(bcatcomgenes,bcatlitcomgenes))$Var1,
      filter="external_gene_name",
      mart=hg)

# genes exclusive bcat lithium
exc_lit<-bcatlitcomgenes[bcatlitcomgenes$Var1 %in% setdiff(bcatlitcomgenes$Var1,bcatcomgenes$Var1),]

bmbcatlit <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
            values= as.character(exc_lit$Var1),
            filter="external_gene_name",
            mart=hg) 

getBM(attributes=c("ensembl_gene_id", "external_gene_name","mmusculus_homolog_ensembl_gene","mmusculus_homolog_associated_gene_name"),
                   values= as.character(exc_lit$Var1),
                   filter="external_gene_name",
                   mart=hg)

# GO TERM table all genes
vas<-toTable(GOTERM)
colnames(vas)[2]<-"repgo_id"
vascc<-vas[vas$Ontology=="CC",]
vasbp<-vas[vas$Ontology=="BP",]
vasmf<-vas[vas$Ontology=="MF",]

# peaks bcat control and bcat lithium or exclusive bcat lithium
geneGObcat<-merge(rbind(bcatcomgenes,bcatlitcomgenes),bmbcat,by.x="Var1",by.y="external_gene_name",all.x=TRUE)
geneGObcatlit<-merge(exc_lit,bmbcatlit,by.x="Var1",by.y="external_gene_name",all.x=TRUE)


#empty values as NA
geneGObcat$go_id<-gsub("^$", NA,geneGObcat$go_id)
geneGObcatlit$go_id<-gsub("^$", NA,geneGObcatlit$go_id)

# all terms
#geneGO<-merge(geneGO,vas,by="go_id",all.x=TRUE)

#only CC term
geneGOCCbcat<-merge(geneGObcat,vascc,by="go_id",all.x=TRUE)
geneGOBPbcat<-merge(geneGObcat,vasbp,by="go_id",all.x=TRUE)
geneGOMFbcat<-merge(geneGObcat,vasmf,by="go_id",all.x=TRUE)

geneGOCCbcatlit<-merge(geneGObcatlit,vascc,by="go_id",all.x=TRUE)
geneGOBPbcatlit<-merge(geneGObcatlit,vasbp,by="go_id",all.x=TRUE)
geneGOMFbcatlit<-merge(geneGObcatlit,vasmf,by="go_id",all.x=TRUE)

sort(unique((geneGOBPbcat[grepl('splicing',geneGOBPbcat$Term) | grepl('RNA processing',geneGOBPbcat$Term),])$Var1))
sort(unique((geneGOBPbcatlit[grepl('splicing',geneGOBPbcatlit$Term) | grepl('RNA processing',geneGOBPbcatlit$Term),])$Var1))


### RNASeq shbcat vs. control Jurkat and RPMI

library(DESeq2)
library(biomaRt)
library(geneplotter)
#library(MDplot)
library(lattice)
library(genefilter)
library(plotly)
library(limma)
library(ggrepel)

#Set wd
#setwd("/Volumes/cancer/Gekas_RNAseq/htseq/")
setwd("/Users/yguillen/Desktop/temp/Gekas_RNASeq/htseq/")

#output.dir="/Volumes/cancer/Gekas_RNAseq/htseq/"
output.dir="/Users/yguillen/Desktop/temp/Gekas_RNASeq/htseq/"

# Import metadata files
metadatar<-read.delim("metadata.txt",sep=" ",header = FALSE)
colnames(metadatar)<-c("sampleID","countFile","line","exp","condition")

#For jurkat sh
metadata_j<-subset(metadatar,metadatar$condition=="sh" & metadatar$line=="J")

#For RPMI sh
metadata_r<-subset(metadatar,metadatar$condition=="sh" & metadatar$line=="R")


# For sh experiments, both jurkat and rpmi
metadata_sh<-subset(metadatar,metadatar$condition=="sh")

library(DESeq2)

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata_r,
                                          directory = output.dir,
                                          design = ~ exp)


rowData(DESeq2Table)
mcols(rowData(DESeq2Table))

# Add additional annotation information using biomaRt
ensembl76 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
listAttributes(ensembl76)

bm_at <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","hgnc_symbol",
                         "gene_biotype","chromosome_name","start_position"),
            values= rownames(DESeq2Table),
            #            filter="external_gene_name",
            mart=ensembl76) 


head(bm_at)

#Add description data to gene counts
DESeq2Features <- data.frame(id = rownames(DESeq2Table))
DESeq2Features$id <- as.character(DESeq2Features$id)


### join them together
#if htseq ensemblid
rowData <- merge(DESeq2Features, bm_at, by.x="id",by.y = "ensembl_gene_id")

#if htseqgene symbol
#rowData <- merge(DESeq2Features, bm, by.x="id",by.y = "external_gene_name")

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
colData(DESeq2Table)$exp <- factor(colData(DESeq2Table)$exp, levels = c("Control","shbcat"))

#### estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

# plot densities of counts for the different samples to assess their distributions
multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))


### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Table, blind=TRUE) ## create a distance matrix between the samples

## Create PCA
DESeq2::plotPCA(rld, intgroup=c("exp"))
dev.off()


## Create PCA 
pcaData <- plotPCA(rld, intgroup=c("line","condition","exp"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggrepel)

ggplot(pcaData, aes(PC1, PC2, color=exp,label=name)) +
  geom_point(size=9,aes(shape=condition),color="black")+
  geom_point(size=8,aes(shape=condition)) +
  stat_ellipse(aes(color=exp),level=0.8)+
  #scale_color_manual(values=c("dodgerblue","coral"))+
  #scale_color_manual(values=c("grey","dodgerblue","coral","darkolivegreen2","black","dodgerblue4","coral4"))+
  geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.5)+
  #geom_text(aes(label=name),hjust=0.5, vjust=1,size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #facet_grid(~gender)+
  #scale_color_manual(values=c("dodgerblue2","black","black","cadetblue1"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed()


ggplot(pcaData, aes(y=PC1, x=exp,label=name)) +
  geom_point(aes(color=exp),size=8) +
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



### Differential EXPRESSION ANALYSIS ###
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

# Statistical testing of DE genes
levels(DESeq2Table$exp)

design(DESeq2Table)

dds<-DESeq(DESeq2Table)
resultsNames(dds)


#results DEGs between condition
rescondition<-results(dds,alpha=0.1)


#Example differences genes between groups
plotCounts(dds, gene=which(rownames(rescondition)=="ENSG00000168214"), intgroup="exp")

head(rescondition)

summary(rescondition)

#genes with < 0.05
sum(rescondition$padj < 0.1,na.rm=TRUE)

plotMA(rescondition)

#Heatmap most differentially expressed genes (only ranked by DEGs)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )

library(gplots)
library(RColorBrewer)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", dendrogram="column",trace = "none", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="grey", shbcat="black" )[
             colData(rld)$exp ],cexCol=0.3)

# Selecting specific DEGs
res_df_rna<-as.data.frame(rescondition)
class(rowData)
rowData<-as.data.frame(rowData)
res_df_rna<-merge(rowData,res_df_rna,by.x="id",by.y="row.names",all.y=TRUE)
res_df_sig_rna<-subset(res_df_rna,res_df_rna$padj<=0.1)
res_df_sig_rna <- res_df_sig_rna[order(res_df_sig_rna$pvalue),]



# UP and DOWN genes
UPgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange>0))$external_gene_name)
DOWNgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange<0))$external_gene_name)
length(UPgenes)
length(DOWNgenes)


#Manual heatmap from log trasformed data, using all samples
rld_df<-as.data.frame(assay(rld))


#only degs
deggenes<-unique(res_df_sig_rna$id)

# for condition
#all genes
deggenes<-unique(res_df_rna$id)


deg_rld<-which(rownames(rld_df) %in% deggenes)
deg_rld_df<-rld_df[deg_rld,]

deg_rld<-merge(deg_rld_df,res_df_rna,by.x="row.names",by.y="id")

deg_rld<- deg_rld[order(deg_rld$log2FoldChange),]
deg_rld$Row.names

heatmap.2( as.matrix(deg_rld[,2:5]), scale="row", dendrogram="column",trace="none", Rowv=FALSE,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           cexCol=0.3)



DOWNgenes
UPgenes

DEGsh<-unique(c(DOWNgenes,UPgenes))


## VENN DIAGRAM RNASEQ CHIPSEQ
# Control bcat genes vs lithium bcat
geneLists<-list(CHIPbcat = as.character(unique((rbind(bcatcomgenes,bcatlitcomgenes))$Var1)),
                unique(lefgenes),
                unique(tcfgenes),DEGsh)

length(unique(lefgenes))
length(unique(tcfgenes))
length(as.character(unique((rbind(bcatcomgenes,bcatlitcomgenes))$Var1)))
length(DEGsh)

geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off() 

# FOR four GROUPS
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue","coral3","grey"), 
                          alpha=c(0.3,0.3,0.3,0.3), 
                          cex = 2, cat.fontface=4, 
                          category.names=c("bcat","LEF","TCF","DEG_sh_rpmi"), main="Genes ChIP bcat")

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
length(inters[[8]])
length(inters[[9]])
length(inters[[10]])
length(inters[[11]])
length(inters[[12]])
length(inters[[13]])

### Plot RNA-Seq vs. peak
#RPMI
res_df_rna_rpmi<-res_df_rna[!duplicated(res_df_rna$external_gene_name), ]
write.table(res_df_rna_rpmi,"/Users/yguillen/Desktop/temp/Gekas_RNASeq/DESeq2/RPMI/res_df_rna_rpmi.txt",row.names=FALSE,quote = FALSE,sep="\t")

#JURKAT
res_df_rna_jurkat<-res_df_rna[!duplicated(res_df_rna$external_gene_name), ]
write.table(res_df_rna_jurkat,"/Users/yguillen/Desktop/temp/Gekas_RNASeq/DESeq2/RPMI/res_df_rna_jurkat.txt",row.names=FALSE,quote = FALSE,sep="\t")




#peaks bcat genes, including lithium
bcat<-c(BRC1_genes,BRC2_genes,DB1_genes,DB2_genes,DB3_genes,BRL2_genes,RBL2_genes)
bcat<-as.data.frame(table(sort(bcat)))
bcat<-bcat[order(bcat$Freq,decreasing = TRUE),]

#peaks lef genes
unique(lefcomgenes$Var1)

#peaks tcf genes
unique(tcfcomgenes$Var1)

#bcatcomgenes --> only controls, no lithium
#bcat -> all genes lithium and controls

#JURKAT
rnapeak<-merge(res_df_rna_jurkat,bcat,by.x="external_gene_name",by.y="Var1",all=TRUE)
rnapeak<-merge(rnapeak,lefcomgenes,by.x="external_gene_name",by.y="Var1",all=TRUE)
rnapeak<-merge(rnapeak,tcfcomgenes,by.x="external_gene_name",by.y="Var1",all=TRUE)

rnapeak$bcat<-ifelse(is.na(rnapeak$Freq.x), 
                        c("nobcat"), c("bcat")) 
rnapeak$LEF<-ifelse(is.na(rnapeak$Freq.y), 
                      c("nolef"), c("lef"))
rnapeak$TCF<-ifelse(is.na(rnapeak$Freq), 
                    c("notcf"), c("tcf"))


rnapeak$sig<-ifelse(rnapeak$pvalue<0.05, 
                      c("sig"), c("nosig"))

rnapeak$peaks<-paste(rnapeak$bcat,rnapeak$LEF,rnapeak$TCF,sep="_")

table(rnapeak$sig,rnapeak$gene_biotype)

#rnapeak<-subset(rnapeak,select=c(external_gene_name,log2FoldChange,pvalue,padj,peaks,sig))
rnapeak<-rnapeak[order(rnapeak$external_gene_name),]

rnapeak_sig<-subset(rnapeak,rnapeak$sig=="sig")
length(unique(rnapeak_sig$external_gene_name))

ggplot(rnapeak,aes(x=peaks,y=log2FoldChange))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_jitter(aes(color=sig))+
  geom_violin(data = subset(rnapeak,rnapeak$sig=="sig"),color="red",alpha=0.2)+
  geom_violin(data = subset(rnapeak,rnapeak$sig!="sig"),color="grey",alpha=0.2)+
  scale_color_manual(values=c("grey","red"))+
  geom_label_repel(aes(label = external_gene_name),color="darkolivegreen", data = subset(rnapeak,rnapeak$bcat=="bcat" & rnapeak$sig=="sig"), fontface = 4,size=2)+
  geom_label_repel(aes(label = external_gene_name),color="darkolivegreen", data = subset(rnapeak,abs(rnapeak$log2FoldChange)>2 & rnapeak$sig=="sig"), fontface = "italic",size=2)+
  geom_label_repel(aes(label = external_gene_name),color="black", data = subset(rnapeak,rnapeak$external_gene_name=="DROSHA"), fontface = "italic",size=2)+
  facet_grid(~bcat,scales="free")+
  ylab("log2FC bcat vs nobcat")+
  theme_bw()+
  theme(legend.position="none")

#RPMI
rnapeak<-merge(res_df_rna_rpmi,bcat,by.x="external_gene_name",by.y="Var1",all=TRUE)
rnapeak<-merge(rnapeak,lefcomgenes,by.x="external_gene_name",by.y="Var1",all=TRUE)
rnapeak<-merge(rnapeak,tcfcomgenes,by.x="external_gene_name",by.y="Var1",all=TRUE)

rnapeak$bcat<-ifelse(is.na(rnapeak$Freq.x), 
                     c("nobcat"), c("bcat")) 
rnapeak$LEF<-ifelse(is.na(rnapeak$Freq.y), 
                    c("nolef"), c("lef"))
rnapeak$TCF<-ifelse(is.na(rnapeak$Freq), 
                    c("notcf"), c("tcf"))


rnapeak$sig<-ifelse(rnapeak$pvalue<0.05, 
                    c("sig"), c("nosig"))

rnapeak$peaks<-paste(rnapeak$bcat,rnapeak$LEF,rnapeak$TCF,sep="_")

table(rnapeak$sig,rnapeak$gene_biotype)

#rnapeak<-subset(rnapeak,select=c(external_gene_name,log2FoldChange,pvalue,padj,peaks,sig))
rnapeak<-rnapeak[order(rnapeak$external_gene_name),]

rnapeak_sig<-subset(rnapeak,rnapeak$sig=="sig")
length(unique(rnapeak_sig$external_gene_name))

ggplot(rnapeak,aes(x=peaks,y=log2FoldChange))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_jitter(aes(color=sig))+
  geom_violin(data = subset(rnapeak,rnapeak$sig=="sig"),color="red",alpha=0.2)+
  geom_violin(data = subset(rnapeak,rnapeak$sig!="sig"),color="grey",alpha=0.2)+
  scale_color_manual(values=c("grey","red"))+
  geom_label_repel(aes(label = external_gene_name),color="darkolivegreen", data = subset(rnapeak,rnapeak$bcat=="bcat" & rnapeak$sig=="sig"), fontface = 4,size=2)+
  geom_label_repel(aes(label = external_gene_name),color="darkolivegreen", data = subset(rnapeak,abs(rnapeak$log2FoldChange)>2 & rnapeak$sig=="sig"), fontface = "italic",size=2)+
  geom_label_repel(aes(label = external_gene_name),color="black", data = subset(rnapeak,rnapeak$external_gene_name=="DROSHA"), fontface = "italic",size=2)+
  facet_grid(~bcat,scales="free")+
  ylab("log2FC bcat vs nobcat")+
  theme_bw()+
  theme(legend.position="none")


## Annotation of differentially expressed genes (RNA binding proteins?)




#### RBP aifantis, compare RBP expression patterns
#list of RBPs from Aifantis
RBP_genes<-read.delim("/Volumes/grcmc/YGUILLEN/beta_catenin_project/gene_lists/RBP_Aifantis.txt",header = FALSE)
RBP_genes<-RBP_genes$V1

RBP_DEGs_rpmi<-res_df_rna_rpmi[res_df_rna_rpmi$hgnc_symbol %in% RBP_genes,]
RBP_DEGs_rpmi<-subset(RBP_DEGs_rpmi,RBP_DEGs_rpmi$pvalue<0.05)

RBP_DEGs_jurkat<-res_df_rna_jurkat[res_df_rna_jurkat$hgnc_symbol %in% RBP_genes,]
RBP_DEGs_jurkat<-subset(RBP_DEGs_jurkat,RBP_DEGs_jurkat$pvalue<0.05)





## List of RNA binding proteins, target bcat chips
RNAbind<-read.delim("/Volumes/grcmc/YGUILLEN/beta_catenin_project/gene_lists/RNA_binding_bcat_peaks.txt",header = TRUE)
RNAbind<-RNAbind$Genes

mRNAbind<-read.delim("/Volumes/grcmc/YGUILLEN/beta_catenin_project/gene_lists/mRNA_binding_bcat_peaks.txt",header = TRUE)
mRNAbind<-mRNAbind$Genes

RNA_DEGs_rpmi<-res_df_rna_rpmi[res_df_rna_rpmi$hgnc_symbol %in% mRNAbind,]
RNA_DEGs_rpmi<-subset(RNA_DEGs_rpmi,RNA_DEGs_rpmi$pvalue<0.05)




### Enhancer peaks as in LOUCY AND DND41 ENCODE
enh_enc<-read.delim('/Volumes/grcmc/YGUILLEN/beta_catenin_project/bed_files_chips/DND_41_LOUCY_H3K4me1/sets/peaks_enhancers.txt',header = FALSE)

enhenc_Li<-read.delim('/Volumes/grcmc/YGUILLEN/beta_catenin_project/bed_files_chips/DND_41_LOUCY_H3K4me1/sets/peaks_enhancers_Liexc.txt',header = FALSE)

enhenc_Cont<-read.delim('/Volumes/grcmc/YGUILLEN/beta_catenin_project/bed_files_chips/DND_41_LOUCY_H3K4me1/sets/peaks_enhancers_Ctexc.txt',header = FALSE)

allbcatpeaks<-rbind(dataf_peakannoBRL2,dataf_peakannoRBL2,
      dataf_peakannoBRC1,dataf_peakannoBRC2,
      dataf_peakannoDB1,dataf_peakannoDB2,dataf_peakannoDB3)
allbcatpeaks$V4<-gsub('^.*peakcall.','',allbcatpeaks$V4)

enhpeaks<-allbcatpeaks[allbcatpeaks$V4 %in% enh_enc$V1,]
enhpeaksLi<-allbcatpeaks[allbcatpeaks$V4 %in% enhenc_Li$V1,]
enhpeaksCont<-allbcatpeaks[allbcatpeaks$V4 %in% enhenc_Cont$V1,]

# peaks in promoters
promenh<-enhpeaks[enhpeaks$annotation=="Promoter",]
table(gsub('_peaks.*','',promenh$V4),promenh$SYMBOL)

# peaks in promoters Li
promenh<-enhpeaksLi[enhpeaksLi$annotation=="Promoter",]
table(gsub('_peaks.*','',promenh$V4),promenh$SYMBOL)

# peaks in promoters Control
promenh<-enhpeaksCont[enhpeaksCont$annotation=="Promoter",]
table(gsub('_peaks.*','',promenh$V4),promenh$SYMBOL)








#table with bcat peaks
allbcatpeaks

unique(allbcatpeaks$ENSEMBL)
unique(allbcatpeaks$SYMBOL)

ensembl76 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

anotense <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
            values=unique(allbcatpeaks$ENSEMBL),
            filter="ensembl_gene_id",
            mart=ensembl76) 

#empty values as NA
anotense$go_id<-gsub("^$", NA, anotense$go_id)
anotense<-anotense[!is.na(anotense$go_id),]


geneGObcat<-merge(allbcatpeaks,anotense,by.x="ENSEMBL",by.y="ensembl_gene_id",all.x=TRUE)


geneGObcat<-merge(geneGObcat,vasbp,by="go_id",all.x=TRUE)

### Select peaks with keywords mRNA and splicing
subt<-geneGObcat[!grepl('splic',geneGObcat$Term) & !grepl('mRNA',geneGObcat$Term) & !grepl('regulation of transcription',geneGObcat$Term),]
unique(subt$Term)
unique(subt$SYMBOL)

RNApeaks<-allbcatpeaks[allbcatpeaks$ENSEMBL %in% unique(subt$ENSEMBL),]


# get peak regions from merged narrowPeak beds
bedbcat<-read.delim('/Volumes/cancer/bcat_Project/peakcall/allbcat_peaks.narrowPeak',header = FALSE)
bedbcat<-bedbcat[bedbcat$V4 %in% RNApeaks$V4,]
write.table(bedbcat,'/Volumes/cancer/bcat_Project/peakcall/bed_inter/RNApeaks.bed',row.names = FALSE,quote = FALSE,sep="\t")




## RNASeq shbcat GO annotation
# JURKAT

sig_deg_jurkat<-res_df_rna_jurkat[res_df_rna_jurkat$padj<0.1 & !is.na(res_df_rna_jurkat$padj),]

anojurk <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                  values=unique(sig_deg_jurkat$id),
                  filter="ensembl_gene_id",
                  mart=ensembl76) 

#empty values as NA
anojurk$go_id<-gsub("^$", NA, anojurk$go_id)
anojurk<-anojurk[!is.na(anojurk$go_id),]

sig_deg_jurkat<-merge(sig_deg_jurkat,anojurk,by.x="id",by.y="ensembl_gene_id",all.x=TRUE)
sig_deg_jurkat<-sig_deg_jurkat[!duplicated(sig_deg_jurkat), ] 



sig_deg_jurkat<-merge(sig_deg_jurkat,vasbp,by="go_id",all.x=TRUE)
sig_deg_jurkat_GO<-as.data.frame(table(sig_deg_jurkat$Term))
sig_deg_jurkat_GO<-sig_deg_jurkat_GO[order(sig_deg_jurkat_GO$Freq,decreasing = TRUE),]


#### MOTIF DISCOVERY
BiocManager::install('JASPAR2018')
#BiocManager::install('TFBSTools')

require(JASPAR2018)
require(TFBSTools)

pfm <- getMatrixSet(JASPAR2018, list(species=9606,name="HOXA4"))[[1]]
pfm

#ICM --> information content matrix
icm <- toICM(pfm)
seqLogo(icm)
dev.off()
