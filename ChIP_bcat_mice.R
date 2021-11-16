#Load libraries
library(data.table)
library(GenomicAlignments)
library(GO.db)
library(DO.db)
library(DiffBind)
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
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ReactomePA)
library(trackViewer)
library(Gviz)
library(rtracklayer)
library(Sushi)
library(venneuler)
library(VennDiagram)
library(readxl)
library(biomaRt)
library(viridis)
library(vegan)
library(ggrepel)

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#setwd("/Volumes/cancer/bcat_Project/peakcall/")
#setwd("/Users/yolanda_guillen/Desktop/IMIM/bcat_Project/peakcall_mod/")
setwd("/Users/yguillen/Desktop/temp/beta_catenin_project/peakcall_mod/")

#When reading narrowPeak, header must be included.

## bcat mice
# bcat fetal liver NOTCH
peakFIR1 <- readPeakFile("FIR1_peaks_peaks.narrowPeak")
peakFISC <- readPeakFile("FISC_peaks_peaks.narrowPeak")
peakFIXP <- readPeakFile("FIXP_peaks_peaks.narrowPeak")

# FL NOTCH bcat, excluding black listed regions and contamination, merged peaks
peakFI_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_FL_NOTCHmerge_nearest.genes.txt")

peakFWR1 <- readPeakFile("FWR1_peaks_peaks.narrowPeak")
peakFWR2 <- readPeakFile("FWR2_peaks_peaks.narrowPeak")
peakFWSC <- readPeakFile("FWSC_peaks_peaks.narrowPeak")
peakFWXP <- readPeakFile("FWXP_peaks_peaks.narrowPeak")

# FL WT bcat, excluding black listed regions and contamination, merged peaks
peakFW_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_FL_WTmerge_nearest.genes.txt")


# bcat leukemic mice
peakBMLB1 <- readPeakFile("BMLB1_peaks_peaks.narrowPeak")
peakBMLB2 <- readPeakFile("BMLB2_peaks_peaks.narrowPeak") 
peakBMLB3 <- readPeakFile("BMLB3_peaks_peaks.narrowPeak")


peakBMLT1 <- readPeakFile("BMLT1_peaks_peaks.narrowPeak") 
peakBMLT2 <- readPeakFile("BMLT2_peaks_peaks.narrowPeak") 
peakBMLT3 <- readPeakFile("BMLT3_peaks_peaks.narrowPeak") 

# Leukemic bcat and TCF BM, excluding black listed regions and contamination, merged peaks
peakBMLB_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_BMLB_merge_nearest.genes.txt")
peakBMLT_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_BMLT_merge_nearest.genes.txt")



## RPMI cell line peaks, after filter and merging (min one replicate, not considering RCB1)
peakDB_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_DBmerge_nearest.genes.txt")
peakDBL_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_DBLmerge_nearest.genes.txt")
peakTCF_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_TCFmerge_nearest.genes.txt")
peakLEF_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_LEFmerge_nearest.genes.txt")

# RPMI SANTA CRUZ ANTIBODIES
peakRBCS_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_RBCSmerge_nearest.genes.txt")
peakRBLS_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_RBLSmerge_nearest.genes.txt")


# DKO
peakRKB_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_RKB_merge_nearest.genes.txt")
peakRKT_black <- readPeakFile("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_RKT_merge_nearest.genes.txt")


covplot(peakRKT_black, weightCol="Score")

# mice genes database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdb_human <- TxDb.Hsapiens.UCSC.hg38.knownGene

### PEAK ANNOTATION
peakAnnoFIR1 <- annotatePeak(peakFIR1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFISC <- annotatePeak(peakFISC, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFIXP <- annotatePeak(peakFIXP, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnnoFI_black <- annotatePeak(peakFI_black, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnnoFWR1 <- annotatePeak(peakFWR1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFWR2 <- annotatePeak(peakFWR2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFWSC <- annotatePeak(peakFWSC, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFWXP <- annotatePeak(peakFWXP, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnnoFW_black <- annotatePeak(peakFW_black, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

#leukemic bcat
peakAnnoBMLB1 <- annotatePeak(peakBMLB1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoBMLB2 <- annotatePeak(peakBMLB2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoBMLB3 <- annotatePeak(peakBMLB3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnnoBMLB_black <- annotatePeak(peakBMLB_black, tssRegion = c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnnoBMLT1 <- annotatePeak(peakBMLT1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoBMLT2 <- annotatePeak(peakBMLT2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoBMLT3 <- annotatePeak(peakBMLT3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnnoBMLT_black <- annotatePeak(peakBMLT_black, tssRegion = c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")


# RPMI cell lines
peakAnnoDB_black <- annotatePeak(peakDB_black, tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnnoDBL_black <- annotatePeak(peakDBL_black, tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnnoLEF_black <- annotatePeak(peakLEF_black, tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnnoTCF_black <- annotatePeak(peakTCF_black, tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")

peakAnnoRBCS_black <- annotatePeak(peakRBCS_black, tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnnoRBLS_black <- annotatePeak(peakRBLS_black, tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")

peakAnnoRKB_black <- annotatePeak(peakRKB_black, tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnnoRKT_black <- annotatePeak(peakRKT_black, tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")


# Visualize genomic annotation
plotAnnoPie(peakAnnoBMLB_black)
dev.off()
plotAnnoBar(peakAnnoFW_black)
dev.off()
vennpie(peakAnnoFW_black)
dev.off()


plotDistToTSS(peakAnnoBMLB_black,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


##DATAFRAME WITH genes id linked to peaks

#bcatenin notch1 RD Ab
dataf_peakannoFIR1<-as.data.frame(peakAnnoFIR1)
dataf_peakannoFIR1$annotation<-gsub('Intron.*','Intron',dataf_peakannoFIR1$annotation)
dataf_peakannoFIR1$annotation<-gsub('Exon.*','Exon',dataf_peakannoFIR1$annotation)
dataf_peakannoFIR1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFIR1$annotation)
dataf_peakannoFIR1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFIR1$annotation)
dataf_peakannoFIR1$peaksite<-paste(dataf_peakannoFIR1$SYMBOL,dataf_peakannoFIR1$annotation,sep='_')
table(dataf_peakannoFIR1$annotation)

#bcatenin notch1 SC Ab
dataf_peakannoFISC<-as.data.frame(peakAnnoFISC)
dataf_peakannoFISC$annotation<-gsub('Intron.*','Intron',dataf_peakannoFISC$annotation)
dataf_peakannoFISC$annotation<-gsub('Exon.*','Exon',dataf_peakannoFISC$annotation)
dataf_peakannoFISC$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFISC$annotation)
dataf_peakannoFISC$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFISC$annotation)
dataf_peakannoFISC$peaksite<-paste(dataf_peakannoFISC$SYMBOL,dataf_peakannoFISC$annotation,sep='_')
table(dataf_peakannoFISC$annotation)

#bcatenin notch1 XP Ab
dataf_peakannoFIXP<-as.data.frame(peakAnnoFIXP)
dataf_peakannoFIXP$annotation<-gsub('Intron.*','Intron',dataf_peakannoFIXP$annotation)
dataf_peakannoFIXP$annotation<-gsub('Exon.*','Exon',dataf_peakannoFIXP$annotation)
dataf_peakannoFIXP$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFIXP$annotation)
dataf_peakannoFIXP$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFIXP$annotation)
dataf_peakannoFIXP$peaksite<-paste(dataf_peakannoFIXP$SYMBOL,dataf_peakannoFIXP$annotation,sep='_')
table(dataf_peakannoFIXP$annotation)

#bcatenin merge FL NOTCH merge
dataf_peakannoFI_black<-as.data.frame(peakAnnoFI_black)
dataf_peakannoFI_black$annotation<-gsub('Intron.*','Intron',dataf_peakannoFI_black$annotation)
dataf_peakannoFI_black$annotation<-gsub('Exon.*','Exon',dataf_peakannoFI_black$annotation)
dataf_peakannoFI_black$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFI_black$annotation)
dataf_peakannoFI_black$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFI_black$annotation)
dataf_peakannoFI_black$peaksite<-paste(dataf_peakannoFI_black$SYMBOL,dataf_peakannoFI_black$annotation,sep='_')
table(dataf_peakannoFI_black$annotation)


#bcatenin WT RD Ab rep 1
dataf_peakannoFWR1<-as.data.frame(peakAnnoFWR1)
dataf_peakannoFWR1$annotation<-gsub('Intron.*','Intron',dataf_peakannoFWR1$annotation)
dataf_peakannoFWR1$annotation<-gsub('Exon.*','Exon',dataf_peakannoFWR1$annotation)
dataf_peakannoFWR1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFWR1$annotation)
dataf_peakannoFWR1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFWR1$annotation)
dataf_peakannoFWR1$peaksite<-paste(dataf_peakannoFWR1$SYMBOL,dataf_peakannoFWR1$annotation,sep='_')
table(dataf_peakannoFWR1$annotation)

#bcatenin WT RD Ab rep 2
dataf_peakannoFWR2<-as.data.frame(peakAnnoFWR2)
dataf_peakannoFWR2$annotation<-gsub('Intron.*','Intron',dataf_peakannoFWR2$annotation)
dataf_peakannoFWR2$annotation<-gsub('Exon.*','Exon',dataf_peakannoFWR2$annotation)
dataf_peakannoFWR2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFWR2$annotation)
dataf_peakannoFWR2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFWR2$annotation)
dataf_peakannoFWR2$peaksite<-paste(dataf_peakannoFWR2$SYMBOL,dataf_peakannoFWR2$annotation,sep='_')
table(dataf_peakannoFWR2$annotation)

#bcatenin WT SC Ab 
dataf_peakannoFWSC<-as.data.frame(peakAnnoFWSC)
dataf_peakannoFWSC$annotation<-gsub('Intron.*','Intron',dataf_peakannoFWSC$annotation)
dataf_peakannoFWSC$annotation<-gsub('Exon.*','Exon',dataf_peakannoFWSC$annotation)
dataf_peakannoFWSC$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFWSC$annotation)
dataf_peakannoFWSC$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFWSC$annotation)
dataf_peakannoFWSC$peaksite<-paste(dataf_peakannoFWSC$SYMBOL,dataf_peakannoFWSC$annotation,sep='_')
table(dataf_peakannoFWSC$annotation)

#bcatenin WT XP Ab 
dataf_peakannoFWXP<-as.data.frame(peakAnnoFWXP)
dataf_peakannoFWXP$annotation<-gsub('Intron.*','Intron',dataf_peakannoFWXP$annotation)
dataf_peakannoFWXP$annotation<-gsub('Exon.*','Exon',dataf_peakannoFWXP$annotation)
dataf_peakannoFWXP$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFWXP$annotation)
dataf_peakannoFWXP$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFWXP$annotation)
dataf_peakannoFWXP$peaksite<-paste(dataf_peakannoFWXP$SYMBOL,dataf_peakannoFWXP$annotation,sep='_')
table(dataf_peakannoFWXP$annotation)


#bcatenin merge FL WT merge
dataf_peakannoFW_black<-as.data.frame(peakAnnoFW_black)
dataf_peakannoFW_black$annotation<-gsub('Intron.*','Intron',dataf_peakannoFW_black$annotation)
dataf_peakannoFW_black$annotation<-gsub('Exon.*','Exon',dataf_peakannoFW_black$annotation)
dataf_peakannoFW_black$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFW_black$annotation)
dataf_peakannoFW_black$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFW_black$annotation)
dataf_peakannoFW_black$peaksite<-paste(dataf_peakannoFW_black$SYMBOL,dataf_peakannoFW_black$annotation,sep='_')
table(dataf_peakannoFW_black$annotation)


#bcatenin leukemia mice 
dataf_peakannoBMLB1<-as.data.frame(peakAnnoBMLB1)
dataf_peakannoBMLB1$annotation<-gsub('Intron.*','Intron',dataf_peakannoBMLB1$annotation)
dataf_peakannoBMLB1$annotation<-gsub('Exon.*','Exon',dataf_peakannoBMLB1$annotation)
dataf_peakannoBMLB1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBMLB1$annotation)
dataf_peakannoBMLB1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBMLB1$annotation)
dataf_peakannoBMLB1$peaksite<-paste(dataf_peakannoBMLB1$SYMBOL,dataf_peakannoBMLB1$annotation,sep='_')
table(dataf_peakannoBMLB1$annotation)

dataf_peakannoBMLB2<-as.data.frame(peakAnnoBMLB2)
dataf_peakannoBMLB2$annotation<-gsub('Intron.*','Intron',dataf_peakannoBMLB2$annotation)
dataf_peakannoBMLB2$annotation<-gsub('Exon.*','Exon',dataf_peakannoBMLB2$annotation)
dataf_peakannoBMLB2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBMLB2$annotation)
dataf_peakannoBMLB2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBMLB2$annotation)
dataf_peakannoBMLB2$peaksite<-paste(dataf_peakannoBMLB2$SYMBOL,dataf_peakannoBMLB2$annotation,sep='_')
table(dataf_peakannoBMLB1$annotation)

dataf_peakannoBMLB3<-as.data.frame(peakAnnoBMLB3)
dataf_peakannoBMLB3$annotation<-gsub('Intron.*','Intron',dataf_peakannoBMLB3$annotation)
dataf_peakannoBMLB3$annotation<-gsub('Exon.*','Exon',dataf_peakannoBMLB3$annotation)
dataf_peakannoBMLB3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBMLB3$annotation)
dataf_peakannoBMLB3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBMLB3$annotation)
dataf_peakannoBMLB3$peaksite<-paste(dataf_peakannoBMLB3$SYMBOL,dataf_peakannoBMLB3$annotation,sep='_')
table(dataf_peakannoBMLB3$annotation)


dataf_peakannoBMLB_black<-as.data.frame(peakAnnoBMLB_black)
dataf_peakannoBMLB_black$annotation<-gsub('Intron.*','Intron',dataf_peakannoBMLB_black$annotation)
dataf_peakannoBMLB_black$annotation<-gsub('Exon.*','Exon',dataf_peakannoBMLB_black$annotation)
dataf_peakannoBMLB_black$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBMLB_black$annotation)
dataf_peakannoBMLB_black$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBMLB_black$annotation)
dataf_peakannoBMLB_black$peaksite<-paste(dataf_peakannoBMLB_black$SYMBOL,dataf_peakannoBMLB_black$annotation,sep='_')
table(dataf_peakannoBMLB_black$annotation)


# TCF leukemia mice 
dataf_peakannoBMLT1<-as.data.frame(peakAnnoBMLT1)
dataf_peakannoBMLT1$annotation<-gsub('Intron.*','Intron',dataf_peakannoBMLT1$annotation)
dataf_peakannoBMLT1$annotation<-gsub('Exon.*','Exon',dataf_peakannoBMLT1$annotation)
dataf_peakannoBMLT1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBMLT1$annotation)
dataf_peakannoBMLT1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBMLT1$annotation)
dataf_peakannoBMLT1$peaksite<-paste(dataf_peakannoBMLT1$SYMBOL,dataf_peakannoBMLT1$annotation,sep='_')
table(dataf_peakannoBMLT1$annotation)

dataf_peakannoBMLT2<-as.data.frame(peakAnnoBMLT2)
dataf_peakannoBMLT2$annotation<-gsub('Intron.*','Intron',dataf_peakannoBMLT2$annotation)
dataf_peakannoBMLT2$annotation<-gsub('Exon.*','Exon',dataf_peakannoBMLT2$annotation)
dataf_peakannoBMLT2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBMLT2$annotation)
dataf_peakannoBMLT2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBMLT2$annotation)
dataf_peakannoBMLT2$peaksite<-paste(dataf_peakannoBMLT2$SYMBOL,dataf_peakannoBMLT2$annotation,sep='_')
table(dataf_peakannoBMLT2$annotation)

dataf_peakannoBMLT3<-as.data.frame(peakAnnoBMLT3)
dataf_peakannoBMLT3$annotation<-gsub('Intron.*','Intron',dataf_peakannoBMLT3$annotation)
dataf_peakannoBMLT3$annotation<-gsub('Exon.*','Exon',dataf_peakannoBMLT3$annotation)
dataf_peakannoBMLT3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBMLT3$annotation)
dataf_peakannoBMLT3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBMLT3$annotation)
dataf_peakannoBMLT3$peaksite<-paste(dataf_peakannoBMLT3$SYMBOL,dataf_peakannoBMLT3$annotation,sep='_')
table(dataf_peakannoBMLT3$annotation)

dataf_peakannoBMLT_black<-as.data.frame(peakAnnoBMLT_black)
dataf_peakannoBMLT_black$annotation<-gsub('Intron.*','Intron',dataf_peakannoBMLT_black$annotation)
dataf_peakannoBMLT_black$annotation<-gsub('Exon.*','Exon',dataf_peakannoBMLT_black$annotation)
dataf_peakannoBMLT_black$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoBMLT_black$annotation)
dataf_peakannoBMLT_black$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoBMLT_black$annotation)
dataf_peakannoBMLT_black$peaksite<-paste(dataf_peakannoBMLT_black$SYMBOL,dataf_peakannoBMLT_black$annotation,sep='_')
table(dataf_peakannoBMLT_black$annotation)

# List of genes
FIR1_genes<-unique(dataf_peakannoFIR1$SYMBOL)
FISC_genes<-unique(dataf_peakannoFISC$SYMBOL)
FIXP_genes<-unique(dataf_peakannoFIXP$SYMBOL)

FI_black_genes<-unique(dataf_peakannoFI_black$gene)

FWR1_genes<-unique(dataf_peakannoFWR1$SYMBOL)
FWR2_genes<-unique(dataf_peakannoFWR2$SYMBOL)
FWSC_genes<-unique(dataf_peakannoFWSC$SYMBOL)
FWXP_genes<-unique(dataf_peakannoFWXP$SYMBOL)

FW_black_genes<-unique(dataf_peakannoFW_black$gene)

BMLB1_genes<-sort(unique(dataf_peakannoBMLB1$SYMBOL))
BMLB2_genes<-sort(unique(dataf_peakannoBMLB2$SYMBOL))
BMLB3_genes<-sort(unique(dataf_peakannoBMLB3$SYMBOL))

BMLB_black_genes<-sort(unique(dataf_peakannoBMLB_black$gene))

BMLT1_genes<-sort(unique(dataf_peakannoBMLT1$SYMBOL))
BMLT2_genes<-sort(unique(dataf_peakannoBMLT2$SYMBOL))
BMLT3_genes<-sort(unique(dataf_peakannoBMLT3$SYMBOL))

BMLT_black_genes<-sort(unique(dataf_peakannoBMLT_black$gene))

## Output list of genes for each sample
#write.table(dataf_peakannoFIR1,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/FIR1.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
#write.table(dataf_peakannoFISC,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/FISC.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
#write.table(dataf_peakannoFIXP,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/FIXP.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
#write.table(dataf_peakannoFWR1,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/FWR1.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
#write.table(dataf_peakannoFWR2,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/FWR2.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
#write.table(dataf_peakannoFWSC,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/FWSC.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
#write.table(dataf_peakannoFWXP,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/FWXP.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
#write.table(bcat_Notch_genes,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/bcat_Notch_genes.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
#write.table(bcat_WT_genes,"/Volumes/grcmc/YGUILLEN/beta_catenin_project/mice_chips/bcat_WT_genes.csv",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")

## Common genes bcat peaks RPMI and leukemic mice

# peaks bcat control RPMI from script ChIPSeq.R minimum 2 replicates
bcatcomgenes[bcatcomgenes$Freq>=2,]$Var1

# peaks bcat Lithium RPMI from script ChIPSeq.R
bcatlitcomgenes[bcatlitcomgenes$Freq>=1,]$Var1


# After substracting blacklisted regions and contamination, genes with peak in at least one sample
BMLB_black_genes
BMLT_black_genes
FW_black_genes
FI_black_genes


# venn
geneLists<-list(Mice_BM = BMLB_black_genes,
                Mice_FL_NICD = FI_black_genes,
                Mice_FL_WT = FW_black_genes)

geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off() 

# FOR THREE GROUPS
venn.plot <- venn.diagram(VENN.LIST , NULL, 
                          fill=c("darkolivegreen", "blue","darkolivegreen2"),
                          alpha=c(0.3,0.3,0.3), 
                          cex = 2, 
                          cat.fontface=4, 
                          category.names=c("mice_BM","mice_FL_NICD","mice_FL_WT"), 
                          main="Genes ChIP bcat min 1 replicate, filtered")


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

# Peak distribution

#FWR1
distFWR1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFWR1$annotation)))))
colnames(distFWR1) <- as.character(unlist(distFWR1[1,]))
distFWR1<-distFWR1[-1,]
distFWR1$Sample<-"FWR1"

colnames(distFWR1)
colnames(distFWR1)[1]<-c("Three_UTR")
colnames(distFWR1)[2]<-c("Five_UTR")
colnames(distFWR1)[3]<-c("Distal_intergenic")

#FWR2
distFWR2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFWR2$annotation)))))
colnames(distFWR2) <- as.character(unlist(distFWR2[1,]))
distFWR2<-distFWR2[-1,]
distFWR2$Sample<-"FWR2"

colnames(distFWR2)
distFWR2$Three_UTR<-"0"
distFWR2$Five_UTR<-"0"
distFWR2$Downstream<-"0"
distFWR2<-distFWR2[,c(6,7,1,8,2,3,4,5)]
colnames(distFWR2)[1]<-c("Three_UTR")
colnames(distFWR2)[2]<-c("Five_UTR")
colnames(distFWR2)[3]<-c("Distal_intergenic")

#FWSC
distFWSC<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFWSC$annotation)))))
colnames(distFWSC) <- as.character(unlist(distFWSC[1,]))
distFWSC<-distFWSC[-1,]
distFWSC$Sample<-"FWSC"

colnames(distFWSC)
distFWSC$Three_UTR<-"0"
distFWSC$Five_UTR<-"0"
distFWSC$Downstream<-"0"
distFWSC<-distFWSC[,c(6,7,1,8,2,3,4,5)]
colnames(distFWSC)[1]<-c("Three_UTR")
colnames(distFWSC)[2]<-c("Five_UTR")
colnames(distFWSC)[3]<-c("Distal_intergenic")



#FWXP
distFWXP<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFWXP$annotation)))))
colnames(distFWXP) <- as.character(unlist(distFWXP[1,]))
distFWXP<-distFWXP[-1,]
distFWXP$Sample<-"FWXP"

colnames(distFWXP)
distFWXP$Three_UTR<-"0"
distFWXP$Five_UTR<-"0"
distFWXP$Downstream<-"0"
distFWXP<-distFWXP[,c(6,7,1,8,2,3,4,5)]
colnames(distFWXP)[1]<-c("Three_UTR")
colnames(distFWXP)[2]<-c("Five_UTR")
colnames(distFWXP)[3]<-c("Distal_intergenic")


#FW_merge filtered
distFW_black<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFW_black$annotation)))))
colnames(distFW_black) <- as.character(unlist(distFW_black[1,]))
distFW_black<-distFW_black[-1,]
distFW_black$Sample<-"FW_black"

colnames(distFW_black)
colnames(distFW_black)[1]<-c("Three_UTR")
colnames(distFW_black)[2]<-c("Five_UTR")
colnames(distFW_black)[3]<-c("Distal_intergenic")

#FIR1
distFIR1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFIR1$annotation)))))
colnames(distFIR1) <- as.character(unlist(distFIR1[1,]))
distFIR1<-distFIR1[-1,]
distFIR1$Sample<-"FIR1"

colnames(distFIR1)
distFIR1$Three_UTR<-"0"
distFIR1$Five_UTR<-"0"
distFIR1$Downstream<-"0"
distFIR1<-distFIR1[,c(6,7,1,8,2,3,4,5)]
colnames(distFIR1)[1]<-c("Three_UTR")
colnames(distFIR1)[2]<-c("Five_UTR")
colnames(distFIR1)[3]<-c("Distal_intergenic")

#FISC
distFISC<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFISC$annotation)))))
colnames(distFISC) <- as.character(unlist(distFISC[1,]))
distFISC<-distFISC[-1,]
distFISC$Sample<-"FISC"

colnames(distFISC)
colnames(distFISC)[1]<-c("Three_UTR")
colnames(distFISC)[2]<-c("Five_UTR")
colnames(distFISC)[3]<-c("Distal_intergenic")


#FIXP
distFIXP<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFIXP$annotation)))))
colnames(distFIXP) <- as.character(unlist(distFIXP[1,]))
distFIXP<-distFIXP[-1,]
distFIXP$Sample<-"FIXP"

colnames(distFIXP)
distFIXP$Five_UTR<-"0"
distFIXP<-distFIXP[,c(1,8,2,3,4,5,6,7)]
colnames(distFIXP)[1]<-c("Three_UTR")
colnames(distFIXP)[2]<-c("Five_UTR")
colnames(distFIXP)[3]<-c("Distal_intergenic")

#FI_merge filtered
distFI_black<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFI_black$annotation)))))
colnames(distFI_black) <- as.character(unlist(distFI_black[1,]))
distFI_black<-distFI_black[-1,]
distFI_black$Sample<-"FI_black"

colnames(distFI_black)
colnames(distFI_black)[1]<-c("Three_UTR")
colnames(distFI_black)[2]<-c("Five_UTR")
colnames(distFI_black)[3]<-c("Distal_intergenic")


#BMLB1
distBMLB1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBMLB1$annotation)))))
colnames(distBMLB1) <- as.character(unlist(distBMLB1[1,]))
distBMLB1<-distBMLB1[-1,]
distBMLB1$Sample<-"BMLB1"

colnames(distBMLB1)
distBMLB1$Five_UTR<-"0"
distBMLB1$Downstream<-"0"
distBMLB1<-distBMLB1[,c(1,7,2,8,3,4,5,6)]
colnames(distBMLB1)[1]<-c("Three_UTR")
colnames(distBMLB1)[2]<-c("Five_UTR")
colnames(distBMLB1)[3]<-c("Distal_intergenic")


#BMLB2
distBMLB2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBMLB2$annotation)))))
colnames(distBMLB2) <- as.character(unlist(distBMLB2[1,]))
distBMLB2<-distBMLB2[-1,]
distBMLB2$Sample<-"BMLB2"

colnames(distBMLB2)
colnames(distBMLB2)[1]<-c("Three_UTR")
colnames(distBMLB2)[2]<-c("Five_UTR")
colnames(distBMLB2)[3]<-c("Distal_intergenic")

#BMLB3
distBMLB3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBMLB3$annotation)))))
colnames(distBMLB3) <- as.character(unlist(distBMLB3[1,]))
distBMLB3<-distBMLB3[-1,]
distBMLB3$Sample<-"BMLB3"

colnames(distBMLB3)
distBMLB3$Five_UTR<-"0"
distBMLB3$Downstream<-"0"
distBMLB3<-distBMLB3[,c(1,7,2,8,3,4,5,6)]
colnames(distBMLB3)[1]<-c("Three_UTR")
colnames(distBMLB3)[2]<-c("Five_UTR")
colnames(distBMLB3)[3]<-c("Distal_intergenic")


#BMLB merged and filtered black
distBMLB_black<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBMLB_black$annotation)))))
colnames(distBMLB_black) <- as.character(unlist(distBMLB_black[1,]))
distBMLB_black<-distBMLB_black[-1,]
distBMLB_black$Sample<-"BMLB_black"

colnames(distBMLB_black)
colnames(distBMLB_black)[1]<-c("Three_UTR")
colnames(distBMLB_black)[2]<-c("Five_UTR")
colnames(distBMLB_black)[3]<-c("Distal_intergenic")


#BMLT1
distBMLT1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBMLT1$annotation)))))
colnames(distBMLT1) <- as.character(unlist(distBMLT1[1,]))
distBMLT1<-distBMLT1[-1,]
distBMLT1$Sample<-"BMLT1"

colnames(distBMLT1)
colnames(distBMLT1)[1]<-c("Three_UTR")
colnames(distBMLT1)[2]<-c("Five_UTR")
colnames(distBMLT1)[3]<-c("Distal_intergenic")


#BMLT2
distBMLT2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBMLT2$annotation)))))
colnames(distBMLT2) <- as.character(unlist(distBMLT2[1,]))
distBMLT2<-distBMLT2[-1,]
distBMLT2$Sample<-"BMLT2"

colnames(distBMLT2)
colnames(distBMLT2)[1]<-c("Three_UTR")
colnames(distBMLT2)[2]<-c("Five_UTR")
colnames(distBMLT2)[3]<-c("Distal_intergenic")

#BMLT3
distBMLT3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBMLT3$annotation)))))
colnames(distBMLT3) <- as.character(unlist(distBMLT3[1,]))
distBMLT3<-distBMLT3[-1,]
distBMLT3$Sample<-"BMLT3"

colnames(distBMLT3)
colnames(distBMLT3)[1]<-c("Three_UTR")
colnames(distBMLT3)[2]<-c("Five_UTR")
colnames(distBMLT3)[3]<-c("Distal_intergenic")

#BMLT merged and filtered black
distBMLT_black<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoBMLT_black$annotation)))))
colnames(distBMLT_black) <- as.character(unlist(distBMLT_black[1,]))
distBMLT_black<-distBMLT_black[-1,]
distBMLT_black$Sample<-"BMLT_black"

colnames(distBMLT_black)
colnames(distBMLT_black)[1]<-c("Three_UTR")
colnames(distBMLT_black)[2]<-c("Five_UTR")
colnames(distBMLT_black)[3]<-c("Distal_intergenic")


# plot dist
distmice<-rbind(distBMLB1,distBMLB2,distBMLB3,distBMLB_black,
                distBMLT1,distBMLT2,distBMLT3,distBMLT_black)


distmicechip<- data.frame(sapply(distmice[,1:7], function(x) as.numeric(as.character(x))*100))
distmicechip<-cbind(distmice$Sample,distmicechip)
colnames(distmicechip)[1]<-"Sample"
distmicechip$condition<-c("Bcat","Bcat","Bcat","Bcat_filter",
                          "Tcf","Tcf","Tcf","Tcf_filter")

distmicechip_melt<-melt(distmicechip,id.vars = c("condition","Sample"))

ggplot(distmicechip_melt,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat="identity",color="black") +
  facet_grid(~condition,scale="free",space="free")+
  scale_fill_brewer(palette="Pastel1")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")


distmicefl<-rbind(distFWR1,distFWR2,distFWXP,distFWSC,distFW_black,
                  distFIR1,distFISC,distFIXP,distFI_black)

distmicechip<- data.frame(sapply(distmicefl[,1:7], function(x) as.numeric(as.character(x))*100))
distmicechip<-cbind(distmicefl$Sample,distmicechip)
colnames(distmicechip)[1]<-"Sample"
distmicechip$condition<-c("WT","WT","WT","WT","WT_black",
                          "NICD","NICD","NICD","NICD_black")

distmicechip_melt<-melt(distmicechip,id.vars = c("condition","Sample"))

ggplot(distmicechip_melt,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat="identity",color="black") +
  facet_grid(~condition,scale="free",space="free")+
  scale_fill_brewer(palette="Pastel1")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")



### Distribution score peaks to detect contamination

# DKO bcat
dataf_peakannoBMLB1$sample<-"BMLB1"
dataf_peakannoBMLB2$sample<-"BMLB2"
dataf_peakannoBMLB3$sample<-"BMLB3"

dataf_peakannoBMLT1$sample<-"BMLT1"
dataf_peakannoBMLT2$sample<-"BMLT2"
dataf_peakannoBMLT3$sample<-"BMLT3"


micepeaks<-rbind(dataf_peakannoBMLB1,dataf_peakannoBMLB2,dataf_peakannoBMLB3,
                 dataf_peakannoBMLT1,dataf_peakannoBMLT2,dataf_peakannoBMLT3)

ggplot(micepeaks,aes(x=sample,y=score,label=SYMBOL))+
  geom_boxplot()+
  geom_point()+
  geom_label_repel(data=micepeaks[micepeaks$score>1000,],aes(label=SYMBOL), size=3,fontface = "italic",alpha=0.5)+
  facet_grid(~sample,scales = "free_x")+
  theme_bw()


# GO annotation of mice chips genes
mm<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

BMLB_gene_anot<-getBM(attributes=c("ensembl_gene_id", "description","go_id","mgi_symbol"),
      values= BMLB_black_genes,
      filter="mgi_symbol",
      mart=mm) 

BMLT_gene_anot<-getBM(attributes=c("ensembl_gene_id", "description","go_id","mgi_symbol"),
                      values= BMLT_black_genes,
                      filter="mgi_symbol",
                      mart=mm) 

FI_gene_anot<-getBM(attributes=c("ensembl_gene_id", "description","go_id","mgi_symbol"),
                      values= FI_black_genes,
                      filter="mgi_symbol",
                      mart=mm) 

FW_gene_anot<-getBM(attributes=c("ensembl_gene_id", "description","go_id","mgi_symbol"),
                      values= FW_black_genes,
                      filter="mgi_symbol",
                      mart=mm) 

# GO TERM table all genes
vas<-toTable(GOTERM)
colnames(vas)[2]<-"repgo_id"
vascc<-vas[vas$Ontology=="CC",]
vasbp<-vas[vas$Ontology=="BP",]
vasmf<-vas[vas$Ontology=="MF",]

# peaks
length(unique(BMLB_gene_anot$mgi_symbol))
length(unique(BMLB_black_genes))

length(unique(BMLT_gene_anot$mgi_symbol))
length(unique(BMLT_black_genes))

length(unique(FI_gene_anot$mgi_symbol))
length(unique(FI_black_genes))

length(unique(FW_gene_anot$mgi_symbol))
length(unique(FW_black_genes))

setdiff(unique(BMLB_black_genes),unique(BMLB_gene_anot$mgi_symbol))
setdiff(unique(BMLT_black_genes),unique(BMLT_gene_anot$mgi_symbol))
setdiff(unique(FI_black_genes),unique(FI_gene_anot$mgi_symbol))
setdiff(unique(FW_black_genes),unique(FW_gene_anot$mgi_symbol))

#empty values as NA
BMLB_gene_anot$go_id<-gsub("^$", NA,BMLB_gene_anot$go_id)
BMLB_gene_anot$go_id<-gsub("^$", NA,BMLB_gene_anot$go_id)

BMLT_gene_anot$go_id<-gsub("^$", NA,BMLT_gene_anot$go_id)
BMLT_gene_anot$go_id<-gsub("^$", NA,BMLT_gene_anot$go_id)

FI_gene_anot$go_id<-gsub("^$", NA,FI_gene_anot$go_id)
FI_gene_anot$go_id<-gsub("^$", NA,FI_gene_anot$go_id)

FW_gene_anot$go_id<-gsub("^$", NA,FW_gene_anot$go_id)
FW_gene_anot$go_id<-gsub("^$", NA,FW_gene_anot$go_id)

# all terms
#geneGO<-merge(geneGO,vas,by="go_id",all.x=TRUE)

#only CC, BP, MF term
geneGOCCBMLB<-merge(BMLB_gene_anot,vascc,by="go_id",all.x=TRUE)
geneGOBPBMLB<-merge(BMLB_gene_anot,vasbp,by="go_id",all.x=TRUE)
geneGOMFBMLB<-merge(BMLB_gene_anot,vasmf,by="go_id",all.x=TRUE)

geneGOCCBMLT<-merge(BMLT_gene_anot,vascc,by="go_id",all.x=TRUE)
geneGOBPBMLT<-merge(BMLT_gene_anot,vasbp,by="go_id",all.x=TRUE)
geneGOMFBMLT<-merge(BMLT_gene_anot,vasmf,by="go_id",all.x=TRUE)

geneGOCCFI<-merge(FI_gene_anot,vascc,by="go_id",all.x=TRUE)
geneGOBPFI<-merge(FI_gene_anot,vasbp,by="go_id",all.x=TRUE)
geneGOMFFI<-merge(FI_gene_anot,vasmf,by="go_id",all.x=TRUE)

geneGOCCFW<-merge(FW_gene_anot,vascc,by="go_id",all.x=TRUE)
geneGOBPFW<-merge(FW_gene_anot,vasbp,by="go_id",all.x=TRUE)
geneGOMFFW<-merge(FW_gene_anot,vasmf,by="go_id",all.x=TRUE)

sort(unique((geneGOBPFI[grepl('splicing',geneGOBPFI$Term) | grepl('RNA processing',geneGOBPFI$Term),])$mgi_symbol))


# Functional enrichment
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "InterPro_Domains_2019",
         "GO_Cellular_Component_2017",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2016")


geneBMLB <- enrichr(as.character(BMLB_black_genes), dbs)
geneBMLB <- geneBMLB[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
geneBMLB$group<-"BMLB"
geneBMLT <- enrichr(as.character(BMLT_black_genes), dbs)
geneBMLT <- geneBMLT[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
geneBMLT$group<-"BMLT"
geneFI <- enrichr(as.character(FI_black_genes), dbs)
geneFI <- geneFI[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
geneFI$group<-"FI"
geneFW <- enrichr(as.character(FW_black_genes), dbs)
geneFW <- geneFW[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
geneFW$group<-"FW"

allGO<-rbind(geneBMLB,geneBMLT,geneFI,geneFW)


bpsub<-subset(allGO,allGO$Adjusted.P.value<=0.1)
bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$group,bpsub$P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))
bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$group,bpsub$P.value)), ]$Term)

ggplot(bpsub,aes(y=-log10(P.value),x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group),color="black")+
  #geom_text(aes(label=tolower(Genes)), size=4,hjust=1,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Spectral")+
  #facet_wrap(~group,scales="free")+
  theme_minimal()+
  coord_flip()+
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=8, face="bold",hjust = 1),
        axis.text.y = element_text(size=10, hjust = 1),
        legend.position = "bottom")+
  facet_grid(~group,space = "free")


### COMPARING CELL LINE AND MICE CHIPS

# Heatmap distribution annotations

DB_freq<-peakAnnoDB_black@annoStat
DB_freq$group="bcat"
DB_freq$model="human_control"

DBL_freq<-peakAnnoDBL_black@annoStat
DBL_freq$group="bcat"
DBL_freq$model="human_Li"

RBCS_freq<-peakAnnoRBCS_black@annoStat
RBCS_freq$group="bcat_SC"
RBCS_freq$model="human_control"

RBLS_freq<-peakAnnoRBLS_black@annoStat
RBLS_freq$group="bcat_SC"
RBLS_freq$model="human_Li"

TCF_freq<-peakAnnoTCF_black@annoStat
TCF_freq$group="TCF"
TCF_freq$model="human_control"

LEF_freq<-peakAnnoLEF_black@annoStat
LEF_freq$group="LEF"
LEF_freq$model="human_control"

RKT_freq<-peakAnnoRKT_black@annoStat
RKT_freq$group="TCF"
RKT_freq$model="human_DKO"

RKB_freq<-peakAnnoRKB_black@annoStat
RKB_freq$group="bcat"
RKB_freq$model="human_DKO"

BMLB_freq<-peakAnnoBMLB_black@annoStat
BMLB_freq$group="bcat"
BMLB_freq$model="mouse_BM"

BMLT_freq<-peakAnnoBMLT_black@annoStat
BMLT_freq$group="TCF"
BMLT_freq$model="mouse_BM"

FI_freq<-peakAnnoFI_black@annoStat
FI_freq$group="bcat"
FI_freq$model="mouse_FLNI"

FW_freq<-peakAnnoFW_black@annoStat
FW_freq$group="bcat"
FW_freq$model="mouse_FLWT"

allfreq<-rbind(DB_freq,DBL_freq,RBCS_freq,RBLS_freq,TCF_freq,LEF_freq,RKT_freq,RKB_freq,
      BMLB_freq,BMLT_freq,FI_freq,FW_freq)

class(allfreq$group)

ggplot(allfreq, aes(x=Feature, y=group)) + 
  geom_tile(aes(fill = Frequency),colour = "white")+
  scale_fill_gradient2(low = "#FFFFFF",
                       high = "red") +
  #scale_fill_gradient(low = "blue", high = "red")+
  theme_dark()+
  theme(axis.text.x = element_text(size=6, hjust = 1,angle=45),
        axis.text.y = element_text(size=5, hjust = 1,face="italic"))+
  facet_grid(~ model, switch = "x", scales = "free_y", space = "free_x")+
  ggtitle("Genome feature distribution")

# Genomics distribution of peaks
distfreqbcat<-rbind(DB_freq,DBL_freq,RBCS_freq,RBLS_freq)

distfreqbcat$Feature<-gsub('Promoter.*','Promoter',distfreqbcat$Feature)
distfreqbcat$Feature<-gsub('.*Exon.*','Exon',distfreqbcat$Feature)
distfreqbcat$Feature<-gsub('.*Intron.*','Intron',distfreqbcat$Feature)

distfreqbcat$groups<-paste(distfreqbcat$group,distfreqbcat$model,sep="_")

distfreqbcat<-data.frame(aggregate(Frequency ~ Feature + groups, data=distfreqbcat, FUN=sum))

distfreqbcat$groups<-factor(distfreqbcat$groups,levels = c("bcat_human_control","bcat_SC_human_control","bcat_human_Li","bcat_SC_human_Li"))



ggplot(distfreqbcat,aes(x=groups,y=Frequency,fill=Feature))+
  geom_bar(stat="identity",color="black") +
  scale_fill_brewer(palette="Pastel1")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")


## Common genes

# Cell lines
DB_TSS<-as.data.frame(cbind(as.vector(peakAnnoDB_black@anno$gene),peakAnnoDB_black@anno$distanceToTSS))
colnames(DB_TSS)<-c("Gene","DB_distance_to_TSS")

DBL_TSS<-as.data.frame(cbind(as.vector(peakAnnoDBL_black@anno$gene),peakAnnoDBL_black@anno$distanceToTSS))
colnames(DBL_TSS)<-c("Gene","DBL_distance_to_TSS")

RBCS_TSS<-as.data.frame(cbind(as.vector(peakAnnoRBCS_black@anno$gene),peakAnnoRBCS_black@anno$distanceToTSS))
colnames(RBCS_TSS)<-c("Gene","RBCS_distance_to_TSS")

RBLS_TSS<-as.data.frame(cbind(as.vector(peakAnnoRBLS_black@anno$gene),peakAnnoRBLS_black@anno$distanceToTSS))
colnames(RBLS_TSS)<-c("Gene","RBLS_distance_to_TSS")

TCF_TSS<-as.data.frame(cbind(as.vector(peakAnnoTCF_black@anno$gene),peakAnnoTCF_black@anno$distanceToTSS))
colnames(TCF_TSS)<-c("Gene","TCF_distance_to_TSS")

LEF_TSS<-as.data.frame(cbind(as.vector(peakAnnoLEF_black@anno$gene),peakAnnoLEF_black@anno$distanceToTSS))
colnames(LEF_TSS)<-c("Gene","LEF_distance_to_TSS")

RKT_TSS<-as.data.frame(cbind(as.vector(peakAnnoRKT_black@anno$gene),peakAnnoRKT_black@anno$distanceToTSS))
colnames(RKT_TSS)<-c("Gene","RKT_distance_to_TSS")

RKB_TSS<-as.data.frame(cbind(as.vector(peakAnnoRKB_black@anno$gene),peakAnnoRKB_black@anno$distanceToTSS))
colnames(RKB_TSS)<-c("Gene","RKB_distance_to_TSS")

human_tss<-merge(DB_TSS,DBL_TSS,by="Gene",all=TRUE)
human_tss<-merge(human_tss,RBCS_TSS,by="Gene",all=TRUE)
human_tss<-merge(human_tss,RBLS_TSS,by="Gene",all=TRUE)
human_tss<-merge(human_tss,TCF_TSS,by="Gene",all=TRUE)
human_tss<-merge(human_tss,LEF_TSS,by="Gene",all=TRUE)
human_tss<-merge(human_tss,RKT_TSS,by="Gene",all=TRUE)
human_tss<-merge(human_tss,RKB_TSS,by="Gene",all=TRUE)



# Mouse chips
BMLB_TSS<-as.data.frame(cbind(as.vector(peakAnnoBMLB_black@anno$gene),peakAnnoBMLB_black@anno$distanceToTSS))
colnames(BMLB_TSS)<-c("Gene","BMLB_distance_to_TSS")
BMLT_TSS<-as.data.frame(cbind(as.vector(peakAnnoBMLT_black@anno$gene),peakAnnoBMLT_black@anno$distanceToTSS))
colnames(BMLT_TSS)<-c("Gene","BMLT_distance_to_TSS")
FI_TSS<-as.data.frame(cbind(as.vector(peakAnnoFI_black@anno$gene),peakAnnoFI_black@anno$distanceToTSS))
colnames(FI_TSS)<-c("Gene","FI_distance_to_TSS")
FW_TSS<-as.data.frame(cbind(as.vector(peakAnnoFW_black@anno$gene),peakAnnoFW_black@anno$distanceToTSS))
colnames(FW_TSS)<-c("Gene","FW_distance_to_TSS")

# Human orthologs
BMLB_human_ort<-getBM(attributes=c("hsapiens_homolog_associated_gene_name", "external_gene_name","description"),
                      values= unique(BMLB_TSS$Gene),
                      filter="external_gene_name",
                      mart=mm)

BMLT_human_ort<-getBM(attributes=c("hsapiens_homolog_associated_gene_name", "external_gene_name","description"),
                      values= unique(BMLT_TSS$Gene),
                      filter="external_gene_name",
                      mart=mm)

FI_human_ort<-getBM(attributes=c("hsapiens_homolog_associated_gene_name", "external_gene_name","description"),
                    values= unique(FI_TSS$Gene),
                    filter="external_gene_name",
                    mart=mm)

FW_human_ort<-getBM(attributes=c("hsapiens_homolog_associated_gene_name", "external_gene_name","description"),
                    values= unique(FW_TSS$Gene),
                    filter="external_gene_name",
                    mart=mm)

BMLB_TSS<-merge(BMLB_TSS,BMLB_human_ort,by.x="Gene",by.y="external_gene_name")
BMLB_TSS<-BMLB_TSS[,c(3,2)]
BMLB_TSS<-BMLB_TSS[BMLB_TSS$hsapiens_homolog_associated_gene_name!="",]
colnames(BMLB_TSS)<-c("Gene","BMLB_distance_to_TSS")

BMLT_TSS<-merge(BMLT_TSS,BMLT_human_ort,by.x="Gene",by.y="external_gene_name")
BMLT_TSS<-BMLT_TSS[,c(3,2)]
BMLT_TSS<-BMLT_TSS[BMLT_TSS$hsapiens_homolog_associated_gene_name!="",]
colnames(BMLT_TSS)<-c("Gene","BMLT_distance_to_TSS")

FI_TSS<-merge(FI_TSS,FI_human_ort,by.x="Gene",by.y="external_gene_name")
FI_TSS<-FI_TSS[,c(3,2)]
FI_TSS<-FI_TSS[FI_TSS$hsapiens_homolog_associated_gene_name!="",]
colnames(FI_TSS)<-c("Gene","FI_distance_to_TSS")

FW_TSS<-merge(FW_TSS,FW_human_ort,by.x="Gene",by.y="external_gene_name")
FW_TSS<-FW_TSS[,c(3,2)]
FW_TSS<-FW_TSS[FW_TSS$hsapiens_homolog_associated_gene_name!="",]
colnames(FW_TSS)<-c("Gene","FW_distance_to_TSS")

## Merge human and mice distance info
mam_tss<-merge(human_tss,BMLB_TSS,by="Gene",all=TRUE)
mam_tss<-merge(mam_tss,BMLT_TSS,by="Gene",all=TRUE)
mam_tss<-merge(mam_tss,FI_TSS,by="Gene",all=TRUE)
mam_tss<-merge(mam_tss,FW_TSS,by="Gene",all=TRUE)

# Common genes
intersect(unique(DB_TSS$Gene),unique(BMLB_TSS$Gene))
intersect(unique(RBCS_TSS$Gene),unique(BMLB_TSS$Gene))

intersect(unique(DBL_TSS$Gene),unique(BMLB_TSS$Gene))
intersect(unique(RBLS_TSS$Gene),unique(BMLB_TSS$Gene))

totli<-c(intersect(unique(DB_TSS$Gene),unique(BMLB_TSS$Gene)),intersect(unique(DBL_TSS$Gene),unique(BMLB_TSS$Gene)))
unique(totli)

intersect(unique(BMLT_TSS$Gene),unique(BMLB_TSS$Gene))

intersect(unique(DBL_TSS$Gene),unique(BMLB_TSS$Gene))


## Subset data to correlation
#chen NAs and select columns to be represented
subpoint<-mam_tss[!is.na(mam_tss$BMLB_distance_to_TSS) & !is.na(mam_tss$RBCS_distance_to_TSS) ,]

#select columns
subpoint<-subpoint[,c(1,4,10)][!duplicated(subpoint[,c(1,4,10)]), ]



ggplot(subpoint, aes(x=as.numeric(as.character(RBCS_distance_to_TSS))/1000,y=as.numeric(as.character(BMLB_distance_to_TSS))/1000,label=Gene))+
  geom_vline(xintercept = 0,color="dodgerblue")+
  geom_hline(yintercept = 0,color="dodgerblue")+
  geom_point(size=3) +
  #geom_density_2d()+
  geom_label_repel(data=subpoint[abs(as.numeric(as.character(subpoint$BMLB_distance_to_TSS)))<20 | abs(as.numeric(as.character(subpoint$RBCS_distance_to_TSS)))<20,],aes(label=Gene), size=3,fontface = "italic",alpha=0.5)+
  theme_bw()




### MEMEChip analysis ##
summary_all<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/MEMEChip/summary_all.tsv",row.names = NULL)
#remove lines with # and NAs and EMPTY
summary_all<-summary_all[!is.na(summary_all$E.VALUE),]
summary_all<-summary_all[!grepl('#',summary_all$MOTIF_INDEX),]
colnames(summary_all)[1]<-"SAMPLE"

summary_all$SAMPLE<-factor(summary_all$SAMPLE,levels = c("BMLB","BMLT","FL_WT","FL_NOTCH","RBCS","RBLS","DB","DBL","TCF","LEF","RKB","RKT"))
summary_all$IDmotif<-paste(summary_all$MOTIF_ID,summary_all$MOST_SIMILAR_MOTIF,sep="_")

# motif ids identifed with jolma
ggplot(summary_all[summary_all$MOTIF_SOURCE=="db/EUKARYOTE/jolma2013.meme",],aes(x=SAMPLE,y=MOTIF_ID))+
  geom_point(aes(color=-log10(E.VALUE)),size=5)+
  facet_grid(~MOTIF_SOURCE)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=5),
        axis.text.y = element_text(size=5))

# motif ids JASPAR and MOUSE
ggplot(summary_all[summary_all$MOTIF_SOURCE=="db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme" | summary_all$MOTIF_SOURCE=="db/MOUSE/uniprobe_mouse.meme",],aes(x=SAMPLE,y=ALT_ID))+
  geom_point(aes(color=-log10(E.VALUE)),size=5)+
  facet_grid(~MOTIF_SOURCE)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=5),
        axis.text.y = element_text(size=5))

# motif ids SIMILAR identifed MEME and DREME
ggplot(summary_all[(summary_all$MOTIF_SOURCE=="MEME" | summary_all$MOTIF_SOURCE=="DREME") & summary_all$MOST_SIMILAR_MOTIF!=" ",],aes(x=SAMPLE,y=MOST_SIMILAR_MOTIF))+
  geom_point(aes(color=-log10(E.VALUE)),size=5)+
  facet_grid(~MOTIF_SOURCE)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=5),
        axis.text.y = element_text(size=5))

# motif ids SIMILAR identifed MEME and DREME not mice, not KOs
summary_all_noKO<-summary_all[summary_all$SAMPLE!="BMLB" & summary_all$SAMPLE!="BMLT" &
                              summary_all$SAMPLE!="FL_WT" & summary_all$SAMPLE!="FL_NOTCH" &
                              summary_all$SAMPLE!="RKB" & summary_all$SAMPLE!="RKT", ]
ggplot(summary_all_noKO[(summary_all_noKO$MOTIF_SOURCE=="MEME" | summary_all_noKO$MOTIF_SOURCE=="DREME") & summary_all_noKO$MOST_SIMILAR_MOTIF!=" ",],aes(x=SAMPLE,y=MOST_SIMILAR_MOTIF))+
  geom_point(color="black",size=7)+
  geom_point(data=summary_all_noKO[(summary_all_noKO$MOTIF_SOURCE=="MEME" | summary_all_noKO$MOTIF_SOURCE=="DREME") & summary_all_noKO$MOST_SIMILAR_MOTIF!=" " & summary_all_noKO$E.VALUE>0,],aes(color=-log10(E.VALUE)),size=6)+
  scale_color_gradient2(low="white",high="dodgerblue")+
  geom_point(data=summary_all_noKO[(summary_all_noKO$MOTIF_SOURCE=="MEME" | summary_all_noKO$MOTIF_SOURCE=="DREME") & summary_all_noKO$MOST_SIMILAR_MOTIF!=" " & summary_all_noKO$E.VALUE==0,],color="dodgerblue",size=6)+
  facet_grid(~MOTIF_SOURCE)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=9),
        axis.text.y = element_text(size=7),legend.position = "bottom")

# motif ids SIMILAR identifed MEME and DREME not mice, not KOs, only RD bcat
summary_all_noKO<-summary_all[summary_all$SAMPLE!="BMLB" & summary_all$SAMPLE!="BMLT" &
                                summary_all$SAMPLE!="FL_WT" & summary_all$SAMPLE!="FL_NOTCH" &
                                summary_all$SAMPLE!="RKB" & summary_all$SAMPLE!="RKT" &
                                summary_all$SAMPLE!="RBCS" & summary_all$SAMPLE!="RBLS" &
                                summary_all$SAMPLE!="TCF" & summary_all$SAMPLE!="LEF", ]

ggplot(summary_all_noKO[(summary_all_noKO$MOTIF_SOURCE=="MEME" | summary_all_noKO$MOTIF_SOURCE=="DREME") & summary_all_noKO$MOST_SIMILAR_MOTIF!=" ",],aes(x=SAMPLE,y=MOST_SIMILAR_MOTIF))+
  geom_point(color="black",size=7)+
  geom_point(data=summary_all_noKO[(summary_all_noKO$MOTIF_SOURCE=="MEME" | summary_all_noKO$MOTIF_SOURCE=="DREME") & summary_all_noKO$MOST_SIMILAR_MOTIF!=" " & summary_all_noKO$E.VALUE>0,],aes(color=-log10(E.VALUE)),size=6)+
  scale_color_gradient2(low="white",high="dodgerblue")+
  geom_point(data=summary_all_noKO[(summary_all_noKO$MOTIF_SOURCE=="MEME" | summary_all_noKO$MOTIF_SOURCE=="DREME") & summary_all_noKO$MOST_SIMILAR_MOTIF!=" " & summary_all_noKO$E.VALUE==0,],color="dodgerblue",size=6)+
  facet_grid(~MOTIF_SOURCE)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=9),
        axis.text.y = element_text(size=7),legend.position = "bottom")


# motif ids not SIMILAR identifed MEME and DREME
ggplot(summary_all[(summary_all$MOTIF_SOURCE=="MEME" | summary_all$MOTIF_SOURCE=="DREME") & summary_all$MOST_SIMILAR_MOTIF==" ",],aes(x=SAMPLE,y=MOTIF_ID))+
  geom_point(aes(color=-log10(E.VALUE)),size=5)+
  facet_grid(~MOTIF_SOURCE)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=5),
        axis.text.y = element_text(size=5))


## ENRICHR
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "InterPro_Domains_2019",
         "GO_Cellular_Component_2017",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2016")



geneRBCS <- enrichr(as.character(unique(peakAnnoRBCS_black@anno$SYMBOL)), dbs)
geneRBCS<- geneRBCS[["GO_Biological_Process_2018"]]
geneRBCS$group<-"bcat_RBCS"

geneRBLS <- enrichr(as.character(unique(peakAnnoRBLS_black@anno$SYMBOL)), dbs)
geneRBLS<- geneRBLS[["GO_Biological_Process_2018"]]
geneRBLS$group<-"bcat_RBLS"

geneDB <- enrichr(as.character(unique(peakAnnoDB_black@anno$SYMBOL)), dbs)
geneDB<- geneDB[["GO_Biological_Process_2018"]]
geneDB$group<-"bcat_DB"

geneDBL <- enrichr(as.character(unique(peakAnnoDBL_black@anno$SYMBOL)), dbs)
geneDBL<- geneDBL[["GO_Biological_Process_2018"]]
geneDBL$group<-"bcat_DBL"

common_RBCS_DB<-enrichr(as.character(unique(com_RBCS_DB$SYMBOL)), dbs)
common_RBCS_DB<- common_RBCS_DB[["GO_Biological_Process_2018"]]
common_RBCS_DB$group<-"common_RBCS_DB"

common_RBLS_DBL<-enrichr(as.character(unique(com_RBLS_DBL$SYMBOL)), dbs)
common_RBLS_DBL<- common_RBLS_DBL[["GO_Biological_Process_2018"]]
common_RBLS_DBL$group<-"common_RBLS_DBL"


allGO<-rbind(common_RBCS_DB,common_RBLS_DBL)

allGO<-geneDB

allGO<-rbind(geneDB,geneDBL,geneRBCS,geneRBLS)

bpsub<-subset(allGO,allGO$Adjusted.P.value<=0.01)

bpsub<-subset(allGO,allGO$P.value<=0.01)

bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=unique((bpsub$Term)[rev(order(bpsub$group,bpsub$P.value))]))

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$group,bpsub$P.value)), ]$Term))

ggplot(bpsub,aes(y=-log10(P.value),x=Term,fill=group),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=group),color="black")+
  #geom_text(aes(label=tolower(Genes)), size=3,hjust=1,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Spectral")+
 # facet_wrap(~group,scales="free")+
  theme_minimal()+
  coord_flip()+
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=10, face="bold",hjust = 1),
        axis.text.y = element_text(size=10, hjust = 1),
        legend.position = "bottom")+
  facet_grid(~group,space = "free")

write.table(bpsub,'/Users/yguillen/Desktop/table_BP_enrich_SC_RD_bcat_padj005.txt',quote = FALSE,sep = "\t",row.names = FALSE)

#### MEME ChIP motif enrichment for genes involved in RNA processing or splicing in bcat genes from SC antibody

hg<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

#merging lithium and control
SCabbm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                values= unique(c(unique(peakAnnoRBCS_black@anno$SYMBOL),unique(peakAnnoRBLS_black@anno$SYMBOL))),
                filter="external_gene_name",
                mart=hg)

SCabbm$go_id<-gsub("^$", NA,SCabbm$go_id)
SCabbm<-merge(SCabbm,vasbp,by="go_id",all.x=TRUE)

SCabbm_rna<-SCabbm[grepl('splic',SCabbm$Term) | grepl('RNA processing',SCabbm$Term) | grepl('RNA metabolic process',SCabbm$Term),]
sort(unique(SCabbm_rna$external_gene_name))

# Select peaks annotated close to these genes
rnapeaksRBCS<-as.data.frame(peakAnnoRBCS_black@anno[peakAnnoRBCS_black@anno$SYMBOL %in% unique(SCabbm_rna$external_gene_name),])
write.table(rnapeaksRBCS,'/Users/yguillen/Desktop/temp/beta_catenin_project/MEMEChip/RNA_GO_peaks_RBCS.tab',quote = FALSE,sep = "\t",row.names = FALSE)

rnapeaksRBLS<-as.data.frame(peakAnnoRBLS_black@anno[peakAnnoRBLS_black@anno$SYMBOL %in% unique(SCabbm_rna$external_gene_name),])
write.table(rnapeaksRBLS,'/Users/yguillen/Desktop/temp/beta_catenin_project/MEMEChip/RNA_GO_peaks_RBLS.tab',quote = FALSE,sep = "\t",row.names = FALSE)


#Comparing distribution of RNA peaks and all peaks
#RBCS
rnapeaksRBCS$annotation<-gsub(' \\(u.*','',rnapeaksRBCS$annotation)
brewer.pal(9, "Pastel1")

dataf_peakannoRBCS<-as.data.frame(peakAnnoRBCS_black@anno)
dataf_peakannoRBCS$annotation<-gsub(' \\(u.*','',dataf_peakannoRBCS$annotation)


#RBLS
rnapeaksRBLS$annotation<-gsub(' \\(u.*','',rnapeaksRBLS$annotation)

dataf_peakannoRBLS<-as.data.frame(peakAnnoRBLS_black@anno)
dataf_peakannoRBLS$annotation<-gsub(' \\(u.*','',dataf_peakannoRBLS$annotation)

rna_all<-rbind(
  data.frame(prop.table(table(rnapeaksRBCS$annotation)),samp="RNA_RBCS"),
  data.frame(prop.table(table(dataf_peakannoRBCS$annotation)),samp="RBCS"),
  data.frame(prop.table(table(rnapeaksRBLS$annotation)),samp="RNA_RBLS"),
  data.frame(prop.table(table(dataf_peakannoRBLS$annotation)),samp="RBLS"))

cols=(c(brewer.pal(9, "Pastel1"),"grey","dodgerblue"))

ggplot(rna_all,aes(x=samp,y=Freq*100))+
  geom_bar(stat="identity",aes(fill=Var1))+
  scale_fill_manual(values =cols)+
  theme_bw()

## Overall distribution of peaks
peakall<-rbind(data.frame(peakAnnoDB1@annoStat,samp="DB1"),
      data.frame(peakAnnoRBL2@annoStat,samp="RBL2"),
      data.frame(peakAnnoRBCS1@annoStat,samp="RBCS1"),
      data.frame(peakAnnoRBLS1@annoStat,samp="RBLS1"))

table(unique(peakall$Feature))

cols=(c(brewer.pal(9, "Pastel1"),"grey","darkseagreen1"))

ggplot(peakall,aes(x=samp,y=Frequency))+
  geom_bar(stat="identity",aes(fill=Feature))+
  scale_fill_manual(values =cols)+
  theme_bw()



## Annotate common peaks SC and RD obtained with script random_peaks.sh
# Control
com_RBCS_DB<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/RBCS_DB_common_sets_d100.bed",header = FALSE)
dataf_peakannoDB_black<-as.data.frame(peakAnnoDB_black@anno)
com_RBCS_DB<-merge(com_RBCS_DB,dataf_peakannoDB_black,by.x="V8",by.y="Name",all.x=TRUE)

prop.table(table(com_RBCS_DB$annotation))

# Lithium
com_RBLS_DBL<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/RBLS_DBL_common_sets_d100.bed",header = FALSE)
dataf_peakannoDBL_black<-as.data.frame(peakAnnoDBL_black@anno)
com_RBLS_DBL<-merge(com_RBLS_DBL,dataf_peakannoDBL_black,by.x="V8",by.y="Name",all.x=TRUE)

prop.table(table(com_RBLS_DBL$annotation))



#### multibigwig correlations between RD and SC CHIPS

RDSC_comp<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/multibigwig/scores_DB1_RBCS1_RBL2_RBLS1_compare.tab")
colnames(RDSC_comp)<-gsub('X\\.','',colnames(RDSC_comp))
colnames(RDSC_comp)<-gsub('\\.','',colnames(RDSC_comp))
colnames(RDSC_comp)



#Remove all zeros
RDSC_filt<-row.names(RDSC_comp[RDSC_comp$DB1==0 & RDSC_comp$RBCS1==0 & RDSC_comp$RBL2==0 & RDSC_comp$RBLS1==0,])
RDSC_comp<-RDSC_comp[!row.names(RDSC_comp) %in% RDSC_filt,]

summary(RDSC_comp$DB1)
summary(RDSC_comp$RBCS1)
summary(RDSC_comp$RBL2)
summary(RDSC_comp$RBLS1)

## filter outliers that are input peaks likely too
RDSC_comp_filt<-RDSC_comp[RDSC_comp$DB1<=100 & RDSC_comp$RBCS1<=100 & RDSC_comp$RBL2<=100 & RDSC_comp$RBLS1<=100,]


ggplot(data=RDSC_comp_filt,aes(x=DB1,y=RBCS1))+
  geom_point(color="black",size=1,alpha=0.3)+
  #geom_point(data=cormul_filt_model[cormul_filt_model$CWT1>cormul_filt_model$upr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO",],color="dodgerblue",size=1,alpha=0.8)+
  #geom_point(data=cormul_filt_model[cormul_filt_model$CWT1<cormul_filt_model$lwr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO",],color="coral",size=1,alpha=0.8)+
  theme_bw()
