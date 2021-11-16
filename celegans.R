source("https://bioconductor.org/biocLite.R")

# Reference genome C. elegans
biocLite("BSgenome.Celegans.UCSC.ce10")
biocLite("TxDb.Celegans.UCSC.ce11.refGene")
BiocManager::install("TxDb.Celegans.UCSC.ce11.ensGene")
biocLite("TxDb.Celegans.UCSC.ce6.refGene")
biocLite("org.Ce.eg.db")
biocLite("trackViewer")

#Load libraries
library(data.table)
library(GenomicAlignments)
library(GO.db)
library(DiffBind)
library(ChIPQC)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(BayesPeak)
library(parallel)
library(ChIPseeker)
library(clusterProfiler)
library(BSgenome.Celegans.UCSC.ce10)
library(ReactomePA)
library(org.Ce.eg.db)
library(trackViewer)
library(Gviz)
library(rtracklayer)
library(Sushi)
library(biomaRt)
library(TxDb.Celegans.UCSC.ce11.refGene)
library(TxDb.Celegans.UCSC.ce11.ensGene)

setwd("/Volumes/cancer/Celegans/ChIPSeq/peakcall_new/")

################## MACS2 PEAK CALLING with ChIPseeker ###############

## Gene germ classification
gaydos<-read_xlsx("/Volumes/grcmc/YGUILLEN/Celegans/gene_lists_ref/germline_vs_somatic_vs_ubiquitous_genes_Gaydos.xlsx")


# L4 stage
peakFMM1 <- readPeakFile("FMM1_peaks_peaks.broadPeak", header=FALSE)
peakFMM2 <- readPeakFile("FMM2_peaks_peaks.broadPeak", header=FALSE)
peakFMM3 <- readPeakFile("FMM3_peaks_peaks.broadPeak", header=FALSE)

peakFMF1<-readPeakFile("FMF1_peaks_peaks.broadPeak",header=FALSE)
peakFMF2<-readPeakFile("FMF2_peaks_peaks.broadPeak",header=FALSE)

peakHK27DKO <- readPeakFile("H3K27_DKO_peaks_peaks.broadPeak", header=FALSE)
peakHK27ikb <- readPeakFile("H3K27_ikb-1KO_peaks_peaks.broadPeak", header=FALSE)
peakHK27nfk <- readPeakFile("H3K27_nfki-1KO_peaks_peaks.broadPeak",header=FALSE)
peakHK27WT <- readPeakFile("H3K27_WT_peaks_peaks.broadPeak",header=FALSE)

peakHK36DKO <- readPeakFile("H3K36_DKO_peaks_peaks.broadPeak",header=FALSE)
peakHK36ikb <- readPeakFile("H3K36_ikb-1KO_peaks_peaks.broadPeak",header=FALSE)
peakHK36nfk <- readPeakFile("H3K36_nfki-1KO_peaks_peaks.broadPeak",header=FALSE)
peakHK36WT <- readPeakFile("H3K36_WT_peaks_peaks.broadPeak",header=FALSE)

# L1 stage
peakIK27 <- readPeakFile("IK27_peaks_peaks.broadPeak",header=FALSE)
peakIK36 <- readPeakFile("IK36_peaks_peaks.broadPeak",header=FALSE)
peakNF27 <- readPeakFile("NF27_peaks_peaks.broadPeak",header=FALSE)
peakNF36 <- readPeakFile("NF36_peaks_peaks.broadPeak",header=FALSE)
peakWT27 <- readPeakFile("WT27_peaks_peaks.broadPeak",header=FALSE)
#no peaks for WT36
#peakWT36 <- readPeakFile("WT36_peaks_peaks.broadPeak",header=FALSE)

#Overall all chromosmes peak coverage
covplot(peakHK36DKO, weightCol="V5")
covplot(peakFMM2, weightCol="V5")
covplot(peakFMM3, weightCol="V5")
covplot(peakWT27, weightCol="V5")

#genome <- BSgenome.Celegans.UCSC.ce10
#seqlengths(genome)

#Annotation form GFF

txdb<-makeTxDbFromGFF(file = "/Volumes/cancer/db_files/Caenorhabditis_elegans/UCSC/ce10/Annotation/Archives/archive-2015-07-17-14-29-29/Genes/genes.gff3",
                      organism = "Caenorhabditis elegans")


## annotation from UCSC

#library(rtracklayer)
#ucscGenomes()[ , "db"]
#supportedUCSCtables()
#txdb <- makeTxDbFromUCSC(genome = "ce10", tablename = "ensGene")


## Profile of ChIP peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

#tagMatrixFMM <- getTagMatrix(peakFMM1, windows=promoter)
#tagMatrixFMF <- getTagMatrix(peakFMF1, windows=promoter)

## to speed up the compilation of this vignettes, we use a precalculated tagMatrix
#data("tagMatrixFMM")
#tagMatrix <- tagMatrixList[[4]]

## Heatmap TSS 
#tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

## Average profile of ChIP peaks binding to TSS region ##

#plotAvgProf(tagMatrix, xlim=c(-3000, 3000),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

### PEAK ANNOTATION
#FMM
peakAnnoFMM1 <- annotatePeak(peakFMM1, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoFMM2 <- annotatePeak(peakFMM2, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoFMM3 <- annotatePeak(peakFMM3, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")

#FMF
peakAnnoFMF1 <- annotatePeak(peakFMF1, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoFMF2 <- annotatePeak(peakFMF2, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")

#H3K27
peakAnnoHK27DKO <- annotatePeak(peakHK27DKO, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoHK27ikb <- annotatePeak(peakHK27ikb, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoHK27nfk <- annotatePeak(peakHK27nfk, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoHK27WT <- annotatePeak(peakHK27WT, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")


#H3K36
peakAnnoHK36DKO <- annotatePeak(peakHK36DKO, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoHK36ikb <- annotatePeak(peakHK36ikb, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoHK36nfk <- annotatePeak(peakHK36nfk, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoHK36WT <- annotatePeak(peakHK36WT, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")

#L1 H3K27
peakAnnoIK27 <- annotatePeak(peakIK27, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoNF27 <- annotatePeak(peakNF27, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoWT27 <- annotatePeak(peakWT27, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")

#L1 H3K36
peakAnnoIK36 <- annotatePeak(peakIK36, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
peakAnnoNF36 <- annotatePeak(peakNF36, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")
# no peaks for WT36
#peakAnnoWT36 <- annotatePeak(peakWT36, tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Ce.eg.db")


# Visualize genomic annotation
plotAnnoPie(peakAnnoNF36)
dev.off()

plotAnnoBar(peakAnnoFMM1)
dev.off()

vennpie(peakAnnoFMM1)
dev.off()


plotDistToTSS(peakAnnoNF36,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


##DATAFRAME WITH genes id lnked to peaks
dataf_peakannoFMM1<-as.data.frame(peakAnnoFMM1)
dataf_peakannoFMM2<-as.data.frame(peakAnnoFMM2)
dataf_peakannoFMM3<-as.data.frame(peakAnnoFMM3)

dataf_peakannoFMF1<-as.data.frame(peakAnnoFMF1)
dataf_peakannoFMF2<-as.data.frame(peakAnnoFMF2)

dataf_peakannoHK27DKO<-as.data.frame(peakAnnoHK27DKO)
dataf_peakannoHK27ikb<-as.data.frame(peakAnnoHK27ikb)
dataf_peakannoHK27nfk<-as.data.frame(peakAnnoHK27nfk)
dataf_peakannoHK27WT<-as.data.frame(peakAnnoHK27WT)

dataf_peakannoHK36DKO<-as.data.frame(peakAnnoHK36DKO)
dataf_peakannoHK36ikb<-as.data.frame(peakAnnoHK36ikb)
dataf_peakannoHK36nfk<-as.data.frame(peakAnnoHK36nfk)
dataf_peakannoHK36WT<-as.data.frame(peakAnnoHK36WT)

#L1
dataf_peakannoIK27<-as.data.frame(peakAnnoIK27)
dataf_peakannoNF27<-as.data.frame(peakAnnoNF27)
dataf_peakannoWT27<-as.data.frame(peakAnnoWT27)

dataf_peakannoIK36<-as.data.frame(peakAnnoIK36)
dataf_peakannoNF36<-as.data.frame(peakAnnoNF36)
#dataf_peakannoWT36<-as.data.frame(peakAnnoWT36)


write.table(dataf_peakannoIK27,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_L1_IKB_H27.tab",quote = FALSE,row.names = FALSE)
write.table(dataf_peakannoNF27,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_L1_NFKI_H27.tab",quote = FALSE,row.names = FALSE)
write.table(dataf_peakannoWT27,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_L1_WT_H27.tab",quote = FALSE,row.names = FALSE)

write.table(dataf_peakannoIK36,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_L1_IKB_H36.tab",quote = FALSE,row.names = FALSE)
write.table(dataf_peakannoNF36,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_L1_NFKI_H36.tab",quote = FALSE,row.names = FALSE)

#list genes

FMM1_genes<-unique(dataf_peakannoFMM1$geneId)
FMM2_genes<-unique(dataf_peakannoFMM2$geneId)
FMM3_genes<-unique(dataf_peakannoFMM3$geneId)

FMF1_genes<-unique(dataf_peakannoFMF1$geneId)
FMF2_genes<-unique(dataf_peakannoFMF2$geneId)

HK27DKO_genes<-unique(dataf_peakannoHK27DKO$geneId)
HK27Dikb_genes<-unique(dataf_peakannoHK27ikb$geneId)
HK27Dnfk_genes<-unique(dataf_peakannoHK27nfk$geneId)
HK27DWT_genes<-unique(dataf_peakannoHK27WT$geneId)

IK27_genes<-unique(dataf_peakannoIK27$geneId)
IK36_genes<-unique(dataf_peakannoIK36$geneId)
NF27_genes<-unique(dataf_peakannoNF27$geneId)
NF36_genes<-unique(dataf_peakannoNF36$geneId)
WT27_genes<-unique(dataf_peakannoWT27$geneId)
#WT36_genes<-unique(dataf_peakannoWT36$geneId)


HK36DKO_genes<-unique(dataf_peakannoHK36DKO$geneId)
HK36Dikb_genes<-unique(dataf_peakannoHK36ikb$geneId)
HK36Dnfk_genes<-unique(dataf_peakannoHK36nfk$geneId)
HK36DWT_genes<-unique(dataf_peakannoHK36WT$geneId)


comgenesFMF<-as.data.frame(table(sort(c(FMF1_genes,FMF2_genes))))
comgenesFMF<-comgenesFMF[order(comgenesFMF$Freq,decreasing = TRUE),]

comgenesFMM<-as.data.frame(table(sort(c(FMM1_genes,FMM2_genes,FMM3_genes))))
comgenesFMM<-comgenesFMM[order(comgenesFMM$Freq,decreasing = TRUE),]


#Tables to plot
tabgenesFMF<-as.data.frame(table(comgenesFMF$Freq))
tabgenesFMM<-as.data.frame(table(comgenesFMM$Freq))


write.table(comgenesFMF$Var1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/all_peaks_FMF.txt",quote = FALSE,row.names = FALSE)
write.table(comgenesFMM$Var1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/all_peaks_FMM.txt",quote = FALSE,row.names = FALSE)



library(ggpubr)

ggbarplot(tabgenesFMF, y = "Freq", x = "Var1",
          fill="steelblue",alpha=0.4,
          label = TRUE, lab.pos = "in", lab.col = "white",lab.size = 9,
          xlab = "Number of replicates", ylab="Number of genes")

ggbarplot(tabgenesFMM, y = "Freq", x = "Var1",
          fill="steelblue",alpha=0.4,
          label = TRUE, lab.pos = "in", lab.col = "white",lab.size = 9,
          xlab = "Number of replicates", ylab="Number of genes")


## PEAK DISTRIBUTION


dataf_peakannoFMF1<-as.data.frame(peakAnnoFMF1)
dataf_peakannoFMF1$annotation<-gsub('Intron.*','Intron',dataf_peakannoFMF1$annotation)
dataf_peakannoFMF1$annotation<-gsub('Exon.*','Exon',dataf_peakannoFMF1$annotation)
dataf_peakannoFMF1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFMF1$annotation)
dataf_peakannoFMF1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFMF1$annotation)
dataf_peakannoFMF1$peaksite<-paste(dataf_peakannoFMF1$geneId,dataf_peakannoFMF1$annotation,sep='_')
table(dataf_peakannoFMF1$annotation)


dataf_peakannoFMF2<-as.data.frame(peakAnnoFMF2)
dataf_peakannoFMF2$annotation<-gsub('Intron.*','Intron',dataf_peakannoFMF2$annotation)
dataf_peakannoFMF2$annotation<-gsub('Exon.*','Exon',dataf_peakannoFMF2$annotation)
dataf_peakannoFMF2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFMF2$annotation)
dataf_peakannoFMF2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFMF2$annotation)
dataf_peakannoFMF2$peaksite<-paste(dataf_peakannoFMF2$geneId,dataf_peakannoFMF2$annotation,sep='_')
table(dataf_peakannoFMF2$annotation)

dataf_peakannoFMM1<-as.data.frame(peakAnnoFMM1)
dataf_peakannoFMM1$annotation<-gsub('Intron.*','Intron',dataf_peakannoFMM1$annotation)
dataf_peakannoFMM1$annotation<-gsub('Exon.*','Exon',dataf_peakannoFMM1$annotation)
dataf_peakannoFMM1$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFMM1$annotation)
dataf_peakannoFMM1$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFMM1$annotation)
dataf_peakannoFMM1$peaksite<-paste(dataf_peakannoFMM1$geneId,dataf_peakannoFMM1$annotation,sep='_')
table(dataf_peakannoFMM1$annotation)


dataf_peakannoFMM2<-as.data.frame(peakAnnoFMM2)
dataf_peakannoFMM2$annotation<-gsub('Intron.*','Intron',dataf_peakannoFMM2$annotation)
dataf_peakannoFMM2$annotation<-gsub('Exon.*','Exon',dataf_peakannoFMM2$annotation)
dataf_peakannoFMM2$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFMM2$annotation)
dataf_peakannoFMM2$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFMM2$annotation)
dataf_peakannoFMM2$peaksite<-paste(dataf_peakannoFMM2$geneId,dataf_peakannoFMM2$annotation,sep='_')
table(dataf_peakannoFMM2$annotation)

dataf_peakannoFMM3<-as.data.frame(peakAnnoFMM3)
dataf_peakannoFMM3$annotation<-gsub('Intron.*','Intron',dataf_peakannoFMM3$annotation)
dataf_peakannoFMM3$annotation<-gsub('Exon.*','Exon',dataf_peakannoFMM3$annotation)
dataf_peakannoFMM3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFMM3$annotation)
dataf_peakannoFMM3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFMM3$annotation)
dataf_peakannoFMM3$peaksite<-paste(dataf_peakannoFMM3$geneId,dataf_peakannoFMM3$annotation,sep='_')
table(dataf_peakannoFMM3$annotation)


dataf_peakannoHK27DKO<-as.data.frame(peakAnnoHK27DKO)
dataf_peakannoHK27DKO$annotation<-gsub('Intron.*','Intron',dataf_peakannoHK27DKO$annotation)
dataf_peakannoHK27DKO$annotation<-gsub('Exon.*','Exon',dataf_peakannoHK27DKO$annotation)
dataf_peakannoHK27DKO$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoHK27DKO$annotation)
dataf_peakannoHK27DKO$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoHK27DKO$annotation)
dataf_peakannoHK27DKO$peaksite<-paste(dataf_peakannoHK27DKO$geneId,dataf_peakannoHK27DKO$annotation,sep='_')
table(dataf_peakannoHK27DKO$annotation)

dataf_peakannoHK27ikb<-as.data.frame(peakAnnoHK27ikb)
dataf_peakannoHK27ikb$annotation<-gsub('Intron.*','Intron',dataf_peakannoHK27ikb$annotation)
dataf_peakannoHK27ikb$annotation<-gsub('Exon.*','Exon',dataf_peakannoHK27ikb$annotation)
dataf_peakannoHK27ikb$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoHK27ikb$annotation)
dataf_peakannoHK27ikb$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoHK27ikb$annotation)
dataf_peakannoHK27ikb$peaksite<-paste(dataf_peakannoHK27ikb$geneId,dataf_peakannoHK27ikb$annotation,sep='_')
table(dataf_peakannoHK27ikb$annotation)

dataf_peakannoHK27nfk<-as.data.frame(peakAnnoHK27nfk)
dataf_peakannoHK27nfk$annotation<-gsub('Intron.*','Intron',dataf_peakannoHK27nfk$annotation)
dataf_peakannoHK27nfk$annotation<-gsub('Exon.*','Exon',dataf_peakannoHK27nfk$annotation)
dataf_peakannoHK27nfk$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoHK27nfk$annotation)
dataf_peakannoHK27nfk$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoHK27nfk$annotation)
dataf_peakannoHK27nfk$peaksite<-paste(dataf_peakannoHK27nfk$geneId,dataf_peakannoHK27nfk$annotation,sep='_')
table(dataf_peakannoHK27nfk$annotation)

dataf_peakannoHK27WT<-as.data.frame(peakAnnoHK27WT)
dataf_peakannoHK27WT$annotation<-gsub('Intron.*','Intron',dataf_peakannoHK27WT$annotation)
dataf_peakannoHK27WT$annotation<-gsub('Exon.*','Exon',dataf_peakannoHK27WT$annotation)
dataf_peakannoHK27WT$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoHK27WT$annotation)
dataf_peakannoHK27WT$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoHK27WT$annotation)
dataf_peakannoHK27WT$peaksite<-paste(dataf_peakannoHK27WT$geneId,dataf_peakannoHK27WT$annotation,sep='_')
table(dataf_peakannoHK27WT$annotation)

dataf_peakannoHK36DKO<-as.data.frame(peakAnnoHK36DKO)
dataf_peakannoHK36DKO$annotation<-gsub('Intron.*','Intron',dataf_peakannoHK36DKO$annotation)
dataf_peakannoHK36DKO$annotation<-gsub('Exon.*','Exon',dataf_peakannoHK36DKO$annotation)
dataf_peakannoHK36DKO$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoHK36DKO$annotation)
dataf_peakannoHK36DKO$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoHK36DKO$annotation)
dataf_peakannoHK36DKO$peaksite<-paste(dataf_peakannoHK36DKO$geneId,dataf_peakannoHK36DKO$annotation,sep='_')
table(dataf_peakannoHK36DKO$annotation)

dataf_peakannoHK36ikb<-as.data.frame(peakAnnoHK36ikb)
dataf_peakannoHK36ikb$annotation<-gsub('Intron.*','Intron',dataf_peakannoHK36ikb$annotation)
dataf_peakannoHK36ikb$annotation<-gsub('Exon.*','Exon',dataf_peakannoHK36ikb$annotation)
dataf_peakannoHK36ikb$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoHK36ikb$annotation)
dataf_peakannoHK36ikb$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoHK36ikb$annotation)
dataf_peakannoHK36ikb$peaksite<-paste(dataf_peakannoHK36ikb$geneId,dataf_peakannoHK36ikb$annotation,sep='_')
table(dataf_peakannoHK36ikb$annotation)

dataf_peakannoHK36nfk<-as.data.frame(peakAnnoHK36nfk)
dataf_peakannoHK36nfk$annotation<-gsub('Intron.*','Intron',dataf_peakannoHK36nfk$annotation)
dataf_peakannoHK36nfk$annotation<-gsub('Exon.*','Exon',dataf_peakannoHK36nfk$annotation)
dataf_peakannoHK36nfk$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoHK36nfk$annotation)
dataf_peakannoHK36nfk$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoHK36nfk$annotation)
dataf_peakannoHK36nfk$peaksite<-paste(dataf_peakannoHK36nfk$geneId,dataf_peakannoHK36nfk$annotation,sep='_')
table(dataf_peakannoHK36nfk$annotation)

dataf_peakannoHK36WT<-as.data.frame(peakAnnoHK36WT)
dataf_peakannoHK36WT$annotation<-gsub('Intron.*','Intron',dataf_peakannoHK36WT$annotation)
dataf_peakannoHK36WT$annotation<-gsub('Exon.*','Exon',dataf_peakannoHK36WT$annotation)
dataf_peakannoHK36WT$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoHK36WT$annotation)
dataf_peakannoHK36WT$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoHK36WT$annotation)
dataf_peakannoHK36WT$peaksite<-paste(dataf_peakannoHK36WT$geneId,dataf_peakannoHK36WT$annotation,sep='_')
table(dataf_peakannoHK36WT$annotation)


# L1 stages 
dataf_peakannoIK27<-as.data.frame(peakAnnoIK27)
dataf_peakannoIK27$annotation<-gsub('Intron.*','Intron',dataf_peakannoIK27$annotation)
dataf_peakannoIK27$annotation<-gsub('Exon.*','Exon',dataf_peakannoIK27$annotation)
dataf_peakannoIK27$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoIK27$annotation)
dataf_peakannoIK27$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoIK27$annotation)
dataf_peakannoIK27$peaksite<-paste(dataf_peakannoIK27$geneId,dataf_peakannoIK27$annotation,sep='_')
table(dataf_peakannoIK27$annotation)

dataf_peakannoNF27<-as.data.frame(peakAnnoNF27)
dataf_peakannoNF27$annotation<-gsub('Intron.*','Intron',dataf_peakannoNF27$annotation)
dataf_peakannoNF27$annotation<-gsub('Exon.*','Exon',dataf_peakannoNF27$annotation)
dataf_peakannoNF27$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoNF27$annotation)
dataf_peakannoNF27$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoNF27$annotation)
dataf_peakannoNF27$peaksite<-paste(dataf_peakannoNF27$geneId,dataf_peakannoNF27$annotation,sep='_')
table(dataf_peakannoNF27$annotation)


dataf_peakannoWT27<-as.data.frame(peakAnnoWT27)
dataf_peakannoWT27$annotation<-gsub('Intron.*','Intron',dataf_peakannoWT27$annotation)
dataf_peakannoWT27$annotation<-gsub('Exon.*','Exon',dataf_peakannoWT27$annotation)
dataf_peakannoWT27$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoWT27$annotation)
dataf_peakannoWT27$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoWT27$annotation)
dataf_peakannoWT27$peaksite<-paste(dataf_peakannoWT27$geneId,dataf_peakannoWT27$annotation,sep='_')
table(dataf_peakannoWT27$annotation)


dataf_peakannoIK36<-as.data.frame(peakAnnoIK36)
dataf_peakannoIK36$annotation<-gsub('Intron.*','Intron',dataf_peakannoIK36$annotation)
dataf_peakannoIK36$annotation<-gsub('Exon.*','Exon',dataf_peakannoIK36$annotation)
dataf_peakannoIK36$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoIK36$annotation)
dataf_peakannoIK36$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoIK36$annotation)
dataf_peakannoIK36$peaksite<-paste(dataf_peakannoIK36$geneId,dataf_peakannoIK36$annotation,sep='_')
table(dataf_peakannoIK36$annotation)

dataf_peakannoNF36<-as.data.frame(peakAnnoNF36)
dataf_peakannoNF36$annotation<-gsub('Intron.*','Intron',dataf_peakannoNF36$annotation)
dataf_peakannoNF36$annotation<-gsub('Exon.*','Exon',dataf_peakannoNF36$annotation)
dataf_peakannoNF36$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoNF36$annotation)
dataf_peakannoNF36$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoNF36$annotation)
dataf_peakannoNF36$peaksite<-paste(dataf_peakannoNF36$geneId,dataf_peakannoNF36$annotation,sep='_')
table(dataf_peakannoNF36$annotation)

# No WT36

write.table(dataf_peakannoIK27,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/dataf_peakanno_IK_H27.tab",quote = FALSE,row.names = FALSE,sep="\t")
write.table(dataf_peakannoNF27,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/dataf_peakanno_NFKI_H27.tab",quote = FALSE,row.names = FALSE,sep="\t")
write.table(dataf_peakannoWT27,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/dataf_peakanno_WT_H27.tab.tab",quote = FALSE,row.names = FALSE,sep="\t")
write.table(dataf_peakannoIK36,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/dataf_peakanno_IK_H36.tab.tab",quote = FALSE,row.names = FALSE,sep="\t")
write.table(dataf_peakannoNF36,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/dataf_peakanno_NFKI_H36.tab.tab",quote = FALSE,row.names = FALSE,sep="\t")


# peaksites
# unique sites
HK27DKO_sites<-unique(dataf_peakannoHK27DKO$peaksite)
HK27Dikb_sites<-unique(dataf_peakannoHK27ikb$peaksite)
HK27Dnfk_sites<-unique(dataf_peakannoHK27nfk$peaksite)
HK27DWT_sites<-unique(dataf_peakannoHK27WT$peaksite)



distFMF1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFMF1$annotation)))))
colnames(distFMF1) <- as.character(unlist(distFMF1[1,]))
distFMF1<-distFMF1[-1,]
distFMF1$Sample<-"FMF1"


distFMF2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFMF2$annotation)))))
colnames(distFMF2) <- as.character(unlist(distFMF2[1,]))
distFMF2<-distFMF2[-1,]
distFMF2$Sample<-"FMF2"


distFMM1<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFMM1$annotation)))))
colnames(distFMM1) <- as.character(unlist(distFMM1[1,]))
distFMM1<-distFMM1[-1,]
distFMM1$Sample<-"FMM1"


distFMM2<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFMM2$annotation)))))
colnames(distFMM2) <- as.character(unlist(distFMM2[1,]))
distFMM2<-distFMM2[-1,]
distFMM2$Sample<-"FMM2"

distFMM3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFMM3$annotation)))))
colnames(distFMM3) <- as.character(unlist(distFMM3[1,]))
distFMM3<-distFMM3[-1,]
distFMM3$Sample<-"FMM3"


distHK27DKO<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoHK27DKO$annotation)))))
colnames(distHK27DKO) <- as.character(unlist(distHK27DKO[1,]))
distHK27DKO<-distHK27DKO[-1,]
distHK27DKO$Sample<-"HK27DKO"

distHK27ikb<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoHK27ikb$annotation)))))
colnames(distHK27ikb) <- as.character(unlist(distHK27ikb[1,]))
distHK27ikb<-distHK27ikb[-1,]
distHK27ikb$Sample<-"HK27ikb"


distHK27nfk<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoHK27nfk$annotation)))))
colnames(distHK27nfk) <- as.character(unlist(distHK27nfk[1,]))
distHK27nfk<-distHK27nfk[-1,]
distHK27nfk$Sample<-"HK27nfk"


distHK27wt<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoHK27WT$annotation)))))
colnames(distHK27wt) <- as.character(unlist(distHK27wt[1,]))
distHK27wt<-distHK27wt[-1,]
distHK27wt$Sample<-"HK27WT"


distHK36DKO<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoHK36DKO$annotation)))))
colnames(distHK36DKO) <- as.character(unlist(distHK36DKO[1,]))
distHK36DKO<-distHK36DKO[-1,]
distHK36DKO$Sample<-"HK36DKO"

distHK36ikb<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoHK36ikb$annotation)))))
colnames(distHK36ikb) <- as.character(unlist(distHK36ikb[1,]))
distHK36ikb<-distHK36ikb[-1,]
distHK36ikb$Sample<-"HK36ikb"


distHK36nfk<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoHK36nfk$annotation)))))
colnames(distHK36nfk) <- as.character(unlist(distHK36nfk[1,]))
distHK36nfk<-distHK36nfk[-1,]
distHK36nfk$Sample<-"HK36nfk"


distHK36wt<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoHK36WT$annotation)))))
colnames(distHK36wt) <- as.character(unlist(distHK36wt[1,]))
distHK36wt<-distHK36wt[-1,]
distHK36wt$Sample<-"HK36WT"


# L1
distIK27<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoIK27$annotation)))))
colnames(distIK27) <- as.character(unlist(distIK27[1,]))
distIK27<-distIK27[-1,]
distIK27$Sample<-"IK27_L1"

distNF27<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoNF27$annotation)))))
colnames(distNF27) <- as.character(unlist(distNF27[1,]))
distNF27<-distNF27[-1,]
distNF27$Sample<-"NF27_L1"

distWT27<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoWT27$annotation)))))
colnames(distWT27) <- as.character(unlist(distWT27[1,]))
distWT27<-distWT27[-1,]
distWT27$Sample<-"WT27_L1"


distIK36<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoIK36$annotation)))))
colnames(distIK36) <- as.character(unlist(distIK36[1,]))
distIK36<-distIK36[-1,]
distIK36$Sample<-"IK36_L1"

distNF36<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoNF36$annotation)))))
colnames(distNF36) <- as.character(unlist(distNF36[1,]))
distNF36<-distNF36[-1,]
distNF36$Sample<-"NF36_L1"

## Cols for FMM and FMF

colnames(distFMF1)
distFMF1$Exon<-"0"
#distFMF1$Three_UTR<-"0"
#distFMF1$Exon<-"0"

colnames(distFMF2)
#colnames(distFMF2)[1]<-"Three_UTR"
distFMF2$Exon<-"0"
#distFMF2<-distFMF2[c(2,3,4,5,6,1,7)]
distFMF2$Intron<-"0"
distFMF2<-distFMF2[c(1,2,6,3,4,5)]

colnames(distFMM1)
#distFMM1$Three_UTR<-"0"
#distFMM1$Exon<-"0"
distFMM1$Intron<-"0"
#distFMM1<-distFMM1[c(1,2,6,3,4,5,7)]
distFMM1<-distFMM1[c(1,2,5,3,4)]
distFMM1$Exon<-"0"

colnames(distFMM2)
#distFMM2$Three_UTR<-"0"
#distFMM2$Exon<-"0"
distFMM2$`Distal Intergenic`<-"0"
#distFMM2<-distFMM2[c(1,2,6,3,4,5,7)]
distFMM2<-distFMM2[c(5,1,2,3,4)]
distFMM2$Exon<-"0"

colnames(distFMM3)
#distFMM3$Three_UTR<-"0"
#distFMM3<-distFMM3[c(1,2,4,5,6,7,3)]
distFMM3$`Distal Intergenic`<-"0"
distFMM3<-distFMM3[c(6,1,3,4,5,2)]

distFMF<-rbind(distFMF1,distFMF2)
distFMF$group<-"FMF"
colnames(distFMF)


distFMM<-rbind(distFMM1,distFMM2,distFMM3)
distFMM$group<-"FMM"
colnames(distFMM)

distchip<-rbind(distFMF,distFMM)

distchip_melt<-melt(distchip,id.vars = c("group","Sample"))

distchip_melt$value<-as.numeric(distchip_melt$value)

class(distchip_melt$variable)
distchip_melt$variable <- factor(distchip_melt$variable, levels = c("Exon","Distal Intergenic","Intron","Downstream","Promoter"))

ggplot(distchip_melt,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat="identity",color="black") +
  facet_wrap(~group,scale="free",ncol=4)+
  scale_fill_brewer(palette="Pastel1")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")


# Cols for histone chips L4

colnames(distHK27DKO)
colnames(distHK27ikb)
colnames(distHK27nfk)
colnames(distHK27wt)

colnames(distHK36DKO)
colnames(distHK36ikb)
colnames(distHK36nfk)
colnames(distHK36wt)


# L4 H27
distchipHK27<-rbind(distHK27DKO,distHK27ikb,distHK27nfk,distHK27wt)

distchipHK27_melt<-melt(distchipHK27,id.vars = c("Sample"))

distchipHK27_melt$value<-as.numeric(distchipHK27_melt$value)
distchipHK27_melt$group<-"HK27"

#L4 H36
distchipHK36<-rbind(distHK36DKO,distHK36ikb,distHK36nfk,distHK36wt)

distchipHK36_melt<-melt(distchipHK36,id.vars = c("Sample"))

distchipHK36_melt$value<-as.numeric(distchipHK36_melt$value)
distchipHK36_melt$group<-"HK36"

# merge both H36 and H27 L4
distchipHist_melt<-rbind(distchipHK27_melt,distchipHK36_melt)


ggplot(distchipHist_melt,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat="identity") +
  facet_wrap(~group,scale="free",ncol=4)+
  scale_fill_brewer(palette="Pastel1")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")


#L1 H27

#histone chips L1
colnames(distIK27)
colnames(distNF27)
colnames(distWT27)

distIK27$`5' UTR`<-"0"
distIK27<-distIK27[c(1,8,2,3,4,5,6,7)]

distNF27$`5' UTR`<-"0"
distNF27<-distNF27[c(1,8,2,3,4,5,6,7)]


distWT27$`3' UTR`<-"0"
distWT27$Downstream<-"0"
distWT27$Intron<-"0"
distWT27$`5' UTR`<-"0"
distWT27<-distWT27[c(5,8,1,6,2,7,3,4)]


distchiL127<-rbind(distIK27,distNF27,distWT27)

distchiL127_melt<-melt(distchiL127,id.vars = c("Sample"))

distchiL127_melt$value<-as.numeric(distchiL127_melt$value)
distchiL127_melt$group<-"L1_H27"


colnames(distIK36)
colnames(distNF36)


distIK36$`5' UTR`<-"0"
distIK36$`Distal Intergenic`<-"0"
distIK36<-distIK36[c(1,7,8,2,3,4,5,6)]
colnames(distIK36)

distNF36$`5' UTR`<-"0"
distNF36<-distNF36[c(1,8,2,3,4,5,6,7)]
colnames(distNF36)


distchiL136<-rbind(distIK36,distNF36)

distchiL136_melt<-melt(distchiL136,id.vars = c("Sample"))

distchiL136_melt$value<-as.numeric(distchiL136_melt$value)
distchiL136_melt$group<-"L1_H36"



# merge both H36 and H27 L4
distchipL1_melt<-rbind(distchiL127_melt,distchiL136_melt)


ggplot(distchipL1_melt,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat="identity") +
  facet_wrap(~group,scale="free",ncol=4)+
  scale_fill_brewer(palette="Pastel1")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")


# TABLE WITH NUMBER OF PEAKS PER EXPERIMENT
peaksall <- data.frame("Sample" = c("WT","IKB_KO","NFKI_KO","DKO","WT","IKB_KO","NFKI_KO","DKO","WT","IKB_KO","NFKI_KO","WT","IKB_KO","NFKI_KO"), 
                       "Histone_mark" = c("HK36","HK36","HK36","HK36","HK27","HK27","HK27","HK27","HK36","HK36","HK36","HK27","HK27","HK27"),
                       "Stage"=c("L4","L4","L4","L4","L4","L4","L4","L4","L1","L1","L1","L1","L1","L1"),
                       "Number_peaks"=c(5684,5422,5480,4992,11319,10817,11505,11848,0,510,2372,79,1319,592))


peaksall$Sample <- factor(peaksall$Sample, levels = c("WT","IKB_KO","NFKI_KO","DKO"))
peaksall$Histone_mark <- factor(peaksall$Histone_mark , levels = c("HK27","HK36"))


ggplot(peaksall,aes(x=Sample,y=Number_peaks))+
  geom_point(size=5,aes(color=Stage),color="black")+
  geom_point(size=4,aes(color=Stage))+
  scale_color_brewer(palette="Set2")+
  facet_wrap(~Histone_mark,ncol=4,scales="free_x")+
  theme_bw()+
  theme(legend.position="bottom",
      axis.text.x = element_text(size=12,angle = 45, hjust = 1),
      axis.text.y = element_text(size=12, hjust = 1))


## VENN DIAGRAM  
library(VennDiagram)

geneLists<-list(comgenesFMF$Var1,comgenesFMM$Var1)


# And for histones
geneLists<-list(UPNFKI,HK27DWT_genes,HK27Dnfk_genes)

geneLists<-list(HK27DWT_genes,HK27Dikb_genes)

geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off() 


library(VennDiagram)
# FOR four GROUPS
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue","pink"), alpha=c(0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("UPKO","peakWT","peakKOnfki"), main="Genes ChIP")

venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue"), alpha=c(0.3,0.3), cex = 2, cat.fontface=4, category.names=c("genesWT","genesKOikb"), main="Genes ChIP HK27 L4")

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

# To get the list of gene present in each Venn compartment we can use the gplots package
require("gplots")


a <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(a,"intersections")
length(inters[[1]])
length(inters[[2]])
length(inters[[3]])
#length(inters[[4]])

#Genes upregulated in KO, with peak only in WT

uniqueWTikb<-dataf_peakannoHK27WT[which(dataf_peakannoHK27WT$geneId %in% inters[[1]]),]
uniqueWTikb$state<-"uniqueWT_ikb"

uniqueKOikb<-dataf_peakannoHK27ikb[which(dataf_peakannoHK27ikb$geneId %in% inters[[2]]),]
uniqueKOikb$state<-"uniqueKO_ikb"

commonikbWT<-dataf_peakannoHK27WT[which(dataf_peakannoHK27WT$geneId %in% inters[[3]]),]
commonikbWT$state<-"common_WT_ikb"

ikbpeaks_state<-rbind(uniqueWTikb,uniqueKOikb,commonikbWT)
ikbpeaks_state<-subset(ikbpeaks_state,select=c(geneId,state))
colnames(ikbpeaks_state)[1]<-"id"

#dataset$V5<-as.numeric(dataset$V5)
#dataset<-dataset[order(-(dataset$V5)),] 

write.table(inters[[1]],"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_ikbKO_unique.txt",quote = FALSE,row.names = FALSE)
write.table(inters[[2]],"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_ikbWT_unique.txt",quote = FALSE,row.names = FALSE)
write.table(inters[[1]],"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_nfkiKO_unique.txt",quote = FALSE,row.names = FALSE)
write.table(inters[[2]],"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_nfkiWT_unique.txt",quote = FALSE,row.names = FALSE)
write.table(inters[[1]],"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_DKO_unique.txt",quote = FALSE,row.names = FALSE)
write.table(inters[[2]],"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/peaks_DWT_unique.txt",quote = FALSE,row.names = FALSE)

### Peak distribution across genome

chromchipFMM<-as.data.frame(table(sort(c(FMM1_genes,FMM2_genes,FMM3_genes))))
colnames(chromchipFMM)[1]<-"geneId"
chromchipFMM$group<-"FMM"

chromchipFMF<-as.data.frame(table(sort(c(FMF1_genes,FMF2_genes))))
colnames(chromchipFMF)[1]<-"geneId"
chromchipFMF$group<-"FMF"


chromchipHK27DKO<-as.data.frame(table(sort(HK27DKO_genes)))
colnames(chromchipHK27DKO)[1]<-"geneId"
chromchipHK27DKO$group<-"HK27DKO"

chromchipHK27ikb<-as.data.frame(table(sort(HK27Dikb_genes)))
colnames(chromchipHK27ikb)[1]<-"geneId"
chromchipHK27ikb$group<-"HK27ikb"

chromchipHK27nfk<-as.data.frame(table(sort(HK27Dnfk_genes)))
colnames(chromchipHK27nfk)[1]<-"geneId"
chromchipHK27nfk$group<-"HK27nfk"

chromchipHK27WT<-as.data.frame(table(sort(HK27DWT_genes)))
colnames(chromchipHK27WT)[1]<-"geneId"
chromchipHK27WT$group<-"HK27wt"

chromchipHK36DKO<-as.data.frame(table(sort(HK36DKO_genes)))
colnames(chromchipHK36DKO)[1]<-"geneId"
chromchipHK36DKO$group<-"HK36DKO"

chromchipHK36ikb<-as.data.frame(table(sort(HK36Dikb_genes)))
colnames(chromchipHK36ikb)[1]<-"geneId"
chromchipHK36ikb$group<-"HK36ikb"

chromchipHK36nfk<-as.data.frame(table(sort(HK36Dnfk_genes)))
colnames(chromchipHK36nfk)[1]<-"geneId"
chromchipHK36nfk$group<-"HK36nfk"

chromchipHK36WT<-as.data.frame(table(sort(HK36DWT_genes)))
colnames(chromchipHK36WT)[1]<-"geneId"
chromchipHK36WT$group<-"HK36wt"


## L1
chromchipHK36ikb_L1<-as.data.frame(table(sort(IK36_genes)))
colnames(chromchipHK36ikb_L1)[1]<-"geneId"
chromchipHK36ikb_L1$group<-"HK36ikb_L1"

chromchipHK36nfki_L1<-as.data.frame(table(sort(NF36_genes)))
colnames(chromchipHK36nfki_L1)[1]<-"geneId"
chromchipHK36nfki_L1$group<-"HK36nfki_L1"

chromchipHK27ikb_L1<-as.data.frame(table(sort(IK27_genes)))
colnames(chromchipHK27ikb_L1)[1]<-"geneId"
chromchipHK27ikb_L1$group<-"HK27ikb_L1"

chromchipHK27nfki_L1<-as.data.frame(table(sort(NF27_genes)))
colnames(chromchipHK27nfki_L1)[1]<-"geneId"
chromchipHK27nfki_L1$group<-"HK27nfki_L1"

chromchipHK27wt_L1<-as.data.frame(table(sort(WT27_genes)))
colnames(chromchipHK27wt_L1)[1]<-"geneId"
chromchipHK27wt_L1$group<-"HK27wt_L1"


genesinf1<-subset(dataf_peakannoFMF1,select=c(geneChr,start,geneId,V5,annotation))
genesinf2<-subset(dataf_peakannoFMF2,select=c(geneChr,start,geneId,V5,annotation))
genesinfFMF<-rbind(genesinf1,genesinf2)
genesinfFMF<-unique(genesinfFMF)

genesinf3<-subset(dataf_peakannoFMM1,select=c(geneChr,start,geneId,V5,annotation))
genesinf4<-subset(dataf_peakannoFMM2,select=c(geneChr,start,geneId,V5,annotation))
genesinf5<-subset(dataf_peakannoFMM3,select=c(geneChr,start,geneId,V5,annotation))

genesinfFMM<-rbind(genesinf3,genesinf4,genesinf5)
genesinfFMM<-unique(genesinfFMM)

chromFMF<-merge(chromchipFMF,genesinfFMF,by="geneId",all.x=TRUE)
chromFMF$group

chromFMM<-merge(chromchipFMM,genesinfFMM,by="geneId",all.x=TRUE)
chromFMM$group

#L4
genesinf6<-subset(dataf_peakannoHK27DKO,select=c(geneChr,start,geneId,V5,annotation))
chromKK27DKO<-merge(chromchipHK27DKO,genesinf6,by="geneId",all.x="TRUE")
genesinf7<-subset(dataf_peakannoHK27ikb,select=c(geneChr,start,geneId,V5,annotation))
chromKK27ikb<-merge(chromchipHK27ikb,genesinf7,by="geneId",all.x="TRUE")
genesinf8<-subset(dataf_peakannoHK27nfk,select=c(geneChr,start,geneId,V5,annotation))
chromKK27nfk<-merge(chromchipHK27nfk,genesinf8,by="geneId",all.x="TRUE")
genesinf9<-subset(dataf_peakannoHK27WT,select=c(geneChr,start,geneId,V5,annotation))
chromKK27WT<-merge(chromchipHK27WT,genesinf9,by="geneId",all.x="TRUE")

genesinfA<-subset(dataf_peakannoHK36DKO,select=c(geneChr,start,geneId,V5,annotation))
chromKK36DKO<-merge(chromchipHK36DKO,genesinfA,by="geneId",all.x="TRUE")
genesinfB<-subset(dataf_peakannoHK36ikb,select=c(geneChr,start,geneId,V5,annotation))
chromKK36ikb<-merge(chromchipHK36ikb,genesinfB,by="geneId",all.x="TRUE")
genesinfC<-subset(dataf_peakannoHK36nfk,select=c(geneChr,start,geneId,V5,annotation))
chromKK36nfk<-merge(chromchipHK36nfk,genesinfC,by="geneId",all.x="TRUE")
genesinfD<-subset(dataf_peakannoHK36WT,select=c(geneChr,start,geneId,V5,annotation))
chromKK36WT<-merge(chromchipHK36WT,genesinfD,by="geneId",all.x="TRUE")


# L1
genesinfE<-subset(dataf_peakannoIK27,select=c(geneChr,start,geneId,V5,annotation))
chromKK27ikb_L1<-merge(chromchipHK27ikb_L1,genesinfE,by="geneId",all.x="TRUE")
genesinfF<-subset(dataf_peakannoNF27,select=c(geneChr,start,geneId,V5,annotation))
chromKK27nfk_L1<-merge(chromchipHK27nfki_L1,genesinfF,by="geneId",all.x="TRUE")
genesinfG<-subset(dataf_peakannoIK36,select=c(geneChr,start,geneId,V5,annotation))
chromKK36ikb_L1<-merge(chromchipHK36ikb_L1,genesinfG,by="geneId",all.x="TRUE")
genesinfH<-subset(dataf_peakannoNF36,select=c(geneChr,start,geneId,V5,annotation))
chromKK36nfk_L1<-merge(chromchipHK36nfki_L1,genesinfH,by="geneId",all.x="TRUE")

genesinfI<-subset(dataf_peakannoWT27,select=c(geneChr,start,geneId,V5,annotation))
chromKK27WT_L1<-merge(chromchipHK27wt_L1,genesinfI,by="geneId",all.x="TRUE")

chromchip<-rbind(chromFMF,chromFMM,
                 chromKK27DKO,chromKK27ikb,chromKK27nfk,chromKK27WT,
                 chromKK36DKO,chromKK36ikb,chromKK36nfk,chromKK36WT,
                 chromKK27ikb_L1,chromKK27nfk_L1,chromKK27WT_L1,
                 chromKK36ikb_L1,chromKK36nfk_L1)
chromchip$geneChr<-gsub('7','M',chromchip$geneChr)
chromchip$geneChr<-gsub('6','X',chromchip$geneChr)

options(scipen=10000)

library(ggforce)

chromchip$annotation <- factor(chromchip$annotation, levels = c("Distal Intergenic","Downstream","Intron","Promoter","Exon","3' UTR","5' UTR"))


ggplot(subset(chromchip,chromchip$group=="HK27ikb"),aes(x=start,y=V5,label=geneId))+
  geom_point(aes(color=annotation))+
  geom_linerange(aes(x=start, ymax=V5, ymin=0,color=annotation))+
  #geom_text(aes(start, cov, label = Gene), data = subset(uniquegenes_KH3_WTH3_noX,uniquegenes_KH3_WTH3_noX$source=="SC"),check_overlap = TRUE,vjust=-0.5)+
  #geom_label_repel(aes(fill = type),colour = "black", fontface = "italic",size=2.5)+
  scale_color_manual(values=c("grey","coral3","salmon","firebrick1","dodgerblue","cornflowerblue","darkolivegreen"))+
  #scale_color_manual(values=c("coral3","cornflowerblue"))+
  #scale_fill_manual(values=c("darksalmon","firebrick1","darkslateblue"))+
  #facet_zoom(x = geneChr=="1", ylim = c(0, 200))+
  facet_wrap(group~geneChr,scales="free",ncol=7)+
  theme_bw()+
  #ylim(c(0,50))+
  theme(axis.text.x = element_blank())+
  ylab("Score")


#Add differentially expressed genes ikb and nkfi KOs vs WT
chromchip_gene<-chromchip
chromchip_gene$geneId<-as.character(chromchip$geneId)
colnames(chromchip_gene)[1]<-"id"


## NEED rnatab from later analysis RNA-seq

#chromchip_gene<-merge(chromchip_gene,rnatab,by="id",all.x=TRUE)
# For L4
chromchip_gene<-merge(chromchip_gene,rnatab,by="id")

#For L1
chromchip_gene<-merge(chromchip_gene,rnatab_L1,by="id")


ikb<-unique((chromchip_gene[which (chromchip_gene$ikbsig=="sig_ikb"),])$id)
nfki<-unique((chromchip_gene[which (chromchip_gene$nfksig=="sig_nfki"),])$id)


#Select genes 
p1<-ggplot(subset(chromchip,(chromchip$group=="FMM")),aes(x=start,y=V5,label=geneId))+
  geom_point(aes(color=annotation,shape=group),size=3)+
  geom_point(aes(color=annotation,shape=group),size=3)+
  geom_linerange(aes(x=start, ymax=V5, ymin=0,color=annotation))+
  #geom_text(aes(start, V5, label = geneId), data = subset(chromchip,(chromchip$group=="FMM" | chromchip$group=="FMF") & chromchip$geneChr=="1" | chromchip$V5 >100),check_overlap = TRUE,vjust=-0.5)+
  #geom_label_repel(data = chromchip[(chromchip$geneId %in% ikb) & chromchip$group =="FMM",], fontface = "italic",size=2.5)+
  scale_color_manual(values=c("grey","coral3","green","darkolivegreen","dodgerblue"))+
  #scale_color_manual(values=c("coral3","cornflowerblue"))+
  #scale_fill_manual(values=c("darksalmon","firebrick1","darkslateblue"))+
  #facet_zoom( ylim = c(0, 75))+
  facet_wrap(~geneChr,scales="free",ncol=6)+
  theme_bw()+
  ylim(c(0,125))+
  theme(axis.text.x = element_blank())+
  ylab("Score")

p2<-ggplot(subset(chromchip,(chromchip$group=="FMF")),aes(x=start,y=V5,label=geneId))+
  geom_point(aes(color=annotation,shape=group),size=3)+
  geom_linerange(aes(x=start, ymax=V5, ymin=0,color=annotation))+
  #geom_text(aes(start, V5, label = geneId), data = subset(chromchip,(chromchip$group=="FMM" | chromchip$group=="FMF") & chromchip$geneChr=="1" | chromchip$V5 >100),check_overlap = TRUE,vjust=-0.5)+
  #geom_label_repel(data = chromchip[(chromchip$geneId %in% ikb) & chromchip$group =="FMF",], fontface = "italic",size=2.5)+
  scale_color_manual(values=c("grey","coral3","green","darkolivegreen","dodgerblue"))+
  #scale_color_manual(values=c("coral3","cornflowerblue"))+
  #scale_fill_manual(values=c("darksalmon","firebrick1","darkslateblue"))+
  #facet_zoom( ylim = c(0, 75))+
  facet_wrap(~geneChr,scales="free",ncol=6)+
  theme_bw()+
  ylim(c(0,80))+
  theme(axis.text.x = element_blank())+
  ylab("Score")

grid.arrange(p1,p2,ncol=1)


### Add info from Gaydos
gaydos$gene
chromchip_gene<-merge(gaydos,chromchip_gene,by.x="gene",by.y="ensembl_gene_id",all.y=TRUE)

# merge FMF and FMM gene names with gaydos
chromchip$geneId

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl")
listAttributes(ensembl)
chr_genes <- getBM(attributes=c('ensembl_gene_id',
                                'wormbase_gene','external_gene_name',
                                'chromosome_name','start_position','end_position'), mart = ensembl)

chromchip_gaydos<-merge(chromchip,chr_genes,by.x="geneId",by.y="external_gene_name",all=TRUE)
chromchip_gaydos<-merge(chromchip_gaydos,gaydos,by.x="wormbase_gene",by.y="gene",all=TRUE)

table(chromchip_gaydos$type_exp)

chromchip_gaydos$type_exp[is.na(chromchip_gaydos$type_exp)] <- "NA"

type_fmm<-ggplot(subset(chromchip_gaydos,(chromchip_gaydos$group=="FMM")),aes(x=start,y=V5,label=geneId))+
  geom_point(aes(color=type_exp,shape=group),size=3)+
  geom_point(aes(color=type_exp,shape=group),size=3)+
  geom_linerange(aes(x=start, ymax=V5, ymin=0,color=type_exp))+
  #geom_text(aes(start, V5, label = geneId), data = subset(chromchip,(chromchip$group=="FMM" | chromchip$group=="FMF") & chromchip$geneChr=="1" | chromchip$V5 >100),check_overlap = TRUE,vjust=-0.5)+
  #geom_label_repel(data = chromchip[(chromchip$geneId %in% ikb) & chromchip$group =="FMM",], fontface = "italic",size=2.5)+
  scale_color_manual(values=c("coral3","darkolivegreen","dodgerblue","grey"))+
  #scale_color_manual(values=c("coral3","cornflowerblue"))+
  #scale_fill_manual(values=c("darksalmon","firebrick1","darkslateblue"))+
  #facet_zoom( ylim = c(0, 75))+
  facet_wrap(~geneChr,scales="free",ncol=6)+
  theme_bw()+
  ylim(c(0,125))+
  theme(axis.text.x = element_blank())+
  ylab("Score")


type_fmf<-ggplot(subset(chromchip_gaydos,(chromchip_gaydos$group=="FMF")),aes(x=start,y=V5,label=geneId))+
  geom_point(aes(color=type_exp,shape=group),size=3)+
  geom_point(aes(color=type_exp,shape=group),size=3)+
  geom_linerange(aes(x=start, ymax=V5, ymin=0,color=type_exp))+
  #geom_text(aes(start, V5, label = geneId), data = subset(chromchip,(chromchip$group=="FMM" | chromchip$group=="FMF") & chromchip$geneChr=="1" | chromchip$V5 >100),check_overlap = TRUE,vjust=-0.5)+
  #geom_label_repel(data = chromchip[(chromchip$geneId %in% ikb) & chromchip$group =="FMM",], fontface = "italic",size=2.5)+
  scale_color_manual(values=c("coral3","darkolivegreen","dodgerblue","grey"))+
  #scale_color_manual(values=c("coral3","cornflowerblue"))+
  #scale_fill_manual(values=c("darksalmon","firebrick1","darkslateblue"))+
  #facet_zoom( ylim = c(0, 75))+
  facet_wrap(~geneChr,scales="free",ncol=6)+
  theme_bw()+
  ylim(c(0,125))+
  theme(axis.text.x = element_blank())+
  ylab("Score")

grid.arrange(type_fmm,type_fmf)


table(chromchip_gaydos$type_exp,chromchip_gaydos$group)

#subset without FMM, FMF and wt histones
subchromchip_gene_HK36<-subset(chromchip_gene,chromchip_gene$group=="HK36nfk" |
                            chromchip_gene$group=="HK36ikb" |
                          chromchip_gene$group=="HK36DKO")


subchromchip_gene_HK36ikb<-subset(chromchip_gene,chromchip_gene$group=="HK36ikb")
subchromchip_gene_HK36ikb$annotation <- factor(subchromchip_gene_HK36ikb$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))

subchromchip_gene_HK36nfki<-subset(chromchip_gene,chromchip_gene$group=="HK36nfk")
subchromchip_gene_HK36nfki$annotation <- factor(subchromchip_gene_HK36nfki$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))

subchromchip_gene_HK36DKO<-subset(chromchip_gene,chromchip_gene$group=="HK36DKO")
subchromchip_gene_HK36DKO$annotation <- factor(subchromchip_gene_HK36DKO$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))


subchromchip_gene_HK27<-subset(chromchip_gene,chromchip_gene$group=="HK27nfk" |
                                 chromchip_gene$group=="HK27ikb" |
                                 chromchip_gene$group=="HK27DKO")

subchromchip_gene_HK27ikb<-subset(chromchip_gene,chromchip_gene$group=="HK27ikb")
subchromchip_gene_HK27ikb$annotation <- factor(subchromchip_gene_HK27ikb$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))

subchromchip_gene_HK27nfki<-subset(chromchip_gene,chromchip_gene$group=="HK27nfk")
subchromchip_gene_HK27nfki$annotation <- factor(subchromchip_gene_HK27nfki$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))

subchromchip_gene_HK27DKO<-subset(chromchip_gene,chromchip_gene$group=="HK27DKO")
subchromchip_gene_HK27DKO$annotation <- factor(subchromchip_gene_HK27DKO$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))

# L1
subchromchip_gene_HK27_L1<-subset(chromchip_gene,chromchip_gene$group=="HK27nfk_L1" |
                                 chromchip_gene$group=="HK27ikb_L1")

subchromchip_gene_HK27ikb_L1<-subset(chromchip_gene,chromchip_gene$group=="HK27ikb_L1")
subchromchip_gene_HK27ikb_L1$annotation <- factor(subchromchip_gene_HK27ikb_L1$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))

subchromchip_gene_HK27nfki_L1<-subset(chromchip_gene,chromchip_gene$group=="HK27nfki_L1")
subchromchip_gene_HK27nfki_L1$annotation <- factor(subchromchip_gene_HK27nfki_L1$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))


subchromchip_gene_HK36_L1<-subset(chromchip_gene,chromchip_gene$group=="HK36nfk_L1" |
                                    chromchip_gene$group=="HK36ikb_L1")

subchromchip_gene_HK36ikb_L1<-subset(chromchip_gene,chromchip_gene$group=="HK36ikb_L1")
subchromchip_gene_HK36ikb_L1$annotation <- factor(subchromchip_gene_HK36ikb_L1$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))

subchromchip_gene_HK36nfki_L1<-subset(chromchip_gene,chromchip_gene$group=="HK36nfki_L1")
subchromchip_gene_HK36nfki_L1$annotation <- factor(subchromchip_gene_HK36nfki_L1$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter","3'UTR","5'UTR"))


# FOR L1 or F4, (L1 ADD FACET_WRAP PER ANNOTATION BECAUSE THERE ARE A FEW  PEAKS SIGNIFICANT WITH RNASEQ)

p3<-ggplot(subchromchip_gene_HK27ikb,aes(x=group,y=log2FoldChange.x))+
  geom_jitter(data=subset(subchromchip_gene_HK27ikb,subchromchip_gene_HK27ikb$ikbsig=="sig_ikb" & !is.na(subchromchip_gene_HK27ikb$annotation)),aes(color=annotation),alpha=0.5)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK27ikb,subchromchip_gene_HK27ikb$ikbsig=="sig_ikb" & !is.na(subchromchip_gene_HK27ikb$annotation)),aes(fill=annotation),alpha=0.5,colour="black")+
  #geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK27ikb,subchromchip_gene_HK27ikb$ikbsig!="sig_ikb" & !is.na(subchromchip_gene_HK27ikb$annotation)),aes(colour=annotation),alpha=0,lwd=1)+
  #scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0,color="grey")+
  facet_grid(type_exp~.,space="free")+
  scale_color_brewer(palette='Set1')+
  scale_fill_brewer(palette='Set1')+
  ylab("log2FC ikb KO vs WT")+
  theme_bw()+
  theme(legend.position = "none")+
  #ylim(-7.5,8)+
  coord_flip()
p3

p3a<-ggplot(subchromchip_gene_HK27nfki,aes(x=group,y=log2FoldChange.y))+
  geom_jitter(data=subset(subchromchip_gene_HK27nfki,subchromchip_gene_HK27nfki$nfksig=="sig_nfki" & !is.na(subchromchip_gene_HK27nfki$annotation)),aes(color=annotation),alpha=0.5)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK27nfki,subchromchip_gene_HK27nfki$nfksig=="sig_nfki" & !is.na(subchromchip_gene_HK27nfki$annotation)),aes(fill=annotation),alpha=0.5,colour="black")+
  #geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK27nfki,subchromchip_gene_HK27nfki$nfksig!="sig_nfki" & !is.na(subchromchip_gene_HK27nfki$annotation)),aes(colour=annotation),alpha=0,lwd=1)+
  #scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  facet_grid(type_exp~.,space="free")+
  scale_fill_brewer(palette='Set1')+
  scale_color_brewer(palette='Set1')+
  ylab("log2FC nfki KO vs WT")+
  theme_bw()+
  theme(legend.position = "none")+
  #ylim(-7.5,8)+
  coord_flip()
p3a

# DKO only for L4
p3b<-ggplot(subchromchip_gene_HK27DKO,aes(x=group,y=log2FoldChange))+
  geom_jitter(data=subset(subchromchip_gene_HK27DKO,subchromchip_gene_HK27DKO$dkosig=="sig_dko" & !is.na(subchromchip_gene_HK27DKO$annotation)),aes(color=annotation),alpha=0.5)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK27DKO,subchromchip_gene_HK27DKO$dkosig=="sig_dko" & !is.na(subchromchip_gene_HK27DKO$annotation)),aes(fill=annotation),alpha=0.5,colour="black")+
  #geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK27DKO,subchromchip_gene_HK27DKO$dkosig!="sig_dko" & !is.na(subchromchip_gene_HK27DKO$annotation)),aes(colour=annotation),alpha=0,lwd=1)+
  #scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  facet_grid(type_exp~.,space="free")+
  scale_fill_brewer(palette='Set1')+
  scale_color_brewer(palette='Set1')+
  ylab("log2FC DKO KO vs WT")+
  theme_bw()+
  theme(legend.position = "none")+
  #ylim(-7.5,8)+
  coord_flip()

p3b

grid.arrange(p3,p3a,p3b,ncol=3)

subchromchip_gene_HK36$annotation <- factor(subchromchip_gene_HK36$annotation, levels = c("Exon","Downstream","Intron","Distal Intergenic","Promoter"))

p4<-ggplot(subchromchip_gene_HK36ikb,aes(x=group,y=log2FoldChange.x))+
  geom_jitter(data=subset(subchromchip_gene_HK36ikb,subchromchip_gene_HK36ikb$ikbsig=="sig_ikb" & !is.na(subchromchip_gene_HK36ikb$annotation)),aes(color=annotation),alpha=0.5)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK36ikb,subchromchip_gene_HK36ikb$ikbsig=="sig_ikb" & !is.na(subchromchip_gene_HK36ikb$annotation)),aes(fill=annotation),alpha=0.5,colour="black")+
  #geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK36ikb,subchromchip_gene_HK36ikb$ikbsig!="sig_ikb" & !is.na(subchromchip_gene_HK36ikb$annotation)),aes(colour=annotation),alpha=0,lwd=1)+
  #scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  scale_fill_brewer(palette='Set1')+
  scale_color_brewer(palette='Set1')+
  #facet_wrap(~annotation,nrow=1)+
  #scale_color_manual(values=c("red","darkolivegreen","purple","orange"))+
  #scale_fill_manual(values=c("red","darkolivegreen","purple","orange"))+
  facet_grid(type_exp~.,space="free")+
  ylab("log2FC ikb KO vs WT")+
  theme_bw()+
  theme(legend.position = "none")+
  #ylim(-7.5,8)+
  coord_flip()

p4


#Remove distal intergenic annotations and downstream

subchromchip_gene_HK36nfki<-subset(subchromchip_gene_HK36nfki,subchromchip_gene_HK36nfki$annotation!="Distal Intergenic" & subchromchip_gene_HK36nfki$annotation!="Downstream")

p4a<-ggplot(subchromchip_gene_HK36nfki,aes(x=group,y=log2FoldChange.y))+
  geom_jitter(data=subset(subchromchip_gene_HK36nfki,subchromchip_gene_HK36nfki$nfksig=="sig_nfki" & !is.na(subchromchip_gene_HK36nfki$annotation)),aes(color=annotation),alpha=0.5)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK36nfki,subchromchip_gene_HK36nfki$nfksig=="sig_nfki" & !is.na(subchromchip_gene_HK36nfki$annotation)),aes(fill=annotation),alpha=0.5,colour="black")+
  #geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK36nfki,subchromchip_gene_HK36nfki$nfksig!="sig_nfki" & !is.na(subchromchip_gene_HK36nfki$annotation)),aes(colour=annotation),alpha=0,lwd=1)+
  #scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  #scale_color_brewer(palette = 'Set1')+
  #facet_wrap(~annotation,nrow=1)+
  facet_grid(type_exp~.,space="free")+
  scale_color_manual(values=c("red","darkolivegreen","orange"))+
  scale_fill_manual(values=c("red","darkolivegreen","orange"))+
  ylab("log2FC nfki KO vs WT")+
  theme_bw()+
  theme(legend.position = "none")+
  #ylim(-7.5,8)+
  coord_flip()

p4a

#Remove downstream
subchromchip_gene_HK36DKO<-subset(subchromchip_gene_HK36DKO,subchromchip_gene_HK36DKO$annotation!="Downstream")


p4b<-ggplot(subchromchip_gene_HK36DKO,aes(x=group,y=log2FoldChange))+
  geom_jitter(data=subset(subchromchip_gene_HK36DKO,subchromchip_gene_HK36DKO$dkosig=="sig_dko" & !is.na(subchromchip_gene_HK36DKO$annotation)),aes(color=annotation),alpha=0.5)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK36DKO,subchromchip_gene_HK36DKO$dkosig=="sig_dko" & !is.na(subchromchip_gene_HK36DKO$annotation)),aes(fill=annotation),alpha=0.5,colour="black")+
  #geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK36DKO,subchromchip_gene_HK36DKO$dkosig!="sig_dko" & !is.na(subchromchip_gene_HK36DKO$annotation)),aes(colour=annotation),alpha=0,lwd=1)+
  #scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  #scale_color_brewer(palette = 'Set1')+
  facet_grid(type_exp~.,space="free")+
  scale_color_manual(values=c("red","darkolivegreen","purple","orange"))+
  scale_fill_manual(values=c("red","darkolivegreen","purple","orange"))+
  ylab("log2FC DKO vs WT")+
  theme_bw()+
  theme(legend.position = "none")+
  #ylim(-7.5,8)+
  coord_flip()

p4b

grid.arrange(p4,p4a,p4b,ncol=3)

grid.arrange(p3,p3a,p3b,p4,p4a,p4b,ncol=3)


p5<-ggplot(subchromchip_gene_HK27_L1,aes(x=group,y=log2FoldChange.y))+
  #geom_point(aes(color=sig))+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK27_L1,subchromchip_gene_HK27_L1$nfksig=="sig_nfki"),aes(color=annotation),alpha=0.2,lwd=1)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK27_L1,subchromchip_gene_HK27_L1$nfksig!="sig_nfki"),aes(color=annotation),linetype="dotted",alpha=0.2,lwd=1)+
  #scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  scale_color_brewer(type='div')+
  ylab("log2FC nfki KO vs WT")+
  theme_bw()+
  #theme(legend.position = "none")+
  ylim(-7.5,8)+
  coord_flip()
p5

p6<-ggplot(subchromchip_gene_HK36,aes(x=group,y=log2FoldChange.y))+
  #geom_point(aes(color=sig))+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK36,subchromchip_gene_HK36$nfksig=="sig_nfki"),aes(color=annotation),alpha=0.2,lwd=1)+
  geom_violin(draw_quantiles = c(0.5),data=subset(subchromchip_gene_HK36,subchromchip_gene_HK36$nfksig!="sig_nfki"),aes(color=annotation),linetype="dotted",alpha=0.2,lwd=1)+
  #scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  scale_color_brewer(type='div')+
  ylab("log2FC nfki KO vs WT")+
  theme_bw()+
  #theme(legend.position = "none")+
  ylim(-7.5,8)+
  coord_flip()
p6

grid.arrange(p5,p6,ncol=2)

grid.arrange(p3,p4,p5,p6)


#subset FMM and FMF
apl1<-ggplot(data=subset(chromchip_gene,chromchip_gene$group=="FMM"),aes(x=group,y=log2FoldChange.x))+
  geom_point(aes(color=sig),size=4)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMM" & chromchip_gene$sig=="sig_all")),color="red",alpha=0.2,lwd=1)+
  geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMM" & chromchip_gene$sig!="sig_all")),color="grey",alpha=0.2,lwd=1)+
  scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  ylab("log2FC ikb KO vs WT")+
  coord_flip()+
  theme_bw()

apl2<-ggplot(data=subset(chromchip_gene,chromchip_gene$group=="FMF"),aes(x=group,y=log2FoldChange.y))+
  geom_point(aes(color=sig),size=4)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMF" & chromchip_gene$sig=="sig_all")),color="red",alpha=0.2,lwd=1)+
  geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMF" & chromchip_gene$sig!="sig_all")),color="grey",alpha=0.2,lwd=1)+
  scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  ylab("log2FC nfki KO vs WT")+
  coord_flip()+
  theme_bw()

grid.arrange(apl1,apl2)



p4<-ggplot(chromchip_gene,aes(x=group,y=log2FoldChange.y))+
  geom_point(aes(color=sig))+
  geom_violin(data=(subset(chromchip_gene,chromchip_gene$nfksig=="sig_nfki")),color="red",alpha=0.2)+
  geom_violin(data=(subset(chromchip_gene,chromchip_gene$nfksig!="sig_nfki")),color="grey",alpha=0.2)+
  scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  geom_hline(yintercept = 0)+
  ylab("log2FC nfki KO vs WT")+
  theme_bw()

grid.arrange(p3,p4,ncol=1)

p5<-ggplot(chromchip_gene,aes(x=group,y=log2FoldChange.y))+
  geom_point(aes(color=nfksig))+
  geom_violin(data=(subset(chromchip_gene,chromchip_gene$nfksig!="sig_nfki")),alpha=0.2)+
  scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  geom_hline(yintercept = 0)+
  ylab("log2FC nfki KO vs WT")+
  theme_bw()

p6<-ggplot(chromchip_gene,aes(x=group,y=log2FoldChange.x))+
  geom_point(aes(color=ikbsig))+
  geom_violin(alpha=0.2)+
  scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  geom_hline(yintercept = 0)+
  ylab("log2FC ikb KO vs WT")+
  theme_bw()


grid.arrange(p3,p6,p4,p5,ncol=1)


#subset FMM and FMF (for L4  use sig=="sig_all", for L1 use sig_ikb or sig_nkfi)
apl1<-ggplot(data=subset(chromchip_gene,chromchip_gene$group=="FMM"),aes(x=group,y=log2FoldChange.x))+
  geom_point(aes(color=sig),size=4)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  #geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMM" & chromchip_gene$sig=="sig_all")),color="red",alpha=0.2,lwd=1)+
  geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMM" & chromchip_gene$sig!="sig_all")),color="grey",alpha=0.2,lwd=1)+
  geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMM" & chromchip_gene$ikbsig=="sig_ikb")),color="darkolivegreen",alpha=0.2,lwd=1)+
  scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  ylab("log2FC ikb KO vs WT")+
  #ylim(-6,4)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")

apl2<-ggplot(data=subset(chromchip_gene,chromchip_gene$group=="FMF"),aes(x=group,y=log2FoldChange.y))+
  geom_point(aes(color=sig),size=4)+
  #geom_label_repel(aes(label=id),data = chromchip_gene[(chromchip_gene$id %in% intgenes) & chromchip_gene$group=="FMM",], fontface = "italic",size=2.5)+
  geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMF" & chromchip_gene$sig=="sig_all")),color="red",alpha=0.2,lwd=1)+
  geom_violin(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMF" & chromchip_gene$sig!="sig_all")),color="grey",alpha=0.2,lwd=1)+
  geom_point(draw_quantiles = c(0.5),data=subset(chromchip_gene,(chromchip_gene$group=="FMF" & chromchip_gene$nfksig=="sig_nfki")),color="darkolivegreen",alpha=0.2,lwd=1)+
  scale_color_manual(values=c("grey","red"))+
  #geom_text(aes(label = SYMBOL),check_overlap = TRUE,vjust=-0.5)+
  #geom_hline(yintercept = 0)+
  ylab("log2FC nfki KO vs WT")+
  #ylim(-6,4)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")

grid.arrange(apl1,apl2,ncol=1)


## Hox genes from enriched functions (end script)
uphox<-c("ceh-13","pal-1","lin-11","lim-6","lim-7","ceh-37","ceh-38","ceh-39","ceh-49","ceh-44","ceh-83","ceh-91")
downhox<-c("mab-5","egl-5","nob-1","php3","vab-7","ceh-9","vab-15","ceh-43","ceh-6","zag-1","zfh-2","hmbx-1","ceh-20")



# peaks from wormenrichR rnaseq
## Importing genes from Wormbase enrichment up and down
upfunctions<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_UP/some_funct_develop.txt",sep="\t",header = FALSE)
colnames(upfunctions)<-c("id","func")

upfunctchip<-merge(chromchip_gene,upfunctions,by="id",all=TRUE)




### PATHWAYS ENRICHMENT ## Check C. elegans Genome!

#Up pathways from wormEnrichR
anam<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_UP/padj001/Anatomic_Associations_WormBase_2018_table.txt")
GOBP<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_UP/padj001/GO_Biological_Process_2018_table.txt")
wiki<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_UP/padj001/WikiPathways_2018_table.txt")

anam$group<-"anam_assoc"
GOBP$group<-"GO_BP"
wiki$group<-"WikiPath"

up_tab<-rbind(anam,GOBP,wiki)
up_tab$state<-"UP"

#Down pathways
#anam_down<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_DOWN/Anatomic_Associations_WormBase_2018_table.txt")
#GOBP_down<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_DOWN/GO_Biological_Process_2018_table.txt")

#anam_down$group<-"anam_assoc"
#GOBP_down$group<-"GO_BP"

#down_tab<-rbind(anam_down,GOBP_down)
#down_tab$state<-"DOWN"


#allGO<-rbind(up_tab,down_tab)
allGO<-up_tab


bpsub<-subset(allGO,allGO$Adjusted.P.value<=0.01 & allGO$Combined.Score>=50)
bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
#bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
#bpsub$GO<-gsub(':','',bpsub$GO)
#bpsub$GO<-gsub(')','',bpsub$GO)
#bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))

row.names(bpsub)<- factor(row.names(bpsub), levels=row.names(bpsub[rev(order(bpsub$Z.score)), ]))

bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$Z.score)), ]$Term)


ggplot(bpsub,aes(y=-(Z.score),x=Term,fill=state),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=state),color="black")+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Blues")+
  #facet_wrap(~state,scales="free",ncol=1)+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=6, hjust = 1),
        axis.text.y = element_text(size=9, hjust = 1))





#pathway1 <- enrichPathway(as.data.frame(peakAnnoFMM)$geneId,organism="celegans")
#head(pathway1, 2)
#dotplot(pathway1)


#gene <- seq2gene(peakFMF, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
#pathway2 <- enrichPathway(gene,organism = "celegans")
#head(pathway2, 2)

#dotplot(pathway2)




######### RNA-Seq #########

library(data.table)
library(ggplot2)
library(nlme)
library(readxl)
#library(pasilla)

## Get ids from ensembl all genes
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl")
chr_genes <- getBM(attributes=c('ensembl_gene_id',
                                'ensembl_transcript_id','external_gene_name',
                                'chromosome_name','start_position','end_position'), mart = ensembl)

library(htmlwidgets)
library(DESeq2)
library(cummeRbund)
library(biomaRt)
library(geneplotter)
library(MDplot)
library(lattice)
library(genefilter)
library(plotly)
library(limma)
library(ggrepel)


#Set wd
setwd("/Volumes/cancer/Celegans/RNASeq/htseq/")

output.dir="/Volumes/cancer/Celegans/RNASeq/htseq/"

# Import metadata files
#metadata<-read.delim("metadata.txt",sep="\t",header = FALSE)

#Including L1 rnaseq
metadata<-read.delim("metadata_2.txt",sep="\t",header = FALSE)

colnames(metadata)<-c("sampleID","countFile","condition","exp","stage")

metadata$group<-paste(metadata$condition,metadata$exp,sep="_")



#without external GEO
metadata<-subset(metadata,metadata$group!="WT_GEO")

#For KO ikb L4
metadata_ikb<-subset(metadata,metadata$stage=="L4" & (metadata$condition=="IKB" | metadata$group=="WT_KO"))

#For KO ikb L1
metadata_ikb_L1<-subset(metadata,metadata$stage=="L1" & (metadata$condition=="IKB" | metadata$group=="WT_KO" ))


#For KO nfkib L4
metadata_nfkb_L4<-subset(metadata,metadata$stage=="L4" & (metadata$condition=="nfkb"  | metadata$group=="WT_KO"))

#For KO nfkib L1
metadata_nfkb_L1<-subset(metadata,metadata$stage=="L1" & (metadata$condition=="nfkb" | metadata$group=="WT_KO" ))


#For DKO L4
metadata_dko<-subset(metadata,metadata$stage=="L4" & (metadata$condition=="DKO"  | metadata$group=="WT_KO"))


#For iRNA
metadata_i<-subset(metadata,metadata$exp=="iRNA" | metadata$group=="WT_iRNA")


#For WT L1 vs L4
metadata_WT<-subset(metadata,metadata$exp=="KO" & metadata$condition=="WT")

#For IKB KO L1 vs L4
metadata_KOL<-subset(metadata,metadata$exp=="KO" & metadata$condition=="IKB")

# excluding geo and iRNA
metadata_noi<-subset(metadata,metadata$exp!="GEO" & metadata$exp!="iRNA")


DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata_dko,
                                          directory = output.dir,
                                          design = ~ condition)


rowData(DESeq2Table)
mcols(rowData(DESeq2Table))

# Add additional annotation information using biomaRt
ensembl76 <- useMart("ensembl", dataset="celegans_gene_ensembl")


bm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","gene_biotype","chromosome_name","wormbase_gene"),
            values= rownames(DESeq2Table),
            #            filter="external_gene_name",
            mart=ensembl76) 

head(bm)

#Add description data to gene counts
DESeq2Features <- data.frame(id = rownames(DESeq2Table))
DESeq2Features$id <- as.character(DESeq2Features$id)


### join them together
#if htseq ensemblid
#rowData <- merge(DESeq2Features, bm, by.x="id",by.y = "ensembl_gene_id")

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
#colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$exp, levels = c("Control","shbcat"))
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","DKO"))
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","nfkb"))
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("WT","IKB"))
colData(DESeq2Table)$stage <- factor(colData(DESeq2Table)$stage, levels = c("L4","L1"))

#### estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

# plot densities of counts for the different samples to assess their distributions

multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))


### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Table, blind=TRUE) ## create a distance matrix between the samples

## Create PCA
DESeq2::plotPCA(rld, intgroup=c("exp","condition"))
dev.off()



## Create PCA 
pcaData <- plotPCA(rld, intgroup=c("condition","group","stage"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition,label=name)) +
  geom_point(size=9,color="black")+
  geom_point(size=8,aes(color=condition))+
  #geom_point(size=8,aes(shape=condition)) +
  #stat_ellipse(aes(color=exp),level=0.8)+
  #scale_color_manual(values=c("dodgerblue","coral"))+
  scale_color_manual(values=c("dodgerblue","coral"))+
  geom_label_repel(aes(label=name), size=3,fontface = "italic",alpha=0.8)+
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


ggplot(pcaData, aes(y=PC1, x=condition,label=name)) +
  geom_point(aes(color=group),size=8) +
  geom_boxplot(alpha=0)+
  #stat_ellipse(aes(color=condition),level=0.8)+
  #scale_color_manual(values=c("dodgerblue2","red"))+
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
levels(DESeq2Table$condition)
levels(DESeq2Table$stage)

design(DESeq2Table)

dds<-DESeq(DESeq2Table)
resultsNames(dds)


#results DEGs between condition
rescondition<-results(dds,alpha=0.05)

#Example differences genes between groups
plotCounts(dds, gene=which(rownames(rescondition)=="ikb-1"), intgroup="condition")

head(rescondition)

summary(rescondition)

#genes with < 0.05
sum(rescondition$padj < 0.05,na.rm=TRUE)

plotMA(rescondition)

#Heatmap most differentially expressed genes (only ranked by DEGs)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 300 )
library(gplots)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", dendrogram="column",trace = "none", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="grey", shbcat="black" )[
             colData(rld)$condition ],cexCol=0.3)

# Selecting specific DEGs
res_df_rna<-as.data.frame(rescondition)
class(rowData)
rowData<-as.data.frame(rowData)
res_df_rna<-merge(rowData,res_df_rna,by.x="id",by.y="row.names",all.y=TRUE)

res_df_sig_rna<-subset(res_df_rna,res_df_rna$pvalue<=0.01)
res_df_sig_rna <- res_df_sig_rna[order(res_df_sig_rna$pvalue),]


## DEGs ikb WT vs KO L1
IKB_L1<-res_df_rna

# DEGs nfki WT vs KO L1
NFKI_L1<-res_df_rna

# DEGs ikb WT vs KO L4
IKB_L4<-res_df_rna

# DEGs nkfi WT vs KO L4
NFKI_L4<-res_df_rna

# DEGs WT vs DKO L4
DKO_L4<-res_df_rna

## DEGs WT L1 vs WT L4
WTL1_L4<-res_df_rna

# DEGs ikb KO L1 vs IKB KO L4
IKB_L1_L4<-res_df_rna

write.table(IKB_L1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/L1_IKB_RNASeq.tab",quote = FALSE,row.names = FALSE,sep="\t")
write.table(NFKI_L1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/L1_NFKI_RNASeq.tab",quote = FALSE,row.names = FALSE,sep="\t")


FC_sign<-merge(IKB_L1,NFKI_L1,by="id")
# parameters for L4 analysis
#FC_sign<-subset(FC_sign,FC_sign$padj.x<=0.05 & FC_sign$padj.y<=0.05 & abs(FC_sign$log2FoldChange.x)>=2 & abs(FC_sign$log2FoldChange.y)>=2)

# parameters for L1 analysis
FC_sign<-subset(FC_sign,(FC_sign$pvalue.x<=0.05 | FC_sign$pvalue.y<=0.05))

ggplot(FC_sign,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  geom_point(color="black",size=3,alpha=0.5)+
  geom_density_2d(color="green")+
  geom_smooth(method="lm",color="red",se = FALSE)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = id),color="black", data=FC_sign[(FC_sign$log2FoldChange.x>3 & FC_sign$log2FoldChange.y>3),], fontface = "italic",size=5)+
  #xlim(-5,5)+
  #ylim(-5,5)+
  xlab("log2FC WT L1 vs IKB KO L1")+
  ylab("log2FC WT L1 vs NFKI KO L1")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1)) 


cor.test(FC_sign$log2FoldChange.x,FC_sign$log2FoldChange.y)

# UP and DOWN genes
UPgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange>0))$id)
DOWNgenes<-unique((subset(res_df_sig_rna,res_df_sig_rna$log2FoldChange<0))$id)
length(UPgenes)
length(DOWNgenes)


#Manual heatmap from log trasformed data, using all samples
rld_df<-as.data.frame(assay(rld))


#only degs
deggenes<-row.names(res_df_sig_rna)

# for condition
#all genes
deggenes<-row.names(res_df_rna)


deg_rld<-which(rownames(rld_df) %in% deggenes)
deg_rld_df<-rld_df[deg_rld,]

deg_rld<-merge(deg_rld_df,res_df_rna,by.x="row.names",by.y="row.names")

deg_rld<- deg_rld[order(deg_rld$log2FoldChange),]
deg_rld$Row.names

#CAREFUL!
#heatmap.2( as.matrix(deg_rld[,2:5]), scale="row", dendrogram="column",trace="none", Rowv=FALSE,
#           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
#           cexCol=0.3)



ikbtable<-IKB_L4
nfkitable<-NFKI_L4
dkotable<-DKO_L4

ikbtable_L1<-IKB_L1
nfkitable_L1<-NFKI_L1


write.table(ikbtable,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/IKB_L4_RNAseq_FC_genes.tab",quote = FALSE,sep="\t",row.names = FALSE)
write.table(nfkitable,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/NFKI_L4_RNAseq_FC_genes.tab",quote = FALSE,sep="\t",row.names = FALSE)
write.table(dkotable,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/DKO_L4_RNAseq_FC_genes.tab",quote = FALSE,sep="\t",row.names = FALSE)

write.table(ikbtable_L1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/IKB_L1_RNAseq_FC_genes.tab",quote = FALSE,sep="\t",row.names = FALSE)
write.table(nfkitable_L1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/NFKI_L1_RNAseq_FC_genes.tab",quote = FALSE,sep="\t",row.names = FALSE)

## Genes up and down in ikb, nfki and DKO KO (padjusted 0.01 para L4, pval 0.05 para L1)


UPIKB_L1<-unique((subset(ikbtable_L1,ikbtable_L1$pvalue<0.05 & ikbtable_L1$log2FoldChange>0))$id)
DOWNIKB_L1<-unique((subset(ikbtable_L1,ikbtable_L1$pvalue<0.05 & ikbtable_L1$log2FoldChange<0))$id)

UPNFKI_L1<-unique((subset(nfkitable_L1,nfkitable_L1$pvalue<0.05 & nfkitable_L1$log2FoldChange>0))$id)
DOWNNFKI_L1<-unique((subset(nfkitable_L1,nfkitable_L1$pvalue<0.05 & nfkitable_L1$log2FoldChange<0))$id)

UPDKO_L4<-unique((subset(dkotable,dkotable$padj<0.01 & dkotable$log2FoldChange>0))$id)
DOWNDKO_L4<-unique((subset(dkotable,dkotable$padj<0.01 & dkotable$log2FoldChange<0))$id)

UPIKB_L4<-unique((subset(ikbtable,ikbtable$padj<0.01 & ikbtable$log2FoldChange>0))$id)
DOWNIKB_L4<-unique((subset(ikbtable,ikbtable$padj<0.01 & ikbtable$log2FoldChange<0))$id)

UPNFKI_L4<-unique((subset(nfkitable,nfkitable$padj<0.01 & nfkitable$log2FoldChange>0))$id)
DOWNNFKI_L4<-unique((subset(nfkitable,nfkitable$padj<0.01 & nfkitable$log2FoldChange<0))$id)

write.table(UPIKB_L1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Upgenes_L1_ikb.txt",quote = FALSE,row.names = FALSE)
write.table(DOWNIKB_L1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Downgenes_L1_ikb.txt",quote = FALSE,row.names = FALSE)

write.table(UPNFKI_L1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Upgenes_L1_nfkb.txt",quote = FALSE,row.names = FALSE)
write.table(DOWNNFKI_L1,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Downgenes_L1_nfkb.txt",quote = FALSE,row.names = FALSE)


write.table(UPIKB,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Upgenes_ikb.txt",quote = FALSE,row.names = FALSE)
write.table(DOWNIKB,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Downgenes_ikb.txt",quote = FALSE,row.names = FALSE)

write.table(UPNFKI,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Upgenes_nfkb.txt",quote = FALSE,row.names = FALSE)
write.table(DOWNNFKI,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Downgenes_nfkb.txt",quote = FALSE,row.names = FALSE)

write.table(UPDKO,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Upgenes_DKO.txt",quote = FALSE,row.names = FALSE)
write.table(DOWNDKO,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/Downgenes_DKO.txt",quote = FALSE,row.names = FALSE)



## VENN DIAGRAM RNASEQ 
library(VennDiagram)
# Control bcat genes vs lithium bcat
geneLists<-list(UPIKB,UPDKO,UPNFKI)

geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off() 

# FOR four GROUPS
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "darkblue","coral3"), alpha=c(0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("UPIKB","UPDKO","UPNFKI"), main="Genes RNASeq")

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

# To get the list of gene present in each Venn compartment we can use the gplots package
require("gplots")

a <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(a,"intersections")
length(inters[[1]])
length(inters[[2]])
length(inters[[3]])


## Plot RNAseq two conditions correlation
ikbtable$sample<-"IKB"
nfkitable$sample<-"nfki"
dkotable$sample<-"DKO"

ikbtable_L1$sample<-"IKB_L1"
nfkitable_L1$sample<-"nfki_L1"

#L4
rnatab<-merge(ikbtable,nfkitable,by="id",all=TRUE)
rnatab<-merge(rnatab,dkotable,by="id")

rnatab$ikbsig<-ifelse(rnatab$padj.x < 0.01, 
                                   c("sig_ikb"), c("nosig_ikb"))  
rnatab$nfksig<-ifelse(rnatab$padj.y < 0.01, 
                      c("sig_nfki"), c("nosig_nfki"))  

rnatab$dkosig<-ifelse(rnatab$padj < 0.01, 
                      c("sig_dko"), c("nosig_dko")) 

rnatab$sig<-ifelse(rnatab$padj.y < 0.01 & rnatab$padj.x<0.01 & rnatab$padj <0.01, 
                      c("sig_all"), c("nosig_all")) 


rnatab$sigtype<-paste(rnatab$ikbsig,rnatab$nfksig,rnatab$dkosig)



ggplot(subset(rnatab,rnatab$pvalue.x < 0.05 | rnatab$pvalue.y < 0.05),aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  #annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "blue",alpha=0.1)+
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "red",alpha=0.1)+
  geom_point(color="black",size=3)+
  geom_point(aes(color=sig),size=0.5)+
  scale_color_manual(values=c("white","coral"))+
  geom_density_2d(color="black",bins=15)+
  #geom_density_2d(data=(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$sigtype=="sig_all sig_HSC")),color="purple")+
  #scale_color_manual(values=c("black","grey"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = id),color="coral", data = subset(rnatab,rnatab$id=="C33A11.1" | rnatab$id=="ikb-1"), fontface = "italic",size=2.5)+
  #geom_label_repel(aes(label = mgi_symbol.x), color="red",data = subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$significativeall=="sig_all" & B1B2B3_all_NP_comp$state.x=="KOdown" & B1B2B3_all_NP_comp$significativeHSC=="sig_HSC" & B1B2B3_all_NP_comp$state=="KOdown"), fontface = "italic",size=2.5)+
  #xlim(-10,10)+
  #ylim(-10,10)+
  xlab("log2FC ikb KO vs. WT")+
  ylab("log2FC nfkb vs. WT")+
  theme_bw()+
  #facet_wrap(~chromosome_name.y,ncol=2)+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1))+
  xlab("log2FC ikb-1 KO vs WT")+
  #xlab("log2FC nfki-1 (C33A11.1) vs. WT")+
  #ylab("log2FC DKO vs. WT")+
  ylab("log2FC nfki-1 vs. WT")


cor.test(rnatab$log2FoldChange.x,rnatab$log2FoldChange.y)


#L1
rnatab_L1<-merge(ikbtable_L1,nfkitable_L1,by="id",all=TRUE)

rnatab_L1$ikbsig<-ifelse(rnatab_L1$pvalue.x < 0.05, 
                      c("sig_ikb"), c("nosig_ikb"))  
rnatab_L1$nfksig<-ifelse(rnatab_L1$pvalue.y< 0.05, 
                      c("sig_nfki"), c("nosig_nfki"))  

rnatab_L1$sig<-ifelse(rnatab_L1$pvalue.x < 0.05 & rnatab_L1$pvalue.y<0.05, 
                   c("sig_all"), c("nosig_all")) 


rnatab_L1$sigtype<-paste(rnatab_L1$ikbsig,rnatab_L1$nfksig,rnatab$dkosig)



ggplot(subset(rnatab_L1,rnatab_L1$pvalue.x < 0.05 | rnatab_L1$pvalue.y < 0.05),aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  #annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "blue",alpha=0.1)+
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "red",alpha=0.1)+
  geom_point(color="black",size=3)+
  geom_point(aes(color=sig),size=0.5)+
  scale_color_manual(values=c("white","coral"))+
  geom_density_2d(color="black",bins=15)+
  #geom_density_2d(data=(subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$sigtype=="sig_all sig_HSC")),color="purple")+
  #scale_color_manual(values=c("black","grey"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = id),color="coral", data = subset(rnatab,rnatab$id=="C33A11.1" | rnatab$id=="ikb-1"), fontface = "italic",size=2.5)+
  #geom_label_repel(aes(label = mgi_symbol.x), color="red",data = subset(B1B2B3_all_NP_comp,B1B2B3_all_NP_comp$significativeall=="sig_all" & B1B2B3_all_NP_comp$state.x=="KOdown" & B1B2B3_all_NP_comp$significativeHSC=="sig_HSC" & B1B2B3_all_NP_comp$state=="KOdown"), fontface = "italic",size=2.5)+
  #xlim(-10,10)+
  #ylim(-10,10)+
  xlab("log2FC ikb KO vs. WT")+
  ylab("log2FC nfkb vs. WT")+
  theme_bw()+
  #facet_wrap(~chromosome_name.y,ncol=2)+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1))+
  xlab("log2FC ikb-1 KO vs WT")+
  #xlab("log2FC nfki-1 (C33A11.1) vs. WT")+
  #ylab("log2FC DKO vs. WT")+
  ylab("log2FC nfki-1 vs. WT")


cor.test(rnatab_L1$log2FoldChange.x,rnatab_L1$log2FoldChange.y)


#3D plots
p <- plot_ly(rnatab, x = ~log2FoldChange.x, y = ~log2FoldChange.y, z = ~log2FoldChange, color = ~sig,marker=list(size=5)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'ikb'),
                      yaxis = list(title = 'nfki'),
                      zaxis = list(title = 'DKO')))
p


## heatmap FC three experiments

rnatab_heat<-subset(rnatab,select=c(id,log2FoldChange.x,log2FoldChange.y,log2FoldChange,sample.x,sample.y,sample,pvalue.x,pvalue.y,pvalue))
rnatab_heat<-melt(rnatab_heat,id.vars = c("id","sample.x","sample.y","sample","pvalue.x","pvalue.y","pvalue"))
rnatab_heat<-melt(rnatab_heat,id.vars = c("id","variable","value","pvalue.x","pvalue.y","pvalue"))
rnatab_heat<-melt(rnatab_heat,id.vars = c("id","variable","value","pvalue.x","pvalue.y","pvalue"))
colnames(rnatab_heat)<-c("id","variable","value","pvalue.x","pvalue.y","pvalue","sample","samp")
rnatab_heat$sample<-NULL

rnatab_heat<-melt(rnatab_heat,id.vars = c("id","variable","value","samp"))
colnames(rnatab_heat)<-c("id","variable","value","sample","pvalue","val")

rnatab_heat$variable<-gsub('\\..*','',rnatab_heat$variable)
rnatab_heat$pvalue<-NULL

#SELECT ONLY GENES WITH p-val 

## TOO MANY GENES FOR A HEATMAP ##
##rnatab_heat<-subset(rnatab_heat,rnatab_heat$val<=0.00000000000000001)

ggplot(rnatab_heat, aes(x=gene,y=variable)) +
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





### Plot RNA-Seq vs. peak
res_df_rna_rpmi<-res_df_rna
res_df_rna_jurkat<-res_df_rna

#peaks bcat genes, including lithium
bcat<-c(BRC1_genes,BRC2_genes,DB1_genes,DB2_genes,DB3_genes,DL2_genes)
bcat<-as.data.frame(table(sort(bcat)))
bcat<-bcat[order(bcat$Freq,decreasing = TRUE),]

#peaks lef genes
unique(lefcomgenes$Var1)

#peaks tcf genes
tcfcomgenes$Var1


rnapeak<-merge(res_df_rna_rpmi,bcatcomgenes,by.x="external_gene_name",by.y="Var1",all=TRUE)
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
  geom_jitter(aes(color=sig))+
  geom_violin(data = subset(rnapeak,rnapeak$sig=="sig"),alpha=0.2)+
  scale_color_manual(values=c("grey","red"))+
  geom_label_repel(aes(label = external_gene_name),color="darkolivegreen", data = subset(rnapeak,rnapeak$bcat=="bcat" & rnapeak$sig=="sig"), fontface = "italic",size=2.5)+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("log2FC bcat vs nobcat")+
  theme_bw()



##### Orthologs C.elegans vs. H sapiens ###

ensembl76_h <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# import ortholog list from different categories
germcell<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_UP/orthologs/germ_cell_orthologlist.txt",sep=",")
germcell$X=NULL
germcell$func<-"germ_cell"

gonad<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_UP/orthologs/gonad_orthologlist.csv",sep=",")
gonad$X=NULL
gonad$func<-"gonad"


nemlarv<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_UP/orthologs/nematod_larva_dev_orthologlist.csv",sep=",")
nemlarv$X=NULL
nemlarv$func<-"nematode_larval_development"

regemb<-read.delim("/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/common_UP/orthologs/regulation_embr_develop.csv",sep=",")
regemb$X=NULL
regemb$func<-"regulation_emrbyonic_development"


uport<-rbind(germcell,gonad,nemlarv,regemb)


library(DataCombine)
uport<-grepl.sub(data=uport,pattern="ENS",Var="Ensembl.Compara")
uport$Ensembl.Compara<-gsub('^ ','',uport$Ensembl.Compara)
#select only one ensembl per gene
uport$Ensembl.Compara<-gsub(' ENS.*','',uport$Ensembl.Compara)


uportlist<-uport$Ensembl.Compara
class(uportlist)

ortup<-getBM(attributes = c("external_gene_name", 'chromosome_name',
                     'start_position', 'end_position','ensembl_gene_id'),
      filters = 'ensembl_gene_id', 
      values = uportlist, 
      mart = ensembl76)

ortab<-merge(ortup,uport,by.x="ensembl_gene_id",by.y="Ensembl.Compara")

write.table(ortab,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/wormenrich_tabs/orhtologs_table.txt",row.names = FALSE, quote = FALSE,sep="\t")



##### L1 ChIPSeq annotations gene ids #####

L1_peaks<-rbind(dataf_peakannoIK27,dataf_peakannoNF27,dataf_peakannoWT27,
                dataf_peakannoIK36,dataf_peakannoNF36)

L1_peaks_genes<-unique(L1_peaks$geneId)

L1_bm_genes_peaks<-getBM(attributes = c('external_gene_name','wormbase_cds', 'chromosome_name',
                            'start_position', 'end_position','ensembl_gene_id'),
             filters = 'wormbase_cds', 
             values = L1_peaks_genes, 
             mart = ensembl76)

L1_peaks<-merge(L1_peaks,L1_bm_genes_peaks,by.x="geneId",by.y = "wormbase_cds",all.x=TRUE)
L1_peaks$V4<-gsub('.projects.*\\/','',L1_peaks$V4)

write.table(L1_peaks,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/L1_peaks_all_chips.tab",quote = FALSE,row.names = FALSE,sep="\t")

##### L4 ChIPSeq annotations gene ids #####

L4_peaks<-rbind(dataf_peakannoHK27ikb,dataf_peakannoHK27nfk,dataf_peakannoHK27WT,
                dataf_peakannoHK36ikb,dataf_peakannoHK36nfk,dataf_peakannoHK36WT)

L4_peaks_genes<-unique(L4_peaks$geneId)

L4_bm_genes_peaks<-getBM(attributes = c('external_gene_name','wormbase_cds', 'chromosome_name',
                                        'start_position', 'end_position','ensembl_gene_id'),
                         filters = 'wormbase_cds', 
                         values = L4_peaks_genes, 
                         mart = ensembl76)

L4_peaks<-merge(L4_peaks,L4_bm_genes_peaks,by.x="geneId",by.y = "wormbase_cds",all.x=TRUE)
L4_peaks$V4<-gsub('.projects.*\\/','',L4_peaks$V4)

write.table(L4_peaks,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/L4_peaks_all_chips.tab",quote = FALSE,row.names = FALSE,sep="\t")


#### FMF and FMM Peaks #####


FMF_FMM_peaks<-rbind(dataf_peakannoFMF1,dataf_peakannoFMF2,
                     dataf_peakannoFMM1,dataf_peakannoFMM2,dataf_peakannoFMM3)

FMF_FMM_peaks_genes<-unique(FMF_FMM_peaks$geneId)

FMF_FMM_peaks_genes_peaks<-getBM(attributes = c('external_gene_name','wormbase_cds', 'chromosome_name',
                                        'start_position', 'end_position','ensembl_gene_id'),
                         filters = 'wormbase_cds', 
                         values = FMF_FMM_peaks_genes, 
                         mart = ensembl76)

FMF_FMM_peaks<-merge(FMF_FMM_peaks,FMF_FMM_peaks_genes_peaks,by.x="geneId",by.y = "wormbase_cds",all.x=TRUE)
FMF_FMM_peaks$V4<-gsub('.projects.*\\/','',FMF_FMM_peaks$V4)


write.table(FMF_FMM_peaks,"/Volumes/grcmc/YGUILLEN/Celegans/Data_R/FMF_FMM_peaks_all_chips.tab",quote = FALSE,row.names = FALSE,sep="\t")

