#Load libraries
library(data.table)
library(GenomicAlignments)
library(GO.db)
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
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(ReactomePA)
#library(trackViewer)
#library(Gviz)
#library(rtracklayer)
#library(Sushi)
library(ChIPpeakAnno)
#library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(DO.db)
library(VennDiagram)
library(readxl)
library(enrichR)
library(scatterpie)
library(Gviz)

source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")

install.packages("rJava",dependencies=TRUE,type="source")
library(xlsx)
library(rJava)
library(venneuler)


setwd("/Volumes/grcmc/YGUILLEN/ATAC-Seq/")

## Blacklist regions of the genome prone to represent artifacts (centromeres, telomers and satellite repeats)
blacklist<-read.delim("db/ENCFF547MET.bed",header = FALSE)
colnames(blacklist)<-c("Chr","start","end")
blacklist$region<-(blacklist$end)-(blacklist$start)

## Artifacts paired reads
samp1<-read.delim("fragment_align/artifacts_2_1_29876_TAAGGC.txt",header = FALSE,sep="")
colnames(samp1)<-c("Mate1","start1","Mate2","start2")

ggplot(samp1,aes(x=start1,y=start2,color=Mate2))+
  geom_point()+
  facet_wrap(~Mate1,scales="free")+
  theme_bw()

## Paired fragment distribution per sample
psamp1<-read.delim("fragment_align/hist_paired_random_2_1_29876_TAAGGC.txt",header = FALSE,sep="")
colnames(psamp1)<-c("Freq","Chr","size")
psamp1$sample<-"Samp_2_1"
psamp1$gender<-"female"

psamp2<-read.delim("fragment_align/hist_paired_random_2_2_29877_CGTACT.txt",header = FALSE,sep="")
colnames(psamp2)<-c("Freq","Chr","size")
psamp2$sample<-"Samp_2_2"
psamp2$gender<-"male"

psamp3<-read.delim("fragment_align/hist_paired_random_2_3_29878_AGGCAG.txt",header = FALSE,sep="")
colnames(psamp3)<-c("Freq","Chr","size")
psamp3$sample<-"Samp_2_3"
psamp3$gender<-"female"

psamp5<-read.delim("fragment_align/hist_paired_random_2_5_29880_GGACTC.txt",header = FALSE,sep="")
colnames(psamp5)<-c("Freq","Chr","size")
psamp5$sample<-"Samp_2_5"
psamp5$gender<-"male"

psamp6<-read.delim("fragment_align/hist_paired_random_2_6_29881_TAGGCA.txt",header = FALSE,sep="")
colnames(psamp6)<-c("Freq","Chr","size")
psamp6$sample<-"Samp_2_6"
psamp6$gender<-"male"

allpsamp<-rbind(psamp1,psamp2,psamp3,psamp5,psamp6)
allpsamp<-subset(allpsamp,allpsamp$size!=0 & allpsamp$size<=1000)

#allpsamp$Chr<-gsub('chr','',allpsamp$Chr)
allpsamp$Chr <- ordered(allpsamp$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrM","chrY"))
allpsamp$Chr<-as.factor(allpsamp$Chr)

allpsamp$Freq<-as.numeric(allpsamp$Freq)

allpsamp_filt<-subset(allpsamp,allpsamp$Chr!="chrM" & allpsamp$Chr!="chrY")


ggplot(allpsamp_filt,aes(x=size,y=Freq))+
  geom_line(aes(color=sample,linetype=gender))+
  #geom_linerange(aes(color=sample,x=size, ymax=Freq, ymin=0),position = position_jitter(height = 0L, seed = 1L),alpha=0.3)+
  scale_color_brewer(palette="Set1")+
  facet_wrap(~Chr,scales="free",ncol=4)+
  theme(axis.text.x = element_blank())+
  theme_bw()

allpsamp_filt<-as.data.table(allpsamp_filt)
distfrag<-allpsamp_filt[ , .SD[which.max(Freq)], by=c("Chr","sample")]

distfrag$cond <- ifelse(distfrag$sample != "Samp_2_5" & distfrag$sample != "Samp_2_6", 
                             c("KO"), c("WT")) 

#Distribution maximum fragment size per chromosome
ggplot(distfrag,aes(x=Chr,y=size,group=sample))+
  geom_point(aes(color=sample,shape=cond),size=7)+
  geom_line(aes(group=sample,color=sample,linetype=gender))+
  theme(axis.text.x = element_blank())+
  theme_bw()

#Overall chromosomes, differences per sample
ggplot(distfrag,aes(x=sample,y=size))+
  geom_point(aes(shape=Chr),size=6)+
  scale_shape_manual(values=1:length(unique(distfrag$Chr)))+
  geom_boxplot(aes(color=cond,linetype=gender),alpha=0)+
  scale_color_manual(values=c("dodgerblue","coral3"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=0.5))

#Overall chromosomes, differences per sample
ggplot(distfrag,aes(x=cond,y=size))+
  geom_point(aes(shape=Chr,color=sample),size=6)+
  scale_shape_manual(values=1:length(unique(distfrag$Chr)))+
  geom_boxplot(alpha=0)+
  #scale_color_manual(values=c("dodgerblue","coral3"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=0.5))

wilcox.test(distfrag$size~distfrag$cond)


## We will consider Broad peaks for merged bam files.
#Peaks from random 30 M reads
peak1 <- readPeakFile("peak_calling_random/random_2_1_29876_TAAGGC_peaks_peaks.broadPeak")
peak2 <- readPeakFile("peak_calling_random/random_2_2_29877_CGTACT_peaks_peaks.broadPeak")
peak3 <- readPeakFile("peak_calling_random/random_2_3_29878_AGGCAG_peaks_peaks.broadPeak")
peakKO <- readPeakFile("peak_calling_random/random_KO_peaks_peaks.broadPeak")

peak5 <- readPeakFile("peak_calling_random/random_2_5_29880_GGACTC_peaks_peaks.broadPeak")
peak6 <- readPeakFile("peak_calling_random/random_2_6_29881_TAGGCA_peaks_peaks.broadPeak")
peakWT<- readPeakFile("peak_calling_random/random_WT_peaks_peaks.broadPeak")

#Peaks from random 14 M reads including ENCODE experiments
peak1 <- readPeakFile("/Volumes/grcmc/YGUILLEN/ATAC-Seq/ENCODE/peakcall/random_ENCODE_2_1_29876_TAAGGC_peaks_peaks.broadPeak")
peak2 <- readPeakFile("/Volumes/grcmc/YGUILLEN/ATAC-Seq/ENCODE/peakcall/random_ENCODE_2_2_29877_CGTACT_peaks_peaks.broadPeak")
peak3 <- readPeakFile("/Volumes/grcmc/YGUILLEN/ATAC-Seq/ENCODE/peakcall/random_ENCODE_2_3_29878_AGGCAG_peaks_peaks.broadPeak")
peak5 <- readPeakFile("/Volumes/grcmc/YGUILLEN/ATAC-Seq/ENCODE/peakcall/random_ENCODE_2_5_29880_GGACTC_peaks_peaks.broadPeak")
peak6 <- readPeakFile("/Volumes/grcmc/YGUILLEN/ATAC-Seq/ENCODE/peakcall/random_ENCODE_2_6_29881_TAGGCA_peaks_peaks.broadPeak")
peakENC1 <- readPeakFile("/Volumes/grcmc/YGUILLEN/ATAC-Seq/ENCODE/peakcall/random_ENCODE_ENCFF202OJI_peaks_peaks.broadPeak")
peakENC2 <- readPeakFile("/Volumes/grcmc/YGUILLEN/ATAC-Seq/ENCODE/peakcall/ENCODE_ENCFF680ZRI_peaks_peaks.broadPeak")


#genes database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#peak annotation (CHECK PEAKS!!!!)

peakAnno1 <- annotatePeak(peak1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno2 <- annotatePeak(peak2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno3 <- annotatePeak(peak3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoKO <- annotatePeak(peakKO,tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnno5 <- annotatePeak(peak5, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno6 <- annotatePeak(peak6, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWT <- annotatePeak(peakWT, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnnoENC1 <- annotatePeak(peakENC1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoENC2 <- annotatePeak(peakENC2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# Visualize genomic annotation

plotAnnoBar(peakAnno1,title = "Peak distribution KO_1")
plotAnnoBar(peakAnno2,title = "Peak distribution KO_2")
plotAnnoBar(peakAnno3,title = "Peak distribution KO_3")
plotAnnoBar(peakAnnoKO,title = "Peak distribution merged KO")

plotAnnoBar(peakAnno5,title = "Peak distribution WT_5")
plotAnnoBar(peakAnno6,title = "Peak distribution WT_6")
plotAnnoBar(peakAnnoWT,title = "Peak distribution merged WT")

plotAnnoBar(peakAnnoENC1,title = "Peak distribution ENC_1")
plotAnnoBar(peakAnnoENC2,title = "Peak distribution ENC_2")
dev.off()


##DATAFRAME WITH genes id linked to peaks
KO_1_peakanno<-as.data.frame(peakAnno1)
KO_1_peakanno$annotation<-gsub('Intron.*','Intron',KO_1_peakanno$annotation)
KO_1_peakanno$annotation<-gsub('Exon.*','Exon',KO_1_peakanno$annotation)
KO_1_peakanno$annotation<-gsub('Downstream.*','Downstream',KO_1_peakanno$annotation)
KO_1_peakanno$annotation<-gsub('Promoter.*','Promoter',KO_1_peakanno$annotation)

KO_2_peakanno<-as.data.frame(peakAnno2)
KO_2_peakanno$annotation<-gsub('Intron.*','Intron',KO_2_peakanno$annotation)
KO_2_peakanno$annotation<-gsub('Exon.*','Exon',KO_2_peakanno$annotation)
KO_2_peakanno$annotation<-gsub('Downstream.*','Downstream',KO_2_peakanno$annotation)
KO_2_peakanno$annotation<-gsub('Promoter.*','Promoter',KO_2_peakanno$annotation)

KO_3_peakanno<-as.data.frame(peakAnno3)
KO_3_peakanno$annotation<-gsub('Intron.*','Intron',KO_3_peakanno$annotation)
KO_3_peakanno$annotation<-gsub('Exon.*','Exon',KO_3_peakanno$annotation)
KO_3_peakanno$annotation<-gsub('Downstream.*','Downstream',KO_3_peakanno$annotation)
KO_3_peakanno$annotation<-gsub('Promoter.*','Promoter',KO_3_peakanno$annotation)


KO_merge_peakanno<-as.data.frame(peakAnnoKO)
KO_merge_peakanno$annotation<-gsub('Intron.*','Intron',KO_merge_peakanno$annotation)
KO_merge_peakanno$annotation<-gsub('Exon.*','Exon',KO_merge_peakanno$annotation)
KO_merge_peakanno$annotation<-gsub('Downstream.*','Downstream',KO_merge_peakanno$annotation)
KO_merge_peakanno$annotation<-gsub('Promoter.*','Promoter',KO_merge_peakanno$annotation)


WT_5_peakanno<-as.data.frame(peakAnno5)
WT_5_peakanno$annotation<-gsub('Intron.*','Intron',WT_5_peakanno$annotation)
WT_5_peakanno$annotation<-gsub('Exon.*','Exon',WT_5_peakanno$annotation)
WT_5_peakanno$annotation<-gsub('Downstream.*','Downstream',WT_5_peakanno$annotation)
WT_5_peakanno$annotation<-gsub('Promoter.*','Promoter',WT_5_peakanno$annotation)

WT_6_peakanno<-as.data.frame(peakAnno6)
WT_6_peakanno$annotation<-gsub('Intron.*','Intron',WT_6_peakanno$annotation)
WT_6_peakanno$annotation<-gsub('Exon.*','Exon',WT_6_peakanno$annotation)
WT_6_peakanno$annotation<-gsub('Downstream.*','Downstream',WT_6_peakanno$annotation)
WT_6_peakanno$annotation<-gsub('Promoter.*','Promoter',WT_6_peakanno$annotation)

WT_merge_peakanno<-as.data.frame(peakAnnoWT)
WT_merge_peakanno$annotation<-gsub('Intron.*','Intron',WT_merge_peakanno$annotation)
WT_merge_peakanno$annotation<-gsub('Exon.*','Exon',WT_merge_peakanno$annotation)
WT_merge_peakanno$annotation<-gsub('Downstream.*','Downstream',WT_merge_peakanno$annotation)
WT_merge_peakanno$annotation<-gsub('Promoter.*','Promoter',WT_merge_peakanno$annotation)


ENC_1_peakanno<-as.data.frame(peakAnnoENC1)
ENC_1_peakanno$annotation<-gsub('Intron.*','Intron',ENC_1_peakanno$annotation)
ENC_1_peakanno$annotation<-gsub('Exon.*','Exon',ENC_1_peakanno$annotation)
ENC_1_peakanno$annotation<-gsub('Downstream.*','Downstream',ENC_1_peakanno$annotation)
ENC_1_peakanno$annotation<-gsub('Promoter.*','Promoter',ENC_1_peakanno$annotation)

ENC_2_peakanno<-as.data.frame(peakAnnoENC2)
ENC_2_peakanno$annotation<-gsub('Intron.*','Intron',ENC_2_peakanno$annotation)
ENC_2_peakanno$annotation<-gsub('Exon.*','Exon',ENC_2_peakanno$annotation)
ENC_2_peakanno$annotation<-gsub('Downstream.*','Downstream',ENC_2_peakanno$annotation)
ENC_2_peakanno$annotation<-gsub('Promoter.*','Promoter',ENC_2_peakanno$annotation)

#unify column names
colnames(WT_merge_peakanno)[6]<-c("Name")
colnames(WT_merge_peakanno)[7]<-c("Score")
colnames(WT_merge_peakanno)[9]<-c("FE")
colnames(WT_merge_peakanno)[10]<-c("log10pval")
colnames(WT_merge_peakanno)[11]<-c("log10qval")

#unify column names
colnames(WT_5_peakanno)[6]<-c("Name")
colnames(WT_5_peakanno)[7]<-c("Score")
colnames(WT_5_peakanno)[9]<-c("FE")
colnames(WT_5_peakanno)[10]<-c("log10pval")
colnames(WT_5_peakanno)[11]<-c("log10qval")

#unify column names
colnames(WT_6_peakanno)[6]<-c("Name")
colnames(WT_6_peakanno)[7]<-c("Score")
colnames(WT_6_peakanno)[9]<-c("FE")
colnames(WT_6_peakanno)[10]<-c("log10pval")
colnames(WT_6_peakanno)[11]<-c("log10qval")


colnames(KO_1_peakanno)[6]<-c("Name")
colnames(KO_1_peakanno)[7]<-c("Score")
colnames(KO_1_peakanno)[9]<-c("FE")
colnames(KO_1_peakanno)[10]<-c("log10pval")
colnames(KO_1_peakanno)[11]<-c("log10qval")

colnames(KO_2_peakanno)[6]<-c("Name")
colnames(KO_2_peakanno)[7]<-c("Score")
colnames(KO_2_peakanno)[9]<-c("FE")
colnames(KO_2_peakanno)[10]<-c("log10pval")
colnames(KO_2_peakanno)[11]<-c("log10qval")

colnames(KO_3_peakanno)[6]<-c("Name")
colnames(KO_3_peakanno)[7]<-c("Score")
colnames(KO_3_peakanno)[9]<-c("FE")
colnames(KO_3_peakanno)[10]<-c("log10pval")
colnames(KO_3_peakanno)[11]<-c("log10qval")

colnames(KO_merge_peakanno)[6]<-c("Name")
colnames(KO_merge_peakanno)[7]<-c("Score")
colnames(KO_merge_peakanno)[9]<-c("FE")
colnames(KO_merge_peakanno)[10]<-c("log10pval")
colnames(KO_merge_peakanno)[11]<-c("log10qval")


#Quality filter peaks, log10qval > 2 (qval < 0.01)
WT_merge_peakanno<-subset(WT_merge_peakanno,WT_merge_peakanno$log10qval>2)
WT_5_peakanno<-subset(WT_5_peakanno,WT_5_peakanno$log10qval>2)
WT_6_peakanno<-subset(WT_6_peakanno,WT_6_peakanno$log10qval>2)
KO_merge_peakanno<-subset(KO_merge_peakanno,KO_merge_peakanno$log10qval>2)
KO_1_peakanno<-subset(KO_1_peakanno,KO_1_peakanno$log10qval>2)
KO_2_peakanno<-subset(KO_2_peakanno,KO_2_peakanno$log10qval>2)
KO_3_peakanno<-subset(KO_3_peakanno,KO_3_peakanno$log10qval>2)

#All chromosomes
length(unique(KO_1_peakanno$SYMBOL))
length(unique(KO_2_peakanno$SYMBOL))
length(unique(KO_3_peakanno$SYMBOL))
length(unique(KO_merge_peakanno$SYMBOL))

length(unique(WT_5_peakanno$SYMBOL))
length(unique(WT_6_peakanno$SYMBOL))
length(unique(WT_merge_peakanno$SYMBOL))

length(unique(ENC_1_peakanno$SYMBOL))
length(unique(ENC_2_peakanno$SYMBOL))


#Discarding X,Y and Mitochondrial, number of genes
length(unique(subset(KO_1_peakanno,KO_1_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))
length(unique(subset(KO_2_peakanno,KO_2_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))
length(unique(subset(KO_3_peakanno,KO_3_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))
length(unique(subset(KO_merge_peakanno,KO_merge_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))

length(unique(subset(WT_5_peakanno,WT_5_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))
length(unique(subset(WT_6_peakanno,WT_6_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))
length(unique(subset(WT_merge_peakanno,WT_merge_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))

length(unique(subset(ENC_1_peakanno,ENC_1_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))
length(unique(subset(ENC_2_peakanno,ENC_2_peakanno$seqnames!=c("chrM","chrY","chrX"))$SYMBOL))

#Discarding X,Y and Mitochondrial, number of peaks
nrow(subset(KO_1_peakanno,KO_1_peakanno$seqnames!=c("chrM","chrY","chrX")))
nrow(subset(KO_2_peakanno,KO_2_peakanno$seqnames!=c("chrM","chrY","chrX")))
nrow(subset(KO_3_peakanno,KO_3_peakanno$seqnames!=c("chrM","chrY","chrX")))
nrow(subset(KO_merge_peakanno,KO_merge_peakanno$seqnames!=c("chrM","chrY","chrX")))

nrow(subset(WT_5_peakanno,WT_5_peakanno$seqnames!=c("chrM","chrY","chrX")))
nrow(subset(WT_6_peakanno,WT_6_peakanno$seqnames!=c("chrM","chrY","chrX")))
nrow(subset(WT_merge_peakanno,WT_merge_peakanno$seqnames!=c("chrM","chrY","chrX")))

nrow(subset(ENC_1_peakanno,ENC_1_peakanno$seqnames!=c("chrM","chrY","chrX")))
nrow(subset(ENC_2_peakanno,ENC_2_peakanno$seqnames!=c("chrM","chrY","chrX")))

#All chromosomes (30M subs, fitering X, Y and M chromosomes)
peakscomp <- data.frame("Sample" = c("KO_1","KO_2","KO_3","KO_merge","WT_5","WT_6","WT_merge"), "Gender"=c("Female","Male","Female","Merge","Male","Male","Merge"),"Condition" = c("KO","KO","KO","KO","WT","WT","WT"),"Unique_genes_with_peaks" = c(22342,17009,22627,22570,14150,11561,1115),"Call"=c("ind","ind","ind","merge","ind","ind","merge"),"Unique_peaks" = c(385137,90768,454973,478609,52885,33848,29952))

#All chromosomes (14M subs ENCODE)
#peakscomp_ENC <- data.frame("Sample" = c("KO_1","KO_2","KO_3","WT_5","WT_6","ENC_1","ENC_2"), "Condition" = c("KO","KO","KO","WT","WT","ENC","ENC"),"Unique_genes_with_peaks" = c(22195,17751,22462,16719,9041,18251,18822),"Call"=c("ind","ind","ind","ind","ind","ind","ind"),"Unique_peaks" = c(366291,104766,411374,82583,20237,61119,71906))


#Melting (check dataframe)
#peakscomp_melt<-melt(peakscomp_ENC,id.vars = c("Sample","Condition","Call"))
peakscomp_melt<-melt(peakscomp,id.vars = c("Sample","Condition","Call","Gender"))

ggplot(peakscomp_melt,aes(x=Condition,y=value,label=Sample))+
  geom_point(color="black",size=7,aes(shape=Gender))+
  geom_point(aes(color = Call,shape=Gender),size=5,alpha=0.8)+
  geom_boxplot(alpha=0)+
  facet_wrap(~variable,scales="free_y")+
  #geom_text(aes(label=Sample), hjust=0,vjust=-1.2,angle=30,colour = "black", fontface = "italic",size=3)+
  scale_color_manual(values=c("dodgerblue","coral3","darkolivegreen"))+
  #scale_color_brewer()+
  geom_label_repel(aes(label=Sample), size=3,hjust=1.5,colour = "black", fontface = "italic")+
  theme_light()+
  expand_limits(x = 0, y = 0)+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  xlab("")+ylab("Number of peaks")+
  ggtitle("Broad peaks in KO and WT excluding mtDNA and sex chromosomes")

kruskal.test(peakscomp$Unique_peaks,peakscomp$Condition)
#pairwise.wilcox.test(peakscomp_ENC$Unique_peaks,peakscomp_ENC$Condition)

#distribution of peak regions
distKO1<-as.data.frame(t(as.data.frame(prop.table(table(KO_1_peakanno$annotation)))))
colnames(distKO1) <- as.character(unlist(distKO1[1,]))
distKO1<-distKO1[-1,]
distKO1$Sample<-"KO_1"
#check ENCODE dataframe! peakscomp (30M) or peakscomp_ENC
distKO1<-merge(distKO1,peakscomp,by="Sample")


distKO2<-as.data.frame(t(as.data.frame(prop.table(table(KO_2_peakanno$annotation)))))
colnames(distKO2) <- as.character(unlist(distKO2[1,]))
distKO2<-distKO2[-1,]
distKO2$Sample<-"KO_2"
#check ENCODE dataframe! peakscomp (30M) or peakscomp_ENC
distKO2<-merge(distKO2,peakscomp,by="Sample")

distKO3<-as.data.frame(t(as.data.frame(prop.table(table(KO_3_peakanno$annotation)))))
colnames(distKO3) <- as.character(unlist(distKO3[1,]))
distKO3<-distKO3[-1,]
distKO3$Sample<-"KO_3"
#check ENCODE dataframe! peakscomp (30M) or peakscomp_ENC
distKO3<-merge(distKO3,peakscomp,by="Sample")

distKOmerg<-as.data.frame(t(as.data.frame(prop.table(table(KO_merge_peakanno$annotation)))))
colnames(distKOmerg) <- as.character(unlist(distKOmerg[1,]))
distKOmerg<-distKOmerg[-1,]
distKOmerg$Sample<-"KO_merge"
distKOmerg<-merge(distKOmerg,peakscomp,by="Sample")

distWT5<-as.data.frame(t(as.data.frame(prop.table(table(WT_5_peakanno$annotation)))))
colnames(distWT5) <- as.character(unlist(distWT5[1,]))
distWT5<-distWT5[-1,]
distWT5$Sample<-"WT_5"
#check ENCODE dataframe! peakscomp (30M) or peakscomp_ENC
distWT5<-merge(distWT5,peakscomp,by="Sample")

distWT6<-as.data.frame(t(as.data.frame(prop.table(table(WT_6_peakanno$annotation)))))
colnames(distWT6) <- as.character(unlist(distWT6[1,]))
distWT6<-distWT6[-1,]
distWT6$Sample<-"WT_6"
#check ENCODE dataframe! peakscomp (30M) or peakscomp_ENC
distWT6<-merge(distWT6,peakscomp,by="Sample")

distWTmerg<-as.data.frame(t(as.data.frame(prop.table(table(WT_merge_peakanno$annotation)))))
colnames(distWTmerg) <- as.character(unlist(distWTmerg[1,]))
distWTmerg<-distWTmerg[-1,]
distWTmerg$Sample<-"WT_merge"
distWTmerg<-merge(distWTmerg,peakscomp,by="Sample")

distENC1<-as.data.frame(t(as.data.frame(prop.table(table(ENC_1_peakanno$annotation)))))
colnames(distENC1) <- as.character(unlist(distENC1[1,]))
distENC1<-distENC1[-1,]
distENC1$Sample<-"ENC_1"
#check ENCODE dataframe! peakscomp (30M) or peakscomp_ENC
distENC1<-merge(distENC1,peakscomp_ENC,by="Sample")

distENC2<-as.data.frame(t(as.data.frame(prop.table(table(ENC_2_peakanno$annotation)))))
colnames(distENC2) <- as.character(unlist(distENC2[1,]))
distENC2<-distENC2[-1,]
distENC2$Sample<-"ENC_2"
#check ENCODE dataframe! peakscomp (30M) or peakscomp_ENC
distENC2<-merge(distENC2,peakscomp_ENC,by="Sample")


distall<-rbind(distKO1,distKO2,distKO3,distKOmerg,distWT5,distWT6,distWTmerg)

#For ENCODE
distall<-rbind(distKO1,distKO2,distKO3,distWT5,distWT6,distENC1,distENC2)

colnames(distall)[2]<-c("Three_UTR")
colnames(distall)[3]<-c("Five_UTR")
#colnames(distall)[4]<-c("Distal_intergenic")
#colnames(distall)[5]<-c("Downstream_1kb")
#colnames(distall)[6]<-c("Downstream_1_2_kb")
#colnames(distall)[7]<-c("Downstream_2_3_kb")
#colnames(distall)[10]<-c("Promoter_1kb")
#colnames(distall)[11]<-c("Promoter_1_2kb")
#colnames(distall)[12]<-c("Promoter_2_3kb")

distall_num <- data.frame(sapply(distall[,2:8], function(x) as.numeric(as.character(x))*100))

distall_num<-cbind(distall$Sample,distall_num,distall$Condition,distall$Unique_genes_with_peaks,distall$Call,distall$Unique_peaks)

colnames(distall_num)[1]<-"Sample"
colnames(distall_num)[9]<-"Condition"
colnames(distall_num)[10]<-"Unique_genes_with_peaks"
colnames(distall_num)[11]<-"Call"
colnames(distall_num)[12]<-"Unique_peaks"

names<-colnames(distall_num[c(2:8)])

distall_num_melt<-melt(distall_num,id.vars = c("Sample","Condition","Unique_genes_with_peaks","Call","Unique_peaks"))

ggplot(distall_num_melt,aes(x=variable,y=Unique_peaks,label=round(value,2)))+
  geom_point(aes(size = value*3),color="black") +
  geom_point(aes(size=value,color = variable))+
  #geom_boxplot(alpha=0)+
  facet_wrap(~Condition)+
  geom_text(aes(label=round(value,2)), hjust=0,vjust=-1.2,angle=30,colour = "black", fontface = "italic",size=3)+
  #scale_color_manual(values=c("dodgerblue","coral3"))+
  scale_color_brewer()+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  xlab("")+ylab("Number of peaks")

ggplot(distall_num_melt,aes(x=Sample,y=value,fill=variable))+
  geom_bar(stat="identity") +
  #geom_point(aes(size=value,color = Condition))+10^
  #geom_boxplot(alpha=0)+
  facet_wrap(~Condition,scale="free")+
  #geom_text(aes(label=round(value,2)), hjust=0,vjust=-1.2,angle=30,colour = "black", fontface = "italic",size=3)+
  #scale_color_manual(values=c("dodgerblue","coral3","darkolivegreen"))+
  scale_fill_brewer(palette="Blues")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  xlab("")+ylab("Fraction of peaks")



#Unique peaks and genes in WT and KO regions

WTmerge_genes<-unique(WT_merge_peakanno$SYMBOL)
WT_5_genes<-unique(WT_5_peakanno$SYMBOL)
WT_6_genes<-unique(WT_6_peakanno$SYMBOL)
KOmerge_genes<-unique(KO_merge_peakanno$SYMBOL)
KO_1_genes<-unique(KO_1_peakanno$SYMBOL)
KO_2_genes<-unique(KO_2_peakanno$SYMBOL)
KO_3_genes<-unique(KO_3_peakanno$SYMBOL)

#genes
comgenes_WT<-as.data.frame(table(sort(c(WT_5_genes,WT_6_genes,WTmerge_genes))))
comgenes_KO<-as.data.frame(table(sort(c(KO_1_genes,KO_2_genes,KO_3_genes,KOmerge_genes))))


#venn diagram common peaks in KO replicates
### venn diagram ##
geneLists <- list(KO1 = KO_1_genes, KO2 = KO_2_genes, KO3 = KO_3_genes, KO_merge = KOmerge_genes)
geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off()
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkblue","red","darkgreen","grey"), alpha=c(0.3,0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("KO_1", "KO_2","KO_3","KO_merg"), main="Genes")

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
length(inters[[5]])
length(inters[[6]])
length(inters[[7]])
length(inters[[8]])
length(inters[[9]])
length(inters[[10]])
length(inters[[11]])
length(inters[[12]])
length(inters[[13]])


femgenepeaks<-inters[[6]]
malegenepeaks<-inters[[3]]

genebias<-c(femgenepeaks,malegenepeaks)
sort(genebias)

#select genebias genes from peaks dataframesf
KO_1_bias<-as.character((KO_1_peakanno[with(KO_1_peakanno, (KO_1_peakanno$SYMBOL) %in% genebias),])$Name)
KO_2_bias<-as.character((KO_2_peakanno[with(KO_2_peakanno, (KO_2_peakanno$SYMBOL) %in% genebias),])$Name)
KO_3_bias<-as.character((KO_3_peakanno[with(KO_3_peakanno, (KO_3_peakanno$SYMBOL) %in% genebias),])$Name)
KO_bias<-as.character((KO_merge_peakanno[with(KO_merge_peakanno, (KO_merge_peakanno$SYMBOL) %in% genebias),])$Name)
WT_5_bias<-as.character((WT_5_peakanno[with(WT_5_peakanno, (WT_5_peakanno$SYMBOL) %in% genebias),])$Name)
WT_6_bias<-as.character((WT_6_peakanno[with(WT_6_peakanno, (WT_6_peakanno$SYMBOL) %in% genebias),])$Name)
WT_bias<-as.character((WT_merge_peakanno[with(WT_merge_peakanno, (WT_merge_peakanno$SYMBOL) %in% genebias),])$Name)

peaks_gender_bias<-c(KO_1_bias,KO_2_bias,KO_3_bias,KO_bias,WT_5_bias,WT_6_bias,WT_5_bias)
peaks_gender_bias<-gsub('.projects.*random','random',peaks_gender_bias)
peaks_gender_bias<-unique(peaks_gender_bias)
length(peaks_gender_bias)

#write output list
write(peaks_gender_bias,"peaks_gender_bias.txt")

#venn diagram common peaks in WT replicates
### venn diagram ##
geneLists <- list(WT_5 = WT_5_genes, WT_6 = WT_6_genes, WT_merg = WTmerge_genes)
geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off()
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkblue","red","darkgreen"), alpha=c(0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("WT_5","WT_6","WT_merge"), main="Genes")

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
length(inters[[5]])
length(inters[[6]])
length(inters[[7]])


comgenes<-as.data.frame(table(sort(c(WTmerge_genes,KOmerge_genes))))

comgenes<-comgenes[order(comgenes$Freq,decreasing = TRUE),]
comgenes_all<-subset(comgenes,comgenes$Freq==2)
uniqueenes_all<-subset(comgenes,comgenes$Freq==1)

#genes
uniquegenes_WT<-setdiff(WTmerge_genes,comgenes_all$Var1)
uniquegenes_KO<-setdiff(KOmerge_genes,comgenes_all$Var1)


#substract gender biases peaks from KO
KO_merge_peakanno_nogen<-KO_merge_peakanno[with(KO_merge_peakanno, !(KO_merge_peakanno$SYMBOL) %in% genebias),]
WT_merge_peakanno_nogen<-WT_merge_peakanno[with(WT_merge_peakanno, !(WT_merge_peakanno$SYMBOL) %in% genebias),]

#genes
uniquegenes_WT_sites<-WT_merge_peakanno_nogen[with(WT_merge_peakanno_nogen,WT_merge_peakanno_nogen$SYMBOL %in% uniquegenes_WT),]
uniquegenes_KO_sites<-KO_merge_peakanno_nogen[with(KO_merge_peakanno_nogen, KO_merge_peakanno_nogen$SYMBOL %in% uniquegenes_KO),]

uniqueWT<-unique(uniquegenes_WT_sites$SYMBOL)
uniqueKO<-unique(uniquegenes_KO_sites$SYMBOL)

colnames(uniquegenes_WT_sites)[6]<-c("Name")
colnames(uniquegenes_WT_sites)[7]<-c("Score")
colnames(uniquegenes_WT_sites)[9]<-c("FE")
colnames(uniquegenes_WT_sites)[10]<-c("log10pval")
colnames(uniquegenes_WT_sites)[11]<-c("log10qval")

colnames(uniquegenes_KO_sites)[6]<-c("Name")
colnames(uniquegenes_KO_sites)[7]<-c("Score")
colnames(uniquegenes_KO_sites)[9]<-c("FE")
colnames(uniquegenes_KO_sites)[10]<-c("log10pval")
colnames(uniquegenes_KO_sites)[11]<-c("log10qval")

uniquegenes_WT_sites<-uniquegenes_WT_sites[rev(order(uniquegenes_WT_sites$log10pval)),]
uniquegenes_KO_sites<-uniquegenes_KO_sites[rev(order(uniquegenes_KO_sites$log10pval)),]

WT_nobias<-unique(WT_merge_peakanno_nogen$SYMBOL)
KO_nobias<-unique(KO_merge_peakanno_nogen$SYMBOL)

# Distribution unique peaks
nrow(uniquegenes_WT_sites)
table(uniquegenes_WT_sites$annotation)

nrow(uniquegenes_KO_sites)
table(uniquegenes_KO_sites$annotation)

#venn diagram common peaks in KO replicates
### venn diagram ##
geneLists <- list( WT = WT_nobias, KO =KO_nobias)
geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off()
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkblue","grey"), alpha=c(0.3,0.3), cex = 2, cat.fontface=4, category.names=c("WT","KO"), main="Genes")

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

KO_exclusive<-inters[[1]]
KO_exclusive<-data.frame(KO_exclusive)
KO_exclusive$peak<-"KO"
names(KO_exclusive)<-c("Row.names","peak")

WT_exclusive<-inters[[2]]
WT_exclusive<-data.frame(WT_exclusive)
WT_exclusive$peak<-"WT"
names(WT_exclusive)<-c("Row.names","peak")

common_KOWT<-inters[[3]]
common_KOWT<-data.frame(common_KOWT)
common_KOWT$peak<-"common"
names(common_KOWT)<-c("Row.names","peak")

#Cross with RNAseq previous script
B12_comp_atac<-merge(B12_comp,KO_exclusive,by="Row.names",all=TRUE)
B12_comp_atac<-merge(B12_comp_atac,WT_exclusive,by="Row.names",all=TRUE)
B12_comp_atac<-merge(B12_comp_atac,common_KOWT,by="Row.names",all=TRUE)
B12_comp_atac$peak<-paste(B12_comp_atac$peak,B12_comp_atac$peak.x,B12_comp_atac$peak.y)


#Gata1 targets
Gata1<-read.delim("Gata1_targets.txt",header = FALSE)
Gata1$V1<-as.character(Gata1$V1)
Gata1$V1<-paste0(toupper(substr(Gata1$V1, 1, 1)), tolower(substr(Gata1$V1, 2, nchar(Gata1$V1))))
colnames(Gata1)<-c("Row.names")
Gata1$target<-"Gata1"

B12_comp_atac<-merge(B12_comp_atac,Gata1,by="Row.names",all=TRUE)

ggplot(B12_comp_atac,aes(x=peak,y=log2FoldChange_B1,label=Row.names))+
  geom_jitter(aes(color=state))+
  geom_violin(data = subset(B12_comp_atac,B12_comp_atac$pvalue.x<=0.05),alpha=0.0)+
  #geom_label_repel(data = subset(B12_comp_atac,B12_comp_atac$target=="Gata1"),colour = "black", fontface = "italic",size=2.5)+
  #geom_violin(alpha=0.2)+
  #scale_color_manual(values=c("grey","purple","grey","chocolate2","grey","cyan4"))+
  ylab("log2FC KO vs WT")+
  theme_bw()

ggplot(B12_comp_atac,aes(x=peak,y=log2FoldChange_NP,label=Row.names))+
  geom_jitter()+
  #geom_violin(data = subset(B12_comp_atac,B12_comp_atac$padj.x<=0.1),alpha=0.0)+
  geom_violin(alpha=0.0)+
  geom_label_repel(data = subset(B12_comp_atac,B12_comp_atac$peak=="NA KO NA"),colour = "black", fontface = "italic",size=2.5)+
  #geom_violin(alpha=0.2)+
  scale_color_manual(values=c("red","cyan4"))+
  ylab("log2FC KO vs WT")+
  theme_bw()


##unique genes in KO functional enrichment
## GO terms and pathways
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018","TF_Perturbations_Followed_by_Expression")
enrichedWT <- enrichr(uniqueWT, dbs)
enrichedKO <- enrichr(uniqueKO,dbs)

WTGO <- enrichedWT[["GO_Biological_Process_2018"]]
WTGO$peak<-c("WT")

KOGO <- enrichedKO[["GO_Biological_Process_2018"]]
KOGO$peak<-c("KO")

allGO<-rbind(WTGO,KOGO)
allGO$db<-c("GO")

library(tidyr)
allGO<-separate(data = allGO, col = Overlap, into = c("counts", "pathway"), sep = "/")
allGO$counts<-as.numeric(allGO$counts)
allGO<-subset(allGO,allGO$counts>=3)

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

KOTF <- enrichedKO[["TF_Perturbations_Followed_by_Expression"]]
KOTF$peak<-c("KO")

allTF<-rbind(WTTF,KOTF)
allTF$db<-c("TF")


bpsub<-subset(allTF,allTF$P.value<=0.05)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)
class(bpsub$Combined.Score)

bpsub$Genes<-gsub(';','/',bpsub$Genes)
bpsub$Term<-gsub('.GO.*','',bpsub$Term)

bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$peak, bpsub$P.value)), ]$Term)

ggplot(bpsub,aes(x=Term,y=-(log(P.value)),fill=peak))+
  geom_bar(position="dodge",stat = "identity",color="black")+
  #geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Blues")+
  theme(axis.text.y =element_text(size=9),axis.title=element_text(size=14,face="bold"),axis.text.x = element_text(size=12))+
  coord_flip()

dev.off()


## Number of peaks per gene ##
peak_unique_KO<-as.data.frame(table(uniquegenes_KO_sites$SYMBOL))
peak_unique_KO<- peak_unique_KO[rev(order(peak_unique_KO$Freq)),]
peak_unique_KO$condition<-"KO"

peak_unique_WT<-as.data.frame(table(uniquegenes_WT_sites$SYMBOL))
peak_unique_WT<- peak_unique_WT[rev(order(peak_unique_WT$Freq)),]
peak_unique_WT$condition<-"WT"

peak_dist_gene<-rbind(peak_unique_KO,peak_unique_WT)


#merging from RNAseq
peak_dist_gene<-merge(peak_dist_gene,alltab,by.x="Var1",by.y="Row.names",all=TRUE)

peak_dist_gene$significative[is.na(peak_dist_gene$significative)]<-"noDEG"

ggplot(peak_dist_gene,aes(x=condition,y=log(Freq),label=Var1))+
  geom_jitter(aes(color=significative))+
  geom_boxplot(alpha=0)+
  geom_label_repel(aes(label = Var1,color=significative), data = subset(peak_dist_gene,peak_dist_gene$significative!="noDEG" & peak_dist_gene$Freq>12), fontface = "italic",size=2.5)+
  theme_bw()



### DIFFERENTIAL ACCESSIBLE CHROMATIN SITES WITH DESEQ2 ###

##  merged peaks KO WT macs2 annotation
setwd("/Volumes/grcmc/YGUILLEN/ATAC-Seq/HTSEQ/peak_merged/")
output.dir="/Volumes/grcmc/YGUILLEN/ATAC-Seq/HTSEQ/peak_merged/"

# Import metadata files
#metadata does not include bcat KO
#metadata.txt = all peaks
#metadata.txt2 = files substracting from htseq those peaks corresponding to gender-biased genes and high p-value.
#metadata_filter.txt = htseq calling extracting low threshold peaks and sex-biased

#metadata<-read.delim("metadata.txt",header = FALSE)
metadata<-read.delim("metadata_filter.txt",header = FALSE)
colnames(metadata)<-c("sampleID","countFile","condition","gender")


## GTF mm10 gene CDS annotations ###
setwd("/Volumes/grcmc/YGUILLEN/ATAC-Seq/HTSEQ/mm10_gene/")
output.dir="/Volumes/grcmc/YGUILLEN/ATAC-Seq/HTSEQ/mm10_gene"

# Import metadata files
#metadata does not include bcat KO
metadata<-read.delim("mm10_gene/metadata.txt",header = FALSE)
colnames(metadata)<-c("sampleID","countFile","condition","gender")


metadata$group <- factor(paste0(metadata$gender, metadata$condition))




## common steps
#change metadata or metadata2

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = metadata,
                                          directory = output.dir,
                                          design = ~ condition)

rowData(DESeq2Table)
mcols(rowData(DESeq2Table))

# Add additional annotation information using biomaRt
ensembl76 <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "description","external_gene_name"),
            values= rownames(DESeq2Table),
            filter="mgi_symbol",
            mart=ensembl76) 
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
keep <- rowSums(counts(DESeq2Table)) >=10
DESeq2Table <- DESeq2Table[keep,]

GeneCounts <- counts(DESeq2Table)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

### make sure to get fold change WT-deletion
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$condition, levels = c("KO", "WT"))
colData(DESeq2Table)$group <- factor(colData(DESeq2Table)$group, levels = c("FKO","MWT","MKO"))
colData(DESeq2Table)$gender <- factor(colData(DESeq2Table)$gender, levels = c("F", "M"))


#### estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)


### produce rlog-transformed data
rld <- rlogTransformation(DESeq2Table, blind=TRUE) ## create a distance matrix between the samples

## Create PCA
DESeq2::plotPCA(rld, intgroup=c("group"))
dev.off()

## Create heatmap
pdf("HeatmapPlots_filter.pdf")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

dev.off()

## Create PCA 
pcaData <- plotPCA(rld, intgroup=c("group"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

#pcaData$type<-paste(pcaData$condition,pcaData$population,sep="_")
ggplot(pcaData, aes(PC1, PC2,label=name)) +
  #geom_point(size=7,color="black",alpha=0.5)+
  geom_point(size=5,aes(color=group)) +
  geom_text_repel(aes(label=name),hjust=1, vjust=0)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #scale_color_manual(values=c("darkred","darkolivegreen4"))+
  theme_bw()+
  coord_fixed()

dev.off()



### Differential EXPRESSION ANALYSIS ###
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

# Statistical testing of DE genes
levels(DESeq2Table$group)
#levels(DESeq2Table$population)

design(DESeq2Table)
#design(DESeq2Table) <- formula(~ population + condition)

dds<-DESeq(DESeq2Table)

resultsNames(dds)

#Example differences genes between groups
plotCounts(dds, gene=which.min(resatac$pvalue), intgroup="group")

#results considering only males or females, DEGs KO vs WT
#resmales<-results(dds, contrast=c("group", "MWT", "MKO"),alpha=0.1)
#resfemales<-results(dds,contrast=c("group","FWT","FKO"),alpha=0.1)
resatac<-results(dds, contrast=c("condition", "KO", "WT"),alpha=0.1)

#results DEGs between males and females
#resgender<-results(dds,contrast=c("group","FWT","MWT"),alpha=0.1)

#head(resmales)
#head(resfemales)
#head(resgender)
head(resatac)

summary(resgender)
summary(resmales)
summary(resfemales)
summary(resatac)

#genes with < 0.05
sum(resgender$padj < 0.1, na.rm=TRUE)
sum(resatac$pvalue < 0.05, na.rm=TRUE)
sum(resfemales$padj < 0.1, na.rm=TRUE)

plotMA(resgender)
plotMA(resfemales)
plotMA(resmales)
plotMA(resatac)

#Heatmap most differentially expressed genes (only ranked by DEGs)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", dendrogram="column",trace = "none", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( FKO="darkgreen", MWT="blue",MKO="black" )[
             colData(rld)$group ],cexCol=0.3)

# Selecting specific DEGs
#To dataframe (ckeck females, males or gender)
res_df_atac<-as.data.frame(resatac)
res_df_sig_atac<-subset(res_df_atac,res_df_atac$pvalue<0.01)
res_df_sig_atac<- res_df_sig_atac[order(res_df_sig_atac$pvalue),]

UPgenes<-row.names(subset(res_df_sig,res_df_sig$log2FoldChange>0))
DOWNgenes<-row.names(subset(res_df_sig,res_df_sig$log2FoldChange<0))

#Manual heatmap from log trasformed data, using all samples
rld_df<-as.data.frame(assay(rld))

# Selecting only WT samples
#WTsamples<-(subset(metadata_A,metadata_A$condition=="WT"))$sampleID

# Selecting only male samples
Malesamples<-(subset(metadata,metadata$gender=="M"))$sampleID

# Selecting only female samples
#Femalesamples<-(subset(metadata_A,metadata_A$gender=="F"))$sampleID

wt_rld<-which(colnames(rld_df) %in% Malesamples)
wt_rld<-rld_df[,wt_rld]
#colnames(wt_rld)<-c("F","M","M","M")
#colnames(wt_rld)<-c("KO","KO","KO","WT","WT","WT")
colnames(wt_rld)<-c("KO","WT","WT")

#Selecting DEGs gender from wt_rld
deggenes<-row.names(res_df_sig)
deg_rld<-which(rownames(wt_rld) %in% deggenes)
deg_rld<-wt_rld[deg_rld,]

deg_rld<-merge(deg_rld,res_df_sig,by="row.names")
deg_rld<- deg_rld[order(deg_rld$pvalue),]
row.names(deg_rld)<-deg_rld$Row.names

heatmap.2( as.matrix(deg_rld[,2:4]), scale="row", dendrogram="column",trace="none", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( KO="grey", WT="purple",WT.1="purple" )[
             colnames(wt_rld) ],cexCol=0.3)

DEGsfemales<-deggenes
DEGsmales<-deggenes








## Annotate peaks
allpeaks<-rbind(WT_merge_peakanno,KO_merge_peakanno)
allpeaks<-subset(allpeaks,select=c(Name,SYMBOL,GENENAME,annotation))
allpeaks$Name<-gsub('.projects.*peakcall.','',allpeaks$Name)


row.names(res_df)<-gsub('_random.*','',row.names(res_df))

#merge annotations
res_df<-merge(allpeaks,res_df,by.x="Name",by.y="row.names")
res_df <- res_df[order(res_df$pvalue),]

res_df_sig<-subset(res_df,res_df$pvalue<0.05)
res_df_sig <- res_df_sig[order(res_df_sig$SYMBOL,res_df_sig$annotation),]
length(unique(res_df_sig$SYMBOL))

## Number of peaks per gene ##
peak_freq<-as.data.frame(table(res_df_sig$SYMBOL))
peak_freq<- peak_freq[rev(order(peak_freq$Freq)),]

### Intersect with RNAseq table
res_df_sig_rnaD<-res_df_sig[with(res_df_sig, res_df_sig$SYMBOL %in% D_KO_genes),]
upatac_downrna<-subset(res_df_sig_rnaD,res_df_sig_rnaD$log2FoldChange>0)
nrow(upatac_downrna)
downatac_downrna<-subset(res_df_sig_rnaD,res_df_sig_rnaD$log2FoldChange<0)
nrow(downatac_downrna)

res_df_sig_rnaU<-res_df_sig[with(res_df_sig, res_df_sig$SYMBOL %in% U_KO_genes),]
upatac_uprna<-subset(res_df_sig_rnaU,res_df_sig_rnaU$log2FoldChange>0)
nrow(upatac_uprna)
downatac_uprna<-subset(res_df_sig_rnaU,res_df_sig_rnaU$log2FoldChange<0)
nrow(downatac_uprna)

table(upatac_downrna$annotation)
table(upatac_uprna$annotation)
table(downatac_downrna$annotation)
table(downatac_uprna$annotation)

# Functional for atacseq
upgenes<-unique((subset(res_df_sig,res_df_sig$log2FoldChange>0))$SYMBOL)
downgenes<-unique((subset(res_df_sig,res_df_sig$log2FoldChange<0))$SYMBOL)

upatac_uprna$SYMBOL


#functional for combination atacseq and rnaseq
unique(res_df_sig_rnaD$SYMBOL)
unique(res_df_sig_rnaU$SYMBOL)

## GO terms and pathways
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018","GO_Molecular_Function_2018","TF_Perturbations_Followed_by_Expression","Chromosome_Location")
downKO <- enrichr(unique(res_df_sig_rnaD$SYMBOL), dbs)
upKO <- enrichr(unique(res_df_sig_rnaU$SYMBOL),dbs)

KOdown <- downKO[["GO_Biological_Process_2018"]]
KOdown$source<-c("downKO")

KOup <- upKO[["GO_Biological_Process_2018"]]
KOup$source<-c("upKO")

allGO<-rbind(KOdown,KOup)
allGO$db<-c("GO")

allGO<-separate(data = allGO, col = Overlap, into = c("counts", "pathway"), sep = "/")
#allGO<-subset(allGO,allGO$counts>=5)

bpsub<-subset(allGO,allGO$P.value<=0.01)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)
class(bpsub$Combined.Score)

bpsub$Genes<-gsub(';','/',bpsub$Genes)
bpsub$Term<-gsub('.GO.*','',bpsub$Term)

bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$source, bpsub$P.value)), ]$Term)
ggplot(bpsub,aes(x=Term,y=-(log(P.value)),fill=source))+
  geom_bar(position="dodge",stat = "identity",color="black")+
  geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Blues")+
  theme(axis.text.y =element_text(size=9),axis.title=element_text(size=14,face="bold"),axis.text.x = element_text(size=12))+
  coord_flip()





