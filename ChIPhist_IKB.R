
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

#setwd("/Volumes/grcmc/YGUILLEN/LMarrueChIP/")
setwd("/Volumes/cancer/IKB/peakcall/")

## Table summary input reads (WT and KO in adults, in LMarrueChIP/)
#projtab<-read.delim("samples_tab.txt",header = TRUE)

## List of genes intestinal stem cell signature
scsignature<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/stem_cell_signature_list.txt")
fetalsignature_up<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/gene_signature_fetal_UP.txt")
fetalsignature_up$dir<-c("up")
fetalsignature_down<-read.delim("/Volumes/grcmc/YGUILLEN/LMarrueChIP/gene_signature_fetal_DOWN.txt")
fetalsignature_down$dir<-c("down")
fetalsignature<-rbind(fetalsignature_down,fetalsignature_up)

################## MACS2 PEAK CALLING with ChIPseeker ###############

##### ADULTS ###### 

## ADULT K27H3 WT
peakWT1H3 <- readPeakFile("W1H3_peaks_summits.bed")
peakWT2H3 <- readPeakFile("W2H3_peaks_summits.bed")

# We will consider Broad peaks for merged bam files
peakWTH3 <- readPeakFile("WH3_peaks_broad_peaks.broadPeak")


## K27H3 KO
peakK8H3 <- readPeakFile("K8H3_peaks_summits.bed")
peakK9H3 <- readPeakFile("K9H3_peaks_summits.bed")

## We will consider Broad peaks for merged bam files
peakKH3 <- readPeakFile("KH3_peaks_broad_peaks.broadPeak")


######### FETAL #######

# WT
peakFWT1 <- readPeakFile("CWT1_33333_peaks_peaks.broadPeak")
peakFWT2 <- readPeakFile("CWT2_33334_peaks_peaks.broadPeak")

# We will consider Broad peaks for merged bam files
peakFWT <- readPeakFile("CWT_3333_peaks_peaks.broadPeak")

## K27H3 KO
peakFKO1 <- readPeakFile("CKO1_33335_peaks_peaks.broadPeak")
peakFKO2 <- readPeakFile("CKO2_33336_peaks_peaks.broadPeak")

## We will consider Broad peaks for merged bam files
peakFKO <- readPeakFile("CKO_3333_peaks_peaks.broadPeak")

#Overall all chromosmes peak coverage
covplot(peakWT1H3, weightCol="V5")
covplot(peakWT2H3, weightCol="V5")
covplot(peakWTH3)


#Overall all chromosmes peak coverage
covplot(peakK81H3, weightCol="V5")
covplot(peakK91H3, weightCol="V5")
covplot(peakKH3)

peak=GenomicRanges::GRangesList(KF=peakFKO,WF=peakFWT)
covplot(peak)

dev.off()

# peaks to consider
#peakWTH3
#peakKH3
#peakFWT
#peakFKO



#genes database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

### PEAK ANNOTATION WT ADULTS
peakAnnoWT1H3 <- annotatePeak(peakWT1H3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWT2H3 <- annotatePeak(peakWT2H3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoWTH3 <- annotatePeak(peakWTH3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

### PEAK ANNOTATION KO ADULTS
peakAnnoK8H3 <- annotatePeak(peakK8H3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoK9H3 <- annotatePeak(peakK8H3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoKH3 <- annotatePeak(peakKH3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

### PEAK ANNOTATION WT FETAL
peakAnnoFWT1 <- annotatePeak(peakFWT1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFWT2 <- annotatePeak(peakFWT2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFWT <- annotatePeak(peakFWT, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

### PEAK ANNOTATION KO FETAL
peakAnnoFKO1 <- annotatePeak(peakFKO1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFKO2 <- annotatePeak(peakFKO2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoFKO <- annotatePeak(peakFKO, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")


# Visualize genomic annotation
plotAnnoPie(peakAnnoWTH3)
dev.off()
plotAnnoBar(peakAnnoWTH3,title = "Peak distribution WTH3")
dev.off()
vennpie(peakAnnoWTH3)
dev.off()

plotAnnoPie(peakAnnoKH3)
dev.off()
plotAnnoBar(peakAnnoKH3,title = "Peak distribution KOH3")
dev.off()
vennpie(peakAnnoKH3)
dev.off()

upsetplot(peakAnnoWTH3, vennpie=TRUE)
dev.off()
upsetplot(peakAnnoKH3, vennpie=TRUE)
dev.off()

plotDistToTSS(peakAnnoWTH3,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

plotDistToTSS(peakAnnoKH3,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


##DATAFRAME WITH genes id linked to peaks
dataf_peakannoWTH3<-as.data.frame(peakAnnoWTH3)
dataf_peakannoWTH3$peaksite<-paste(dataf_peakannoWTH3$SYMBOL,dataf_peakannoWTH3$annotation,sep='_')
colnames(dataf_peakannoWTH3)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakannoWTH3<- dataf_peakannoWTH3[rev(order(dataf_peakannoWTH3$FE)),]
dataf_peakannoWTH3$annotation<-gsub('Intron.*','Intron',dataf_peakannoWTH3$annotation)
dataf_peakannoWTH3$annotation<-gsub('Exon.*','Exon',dataf_peakannoWTH3$annotation)
dataf_peakannoWTH3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoWTH3$annotation)
dataf_peakannoWTH3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoWTH3$annotation)

WTH3_genes<-unique(dataf_peakannoWTH3$SYMBOL)
WTH3_peaks<-unique(dataf_peakannoWTH3$peaksite)


dataf_peakannoKH3<-as.data.frame(peakAnnoKH3)
dataf_peakannoKH3$peaksite<-paste(dataf_peakannoKH3$SYMBOL,dataf_peakannoKH3$annotation,sep='_')
colnames(dataf_peakannoKH3)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakannoKH3<- dataf_peakannoKH3[rev(order(dataf_peakannoKH3$FE)),]
dataf_peakannoKH3$annotation<-gsub('Intron.*','Intron',dataf_peakannoKH3$annotation)
dataf_peakannoKH3$annotation<-gsub('Exon.*','Exon',dataf_peakannoKH3$annotation)
dataf_peakannoKH3$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoKH3$annotation)
dataf_peakannoKH3$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoKH3$annotation)

KH3_genes<-unique(dataf_peakannoKH3$SYMBOL)
KH3_peaks<-unique(dataf_peakannoKH3$peaksite)



##DATAFRAME WITH genes id linked to peaks
dataf_peakannoFWT<-as.data.frame(peakAnnoFWT)
dataf_peakannoFWT$peaksite<-paste(dataf_peakannoFWT$SYMBOL,dataf_peakannoFWT$annotation,sep='_')
colnames(dataf_peakannoFWT)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakannoFWT<- dataf_peakannoFWT[rev(order(dataf_peakannoFWT$FE)),]
dataf_peakannoFWT$annotation<-gsub('Intron.*','Intron',dataf_peakannoFWT$annotation)
dataf_peakannoFWT$annotation<-gsub('Exon.*','Exon',dataf_peakannoFWT$annotation)
dataf_peakannoFWT$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFWT$annotation)
dataf_peakannoFWT$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFWT$annotation)

FWT_genes<-unique(dataf_peakannoFWT$SYMBOL)
FWT_peaks<-unique(dataf_peakannoFWT$peaksite)


dataf_peakannoFKO<-as.data.frame(peakAnnoFKO)
dataf_peakannoFKO$peaksite<-paste(dataf_peakannoFKO$SYMBOL,dataf_peakannoFKO$annotation,sep='_')
colnames(dataf_peakannoFKO)[c(6,7,9,10,11)]<-c("Name","Score","FE","log10pval","log10qval")
dataf_peakannoFKO<- dataf_peakannoFKO[rev(order(dataf_peakannoFKO$FE)),]
dataf_peakannoFKO$annotation<-gsub('Intron.*','Intron',dataf_peakannoFKO$annotation)
dataf_peakannoFKO$annotation<-gsub('Exon.*','Exon',dataf_peakannoFKO$annotation)
dataf_peakannoFKO$annotation<-gsub('Downstream.*','Downstream',dataf_peakannoFKO$annotation)
dataf_peakannoFKO$annotation<-gsub('Promoter.*','Promoter',dataf_peakannoFKO$annotation)

FKO_genes<-unique(dataf_peakannoFKO$SYMBOL)
FKO_peaks<-unique(dataf_peakannoFKO$peaksite)

length(FKO_genes)
length(KH3_genes)
length(WTH3_genes)
length(KH3_genes)


distWtH3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoWTH3$annotation)))))
colnames(distWtH3) <- as.character(unlist(distWtH3[1,]))
distWtH3<-distWtH3[-1,]
distWtH3$Sample<-"WTH3"

distFWT<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFWT$annotation)))))
colnames(distFWT) <- as.character(unlist(distFWT[1,]))
distFWT<-distFWT[-1,]
distFWT$Sample<-"FWT"


distKH3<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoKH3$annotation)))))
colnames(distKH3) <- as.character(unlist(distKH3[1,]))
distKH3<-distKH3[-1,]
distKH3$Sample<-"KH3"


distFKO<-as.data.frame(t(as.data.frame(prop.table(table(dataf_peakannoFKO$annotation)))))
colnames(distFKO) <- as.character(unlist(distFKO[1,]))
distFKO<-distFKO[-1,]
distFKO$Sample<-"FKO"

distWTKH3<-rbind(distKH3,distWtH3,distFWT,distFKO)
colnames(distWTKH3)[1]<-c("Three_UTR")
colnames(distWTKH3)[2]<-c("Five_UTR")

distall_WTKH3 <- data.frame(sapply(distWTKH3[,1:7], function(x) as.numeric(as.character(x))*100))
distall_WTKH3<-cbind(distWTKH3$Sample,distall_WTKH3)
colnames(distall_WTKH3)[1]<-"Condition"



distall_melt<-melt(distall_WTKH3,id.vars = c("Condition"))

ggplot(distall_melt,aes(x=Condition,y=value,fill=variable))+
  geom_bar(stat="identity") +
  facet_wrap(~Condition,scale="free",ncol=4)+
  scale_fill_brewer(palette="Accent")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  xlab("")+ylab("Fraction of peaks")








#genes
comgenes<-as.data.frame(table(sort(c(WTH3_genes,KH3_genes))))

#peaks
comgenes<-as.data.frame(table(sort(c(WTH3_peaks,KH3_peaks))))

comgenes<-comgenes[order(comgenes$Freq,decreasing = TRUE),]
comgenes_all<-subset(comgenes,comgenes$Freq==2)
uniqueenes_all<-subset(comgenes,comgenes$Freq==1)

#genes
uniquegenes_KH3<-setdiff(KH3_genes,comgenes_all$Var1)
uniquegenes_WTH3<-setdiff(WTH3_genes,comgenes_all$Var1)

#peaks
uniquepeaks_KH3<-setdiff(KH3_peaks,comgenes_all$Var1)
uniquepeaks_WTH3<-setdiff(WTH3_peaks,comgenes_all$Var1)

#genes
uniquegenes_KH3_sites<-dataf_peakannoKH3[with(dataf_peakannoKH3, dataf_peakannoKH3$SYMBOL %in% uniquegenes_KH3),]
uniquegenes_WTH3_sites<-dataf_peakannoWTH3[with(dataf_peakannoWTH3, dataf_peakannoWTH3$SYMBOL %in% uniquegenes_WTH3),]


uniquegenes_KH3_sites$type=c("KH3")
uniquegenes_KH3_sites_sub<-subset(uniquegenes_KH3_sites,select=c(start,seqnames,width,type,X90,SYMBOL))
colnames(uniquegenes_KH3_sites_sub)<-c("start","seqnames","width","type","cov","Gene")
uniquegenes_WTH3_sites$type=c("WTH3")
uniquegenes_WTH3_sites_sub<-subset(uniquegenes_WTH3_sites,select=c(start,seqnames,width,type,X62,SYMBOL))
colnames(uniquegenes_WTH3_sites_sub)<-c("start","seqnames","width","type","cov","Gene")

uniquegenes_KH3_WTH3<-rbind(uniquegenes_KH3_sites_sub,uniquegenes_WTH3_sites_sub)
uniquegenes_KH3_WTH3$seqnames<-gsub('chr','',uniquegenes_KH3_WTH3$seqnames)
uniquegenes_KH3_WTH3$seqnames<-as.factor(uniquegenes_KH3_WTH3$seqnames)

uniquegenes_KH3_WTH3$seqnames <- ordered(uniquegenes_KH3_WTH3$seqnames,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X"))

#merge with ISC signature
scsignature$source<-"SC"

#merge with fetal signature
fetalsignature$source<-"Fetal"

#uniquegenes_KH3_WTH3<-merge(scsignature,uniquegenes_KH3_WTH3,all.y="TRUE")
uniquegenes_KH3_WTH3<-merge(fetalsignature,uniquegenes_KH3_WTH3,all.y="TRUE")

uniquegenes_KH3_WTH3$source[is.na(uniquegenes_KH3_WTH3$source)] <- "noSC"
uniquegenes_KH3_WTH3$dir[is.na(uniquegenes_KH3_WTH3$dir)] <- "NULL"

uniquegenes_KH3_WTH3$class <- paste(uniquegenes_KH3_WTH3$dir,uniquegenes_KH3_WTH3$type,sep="_")

uniquegenes_KH3_WTH3_noX<-subset(uniquegenes_KH3_WTH3,uniquegenes_KH3_WTH3$seqnames!="X")

ggplot(uniquegenes_KH3_WTH3_noX,aes(x=start,y=cov,label=Gene))+
  geom_point(aes(color=type))+
  geom_linerange(aes(x=start, ymax=cov, ymin=0,color=type),position = position_jitter(height = 0L, seed = 1L))+
  #geom_text(aes(start, cov, label = Gene), data = subset(uniquegenes_KH3_WTH3_noX,uniquegenes_KH3_WTH3_noX$source=="SC"),check_overlap = TRUE,vjust=-0.5)+
  geom_label_repel(aes(fill = class), data = subset(uniquegenes_KH3_WTH3_noX,uniquegenes_KH3_WTH3_noX$source!="noSC"),colour = "black", fontface = "italic",size=2.5)+
  #scale_color_manual(values=c("coral3","cornflowerblue","firebrick1","cyan"))+
  scale_color_manual(values=c("coral3","cornflowerblue"))+
  scale_fill_manual(values=c("darksalmon","firebrick1","darkslateblue"))+
  facet_wrap(~seqnames,scales="free",ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_blank())

  

#peaks
uniquepeaks_KH3_sites<-dataf_peakannoKH3[with(dataf_peakannoKH3, dataf_peakannoKH3$peaksite %in% uniquepeaks_KH3),]
uniquepeaks_WTH3_sites<-dataf_peakannoWTH3[with(dataf_peakannoWTH3, dataf_peakannoWTH3$peaksite %in% uniquepeaks_WTH3),]


## Create excel with table KH3 genes, WTH3 genes and ISC genes
write.table(uniquegenes_KH3_sites,"KH3_genes_sites.csv",quote = FALSE, row.names = FALSE,sep='\t')
write.table(uniquegenes_WTH3_sites,"WTH3_genes_sites.csv",quote = FALSE, row.names = FALSE,sep='\t')

## Create excel with table KH3 peaks, WTH3 peaks 
write.table(uniquepeaks_KH3_sites,"KH3_uniquepeaks_sites.csv",quote = FALSE, row.names = FALSE,sep='\t')
write.table(uniquepeaks_WTH3_sites,"WTH3_uniquepeaks_sites.csv",quote = FALSE, row.names = FALSE,sep='\t')



#venn diagram common peaks or genes KO and WT
### venn diagram ##
geneLists <- list(KH3 = KH3_genes, FKO = FKO_genes, FWT = FWT_genes, WTH3 = WTH3_genes)
#geneLists <- list(WH3 = WTH3_genes,FWT = FWT_genes)
#geneLists <- list(KH3 = KH3_genes, WTH3 = WTH3_genes)
geneLists <- lapply(geneLists, function(x) x[!is.na(x)])
head(geneLists)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneLists)

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
VENN.LIST <- geneLists
dev.off()
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkolivegreen", "coral","green","orange"), alpha=c(0.3,0.3,0.3,0.3), cex = 2, cat.fontface=4, category.names=c("KH3","FKO","FWT","WH3"), main="Genes")

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

common_KO_WT<-inters[[3]]

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
