
#BiocManager::install("WGCNA")
BiocManager::install("GGally")
BiocManager::install("randomcoloR")
BiocManager::install("annotatr")
BiocManager::install("regioneR")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

library(dplyr)
library(reshape2)
library(ggnewscale)
library(ggplot2)
library(gridExtra)
#library(WGCNA)
#library(GGally)
library(randomcoloR)
library(annotatr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggnewscale)

### Scores files in temp (teleworking) are in: temp/CKC_Irene/peakcall/

### Multibigwig correlation ###
#cormul<-read.delim("/Users/yolanda_guillen/Desktop/IMIM/CR_VI_H4/peakcall/scores_cond_compare.tab")
cormul<-read.delim("/Volumes/cancer/CR_VI_H4/peakcall/norm_bigwig/scores_cond_compare.tab")
colnames(cormul)<-gsub('X.','',colnames(cormul))
colnames(cormul)<-gsub('\\.','',colnames(cormul))

#Remove all zeros
remwin<-row.names(cormul[cormul$CHC1==0 & cormul$CHC2==0 &
                           cormul$CHV1==0 & cormul$CHV2==0 &
                           cormul$H4C1==0 & cormul$H4C2==0 &
                           cormul$H4V1==0 & cormul$H4V2==0 &
                           cormul$INCR<100 & cormul$INVI<100,])
cormul<-cormul[!row.names(cormul) %in% remwin,]

summary(cormul$CHC2)
summary(cormul$H4C1)
summary(cormul$H4C2)
summary(cormul$INVI)

## filter outliers that are input peaks likely too
cormul_filt<-cormul[cormul$INVI<50,]
#Asign color value to chromosome
# Create random colors
library(randomcoloR)
n <- dim(table(cormul_filt$chr))
palette <- distinctColorPalette(n)

colchr<-data.frame(chr=table(cormul_filt$chr),
                   colored=palette)

cormul_filt<-merge(cormul_filt,colchr,by.x="chr",by.y="chr.Var1")


# Plot correlations

library(plotly)
library(psych)
library(corrplot)
library(Hmisc)
library(ggplot2)
library(GGally)

#Select random windows
ranchr<-cormul_filt[,c(4,6,8,9,1)][sample(nrow(cormul_filt[,c(4,6,8,9,1)]), 15000), ]

# Or select one or multiple chromosome
chrcormul<-cormul_filt[cormul_filt$chr=="chr11" | cormul_filt$chr=="chr18",c(4,6,8,9,1)]

# Correlation matrix corrected for multiple values
corrplot(cor(cormul_filt[,c(4:11)]))


## predictive model for new peaks in H4 or CH

#model Crypt CH
model_CHC<-lm(CHC1~CHC2,data=cormul_filt)
model_CHC

#model Crypt H4
model_H4C<-lm(H4C1~H4C2,data=cormul_filt)
model_H4C

# model Villi CH
model_CHV<-lm(CHV1~CHV2,data=cormul_filt)
model_CHV

# model Villi H4
model_H4V<-lm(H4V1~H4V2,data=cormul_filt)
model_H4V

#observed H4C1 in model CHC
obs_H4C1<-data.frame(
  H4C1=c(cormul_filt[,colnames(cormul_filt) == "H4C1"]))

row.names(obs_H4C1)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
colnames(obs_H4C1)<-"CHC2"

#observed H4C2 in model CHC
obs_H4C2<-data.frame(
  H4C2=c(cormul_filt[,colnames(cormul_filt) == "H4C2"]))

row.names(obs_H4C2)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
colnames(obs_H4C2)<-"CHC2"


#observed CHC1 in model H4C
obs_CHC1<-data.frame(
  CHC1=c(cormul_filt[,colnames(cormul_filt) == "CHC1"]))

row.names(obs_CHC1)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
colnames(obs_CHC1)<-"H4C2"

#observed CHC2 in model H4C
obs_CHC2<-data.frame(
  CHC2=c(cormul_filt[,colnames(cormul_filt) == "CHC2"]))

row.names(obs_CHC2)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
colnames(obs_CHC2)<-"H4C2"

# prediction for H4C1 from model CHC
pre_H4C1<-data.frame(predict(model_CHC,
                             level=0.99,
                   newdata = obs_H4C1,
                   interval="prediction"))

colnames(pre_H4C1)<-c("fit_H4C1_from_CHCmodel",
                     "lwr_H4C1_from_CHCmodel",
                     "upr_H4C1_from_CHCmodel")

# prediction for H4C2 from model CHC
pre_H4C2<-data.frame(predict(model_CHC,
                             level=0.99,
                             newdata = obs_H4C2,
                             interval="prediction"))

colnames(pre_H4C2)<-c("fit_H4C2_from_CHCmodel",
                      "lwr_H4C2_from_CHCmodel",
                      "upr_H4C2_from_CHCmodel")


# prediction for CHC1 from model H4C
pre_CHC1<-data.frame(predict(model_H4C,
                             level=0.99,
                             newdata = obs_CHC1,
                             interval="prediction"))

colnames(pre_CHC1)<-c("fit_CHC1_from_H4Cmodel",
                      "lwr_CHC1_from_H4Cmodel",
                      "upr_CHC1_from_H4Cmodel")

# prediction for CHC2 from model H4C
pre_CHC2<-data.frame(predict(model_H4C,
                             level=0.99,
                             newdata = obs_CHC2,
                             interval="prediction"))

colnames(pre_CHC2)<-c("fit_CHC2_from_H4Cmodel",
                      "lwr_CHC2_from_H4Cmodel",
                      "upr_CHC2_from_H4Cmodel")



## Merge predictive CHV from model using observed H4V1 with observed values H4V2
cormul_filt_model<-cbind(cormul_filt,pre_H4C1,pre_H4C2,pre_CHC1,pre_CHC2)


# Genes that do no fit the model CHV1 from H4C1
cormul_filt_model$model_CHC_H4C1<-ifelse((cormul_filt_model$lwr_H4C1_from_CHCmodel < cormul_filt_model$H4C1 & 
                                      cormul_filt_model$H4C1 < cormul_filt_model$upr_H4C1_from_CHCmodel), 
                             c("YES"),c("NO")) 

cormul_filt_model$model_CHC_H4C2<-ifelse((cormul_filt_model$lwr_H4C2_from_CHCmodel < cormul_filt_model$H4C2 & 
                                            cormul_filt_model$H4C2 < cormul_filt_model$upr_H4C2_from_CHCmodel), 
                                         c("YES"),c("NO")) 


cormul_filt_model$model_CHC_H4C<-paste(cormul_filt_model$model_CHC_H4C1,cormul_filt_model$model_CHC_H4C2,sep="_")

# Genes that do no fit the model H4C
cormul_filt_model$model_H4C_CHC1<-ifelse((cormul_filt_model$lwr_CHC1_from_H4Cmodel < cormul_filt_model$CHC1 & 
                                       cormul_filt_model$CHC1 < cormul_filt_model$upr_CHC1_from_H4Cmodel), 
                                    c("YES"),c("NO")) 

cormul_filt_model$model_H4C_CHC2<-ifelse((cormul_filt_model$lwr_CHC2_from_H4Cmodel < cormul_filt_model$CHC2 & 
                                         cormul_filt_model$CHC2 < cormul_filt_model$upr_CHC2_from_H4Cmodel), 
                                    c("YES"),c("NO")) 


cormul_filt_model$model_H4C_CHC<-paste(cormul_filt_model$model_H4C_CHC1,cormul_filt_model$model_H4C_CHC2,sep="_")


cormul_filt_model$mergemodel<-paste(cormul_filt_model$model_CHC_H4C,
                                    cormul_filt_model$model_H4C_CHC,sep="_")


table(cormul_filt_model$mergemodel)

cormul_filt_model$dif_1<-cormul_filt_model$CHV1-cormul_filt_model$CHC1
cormul_filt_model$dif_2<-cormul_filt_model$CHV2-cormul_filt_model$CHC1



# model CHC 

ggplot(data=cormul_filt_model,aes(x=CHC1,y=H4C1))+
  #  geom_point(color="black",size=1,alpha=0.3)+
  geom_point(data=cormul_filt_model[cormul_filt_model$model_CHC_H4C!="YES_YES",],aes(color=model_CHC_H4C),size=2,alpha=0.8)+
  scale_colour_brewer(palette = "Rainbow")+
  geom_point(data=cormul_filt_model[cormul_filt_model$model_H4C_CHC!="YES_YES",],aes(shape=model_H4C_CHC),size=1,alpha=0.8)+
#  geom_point(data=cormul_filt_model[cormul_filt_model$model_CHC=="NO",],color="coral",size=1,alpha=0.8)+
#  geom_smooth(data=cormul_filt_model[cormul_filt_model$model_CHC=="YES",][sample(nrow(cormul_filt_model[cormul_filt_model$model_CHC=="YES",]),500),],method = "lm",se=TRUE)+
  theme_bw()+
  coord_fixed(ratio=1)



cor.test(cormul_filt_model$CHC1,cormul_filt_model$CHC2)
cor.test(cormul_filt_model$CHC2,cormul_filt_model$CHV1)
cor.test(cormul_filt_model$CHC1,cormul_filt_model$CHC2)
cor.test(cormul_filt_model$CHV1,cormul_filt_model$CHV2)


# Write data.table with all regions
write.table(cormul_filt_model[,c(1:3,34)],"/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/allregions_crypt.bed",quote = FALSE,row.names = FALSE,sep="\t")

# Write data table with peaks discordance
write.table(cormul_filt_model[cormul_filt_model$mergemodel!="YES_YES_YES_YES",c(1:3,34)],"/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/regions_no_model_crypt.bed",quote = FALSE,row.names = FALSE,sep="\t")


## Annotation of regions with annotatr
library(annotatr)
vignette("annotatr-vignette")

annots = c('mm10_cpgs', 'mm10_basicgenes', 'mm10_genes_intergenic',
           'mm10_genes_intronexonboundaries')


annots = c('mm10_basicgenes','mm10_genes_intergenic')

annotations = build_annotations(genome = 'mm10', annotations = annots)

annots_order = c(
  'mm10_genes_intergenic',
  'mm10_cpg_inter',
  'mm10_cpg_islands',
  'mm10_cpg_shelves',
  'mm10_cpg_shores',
  'mm10_genes_1to5kb',
  'mm10_genes_promoters',
  'mm10_genes_5UTRs',
  'mm10_genes_exons',
  'mm10_genes_introns',
  'mm10_genes_intronexonboundaries',
  'mm10_genes_3UTRs')

# Annotate all regions
CHC_all<-read_regions(con="/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/allregions_crypt.bed", genome = 'mm10', format = 'bed')

all_annotated = annotate_regions(
  regions = CHC_all,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## Regions not fitting the model
CHC_nomodel<-read_regions(con = "/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/regions_no_model_crypt.bed", genome = 'mm10', format = 'bed')

CHC_annotated = annotate_regions(
  regions = CHC_nomodel,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)



# Compare with random regions
# create random regions
library(regioneR)
library(BSgenome.Mmusculus.UCSC.mm10)

random_regions<-createRandomRegions(nregions=length(CHC_annotated$name),
                    length.mean=1000,
                    length.sd=20,
                    genome="mm10", mask=NULL, non.overlapping=TRUE)

random_annotated = annotate_regions(
  regions = random_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


CHC_annotated_wrandom = plot_annotation(
  annotated_regions = CHC_annotated,
  annotated_random = random_annotated,
  annotation_order = annots_order,
  plot_title = 'CHC_annotated vs random',
  x_label = 'Annotations',
  y_label = 'Count')

comp_ran<-print(CHC_annotated_wrandom)




# Regions in pairs of annotations
CHC_coannotations = plot_coannotations(
  annotated_regions = CHC_annotated,
  annotation_order = annots_order,
  axes_label = 'Annotations',
  plot_title = 'Regions in Pairs of Annotations')

print(CHC_coannotations)


# proportions

all_annotated<-data.frame(all_annotated)
propall_annot<-data.frame(prop.table(table(all_annotated$annot.type)))
propall_annot$sam<-"all_CHC"

CHC_annotated<-data.frame(CHC_annotated)
propCHC<-data.frame(prop.table(table(CHC_annotated$annot.type)))
propCHC$sam<-"CHC_nomodel"

random_annotated<-data.frame(random_annotated)
proprandom<-data.frame(prop.table(table(random_annotated$annot.type)))
proprandom$sam<-"random"


propall<-rbind(propall_annot,propCHC)
propall<-rbind(propall,proprandom)

propall$sam<-factor(propall$sam,levels = c("random","all_CHC","CHC_nomodel"))

ggplot(propall,aes(x=sam,y=Freq,fill=Var1))+
  geom_bar(stat="identity",color="black") +
  scale_fill_brewer(palette="Paired")+
  
  
      
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Frequency annotation")


# merge with scores models

cormul_chr<-subset(cormul_filt_model,select=c("chr","start","CHV1","CHV2","CHC1","CHC2","mergemodel","Row.names"))
cormul_chr<-melt(cormul_chr,id.vars=c("Row.names","chr","start","mergemodel"))

cormul_chr$chr<-gsub('chr','',cormul_chr$chr)
cormul_chr$chr<-as.factor(cormul_chr$chr)

cormul_chr$chr<-factor(cormul_chr$chr,levels = c("1","2","3","4","5","6","7","8","9","10","11","12",
                                                 "13","14","15","16","17","18","19","X","Y","M"))
  
  
all_annotated$start<-all_annotated$start-1
all_annotated$name<-paste(all_annotated$seqnames,all_annotated$start,sep="_")

all_annotated_model<-merge(cormul_chr,all_annotated[,c(6,15,16)],by.x="Row.names",by.y="name")
all_annotated_model<-all_annotated_model[!duplicated(all_annotated_model),]

all_annotated_model$annot.type<-gsub('mm10_genes_','',all_annotated_model$annot.type)
all_annotated_model$annot.type<-factor(all_annotated_model$annot.type,levels=c("promoters","5UTRs","exons","introns","3UTRs","1to5kb"))
all_annotated_model$Row.names<-NULL



######## Annotation of intergenic regions with chipseeker
# mice genes database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


library(ChIPseeker)

# ALL REGIONS BRD4
Brd4_all_chiptab<-annotatePeak(Brd4_all, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
Brd4_all_chip<-as.data.frame(Brd4_all_chiptab@anno)

Brd4_all_chip$annotation<-gsub('\\(.*,.intron.*','',Brd4_all_chip$annotation)
Brd4_all_chip$annotation<-gsub('\\(.*, exon.*','',Brd4_all_chip$annotation)

Brd4_all_prop<-as.data.frame(prop.table(table(Brd4_all_chip$annotation)))
Brd4_all_prop$sample<-"Brd4_all"

# WTIR and KOIR regions
WTIR_KOIR_chiptab<-annotatePeak(WTIR_KOIR_regions, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
WTIR_KOIR_chip<-as.data.frame(WTIR_KOIR_chiptab@anno)

WTIR_KOIR_chip$annotation<-gsub('\\(.*,.intron.*','',WTIR_KOIR_chip$annotation)
WTIR_KOIR_chip$annotation<-gsub('\\(.*, exon.*','',WTIR_KOIR_chip$annotation)

WTIR_KOIR_prop<-as.data.frame(prop.table(table(WTIR_KOIR_chip$annotation)))
WTIR_KOIR_prop$sample<-"WTIR_KOIR"

# CWT1 regions
CWT1_chiptab<-annotatePeak(CWT1_regions, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
CWT1_chip<-as.data.frame(CWT1_chiptab@anno)

CWT1_chip$annotation<-gsub('\\(.*,.intron.*','',CWT1_chip$annotation)
CWT1_chip$annotation<-gsub('\\(.*, exon.*','',CWT1_chip$annotation)

CWT1_prop<-as.data.frame(prop.table(table(CWT1_chip$annotation)))
CWT1_prop$sample<-"CWT1"

# CKO1 regions
CKO1_chiptab<-annotatePeak(CKO1_regions, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
CKO1_chip<-as.data.frame(CKO1_chiptab@anno)

CKO1_chip$annotation<-gsub('\\(.*,.intron.*','',CKO1_chip$annotation)
CKO1_chip$annotation<-gsub('\\(.*, exon.*','',CKO1_chip$annotation)

CKO1_prop<-as.data.frame(prop.table(table(CKO1_chip$annotation)))
CKO1_prop$sample<-"CKO1"


# WTIR regions
WTIR_chiptab<-annotatePeak(WTIR_regions, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
WTIR_chip<-as.data.frame(WTIR_chiptab@anno)

WTIR_chip$annotation<-gsub('\\(.*,.intron.*','',WTIR_chip$annotation)
WTIR_chip$annotation<-gsub('\\(.*, exon.*','',WTIR_chip$annotation)

WTIR_prop<-as.data.frame(prop.table(table(WTIR_chip$annotation)))
WTIR_prop$sample<-"WTIR"

# KOIR regions
KOIR_chiptab<-annotatePeak(KOIR_regions, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
KOIR_chip<-as.data.frame(KOIR_chiptab@anno)

KOIR_chip$annotation<-gsub('\\(.*,.intron.*','',KOIR_chip$annotation)
KOIR_chip$annotation<-gsub('\\(.*, exon.*','',KOIR_chip$annotation)

KOIR_prop<-as.data.frame(prop.table(table(KOIR_chip$annotation)))
KOIR_prop$sample<-"KOIR"

# random regions
random_chiptab<-annotatePeak(random_regions, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
random_chip<-as.data.frame(random_chiptab@anno)

random_chip$annotation<-gsub('\\(.*,.intron.*','',random_chip$annotation)
random_chip$annotation<-gsub('\\(.*, exon.*','',random_chip$annotation)

random_prop<-as.data.frame(prop.table(table(random_chip$annotation)))
random_prop$sample<-"random"

annot_chip<-rbind(Brd4_all_prop,WTIR_KOIR_prop)
annot_chip<-rbind(annot_chip,CWT1_prop)
annot_chip<-rbind(annot_chip,KOIR_prop)
annot_chip<-rbind(annot_chip,random_prop)
annot_chip<-rbind(annot_chip,WTIR_prop)
annot_chip<-rbind(annot_chip,CKO1_prop)

annot_chip$sample<-factor(annot_chip$sample,levels = c("random","Brd4_all","CWT1","CKO1","WTIR_KOIR","WTIR","KOIR"))

ggplot(annot_chip,aes(x=sample,y=Freq,fill=Var1))+
  geom_bar(stat="identity",color="black") +
  facet_grid(~sample,scale="free",space = "free_x")+
  scale_fill_brewer(palette="Paired")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Fraction of peaks")


chipbind<-rbind(CWT1_chip,CKO1_chip,WTIR_KOIR_chip,WTIR_chip,KOIR_chip)
chipbind$Row.names<-paste(chipbind$seqnames,chipbind$start-1,sep='_')
namech<-cormul_chr[,c(1,5:6)]
chipbind<-merge(chipbind,namech,by="Row.names")

ggplot(chipbind,aes(x=start,y=value))+
  geom_point(data=chipbind[chipbind$name=="WTIR_CWT1_KOIR_noCWT1",],aes(color=variable,shape=annotation),size=1)+
  scale_shape_manual(values=c(0,1,2,2,2,3,4,5,5,5))+
  scale_color_brewer(palette = "Dark2",na.translate=F)+
  facet_wrap(~geneChr,nrow=3,scales="free_x")+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  ylab("bigwig score")+
  ggtitle("Higher in KOIR vs other")

write.table(chipbind,"/Users/yguillen/Desktop/temp/CKC_Irene/ChIPs_August_2020/genes_peaks/all_peaks/annotation_chipseeker_peaks.txt",quote = FALSE,row.names = FALSE,sep="\t")

############ Functional Enrichment

## GO terms and pathways
library(devtools)
library(enrichR)

dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "GO_Molecular_Function_2018",
         "Chromosome_Location",
         "Epigenomics_Roadmap_HM_ChIP-seq",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

enriched_Brd4 <- enrichr(unique(Brd4_all_chiptab@anno$SYMBOL), dbs)
enr_Brd4 <- enriched_Brd4[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
enr_Brd4$cat<-"Brd4"

enriched_random <- enrichr(unique(random_chiptab@anno$SYMBOL), dbs)
enr_random <- enriched_random[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
enr_random$cat<-"random"

enriched_CWT1 <- enrichr(unique(CWT1_chiptab@anno$SYMBOL), dbs)
enr_CWT1 <- enriched_CWT1[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
enr_CWT1$cat<-"CWT1"

enriched_CKO1 <- enrichr(unique(CKO1_chiptab@anno$SYMBOL), dbs)
enr_CKO1 <- enriched_CKO1[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
enr_CKO1$cat<-"CKO1"


enriched_WTIR_KOIR <- enrichr(unique(WTIR_KOIR_chiptab@anno$SYMBOL), dbs)
enr_WTIR_KOIR <- enriched_WTIR_KOIR[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
enr_WTIR_KOIR$cat<-"WTIR_KOIR"


enriched_WTIR <- enrichr(unique(WTIR_chiptab@anno$SYMBOL), dbs)
enr_WTIR <- enriched_WTIR[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
enr_WTIR$cat<-"WTIR"

enriched_KOIR <- enrichr(unique(KOIR_chiptab@anno$SYMBOL), dbs)
enr_KOIR <- enriched_KOIR[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
enr_KOIR$cat<-"KOIR"



allGO<-rbind(enr_random,enr_CWT1,enr_CKO1,enr_WTIR_KOIR,enr_WTIR,enr_KOIR)
allGO$db<-c("GO_BP")

bpsub<-subset(allGO,allGO$Adjusted.P.value<0.056)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)

bpsub$cat<-as.factor(bpsub$cat)
bpsub$cat<-factor(bpsub$cat,levels=c("CWT1","CKO1","WTIR_KOIR","WTIR","KOIR","random"))

bpsub$Genes<-gsub(';','/',bpsub$Genes)
bpsub$Term<-gsub('.GO.*','',bpsub$Term)

bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$cat, bpsub$P.value)), ]$Term))

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


# Heatmap
ggplot(bpsub, aes(x=cat, y=Term)) + 
  geom_tile(aes(fill = -(log(P.value))),colour = "black")+
  scale_fill_gradient2(low = "white",
                       high = "red") +
  #scale_fill_gradient(low = "blue", high = "red")+
  theme_dark()+
  theme(axis.text.x = element_text(size=10, hjust = 1,angle=45),
        axis.text.y = element_text(size=7, hjust = 1,face="italic"))+
  ggtitle("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

# tables before january 9th
#write.table(bpsub,"/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/Functional_enrichment/BP_enrichR_padj01.txt",sep="\t",row.names = FALSE,quote = FALSE)
#write.table(bpsub,"/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/Functional_enrichment/TF_enrichR_padj005.txt",sep="\t",row.names = FALSE,quote = FALSE)

# Heatmap with WTIR significant
bpsub_sel<-bpsub[bpsub$Term %in% c("UBTF ENCODE","SMAD4 CHEA","STAT3 CHEA","REST CHEA","SUZ12 CHEA","TRIM28 CHEA"),]$Term
ggplot(bpsub[bpsub$Term %in% bpsub_sel,], aes(x=cat, y=Term)) + 
  geom_tile(aes(fill = -(log(P.value))),colour = "black")+
  scale_fill_gradient2(low = "white",
                       high = "red") +
  #scale_fill_gradient(low = "blue", high = "red")+
  theme_dark()+
  theme(axis.text.x = element_text(size=14, hjust = 1,angle=45),
        axis.text.y = element_text(size=14, hjust = 1,face="italic"))+
  ggtitle("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")


