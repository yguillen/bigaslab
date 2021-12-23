
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

summary(cormul$CHC1)
summary(cormul$CHC2)
summary(cormul$CHV1)
summary(cormul$CHV2)

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


# Distribution of intensities
cormul_density<-cormul_filt[,4:7]
cormul_density<-melt(cormul_density)

p<-ggplot(data=cormul_density[cormul_density$value<20,],aes(value))+
  geom_density(aes(color=variable))+
  theme_bw()
pdf("density_exp.pdf")
print(p)
dev.off()


## predictive model for new peaks in CHC

#model Crypt CH, taking values within normal distribution, maximum value 10
boxplot(cormul_filt[cormul_filt$CHC2<10,]$CHC2)

p<-ggplot(data=cormul_filt,aes(x=CHC1,y=CHC2))+
  geom_point(data=cormul_filt[sample(nrow(cormul_filt),5000),],color="black",size=1,alpha=0.3)+
  geom_smooth(data=cormul_filt,method = "lm",se=TRUE)+
  theme_bw()+
  coord_fixed(ratio=1)
pdf("CH1_CHC2_random5000.pdf")
print(p)
dev.off()

model_CHC<-lm(CHC1~CHC2,data=cormul_filt[cormul_filt$CHC1<10 & cormul_filt$CHC2<10,])
model_CHC


## predictive model for new peaks in CHV

p<-ggplot(data=cormul_filt,aes(x=CHV1,y=CHV2))+
  geom_point(data=cormul_filt[sample(nrow(cormul_filt),5000),],color="black",size=1,alpha=0.3)+
#  geom_smooth(data=cormul_filt,method = "lm",se=TRUE)+
  theme_bw()+
  coord_fixed(ratio=1)
pdf("CHV1_CHV2_random5000.pdf")
print(p)
dev.off()

#model Villi CH, taking values within normal distribution, maximum value 10
boxplot(cormul_filt[cormul_filt$CHV2<10,]$CHV2)

model_CHV<-lm(CHV1~CHV2,data=cormul_filt[cormul_filt$CHV1<10 & cormul_filt$CHV2<10,])
model_CHV


#observed CHV1 in model CHC
obs_CHV1<-data.frame(
  CHV1=c(cormul_filt[,colnames(cormul_filt) == "CHV1"]))

row.names(obs_CHV1)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
colnames(obs_CHV1)<-"CHC2"
obs_CHV1$CHC2=obs_CHV1$CHC2

#observed CHV2 in model CHCC
obs_CHV2<-data.frame(
  CHV2=c(cormul_filt[,colnames(cormul_filt) == "CHV2"]))

row.names(obs_CHV2)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
colnames(obs_CHV2)<-"CHC2"
obs_CHV2$CHC2=obs_CHV2$CHC2/2

# prediction for CHV1 from model CHC
pre_CHV1<-data.frame(predict(model_CHC,
                             level=0.99,
                             newdata = obs_CHV1,
                             interval="prediction"))

colnames(pre_CHV1)<-c("fit_CHV1_from_CHCmodel",
                      "lwr_CHV1_from_CHCmodel",
                      "upr_CHV1_from_CHCmodel")

# prediction for CHV2 from model CHC
pre_CHV2<-data.frame(predict(model_CHC,
                             level=0.99,
                             newdata = obs_CHV2,
                             interval="prediction"))

colnames(pre_CHV2)<-c("fit_CHV2_from_CHCmodel",
                      "lwr_CHV2_from_CHCmodel",
                      "upr_CHV2_from_CHCmodel")



## Merge predictive CHV from model using observed H4V1 with observed values H4V2
cormul_filt_model<-cbind(cormul_filt,pre_CHV1,pre_CHV2)


# Genes that do no fit the model CHV1 from H4C1
cormul_filt_model$model_CHC_CHV1<-ifelse((cormul_filt_model$lwr_CHV1_from_CHCmodel < cormul_filt_model$CHV1 & 
                                            cormul_filt_model$CHV1 < cormul_filt_model$upr_CHV1_from_CHCmodel), 
                                         c("YES"),c("NO")) 

cormul_filt_model$model_CHC_CHV2<-ifelse((cormul_filt_model$lwr_CHV2_from_CHCmodel < cormul_filt_model$CHV2 & 
                                            cormul_filt_model$CHV2 < cormul_filt_model$upr_CHV2_from_CHCmodel), 
                                         c("YES"),c("NO")) 


cormul_filt_model$model_CHC_CHV<-paste(cormul_filt_model$model_CHC_CHV1,cormul_filt_model$model_CHC_CHV2,sep="_")


table(cormul_filt_model$model_CHC_CHV)

cormul_filt_model$ratio_1<-cormul_filt_model$CHC1/(cormul_filt_model$CHV1)
cormul_filt_model$ratio_2<-cormul_filt_model$CHC1/(cormul_filt_model$CHV2/2)


# model CHC 

p<-ggplot(data=cormul_filt_model,aes(x=CHC1,y=CHV2))+
  geom_point(data=cormul_filt_model[cormul_filt_model$model_CHC_CHV=="YES_YES",][sample(nrow(cormul_filt_model[cormul_filt_model$model_CHC_CHV=="YES_YES",]),500),],color="black",size=1,alpha=0.3)+
  geom_point(data=cormul_filt_model[cormul_filt_model$model_CHC_CHV!="YES_YES",],aes(color=model_CHC_CHV),size=1,alpha=0.3)+
  scale_colour_brewer(palette = "Set1")+
  #  geom_smooth(data=cormul_filt_model[cormul_filt_model$model_CHC=="YES",][sample(nrow(cormul_filt_model[cormul_filt_model$model_CHC=="YES",]),500),],method = "lm",se=TRUE)+
  theme_bw()+
  coord_fixed()
p

pdf("model_CR_VI_cor.pdf")
print(p)
dev.off()

cor.test(cormul_filt_model$CHC1,cormul_filt_model$CHC2)
cor.test(cormul_filt_model$CHC2,cormul_filt_model$CHV1)
cor.test(cormul_filt_model$CHC1,cormul_filt_model$CHC2)
cor.test(cormul_filt_model$CHV1,cormul_filt_model$CHV2)


# Write data.table with all regions
write.table(cormul_filt_model[,c(1:3,24)],"/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/allregions_crypt_villi_CH.bed",quote = FALSE,row.names = FALSE,sep="\t")
#write.table(cormul_filt_model[,c(1:3,24)],"/Users/yolanda_guillen/Desktop/IMIM/CR_VI_H4/allregions_crypt_villi_CH.bed",quote = FALSE,row.names = FALSE,sep="\t")


# Write data table with peaks discordance
write.table(cormul_filt_model[cormul_filt_model$model_CHC_CHV=="NO_NO",c(1:3,24)],"/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/regions_doble_no_model_crypt_villi_CH.bed",quote = FALSE,row.names = FALSE,sep="\t")
write.table(cormul_filt_model[cormul_filt_model$model_CHC_CHV=="NO_YES" | cormul_filt_model$model_CHC_CHV=="YES_NO",c(1:3,24)],"/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/regions_single_no_model_crypt_villi_CH.bed",quote = FALSE,row.names = FALSE,sep="\t")

write.table(cormul_filt_model[cormul_filt_model$model_CHC_CHV=="NO_NO",c(1:3,24)],"/Users/yolanda_guillen/Desktop/IMIM/CR_VI_H4/regions_doble_no_model_crypt_villi_CH.bed",quote = FALSE,row.names = FALSE,sep="\t")
write.table(cormul_filt_model[cormul_filt_model$model_CHC_CHV=="NO_YES" | cormul_filt_model$model_CHC_CHV=="YES_NO",c(1:3,24)],"/Users/yolanda_guillen/Desktop/IMIM/CR_VI_H4/regions_single_no_model_crypt_villi_CH.bed",quote = FALSE,row.names = FALSE,sep="\t")


## Annotation of regions with annotatr
library(annotatr)

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
CHC_all<-read_regions(con="/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/allregions_crypt_villi_CH.bed", genome = 'mm10', format = 'bed')
#CHC_all<-read_regions(con="/Users/yolanda_guillen/Desktop/IMIM/CR_VI_H4/allregions_crypt_villi_CH.bed", genome = 'mm10', format = 'bed')


all_annotated = annotate_regions(
  regions = CHC_all,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## Regions not fitting the model SINGLE
CHC_nomodel_single<-read_regions(con = "/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/regions_single_no_model_crypt_villi_CH.bed", genome = 'mm10', format = 'bed')
#CHC_nomodel_single<-read_regions(con = "/Users/yolanda_guillen/Desktop/IMIM/CR_VI_H4/regions_single_no_model_crypt_villi_CH.bed", genome = 'mm10', format = 'bed')


CHC_nomodel_single_annotated = annotate_regions(
  regions = CHC_nomodel_single,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## Regions not fitting the model DOBLE
CHC_nomodel_doble<-read_regions(con = "/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/regions_doble_no_model_crypt_villi_CH.bed", genome = 'mm10', format = 'bed')
#CHC_nomodel_doble<-read_regions(con = "/Users/yolanda_guillen/Desktop/IMIM/CR_VI_H4/regions_doble_no_model_crypt_villi_CH.bed", genome = 'mm10', format = 'bed')


CHC_nomodel_doble_annotated = annotate_regions(
  regions = CHC_nomodel_doble,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

# Compare with random regions
# create random regions
library(regioneR)
library(BSgenome.Mmusculus.UCSC.mm10)

random_regions<-createRandomRegions(nregions=length(c(CHC_nomodel_doble$name,CHC_nomodel_single$name)),
                                    length.mean=1000,
                                    length.sd=20,
                                    genome="mm10", mask=NULL, non.overlapping=TRUE)

random_annotated = annotate_regions(
  regions = random_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


CHC_annotated_wrandom = plot_annotation(
  annotated_regions = CHC_nomodel_doble_annotated,
  annotated_random = random_annotated,
  annotation_order = annots_order,
  plot_title = 'CHC_annotated vs random',
  x_label = 'Annotations',
  y_label = 'Count')

comp_ran<-print(CHC_annotated_wrandom)
pdf("comparative_doble.pdf")
print(comp_ran)
dev.off()

# Regions in pairs of annotations
CHC_coannotations = plot_coannotations(
  annotated_regions = CHC_nomodel_doble,
  annotation_order = annots_order,
  axes_label = 'Annotations',
  plot_title = 'Regions in Pairs of Annotations')

print(CHC_coannotations)


# proportions

all_annotated<-data.frame(all_annotated)
propall_annot<-data.frame(prop.table(table(all_annotated$annot.type)))
propall_annot$sam<-"all_CHC"

CHC_nomodel_single_annotated<-data.frame(CHC_nomodel_single_annotated)
propCHC_nomodel_single<-data.frame(prop.table(table(CHC_nomodel_single_annotated$annot.type)))
propCHC_nomodel_single$sam<-"CHC_nomodel_single"

CHC_nomodel_doble_annotated<-data.frame(CHC_nomodel_doble_annotated)
propCHC_nomodel_doble<-data.frame(prop.table(table(CHC_nomodel_doble_annotated$annot.type)))
propCHC_nomodel_doble$sam<-"CHC_nomodel_doble"

random_annotated<-data.frame(random_annotated)
proprandom<-data.frame(prop.table(table(random_annotated$annot.type)))
proprandom$sam<-"random"


propall<-rbind(propall_annot,propCHC_nomodel_single)
propall<-rbind(propall,propCHC_nomodel_doble)
propall<-rbind(propall,proprandom)

propall$sam<-factor(propall$sam,levels = c("random","all_CHC","CHC_nomodel_single","CHC_nomodel_doble"))

p<-ggplot(propall,aes(x=sam,y=Freq,fill=Var1))+
  geom_bar(stat="identity",color="black") +
  scale_fill_brewer(palette="Paired")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Frequency annotation")
p
pdf("Frequency_annotations.pdf")
print(p)
dev.off()

# merge with scores models

cormul_chr<-subset(cormul_filt_model,select=c("chr","start","CHV1","CHV2","CHC1","CHC2","model_CHC_CHV"))
cormul_chr$name<-row.names(cormul_chr)
cormul_chr<-melt(cormul_chr,id.vars=c("chr","start","model_CHC_CHV","name"))

#cormul_chr$chr<-gsub('chr','',cormul_chr$chr)
#cormul_chr$chr<-as.factor(cormul_chr$chr)

#cormul_chr$chr<-factor(cormul_chr$chr,levels = c("1","2","3","4","5","6","7","8","9","10","11","12",
     #                                            "13","14","15","16","17","18","19","X","Y","M"))


all_annotated$start<-all_annotated$start-1
all_annotated$name<-paste(all_annotated$seqnames,all_annotated$start,sep="_")

all_annotated_model<-merge(cormul_chr,all_annotated[,c(6,15,16)],by="name")
all_annotated_model<-all_annotated_model[!duplicated(all_annotated_model$name),]

all_annotated_model$annot.type<-gsub('mm10_genes_','',all_annotated_model$annot.type)
all_annotated_model$annot.type<-factor(all_annotated_model$annot.type,levels=c("promoters","5UTRs","exons","introns","3UTRs","1to5kb"))


test<-all_annotated_model[all_annotated_model$model_CHC_CHV=="NO_NO",]
unique(test$annot.symbol)

test<-merge(test,cormul_filt_model,by.x="name",by.y="row.names")
test<-test[order(test$ratio_1),]

#remove discordant CHC replicates
test$rat_CHV<-test$CHV1/(test$CHV2/2)

test$rat_CHC<-test$CHC1/test$CHC2

boxplot(test$rat_CHC)
test<-test[test$rat_CHV<summary(test$rat_CHV)[5],]

plot(test$ratio_1,test$ratio_2)
boxplot(test$rat_CHC)

test_filt<-test[test$rat_CHC<2,]
boxplot(test_filt$rat_CHC)
test_filt<-test_filt[test_filt$rat_CHC>0.5,]
sort(unique(test_filt$annot.symbol))

test_filt<-test_filt[test_filt$ratio_1<0.23 | test_filt$ratio_1>1.25,]
sort(unique(test_filt$annot.symbol))

write.table(test_filt,"/Volumes/grcmc/YGUILLEN/CR_VI_H4/cor_tests/candidate_genes_023_125.tsv",sep="\t",quote=FALSE)

ggplot(data=test[test$dif_1>13,],aes(x=CHC1,y=CHV1,label=annot.symbol))+
  geom_point(data=test[test$model_CHC_CHV.x=="YES_YES"& test$dif_1>13,],color="black",size=1,alpha=0.3)+
  geom_point(data=test[test$model_CHC_CHV.x!="YES_YES"& test$dif_1>13,],aes(color=model_CHC_CHV.x),size=1,alpha=0.3)+
  geom_label_repel(data=test[(test$CHC1>20 | test$CHV1>30),],aes(label=annot.symbol),size=3)+
  scale_colour_brewer(palette = "Set1")+
  #  geom_smooth(data=cormul_filt_model[cormul_filt_model$model_CHC=="YES",][sample(nrow(cormul_filt_model[cormul_filt_model$model_CHC=="YES",]),500),],method = "lm",se=TRUE)+
  theme_bw()+
  coord_fixed()


write.table(test_filt,'/Users/yolanda_guillen/Desktop/IMIM/CR_VI_H4/candidates_dif_1_20.tsv',quote = FALSE,sep="\t")
write.table(test_filt,'/Users/yguillen/Downloads/candidates_dif_1_15.tsv',quote = FALSE,sep="\t")


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

