
#BiocManager::install("WGCNA")
#BiocManager::install("GGally")
BiocManager::install("randomcoloR")
BiocManager::install("annotatr")
BiocManager::install("regioneR")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

library(reshape2)
library(ggnewscale)
library(ggplot2)
library(gridExtra)
library(WGCNA)
library(GGally)
library(randomcoloR)
library(annotatr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

### Scores files in temp (teleworking) are in: temp/CKC_Irene/peakcall/

#compIKK<-read.delim("/Volumes/cancer/CKC/peakcall/norm_bigwig/scores_WTIR_KOIR_CKO1_compare.tab")
compIKK<-read.delim("/Users/yguillen/Desktop/temp/CKC_Irene/peakcall/scores_WTIR_KOIR_CKO1_compare.tab")

colnames(compIKK)<-gsub('X\\.','',colnames(compIKK))
colnames(compIKK)<-gsub('\\.','',colnames(compIKK))
compIKK<-melt(compIKK,id.vars = c("chr","start","end"))
compIKK$group<-"common_genes"

#exIKK<-read.delim("/Volumes/cancer/CKC/peakcall/norm_bigwig/scores_excKOIR_compare.tab")
exIKK<-read.delim("/Users/yguillen/Desktop/temp/CKC_Irene/peakcall/scores_excKOIR_compare.tab")
colnames(exIKK)<-gsub('X\\.','',colnames(exIKK))
colnames(exIKK)<-gsub('\\.','',colnames(exIKK))
exIKK<-melt(exIKK,id.vars = c("chr","start","end"))
exIKK$group<-"exc_KOIR"

exWT<-read.delim("/Volumes/cancer/CKC/peakcall/norm_bigwig/scores_excWTIR_compare.tab")
colnames(exWT)<-gsub('X\\.','',colnames(exWT))
colnames(exWT)<-gsub('\\.','',colnames(exWT))
exWT<-melt(exWT,id.vars = c("chr","start","end"))
exWT$group<-"exc_WTIR"


fcIKK<-rbind(compIKK,exIKK,exWT)

fcIKK<-fcIKK[fcIKK$value<50,]

summary(fcIKK[fcIKK$variable=="CKO1" & fcIKK$group=="common_genes",]$value)
summary(fcIKK[fcIKK$variable=="CWT1" & fcIKK$group=="common_genes",]$value)
summary(fcIKK[fcIKK$variable=="WTIR" & fcIKK$group=="common_genes",]$value)
summary(fcIKK[fcIKK$variable=="KOIR1" & fcIKK$group=="common_genes",]$value)

summary(fcIKK[fcIKK$variable=="CKO1" & fcIKK$group=="exc_KOIR",]$value)
summary(fcIKK[fcIKK$variable=="CWT1" & fcIKK$group=="exc_KOIR",]$value)
summary(fcIKK[fcIKK$variable=="WTIR" & fcIKK$group=="exc_KOIR",]$value)
summary(fcIKK[fcIKK$variable=="KOIR1" & fcIKK$group=="exc_KOIR",]$value)

summary(fcIKK[fcIKK$variable=="CKO1" & fcIKK$group=="exc_WTIR",]$value)
summary(fcIKK[fcIKK$variable=="CWT1" & fcIKK$group=="exc_WTIR",]$value)
summary(fcIKK[fcIKK$variable=="WTIR" & fcIKK$group=="exc_WTIR",]$value)
summary(fcIKK[fcIKK$variable=="KOIR1" & fcIKK$group=="exc_WTIR",]$value)

wilcox.test(fcIKK[fcIKK$variable=="CKO1" & fcIKK$group=="common_genes",]$value,fcIKK[fcIKK$variable=="CWT1" & fcIKK$group=="common_genes",]$value)

wilcox.test(fcIKK[fcIKK$variable=="CKO1" & fcIKK$group=="exc_KOIR",]$value,fcIKK[fcIKK$variable=="CWT1" & fcIKK$group=="exc_KOIR",]$value)

wilcox.test(fcIKK[fcIKK$variable=="CKO1" & fcIKK$group=="exc_WTIR",]$value,fcIKK[fcIKK$variable=="CWT1" & fcIKK$group=="exc_WTIR",]$value)

wilcox.test(fcIKK[fcIKK$variable=="WTIR" & fcIKK$group=="common_genes",]$value,fcIKK[fcIKK$variable=="KOIR1" & fcIKK$group=="common_genes",]$value)

wilcox.test(fcIKK[fcIKK$variable=="WTIR" & fcIKK$group=="exc_KOIR",]$value,fcIKK[fcIKK$variable=="KOIR1" & fcIKK$group=="exc_KOIR",]$value)

wilcox.test(fcIKK[fcIKK$variable=="WTIR" & fcIKK$group=="exc_WTIR",]$value,fcIKK[fcIKK$variable=="KOIR1" & fcIKK$group=="exc_WTIR",]$value)


ggplot(fcIKK,aes(x=variable,y=value))+
  geom_point(aes(color=variable))+
  geom_violin(aes(color=variable),alpha=0.2)+
  facet_wrap(~group,scales="free")+
  scale_color_brewer(palette="Set2")+
  theme_light()
  #scale_y_continuous(limits=c(0,25))


### Multibigwig correlation ###
cormul<-read.delim("/Volumes/cancer/CKC/peakcall/norm_bigwig/scores_CWT1_CKO1_WTIR_KOIR_compare.tab")
colnames(cormul)<-c("chr","start","end","CWT1","CWC2","CKO1","CKC2","WTIR","KOIR1")

#Remove all zeros
remwin<-row.names(cormul[cormul$CWT1==0 & cormul$CWC2==0 & cormul$CKO1==0 & cormul$CKC2==0 & cormul$WTIR==0 & cormul$KOIR1==0,])
cormul<-cormul[!row.names(cormul) %in% remwin,]

summary(cormul$CWT1)
summary(cormul$CKO1)
summary(cormul$WTIR)
summary(cormul$KOIR1)

## filter outliers that are input peaks likely too
cormul_filt<-cormul[cormul$CWT1>=5 & cormul$CWT1<=40 & cormul$CKO1<=40 & cormul$WTIR<=40 & cormul$KOIR<=40,]

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

ggpairs(
  data=ranchr,mapping=aes(color = chr),
  upper = list(continuous = wrap("cor", alpha = 1,size=2), combo = wrap("box_no_facet",alpha=0.4)),
  lower = list(continuous = wrap("points", alpha = 0.3,size=0.1), 
               combo = wrap("dot", alpha = 0.4,size=0.2)),cardinality_threshold = 22)+theme_bw()


#Considering all windows independently of chromosomes
ggpairs(
  data=cormul_filt[,c(4,6,8,9)],
  upper = list(continuous = wrap("cor", alpha = 1,size=4), combo = wrap("box_no_facet",alpha=0.4)),
  lower = list(continuous = wrap("points", alpha = 0.3,size=0.1), 
               combo = wrap("dot", alpha = 0.4,size=0.2)),cardinality_threshold = 22)+theme_bw()



# Plot pair panels
pairs.panels(cormul_filt[,c(4,6,8,9)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             pch=".",
             cex.cor = 0.5,
             points.false=TRUE
)


# Correlation matrix corrected for multiple values
cor(cormul_filt[,c(4,6,8,9)])
corPvalueFisher(cor(cormul_filt[,c(4,6,8,9)]),4)


# 3D plot
plot_ly(x = cormul_filt$CWT1, y = cormul_filt$WTIR, z = cormul_filt$KOIR,marker=list(size=2)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'CWT1'),
                      yaxis = list(title = 'WTIR'),
                      zaxis = list(title = 'KOIR')))


## predictive model for new peaks in WTIR or KOIR

#model
model_WT_KO<-lm(CWT1~CKO1,data=cormul_filt)
model_WT_KO

#observed WTIR (constant). Change name to CKO1 to match with the model
obs_WTIR<-data.frame(
  WTIR=c(cormul_filt[,colnames(cormul_filt) == "WTIR"]))

row.names(obs_WTIR)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
colnames(obs_WTIR)<-"CKO1"

# prediction for WTIR
pre_WTIR<-data.frame(predict(model_WT_KO,
                             level=0.99,
                   newdata = obs_WTIR,
                   interval="prediction"))

colnames(pre_WTIR)<-c("fit_CWT1_from_WTIR","lwr_CWT1_from_WTIR","upr_CWT1_from_WTIR")

## Merge predictive CWT1 from model using observe WTIR with observed values CWT1
row.names(cormul_filt)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
cormul_filt_model<-merge(cormul_filt,pre_WTIR,by=0)

row.names(cormul_filt_model)<-cormul_filt_model$Row.names
cormul_filt_model$Row.names<-NULL


#observed KOIR (constant). Change name to CKO1 to match with the model
obs_KOIR<-data.frame(
  KOIR=c(cormul_filt[,colnames(cormul_filt) == "KOIR1"]))

row.names(obs_KOIR)<-paste(cormul_filt$chr,cormul_filt$start,sep="_")
colnames(obs_KOIR)<-"CKO1"

# prediction for KOIR
pre_KOIR<-data.frame(predict(model_WT_KO,
                             level=0.99,
                             newdata = obs_KOIR,
                             interval="prediction"))

colnames(pre_KOIR)<-c("fit_CWT1_from_KOIR","lwr_CWT1_from_KOIR","upr_CWT1_from_KOIR")
row.names(pre_KOIR)


## Merge models WTIR and KOIR
cormul_filt_model<-merge(cormul_filt_model,pre_KOIR,by=0)



# Genes that fit the model WTIR
cormul_filt_model$model_CWT1_from_WTIR<-ifelse(cormul_filt_model$lwr_CWT1_from_WTIR < cormul_filt_model$CWT1 &  cormul_filt_model$CWT1 < cormul_filt_model$upr_CWT1_from_WTIR, 
                             c("YES"),c("NO")) 

# Genes that have high score in WTIR but disseappear in CWT1
cormul_filt_model$fitmod_yesWTIR_noCWT1<-ifelse(cormul_filt_model$model_CWT1_from_WTIR =="NO" & cormul_filt_model$CWT1<cormul_filt_model$lwr_CWT1_from_WTIR, 
                                                c("WTIR_noCWT1"),c("WTIR_CWT1")) 


# Genes that fit the model KOIR
cormul_filt_model$model_CWT1_from_KOIR<-ifelse(cormul_filt_model$lwr_CWT1_from_KOIR < cormul_filt_model$CWT1 &  cormul_filt_model$CWT1 < cormul_filt_model$upr_CWT1_from_KOIR, 
                                               c("YES"),c("NO")) 

# Genes that have high score in KOIR but disseappear in CWT1
cormul_filt_model$fitmod_yesKOIR_noCWT1<-ifelse(cormul_filt_model$model_CWT1_from_KOIR =="NO" & cormul_filt_model$CWT1<cormul_filt_model$lwr_CWT1_from_KOIR, 
                                                c("KOIR_noCWT1"),c("KOIR_CWT1")) 

# Regions no model CWT1_WTR1, CWT1 < lwr and no model CWT1_KOIR, CWT1 < lwr
cormul_filt_model$mergemodel<-paste(cormul_filt_model$fitmod_yesWTIR_noCWT1,cormul_filt_model$fitmod_yesKOIR_noCWT1,sep="_")
table(cormul_filt_model$mergemodel)

ggpairs(
  data=cormul_filt_model[,c(5,7,9,10,23)],mapping=aes(color = mergemodel),
  upper = list(continuous = wrap("cor", alpha = 1,size=2), combo = wrap("box_no_facet",alpha=0.4)),
  lower = list(continuous = wrap("points", alpha = 0.3,size=0.1), 
               combo = wrap("dot", alpha = 0.4,size=0.2)),cardinality_threshold = 22)+theme_bw()

table(cormul_filt_model$mergemodel)
cormul_filt_model$dif<-cormul_filt_model$WTIR-cormul_filt_model$KOIR1


# In blue, the score or WTIR peaks is lower than that expected by CWT1 vs CKO1 model.
# In coral, the score of WTIR peaks is higher than that expected by CWT1 vs CKO1 model 
p1cor<-ggplot(data=cormul_filt_model,aes(x=CWT1,y=CKO1))+
  geom_point(color="black",size=1,alpha=0.3)+
  geom_point(data=cormul_filt_model[cormul_filt_model$CWT1>cormul_filt_model$upr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO",],color="dodgerblue",size=1,alpha=0.8)+
  geom_point(data=cormul_filt_model[cormul_filt_model$CWT1<cormul_filt_model$lwr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO",],color="coral",size=1,alpha=0.8)+
  theme_bw()

cor.test(cormul_filt_model$CWT1,cormul_filt_model$CKO1)

p2cor<-ggplot(data=cormul_filt_model,aes(x=CWT1,y=WTIR))+
  geom_point(color="black",size=1,alpha=0.3)+
  geom_point(data=cormul_filt_model[cormul_filt_model$CWT1>cormul_filt_model$upr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO",],color="dodgerblue",size=1,alpha=0.8)+
  geom_point(data=cormul_filt_model[cormul_filt_model$CWT1<cormul_filt_model$lwr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO",],color="coral",size=1,alpha=0.8)+
  geom_hline(yintercept = min(cormul_filt_model[cormul_filt_model$CWT1<cormul_filt_model$lwr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO",]$WTIR),color="white")+
  theme_bw()

cor.test(cormul_filt_model$CWT1,cormul_filt_model$WTIR)

# In blue peaks with lower score in KOIR than expected in CWT or CKO
# In purple peaks with higher score in KOIR than expected in CWT or CKO
p2acor<-ggplot(data=cormul_filt_model,aes(x=CKO1,y=KOIR1))+
  geom_point(color="black",size=1,alpha=0.3)+
  geom_point(data=cormul_filt_model[cormul_filt_model$CWT1>cormul_filt_model$upr_CWT1_from_KOIR & cormul_filt_model$model_CWT1_from_KOIR=="NO",],color="dodgerblue",size=1,alpha=0.8)+
  geom_point(data=cormul_filt_model[cormul_filt_model$CWT1<cormul_filt_model$lwr_CWT1_from_KOIR & cormul_filt_model$model_CWT1_from_KOIR=="NO",],color="darkorchid1",size=1,alpha=0.8)+
  geom_hline(yintercept = min(cormul_filt_model[cormul_filt_model$CWT1<cormul_filt_model$lwr_CWT1_from_KOIR & cormul_filt_model$model_CWT1_from_KOIR=="NO",]$KOIR),color="white")+
  theme_bw()

cor.test(cormul_filt_model$CKO1,cormul_filt_model$KOIR)

p3cor<-ggplot(data=cormul_filt_model,aes(x=KOIR1,y=WTIR))+
  geom_point(color="black",size=1,alpha=0.3)+
  geom_point(data=cormul_filt_model[cormul_filt_model$mergemodel=="WTIR_noCWT1_KOIR_CWT1",],color="darkslateblue",size=1,alpha=0.8)+
  geom_point(data=cormul_filt_model[cormul_filt_model$mergemodel=="WTIR_CWT1_KOIR_noCWT1",],color="goldenrod3",size=1,alpha=0.8)+
  theme_bw()

cor.test(cormul_filt_model$WTIR,cormul_filt_model$KOIR)

p3cor


meltdat<-cormul_filt_model[cormul_filt_model$mergemodel=="WTIR_noCWT1_KOIR_CWT1",c(9,10)]
meltdat<-melt(meltdat)

pden<-ggplot(meltdat,aes(x=value))+
  geom_density(aes(color=variable),size=1,alpha=0.3)+
  scale_color_manual(values=c("coral","gray"))+
  geom_rug(aes(color=variable),size=0.7)+
  theme_bw()
pden

grid.arrange(p1cor,p2cor,p2acor,p3cor,ncol=2)




# Write data.table with all regions
write.table(cormul_filt_model[,c(2:4,23)],"/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/allregions.bed",quote = FALSE,row.names = FALSE,sep="\t")

# Write data table with peaks in WTIR, no CWT1 or CKO1, and no KOIR to annotate later
write.table(cormul_filt_model[cormul_filt_model$mergemodel=="WTIR_noCWT1_KOIR_CWT1",c(2:4,23)],"/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_WTIR.bed",quote = FALSE,row.names = FALSE,sep="\t")

# Write data table with peaks in WTIR and KOIR, but no CWT1 or CKO1
write.table(cormul_filt_model[cormul_filt_model$mergemodel=="WTIR_noCWT1_KOIR_noCWT1",c(2:4,23)],"/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_WTIR_KOIR.bed",quote = FALSE,row.names = FALSE,sep="\t")

# Write data table with peaks in CWT1 but no in WTIR
write.table(cormul_filt_model[cormul_filt_model$CWT1>cormul_filt_model$upr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO",c(2:4,23)],"/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_CWT1_noWTRmodel.bed",quote = FALSE,row.names = FALSE,sep="\t")

# Write data table with peaks in CKO1 but no in KOIR or WTIR
#write.table(cormul_filt_model[cormul_filt_model$CKO1>cormul_filt_model$upr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO" & cormul_filt_model$model_CWT1_from_KOIR=="NO",c(2:4,23)],"/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_CWT1_noWTRmodel.bed",quote = FALSE,row.names = FALSE,sep="\t")
write.table(cormul_filt_model[cormul_filt_model$CKO1>cormul_filt_model$upr_CWT1_from_WTIR & cormul_filt_model$model_CWT1_from_WTIR=="NO" & cormul_filt_model$model_CWT1_from_KOIR=="NO",c(2:4,23)],"/Users/yguillen/Desktop/temp/CKC_Irene/ChIPs_August_2020/genes_peaks/regions_no_model_CKO1_noWTRmodel_noKOIRmodel.bed",quote = FALSE,row.names = FALSE,sep="\t")

# Write data table with peaks in KOIR but no in WTIR, and not WT or KO
write.table(cormul_filt_model[cormul_filt_model$mergemodel=="WTIR_CWT1_KOIR_noCWT1",c(2:4,23)],"/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_KOIR.bed",quote = FALSE,row.names = FALSE,sep="\t")



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
Brd4_all<-read_regions(con="/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/allregions.bed", genome = 'mm10', format = 'bed')

all_annotated = annotate_regions(
  regions = Brd4_all,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## Regions exclusive WTIR, not KOIR, according to model
WTIR_regions<-read_regions(con = "/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_WTIR.bed", genome = 'mm10', format = 'bed')

WTIR_annotated = annotate_regions(
  regions = WTIR_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## Regions exclusive WTIR and KOIR, not CWT1 or CKO1, according to model
WTIR_KOIR_regions<-read_regions(con = "/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_WTIR_KOIR.bed", genome = 'mm10', format = 'bed')

WTIR_KOIR_annotated = annotate_regions(
  regions = WTIR_KOIR_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## Regions exclusive in CWT1, but no WTIR according to the model (plot 2, blue dots)
CWT1_regions<-read_regions(con = "/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_CWT1_noWTRmodel.bed", genome = 'mm10', format = 'bed')

CWT1_annotated = annotate_regions(
  regions = CWT1_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## Regions exclusive in CKO1, but no WTIR or KOIR according to the model
CKO1_regions<-read_regions(con = "/Users/yguillen/Desktop/temp/CKC_Irene/ChIPs_August_2020/genes_peaks/regions_no_model_CKO1_noWTRmodel_noKOIRmodel.bed", genome = 'mm10', format = 'bed')

CKO1_annotated = annotate_regions(
  regions = CKO1_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


## Regions exclusive in KOIR, but no WTIR according to the model (plot 2, blue dots) and no CWT1 or CKO1
KOIR_regions<-read_regions(con = "/Volumes/grcmc/IRENE PECHARROMAN/YGuillen/ChIPs_August_2020/genes_peaks/regions_no_model_KOIR.bed", genome = 'mm10', format = 'bed')

KOIR_annotated = annotate_regions(
  regions = KOIR_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)




# Compare with random regions
# create random regions
library(regioneR)
library(BSgenome.Mmusculus.UCSC.mm10)

random_regions<-createRandomRegions(nregions=length(WTIR_regions$name),
                    length.mean=1000,
                    length.sd=20,
                    genome="mm10", mask=NULL, non.overlapping=TRUE)

random_annotated = annotate_regions(
  regions = random_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


WTIR_annotated_wrandom = plot_annotation(
  annotated_regions = WTIR_annotated,
  annotated_random = random_annotated,
  annotation_order = annots_order,
  plot_title = 'WTIR_annotated vs random',
  x_label = 'Annotations',
  y_label = 'Count')

comp_ran<-print(WTIR_annotated_wrandom)


WTIR_annotated_WTIR_KOIR = plot_annotation(
  annotated_regions = WTIR_annotated,
  annotated_random = WTIR_KOIR_annotated,
  annotation_order = annots_order,
  plot_title = 'WTIR_annotated vs common_WTIR_KOIR',
  x_label = 'Annotations',
  y_label = 'Count')

comp_WTIRKOIR<-print(WTIR_annotated_WTIR_KOIR)


WTIR_annotated_CWT1 = plot_annotation(
  annotated_regions = WTIR_annotated,
  annotated_random = CWT1_annotated,
  annotation_order = annots_order,
  plot_title = 'WTIR_annotated vs CWT1_no_WTIR',
  x_label = 'Annotations',
  y_label = 'Count')

comp_WTIR_CWT1<-print(WTIR_annotated_CWT1)


grid.arrange(comp_ran,comp_WTIRKOIR,comp_WTIR_CWT1,ncol=2)


# Regions in pairs of annotations
WTIR_coannotations = plot_coannotations(
  annotated_regions = WTIR_annotated,
  annotation_order = annots_order,
  axes_label = 'Annotations',
  plot_title = 'Regions in Pairs of Annotations')

print(WTIR_coannotations)


# proportions

all_annotated<-data.frame(all_annotated)
propall_annot<-data.frame(prop.table(table(all_annotated$annot.type)))
propall_annot$sam<-"all_bigwig"

WTIR_annotated<-data.frame(WTIR_annotated)
propWTIR<-data.frame(prop.table(table(WTIR_annotated$annot.type)))
propWTIR$sam<-"WTIR"

WTIR_KOIR_annotated<-data.frame(WTIR_KOIR_annotated)
propWTIR_KOIR<-data.frame(prop.table(table(WTIR_KOIR_annotated$annot.type)))
propWTIR_KOIR$sam<-"WTIR_KOIR"

CWT1_annotated<-data.frame(CWT1_annotated)
propCWT1<-data.frame(prop.table(table(CWT1_annotated$annot.type)))
propCWT1$sam<-"CWT1"

CKO1_annotated<-data.frame(CKO1_annotated)
propCKO1<-data.frame(prop.table(table(CKO1_annotated$annot.type)))
propCKO1$sam<-"CKO1"

KOIR_annotated<-data.frame(KOIR_annotated)
propKOIR<-data.frame(prop.table(table(KOIR_annotated$annot.type)))
propKOIR$sam<-"KOIR"

random_annotated<-data.frame(random_annotated)
proprandom<-data.frame(prop.table(table(random_annotated$annot.type)))
proprandom$sam<-"random"


propall<-rbind(propWTIR,propWTIR_KOIR)
propall<-rbind(propall,propCWT1)
propall<-rbind(propall,proprandom)
propall<-rbind(propall,propKOIR)
propall<-rbind(propall,propall_annot)
propall<-rbind(propall,propCKO1)

propall$sam<-factor(propall$sam,levels = c("random","all_bigwig","CWT1","CKO1","WTIR_KOIR","WTIR","KOIR"))

ggplot(propall,aes(x=sam,y=Freq,fill=Var1))+
  geom_bar(stat="identity",color="black") +
  scale_fill_brewer(palette="Paired")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Frequency annotation")


# merge with scores models

cormul_chr<-subset(cormul_filt_model,select=c("chr","start","CWT1","CKO1","WTIR","KOIR1","mergemodel","Row.names"))
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

## peaks called with chipseeker in ChIPhist_CKC.R script
allmacs<-rbind(CWT1_aug_df,CKO1_aug_df,WTIR_aug_df[,-c(20)],KOIR_aug_df)
allmacs<-subset(allmacs,select=c("geneChr","start","Name","Score","annotation","SYMBOL","type"))
colnames(allmacs)[7]<-c("mergemodel")
allmacs$mergemodel<-gsub('^','MACS2',allmacs$mergemodel)
table(allmacs$mergemodel)

allmacs$Name<-gsub('_.*','',allmacs$Name)
allmacs<-allmacs[,c(1,2,7,3,4,6,5)]
colnames(allmacs)=colnames(all_annotated_model)

allmacs$chr<-as.factor(allmacs$chr)

# change scale of score to represent with bigwig values
allmacs$value<-allmacs$value/2

all_annotated_model_macs<-rbind(all_annotated_model,allmacs)


# Choose all_annotated_model_macs or all_annotated_model with no MACS peaks info
ggplot(all_annotated_model,aes(x=start,y=value))+
  geom_point(data=all_annotated_model[all_annotated_model$mergemodel=="WTIR_noCWT1_KOIR_CWT1",],aes(color=variable,shape=annot.type),size=1)+
  scale_color_brewer(palette = "Dark2",na.translate=F)+
  facet_wrap(~chr,nrow=3,scales="free_x")+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  ylab("bigwig score")+
  ggtitle("Higher in WTIR vs other")


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
enr_Brd4 <- enriched_Brd4[["GO_Biological_Process_2018"]]
enr_Brd4$cat<-"Brd4"

enriched_random <- enrichr(unique(random_chiptab@anno$SYMBOL), dbs)
enr_random <- enriched_random[["GO_Biological_Process_2018"]]
enr_random$cat<-"random"

enriched_CWT1 <- enrichr(unique(CWT1_chiptab@anno$SYMBOL), dbs)
enr_CWT1 <- enriched_CWT1[["GO_Biological_Process_2018"]]
enr_CWT1$cat<-"CWT1"

enriched_CKO1 <- enrichr(unique(CKO1_chiptab@anno$SYMBOL), dbs)
enr_CKO1 <- enriched_CKO1[["GO_Biological_Process_2018"]]
enr_CKO1$cat<-"CKO1"


enriched_WTIR_KOIR <- enrichr(unique(WTIR_KOIR_chiptab@anno$SYMBOL), dbs)
enr_WTIR_KOIR <- enriched_WTIR_KOIR[["GO_Biological_Process_2018"]]
enr_WTIR_KOIR$cat<-"WTIR_KOIR"


enriched_WTIR <- enrichr(unique(WTIR_chiptab@anno$SYMBOL), dbs)
enr_WTIR <- enriched_WTIR[["GO_Biological_Process_2018"]]
enr_WTIR$cat<-"WTIR"

enriched_KOIR <- enrichr(unique(KOIR_chiptab@anno$SYMBOL), dbs)
enr_KOIR <- enriched_KOIR[["GO_Biological_Process_2018"]]
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


# Heatmap with only ChEA
ggplot(bpsub[grepl('CHEA',bpsub$Term),], aes(x=cat, y=Term)) + 
  geom_tile(aes(fill = -(log(P.value))),colour = "black")+
  scale_fill_gradient2(low = "white",
                       high = "red") +
  #scale_fill_gradient(low = "blue", high = "red")+
  theme_dark()+
  theme(axis.text.x = element_text(size=10, hjust = 1,angle=45),
        axis.text.y = element_text(size=10, hjust = 1,face="italic"))+
  ggtitle("ChEA_Consensus_TFs_from_ChIP-X")

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


