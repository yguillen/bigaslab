### NOVEMBER 2020 ### Yolanda Guillén
# Script for exploring T-ALL vs thymus dysregulated splicing events candidates.


#basic libraries
library(data.table)
library(ggplot2)
library(nlme)
library(readxl)
library(ggrepel)
library(ggnewscale)
library(GO.db)
library(topGO)
library(gridExtra)
library(Hmisc)
library(corrplot)
library(DESeq2)
library(gplots)
library(genefilter)
library(gRbase)
library(plyr)
library(pheatmap)
library(heatmaply)
library(limma)
library(remotes)
library(dplyr)
library(stringr)
library(tidyr)


### EGA_TLE
#####
# vastools diff alt splicing EGA_TLE
EGATLE_dev<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/EGA_TLE_expr_out/isoforms_diff/EGA_TLE_diffisoforms.txt',header = TRUE)
colnames(EGATLE_dev)<-gsub('EGAR.*_','',colnames(EGATLE_dev))

columnnames<-colnames(EGATLE_dev)
EGATLE_dev$cohort<-"EGA_TLE"

#select only EGA TLE samples
egasamp<-metadata[metadata$GSEA=="EGA_TLE",]$SampID
egasamp<-which(colnames(EGATLE_dev) %in% egasamp)
EGATLE_dev<-EGATLE_dev[,c(1:6,egasamp)]

EGATLE_dev_med<-EGATLE_dev[,c(1:2,6,11,7:10,12:ncol(EGATLE_dev))]
EGATLE_dev_med$med<-apply(EGATLE_dev_med[-c(1:4)], 1, FUN=median, na.rm=TRUE)
EGATLE_dev_med$dif<-EGATLE_dev_med$Thymus-EGATLE_dev_med$med
EGATLE_dev_med$state<-ifelse(EGATLE_dev_med$dif>0,c("TALL_down"),c("TALL_UP"))

table(EGATLE_dev_med[grepl("EX",EGATLE_dev_med$EVENT),]$state)
table(EGATLE_dev_med[grepl("INT",EGATLE_dev_med$EVENT),]$state)
table(EGATLE_dev_med[grepl("ALT",EGATLE_dev_med$EVENT),]$state)

table(EGATLE_dev_med$COMPLEX,EGATLE_dev_med$state)



## vastools Background splicing events that do not change between conditions in RPMI
EGATLE_BG<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/EGA_TLE_expr_out/isoforms_diff/AS_NC-EGA_TLE_diffisoforms-Max_dPSI2.tab',header = FALSE)
colnames(EGATLE_BG)<-columnnames

table(EGATLE_BG$COMPLEX)


### Generate output for rmap2
#UP IR and EX in TALL
EGATLE_UPrmap<-EGATLE_dev[which(EGATLE_dev_med$state=="TALL_UP"),]

EGATLE_UP_IRrmap<-EGATLE_UPrmap[EGATLE_UPrmap$COMPLEX=="IR-C" | EGATLE_UPrmap$COMPLEX=="IR-S",]
EGATLE_UP_IRrmap<-as.data.frame(EGATLE_UP_IRrmap$FullCO)
colnames(EGATLE_UP_IRrmap)<-"UP_IR_TALL_EGATLE"

write.table(EGATLE_UP_IRrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/IR/up_IR_TALL_EGATLE.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

EGATLE_UP_EXrmap<-EGATLE_UPrmap[EGATLE_UPrmap$COMPLEX=="C1" | EGATLE_UPrmap$COMPLEX=="C2" | EGATLE_UPrmap$COMPLEX=="C3" | EGATLE_UPrmap$COMPLEX=="S",]

EGATLE_UP_EXrmap<-as.data.frame(EGATLE_UP_EXrmap$FullCO)
colnames(EGATLE_UP_EXrmap)<-"UP_EX_TALL_EGATLE"

write.table(EGATLE_UP_EXrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/ES/up_EX_TALL_EGATLE.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)



#DOWN IR and EX in TALL
EGATLE_DOWNrmap<-EGATLE_dev[which(EGATLE_dev_med$state=="TALL_down"),]

EGATLE_DOWN_IRrmap<-EGATLE_DOWNrmap[EGATLE_DOWNrmap$COMPLEX=="IR-S" | EGATLE_DOWNrmap$COMPLEX=="IR-C",]

EGATLE_DOWN_IRrmap<-as.data.frame(EGATLE_DOWN_IRrmap$FullCO)
colnames(EGATLE_DOWN_IRrmap)<-"DOWN_IR_TALL_EGATLE"

write.table(EGATLE_DOWN_IRrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/IR/down_IR_TALL_EGATLE.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

EGATLE_DOWN_EXrmap<-EGATLE_DOWNrmap[EGATLE_DOWNrmap$COMPLEX=="C1" | EGATLE_DOWNrmap$COMPLEX=="C2" | EGATLE_DOWNrmap$COMPLEX=="C3" | EGATLE_DOWNrmap$COMPLEX=="S",]

EGATLE_DOWN_EXrmap<-as.data.frame(EGATLE_DOWN_EXrmap$FullCO)
colnames(EGATLE_DOWN_EXrmap)<-"DOWN_EX_TALL_EGATLE"

write.table(EGATLE_DOWN_EXrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/ES/down_IR_TALL_EGATLE.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)


#No changing AS in EGA_TLE
# Intron
EGATLE_BG_INT<-EGATLE_BG[EGATLE_BG$COMPLEX=="IR-S" | EGATLE_BG$COMPLEX=="IR-C",]

EGATLE_BG_INT<-as.data.frame(EGATLE_BG_INT$FullCO)
colnames(EGATLE_BG_INT)<-"BG_IR_EGATLE"
write.table(EGATLE_BG_INT,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/BG_IR_EGATLE.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

# Exon
EGATLE_BG_EX<-EGATLE_BG[EGATLE_BG$COMPLEX=="C1" | EGATLE_BG$COMPLEX=="C2" | EGATLE_BG$COMPLEX=="C3" | EGATLE_BG$COMPLEX=="S",]

EGATLE_BG_EX<-as.data.frame(EGATLE_BG_EX$FullCO)
colnames(EGATLE_BG_EX)<-"BG_EX_EGATLE"
write.table(EGATLE_BG_EX,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/BG_EX_EGATLE.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)



#####

### GSE57982
#####
# vastools diff alt splicing GSE57982
GSE57982_dev<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/GSE57982_expr_out/isoforms_diff/GSE57982_minPSI10.txt',header = TRUE)
colnames(GSE57982_dev)<-gsub('EGAR.*_','',colnames(GSE57982_dev))
colnames(GSE57982_dev)<-gsub('_1','',colnames(GSE57982_dev))

columnnames<-colnames(GSE57982_dev)
GSE57982_dev$cohort<-"GSE57982"

#select only GSE57982 samples
GSE57982samp<-metadata[metadata$GSEA=="GSE57982",]$SampID
GSE57982samp<-which(colnames(GSE57982_dev) %in% GSE57982samp)
GSE57982_dev<-GSE57982_dev[,c(1:6,GSE57982samp)]

GSE57982_dev_med<-GSE57982_dev[,c(1:2,6:ncol(GSE57982_dev))]
GSE57982_dev_med$medthym<-apply(GSE57982_dev_med[-c(1:3,6:ncol(GSE57982_dev_med))], 1, FUN=median, na.rm=TRUE)
GSE57982_dev_med$medtall<-apply(GSE57982_dev_med[-c(1:5,ncol(GSE57982_dev_med))], 1, FUN=median, na.rm=TRUE)
GSE57982_dev_med$dif<-GSE57982_dev_med$medthym-GSE57982_dev_med$medtall
GSE57982_dev_med$state<-ifelse(GSE57982_dev_med$dif>0,c("TALL_down"),c("TALL_UP"))

table(GSE57982_dev_med[grepl("EX",GSE57982_dev_med$EVENT),]$state)
table(GSE57982_dev_med[grepl("INT",GSE57982_dev_med$EVENT),]$state)
table(GSE57982_dev_med[grepl("ALT",GSE57982_dev_med$EVENT),]$state)

table(GSE57982_dev_med$COMPLEX,GSE57982_dev_med$state)



## vastools Background splicing events that do not change between conditions in RPMI
GSE57982_BG<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/GSE57982_expr_out/isoforms_diff/AS_NC-GSE57982_minPSI10-Max_dPSI2.tab',header = FALSE)
colnames(GSE57982_BG)<-columnnames

table(GSE57982_BG$COMPLEX)


### Generate output for rmap2
#UP IR in TALL
GSE57982_UP_IRrmap<-GSE57982_dev[which(GSE57982_dev_med$state=="TALL_UP"),]
GSE57982_UP_IRrmap<-GSE57982_UP_IRrmap[GSE57982_UP_IRrmap$COMPLEX=="IR-C" | GSE57982_UP_IRrmap$COMPLEX=="IR-S",]

GSE57982_UP_IRrmap<-as.data.frame(GSE57982_UP_IRrmap$FullCO)
colnames(GSE57982_UP_IRrmap)<-"UP_IR_TALL_GSE57982"

write.table(GSE57982_UP_IRrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/IR/up_IR_TALL_GSE57982.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

#DOWN IR in TALL
GSE57982_DOWN_IRrmap<-GSE57982_dev[which(GSE57982_dev_med$state=="TALL_down"),]
GSE57982_DOWN_IRrmap<-GSE57982_DOWN_IRrmap[GSE57982_DOWN_IRrmap$COMPLEX=="IR-S" | GSE57982_DOWN_IRrmap$COMPLEX=="IR-C",]

GSE57982_DOWN_IRrmap<-as.data.frame(GSE57982_DOWN_IRrmap$FullCO)
colnames(GSE57982_DOWN_IRrmap)<-"DOWN_IR_TALL_GSE57982"

write.table(GSE57982_DOWN_IRrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/IR/down_IR_TALL_GSE57982.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)


#No changing AS 
GSE57982_BG_INT<-GSE57982_BG[GSE57982_BG$COMPLEX=="IR-S" | GSE57982_BG$COMPLEX=="IR-C",]

GSE57982_BG_INT<-as.data.frame(GSE57982_BG_INT$FullCO)
colnames(GSE57982_BG_INT)<-"BG_IR_GSE57982"
write.table(GSE57982_BG_INT,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/BG_IR_GSE57982.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

#####

### GSE109231
#####
# vastools diff alt splicing GSE109231
GSE109231_dev<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/GSE109231_expr_out/isoforms_diff/GSE109231_minPSI10.txt',header = TRUE)
colnames(GSE109231_dev)<-gsub('EGAR.*_','',colnames(GSE109231_dev))
colnames(GSE109231_dev)<-gsub('_1','',colnames(GSE109231_dev))

columnnames<-colnames(GSE109231_dev)
GSE109231_dev$cohort<-"GSE109231"

#select only EGA TLE samples
GSE109231samp<-metadata[metadata$GSEA=="GSE109231",]$SampID
GSE109231samp<-which(colnames(GSE109231_dev) %in% GSE109231samp)
GSE109231_dev<-GSE109231_dev[,c(1:6,GSE109231samp)]

GSE109231_dev_med<-GSE109231_dev[,c(1:2,6,10,11,7:9,12:ncol(GSE109231_dev))]
GSE109231_dev_med$medthym<-apply(GSE109231_dev_med[-c(1:3,6:ncol(GSE109231_dev_med))], 1, FUN=median, na.rm=TRUE)
GSE109231_dev_med$medtall<-apply(GSE109231_dev_med[-c(1:5,(ncol(GSE109231_dev_med)-1),ncol(GSE109231_dev_med))], 1, FUN=median, na.rm=TRUE)
GSE109231_dev_med$dif<-GSE109231_dev_med$medthym-GSE109231_dev_med$medtall
GSE109231_dev_med$state<-ifelse(GSE109231_dev_med$dif>0,c("TALL_down"),c("TALL_UP"))

table(GSE109231_dev_med[grepl("EX",GSE109231_dev_med$EVENT),]$state)
table(GSE109231_dev_med[grepl("INT",GSE109231_dev_med$EVENT),]$state)
table(GSE109231_dev_med[grepl("ALT",GSE109231_dev_med$EVENT),]$state)

table(GSE109231_dev_med$COMPLEX,GSE109231_dev_med$state)



## vastools Background splicing events that do not change between conditions in RPMI
GSE109231_BG<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/GSE109231_expr_out/isoforms_diff/AS_NC-GSE109231_minPSI10-Max_dPSI2.tab',header = FALSE)
colnames(GSE109231_BG)<-columnnames

table(GSE109231_BG$COMPLEX)


### Generate output for rmap2
#UP IR in TALL
GSE109231_UP_IRrmap<-GSE109231_dev[which(GSE109231_dev_med$state=="TALL_UP"),]
GSE109231_UP_IRrmap<-GSE109231_UP_IRrmap[GSE109231_UP_IRrmap$COMPLEX=="IR-C" | GSE109231_UP_IRrmap$COMPLEX=="IR-S",]

GSE109231_UP_IRrmap<-as.data.frame(GSE109231_UP_IRrmap$FullCO)
colnames(GSE109231_UP_IRrmap)<-"UP_IR_TALL_GSE109231"

write.table(GSE109231_UP_IRrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/IR/up_IR_TALL_GSE109231.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

#DOWN IR in TALL
GSE109231_DOWN_IRrmap<-GSE109231_dev[which(GSE109231_dev_med$state=="TALL_down"),]
GSE109231_DOWN_IRrmap<-GSE109231_DOWN_IRrmap[GSE109231_DOWN_IRrmap$COMPLEX=="IR-S" | GSE109231_DOWN_IRrmap$COMPLEX=="IR-C",]

GSE109231_DOWN_IRrmap<-as.data.frame(GSE109231_DOWN_IRrmap$FullCO)
colnames(GSE109231_DOWN_IRrmap)<-"DOWN_IR_TALL_GSE109231"

write.table(GSE109231_DOWN_IRrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/IR/down_IR_TALL_GSE109231.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)


#No changing AS in EGA_TLE
GSE109231_BG_INT<-GSE109231_BG[GSE109231_BG$COMPLEX=="IR-S" | GSE109231_BG$COMPLEX=="IR-C",]

GSE109231_BG_INT<-as.data.frame(GSE109231_BG_INT$FullCO)
colnames(GSE109231_BG_INT)<-"BG_IR_GSE109231"
write.table(GSE109231_BG_INT,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/BG_IR_GSE109231.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)
#####

## DP vs TALL
#####
# vastools diff alt splicing GSE109231
DPTALL_dev<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/TALL_DP_expr_out/isoforms_diff/TALL_DP_minPSI10.txt',header = TRUE)
colnames(DPTALL_dev)<-gsub('EGAR.*_','',colnames(DPTALL_dev))
colnames(DPTALL_dev)<-gsub('_1','',colnames(DPTALL_dev))

columnnames<-colnames(DPTALL_dev)
DPTALL_dev$cohort<-"DPTALL"

#select only EGA TLE samples
DPTALLsamp<-metadata[((metadata$GSEA=="GSE109231" | metadata$GSEA=="EGA_TLE" | metadata$GSEA=="GSE57982") & metadata$Source!="Thymus") | metadata$SampID=="DP",]$SampID
DPTALLsamp<-which(colnames(DPTALL_dev) %in% DPTALLsamp)
DPTALL_dev<-DPTALL_dev[,c(1:6,DPTALLsamp)]

DPTALL_dev_med<-DPTALL_dev[,c(1:2,6,7:ncol(DPTALL_dev))]
DPTALL_dev_med$medtall<-apply(DPTALL_dev_med[-c(1:4)], 1, FUN=median, na.rm=TRUE)
DPTALL_dev_med$dif<-DPTALL_dev_med$DP-DPTALL_dev_med$medtall
DPTALL_dev_med$state<-ifelse(DPTALL_dev_med$dif>0,c("TALL_down"),c("TALL_UP"))

table(DPTALL_dev_med[grepl("EX",DPTALL_dev_med$EVENT),]$state)
table(DPTALL_dev_med[grepl("INT",DPTALL_dev_med$EVENT),]$state)
table(DPTALL_dev_med[grepl("ALT",DPTALL_dev_med$EVENT),]$state)

table(DPTALL_dev_med$COMPLEX,DPTALL_dev_med$state)



## vastools Background splicing events that do not change between conditions in RPMI
DPTALL_BG<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/TALL_DP_expr_out/isoforms_diff/AS_NC-TALL_DP_minPSI10-Max_dPSI2.tab',header = FALSE)
colnames(DPTALL_BG)<-columnnames

table(DPTALL_BG$COMPLEX)


### Generate output for rmap2
#UP IR in TALL
DPTALL_UP_IRrmap<-DPTALL_dev[which(DPTALL_dev_med$state=="TALL_UP"),]
DPTALL_UP_IRrmap<-DPTALL_UP_IRrmap[DPTALL_UP_IRrmap$COMPLEX=="IR-C" | DPTALL_UP_IRrmap$COMPLEX=="IR-S",]

DPTALL_UP_IRrmap<-as.data.frame(DPTALL_UP_IRrmap$FullCO)
colnames(DPTALL_UP_IRrmap)<-"UP_IR_TALL_DPTALL"

write.table(DPTALL_UP_IRrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/IR/up_IR_TALL_DPTALL.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

#DOWN IR in TALL
DPTALL_DOWN_IRrmap<-DPTALL_dev[which(DPTALL_dev_med$state=="TALL_down"),]
DPTALL_DOWN_IRrmap<-DPTALL_DOWN_IRrmap[DPTALL_DOWN_IRrmap$COMPLEX=="IR-S" | DPTALL_DOWN_IRrmap$COMPLEX=="IR-C",]

DPTALL_DOWN_IRrmap<-as.data.frame(DPTALL_DOWN_IRrmap$FullCO)
colnames(DPTALL_DOWN_IRrmap)<-"DOWN_IR_TALL_DPTALL"

write.table(DPTALL_DOWN_IRrmap,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/IR/down_IR_TALL_DPTALL.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)


#No changing AS in EGA_TLE
DPTALL_BG_INT<-DPTALL_BG[DPTALL_BG$COMPLEX=="IR-S" | DPTALL_BG$COMPLEX=="IR-C",]

DPTALL_BG_INT<-as.data.frame(DPTALL_BG_INT$FullCO)
colnames(DPTALL_BG_INT)<-"BG_IR_DPTALL"
write.table(DPTALL_BG_INT,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/rMATTs_TALL/BG_IR_DPTALL.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)
#####

# Common AS events dysregulated
EGATLE_dev_med$cohort<-"EGA_TLE"
GSE109231_dev_med$cohort<-"GSE109231"
GSE57982_dev_med$cohort<-"GSE57982"
DPTALL_dev_med$cohort<-"DPTALL"

AS_TALL<-rbind(EGATLE_dev_med[,c(1:2,(ncol(EGATLE_dev_med)-2):ncol(EGATLE_dev_med))],
      GSE109231_dev_med[,c(1:2,(ncol(GSE109231_dev_med)-2):ncol(GSE109231_dev_med))],
      GSE57982_dev_med[,c(1:2,(ncol(GSE57982_dev_med)-2):ncol(GSE57982_dev_med))],
      DPTALL_dev_med[,c(1:2,(ncol(DPTALL_dev_med)-2):ncol(DPTALL_dev_med))])

as<-data.frame(table(AS_TALL$EVENT))
colnames(as)[1]<-"EVENT"

AS_TALL<-merge(as,AS_TALL,by="EVENT")
genefreq<-data.frame(table(AS_TALL$GENE))
colnames(genefreq)<-c("GENE","Freq_gene")
AS_TALL<-merge(AS_TALL,genefreq,by="GENE")

# Cross with differentially expressed genes
# DEGs_cohorts from Clinical_TCell.R
colnames(DEGs_cohorts)
colnames(DEGs_cohorts)[1]<-"GENE"
unique(DEGs_cohorts$GENE)

# AS in differentially expressed genes
Deg_AS<-AS_TALL[AS_TALL$GENE %in% unique(DEGs_cohorts$GENE),]
write.table(Deg_AS,"/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Deg_AS.tab",quote = FALSE,row.names = FALSE,sep="\t")


## FUNCTIONAL ENRICHMENT
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "Chromosome_Location",
         "GO_Molecular_Function_2018",
         "InterPro_Domains_2019",
         "GO_Cellular_Component_2017",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "WikiPathways_2016")

upskip <- enrichr(as.character(uprmap$GENE), dbs)
downskip <- enrichr(as.character(downrmap$GENE), dbs)

upskip_en <- upskip[["GO_Biological_Process_2018"]]
downskip_en <- downskip[["GO_Biological_Process_2018"]]

bpsub<-subset(downskip_en,downskip_en$P.value<0.01)
bpsub<- bpsub[order(bpsub$P.value),]


bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=unique((bpsub$Term))[rev(order(bpsub$P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))

# for one group
bpsub$Term <- factor(bpsub$Term, levels=unique(bpsub[rev(order(bpsub$P.value)), ]$Term))


ggplot(bpsub,aes(y=-log10(P.value),x=Term),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity")+
  #geom_text(aes(label=tolower(Genes)), size=2,hjust=0,color="black")+
  #scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=7),axis.title=element_text(size=14,face="bold"))+
  coord_flip()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=8, hjust = 1),
        axis.text.y = element_text(size=10, hjust = 1))





