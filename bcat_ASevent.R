### SEPTEMBER 2020 ### Yolanda Guillén
# Script for exploring b-cat regulated splicing events candidates.


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

# vastools diff alt splicing RPMI
#RPMI_dev<-read.delim('/Volumes/cancer/Gekas_RNAseq/vast_output/RPMI_PSI20/RPMI_shbcat_minpSI20.txt',header = TRUE)
RPMI_dev<-read.delim('/Users/yguillen/Desktop/temp/Gekas_RNASeq/vast_output/RPMI_PSI20/RPMI_shbcat_minpSI20.txt',header = TRUE)

columnnames<-colnames(RPMI_dev)
RPMI_dev$type<-"RPMI"

RPMI_dev$dif_rpmi<-(RPMI_dev$shbcat_RPMI)-(RPMI_dev$control_RPMI)
RPMI_dev$dif_Jur<-(RPMI_dev$shbcat_Jurkat)-(RPMI_dev$control_Jurkat)
RPMI_dev$sense<-RPMI_dev$dif_Jur*RPMI_dev$dif_rpmi



# Selecting same trend RPMI and Jurkat
#RPMI_dev<-RPMI_dev[!is.na(RPMI_dev$sense) & RPMI_dev$sense>0,]
#RPMI_dev<-RPMI_dev[with(RPMI_dev,rev(order(RPMI_dev$sense))),]

RPMI_dev$state<-ifelse(RPMI_dev$dif_rpmi<0,c("shbcat_down"),c("shbcat_up"))

table(RPMI_dev$COMPLEX,RPMI_dev$state)

table(RPMI_dev[grepl("EX",RPMI_dev$EVENT),]$state)
table(RPMI_dev[grepl("INT",RPMI_dev$EVENT),]$state)
table(RPMI_dev[grepl("ALT",RPMI_dev$EVENT),]$state)


## vastools Background splicing events that do not change between conditions in RPMI
#RPMI_BG<-read.delim('/Volumes/cancer/Gekas_RNAseq/vast_output/RPMI_PSI20/AS_NC-RPMI_minPSI20-Max_dPSI4.tab',header = FALSE)
RPMI_BG<-read.delim('/Users/yguillen/Desktop/temp/Gekas_RNASeq/vast_output/RPMI_PSI20/AS_NC-RPMI_minPSI20-Max_dPSI4.tab',header = FALSE)

colnames(RPMI_BG)<-columnnames

table(RPMI_BG$COMPLEX)


### Generate output for rmap2
#UP AS shbcat in RPMI
uprmap<-RPMI_dev[RPMI_dev$state=="shbcat_up",]
uprmap_exskip<-uprmap[uprmap$COMPLEX=="C1" | uprmap$COMPLEX=="C2" | uprmap$COMPLEX=="C3" | uprmap$COMPLEX=="S",]

up_rpmi_exsk<-as.data.frame(uprmap_exskip$FullCO)
colnames(up_rpmi_exsk)<-"UP_AS_RPMI"

write.table(up_rpmi_exsk,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/rMATTs_bcat/RPMIsig/up_rpmiall_exsp.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

#DOWN AS shbcat in RPMI
downrmap<-RPMI_dev[RPMI_dev$state=="shbcat_down",]
downrmap_exskip<-downrmap[downrmap$COMPLEX=="C1" | downrmap$COMPLEX=="C2" | downrmap$COMPLEX=="C3" | downrmap$COMPLEX=="S",]

down_rpmi_exsk<-as.data.frame(downrmap_exskip$FullCO)
colnames(down_rpmi_exsk)<-"DOWN_AS_RPMI"
write.table(down_rpmi_exsk,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/rMATTs_bcat/RPMIsig/down_rpmiall_exsp.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)


#No changing AS shbcat in RPMI
RPMI_BG_exskip<-RPMI_BG[RPMI_BG$COMPLEX=="C1" | RPMI_BG$COMPLEX=="C2" | RPMI_BG$COMPLEX=="C3" | RPMI_BG$COMPLEX=="S",]

BG_rpmi_exsk<-as.data.frame(RPMI_BG_exskip$FullCO)
colnames(BG_rpmi_exsk)<-"BG_AS_RPMI"
write.table(BG_rpmi_exsk,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/rMATTs/BG_AS_RPMI.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)



# vastools diff alt splicing Jurkat
Jurkat_dev<-read.delim('/Volumes/cancer/Gekas_RNAseq/vast_output/Jurkat_PSI20/Jurkat_minPSI20.txt',header = TRUE)
columnnames<-colnames(Jurkat_dev)
Jurkat_dev$type<-"Jurkat"

Jurkat_dev$dif_rpmi<-(Jurkat_dev$shbcat_RPMI)-(Jurkat_dev$control_RPMI)
Jurkat_dev$dif_Jur<-(Jurkat_dev$shbcat_Jurkat)-(Jurkat_dev$control_Jurkat)
Jurkat_dev$sense<-Jurkat_dev$dif_Jur*Jurkat_dev$dif_rpmi



# Selecting same trend RPMI and Jurkat
#Jurkat_dev<-Jurkat_dev[!is.na(Jurkat_dev$sense) & Jurkat_dev$sense>0,]
#Jurkat_dev<-Jurkat_dev[with(Jurkat_dev,rev(order(Jurkat_dev$sense))),]

Jurkat_dev$state<-ifelse(Jurkat_dev$dif_Jur<0,c("shbcat_down"),c("shbcat_up"))

table(Jurkat_dev$COMPLEX,Jurkat_dev$state)

table(Jurkat_dev[grepl("EX",Jurkat_dev$EVENT),]$state)
table(Jurkat_dev[grepl("INT",Jurkat_dev$EVENT),]$state)
table(Jurkat_dev[grepl("ALT",Jurkat_dev$EVENT),]$state)


## vastools Background splicing events that do not change between conditions in RPMI
Jurkat_BG<-read.delim('/Volumes/cancer/Gekas_RNAseq/vast_output/Jurkat_PSI20/AS_NC-Jurkat_minPSI20-Max_dPSI4.tab',header = FALSE)
colnames(Jurkat_BG)<-columnnames

table(Jurkat_BG$COMPLEX)


### Generate output for rmap2
#UP AS shbcat in Jurkat
uprmap<-Jurkat_dev[Jurkat_dev$state=="shbcat_up",]
uprmap_exskip<-uprmap[uprmap$COMPLEX=="C1" | uprmap$COMPLEX=="C2" | uprmap$COMPLEX=="C3" | uprmap$COMPLEX=="S",]

up_Jurkat_exsk<-as.data.frame(uprmap_exskip$FullCO)
colnames(up_Jurkat_exsk)<-"UP_AS_Jurkat"

write.table(up_Jurkat_exsk,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/rMATTs_bcat/Jurkatsig/up_Jurkatall_exsp.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)

#DOWN AS shbcat in RPMI
downrmap<-Jurkat_dev[Jurkat_dev$state=="shbcat_down",]
downrmap_exskip<-downrmap[downrmap$COMPLEX=="C1" | downrmap$COMPLEX=="C2" | downrmap$COMPLEX=="C3" | downrmap$COMPLEX=="S",]

down_Jurkat_exsk<-as.data.frame(downrmap_exskip$FullCO)
colnames(down_Jurkat_exsk)<-"DOWN_AS_RPMI"
write.table(down_Jurkat_exsk,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/rMATTs_bcat/Jurkatsig/down_Jurkatall_exsp.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)


#No changing AS shbcat in RPMI
Jurkat_BG_exskip<-Jurkat_BG[Jurkat_BG$COMPLEX=="C1" | Jurkat_BG$COMPLEX=="C2" | Jurkat_BG$COMPLEX=="C3" | Jurkat_BG$COMPLEX=="S",]

BG_Jurkat_exsk<-as.data.frame(Jurkat_BG_exskip$FullCO)
colnames(BG_Jurkat_exsk)<-"BG_AS_Jurkat"
write.table(BG_Jurkat_exsk,"/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/rMATTs_bcat/Jurkatsig/BG_AS_Jurkat.tab",quote = FALSE,row.names = FALSE,col.names = FALSE)






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





