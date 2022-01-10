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
library(scatterpie)
library(biomaRt)

set.seed(100)


#latest version
hg<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")


# Import datasets from TCell splicing project
##### 
setwd("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/") 
#setwd("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Public_data/") 

# IntoGen cancer driver genes
candriver<-read.delim("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Drivers_intoGens/Compendium_Cancer_Genes.tsv",header = TRUE)
#candriver<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Drivers_intoGens/Compendium_Cancer_Genes.tsv",header = TRUE)
unique(candriver$SYMBOL)

# metadata for EGA
#progega<-read_xlsx("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Public_data/EGA/EGAD00001000849/metadata_EGA.xlsx")
#row.names(progega)<-progega$Sample
#colnames(progega)<-c("SampID","Prog")

## Structural and mutational data from TARGET
somatic_target<-read.delim("/Users/yguillen/Documents/Projects/TARGET/stjude.org_TARGET_TALL_WXS_Diagnosis_IlluminaHiSeq_somatic.maf.txt")
structural_target<-read.delim("/Users/yguillen/Documents/Projects/TARGET/TARGET_TALL_mRNA-seq_sv.tsv")

# Selecting RBP genes form different sources
Huang_RBP<-as.data.frame(read_xlsx('/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/RBP_gene_lists/Huang_et_2018_PNAS.xlsx',sheet = "Sheet1"))
Wang_RBP<-as.data.frame(read_xlsx('/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/RBP_gene_lists/Wang_Aifantis_et_2019_CancerCell.xlsx',sheet = "Sheet1"))
Sebest_RBP<-as.data.frame(read_xlsx('/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/RBP_gene_lists/Sebestyen_et_2016_GenRes.xlsx',sheet = "Sheet1"))
rMATs_RBP<-read.delim('/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/RBP_gene_lists/rMAPS2_list.txt')


# Genes id and ENS id
geneid<-read.delim('/Volumes/cancer/TCell/vast_output/genes_ENS.txt',header = TRUE,sep=' ')
#geneid<-read.delim('/Users/yguillen/Desktop/temp/TCell/vast_output/genes_ENS.txt',header = TRUE,sep=' ')

## For all samples
metadata<-read_xlsx("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/metadata_all.xlsx")
#metadata<-read_xlsx("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/Public_data/metadata_all.xlsx")

metadata_target<-as.data.frame(read_xlsx("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/metadata_gender_TARGET.xlsx"))
row.names(metadata_target)<-metadata_target$GSM
metadata_target$GSM<-NULL

#####

# Metadata summary
#####

metadata_sub<-subset(metadata,select=c(SampID,total_reads))
metadata_sub<-aggregate(x = metadata_sub$total_reads, by = list(metadata_sub$SampID), FUN = "sum")
metadata_sub$thres<-ifelse(metadata_sub$x>=65000000,c("validvast"),c("novalidvast"))

metadata<-merge(metadata,metadata_sub,by.x="SampID",by.y="Group.1")

metadata$Event_free_surv_days<-as.numeric(metadata$Event_free_surv_days)

metadata$Reads_pair1<-NULL
metadata$total_reads<-NULL
metadata<-metadata[!duplicated(metadata$SampID),]
names(metadata)[16]<-"Reads"

ggplot(metadata,aes(x=SampID,y=Reads))+
  geom_segment(aes(yend = Reads, x = SampID, xend = SampID, y = 0),color="grey",size=3)+
  geom_segment(data=metadata[metadata$First_event!="NA",],aes(yend = Reads, x = SampID, xend = SampID, y = 0,color=First_event),size=2)+
  geom_point(aes(y=Reads+2000000,shape=Stage),alpha=0.5,size=4)+
  scale_shape_manual(values=c(13, 1,9,8))+
  #geom_jitter(aes(color=Pheno,shape=thres),size=3)+
  #geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_hline(yintercept=65000000)+
  scale_color_brewer(palette = "Set3")+
  facet_grid(~GSEA,scales="free",space="free")+
  theme_light()+
  #ylim(c(0,300000000))+
  xlab("")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, hjust = 1),
        plot.title=element_text(family="Tahoma",angle=45,size=5),
        text=element_text(family="Tahoma"),
        legend.position = "bottom")

table(metadata$First_event,metadata$Vital_status)

#####


# CIBERSORT adding data
#####

## Expression profile CIBERSORT from above vastout
#cibersort<-read.delim("/Users/yguillen/Desktop/temp/IsoTALL_YGuillen/CIBERSORT/CIBERSORT.Output_Job5.txt",header = TRUE)
cibersort<-read.delim("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/CIBERSORT/CIBERSORT.Output_Job5.txt",header = TRUE)
row.names(cibersort)<-cibersort$Input.Sample
colnames(cibersort)[1]<-"SampID"

ciberbar<-cibersort[,c(1:(ncol(cibersort)-3))]
corciber<-merge(ciberbar,metadata,by="SampID")

corciber$Blasts<-as.numeric(corciber$Blasts)
cor.test(corciber$Blasts,corciber$Thy1)
plot(corciber$Blasts,corciber$Thy1)

ciberbar<-melt(ciberbar,id.vars="SampID")

ciberbar<-merge(metadata,ciberbar,by="SampID")

ciberbar$SampID <- factor(ciberbar$SampID, levels=unique(ciberbar[rev(order(ciberbar$value)), ]$SampID))

sampord<-ciberbar[ciberbar$variable=="Thy2",][rev(order(ciberbar[ciberbar$variable=="Thy2",]$value)),]$SampID
sampord
ciberbar$SampID <- factor(ciberbar$SampID, levels=sampord)

ciberbar$Blasts<-as.numeric(ciberbar$Blasts)

ggplot(ciberbar, aes(x=as.factor(SampID), y=value,fill=variable))+
  geom_bar(stat="identity", fill="black",width = 1) +
  geom_bar(stat="identity", aes(fill=variable),width = 0.8) +
  scale_fill_brewer(palette="Paired")+
  facet_grid(~First_event,scales="free",space="free")+
  #geom_point(aes(y=genes+500,color=Reads),size=5)+
  new_scale_color()+
  geom_point(data=ciberbar[!is.na(ciberbar$First_event) & ciberbar$variable=="Thy2",],aes(y=1.3),color="black",size=3,shape=17)+
  geom_point(data=ciberbar[!is.na(ciberbar$First_event) & ciberbar$variable=="Thy2",],aes(y=1.3,color=First_event),size=2,shape=17)+
  scale_color_brewer(palette="Set1")+
  new_scale_color()+
  geom_point(data=ciberbar[ciberbar$variable=="Thy2",],aes(y=1.1),color="black",size=3)+
  geom_point(data=ciberbar[ciberbar$GSEA!="line" & ciberbar$variable=="Thy2",],aes(y=1.1,color=Pheno),size=2)+
  scale_color_brewer(palette="Set3")+
  new_scale_color()+
  geom_point(data=ciberbar[ciberbar$variable=="Thy2",],aes(y=1.2),color="black",size=3)+
  geom_point(data=ciberbar[ciberbar$GSEA!="line" & ciberbar$variable=="Thy2",],aes(y=1.2,color=GSEA),size=2)+
  scale_color_brewer(palette="Dark2")+
  geom_point(data=ciberbar[ciberbar$GSEA!="line" & ciberbar$variable=="Thy2" & !is.na(ciberbar$Blasts),],aes(y=1.4,size=Blasts/10),color="grey")+
  #  scale_colour_gradient(low = "white", high = "black")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom")


# Summary of variable by group
tapply(ciberbar[!is.na(ciberbar$Blasts),]$Blasts, ciberbar[!is.na(ciberbar$Blasts),]$First_event, summary)

#####


## B-CAT CORRELATION NETWORWK MODELS
#####
###### CORRELATION MATRIX 

#Bcat gene = ENSG00000168036

## Correlation between one target genes and others!

## Dataframes categorized by T-ALL, or normal T-cell development (DSmat is not corrected)
# DO NOT INCLUDE TLE samples
# DSmat_noTLE comes from Clinical_TCell.R code

Tdev<-subset(DSmat_noTLE,DSmat_noTLE$GSEA=="GSE69239")
Tall<-subset(DSmat_noTLE,DSmat_noTLE$Source=="TALL")
Tlbl<-subset(DSmat_noTLE,DSmat_noTLE$Source=="TLBL")
Tleuk<-subset(DSmat_noTLE,DSmat_noTLE$Source=="TALL" | DSmat_noTLE$Source=="TLBL")
Ttarget<-subset(DSmat_noTLE,DSmat_noTLE$GSEA=="TARGET")

# Cell lines to predict later
Tline<-subset(DSmat_noTLE,DSmat_noTLE$GSEA=="line")

# In TALL and T-LBL databases, bcat coexpressed genes or kaiso or brca1
corout<-apply(Tleuk[,18:ncol(Tleuk)], 2, cor.test,Tleuk$ENSG00000168036,method="spearman")

#For only target
corout<-apply(Ttarget[,18:ncol(Ttarget)], 2, cor.test,Ttarget$ENSG00000168036,method="spearman")


corout1<-as.data.frame(sapply(corout, "[[", "p.value"))
corout2<-as.data.frame(sapply(corout, "[[", "estimate"))
corout<-cbind(corout1,corout2)
names(corout)<-c("Pval","Rho")
corout$Gene<-row.names(corout)
corout<-merge(geneid,corout,by="Gene",all.y=TRUE)


corout_sig<-subset(corout,corout$Rho>=0.4)
corout_sig<-corout_sig[with(corout_sig,order(corout_sig$Pval, corout_sig$Rho)),]

negcorctnnb1<-(corout_sig[corout_sig$Rho<0,])$NAME
length(negcorctnnb1)
poscorctnnb1<-(corout_sig[corout_sig$Rho>0,])$NAME
length(poscorctnnb1)

# vector for functional enrichment
(corout_sig[corout_sig$Rho>0.4,])$NAME

# KAISO for all
corout_kai<-apply(Tleuk[,18:ncol(Tleuk)], 2, cor.test,Tleuk$ENSG00000177485,method="spearman")

# KAISO for TARGET
corout_kai<-apply(Ttarget[,18:ncol(Ttarget)], 2, cor.test,Ttarget$ENSG00000177485,method="spearman")

corout1_kai<-as.data.frame(sapply(corout_kai, "[[", "p.value"))
corout2_kai<-as.data.frame(sapply(corout_kai, "[[", "estimate"))
corout_kai<-cbind(corout1_kai,corout2_kai)
names(corout_kai)<-c("Pval","Rho")
corout_kai$Gene<-row.names(corout_kai)
corout_kai<-merge(geneid,corout_kai,by="Gene",all.y=TRUE)


# BRCA1 for all
corout_brca<-apply(Tleuk[,18:ncol(Tleuk)], 2, cor.test,Tleuk$ENSG00000012048,method="spearman")

# BRCA1 for TARGET
corout_brca<-apply(Ttarget[,18:ncol(Ttarget)], 2, cor.test,Ttarget$ENSG00000012048,method="spearman")

corout1_brca<-as.data.frame(sapply(corout_brca, "[[", "p.value"))
corout2_brca<-as.data.frame(sapply(corout_brca, "[[", "estimate"))
corout_brca<-cbind(corout1_brca,corout2_brca)
names(corout_brca)<-c("Pval","Rho")
corout_brca$Gene<-row.names(corout_brca)
corout_brca<-merge(geneid,corout_brca,by="Gene",all.y=TRUE)


# Model correlation example RBM39
as.character((geneid[geneid$NAME=="RBM39",])$Gene)
as.character((geneid[geneid$NAME=="CTNNB1",])$Gene)

#model in Tall + tlbl
model_RBM39<-lm(ENSG00000131051~ENSG00000168036,data=Tleuk)

#observed CTNNB1 (constant)
obs_CTNNB1<-data.frame(
  ENSG00000168036=c(Tline[,colnames(Tline) == "ENSG00000168036"]))

row.names(obs_CTNNB1)<-row.names(Tline[,colnames(Tline) == "ENSG00000168036" | colnames(Tline) == "ENSG00000131051" ,])

# Observed x gene (specific for each gene)
obs_RBM39<-data.frame(
  ENSG00000131051=c(Tline[,colnames(Tline) == "ENSG00000131051" ,]))

row.names(obs_RBM39)<-row.names(Tline[,colnames(Tline) == "ENSG00000168036" | colnames(Tline) == "ENSG00000131051" ,])

# prediction
pre_RBM39<-predict(lm(ENSG00000131051~ENSG00000168036,data=Tleuk),
                   newdata = obs_CTNNB1,
                   interval="prediction")
pre_RBM39


# Save data models for each correlated gene and prediction
#Select matrix with expression of correlated genes with bcat
modcor<-Tleuk[,colnames(Tleuk) %in% c(as.character(corout_sig$Gene),"ENSG00000168036")]
dim(modcor)

# Output matrix with expected expression based on model
exp_fit<-data.frame(t(apply(modcor, 2, function(x) predict(lm(x~ENSG00000168036,data=modcor),
                                                           newdata = obs_CTNNB1,
                                                           interval="prediction")[,1])))

colnames(exp_fit)<-gsub('^',"expected_",colnames(exp_fit))

exp_lwr<-data.frame(t(apply(modcor, 2, function(x) predict(lm(x~ENSG00000168036,data=modcor),
                                                           newdata = obs_CTNNB1,
                                                           interval="prediction")[,2])))
colnames(exp_lwr)<-gsub('^',"lwr_",colnames(exp_lwr))

exp_upr<-data.frame(t(apply(modcor, 2, function(x) predict(lm(x~ENSG00000168036,data=modcor),
                                                           newdata = obs_CTNNB1,
                                                           interval="prediction")[,3])))
colnames(exp_upr)<-gsub('^',"upr_",colnames(exp_upr))

# Observed x gene (specific for each gene)
obs_gene<-data.frame(t(Tline[,colnames(Tline) %in% c(as.character(corout_sig$Gene),"ENSG00000168036")]))
colnames(obs_gene)<-gsub('^','observed_',colnames(obs_gene))

# bind observed and expected
est_mod<-cbind(exp_fit,exp_lwr,exp_upr,obs_gene)

est_mod$as_cont_jurk<-ifelse(est_mod$lwr_control_Jurkat < est_mod$observed_control_Jurkat &  est_mod$observed_control_Jurkat < est_mod$upr_control_Jurkat, 
                             c("YES"),c("NO")) 

est_mod$as_cont_RPMI<-ifelse(est_mod$lwr_control_RPMI < est_mod$observed_control_RPMI &  est_mod$observed_control_RPMI < est_mod$upr_control_RPMI, 
                             c("YES"),c("NO")) 

est_mod$as_sh_jurk<-ifelse(est_mod$lwr_shbcat_Jurkat < est_mod$observed_shbcat_Jurkat &  est_mod$observed_shbcat_Jurkat < est_mod$upr_shbcat_Jurkat, 
                           c("YES"),c("NO")) 

est_mod$as_sh_RPMI<-ifelse(est_mod$lwr_shbcat_RPMI < est_mod$observed_shbcat_RPMI &  est_mod$observed_shbcat_RPMI < est_mod$upr_shbcat_RPMI, 
                           c("YES"),c("NO")) 


est_mod<-merge(geneid,est_mod,by.x="Gene",by.y="row.names")

coreg_bcat<-est_mod[est_mod$as_sh_jurk =="YES" & est_mod$as_sh_RPMI=="YES" & est_mod$as_cont_jurk=="YES" & est_mod$as_cont_RPMI=="YES",]
nocont<-est_mod[est_mod$as_cont_jurk =="NO" & est_mod$as_cont_RPMI=="NO" & est_mod$as_sh_jurk =="NO" & est_mod$as_sh_RPMI=="NO",]

coreg_bcat$NAME

gene<-as.character((geneid[geneid$NAME=="ZNF76",])$Gene)
bcat<-as.character((geneid[geneid$NAME=="CTNNB1",])$Gene)

ggplot(rbind(Tleuk,Tline,Tdev),aes_string(x=bcat,y=gene))+
  geom_point(size=5)+
  geom_point(data=Tline,aes(color=Pheno,shape=Sample),size=7)+
  scale_color_brewer("Spectral")+
  new_scale_color()+
  geom_point(aes(color=Type),size=3)+
  scale_color_brewer(palette="Paired")+
  geom_smooth(data=rbind(Tleuk,Tdev),method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,aes(color=GSEA))+
  #geom_smooth(data=DSmat[DSmat$Source!="TALL" ,],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="purple")+
  #geom_smooth(data=gene_sam[gene_sam$RBPJ<75 & gene_sam$Source!="TALL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="red")+
  #geom_label_repel(data=DSmat[DSmat$GSEA=="line" ,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  facet_wrap(~GSEA,ncol=2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))
#coord_fixed()


# Plots
gene<-as.character((geneid[geneid$NAME=="BRCA1",])$Gene)
bcat<-as.character((geneid[geneid$NAME=="CTNNB1",])$Gene)
kaiso<-as.character((geneid[geneid$NAME=="ZBTB3",])$Gene)

DSmat_noTLE$Source<-factor(DSmat_noTLE$Source,levels = c("HSC","LMPP","CLP","BCP","Thy1","Thy2","DN","DP","CD4","CD8","Thymus","TALL","TLBL","Jurkat","RPMI"))

ggplot(DSmat_noTLE[DSmat_noTLE$GSEA=="TARGET",],aes_string(x=bcat,y=kaiso))+
  geom_point(size=5)+
  geom_point(aes(color=First_event),size=4)+
  scale_color_brewer(palette="Paired")+
  new_scale_color()+
  #  geom_point(data=DSmat_noTLE[DSmat_noTLE$GSEA=="GSE69239",],size=3.5)+
  #  geom_point(data=DSmat_noTLE[DSmat_noTLE$GSEA=="GSE69239",],size=2)+
  #  geom_smooth(data=DSmat_noTLE[DSmat_noTLE$Source=="TALL" | DSmat_noTLE$Source=="TLBL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,aes(color=GSEA))+
  #geom_smooth(data=DSmat[DSmat$Source!="TALL" ,],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="purple")+
  #geom_smooth(data=gene_sam[gene_sam$RBPJ<75 & gene_sam$Source!="TALL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="red")+
  #  geom_label_repel(data=DSmat_noTLE[DSmat_noTLE$GSEA=="line" ,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  geom_smooth(method="lm")+
  facet_wrap(~GSEA,ncol=2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))
#coord_fixed()

cor.test(DSmat_noTLE[DSmat_noTLE$GSEA=="TARGET",]$ENSG00000185670,DSmat_noTLE[DSmat_noTLE$GSEA=="TARGET",]$ENSG00000168036)

#Plot only normal T-cell cohort
ggplot(DSmat_noTLE[DSmat_noTLE$Source!="TALL" & DSmat_noTLE$Source!="TLBL" & DSmat_noTLE$GSEA!="line",],aes_string(x=gene,y=kaiso))+
  geom_point(size=7)+
  geom_point(aes(color=Type),size=6)+
  scale_color_brewer(palette="Paired")+
  geom_smooth(data=DSmat_noTLE[DSmat_noTLE$Source!="TALL" & DSmat_noTLE$Source!="TLBL" & DSmat_noTLE$GSEA!="line",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="blue")+
  #geom_smooth(data=gene_sam[gene_sam$RBPJ<75 & gene_sam$Source!="TALL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="red")+
  #geom_label_repel(data=gene_sam[gene_sam$DDX3X>350,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

gene
bcat

#Plot only T-ALL cohort
ggplot(DSmat_noTLE[DSmat_noTLE$Source=="TALL" | DSmat_noTLE$Source=="TLBL" | DSmat_noTLE$GSEA=="line",],aes(x=log(ENSG00000012048),y=log(ENSG00000185670)))+
  geom_point(size=7)+
  geom_point(aes(color=Type),size=6)+
  scale_color_brewer(palette="Paired")+
  geom_smooth(data=DSmat_noTLE[DSmat_noTLE$Source=="TALL" | DSmat_noTLE$Source=="TLBL" | DSmat_noTLE$GSEA=="line",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="purple")+
  #geom_smooth(data=gene_sam[gene_sam$RBPJ<75 & gene_sam$Source!="TALL",],method="lm",se = FALSE,linetype="dotted",fullrange=TRUE,color="red")+
  #geom_label_repel(data=gene_sam[gene_sam$DDX3X>350,],aes(label=Pheno), size=3,fontface = "italic",alpha=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

#####


# list of DEGs shbcat in RPMI from chip.R script cross with b-cat chips and correlation models (INCLUDES VIOLIN PLOT FOR RNASEQ AND CHIPS)
#####

res_df_rna_rpmi<-read.delim('/Users/yguillen/Desktop/temp/Gekas_RNASeq/DESeq2/RPMI/res_df_rna_rpmi.txt')
res_df_rna_rpmi

res_df_rna_jurkat<-read.delim('/Users/yguillen/Desktop/temp/Gekas_RNASeq/DESeq2/Jurkat/res_df_rna_jurkat.txt')
res_df_rna_jurkat

# Kaiso chips
kaiso_peaks<-read_xlsx("/Users/yguillen/Desktop/temp/beta_catenin_project/KAISO_ChIPs/kaiso_in4.xlsx")
kaiso_peaks_tab<-kaiso_peaks[,colnames(kaiso_peaks) %in% c("annotation","ENSEMBL")]
kaiso_peaks_tab<-kaiso_peaks_tab[!duplicated(kaiso_peaks_tab),]
kaiso_peaks_tab$annotation<-gsub('Intron\\s(.*)','Intron',kaiso_peaks_tab$annotation)
kaiso_peaks_tab$annotation<-gsub('Exon\\s(.*)','Exon',kaiso_peaks_tab$annotation)
kaiso_target<-data.frame(table(kaiso_peaks_tab$ENSEMBL))
row.names(kaiso_target)<-kaiso_target$Var1
kaiso_target$Var1<-NULL
colnames(kaiso_target)<-("kaiso")



# b-cat chips
# genes with peaks in control 
# RD antibody
bcat_cpeaks<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_DBmerge_nearest.genes.txt",header = TRUE)
dim(bcat_cpeaks)
length(unique(bcat_cpeaks$gene))
length(unique(bcat_cpeaks$Name))

# SC antibody
#bcat_cpeaks<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_RBCSmerge_nearest.genes.txt",header = TRUE)
#dim(bcat_cpeaks)
#length(unique(bcat_cpeaks$Name))

# if a peak has two annotated genes, they are all taken into account
bcat_cpeaks_ens<-getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                       values= unique(bcat_cpeaks$gene),
                       filter="external_gene_name",
                       mart=hg) 

length(unique(bcat_cpeaks_ens$external_gene_name))

bcat_cpeaks<-merge(bcat_cpeaks_ens,bcat_cpeaks,by.x="external_gene_name",by.y="gene",all.y=TRUE)
dim(bcat_cpeaks)
length(unique(bcat_cpeaks$Name))
length(unique(bcat_cpeaks$external_gene_name))

# genes with peaks in lithium
#RD Antibody
bcat_lpeaks<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_DBLmerge_nearest.genes.txt",header = TRUE)
dim(bcat_lpeaks)
length(unique(bcat_lpeaks$Name))

#SC Antibody
#bcat_lpeaks<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_RBLSmerge_nearest.genes.txt",header = TRUE)


bcat_lpeaks_ens<-getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                       values= unique(bcat_lpeaks$gene),
                       filter="external_gene_name",
                       mart=hg) 

bcat_lpeaks<-merge(bcat_lpeaks_ens,bcat_lpeaks,by.x="external_gene_name",by.y="gene",all.y=TRUE)
dim(bcat_lpeaks)
length(unique(bcat_lpeaks$external_gene_name))

# genes with peaks in TCF 
tcf_peaks_dt<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_TCFmerge_nearest.genes.txt",header = TRUE)

tcf_peaks<-getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                 values= unique(tcf_peaks_dt$gene),
                 filter="external_gene_name",
                 mart=hg)

tcf_peaks<-merge(tcf_peaks_dt,tcf_peaks,by.y="external_gene_name",by.x="gene",all.x = TRUE)

length(unique(tcf_peaks$gene))

tcf_peaks_tab<-data.frame(TCF=unique(tcf_peaks[,colnames(tcf_peaks) %in% c("ensembl_gene_id")]))
tcf_target<-data.frame(table(tcf_peaks_tab$TCF))
row.names(tcf_target)<-tcf_target$Var1
tcf_target$Var1<-NULL
colnames(tcf_target)<-("TCF")


# genes with peaks in LEF
lef_peaks_dt<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/bed_merge_min_1_rep_decont/R_decont_LEFmerge_nearest.genes.txt",header = TRUE)
lef_peaks<-getBM(attributes=c("ensembl_gene_id", "description","external_gene_name","go_id"),
                 values= unique(lef_peaks_dt$gene),
                 filter="external_gene_name",
                 mart=hg)

lef_peaks<-merge(lef_peaks_dt,lef_peaks,by.y="external_gene_name",by.x="gene",all.x = TRUE)

lef_peaks_tab<-data.frame(LEF=unique(lef_peaks[,colnames(lef_peaks) %in% c("ensembl_gene_id")]))
lef_target<-data.frame(table(lef_peaks_tab$LEF))
row.names(lef_target)<-lef_target$Var1
lef_target$Var1<-NULL
colnames(lef_target)<-("LEF")

# RPMI
bcat_allpeaks<-rbind(bcat_cpeaks,bcat_lpeaks)
dim(bcat_allpeaks)
length(unique(bcat_allpeaks$Name))
length(unique(bcat_allpeaks$external_gene_name))

# filter genes which peak is found in > replicates
rep_bcat<-read.delim("/Volumes/cancer/bcat_Project/peakcall/bed_inter/unique_repl_RD.txt",header = FALSE)
colnames(rep_bcat)<-c("chr","Start","End","Name","Score")

bcat_allpeaks<-bcat_allpeaks[!(bcat_allpeaks$Name %in% rep_bcat$Name),]

## Annotate selection of peaks >2 reps (one of them RD)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb_human <- TxDb.Hsapiens.UCSC.hg38.knownGene

annot_bcat_allpeaks<-as.data.frame(bcat_allpeaks[,5:ncol(bcat_allpeaks)])
annot_bcat_allpeaks<-makeGRangesFromDataFrame(annot_bcat_allpeaks)
annot_bcat_allpeaks<-annotatePeak(annot_bcat_allpeaks,tssRegion = c(-3000, 3000), TxDb=txdb_human, annoDb="org.Hs.eg.db")

annostat_min2rep<-annot_bcat_allpeaks@annoStat
annostat_min2rep$group<-"peaks >=2 replicates"

ggplot(annostat_min2rep,aes(x=group,y=Frequency,fill=Feature))+
  geom_bar(stat="identity",color="black") +
  scale_fill_brewer(palette="Paired")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face = "bold"))+
  xlab("")+ylab("Percentage of peaks (%)")+
  coord_flip()

## functional renrichment enrichr output web genes min 2 reps
bp_genes_min2rep<-read.delim('/Volumes/grcmc/YGUILLEN/Projects_bioinfo_data_Yolanda/bcat/ChIPs_bcat_peaks_selection/BP_function_enrichment_peaks_min_2reps.txt',header = TRUE)
bpsub<-subset(bp_genes_min2rep,bp_genes_min2rep$Adjusted.P.value<=0.05)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$Adjusted.P.value))])

library(tidyr)
#For wiki pathways
bpsub<-separate(data = bpsub, col = Overlap, into = c("counts", "pathway"), sep = "/")
bpsub<-separate(data = bpsub, col = Term, into = c("Term", "GO"), sep = "GO")
bpsub$GO<-gsub(':','',bpsub$GO)
bpsub$GO<-gsub(')','',bpsub$GO)
bpsub$Term<-gsub(' \\(','',bpsub$Term)

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))
bpsub$Term <- factor(bpsub$Term, levels=bpsub[rev(order(bpsub$Adjusted.P.value)), ]$Term)

ggplot(bpsub[1:20,],aes(y=-log10(P.value),x=Term),alpha = 0.5)+
  geom_bar(position="dodge",stat = "identity",aes(fill=Combined.Score))+
  #facet_wrap(~group,scales="free")+
  theme_minimal()+
  coord_flip()+
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=8, face="bold",hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1),
        legend.position = "bottom")+
  ggtitle("Top 20 Biological Processes")

ggplot(bpsub[1:20,], aes(x = Term, y = -log10(P.value))) +
  geom_segment(aes(x = Term, xend = Term, y = 0, yend = -log10(P.value)), 
    color = "lightblue",size=3) + 
  geom_point(color="darkblue", aes(size = Combined.Score))+
  theme_void() +
  coord_flip()+
  theme(axis.title=element_text(size=10,face="bold"),
        axis.text.x = element_text(size=10, face="bold",hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1),
        legend.position = "bottom")+
  ggtitle("Top 20 Enriched Biological Processes")


## merge with OLD CHRISTOS GEKAS diff expression analysis

rpmi_al<-merge(bcat_allpeaks[!duplicated(bcat_allpeaks$ensembl_gene_id),],res_df_rna_rpmi,by.x="ensembl_gene_id",by.y="id",all=TRUE)
rpmi_al<-merge(rpmi_al,est_mod,by.x="ensembl_gene_id",by.y="Gene",all=TRUE)
rpmi_al<-merge(rpmi_al,corout,by.x="ensembl_gene_id",by.y="Gene",all=TRUE)

#variables
rpmi_al$peak<-ifelse((!is.na(rpmi_al$external_gene_name.x)),c("peak"),c("nopeak"))
rpmi_al$deg<-ifelse(!is.na(rpmi_al$pvalue) & rpmi_al$pvalue<0.05,c("DEG"),c("noDEG"))
rpmi_al$cor<-ifelse(!is.na(rpmi_al$Rho) & abs(rpmi_al$Rho)>=0.4 & rpmi_al$Pval<0.01, c("cor"),c("noncor"))
rpmi_al$fit<-ifelse(!is.na(rpmi_al$as_cont_RPMI) & rpmi_al$as_cont_RPMI=="YES" & rpmi_al$as_sh_RPMI=="YES",c("fit"),c("nofit"))

rpmi_al$comb<-paste(rpmi_al$peak,rpmi_al$deg,sep="_")
rpmi_al$comb<-paste(rpmi_al$comb,rpmi_al$cor,sep="_")
rpmi_al$comb<-paste(rpmi_al$comb,rpmi_al$fit,sep="_")

row.names(rpmi_al)<-rpmi_al$ensembl_gene_id

table(rpmi_al$comb)

## number of genes per each category
clas_gene<-data.frame(class=c("Peaks",
                              "DEGs",
                              "coexpression",
                              "fit_model"),
                      frequency=c(dim(rpmi_al[rpmi_al$peak=="peak",])[1],
                                  dim(rpmi_al[rpmi_al$deg=="DEG",])[1],
                                  dim(rpmi_al[rpmi_al$cor=="cor",])[1],
                                  dim(rpmi_al[rpmi_al$fit=="fit",])[1]))

clas_gene$class<- factor(clas_gene$class, levels = clas_gene$class[order(clas_gene$frequency, decreasing = TRUE)])


ggplot(clas_gene, aes(x=class, y=frequency)) +
  geom_bar(stat="identity", color="black",fill="orange",alpha=0.5)+
  geom_text(aes(label=frequency), vjust=-0.1, color="black", size=10)+
  theme_minimal()+
  theme(axis.text.x=element_text(size=14,angle=90),
        axis.text.y=element_text(size=14),
        axis.title=element_blank())



bardata<-data.frame(table(rpmi_al$comb))
colnames(bardata)<-c("combination","Freq")

bardata$combination<- factor(bardata$combination, levels = bardata$combination[order(bardata$Freq, decreasing = FALSE)])

ggplot(bardata[bardata$combination!="nopeak_noDEG_noncor_nofit",], aes(x=combination, y=Freq)) +
  geom_bar(stat="identity", color="black",fill="#56B4E9",alpha=0.5)+
  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=6)+
  theme_minimal()+
  theme(axis.text.x=element_text(size=10,angle=90,face="bold.italic"),
        axis.text.y=element_text(size=14,face="bold.italic"),
        axis.title=element_blank())

rpmi_al[rpmi_al$comb=="peak_DEG_cor_fit",]$external_gene_name.x
rpmi_al[rpmi_al$comb=="peak_DEG_cor_nofit",]$external_gene_name.x


## VIOLIN PLOT WITH BCAT,TCF AND LEF PEAKS AND RNASEQ FC shBCAT
unique(tcf_peaks$ensembl_gene_id)

tcf_peaksdup<-tcf_peaks[!duplicated(tcf_peaks$ensembl_gene_id), ]
#tcf_peaksdup<-tcf_peaks[!duplicated(tcf_peaks$gene), ]

tcf_peaksdup<-tcf_peaksdup[,c(1,6,12)]
tcf_peaksdup$TCF<-"TCF"

lef_peaksdup<-lef_peaks[!duplicated(lef_peaks$ensembl_gene_id), ]
#lef_peaksdup<-lef_peaks[!duplicated(lef_peaks$gene), ]

lef_peaksdup<-lef_peaksdup[,c(1,6,12)]
lef_peaksdup$LEF<-"LEF"


#rpmi_al<-merge(rpmi_al,tcf_peaksdup,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.x="TRUE")
#rpmi_al<-merge(rpmi_al,lef_peaksdup,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.x="TRUE")

#rpmi_al$peaks<-paste(rpmi_al$peak,rpmi_al$TCF,sep="_")
#rpmi_al$peaks<-paste(rpmi_al$peaks,rpmi_al$LEF,sep="_")

#ggplot(rpmi_al,aes(x=peaks,y=log2FoldChange,label=gene.x))+
#  geom_jitter(data=rpmi_al[rpmi_al$padj>0.1 & !is.na(rpmi_al$padj),],color="grey")+
#  geom_jitter(data=rpmi_al[rpmi_al$padj<=0.1 & !is.na(rpmi_al$padj),],color="red")+
#  #geom_label_repel(data=rpmi_al[rpmi_al$padj<=0.05 & !is.na(rpmi_al$padj),],aes(label=gene.y),size=2)+
#  facet_grid(~peak,scales= "free_x")+
#  geom_violin(data=rpmi_al[rpmi_al$padj>0.1 & !is.na(rpmi_al$padj),],color="grey",alpha=0.3)+
#  geom_violin(data=rpmi_al[rpmi_al$padj<=0.1 & !is.na(rpmi_al$padj),],color="red",alpha=0.3)+
#  theme_bw()


#####


## b-cat signature
#####

## Genes with peak b-cat

#rpmi_al has all the data peaks, DEGs, cor, fit
rpmi_al
length(unique(rpmi_al$external_gene_name.y))
length(unique(rpmi_al$ensembl_gene_id))

# output for ssGSEA all data
ssGSEA<-as.data.frame(t(merge(metadata,vastout_t,by="row.names")))
colnames(ssGSEA)<-ssGSEA[1,]
ssGSEA<-ssGSEA[19:nrow(ssGSEA),]
ssGSEA<-merge(geneid,ssGSEA,by.x="Gene",by.y="row.names")
ssGSEA$Gene<-NULL
ssGSEA$Description<-'NA'
ssGSEA<-ssGSEA[,c(1,ncol(ssGSEA),2:(ncol(ssGSEA)-1))]
dim(ssGSEA)
write.table(ssGSEA,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/ggsea_bcatexp_allgenes_incTLE.gct",sep="\t",row.names=FALSE,quote = FALSE)

#only b-cat targets, from clusters unsupervised microarrays all b-cat targets
ssGSEA_bcat<-ssGSEA[ssGSEA$NAME %in% clusters$Gene,]
write.table(ssGSEA_bcat,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/unsupervised_microarrays_to_TARGET/ggsea_bcatexp_targets.gct",sep="\t",row.names=FALSE,quote = FALSE)



# All Samples
# excluding genes with no conditions
bcatexp<-vastout_t[,colnames(vastout_t) %in% rpmi_al[rpmi_al$comb!="nopeak_noDEG_noncor_nofit",]$ensembl_gene_id,]

#including genes with condition peak
bcatexp<-vastout_t[,colnames(vastout_t) %in% rpmi_al[rpmi_al$peak=="peak",]$ensembl_gene_id,]


dim(bcatexp)
colcol<-(rpmi_al[rpmi_al$ensembl_gene_id %in% colnames(bcatexp),c(6,11,9,16,49,50,52:54)])
row.names(colcol)<-colnames(bcatexp)
dim(colcol)
colcol$dss<-colcol$Start-colcol$genestart

#merge with TARGET correlation networwk, not all
row.names(corout)<-corout$Gene
colnames(corout)<-c("Gene","NAME","Pval_target","Rho_target")
colcol<-merge(colcol,corout,by="row.names")
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL
colcol$cortarget<-ifelse(!is.na(colcol$Pval_target) & abs(colcol$Rho_target)>=0.4 & colcol$Pval_target<0.01, c("cortar"),c("noncortar"))

colcol$corsen<-ifelse(colcol$Rho_target>0,'pos','neg')
colcol$corsen<-paste(colcol$cortarget,colcol$corsen,sep="_")


# KAISO, TCF and LEF targets

row.names(colcol)
row.names(kaiso_target)
dim(colcol)

colcol<-merge(colcol,kaiso_target,by="row.names",all.x=TRUE)
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL
colcol<-merge(colcol,tcf_target,by="row.names",all.x=TRUE)
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL
colcol<-merge(colcol,lef_target,by="row.names",all.x=TRUE)
row.names(colcol)<-colcol$Row.names
colcol$Row.names<-NULL

colcol$kaiso[is.na(colcol$kaiso)]<-0
colcol$kaiso[colcol$kaiso>0]<-1
colcol$kaiso<-as.factor(colcol$kaiso)

colcol$TCF[is.na(colcol$TCF)]<-0
colcol$TCF[colcol$TCF>0]<-1
colcol$TCF<-as.factor(colcol$TCF)

colcol$LEF[is.na(colcol$LEF)]<-0
colcol$LEF[colcol$LEF>0]<-1
colcol$LEF<-as.factor(colcol$LEF)

# markers matrix with no NAs
complete.cases(bcatexp)
is.na(bcatexp)
bcatexp<-bcatexp[complete.cases(bcatexp),]
dim(bcatexp)

length(colnames(bcatexp))
length((geneid[geneid$Gene %in% colnames(bcatexp),])$NAME)

# targets with no expression
noexp<-colnames(bcatexp[, which(numcolwise(sum)(bcatexp) ==0)])
noexp

# targets with expressed,Remove genes with sum 0 in all samples
exp<-colnames(bcatexp[, which(numcolwise(sum)(bcatexp) !=0)])
exp

bcatexp<-bcatexp[,colnames(bcatexp) %in% exp]

row.names(metadata)<-metadata$SampID
bcatexp<-merge(subset(metadata,select=c(Source,Pheno,Type,GSEA,First_event,Event_free_surv_days,Vital_status,Lesion,Age_Dx_days)),bcatexp,by="row.names",sort=FALSE)
dim(metadata)
dim(bcatexp)

## Add mutations metadata info
somatic_target_bcat<-somatic_target
somatic_target_bcat$Tumor_Sample_Barcode<-gsub('-0.*','',somatic_target_bcat$Tumor_Sample_Barcode)
somatic_target_bcat<-merge(somatic_target_bcat,metadata,by.x="Tumor_Sample_Barcode",by.y="Sample")

somatic_target_bcat_notch<-as.data.frame(table(somatic_target_bcat[grepl('NOTCH1',somatic_target_bcat$Hugo_Symbol),]$SampID))
row.names(somatic_target_bcat_notch)<-somatic_target_bcat_notch$Var1
colnames(somatic_target_bcat_notch)<-c('Samplemut','NOTCH1_freq')

somatic_target_bcat_fbxw7<-as.data.frame(table(somatic_target_bcat[grepl('FBXW7',somatic_target_bcat$Hugo_Symbol),]$SampID))
row.names(somatic_target_bcat_fbxw7)<-somatic_target_bcat_fbxw7$Var1
colnames(somatic_target_bcat_fbxw7)<-c('Samplemut','FBXW7_freq')

somatic_target_bcat_lef<-as.data.frame(table(somatic_target_bcat[grepl('LEF1',somatic_target_bcat$Hugo_Symbol),]$SampID))
row.names(somatic_target_bcat_lef)<-somatic_target_bcat_lef$Var1
colnames(somatic_target_bcat_lef)<-c('Samplemut','LEF1_freq')

somatic_target_mut<-merge(somatic_target_bcat_notch,somatic_target_bcat_fbxw7,by="Samplemut",all=TRUE)
somatic_target_mut<-merge(somatic_target_mut,somatic_target_bcat_lef,by="Samplemut",all=TRUE)
row.names(somatic_target_mut)<-somatic_target_mut$Samplemut

row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL


somatic_target_mut$NOTCH1_freq[is.na(somatic_target_mut$NOTCH1_freq)] = "0"
somatic_target_mut$FBXW7_freq[is.na(somatic_target_mut$FBXW7_freq)] = "0"
somatic_target_mut$LEF1_freq[is.na(somatic_target_mut$LEF1_freq)] = "0"


somatic_target_mut$NOTCH1_freq<-as.numeric(somatic_target_mut$NOTCH1_freq)
somatic_target_mut$LEF1_freq<-as.numeric(somatic_target_mut$LEF1_freq)
somatic_target_mut$FBXW7_freq<-as.numeric(somatic_target_mut$FBXW7_freq)

bcatexp<-merge(somatic_target_mut,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)

row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

bcatexp$LEF1_freq<-as.factor(bcatexp$LEF1_freq)

### Add bcat expression variable
bcat_lev<-data.frame(row.names = row.names(vastout_t_TLE),
                     bcat = vastout_t_TLE[,colnames(vastout_t_TLE) == "ENSG00000168036"])

bcatexp<-merge(bcat_lev,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

### Add kaiso expression variable
kaiso_lev<-data.frame(row.names = row.names(vastout_t_TLE),
                      kaiso_exp = vastout_t_TLE[,colnames(vastout_t_TLE) == "ENSG00000177485"])

bcatexp<-merge(kaiso_lev,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

### Add MYC expression variable
myc_lev<-data.frame(row.names = row.names(vastout_t_TLE),
                    myc_exp = vastout_t_TLE[,colnames(vastout_t_TLE) == "ENSG00000136997"])

bcatexp<-merge(myc_lev,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL


### Add BRCA1 expression variable
brca1_lev<-data.frame(row.names = row.names(vastout_t_TLE),
                      brca1_exp = vastout_t_TLE[,colnames(vastout_t_TLE) == "ENSG00000012048"])

bcatexp<-merge(brca1_lev,bcatexp,by="row.names",sort=FALSE,all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

# add gender info and overall survival
bcatexp<-merge(metadata_target[,c(2:5)],bcatexp,by="row.names",all.y=TRUE)
row.names(bcatexp)<-bcatexp$Row.names
bcatexp$Row.names<-NULL

#write.table(bcatexp,"/Users/yguillen/Downloads/bcat_RD_target_genes_expression_TALL.txt",quote = FALSE,row.names = FALSE,sep="\t")

### Using all patients, no matter the first event status
library(RColorBrewer)
cols <- colorRampPalette(rev(brewer.pal(10,name="RdBu")))(50)

paletteLength<-length(cols)

#cols <- colorRampPalette(c("red", "green", "blue"))(50)

## Add or not GSE69239 ( | bcatexp$GSEA=="GSE69239")

myBreaks <- c(seq(min(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])]))),
               annotation_col = as.data.frame(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,c(1,2,3,5,7,8,10,11,12,14,17,19)]),
               annotation_row = colcol[,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 3,
               cutree_rows = 10,
               color=cols, breaks=myBreaks)

clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=10)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=3)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==4,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,22:ncol(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",])])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(bcatexp,patients,by="row.names")
ggplot(patinfo,aes(x=as.factor(cluster),y=bcat))+
  geom_boxplot(aes(color=as.factor(cluster)))+
  theme_bw()

# boxplot expression bcat, tcf, lef and kaiso and brca1 clusters patients
clust_bcat<-data.frame(cluster="cluster1",
                       exp_bcat=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000168036,
                       exp_kaiso=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000177485,
                       exp_tcf=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000081059,
                       exp_lef=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000138795,
                       exp_brca1=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==1)),]$ENSG00000012048)

clust_bcat<-rbind(clust_bcat,
                  data.frame(cluster="cluster2",
                             exp_bcat=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000168036,
                             exp_kaiso=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000177485,
                             exp_tcf=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000081059,
                             exp_lef=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000138795,
                             exp_brca1=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==2)),]$ENSG00000012048))

clust_bcat<-rbind(clust_bcat,
                  data.frame(cluster="cluster3",
                             exp_bcat=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000168036,
                             exp_kaiso=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000177485,
                             exp_tcf=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000081059,
                             exp_lef=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000138795,
                             exp_brca1=vastout_t[row.names(vastout_t) %in% row.names(subset(patients,patients$cluster==3)),]$ENSG00000012048))

clust_bcat<-melt(clust_bcat)

ggplot(clust_bcat,aes(x=cluster,y=value,fill=value))+
  geom_boxplot(aes(fill=cluster),alpha=0.3)+
  facet_wrap(~variable,scales="free",nrow = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle=45,hjust = 1),
        axis.text.y = element_text(size=14),
        strip.text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set3")


write.table(clusters,"/Users/yguillen/Downloads/clusters_TARGET_bcatpeaks.txt",sep="\t",quote = FALSE,row.names = FALSE)


## Create correlation across genes using only expression data in patients with no NAs
bcatexp_survival <-bcatexp[bcatexp$GSEA=="TARGET",]
table(bcatexp_survival$First_event)
dim(bcatexp_survival)

#Discard NAs and Censored
bcatexp_survival<-subset(bcatexp_survival,bcatexp_survival$First_event!="NA")
bcatexp_survival<-subset(bcatexp_survival,bcatexp_survival$First_event!="Censored")
names<-row.names(bcatexp_survival)
length(names)

# Remission = 0, no remission = 1
bcatexp_survival<-bcatexp_survival %>% mutate(disease_free =  ifelse(First_event == "Remission", 0, 1))
bcatexp_survival$disease_free<-as.factor(bcatexp_survival$disease_free)
dim(bcatexp_survival)
str(bcatexp_survival)

table(bcatexp_survival$disease_free)
row.names(bcatexp_survival)<-names

# Create variable brca1 expression
bcatexp_survival$brcagroup<-NA
bcatexp_survival$brcagroup<-ifelse(bcatexp_survival$brca1_exp>median(bcatexp_survival$brca1_exp),
                                   "high_brca1",
                                   "low_brca1")


# Create variable kaiso expression
bcatexp_survival$kaisogroup<-NA
bcatexp_survival$kaisogroup<-ifelse(bcatexp_survival$kaiso_exp>median(bcatexp_survival$kaiso_exp),
                                    "high_kaiso",
                                    "low_kaiso")

# Create variable bcat expression
bcatexp_survival$bcatgroup<-NA
bcatexp_survival$bcatgroup<-ifelse(bcatexp_survival$bcat>median(bcatexp_survival$bcat),
                                   "high_bcat",
                                   "low_bcat")


library(RColorBrewer)
cols <- colorRampPalette(rev(brewer.pal(10,name="RdBu")))(50)

paletteLength<-length(cols)

bcatexp_survival$Age_Dx_days<-as.numeric(bcatexp_survival$Age_Dx_days)
bcatexp_survival$Age_Dx_days<-bcatexp_survival$Age_Dx_days/365
colnames(bcatexp_survival)[21]<-c("Age_Dx_years")

#cols <- colorRampPalette(c("red", "green", "blue"))(50)
myBreaks <- c(seq(min(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)]))),
               annotation_col = as.data.frame(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",c(7,8,14,17,19)]),
               annotation_row = colcol[,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 5,
               cutree_rows = 7,
               color=cols, breaks=myBreaks)

clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=7)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=5)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==2,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",22:(ncol(bcatexp_survival)-4)])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(bcatexp_survival,patients,by="row.names")

patinfo$group<-patinfo$cluster
patinfo$group<-gsub('4','1',patinfo$group)
patinfo$group<-gsub('5','3',patinfo$group)

ggplot(patinfo[patinfo$First_event!="Second Malignant Neoplasm",],aes(x=reorder(as.factor(group),bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3)+
  geom_boxplot(alpha=0.2,lwd=1)+
  #  geom_hline(yintercept = median(scale(bcatexp_survival$brca1_exp))+0.1,color="grey",lwd=1,linetype = "dotted")+
  #  geom_hline(yintercept = median(scale(bcatexp_survival$brca1_exp))-0.1,color="grey",lwd=1,linetype = "dotted")+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[4],color="grey",lwd=2)+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))

patinfo$scaled_brca1<-as.numeric(scale(patinfo$brca1_exp))
median(patinfo$scaled_brca1)
median(patinfo[patinfo$cluster=="1",]$scaled_brca1)

## survival curves with all patients, REDO PATIENTS WITH HEATMAP

bcatexp_survival_cluster<-merge(patients,bcatexp_survival,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Row.names<-NULL


library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)
bcatexp_survival_cluster$Vital_status<-ifelse(bcatexp_survival_cluster$Vital_status=="Alive",0,1)
bcatexp_survival_cluster$Vital_status<-as.integer(bcatexp_survival_cluster$Vital_status)

# Create variable group of clusters
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 1,
                                cluster == 2 ~ 2,
                                cluster == 3 ~ 3,
                                cluster == 4 ~ 1,
                                cluster == 5 ~ 3
                                #                                cluster == 6 ~ 23679,
                                #                                cluster == 7 ~ 23679,
                                #                                cluster == 8 ~ 1854,
                                #                                cluster == 9 ~ 23679
  ))

# Use vital_status or diseas_free and overall_survival_time_days or Event_Free_survival_days
surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])
surv_pvalue(surv_C_fit_TALL)
ggsurvplot(surv_C_fit_TALL, data = bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",], 
           pval = TRUE,
           size=2,
           palette = c("red", "gray87","gray60"),
           ggtheme = theme_bw()+theme(legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      text = element_text(size=20),
                                      axis.text.x = element_text(angle=45, size=20,hjust=1),
                                      axis.text.y = element_text(size=20,hjust=1)),
           font.x = c(20),
           font.y = c(20))

pairwise_survdiff(Surv(Event_free_surv_days, disease_free) ~ cluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])


library(peperr)
bcatexp_survival_cluster$Gender_surv<-ifelse(bcatexp_survival_cluster$Gender=="Female",1,2)
bcatexp_survival_cluster$newcluster_surv<-gsub('235','high_brca1_clusters',bcatexp_survival_cluster$newcluster)
bcatexp_survival_cluster$newcluster_surv<-gsub('1749','medium_brca1_clusters',bcatexp_survival_cluster$newcluster_surv)
bcatexp_survival_cluster$newcluster_surv<-gsub('6810','low_brca1_clusters',bcatexp_survival_cluster$newcluster_surv)

bcatexp_survival_cluster$newcluster_surv<-ifelse(bcatexp_survival_cluster$newcluster=="235","high_brca1_cl","medium/low_brca1_cl")


fit.coxph<- coxph(Surv(Event_free_surv_days, disease_free) ~ newcluster_surv+bcatgroup+kaisogroup+brcagroup+Pheno, 
                  data = bcatexp_survival_cluster[,bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm"],
                  ties = 'exact')
summary(fit.coxph)

ggforest(fit.coxph)


## DEGs remission vs relapse

#Groups
# T-ALL TARGET
# using REMISSION vs RELAPSE without second malignant neoplasm
codata<-as.data.frame(bcatexp_survival[rownames(bcatexp_survival) %in% rownames(patients),][,c(1,5,7,8,10:12,14,17,19,(ncol(bcatexp_survival)-4):ncol(bcatexp_survival))])
codata<-codata[codata$First_event!="Second Malignant Neoplasm",]
dim(codata)
#row.names(codata)<-codata$Row.names
#codata$Row.names<-NULL
codata$disease_free<-as.factor(codata$disease_free)

table(codata$disease_free,codata$First_event)


#make deseq object for disease free REMISSION vs RELAPSE TARGET
TARGETCountTable <- DESeqDataSetFromMatrix(
  countData=countdata[,colnames(countdata) %in% row.names(codata)],
  colData = codata,
  ~ disease_free)


colnames(TARGETCountTable)

#how many genes we capture, counting the number of genes that have non–zero counts in all samples.
GeneCounts <- counts(TARGETCountTable)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)}) 
sum(idx.nz)

### random sample from the count matrix
nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]


### NORMALIZATION
# Remove genes with low expression levels ()
keep <- rowSums(counts(TARGETCountTable)) >= 10
TARGETCountTable <- TARGETCountTable[keep,]

GeneCounts <- counts(TARGETCountTable)

#colData(TARGETCountTable)$cluster <- factor(colData(TARGETCountTable)$cluster, levels = c("1","2"))

#### estimate size factors
TARGETCountTable <- estimateSizeFactors(TARGETCountTable)
sizeFactors(TARGETCountTable)


### produce rlog-transformed data
rldTARGET <- varianceStabilizingTransformation(TARGETCountTable, blind=TRUE) ## create a distance matrix between the samples


###### Differential EXPRESSION ANALYSIS ###
TARGETCountTable <- estimateDispersions(TARGETCountTable)
plotDispEsts(TARGETCountTable)


# Statistical testing of DE genes
levels(TARGETCountTable$disease_free)
design(TARGETCountTable)

dds<-DESeq(TARGETCountTable)
resultsNames(dds)
rescondition<-results(dds)

# Example plot counts
plotCounts(dds, gene=which(rownames(rescondition)=="ENSG00000165323"), intgroup="disease_free")

head(rescondition)
summary(rescondition)
sum(rescondition$pvalue < 0.05,na.rm=TRUE)
sum(rescondition$padj < 0.1,na.rm=TRUE)
plotMA(rescondition)

# Selecting specific DEGs
target_prog_degs<-as.data.frame(rescondition)
target_prog_degs<-merge(target_prog_degs,geneid,by.x=0,by.y="Gene",all.x=TRUE)

countab<-as.data.frame(counts(dds, normalized=TRUE))
countab<-t(countab)
countab<-melt(countab)
colnames(countab)<-c("SampID","Gene","value")
countab<-merge(countab,geneid,by="Gene")

countab<-merge(countab,codata,by.x="SampID",by.y="row.names")
genesdown<-subset(target_prog_degs,target_prog_degs$pvalue<=0.05 & target_prog_degs$log2FoldChange<0)
genesdown<-genesdown[order(genesdown$pvalue),]
genesup<-subset(target_prog_degs,target_prog_degs$pvalue<=0.05 & target_prog_degs$log2FoldChange>0)
genesup<-genesup[order(genesup$pvalue),]

genesdown[genesdown$Row.names %in% colnames(bcatexp_survival[,c(21:(ncol(bcatexp_survival)-4))]),]
genesup[genesup$Row.names %in% colnames(bcatexp_survival[,c(21:(ncol(bcatexp_survival)-4))]),]

countab$group<-ifelse(countab$First_event=="Progression" | countab$First_event=="Death" | countab$First_event=="Relapse",
                      "Relapse","Remission")

## Differentially expressed genes relapse vs remission
write.table(target_prog_degs,"/Volumes/grcmc/YGUILLEN/Projects_bioinfo_data_Yolanda/bcat/Patients_transcriptomes/TARGET/supervised/DEGs_rem_relapse_TARGET.tsv",sep="\t",quote = FALSE,row.names = FALSE)

target_prog_degs_sig<-(subset(target_prog_degs,target_prog_degs$pvalue<0.05))
row.names(target_prog_degs_sig)<-target_prog_degs_sig$Row.names
target_prog_degs_sig$Row.names<-NULL

deg_bcattarg<-row.names(target_prog_degs_sig[row.names(target_prog_degs_sig) %in% colnames(bcatexp_survival[,21:(ncol(bcatexp_survival)-4)]),])


# Gene signature
ggplot(countab[countab$Gene %in% c(deg_bcattarg),],aes(x=group,y=log(value)))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=First_event),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_shape_manual(values=c(15,16,16,17,18,19,19,20,20,0))+
  scale_color_brewer(palette = "Spectral")+
  ylab("Normalized counts")+
  xlab('')+
  facet_wrap(~NAME,scales = "free",ncol=4)+
  theme_classic()+
  theme(axis.text.x = element_blank())

# Add b-cat, kaiso, lef and tcf
ggplot(countab[countab$Gene %in% c("ENSG00000168036","ENSG00000177485","ENSG00000081059","ENSG00000138795"),],aes(x=group,y=log(value)))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=First_event),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_shape_manual(values=c(15,16,16,17,18,19,19,20,20,0))+
  scale_color_brewer(palette = "Spectral")+
  ylab("Normalized counts")+
  xlab('')+
  facet_wrap(~NAME,scales = "free",ncol=4)+
  theme_classic()


rm(TCeltarg)
rm(nz.counts)

## Common targets bcat and brca1
brca_bcat<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/gene_lists/common_brca1_bcat.txt",header = FALSE)
colnames(brca_bcat)<-"GENE"

brca_bcat_genes<-rpmi_al[rpmi_al$external_gene_name.y %in% brca_bcat$GENE & !is.na(rpmi_al$external_gene_name.y),]$ensembl_gene_id

myBreaks <- c(seq(min(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcattarg]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% deg_bcattarg]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% deg_bcattarg]))),
               annotation_col = as.data.frame(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",c(7,8,14,17,19,21)]),
               annotation_row = colcol[row.names(colcol) %in% deg_bcattarg,c(17,18,19)],
               show_colnames=F,
               show_rownames = T,
               fontsize = 6,
               cutree_cols = 4,
               cutree_rows = 4,
               color=cols, breaks=myBreaks)


listbcat<-c("ENSG00000008300","ENSG00000064703","ENSG00000065029","ENSG00000093183","ENSG00000096401",
            "ENSG00000109805","ENSG00000124574","ENSG00000125871","ENSG00000133119","ENSG00000139629",
            "ENSG00000164068","ENSG00000175463","ENSG00000105810")


myBreaks <- c(seq(min(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% listbcat]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% listbcat]))), length.out=round(floor(paletteLength*0.10)))) 


phet<-pheatmap(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% listbcat]))),
               annotation_col = as.data.frame(bcatexp_survival[,c(7,8,14,17,19,21)]),
               #               annotation_row = colcol[,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = T,
               fontsize = 6,
               cutree_cols = 6,
               cutree_rows = 2,
               color=cols, breaks=myBreaks)



myBreaks <- c(seq(min(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% brca_bcat_genes]))),
               annotation_col = as.data.frame(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",c(5,7,8,10:12,14,17,19)]),
               annotation_row = colcol[row.names(colcol) %in% brca_bcat_genes,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 8,
               cutree_rows = 5,
               color=cols, breaks=myBreaks)

# random selection within DEGs relapse vs remission
listbcat<-sample(row.names(target_prog_degs_sig),23)

myBreaks <- c(seq(min(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% listbcat]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% listbcat]))), length.out=round(floor(paletteLength*0.10)))) 

phet<-pheatmap(as.data.frame(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% listbcat]))),
               annotation_col = as.data.frame(bcatexp_survival[,c(5,7,8,10:12,14,17,19)]),
               #               annotation_row = colcol[,c(14,16,17,18,19)],
               show_colnames=F,
               show_rownames = F,
               fontsize = 6,
               cutree_cols = 15,
               cutree_rows = 3,
               color=cols, breaks=myBreaks)


clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=4)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=4)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==2,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% deg_bcattarg])))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% row.names(target_prog_degs_sig)])))[phet$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes])))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% deg_bcattarg])))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(vastout_t_TLE[row.names(vastout_t_TLE) %in% row.names(bcatexp_survival),colnames(vastout_t_TLE) %in% row.names(target_prog_degs_sig)])))[,phet$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% brca_bcat_genes])))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(bcatexp_survival,patients,by="row.names")

patinfo$groups<-ifelse(patinfo$cluster=="1" | patinfo$cluster=="4", "P1","P2")


### Add TCF7 and LEF1 expression variable
tcf_lev<-data.frame(row.names = row.names(vastout_t),
                    tcf_exp = vastout_t[,colnames(vastout_t) == "ENSG00000081059"])

lef_lev<-data.frame(row.names = row.names(vastout_t),
                    lef_exp = vastout_t[,colnames(vastout_t) == "ENSG00000138795"])

tcflefexp<-merge(tcf_lev,lef_lev,by="row.names",sort=FALSE,all=TRUE)
row.names(tcflefexp)<-tcflefexp$Row.names
tcflefexp$Row.names<-NULL

row.names(patinfo)<-patinfo$Row.names
patinfo<-merge(patinfo,tcflefexp,by="row.names")
patinfo<-patinfo[,-2]

patinfo$outcome_group<-ifelse(patinfo$First_event=="Progression" | patinfo$First_event=="Death" | patinfo$First_event=="Relapse",
                              "Relapse","Remission")


ggplot(patinfo[patinfo$First_event!="Second Malignant Neoplasm",],aes(x=reorder(outcome_group,bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3,alpha=0.8)+
  geom_boxplot(alpha=0.2,lwd=1)+
  #  geom_hline(yintercept = median(scale(bcatexp_survival$bcat))+0.1,color="grey",lwd=1,linetype = "dotted")+
  # geom_hline(yintercept = median(scale(bcatexp_survival$brca1_exp))-0.1,color="grey",lwd=1,linetype = "dotted")+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[4],color="grey",lwd=2)+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))

patinfo_melt<-subset(patinfo,select=c("Row.names","First_event","outcome_group","groups","bcat","kaiso_exp","tcf_exp","lef_exp"))
patinfo_melt<-melt(patinfo_melt,id.vars=c("Row.names","First_event","outcome_group","groups"))

ggplot(patinfo_melt,aes(x=outcome_group,y=log(value)))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=First_event),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_color_brewer(palette = "Spectral")+
  ylab("log normalized expression")+
  xlab('')+
  facet_wrap(~variable,scales = "free",ncol=4)+
  theme_classic()

wilcox.test(log(patinfo[patinfo$First_event!="Second Malignant Neoplasm",]$lef_exp)~patinfo[patinfo$First_event!="Second Malignant Neoplasm",]$outcome_group)

patinfo$scaled_brca1<-as.numeric(scale(patinfo$brca1_exp))
median(patinfo$scaled_brca1)
median(patinfo[patinfo$cluster=="1",]$scaled_brca1)

## survival curves with all patients, REDO PATIENTS WITH HEATMAP

bcatexp_survival_cluster<-merge(patients,bcatexp_survival,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Sample<-row.names(bcatexp_survival_cluster)
bcatexp_survival_cluster$Row.names<-NULL

library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)
bcatexp_survival_cluster$Vital_status<-ifelse(bcatexp_survival_cluster$Vital_status=="Alive",0,1)
bcatexp_survival_cluster$Vital_status<-as.integer(bcatexp_survival_cluster$Vital_status)

# Create variable group of clusters
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 14,
                                cluster == 2 ~ 23,
                                cluster == 3 ~ 23,
                                cluster == 4 ~ 14
  ))

# Use vital_status or diseas_free and overall_survival_time_days or Event_Free_survival_days
surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])

surv_C_fit_TALL <- survfit(Surv(Overall_survival_time_days, Vital_status) ~ newcluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])

surv_pvalue(surv_C_fit_TALL)
ggsurvplot(surv_C_fit_TALL, 
           pval = TRUE,
           size=2,
           palette = c("red", "gray60"),
           ggtheme = theme_bw()+theme(legend.position = "bottom",
                                      legend.text = element_text(size=14),
                                      text = element_text(size=20),
                                      axis.text.x = element_text(angle=45, size=20,hjust=1),
                                      axis.text.y = element_text(size=20,hjust=1)),
           font.x = c(20),
           font.y = c(20))


# Add metadata (Age, WBC count, CNS, phenotype, lession)
metadata_target_add<-read_excel("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/metadata_all.xlsx",sheet="TARGET")
row.names(metadata_target_add)<-metadata_target_add$Sample

datat<-merge(bcatexp_survival_cluster,metadata_target_add,by.x="Sample",by.y="GSM")
row.names(datat)<-datat$Sample
datat$Row.names<-NULL  

datat<-subset(datat,datat$First_event!="Second Malignant Neoplasm")
datat$newcluster_surv<-ifelse(datat$newcluster=="14","C_14","C_23")

datat$Age_Dx<-ifelse(datat$`Age at Diagnosis in years`<10,"<10",">10")
datat$WBC_Dx<-ifelse(datat$`WBC at Diagnosis`<200,"<200",">200")
datat$MRD_29<-ifelse(datat$`MRD Day 29`>0,">0","0")

datat$newcluster_surv <- factor(datat$newcluster_surv, levels = c("C_23","C_14"))
datat$Age_Dx <- factor(datat$Age_Dx, levels = c(">10","<10"))
datat$MRD_29 <- factor(datat$MRD_29, levels = c("0",">0"))

datat$CNS_hr <- gsub('2.*','2',datat$`CNS Status at Diagnosis`)
datat$CNS_hr <- gsub('3.*','3',datat$CNS_hr)

fit.coxph<- coxph(Surv(Event_free_surv_days, disease_free) ~ newcluster_surv+Gender.x+Age_Dx+WBC_Dx+MRD_29+CNS_hr, 
                  data = datat[datat$CNS_hr!=".",],
                  ties = 'exact')
summary(fit.coxph)

ggforest(fit.coxph)
##

write.table(clusters,"/Users/yguillen/Downloads/clusters_TARGET_bcatpeaks.txt",sep="\t",quote = FALSE,row.names = FALSE)


# Display genes in each cluster and order from the heatmap and heatmap only these
# clust1_microarray from script T_ALL_GEOquery.R includes all genes in cluster1 separating patients 145 vs 2367

# targdeg are b-cat targets differentially expressed between 145 and 2367 (deg_clust_GSE14618)
#deg_bcat_prog<-unique(deg_clust_GSE14618[deg_clust_GSE14618$external_gene_name %in% targdeg,]$ensembl_gene_id)
unique(comp_all_tot$ensembl_gene_id)



# corout_brca matrix with p and R for brca1
corout_brca[corout_brca$Pval<0.05,]$Gene

#brca and kaiso
corout_brca[corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,]$Gene
corout_kai[corout_kai$Pval<0.01 & abs(corout_kai$Rho)>0.2,]$Gene

kaibrcacor<-unique(c(corout_brca[corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,]$Gene,
                     corout_kai[corout_kai$Pval<0.01 & abs(corout_kai$Rho)>0.2,]$Gene))

myBreaks <- c(seq(min(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)]))), length.out=round(floor(paletteLength*0.10)))) 

#row.names(corout_brca)<-corout_brca$Gene

colgenes<-comp_all_tot[comp_all_tot$ensembl_gene_id %in% colnames(bcatexp_survival),]
colgenes<-colgenes[!duplicated(colgenes$ensembl_gene_id),]
row.names(colgenes)<-colgenes$ensembl_gene_id

phetsel<-pheatmap(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)]))),
                  annotation_col = as.data.frame(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",c(1,5,7,8,10:12,14,17,19,(ncol(bcatexp_survival)-3):ncol(bcatexp_survival))]),
                  #         annotation_row = corout_brca[corout_brca$Gene %in% colnames(bcatexp_survival) & corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,][,3:4],
                  annotation_row = colgenes[colgenes$comp=="Rem_indfail",c(14,16)],
                  show_colnames=F,
                  show_rownames = F,
                  fontsize = 6,
                  cutree_cols =6,
                  cutree_rows = 3,
                  color=cols, breaks=myBreaks)


# deg_clust_GSE14618 has all genes DEGs between the previous two groups
deg_micro<-vastout_t_TLE[,colnames(vastout_t_TLE) %in% unique(genedeg$ensembl_gene_id),]
deg_micro<-deg_micro[row.names(deg_micro) %in% row.names(bcatexp_survival),]

myBreaks <- c(seq(min(t(scale(deg_micro))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(deg_micro))), length.out=round(floor(paletteLength*0.10)))) 


phetsel<-pheatmap(as.data.frame(t(scale(deg_micro))),
                  annotation_col = as.data.frame(bcatexp_survival[,c(5,7,8,10:12,14,17,19)]),
                  annotation_row = genedeg[,c(1,3)],
                  show_colnames=F,
                  show_rownames = F,
                  fontsize = 6,
                  cutree_cols =8,
                  cutree_rows = 6,
                  color=cols, breaks=myBreaks)


clusters<-data.frame(cluster=sort(cutree(phetsel$tree_row, k=6)))
patients<-data.frame(cluster=sort(cutree(phetsel$tree_col, k=6)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==2,]$NAME

# Order genes (clusters)
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(deg_micro)))[phetsel$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% clusters[clusters$cluster==6 | clusters$cluster==3,]$Gene])))[phetsel$tree_row[["order"]],]))
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail,"]$ensembl_gene_id)])))[phetsel$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% corout_brca[corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,]$Gene])))[phetsel$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(deg_micro)))[,phetsel$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% clusters[clusters$cluster==6 | clusters$cluster==3,]$Gene])))[,phetsel$tree_col[["order"]]]))
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[bcatexp_survival$First_event!="Second Malignant Neoplasm",colnames(bcatexp_survival) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail,"]$ensembl_gene_id)])))[,phetsel$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcatexp_survival[,colnames(bcatexp_survival) %in% corout_brca[corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,]$Gene])))[,phetsel$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(bcatexp_survival,patients,by="row.names")

ggplot(patinfo[patinfo$First_event!="Second Malignant Neoplasm",],aes(x=reorder(as.factor(cluster), bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3)+
  geom_boxplot(alpha=0.2,lwd=1)+
  geom_hline(yintercept = median(scale(bcatexp_survival$bcat)),color="red",lwd=2)+
  geom_hline(yintercept = quantile(scale(bcatexp_survival$bcat))[4],color="grey",lwd=2)+
  geom_hline(yintercept = quantile(scale(bcatexp_survival$bcat))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))



## survival curves with all patients, REDO PATIENTS WITH HEATMAP

bcatexp_survival_cluster<-merge(patients,bcatexp_survival,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Row.names<-NULL


library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)
bcatexp_survival_cluster$Vital_status<-ifelse(bcatexp_survival_cluster$Vital_status=="Alive",0,1)
bcatexp_survival_cluster$Vital_status<-as.integer(bcatexp_survival_cluster$Vital_status)

# Create variable group of clusters
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 1,
                                cluster == 2 ~ 256,
                                cluster == 3 ~ 3,
                                cluster == 4 ~ 4,
                                cluster == 5 ~ 256,
                                cluster == 6 ~ 256
              #                  cluster == 7 ~ 247,
              #                  cluster == 8 ~ 13658
                                #                                cluster == 9 ~ 49731112,
                                #                                cluster == 10 ~ 2810,
                                #                                cluster == 11 ~ 49731112,
                                #                                cluster == 12 ~ 49731112
  ))


#surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster)

surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days/365, disease_free) ~ newcluster, data=bcatexp_survival_cluster[bcatexp_survival_cluster$First_event!="Second Malignant Neoplasm",])

ggsurvplot(surv_C_fit_TALL, pval = TRUE,
           font.x = c(20),
           font.y = c(20))
#           palette = c("green", "indianred","lightseagreen","purple"))

pairwise_survdiff(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster)

surv_C_fit_TALL$time




#### Reducing b-cat targets signature

flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

#correlation matrix with all targets
rcorTALL<-rcorr(as.matrix(bcatexp_survival[,20:(ncol(bcatexp_survival)-4)]), type = "spearman")

#correlation matrix with targets cluster 1 microarray data
rcorTALL<-rcorr(as.matrix(bcatexp_survival[,colnames(bcatexp_survival) %in% clust1_microarray]), type = "spearman")

# correlation matrix
corrplot(rcorTALL$r, order="hclust",tl.col="black", tl.cex = 0.1, type = "lower",
         p.mat = rcorTALL$P, sig.level = 0.05, insig = "blank", method = "color")


ordercorTALL <- flat_cor_mat(rcorTALL$r, rcorTALL$P)
ordercorTALL<-ordercorTALL[!is.na(ordercorTALL$p),]

ordercorTALL<-merge(ordercorTALL,geneid,by.x="row",by.y="Gene",all.x=TRUE)

# selection by p-value (0.85 FOR ALL, 0.8 FOR CLUST1_MICROARRAY SUBSET)
pair1<-as.vector(unique(ordercorTALL[ordercorTALL$p<0.01 & abs(ordercorTALL$cor)>0.8,]$row))
pair2<-as.vector(unique(ordercorTALL[ordercorTALL$p<0.01 & abs(ordercorTALL$cor)>0.8,]$column))

genes<-unique(c(pair1,pair2))
length(genes)


clusterTALL <-bcatexp_survival[,genes]
dim(clusterTALL)

clusterTALL<-t(clusterTALL)
clusterTALL<-merge(geneid,clusterTALL,by.y="row.names",by.x="Gene",all.y=TRUE)
row.names(clusterTALL)<-clusterTALL$NAME
clusterTALL$Gene<-NULL
clusterTALL$NAME<-NULL
clusterTALL<-t(clusterTALL)

rclusterTALL<-rcorr(as.matrix(clusterTALL), type = "pearson")
corrplot(rclusterTALL$r, order="hclust",tl.col="black", tl.cex = 0.5,
         p.mat = rclusterTALL$P, sig.level = 0.01, insig = "blank",addrect = 2)



genes

# remove two elements from genes list x[x != "b"], THESE ARE TWO HEMOGLOBIN SUBUNITS, THEY HAVE A HIGH CORRELATION, BUT ONLY AMONG THEM
genes<-genes[genes != "ENSG00000244734"]
genes<-genes[genes != "ENSG00000086506"]

#clusterTALL <-bcatexp_survival[,genes]
clusterTALL <-bcatexp_survival[,deg_bcattarg]


dim(clusterTALL)

rclusterTALL<-rcorr(as.matrix(clusterTALL), type = "pearson")

cor(as.matrix(clusterTALL))
corrplot(rclusterTALL$r, order="hclust",tl.col="black", tl.cex = 0.5,
         p.mat = rclusterTALL$P, sig.level = 0.05, insig = "blank",addrect = 4,tl.pos = "n")


# Using the reduced gene list from correlation networks
bcat_signature<-genes
sort(geneid[geneid$Gene %in% bcat_signature,]$NAME)

bcat_matrix<-bcatexp_survival[,colnames(bcatexp_survival) %in% bcat_signature,]

# or using gene liste crossing DEGs leukemia vs thymus in scRNASeq data
scsig<-c("ROMO1","SERINC5","TCF7","FAM40A","BCLAF1",
         "RBM39","RPS28","DDX3X","NDUFA11","SLC38A2","UBE2D3")

scsig<-geneid[geneid$NAME %in% scsig,]$Gene

bcat_matrix<-bcatexp_survival[,colnames(bcatexp_survival) %in% scsig,]


bcat_matrix<-merge(subset(bcatexp_survival,select=c(bcat,kaiso_exp,brca1_exp,
                                                    bcatgroup,kaisogroup,brcagroup,disease_free,Vital_status,First_event,Pheno,
                                                    NOTCH1_freq,FBXW7_freq,LEF1_freq,Event_free_surv_days)),bcat_matrix,by="row.names",sort=FALSE)
dim(metadata)
dim(bcat_matrix)

row.names(bcat_matrix)<-bcat_matrix$Row.names
bcat_matrix$Row.names<-NULL


myBreaks1 <- c(seq(min(t(scale(bcat_matrix[,15:ncol(bcat_matrix)]))), 2.5, length.out=round(ceiling(paletteLength*0.90))),
               seq(2.6,max(t(scale(bcat_matrix[,15:ncol(bcat_matrix)]))),length.out = round(ceiling(paletteLength*0.10)))) 

colcolmat<-colcol[row.names(colcol) %in% colnames(as.data.frame(bcat_matrix[,15:ncol(bcat_matrix)])),]

phetsel<-pheatmap(as.data.frame(t(scale(bcat_matrix[,15:ncol(bcat_matrix)]))),
                  annotation_col = as.data.frame(bcat_matrix[,c(1:13)]),
                  annotation_row = colcolmat[,c(14,16,17,18,19)],
                  show_colnames=F,
                  show_rownames = F,
                  fontsize = 6,
                  cutree_cols = 6,
                  cutree_rows = 2,
                  color=cols, breaks=myBreaks1)

# order by bcat group
matplot<-as.data.frame(t(scale(bcat_matrix[,15:ncol(bcat_matrix)])))
matplot[, order(bcat_matrix$bcatgroup)]
matplot<-matplot[,row.names(bcat_matrix[order(bcat_matrix$bcatgroup),])]

pheatmap(matplot,
         annotation_col = as.data.frame(bcat_matrix[,c(1:13)]),
         annotation_row = colcolmat[,c(14,16,17,18,19)],
         show_colnames=F,
         show_rownames = F,
         fontsize = 6,
         cutree_cols = FALSE,
         cluster_cols = FALSE,
         cutree_rows = 5,
         color=cols, breaks=myBreaks1)

##

clusters<-data.frame(cluster=sort(cutree(phetsel$tree_row, k=2)))
patients<-data.frame(cluster=sort(cutree(phetsel$tree_col, k=6)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==1,]$NAME

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(bcat_matrix[,15:ncol(bcat_matrix)])))[phetsel$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(bcat_matrix[,15:ncol(bcat_matrix)])))[,phetsel$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)




## survival curves

bcatexp_survival_cluster<-merge(patients,bcat_matrix,by=0)
row.names(bcatexp_survival_cluster)<-bcatexp_survival_cluster$Row.names
bcatexp_survival_cluster$Row.names<-NULL

table(bcatexp_survival_cluster$cluster,bcatexp_survival_cluster$Vital_status)

# Create variable group of clusters
bcatexp_survival_cluster$newcluster<-NA
bcatexp_survival_cluster<-bcatexp_survival_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 13,
                                cluster == 2 ~ 2456,
                                cluster == 3 ~ 13,
                                cluster == 4 ~ 2456,
                                cluster == 5 ~ 2456,
                                cluster == 6 ~ 2456))


library(survminer)
bcatexp_survival_cluster$disease_free<-as.integer(bcatexp_survival_cluster$disease_free)
surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, disease_free) ~ newcluster, data=bcatexp_survival_cluster)
ggsurvplot(surv_C_fit_TALL, data = bcatexp_survival_cluster, pval = TRUE,
           font.x = c(20),
           font.y = c(20))
#           palette = c("green", "indianred","lightseagreen","purple"))

pairwise_survdiff(Surv(Event_free_surv_days, disease_free) ~ cluster, data=bcatexp_survival_cluster)

surv_C_fit_TALL$time


# Display genes in each cluster and order from the heatmap and heatmap only these

myBreaks <- c(seq(min(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene]))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene]))), length.out=round(floor(paletteLength*0.10)))) 

# annotate genes from microarrays unsupervised T_ALL_GEOquary.R --> colgene


phetsel<-pheatmap(as.data.frame(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene]))),
                  annotation_col = as.data.frame(bcat_matrix[,c(1:13)]),
                  annotation_row = colcolmat[,c(14,16,17,18,19)],
                  show_colnames=F,
                  show_rownames = T,
                  fontsize = 6,
                  cutree_cols = 4,
                  #  cutree_rows = 3,
                  color=cols, breaks=myBreaks)


patients<-data.frame(cluster=sort(cutree(phetsel$tree_col, k=4)))
matsel<-as.data.frame(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene])))

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(matsel[,phetsel$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)



patinfosel<-merge(t(as.data.frame(t(scale(bcat_matrix[,colnames(bcat_matrix) %in% clusters[clusters$cluster==4 | clusters$cluster==3 | clusters$cluster==1,]$Gene])))),patients,by="row.names")
row.names(patinfosel)<-patinfosel$Row.names
patinfosel$Row.names<-NULL
patinfosel<-merge(bcat_matrix[,c(1:13)],patinfosel,by="row.names")
row.names(patinfosel)<-patinfosel$Row.names
patinfosel$Row.names<-NULL

patinfosel$cluster<-as.factor(patinfosel$cluster)
table(patinfosel$cluster,patinfosel$Vital_status)

ggplot(patinfosel,aes(x=reorder(as.factor(cluster), bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3)+
  geom_boxplot(alpha=0.2,lwd=1)+
  geom_hline(yintercept = median(scale(bcat_matrix$bcat)),color="red",lwd=2)+
  #  geom_hline(yintercept = quantile(scale(bcatexp_survival$brca1_exp))[4],color="grey",lwd=2)+
  geom_hline(yintercept = quantile(scale(bcat_matrix$bcat))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))

# output for ssGSEA
write.table(clusters[clusters$cluster==4,]$NAME,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/cluster4.txt",row.names=FALSE,quote = FALSE)
write.table(clusters[clusters$cluster==5 | clusters$cluster==2,]$NAME,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/cluster2_5.txt",row.names=FALSE,quote = FALSE)
write.table(clusters[clusters$cluster==1 | clusters$cluster==3,]$NAME,"/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/cluster1_3.txt",row.names=FALSE,quote = FALSE)

# Results ssGSEA
ssGSEAout<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/output/TALL_clusters_mod.gct")
ssGSEAout<-as.data.frame(t(ssGSEAout))
colnames(ssGSEAout)<-ssGSEAout[1,]
ssGSEAout<-ssGSEAout[-c(1:2),]

## For other cohorts, first remove genes with sum zero values
i <- (colSums(bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",22:ncol(bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",])], na.rm=T) != 0)
matnonzero <- bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",22:ncol(bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",])][, i]

matnonzero<-cbind(bcatexp[bcatexp$GSEA=="EGA_R" | bcatexp$GSEA=="EGA_TLE",1:21],matnonzero)

matnonzero<-matnonzero[row.names(matnonzero)!="Thymus" & row.names(matnonzero)!="R3",]

matnonzero<-subset(matnonzero,matnonzero$First_event!="Secondary leukemia" & matnonzero$First_event!="Disease" & matnonzero$First_event!="Death")


myBreaks1 <- c(seq(min(t(scale(matnonzero[,22:ncol(matnonzero)]))), 0, length.out=round(ceiling(paletteLength*0.40)) + 1),
               seq(0.1,2,length.out = round(ceiling(paletteLength*0.40))),
               seq(2.1, max(t(scale(matnonzero[,22:ncol(matnonzero)]))), length.out=round(floor(paletteLength*0.20)-1))) 

# annotate genes from unsupervised groups microarrays G1, G2, G3
clusters_micro_uns_ensembl<-merge(geneid,clusters_micro_uns,by.x="NAME",by.y="Gene",all.y=TRUE)

colgene_ega<-merge(as.data.frame(t(scale(matnonzero[,22:ncol(matnonzero)]))),clusters_micro_uns_ensembl,by.x=0,by.y="Gene",all.x=TRUE)
row.names(colgene_ega)<-colgene_ega$Row.names
colgene_ega$Row.names<-NULL

phet1<-pheatmap(as.data.frame(t(scale(matnonzero[,22:ncol(matnonzero)]))),
                annotation_col = as.data.frame(matnonzero[,c(14,17,19,21)]),
                annotation_row = data.frame(row.names=row.names(colgene_ega),
                                            cluster=paste("cluster",as.factor(colgene_ega[,(ncol(colgene_ega))]))),
                show_colnames=T,
                show_rownames = F,
                fontsize = 6,
                cutree_cols = 3,
                cutree_rows = 3,
                color=cols,
                breaks = myBreaks1)


myBreaks1 <- c(seq(min(t(scale(matnonzero[,colnames(matnonzero) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)]))), 0, length.out=round(ceiling(paletteLength*0.40)) + 1),
               seq(0.1,2,length.out = round(ceiling(paletteLength*0.40))),
               seq(2.1, max(t(scale(matnonzero[,colnames(matnonzero) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)]))), length.out=round(floor(paletteLength*0.20)-1))) 


matnonzero$Age_Dx_days<-(as.numeric(matnonzero$Age_Dx_days))/365
colnames(matnonzero)[21]<-"Age_Dx"

phet1<-pheatmap(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)]))),
                annotation_col = as.data.frame(matnonzero[,c(14,17,19,21)]),
                #                annotation_row = colcol[row.names(colcol) %in% colnames(matnonzero[,colnames(matnonzero) %in% deg_bcattarg]),c(6,16)],
                show_colnames=F,
                show_rownames = F,
                fontsize = 6,
                cutree_cols = 2,
                cutree_rows = 4,
                color=cols,
                breaks = myBreaks1)

phet1<-pheatmap(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% colnames(bcatexp[,22:ncol(bcatexp)])]))),
                annotation_col = as.data.frame(matnonzero[,c(14,17,19,21)]),
                annotation_row = colcol[row.names(colcol) %in% colnames(matnonzero[,colnames(matnonzero) %in% deg_bcattarg]),c(6,16)],
                show_colnames=F,
                show_rownames = F,
                fontsize = 6,
                cutree_cols = 3,
                cutree_rows = 4,
                color=cols,
                breaks = myBreaks1)

phet1<-pheatmap(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% geneid[geneid$NAME %in% targdeg,]$Gene]))),
                annotation_col = as.data.frame(matnonzero[,c(14,17,19,21)]),
                #                annotation_row = colcol[row.names(colcol) %in% colnames(matnonzero[,colnames(matnonzero) %in% deg_bcattarg]),c(6,16)],
                show_colnames=F,
                show_rownames = F,
                fontsize = 6,
                cutree_cols = 3,
                cutree_rows = 4,
                color=cols,
                breaks = myBreaks1)

colnames(matnonzero[,colnames(matnonzero) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)])
colgenesmat<-comp_all_tot[comp_all_tot$ensembl_gene_id %in% colnames(matnonzero[,colnames(matnonzero) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)]),]
colgenesmat<-colgenesmat[!duplicated(colgenesmat$ensembl_gene_id),]
row.names(colgenesmat)<-colgenesmat$ensembl_gene_id

phet1<-pheatmap(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)]))),
                annotation_col = as.data.frame(matnonzero[,c(5,8,14,17,19,21)]),
                annotation_row = colgenesmat[,c(14,16)],
                show_colnames=T,
                show_rownames = F,
                fontsize = 6,
                cutree_cols = 3,
                cutree_rows = 3,
                color=cols,
                breaks = myBreaks1)


clusters1<-data.frame(cluster=sort(cutree(phet1$tree_row, k=3))) 
patients1<-data.frame(cluster=sort(cutree(phet1$tree_col, k=3)))




table(patients1$cluster)
table(clusters1$cluster)
clusters1$Gene<-row.names(clusters1)
clusters1[clusters1$cluster==1,]$Gene

# Order genes (clusters)
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% deg_bcattarg])))[phet1$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA"  & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% deg_bcattarg])))[phet1$tree_row[["order"]],]))
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(matnonzero[,22:ncol(matnonzero)])))[phet1$tree_row[["order"]],]))
#geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)])))[phet1$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters1<-merge(clusters1,geneorder,by="Gene")
clusters1$order<-as.numeric(clusters1$order)

clusters1<-merge(clusters1,geneid,by="Gene",all.x=TRUE)

# Order patients (patients)
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% deg_bcattarg])))[,phet1$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(matnonzero[row.names(matnonzero)!="Thymus" & matnonzero$First_event!="NA" & matnonzero$First_event!="Disease" & matnonzero$First_event!="Secondary leukemia",colnames(matnonzero) %in% deg_bcattarg])))[,phet1$tree_col[["order"]]]))
#patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(matnonzero[,colnames(matnonzero) %in% unique(comp_all_tot[comp_all_tot$comp=="Rem_indfail",]$ensembl_gene_id)])))[,phet1$tree_col[["order"]]]))
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(matnonzero[,22:ncol(matnonzero)])))[,phet1$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients1<-merge(patients1,patientorder,by=0)
row.names(patients1)<-patients1$Row.names
patients1$Row.names<-NULL
patients1$order<-as.numeric(patients1$order)


patinfo1<-merge(matnonzero,patients1,by="row.names")


### Add TCF7 and LEF1 and bcat and kaiso in EGA_TLE expression variable
tcf_lev<-data.frame(row.names = row.names(vastout_t),
                    tcf_exp_EGA = vastout_t[,colnames(vastout_t) == "ENSG00000081059"])

lef_lev<-data.frame(row.names = row.names(vastout_t),
                    lef_exp_EGA = vastout_t[,colnames(vastout_t) == "ENSG00000138795"])

bcat_lev<-data.frame(row.names = row.names(vastout_t),
                     bcat_exp_EGA = vastout_t[,colnames(vastout_t) == "ENSG00000168036"])

kaiso_lev<-data.frame(row.names = row.names(vastout_t),
                      kaiso_exp_EGA = vastout_t[,colnames(vastout_t) == "ENSG00000177485"])

tf_expTLE<-merge(tcf_lev,lef_lev,by="row.names",all=TRUE)
row.names(tf_expTLE)<-tf_expTLE$Row.names
tf_expTLE$Row.names<-NULL
tf_expTLE<-merge(tf_expTLE,bcat_lev,by="row.names",all=TRUE)
row.names(tf_expTLE)<-tf_expTLE$Row.names
tf_expTLE$Row.names<-NULL
tf_expTLE<-merge(tf_expTLE,kaiso_lev,by="row.names",all=TRUE)
row.names(tf_expTLE)<-tf_expTLE$Row.names
tf_expTLE$Row.names<-NULL

row.names(patinfo1)<-patinfo1$Row.names
patinfo1<-merge(patinfo1,tf_expTLE,by="row.names")
patinfo1<-patinfo1[,-2]

#patinfo1$outcome_group<-ifelse(patinfo1$First_event=="Progression" | patinfo1$First_event=="Death" | patinfo1$First_event=="Relapse",
#                               "Relapse","Remission")

patinfo1$cluster<- factor(patinfo1$cluster, levels=c("2","1","3"))
ggplot(patinfo1,aes(x=reorder(as.factor(cluster), bcat_exp_EGA, FUN = median),y=bcat_exp_EGA))+
  geom_point(aes(color=First_event),size=5)+
  geom_boxplot(alpha=0.2,lwd=1)+
  #  scale_color_manual(values=c("chartreuse3","hotpink","deepskyblue"))+
  #  geom_hline(yintercept = median(scale(patinfo$bcat_exp_prob1)),color="red",lwd=2)+
  #  geom_hline(yintercept = median(scale(patinfo$bcat_exp_prob1))-0.2,color="grey",lwd=2)+
  #  geom_hline(yintercept = median(scale(patinfo$kaiso_exp_prob2))+0.25,color="grey",lwd=2)+
  # geom_hline(yintercept = 0,color="grey",lwd=1)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_text(angle=45, size=16,hjust=1),
        axis.title.x = element_blank(),
        panel.border = element_rect(size = 2))

patinfo1$bcat_scale<-as.numeric(scale(patinfo1$bcat_exp_EGA))
patinfo1$kaiso_scale<-as.numeric(scale(patinfo1$kaiso_exp_EGA))
patinfo1$tcf_scale<-as.numeric(scale(patinfo1$tcf_exp_EGA))
patinfo1$lef_scale<-as.numeric(scale(patinfo1$lef_exp_EGA))

wilcox.test(patinfo1$bcat_exp_EGA~as.factor(patinfo1$cluster))

kruskal.test(patinfo1$bcat_exp_EGA~as.factor(patinfo1$cluster))

patinfo1_melt<-subset(patinfo1,select=c("Row.names","First_event","cluster","bcat_scale","kaiso_scale","tcf_scale","lef_scale"))
patinfo1_melt<-melt(patinfo1_melt,id.vars=c("Row.names","First_event","cluster"))

class(patinfo1_melt$value)<-as.numeric(patinfo1_melt$value)
ggplot(patinfo1_melt,aes(x=First_event,y=value))+
  # geom_point(color="black",size=5,alpha=0.5)+
  geom_jitter(aes(color=First_event),size=2,alpha=0.8)+
  geom_boxplot(alpha=0.2)+
  scale_color_brewer(palette = "Spectral")+
  ylab("scaled expression")+
  xlab('')+
  facet_wrap(~variable,scales = "free",ncol=4)+
  theme_classic()


prop_clust<-data.frame(prop.table(table(patinfo1$cluster,patinfo1$First_event),1))
colnames(prop_clust)<-c("cluster","First_event","Proportion")
prop_clust$cluster<-as.factor(prop_clust$cluster)

prop_clust$cluster<- factor(prop_clust$cluster, levels=c("2","3","1"))

ggplot(prop_clust,aes(x=cluster,y=Proportion,fill=First_event))+
  geom_col(color="black")+
  scale_fill_manual(values=c("dodgerblue","green","darkolivegreen","cyan"))+
  theme_bw()+
  theme(legend.position="right",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=16))



## Read ssGSEA from b-cat signature
#UP AND DOWN SEPARATE all genes
ssbcat<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/ssGSEA_output/TALL_ALL.gct")

#UP AND DOWN SEPARATE excluding SAP18,HLA-DOA and KRT5
ssbcat<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/ssGSEA_output/TALL_reduced.gct")

#clusters unsupervised microarrays in all TARGETs
ssbcat<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/unsupervised_microarrays_to_TARGET/TALL_bcat.gct")

ssbcat<-as.data.frame(t(ssbcat))
colnames(ssbcat)<-ssbcat[1,]
ssbcat<-ssbcat[-c(1,2),]

metadata_ssbcat<-merge(ssbcat,metadata,by=0,all.y=TRUE)
row.names(metadata_ssbcat)<-metadata_ssbcat$Row.names
metadata_ssbcat$Row.names<-NULL


metadata_ssbcat<-merge(metadata_ssbcat,bcat_lev,by=0)
#metadata_ssbcat$deg_up_relapse<-as.numeric(metadata_ssbcat$deg_up_relapse)
#metadata_ssbcat$deg_down_relapse<-as.numeric(metadata_ssbcat$deg_down_relapse)

metadata_ssbcat$cluster1<-as.numeric(metadata_ssbcat$cluster1)
metadata_ssbcat$cluster2<-as.numeric(metadata_ssbcat$cluster2)
metadata_ssbcat$cluster3<-as.numeric(metadata_ssbcat$cluster3)

cohorts<-c("TARGET","EGA_R","EGA_TLE")
cohorts<-c("TARGET")
cohorts<-c("EGA_R","EGA_TLE")

metadata_ssbcat<-subset(metadata_ssbcat,(metadata_ssbcat$First_event!="Second Malignant Neoplasm" & 
                                           metadata_ssbcat$First_event!="Secondary leukemia" &
                                           metadata_ssbcat$First_event!="Disease"))

summary(metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts , ]$cluster1)
boxplot(metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts , ]$cluster1)

ggplot(metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts , ],aes(x=cluster1,y=cluster3))+
  geom_point(data=metadata_ssbcat[metadata_ssbcat$First_event!="NA" & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],color="black",size=3)+
  geom_point(data=metadata_ssbcat[metadata_ssbcat$First_event!="NA" & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event),size=2.5)+
  scale_color_manual(values=c("red","grey","darkolivegreen","green","orange"))+
#  geom_density_2d(data=metadata_ssbcat[metadata_ssbcat$First_event!="Remission",],color="darkolivegreen",bins=20)+
#  geom_vline(xintercept = 130,linetype="dashed")+
#  geom_vline(xintercept = 100,linetype="dashed")+
  geom_vline(xintercept = 93,linetype="dashed")+
  geom_vline(xintercept = 115,linetype="dashed")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.text = element_text(size=12),
        text = element_text(size=10),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16))+
  coord_fixed(ratio=1)



boxplot(metadata_ssbcat[metadata_ssbcat$First_event!="NA" & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts & metadata_ssbcat$cluster1>130,]$bcat_exp_EGA)
boxplot(metadata_ssbcat[metadata_ssbcat$First_event!="NA" & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts & metadata_ssbcat$cluster1<130,]$bcat_exp_EGA)

clust1_clas<-rbind(data.frame(cluster1_lev="highest",metadata_ssbcat[metadata_ssbcat$cluster1>130,]),
      data.frame(cluster1_lev="lowest",metadata_ssbcat[metadata_ssbcat$cluster1<100,]),
      data.frame(cluster1_lev="medium",metadata_ssbcat[metadata_ssbcat$cluster1<130 & metadata_ssbcat$cluster1>100,]))



ggplot(clust1_clas[clust1_clas$GSEA %in% cohorts, ],aes(x=cluster1_lev,y=bcat_exp_EGA))+
  geom_jitter(data=clust1_clas[!is.na(clust1_clas$First_event) & clust1_clas$First_event!="Censored" & clust1_clas$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),size=3)+
  geom_boxplot(data=clust1_clas[!is.na(clust1_clas$First_event) & clust1_clas$First_event!="Censored" & clust1_clas$GSEA %in% cohorts,],alpha=0.1)+
  scale_color_manual(values=c("red","grey","darkolivegreen","green"))+
  theme_bw()+
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=16))

kruskal.test(clust1_clas$bcat_exp_EGA,clust1_clas$cluster1_lev)


library(survminer)

clust1_clas<-clust1_clas %>% mutate(disease_free =  ifelse(First_event == "Remission", 0, 1))
clust1_clas$disease_free<-as.factor(clust1_clas$disease_free)
dim(clust1_clas)

clust1_clas$disease_free<-as.integer(clust1_clas$disease_free)
clust1_clas$Vital_status<-ifelse(clust1_clas$Vital_status=="Alive",0,1)
clust1_clas$Vital_status<-as.integer(clust1_clas$Vital_status)



surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days/365, disease_free) ~ cluster1_lev, data=clust1_clas[clust1_clas$First_event!="Second Malignant Neoplasm",])

ggsurvplot(surv_C_fit_TALL, pval = TRUE,
           font.x = c(20),
           font.y = c(20))
#           palette = c("green", "indianred","lightseagreen","purple"))

pairwise_survdiff(Surv(Event_free_surv_days, disease_free) ~ cluster, data=data_boot_cluster)

surv_C_fit_TALL$time





ggplot(metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts, ],aes(x=First_event,y=cluster1))+
  geom_point(data=metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),color="black",size=4)+
  geom_point(data=metadata_ssbcat[!is.na(metadata_ssbcat$First_event) & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),size=3)+
  geom_boxplot(data=metadata_ssbcat[!is.na(metadata_ssbcat$First_event) & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],alpha=0.1)+
  scale_color_brewer(palette="Set2")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=16))

#UNIQUE
ssbcat_unique<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/ssGSEA_genepattern/ssGSEA_output/TALL_ALL_unique.gct")
ssbcat_unique<-as.data.frame(t(ssbcat_unique))
colnames(ssbcat_unique)<-ssbcat_unique[1,]
ssbcat_unique$method<-"ssGSEA"
ssbcat_unique<-ssbcat_unique[-c(1,2),]

metadata_ssbcat<-merge(ssbcat_unique,metadata,by=0,all.y=TRUE)
row.names(metadata_ssbcat)<-metadata_ssbcat$Row.names
metadata_ssbcat$Row.names<-NULL

metadata_ssbcat$bcat_relapse<-as.numeric(metadata_ssbcat$bcat_relapse)

cohorts<-c("TARGET","EGA_R","EGA_TLE")

ggplot(metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts, ],aes(x=First_event,y=bcat_relapse))+
  geom_point(data=metadata_ssbcat[metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),color="black",size=4)+
  geom_point(data=metadata_ssbcat[!is.na(metadata_ssbcat$First_event) & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],aes(color=First_event,shape=GSEA),size=3)+
  geom_boxplot(data=metadata_ssbcat[!is.na(metadata_ssbcat$First_event) & metadata_ssbcat$First_event!="Censored" & metadata_ssbcat$GSEA %in% cohorts,],alpha=0.1)+
  scale_color_brewer(palette="Set2")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.text = element_text(size=9),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=16,hjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=16))


# Empty space for subsequent data
rm(ty2)
rm(res.pca)
rm(res.pcaTLE)
rm(c1)

#####



##### bootstrapping ####
# dataset b-cat targets
#it8<-data_boot

relapses<-bcatexp_survival[bcatexp_survival$First_event=="Relapse" | bcatexp_survival$First_event=="Progression" | bcatexp_survival$First_event=="Death",]
remission<-bcatexp_survival[bcatexp_survival$First_event=="Remission",][sample(nrow(bcatexp_survival[bcatexp_survival$First_event=="Remission",]),50),]

data_boot<-rbind(relapses,remission)

# with the signature
phetsel<-pheatmap(as.data.frame(t(scale(data_boot[data_boot$First_event!="Second Malignant Neoplasm",colnames(data_boot) %in% unique(comp_all_tot$ensembl_gene_id)]))),
                  annotation_col = as.data.frame(data_boot[data_boot$First_event!="Second Malignant Neoplasm",c(1,5,7,8,14,17,19,(ncol(data_boot)-3):ncol(data_boot))]),
                  #         annotation_row = corout_brca[corout_brca$Gene %in% colnames(bcatexp_survival) & corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,][,3:4],
                  annotation_row = colgenes[,c(14,16)],
                  show_colnames=F,
                  show_rownames = F,
                  fontsize = 6,
                  cutree_cols =3,
                  cutree_rows = 6,
                  color=cols, breaks=myBreaks)

# all targets
phetsel<-pheatmap(as.data.frame(t(scale(data_boot[data_boot$First_event!="Second Malignant Neoplasm",22:(ncol(data_boot)-4)]))),
                  annotation_col = as.data.frame(data_boot[data_boot$First_event!="Second Malignant Neoplasm",c(1,5,7,8,14,17,19,(ncol(data_boot)-3):ncol(data_boot))]),
                  #         annotation_row = corout_brca[corout_brca$Gene %in% colnames(bcatexp_survival) & corout_brca$Pval<0.01 & abs(corout_brca$Rho)>0.2,][,3:4],
                  annotation_row = colgenes[,c(14,16)],
                  show_colnames=F,
                  show_rownames = F,
                  fontsize = 6,
                  cutree_cols =2,
                  cutree_rows = 6,
                  color=cols, breaks=myBreaks)

clusters<-data.frame(cluster=sort(cutree(phetsel$tree_row, k=6)))
patients<-data.frame(cluster=sort(cutree(phetsel$tree_col, k=2)))

clusters$Gene<-row.names(clusters)
clusters<-merge(clusters,geneid,by='Gene',all.x=T)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==2,]$NAME

# Order genes (clusters)
# all targets
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(data_boot[data_boot$First_event!="Second Malignant Neoplasm",22:(ncol(data_boot)-4)])))[phetsel$tree_row[["order"]],]))
#signature
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(data_boot[data_boot$First_event!="Second Malignant Neoplasm",colnames(data_boot) %in% unique(comp_all_tot$ensembl_gene_id)])))[phetsel$tree_row[["order"]],]))

geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
# all targets
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(data_boot[data_boot$First_event!="Second Malignant Neoplasm",22:(ncol(data_boot)-4)])))[,phetsel$tree_col[["order"]]]))
#signature
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(data_boot[data_boot$First_event!="Second Malignant Neoplasm",colnames(data_boot) %in% unique(comp_all_tot$ensembl_gene_id)])))[,phetsel$tree_col[["order"]]]))

patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(data_boot,patients,by="row.names")

ggplot(patinfo[patinfo$First_event!="Second Malignant Neoplasm",],aes(x=reorder(as.factor(cluster), bcat, FUN = median),y=scale(bcat)))+
  geom_point(aes(color=as.factor(First_event)),size=3)+
  geom_boxplot(alpha=0.2,lwd=1)+
  geom_hline(yintercept = median(scale(bcatexp_survival$bcat)),color="red",lwd=2)+
  geom_hline(yintercept = quantile(scale(bcatexp_survival$bcat))[4],color="grey",lwd=2)+
  geom_hline(yintercept = quantile(scale(bcatexp_survival$bcat))[2],color="grey",lwd=2)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, size=14,hjust=1),
        axis.text.y = element_text(angle=45, size=14,hjust=1))



## survival curves with all patients, REDO PATIENTS WITH HEATMAP

data_boot_cluster<-merge(patients,data_boot,by=0)
row.names(data_boot_cluster)<-data_boot_cluster$Row.names
data_boot_cluster$Row.names<-NULL


library(survminer)
data_boot_cluster$disease_free<-as.integer(data_boot_cluster$disease_free)
data_boot_cluster$Vital_status<-ifelse(data_boot_cluster$Vital_status=="Alive",0,1)
data_boot_cluster$Vital_status<-as.integer(data_boot_cluster$Vital_status)

# Create variable group of clusters
data_boot_cluster$newcluster<-NA
data_boot_cluster<-data_boot_cluster %>% 
  mutate(newcluster = case_when(cluster == 1 ~ 12,
                                cluster == 2 ~ 12,
                                cluster == 3 ~ 34,
                                cluster == 4 ~ 34
  ))


surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days/365, disease_free) ~ cluster, data=data_boot_cluster[data_boot_cluster$First_event!="Second Malignant Neoplasm",])

ggsurvplot(surv_C_fit_TALL, pval = TRUE,
           font.x = c(20),
           font.y = c(20))
#           palette = c("green", "indianred","lightseagreen","purple"))

pairwise_survdiff(Surv(Event_free_surv_days, disease_free) ~ cluster, data=data_boot_cluster)

surv_C_fit_TALL$time
#+Gender.x+Age_Dx+WBC_Dx+MRD_29+CNS_hr
fit.coxph<- coxph(Surv(Event_free_surv_days/365, disease_free) ~ cluster, 
                  data = data_boot_cluster,
                  ties = 'exact')
summary(fit.coxph)
 
ggforest(fit.coxph)


# bootstraping results 10 iterations
BS_pval<-data.frame(exp=c(rep("all",10),rep("microarrays_signature",10)),
           pval=c(0.75,0.25,0.97,0.4,0.64,0.76,0.038,0.37,0.53,0.72,0.76,0.80,0.48,0.25,0.16,0.6,0.84,0.74,0.94,0.79))

BS_pval_HR<-data.frame(exp=c(rep("all",10)),
                    pval=c(0.545,0.013,0.966,0.900,0.895,0.893,0.799,0.848,0.262,0.693),
                    HR=c(0.77,0.25,1,1.1,1.1,1.1,0.9,1.1,1.7,0.84))

BS_pval_HR<-melt(BS_pval_HR,id.vars="exp")

ggplot(BS_pval_HR,aes(x=value))+
  geom_density(aes(color=variable),size=2)+
  scale_color_brewer(palette="Set2")+
  geom_vline(xintercept = 0.05,linetype="dotted",color="red",size=1)+
  theme_bw()

#####