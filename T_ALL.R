## Expression profile of T-ALL human samples extracted from GEO ##

source("https://bioconductor.org/biocLite.R")
biocLite("affy")

library(biomaRt)
library(Biobase)
library(affy)
library(ggrepel)
library(dplyr)
library(data.table)

# Import datasets:

setwd("/Users/yguillen/Desktop/temp/beta_catenin_project/TALL_bcatenin_data/GEO_Expression/") 

# DS1: GSE62156
# DS2: GSE28703
# DS3: GSE26713
# DS4: GSE8879
# DS5: GSE110636
# Already have it in RNASeq DS6: GSE110633

DS1<-read.delim("GSE62156/GSE62156_series_matrix_R.txt",sep="\t")
DS1<-DS1[c(1,9,10,30:nrow(DS1)),]
rownames<-colnames(DS1)
row.names(DS1)<-DS1$Sample_title

metadata<-DS1[c(1:4),]
metadata<-as.data.frame(t(metadata))
metadata<-metadata[-1,]

DS1<-DS1[c(5:nrow(DS1)),]
DS1<-as.data.frame(t(DS1))
colnames(DS1)<-gsub('^','gene_',colnames(DS1))


DS1merge<-merge(metadata,DS1,by="row.names")
DS1merge$Sample_group<-gsub('molecular subgroup: ','',DS1merge$Sample_group)
DS1merge$Sample_group<-gsub(' subgroup','',DS1merge$Sample_group)

ensembl76_h <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
listgenes=c("CTNNB1","LEF1","TCF7","AXIN2","MYC")
getBM(attributes = c("external_gene_name",'hgnc_symbol', 'chromosome_name',
                     'start_position', 'end_position','affy_hg_u133_plus_2','ensembl_gene_id','entrezgene','ensembl_transcript_id_version'),
      filters = 'external_gene_name', 
      values = listgenes, 
      mart = ensembl76_h)

DS1merge$Sample_group[63]<-c("immature, HOXA13-t")

DS1merge<-subset(DS1merge,DS1merge$Sample_group!="CALM-AF10")
DS1merge<-subset(DS1merge,DS1merge$Sample_group!="HOX")
DS1merge<-subset(DS1merge,DS1merge$Sample_group!="HOXA of unknown mechanism")
DS1merge<-subset(DS1merge,DS1merge$Sample_group!="unknown, clusters with MLL but does not express HOXA genes")


DS1merge$sgroup<-ifelse(DS1merge$Sample_group == "immature" | DS1merge$Sample_group == "immature, HOXA13-t", 
                                           c("immature"), c("non_immature"))  


#Plot by class
DS1merge_test<-subset(DS1merge,select=c("Sample_group","sgroup","gene_1554411_at","gene_201533_at","gene_223679_at",
                                        "gene_210948_s_at","gene_221557_s_at","gene_221558_s_at",
                                        "gene_205254_x_at","gene_205255_x_at",
                                        "gene_222695_s_at","gene_222696_at","gene_224176_s_at","gene_224498_x_at",
                                        "gene_202431_s_at","gene_239931_at"))
colnames(DS1merge_test)<-c("Sample_group","sgroup","CTNNB1","CTNNB1","CTNNB1",
                           "LEF1","LEF1","LEF1",
                           "TCF7","TCF7",
                           "AXIN2","AXIN2","AXIN2","AXIN2",
                           "MYC","MYC")
DS1merge_test<-melt(DS1merge_test,id.vars=c("Sample_group","sgroup"))

ggplot(DS1merge_test,aes(x=sgroup,y=log(as.numeric(value))))+
  geom_violin(aes(color=sgroup),alpha=0,trim = FALSE)+
  geom_point()+
  geom_boxplot(width=0.3,alpha=0)+
  scale_color_manual(values=c("darkolivegreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "right")+
  facet_wrap(~variable)+
  ggtitle("DS1")+
  ylab("log norm exp")



## Download statistic test GEO2R
DS1_geo2R<-read.delim("GSE62156/GSE62156_GEO2R_padj01.txt",sep="\t")
DS1_geo2R_sel<-subset(DS1_geo2R,DS1_geo2R$Gene.symbol=="CTNNB1" | DS1_geo2R$Gene.symbol=="TCF7" | DS1_geo2R$Gene.symbol=="LEF1"|  DS1_geo2R$Gene.symbol=="AXIN2" | DS1_geo2R$Gene.symbol=="MYC")



# 2nd dataset
DS2<-read.delim("GSE28703/GSE28703_series_matrix_R.txt",sep="\t")
rownames<-colnames(DS2)
row.names(DS2)<-DS2$Sample_title

metadata<-DS2[c(1:4),]
metadata<-as.data.frame(t(metadata))
metadata<-metadata[-1,]

DS2<-DS2[c(5:nrow(DS2)),]
DS2<-as.data.frame(t(DS2))
colnames(DS2)<-gsub('^','gene_',colnames(DS2))


DS2merge<-merge(metadata,DS2,by="row.names")
DS2merge$Sample_source_name_ch1
DS2merge$Sample_group<-DS2merge$Sample_source_name_ch1

DS2merge$sgroup<-ifelse(DS2merge$Sample_group == "early T-cell precursor acute lymphoblastic leukemia", 
                        c("ETP"), c("non_ETP"))  

#CTNNB1: 201533_PM_at
#LEF1: 210948_PM_s_at, 221557_PM_s_at, 221558_PM_s_at
#TCF7: 205254_PM_x_at, 205255_PM_x_at
#AXIN2: gene_222695_PM_s_at, gene_222696_PM_s_at, gene_224176_PM_s_at,gene_224498_PM_x_at
#MYC: gene_202431_PM_s_at, gene_239931_PM_s_at

#Plot by class
DS2merge_test<-subset(DS2merge,select=c("sgroup","gene_1554411_PM_at","gene_201533_PM_at","gene_223679_PM_at",
                                        "gene_210948_PM_s_at","gene_221557_PM_s_at","gene_221558_PM_s_at",
                                        "gene_205254_PM_x_at","gene_205255_PM_x_at",
                                        "gene_222695_PM_s_at", "gene_222696_PM_at", "gene_224176_PM_s_at","gene_224498_PM_x_at",
                                        "gene_202431_PM_s_at","gene_239931_PM_at"
                                        ))
colnames(DS2merge_test)<-c("sgroup","CTNNB1","CTNNB1","CTNNB1",
                           "LEF1","LEF1","LEF1",
                           "TCF7","TCF7",
                           "AXIN2","AXIN2","AXIN2","AXIN2",
                           "MYC","MYC")
DS2merge_test<-melt(DS2merge_test,id.vars="sgroup")

ggplot(DS2merge_test,aes(x=sgroup,y=log(as.numeric(value))))+
  geom_violin(aes(color=sgroup),alpha=0,trim = FALSE)+
  geom_point()+
  geom_boxplot(width=0.3,alpha=0)+
  scale_color_manual(values=c("darkolivegreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "none")+
  facet_wrap(~variable)+
  ggtitle("DS2")+
  ylab("log norm exp")


## Download statistic test GEO2R
DS2_geo2R<-read.delim("GSE28703/GSE28703_GEO2R.txt",sep="\t")
DS2_geo2R_sel<-subset(DS2_geo2R,DS2_geo2R$Gene.symbol=="CTNNB1" | DS2_geo2R$Gene.symbol=="TCF7" | DS2_geo2R$Gene.symbol=="LEF1" |DS2_geo2R$Gene.symbol=="AXIN2" | DS2_geo2R$Gene.symbol=="MYC")


# 3rd dataset, two parts
DS3<-read.delim("GSE14618/GSE14618_series_matrix_R.txt",sep="\t")
rownames<-colnames(DS3)
row.names(DS3)<-DS3$Sample_geo_accession

metadata<-DS3[c(1:5),]
metadata<-as.data.frame(t(metadata))
metadata<-metadata[-1,]

DS3<-DS3[c(6:nrow(DS3)),]
DS3<-as.data.frame(t(DS3))
colnames(DS3)<-gsub('^','gene_',colnames(DS3))


DS3merge<-merge(metadata,DS3,by="row.names")
table(DS3merge$Sample_source_name_ch1)
DS3merge$Sample_group<-DS3merge$Sample_source_name_ch1


#Plot by class
DS3merge_test<-subset(DS3merge,select=c("Sample_group","gene_1554411_at","gene_201533_at","gene_223679_at",
                                        "gene_210948_s_at","gene_221557_s_at","gene_221558_s_at",
                                        "gene_205254_x_at","gene_205255_x_at",
                                        "gene_222695_s_at","gene_222696_at","gene_224176_s_at","gene_224498_x_at",
                                        "gene_202431_s_at","gene_239931_at" ))
colnames(DS3merge_test)<-c("sgroup","CTNNB1","CTNNB1","CTNNB1",
                           "LEF1","LEF1","LEF1",
                           "TCF7","TCF7",
                           "AXIN2","AXIN2","AXIN2","AXIN2",
                           "MYC","MYC")
DS3merge_test<-melt(DS3merge_test,id.vars="sgroup")

DS3merge_test$sgroup<-gsub('Primary T-ALL cells; newly diagnosed patient; ','',DS3merge_test$sgroup)

ggplot(DS3merge_test,aes(x=sgroup,y=log(as.numeric(value))))+
  geom_violin(aes(color=sgroup),alpha=0,trim = FALSE)+
  geom_point()+
  geom_boxplot(width=0.3,alpha=0)+
  #scale_color_manual(values=c("darkolivegreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "none")+
  facet_wrap(~variable)+
  ggtitle("DS3")+
  ylab("log norm exp")

## Download statistic test GEO2R
DS3_geo2R<-read.delim("GSE28703/GSE28703_GEO2R.txt",sep="\t")
DS3_geo2R_sel<-subset(DS3_geo2R,DS3_geo2R$Gene.symbol=="CTNNB1" | DS3_geo2R$Gene.symbol=="TCF7" | DS3_geo2R$Gene.symbol=="LEF1" | DS3_geo2R$Gene.symbol=="AXIN2" | DS3_geo2R$Gene.symbol=="MYC")


# 3rd dataset, two parts
DS3b<-read.delim("GSE14618/GSE14618_GPL96series_matrix_R.txt",sep="\t")
rownames<-colnames(DS3b)
row.names(DS3b)<-DS3b$Sample_geo_accession
  
metadata<-DS3b[c(1:3),]
metadata<-as.data.frame(t(metadata))
metadata<-metadata[-1,]

DS3b<-DS3b[c(4:nrow(DS3b)),]
DS3b<-as.data.frame(t(DS3b))
colnames(DS3b)<-gsub('^','gene_',colnames(DS3b))


DS3bmerge<-merge(metadata,DS3b,by="row.names")
table(DS3bmerge$Sample_source_name_ch1)
DS3bmerge$Sample_group<-DS3bmerge$Sample_source_name_ch1


#Plot by class
DS3bmerge_test<-subset(DS3bmerge,select=c("Sample_group","gene_201533_at",
                                          "gene_210948_s_at","gene_221557_s_at","gene_221558_s_at",
                                          "gene_205254_x_at","gene_205255_x_at",
                                          "gene_202431_s_at"))
colnames(DS3bmerge_test)<-c("sgroup","CTNNB1","LEF1","LEF1",
                            "LEF1",
                            "TCF7","TCF7",
                            "MYC")
DS3bmerge_test<-melt(DS3bmerge_test,id.vars="sgroup")

DS3bmerge_test$sgroup<-gsub('Primary T-ALL cells; newly diagnosed patient; ','',DS3bmerge_test$sgroup)

ggplot(DS3bmerge_test,aes(x=sgroup,y=log(as.numeric(value))))+
  geom_violin(aes(color=sgroup),alpha=0,trim = FALSE)+
  geom_point()+
  geom_boxplot(width=0.3,alpha=0)+
  scale_color_manual(values=c("coral","darkgreen"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "none")+
  facet_wrap(~variable)+
  ggtitle("DS3b")+
  ylab("log norm exp")



## Download statistic test GEO2R
DS3_geo2R<-read.delim("GSE28703/",sep="\t")
DS3_geo2R_sel<-subset(DS3_geo2R,DS3_geo2R$Gene.symbol=="CTNNB1" | DS3_geo2R$Gene.symbol=="TCF7" | DS3_geo2R$Gene.symbol=="LEF1" | DS3_geo2R$Gene.symbol=="MYC")

# 4rd dataset
DS4<-read.delim("GSE26713/GSE26713_series_matrix_R.txt",sep="\t")
rownames<-colnames(DS4)
row.names(DS4)<-DS4$Sample_geo_accession

metadata<-DS4[c(1:4),]
metadata<-as.data.frame(t(metadata))
metadata<-metadata[-1,]

DS4<-DS4[c(5:nrow(DS4)),]
DS4<-as.data.frame(t(DS4))
colnames(DS4)<-gsub('^','gene_',colnames(DS4))


DS4merge<-merge(metadata,DS4,by="row.names")
table(DS4merge$Sample_group)
DS4merge$Sample_group<-DS4merge$Sample_type

DS4merge$Sample_characteristics_ch1<-gsub('cell type: ','',DS4merge$Sample_characteristics_ch1)



#Plot by class
DS4merge_test<-subset(DS4merge,select=c("Sample_characteristics_ch1","Sample_group","gene_201533_at","gene_223679_at","gene_1554411_at",
                                        "gene_210948_s_at","gene_221557_s_at","gene_221558_s_at",
                                        "gene_205254_x_at","gene_205255_x_at",
                                        "gene_222695_s_at","gene_222696_at","gene_224176_s_at","gene_224498_x_at",
                                        "gene_202431_s_at","gene_239931_at"))
colnames(DS4merge_test)<-c("Type","sgroup","CTNNB1","CTNNB1","CTNNB1",
                           "LEF1","LEF1","LEF1",
                           "TCF7","TCF7",
                           "AXIN2","AXIN2","AXIN2","AXIN2",
                           "MYC","MYC")
DS4merge_test<-melt(DS4merge_test,id.vars=c("Type","sgroup"))

DS4merge_test$sgroup<-gsub('cytogenetics: ','',DS4merge_test$sgroup)

ggplot(DS4merge_test,aes(x=Type,y=log(as.numeric(value))))+
  geom_violin(aes(color=Type),alpha=0,trim = FALSE)+
  geom_point()+
  geom_boxplot(width=0.3,alpha=0)+
  scale_color_manual(values=c("darkgreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "none")+
  facet_wrap(~variable)+
  ggtitle("DS4")+
  ylab("log norm exp")


## Download statistic test GEO2R
DS4_geo2R<-read.delim("GSE26713/GSE26713_geo2R.txt",sep="\t")
DS4_geo2R_sel<-subset(DS4_geo2R,DS4_geo2R$Gene.symbol=="CTNNB1" | DS4_geo2R$Gene.symbol=="TCF7" | DS4_geo2R$Gene.symbol=="LEF1" | DS4_geo2R$Gene.symbol=="AXIN2" | DS4_geo2R$Gene.symbol=="MYC")

subset(DS4_geo2R,DS4_geo2R$Gene.symbol=="CTNNB1")
subset(DS4_geo2R,DS4_geo2R$Gene.symbol=="LEF1")
subset(DS4_geo2R,DS4_geo2R$Gene.symbol=="TCF7")
subset(DS4_geo2R,DS4_geo2R$Gene.symbol=="AXIN2")
subset(DS4_geo2R,DS4_geo2R$Gene.symbol=="MYC")

# 5th dataset
DS5<-read.delim("GSE8879/GSE8879_series_matrix_R.txt",sep="\t")
rownames<-colnames(DS5)
row.names(DS5)<-DS5$Sample_geo_accession

metadata<-DS5[c(1:2),]
metadata<-as.data.frame(t(metadata))
metadata<-metadata[-1,]

DS5<-DS5[c(3:nrow(DS5)),]
DS5<-as.data.frame(t(DS5))
colnames(DS5)<-gsub('^','gene_',colnames(DS5))


DS5merge<-merge(metadata,DS5,by="row.names")
table(DS5merge$Sample_characteristics_ch1)
DS5merge$Sample_group<-DS5merge$Sample_characteristics_ch1

DS5merge$Sample_group<-gsub('cell type: diagnostic leukemic blasts of early T-cell precursor acute lymphoblastic leukemia ','T-ALL',DS5merge$Sample_group)
DS5merge$Sample_group<-gsub('cell type: diagnostic leukemic blasts of T-cell precursor acute lymphoblastic leukemia','T-ALL(nonETP)',DS5merge$Sample_group)
table(DS5merge$Sample_group)



#Plot by class
DS5merge_test<-subset(DS5merge,select=c("Sample_group","gene_201533_at",
                                        "gene_210948_s_at","gene_221557_s_at","gene_221558_s_at",
                                        "gene_205254_x_at","gene_205255_x_at",
                                        "gene_202431_s_at"))
colnames(DS5merge_test)<-c("sgroup","CTNNB1",
                           "LEF1","LEF1","LEF1",
                           "TCF7","TCF7",
                           "MYC")
DS5merge_test<-melt(DS5merge_test,id.vars=c("sgroup"))

ggplot(DS5merge_test,aes(x=sgroup,y=log(as.numeric(value))))+
  geom_violin(aes(color=sgroup),alpha=0,trim = FALSE)+
  geom_point()+
  geom_boxplot(width=0.3,alpha=0)+
  scale_color_manual(values=c("darkgreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "none")+
  facet_wrap(~variable)+
  ggtitle("DS5")+
  ylab("log norm exp")



## Download statistic test GEO2R
DS5_geo2R<-read.delim("GSE8879/GSE8879_GEO2R.txt",sep="\t")
DS5_geo2R_sel<-subset(DS5_geo2R,DS5_geo2R$Gene.symbol=="CTNNB1" | DS5_geo2R$Gene.symbol=="TCF7" | DS5_geo2R$Gene.symbol=="LEF1" | DS5_geo2R$Gene.symbol=="AXIN2" | DS5_geo2R$Gene.symbol=="MYC")

subset(DS5_geo2R,DS5_geo2R$Gene.symbol=="CTNNB1")
subset(DS5_geo2R,DS5_geo2R$Gene.symbol=="LEF1")
subset(DS5_geo2R,DS5_geo2R$Gene.symbol=="TCF7")
subset(DS5_geo2R,DS5_geo2R$Gene.symbol=="MYC")

# E-MTAB-604

load("E-MTAB-604/E-MTAB-604.eSet.r")
#Normalize expression data
eset <- rma(study)
write.exprs(eset,file="data.txt")

featureNames(eset)
dim(eset)

# Read data and metadata
metadata<-read.delim("E-MTAB-604/samples.txt",sep="\t")
expdata<-read.delim("E-MTAB-604/data.txt",sep="\t")
row.names(expdata)<-expdata$X
expdata$X<-NULL

expdata<-as.data.frame(t(expdata))

row.names(expdata)
colnames(expdata)<-gsub('^','gene_',colnames(expdata))

row.names(expdata)
row.names(metadata)<-metadata$Hybridization.Name
DS6merge<-merge(metadata,expdata,by="row.names")

table(DS6merge$Characteristics..ClinicalInformation.)
DS6merge$Sample_group<-paste(DS6merge$Material.Type,DS6merge$Characteristics..OrganismPart.,DS6merge$Characteristics..DiseaseState.)
DS6merge$Sample_group<-gsub('^cell.*','Cell line T-ALL',DS6merge$Sample_group)
DS6merge$Sample_group<-gsub('organism_part blood.*','blood T-ALL',DS6merge$Sample_group)
DS6merge$Sample_group<-gsub('organism_part bone marrow normal','BM normal',DS6merge$Sample_group)
DS6merge$Sample_group<-gsub('organism_part bone marrow T.*','BM T-ALL',DS6merge$Sample_group)
DS6merge$Sample_group<-gsub('organism_part xenograft.*','Xenograft T-ALL',DS6merge$Sample_group)

table(DS6merge$Sample_group)

DS6merge$Characteristics..CellLine.<-gsub('  ','',DS6merge$Characteristics..CellLine.)
DS6merge$Characteristics..CellLine.




#Plot by class
DS6merge_test<-subset(DS6merge,select=c("Row.names","Characteristics..CellLine.","Sample_group","gene_201533_at",
                                        "gene_210948_s_at","gene_221557_s_at","gene_221558_s_at",
                                        "gene_205254_x_at","gene_205255_x_at",
                                        "gene_222695_s_at","gene_222696_at","gene_224176_s_at","gene_224498_x_at",
                                        "gene_202431_s_at"))
colnames(DS6merge_test)<-c("Row.names","Cell_line","sgroup","CTNNB1",
                           "LEF1","LEF1","LEF1",
                           "TCF7","TCF7",
                           "AXIN2","AXIN2","AXIN2","AXIN2",
                           "MYC")
DS6merge_test<-melt(DS6merge_test,id.vars=c("Row.names","Cell_line","sgroup"))

DS6merge_test$sgroup<-gsub('^organism_part.*','NA',DS6merge_test$sgroup)
DS6merge_test<-subset(DS6merge_test,DS6merge_test$sgroup!="NA")

ggplot(DS6merge_test,aes(x=sgroup,y=log(as.numeric(value))))+
  geom_violin(aes(color=sgroup),alpha=0,trim = FALSE)+
  geom_point()+
  #geom_label_repel(aes(label=Cell_line),color="grey", fontface = "italic",size=2.5)+
  geom_boxplot(width=0.3,alpha=0)+
  #scale_color_manual(values=c("darkgreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "none")+
  facet_wrap(~variable)+
  ggtitle("DS6")+
  ylab("log norm exp")

tabcom<-compare_means(value ~ sgroup, data = DS6merge_test, 
              group.by = "variable")


ggboxplot(DS6merge_test, x = "sgroup", y = "value",
          color = "sgroup", palette = "jco",facet.by = "variable")+
  stat_compare_means(method = "anova")+stat_compare_means(label = "p.signif", method = "t.test",ref.group = "BM T-ALL") 
  


# Correlation bcat, LEF1, TCF7, MYC, AXIN2

#corgenes<-subset(DS6merge,select=c(Row.names,Sample_group,Characteristics..CellLine.,gene_201533_at,gene_210948_s_at,gene_221558_s_at,gene_221557_s_at,gene_205254_x_at,gene_205255_x_at))
#colnames(corgenes)[4:ncol(corgenes)]<-c("CTNNB1","LEF1","LEF1","LEF1","TCF7","TCF7")

colnames(DS6merge_test)
DS6merge_test$Cell_line<-gsub('^$','Primary T-ALL',DS6merge_test$Cell_line)
DS6merge_test$group<-paste(DS6merge_test$Cell_line,DS6merge$Sample_group)

DS6merge_test$id<-paste(DS6merge_test$Row.names,DS6merge_test$Cell_line)
DS6merge_test$id<-gsub(' Primary.*','',DS6merge_test$id)
DS6merge_test$id<-gsub('_U133_2','',DS6merge_test$id)

DS6merge_test$value<-as.numeric(DS6merge_test$value)

ggplot(DS6merge_test, aes(id,as.character(variable))) +
  geom_tile(aes(fill = log(value)), color = "white") +
  #scale_fill_manual(values = c("darkolivegreen","coral3")) +
  theme_minimal()+
  scale_fill_gradient2(low="white", mid="grey", high="blue",midpoint=1.5,limits=range(log(DS6merge_test$value))) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=8,face="bold"),
        axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=4, hjust = 1)) +
  labs(fill = "Expression level norm")+
  facet_wrap(~sgroup,scales="free",ncol=3)+
  coord_flip()


ggplot(subset(DS6merge_test,DS6merge_test$sgroup!="BM normal"),aes(log(value),colour=variable))+
  geom_density()+
  theme_bw()+
  facet_grid(~sgroup)


# select only paired T-ALL and xenografted
cogenes_melt_sel<-as.data.frame(subset(cogenes_melt,
                            cogenes_melt$Row.names=="JSR_129_U133_2" |
                            cogenes_melt$Row.names=="JSR_130_U133_2" | 
                            cogenes_melt$Row.names=="JSR_131_U133_2" |
                            cogenes_melt$Row.names=="JSR_138_U133_2" |
                            cogenes_melt$Row.names=="JSR_139_U133_2" |
                            cogenes_melt$Row.names=="JSR_140_U133_2" |
                            cogenes_melt$Row.names=="JSR_141_U133_2" |
                            cogenes_melt$Row.names=="JSR_142_U133_2" |
                            cogenes_melt$Row.names=="JSR_143_U133_2" |
                            cogenes_melt$Row.names=="JSR_112_U133_2" |
                            cogenes_melt$Row.names=="JSR_115_U133_2" |
                            cogenes_melt$Row.names=="JSR_127_U133_2" |
                            cogenes_melt$Row.names=="JSR_036_U133_2" |
                            cogenes_melt$Row.names=="JSR_034_U133_2" |
                            cogenes_melt$Row.names=="JSR_077_U133_2" |
                            cogenes_melt$Row.names=="JSR_079_U133_2" |
                            cogenes_melt$Row.names=="JSR_088_U133_2" |
                            cogenes_melt$Row.names=="JSR_051_U133_2"))

cogenes_melt_sel$Row.names
cogenes_melt_sel$Pair<-c("Pair5","Pair4","Pair9","Pair6","Pair7","Pair8","Pair1","Pair2","Pair3","Pair1","Pair2","Pair3","Pair4",
                         "Pair5","Pair6","Pair7","Pair8","Pair9","Pair5","Pair4","Pair9","Pair6","Pair7","Pair8","Pair1","Pair2","Pair3","Pair1","Pair2","Pair3","Pair4",
                         "Pair5","Pair6","Pair7","Pair8","Pair9", "Pair5","Pair4","Pair9","Pair6","Pair7","Pair8","Pair1","Pair2","Pair3","Pair1","Pair2","Pair3","Pair4",
                         "Pair5","Pair6","Pair7","Pair8","Pair9")

cogenes_melt_sel$group<-gsub(' blood.*','',cogenes_melt_sel$group)
cogenes_melt_sel$group<-gsub(' BM.*','',cogenes_melt_sel$group)
cogenes_melt_sel$group<-gsub('^.*Xenograft','Xenograft',cogenes_melt_sel$group)
#change name
cogenes_melt_sel$type<-cogenes_melt_sel$group


cogenes_melt_sel$gene<-paste(cogenes_melt_sel$Pair,cogenes_melt_sel$variable)
cogenes_melt_sel$class<-paste(cogenes_melt_sel$Pair,cogenes_melt_sel$group)

ggplot(cogenes_melt_sel,aes(x=type,y=value,group=Pair))+
  #geom_violin(aes(linetype=type),alpha=0,trim = FALSE)+
  geom_point(aes(color=variable))+
  geom_boxplot(aes(linetype=class),width=0.3,alpha=0)+
  geom_path(aes(group=gene))+
  #scale_color_manual(values=c("darkolivegreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1))+
  facet_grid(~Pair)
  #ggtitle("TCF7")+
  #ylab("log norm exp TCF7")


# RNASeq data

#GSE110633

metadata<-read.delim('GSE110633/GSE110633_RAW/input_R/info_samples.txt',sep="\t")
metadata<-as.data.frame(t(metadata))
colnames<-metadata[1,]
metadata<-as.data.frame(metadata[-1,])
colnames(metadata)<-colnames
row.names(metadata)

ff <- list.files( path = "GSE110633/GSE110633_RAW/raw_counts/", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
#Read each file, skipping first 4 lines
counts.files <- lapply( ff, read.table, skip = 4 )
#Select 2nd column as counts
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )
ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/counts/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1

colnames(counts)<-gsub('^GSE.*GSM','GSM',colnames(counts))
colnames(counts)<-gsub('_Reads.*','',colnames(counts))
colnames(counts)<-gsub('_sample.*','',colnames(counts))
counts<-as.matrix(counts)

normcounts <- rlogTransformation(counts, blind=TRUE) ## create a distance matrix between the samples
colnames(normcounts)
row.names(normcounts)<-row.names(counts)

t_normcounts<-t(normcounts)
datanorm<-merge(metadata,t_normcounts,by="row.names")

library(DataCombine)

datanorm<-melt(datanorm,id.vars=c("Row.names","Sample_ID"))
datanorm$variable<-as.character(datanorm$variable)

# Search pattern within variable
#grepl.sub(data=datanorm,pattern="ENSG00000168646",Var="variable")
          
subdata<-subset(datanorm,datanorm$variable=="ENSG00000168036.16" | 
                  datanorm$variable=="ENSG00000138795.9" | 
                  datanorm$variable=="ENSG00000081059.19" | 
                  datanorm$variable=="ENSG00000136997.15" |
                  datanorm$variable=="ENSG00000168646.12")



subdata$variable<-gsub('ENSG00000168036.16','CTNNB1',subdata$variable)
subdata$variable<-gsub('ENSG00000138795.9','LEF1',subdata$variable)
subdata$variable<-gsub('ENSG00000081059.19','TCF7',subdata$variable)
subdata$variable<-gsub('ENSG00000136997.15','MYC',subdata$variable)
subdata$variable<-gsub('ENSG00000168646.12','AXIN2',subdata$variable)

ggplot(subdata,aes(x=Sample_ID,y=value))+
  geom_violin(aes(color=Sample_ID),alpha=0,trim = FALSE)+
  geom_point()+
  #geom_label_repel(aes(label=Cell_line),color="grey", fontface = "italic",size=2.5)+
  geom_boxplot(width=0.3,alpha=0)+
  #scale_color_manual(values=c("darkgreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "none")+
  facet_wrap(~variable)+
  ggtitle("GSE110633")+
  ylab("log norm exp")


tabcom2<-compare_means(value ~ Sample_ID, data = subdata, 
                      group.by = "variable")


ggboxplot(subdata, x = "Sample_ID", y = "value",
          color = "Sample_ID", palette = "jco",facet.by = "variable")+
  stat_compare_means(method = "anova")+stat_compare_means(label = "p.signif", method = "t.test",ref.group = "IMM") 



#GSE110636

metadata<-read.delim('GSE110636/GSE110636_RAW/input_R/info_samples.txt',sep="\t")
metadata<-as.data.frame(t(metadata))
colnames<-metadata[1,]
metadata<-as.data.frame(metadata[-1,])
colnames(metadata)<-colnames
row.names(metadata)

ff <- list.files( path = "GSE110636/GSE110636_RAW/raw_counts/", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
ff
#Read each file, skipping first 4 lines
counts.files <- lapply( ff, read.table, skip = 4 )
#Select 2nd column as counts
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )
ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/counts/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1

colnames(counts)<-gsub('^GSE.*GSM','GSM',colnames(counts))
colnames(counts)<-gsub('_Reads.*','',colnames(counts))
colnames(counts)<-gsub('_sample.*','',colnames(counts))
counts<-as.matrix(counts)

normcounts <- rlogTransformation(counts, blind=TRUE) ## create a distance matrix between the samples
colnames(normcounts)
row.names(normcounts)<-row.names(counts)

t_normcounts<-t(normcounts)
datanorm<-merge(metadata,t_normcounts,by="row.names")

datanorm<-melt(datanorm,id.vars=c("Row.names","Sample_ID"))
datanorm$variable<-as.character(datanorm$variable)
subdata<-subset(datanorm,datanorm$variable=="ENSG00000168036.16" | 
                  datanorm$variable=="ENSG00000138795.9" | 
                  datanorm$variable=="ENSG00000081059.19" |
                  datanorm$variable=="ENSG00000136997.15" |
                  datanorm$variable=="ENSG00000168646.12")


subdata$variable<-gsub('ENSG00000168036.16','CTNNB1',subdata$variable)
subdata$variable<-gsub('ENSG00000138795.9','LEF1',subdata$variable)
subdata$variable<-gsub('ENSG00000081059.19','TCF7',subdata$variable)
subdata$variable<-gsub('ENSG00000136997.15','MYC',subdata$variable)
subdata$variable<-gsub('ENSG00000168646.12','AXIN2',subdata$variable)



subdata$sgroup<-ifelse(subdata$Sample_ID == "IMM" , 
                        c("immature"), c("non_immature")) 

ggplot(subdata,aes(x=sgroup,y=value))+
  geom_violin(aes(color=sgroup),alpha=0,trim = FALSE)+
  geom_point()+
  #geom_label_repel(aes(label=Cell_line),color="grey", fontface = "italic",size=2.5)+
  geom_boxplot(width=0.3,alpha=0)+
  #scale_color_manual(values=c("darkgreen","coral"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=8, hjust = 1),
        legend.position = "none")+
  facet_wrap(~variable)+
  ggtitle("GSE110636")+
  ylab("log norm exp")


tabcom3<-compare_means(value ~ Sample_ID, data = subdata, 
                       group.by = "variable")


ggboxplot(subdata, x = "Sample_ID", y = "value",
          color = "Sample_ID", palette = "jco",facet.by = "variable")+
  stat_compare_means(method = "anova")+stat_compare_means(label = "p.signif", method = "t.test",ref.group = "IMM") 


ggboxplot(subdata, x = "sgroup", y = "value",
          color = "sgroup", palette = "jco",facet.by = "variable")+
  stat_compare_means()+stat_compare_means(label = "p.signif", method = "t.test") 


