### Correlation DEGs

library(data.table)
library(ggplot2)
library(nlme)
library(readxl)
library(ggrepel)
library(ggnewscale)

setwd("/Volumes/grcmc/TERESA LJ/Experiments/Resistencias PDT005/2. RNAseq/") 

CTvsT<-read_xlsx("CvsT_GO_Rscript results/res_df_sig_rna.xlsx")
CTvsIC20<-read_xlsx("CvsIC20_GO_Rscript results/res_df_sig_rna20.xlsx")
CTvsIC30<-read_xls("CvsIC30_GO_Rscript results/res_df_sig_rna30.xls")

fetal<-read_xlsx("GSEA/Gene dataset signatures/Fetal/total_fetal.xlsx")

colnames(CTvsIC20)
CTvsIC20$FoldChange<-NULL

colnames(CTvsIC20)
CTvsIC20$exp<-"CTvsIC20"
colnames(CTvsIC30)
CTvsIC30$exp<-"CTvsIC30"
colnames(CTvsT)
CTvsT$exp<-"CTvsT"

colnames(fetal)<-c("log2FoldChange","external_gene_name","description","state")
fetal$exp<-"fetal"


# Merge organoids with fetal
IC20fet<-merge(CTvsIC20,fetal,by="external_gene_name")
IC30fet<-merge(CTvsIC30,fetal,by="external_gene_name")
Tfet<-merge(CTvsT,fetal,by="external_gene_name")

#merge IC20 vs IC30
ICsmer<-merge(CTvsIC20,CTvsIC30,by="id")


CTplot<-ggplot(Tfet,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  geom_point(color="black",size=3,alpha=0.5)+
  geom_density_2d(color="green")+
  geom_smooth(method="lm",color="red",se = FALSE)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = external_gene_name),color="black", data=Tfet[(Tfet$state=="UP" & Tfet$log2FoldChange.x>0),], fontface = "italic",size=5)+
  #xlim(-5,5)+
  #ylim(-5,5)+
  xlab("log2FC CT vs T")+
  ylab("log2FC fetal vs adult")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1)) 

CTplot

cor.test(Tfet$log2FoldChange.x,Tfet$log2FoldChange.y)


CIC20plot<-ggplot(IC20fet,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  geom_point(color="black",size=3,alpha=0.5)+
  geom_density_2d(color="green")+
  geom_smooth(method="lm",color="red",se = FALSE)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_label_repel(aes(label = external_gene_name),color="black", data=IC20fet[(IC20fet$state=="UP" & IC20fet$log2FoldChange.x>0),], fontface = "italic",size=5)+
  #xlim(-5,5)+
  #ylim(-5,5)+
  xlab("log2FC CT vs IC20")+
  ylab("log2FC fetal vs adult")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1)) 


cor.test(IC20fet$log2FoldChange.x,IC20fet$log2FoldChange.y)

CIC30plot<-ggplot(IC30fet,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  geom_point(color="black",size=3,alpha=0.5)+
  geom_density_2d(color="green")+
  geom_smooth(method="lm",color="red",se = FALSE)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_label_repel(aes(label = external_gene_name),color="black", data=IC30fet[(IC30fet$state=="UP" & IC30fet$log2FoldChange.x>0),], fontface = "italic",size=5)+
  #xlim(-5,5)+
  #ylim(-5,5)+
  xlab("log2FC CT vs IC30")+
  ylab("log2FC fetal vs adult")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1)) 


cor.test(IC30fet$log2FoldChange.x,IC30fet$log2FoldChange.y)

ICsplot<-ggplot(ICsmer,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  geom_point(color="black",size=3,alpha=0.5)+
  geom_density_2d(color="green")+
  geom_smooth(method="lm",color="red",se = FALSE)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_label_repel(aes(label = external_gene_name),color="black", data=ICsmer, fontface = "italic",size=5)+
  xlim(-7,7)+
  ylim(-7,7)+
  xlab("log2FC CT vs IC20")+
  ylab("log2FC CT vs IC30")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, hjust = 1))

ICsplot

cor.test(IC30fet$log2FoldChange.x,IC30fet$log2FoldChange.y)


library(gridExtra)
grid.arrange(CTplot,CIC30plot,CIC20plot)


### Senescence gene list from CellAge (https://genomics.senescence.info/cells/query.php?cell_type=&cell_line=&method=&cancer_type=&senescence_type=&search=&show=5&sort=1&page=1)
senes<-read.delim("/Volumes/grcmc/YGUILLEN/organ_treat/gene_lists/senescence_genes.txt")
colnames(senes)[1]<-"external_gene_name"

# senescence signature GSEA Teresa
SETE<-read_xlsx("GSEA/Gene dataset signatures/genesetpropio.xlsx")
SETE<-SETE[!is.na(SETE$senescence),]

seteme<-merge(SETE,senes,by="external_gene_name",all = TRUE)

