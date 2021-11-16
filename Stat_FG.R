## Expression profile of T-ALL human samples extracted from GEO ##

source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("factoextra")

library(biomaRt)
library(Biobase)
library(affy)
library(ggrepel)
library(dplyr)
library(data.table)
library(factoextra)

# Import datasets:

setwd("/Volumes/grcmc/YGUILLEN/STAT_paper_Gallardo/") 

DS<-read.delim("GSEA/CT_mi124_all.txt")
rownames<-colnames(DS)

DSt<-as.data.frame(t(DS))
colnames(DSt)

colnames(DSt) <- as.character(unlist(DSt[1,]))
DSt=DSt[-1, ]

ensembl76 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#genelist<-colnames(DSt[,1:(ncol(DSt)-1)])
#genenames<-getBM(attributes = c("external_gene_name",'hgnc_symbol', 'chromosome_name',
#                     'start_position', 'end_position','affy_hg_u133_plus_2','ensembl_gene_id','entrezgene'),
      #filters = 'external_gene_name', 
#      values = genelist, 
#      mart = ensembl76)


#mergenames
#genemat<-merge(genenames,DS,by.x="affy_hg_u133_plus_2",by.y="row.names")

DSmat <- data.frame(sapply(DSt, function(x) as.numeric(as.character(x))))
row.names(DSmat)<-row.names(DSt)
class(DSmat$ABCA1.1)


#select those with variance != 0 (error in prcomp otherwise)
cond<-(apply(DSmat, 2, var)!=0)
DSmat<-DSmat[, cond, drop = FALSE]

#PCA
res.pca<-prcomp(DSmat,scale=TRUE)
summary(res.pca)

#Proportion of variance:
# PC1: 49.53 %
# PC2: 29.08 %

pcatab<-as.data.frame(res.pca$x)
pcatab$condition<-c("Control","Control","miR124","miR124")
pcatab$names<-row.names(pcatab)

ggplot(pcatab, aes(PC1, PC2, color=condition,label=names)) +
  geom_point(size=8) +
  #stat_ellipse(aes(color=condition),level=0.8)+
  scale_color_manual(values=c("dodgerblue","coral"))+
  #scale_color_manual(values=c("grey","dodgerblue","coral","darkolivegreen2","black","dodgerblue4","coral4"))+
  geom_label_repel(aes(label=names), size=4,fontface = "italic",alpha=0.7)+
  #scale_color_manual(values=c("purple","red"))+
  #ylim(-50,50)+
  ylab("PC2 (29.08% variance)")+
  xlab("PC1 (49.53% variance)")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  coord_fixed()


### FUNCTIONAL ENRICHMENT ###
FEn<-read_xlsx("CtrlvsmiR124.p05%20Filtrada.290618.xlsx")
FEn<-subset(FEn,select=c("Symbol","logFC_CT_miR124","Pval"))

FEn_filt<-subset(FEn,FEn$Pval<=0.05)
FEn_filt<-subset(FEn_filt,abs(FEn_filt$logFC_CT_miR124)>=0.6)
FEn_filt<-subset(FEn_filt,FEn_filt$Symbol!="---")

Up_CT_genes<-(subset(FEn_filt,FEn_filt$logFC_CT_miR124>0))$Symbol
Down_CT_genes<-(subset(FEn_filt,FEn_filt$logFC_CT_miR124<0))$Symbol

## GO terms and pathways
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018","TF_Perturbations_Followed_by_Expression","Reactome_2015")

enrichedUp_CT <- enrichr(Up_CT_genes, dbs)
enrichedDown_CT <- enrichr(Down_CT_genes, dbs)

Up_CT_enrich <- enrichedUp_CT[["GO_Biological_Process_2018"]]
Up_CT_enrich$group<-c("Up_CT")

Down_CT_enrich <- enrichedDown_CT[["GO_Biological_Process_2018"]]
Down_CT_enrich$group<-c("Down_CT")

allGO<-rbind(Up_CT_enrich,Down_CT_enrich )

bpsub<-subset(allGO,allGO$Adjusted.P.value<0.1)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term<-as.factor(bpsub$Term)
bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$P.value))])

bpsub$Term <-reorder(bpsub$Term, rev(bpsub$P.value))

bpsub$Overlap<-gsub('/.*','',bpsub$Overlap)
bpsub<-subset(bpsub,bpsub$Overlap>=2)

ggplot(bpsub,aes(x=Term,y=-(log(P.value)),fill=group))+
  geom_bar(position="dodge",stat = "identity",color="black")+
  geom_text(aes(label=tolower(Genes)), size=3,hjust=1)+
  theme_minimal()+
  scale_fill_brewer(palette = "Blues")+
  xlab("")+
  ylab("-log (P-Value)")+
  theme(axis.text.y =element_text(size=6),axis.title=element_text(size=8,face="bold"),axis.text.x = element_text(size=12))+
  coord_flip()

### GSEA heatmaps ###
## Hallmark databse
myc1<-read.delim("GSEA/R_GSEA/HALLMARK_MYC_TARGETS_V1.txt")
myc1<-subset(myc1,select=c(PROBE,CORE.ENRICHMENT))
myc1$PATH<-"MYC1_targets"
colnames(myc1)<-c("NAME","CORE","PATH")

length(unique(DS$NAME))

genesel<-merge(myc1,DS,by="NAME")
genesel<-genesel[,c(1,3,4,5,6,7)]

require(scales)
genesel_path<-as.data.frame(t(apply(genesel[, -c(1,2)], 1, rescale)))


genesel_path$NAME<-genesel$NAME
genesel_path$PATH<-genesel$PATH

genesel_path<-melt(genesel_path,id.vars=c("NAME","PATH"))
genesel_path$NAME<-as.character(genesel_path$NAME)
genesel_path$gene<-paste(genesel_path$NAME,row.names(genesel_path),sep = "_")

ggplot(genesel_path, aes(x=gene,y=variable)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="grey60", high="red",mid="white",midpoint = 0.5)+
  xlab("Genes") +
  #ylab("Batch") +
  theme_minimal()+
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=8,face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6, hjust = 1,face="italic")) +
  labs(fill = "Expression level fold change")+
  facet_wrap(~PATH,scales="free",ncol=2)+
  xlab("")+
  ylab("")+
  coord_flip()


