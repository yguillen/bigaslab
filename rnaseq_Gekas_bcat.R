## Yolanda Guilén ##
## October 2018 ##

##################### RNA-seq pipeline ####################
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("DESeq2")
biocLite("cummeRbund")
biocLite("pasilla")

install.packages("htmlwidgets")
install.packages("readxl")

#basic libraries
library(data.table)
library(ggplot2)
library(nlme)
library(readxl)
#library(pasilla)

## Get ids from ensembl all genes
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version=77)
chr_genes <- getBM(attributes=c('ensembl_gene_id',
                                 'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)

library(htmlwidgets)
library(DESeq2)
library(cummeRbund)

##################### Pipeline bash form tophat alignments Christos data ############

### Tophat alignments at:
# /Users/instalar/Desktop/RNA_seq/Christos/tophat/bam_files/bamfiles_Tophat
# parameters used: 
# tophat 
#-p 4 num threads
#-a 8 anchor length
# -m 0 splice-missmatches
# -i 70 min anchor length
# -I 500000 MAx intron length
#-g 2 max multihits
#--bowtie-n 
#--library-type fr-unstranded library unstranded
#--max-insertion-length 3 
#--max-deletion-length 3 
#--coverage-search coverage based search for junctions
#--min-coverage-intron 50 
#--max-coverage-intron 20000 
#-o /users/fk/cbello/CGekas_RNAseq/tophat_out/R_Samples/2 
# /software/fk/el6.3/bowtie-1.0.1/indexes/Homo_sapiens.GRCh38.dna.ebwt 
#/users/fk/cbello/CGekas_RNAseq/R_8204_ACTTGA.fastq.trimmed.qual_filtered #fastq file


## Summary alignments
sum_align<-read.delim('/Users/instalar/Desktop/RNA_seq/Christos/tophat/bam_files/summary_tophat_alignments/summary_alignments.txt',sep="\t")
sum_align$Exp<-NULL


## cufflinks on each separated line and condition

## cuffmerge of all transcripts from each line
#cuffmerge 
#-p 8 
#-o /users/fk/cbello/CGekas_RNAseq/tophat_out/cuffmerge/assembly_list.merged 
#-g /users/fk/cbello/CGekas_RNAseq/Human_annotation_files/Homo_sapiens.GRCh38.77.gtf 
#-s /users/fk/cbello/CGekas_RNAseq/Human_annotation_files/Homo_sapiens.GRCh38.dna.fa 
#/users/fk/cbello/CGekas_RNAseq/tophat_out/cuffmerge/assembly_list.txt

## cuffdiff for each paired condition line
#cuffdiff 
#-p 16 
#-o /users/fk/cbello/CGekas_RNAseq/tophat_out/cuffdiff/ 
#-L C_Samples_1,C_Samples_2 
#-b /users/fk/cbello/CGekas_RNAseq/Human_annotation_files/Homo_sapiens.GRCh38.dna.fa 
#-u --library-type fr-unstranded 
#/users/fk/cbello/CGekas_RNAseq/tophat_out/cuffmerge/assembly_list.merged/merged.gtf 
#/users/fk/cbello/CGekas_RNAseq/tophat_out/C_Samples/1/accepted_hits.bam 
#/users/fk/cbello/CGekas_RNAseq/tophat_out/C_Samples/2/accepted_hits.bam 



## importing cuffdiff results from Gekas analyses DMSO vs PFK##
# LINE RPMI
RPMI_dif<-read.delim('/Volumes/grcmc/YGUILLEN/Christos_RNASeq/Christos_pkf/cuffdiff/R_Samples/gene_exp.diff',sep="\t")
RPMI_dif$line<-c("RPMI")
#line RPMI conditions 1 and 2 need to be exchange to match with the rest of lines. Change signal log2foldchange
RPMI_dif$log2.fold_change.<-RPMI_dif$log2.fold_change.*(-1)
RPMI_dif <- RPMI_dif[ , c("test_id", "gene_id", "gene", "locus", "sample_2", "sample_1", "status", "value_2", "value_1", "log2.fold_change.", "test_stat", "p_value", "q_value", "significant", "line")]
names(RPMI_dif)[8:9]<-c("value_1","value_2")
names(RPMI_dif)[5:6]<-c("sample_1","sample_2")

# LINE CEM
CEM_dif<-read.delim('/Volumes/grcmc/YGUILLEN/Christos_RNASeq/Christos_pkf/cuffdiff/C_Samples/gene_exp.diff',sep="\t")
CEM_dif$line<-c("CEM")

# LINE HBP
HBP_dif<-read.delim('/Volumes/grcmc/YGUILLEN/Christos_RNASeq/Christos_pkf/cuffdiff/H_Samples/gene_exp.diff',sep="\t")
HBP_dif$line<-c("HBP")

# LINE JURKAT
JUR_dif<-read.delim('/Volumes/grcmc/YGUILLEN/Christos_RNASeq/Christos_pkf/cuffdiff/J_Samples/gene_exp.diff',sep="\t")
JUR_dif$line<-c("JUR")

# merge all lines
all_bind<-rbind(RPMI_dif,CEM_dif,HBP_dif,JUR_dif)
all_bind$chromosome<-gsub(':.*','',all_bind$locus)


### Overall expression levels ##
# exclude expression levels higher than 500
all_bind_fpkm<-subset(all_bind,all_bind$value_1<=600 & all_bind$value_2<=600)

ggplot(all_bind_fpkm,aes(x=log(value_1+0.001),y=log(value_2+0.001)))+
  geom_point(size=1)+
  facet_wrap(~line,scales="free")+
  geom_smooth(method='lm',color="blue")+
  geom_abline(intercept = 0, slope = 1,color="red")+
  theme_minimal()+
  xlab("DMSO")+
  ylab("PKF")

lmList(value_2 ~ value_1 | gene , data=all_bind_fpkm, na.action=na.exclude)


### subset significant FDR <= 0.1
sig_exp<-subset(all_bind,all_bind$q_value<=0.1)
sig_exp$linegene<-paste(sig_exp$line,sig_exp$gene,sep="_")
sig_exp<-sig_exp[order(sig_exp$gene),]

#common matches
common<-as.data.frame(table(sig_exp$gene))
names(common)<-c("gene","Freq")
common<-subset(common,common$Freq>=1)
common<-common[rev(order(common$Freq)),]
list_genes<-common$gene

#extract common genes from all 4 lines datasets

all_sel<-subset(all_bind, gene %in% common$gene)
all_sel<-subset(all_sel,select=c(line,gene,log2.fold_change.,q_value,p_value))
all_sel<-merge(all_sel,common,by="gene")


### subset RPMI FDR <= 0.1
sig_exp<-subset(RPMI_dif,RPMI_dif$q_value<=0.1)
sig_exp<-sig_exp[order(sig_exp$gene),]
sig_exp[,10][(sig_exp[,10]== -Inf)] <- min(sig_exp[,10][sig_exp[,10] != -Inf])
sig_exp[,10][(sig_exp[,10]== Inf)] <- max(sig_exp[,10][sig_exp[,10] != Inf])
all_sel<-sig_exp


#substitute infinite values by maximum and minimum
all_sel[,10][(all_sel[,10]== -Inf)] <- min(all_sel[,10][all_sel[,10] != -Inf])
all_sel[,10][(all_sel[,10]== Inf)] <- max(all_sel[,10][all_sel[,10] != Inf])

#heatmap
all_sel<-all_sel[order(all_sel$log2.fold_change.), ]
all_sel$gene <-reorder(all_sel$gene, all_sel$log2.fold_change.)

ggplot(subset(all_sel,all_sel$gene!="Y_RNA"), aes(gene, line )) +
  geom_tile(aes(fill = log2.fold_change.), color = "white") +
  geom_point(aes(size=(-log(q_value))),data = subset(all_sel,all_sel$q_value<=0.1 & all_sel$gene!="Y_RNA"))+
  #scale_fill_gradient(low = "red", high = "green") +
  scale_fill_gradient2(low = "red", mid = "white",high = "blue", midpoint = 0)+
  xlab("Genes ") +
  ylab("Lines") +
  theme_minimal()+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=8,face="bold"),
        axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=4, hjust = 1)) +
  labs(fill = "Expression level fold change")+
  #facet_wrap(~Freq,scales="free",ncol=4)+
  coord_flip()
 

write.table(all_sel,'/Users/instalar/Desktop/RNA_seq/all.sel.txt',row.names = FALSE)


#volcano plot (select line, or merged)
ggplot(RPMI_dif,aes(log2.fold_change.,-log(p_value)))+
  geom_point()+
  geom_vline(xintercept = c(-1.5,1.5),color="red")+
  geom_hline(yintercept = 2,color="blue")+
  theme_minimal()



###### importing cuffdiff results from Gekas analyses shBcat #####

# LINE RPMI shbcat
shbcat_rpmi<-read.delim('/Volumes/grcmc/YGUILLEN/Christos_RNASeq/Christos_shbcat/Christos_shbcat/RPMI8402_per-cond_classic/gene_exp.diff',sep="\t")
shbcat_rpmi$line<-c("RPMI")

#volcano plot (select line)
ggplot(shbcat_rpmi,aes(log2.fold_change.,-log(p_value)))+
  geom_point()+
  geom_vline(xintercept = c(-1.5,1.5),color="red")+
  geom_hline(yintercept = 2,color="blue")+
  theme_minimal()

### Overall expression levels ### exclude expression levels higher than 500
shbcat_rpmi_fpkm<-subset(shbcat_rpmi,shbcat_rpmi$value_1<=600 & shbcat_rpmi$value_2<=600)

ggplot(shbcat_rpmi_fpkm,aes(x=log(value_1+0.1),y=log(value_2+0.1)))+
  geom_point(size=1)+
  #facet_wrap(~line,scales="free")+
  geom_smooth(method='lm',color="blue")+
  geom_abline(intercept = 0, slope = 1,color="red")+
  theme_minimal()+
  xlab("Control")+
  ylab("shbcat")


## significant DEG in shBcat FDR 0.1
sig_exp_rpmi_bcat<-subset(shbcat_rpmi,shbcat_rpmi$p_value<=0.05)
genes_sh<-sig_exp_rpmi_bcat$gene

#heatmap
sig_exp_rpmi_bcat<-sig_exp_rpmi_bcat[order(sig_exp_rpmi_bcat$log2.fold_change.), ]
sig_exp_rpmi_bcat$gene <-reorder(sig_exp_rpmi_bcat$gene, sig_exp_rpmi_bcat$log2.fold_change.)
#substitute infinite values by maximum and minimum
sig_exp_rpmi_bcat[,10][(sig_exp_rpmi_bcat[,10]== -Inf)] <- min(sig_exp_rpmi_bcat[,10][sig_exp_rpmi_bcat[,10] != -Inf])
sig_exp_rpmi_bcat[,10][(sig_exp_rpmi_bcat[,10]== Inf)] <- max(sig_exp_rpmi_bcat[,10][sig_exp_rpmi_bcat[,10] != Inf])


ggplot(sig_exp_rpmi_bcat, aes(gene, line )) +
  geom_tile(aes(fill = log2.fold_change.), color = "white") +
  #scale_fill_gradient(low = "red", high = "green") +
  scale_fill_gradient2(low = "red", mid = "white",high = "blue", midpoint = 0)+
  xlab("Genes ") +
  ylab("Lines") +
  theme_minimal()+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=8,face="bold"),
        axis.text.x = element_text(size=8,angle = 45, hjust = 1),
        axis.text.y = element_text(size=5, hjust = 1)) +
  labs(fill = "Expression level fold change")+
  coord_flip()

write.table(sig_exp_rpmi_bcat,'/Volumes/grcmc/YGUILLEN/bcatenin_ChIPSeq/RNA_Seq_Gekas_R/sh_DEG_fdr01.csv',row.names = FALSE,quote = FALSE,sep="\t")

# rbind sh and PKF datasets
RPMI_dif$source<-"PKF"
shbcat_rpmi$source<-"shbcat"
sh_pkf_rpmi<-rbind(shbcat_rpmi,RPMI_dif)

#select significant shbcat DEG in bind datasets
shcat_bind<-subset(sh_pkf_rpmi, gene %in% genes_sh)
shcat_bind$dataset<-c("shbcat-RNAseq")

#heatmap

shcat_bind[,10][(shcat_bind[,10]== -Inf)] <- min(shcat_bind[,10][shcat_bind[,10] != -Inf])
shcat_bind[,10][(shcat_bind[,10]== Inf)] <- max(shcat_bind[,10][shcat_bind[,10] != Inf])

shcat_bind<-shcat_bind[order(shcat_bind$log2.fold_change.), ]
shcat_bind$gene <-reorder(shcat_bind$gene, shcat_bind$log2.fold_change.)
ggplot(shcat_bind, aes(source, gene )) +
  geom_tile(aes(fill = log2.fold_change.), color = "white") +
  geom_point(aes(size=(-log(q_value))),data = subset(shcat_bind,shcat_bind$p_value<=0.05))+
  scale_fill_gradient2(low = "red", mid = "white",high = "blue", midpoint = 0)+
  xlab("Genes ") +
  ylab("Experiment") +
  theme_minimal()+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=6,hjust = 1),
        axis.text.x = element_text(size=12,angle = 45, hjust = 1)) +
  labs(fill = "Expression level fold change")



#subset significant FDR <= 0.1 in RPMI PKF
sig_exp_rpmi_pkf<-subset(RPMI_dif,RPMI_dif$q_value<=0.1)
sig_exp_rpmi_pkf<-sig_exp_rpmi_pkf[order(sig_exp_rpmi_pkf$p_value), ]
genes_pkf<-sig_exp_rpmi_pkf$gene
#select significant in bind datasets
pkf_bind<-subset(sh_pkf_rpmi, gene %in% genes_pkf)
pkf_bind$dataset<-c("PKF-RNAseq")

sh_rpmi_merge<-rbind(shcat_bind,pkf_bind)


sh_rpmi_merge[,10][(sh_rpmi_merge[,10]== -Inf)] <- min(sh_rpmi_merge[,10][sh_rpmi_merge[,10] != -Inf])
sh_rpmi_merge[,10][(sh_rpmi_merge[,10]== Inf)] <- max(sh_rpmi_merge[,10][sh_rpmi_merge[,10] != Inf])

ggplot(sh_rpmi_merge, aes(sample_2, gene )) +
  geom_tile(aes(fill = log2.fold_change.), color = "white") +
  #scale_fill_gradient(low = "red", high = "green") +
  scale_fill_gradient2(low = "red", mid = "white",high = "blue", midpoint = 0)+
  xlab("Genes ") +
  ylab("Lines") +
  theme_minimal()+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=6,hjust = 1),
        axis.text.x = element_text(size=6,angle = 45, hjust = 1)) +
  labs(fill = "Expression level fold change")+
  facet_wrap(~dataset,scales="free",ncol=2)



##### EnrichR ####
## GO terms and pathways
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018","GO_Biological_Process_2018","TF_Perturbations_Followed_by_Expression")


genes_up<-(subset(sig_exp_rpmi_bcat,sig_exp_rpmi_bcat$log2.fold_change.>0))$gene
genes_down<-(subset(sig_exp_rpmi_bcat,sig_exp_rpmi_bcat$log2.fold_change.<0))$gene

genes_up<-tolower(genes_up)
genes_down<-tolower(genes_down)

enrichUP <- enrichr(genes_up, dbs)
enrichDOWN <- enrichr(genes_down, dbs)

UP <- enrichUP[["TF_Perturbations_Followed_by_Expression"]]
UP$state<-c("UP")

DOWN <- enrichDOWN[["TF_Perturbations_Followed_by_Expression"]]
DOWN$state<-c("DOWN")


allGO<-rbind(UP,DOWN)


bpsub<-subset(allGO,allGO$P.value<=0.05)
bpsub<- bpsub[order(bpsub$P.value),]

bpsub$Term <- factor(bpsub$Term, levels=(bpsub$Term)[rev(order(bpsub$P.value))])
bpsub<-bpsub[order(bpsub$P.value), ]

bpsub$db<-c("TF")

write.table(bpsub,"/Volumes/grcmc/YGUILLEN/Christos_RNASeq/Leonie_GO/GO_shbcat_TF_pval005.csv",quote = FALSE,row.names = FALSE,sep="\t")

ggplot(bpsub,aes(x=state,y=Term,fill=P.value),alpha = 0.5)+
  geom_tile(color="black",width=1,height=0.8)+
  geom_text(aes(label=tolower(Genes)), size=3,hjust=1,color="white")+
  scale_fill_gradient(low="darkblue", high="cadetblue1",limits=range(bpsub$P.value),na.value="white")+
  #scale_fill_viridis()+
  #facet_wrap(~db,scales="free")+
  theme_minimal()+
  theme(axis.text.y =element_text(size=9),axis.title=element_text(size=14,face="bold"))



### Import ChIPSeq results RPMI line Bcatenin ### own peaks cluster!

chip_bcat_rpmi_sarah<-read_xlsx('/Volumes/grcmc/YGUILLEN/bcatenin_ChIPSeq/ChIPSeq_Sarah/merged_DB_annotation_fdr005_final.xlsx')

#chip_bcat_rpmi<-read_xlsx('/Users/instalar/Desktop/bcatenin_ChIPSeq/ChIPSeq_Sarah/merged_DB1DB2_annotation_fdr001_final.xlsx')

#chip_bcat_rpmi<-read_xlsx('/Users/instalar/Desktop/bcatenin_ChIPSeq/ChIPSeq_Sarah/Bcatenin_AbRD_RPMI_DA3vsDA1_annotation_final.xlsx')

#chip_bcat_rpmi<-read_xlsx('/Users/instalar/Desktop/bcatenin_ChIPSeq/ChIPSeq_Sarah/Bcatenin_AbSC_RPMI_DA5vsDA1_annotation_final.xlsx')


#list of genes next to significat peaks in chipseq
genecat<-chip_bcat_rpmi_sarah$external_gene_name
genes_sh

#Common genes peakds Bcatenin ChIPseq and RNAseq shbcat
peak_rpmi_rna<-subset(RPMI_dif, gene %in% genecat)
peak_rpmi_rna_sig<-subset(peak_rpmi_rna,peak_rpmi_rna$p_value<=0.05)




#Common genes peaks Bcatenin ChIPseq and RNAseq RPMI sh
peak_rpmi_sh<-subset(shbcat_rpmi, gene %in% genecat)
peak_rpmi_sh_sig<-subset(peak_rpmi_sh,peak_rpmi_sh$p_value<=0.05)


peak_rpmi_sh_sig$gene<- factor(peak_rpmi_sh_sig$gene, levels=(peak_rpmi_sh_sig$gene)[rev(order(peak_rpmi_sh_sig$log2.fold_change.))])
peak_rpmi_sh_sig<-peak_rpmi_sh_sig[order(peak_rpmi_sh_sig$log2.fold_change.), ]

ggplot(peak_rpmi_sh_sig, aes(source, gene )) +
  geom_tile(aes(fill = log2.fold_change.), color = "white") +
  geom_point(aes(size=(q_value)))+
  #scale_fill_gradient(low = "red", high = "green") +
  scale_fill_gradient2(low = "red", mid = "white",high = "blue", midpoint = 0)+
  xlab("Genes ") +
  ylab("Lines") +
  theme_minimal()+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y = element_text(size=10,hjust = 1),
        axis.text.x = element_text(size=6,angle = 45, hjust = 1)) +
  labs(fill = "Expression level fold change")





###


######## RUN cummeRbund ########

library(cummeRbund)
cuff_R<-readCufflinks('/Users/instalar/Desktop/RNA_seq/Christos_PKF/cuffdiff/R_Samples/')
cuff_R

#Dispersion plot
disp<-dispersionPlot(genes(cuff_R))
disp

#Density plot
dens<-csDensity(genes(cuff_R))
dens

#Scatterplot
s<-csScatterMatrix(genes(cuff_R))
s

#boxplot
b<-csBoxplot(genes(cuff_R))
b

#volcano plot
v<-csVolcano(genes(cuff_R),"R_Samples_1","R_Samples_2")
v

#significant genes
mySigMat<-sigMatrix(cuff_R,level='genes',alpha=0.05)
mySigMat

mySigGeneIds<-getSig(cuff_R,alpha=0.05,level='genes')
mySigGenes<-getGenes(cuff_R,mySigGeneIds)

#Find similar genes (up to 20 most similar) (it makes sense with more than one line)
mySimilar<-findSimilar(cuff_R,"CTNNB1",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)



