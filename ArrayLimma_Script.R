#### Yolanda Guillén ###
## May 2020 Barcelona ##


# MANUALS
# https://www.bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html
# https://wikibits.ugent.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor


#BiocManager::install("oligo")
#BiocManager::install("mogene10sttranscriptcluster.db")
#BiocManager::install("mogene20sttranscriptcluster.db")

library("oligo")
library("tidyverse")
library("limma")
library("affy")
library("ggrepel")
library("mogene10sttranscriptcluster.db")
library("mogene20sttranscriptcluster.db")
library(pheatmap)
library(RColorBrewer)
library(reshape2)

head(ls("package:mogene10sttranscriptcluster.db"))

setwd("/Users/yguillen/Desktop/temp/beta_catenin_project/microarrays/FL")


# Read old annotations Christos
anotchrist<-read.delim('data_ann_Gekas.txt',header = TRUE)
anotchrist$PROBEID<-anotchrist$Row.names
anotchrist$Row.names<-NULL

### Add phenotype
targets <- readTargets()
row.names(targets)
targets$cond <- paste(targets$Genotype, targets$Infection, sep=".")
targets$cond <- as.factor(targets$cond)
levels(targets$cond)

targets<-AnnotatedDataFrame(targets)
row.names(targets)
pData(eset)

## input matrix expression
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles,phenoData = targets)

eset <- oligo::rma(affyRaw,target="core")
pData(eset)



## QUality control

# All samples
PCA_norm <- prcomp(t(exprs(eset)), scale = FALSE)

summary(PCA_norm)

percentVar <- round(100*PCA_norm$sdev^2/sum(PCA_norm$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_norm$x[,1], PC2 = PCA_norm$x[,2],
                     Genotype = pData(eset)$Genotype,
                     Infection = pData(eset)$Infection,
                     cond = pData(eset)$cond,
                     label = gsub('_\\(.*','',colnames(eset)))

ggplot(dataGG, aes(PC1, PC2,label)) +
  geom_point(aes(colour = cond),size=5) +
  geom_label_repel(aes(label=label,fill=cond),alpha=0.2,size=3)+
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  theme_bw()


oligo::boxplot(eset, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")


### HEATMAP
annotation_for_heatmap <- 
  data.frame(cond = eset$cond,
             sampID = row.names(pData(eset)))
row.names(annotation_for_heatmap)<-row.names(pData(eset))
annotation_for_heatmap$sampID<-("sampID")

dists <- as.matrix(dist(t(exprs(eset)), method = "manhattan"))

rownames(dists) <- row.names(pData(eset))

diag(dists) <- NA

myColors = list(
  sampID = c(sampID=brewer.pal(2,"BrBG")[1]),
  cond = c("KO.C"=brewer.pal(4,"Paired")[1],
            "KO.N"=brewer.pal(4,"Paired")[2],
            "WT.C"=brewer.pal(4,"Paired")[3],
            "WT.N"=brewer.pal(4,"Paired")[4]))


pheatmap(dists,
         annotation_row = annotation_for_heatmap,
         fontsize = 9,
         fontsize_col = 9,
         legend = FALSE,
         legend_breaks = c(min(dists, na.rm = TRUE), 
                           max(dists, na.rm = TRUE)),
         annotation_colors = myColors,
         main = "Clustering heatmap for the calibrated samples")





# Excluding probes with low expression
eset_medians <- rowMedians(Biobase::exprs(eset))

hist_res <- hist(eset_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
# 4 for BM, 3 for FL
man_threshold=3

abline(v = man_threshold, col = "coral4", lwd = 2)


no_of_samples <- table(pData(eset)$cond)
no_of_samples 

# DISCARD FILTERS TO RESEMBLE CHRISTOS DATA (eset_manfiltered)

#####
#minimum in one sample
#samples_cutoff <- 1

#idx_man_threshold <- apply(Biobase::exprs(eset), 1,
#                           function(x){
#                             sum(x > man_threshold) >= samples_cutoff})
#table(idx_man_threshold)

#eset_manfiltered <- subset(eset, idx_man_threshold)

#####


write.exprs(eset,file="data.txt")

## ANNOTATION BM is 10, FL is 20

Annot_bioc <- AnnotationDbi::select(mogene20sttranscriptcluster.db,
                                       keys = (featureNames(eset)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

Annot_bioc <- subset(Annot_bioc, !is.na(SYMBOL))

# highlighting multiple mappings, one probe multiple genes

anno_grouped <- group_by(Annot_bioc, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

head(anno_summarized)

anno_filtered <- filter(anno_summarized, no_of_matches >= 1)

head(anno_filtered)
probe_stats <- anno_filtered 

nrow(probe_stats)
Annot_bioc<-merge(Annot_bioc,probe_stats,by="PROBEID",all.x=TRUE)

### DON'T EXCLUDE MULTIPLE MAPPINGS, JUST HIGHLIGHT THEM IN ANNOTATION FILE
#ids_to_exlude <- (featureNames(eset) %in% probe_stats$PROBEID)
#table(ids_to_exlude)
#eset_final <- subset(eset, !ids_to_exlude)

validObject(eset_final)
head (eset_final)
head(Annot_bioc)

eset_final<-eset

fData(eset_final)$PROBEID <- rownames(fData(eset_final))

validObject(eset_final)


### LINEAR MODELS FOR MICROARRAYS

geno<-pData(eset_final)$Genotype
infection<-pData(eset_final)$Infection


##### Comparison bcat WT vs bcat KO in NICD+ context
i_NICD <- gsub('_\\(.*','',row.names(pData(eset_final)))
i_NICD<-i_NICD[infection == "N"]

design_NICD <- model.matrix(~ 0 + geno[infection == "N"])
colnames(design_NICD)[1:2] <- c("KO_NICD", "WT_NICD")
rownames(design_NICD) <- i_NICD 

##### Comparison bcat WT vs bcat KO in NICD WT context
w_NICD <- gsub('_\\(.*','',row.names(pData(eset_final)))
w_NICD<-w_NICD[infection == "C"]

design_wNICD <- model.matrix(~ 0 + geno[infection == "C"])
colnames(design_wNICD)[1:2] <- c("KO_bcat", "WT_bcat")
rownames(design_wNICD) <- w_NICD 

##### Comparison NICD+ vs Control in bcat WT context
bw_NICD <- gsub('_\\(.*','',row.names(pData(eset_final)))
bw_NICD<-bw_NICD[geno == "WT"]

design_bwNICD <- model.matrix(~ 0 + infection[geno == "WT"])
colnames(design_bwNICD)[1:2] <- c("C_bcatWT", "NICD_bcatWT")
rownames(design_bwNICD) <- bw_NICD 


##### Comparison NICD+ vs Control in bcat KO context
bko_NICD <- gsub('_\\(.*','',row.names(pData(eset_final)))
bko_NICD<-bko_NICD[geno == "KO"]

design_bkoNICD <- model.matrix(~ 0 + infection[geno == "KO"])
colnames(design_bkoNICD)[1:2] <- c("C_bcatKO", "NICD_bcatKO")
rownames(design_bkoNICD) <- bko_NICD 


# Example checking expression of one gene
#i_NICD<-i_NICD[infection == "N"]
crat_expr <- exprs(eset_final)[row.names(exprs(eset_final)) %in% c(Annot_bioc[Annot_bioc$SYMBOL=="Ctnnb1",]$PROBEID), ]
#crat_expr <- exprs(eset_final)[row.names(exprs(eset_final)) == "17523257", ]

crat_data <- as.data.frame(crat_expr)
#crat_data<-melt(crat_data)
crat_data$samp<-row.names(crat_data)

ggplot(data = crat_data, aes(x = samp, y = crat_expr)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))+
  ggtitle("Expression changes for gene")


### MAKE CONTRASTS
# NICD context, bcat KO vs bcat WT
cont.matrix.NICD <- makeContrasts(WT_NICD-KO_NICD,
                             levels=design_NICD)

fit_NICD <- eBayes(contrasts.fit(lmFit(eset_final[,infection == "N"],
                                              design = design_NICD),
                                        cont.matrix.NICD))

FCtab.NICD<-topTable(fit_NICD,number=Inf)

FCtab.NICD<-merge(FCtab.NICD,Annot_bioc,by="PROBEID",all.x=TRUE)
dim(FCtab.NICD[!is.na(FCtab.NICD$SYMBOL),])

write.table(FCtab.NICD,file="FCtab_NICD.tab",sep="\t")


# NICD WT context, bcat KO vs bcat WT
cont.matrix.wNICD <- makeContrasts(WT_bcat-KO_bcat,
                                  levels=design_wNICD)

fit_wNICD <- eBayes(contrasts.fit(lmFit(eset_final[,infection == "C"],
                                       design = design_wNICD),
                                 cont.matrix.wNICD))

FCtab.wNICD<-topTable(fit_wNICD,number=Inf)

FCtab.wNICD<-merge(FCtab.wNICD,Annot_bioc,by="PROBEID",all.x=TRUE)
dim(FCtab.wNICD[!is.na(FCtab.wNICD$SYMBOL),])

write.table(FCtab.wNICD,file="FCtab_wNICD.tab",sep="\t")

# bcat WT context, NICD CT vs NICD+
cont.matrix.bwNICD <- makeContrasts(C_bcatWT-NICD_bcatWT,
                                   levels=design_bwNICD)

fit_bwNICD <- eBayes(contrasts.fit(lmFit(eset_final[,geno == "WT"],
                                        design = design_bwNICD),
                                  cont.matrix.bwNICD))

FCtab.bwNICD<-topTable(fit_bwNICD,number=Inf)

FCtab.bwNICD<-merge(FCtab.bwNICD,Annot_bioc,by="PROBEID",all.x=TRUE)
dim(FCtab.bwNICD[!is.na(FCtab.bwNICD$SYMBOL),])

write.table(FCtab.bwNICD,file="FCtab_bwNICD.tab",sep="\t")

# bcat KO context, NICD CT vs NICD+
cont.matrix.bkoNICD <- makeContrasts(C_bcatKO-NICD_bcatKO,
                                    levels=design_bkoNICD)

fit_bkoNICD <- eBayes(contrasts.fit(lmFit(eset_final[,geno == "KO"],
                                         design = design_bkoNICD),
                                   cont.matrix.bkoNICD))

FCtab.bkoNICD<-topTable(fit_bkoNICD,number=Inf)

FCtab.bkoNICD<-merge(FCtab.bkoNICD,Annot_bioc,by="PROBEID",all.x=TRUE)
dim(FCtab.bkoNICD[!is.na(FCtab.bkoNICD$SYMBOL),])

write.table(FCtab.bkoNICD,file="FCtab_bkoNICD.tab",sep="\t")



# MERGE DATASETS
colnames(FCtab.NICD)[2]="logFC_NICD_pos"
colnames(FCtab.wNICD)[2]="logFC_NICD_C"
colnames(FCtab.bwNICD)[2]="logFC_bcat_WT"
colnames(FCtab.bkoNICD)[2]="logFC_bcat_KO"

FC.tab<-merge(FCtab.NICD[,c(1:7)],FCtab.wNICD[,c(1:7)],by="PROBEID")
FC.tab<-FC.tab[,c(1,2,5,6,8,11,12)]
FC.tab<-merge(FC.tab,FCtab.bwNICD[,c(1:7)],by="PROBEID")
FC.tab<-FC.tab[,c(1:8,11,12)]
FC.tab<-merge(FC.tab,FCtab.bkoNICD[,c(1:7)],by="PROBEID")
FC.tab<-FC.tab[,c(1:10,15,18,19)]
colnames(FC.tab)<-c("PROBEID","SYMBOL","GENENAME",
                    "logFC_NICD_pos","P.val.NICD_pos","adj.P.val.NICD_pos",
                    "logFC_NICD_C","P.val.NICD_C","adj.P.val.NICD_C",
                    "logFC_bcat_WT","P.val.bcat_WT","adj.P.val.bcat_WT",
                    "logFC_bcat_KO","P.val.bcat_KO","adj.P.val.bcat_KO")

FC.tab[FC.tab$P.val.NICD_pos<0.01,]$PROBEID
FC.tab[FC.tab$P.val.NICD_C<0.01,]$PROBEID


ggplot(data=FC.tab[FC.tab$P.val.bcat_WT<0.01 | FC.tab$P.val.bcat_KO<0.01,],aes(x=logFC_bcat_WT,y=logFC_bcat_KO))+
  geom_point(size=1,color="grey")+
  geom_point(data=FC.tab[FC.tab$P.val.NICD_pos<0.01,],color="red")+
  geom_density_2d(data=FC.tab[FC.tab$P.val.NICD_pos<0.01,],color="red")+
  #geom_point(data=FC.tab[FC.tab$P.val.NICD_C<0.01,],color="darkolivegreen")+
  geom_vline(xintercept=0)+
  geom_hline(yintercept = 0)+
  theme_bw()



######## Christos script ###

cont.matrix <- makeContrasts(NvsCinWT=WT.N-WT.C,
                             NvsCinKO=KO.N-KO.C,
                             #Diff=(KO.N-KO.C)-(WT.N-WT.C),
                             levels=design)


my_frame <- data.frame(exprs(eset))

design <- model.matrix(~0+TS)
View(design)

colnames(design) <- levels(TS)
fit <- lmFit(eset, design)
summary(fit)


fit2 <- contrasts.fit(fit, cont.matrix)
summary(fit2)


write.table(fit2,file="fit2",sep="\t")

results <- decideTests(fit2)
vennDiagram(results)
