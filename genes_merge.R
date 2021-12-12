
Dx_genecounts<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/data/TALL_sc/TALL_Bigas/Dx/Dx_genecounts.txt')
row.names(Dx_genecounts)<-Dx_genecounts$gene
Rx_genecounts<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/data/TALL_sc/TALL_Bigas/Rx/Rx_genecounts.txt')
row.names(Rx_genecounts)<-Rx_genecounts$gene
gene_annot<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/data/TALL_sc/gene_names.txt',sep=" ",header = FALSE)
colnames(gene_annot)<-c("gene","symbol")
row.names(gene_annot)<-gene_annot$gene

index_Dx<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/data/TALL_sc/TALL_Bigas/Dx/indexes.tsv',header = FALSE)

index_Rx<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/data/TALL_sc/TALL_Bigas/Rx/indexes.tsv',header = FALSE)
index_Rx<-index_Rx[-1,]
index_Rx$V3<-NULL

index_Dx$sample<-"Dx"
index_Rx$sample<-"Rx"

index_Dx_Rx<-rbind(index_Dx,index_Rx)
index_Dx_Rx$V1<-as.character(index_Dx_Rx$V1)
index_Dx_Rx$cell<-gsub('^P4.*','CD7+',index_Dx_Rx$V2)
index_Dx_Rx$cell<-gsub('^P.*CD7','CD7',index_Dx_Rx$cell)
index_Dx_Rx$cell<-gsub('_.*','',index_Dx_Rx$cell)

index_Dx_Rx$time_cell<-paste(index_Dx_Rx$sample,index_Dx_Rx$cell,sep="_")

index_Dx_Rx$batch<-index_Dx_Rx$time_cell
index_Dx_Rx$batch<-gsub('Dx_CD7\\+','0',index_Dx_Rx$batch)
index_Dx_Rx$batch<-gsub('Dx_.*','1',index_Dx_Rx$batch)
index_Dx_Rx$batch<-gsub('Rx_.*','2',index_Dx_Rx$batch)

gene_annot$gene<-as.character(gene_annot$gene)

Dx_Rx_genecounts<-merge(Dx_genecounts,Rx_genecounts,by=0)
row.names(Dx_Rx_genecounts)<-Dx_Rx_genecounts$Row.names
Dx_Rx_genecounts$Row.names<-NULL
Dx_Rx_genecounts$gene.x<-NULL

Dx_Rx_genecounts<-merge(gene_annot,Dx_Rx_genecounts,by=0)
row.names(Dx_Rx_genecounts)<-Dx_Rx_genecounts$gene
Dx_Rx_genecounts$Row.names<-NULL
Dx_Rx_genecounts$gene.y<-NULL
Dx_Rx_genecounts$gene<-NULL

row.names(Dx_Rx_genecounts)<-gsub('\\.','',row.names(Dx_Rx_genecounts))
Dx_Rx_genecounts$gene<-gsub('\\.','',Dx_Rx_genecounts$gene)

gene_names<-subset(Dx_Rx_genecounts,select=c("gene","symbol"))
Dx_Rx_genecounts$gene<-NULL
Dx_Rx_genecounts$symbol<-NULL

write.table(Dx_Rx_genecounts,'/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/data/TALL_sc/Dx_Rx_symbols_genecounts.txt',quote = FALSE,sep="\t")
write.table(gene_names,'/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/data/TALL_sc/gene_names_R.txt',quote = FALSE,sep="\t",row.names=FALSE)
write.table(index_Dx_Rx,'/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/data/TALL_sc/index.txt',quote = FALSE,sep="\t",row.names=FALSE)


## IMPORT scanpy data and b-cat signature
comp_all<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/beta_catenin_project/gene_lists/Dec_2021/comp_all.txt')

exp_mat<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/TALL_bcat_project/cell_browser/exprMatrix.tsv')
exp_mat$symbol<-exp_mat$gene
exp_mat$symbol<-gsub('ENS.*\\|','',exp_mat$symbol)

meta<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/TALL_bcat_project/cell_browser/meta.tsv')
colnames(meta)[1]<-"V1"
meta<-merge(meta,index_Dx_Rx,by="V1")
row.names(meta)<-meta$V1
meta$Louvain.Cluster<-as.factor(meta$Louvain.Cluster)

umap<-read.delim('/Users/yolanda_guillen/Desktop/IMIM/bigas_collab_scDx_Rel/TALL_bcat_project/cell_browser/umap_coords.tsv')
row.names(umap)<-umap$X
umap$X<-NULL
colnames(umap)<-c("UMAP1","UMAP2")

comp_sc<-exp_mat[exp_mat$symbol %in% c(as.character(comp_all$external_gene_name),"CTNNB1"),]
row.names(comp_sc)<-comp_sc$gene
comp_sc$gene<-NULL
comp_sc$symbol<-NULL

dim(comp_sc)
dim(meta)

library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(rev(brewer.pal(10,name="RdBu")))(50)
paletteLength<-length(cols)

comp_sct<-t(comp_sc)

myBreaks <- c(seq(min(t(scale(comp_sct))), -2.5, length.out=round(ceiling(paletteLength*0.10))),
              seq(-2.6,2.5,length.out = round(ceiling(paletteLength*0.80))),
              seq(2.6, max(t(scale(comp_sct))), length.out=round(floor(paletteLength*0.10)))) 

colnames(comp_sct)<-gsub('ENS.*\\|','',colnames(comp_sct))

phet<-pheatmap(as.data.frame(t(scale(comp_sct))),
        # annotation_col = as.data.frame(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,c(1,2,3,5,7,8,10,11,12,14,17,19)]),
         annotation_col = meta[,c(3,4,7,8,9)],
         show_colnames=F,
         show_rownames = T,
         fontsize = 4,
        # cutree_cols = 12,
        color=cols,
        breaks=myBreaks,
         cutree_rows = 10)

# Using only targets in more than 3 cells...
comp_bin<-comp_sct
comp_bin[comp_bin>0]<-1
# minimum 20 cells
sel<-colnames(comp_bin[,colSums(comp_bin)>3])

pheatmap(as.data.frame(t(scale(comp_sct[,colnames(comp_sct) %in% sel]))),
         # annotation_col = as.data.frame(bcatexp[bcatexp$GSEA=="TARGET" | bcatexp$GSEA=="line",][,c(1,2,3,5,7,8,10,11,12,14,17,19)]),
         annotation_col = meta[,c(7,8,9)],
         show_colnames=F,
         show_rownames = T,
         fontsize = 4,
         #        cutree_cols = 12,
         color=cols,
         breaks=myBreaks)


clusters<-data.frame(cluster=sort(cutree(phet$tree_row, k=10)))
patients<-data.frame(cluster=sort(cutree(phet$tree_col, k=12)))

clusters$Gene<-row.names(clusters)

table(patients$cluster)
table(clusters$cluster)
clusters[clusters$cluster==4,]$Gene

# Order genes (clusters)
geneorder<-data.frame(Gene=rownames(as.data.frame(t(scale(comp_sct)))[phet$tree_row[["order"]],]))
geneorder$order<-row.names(geneorder)
row.names(geneorder)<-geneorder$Gene

clusters<-merge(clusters,geneorder,by="Gene")
clusters$order<-as.numeric(clusters$order)

# Order patients (patients)
patientorder<-data.frame(Sample=colnames(as.data.frame(t(scale(comp_sct)))[,phet$tree_col[["order"]]]))
patientorder$order<-row.names(patientorder)
row.names(patientorder)<-patientorder$Sample
patientorder$Sample<-NULL

patients<-merge(patients,patientorder,by="row.names")
row.names(patients)<-patients$Row.names
patients$Row.names<-NULL
patients$order<-as.numeric(patients$order)


patinfo<-merge(comp_sct,patients,by="row.names")

library(ggplot2)
ggplot(patinfo,aes(x=as.factor(cluster),y=CTNNB1))+
  geom_boxplot(aes(color=as.factor(cluster)))+
  theme_bw()

patinfo<-merge(index_Dx_Rx,patinfo,by.y="Row.names",by.x="V1")
patinfo$order

ggplot(patinfo,aes(x=order,y=1,color=sample))+
  geom_point(aes(color=time_cell),size=1,alpha=0.3)+
  theme_bw()

#### UMAP COORDS AND expression matrix
exp_meta<-merge(comp_sct,umap,by=0)
row.names(exp_meta)<-exp_meta$Row.names
exp_meta$Row.names<-NULL

exp_meta<-merge(meta,exp_meta,by=0)

library(ggnewscale)

ggplot(exp_meta,aes(x=UMAP1,y=UMAP2))+
  geom_point(color="black",size=4)+
  scale_color_brewer(palette="Spectral")+
  new_scale("color")+
  geom_point(aes(color=as.factor(Louvain.Cluster)),size=3)+
  scale_color_brewer(palette="Spectral")+
  theme_bw()+
  theme(legend.position = "bottom")


p1<-ggplot(exp_meta,aes(x=UMAP1,y=UMAP2))+
  geom_point(color="grey",size=4)+
  new_scale("color")+
  geom_point(data=exp_meta[exp_meta$time_cell=="Rx_CD7+",],color="black",size=4)+
  geom_point(data=exp_meta[exp_meta$time_cell=="Rx_CD7+",],aes(color=TMED1),size=3)+
  scale_color_gradient(low="white",high="red")+
  theme_bw()+
  theme(legend.position = "bottom")

p2<-ggplot(exp_meta,aes(x=UMAP1,y=UMAP2))+
  geom_point(color="grey",size=4)+
  new_scale("color")+
  geom_point(data=exp_meta[exp_meta$time_cell=="Dx_CD7+",],color="black",size=4)+
  geom_point(data=exp_meta[exp_meta$time_cell=="Dx_CD7+",],aes(color=TMED1),size=3)+
  scale_color_gradient(low="white",high="red")+
  theme_bw()+
  theme(legend.position = "bottom")

p3<-ggplot(exp_meta,aes(x=UMAP1,y=UMAP2))+
  geom_point(color="grey",size=4)+
  new_scale("color")+
  geom_point(data=exp_meta[exp_meta$time_cell=="Dx_CD7-.CD34+",],color="black",size=4)+
  geom_point(data=exp_meta[exp_meta$time_cell=="Dx_CD7-.CD34+",],aes(color=TMED1),size=3)+
  scale_color_gradient(low="white",high="red")+
  theme_bw()+
  theme(legend.position = "bottom")

library(gridExtra)
grid.arrange(p1,p2,p3,ncol=2)

ggplot(exp_meta,aes(x=time_cell,y=CTNNB1))+
  geom_jitter(aes(color=time_cell))+
  theme_bw()


p4<-ggplot(exp_meta,aes(x=UMAP1,y=UMAP2))+
  geom_point(color="grey",size=4)+
  new_scale("color")+
  geom_point(data=exp_meta[exp_meta$time_cell=="Dx_CD7-.CD34+" & exp_meta$CTNNB1>=2,],color="black",size=4)+
  geom_point(data=exp_meta[exp_meta$time_cell=="Dx_CD7-.CD34+" & exp_meta$CTNNB1>=2,],aes(color=CTNNB1),size=3)+
  scale_color_gradient(low="white",high="red")+
  theme_bw()+
  theme(legend.position = "bottom")

p5<-ggplot(exp_meta,aes(x=UMAP1,y=UMAP2))+
  geom_point(color="grey",size=4)+
  new_scale("color")+
  geom_point(data=exp_meta[exp_meta$time_cell=="Rx_CD7+" & exp_meta$CTNNB1>=2,],color="black",size=4)+
  geom_point(data=exp_meta[exp_meta$time_cell=="Rx_CD7+" & exp_meta$CTNNB1>=2,],aes(color=CTNNB1),size=3)+
  scale_color_gradient(low="white",high="red")+
  theme_bw()+
  theme(legend.position = "bottom")

p6<-ggplot(exp_meta,aes(x=UMAP1,y=UMAP2))+
  geom_point(color="grey",size=4)+
  new_scale("color")+
  geom_point(data=exp_meta[exp_meta$time_cell=="Dx_CD7+" & exp_meta$CTNNB1>=2,],color="black",size=4)+
  geom_point(data=exp_meta[exp_meta$time_cell=="Dx_CD7+" & exp_meta$CTNNB1>=2,],aes(color=CTNNB1),size=3)+
  scale_color_gradient(low="white",high="red")+
  theme_bw()+
  theme(legend.position = "bottom")

grid.arrange(p4,p5,p6,ncol=2)
