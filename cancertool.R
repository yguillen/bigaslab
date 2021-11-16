### cancer tools ####

# Import patients months
timesurvival<-read.delim("/Volumes/grcmc/TERESA LJ/Experiments/Resistencias PDT005/2. RNAseq/Cancertool/TCGA_p53core/time_survival.txt")

normcounts_TCGA<-read.delim("/Volumes/grcmc/TERESA LJ/Experiments/Resistencias PDT005/2. RNAseq/Cancertool/TCGA_p53core/normcounts_TCGA.txt")
row.names(normcounts_TCGA)<-normcounts_TCGA$Sample
normcounts_TCGA$Sample=NULL

# matrix
normcounts_TCGA<-as.matrix(normcounts_TCGA)

heatmap.2(normcounts_TCGA,col=rev(morecols(50)),trace="none", main="TCGA p53 genes",scale="row")


library(heatmaply)
heatmaply(normcounts_TCGA,
          plot_method = "plotly",
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", 
          high = "red"), 
          k_row = 2,
          k_col = 2)
