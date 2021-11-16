library(reshape2)
library(ggnewscale)
library(ggplot2)
library(gridExtra)

ran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/kaiso_RD_match_plot_distance_dist.txt")
ran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/kaiso_RD_Li_match_plot_distance_dist.txt")
ran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/kaiso_SC_match_plot_distance_dist.txt")
ran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/kaiso_SC_Li_match_plot_distance_dist.txt")


ran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/RBLS_RBCS_plot_distance_dist.txt")


ran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HCT116/random_peaks/IKB_pIKB_H4K12_plot_distance_dist.txt",header = TRUE)
ran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HCT116/random_peaks/IKB_pIKB_panac_plot_distance_dist.txt",header = TRUE)
ran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HCT116/random_peaks/IKB_H4K12_HCT116.bed_plot_distance_dist.txt",header = TRUE)

ran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HT29/random_peaks/IKB_pIKB_HT29_H4K12_plot_distance_dist.txt",header = TRUE)
ran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HT29/random_peaks/IKB_pIKB_HT29_panacet_plot_distance_dist.txt",header = TRUE)

ran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/Caco2/random_peaks/IKB_pIKB_Caco2_H4K12_plot_distance_dist.txt",header = TRUE)
ran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/Caco2/random_peaks/IKB_pIKB_Caco2_panacet_plot_distance_dist.txt",header = TRUE)


distplot<-ggplot(ran,aes(x=reldist,y=fraction*100,group=data))+
  geom_point(aes(color=data))+
  geom_line(aes(group=data,color=data))+
  scale_color_brewer(palette="Set2")+
  geom_ribbon(aes(group=data,fill=data,ymax=fraction,ymin=0),alpha=0.2)+
  scale_fill_brewer(palette="Set2")+
  theme_light()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position = "top")

distplot

distran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/kaiso_RD_match_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/kaiso_RD_Li_match_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/kaiso_SC_match_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/kaiso_SC_Li_match_dist_iter_peaks.txt")



distran<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/random_test/RBLS_RBCS_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HCT116/random_peaks/IKB_pIKB_H4K12_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HCT116/random_peaks/IKB_pIKB_panac_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HCT116/random_peaks/IKB_H4K12_HCT116.bed_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HCT116/random_peaks/IKB_panacet_HCT116_dist_iter_peaks.txt")


distran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HT29/random_peaks/IKB_pIKB_HT29_H4K12_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HT29/random_peaks/IKB_pIKB_HT29_panacet_dist_iter_peaks.txt")

distran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/Caco2/random_peaks/IKB_pIKB_Caco2_H4K12_dist_iter_peaks.txt")
distran<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/Caco2/random_peaks/IKB_pIKB_Caco2_panacet_dist_iter_peaks.txt")


boxplot(distran$intersect_100bp)

distran<-melt(distran)

meand<-mean(distran[distran$variable=="intersect_100bp",]$value)
estd<-sd(distran[distran$variable=="intersect_100bp",]$value)
n<-nrow(distran)/2


left<-meand + qt(1 - 0.01 / 2, n - 1) * estd / sqrt(n) * c(-1, 1)[1]
right<-meand + qt(1 - 0.01 / 2, n - 1) * estd / sqrt(n) * c(-1, 1)[2]


left
right

pvalplot<-ggplot(distran,aes(x=value))+
  geom_density(aes(color=variable),adjust=1,size=2)+
  scale_color_brewer(palette="Paired")+
  geom_vline(xintercept = right,color="red",size=2,linetype="dotted")+
  xlim(c(0,max(distran$value)+median(distran$value)))+
#  ylim(c(0,5))+
  theme_light()+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.y = element_blank(),
        legend.position = "top")


grid.arrange(distplot,pvalplot,ncol=2)


## For IKB score distribution
distscore<-read.delim("/Users/yguillen/Desktop/temp/IKB_Hist_2020/Colon_cell_lines_ChIP/HCT116/random_peaks/IKB_PIKB_common_sets_d1.bed",header = FALSE)
colnames(distscore)<-c("chr_IKB","start_IKB","end_IKB","name_IKB","score_IKB",
                       "chr_pIKB","start_pIKB","end_pIKB","name_pIKB","score_pIKB")


wilcox.test(distscore$score_IKB,distscore$score_pIKB)

distscore<-subset(distscore,select=c(chr_IKB,score_IKB,score_pIKB))
distscore<-melt(distscore)

ggplot(distscore,aes(x=variable,y=value))+
  geom_point(aes(color=variable))+
  geom_boxplot(aes(color=variable),alpha=0.2)+
  scale_color_brewer(palette="Set2")+
  theme_light()

