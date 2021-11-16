
library(reshape2)
library(ggnewscale)
library(ggplot2)
library(gridExtra)

fcbcathist<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/multibigwig/scores_RBCSpeaks_bw_compare.tab")
fcbcathist<-melt(fcbcathist,id.vars = c("chr","start","end"))
fcbcathist$group<-"hist_target_bcat_SC"

fcbcat<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/multibigwig/scores_bcat_RBCS_bw_compare.tab")
fcbcat<-melt(fcbcat,id.vars = c("chr","start","end"))
fcbcat$group<-"bcat_SC_target"


fcrandhist<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/multibigwig/scores_randomnoRBCS_peaks_bw_compare.tab")
fcrandhist<-melt(fcrandhist,id.vars = c("chr","start","end"))
fcrandhist$group<-"hist_random"

fcrand<-read.delim("/Users/yguillen/Desktop/temp/beta_catenin_project/bed_files_chips/multibigwig/scores_bcat_randomRBCS_bw_compare.tab")
fcrand<-melt(fcrand,id.vars = c("chr","start","end"))
fcrand$group<-"bcat_SC_random"

fcdata_bcat<-rbind(fcbcat,fcrand)

fcdata_hist<-rbind(fcbcathist,fcrandhist)

summary(log2(fcrand[fcrand$variable=="bcat_SC_Ctl",]$value))
summary(log2(fcrand[fcrand$variable=="bcat_SC_Li",]$value))

summary(log2(fcbcat[fcbcat$variable=="bcat_SC_Ctl",]$value))
summary(log2(fcbcat[fcbcat$variable=="bcat_SC_Li",]$value))

summary(log2(fcrandhist[fcrandhist$variable=="H3Ac_Ct",]$value))
summary(log2(fcrandhist[fcrandhist$variable=="H3Ac_Li",]$value))

summary(log2(fcbcathist[fcbcathist$variable=="H3Ac_Ct",]$value))
summary(log2(fcbcathist[fcbcathist$variable=="H3Ac_Li",]$value))


wilcox.test(fcbcat[fcbcat$variable=="bcat_SC_Ctl",]$value,fcbcat[fcbcat$variable=="bcat_SC_Li",]$value)
wilcox.test(fcrand[fcrand$variable=="bcat_SC_Ctl",]$value,fcrand[fcrand$variable=="bcat_SC_Li",]$value)

wilcox.test(fcbcat[fcbcat$variable=="bcat_SC_Ctl",]$value,fcrand[fcrand$variable=="bcat_SC_Ctl",]$value)
wilcox.test(fcbcat[fcbcat$variable=="bcat_SC_Li",]$value,fcrand[fcrand$variable=="bcat_SC_Li",]$value)


wilcox.test(fcbcathist[fcbcathist$variable=="H3Ac_Ct",]$value,fcbcathist[fcbcathist$variable=="H3Ac_Li",]$value)
wilcox.test(fcrandhist[fcrandhist$variable=="H3Ac_Ct",]$value,fcrandhist[fcrandhist$variable=="H3Ac_Li",]$value)

wilcox.test(fcbcathist[fcbcathist$variable=="H3Ac_Ct",]$value,fcrandhist[fcrandhist$variable=="H3Ac_Ct",]$value)
wilcox.test(fcbcathist[fcbcathist$variable=="H3Ac_Li",]$value,fcrandhist[fcrandhist$variable=="H3Ac_Li",]$value)


bc<-ggplot(fcdata_bcat,aes(x=variable,y=log2(value)))+
  geom_point(aes(color=variable))+
  geom_boxplot(aes(color=variable),alpha=0.2)+
  facet_wrap(~group)+
  scale_color_brewer(palette="Set2")+
  theme_light()+
  scale_y_continuous(limits=c(-6,6))
bc

hc<-ggplot(fcdata_hist,aes(x=variable,y=log2(value)))+
  geom_point(aes(color=variable))+
  geom_boxplot(aes(color=variable),alpha=0.2)+
  facet_wrap(~group)+
  scale_color_brewer(palette="Set2")+
  theme_light()+
  scale_y_continuous(limits=c(-6,8))

hc

grid.arrange(bc,hc,ncol=1)
