#install.packages('pcr')
library(pcr)

# Read two column file with Ct values of target gene and housekeeping gene
#fl <- system.file('extdata','ct1.csv',package = 'pcr')
# csv file data from Zara excel file google drive IMIM/Zara_Lef1/11Feb20...
ct1<-read.csv("/Users/yolanda_guillen/Desktop/IMIM/beta_catenin_project/Zara_RTPCR/qPCR_exp_LEF1_R.csv",dec = ",",stringsAsFactors=FALSE)
ct1$Condition<-factor(ct1$Condition,levels = c("Control","CHIR","LiCl","FH-535","IGC","XAV"))

#PLOT VARIATIONS
ct1metl<-melt(ct1,id.vars=c("Cell_line","Condition","Primers"))
ggplot(ct1metl,aes(x=Condition,y=value))+
  geom_point(aes(color=variable))+
  scale_color_brewer(palette="Paired")+
  facet_wrap(Cell_line~Primers)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
        axis.text.y = element_text(size=12))

ct1_Jurkat_p4_5<-ct1[!is.na(ct1$Ct_LEF1) & !is.na(ct1$Ct_TBP) & ct1$Cell_line=="Jurkat" & ct1$Primers=="P4-5",]
ct1_Jurkat_p4_5$Condition<-factor(ct1_Jurkat_p4_5$Condition,levels = c("Control","CHIR","LiCl","FH-535","IGC","XAV"))

ct1_Jurkat_p6<-ct1[!is.na(ct1$Ct_LEF1) & !is.na(ct1$Ct_TBP) & ct1$Cell_line=="Jurkat" & ct1$Primers=="P6",]
ct1_Jurkat_p6$Condition<-factor(ct1_Jurkat_p6$Condition,levels = c("Control","CHIR","LiCl","FH-535","IGC","XAV"))

ct1_Jurkat_p5_7<-ct1[!is.na(ct1$Ct_LEF1) & !is.na(ct1$Ct_TBP) & ct1$Cell_line=="Jurkat" & ct1$Primers=="P5_7",]
ct1_Jurkat_p5_7$Condition<-factor(ct1_Jurkat_p5_7$Condition,levels = c("Control","CHIR","LiCl","FH-535","IGC","XAV"))

ct1_RPMI_p4_5<-ct1[!is.na(ct1$Ct_LEF1) & !is.na(ct1$Ct_TBP) & ct1$Cell_line=="RPMI" & ct1$Primers=="P4-5",]
ct1_RPMI_p4_5$Condition<-factor(ct1_RPMI_p4_5$Condition,levels = c("Control","CHIR","LiCl","FH-535","IGC","XAV"))

ct1_RPMI_p6<-ct1[!is.na(ct1$Ct_LEF1) & !is.na(ct1$Ct_TBP) & ct1$Cell_line=="RPMI" & ct1$Primers=="P6",]
ct1_RPMI_p6$Condition<-factor(ct1_RPMI_p6$Condition,levels = c("Control","CHIR","LiCl","FH-535","IGC","XAV"))

ct1_RPMI_p5_7<-ct1[!is.na(ct1$Ct_LEF1) & !is.na(ct1$Ct_TBP) & ct1$Cell_line=="RPMI" & ct1$Primers=="P5_7",]
ct1_RPMI_p5_7$Condition<-factor(ct1_RPMI_p5_7$Condition,levels = c("Control","CHIR","LiCl","FH-535","IGC","XAV"))



res <- pcr_analyze(ct1_RPMI_p6[,c(4,5)],
                   group_var = ct1_RPMI_p6$Condition,
                   reference_gene = 'Ct_TBP',
                   reference_group = 'Control')

res


gg6 <- pcr_analyze(ct1_RPMI_p6[,c(4,5)],
                   group_var = ct1_RPMI_p6$Condition,
                   reference_gene = 'Ct_TBP',
                   reference_group = 'Control',
                   plot = TRUE) +
  labs(x = '', y = 'Relative mRNA expression') +
  ggtitle(label = 'RPMI_P6')

gg1
gg2
gg3

gg4
gg5
gg6

grid.arrange(gg1,gg4,
             gg2,gg5,
             gg3,gg6,ncol=2)


# FOLD CHANGE METHOD

# Jurkat
p4_5=ct1[ct1$Cell_line=="Jurkat" & ct1$Primers=="P4-5",]
#Remove one control and one CHIR from p6 to merge
p6=ct1[ct1$Cell_line=="Jurkat" & ct1$Primers=="P6",]
p6=p6[-c(1,6),]

dmet<-data.frame(p4_5$Condition,p4_5$Ct_LEF1,p6$Ct_LEF1)
dmet<-dmet[!is.na(dmet$p4_5.Ct_LEF1) & !is.na(dmet$p6.Ct_LEF1),]


# RPMI
p4_5=ct1[ct1$Cell_line=="RPMI" & ct1$Primers=="P4-5",]
#Remove one FH from p6 to merge
p6=ct1[ct1$Cell_line=="RPMI" & ct1$Primers=="P6",]
p6=p6[-c(10),]

dmet<-data.frame(p4_5$Condition,p4_5$Ct_LEF1,p6$Ct_LEF1)
dmet<-dmet[!is.na(dmet$p4_5.Ct_LEF1) & !is.na(dmet$p6.Ct_LEF1),]

# delta_ct method

# delta_ct method
## calculate caliberation
res <- pcr_analyze(dmet[,c(2,3)],
                   group_var = dmet$p4_5.Condition,
                   reference_group = 'Control',
                   method = 'delta_ct')

R1<-pcr_analyze(dmet[,c(2,3)],
            group_var = dmet$p4_5.Condition,
            reference_group = 'Control',
            method = 'delta_ct', 
            plot = TRUE) +
  scale_fill_brewer(palette = "Paired")+
  theme(legend.position = 'top',
        legend.direction = 'horizontal') +
  theme_bw()+
  labs(x = '', y = 'Relative fold change')+
  ggtitle("RPMI - Delta fold change method")

j1
R1

grid.arrange(j1,R1)
