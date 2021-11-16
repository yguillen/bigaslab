library(corrplot)
library(Hmisc)

## T-ALL and b-cat ##

df_tall <-read.delim("/Volumes/grcmc/TERESA LJ/Experimentos/Proyecto Quiescencia PDT005/2. RNAseq/bcat_RD_target_genes_expression_TALL.txt")

#df_tall <-read.delim("bcat_RD_target_genes_expression_TALL.txt")


df_target <-subset(df_tall,df_tall$GSEA=="TARGET")
dim(df_target) #265 pacientes


#creo una nueva variable para identificar los pacientes
Patient <-as.factor(rep(1:265,each=1))

df_target <-cbind(Patient,df_target)
df_target$Patient<-interaction('A', df_target$Patient)
str(df_target)
colnames(df_target)
row.names(df_target)

#solo survival data
df_target_survival <-df_target[,c(1,7,8)]
df_target_survival$Vital_status[df_target_survival$Vital_status== "Alive"] <- "0"
df_target_survival$Vital_status[df_target_survival$Vital_status== "Dead"] <- "1"
df_target_survival$Vital_status <-as.integer(df_target_survival$Vital_status)
df_target_survival <-na.omit(df_target_survival)
dim(df_target_survival)
str(df_target_survival)

#solo genes data
dim(df_target)
df_target_genes <-df_target[,c(1,9:ncol(df_target))]
row.names(df_target_genes)<-df_target_genes$Patient
df_target_genes<-df_target_genes[row.names(df_target_genes) %in% df_target_survival$Patient,]
dim(df_target_genes) #262 pacientes y 602 genes 

df_target_genes$Patient <-NULL

#Corrplot, 
#col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F","cyan", "#007FFF", "blue", "#00007F"))

rcorTALL<-rcorr(as.matrix(df_target_genes), type = "spearman")
corrplot(rcorTALL$r, order="hclust",tl.col="black", tl.cex = 0.1, type = "lower",
         p.mat = rcorTALL$P, sig.level = 0.05, insig = "blank", method = "color")

rcorTALL_r <-as.data.frame(rcorTALL$r)
#View(rcorTALL_r)


#replico:
ordercorTALL <-cor(as.matrix(df_target_genes))
ordercorTALL[lower.tri(ordercorTALL,diag=TRUE)] <-NA  #selecciono solo un triángulo porque la matriz es igual
ordercorTALL <-as.data.frame(as.table(ordercorTALL))  #Tlo convierto en tabla
ordercorTALL <-na.omit(ordercorTALL)  #elimino los valores donde hay NA
ordercorTALL <-ordercorTALL[order(abs(ordercorTALL$Freq),decreasing = T),] #ordeno de mayor a menor correlación

dim(ordercorTALL)
seleccionTALL <-ordercorTALL[c(1:150),] #para 68 150
hist(ordercorTALL$Freq)
selectionTALL<-subset(ordercorTALL,abs(ordercorTALL$Freq)>0.8)

var1 <-as.vector(unique(seleccionTALL$Var1))
length(var1)
var2 <-as.vector(unique(seleccionTALL$Var2))
length(var2)

genes <-unique(c(var1,var2))
length(genes)
#genes <-genes[-c(14,24,53,46)]
#genes

#dfGenes <-as.data.frame(genes) 

#write.table(dfGenes, "/Volumes/grcmc/TERESA LJ/Experimentos/Proyecto Quiescencia PDT005/2. RNAseq/T-ALL and Bcat/genes.txt", sep="\t",col.names = NA,  row.names = T)

clusterTALL <-df_target_genes[,genes]
dim(clusterTALL)

rclusterTALL<-rcorr(as.matrix(clusterTALL), type = "pearson")
corrplot(rclusterTALL$r, order="hclust",tl.col="black", tl.cex = 0.5, type = "lower",
         p.mat = rclusterTALL$P, sig.level = 0.05, insig = "blank", method = "color")


dendpacTALL <- hclust(dist(clusterTALL, method = "euclidean"), method = "complete") #filas/pacientes
dendTALL <-as.dendrogram(dendpacTALL)

library(dendextend)
dendTALL <-color_branches(dendTALL, k=2)
pacTALL <-as.data.frame(cutree(dendTALL, k=2))
bypal <-c("#E31A1C", "gold2","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C" )
pal<-colorRampPalette(bypal)
library(heatmaply)
heatmaply(clusterTALL,
          plot_method = "plotly",
          showticklabels = c(TRUE, FALSE),
          fontsize_col = 4,
          col = gplots::bluered(50),
          scale="column",
          Rowv = dendTALL)

# CLUSTERS DE PACIENTES
pac1_TALL <-subset(pacTALL,pacTALL$`cutree(dendTALL, k = 2)`=="1")
dim(pac1_TALL) 

pac2_TALL<-subset(pacTALL,pacTALL$`cutree(dendTALL, k = 2)`=="2")
dim(pac2_TALL) 

pac3_TALL<-subset(pacTALL,pacTALL$`cutree(dendTALL, k = 3)`=="3")
dim(pac3_TALL) 

#pac4_TALL<-subset(pacTALL,pacTALL$`cutree(dendTALL, k = 4)`=="4")
#dim(pac4_TALL) 

#Para el cluster de pacientes 1
survpac1_TALL<- as.data.frame(df_target_survival[df_target_survival$Patient %in% rownames(pac1_TALL),])
dim(survpac1_TALL) 
survpac1_TALL$cluster <- "1"
colnames(survpac1_TALL)

#Para el cluster de pacientes 2
survpac2_TALL<- as.data.frame(df_target_survival[df_target_survival$Patient %in% rownames(pac2_TALL),])
dim(survpac2_TALL) 
survpac2_TALL$cluster <- "2"
colnames(survpac2_TALL)

#Para el cluster de pacientes 3
survpac3_TALL<- as.data.frame(df_target_survival[df_target_survival$Patient %in% rownames(pac3_TALL),])
dim(survpac3_TALL) 
survpac3_TALL$cluster <- "3"
colnames(survpac3_TALL)

#Para el cluster de pacientes 4
survpac4_TALL<- as.data.frame(df_target_survival[df_target_survival$Patient %in% rownames(pac4_TALL),])
dim(survpac4_TALL) 
survpac4_TALL$cluster <- "4"
colnames(survpac4_TALL)

# Junto ambos cluster en una nueva tabla, por ejemplo el cluster 1 de pacientes y el cluster 2 de pacientes
comp_TALL <-rbind(survpac1_TALL,survpac2_TALL)
colnames(comp_TALL)


library(survminer)
surv_C_fit_TALL <- survfit(Surv(Event_free_surv_days, Vital_status) ~ cluster, data=comp_TALL)
ggsurvplot(surv_C_fit_TALL, data = comp_TALL, pval = TRUE,
           font.x = c(20),
           font.y = c(20),
           palette = c("green", "indianred","lightseagreen","purple"))


surv_C_fit_TALL$time
