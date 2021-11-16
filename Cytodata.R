
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
biocLite("openCyto")
biocLite("ggcyto")

library(flowCore)
library(openCyto)
library(ggcyto)
library(flowViz)

setwd("/Volumes/grcmc/YGUILLEN/FlowCore/Mice_CRuiz/ciclo/") # add the appropriate destination for your directory. 

#if more than one file
file_list<-(c("280219_84 ciclo.fcs","280219_85 ciclo.fcs")) 

#multiple fcs at a time
fs <- read.flowSet(files = file_list,
                   alter.names = TRUE,
                   column.pattern = ".A",
                   full=TRUE,
                   phenoData=list(name="SAMPLE ID", Filename="$FIL"))
fs
summary(fs)
pData(phenoData(fs))


#Example median of all markers
fsApply(fs, each_col, median)


#Compensation matrix
#two files
x <- fs[[1]]
y <- fs[[2]]

comp_list_x <- spillover(x) 
comp_list_y <- spillover(y) 

comp_x<-comp_list_x[[1]]
comp_y<-comp_list_y[[1]]


x_comp <- compensate(x, comp_x)
y_comp <- compensate(y, comp_y)

#Compensation effects
library(gridExtra)
#Example overlap between PE.A and PE.Cy7.A
transList_x <- estimateLogicle(x, c("PE.A","PE.Cy7.A"))
transList_y <- estimateLogicle(y, c("PE.A","PE.Cy7.A"))

p1 <- autoplot(transform(x, transList_x), "PE.A", "PE.Cy7.A") +
  ggtitle("Before")
p2 <- autoplot(transform(x_comp, transList_x), "PE.A", "PE.Cy7.A") +
  ggtitle("Compensation")

p3 <- autoplot(transform(y, transList_y), "PE.A","PE.Cy7.A") +
  ggtitle("Before")
p4 <- autoplot(transform(y_comp, transList_y), "PE.A","PE.Cy7.A") +
  ggtitle("Compensation")

grid.arrange(as.ggplot(p1), as.ggplot(p2),as.ggplot(p3), as.ggplot(p4), ncol = 2)

#2D plot
#Compensated first
autoplot(x_comp, "FSC.A", "SSC.A")
#First file
autoplot(fs[[1]], "FSC.A", "SSC.A")
#All files
autoplot(fs, "FSC.A", "SSC.A")
autoplot(fs, "FSC.A")

#Transform Data
P1<-autoplot(transform(fs[[1]],
                   FSC.A=log(FSC.A), SSC.A=log(SSC.A) ),
         "log.FSC.A","log.SSC.A")
P2<-autoplot(fs[[1]],"FSC.A","SSC.A")
grid.arrange(as.ggplot(P1),as.ggplot(P2))


## GATING ##
#BASIC
rectGate <- rectangleGate(filterId="Fluorescence Region", "FSC.A"=c(0,150000), "SSC.A"=c(0, 100000))
result = filter(fs[[1]],rectGate) 
result
summary(result)
summary(result)$n
summary(result)$true
summary(result)$p


summary(filter(fs,rectGate))

#multiple
summary(filter(fs[[1]],
               kmeansFilter("FSC.A"=c("Low", "Medium", "High"),
                            filterId="myKMeans")))
## SUBSETTING
morphGate <- norm2Filter("FSC.A", "SSC.A", filterId="MorphologyGate", scale=2)
smaller <- Subset(fs, morphGate) 
fs[[1]]
autoplot(smaller,"FSC.A","SSC.A")

#Create compensation matrix from controls

keyword(x,c("$P1E", "$P2E", "$P3E", "$P4E","$P5E", "$P6E", "$P7E", "$P8E"))

n <- as.data.frame(exprs(ntl))  # extract expression values and put into a data frame
colnames(n)
#with colours indicating density
colfunc <- colorRampPalette(c("white", "lightblue", "green", "yellow", "red"))
# this colour palette can be changed to your taste 

df_n<-n[,c(1,4,7,10,13,16,19,22,25,28,31)]

colnames(df_n)<-c("FSC_A","SSC_A","CD45","DAPI","CD150","Kit","CD48","Sca1","Lin","CD34","CD135")


ggplot(df_n, aes(x=FSC_A,y=SSC_A)) +
  #geom_point(size=2)+
  #ylim(0, 500000) +
  #xlim(0,5000000) +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  #geom_density2d(colour="black", bins=5) # draws the lines inside
  theme_bw()


ggplot(subset(df_n,df_n$SSC_A<75000 & df_n$FSC_A<125000), aes(x=FSC_A, y=SSC_A)) +
  geom_point(aes(size=DAPI))+
  #ylim(0, 500000) +
  #xlim(0,5000000) +
  #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  #scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  #geom_density2d(colour="black", bins=5) # draws the lines inside
  theme_bw()


ggplot(subset(df_n,df_n$SSC_A<75000 & df_n$FSC_A<100000), aes(x=FSC_A, y=DAPI)) +
  geom_point()+
  #ylim(0, 500000) +
  #xlim(0,5000000) +
  #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  #scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  #geom_density2d(colour="black", bins=5) # draws the lines inside
  theme_bw()
