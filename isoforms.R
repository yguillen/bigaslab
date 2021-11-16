
###### ISOFORM SWITCHING SCRIPT #######

# For Gekas sh and inhibitors in cell lines
# setwd('/Volumes/cancer/Gekas_RNAseq/isoformSwitch/')

# For T-ALL public data
setwd('/Volumes/cancer/TCell/salmon/')

#source("https://bioconductor.org/biocLite.R")
#
#if (!requireNamespace("BiocManager", quietly = TRUE)){
#  install.packages("BiocManager")
#}
#BiocManager::install("IsoformSwitchAnalyzeR")

library(IsoformSwitchAnalyzeR)

#browseVignettes("IsoformSwitchAnalyzeR")

#Example importing data
#salmonQuant<-importIsoformExpression(parentDir = system.file("extdata/",package="IsoformSwitchAnalyzeR"))


# Cell lines data
#salmonQuant<-importIsoformExpression(parentDir = "/Volumes/cancer/Gekas_RNAseq/salmon_sh/")

# TALL data
salmonQuant<-importIsoformExpression(parentDir = "/Volumes/cancer/TCell/salmon/salmon_R/")


# list containing
head(salmonQuant$length)
head(salmonQuant$abundance,2)
head(salmonQuant$counts,2)

#metadata
metadata<-read_xlsx("/Volumes/grcmc/YGUILLEN/IsoTALL_YGuillen/Public_data/metadata_all.xlsx")


salmonmeta<-as.data.frame(colnames(salmonQuant$abundance)[-1])
colnames(salmonmeta)<-"Sample"
salmonmeta$Sample<-gsub('quant_','',salmonmeta$Sample)
salmonmeta$Sample<-gsub('EGAR.*_','',salmonmeta$Sample)

salmonmeta<-merge(salmonmeta,metadata,by="Sample",sort=FALSE)


## We need to manually construct a design matrix
myDesign <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = salmonmeta$Type)

myDesign

# Create switchAnalyzeRlist
aSwitchList <- importRdata( isoformCountMatrix = salmonQuant$counts, 
                            isoformRepExpression = salmonQuant$abundance, designMatrix = myDesign,
                            isoformExonAnnoation = "/Volumes/cancer/db_files/transcriptomes/Homo_sapiens/ENSEMBL/gencode.v29.annotation.gtf",
                            showProgress = FALSE)

aSwitchList 

#Test differential isoforms
aSwitchList  <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchList ,
  reduceToSwitchingGenes=FALSE)


p1<-switchPlotTranscript(aSwitchList, gene = 'LEF1')
p2<-switchPlotGeneExp(aSwitchList, gene = 'LEF1')
p3<-switchPlotIsoExp(aSwitchList, gene = 'LEF1')
p4<-switchPlotIsoUsage(aSwitchList, gene = 'LEF1')

grid.arrange(p1,p2,p3,p4)

#data("exampleSwitchList")
#exampleSwitchList

## Need hsapiens annotations, check versions
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)

BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
library(BSgenome.Hsapiens.NCBI.GRCh38)

### isoformSwitchAnalysisPart1 needs the genomic sequence to predict ORFs
runSwitchList <- isoformSwitchAnalysisPart1( switchAnalyzeRlist = aSwitchList, genomeObject = Hsapiens,
                                                 # the reference to the human BS genome 
                                                 dIFcutoff = 0.3,
                                                 # Cutoff for finding switches - set high for short runtime in example data
                                                 pathToOutput = '/Volumes/cancer/TCell/isoformSwitchresults/',
                                                 outputSequences = FALSE # prevents outputting of the fasta files used for external sequence analysis
)




## Switching features summary
extractSwitchSummary(runSwitchList, 
                     dIFcutoff = 0.3 # supply the same cutoff to the summary function
)

#Swithced isoforms and conditions
extractTopSwitches(runSwitchList,filterForConsequences = TRUE)


## Part2

# Need run previously CPAT, PFAM and SIG

#Step by step
### Add CPAT analysis
exampleSwitchListAnalyzed <- analyzeCPAT(
  switchAnalyzeRlist   = runSwitchList,
  pathToCPATresultFile = "/Volumes/cancer/Gekas_RNAseq/isoformSwitch/cpat_result.txt",
  codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)


### Add PFAM analysis
exampleSwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = exampleSwitchListAnalyzed,
  pathToPFAMresultFile = "/Volumes/cancer/Gekas_RNAseq/isoformSwitch/pfam_results.txt",
  showProgress=FALSE
)

exampleSwitchListAnalyzed <- analyzeSignalP(
  switchAnalyzeRlist       = exampleSwitchListAnalyzed,
  pathToSignalPresultFile  = "/Volumes/cancer/Gekas_RNAseq/isoformSwitch/signalP_results.txt"
)

exampleSwitchListAnalyzed <- analyzeAlternativeSplicing(exampleSwitchListAnalyzed, quiet=TRUE)
table(exampleSwitchListAnalyzed$AlternativeSplicingAnalysis$IR)

switchPlot(exampleSwitchListAnalyzed, gene = 'ZAK')

#all steps together
progSwitchList <- isoformSwitchAnalysisPart2(switchAnalyzeRlist = runSwitchList,
                                                 dIFcutoff = 0.2, # Cutoff for finding switches - set high for short runtime in example data
                                                 #n = 15, # if plotting was enabled, it would only output the top 10 switches
                                                 removeNoncodinORFs = FALSE, # Because ORF was predicted de novo
                                                 pathToCPATresultFile = "/Volumes/cancer/Gekas_RNAseq/isoformSwitch/cpat_result.txt",
                                                 pathToPFAMresultFile = "/Volumes/cancer/Gekas_RNAseq/isoformSwitch/pfam_results.txt",
                                                 pathToSignalPresultFile = "/Volumes/cancer/Gekas_RNAseq/isoformSwitch/signalP_results.txt", 
                                                 codingCutoff= 0.725, # the coding potential cutoff we suggested for human
                                                 outputPlots = TRUE) # keeps the function from outputting the plots from this example )
                                


extractSwitchSummary( aSwitchList, filterForConsequences = FALSE, dIFcutoff = 0.3)
                      # supply the same cutoff to the summary function

# Top switches
extractTopSwitches(aSwitchList, filterForConsequences = TRUE, n=3)

# Visualize genes
switchPlot(runSwitchList, gene='CTNNB1',localTheme = theme_minimal())

# Global splicing
extractSplicingSummary(aSwitchList)
