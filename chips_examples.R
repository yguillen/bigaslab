### After BAM alignment (docker_chipseq_pipeline) example bam managing with Rsamtools

## Here's an example from day4 course
#First we create a BamFilelist object
dataDir<-"/Day4/data_for_practical/"
bamFiles<-dir(file.path(getwd(), dataDir), pattern="*.bam",full.name=T)
names(bamFiles)<-gsub(".bam","",basename(bamFiles))
bfList<-BamFileList(bamFiles)
bfList

path(bfList)
samHeader<-scanBamHeader(path(bfList["TF_1"]))
str(samHeader,max.level=2)

samHeader[[1]]$targets
samHeader[[1]]$text["@HD"]

## If we want extract only information from Chromosome 1:
chr1dat<-samHeader[[1]]$targets["chr1"]
chr1dat
chr1range<-GRanges(seqnames=names(chr1dat),ranges=IRanges(1,chr1dat))
param<-ScanBamParam(which=chr1range)

# Note that in this example, we use the default ScanBamParam ’flag’ argument. 
# However we could use this to build further filters (e.g. to remove duplicates or other potential artifacts) using the scanBamFlag() function. 
# Look at the options for this function for more details using ?scanBamFlag

# Calculate the average read length of the selected reads (in case some reads were trimmed).
alignDat <- readGAlignmentPairs(path(bfList["TF_1"]),param=param)
alignGR<-granges(alignDat)
seqlevels(alignGR)<-"chr1"
median(width(alignGR))


################## BAYESPEAK PEAK CALLING FROM BED ALIGNMENT FILES ###############

setwd("/Users/instalar/Desktop/bcatenin_ChIPSeq/raw_data_ChIPSeq")
raw.output <- bayespeak("output_unique_qual.bed","output_IM_qual.bed")

output<-summarize.peaks(raw.output,method="lowerbound")
write.table(as.data.frame(output),file="peaks.txt",quote=FALSE)

raw.output$peaks

# Example
data(raw.output)
raw.output$call

min.job <- min(raw.output$peaks$job)
max.job<-max(raw.output$peaks$job)

par(mfrow=c(2,2),ask=TRUE)
for(i in min.job:max.job) {plot.PP(raw.output,job=i,ylim=c(0,50))}


### MACS2 peakcalling ###

# After MACS2 peakcalling: docker run -v ${RAWDATA}:/data dceoy/macs2 callpeak -t /data/output_unique.bam -c /data/Input.bam -n /data/mypeaks

# Four files are created: mypeakds_peaks.narrowPeak, mypeaks_summits.bed, mypeaks_peaks.xls and mypeaks_model.r

# Loading peak data:
mypeaks<-read.delim("//DB_peaks_peaks.narrowPeak",header=FALSE)
colnames(mypeaks)<-c("chrom","chromStart","chromEnd","name","score","strand","fold.enrichment","log10.pvalue","log10.qvalue","peak")

#Add +1 to 0 start coordinates
GRanges(mypeaks$chrom,IRanges(mypeaks$chromStart+1,mypeaks$chromEnd))

## Example using ChIPQsample()
#To faster computation, select only peaks from chromosome 1 TF_1.bam
bamFile<-file.path("/Users/instalar/Desktop/bcatenin_ChIPSeq/peakcall_cluster/sorted_unique_output_DB1_21487_ATTCTC_trimmed.fq.gz.sam.bam")
#annotation refers to gene annotation of reference genome, hg19 for human; peaks is the otuput of peak caller (in this case MACS2)
exampleExp = ChIPQCsample(bamFile,peaks=mypeaks,blacklist = blacklist_HUM)

## QC metric summary
QCmetrics(exampleExp)

## plotCC function to calculate cross-coverage scores
plotCC(exampleExp)
#shift size of the maximum should corresponds to fragment size

fragmentlength(exampleExp)

#We can extract RelCC and FragCC scores. They relate the efficiency of ChIP (fragCC, or the maximum of cross coverage score) and efficiency of ChiP compared to artifact signal (RelCC)
#FragCC = CrossCoverageScoremax
#RelCC = CrossCoverageScoremax/CrossCoverageScorereadlength

FragmentLengthCrossCoverage(exampleExp)
ReadLengthCrossCoverage(exampleExp)

RelCC=FragmentLengthCrossCoverage(exampleExp)/ReadLengthCrossCoverage(exampleExp)
RelCC
#or
RelativeCrossCoverage(exampleExp)

## RCC is > 1, so the ChIP is successful

## Ditribution of signal within peaks
#Percentage of reads within peaks:
frip(exampleExp)*100
plotFrip(exampleExp)

## Distribution of signal in annotated regions
#Distribution of signal within genomic interval annotation over that expected given their size
regi(exampleExp)
plotRegi(exampleExp)


path <- system.file("extdata", package="ChIPpeakAnno")
files <- dir(path, "broadPeak")
data <- sapply(file.path(path, files), toGRanges, format="broadPeak")
names(data) <- gsub(".broadPeak", "", files)




