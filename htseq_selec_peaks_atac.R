### From R sscript ATACseq.R extract those peaks that are linked to genes biased by gender, and those peaks that have a selected threshold value

WTbias<-as.character((WT_merge_peakanno[with(WT_merge_peakanno, (WT_merge_peakanno$SYMBOL) %in% genebias),])$Name)
KObias<-as.character((KO_merge_peakanno[with(KO_merge_peakanno, (KO_merge_peakanno$SYMBOL) %in% genebias),])$Name)
WTthres<-as.character(WT_merge_peakanno$Name)
KOthres<-as.character(KO_merge_peakanno$Name)

remnames<-c(WTbias,KObias)
mainnames<-c(WTthres,KOthres)
remnames<-gsub('/projects.*random_','random_',remnames)
mainnames<-gsub('/projects.*random_','random_',mainnames)


# Extract peaks from htseq
setwd("/Volumes/grcmc/YGUILLEN/ATAC-Seq/HTSEQ/peak_merged/")
samp2_1<-read.delim("htseq_mergepeaks_counts_sorted_random_2_1_29876_TAAGGC.bam.txt",header=FALSE)
samp2_2<-read.delim("htseq_mergepeaks_counts_sorted_random_2_2_29877_CGTACT.bam.txt",header=FALSE)
samp2_3<-read.delim("htseq_mergepeaks_counts_sorted_random_2_3_29878_AGGCAG.bam.txt",header=FALSE)
samp2_5<-read.delim("htseq_mergepeaks_counts_sorted_random_2_5_29880_GGACTC.bam.txt",header=FALSE)
samp2_6<-read.delim("htseq_mergepeaks_counts_sorted_random_2_6_29881_TAGGCA.bam.txt",header=FALSE)


samp2_1$V1<-gsub('_random.*','',samp2_1$V1)
samp2_1<-(samp2_1[with(samp2_1, !(samp2_1$V1) %in% remnames),])
samp2_1<-(samp2_1[with(samp2_1, (samp2_1$V1) %in% mainnames),])

samp2_2$V1<-gsub('_random.*','',samp2_2$V1)
samp2_2<-(samp2_2[with(samp2_2, !(samp2_2$V1) %in% remnames),])
samp2_2<-(samp2_2[with(samp2_2, (samp2_2$V1) %in% mainnames),])

samp2_3$V1<-gsub('_random.*','',samp2_3$V1)
samp2_3<-(samp2_3[with(samp2_3, !(samp2_3$V1) %in% remnames),])
samp2_3<-(samp2_3[with(samp2_3, (samp2_3$V1) %in% mainnames),])

samp2_5$V1<-gsub('_random.*','',samp2_5$V1)
samp2_5<-(samp2_5[with(samp2_5, !(samp2_5$V1) %in% remnames),])
samp2_5<-(samp2_5[with(samp2_5, (samp2_5$V1) %in% mainnames),])

samp2_6$V1<-gsub('_random.*','',samp2_6$V1)
samp2_6<-(samp2_6[with(samp2_6, !(samp2_6$V1) %in% remnames),])
samp2_6<-(samp2_6[with(samp2_6, (samp2_6$V1) %in% mainnames),])

# Write output
write.table(samp2_1,"htseq_mergepeaks_counts_sorted_random_2_1_29876_TAAGGC.bam.txt2",quote = FALSE,row.names = F)
write.table(samp2_2,"htseq_mergepeaks_counts_sorted_random_2_2_29877_CGTACT.bam.txt2",quote = FALSE, row.names = F)
write.table(samp2_3,"htseq_mergepeaks_counts_sorted_random_2_3_29878_AGGCAG.bam.txt2",quote = FALSE,row.names = F)
write.table(samp2_5,"htseq_mergepeaks_counts_sorted_random_2_5_29880_GGACTC.bam.txt2",quote = FALSE, row.names = F)
write.table(samp2_6,"htseq_mergepeaks_counts_sorted_random_2_6_29881_TAGGCA.bam.txt2",quote = FALSE, row.names = F)


