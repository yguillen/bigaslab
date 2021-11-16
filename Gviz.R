library(Gviz)
library(GenomicRanges)

data("cpgIslands")
class(cpgIslands)

chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands)
atrack <- AnnotationTrack(cpgIslands, name = "CpG")

plotTracks(atrack)
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))

itrack <- IdeogramTrack(genome = gen, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack))

data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = gen,chromosome = chr, name = "Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))

biocLite('BSgenome.Hsapiens.UCSC.hg19')
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack, grtrack,strack), from = 26591822, to = 26591852, cex = 0.8)
