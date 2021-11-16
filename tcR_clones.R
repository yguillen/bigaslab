#### Install packages (do only once)


BiocManager::install("tcR")

library(tcR)
data(twa)
head(twa[[1]])

library(stringdist)

# Gene alphabets - character vectors with names of genes for TCR and Ig.
?genesegments
data(genesegments)

# TCRalpha
TCRalpha_1_P<-parse.mixcr("/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRalpha_1_P.clones.txt")
TCRalpha_2_P<-parse.mixcr("/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRalpha_2_P.clones.txt")
TCRalpha_1_R<-parse.mixcr("/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRalpha_1_R.clones.txt")
TCRalpha_2_R<-parse.mixcr("/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRalpha_2_R.clones.txt")

#TCRbeta
TCRbeta_1_P<-parse.mixcr("/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRbeta_1_P.clones.txt")
TCRbeta_1_R<-parse.mixcr("/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRbeta_1_R.clones.txt")
TCRbeta_2_P<-parse.mixcr("/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRbeta_2_P.clones.txt")
TCRbeta_2_R<-parse.mixcr("/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRbeta_2_R.clones.txt")




# list of dataframes for TCR alpha
TCRalpha <- list(P1 = TCRalpha_1_P, 
                   R1 = TCRalpha_1_R,
                   P2 = TCRalpha_2_P,
                   R2 = TCRalpha_2_R)  ## create a named list of your data frames


### stats ###
TCRalpha_rep<-repseq.stats(TCRalpha)
write.table(TCRalpha_rep,"/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRalpha_rep.txt",sep="\t",quote = FALSE)

### Most abundant clonotypes (accounting for 25% of all reads)
clonal.proportion(TCRalpha, 25)

## TOP 10 clonotypes (10 most abundant clonotypes)
top.proportion(TCRalpha, 10)
# proportion of grouped top (top10, top11-100, etc)
vis.top.proportions(TCRalpha,c(10, 500, 3000, 10000))


## Search for target CDR3 sequence patient 1 -- ratón 1. Extracting TOP500 from PATIENT 1

tcra_1_p <- data.frame(CDR3.amino.acid.sequence = TCRalpha_1_P$CDR3.amino.acid.sequence,
                  V.genes = TCRalpha_1_P$V.gene, stringsAsFactors = F)

# All 
tcra_1_p <- data.frame(CDR3.amino.acid.sequence = TCRalpha_1_P$CDR3.amino.acid.sequence,
                       V.genes = TCRalpha_1_P$V.gene, stringsAsFactors = F)

tcra_1_p

write.table(tcra_1_p,"/Volumes/grcmc/YGUILLEN/TCR_mice/output/tcra_1_p.txt",sep="\t",quote = FALSE)

TCRalpha <- set.rank(TCRalpha)

# Search these top10 clonotypes from patient 1 in mouse 1
# exact match
cmv.imm.ex <- 
  find.clonotypes(.data = TCRalpha[c(1,2)], .targets = tcra_1_p[,1], .method = 'exact',
                  .col.name = c('Read.count', 'Umi.proportion'),
                  .verbose = F)

cmv.imm.ex<-cmv.imm.ex[!is.na(cmv.imm.ex$Read.count.R1) & !is.na(cmv.imm.ex$Read.count.P1),]

write.table(cmv.imm.ex,"/Volumes/grcmc/YGUILLEN/TCR_mice/output/cmv.imm.ex.txt",sep="\t",quote = FALSE)



#hammer distance
cmv.imm.ex <- 
  find.clonotypes(.data = TCRalpha[c(1,2)], tcra_1_p, 'hamm', 
                  .target.col = c('CDR3.amino.acid.sequence', 'V.gene'),
                  .col.name = c('Read.count', 'Umi.proportion'),
                  .verbose = T)

cmv.imm.ex[!is.na(cmv.imm.ex$Read.count.R1) & !is.na(cmv.imm.ex$Read.count.P1),]


## Search for target CDR3 sequence patient 2 -- ratón 2. Extracting TOP100 from PATIENT 2

tcra_2_p <- data.frame(CDR3.amino.acid.sequence = TCRalpha_2_P$CDR3.amino.acid.sequence,
                       V.genes = TCRalpha_2_P$V.gene, stringsAsFactors = F)

tcra_2_p

TCRalpha <- set.rank(TCRalpha)

# Search these top100 clonotypes from patient 2 in mouse 2
# exact match
cmv.imm.ex <- 
  find.clonotypes(.data = TCRalpha[c(3,4)], .targets = tcra_2_p[,1], .method = 'exact',
                  .col.name = c('Read.count', 'Umi.proportion'),
                  .verbose = F)

cmv.imm.ex[!is.na(cmv.imm.ex$Read.count.R2),]



#levenshtein distance
cmv.imm.ex <- 
  find.clonotypes(.data = TCRalpha[c(3,4)], tcra_2_p, 'hamm', 
                  .target.col = c('CDR3.amino.acid.sequence', 'V.gene'),
                  .col.name = c('Read.count', 'Umi.proportion'),
                  .verbose = T)

cmv.imm.ex<-cmv.imm.ex[!is.na(cmv.imm.ex$Read.count.R2),]

write.table(cmv.imm.ex,"/Volumes/grcmc/YGUILLEN/TCR_mice/output/cmv.imm.ex_pat2.txt",sep="\t",quote = FALSE)


# Similarities
repOverlap(TCRalpha, 'hamm', 'nuc', .norm = T, .verbose = T)
vis.heatmap(repOverlap(TCRalpha, 'hamm', 'aa', .vgene = T, .verbose = T), .title = 'TCRalpha - (ave)-intersection', .labs = '')



TCRalpha.top <- top.cross(.data = TCRalpha, .n = seq(500, 10000, 500), .verbose = F, .norm = T)
top.cross.plot(TCRalpha.top)




# list of dataframes for TCR beta
TCRbeta <- list(P1 = TCRbeta_1_P, 
                 R1 = TCRbeta_1_R,
                 P2 = TCRbeta_2_P,
                 R2 = TCRbeta_2_R)  ## create a named list of your data frames

### stats ###
TCRbeta_rep<-repseq.stats(TCRbeta)

write.table(TCRbeta_rep,"/Volumes/grcmc/YGUILLEN/TCR_mice/output/TCRbeta_rep.txt",sep="\t",quote = FALSE)


### Most abundant clonotypes (accounting for 25% of all reads)
clonal.proportion(TCRbeta, 25)

## TOP 10 clonotypes (10 most abundant clonotypes)
top.proportion(TCRbeta, 10)
# proportion of grouped top (top10, top11-100, etc)
vis.top.proportions(TCRbeta,c(10, 500, 3000, 10000))


## Search for target CDR3 sequence patient 1 -- ratón 1. Extracting TOP500 from PATIENT 1

tcrb_1_p <- data.frame(CDR3.amino.acid.sequence = TCRbeta_1_P$CDR3.amino.acid.sequence,
                       V.genes = TCRbeta_1_P$V.gene, stringsAsFactors = F)

tcrb_1_p

TCRbeta <- set.rank(TCRbeta)
# Search these top10 clonotypes from patient 1 in mouse 1
# exact match
cmv.imm.ex <- 
  find.clonotypes(.data = TCRbeta[c(1,2)], .targets = tcrb_1_p[,1], .method = 'exact',
                  .col.name = c('Read.count', 'Umi.proportion'),
                  .verbose = F)

cmv.imm.ex[!is.na(cmv.imm.ex$Read.count.R1),]

#levenshtein distance
cmv.imm.ex <- 
  find.clonotypes(.data = TCRbeta[c(1,2)], tcrb_1_p, 'hamm', 
                  .target.col = c('CDR3.amino.acid.sequence', 'V.gene'),
                  .col.name = c('Read.count', 'Umi.proportion'),
                  .verbose = T)

cmv.imm.ex[!is.na(cmv.imm.ex$Read.count.R1) & !is.na(cmv.imm.ex$Read.count.P1),]


## Search for target CDR3 sequence patient 2 -- ratón 2. Extracting TOP100 from PATIENT 2

tcrb_2_p <- data.frame(CDR3.amino.acid.sequence = TCRbeta_2_P$CDR3.amino.acid.sequence,
                       V.genes = TCRbeta_2_P$V.gene, stringsAsFactors = F)

tcrb_2_p

TCRbeta <- set.rank(TCRbeta)
# Search these top100 clonotypes from patient 2 in mouse 2
# exact match
cmv.imm.ex <- 
  find.clonotypes(.data = TCRbeta[c(3,4)], .targets = tcrb_2_p[,1], .method = 'exact',
                  .col.name = c('Read.count', 'Umi.proportion'),
                  .verbose = F)

cmv.imm.ex[!is.na(cmv.imm.ex$Read.count.R2),]

#hamming distance
cmv.imm.ex <- 
  find.clonotypes(.data = TCRbeta[c(3,4)], tcrb_2_p, 'hamm', 
                  .target.col = c('CDR3.amino.acid.sequence', 'V.gene'),
                  .col.name = c('Read.count', 'Umi.proportion'),
                  .verbose = T)

cmv.imm.ex[!is.na(cmv.imm.ex$Read.count.R2),]


## Similarities
repOverlap(TCRbeta[1:2], 'hamm', 'nuc', .verbose = T)

repOverlap(TCRbeta, 'hamm', 'nuc', .norm = T, .verbose = T)

vis.heatmap(repOverlap(TCRbeta, 'hamm', 'aa', .vgene = T, .verbose = T), .title = 'TCRbeta - (ave)-intersection', .labs = '')



TCRbeta.top <- top.cross(.data = TCRbeta, .n = seq(500, 10000, 500), .verbose = F, .norm = T)
top.cross.plot(TCRbeta.top)



#### GENE USAGE
TCRalpha.jusage <- geneUsage(TCRalpha, HUMAN_TRBJ)
vis.gene.usage(TCRalpha.jusage, .main = 'Samples J-usage TCR alpha', .dodge = T)

TCRalpha.vusage <- geneUsage(TCRalpha, HUMAN_TRBV)
vis.gene.usage(TCRalpha.vusage, .main = 'Samples V-usage', .dodge = T)

TCRbeta.jusage <- geneUsage(TCRbeta, HUMAN_TRBJ)
vis.gene.usage(TCRbeta.jusage, .main = 'Samples J-usage TCR beta', .dodge = T)

TCRbeta.vusage <- geneUsage(TCRbeta, HUMAN_TRBV)
vis.gene.usage(TCRbeta.vusage, .main = 'Samples V-usage', .dodge = T)


### PCA
TCRalpha.pca <- pca.segments(TCRalpha, .do.plot = F) 
vis.pca(pca.segments(TCRalpha, .do.plot = F, .genes = HUMAN_TRBV), .groups = list(GroupA = c(1,2), GroupB = c(3,4)))


