## November 2020 ##
## Script to classify cancer samples in different subgroups following RANDOM FOREST MODEL
# https://www.bioconductor.org/packages/devel/bioc/vignettes/MLSeq/inst/doc/MLSeq.pdf

install.packages("BiocManager")
library(BiocManager)

BiocManager::install("ROCR")
#BiocManager::install("DeSousa2013")
BiocManager::install("MLSeq")
BiocManager::install("e1071")

#library(ROCR)
#library(DeSousa2013)
library(MLSeq)
library(S4Vectors)
library(DESeq2)
library(e1071)

set.seed(100)

# Example
#COUNT TABLE
filepath <- system.file("extdata/cervical.txt", package = "MLSeq")
cervical <- read.table(filepath, header=TRUE)

#USE TARGET SAMPLES relapse or remission
colnames(tallcount)

tallml_rem<-tallcount[,colnames(tallcount) %in% metadata[metadata$GSEA=="TARGET" & metadata$First_event=="Remission",]$SampID]
dim(tallml_rem)[2]
#230 patients remission

#subset 23 patients
#tallml_rem<-sample(tallml_rem, 23)
#dim(tallml_rem)

tallml_rel<-tallcount[,colnames(tallcount) %in% metadata[metadata$GSEA=="TARGET" & metadata$First_event!="Remission" & metadata$First_event!="NA" & metadata$First_event!="Censored" & metadata$First_event!="Second Malignant Neoplasm",]$SampID]
dim(tallml_rel)[2]
#23 patients no remission, no NA, no censored, no second malignant neoplasm


tallml<-cbind(tallml_rem,tallml_rel)
dim(tallml)
row.names(tallml)<-tallcount$ID


### Remove Genes with zero counts in all samples and genes with NA's

tallml<-tallml[complete.cases(tallml), ]
dim(tallml)
summary(tallml$SRR3162149)

tallml<-tallml[rowSums(tallml) > 0, ]
dim(tallml)


# METADATA TABLE
#class<-DataFrame(condition = factor(rep(c("N","T"), c(29, 29))))
#class

metadata_class<-DataFrame(condition = factor(rep(c("Remission","Non_Remission"), c(dim(tallml_rem)[2], dim(tallml_rel)[2]))))
metadata_class

# We do not perform a differential expression analysis to select differentially
# expressed genes. However, in practice, DE analysis might be performed before
# fitting classifiers. Here, we selected bcat targets within the DEGs remission vs relapse.


## Using genes up and down Remission vs. non_Remission from Clinical_TCell.R script
deg_bcattarg 
data <- tallml[c(deg_bcattarg), ]
dim(data)



#10% of data will be used as test dataset, 90% used as training dataset
nTest <- ceiling(ncol(data) * 0.5)
nTest

#random selection of 10% of samples
ind <- sample(ncol(data), nTest)
ind

# Minimum count is set to 1 in order to prevent 0 division problem within
# classification models.
data.train <- as.matrix(data[ ,-ind])
data.test <- as.matrix(data[ ,ind])

dim(data.train)
dim(data.test)

classtr <- DataFrame(condition = metadata_class[-ind, ])
classts <- DataFrame(condition = metadata_class[ind, ])

#classtr <- DataFrame(condition = class[-ind, ])
#classts <- DataFrame(condition = class[ind, ])


# Now, we have 25 samples which will be used to train the classification models and have
#remaining 18 samples to be used to test the model performances. The training and test sets
#are stored in a DESeqDataSet using related functions from DESeq2 [7]. This object is then
#used as input for MLSeq.

data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                      design = formula(~condition))
table(data.trainS4$condition)

data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,
                                     design = formula(~condition))
table(data.testS4$condition)

availableMethods()


# Consider fitting a voomNSC model on deseq normalization vst
fit <- classify(data = data.trainS4, method = "voomNSC",
                preProcessing = "deseq-vst", ref = "Remission",
                control = voomControl(tuneLength = 20),
                number = 5,
                repeats = 5, classProbs = TRUE)

show(fit)
plot(fit)

trained(fit)

# The model were trained using 5-fold cross validation repeated 10 times. The number of levels
#for tuning parameter is set to 10. The length of tuning parameter space, tuneLength, may be
#increased to be more sensitive while searching optimal value of the parameters. However, this
#may drastically increase the total computation time. 


# PREDICTION
pred.svm <- predict(fit, data.testS4)
pred.svm

pred.svm <- relevel(pred.svm, ref = "Remission")
actual <- relevel(classts$condition, ref = "Remission")
tbl <- table(Predicted = pred.svm, Actual = actual)
tbl

confusionMatrix(tbl, positive = "Remission")


selectedGenes(fit)


# COMPARING THE PERFORMANCE OF CLASSIFIERS

# We selected SVM, voomDLDA and NBLDA as non-sparse classifiers and PLDA with power
# transformation, voomNSC and NSC as sparse classifiers for the comparison of fitted models.
# Raw counts are normalized using deseq method and vst transformation is used for continuous
# classifiers (NSC and SVM).

# Define control lists.
ctrl.continuous <- trainControl(method = "repeatedcv", number = 5, repeats = 10)

ctrl.discrete <- discreteControl(method = "repeatedcv", number = 5, repeats = 10,tuneLength = 10)

ctrl.voom <- voomControl(method = "repeatedcv", number = 5, repeats = 10,tuneLength = 10)


# 1. Continuous classifiers, SVM and NSC
fit.svm <- classify(data = data.trainS4, method = "svmRadial",
                    preProcessing = "deseq-vst", ref = "Remission", tuneLength = 10,
                    control = ctrl.continuous)

fit.NSC <- classify(data = data.trainS4, method = "pam",
                    preProcessing = "deseq-vst", ref = "Remission", tuneLength = 10,
                    control = ctrl.continuous)

# 2. Discrete classifiers
fit.plda <- classify(data = data.trainS4, method = "PLDA", normalize = "deseq",
                     ref = "Remission", control = ctrl.discrete)
selectedGenes(fit.plda)


fit.plda2 <- classify(data = data.trainS4, method = "PLDA2", normalize = "deseq",
                      ref = "Remission", control = ctrl.discrete)
selectedGenes(fit.plda2)

fit.nblda <- classify(data = data.trainS4, method = "NBLDA", normalize = "deseq",
                      ref = "Remission", control = ctrl.discrete)
selectedGenes(fit.nblda)

# 3. voom-based classifiers
fit.voomDLDA <- classify(data = data.trainS4, method = "voomDLDA",
                         normalize = "deseq", ref = "Remission", control = ctrl.voom)
selectedGenes(fit.voomDLDA)

fit.voomNSC <- classify(data = data.trainS4, method = "voomNSC",
                        normalize = "deseq", ref = "Remission", control = ctrl.voom)
selectedGenes(fit.voomNSC)




# 4. Predictions
# Continuous
#pred.svm <- predict(fit.svm, data.testS4)
#pred.NSC <- predict(fit.NSC, data.testS4)

# Discrete 
# PLDA
pred.plda <- predict(fit.plda, data.testS4)
pred.plda <- relevel(pred.plda, ref = "Remission")
actual <- relevel(classts$condition, ref = "Remission")
tbl <- table(Predicted = pred.plda, Actual = actual)
tbl

confusionMatrix(tbl, positive = "Remission")

#PLDA2
pred.plda2 <- predict(fit.plda2, data.testS4)
pred.plda2 <- relevel(pred.plda2, ref = "Remission")
actual <- relevel(classts$condition, ref = "Remission")
tbl <- table(Predicted = pred.plda2, Actual = actual)
tbl

confusionMatrix(tbl, positive = "Remission")

#NBLDA
pred.nblda <- predict(fit.nblda, data.testS4)
pred.nblda <- relevel(pred.nblda, ref = "Remission")
actual <- relevel(classts$condition, ref = "Remission")
tbl <- table(Predicted = pred.nblda, Actual = actual)
tbl

confusionMatrix(tbl, positive = "Remission")

# voom-based
# voomDLDA
pred.voomDLDA <- predict(fit.voomDLDA, data.testS4)
pred.voomDLDA <- relevel(pred.voomDLDA, ref = "Remission")
actual <- relevel(classts$condition, ref = "Remission")
tbl <- table(Predicted = pred.voomDLDA, Actual = actual)
tbl

confusionMatrix(tbl, positive = "Remission")

# voomNSC
pred.voomNSC <- predict(fit.voomNSC, data.testS4)
pred.voomNSC <- relevel(pred.voomNSC, ref = "Remission")
actual <- relevel(classts$condition, ref = "Remission")
tbl <- table(Predicted = pred.voomNSC, Actual = actual)
tbl

confusionMatrix(tbl, positive = "Remission")



