install.packages("supclust")
library(supclust)
data(leukemia, package="supclust")

class(leukemia)  
dim(leukemia.x)
length(leukemia.y)


## Generating random test data: 3 observations and 250 variables (genes)
set.seed(724)
xN <- matrix(rnorm(750), nrow = 3, ncol = 250)

## Fitting Wilma
fit  <- wilma(leukemia.x, leukemia.y, noc = 3, trace = 1)

## Working with the output
fit
summary(fit)
plot(fit)
fitted(fit)

## Fitted values and class predictions for the training data
predict(fit, type = "cla")
predict(fit, type = "fitt")

## Predicting fitted values and class labels for test data
predict(fit, newdata = xN)
predict(fit, newdata = xN, type = "cla", classifier = "nnr", noc = c(1,2,3))
predict(fit, newdata = xN, type = "cla", classifier = "dlda", noc = c(1,3))
predict(fit, newdata = xN, type = "cla", classifier = "logreg")
predict(fit, newdata = xN, type = "cla", classifier = "aggtrees")