# Required packages

library(MASS)   # for LDA
library(abcrf)  # for data reading

########## Population genetics example ##########

set.seed(1)

refTable <- readRefTable("reftable.bin", header="header.txt")

data("snp.obs")  # Two real observed data (from the abcrf package)

### Data formulation

N <- length(refTable$scenarios) #70100 total data

N.train <- 10000 # number of training data

# We would like to check the prediction accuracy in a difficult prediction area
# (thanks to the LDA axes), inside the ]-1,1[x]-1,1[ box

N.test <- N - N.train

indicesTrain <- sample(1:N, N.train, replace=FALSE)

indicesTestRemaining <- c(1:N)[-indicesTrain]

ySNP <- refTable$scenarios[indicesTrain]
xSNP <- as.data.frame(refTable$stats[indicesTrain,])

ytestSNP <- refTable$scenarios[indicesTestRemaining]
xtestSNP <- as.data.frame(refTable$stats[indicesTestRemaining,])

trainSNP <- data.frame(mod = ySNP, xSNP)

table(ySNP)

### We use the LDA axes to find the closest data to the difficult observed one

model.ldaSNP <- lda(mod~. , trainSNP)

xSNP.lda <- predict(model.ldaSNP)$x

colnames(xtestSNP) <- colnames(model.ldaSNP$means)

xtestSNP.lda <- predict(model.ldaSNP, xtestSNP)$x

colnames(snp.obs) <- colnames(model.ldaSNP$means)

xObs.lda <- predict(model.ldaSNP, snp.obs)$x

train.lda <- data.frame(mod=ySNP, xSNP.lda)

### Data visualization

indMel <- sample(c(1:length(ySNP)), length(ySNP), replace=FALSE)

plot(xSNP.lda[indMel,1], xSNP.lda[indMel,2], col=as.numeric(ySNP[indMel])+1, pch="*", cex=2)

points(xObs.lda[,1], xObs.lda[,2], col=c("orange","magenta"), lwd=5)

### Testing set (nTest sampled inside the desired box)

nTest <- 500

nTestScen1 <- ceiling(nTest/3)
nTestScen2 <- ceiling(nTest/3)
nTestScen3 <- nTest - 2*ceiling(nTest/3)

# Keep only test data inside the ]-1,1[x]-1,1[ box
idxWindows <- c(1:N.test)[-1<xtestSNP.lda[,1] & xtestSNP.lda[,1]<1 & -1<xtestSNP.lda[,2] & xtestSNP.lda[,2]<1]
length(idxWindows)

# Sample 500 test data with equal proportions in each class

indiceTestScen1 <- sample(idxWindows[ytestSNP[idxWindows] == 1], nTestScen1, replace = FALSE)
indiceTestScen2 <- sample(idxWindows[ytestSNP[idxWindows] == 2], nTestScen2, replace = FALSE)
indiceTestScen3 <- sample(idxWindows[ytestSNP[idxWindows] == 3], nTestScen3, replace = FALSE)

indiceTest <- c(indiceTestScen1, indiceTestScen2, indiceTestScen3)

x.test <- xtestSNP[indiceTest,]
x.testlda <- xtestSNP.lda[indiceTest,]
y.test <- ytestSNP[indiceTest]

testSNP <- data.frame(mod = y.test, x.test)

points(x = x.testlda[,1], y = x.testlda[,2], col=as.numeric(y.test)+1)

table(y.test)

### Notatons to be in agreement with the other examples

n <- N.train
classe <- ySNP
x.train <- xSNP
data.train <- trainSNP
data.test <- testSNP
classeTest <- y.test

# Name consistency
colnames(x.train) <- colnames(data.train)[-1]
colnames(x.test) <- colnames(data.train)[-1]



###################################################################
##################### Performance comparisons #####################


# Number of cores
ncores <- 7


###################################
########## Breiman's random forests

library(ranger)

### Bagging

baggedRf <- ranger(formula = mod~., data = data.train, num.trees = 100, 
                   mtry = dim(x.train)[2], num.threads = ncores)

resBagging <- predict(object = baggedRf, data = x.test, num.threads = ncores)

mean(resBagging$predictions != classeTest)


###################################
############### Random forest

classicRF <- ranger(formula = mod~., data = data.train,
                    num.trees = 100, num.threads = ncores)

resRF <- predict(object = classicRF, data = x.test, num.threads = ncores)

mean(resRF$predictions != classeTest)


##############################################
############### Local variable importance RF

source("LocalVarImpRF.R")


### When bagged trees are used for the first forest

# Step 1: train a tree ensemble
rf.ranger <- ranger(mod ~ ., data = data.train, num.trees = 100,
                    num.threads = ncores, mtry = dim(x.train)[2])

# Step 2: Compute local importance of covariates for each observed data
impxStd <- matrix(NA, nrow = nTest, ncol=dim(x.train)[2])
for(i in 1:nTest){
  impxStd[i,] <- LocalVarImp(rf.ranger, x.test[i,,drop=FALSE])
}

# Step 3 : Use these weights during the sampling of covariates
predLVIRF1 <- factor(c(), levels=levels(classe))

for(i in 1:nTest){
  rf.local.ranger <- ranger(mod ~ ., data = data.train, num.trees = 100,
                            split.select.weights = impxStd[i,], num.threads = ncores)
  predLVIRF1[i] <- predict(rf.local.ranger, data=x.test[i,,drop=FALSE])$predictions
}

mean(predLVIRF1 != classeTest)


### When classic RF are used for the first forest

# Step 1: train a tree ensemble
rf.ranger <- ranger(mod ~ ., data = data.train, num.trees = 100,
                    num.threads = ncores)

# Step 2: Compute local importance of covariates for each observed data
impxStd <- matrix(NA, nrow = nTest, ncol=dim(x.train)[2])
for(i in 1:nTest){
  impxStd[i,] <- LocalVarImp(rf.ranger, x.test[i,,drop=FALSE])
}

# Step 3 : Use these weights during the sampling of covariates
predLVIRF2 <- factor(c(), levels=levels(classe))

for(i in 1:nTest){
  rf.local.ranger <- ranger(mod ~ ., data = data.train, num.trees = 100, 
                            split.select.weights = impxStd[i,], num.threads = ncores)
  predLVIRF2[i] <- predict(rf.local.ranger, data=x.test[i,,drop=FALSE])$predictions
}

mean(predLVIRF2 != classeTest)

###############################################################################
############### Case Specific Random Forest (function csrf from ranger package)

# Nmin=5
predCsrf5 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test),
                  params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 5),
                  params2 = list(num.trees=100))

mean(predCsrf5 != classeTest)

# Nmin=10
predCsrf10 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test),
                   params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 10),
                   params2 = list(num.trees=100))

mean(predCsrf10 != classeTest)

# Nmin=50
predCsrf50 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test),
                   params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 50),
                   params2 = list(num.trees=100))

mean(predCsrf50 != classeTest)

# Nmin=150
predCsrf150 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test),
                    params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 150),
                    params2 = list(num.trees=100))

mean(predCsrf150 != classeTest)

# Nmin=250
predCsrf250 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test),
                    params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 250),
                    params2 = list(num.trees=100))

mean(predCsrf250 != classeTest)

# Nmin=350
predCsrf350 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test),
                    params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 350),
                    params2 = list(num.trees=100))

mean(predCsrf350 != classeTest)


########################################################
############### Kernel voting

source("KernelVotingRF.R")

# alpha = 1
resKVRF1 <- kernelVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), 
                         ntree = 100, ncores = ncores, rule = "quantile", alpha = 1)

mean(resKVRF1$prediction != classeTest)

# alpha = 0.75
resKVRF2 <- kernelVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), 
                         ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.75)

mean(resKVRF2$prediction != classeTest)

# alpha = 0.5
resKVRF3 <- kernelVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), 
                         ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.5)

mean(resKVRF3$prediction != classeTest)

# alpha = 0.25
resKVRF4 <- kernelVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), 
                         ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.25)

mean(resKVRF4$prediction != classeTest)


########################################################
############### Nearest-neighbors followed by classic RF

# Compute the mad for the standardized Euclidean distance computation
madInit <- apply(X = x.train, 2, mad)

x.train.mat <- as.matrix(x.train) # to make feasible the computation of the distances
x.test.mat <- as.matrix(x.test)

# 1000 nearest neighbors
K <- 1000

predNNRF1 <- factor(c(), levels=levels(classe))

for(i in 1:nTest){
  
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.train.mat[X,]-x.test.mat[i,])/madInit )^2) ) )
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.train[toKeep,]
  
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  
  predNNRF1[i] <- predict(rfNN, data=data.frame(x.test[i,,drop=FALSE]), num.threads=ncores)$predictions
  
}

mean(predNNRF1 != classeTest)

# 1500 nearest neighbors
K <- 1500

predNNRF2 <- factor(c(), levels=levels(classe))

for(i in 1:nTest){
  
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.train.mat[X,]-x.test.mat[i,])/madInit )^2) ) )
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.train[toKeep,]
  
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  
  predNNRF2[i] <- predict(rfNN, data=data.frame(x.test[i,,drop=FALSE]), num.threads=ncores)$predictions
  
}

mean(predNNRF2 != classeTest)

# 2500 nearest neighbors
K <- 2500

predNNRF3 <- factor(c(), levels=levels(classe))

for(i in 1:nTest){
  
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.train.mat[X,]-x.test.mat[i,])/madInit )^2) ) )
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.train[toKeep,]
  
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  
  predNNRF3[i] <- predict(rfNN, data=data.frame(x.test[i,,drop=FALSE]), num.threads=ncores)$predictions
  
}

mean(predNNRF3 != classeTest)