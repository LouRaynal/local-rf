# Required packages

library(mvtnorm)
library(doParallel)

########## Unbalanced Gaussian example ##########

set.seed(1)

### Training data generation

# Model probabilities (4 balanced classes)
pi0 <- 0.40
pi1 <- 0.40
pi2 <- 0.10
pi3 <- 0.10

# Gaussian dimension
l <- 20

# Gaussian parameters
mu0 <- c(c(0.8,3), rep(c(1,2.5), l/2-1)) ; Sigma0 <- diag(c(c(2,1), rep(c(3,1),l/2-1)) )
mu1 <- c(c(3.2,3), rep(c(2.5,2.5), l/2-1)) ; Sigma1 <- diag(c(c(2,1), rep(c(3,5),l/2-1)) )
mu2 <- c(c(2,1), rep(c(2,2.3),l/2-1)) ; Sigma2 <- diag( (rep(c(4,1),l/2) ) )
mu3 <- c(c(2,0), rep(c(2,1.8),l/2-1)) ; Sigma3 <- diag( (rep(c(2.5,1),l/2) ) )

# Number of training data
n <- 3000
# Sample class label
classe <- sample(x = c(0,1,2,3), size = n, replace = TRUE, prob = c(pi0,pi1,pi2,pi3))
classe <- sort(classe)

couleur <- rep("orange", n)
couleur[classe==1] <- "cyan"
couleur[classe==2] <- "purple"
couleur[classe==3] <- "green2"

n0 <- sum(classe==0)
n1 <- sum(classe==1)
n2 <- sum(classe==2)
n3 <- sum(classe==3)

# Sample from the Gaussians
x.train <- rbind(rmvnorm(n0, mu0, Sigma0),
                 rmvnorm(n1, mu1, Sigma1),
                 rmvnorm(n2, mu2, Sigma2),
                 rmvnorm(n3, mu3, Sigma3))

plot(x.train[,1], x.train[,2], col=couleur)

# Plot the Bayes frontier
len <- 100
xp <- seq(min(x.train[,1]),max(x.train[,1]), length=len)
yp <- seq(min(x.train[,2]),max(x.train[,2]), length=len)

grille <- expand.grid(z1=xp,z2=yp)

Z <- pi0*dmvnorm(grille, mu0[1:2], Sigma0[1:2,1:2])
Z <- cbind(Z,pi1*dmvnorm(grille, mu1[1:2], Sigma1[1:2,1:2]))
Z <- cbind(Z,pi2*dmvnorm(grille, mu2[1:2], Sigma2[1:2,1:2]))
Z <- cbind(Z,pi3*dmvnorm(grille, mu3[1:2], Sigma3[1:2,1:2]))
zp4 <- Z[,4] - pmax(Z[,3], Z[,2], Z[,1])
contour(xp, yp, matrix(zp4, len), add=TRUE, levels=0, drawlabels=FALSE)
zp3 <- Z[,3] - pmax(Z[,1], Z[,2], Z[,4])
contour(xp, yp, matrix(zp3, len), add=TRUE, levels=0, drawlabels=FALSE)
zp2 <- Z[,2] - pmax(Z[,1], Z[,3], Z[,4])
contour(xp, yp, matrix(zp2, len), add=TRUE, levels=0, drawlabels=FALSE)
zp1 <- Z[,1] - pmax(Z[,2], Z[,3], Z[,4])
contour(xp, yp, matrix(zp1, len), add=TRUE, levels=0, drawlabels=FALSE)

# Bayes classifier on training set
BayesClassifieur <- function(x){
  c0 <- pi0*dmvnorm(x,mean=mu0,sigma=Sigma0)
  c1 <- pi1*dmvnorm(x,mean=mu1,sigma=Sigma1)
  c2 <- pi2*dmvnorm(x,mean=mu2,sigma=Sigma2)
  c3 <- pi3*dmvnorm(x,mean=mu3,sigma=Sigma3)
  return(c(0,1,2,3)[which.max(c(c0,c1,c2,c3))])
}

predTrainingBayes <- rep(NA, n)
for(i in 1:n) predTrainingBayes[i] <- BayesClassifieur(x.train[i,])

resBayes <- table(classe, predTrainingBayes)
resBayes
(sum(resBayes)-sum(diag(resBayes))) / n
(sum(diag(resBayes))) / n

### Test data generation

# Number of testing data
nTest <- 500

classeTest <- sample(c(2,3), size=nTest, prob=c(pi2,pi3), replace=TRUE)
classeTest <- sort(classeTest)

# nTest0 <- sum(classeTest==0)
# nTest1 <- sum(classeTest==1)
nTest2 <- sum(classeTest==2)
nTest3 <- sum(classeTest==3)

x.test <- rbind(rmvnorm(nTest2, mu2, Sigma2),
                rmvnorm(nTest3, mu3, Sigma3))


couleurTest <- rep("orange", nTest)
couleurTest[classeTest==1] <- "cyan"
couleurTest[classeTest==2] <- "purple"
couleurTest[classeTest==3] <- "green2"

predTestingBayes <- rep(NA, nTest)
for(i in 1:nTest) predTestingBayes[i] <- BayesClassifieur(x.test[i,])

resTestBayes <- table(predTestingBayes,classeTest)
resTestBayes
sum(classeTest!=predTestingBayes)/nTest
sum(classeTest==predTestingBayes)/nTest

### We add 20 noise explanatory variables simulated from uniform distributions

x.trainNoised <- cbind(x.train, matrix(runif(n*l, 0, 10000), nrow=n))
x.testNoised <- cbind(x.test, matrix(runif(nTest*l, 0, 10000), nrow=nTest))

data.train <- data.frame(mod = as.factor(classe), x.train)
data.test <- data.frame(mod = as.factor(classeTest), x.test)

data.trainNoised <- data.frame(mod = as.factor(classe), x.trainNoised)

colnames(x.testNoised) <- colnames(data.trainNoised)[-1]


###################################################################
##################### Performance comparisons #####################


# Number of cores
ncores <- 7


###################################
########## Breiman's random forests

library(ranger)

### Bagging

baggedRf <- ranger(formula = mod~., data = data.trainNoised, num.trees = 100, 
                   mtry = dim(x.trainNoised)[2], num.threads = ncores)

resBagging <- predict(object = baggedRf, data = x.testNoised, num.threads = ncores)

mean(resBagging$predictions != classeTest)


###################################
############### Random forest

classicRF <- ranger(formula = mod~., data = data.trainNoised,
                    num.trees = 100, num.threads = ncores)

resRF <- predict(object = classicRF, data = x.testNoised, num.threads = ncores)

mean(resRF$predictions != classeTest)


###################################
############### LazyDRF

source("LazyDRF.R")

cl <- makeCluster(ncores)
registerDoParallel(cl)

predLazyDRF <- foreach(i=1:nTest, .combine="c", .packages = c("discretization")) %dopar% {
  
  modelLazyDRF <- LazyOptionForest(x = x.trainNoised, y = classe, obs = x.testNoised[i,], Nmin = 1, pGain = 0.9, bootstrap = TRUE, num.trees = 100)
  
  return(modelLazyDRF$pred)
  
}

stopCluster(cl)

mean(predLazyDRF != classeTest)


##############################################
############### Local variable importance RF

source("LocalVarImpRF.R")


### When bagged trees are used for the first forest

# Step 1: train a tree ensemble
rf.ranger <- ranger(mod ~ ., data = data.trainNoised, num.trees = 100,
                    num.threads = ncores, mtry = dim(x.trainNoised)[2])

# Step 2: Compute local importance of covariates for each observed data
impxStd <- matrix(NA, nrow = nTest, ncol=dim(x.trainNoised)[2])
for(i in 1:nTest){
  impxStd[i,] <- LocalVarImp(rf.ranger, x.testNoised[i,,drop=FALSE])
}

# Step 3 : Use these weights during the sampling of covariates
predLVIRF1 <- rep(NA, nTest)
for(i in 1:nTest){
  rf.local.ranger <- ranger(mod ~ ., data = data.trainNoised, num.trees = 100,
                            split.select.weights = impxStd[i,], num.threads = ncores)
  predLVIRF1[i] <- predict(rf.local.ranger, data=x.testNoised[i,,drop=FALSE])$predictions
}

mean(predLVIRF1-1 != classeTest) # !Different levels are returned by ranger! (hence the -1)


### When classic RF are used for the first forest

# Step 1: train a tree ensemble
rf.ranger <- ranger(mod ~ ., data = data.trainNoised, num.trees = 100,
                    num.threads = ncores)

# Step 2: Compute local importance of covariates for each observed data
impxStd <- matrix(NA, nrow = nTest, ncol=dim(x.trainNoised)[2])
for(i in 1:nTest){
  impxStd[i,] <- LocalVarImp(rf.ranger, x.testNoised[i,,drop=FALSE])
}

# Step 3 : Use these weights during the sampling of covariates
predLVIRF2 <- rep(NA, nTest)
for(i in 1:nTest){
  rf.local.ranger <- ranger(mod ~ ., data = data.trainNoised, num.trees = 100, 
                            split.select.weights = impxStd[i,], num.threads = ncores)
  predLVIRF2[i] <- predict(rf.local.ranger, data=x.testNoised[i,,drop=FALSE])$predictions
}

mean(predLVIRF2-1 != classeTest) # !Different levels are returned by ranger! (hence the -1)


###############################################################################
############### Case Specific Random Forest (function csrf from ranger package)

# Nmin=5
predCsrf5 <- csrf(mod~., training_data = data.trainNoised, test_data = data.frame(x.testNoised),
                  params1 = list(num.trees=100, mtry = dim(x.trainNoised)[2], min.node.size = 5),
                  params2 = list(num.trees=100))

mean(predCsrf5 != classeTest)

# Nmin=10
predCsrf10 <- csrf(mod~., training_data = data.trainNoised, test_data = data.frame(x.testNoised),
                   params1 = list(num.trees=100, mtry = dim(x.trainNoised)[2], min.node.size = 10),
                   params2 = list(num.trees=100))

mean(predCsrf10 != classeTest)

# Nmin=50
predCsrf50 <- csrf(mod~., training_data = data.trainNoised, test_data = data.frame(x.testNoised),
                   params1 = list(num.trees=100, mtry = dim(x.trainNoised)[2], min.node.size = 50),
                   params2 = list(num.trees=100))

mean(predCsrf50 != classeTest)

# Nmin=150
predCsrf150 <- csrf(mod~., training_data = data.trainNoised, test_data = data.frame(x.testNoised),
                    params1 = list(num.trees=100, mtry = dim(x.trainNoised)[2], min.node.size = 150),
                    params2 = list(num.trees=100))

mean(predCsrf150 != classeTest)

# Nmin=250
predCsrf250 <- csrf(mod~., training_data = data.trainNoised, test_data = data.frame(x.testNoised),
                    params1 = list(num.trees=100, mtry = dim(x.trainNoised)[2], min.node.size = 250),
                    params2 = list(num.trees=100))

mean(predCsrf250 != classeTest)

# Nmin=350
predCsrf350 <- csrf(mod~., training_data = data.trainNoised, test_data = data.frame(x.testNoised),
                    params1 = list(num.trees=100, mtry = dim(x.trainNoised)[2], min.node.size = 350),
                    params2 = list(num.trees=100))

mean(predCsrf350 != classeTest)


##############################################
############### Local dynamic selection RF

source("DynamicVotingWithSelectionRF.R")

# 3000 neighbors, we keep 100 best trees (all)

resDVSRF1 <- dynamicVoting(formula = mod~., data = data.trainNoised, dataTest = data.frame(x.testNoised),
                           K = 3000, nTree = 100, nTreeToKeep = 100, ncores = ncores)

mean(resDVSRF1$prediction !=  classeTest)


# 3000 neighbors, we keep 50 best trees

resDVSRF2 <- dynamicVoting(formula = mod~., data = data.trainNoised, dataTest = data.frame(x.testNoised), 
                           K = 3000, nTree = 100, nTreeToKeep = 50, ncores = ncores)

mean(resDVSRF2$prediction !=  classeTest)


########################################################
############### Kernel voting

source("KernelVotingRF.R")

# alpha = 1
resKVRF1 <- kernelVoting(formula = mod~., data = data.trainNoised, dataTest = data.frame(x.testNoised), 
                         nTree = 100, ncores = ncores, rule = "quantile", alpha = 1)

mean(resKVRF1$prediction != classeTest)

# alpha = 0.75
resKVRF2 <- kernelVoting(formula = mod~., data = data.trainNoised, dataTest = data.frame(x.testNoised), 
                         nTree = 100, ncores = ncores, rule = "quantile", alpha = 0.75)

mean(resKVRF2$prediction != classeTest)

# alpha = 0.5
resKVRF3 <- kernelVoting(formula = mod~., data = data.trainNoised, dataTest = data.frame(x.testNoised), 
                         nTree = 100, ncores = ncores, rule = "quantile", alpha = 0.5)

mean(resKVRF3$prediction != classeTest)

# alpha = 0.25
resKVRF4 <- kernelVoting(formula = mod~., data = data.trainNoised, dataTest = data.frame(x.testNoised), 
                         nTree = 100, ncores = ncores, rule = "quantile", alpha = 0.25)

mean(resKVRF4$prediction != classeTest)


########################################################
############### Nearest-neighbors followed by classic RF

# Compute the mad for the standardized Euclidean distance computation
madInit <- apply(X = x.trainNoised, 2, mad)

# 1000 nearest neighbors
K <- 1000

predNNRF1 <- rep(NA, nTest)

for(i in 1:nTest){
  
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.trainNoised[X,]-x.testNoised[i,])/madInit )^2)) )
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.trainNoised[toKeep,]
  
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  
  predNNRF1[i] <- predict(rfNN, data=data.frame(x.testNoised[i,,drop=FALSE]), num.threads=ncores)$predictions
  
}

mean(predNNRF1-1 != classeTest) # !Different levels are returned by ranger! (hence the -1)

# 1500 nearest neighbors
K <- 1500

predNNRF2 <- rep(NA, nTest)

for(i in 1:nTest){
  
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.trainNoised[X,]-x.testNoised[i,])/madInit )^2)) )
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.trainNoised[toKeep,]
  
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  
  predNNRF2[i] <- predict(rfNN, data=data.frame(x.testNoised[i,,drop=FALSE]), num.threads=ncores)$predictions
  
}

mean(predNNRF2-1 != classeTest)

# 2500 nearest neighbors
K <- 2500

predNNRF3 <- rep(NA, nTest)

for(i in 1:nTest){
  
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.trainNoised[X,]-x.testNoised[i,])/madInit )^2)) )
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.trainNoised[toKeep,]
  
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  
  predNNRF3[i] <- predict(rfNN, data=data.frame(x.testNoised[i,,drop=FALSE]), num.threads=ncores)$predictions
  
}

mean(predNNRF3-1 != classeTest)


########################################################
############### Unidimensional kernel RF

source("UniDKernelRF.R")

minNodeSize <- c(1:25) # minNodeSize values


# alpha = 1, hfixe=TRUE

cl <- makeCluster(ncores)
registerDoParallel(cl)

predUKRF1 <- foreach(i=1:nTest, .combine="rbind", .packages = "matrixStats",
                     .export = c("KernelUnif", "KernelGauss", "KernelEpan", "KernelRF")) %dopar% {
                       
                       
                       modelUKRF <- forestLocalMultiple(y = factor(classe), x = x.trainNoised, obs = x.testNoised[i,],
                                                        multiMinNodeSize = minNodeSize, alpha = 1, ntree = 100, 
                                                        bootstrap = TRUE, hfixe = TRUE, whichKernel = "Gauss")
                       return(modelUKRF$prediction)
                       
                     }

stopCluster(cl)

for(i in 1:length(minNodeSize)){
  cat("Min. node size: ", minNodeSize[i], "\n")
  cat("Error: ", mean(predUKRF1[,i]!=classeTest),"\n\n")
}


# alpha = 0.75, hfixe=TRUE

cl <- makeCluster(ncores)
registerDoParallel(cl)

predUKRF2 <- foreach(i=1:nTest, .combine="rbind", .packages = "matrixStats",
                     .export = c("KernelUnif", "KernelGauss", "KernelEpan", "KernelRF")) %dopar% {
                       
                       modelUKRF <- forestLocalMultiple(y = factor(classe), x = x.trainNoised, obs = x.testNoised[i,],
                                                        multiMinNodeSize = minNodeSize, alpha = 0.75, ntree = 100, 
                                                        bootstrap = TRUE, hfixe = TRUE, whichKernel = "Gauss")
                       return(modelUKRF$prediction)
                       
                     }

stopCluster(cl)

for(i in 1:length(minNodeSize)){
  cat("Min. node size: ", minNodeSize[i], "\n")
  cat("Error: ", mean(predUKRF2[,i]!=classeTest),"\n\n")
}


# alpha = 0.50, hfixe=TRUE

cl <- makeCluster(ncores)
registerDoParallel(cl)

predUKRF3 <- foreach(i=1:nTest, .combine="rbind", .packages = "matrixStats",
                     .export = c("KernelUnif", "KernelGauss", "KernelEpan", "KernelRF")) %dopar% {
                       
                       modelUKRF <- forestLocalMultiple(y = factor(classe), x = x.trainNoised, obs = x.testNoised[i,],
                                                        multiMinNodeSize = minNodeSize, alpha = 0.50, ntree = 100, 
                                                        bootstrap = TRUE, hfixe = TRUE, whichKernel = "Gauss")
                       return(modelUKRF$prediction)
                       
                     }

stopCluster(cl)

for(i in 1:length(minNodeSize)){
  cat("Min. node size: ", minNodeSize[i], "\n")
  cat("Error: ", mean(predUKRF3[,i]!=classeTest),"\n\n")
}



########################################################
############### Multidimensional kernel RF

source("MultiDKernelRF.R")

minNodeSize <- c(1:25) # minNodeSize values


# alpha = 1, hfixe=TRUE

cl <- makeCluster(ncores)
registerDoParallel(cl)

predMKRF1 <- foreach(i=1:nTest, .combine="rbind", .packages = "matrixStats",
                     .export = c("KernelMultiGauss", "KernelRF")) %dopar% {
                       
                       
                       modelMKRF <- forestLocalMultiDim(y = factor(classe), x = x.trainNoised, obs = x.testNoised[i,],
                                                        multiMinNodeSize = minNodeSize, alpha = 1, ntree = 100, 
                                                        bootstrap = TRUE, hfixe = TRUE, whichKernel = "MultiGauss")
                       return(modelMKRF$prediction)
                       
                     }

stopCluster(cl)

for(i in 1:length(minNodeSize)){
  cat("Min. node size: ", minNodeSize[i], "\n")
  cat("Error: ", mean(predMKRF1[,i]!=classeTest),"\n\n")
}


# alpha = 0.75, hfixe=TRUE

cl <- makeCluster(ncores)
registerDoParallel(cl)

predMKRF2 <- foreach(i=1:nTest, .combine="rbind", .packages = "matrixStats",
                     .export = c("KernelMultiGauss", "KernelRF")) %dopar% {
                       
                       
                       modelMKRF <- forestLocalMultiDim(y = factor(classe), x = x.trainNoised, obs = x.testNoised[i,],
                                                        multiMinNodeSize = minNodeSize, alpha = 0.75, ntree = 100, 
                                                        bootstrap = TRUE, hfixe = TRUE, whichKernel = "MultiGauss")
                       return(modelMKRF$prediction)
                       
                     }

stopCluster(cl)

for(i in 1:length(minNodeSize)){
  cat("Min. node size: ", minNodeSize[i], "\n")
  cat("Error: ", mean(predMKRF2[,i]!=classeTest),"\n\n")
}


# alpha = 0.5, hfixe=TRUE

cl <- makeCluster(ncores)
registerDoParallel(cl)

predMKRF3 <- foreach(i=1:nTest, .combine="rbind", .packages = "matrixStats",
                     .export = c("KernelMultiGauss", "KernelRF")) %dopar% {
                       
                       
                       modelMKRF <- forestLocalMultiDim(y = factor(classe), x = x.trainNoised, obs = x.testNoised[i,],
                                                        multiMinNodeSize = minNodeSize, alpha = 0.5, ntree = 100, 
                                                        bootstrap = TRUE, hfixe = TRUE, whichKernel = "MultiGauss")
                       return(modelMKRF$prediction)
                       
                     }

stopCluster(cl)

for(i in 1:length(minNodeSize)){
  cat("Min. node size: ", minNodeSize[i], "\n")
  cat("Error: ", mean(predMKRF3[,i]!=classeTest),"\n\n")
}