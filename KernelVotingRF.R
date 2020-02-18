# Required package
library(ranger)

### Functions for the kernel voting RF approach

# General function for the kernel voting RF approach
# A multivariate Gaussian kernel is used.

# formula : Object of class formula
# data : the training data sets (in a data.frame)
# dataTest : the testing data sets (in a data.frame); covariates names must be identical to the training ones
# ntree : number of trees
# ntreeToKeep : number of trees with the highest scores we keep
# ncores : number of cores to use
# rule : which rule to use for the bandwidth computation (values are "quantile" or "Silverman")
# alpha : the quantile order, if the bandwidth is computed thanks to quantiles

kernelVoting <- function(formula, data, dataTest, ntree, ntreeToKeep=ntree,
                        ncores=7, rule="quantile", alpha=1, ...){

  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data"), names(mf))
  mf <- mf[c(1L,m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  responseValues <- model.response(mf)

  covMatrix <- get_all_vars(formula, data=data)[ , names(get_all_vars(formula, data=data)) != as.character(formula[[2]])]
  covMatrixTest <- dataTest[,colnames(covMatrix)]
  
  nTrain <- nrow(data)
  nTest <- nrow(dataTest)

  # classic RF training
  rf.ranger <- ranger(formula, data=data, keep.inbag = TRUE, num.trees = ntree, num.threads = ncores, ...)

  predRF <- predict(object = rf.ranger, data = dataTest)$predictions

  # Recover the inbag matrix per tree
  inbag <- simplify2array(rf.ranger$inbag.counts)
  
  ### Instead of RF similarity we use a multidimensional kernel
  
  pred.Trainresponse <- predict(object = rf.ranger, data = data, predict.all = TRUE, type = "response")
  predTrainingResponse <- pred.Trainresponse$predictions
  
  pred.Testresponse <- predict(object = rf.ranger, data = dataTest, predict.all = TRUE, type = "response")
  predTestResponse <- pred.Testresponse$predictions
  
  rightMatrix <- matrix(0, nrow=nTrain, ncol=ntree)
  
  # For each tree, check whether the prediction on training is good or not
  for(k in 1:ntree){
    
    toChange <- factor(rf.ranger$forest$levels[predTrainingResponse[,k]]) == responseValues
    rightMatrix[toChange,k] <- 1 
    
  }
  
  # Recover the matrix of out-of-bag identifier (0=inbag, 1=out-of-bag)
  matrixOOB <- matrix(0, nrow=nTrain, ncol=ntree)
  
  for(k in 1:ntree){
    matrixOOB[inbag[,k]==0,k] <- 1
  }
  
  ### Tree weight computation for each test instances
  #  w_k = sum_j^nTrain( 1_xjOOB_inN * K(xj, x*) * indicatriceBienPredit ) / sum( 1_xjOOB_inN * K(xj, x*))
  
  # KMatrix will store the distances (kernel values)
  KMatrix <- matrix(NA, nrow=nTest, ncol=nTrain)
  
  # Compute bandwidth matrix
  for(i in 1:nTest){
    
    if( rule == "Silverman" ){
      h <- sapply(1:ncol(covMatrix), function(x) (4/(3*nTrain))^(1/5) * sd(covMatrix[,x]) )
      matH <- diag(h)
    } else if (rule == "quantile") {
      h <- sapply(1:ncol(covMatrix), function(x) quantile(abs(covMatrix[,x]-covMatrixTest[i,x]), alpha))
      matH <- diag(h)
    }
    
    # Center each covariate with the observed data to compute the kernel values
    matXCentered <- sapply(1:ncol(covMatrix), function(x) covMatrix[,x]-covMatrixTest[i,x])
    
    KMatrix[i,] <- KernelMultiGauss(matXCentered, matH)
    
  }
  
  # To store the tree weights for each test instance
  weights <- matrix(NA, nrow=nTest, ncol=ntree)
  
  for(i in 1:nTest){
    
    weightsPredXTest <- rep(NA, ntree)
    
    for(k in 1:ntree){
      
      denomPos <- sum(matrixOOB[,k] * KMatrix[i,])
      
      if(denomPos==0){
        
        weightsPredXTest[k] <- 0
        
      } else{
        
        weightsPredXTest[k] <- sum(matrixOOB[,k] * KMatrix[i,] * rightMatrix[,k]) / denomPos
        
      }
      
    }
    
    weights[i,] <- weightsPredXTest/sum(weightsPredXTest)
    
  }
  
  # We weight each tree prediction depending on the matrix "weights"
  # we compute weighted proportions
  
  matrixPropWeighted <- matrix(NA, nrow=nTest, ncol=length(rf.ranger$forest$levels))
  
  for(i in 1:nTest){
    
    vectorWeightedProp <- rep(NA, nlevels(responseValues))
    
    bestTreeIndex <- order(weights[i,], decreasing = TRUE)[1:ntreeToKeep]
    
    cptr <- 0
    
    for(j in rf.ranger$forest$levels){
      
      cptr <- cptr + 1
      
      if(sum(weights[i,bestTreeIndex])==0){
        
        vectorWeightedProp[cptr] <- 0
        
      } else {
        
        vectorWeightedProp[cptr] <- sum( weights[i,bestTreeIndex] * (rf.ranger$forest$levels[predTestResponse[i,bestTreeIndex]]==j) ) / sum(weights[i,bestTreeIndex])
        
      }
      
    }
    
    matrixPropWeighted[i,] <- vectorWeightedProp
    
  }
  
  # We predict as the weighted majority rule
  predKVtmp <- apply(matrixPropWeighted, 1, which.max)
  
  predKV <- levels(responseValues)[predKVtmp]
  
  return(list(prediction = predKV, weightedPropMatrix = matrixPropWeighted, weightsTreeMatrix = weights, predictionRF = predRF))
  
}


# Function for the kernel computation (where matX is centered in the observed data)
# matX : the covariate matrix centered in the observed data
# matH : the bandwidth matrix

KernelMultiGauss <- function(matX, matH){
  
  # Checkings / Initialisations
  if(is.null(dim(matX))){
    matX <- matrix(matX, nrow=1)
  }
  
  if(any(diag(matH)==0)){
    cat("Zero detected","\n")
    indOfInterest <- which(diag(matH)==0)
    for(k in indOfInterest){
      matH[k,k] <- min(abs(matX[,k][abs(matX[,k])!=0]))
    }
  }
  
  if(any(is.na(matH))){
    warnings("Kernel's bandwidth matrix H has NA !")
  }
  
  d <- nrow(matX)
  nCov <- ncol(matX)
  invMatH <- solve(matH^2)
  return( sapply(1:d, function(x)  exp(-0.5 * matX[x,] %*% invMatH %*% matX[x,])) )
  
}