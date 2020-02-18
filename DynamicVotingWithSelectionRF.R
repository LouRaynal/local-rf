# Required packages
library(Rcpp)
library(ranger)

sourceCpp("similarityRfBreiman.cpp") # To compute the similarity matrix of Breiman

### Function for dynamic voting with selection random forest (DVSRF)

# data : training data set
# dataTest: the observed data (one or more), colnames have to be identical to the training
# K : the number of nearest neightbors to consider for each test data
# ntree : number of trees in the forest
# ntreeToKeep : number of trees with the highest scores we keep
# ncores : number of cores to use

dynamicVoting <- function(formula, data, dataTest, K, ntree, ntreeToKeep, ncores, ...){
  
  if(ntreeToKeep > ntree) stop("ntreeToKeep must be equal or smaller than ntree")
  
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data"), names(mf))
  mf <- mf[c(1L,m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  responseValues <- model.response(mf)
  
  nTrain <- nrow(data)
  nTest <- nrow(dataTest)
  
  if(K > nTrain) stop("Number of nearest neighbors K too large compared to the training size")
  
  # We train a classique RF
  rf.ranger <- ranger(formula, data=data, keep.inbag = TRUE, num.trees = ntree, num.threads = ncores, ...)
  
  predTestingResponseWithRF <- predict(object = rf.ranger, data = dataTest)$predictions
  
  # Recover the inbag matrix per tree
  inbag <- simplify2array(rf.ranger$inbag.counts)
  
  # We want to compute the similarity between training data and test data
  
  pred.NodeIDTrain <- predict(object = rf.ranger, data = data, predict.all = TRUE, type = 'terminalNodes')
  predTraining <- pred.NodeIDTrain$predictions # node id. for training data
  
  pred.NodeIDTest <- predict(object = rf.ranger, data = dataTest, predict.all = TRUE, type = 'terminalNodes')
  predTesting <- pred.NodeIDTest$predictions   # node id. for testing data

  similarityMeasure <- findweights(predTraining = predTraining, predTesting = predTesting, nTrain = nTrain, nTest = nTest, ntree = ntree)
  similarityMeasure3 <- similarityMeasure^3
  # This is a matrix with nTrain rows and nTest columns
  
  # Compute the margin per tree and training data
  # +1 if correct prediction, -1 otherwise
  pred.Trainresponse <- predict(object = rf.ranger, data = data, predict.all = TRUE, type = "response")
  predTrainingResponse <- pred.Trainresponse$predictions
  
  pred.Testresponse <- predict(object = rf.ranger, data = dataTest, predict.all = TRUE, type = "response")
  predTestResponse <- pred.Testresponse$predictions
  
  margineMatrix <- matrix(-1, nrow=nTrain, ncol=ntree)
  
  for(k in 1:ntree){
    toChange <- predTrainingResponse[,k] == as.numeric(responseValues)
    margineMatrix[toChange,k] <- 1
  }
  
  # Recover the matrix of out-of-bag identifier (0=inbag, 1=out-of-bag)
  matrixOOB <- matrix(0, nrow=nTrain, ncol=ntree)
  
  for(k in 1:ntree){
    matrixOOB[inbag[,k]==0,k] <- 1
  }
  
  ### Tree weight computation for each test instance
  
  # First, search for the nearest-neighbors thanks to similarityMeasure3
  # Second, compute the associated weights  
  # w_k = sum_j^nTrain( 1_xjOOB_ink * similarity(xj, x*) * marge_k_j ) / sum( 1_xjOOB_ink * similarity(xj, x*))
  
  localWeightsDV <- matrix(NA, nrow=nTest, ncol=ntree)
  
  for(i in 1:nTest){
    
    KMostSimilarIndex <- order(similarityMeasure3[,i], decreasing = TRUE)[1:K]
    weightsPredXTest <- rep(NA, ntree)
    
    for(k in 1:ntree){
      
      denomPos <- sum(matrixOOB[KMostSimilarIndex,k] * similarityMeasure3[KMostSimilarIndex,i])
      
      if(denomPos==0){
        weightsPredXTest[k] <- 0
      } else{
        weightsPredXTest[k] <- sum(matrixOOB[KMostSimilarIndex,k] * similarityMeasure3[KMostSimilarIndex,i] * margineMatrix[KMostSimilarIndex,k]) / denomPos
      }
      
    }
    
    localWeightsDV[i,] <- (weightsPredXTest+abs(min(weightsPredXTest)))/sum(weightsPredXTest+abs(min(weightsPredXTest)))
    # It contains the tree weights for each test data
  }
  
  # We compute some weighted proportions for final prediction
  matrixPropWeighted <- matrix(NA, nrow=nTest, ncol=length(rf.ranger$forest$levels))
  
  for(i in 1:nTest){
    
    vectorWeightedProp <- rep(NA, nlevels(responseValues))
    
    bestTreeIndex <- order(localWeightsDV[i,], decreasing = TRUE)[1:ntreeToKeep]
    
    for(j in 1:nlevels(responseValues)){
      
      if(sum(localWeightsDV[i,bestTreeIndex])==0){
        vectorWeightedProp[j] <- 0
      } else {
        vectorWeightedProp[j] <- sum( localWeightsDV[i,bestTreeIndex] * (as.numeric(predTestResponse[i,bestTreeIndex])==j) ) / sum(localWeightsDV[i,bestTreeIndex])
      }
      
    }
    
    matrixPropWeighted[i,] <- vectorWeightedProp
    
  }
  
  # We predict as the weighted majority rule
  predDVtmp <- apply(matrixPropWeighted, 1, which.max)
  
  predDV <- levels(responseValues)[predDVtmp]
  
  return(list(prediction = predDV, weightedPropMatrix = matrixPropWeighted, weightsTreeMatrix = localWeightsDV, predictionRFClassic = predTestingResponseWithRF))
  
}