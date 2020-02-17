### Functions for the training of unidimensional kernel RF

### Function for the computation of multidimensional Gaussian kernel
# x : the covariates, assumed centered in the observed data
# h : the bandwidth value

KernelMultiGauss <- function(matX, matH){

  if(is.null(dim(matX))){
    matX <- matrix(matX, nrow=1)
  }
  
  if(any(diag(matH)==0)){
    indOfInterest <- which(diag(matH)==0)
    for(k in indOfInterest){
      matH[k,k] <- min(abs(matX[k,][abs(matX[k,])!=0]))
    }
  }
  
  if(any(is.na(matH))){
    warnings("Kernel's bandwidth matrix H has NA !")
  }
  d <- nrow(matX)
  nCov <- ncol(matX)
  invMatH <- solve(matH^2)
  sapply(1:d, function(x)  exp(-0.5 * matX[x,] %*% invMatH %*% matX[x,] ) )
}

# To retrieve the RF predictions (by giving weights equal to 1 for each training data)
KernelRF <- function(matX, matH){
  res <- rep(1,nrow(matX))
  return(res)
}


### Function to compute weighted proportions

propWeighted <- function(y, weights){
  
  nbrMod <- nlevels(y)
  lvl <- levels(y)
  propW <- rep(NA, nbrMod)
  
  for(k in 1:nbrMod){
    condEg <- y == lvl[k]
    numCondEg <- as.numeric(condEg)
    cat(length(numCondEg),"\n")
    cat(length(weights),"\n\n")
    propW[k] <- sum(numCondEg*weights)/sum(weights)
  }
  
  names(propW) <- levels(y)
  return(propW)
}


### Function to find the best split for a given covariate using the multidim. kernel as weights

# y : the response values
# xInt : the value for a given covariate
# obsInt : the corresponding value for the observed data
# weights : the weights (multidimensional kernel values)
# minNodeSize : the threshold number of data below which the node is terminal

critLocalMultiDim <- function(y, xInt, obsInt, weights, minNodeSize){
  
  # Are we in a terminal node?
  theEnd <- TRUE      # By default, yes
  critValue <- 0      # Gain initialization
  cutValue <- NA      
  posStar <- TRUE     # indicate the position of the observed data toward the cut value, TRUE = x*<= cutValue

  
  N <- length(xInt)
  xOrd <- order(xInt)
  xSort <- xInt[xOrd]
  ySort <- y[xOrd]
  nbrMod <- nlevels(ySort)
  lvl <- levels(ySort)
  
  # If note enough data we stop
  if(N < minNodeSize){
    return( list(critValue = critValue, cutValue= cutValue, posStar= posStar, theEnd = theEnd) )
  }
  
  # If all covariate values are identical
  if( all(xInt == xInt[1] ) ) {
    
    theEnd <- FALSE
    return( list(critValue = critValue, cutValue= cutValue, posStar= posStar, theEnd = theEnd) )
    
  }
  
  # What are the possible cutValues
  v <- unique(xSort)
  v <- v[-length(v)] + diff(v)/2
  
  res <- NA # contiendra la valeur du crit?re temporaire
  
  # Mother node impurity
  
  denomMere <- sum(weights)
  
  nMere <- rep(NA, nbrMod)
  
  for(k in 1:nbrMod){
    
    condEgMere <- ySort == lvl[k]
    numcondEgMere <- as.numeric(condEgMere)
    
    nMere[k] <- sum(numcondEgMere*weights)
    
  }
  
  if(denomMere!=0){
    pMere <- nMere/denomMere
  } else{
    pMere <- table(ySort)/length(ySort)
  }
  
  critMere <- sum(pMere*(1-pMere))
  
  # For the first cut value
  i <- 1
  
  wix <- (obsInt <= v[i])
  indLeft <- xSort <= v[i]
  indRight <- !indLeft
  yLeft <- ySort[indLeft]
  yRight <- ySort[indRight]
  
  denomLeft <- sum(weights[indLeft], na.rm=TRUE)
  denomRight <- sum(weights[indRight], na.rm=TRUE)
  
  nLeft <- rep(NA, nbrMod)
  nRight <- rep(NA, nbrMod)
  
  for(k in 1:nbrMod){
    
    condEgLeft <- yLeft == lvl[k]
    condEgRight <- yRight == lvl[k]
    
    numCondEgLeft <- as.numeric(condEgLeft)
    numCondEgRight <- as.numeric(condEgRight)
    
    nLeft[k] <- sum(numCondEgLeft*weights[indLeft])
    nRight[k] <- sum(numCondEgRight*weights[indRight])
    
  }
  
  if(denomLeft!=0){
    pLeft <- nLeft/denomLeft
  }else{
    pLeft <- rep(0,nbrMod)
  }
  
  if(denomRight!=0){
    pRight <- nRight/denomRight
  }else{
    pRight <- rep(0,nbrMod)
  }
  
  if(denomLeft == 0 && denomRight==0){ # if I have a null kernel, I take the CART criterion
    pLeft <- table(ySort[xSort <= v[i]])/sum(xSort <= v[i])
    pRight <- table(ySort[xSort > v[i]])/sum(xSort > v[i])
    res <- sum(xSort <= v[i])/length(xSort)*sum(pLeft*(1-pLeft)) + sum(xSort > v[i])/length(xSort)*sum(pRight*(1-pRight))
  } else{
    res <- denomLeft/(denomMere)*sum(pLeft*(1-pLeft)) + denomRight/(denomMere)*sum(pRight*(1-pRight))
  }
  
  res <- critMere - res
  
  # Update
  if ( res > critValue ) {
    critValue <- res
    cutValue <- v[i]
    posStar <- wix
    theEnd <- FALSE
  }
  
  i <- i+1
  
  while(i <= length(v)){
    
    wix <- (obsInt <= v[i])
    indLeft <- xSort <= v[i]
    indRight <- !indLeft
    yLeft <- ySort[indLeft]
    yRight <- ySort[indRight]
    
    denomLeft <- sum(weights[indLeft])
    denomRight <- sum(weights[indRight])
    
    nLeft <- rep(NA, nbrMod)
    nRight <- rep(NA, nbrMod)
    
    for(k in 1:nbrMod){
      
      condEgLeft <- yLeft == lvl[k]
      condEgRight <- yRight == lvl[k]
      
      numCondEgLeft <- as.numeric(condEgLeft)
      numCondEgRight <- as.numeric(condEgRight)
      
      nLeft[k] <- sum(numCondEgLeft*weights[indLeft])
      nRight[k] <- sum(numCondEgRight*weights[indRight])
      
    }
    
    if(denomLeft!=0){
      pLeft <- nLeft/denomLeft
    }else{
      pLeft <- rep(0,nbrMod)
    }    
    
    if(denomRight!=0){
      pRight <- nRight/denomRight
    }else{
      pRight <- rep(0,nbrMod)
    }
    
    if(denomLeft == 0 && denomRight==0){ # !!! Si on a un noyau nul partout je reprends le crit?re CART
      pLeft <- table(ySort[xSort <= v[i]])/sum(xSort <= v[i])
      pRight <- table(ySort[xSort > v[i]])/sum(xSort > v[i])
      res <- sum(xSort <= v[i])/length(xSort)*sum(pLeft*(1-pLeft)) + sum(xSort > v[i])/length(xSort)*sum(pRight*(1-pRight))
    } else{
      res <- denomLeft/(denomMere)*sum(pLeft*(1-pLeft)) + denomRight/(denomMere)*sum(pRight*(1-pRight))
    }
    
    res <- critMere - res
    
    # Update
    if ( res > critValue ) {
      critValue <- res
      cutValue <- v[i]
      posStar <- wix
      theEnd <- FALSE
    }
    
    i <- i+1
    
  }
  
  return( list(critValue = critValue, cutValue= cutValue, posStar= posStar, theEnd = theEnd) )
  
}


### Build an unidimensional kernel tree

# x : the training explanatory variables
# y : the training response
# obs : the observed data
# mtry : the number of covariate to sample
# minNodeSize : the threshold below which the node is terminal (1 must be used, because we stock the tree path with Nstock)
# alpha : the quantile order
# bootstrap : do we use bootstrap samples for the tree construction
# Nstock : the value below which we store the tree (path) evolution, so that we can recover any prediction for different minNodeSize value
# hfixe : does the weights are computed only once at the root (TRUE), or are they updated at each internal node (FALSE)?
# whichKernel : which type of kernel to use ("MultiGauss" or "RF" to retrieve the usual RF prediction)
# covWeights : to add covariate weights for their sampling


treeLocalMultiDim <- function(x, y, obs, mtry, minNodeSize, alpha=1, bootstrap=FALSE, Nstock, hfixe,
                      whichKernel, covWeights=NULL){

  covIDsEvol <- c()
  cutValuesEvol <- c()
  supOrInfEvol <- c()
  
  obs <- as.numeric(obs)
  
  # If only one covariate
  if(is.vector(x)) x <- as.matrix(x)
  
  N <- length(y)
  
  if( bootstrap == TRUE ){
    indBoot <- sample(1:N, N, replace=TRUE)
  }
  else if( bootstrap == FALSE ){
    indBoot <- c(1:N)
  }
  
  yBoot <- y[indBoot]
  xBoot <- x[indBoot, , drop=FALSE]
  yOld <- yBoot
  
  stockMatrix <- NULL
  propMatrix <- NULL 
  isLimit <- FALSE

  # the bandwidth is taken as a quantile of order alpha
  h <- sapply(1:ncol(xBoot), function(x) quantile(abs(xBoot[,x]-obs[x]), alpha))
  matH <- diag(h)
  
  matX <- sapply(1:ncol(xBoot), function(x) xBoot[,x]-obs[x])
  weights <- eval(parse(text = paste("Kernel", whichKernel, "(matX, matH)", sep="") ) )
  
  theEnd <- FALSE # Initalisation
  
  if( length(yBoot) < minNodeSize ) theEnd <- TRUE
  
  if(nrow(unique(xBoot, MARGIN=1))==1){
    
    proportionFinale <- propWeighted(y = yBoot, weights = weights)
    allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
    
    return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                 propMatrix = propMatrix, covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
    
  }
  
  critValue <- 0  # Initialisation
  
  while( (!theEnd) && (nlevels(as.factor(as.numeric(yBoot))) != 1) ){
    
    if( !hfixe ){
      
      h <- sapply(1:ncol(xBoot), function(x) quantile(abs(xBoot[,x]-obs[x]), alpha))
      matH <- diag(h)
      
      matX <- sapply(1:ncol(xBoot), function(x) xBoot[,x]-obs[x])
      weights <- eval(parse(text = paste("Kernel", whichKernel, "(matX, matH)", sep="") ) )
      
    }
    
    # Sampling of covariates
    covToTry <- sample(1:ncol(x), mtry, replace = FALSE, prob = covWeights)
    
    # We start with the first covariate
    imax <- covToTry[1]
    
    # If all covariates are identical, stop
    if( all(xBoot[,covToTry[1]] == xBoot[1,covToTry[1]] ) ) {
      
      proportionFinale <- propWeighted(y = yBoot, weights = weights)
      allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
      
      return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                   propMatrix = propMatrix, covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
      
    }
    
    resCrit <- critLocalMultiDim(yBoot, xBoot[,imax], obs[imax], weights = weights, minNodeSize = minNodeSize)
    
    # If this is the end
    if(resCrit$theEnd){
      
      proportionFinale <- propWeighted(y = yBoot, weights = weights)
      allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
      
      return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                   propMatrix = propMatrix, covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
      
    }
    
    critValue <- resCrit$critValue
    cutValue <- resCrit$cutValue
    posStar <- resCrit$posStar
    
    # Let's look at all other covariates
    for(i in covToTry[-1]){
      
      if( all( xBoot[,i] == xBoot[1,i] ) ) {
        
        proportionFinale <- propWeighted(y = yBoot, weights = weights)
        allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
        
        return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                     propMatrix = propMatrix, covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
        
      }
      
      resCrit <- critLocalMultiDim(yBoot, xBoot[,i], obs[i], weights = weights, minNodeSize = minNodeSize)
      
      if(resCrit$theEnd){
        
        proportionFinale <- propWeighted(y = yBoot, weights = weights)
        allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
        
        return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                     propMatrix = propMatrix, covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
        
      }
      
      # Update
      if(resCrit$critValue > critValue){
        critValue <- resCrit$critValue
        cutValue <- resCrit$cutValue
        posStar <- resCrit$posStar
        imax <- i
      }
      
    }
    
    if(critValue==0){
      
      proportionFinale <- propWeighted(y = yBoot, weights = weights)
      allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
      
      return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                   propMatrix = propMatrix, covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
      
    }
    
    covIDsEvol <- c(covIDsEvol, imax)
    cutValuesEvol <- c(cutValuesEvol, cutValue)
    if(posStar) supOrInfEvol <- c(supOrInfEvol, "Inferior")
    if(!posStar) supOrInfEvol <- c(supOrInfEvol, "Superior")
    
    # Data update
    if( posStar == TRUE ){
      yOld <- yBoot
      toKeep <- xBoot[,imax] <= cutValue
      yBoot <- yBoot[toKeep]
      xBoot <- xBoot[toKeep, ,drop=FALSE]
    }
    else if( posStar == FALSE ){
      yOld <- yBoot
      toKeep <- xBoot[,imax] > cutValue
      yBoot <- yBoot[toKeep]
      xBoot <- xBoot[toKeep, ,drop=FALSE]
    }
    
    
    if( nrow(unique(xBoot, MARGIN=1)) == 1 ) {
      
      proportionFinale <- propWeighted(y = yBoot, weights = weights[toKeep])
      allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
      
      return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                   propMatrix = propMatrix, covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
      
    }
    
    weightsOld <- weights
    
    if(hfixe){
      weights <- weights[toKeep]
      matX <- matX[toKeep, ,drop=FALSE]
    }
    
    
    # We would like to stock the tree path, so that we can study influence of minNodeSize thanks to getPredMulti
    if(isLimit == FALSE && length(yBoot)<=Nstock){
      isLimit <- TRUE
      propMatrix <- matrix(NA, ceiling(length(yOld)/minNodeSize)+1, nlevels(yOld))
      colnames(propMatrix) <- levels(yOld)
      stockMatrix <- matrix(NA, ceiling(length(yOld)/minNodeSize)+1, length(yOld))
      k <- 1
    }
    
    if(isLimit){
      propMatrix[k,] <- propWeighted(y = yOld, weights = weightsOld)
      stockMatrix[k,c(1:length(yOld))] <- as.numeric(levels(yOld))[yOld]
      k <- k+1
    }
    
    
  }
  
  if(!hfixe){
    h <- sapply(1:ncol(xBoot), function(x) quantile( abs(xBoot[,x] - obs[x]), alpha) )
    matH <- diag(h)
    
    # Compute data weights
    matX <- sapply(1:ncol(xBoot), function(x) xBoot[,x]-obs[x] )
    weights <- eval(parse(text = paste("Kernel", whichKernel, "(matX, matH)", sep="") ) )
  }
  
  if(isLimit){
    propMatrix[k,] <- propWeighted(y = yBoot, weights = weights)
    stockMatrix[k,c(1:length(yBoot))] <- as.numeric(levels(yBoot))[yBoot]
  }
  
  proportionFinale <- propWeighted(y = yBoot, weights = weights)
  allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
  
  return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
               propMatrix = propMatrix, covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
  
}

### Function to predict thanks to the storage matrix

getPredMulti <- function(resLocalTree, Npred){
  
  if(is.null(resLocalTree$stockMatrix)){
    return(resLocalTree$allocation)
  }
  if(length(resLocalTree$stockMatrix[1,]) < Npred ){
    stop("Npred must be smaller to match the previous local tree")
  }
  
  if(Npred==1){
    res <- resLocalTree$allocation
  } else{

    nombreNonNa <- rowSums(!is.na(resLocalTree$stockMatrix))
    nombreNonNaNo0 <- nombreNonNa[nombreNonNa != 0]
    
    indToComputeTmp <- which( nombreNonNaNo0 < Npred )
    
    if(length(indToComputeTmp) == 0){
      res <- resLocalTree$allocation
    } else{
      indToCompute <- min(indToComputeTmp)
      
      proportionFinale <- resLocalTree$propMatrix[indToCompute,]
      res <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
    }
    
  }
  
  return(res)
}


### Function for the construction of forests
### It can return predictions for more than one minNodeSize value thanks to multiMinNodeSize

# x : the training explanatory variables
# y : the training response
# obs : the observed data
# mtry : the number of covariate to sample
# multiMinNodeSize : to return predictions for different minNodeSize values
# alpha : the quantile order
# ntree : the number of trees
# bootstrap : do we use bootstrap samples for the tree construction
# whichKernel : which type of kernel to use ("MultiGauss" or "RF" allowed)
# hfixe : does the weights are computed only once at the root (TRUE), or are they updated at each internal node (FALSE)?
# covWeights : to add covariate weights for their sampling

forestLocalMultiDim <- function(x, y, obs, mtry=floor(sqrt(ncol(x))), multiMinNodeSize = 1, alpha,
                                ntree = 100, bootstrap = TRUE, whichKernel, hfixe, covWeights=NULL){
  q <- length(multiMinNodeSize)
  allocationTree <- matrix(NA, q , ntree)
  
  for(i in 1:ntree){
    tree.tmp <- treeLocalMultiDim(x = x, y = y, obs = obs, mtry = mtry, minNodeSize = 1, alpha=alpha, bootstrap = bootstrap,
                          Nstock = max(multiMinNodeSize), hfixe = hfixe, whichKernel = whichKernel, covWeights=covWeights)
    l <- 0
    for(j in multiMinNodeSize){
      l <- l+1
      allocationTree[l,i] <- getPredMulti(tree.tmp, j)
    }
  }
  row.names(allocationTree) <- multiMinNodeSize
  predictionForest <- sapply(c(1:q), function(X) sample(names(table(allocationTree[X,])[which(table(allocationTree[X,])==max(table(allocationTree[X,])))]),1) )
  names(predictionForest) <- multiMinNodeSize
  return( list(predictionForest = predictionForest, allocationTree = allocationTree) )
}