### Functions for the training of unidimensional kernel RF

### Function for the computation of various kernels
# x : the covariates, assumed centered in the observed data
# h : the bandwidth value

# Epanechnikov
KernelEpan <- function(x, h){
  if(h==0){
    h <- min(abs(x[abs(x)!=0]))
  }
  if(is.na(h)){
    warnings("Kernel's bandwidth h is NA !")
  }
  ifelse(abs(x/h)<=1, 3/(4 * h) * (1 - (x/h)^2), 0)
}

# Gaussian
KernelGauss <- function(x, h){
  if(h==0){
    h <- min(abs(x[abs(x)!=0]))
  }
  if(is.na(h)){
    warnings("Kernel's bandwidth h is NA !")
  }
  exp(-0.5*(abs(x))^2/(h^2)) *  1 / (sqrt(2*pi)*h)
}

# Uniform
KernelUnif <- function(x, h){
  if(h==0){
    h <- min(abs(x[abs(x)!=0]))
  }
  if(is.na(h)){
    warnings("Kernel's bandwidth h is NA !")
  }
  ifelse(abs(x/h)<=1, 0.5/h, 0)
}

# To retrieve the RF predictions (by giving weights equal to 1 for each training data)
KernelRF <- function(x, h){
  res <- rep(1,length(x))
  return(res)
}

### Function to find the best split for a given covariate using an unidimensional kernel

# y : the response values
# xInt : the value for a given covariate
# obsInt : the corresponding value for the observed data
# minNodeSize : the threshold number of data below which the node is terminal
# h : the bandwidth
# whichKernel : which type of kernel to use, possible values are "Epan", "Gauss", "Unif", "RF"

critLocal <- function(y, xInt, obsInt, minNodeSize, h, whichKernel){
  
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
  
  vectorDist <- xSort - as.numeric(obsInt)
  vectorKernel <- eval(parse(text = paste("Kernel", whichKernel, "(vectorDist, h)", sep="") ) )
  
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
  
  res <- NA # will contain the criterion value
  
  # Mother node impurity
  denomMoth <- sum(vectorKernel)
  
  nMoth <- rep(NA, nbrMod)
  
  for(k in 1:nbrMod){
    
    condEgMoth <- ySort == lvl[k]
    numcondEgMoth <- as.numeric(condEgMoth)
    nMoth[k] <- sum(numcondEgMoth*vectorKernel)
    
  }
  
  if(denomMoth!=0){
    pMoth <- nMoth/denomMoth
  } else{
    pMoth <- table(ySort)/length(ySort)
  }
  
  critMoth <- sum(pMoth*(1-pMoth))
  
  # For the first cut value
  i <- 1
  
  wix <- (obsInt <= v[i])
  indLeft <- xSort <= v[i]
  indRight <- !indLeft
  yLeft <- ySort[indLeft]
  yRight <- ySort[indRight]
  
  denomLeft <- sum(vectorKernel[indLeft], na.rm=TRUE)
  denomRight <- sum(vectorKernel[indRight], na.rm=TRUE)
  
  nLeft <- rep(NA, nbrMod)
  nRight <- rep(NA, nbrMod)
  
  for(k in 1:nbrMod){
    
    condEgLeft <- yLeft == lvl[k]
    condEgRight <- yRight == lvl[k]
    
    numCondEgLeft <- as.numeric(condEgLeft)
    numCondEgRight <- as.numeric(condEgRight)
    
    nLeft[k] <- sum(numCondEgLeft*vectorKernel[indLeft])
    nRight[k] <- sum(numCondEgRight*vectorKernel[indRight])
    
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
    res <- denomLeft/(denomMoth)*sum(pLeft*(1-pLeft)) + denomRight/(denomMoth)*sum(pRight*(1-pRight))
  }
  
  res <- critMoth - res
  
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
    
    denomLeft <- sum(vectorKernel[indLeft])
    denomRight <- sum(vectorKernel[indRight])
    
    nLeft <- rep(NA, nbrMod)
    nRight <- rep(NA, nbrMod)
    
    for(k in 1:nbrMod){
      
      condEgLeft <- yLeft == lvl[k]
      condEgRight <- yRight == lvl[k]
      
      numCondEgLeft <- as.numeric(condEgLeft)
      numCondEgRight <- as.numeric(condEgRight)
      
      nLeft[k] <- sum(numCondEgLeft*vectorKernel[indLeft])
      nRight[k] <- sum(numCondEgRight*vectorKernel[indRight])
      
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
    
    if(denomLeft == 0 && denomRight==0){
      pLeft <- table(ySort[xSort <= v[i]])/sum(xSort <= v[i])
      pRight <- table(ySort[xSort > v[i]])/sum(xSort > v[i])
      res <- sum(xSort <= v[i])/length(xSort)*sum(pLeft*(1-pLeft)) + sum(xSort > v[i])/length(xSort)*sum(pRight*(1-pRight))
    } else{
      res <- denomLeft/(denomMoth)*sum(pLeft*(1-pLeft)) + denomRight/(denomMoth)*sum(pRight*(1-pRight))
    }
    
    res <- critMoth - res
    
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

# y : the training response
# x : the training explanatory variables
# obs : the observed data
# mtry : the number of covariate to sample
# minNodeSize : the threshold below which the node is terminal (1 must be used, because we stock the tree path with Nstock)
# alpha : the quantile order
# bootstrap : do we use bootstrap samples for the tree construction
# Nstock : the value below which we store the tree (path) evolution, so that we can recover any prediction for different minNodeSize value
# hfixe : does the weights are computed only once at the root (TRUE), or are they updated at each internal node (FALSE)?
# whichKernel : which type of kernel to use
# covWeights : to add covariate weights for their sampling

treeLocal <- function(y, x, obs, mtry, minNodeSize, alpha=1, bootstrap=FALSE, Nstock, hfixe,
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
  isLimit <- FALSE
  
  theEnd <- FALSE # Initalisation
  
  if( length(yBoot) < minNodeSize ) theEnd <- TRUE
  
  if(nrow(unique(xBoot, MARGIN=1))==1){
    
    proportionFinale <- table(yBoot)
    allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
    
    return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                 covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
    
  }
  
  critValue <- 0 # Initialisation
  
  # the bandwidth is taken as a quantile of order alpha
  if( hfixe ) h <- sapply(1:ncol(xBoot), function(x) quantile( abs(xBoot[,x] - obs[x]), alpha) )
  
  while( (!theEnd) && (nlevels(as.factor(as.numeric(yBoot))) != 1) ){
    
    if( !hfixe ) h <- sapply(1:ncol(xBoot), function(x) quantile( abs(xBoot[,x] - obs[x]), alpha) )
    
    # Sampling of covariates
    covToTry <- sample(1:ncol(x), mtry, replace = FALSE, prob = covWeights)
    
    # We start with the first covariate
    imax <- covToTry[1]
    
    # If all covariates are identical, stop
    if( all(xBoot[,covToTry[1]] == xBoot[1,covToTry[1]] ) ) {
      
      proportionFinale <- table(yBoot)
      allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
      
      return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                   covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
      
    }
    
    resCrit <- critLocal(yBoot, xBoot[,imax], obs[imax], h = h[imax], minNodeSize = minNodeSize, whichKernel = whichKernel)
    
    # If this is the end
    if(resCrit$theEnd){
      
      proportionFinale <- table(yBoot)
      allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
      
      return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                   covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
      
    }
    
    critValue <- resCrit$critValue
    cutValue <- resCrit$cutValue
    posStar <- resCrit$posStar
    
    # Let's look at all other covariates
    for(i in covToTry[-1]){
      
      if( all( xBoot[,i] == xBoot[1,i] ) ) {
        
        proportionFinale <- table(yBoot)
        allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
        
        return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                     covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
        
      }
      
      resCrit <- critLocal(yBoot, xBoot[,i], obs[i], h = h[i], minNodeSize = minNodeSize, whichKernel = whichKernel)
      
      if(resCrit$theEnd){
        
        proportionFinale <- table(yBoot)
        allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
        
        return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                     covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
        
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
      
      proportionFinale <- table(yBoot)
      allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
      
      return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                   covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
      
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
      
      proportionFinale <- table(yBoot)
      allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
      
      return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
                   covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
      
    }
    
    # We would like to stock the tree path, so that we can study influence of minNodeSize thanks to getPred
    if(isLimit == FALSE && length(yBoot)<=Nstock){
      isLimit <- TRUE
      stockMatrix <- matrix(NA, ceiling(length(yOld)/minNodeSize)+1, length(yOld))
      k <- 1
    }
    
    if(isLimit){
      stockMatrix[k,c(1:length(yOld))] <- as.numeric(levels(yOld))[yOld]
      k <- k+1
    }
    
  }
  
  if(isLimit){
    stockMatrix[k,c(1:length(yBoot))] <- as.numeric(levels(yBoot))[yBoot]
  }
  
  proportionFinale <- table(yBoot)
  allocation <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
  
  return( list(allocation = allocation, repartitionFeuille = proportionFinale, stockMatrix=stockMatrix,
               covIDs = covIDsEvol, cutValues = cutValuesEvol, supOrInf = supOrInfEvol) )
  
}

### Function to predict thanks to the storage matrix

getPred <- function(resLocalTree, Npred){
  
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
      
      proportionFinale <- table(resLocalTree$stockMatrix[indToCompute,])
      res <- sample(names(proportionFinale[which(proportionFinale == max(proportionFinale))]),1)
    }
    
  }
  
  return(res)
}


### Function for the construction of forests
### It can return predictions for more than one minNodeSize value thanks to multiMinNodeSize

# y : the training response
# x : the training explanatory variables
# obs : the observed data
# mtry : the number of covariate to sample
# multiMinNodeSize : to return predictions for different minNodeSize values
# alpha : the quantile order
# ntree : the number of trees
# bootstrap : do we use bootstrap samples for the tree construction
# whichKernel : which type of kernel to use
# hfixe : does the weights are computed only once at the root (TRUE), or are they updated at each internal node (FALSE)?
# covWeights : to add covariate weights for their sampling

forestLocalMultiple <- function(y, x, obs, mtry=floor(sqrt(ncol(x))), multiMinNodeSize = 1, alpha,
                                ntree = 100, bootstrap = TRUE, whichKernel, hfixe, covWeights=NULL){
  q <- length(multiMinNodeSize)
  allocationTree <- matrix(NA, q , ntree) # tree votes for the observed data
  
  for(i in 1:ntree){
    tree.tmp <- treeLocal(y = y, x = x, obs = obs, mtry = mtry, minNodeSize = 1, alpha=alpha, bootstrap = bootstrap,
                          Nstock = max(multiMinNodeSize), hfixe = hfixe, whichKernel = whichKernel, covWeights=covWeights)
    l <- 0
    for(j in multiMinNodeSize){
      l <- l+1
      allocationTree[l,i] <- getPred(tree.tmp, j)
    }
  }
  row.names(allocationTree) <- multiMinNodeSize
  prediction <- sapply(c(1:q), function(X) sample(names(table(allocationTree[X,])[which(table(allocationTree[X,])==max(table(allocationTree[X,])))]),1) )
  names(prediction) <- multiMinNodeSize
  
  return( list(prediction = prediction, allocationTree = allocationTree) )
  
}