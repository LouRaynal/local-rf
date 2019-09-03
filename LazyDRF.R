# Required packages

library(discretization) # For data discretization


### Functions to train some lazy decision RF

# x : explanatory variables
# y : response
# obs : the observed data (explanatory variables)
# Nmin : minimal number of observations in a leaf
# pGain : percentage of the maximal achievable gain above which options are considered
# mtry : covariate sampling parameter

# We assume x is a matrix, y and obs are vectors, 0<pGain<=1, Nmin and mtry are numerical values

LazyOptionTree <- function(x, y, obs, Nmin, pGain, mtry){
  
  if(mtry > length(obs) || mtry<1){
    stop("mtry is either too large or to small")
  }
  
  dataDisc <- discret(y = y, x = x, obs = obs)
  
  result <- explore(dataDisc$x, y, dataDisc$obs, Nmin, pGain, mtry)
  
  return(result)
  
}

### Function to discretize x and obs

discret <- function(y, x, obs){
  
  data <- data.frame(x, mod = as.factor(y))
  res.mdlp <- mdlp(data)
  
  x.disc <- res.mdlp$Disc.data[,-ncol(res.mdlp$Disc.data), drop=FALSE]
  
  obs.disc <- matrix(NA, nrow=1, ncol=ncol(x.disc))
  
  for( i in 1:ncol(x.disc) ){
    if(any(res.mdlp$cutp[[i]]=="All")){
      obs.disc[,i] <- 1
    } else{
      obs.disc[,i] <- cut(obs[i], breaks = c(-Inf,res.mdlp$cutp[[i]],Inf), labels = FALSE)
    }
  }
  
  return(list(x = x.disc, obs = obs.disc))
  
}

### Develop the lazy tree

explore <- function(x, y, obs, Nmin, pGain, mtry){
  
  maxGain <- 0  # Initialize the gain
  done <- FALSE # Is this branch maximaly developed?
  maxSupport <- 0 # Record the number of individuals used to predict
  
  # Preliminary checks (stopping rules)
  
  # Unique class?
  if(equalVector(y)){
    result <- list(pred = y[1], proportion = table(y), nSupport = length(y))
    done <- TRUE
  }
  
  # Identical covariates?
  if(equalData(x)){
    result <- majority(y)
    done <- TRUE
  }
  
  # Too few data at the leaf?
  if(nrow(x) < Nmin){
    result <- majority(y)
    done <- TRUE
  }
  
  if(!done){
    
    nbrCovariates <- length(obs)
    covToTry <- sample(1:nbrCovariates, mtry, replace = FALSE)
    
    # Compute the mother entropy
    
    weights <- classWeights(y)
    entMother <- entropy(x, y, covNumber = NULL, cut = NULL, weights = weights)
    
    # What are the possible cuts?
    
    listOfCuts <- computeCutsToTry(x, obs)
    
    # Store the different gains in a list
    
    listOfGains <- listOfCuts
    
    # Compute the entropy for every possible covariate and cut value

    for(k in covToTry){
      
      if( !is.null(listOfCuts[[k]]) ){
        
        for(j in 1:length(listOfCuts[[k]])){
          
          entChild <- entropy(x, y, covNumber = k, cut = listOfCuts[[k]][j], weights = weights)
          newGain <- entMother - entChild
          
          # Store the gain
          listOfGains[[k]][j] <- newGain
          
          if(newGain > maxGain){
            maxGain <- newGain
          }
          
        }
        
      }
      
    }
  
    optionsList <- listOfGains
    
    # Option selection
    
    for(k in covToTry){
      
      if( !is.null(listOfGains[[k]]) ){
        optionsList[[k]] <- (listOfGains[[k]] >= (pGain*maxGain)) & (listOfGains[[k]]!=0)
      }
      
    }
    
    for(k in covToTry){
      
      if( !is.null(optionsList[[k]]) ){
        
        for(j in 1:length(optionsList[[k]])){
          
          if(optionsList[[k]][j]){
            
            # For each valid option, we develop the corresponding node recursively using explore()
            partial <- tryNewExplore(x, y, obs, covNumber = k, cut = listOfCuts[[k]][j], Nmin, pGain, mtry)
            if(partial$nSupport > maxSupport){
              maxSupport <- partial$nSupport
              result <- partial
            }
            
          }
        }
      }
    }
    
    if(maxGain == 0){
      result <- majority(y)
    }
    
  }
  
  return(result)
  
}

### Check if only one modality in y

equalVector <- function(y){
  
  if(all(y[1]==y)){
    return(TRUE)
  } else{
    return(FALSE)
  }
  
}

### Check if all covariates lines are equal

equalData <- function(x){
  
  if(nrow(unique(x))==1){
    return(TRUE)
  } else{
    return(FALSE)
  }
  
}

### Compute the majority class

majority <- function(y){
  
  proportion <- table(y)
  # If equality we take the majority at random among the equal classes
  pred <- sample(names(proportion[which(proportion == max(proportion))]),1)
  nSupport <- sum(y==pred)
  return( list(pred = pred, proportion = proportion, nSupport = nSupport) )
  
}

### Compute cuts to try
### Reminder, a cut value is the obs modality

computeCutsToTry <- function(x, obs){
  
  listOfCuts <- vector("list", length(obs))
  
  for(k in 1:length(obs)){
    
    whichIsEqual <- obs[k] == x[,k]
    isXStarMod <- sum(whichIsEqual)!=0
    if(isXStarMod){
      listOfCuts[[k]] <- obs[k]
    }
    
  }
  
  return(listOfCuts)
  
}

### Compute the class weights. An instance in class k has weights 1/(n_k*K) if K classes

classWeights <- function(y){
  
  y <- factor(y)
  weights <- rep(NA, length(y))
  nClasse <- nlevels(y)
  lvl <- levels(y)
  
  # On donne les poids aux données, par classe
  for(k in 1:nClasse){
    idx <- y==lvl[k]
    nLvl <- sum(idx)
    weights[idx] <- 1/(nLvl*nClasse)
  }
  
  return(weights)
  
}

### Compute the weighted entropy

entropy <- function(x, y, covNumber, cut, weights){
  
  if(!is.null(covNumber) && !is.null(cut)){
    onWhichToCompute <- x[,covNumber] == cut
    y <- factor(y[onWhichToCompute])
    weights <- weights[onWhichToCompute]
  }
  
  y <- factor(y)
  nClasse <- length(unique(y))
  n_k.w_k <- rep(NA, nClasse)
  p_k <- rep(NA, nClasse)
  lvl <- levels(y)
  
  for(k in 1:nClasse){
    n_k.w_k[k] <- sum(weights[y==lvl[k]])
  }
  
  p_k <- n_k.w_k/sum(n_k.w_k)
  res <- 0
  for(k in 1:nClasse){
    if(p_k[k]!=0){
      res <- res + p_k[k]*log(p_k[k])
    }
  }
  res <- -res
  
  return(res)
  
}

### Explore a selected option

tryNewExplore <- function(x, y, obs, covNumber, cut, Nmin, pGain,  mtry){
  
  whichIsToKeep <- x[,covNumber]==cut
  
  newx <- x[whichIsToKeep,,drop=FALSE]
  newy <- y[whichIsToKeep]
  
  result <- explore(newx, newy, obs, Nmin, pGain, mtry)
  
  return(result)
  
}

### Build a lazyDRF

# bootstrap : should each tree be built on bootstrap samples?
# num.trees : number of trees

LazyOptionForest <- function(x, y, obs, Nmin, pGain,  mtry=floor(sqrt(length(obs))), bootstrap, num.trees){
  
  predictionTrees <- rep(NA, num.trees)
  nSupportPerTree <- rep(NA, num.trees) # personal
  
  for(i in 1:num.trees){
    
    if(bootstrap){
      toBoot <- sample(x = 1:length(y), size = length(y), replace = TRUE)
    } else {
      toBoot <- 1:length(y)
    }
    
    resLazyDT <- LazyOptionTree(x[toBoot,,drop=FALSE], y[toBoot], obs, Nmin, pGain,  mtry)
    predictionTrees[i] <- resLazyDT$pred
    nSupportPerTree[i] <- resLazyDT$nSupport
    
  }
  
  result <- majority(predictionTrees)
  
  return(result)
  
}
