# Required package
library(ranger)

### Function to determine local importance of covariates

LocalVarImp <- function(ranger.rf.object, obs){
  
  nodeIDs <- ranger.rf.object$forest$child.nodeIDs
  varIDs <- ranger.rf.object$forest$split.varIDs
  
  predNodeID <- predict(object = ranger.rf.object, data = obs, type = "terminalNodes")$predictions
  
  # We count the number of times each covariate has been used to allocate the observed data
  
  impx <- rep(0, ranger.rf.object$num.independent.variables)
  ntree <- rf.ranger$forest$num.trees
  for( b in 1:ntree ){
    
    motherNode <- predNodeID[b]
    
    while(motherNode != 0){
      
      if(motherNode%%2 == 1){
        motherNode <- which(nodeIDs[[b]][[1]] == motherNode)-1 # Odd
      }
      else{
        motherNode <- which(nodeIDs[[b]][[2]] == motherNode)-1 # Even
      }
      
      impx[varIDs[[b]][motherNode+1]] <- impx[varIDs[[b]][motherNode+1]] + 1 # Increment
      
    }
    
  }
  
  impxStand <- impx/sum(impx) 
  
  return(stdImp = impxStand)
}