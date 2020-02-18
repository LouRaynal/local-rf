#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix findweights(NumericMatrix predTraining, NumericMatrix predTesting, int nTrain, int nTest, int ntree){
  NumericMatrix result(nTrain, nTest);
  IntegerVector counti(nTrain);
  for(int k=0; k<ntree; k++){
    for(int i=0; i<nTest; i++){
      for(int j=0; j<nTrain; j++){
        if( predTraining(j,k) == predTesting(i,k) ){
          counti[j] = 1;
        } else{
          counti[j] = 0;
        }
      }
      for(int j=0; j<nTrain; j++){
        result(j,i) = result(j,i) + counti[j];
      }
    } 
  }
  return(result);
}