
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mangolassi_types.h"
using namespace Rcpp;

////------------------------------------------------------
// General code

bool columns_equal(MSpMat X, int col_1, int col_2){
  
      InIterMat iter_2(X, col_1);
      for (InIterMat iter_1(X, col_2); iter_1; ++iter_1){
         if(!iter_2){
           //iter_2 is shorter
           return(false);

         }
         if(iter_1.index()!=iter_2.index()){
          //index mismatch
          // Rcout << "index mismatch " << std::endl;
          return(false);
         }
         
         ++iter_2;
      }
      
      if(iter_2){
        //iter_2 is longer
        return(false);
      }
      
      return(true);
}
// [[Rcpp::export]]
LogicalVector dedupe(MSpMat X){
  int p=X.cols();
  int k,j;
  LogicalVector unique_cols(p, true);
  
  for (k=0; k<X.outerSize(); ++k){
    if(!unique_cols[k]){
      // we've already determined this column to be a duplicate
      // therefore, we've already checked its original against all other columns
      // so we can move on
      continue;
    }
    
    for (j=k+1; j<X.outerSize(); ++j){
      
      if(!unique_cols[j]){
        // we've already determined this column to be a duplicate
        // so we can move on
        continue;
      }
      
      unique_cols[j]=!columns_equal(X, k, j);
    }
    
  }
  
  return(unique_cols);
}

// [[Rcpp::export]]
LogicalVector dedupe_math(MSpMat X){
  int p=X.cols();
  int k,j;
  LogicalVector unique_cols(p, true);
  SpVec col_1;
  SpVec col_2;
  bool is_dupe=false;
  
  for (k=0; k<X.outerSize(); ++k){
    if(!unique_cols[k]){
      // we've already determined this column to be a duplicate
      // therefore, we've already checked its original against all other columns
      // so we can move on
      continue;
    }
    
    for (j=k+1; j<X.outerSize(); ++j){
      // start out assuming duplicate
      // Rcout << "Comparing " << j << " and " << k << " ";
      is_dupe=true;
      
      col_1=X.innerVector(k);
      col_2=X.innerVector(j);
      // can't be a dupe if the lengths are different
      if(col_1.nonZeros()!=col_2.nonZeros()){
        // Rcout << "nonzero length mismatch " << std::endl;
        is_dupe=false;
        continue;
      }
      
      InIterVec iter_2(col_2);
      // compare all indices (we know the values are equal)
      for (InIterVec iter_1(col_1); iter_1; ++iter_1){
        if(iter_1.index()!=iter_2.index()){
          //index mismatch
          // Rcout << "index mismatch " << std::endl;
          is_dupe=false;
          break;
        }
        
        ++iter_2;
      }
      // Rcout << is_dupe << std::endl;
      unique_cols[j]=unique_cols[j]&!(is_dupe);
    }
    
  }
  
  return(unique_cols);
}
