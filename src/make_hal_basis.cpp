
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mangolassi_types.h"
using namespace Rcpp;

// [[Rcpp::export]]
SpMat make_hal_basis(NumericMatrix x){
  SpMat x_basis;
  
  return(x_basis);
}

NumericVector get_cutpoints(NumericMatrix::Column col){
  NumericVector unique_vals=unique(col);
  
  // drop first element
  unique_vals.sort();
  unique_vals.erase(0);
  
  return(unique_vals);
}


void eval_cutpoints(NumericMatrix::Column x_col, NumericVector cutpoints, int offset, SpMat & x_basis){
  // initialize all indicators to be 1s
  int n=x_col.size();
  
  //begin by assuming every observation is bigger than the cutpoint
  SpVec indicators=Eigen::VectorXd::Ones(n).sparseView();
  
  //for each cutpoint in the list
  for(NumericVector::iterator cut_it = cutpoints.begin(); cut_it != cutpoints.end(); ++cut_it) {
    //verify that previously true indicators remain true
    double cutpoint=*cut_it;

    for (InIterVec ind_it(indicators); ind_it; ++ind_it){
      //if we exceed the cutpoint
      if(x_col[ind_it.index()]>=cutpoint){
        //we can add a positive indicator for this row, basis
        x_basis.coeffRef(ind_it.index(),offset)=1;
      } else{
        //otherwise we should stop checking this row going forward
        indicators.coeffRef(ind_it.index())=0;
      }
    }
    
    
    offset++;
  }
  
}
// [[Rcpp::export]]
SpMat make_univariate_basis(NumericMatrix X){
  int n=X.rows();
  int p=X.cols();
  int basis_p=0;
  List all_cutpoints(p);
  
  for(int base_col=0; base_col<p; base_col++){
    NumericVector cutpoints=get_cutpoints(X.column(base_col));
    basis_p+=cutpoints.length();
    all_cutpoints[base_col]=cutpoints;
  }
  
  
  SpMat x_basis(n,basis_p);
  int offset=0;
  for(int base_col=0; base_col<p; base_col++){
    NumericMatrix::Column x_col=X.column(base_col);
    NumericVector cutpoints=all_cutpoints[base_col];
    eval_cutpoints(x_col, cutpoints, offset, x_basis);
    offset+=cutpoints.length();
  }
  
  x_basis.prune(0.5,0.1);
  return(x_basis);
}

//define basis function as set of indicators
//write function to cross two basis functions
//can use that via induction I think