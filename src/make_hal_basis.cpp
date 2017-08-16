
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mangolassi_types.h"
using namespace Rcpp;

// [[Rcpp::export]]
SpMat make_hal_basis(NumericMatrix x){
  SpMat x_basis;
  
  return(x_basis);
}

//quick and dirty way of validating comparator
// [[Rcpp::export]]
bool compare_vectors(MatRow v1, MatRow v2){
  BasisMap bmap;
  bmap[v1]=0;
  return(bmap.key_comp()(v1,v2));
}


// generates a basismap with keys the values to test against
// [[Rcpp::export]]
SpMat make_basis(NumericMatrix X, IntegerVector cols){
  int n=X.rows();
  int p=X.cols();
  int sub_p=cols.length();
  BasisMap bmap;
  
  //construct submatrix with only basis cols
  List all_cutpoints(p);
  NumericMatrix X_sub(n,sub_p);
  for(int j=0; j<sub_p; j++){
    X_sub.column(j)=X.column(cols[j]);
  }
  
  //find unique basis functions
  for(int i=0; i<n; i++){
    MatRow row=X_sub.row(i);
    bmap.insert(std::pair<MatRow,int>(row, -1));
  }
  //erase the lowest(always true) basis function
  bmap.erase(bmap.begin());
  
  //now generate an indicator vector for each
  int basis_p=bmap.size();
  SpMat x_basis(n,basis_p);
  
  //begin by assuming every observation is bigger than the smallest basis function
  SpVec indicators=Eigen::VectorXd::Ones(n).sparseView();
  
  //for each basis function
  int offset=0;
  MatRow basis;
  for (BasisMap::iterator it=bmap.begin(); it!=bmap.end(); ++it){
    basis=it->first;
    //verify that previously true indicators remain true
    for (InIterVec ind_it(indicators); ind_it; ++ind_it){
      MatRow current_row=X_sub.row(ind_it.index());
      if(!bmap.key_comp()(current_row,basis)){
        //we can add a positive indicator for this row, basis
        x_basis.coeffRef(ind_it.index(),offset)=1;
      } else{
        //otherwise we should stop checking this row going forward
        indicators.coeffRef(ind_it.index())=0;
      }
    }
    offset++;
  }
  
  x_basis.makeCompressed();
  return(x_basis);
}

