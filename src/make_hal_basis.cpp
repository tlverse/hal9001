
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mangolassi_types.h"

using namespace Rcpp;

// ------------------------------------------------------------------------
// Functions to enumerate basis functions

// populates a map with unique basis functions based on data in xsub
// values are thresholds, keys are column indicies
BasisMap enumerate_basis(const NumericMatrix& X_sub, const NumericVector& cols){
  BasisMap bmap;


  //find unique basis functions
  int n=X_sub.rows();
  for(int i=0; i<n; i++){
    NumericVector cutoffs=X_sub.row(i);
    bmap.insert(std::pair<NumericVector,NumericVector>(cutoffs, cols));
  }

  //erase the lowest(always true) basis function
  //actually, I think we don't want to do this
  //because it might not be true in an OOB prediction set
  //bmap.erase(bmap.begin());


  return(bmap);
}



// returns a sorted list of unique basis functions based on columns in cols (so basis order=cols.length())
// each basis function is a list(cols,cutoffs)
// X_sub is a subset of the columns of X (the original design matrix)
// cols is an index of the columns that were subsetted
// [[Rcpp::export]]
List make_basis_list(const NumericMatrix& X_sub, const NumericVector& cols){
  
  BasisMap bmap=enumerate_basis(X_sub, cols);
  List basis_list(bmap.size());
  int index=0;
  for (BasisMap::iterator it=bmap.begin(); it!=bmap.end(); ++it){
    // List basis(2);
    // basis[0]=it->second;
    // basis[1]=it->first;
    
    // List basis();
    // basis["cols"]=it->second;
    // basis["cutoffs"]=it->first;
    
    List basis=List::create(
      Rcpp::Named("cols") = it->second,
      Rcpp::Named("cutoffs") = it->first
    );
      
    basis_list[index++]=basis;
  }
  
  return(basis_list);
}

// ------------------------------------------------------------------------
// Functions to make a design matrix based on a list of basis functions

// returns the indicator value for the basis described by cols,cutoffs for X[row_num,] 
// X is the original design matrix
// row_num is a row index to evaluate
// cols are the column incides of the basis function
// cutoffs are thresholds
// [[Rcpp::export]]
bool meets_basis(const NumericMatrix& X, const int row_num, const IntegerVector& cols, const NumericVector& cutoffs){
  int p=cols.length();
  for(int i=0; i<p; i++){
    double obs=X(row_num,cols[i]-1); //we're using 1-indexing for basis columns
    if(!(obs>=cutoffs[i])){
      return(false);
    }
  }
  
  return(true);
}


// populates a column (indexed by basis_col) of x_basis with basis indicators
// basis is the basis function
// X is the original design matrix
// x_basis is the hal design matrix
// basis_col indicates which column to populate
// [[Rcpp::export]]
void evaluate_basis(const List& basis, const NumericMatrix& X, SpMat& x_basis, int basis_col){
  int n=X.rows();
  //split basis into x[1] x[-1]
  //find sub-basises
  //intersect
  
  IntegerVector cols=as<IntegerVector>(basis["cols"]);
  NumericVector cutoffs=as<NumericVector>(basis["cutoffs"]);
  for (int row_num=0; row_num<n; row_num++){
    
    if(meets_basis(X,row_num,cols,cutoffs)){
      //we can add a positive indicator for this row, basis
      x_basis.insert(row_num,basis_col)=1;
    }    
  }
  
}

// makes a hal design matrix based on original design matrix X and
// a list of basis functions in blist
// [[Rcpp::export]]
SpMat make_design_matrix(NumericMatrix X, List blist){
  //now generate an indicator vector for each
  int n=X.rows();
  int basis_p=blist.size();


  SpMat x_basis(n,basis_p);  
  x_basis.reserve(0.5*n*basis_p);


  List basis;
  NumericVector cutoffs, current_row; 
  IntegerVector last_cols, cols;
  NumericMatrix X_sub;
  
  //for each basis function  
  for (int basis_col=0; basis_col<basis_p; basis_col++){
    last_cols=cols;
    
    basis=(List) blist[basis_col];

    evaluate_basis(basis, X, x_basis, basis_col);

  }

  x_basis.makeCompressed();
  return(x_basis);
}
