// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "hal9001_types.h"
using namespace Rcpp;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Functions to enumerate basis functions
//------------------------------------------------------------------------------

// populates a map with unique basis functions based on data in xsub values are
// thresholds, keys are column indices
BasisMap enumerate_basis(const NumericMatrix& X_sub,
                         const NumericVector& cols) {
  BasisMap bmap;

  //find unique basis functions
  int n = X_sub.rows();
  for(int i = 0; i < n; i++) {
    NumericVector cutoffs = X_sub.row(i);
    bmap.insert(std::pair<NumericVector, NumericVector>(cutoffs, cols));
  }
  //erase the lowest(always true) basis function
  //actually, I think we don't want to do this
  //because it might not be true in an OOB prediction set
  //bmap.erase(bmap.begin());
  return(bmap);
}

//------------------------------------------------------------------------------

//' Sort Basis Functions
//'
//' Build a sorted list of unique basis functions based on columns, where each
//' basis function is a list
//'
//' @details Note that sorting of columns is performed such that the basis order
//' equals cols.length() and each basis function is a list(cols, cutoffs).
//'
//' @param X_sub A subset of the columns of X, the original design matrix.
//' @param cols An index of the columns that were reduced to by sub-setting.
// [[Rcpp::export]]
List make_basis_list(const NumericMatrix& X_sub, const NumericVector& cols, const IntegerVector& order_map){

  BasisMap bmap = enumerate_basis(X_sub, cols);
  List basis_list(bmap.size());
  int index = 0;
  for (BasisMap::iterator it = bmap.begin(); it != bmap.end(); ++it) {
    // List basis(2);
    // basis[0]=it->second;
    // basis[1]=it->first;

    // List basis();
    // basis["cols"]=it->second;
    // basis["cutoffs"]=it->first;
    NumericVector subCols = it->second;

    IntegerVector order (subCols.length());
    for(int i=0; i < subCols.length(); i++){
      order[i] = order_map[subCols[i]-1];
    }

    List basis = List::create(
      Rcpp::Named("cols") = it->second,
      Rcpp::Named("cutoffs") = it->first,
      Rcpp::Named("orders") = order
    );
    basis_list[index++] = basis;
  }
  return(basis_list);
}

//------------------------------------------------------------------------------

//' Compute Values of Basis Functions
//'
//' Computes and returns the indicator value for the basis described by
//' cols and cutoffs for a given row of X (X[row_num, ])
//'
//' @param X The design matrix, containing the original data.
//' @param row_num Numeri for  a row index over which to evaluate.
//' @param cols Numeric for the column indices of the basis function.
//' @param cutoffs Numeric providing thresholds.
//'
// [[Rcpp::export]]
double meets_basis(const NumericMatrix& X, const int row_num,
                   const IntegerVector& cols, const NumericVector& cutoffs,  const IntegerVector& orders) {
  int p = cols.length();
  double value = 1;


  for (int i = 0; i<p; i++) {
    double obs = X(row_num,cols[i] - 1); // using 1-indexing for basis columns
    int order =  orders[i];
    double cutoff = cutoffs[i];
    if(!(obs > cutoff)) {
      return(0);
    }
    if(order!=0){
      value = value * pow((obs - cutoff),order);
    }

  }

  return(value);
}



//------------------------------------------------------------------------------

//' Generate Basis Functions
//'
//' Populates a column (indexed by basis_col) of x_basis with basis indicators.
//'
//' @param basis The basis function.
//' @param X The design matrix, containing the original data.
//' @param x_basis The HAL design matrix, containing indicator functions.
//' @param basis_col Numeric indicating which column to populate.
//'
// [[Rcpp::export]]
void evaluate_basis(const List& basis, const NumericMatrix& X, SpMat& x_basis,
                    int basis_col) {
  int n = X.rows();


  //split basis into x[1] x[-1]
  //find sub-bases
  //intersect
  IntegerVector cols = as<IntegerVector>(basis["cols"]);
  NumericVector cutoffs = as<NumericVector>(basis["cutoffs"]);
  IntegerVector orders =  as<IntegerVector>(basis["orders"]);
  for (int row_num = 0; row_num < n; row_num++) {
    double value = meets_basis(X, row_num, cols, cutoffs,  orders);
    if (value!=0) {
      //Add value
      x_basis.insert(row_num, basis_col) = value;
    }



  }
}

//------------------------------------------------------------------------------

//' Build HAL Design Matrix
//'
//' Make a HAL design matrix based on original design matrix X and a list of
//' basis functions in argument blist
//'
//' @param X Matrix of covariates containing observed data in the columns.
//' @param blist List of basis functions with which to build HAL design matrix.
//'
//' @export
//'
//' @examples
//' \donttest{
//' gendata <- function(n) {
//'   W1 <- runif(n, -3, 3)
//'   W2 <- rnorm(n)
//'   W3 <- runif(n)
//'   W4 <- rnorm(n)
//'   g0 <- plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4))
//'   A <- rbinom(n, 1, g0)
//'   Q0 <- plogis(0.15 * (2 * A + 2 * A * W1 + 6 * A * W3 * W4 - 3))
//'   Y <- rbinom(n, 1, Q0)
//'   data.frame(A, W1, W2, W3, W4, Y)
//' }
//' set.seed(1234)
//' data <- gendata(100)
//' covars <- setdiff(names(data), "Y")
//' X <- as.matrix(data[, covars, drop = FALSE])
//' basis_list <- enumerate_basis(X)
//' x_basis <- make_design_matrix(X, basis_list)
//' }
//'
//' @return A \code{dgCMatrix} sparse matrix of indicator basis functions
//'  corresponding to the design matrix in a zero-order highly adaptive lasso.
// [[Rcpp::export]]
SpMat make_design_matrix(const NumericMatrix& X, const List& blist) {
  //now generate an indicator vector for each
  int n = X.rows();
  int basis_p = blist.size();

  SpMat x_basis(n, basis_p);
  x_basis.reserve(0.5 * n * basis_p);

  List basis;
  NumericVector cutoffs, current_row;
  IntegerVector last_cols, cols;
  NumericMatrix X_sub;

  //for each basis function
  for (int basis_col = 0; basis_col < basis_p; basis_col++) {
    last_cols = cols;
    basis = (List) blist[basis_col];
    evaluate_basis(basis, X, x_basis, basis_col);
  }
  x_basis.makeCompressed();
  return(x_basis);
}

