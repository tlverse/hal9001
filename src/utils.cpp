// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include "hal9001_types.h"
using namespace Rcpp;
//------------------------------------------------------------------------------

//' Fast Coercion to Sparse Matrix
//'
//' Fast and efficient coercion of standard matrix objects to sparse matrices.
//' Borrowed from http://gallery.rcpp.org/articles/sparse-matrix-coercion/.
//' INTERNAL USE ONLY.
//'
//' @param XX_ An object of class \code{Matrix} that has a sparse structure
//'  suitable for coercion to a sparse matrix format of \code{dgCMatrix}.
//'
//' @return An object of class \code{dgCMatrix}, coerced from input \code{XX_}.
//'
// [[Rcpp::export]]
SEXP as_dgCMatrix( SEXP XX_ ) {
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::Map<Eigen::MatrixXd> MapMatd; // Input: must be double
  MapMatd X(Rcpp::as<MapMatd>(XX_));
  SpMat Xsparse = X.sparseView();              // Output: sparse matrix
  S4 Xout(wrap(Xsparse));                      // Output: as S4 object
  NumericMatrix Xin(XX_);                      // Copy dimnames
  Xout.slot("Dimnames") = clone(List(Xin.attr("dimnames")));
  return(Xout);
}

//------------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector non_zeros(const MSpMat& X) {
  int p = X.cols();
  int j;
  int nz;

  IntegerVector non_zeros(p);

  for (j = 0; j < p; ++j) {
    nz = 0;
    for (MInIterMat i_(X, j); i_; ++i_) {
      nz++;
    }
    non_zeros[j] = nz;
  }
  return(non_zeros);
}

//------------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector calc_pnz(const MSpMat& X) {
  IntegerVector nz = non_zeros(X);
  int n = X.rows();
  NumericVector pnz = as<NumericVector>(nz)/n;

  return(pnz);
}

// [[Rcpp::export]]
NumericVector get_pnz(const MSpMat& X) {
  IntegerVector nz = non_zeros(X);
  int n = X.rows();
  NumericVector pnz = as<NumericVector>(nz)/n;

  return(pnz);
}
//------------------------------------------------------------------------------
NumericVector not_dumb_sqrt(const NumericVector& x){
  NumericVector res(x.length());
  for(int i=0; i<x.length(); i++){
    res[i] = std::sqrt(x[i]);
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector calc_xscale(const MSpMat& X, const NumericVector& xcenter) {
 int n = X.rows();
 NumericVector pnz = get_pnz(X);
 NumericVector xscale = not_dumb_sqrt(pnz);
 double minx = std::sqrt(1.0 / n);

 xscale = not_dumb_sqrt(xscale * xscale - (xcenter * xcenter));
 xscale[xscale == 0 ] = minx;

 return(xscale);
}

// [[Rcpp::export]]
NumericVector get_xscale(const MSpMat& X, const NumericVector& xcenter) {
  int n = X.rows();
  NumericVector pnz = get_pnz(X);
  NumericVector xscale = not_dumb_sqrt(pnz);
  double minx = std::sqrt(1.0 / n);

  xscale = not_dumb_sqrt(xscale * xscale - (xcenter * xcenter));
  xscale[xscale == 0 ] = minx;

  return(xscale);
}


// [[Rcpp::export]]
bool equal_double(double x, double y){
  return(std::abs(x - y) < 1e-16);
}

double soft_max(double beta, double lambda){
  if (beta > lambda) {
    beta -= lambda;
  } else if (beta < -1 * lambda) {
    beta += lambda;
  } else {
    beta = 0;
  }
  return(beta);
}
