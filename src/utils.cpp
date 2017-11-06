// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

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
SEXP asdgCMatrix_( SEXP XX_ ) {
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::Map<Eigen::MatrixXd> MapMatd; // Input: must be double
  MapMatd X(Rcpp::as<MapMatd>(XX_));
  SpMat Xsparse = X.sparseView();              // Output: sparse matrix
  S4 Xout(wrap(Xsparse));                      // Output: as S4 object
  NumericMatrix Xin(XX_);                      // Copy dimnames
  Xout.slot("Dimnames") = clone(List(Xin.attr("dimnames")));
  return(Xout);
}

