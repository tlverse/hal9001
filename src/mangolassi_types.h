
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::MappedSparseMatrix<bool> MSpMat_bool;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<bool> SpMat_bool;
typedef MSpMat::InnerIterator InIterMat;
typedef MSpMat::InnerVectorReturnType InVec;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;
