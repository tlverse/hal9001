
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
typedef Rcpp::NumericVector MatRow;

struct cmpMatrixRow {
  bool operator()(const MatRow& a, const MatRow& b) const {
    int i=0;
    
    int smaller_length=a.size();
    if(b.size()<smaller_length){
      smaller_length=b.size();
    }
    
    for(i=0; i<smaller_length; i++){
      if(a[i]==b[i]){
        //skip anything at the beginning that matches
        continue;
      } else {
        //once there's a mismatch, determine which one is bigger
        return(a[i]<b[i]);
      }
    }
    

    
    return(a.size() < b.size());
  }
};

typedef std::map<MatRow, int, cmpMatrixRow> BasisMap;