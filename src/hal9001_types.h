// [[Rcpp::depends(RcppEigen)]]
#ifndef HAL9001_TYPES_H
#define HAL9001_TYPES_H
#include <RcppEigen.h>
using namespace Rcpp;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<int> IntSpMat;
typedef Eigen::Map<SpMat> MSpMat;
typedef MSpMat::InnerIterator MInIterMat;

typedef SpMat::InnerIterator InIterMat;
typedef SpMat::InnerVectorReturnType InVec;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

struct cmpMatrixRow {
  bool operator()(const NumericVector& a, const NumericVector& b) const {

    int i = 0;

    int smaller_length = a.size();
    if (b.size()<smaller_length) {
      smaller_length = b.size();
    }

    for (i = 0; i<smaller_length; i++) {
      if(a[i] == b[i]) {
        //skip anything at the beginning that matches
        continue;
      } else {
        //once there's a mismatch, determine which one is bigger
        return(a[i] < b[i]);
      }
    }
    return(a.size() < b.size());
  }
};

typedef std::map<NumericVector, NumericVector, cmpMatrixRow> BasisMap;

typedef  std::pair<const MSpMat&,int> MSpMatCol;

struct cmpCol {
  bool operator()(const MSpMatCol& a, const MSpMatCol& b) const {
    //returns true if a is strictly less than b
    const MSpMat& X = a.first;
    int col_a = a.second;
    int col_b = b.second;

    MInIterMat iter_b(X, col_b);
    for (MInIterMat iter_a(X, col_a); iter_a; ++iter_a,++iter_b) {
      if (!iter_b) {
        //we've matched the entirety of b to a, but there's still more...
        //...elements in a, so it comes after
        //iter_b is shorter
        return(false);
      }
      int index_a = iter_a.index();
      int index_b = iter_b.index();

      // Rcout << index_a << " " << index_b << std::endl;
      if (index_a == index_b) {
        //skip anything at the beginning that matches
        continue;
      } else {
        //once there's a mismatch, determine sort order
        //if iter_a has a lower index, it comes later in the sort order
        return(index_a > index_b);
      }
    }
    //we've matched the entirety of a to b
    //if there are more elements in b, it comes after, otherwise they're a match
    return(iter_b);
  }
};

typedef std::map<MSpMatCol, int, cmpCol> ColMap;

#endif //HAL9001_TYPES_H
