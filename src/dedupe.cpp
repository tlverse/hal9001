// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "hal9001_types.h"
using namespace Rcpp;

//' Find Copies of Columns
//'
//' Index vector that, for each column in X, indicates the index of the first
//' copy of that column
//'
//' @param X Sparse matrix containing columns of indicator functions.
//'
// [[Rcpp::export]]
IntegerVector index_first_copy(const MSpMat& X) {
  int p = X.cols();

  ColMap col_map;
  IntegerVector copy_index(p);

  for (int j = 0; j < p; j++) {
    MSpMatCol current_col(X, j);

    //https://stackoverflow.com/questions/97050/stdmap-insert-or-stdmap-find
    ColMap::iterator match = col_map.lower_bound(current_col);
    if (match != col_map.end() &&
        !(col_map.key_comp()(current_col, match->first))) {
      // column already exists
      copy_index[j] = match->second + 1; //use 1-indexing
    } else {
      // column not yet in map
      col_map.insert(match, ColMap::value_type(current_col, j));
      copy_index[j] = j+1; //use 1-indexing
    }
  }
  return(copy_index);
}

//------------------------------------------------------------------------------

//' Apply copy map
//'
//' OR duplicate training set columns together
//' @param X Sparse matrix containing columns of indicator functions.
//' @param copy_map the copy map
// [[Rcpp::export]]
SpMat apply_copy_map(const MSpMat& X, const List& copy_map) {
  int n = X.rows();
  int basis_p = copy_map.size();

  SpMat x_unique(n, basis_p);
  x_unique.reserve(0.5 * n * basis_p);

  for(int j=0; j<basis_p; ++j){
    IntegerVector copy_group=copy_map[j];

    int p_cols = copy_group.size();
    for (int i = 0; i < n; i++) {

      for (int j_original = 0; j_original < p_cols; j_original++) {
        int col = copy_group[j_original] - 1;  //cols is 1-indexed
        if (X.coeff(i, col) == 1) {
          x_unique.insert(i, j) = 1;
          break;
        }
      }
    }
  }
  x_unique.makeCompressed();
  return(x_unique);
}
