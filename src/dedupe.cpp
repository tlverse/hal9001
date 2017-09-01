// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "hal9001_types.h"
using namespace Rcpp;

//' Find Copies of Columns
//'
//' Index vector that, for each column in X, indicates the index of the first
//' copy of that column
//'
// [[Rcpp::export]]
IntegerVector index_first_copy(const MSpMat& X){
  int p = X.cols();

  ColMap col_map;
  IntegerVector copy_index(p);

  for(int j = 0; j < p; j++) {
    MSpMatCol current_col(X, j);

    //https://stackoverflow.com/questions/97050/stdmap-insert-or-stdmap-find
    ColMap::iterator match = col_map.lower_bound(current_col);
    if(match != col_map.end() && !(col_map.key_comp()(current_col, match->first)))
    {
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

//' Compute Column Ranks
//'
//' \code{TRUE} iff col_1 is strictly less than col_2 in the ordering scheme
//'
// [[Rcpp::export]]
bool column_compare(const MSpMat& X, int col_1, int col_2){
  ColMap cmap;
  MSpMatCol X_1(X, col_1);
  MSpMatCol X_2(X, col_2);
  return(cmap.key_comp()(X_1, X_2));
}

//------------------------------------------------------------------------------

//' Reorder Columns of X
//'
//' ORs the columns of X listed in cols and places the result in column col[1]
//'
// [[Rcpp::export]]
void or_duplicate_columns(MSpMat& X, const IntegerVector& cols){
  int first = cols[0] - 1; //cols is 1-indexed
  int p_cols = cols.length();
  int n = X.rows();
  for(int i = 0; i < n; i++){
    if(X.coeffRef(i, first) == 1) {
      //this is already 1
      continue;
    }
    //search remaining columns for 1, inserting into first if found
    for(int j = 1; j < p_cols; j++) {
      int j_col = cols[j] - 1;  //cols is 1-indexed
      if(X.coeffRef(i, j_col) == 1) {
        X.coeffRef(i, j_col) = 1;
        break;
      }
    }
  }
}

